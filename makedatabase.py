#!/usr/bin/env python
#
# Read a list of polygons from stdin in Castryck format and enter them
# in a database.

import os, sys
from sqlite3 import connect
from multiprocessing import Pool

from betti.polygon import normalize_vertices
from sage.all import Polyhedron, ZZ, walltime

NUMPROCESSES=12


def polygon_values(vertices):
    vertices = normalize_vertices(vertices)

    vertices_str = ",".join(repr(pt).replace(" ", "") for pt in vertices)

    polygon = Polyhedron(vertices)
    volume = float(polygon.volume())
    points = polygon.integral_points()
    num_points = len(points)
    num_interior = len([pt for pt in points if polygon.interior_contains(pt)])

    length = int(max(x for x, y in vertices))
    width = int(max(y for x, y in vertices))

    SymmQQ = polygon.restricted_automorphism_group(output="matrixlist")
    if len(SymmQQ) == 1:
        # Fast path for trivial symmetry group
        symm = 1
    else:
        # Restrict symmetries to symmetries in GL(3, ZZ)
        Symm = [s for s in SymmQQ if abs(s.determinant()) == 1 and all(x in ZZ for x in s.list())]
        if all(s.determinant() > 0 for s in Symm):
            symm = len(Symm)
        else:
            symm = -len(Symm)

    return (vertices_str, len(vertices), volume,
        num_points, num_interior, num_points - num_interior,
        width, length, symm)


def process_data(t):
    counter, vertices = t
    v = (counter,) + polygon_values(vertices)

    vstr = "   ".join(repr(x) for x in v)
    t1 = walltime()
    sys.stdout.write("({:5.1f}Hz) {}\n".format(counter/float(t1 - t0), vstr))

    return v


def raw_data(f):
    counter = 0
    for line in f:
        if '<' not in line:
            continue

        line = line.replace("<", "(")
        line = line.replace(">", ")")

        # Drop all characters except these:
        line = "".join(c for c in line if c in "[](),-0123456789")

        # Drop trailing comma
        if line.endswith(","):
            line = line[:-1]

        counter += 1
        yield (counter, eval(line))


def create_database(f, db):
    global t0
    t0 = walltime()

    P = Pool(NUMPROCESSES)
    it = P.imap_unordered(process_data, raw_data(f), chunksize=100)

    con = connect(db)
    con.execute('''
        CREATE TABLE polygons(
            rowid INTEGER PRIMARY KEY,
            vertices TEXT,
            num_vertices INTEGER,
            volume REAL,
            num_points INTEGER,
            num_interior INTEGER,
            num_border INTEGER,
            width INTEGER,
            length INTEGER,
            symm INTEGER)
    ''')

    con.executemany('INSERT INTO polygons VALUES (?,?,?,?,?,?,?,?,?,?)', it)
    con.commit()
    con.close()


db = sys.argv[1]
if os.path.exists(db):
    os.remove(db)

create_database(sys.stdin, db)
