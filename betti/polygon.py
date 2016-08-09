from sage.all import (ZZ, GL, gcd, vector, point, Polyhedron,
    MixedIntegerLinearProgram)
from .fundom import fundamental_domain

GL_3_Z = GL(3, ZZ)


def lattice_width(vertices):
    if not vertices:
        return ZZ(-1)

    miny = min(y for x, y in vertices)
    maxy = max(y for x, y in vertices)

    # Interval of Y values
    lattice_width = maxy - miny

    # Try to do better than lattice_width: we propose an affine
    # transformation matrix
    # [ a b 0 ]
    # [ d e f ]
    # [ 0 0 1 ]
    # with the goal of minimizing the interval of Y values.
    compute_width = MixedIntegerLinearProgram(maximization=False, solver="PPL")
    t = compute_width.new_variable(integer=True)
    lw = compute_width.new_variable(integer=True)[0]
    compute_width.set_objective(lw)

    # We want d >= 1 to ensure that the X-axis gets mapped to some
    # line which is not parallel to the original X-axis. If we keep the
    # X-axis parallel, then we get a lattice width equal to lattice_width
    # which does not give new information. If d >= 1, we are also
    # guaranteed that we can make the transformation matrix invertible.
    compute_width.set_min(t[0], 1)
    for x, y in vertices:
        compute_width.add_constraint(0 <= t[0] * x + t[1] * y + t[2] <= lw)

    return min(compute_width.solve(), lattice_width)


def normalize_vertices(vertices):
    if not vertices:
        return vertices

    lw = lattice_width(vertices)

    # Find all values (d,e,f) of transformation matrices
    # [ a b 0 ]
    # [ d e f ]
    # [ 0 0 1 ]
    # realizing the lattice width.
    compute_width = MixedIntegerLinearProgram(maximization=False, solver="PPL")
    t = compute_width.new_variable(integer=True)

    # There are always 2 choices of sign for (d,e) but we want to pick
    # only one. So we assume d >= 0. We still need to take care of the
    # case d == 0 and e < 0 later.
    compute_width.set_min(t[0], 0)

    for x, y in vertices:
        compute_width.add_constraint(0 <= t[0] * x + t[1] * y + t[2] <= lw)

    # Get list of (a,b,d,e,f) transformation matrices
    if lw:
        transformations = []
        for d, e, f in compute_width.polyhedron().integral_points():
            # Take care of sign of (d,e) when d == 0 and skip
            # trivial solution (d,e) == (0,0).
            if d == 0 and e <= 0:
                continue

            g, a, b = e.xgcd(-d)   # a e - b d = 1
            assert g == 1
            # Two cases: det=1 and det=-1
            transformations += [(a,b,d,e,f), (-a,-b,d,e,f)]
    else:
        # Special case lattice width zero: the "compute_width"
        # polyhedron is a vector line and we want the point on that
        # line with content 1.
        d, e, f = compute_width.polyhedron().lines()[0]
        g = d.gcd(e)
        d, e, f = d/g, e/g, f/g
        g, a, b = e.xgcd(-d)   # a e - b d = 1
        assert g == 1
        transformations = [(a,b,d,e,f)]

    best_llg = None
    candidates = []  # Contains lists of vertices with optimal lattice length
    for a, b, d, e, f in transformations:
        newvertices = [(a * x + b * y, d * x + e * y + f) for x, y in vertices]
        newvertices.sort(key=lambda xy: (xy[1], xy[0]))

        # Assert that the range of Y values is the interval [0,lw]
        assert newvertices[0][1] == 0
        assert newvertices[-1][1] == lw

        # Now move the first vertex to the origin and shear such
        # that all points have positive X-coordinate.
        # Note that this shearing keeps the vertices sorted.
        originx = newvertices[0][0]

        # Now shear by a variable amount. We only need those candidates
        # with minimal lattice length. Remark that lattice length is a
        # convex function of the shear distance.
        if lw:
            maxshear = max((originx - x) / y for x, y in newvertices if y).ceil()
            newvertices = [(x - originx + maxshear * y, y) for x, y in newvertices]
        else:
            # All points have Y-coordinate 0: shearing does nothing
            newvertices = [(x - originx, y) for x, y in newvertices]
        newllg = max(x for x, y in newvertices)
        if best_llg is None or newllg < best_llg:
            # We improved the best length: delete all
            # previous candidates
            candidates = []
            best_llg = newllg
        if newllg == best_llg:
            candidates.append(newvertices)

        if not lw:
            continue

        # Now continue shearing in both directions as long as we do not
        # increase the lattice length.
        for s in [-1, 1]:
            shearvertices = newvertices
            while True:
                originx = min(x + s*y for x, y in shearvertices)
                shearvertices = [(x + s*y - originx, y) for x, y in shearvertices]
                llg = max(x for x, y in shearvertices)
                if llg > newllg:
                    break
                if llg < best_llg:
                    # We improved the best length: delete all
                    # previous candidates
                    candidates = []
                    best_llg = llg
                if llg == best_llg:
                    candidates.append(shearvertices)

    candidates180 = []
    for vertices in candidates:
        # Also consider the polygon rotated over 180 degrees. We need
        # this because we considered only 1 sign for (d,e) in the
        # transformation matrix.
        vertices180 = [(best_llg-x, lw-y) for x, y in vertices[::-1]]
        candidates180.append(vertices180)

    candidates += candidates180

    # Now find the best candidate
    def vertices_key(v):
        k = []
        for x, y in v:
            k += [y, x]
        return k

    return min(candidates, key=vertices_key)


def polygon_normalized(vertices):
    return Polyhedron(normalize_vertices(vertices))


def polygon_symmetries(P):
    AUT = P.restricted_automorphism_group(output="matrix")
    return [g.matrix() for g in AUT.intersection(GL_3_Z)]


def interior_points(P):
    return tuple(pt for pt in P.integral_points() if P.interior_contains(pt))


def fundamental_points(F, AUT):
    return [orbit[0] for orbit in orbits(F.integral_points(), AUT)]


def orbits(L, G):
    """
    Given a collection `L` of points, group them in orbits under the
    action of `G` (a collection of matrices representing elements of
    AGL.).
    """
    from collections import defaultdict
    orbits = defaultdict(list)
    for pt in L:
        v = vector(ZZ, tuple(pt) + (1,))
        orbit = frozenset(tuple(g*v)[:-1] for g in G)
        orbits[orbit] = orbits[orbit] + [pt]
    return orbits.values()


def plot_polygon(Delta, remove_points=[]):
    P = Delta.integral_points()
    AUT = polygon_symmetries(Delta)
    F = fundamental_domain(Delta, AUT)
    fundam = fundamental_points(F, AUT)

    plt = Delta.plot(fill="yellow", point=False, zorder=-10)
    plt += F.plot(fill=(1,0.9,0), point=False, zorder=-9)

    for pt in P:
        style = {'pointsize': 80, 'color': (0.4,0.5,1)}
        if pt in remove_points:
            style['color'] = (0.7,0,0)
            style['marker'] = "x"
        elif Delta.interior_contains(pt):
            if pt in fundam:
                style['color'] = (0.2,0.8,0)
            else:
                style['color'] = (0,1,0)
        plt += point(pt, **style)
    return plt


def point_lists(Delta, remove_points=[]):
    all_pts = list(Delta.integral_points())
    int_pts = list(interior_points(Delta))
    int_pts2 = list(interior_points(ZZ(2) * Delta))

    for pt in remove_points:
        all_pts.remove(pt)
        for z in int_pts:
            int_pts2.remove(pt + z)

    return all_pts, int_pts, int_pts2


def polygon_expand(Delta, len=1):
    ieqs = []
    for ineq in Delta.inequalities_list():
        c = ineq[0]
        ineq = ineq[1:]
        assert c in ZZ
        assert gcd(ineq) == 1
        ieqs.append([c + len] + ineq)
    return Polyhedron(ieqs=ieqs)


def is_interior(Delta):
    Gamma = polygon_expand(Delta)
    return all(c in ZZ for v in Gamma.vertices_list() for c in v)


def normal_vertices(Delta):
    vertices = []
    for ineq in Delta.inequalities_list():
        c = ineq[0]
        ineq = ineq[1:]
        assert gcd(ineq) == 1
        vertices.append(vector(ZZ, ineq))
    return vertices


def is_weak_fano(Delta):
    """
    Check whether Delta is a (weak) Fano polytope.

    OUTPUT:

    - 2 if Fano

    - 1 if weak Fano

    - 0 if not weak Fano
    """
    N = normal_vertices(Delta)
    Normal = Polyhedron(vertices=N)
    if len(interior_points(Normal)) == 1:
        # Weak Fano => Fano iff all elements of N
        # are vertices
        if Normal.n_vertices() == len(N):
            return 2
        else:
            return 1
    return 0


def is_smooth(Delta):
    for v in Delta.vertices():
        # There should be exactly two edges incident with v
        h1, h2 = v.incident()
        _, a, b = h1.vector()
        _, c, d = h2.vector()
        det = a * d - b * c
        if abs(det) != 1:
            assert det != 0
            return False
    return True
