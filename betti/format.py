from sage.all import ZZ, Matrix
from .genf import genf_diagonal_difference


def format_table(typ, *args, **kwds):
    if typ == "b":
        return format_table_b(*args, **kwds)
    elif typ == "c":
        return format_table_c(*args, **kwds)
    else:
        raise ValueError("unknown typ {!r}".format(typ))


def format_table_b(table, l, dual=False, **kwds):
    # Simple table without removed points
    table0 = table.__class__(table.polygon, ZZ, remove_points=[])
    N = table0.N

    # Definition
    A = table0.computer(l, 1)
    Asupp, Ad, Ac = A.bidegree_tables(short=True, flip=dual)

    # Reduce support via dual
    B = table0.dual_computer(N - 3 - l, 2)
    Bsupp, Bd, Bc = B.bidegree_tables(short=True, flip=not dual)
    Mzero = Bsupp.parent()(0)

    # Intersect Asupp and Bsupp
    nr = min(Asupp.nrows(), Bsupp.nrows())
    nc = min(Asupp.ncols(), Bsupp.ncols())
    Supp = Matrix(ZZ, nr, nc)
    for i in range(nr):
        for j in range(nc):
            Supp[i,j] = Asupp[i,j] and Bsupp[i,j]

    # Compute b_l using the corresponding c_cl and the formula for the
    # diagonal differences
    lc = N - 1 - l
    C = table.c_computer(lc)
    if N - 1 - lc > table.hering_schenck_bound:
        Msupp, Md, Mc, Mker, Mim = C.bidegree_tables(flip=not dual)
        Cohom = Mker
    else:
        Msupp, Md, Mc = C.bidegree_tables(short=True, flip=not dual)
        Cohom = Mzero

    # Apply diagonal differences
    dd = genf_diagonal_difference(table.polygon, l)
    for i in range(nr):
        for j in range(nc):
            Cohom[i,j] -= dd[i,j]

    return format_matrices(Cohom, Md, Supp, **kwds)


def format_table_c(table, l, dual=False, **kwds):
    # Simple table without removed points
    table0 = table.__class__(table.polygon, ZZ, remove_points=[])
    N = table0.N

    # Definition
    A = table0.computer(N - 2 - l, 2)
    Asupp, Ad, Ac = A.bidegree_tables(short=True, flip=dual)

    # Reduce support via dual
    B = table0.dual_computer(l - 1, 1)
    Bsupp, Bd, Bc = B.bidegree_tables(short=True, flip=not dual)
    Mzero = Bsupp.parent()(0)

    # Intersect Asupp and Bsupp
    nr = min(Asupp.nrows(), Bsupp.nrows())
    nc = min(Asupp.ncols(), Bsupp.ncols())
    Supp = Matrix(ZZ, nr, nc)
    for i in range(nr):
        for j in range(nc):
            Supp[i,j] = Asupp[i,j] and Bsupp[i,j]

    C = table.c_computer(l)
    if N - 1 - l > table.hering_schenck_bound:
        Msupp, Md, Mc, Mker, Mim = C.bidegree_tables(flip=not dual)
        Cohom = Mker
    else:
        Msupp, Md, Mc = C.bidegree_tables(short=True, flip=not dual)
        Cohom = Mzero

    return format_matrices(Cohom, Md, Supp, **kwds)


def format_table_c_coimage(table, l, **kwds):
    C = table.c_computer(l)
    Msupp, Md, Mc, Mker, Mim = C.bidegree_tables()

    return format_matrices(Mc - (Md - Mker), Md, Msupp, **kwds)


def format_matrices(M, Mnontriv, Msupp, format, checksum=None, **kwds):
    if checksum is not None:
        sum = ZZ()
        for i in range(Msupp.nrows()):
            for j in range(Msupp.ncols()):
                if Msupp[i,j]:
                    sum += M[i,j]
        if sum != checksum:
            raise AssertionError("sum ({}) does not equal checksum ({})".format(sum, checksum))

    if format == "tikz":
        return tikz_table(M, Mnontriv, Msupp, **kwds)
    elif format == "html":
        return html_table(M, Mnontriv, Msupp, **kwds)
    else:
        raise ValueError("unknown format {!r}".format(format))


def html_table(M, Mnontriv, Msupp):
    supp = [(x,y)
        for y in range(Msupp.nrows()) for x in range(Msupp.ncols())
        if Msupp[y,x]]

    if not supp:
        return ""

    minx = min(bideg[0] for bideg in supp)
    maxx = max(bideg[0] for bideg in supp)
    miny = min(bideg[1] for bideg in supp)
    maxy = max(bideg[1] for bideg in supp)

    s = '<table>\n'

    for y in range(maxy, miny-1, -1):
        s += '<tr>'
        s += '<th class="y">{}</th>'.format(y)
        for x in range(minx, maxx+1):
            if Msupp[y,x]:
                if Mnontriv[y,x]:
                    s += '<td class="nontriv"'
                else:
                    s += '<td class="triv"'
                if M[y,x]:
                    s += '>{}</td>'.format(M[y,x])
                else:
                    s += ' style="color:#bbb">{}</td>'.format(M[y,x])
            else:
                s += '<td></td>'
        s += '</tr>\n'

    s += '<tr><th></th>'
    for x in range(minx, maxx+1):
        s += '<th class="x">{}</th>'.format(x)
    s += '</tr>\n'

    s += '</table>'

    return s


def tikz_table(M, Mnontriv, Msupp, options=""):
    supp = [(x,y)
        for y in range(Msupp.nrows()) for x in range(Msupp.ncols())
        if Msupp[y,x]]

    minx = min(bideg[0] for bideg in supp)
    maxx = max(bideg[0] for bideg in supp)
    miny = min(bideg[1] for bideg in supp)
    maxy = max(bideg[1] for bideg in supp)

    s = "\\begin{tikzpicture}["
    s += "trivial/.style={color=yellow!20},"
    s += "nontrivial/.style={color=green!20},"
    s += options
    s += "]\n"

    for x in range(minx, maxx+1):
        X = float(x - minx)
        Y = -1.25
        s += "\\node at (%s, %s) {$\\mathbf{%s}$};\n" % (X, Y, x)

    for y in range(miny, maxy+1):
        X = -1.25
        Y = float(y - miny)
        s += "\\node at (%s, %s) {$\\mathbf{%s}$};\n" % (X, Y, y)

    for x, y in supp:
        X = float(x - minx)
        Y = float(y - miny)

        if Mnontriv[y,x]:
            style = "nontrivial"
        else:
            style = "trivial"
        s += "\\fill[%s] (%s, %s) rectangle (%s, %s);\n" % (style, X-0.5, Y-0.5, X+0.5, Y+0.5)
        s += "\\node at (%s, %s) {$%s$};\n" % (X, Y, M[y,x])

    s += "\\end{tikzpicture}\n"

    return s
