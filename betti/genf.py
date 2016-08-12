from sage.all import (ZZ, PolynomialRing, PowerSeriesRing,
        prod, factorial, cached_function)

ZZ_X_Y = PolynomialRing(ZZ, ('X', 'Y'), order="lex")
ZZ_X_Y_T = PowerSeriesRing(ZZ_X_Y, 'T')


def genf_points_multiples(pts, N):
    X, Y = ZZ_X_Y.gens()

    G = [ZZ_X_Y(1)] + [ZZ_X_Y(0)]*N
    for pt in pts:
        Z = X**pt[0] * Y**pt[1]
        for i in range(1, N+1):
            G[i] += Z**i
    return tuple(G)


def genf_points(pts, t=[1]):
    t = tuple(t)
    G = genf_points_multiples(pts, sum(t))
    return formula_partition_tuple(G, t)


def genf_polygon(P, t=[1]):
    return genf_points(P.integral_points(), t)


def genf_polygon_interior(P, t=[1]):
    from .polygon import interior_points
    return genf_points(interior_points(P), t)


def genf_subsets(pts, k):
    X, Y = ZZ_X_Y.gens()
    T = ZZ_X_Y_T.gen().add_bigoh(k + 1)
    return prod(1 + X**pt[0] * Y**pt[1] * T for pt in pts)[k]


def genf_koszul_dual(P, Q, p):
    if p < 0:
        return ZZ_X_Y()
    f1 = genf_polygon(P, [1]*p) // factorial(p)
    f2 = genf_polygon_interior(Q)
    return f1 * f2


def genf_koszul_points(P, Q, p):
    if p < 0:
        return ZZ_X_Y()
    f1 = genf_points(P, [1]*p) // factorial(p)
    f2 = genf_points(Q)
    return f1 * f2


def genf_diagonal_difference(Delta, l):
    f = ZZ_X_Y()
    for j in range(l+2):
        j = ZZ(j)
        t = genf_koszul_points(Delta.integral_points(), (j*Delta).integral_points(), l+1-j)
        if j % 2:
            f -= t
        else:
            f += t
    return f


def formula_partition(t):
    t = tuple(t)
    R = PolynomialRing(ZZ, ["f"+str(i+1) for i in range(sum(t))])
    G = (R(1),) + R.gens()
    return formula_partition_tuple(G, t)


@cached_function
def formula_partition_tuple(G, t):
    if not t:
        return G[0]
    if len(t) == 1:
        return G[t[0]]

    t0 = t[:1]
    t1 = t[1:]

    f = formula_partition_tuple(G, t0) * formula_partition_tuple(G, t1)

    # Substract doubly-counted tuples
    for i in range(1, len(t)):
        u = list(t)
        a = u.pop(i)
        u[0] += a
        f -= formula_partition_tuple(G, tuple(u))

    return f
