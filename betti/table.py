import sys
import gc
from threading import Thread, Lock
from sage.all import (walltime, cached_method,
        ZZ, vector, binomial, Matrix, diagonal_matrix)

from .polygon import (polygon_symmetries, interior_points,
        fundamental_points, orbits, fundamental_domain)
from .genf import genf_koszul_points, ZZ_X_Y
from .bidegmap import bideg_map_points


def unique_intersection(A, B):
    """
    Return the unique element which is both in A and in B.
    """
    B = list(B)
    L = [x for x in A if x in B]
    if len(L) == 1:
        return L[0]
    elif len(L) < 1:
        raise ValueError("%r and %r have empty intersection" % (A,B))
    else:
        raise ValueError("%r and %r have non-unique intersection %r" % (A,B,L))


def optimal_remove_points(Delta, Fundam):
    l = ZZ(5)

    candidates = [(vector(pt),) for pt in Delta.integral_points()]

    V = Delta.vertices()
    if len(V) == 3:
        candidates.append(tuple(vector(pt) for pt in V))
    elif len(V) == 4:
        VV = list(V)
        WW = []
        V0 = VV[0]
        for edge in Delta.faces(1):
            edge = edge.vertices()
            if edge[0] == V0:
                e = edge[1]
                VV.remove(e); WW.append(e)
            elif edge[1] == V0:
                e = edge[0]
                VV.remove(e); WW.append(e)
        # We should be left with 2 pairs
        assert len(VV) == 2
        assert len(WW) == 2

        candidates.append(tuple(vector(pt) for pt in VV))
        candidates.append(tuple(vector(pt) for pt in WW))

    Delta_pts = Delta.integral_points()
    Delta_int_pts = interior_points(Delta)
    Delta2_int_pts = interior_points(ZZ(2) * Delta)

    # Check monomials in the fundamental domain
    X, Y = ZZ_X_Y.gens()
    AUT = polygon_symmetries(l * Delta)
    lFundam = fundamental_points(l * Fundam, AUT)
    monomials = [X**pt[0] * Y**pt[1] for pt in lFundam]

    def difficulty(remove_points):
        all_pts = list(Delta_pts)
        int_pts = list(Delta_int_pts)
        int_pts2 = list(Delta2_int_pts)
        for pt in remove_points:
            all_pts.remove(pt)
            for z in int_pts:
                int_pts2.remove(pt + z)

        fd_coeff = genf_koszul_points(all_pts, int_pts, l-1).monomial_coefficient
        fc_coeff = genf_koszul_points(all_pts, int_pts2, l-2).monomial_coefficient

        return max(fd_coeff(monom) * fc_coeff(monom) for monom in monomials)

    return min(candidates, key=difficulty)


class ToricSurfaceBettiTable(object):
    """
    Betti table of a toric surface.
    """
    def __init__(self, polygon, base_ring, remove_points=None, use_hering_schenck=True, threads=None, verbose=0):
        self.polygon = polygon
        self.base_ring = base_ring
        self.verbose = verbose
        self.threads = threads

        self.volume = polygon.volume()
        self.polygon_points = self.polygon.integral_points()
        self.N = ZZ(len(self.polygon_points))
        self.symmetry = list(polygon_symmetries(polygon))
        self.fundamental_domain = fundamental_domain(polygon, self.symmetry)
        if remove_points is None:
            self.remove_points = optimal_remove_points(self.polygon, self.fundamental_domain)
        else:
            self.remove_points = [vector(pt) for pt in remove_points]

        if use_hering_schenck:
            R = len([pt for pt in self.polygon_points if not polygon.interior_contains(pt)])
            self.hering_schenck_bound = R - 2
        else:
            self.hering_schenck_bound = 0

    def __repr__(self):
        return "Betti table for %r over %r" % (self.polygon, self.base_ring)

    @cached_method
    def computer(self, p, q):
        return KoszulRankComputer(self, p, q)

    @cached_method
    def dual_computer(self, p, q):
        return DualKoszulRankComputer(self, p, q)

    def b_computer(self, l):
        return self.computer(l, 1)

    def b_computer_im(self, l):
        return self.computer(l+1, 0)

    def b(self, l):
        """
        Compute b_l without diagonal formula
        """
        B = self.b_computer(l)
        B0 = self.b_computer_im(l)
        return B.ker() - B0.dom()

    def b_bidegrees(self, l):
        """
        Compute b_l in each bidegree without diagonal formula
        """
        B = self.b_computer(l)
        B0 = self.b_computer_im(l)
        return B.bidegree_tables()[3] - B0.bidegree_tables(short=True)[1]

    def c_computer(self, l):
        return self.dual_computer(l-1, 1)

    def c_computer_im(self, l):
        return self.dual_computer(l, 0)

    def c(self, l):
        """
        Compute c_l without diagonal formula
        """
        return self.c_computer(l).ker()

    def c_bidegrees(self, l):
        """
        Compute c_l in each bidegree without diagonal formula
        """
        C = self.c_computer(l)
        return C.bidegree_tables()[3]

    def betti_number(self, p, q):
        """
        Return the Betti number b_{p,q}. Either use general theory to
        determine this number or compute it using some b_l or c_l
        (whichever is optimal).
        """
        p = ZZ(p)
        q = ZZ(q)
        if q == 0 and p == 0:
            return ZZ(1)
        elif q == 1:
            return self._betti_number_row1(p)
        elif q == 2:
            b = self._betti_number_row1(p+1)
            N = self.N
            l = N - 2 - p
            diff = (N-1-l)*binomial(N-1, l-1) - 2*self.volume*binomial(N-3, l-1)
            return b - diff
        else:
            return ZZ(0)

    @cached_method
    def _betti_number_row1(self, p):
        if p <= 0 or p > self.N - 3:
            return ZZ(0)

        N = self.N
        l = N - 1 - p

        if p <= self.hering_schenck_bound:
            return (N-1-l)*binomial(N-1, l-1) - 2*self.volume*binomial(N-3, l-1)

        # Compute via b or c?
        B = self.b_computer(p)
        C = self.c_computer(l)

        bdiff = (B.difficulty(), B.pq)
        cdiff = (C.difficulty(), C.pq)

        if bdiff < cdiff:
            B0 = self.b_computer_im(p)
            return B.ker() - B0.dom()
        else:
            diff = (N-1-l)*binomial(N-1, l-1) - 2*self.volume*binomial(N-3, l-1)
            return C.ker() + diff

    def betti_table(self):
        if self.verbose:
            # Print all numbers individually
            for p in range(self.N - 3, 0, -1):
                for q in [2,1]:
                    print("K_{},{} = {}".format(p, q, self.betti_number(p,q)))
                    sys.stdout.flush()
        return Matrix(3, self.N - 2, lambda q,p: self.betti_number(p,q))


class KoszulRankComputer(object):
    """
    Class to compute the dimension of the kernel and image of a
    particular map in the Koszul cohomology complex.
    """
    def __init__(self, table, p, q):
        self.p = p = ZZ(p)
        self.q = q = ZZ(q)
        self.pq = p + q
        self.threads = table.threads
        self.verbose = table.verbose
        self.base_ring = table.base_ring

        self.Delta = table.polygon
        self.Delta_symmetry = table.symmetry
        self.Delta_fundam = table.fundamental_domain
        self.Delta_pts = list(table.polygon_points)

        self.Tensor0_pts = self.tensor_points((q-1) * self.Delta) if q >= 1 else []
        self.Tensor1_pts = self.tensor_points(q * self.Delta)
        self.Tensor2_pts = self.tensor_points((q+1) * self.Delta)

        # Remove points
        for pt in table.remove_points:
            self.Delta_pts.remove(pt)
            for z in self.Tensor1_pts:
                self.Tensor2_pts.remove(pt + z)
            for z in self.Tensor0_pts:
                self.Tensor1_pts.remove(pt + z)

        # Generating function for domain/codomain dimensions
        self.fd = genf_koszul_points(self.Delta_pts, self.Tensor1_pts, p)
        self.fc = genf_koszul_points(self.Delta_pts, self.Tensor2_pts, p-1)

    def __str__(self):
        return "<%i,%i>" % (self.p, self.q)

    def difficulty(self):
        monomials = self.fd.monomials()
        if not monomials:
            return ZZ()
        fd_coeff = self.fd.monomial_coefficient
        fc_coeff = self.fc.monomial_coefficient
        return max(fd_coeff(monom) * fc_coeff(monom) for monom in monomials)

    @staticmethod
    def tensor_points(P):
        return list(P.integral_points())

    @cached_method
    def ker_im(self):
        T0 = walltime()
        S = diagonal_matrix([self.pq, self.pq, 1])
        Si = ~S
        AUT = [S * M * Si for M in self.Delta_symmetry]

        pqDelta = self.pq * self.Delta
        pqDelta_fundam = self.pq * self.Delta_fundam
        fundam = fundamental_points(pqDelta_fundam, AUT)
        self.bidegrees = []
        self.numsymm = {}
        for o in orbits(self.tensor_points(pqDelta), AUT):
            bideg = unique_intersection(o, fundam)
            if self.fd[bideg]:
                # Do not handle bidegrees where the domain is
                # 0-dimensional. These do not contribute to the kernel
                # or image.
                self.bidegrees.append(bideg)
            self.numsymm[bideg] = len(o)

        self.TM = 0
        self.TR = 0
        self.done = 0
        self.results = {}

        # Bidegrees to handle in the computation, sorted according to
        # dimension of codomain (=> roughly by increasing difficulty)
        self._bidegrees_todo = sorted(self.bidegrees, key=lambda bideg: self.fc[bideg])

        self.lock = Lock()
        T1 = walltime()
        if not self.threads:
            self._run()
        else:
            T = [Thread(target=self._run, args=(i,)) for i in range(self.threads)]
            for t in T:
                t.start()
            for t in T:
                t.join()
        TT = walltime() - T0

        dimdomain = ZZ(sum(self.fd[bideg] * self.numsymm[bideg] for bideg in self.bidegrees))
        rank = ZZ(sum(self.results[bideg][1] * self.numsymm[bideg] for bideg in self.bidegrees))
        ker = dimdomain - rank
        if self.verbose:
            sys.stderr.flush()
            print("%-8s[%5i]:%10s   (t = %.2fs + %.2fs + %.2fs = %.2fs)" %
                  (str(self), len(self.bidegrees), ker, T1 - T0, self.TM, self.TR, TT))
            sys.stdout.flush()
        return ker, rank

    def dom(self):
        return self.fd(1,1)

    def codom(self):
        return self.fc(1,1)

    def ker(self):
        return self.ker_im()[0]

    def im(self):
        return self.ker_im()[1]

    def bidegree_tables(self, short=False, flip=False):
        """
        Return a 5-tuple of matrices representing, for each bidegree:

        - [0] whether or not it is in the support
        - [1] the dimension of the domain
        - [2] the dimension of the codomain
        - [3] the dimension of the kernel
        - [4] the dimension of the image

        If ``short=True``, return only a 3-tuple with the first 3
        entries from the above list.

        If ``flip=True``, return the tables for the dual.
        """

        S = diagonal_matrix([self.pq, self.pq, 1])
        Si = ~S
        AUT = [S * M * Si for M in self.Delta_symmetry]

        pqDelta = self.pq * self.Delta
        if flip:
            dx, dy = sum(self.Delta.integral_points())
            xmax = dx
            ymax = dy
        else:
            xmax = max(pt[0] for pt in pqDelta.vertices())
            ymax = max(pt[1] for pt in pqDelta.vertices())
        Msupp = Matrix(ZZ, ymax+1, xmax+1)
        Md = Matrix(ZZ, ymax+1, xmax+1)
        Mc = Matrix(ZZ, ymax+1, xmax+1)

        for bideg in self.tensor_points(pqDelta):
            x, y = bideg
            if flip:
                x = dx - x
                y = dy - y
                if x < 0 or y < 0:
                    continue
            Md[y,x] = self.fd[bideg]
            Mc[y,x] = self.fc[bideg]
            Msupp[y,x] = 1

        if short:
            return Msupp, Md, Mc

        self.ker_im()  # Compute results
        Mker = Matrix(ZZ, ymax+1, xmax+1)
        Mim = Matrix(ZZ, ymax+1, xmax+1)

        for bideg in self.bidegrees:
            v = vector([bideg[0], bideg[1], 1])
            for g in AUT:
                x, y, _ = g * v
                if flip:
                    x = dx - x
                    y = dy - y
                    if x < 0 or y < 0:
                        continue
                Mker[y,x] = self.results[bideg][0]
                Mim[y,x] = self.results[bideg][1]

        return Msupp, Md, Mc, Mker, Mim

    def _run(self, threadnum=0):
        while True:
            try:
                if threadnum < 1:
                    # One thread for the easy cases
                    bideg = self._bidegrees_todo.pop(0)
                else:
                    # Other threads for the hard cases
                    bideg = self._bidegrees_todo.pop()
            except IndexError:
                return

            gc.collect(0)

            t0 = walltime()
            M = bideg_map_points(self.base_ring, self.p,
                                 self.Delta_pts, self.Tensor1_pts, self.Tensor2_pts,
                                 bideg[0], bideg[1])
            # Verify dimensions
            assert M.ncols() == self.fd[bideg]
            assert M.nrows() == self.fc[bideg]
            t1 = walltime()
            bideg_rank = M.rank()
            bideg_dimker = M.ncols() - bideg_rank
            t2 = walltime()
            del M

            self.lock.acquire()
            self.TM += (t1 - t0)
            self.TR += (t2 - t1)
            self.done += 1
            if self.verbose >= 2:
                sys.stderr.write("%-8s%10s:%8i x %-8i ker=%-8iim=%-8i(t = %.2fs) %i/%i\n" % (
                    str(self), bideg, self.fc[bideg], self.fd[bideg],
                    bideg_dimker, bideg_rank,
                    t2 - t0, self.done, len(self.bidegrees)))
                sys.stderr.flush()
            self.results[bideg] = (bideg_dimker, bideg_rank)
            self.lock.release()


class DualKoszulRankComputer(KoszulRankComputer):
    """
    Class to compute the dimension of the kernel and image of a
    particular map in the dual Koszul cohomology complex.
    """
    def __str__(self):
        return "<%i,%i'>" % (self.p, self.q)

    @staticmethod
    def tensor_points(P):
        return [pt for pt in P.integral_points() if P.interior_contains(pt)]
