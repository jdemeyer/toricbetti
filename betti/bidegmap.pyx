# distutils: language = c++
# distutils: libraries = givaro
# distutils: extra_compile_args = -DDISABLE_COMMENTATOR
# cython: cdivision = True

from libc.stdint cimport uint32_t
from cpython.array cimport array
cimport cython
include "cysignals/memory.pxi"

from sage.rings.integer cimport Integer
from sage.all import binomial

cdef ZERO = Integer(0)


# Type for sequence number
ctypedef uint32_t seq_t

cdef extern from "linbox/field/Modular/modular-double.h" nogil:
    cdef cppclass Fp "LinBox::Modular<double>":
        Fp()
        Fp(unsigned long p)

cdef extern from "linbox/vector/vector-traits.h" nogil:
    pass

cdef extern from "linbox/blackbox/sparse.h" nogil:
    cdef cppclass SparseMatrix "LinBox::SparseMatrix<LinBox::Modular<double>, LinBox::Vector<LinBox::Modular<double> >::SparseSeq>":
        SparseMatrix(Fp F, size_t m, size_t n)
        size_t coldim()
        size_t rowdim()
        void setEntry(size_t i, size_t j, double value)

cdef extern from "linbox/solutions/methods.h" nogil:
    cdef cppclass SparseElimination "LinBox::Method::SparseElimination":
        pass

cdef extern from "linbox/solutions/rank.h" nogil:
    unsigned long rank(unsigned long r, SparseMatrix M, SparseElimination m)


@cython.final
cdef class FpKoszulMatrix(object):
    cdef Fp field
    cdef SparseMatrix* m

    def __cinit__(self, unsigned long p, size_t y, size_t x):
        self.field = Fp(p)
        self.m = new SparseMatrix(self.field, y, x)

    def __dealloc__(self):
        del self.m

    def __repr__(self):
        return "Sparse matrix in dimensions {0} -> {0}",format(self.coldim(), self.rowdim())

    def rank(self):
        if self.ncols() == 0 or self.nrows() == 0:
            return 0

        cdef unsigned long r
        cdef SparseElimination m
        with nogil:
            sig_on()
            rank(r, self.m[0], m)
            sig_off()
        return int(r)

    def ncols(self):
        return int(self.m.coldim())

    def nrows(self):
        return int(self.m.rowdim())


assert array("I").itemsize == sizeof(seq_t)


# Precompute all binomial sums up to this limit
DEF BINOMIAL_SUM_CACHE_SIZE = 50

# Size of bideg bitmap (must be multiple of 8)
DEF DEG_SIZE = 256


cdef seq_t binomial_sum_cache[BINOMIAL_SUM_CACHE_SIZE * BINOMIAL_SUM_CACHE_SIZE * BINOMIAL_SUM_CACHE_SIZE]

cpdef inline seq_t binomial_sum(int l, int g, int k):
    """
    Return sum_(1 <= i <= k) binomial(g-i, l)
    """
    return binomial_sum_cache[l*(BINOMIAL_SUM_CACHE_SIZE**2) + g*(BINOMIAL_SUM_CACHE_SIZE) + k]

cdef precompute():
    cdef size_t size = BINOMIAL_SUM_CACHE_SIZE * BINOMIAL_SUM_CACHE_SIZE * BINOMIAL_SUM_CACHE_SIZE

    cdef int l, g, k
    cdef size_t idx = 0
    for l in range(BINOMIAL_SUM_CACHE_SIZE):
        for g in range(BINOMIAL_SUM_CACHE_SIZE):
            N = ZERO
            for k in range(BINOMIAL_SUM_CACHE_SIZE):
                try:
                    binomial_sum_cache[idx] = N
                except OverflowError:
                    binomial_sum_cache[idx] = <seq_t>(-1)
                idx += 1
                if l+k+1 <= g:
                    N += binomial(g-1-k, l)


cdef seq_t sequence_number(int* seq, int L, int N):
    """
    Given an increasing sequence seq[0] < ... < seq[L-1] of numbers in the interval [0, ..., N-1],
    return the index of this sequence in the lexicographic ordering.
    """
    cdef seq_t S = 0
    cdef int i
    cdef int m = 0
    for i in range(L):
        L -= 1
        S += binomial_sum(L, N-m, seq[i]-m)
        m = seq[i] + 1
    return S

cdef seq_t sequence_number_hole(int* seq, int L, int N, int hole):
    """
    Given an increasing sequence seq[0] < ... < seq[L-1] of numbers in the interval [0, ..., N-1],
    return the index of the sub-sequence seq[0] < ... < seq[hole-1] < seq[hole+1] < ... < seq[L-1]
    in the lexicographic ordering.
    """
    cdef seq_t S = 0
    cdef int i
    cdef int m = 0
    L -= 1  # effective length is 1 less
    for i in range(L):
        if i == hole:  # Shift seq
            seq += 1
        L -= 1
        S += binomial_sum(L, N-m, seq[i]-m)
        m = seq[i] + 1
    return S

cdef int sequence_next(int* seq, int L, int N):
    """
    Given an increasing sequence seq[0] < ... < seq[L-1] of numbers in the interval [0, ..., N-1],
    change seq to the next such sequence and return 1. If the given sequence was the last sequence,
    leave seq unchanged and return 0.
    """
    if L <= 0:  # The empty sequence is always the last one
        return 0
    cdef int slack = N - L
    cdef int i = L - 1
    while seq[i] >= i + slack:
        if i == 0:
            return 0
        i -= 1
    seq[i] += 1
    i += 1
    while i < L:
        seq[i] = seq[i-1] + 1
        i += 1
    return 1

def py_sequence_number(seq, N):
    cdef Py_ssize_t L = len(seq)
    cdef array arr = array("i", seq)
    return sequence_number(arr.data.as_ints, L, N)

def py_sequence_number_hole(seq, N, hole):
    cdef Py_ssize_t L = len(seq)
    cdef array arr = array("i", seq)
    return sequence_number_hole(arr.data.as_ints, L, N, hole)

def py_sequence_next(seq, N):
    cdef Py_ssize_t L = len(seq)
    cdef array arr = array("i", seq)
    cdef int r = sequence_next(arr.data.as_ints, L, N)
    if not r:
        raise StopIteration
    return list(arr)


cpdef array bideg_basis(int L, polygon_pts, tensor_pts, long a, long b):
    """
    Let V be the vector space spanned by polygon and W be the vector
    space spanned by tensor.

    Return a basis for Wedge^L(V) consisting of those elements
    which can be tensored with an element of W to obtain an element
    of bidegree (a,b).
    """
    cdef int N = len(polygon_pts)
    if N >= BINOMIAL_SUM_CACHE_SIZE:
        raise OverflowError("number of points must be less than {}".format(BINOMIAL_SUM_CACHE_SIZE))
    cdef seq_t seqmax = binomial(N, L)

    cdef int i
    cdef uint32_t block, bit
    cdef long x, y

    # Allocate empty basis
    cdef array basis = array("I")

    if L < 0:
        return basis

    # Check degree bounds
    cdef long degmax
    if L > 0:
        degmax = (DEG_SIZE - 1) // L
        for x, y in polygon_pts:
            if x < 0 or y < 0:
                raise ValueError("polygon contains points with negative degree")
            if x > degmax or y > degmax:
                raise ValueError("polygon contains points with degree larger than {}".format(degmax))

    # Create bideg bitmap
    cdef size_t* bidegmap = <size_t*>check_calloc(DEG_SIZE, DEG_SIZE//8)
    for x, y in tensor_pts:
        x = a - x
        y = b - y
        if x < 0 or y < 0 or x >= DEG_SIZE or y >= DEG_SIZE:
            continue
        bit = x + y * DEG_SIZE
        block = bit // (sizeof(size_t) * 8)
        bit = bit % (sizeof(size_t) * 8)
        bidegmap[block] |= (<size_t>(1) << bit)

    # Convert bideg of polygon points to bitmap indices
    cdef uint32_t* bideg = <uint32_t*>check_allocarray(N, sizeof(uint32_t))
    for i in range(N):
        x, y = polygon_pts[i]
        bideg[i] = x + y * DEG_SIZE

    # Initial sequence
    cdef int* seq = <int*>check_allocarray(L, sizeof(int))
    for i in range(L):
        seq[i] = i
    cdef seq_t seqnum = 0

    # Loop over all sequences
    while True:
        sig_check()
        bit = 0
        for i in range(L):
            bit += bideg[seq[i]]
        block = bit // (sizeof(size_t) * 8)
        bit = bit % (sizeof(size_t) * 8)
        if bidegmap[block] & (<size_t>(1) << bit):
            basis.append(seqnum)

        if not sequence_next(seq, L, N):
            break
        seqnum += 1

    sig_free(bidegmap)
    sig_free(bideg)
    sig_free(seq)

    return basis


def bideg_map_points(F, int L, basis, tensor, tensor2, long a, long b):
    cdef int i
    cdef seq_t j
    cdef unsigned long p = F.characteristic()

    cdef int N = len(basis)
    domain = bideg_basis(L, basis, tensor, a, b)
    codomain = bideg_basis(L-1, basis, tensor2, a, b)

    cdef seq_t* domainptr = domain.data.as_uints
    cdef seq_t* codomainptr = codomain.data.as_uints

    # Create Koszul matrix
    cdef FpKoszulMatrix M = FpKoszulMatrix(p, len(codomain), len(domain))

    # Trivial case
    if len(codomain) == 0 or len(domain) == 0:
        return M

    # Sequence -> Index map for codomain (-1 means: not in codomain)
    cdef seq_t B = binomial(N, L-1)
    cdef seq_t* codomain_idx = <seq_t*>check_allocarray(B, sizeof(seq_t))
    for j in range(B):
        codomain_idx[j] = <seq_t>(-1)
    for j in range(len(codomain)):
        assert codomain[j] < B
        codomain_idx[codomain[j]] = j

    # Initial sequence
    cdef int* seq = <int*>check_allocarray(L, sizeof(int))
    for i in range(L):
        seq[i] = i
    cdef seq_t seqnum = 0

    # Matrix entries
    cdef double pos1 = 1
    cdef double neg1 = p - 1
    cdef double val

    # Loop over all sequences from the domain
    j = 0
    cdef seq_t idx
    while True:
        sig_check()
        if seqnum == domainptr[j]:
            for i in range(L):
                idx = codomain_idx[sequence_number_hole(seq, L, N, i)]
                if idx == <seq_t>(-1):
                    continue
                if i % 2 == 0:
                    val = pos1
                else:
                    val = neg1
                M.m.setEntry(idx, j, val)
            j += 1
            if j >= len(domain):
                break
        if not sequence_next(seq, L, N):
            raise AssertionError("sequences exhausted")
        seqnum += 1

    sig_free(codomain_idx)
    sig_free(seq)

    return M


precompute()
