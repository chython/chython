from numpy import uint32, empty
from math import sqrt
cimport cython


@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
def _possible_bonds(double[:,::1] xyz, double[::1] radii, double multiplier):
    cdef int size, max_bonds, c=0, n, m
    cdef double r, d, nx, ny, nz, mx, my, mz, rn
    size = xyz.shape[0]
    max_bonds = size * 10  # each atom has less then 10 neighbors approximately
    cdef unsigned int[:,:] nm = empty((max_bonds, 2), dtype=uint32)
    cdef double[:] ds = empty(max_bonds)

    for n in range(size - 1):
        nx, ny, nz = xyz[n]
        rn = radii[n]
        for m in range(n + 1, size):
            mx, my, mz = xyz[m]
            d = sqrt((nx - mx) ** 2 + (ny - my) ** 2 + (nz - mz) ** 2)
            r = (rn + radii[m]) * multiplier
            if d <= r:
                nm[c,0] = n + 1
                nm[c,1] = m + 1
                ds[c] = d
                c += 1
    return nm[:c], ds[:c]
