"""useful things
"""
import numpy as np
from sympy.combinatorics import Permutation
import sympy


###############################
### Standard form of a code ###
###############################


def find_next_pivot(matrix_stab, j, n, rx):
    """find next pivot position
    """
    for kp in range(j, n):
        for jp in range(j, rx):
            if matrix_stab[jp, kp] != 0:
                return (jp, kp)
    return (-1, -1)


def standard_form(matrix_stab_in, n, rx):
    """put matrix_stab in standard form
    """
    matrix_stab = np.array(matrix_stab_in, dtype='uint8')
    qubitpermutation = Permutation([], size=n)

    for j in range(rx):
        jp, kp = find_next_pivot(matrix_stab, j, n, rx)
        if jp == -1:
            break
        if j != jp:
            perm = Permutation([[j, jp]], size=rx)
            matrix_stab = matrix_stab[np.array([i^perm for i in range(rx)]), :]
        if j != kp:
            perm = Permutation([[j, kp]], size=n)
            matrix_stab = matrix_stab[:, np.array([i^perm for i in range(n)])]
            qubitpermutation = qubitpermutation*perm

        for l in range(j+1, rx):
            if matrix_stab[l, j] == 1:
                matrix_stab[l, :] = (matrix_stab[l, :] + matrix_stab[j, :]) % 2

    for l in range(j-1, 0, -1):
        for j in range(l):
            if matrix_stab[j, l] == 1:
                matrix_stab[j, :] = (matrix_stab[j, :] + matrix_stab[l, :]) % 2

    return matrix_stab, qubitpermutation


############################
### Bravyi Haah matrices ###
############################

ZERO = sympy.GF(2).zero
ONE = sympy.GF(2).one

LMAT = sympy.Matrix([[ZERO, ZERO, ZERO, ZERO, ONE, ONE, ONE, ONE],
                     [ZERO, ZERO, ZERO, ZERO, ONE, ONE, ONE, ONE]])

MMAT = sympy.Matrix([[ONE, ONE, ONE, ZERO, ZERO, ZERO],
                     [ZERO, ZERO, ZERO, ONE, ONE, ONE]])

S1MAT = sympy.Matrix([[ZERO, ONE, ZERO, ONE, ZERO, ONE, ZERO, ONE],
                      [ZERO, ZERO, ONE, ONE, ZERO, ZERO, ONE, ONE],
                      [ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE]])

S2MAT = sympy.Matrix([[ONE, ZERO, ONE, ONE, ZERO, ONE],
                      [ZERO, ONE, ONE, ZERO, ONE, ONE],
                      [ZERO, ZERO, ZERO, ZERO, ZERO, ZERO]])

def bravyi_haah_gmat(k):
    """Tri-orthogonal matrix from Bravyi-Haah
    """
    bigl = np.block([[LMAT] for _ in range(k)]+ [[S1MAT]])
    bigm = sympy.diag(*[MMAT for _ in range(k)])
    bigs2 = np.block([S2MAT for _ in range(k)])
    gmatr = np.block([[bigm], [bigs2]])
    gmat = np.block([bigl, gmatr])
    return sympy.Matrix(gmat)

def bravyi_haah_stabx(k):
    """X-stabilizer matrix from Bravyi-Haah
    """
    return np.block([S1MAT] + [S2MAT for _ in range(k)])

def bravyi_haah_stabz(k):
    """Z-stabilizer matrix from Bravyi-Haah
    """
    gmat = bravyi_haah_gmat(k)
    gorth = np.array([g.transpose() for g in gmat.nullspace()], dtype='uint8') % 2
    return gorth
