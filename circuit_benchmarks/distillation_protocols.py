"""Distillation protocols
"""

import numpy as np
from qiskit import ClassicalRegister, QuantumRegister, QuantumCircuit
from distillutils import standard_form, bravyi_haah_stabz, bravyi_haah_stabx


def a_state(qr, j):
    """create a |A> state on qr[j]
    """
    qc = QuantumCircuit(qr, name='a_state')
    qc.h(qr[j])
    qc.t(qr[j])
    return qc


def h_state(qr, j):
    """create a |H> state on qr[j]
    """
    qc = QuantumCircuit(qr, name='a_state')
    qc.h(qr[j])
    qc.t(qr[j])
    return qc
def init_magic_register(magicregister):
    """create |A> states on all qubits of
    magicregister
    """
    qc = QuantumCircuit(magicregister, name='init magic register')
    for j in range(magicregister.size):
        qc += a_state(magicregister, j)

    return qc


def t_gate_injection(magicregister, j, qr, k):
    """perform a T gate on qr[k] using a magic state
    located in magicregister[j]
    """
    cr = ClassicalRegister(1)
    qc = QuantumCircuit(magicregister, qr, cr, name='t_gate_injection')
    qc.cx(magicregister[j], qr[k])
    qc.measure(magicregister[j], cr[0])

    qc.s(qr[k]).c_if(cr, 1)

    return qc



def css_encoding_plus(matrix_stabz, n, rz, name=None):
    """assumes matrix_stab is full rank
    """

    qr = QuantumRegister(n)
    qc = QuantumCircuit(qr, name=name)

    std_form_stb, perm = standard_form(matrix_stabz, n, rz)

    for j in range(rz, n):
        qc.h(qr[j^perm])

    for j in range(rz):
        for k in range(rz, n):
            if std_form_stb[j, k] == 1:
                qc.cx(qr[k^perm], qr[j^perm])
    return qc


def measure_x_syndrome(matrix_stabx, qr, n, rx):
    """measurement circuit for the x syndrome
    given by X-stabilizer matrix matrix_stabx
    """
    cr = ClassicalRegister(rx)
    ancillas = QuantumRegister(rx)
    synextr = QuantumCircuit(qr, ancillas, cr)
    for j in range(rx):
        synextr.h(ancillas[j])
        for k in np.nonzero(matrix_stabx[j, :])[0]:
            synextr.cx(ancillas[j], qr[int(k)])
        synextr.h(ancillas[j])
        synextr.measure(ancillas[j], cr[j])
    return synextr


def standard_distillation(matrix_stabz, matrix_stabx, n, rz, rx, name=None):
    """Structure for the standard distillation circuit
    """
    enc = css_encoding_plus(matrix_stabz, n, rz)
    qr = enc.qregs[0]

    magicregister = QuantumRegister(n)
    circ = QuantumCircuit(qr, magicregister, name=name)
    circ += init_magic_register(magicregister)
    circ += enc
    for j in range(n):
        circ += t_gate_injection(magicregister, j, qr, j)
    # TODO: in principle there is a potential clifford correction missing here
    circ += measure_x_syndrome(matrix_stabx, qr, n, rx)
    # TODO: proper decoding would be required here after

    return circ


def reed_muller_15():
    """Distillation circuit based on
    [[15, 1, 3]]
    """
    matrix_stabz = np.array([[1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0],
                             [1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0],
                             [1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0],
                             [1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1],
                             [1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
                             [1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0],
                             [1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0],
                             [1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0]], dtype='uint8')
    matrix_stabx = np.array([[1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0],
                             [1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0],
                             [1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0],
                             [1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1]], dtype='uint8')
    return standard_distillation(matrix_stabz, matrix_stabx, 15, 10, 4, name='reed_muller_15')


def bravyi_haah_distillation(k):
    """Distillation circuit based on Bravyi-Haah
    """
    matrix_stabz = bravyi_haah_stabz(k)
    matrix_stabx = bravyi_haah_stabx(k)
    rz, n = matrix_stabz.shape
    rx, _ = matrix_stabx.shape
    return standard_distillation(matrix_stabz, matrix_stabx, n, rz, rx, name='bravyi_haah_{}'.format(k))
