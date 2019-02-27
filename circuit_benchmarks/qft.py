# -*- coding: utf-8 -*

# Copyright 2019, IBM.
#
# This source code is licensed under the Apache License, Version 2.0 found in
# the LICENSE.txt file in the root directory of this source tree.


import math

from qiskit import QuantumRegister, QuantumCircuit


def qft(number_qubits: int):
    """Create quantum fourier transform circuit on quantum register qreg."""
    qreg = QuantumRegister(number_qubits)
    circuit = QuantumCircuit(qreg, name="qft")

    for i in range(number_qubits):
        for j in range(i):
            circuit.cu1(math.pi / float(2 ** (i - j)), qreg[i], qreg[j])
        circuit.h(qreg[i])

    return circuit
