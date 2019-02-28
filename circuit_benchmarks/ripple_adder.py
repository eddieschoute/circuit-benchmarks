# -*- coding: utf-8 -*-

# Copyright 2017, IBM.
#
# This source code is licensed under the Apache License, Version 2.0 found in
# the LICENSE.txt file in the root directory of this source tree.

"""
Ripple adder example based on Cuccaro et al., quant-ph/0410184.
"""

from qiskit import QuantumRegister, ClassicalRegister, QuantumCircuit
from qiskit import compile, BasicAer
from qiskit.circuit import Gate, InstructionSet
from qiskit.dagcircuit import DAGCircuit
from qiskit.extensions.standard import *


def ripple_adder(number_qubits: int):
    """
    A quantum adder circuit

    Finds an adder that uses at most number_qubits qubits.
    For odd number_qubits this means that one qubit is not used.
    """
    n = ((number_qubits - 2) // 2) * 2
    cin = QuantumRegister(1)
    a = QuantumRegister(n)
    cout = QuantumRegister(1)

    # Build a temporary subcircuit that adds a to b,
    # storing the result in b
    adder_circuit = QuantumCircuit(cin, a, cout, name="ripple_adder")
    adder_circuit.majority(cin[0], a[0], a[1])
    for j in range((n - 2)//2):
        adder_circuit.majority(a[2*j + 1], a[2*j + 2], a[2*j + 3])
    adder_circuit.cx(a[n - 1], cout[0])
    for j in reversed(range((n-2)//2)):
        adder_circuit.umaj_add(a[2*j + 1], a[2*j + 2], a[2*j + 3])
    adder_circuit.umaj_add(cin[0], a[0], a[1])

    return adder_circuit


class MajorityGate(Gate):
    def __init__(self, c, b, a, circ=None):
        """Create new Toffoli gate."""
        super().__init__("maj", [], [c, b, a], circ)

    def _define_decompositions(self):
        decomposition = DAGCircuit()
        q = QuantumRegister(3)
        decomposition.add_qreg(q)
        c, b, a = q
        rule = [
            CnotGate(a, b),
            CnotGate(a, c),
            ToffoliGate(c, b, a)
        ]
        for inst in rule:
            decomposition.apply_operation_back(inst)
        self._decompositions = [decomposition]

    def inverse(self):
        """Invert this gate."""
        raise NotImplementedError("Inverse not implemented")

    def reapply(self, circ):
        """Reapply this gate to corresponding qubits in circ."""
        self._modifiers(circ.majority(self.qargs[0], self.qargs[1], self.qargs[2]))


def majority(self, c, b, a):
    """Apply a majority function on the given qubits"""
    if isinstance(c, QuantumRegister) and \
            isinstance(a, QuantumRegister) and \
            isinstance(b, QuantumRegister) and \
            len(c) == len(b) and len(b) == len(a):
        instructions = InstructionSet()
        for i in range(c.size):
            instructions.add(self.majority((c, i), (b, i), (a, i)))
        return instructions

    self._check_qubit(c)
    self._check_qubit(b)
    self._check_qubit(a)
    self._check_dups([c, b, a])
    return self._attach(MajorityGate(c, b, a, circ=self))


QuantumCircuit.majority = majority


class UnMajAdd(Gate):
    def __init__(self, a, b, c, circ=None):
        """Create new Unmajority and add (UMA) gate."""
        super().__init__("uma", [], [a, b, c], circ)

    def _define_decompositions(self):
        # Decomposition with minimal nr of CNOTs
        decomposition_cnot = DAGCircuit()
        q = QuantumRegister(3)
        a, b, c = q
        decomposition_cnot.add_qreg(q)
        rule = [
            ToffoliGate(a, b, c),
            CnotGate(c, a),
            CnotGate(a, b)
        ]
        for inst in rule:
            decomposition_cnot.apply_operation_back(inst)

        # Decomposition with minimal depth
        decomposition_par = DAGCircuit()
        q = QuantumRegister(3)
        a, b, c = q
        decomposition_par.add_qreg(q)
        rule = [
            XGate(b),
            CnotGate(a, b),
            ToffoliGate(a, b, c),
            XGate(b),
            CnotGate(c, a),
            CnotGate(c, b)
        ]
        for inst in rule:
            decomposition_par.apply_operation_back(inst)
        self._decompositions = [decomposition_cnot, decomposition_par]

    def inverse(self):
        """Invert this gate."""
        raise NotImplementedError("Inverse not implemented")

    def reapply(self, circ):
        """Reapply this gate to corresponding qubits in circ."""
        self._modifiers(circ.umaj_add(self.qargs[0], self.qargs[1], self.qargs[2]))


def umaj_add(self, a, b, c):
    """Apply a majority function on the given qubits"""
    if isinstance(c, QuantumRegister) and \
            isinstance(a, QuantumRegister) and \
            isinstance(b, QuantumRegister) and \
            len(c) == len(b) and len(b) == len(a):
        instructions = InstructionSet()
        for i in range(a.size):
            instructions.add(self.umaj_add((a, i), (b, i), (c, i)))
        return instructions

    self._check_qubit(c)
    self._check_qubit(b)
    self._check_qubit(a)
    self._check_dups([c, b, a])
    return self._attach(UnMajAdd(a, b, c, circ=self))


QuantumCircuit.umaj_add = umaj_add
