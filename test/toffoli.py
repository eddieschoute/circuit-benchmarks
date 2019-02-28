import unittest
from unittest import TestCase

from circuit_benchmarks.toffoli import *
from qiskit import Aer, execute

class TestToffoli(TestCase):
    def test_3_ntoffoli(self) -> None:
        number_qubits = 3
        toff_circ = toffoli(number_qubits)
        # Construct a simulation circuit
        q = toff_circ.qregs[0] # assert that there is only 1 q.register
        creg = ClassicalRegister(number_qubits)
        circ = QuantumCircuit(q, creg)
        # Apply Hadamards to each control
        for i in range(number_qubits-1):
            circ.h(q[i])
        circ.extend(toff_circ)
        circ.measure(q, creg)

        # simulate the circuit
        simulator = Aer.get_backend('qasm_simulator')
        result = execute(circ, simulator).result()
        counts = result.get_counts(circ)

        # Check that the toffoli was applied.
        self.assertNotIn('011', counts)
        self.assertNotIn('100', counts)
        self.assertNotIn('101', counts)
        self.assertNotIn('110', counts)

