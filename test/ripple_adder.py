import unittest
from unittest import TestCase

from circuit_benchmarks.ripple_adder import *
from qiskit import Aer, execute


class TestRippleAdder(TestCase):
    def test_maj_011(self) -> None:
        qr = QuantumRegister(3)
        cr = ClassicalRegister(3)
        qc = QuantumCircuit(qr, cr)

        qc.x(qr[0])
        qc.x(qr[1])
        qc.majority(qr[0], qr[1], qr[2])
        qc.measure(qr, cr)

        # simulate the circuit
        simulator = Aer.get_backend('qasm_simulator')
        result = execute(qc, simulator).result()
        counts = result.get_counts(qc)

        self.assertEqual({"111"}, set(counts.keys()))

    def test_maj_010(self) -> None:
        qr = QuantumRegister(3)
        cr = ClassicalRegister(3)
        qc = QuantumCircuit(qr, cr)

        qc.x(qr[1])
        qc.majority(qr[0], qr[1], qr[2])
        qc.measure(qr, cr)

        # simulate the circuit
        simulator = Aer.get_backend('qasm_simulator')
        result = execute(qc, simulator).result()
        counts = result.get_counts(qc)

        self.assertEqual({"010"}, set(counts.keys()))

    def test_maj_111(self) -> None:
        qr = QuantumRegister(3)
        cr = ClassicalRegister(3)
        qc = QuantumCircuit(qr, cr)

        qc.x(qr)
        qc.majority(qr[0], qr[1], qr[2])
        qc.measure(qr, cr)
        print(qc)

        # simulate the circuit
        simulator = Aer.get_backend('qasm_simulator')
        result = execute(qc, simulator).result()
        counts = result.get_counts(qc)

        self.assertEqual({"100"}, set(counts.keys()))

    def test_maj_100(self) -> None:
        qr = QuantumRegister(3)
        cr = ClassicalRegister(3)
        qc = QuantumCircuit(qr, cr)

        qc.x(qr[2])
        qc.majority(qr[0], qr[1], qr[2])
        qc.measure(qr, cr)
        print(qc)

        # simulate the circuit
        simulator = Aer.get_backend('qasm_simulator')
        result = execute(qc, simulator).result()
        counts = result.get_counts(qc)

        self.assertEqual({"011"}, set(counts.keys()))

    def test_uma_100(self) -> None:
        qr = QuantumRegister(3)
        cr = ClassicalRegister(3)
        qc = QuantumCircuit(qr, cr)

        qc.x(qr[2])
        qc.umaj_add(qr[0], qr[1], qr[2])
        qc.measure(qr, cr)

        # simulate the circuit
        simulator = Aer.get_backend('qasm_simulator')
        result = execute(qc, simulator).result()
        counts = result.get_counts(qc)
        print(counts)

        self.assertEqual({"111"}, set(counts.keys()))

    def test_uma_000(self) -> None:
        qr = QuantumRegister(3)
        cr = ClassicalRegister(3)
        qc = QuantumCircuit(qr, cr)

        qc.umaj_add(qr[0], qr[1], qr[2])
        qc.measure(qr, cr)
        print(qc)

        # simulate the circuit
        simulator = Aer.get_backend('qasm_simulator')
        result = execute(qc, simulator).result()
        counts = result.get_counts(qc)

        self.assertEqual({"000"}, set(counts.keys()))

    def test_adder_ones(self) -> None:
        nr_qubits = 6
        n = ((nr_qubits - 2) // 2)*2
        adder_circuit = ripple_adder(nr_qubits)
        cin, a, cout = adder_circuit.qregs
        c_cin = ClassicalRegister(1)
        c_a = ClassicalRegister(n)
        c_cout = ClassicalRegister(1)

        qc = QuantumCircuit(cin, a, cout, c_a, c_cin, c_cout)
        # 011110
        # a=11 + b=11 = 110
        qc.x(a)
        qc.extend(adder_circuit)
        qc.measure(cin, c_cin)
        qc.measure(a, c_a)
        qc.measure(cout, c_cout)

        # simulate the circuit
        simulator = Aer.get_backend('qasm_simulator')
        result = execute(qc, simulator).result()
        counts = result.get_counts(qc)

        print(counts)
        # z \xor s_2, ancilla, a1s1a0s0
        self.assertEqual({"1 0 1110"}, set(counts.keys()))

    def test_adder_two(self) -> None:
        nr_qubits = 6
        n = ((nr_qubits - 2) // 2)*2
        adder_circuit = ripple_adder(nr_qubits)
        cin, a, cout = adder_circuit.qregs
        c_cin = ClassicalRegister(1)
        c_a = ClassicalRegister(n)
        c_cout = ClassicalRegister(1)

        qc = QuantumCircuit(cin, a, cout, c_a, c_cin, c_cout)
        # 000110
        # a=01 + b=01 = 010
        qc.x(a[:2])
        qc.extend(adder_circuit)
        qc.measure(cin, c_cin)
        qc.measure(a, c_a)
        qc.measure(cout, c_cout)

        # simulate the circuit
        simulator = Aer.get_backend('qasm_simulator')
        result = execute(qc, simulator).result()
        counts = result.get_counts(qc)

        # z \xor s_2, ancilla, a1s1a0s0
        self.assertEqual({"0 0 0110"}, set(counts.keys()))
