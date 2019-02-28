import itertools
import math

from qiskit import QuantumRegister, QuantumCircuit, ClassicalRegister
from qiskit.circuit import Gate, InstructionSet
from qiskit.dagcircuit import DAGCircuit
from qiskit.extensions.standard import *
from qiskit.qasm import pi


def toffoli(number_qubits: int):
    assert number_qubits >= 2
    q = QuantumRegister(number_qubits)
    qc = QuantumCircuit(q, name="toffoli")
    # for i in range(number_qubits-1):
    #     qc.h(controls[i])
    qc.ntoffoli(q[number_qubits-1], *q[0:number_qubits-1])
    # qc.measure(controls, c_controls)
    # qc.measure(target, c_target)
    return qc

class NcrxGate(Gate):
    """n-controlled x rotation gate."""

    def __init__(self, theta, tgt, *ctls, circ=None):
        """Create new Toffoli gate."""
        assert len(ctls) >= 1
        super().__init__(f"c^{len(ctls)}rx", [theta], [tgt] + list(ctls), circ)

    def _define_decompositions(self):
        decomposition = DAGCircuit()
        nr_qubits = len(self.qargs)
        q = QuantumRegister(nr_qubits)
        last_control = q[1]
        target = q[0]
        decomposition.add_qreg(q)
        if nr_qubits == 2:
            # Equal to crx of theta
            crx_theta = Cu3Gate(self.params[0], -pi/2, pi/2, last_control, target)
            decomposition.apply_operation_back(crx_theta)
        else:
            # Recurse
            rule = [
                # C-sqrt(rx(theta)) gate
                Cu3Gate(self.params[0]/2, -pi/2, pi/2, last_control, target),
                NcrxGate(pi, last_control, *q[2:]), # toffoli
                Cu3Gate(self.params[0]/2, -pi/2, pi/2, last_control, target).inverse(),
                NcrxGate(pi, last_control, *q[2:]), # toffoli
                NcrxGate(self.params[0]/2, target, *q[2:]) # c^nrx(theta/2) gate on n-1 qubits
            ]
            for inst in rule:
                decomposition.apply_operation_back(inst)
            # decomposition.apply_operation_back(ToffoliGate(q[1], q[2], q[0]))
        self._decompositions = [decomposition]

    def inverse(self):
        """Invert this gate."""
        return self  # self-inverse

    def reapply(self, circ):
        """Reapply this gate to corresponding qubits in circ."""
        self._modifiers(circ.ncrx(self.params[0], self.qargs[0], *self.qargs[1:]))

def ncrx(self, theta, tgt, *ctls):
    """Apply n-controlled x-rotation(theta) to target from controls"""
    if all(isinstance(ctl, QuantumRegister) for ctl in ctls) and \
            isinstance(tgt, QuantumRegister) and \
            all(len(ctl) == len(tgt) for ctl in ctls):
        instructions = InstructionSet()
        for i in range(ctls[0].size):
            instructions.add(self.ntoffoli(theta, (tgt, i), *zip(ctls, itertools.repeat(i))))
        return instructions

    for ctl in ctls:
        self._check_qubit(ctl)
    self._check_qubit(tgt)
    self._check_dups(list(ctls) + [tgt])
    return self._attach(NcrxGate(theta, tgt, *ctls, circ=self))

def ntoffoli(self, tgt, *ctls):
    """Apply n-controlled Toffoli to tgt with controls."""
    if all(isinstance(ctl, QuantumRegister) for ctl in ctls) and \
            isinstance(tgt, QuantumRegister) and \
            all(len(ctl) == len(tgt) for ctl in ctls):
        instructions = InstructionSet()
        for i in range(ctls[0].size):
            instructions.add(self.ntoffoli((tgt, i), *zip(ctls, itertools.repeat(i))))
        return instructions

    for ctl in ctls:
        self._check_qubit(ctl)
    self._check_qubit(tgt)
    self._check_dups(list(ctls) + [tgt])
    return self._attach(NcrxGate(pi, tgt, *ctls, circ=self))

QuantumCircuit.ncrx = ncrx
QuantumCircuit.ntoffoli = ntoffoli
