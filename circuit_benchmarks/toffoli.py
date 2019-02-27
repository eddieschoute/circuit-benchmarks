import itertools
import math

from qiskit import QuantumRegister, QuantumCircuit
from qiskit.circuit import Gate, InstructionSet
from qiskit.dagcircuit import DAGCircuit
from qiskit.extensions.standard import *
from qiskit.qasm import pi


def toffoli(number_qubits: int):
    assert number_qubits >= 2
    controls = QuantumRegister(number_qubits-1)
    target = QuantumRegister(1)
    qc = QuantumCircuit(controls, target, name="toffoli")
    qc.ncrx(pi/2, target[0], *controls)
    return qc

class NcrxGate(Gate):
    """n-controlled x rotation gate."""

    def __init__(self, theta, tgt, *ctls, circ=None):
        """Create new Toffoli gate."""
        assert len(ctls) >= 1
        super().__init__(f"c^{len(ctls)}rx", [theta], list(ctls) + [tgt], circ)

    def _define_decompositions(self):
        decomposition = DAGCircuit()
        q = QuantumRegister(3, "q")
        last_control = q[len(q)-2]
        target = q[len(q)-1]
        decomposition.add_qreg(q)
        decomposition.add_basis_element("cu3", 2, 0, 3)
        if len(self.qargs) == 2:
            # Equal to crx
            crx_theta = Cu3Gate(self.param[0], -pi/2, pi/2, last_control, target)
            decomposition.apply_operation_back(crx_theta)
        else:
            # Recurse
            nr_controls = len(q) - 1
            decomposition.add_basis_element(f"c^{nr_controls-1}rx", nr_controls, 0, 1)
            rule = [
                # C-sqrt(rx(theta)) gate
                Cu3Gate(self.param[0]/2, -pi/2, pi/2, last_control, target),
                NcrxGate(pi/2, q[0:nr_controls-1], last_control), # toffoli
                Cu3Gate(self.param[0]/2, -pi/2, pi/2, last_control, target).inverse(),
                NcrxGate(pi/2, q[0:nr_controls-1], last_control), # toffoli
                NcrxGate(self.param[0]/2, q[0:nr_controls-1], q[nr_controls]) # c^nrx(theta/2) gate on n-1 qubits
            ]
            for inst in rule:
                decomposition.apply_operation_back(inst)
        self._decompositions = [decomposition]

    def inverse(self):
        """Invert this gate."""
        return self  # self-inverse

    def reapply(self, circ):
        """Reapply this gate to corresponding qubits in circ."""
        self._modifiers(circ.ntoffoli(self.qargs[0:len(self.qargs)-1], self.qargs[len(self.qargs-1)]))


def ncrx(self, theta, tgt, *ctls):
    """Apply Toffoli to from ctl1 and ctl2 to tgt."""
    if all(isinstance(ctl, QuantumRegister) for ctl in ctls) and \
            isinstance(tgt, QuantumRegister) and \
            all(len(ctl) == len(tgt) for ctl in ctls):
        instructions = InstructionSet()
        for i in range(ctls[0].size):
            instructions.add(self.ntoffoli(theta, zip(ctls, itertools.repeat(i)), (tgt, i)))
        return instructions

    for ctl in ctls:
        self._check_qubit(ctl)
    self._check_qubit(tgt)
    self._check_dups(list(ctls) + [tgt])
    return self._attach(NcrxGate(theta, tgt, *ctls, circ=self))

QuantumCircuit.ncrx = ncrx
