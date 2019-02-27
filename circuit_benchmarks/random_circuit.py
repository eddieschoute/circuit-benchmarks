# -*- coding: utf-8 -*-

# Copyright 2017 IBM RESEARCH. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# =============================================================================

import math

import numpy as np
from qiskit import QuantumCircuit, QuantumRegister
from scipy.linalg import qr
from qiskit.mapper import two_qubit_kak

def random_unitary(n: int, np_random: np.random.RandomState):
    """Return an n x n Haar distributed unitary matrix.
    Return numpy array.
    """
    X = (1.0 / math.sqrt(2.0)) * (np_random.randn(n, n) +
                                  1j * np_random.randn(n, n))
    Q, R = qr(X)
    R = np.diag(np.diag(R) / np.abs(np.diag(R)))
    U = np.dot(Q, R)
    return U


def random_circuit(n: int, depth: int = 20, seed: int = None) -> QuantumCircuit:
    """Create a quantum program containing model circuits.
    The model circuits consist of layers of Haar random
    elements of SU(4) applied between corresponding pairs
    of qubits in a random bipartition.
    name = leading name of circuits
    n = number of qubits
    depth = ideal depth of each model circuit (over SU(4))
    Return a quantum circuit.
    """
    np_random = np.random.RandomState(seed)
    q = QuantumRegister(n)
    qc = QuantumCircuit(q)
    # For each layer
    for j in range(depth):
        # Generate uniformly random permutation Pj of [0...n-1]
        perm = np_random.permutation(n)
        # For each pair p in Pj, generate Haar random SU(4)
        for k in range(math.floor(n / 2)):
            qubits = [int(perm[2 * k]), int(perm[2 * k + 1])]
            U = random_unitary(4, np_random)
            # TODO: Replace this with a general "two-qubit" gate instead of decomposing it.
            for gate in two_qubit_kak(U):
                i0 = qubits[gate["args"][0]]
                if gate["name"] == "cx":
                    i1 = qubits[gate["args"][1]]
                    qc.cx(q[i0], q[i1])
                elif gate["name"] == "u1":
                    qc.u1(gate["params"][2], q[i0])
                elif gate["name"] == "u2":
                    qc.u2(gate["params"][1], gate["params"][2],
                          q[i0])
                elif gate["name"] == "u3":
                    qc.u3(gate["params"][0], gate["params"][1],
                          gate["params"][2], q[i0])
                elif gate["name"] == "id":
                    pass  # do nothing
    return qc

