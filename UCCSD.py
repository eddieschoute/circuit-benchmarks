# -*- coding: utf-8 -*-

# Copyright 2019, IBM.
#
# This source code is licensed under the Apache License, Version 2.0 found in
# the LICENSE.txt file in the root directory of this source tree.
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

import numpy as np
from qiskit_chemistry.aqua_extensions.components.variational_forms import UCCSD


def UCCSD_LiH_sto3g(active_occupied=[], active_unoccupied=[],
                    map_type='parity', two_qubit_reduction=False,
                    depth=1):

    # Set up problem-dependent variables related to molecular system and basis set
    n_orbitals = 12
    n_particles = 4
    n_qubits = 12

    if map_type == 'parity' and two_qubit_reduction:
        n_qubits -= 2

    # Define the variational form
    var_form = UCCSD(n_qubits, depth=depth,
                     num_orbitals=n_orbitals, num_particles=n_particles,
                     active_occupied=active_occupied, active_unoccupied=active_unoccupied,
                     qubit_mapping=map_type,two_qubit_reduction=two_qubit_reduction,
                     num_time_slices=1)

    # Arbitrary list of values for the variational parameters (does not impact circuit structure)
    params = np.ones(var_form._num_parameters)

    # Return the corresponding quantum circuit
    circuit = var_form.construct_circuit(params)

    print(circuit.width())
    print(circuit.depth())
    print(circuit.size())
    print(circuit.count_ops())

    return circuit

UCCSD_LiH_sto3g()
UCCSD_LiH_sto3g(active_occupied=[], active_unoccupied=[],
                map_type='parity', two_qubit_reduction=True)
UCCSD_LiH_sto3g(active_occupied=[0], active_unoccupied=[0, 1],
                map_type='parity', two_qubit_reduction=True)