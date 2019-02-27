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
from qiskit_chemistry.drivers import PySCFDriver, UnitsType
from qiskit_chemistry.aqua_extensions.components.variational_forms import UCCSD


def mol_info_pyscf(atom_coords, basis_set, unit):
    driver = PySCFDriver(atom=atom_coords, unit=unit,
                         charge=0, spin=0, basis=basis_set)
    molecule = driver.run()
    n_particles = molecule.num_alpha + molecule.num_beta
    n_spin_orbitals = molecule.num_orbitals * 2
    return n_particles, n_spin_orbitals


def UCCSD_circuit(n_particles,n_spin_orbitals,
                  map_type='parity', two_qubit_reduction=False,
                  active_occupied=[], active_unoccupied=[],
                  depth=1):

    # Compute number of qubits needed
    n_qubits = n_spin_orbitals
    if map_type == 'parity' and two_qubit_reduction:
        n_qubits -= 2

    # Define the variational form
    var_form = UCCSD(n_qubits, depth=depth,
                     num_orbitals=n_spin_orbitals, num_particles=n_particles,
                     active_occupied=active_occupied, active_unoccupied=active_unoccupied,
                     qubit_mapping=map_type, two_qubit_reduction=two_qubit_reduction,
                     num_time_slices=1)

    # Arbitrary list of values for the variational parameters (does not impact circuit structure)
    params = np.ones(var_form._num_parameters)
    # Return the corresponding quantum circuit, with arbitrary values for parameters
    return var_form.construct_circuit(params)



def UCCSD_H2_sto3g(active_occupied=[], active_unoccupied=[],
                   map_type='parity', two_qubit_reduction=False,
                   depth=1):

    atom_coords = 'H .0 .0 .0; H .0 .0 0.735'
    basis_set = 'sto3g'
    unit = UnitsType.ANGSTROM

    n_particles, n_spin_orbitals = mol_info_pyscf(atom_coords, basis_set, unit)
    circuit = UCCSD_circuit(n_particles, n_spin_orbitals,
                            map_type=map_type, two_qubit_reduction=two_qubit_reduction,
                            active_occupied=active_occupied, active_unoccupied=active_unoccupied,
                            depth=depth)

    return circuit


circuit = UCCSD_H2_sto3g(active_occupied=[], active_unoccupied=[],
               map_type='parity', two_qubit_reduction=False,
               depth=1)

print(circuit.width(), circuit.depth(), circuit.size())
print(circuit.count_ops())


def UCCSD_H2_631g(active_occupied=[], active_unoccupied=[],
                   map_type='parity', two_qubit_reduction=False,
                   depth=1):

    atom_coords = 'H .0 .0 .0; H .0 .0 0.735'
    basis_set = '631g'
    unit = UnitsType.ANGSTROM

    n_particles, n_spin_orbitals = mol_info_pyscf(atom_coords, basis_set, unit)
    circuit = UCCSD_circuit(n_particles, n_spin_orbitals,
                            map_type=map_type, two_qubit_reduction=two_qubit_reduction,
                            active_occupied=active_occupied, active_unoccupied=active_unoccupied,
                            depth=depth)

    return circuit


circuit = UCCSD_H2_631g(active_occupied=[], active_unoccupied=[],
               map_type='parity', two_qubit_reduction=False,
               depth=1)

print(circuit.width(), circuit.depth(), circuit.size())
print(circuit.count_ops())

# def UCCSD_H2_sto3g(active_occupied=[], active_unoccupied=[],
#                     map_type='parity', two_qubit_reduction=False,
#                     depth=1):
#
#     # Set up problem-dependent variables related to molecular system and basis set
#     n_particles = 2
#     n_spin_orbitals = 4
#     n_qubits = n_spin_orbitals
#
#     # Apply qubit_reduction depending on users parameters
#     if map_type == 'parity' and two_qubit_reduction:
#         n_qubits -= 2
#
#     # Define the variational form
#     var_form = UCCSD(n_qubits, depth=depth,
#                      num_orbitals=n_spin_orbitals, num_particles=n_particles,
#                      active_occupied=active_occupied, active_unoccupied=active_unoccupied,
#                      qubit_mapping=map_type,two_qubit_reduction=two_qubit_reduction,
#                      num_time_slices=1)
#
#     # Arbitrary list of values for the variational parameters (does not impact circuit structure)
#     params = np.ones(var_form._num_parameters)
#
#     # Return the corresponding quantum circuit
#     circuit = var_form.construct_circuit(params)
#
#     print(circuit.width())
#     print(circuit.depth())
#     print(circuit.size())
#     print(circuit.count_ops())
#
#     return circuit
#
# UCCSD_H2_sto3g()
# UCCSD_H2_sto3g(active_occupied=[], active_unoccupied=[],
#                 map_type='parity', two_qubit_reduction=True)
#
#
# def UCCSD_H2_631g(active_occupied=[], active_unoccupied=[],
#                     map_type='parity', two_qubit_reduction=False,
#                     depth=1):
#
#     # Set up problem-dependent variables related to molecular system and basis set
#     n_particles = 2
#     n_spin_orbitals = 8
#     n_qubits = n_spin_orbitals
#
#     # Apply qubit_reduction depending on users parameters
#     if map_type == 'parity' and two_qubit_reduction:
#         n_qubits -= 2
#
#     # Define the variational form
#     var_form = UCCSD(n_qubits, depth=depth,
#                      num_orbitals=n_spin_orbitals, num_particles=n_particles,
#                      active_occupied=active_occupied, active_unoccupied=active_unoccupied,
#                      qubit_mapping=map_type,two_qubit_reduction=two_qubit_reduction,
#                      num_time_slices=1)
#
#     # Arbitrary list of values for the variational parameters (does not impact circuit structure)
#     params = np.ones(var_form._num_parameters)
#
#     # Return the corresponding quantum circuit
#     circuit = var_form.construct_circuit(params)
#
#     print(circuit.width())
#     print(circuit.depth())
#     print(circuit.size())
#     print(circuit.count_ops())
#
#     return circuit
#
# UCCSD_H2_631g()
# UCCSD_H2_631g(active_occupied=[], active_unoccupied=[],
#                 map_type='parity', two_qubit_reduction=True)
#
#
# def UCCSD_LiH_sto3g(active_occupied=[], active_unoccupied=[],
#                     map_type='parity', two_qubit_reduction=False,
#                     depth=1):
#
#     # Set up problem-dependent variables related to molecular system and basis set
#     n_particles = 4
#     n_spin_orbitals = 12
#     n_qubits = n_spin_orbitals
#
#     # Apply qubit_reduction depending on users parameters
#     if map_type == 'parity' and two_qubit_reduction:
#         n_qubits -= 2
#
#     # Define the variational form
#     var_form = UCCSD(n_qubits, depth=depth,
#                      num_orbitals=n_spin_orbitals, num_particles=n_particles,
#                      active_occupied=active_occupied, active_unoccupied=active_unoccupied,
#                      qubit_mapping=map_type,two_qubit_reduction=two_qubit_reduction,
#                      num_time_slices=1)
#
#     # Arbitrary list of values for the variational parameters (does not impact circuit structure)
#     params = np.ones(var_form._num_parameters)
#
#     # Return the corresponding quantum circuit
#     circuit = var_form.construct_circuit(params)
#
#     print(circuit.width())
#     print(circuit.depth())
#     print(circuit.size())
#     print(circuit.count_ops())
#
#     return circuit
#
# UCCSD_LiH_sto3g()
# UCCSD_LiH_sto3g(active_occupied=[], active_unoccupied=[],
#                 map_type='parity', two_qubit_reduction=True)
# UCCSD_LiH_sto3g(active_occupied=[0], active_unoccupied=[0, 1],
#                 map_type='parity', two_qubit_reduction=True)
#
#
# def UCCSD_NaH_sto3g(active_occupied=[], active_unoccupied=[],
#                     map_type='parity', two_qubit_reduction=False,
#                     depth=1):
#
#     # Set up problem-dependent variables related to molecular system and basis set
#     n_particles = 12
#     n_spin_orbitals = 20
#     n_qubits = n_spin_orbitals
#
#     # Apply qubit_reduction depending on users parameters
#     if map_type == 'parity' and two_qubit_reduction:
#         n_qubits -= 2
#
#     # Define the variational form
#     var_form = UCCSD(n_qubits, depth=depth,
#                      num_orbitals=n_spin_orbitals, num_particles=n_particles,
#                      active_occupied=active_occupied, active_unoccupied=active_unoccupied,
#                      qubit_mapping=map_type,two_qubit_reduction=two_qubit_reduction,
#                      num_time_slices=1)
#
#     # Arbitrary list of values for the variational parameters (does not impact circuit structure)
#     params = np.ones(var_form._num_parameters)
#
#     # Return the corresponding quantum circuit
#     circuit = var_form.construct_circuit(params)
#
#     print(circuit.width())
#     print(circuit.depth())
#     print(circuit.size())
#     print(circuit.count_ops())
#
#     return circuit
#
# UCCSD_NaH_sto3g()
# UCCSD_NaH_sto3g(active_occupied=[], active_unoccupied=[],
#                 map_type='parity', two_qubit_reduction=True)

