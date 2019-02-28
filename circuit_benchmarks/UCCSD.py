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

"""
    As part of a set of benchmark circuits, this file provides quantum circuits
    used in variational approaches in quantum chemistry, such as VQE.

    It provide the quantum circuits for different molecules and basis sets,
    using the UCCSD ansatz (variational form).

    It also provides a higher-level interface allowing user to define new benchmark circuits
    for their choice of molecule and basis set.

    One could then study how different metrics associated to the circuit (width, depth, runtime...)
    scale with problem size of basis set. The study could also be extended across different variational forms.
"""

import numpy as np
from qiskit_chemistry.drivers import PySCFDriver, UnitsType
from qiskit_chemistry.aqua_extensions.components.variational_forms import UCCSD


def mol_info_pyscf(atom_coords, basis_set, unit):
    """
        Derives the number of particles and spin orbitals of a molecular system
        using the pySCF driver

        Args:
            atom_coords : atom coordinates (string)
            basis_set : basis set (string)
            unit : unit of distance between atoms (qiskit_chemistry.drivers.UnitsType)
        Return:
            n_particles (int)
            n_spin_orbitals (int)
    """

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
    """
        Derives the number of particles and spin orbitals of a molecular system
        using the pySCF driver

        Args:
            n_particles (int)
            n_spin_orbitals (int)
            map_type : the mapping type for symmetry (string)
            two_qubit_reduction : flag to reduce qubit count by 2 (bool)
            active_occupied : active occupied orbitals (list of int)
            active_unoccupied : active unoccupied orbitals (list on int)
            depth : allow the ansatz to reach states further away from initial state by
                    duplicating and concatenating the quantum circuit with itself (int)
        Return:
            n_particles (int)
            n_spin_orbitals (int)
    """

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



def UCCSD_qc(molecule, active_occupied=[], active_unoccupied=[],
            map_type='parity', two_qubit_reduction=False,
            depth=1):
    """
        Generates a UCCSD circuit for H2 in sto-3g basis
        By default, no active orbitals are specified, and there is no 2-qubit reduction
    """

    atom_coords, basis_set, unit = molecule[:]

    n_particles, n_spin_orbitals = mol_info_pyscf(atom_coords, basis_set, unit)
    circuit = UCCSD_circuit(n_particles, n_spin_orbitals,
                            map_type=map_type, two_qubit_reduction=two_qubit_reduction,
                            active_occupied=active_occupied, active_unoccupied=active_unoccupied,
                            depth=depth)
    return circuit

# Some molecules
H2 = 'H .0 .0 .0; H .0 .0 0.735'
LiH = 'Li .0 .0 .0; H .0 .0 1.6'
NaH = 'Na .0 .0 .0; H .0 .0 1.9'
molecules = [H2, LiH, NaH]

# Some bases
bases = ['sto3g', '631g', 'ccpVDZ']

# How the resource required by a quantum circuit scale with molecule size, for a given basis
for m in molecules:
    print("Molecule ---- ")
    qc = UCCSD_qc((m, 'sto3g', UnitsType.ANGSTROM), active_occupied=[], active_unoccupied=[],
                   map_type='parity', two_qubit_reduction=False,
                   depth=1)
    print(qc.width(), qc.depth(), qc.size())
    print(qc.count_ops())

# How the resource required by a quantum circuit scale with basis size, for a given molecule
for b in bases:
    print("Quantum circuit for H2 in ", b, " basis")
    qc = UCCSD_qc((H2, b, UnitsType.ANGSTROM), active_occupied=[], active_unoccupied=[],
                   map_type='parity', two_qubit_reduction=False,
                   depth=1)
    print(qc.width(), qc.depth(), qc.size())
    print(qc.count_ops())


def UCCSD_H2_sto3g(active_occupied=[], active_unoccupied=[],
                   map_type='parity', two_qubit_reduction=False,
                   depth=1):
    """
        Generates a UCCSD circuit for H2 in sto-3g basis
        By default, no active orbitals are specified, and there is no 2-qubit reduction
    """

    atom_coords = 'H .0 .0 .0; H .0 .0 0.735'
    basis_set = 'sto3g'
    unit = UnitsType.ANGSTROM

    n_particles, n_spin_orbitals = mol_info_pyscf(atom_coords, basis_set, unit)
    circuit = UCCSD_circuit(n_particles, n_spin_orbitals,
                            map_type=map_type, two_qubit_reduction=two_qubit_reduction,
                            active_occupied=active_occupied, active_unoccupied=active_unoccupied,
                            depth=depth)
    return circuit


def UCCSD_H2_631g(active_occupied=[], active_unoccupied=[],
                  map_type='parity', two_qubit_reduction=False,
                  depth=1):
    """
        Generates a UCCSD circuit for H2 in 6-31g basis
        By default, no active orbitals are specified, there is no 2-qubit reduction, depth is set to 1
    """

    atom_coords = 'H .0 .0 .0; H .0 .0 0.735'
    basis_set = '631g'
    unit = UnitsType.ANGSTROM

    n_particles, n_spin_orbitals = mol_info_pyscf(atom_coords, basis_set, unit)
    circuit = UCCSD_circuit(n_particles, n_spin_orbitals,
                            map_type=map_type, two_qubit_reduction=two_qubit_reduction,
                            active_occupied=active_occupied, active_unoccupied=active_unoccupied,
                            depth=depth)
    return circuit


def UCCSD_H2_ccpVDZ(active_occupied=[], active_unoccupied=[],
                  map_type='parity', two_qubit_reduction=False,
                  depth=1):
    """
        Generates a UCCSD circuit for H2 in 6-31g basis
        By default, no active orbitals are specified, there is no 2-qubit reduction, depth is set to 1
    """

    atom_coords = 'H .0 .0 .0; H .0 .0 0.735'
    basis_set = 'ccpVDZ'
    unit = UnitsType.ANGSTROM

    n_particles, n_spin_orbitals = mol_info_pyscf(atom_coords, basis_set, unit)
    circuit = UCCSD_circuit(n_particles, n_spin_orbitals,
                            map_type=map_type, two_qubit_reduction=two_qubit_reduction,
                            active_occupied=active_occupied, active_unoccupied=active_unoccupied,
                            depth=depth)
    return circuit


def UCCSD_LiH_sto3g(active_occupied=[], active_unoccupied=[],
                    map_type='parity', two_qubit_reduction=False,
                    depth=1):
    """
        Generates a UCCSD circuit for LiH in sto-3g basis
        By default, no active orbitals are specified, there is no 2-qubit reduction, depth is set to 1
    """

    atom_coords = 'Li .0 .0 .0; H .0 .0 1.6'
    basis_set = 'sto-3g'
    unit = UnitsType.ANGSTROM

    n_particles, n_spin_orbitals = mol_info_pyscf(atom_coords, basis_set, unit)
    circuit = UCCSD_circuit(n_particles, n_spin_orbitals,
                            map_type=map_type, two_qubit_reduction=two_qubit_reduction,
                            active_occupied=active_occupied, active_unoccupied=active_unoccupied,
                            depth=depth)
    return circuit


def UCCSD_NaH_sto3g(active_occupied=[], active_unoccupied=[],
                    map_type='parity', two_qubit_reduction=False,
                    depth=1):
    """
        Generates a UCCSD circuit for NaH in sto-3g basis
        By default, no active orbitals are specified, there is no 2-qubit reduction, depth is set to 1
    """

    atom_coords = 'Na .0 .0 .0; H .0 .0 1.9'
    basis_set = 'sto-3g'
    unit = UnitsType.ANGSTROM

    n_particles, n_spin_orbitals = mol_info_pyscf(atom_coords, basis_set, unit)
    circuit = UCCSD_circuit(n_particles, n_spin_orbitals,
                            map_type=map_type, two_qubit_reduction=two_qubit_reduction,
                            active_occupied=active_occupied, active_unoccupied=active_unoccupied,
                            depth=depth)
    return circuit


# qc = UCCSD_H2_sto3g()
# print(qc.width(), qc.depth(), qc.size())
# print(qc.count_ops())
#
# qc = UCCSD_H2_631g()
# print(qc.width(), qc.depth(), qc.size())
# print(qc.count_ops())
#
# qc = UCCSD_H2_ccpVDZ()
# print(qc.width(), qc.depth(), qc.size())
# print(qc.count_ops())