# -*- coding: utf-8 -*-

"""
A module to find the optimizing vibrational coordinates to reduce intermode coupling
"""

import os
import numpy as np
from copy import deepcopy

import rmgpy.constants as constants


from ape.job.job import Job
from ape.qchem import QChemLog
from ape.common import diagonalize_projected_hessian
from ape.InternalCoordinates import get_RedundantCoords, getXYZ

def optvib(coordinate_system, cart_coords, conformer, hessian, linearity, n_vib, rotors, label, path, ncpus, 
                        charge=None, multiplicity=None, level_of_theory=None, basis=None, unrestricted=None):
    n_rotors = len(rotors)
    internal = get_RedundantCoords(symbols, cart_coords)
    grid_of_hessians = get_grid_of_hessians(symbols, cart_coords, conformer, internal, hessian, linearity, n_vib, rotors, label, path, ncpus, 
                                            charge, multiplicity, level_of_theory, basis, unrestricted):
    return

def get_grid_of_hessians(symbols, cart_coords, conformer, internal, hessian, linearity, n_vib, rotors, label, path, ncpus, 
                        charge=None, multiplicity=None, level_of_theory=None, basis=None, unrestricted=None):
    initial_geometry = cart_coords.copy()
    vib_freq, unweighted_v = diagonalize_projected_hessian(conformer, hessian, linearity, n_vib, rotors, label=label)
    grid_of_hessians = {}
    grid_of_hessians[0] = deepcopy(initial_hessian)
    n_rotors = len(rotors)
    for i in range(nmode):
        if i in range(n_rotors): continue
        mode = i + 1
        vector = unweighted_v[i - n_rotors]
        freq = vib_freq[i - n_rotors]
        magnitude = np.linalg.norm(vector)
        reduced_mass = magnitude ** -2 / constants.amu # in amu
        step_size = np.sqrt(constants.hbar / (reduced_mass * constants.amu) / (freq * 2 * np.pi * constants.c * 100)) * 10 ** 10 # in angstrom
        normalizes_vector = vector / magnitude
        qj = np.matmul(internal.B, normalizes_vector)
        qj = qj.reshape(-1,)
        cart_coords += initial_geometry + internal.transform_int_step((qj * step_size).reshape(-1,))
        xyz = getXYZ(symbols, cart_coords)
        file_name = mode
        job = Job(xyz, path, file_name, jobtype='freq', ncpus=ncpus, charge=charge, multiplicity=multiplicity, \
            level_of_theory=level_of_theory, basis=basis, unrestricted=unrestricted)
        # Write Q-Chem input file
        job.write_input_file()
        # Job submission
        job.submit()
        # Parse output file to get the hessian matrix
        output_file_path = os.path.join(path, '{}.q.out'.format(file_name))
        hessian = QChemLog(output_file_path).load_force_constant_matrix()
        if rotors != []:
            hessian = diagonalize_projected_hessian(conformer, hessian, linearity, n_vib, rotors, get_projected_hessian=True, label=label)
        grid_of_hessians[mode] = hessian
    return grid_of_hessians

def rotation_matrix(angle, axis='z', natoms=2):
    """
    Generate matrices for rotation by some angle around a axis.
    """
    s = np.sin(angle)
    c = np.cos(angle)

    i = 'xyz'.index(axis)
    a1 = (i + 1) % 3
    a2 = (i + 2) % 3
    dim = natoms * 3
    R = np.zeros((dim, dim))
    for n in range(natoms):
        R[..., i + 3 * n, i + 3 * n] = 1.
        R[..., a1 + 3 * n, a1 + 3 * n] = c
        R[..., a1 + 3 * n, a2 + 3 * n] = s
        R[..., a2 + 3 * n, a1 + 3 * n] = -s
        R[..., a2 + 3 * n, a2 + 3 * n] = c

    return R

m = rotation_matrix(angle=np.pi/6, axis='z')
print(m)