# -*- coding: utf-8 -*-

"""
A module to find the optimizing vibrational coordinates to reduce intermode coupling
"""

import os
import logging
import numpy as np
from copy import deepcopy
from scipy import optimize

import rmgpy.constants as constants

from ape.job.job import Job
from ape.qchem import QChemLog
from ape.exceptions import InputError
from ape.common import diagonalize_projected_hessian
from ape.InternalCoordinates import get_RedundantCoords, getXYZ

class OptVib(object):
    def __init__(self, symbols, nmode, coordinate_system, cart_coords, conformer, hessian, linearity, n_vib, rotors, label, path, ncpus, 
                 charge=None, multiplicity=None, level_of_theory=None, basis=None, unrestricted=None, gen_basis="", purecart=None):
        self.symbols = symbols
        self.nmode = nmode
        self.coordinate_system = coordinate_system
        self.cart_coords = cart_coords
        self.conformer = conformer
        self.hessian = hessian
        self.linearity = linearity
        self.n_vib = n_vib
        self.rotors = rotors
        self.label = label
        self.path = path
        self.ncpus = ncpus
        self.charge = charge
        self.multiplicity = multiplicity
        self.level_of_theory = level_of_theory
        self.basis = basis
        self.unrestricted = unrestricted
        self.gen_basis = gen_basis
        self.purecart = purecart
        self.n_rotors = len(self.rotors)
        self.internal_object = get_RedundantCoords(self.symbols, self.cart_coords)
        self.natoms = len(self.symbols)


    def get_optvib(self):
        """
        Algorithms for local and optimal vibrations.
        """
        self.grid_of_hessians = self.get_grid_of_hessians()
        self.weighted_v = diagonalize_projected_hessian(self.conformer, self.hessian, self.linearity, self.n_vib, 
                                                        self.rotors, get_weighted_vectors=True, label=self.label)
        _, unweighted_v = diagonalize_projected_hessian(self.conformer, self.hessian, self.linearity, self.n_vib,
                                                        self.rotors, label=self.label)
        
        b1 = (0, 2*np.pi)
        b2 = (0, 2*np.pi)
        b3 = (0, 2*np.pi)
        bounds = (b1 ,b2 ,b3)
        result = optimize.differential_evolution(self.objectiveFunction, bounds=bounds)

        if not result.success:
            logging.warning('Optimization to find optimal vibrational coordinates fails.')
        else:
            angles = result.x
            logging.info('Optimization vibrational coordinates are found by rotating pairs of eigenvectors around x-axis {:.2f} degree, around y-axis {:.2f} degree and around z-axis {:.2f} degree.'.format(angles[0] / np.pi * 180, angles[1] / np.pi * 180, angles[2] / np.pi * 180))
            U = self.rotation_matrix(angles, natoms=self.natoms)
            unweighted_v = unweighted_v.dot(U.T)
            V = self.weighted_v.dot(U.T)
            H = V.dot(self.grid_of_hessians[0]).dot(V.T)
            vib_freq = []
            for i in range(self.n_vib):
                freq = np.sqrt(H[i][i]) / (2 * np.pi * constants.c * 100)
                vib_freq.append(freq)

        return vib_freq, unweighted_v

    def get_grid_of_hessians(self):
        """
        Hessians are generated on a grid of one point per vibrational mode.
        """
        initial_geometry = self.cart_coords.copy()
        vib_freq, unweighted_v = diagonalize_projected_hessian(self.conformer, self.hessian, self.linearity, self.n_vib, self.rotors, label=self.label)
        grid_of_hessians = {}
        fm = diagonalize_projected_hessian(self.conformer, self.hessian, self.linearity, self.n_vib, self.rotors, get_mass_weighted_hessian=True, label=self.label)
        grid_of_hessians[0] = deepcopy(fm)
        for i in range(self.n_vib):
            internal = deepcopy(self.internal_object)
            mode = i + 1
            vector = unweighted_v[i - self.n_rotors]
            freq = vib_freq[i - self.n_rotors]
            magnitude = np.linalg.norm(vector)
            reduced_mass = magnitude ** -2 / constants.amu # in amu
            step_size = np.sqrt(constants.hbar / (reduced_mass * constants.amu) / (freq * 2 * np.pi * constants.c * 100)) * 10 ** 10 # in angstrom
            normalizes_vector = vector / magnitude
            qj = np.matmul(internal.B, normalizes_vector)
            qj = qj.reshape(-1,)
            cart_coords = initial_geometry + internal.transform_int_step((qj * step_size).reshape(-1,))
            xyz = getXYZ(self.symbols, cart_coords)
            file_name = mode
            job = Job(xyz, self.path, file_name, jobtype='freq', ncpus=self.ncpus, charge=self.charge, multiplicity=self.multiplicity,
                      level_of_theory=self.level_of_theory, basis=self.basis, unrestricted=self.unrestricted, gen_basis=self.gen_basis,
                      purecart=self.purecart)
            # Write Q-Chem input file
            job.write_input_file()
            # Job submission
            job.submit()
            # Parse output file to get the hessian matrix
            output_file_path = os.path.join(self.path, '{}.q.out'.format(file_name))
            hessian = QChemLog(output_file_path).load_force_constant_matrix()
            fm = diagonalize_projected_hessian(self.conformer, hessian, self.linearity, self.n_vib, self.rotors, get_mass_weighted_hessian=True, label=self.label)
            grid_of_hessians[mode] = fm
        return grid_of_hessians

    def objectiveFunction(self, angles):
        """
        To produce optimal coordinates, metrics which quantify off-diagonal couplings
        over a grid of Hessian matrices are minimized through unitary rotations of
        the vibrational basis.
        
        1. coordinate_system == "E-Optimized"
            Return the sum of squared off-diagonal coupling
        2. coordinate_system == "E'-Optimized"
            Return the suml squared change in off-diagonal coupling
        """
        U = self.rotation_matrix(angles, natoms=self.natoms)
        V = self.weighted_v.dot(U.T)
        E = 0
        if self.coordinate_system == "E-Optimized":
            for key in self.grid_of_hessians.keys():
                dim = self.n_vib
                hessian = self.grid_of_hessians[key]
                H = V.dot(hessian).dot(V.T) / ((2 * np.pi * constants.c * 100) ** 2)
                for i in range(dim):
                    for j in range(dim):
                        if i < j:
                            E += (H[i][j]) ** 2
        elif self.coordinate_system == "E'-Optimized":
            H0 = V.dot(self.grid_of_hessians[0]).dot(V.T) / ((2 * np.pi * constants.c * 100) ** 2)
            for key in self.grid_of_hessians.keys():
                if key == 0: continue
                dim = self.n_vib
                hessian = self.grid_of_hessians[key]
                H = V.dot(hessian).dot(V.T) / ((2 * np.pi * constants.c * 100) ** 2)
                for i in range(dim):
                    for j in range(dim):
                        if i < j:
                            E += (H[i][j] - H0[i][j]) ** 2
        else:
            raise InputError("The value of coordinate_system should be E-Optimized or E'-Optimized to produce optimal coordinates.")
        
        return E

    def rotation_matrix(self, angles, axis=None, natoms=1):
        """
        Calculates and returns the rotation matrix defined by three angles of
        rotation about the x, y, and z axes or one angle of rotation about a
        given axis.
        """
        if axis is None:
            Rx = np.array(
                [[1.0, 0.0, 0.0],
                [0.0, np.cos(angles[0]), -np.sin(angles[0])],
                [0.0, np.sin(angles[0]), np.cos(angles[0])]]
            )
            Ry = np.array(
                [[np.cos(angles[1]), 0.0, np.sin(angles[1])],
                [0.0, 1.0, 0.0],
                [-np.sin(angles[1]), 0.0, np.cos(angles[1])]]
            )
            Rz = np.array(
                [[np.cos(angles[2]), -np.sin(angles[2]), 0.0],
                [np.sin(angles[2]), np.cos(angles[2]), 0.0],
                [0.0, 0.0, 1.0]]
            )
            R =  Rx.dot(Ry).dot(Rz)
        else:
            axis = axis/np.sqrt(axis.dot(axis))
            x, y, z = axis[0], axis[1], axis[2]
            angle = angles
            R = np.array(
                [[np.cos(angle) + x ** 2 * (1 - np.cos(angle)),
                x * y * (1 - np.cos(angle)) - z * np.sin(angle),
                x * z * (1 - np.cos(angle))+y * np.sin(angle)],
                [y * x * (1 - np.cos(angle))+z * np.sin(angle),
                np.cos(angle) + y ** 2 * (1 - np.cos(angle)),
                y * z * (1 - np.cos(angle)) - x * np.sin(angle)],
                [z * x * (1 - np.cos(angle)) - y * np.sin(angle),
                z * y * (1 - np.cos(angle)) + x * np.sin(angle),
                np.cos(angle) + z ** 2 * (1 - np.cos(angle))]]
            )
        
        dim = natoms * 3
        R_mat = np.zeros((dim, dim))
        for n in range(natoms):
            R_mat[0 + 3 * n, 0 + 3 * n] = R[0 ,0]
            R_mat[0 + 3 * n, 1 + 3 * n] = R[0 ,1]
            R_mat[0 + 3 * n, 2 + 3 * n] = R[0 ,2]
            R_mat[1 + 3 * n, 0 + 3 * n] = R[1 ,0]
            R_mat[1 + 3 * n, 1 + 3 * n] = R[1 ,1]
            R_mat[1 + 3 * n, 2 + 3 * n] = R[1 ,2]
            R_mat[2 + 3 * n, 0 + 3 * n] = R[2 ,0]
            R_mat[2 + 3 * n, 1 + 3 * n] = R[2 ,1]
            R_mat[2 + 3 * n, 2 + 3 * n] = R[2 ,2]

        return R_mat