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
from ape.exceptions import InputError, ConvergeError
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

    def get_optvib(self):
        """
        Algorithms for local and optimal vibrations.
        """
        logging.info('{0} modes of {1} finding...'.format(self.coordinate_system, label))
        self.grid_of_hessians = self.get_grid_of_hessians()
        self.mwv = diagonalize_projected_hessian(self.conformer, self.hessian, self.linearity, self.n_vib, 
                                                 self.rotors, get_weighted_vectors=True, label=self.label)

        # random rotations by 1° over all pairs of normal modes to break symmetry
        num = int(self.n_vib * (self.n_vib - 1) / 2) # 2-combination of self.n_vib
        angles = np.zeros(num)
        ind = 0
        for i in range(self.n_vib):
            new_raw = True
            for j in range(self.n_vib):
                if i < j:
                    if (i + 1) % 2 == 1 and new_raw:
                        angles[ind] = 1 / 180 * np.pi
                        new_raw = False
                    else:
                        angles[ind] = 0
                    ind += 1
        Ui = self.U(angles)
        self.mwv = self.mwv.T.dot(Ui)

        # Do Jacobi sweeps over all pairs of modes 
        self.Jacobi_sweeps()

        # Calculate anharmonic frequencies
        mwv = self.mwv
        H = mwv.T.dot(self.grid_of_hessians[0]).dot(mwv)
        vib_freq = []
        for i in range(self.n_vib):
            freq = np.sqrt(H[i][i]) / (2 * np.pi * constants.c * 100)
            vib_freq.append(freq)
        
        # Calculate optimal coordinates in terms of not mass-weighted cartesian coordinate
        mass = self.conformer.mass.value_si
        mass_3N_array = np.array([i for i in mass for j in range(3)])
        mass_mat = np.diag(mass_3N_array)
        inv_sq_mass_mat = np.linalg.inv(mass_mat ** 0.5)
        unweighted_v = np.matmul(inv_sq_mass_mat, mwv).T

        # Sort anharmonic frequencies in ascending order
        unweighted_v = [v for _, v in sorted(zip(vib_freq, unweighted_v))]
        vib_freq = sorted(vib_freq)

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

    def objectiveFunction(self, angle):
        """
        To produce optimal coordinates, metrics which quantify off-diagonal couplings
        over a grid of Hessian matrices are minimized through unitary rotations of
        the vibrational basis.
        
        1. coordinate_system == "E-Optimized"
            Return the sum of squared off-diagonal coupling
        2. coordinate_system == "E'-Optimized"
            Return the sum squared change in off-diagonal coupling
        3. coordinate_system == "Pipek_Mezey"
            Return the sum of squares of the atomic contribution to the modes
        """
        angles = self.angles.copy()
        angles[self.n] = angle

        # Rotate the paires of eigenvectors
        U = self.U(angles)
        V = self.mwv.dot(U)
        E = 0
        if self.coordinate_system == "E-Optimized":
            for key in self.grid_of_hessians.keys():
                dim = self.n_vib
                hessian = self.grid_of_hessians[key]
                H = V.T.dot(hessian).dot(V) / ((2 * np.pi * constants.c * 100) ** 2)
                for i in range(dim):
                    for j in range(dim):
                        if i < j:
                            E += (H[i][j]) ** 2
        elif self.coordinate_system == "E'-Optimized":
            H0 = V.T.dot(self.grid_of_hessians[0]).dot(V) / ((2 * np.pi * constants.c * 100) ** 2)
            for key in self.grid_of_hessians.keys():
                if key == 0: continue
                dim = self.n_vib
                hessian = self.grid_of_hessians[key]
                H = V.T.dot(hessian).dot(V) / ((2 * np.pi * constants.c * 100) ** 2)
                for i in range(dim):
                    for j in range(dim):
                        if i < j:
                            E += (H[i][j] - H0[i][j]) ** 2
        elif self.coordinate_system == "Pipek_Mezey":
            modes = V.T
            squared_modes = modes ** 2
            c = squared_modes[:, 0::3] + squared_modes[:, 1::3] + squared_modes[:, 2::3]
            E = -np.linalg.norm(c) ** 2            
        else:
            raise InputError("The value of coordinate_system should be E-Optimized or E'-Optimized to produce optimal coordinates.")

        return E

    def Jacobi_sweeps(self, thresh=1e-6, thresh2=1e-4, printing=True):
        """
        Jacobi sweeps are performed over the angles until minimization was reached.
        """
        num = int(self.n_vib * (self.n_vib - 1) / 2) # 2-combination of a set self.n_vib
        err  = 1e10
        err2 = 1e10
        isweep = 0
        while (err > thresh) or (err2 > thresh2):
            isweep += 1
            self.angles = np.zeros(num)
            if isweep == 1:
                self.n = 0
                old_E = self.objectiveFunction(0)
            for n in range(num):
                self.n = n
                bounds = [[0, 2 * np.pi]]
                result = optimize.minimize_scalar(self.objectiveFunction, bounds=bounds)
                if not result.success:
                    raise ConvergeError('Optimization to find optimal vibrational coordinates fails.')
                else:
                    self.angles[n] = result.x
            
            # Update vectors
            U = self.U(self.angles)
            self.mwv = self.mwv.dot(U)

            # Check convergence
            E = result.fun
            err = E - old_E
            err2 = abs(self.angles.sum())
            old_E = E
            
            if printing:
                logging.info('Normal mode localization: Cycle {:3d}    E: {:8.3f}   change: {:10.7f}  {:10.5f}'\
                             .format(isweep, E, err, err2))
    
    def Ui(self, angle, i, j):
        """
        Ui is a Jacobi rotation matrix of two by two rotation (among mode i and mode j).
        """
        Ui = np.identity(self.n_vib)
        if i > j:
            tmp = i
            i = j
            j = tmp
        c = np.cos(angle)
        s = np.sin(angle)       
        Ui[i][i] = c
        Ui[j][j] = c
        Ui[i][j] = -s
        Ui[j][i] = s
        return Ui

    def U(self, angles):
        """
        Matrix U is expressed as consecutive Jacobi rotations.
        """
        U = np.identity(self.n_vib)
        ind = 0
        for i in range(self.n_vib):
            for j in range(self.n_vib):
                if i < j:
                    angle = angles[ind]
                    Ui = self.Ui(angle, i, j)
                    U = np.matmul(U, Ui)
                    ind += 1
        return U