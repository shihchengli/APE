# -*- coding: utf-8 -*-

"""
A module to find the optimizing vibrational coordinates to reduce intermode coupling
"""

import os
import time
import logging
import numpy as np
from copy import deepcopy
from scipy import optimize
from numba import jit

import rmgpy.constants as constants

from ape.job.job import Job
from ape.qchem import QChemLog
from ape.exceptions import InputError, ConvergeError
from ape.common import diagonalize_projected_hessian
from ape.intcoords.InternalCoordinates import getXYZ
from ape.intcoords.constants import BOHR2ANG

class OptVib(object):
    def __init__(self, symbols, nmode, coordinate_system, cart_coords, internal_object, conformer, hessian, linearity, n_vib, rotors, label, path, ncpus, 
                 charge=None, multiplicity=None, rem_variables_dict=None, gen_basis="", is_QM_MM_INTERFACE=None, QM_USER_CONNECT=None, QM_ATOMS=None,
                 ISOTOPE=None, force_field_params=None, fixed_molecule_string=None, opt=None):
        self.symbols = symbols
        self.nmode = nmode
        self.coordinate_system = coordinate_system
        self.cart_coords = cart_coords
        self.internal_object = internal_object
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
        self.rem_variables_dict = rem_variables_dict
        self.gen_basis = gen_basis
        self.n_rotors = len(self.rotors)

        # QMMM parameters
        self.is_QM_MM_INTERFACE = is_QM_MM_INTERFACE
        self.QM_USER_CONNECT = QM_USER_CONNECT
        self.QM_ATOMS = QM_ATOMS
        self.ISOTOPE = ISOTOPE
        self.force_field_params = force_field_params
        self.fixed_molecule_string = fixed_molecule_string
        self.opt = opt
        
    def get_optvib(self):
        """
        Algorithms for local and optimal vibrations.
        """
        logging.info('{0} modes of {1} finding...'.format(self.coordinate_system, self.label))
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
        logging.info('-------------------------------------------------------------------------------------------------------------------')
        logging.info('                                                  Jacobi sweeps                                                    ')
        logging.info('-------------------------------------------------------------------------------------------------------------------')
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
        logging.info('A grid of Hessians generating...\n')
        vib_freq, unweighted_v = diagonalize_projected_hessian(self.conformer, self.hessian, self.linearity, self.n_vib, self.rotors, label=self.label)
        grid_of_hessians = {}
        fm = diagonalize_projected_hessian(self.conformer, self.hessian, self.linearity, self.n_vib, self.rotors, get_mass_weighted_hessian=True, label=self.label)
        grid_of_hessians[0] = deepcopy(fm)
        for i in range(self.nmode):
            if i in range(self.n_rotors): continue
            logging.info('Sampling Mode {mode}'.format(mode=i+1))
            mode = i + 1
            vector = unweighted_v[i - self.n_rotors]
            freq = vib_freq[i - self.n_rotors]
            magnitude = np.linalg.norm(vector)
            reduced_mass = magnitude ** -2 / constants.amu # in amu
            step_size = np.sqrt(constants.hbar / (reduced_mass * constants.amu) / (freq * 2 * np.pi * constants.c * 100)) * 10 ** 10 # in angstrom
            normalizes_vector = vector / magnitude
            qj = np.matmul(self.internal_object.B, normalizes_vector/BOHR2ANG)
            qj = qj.reshape(-1,)
            
            initial_geometry = self.cart_coords.copy()
            cart_coords = initial_geometry.copy()
            internal = deepcopy(self.internal_object)

            cart_coords += internal.transform_int_step((qj * step_size).reshape(-1,)) * BOHR2ANG
            xyz = getXYZ(self.symbols, cart_coords)
            file_name = mode
            if self.is_QM_MM_INTERFACE:
                QMMM_xyz_string = ''
                for i, xyz in enumerate(xyz.split('\n')):
                    QMMM_xyz_string += " ".join([xyz, self.QM_USER_CONNECT[i]]) + '\n'
                    if i == len(self.QM_ATOMS)-1:
                        break
                QMMM_xyz_string += self.fixed_molecule_string
                job = Job(QMMM_xyz_string, self.path, file_name, jobtype='freq', ncpus=self.ncpus, charge=self.charge, multiplicity=self.multiplicity,
                          rem_variables_dict=self.rem_variables_dict, gen_basis=self.gen_basis, QM_atoms=self.QM_ATOMS, ISOTOPE=self.ISOTOPE,
                          force_field_params=self.force_field_params, opt=self.opt)
            else:
                job = Job(xyz, self.path, file_name, jobtype='freq', ncpus=self.ncpus, charge=self.charge, multiplicity=self.multiplicity,
                          rem_variables_dict=self.rem_variables_dict, gen_basis=self.gen_basis)

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
                E += E_Optimized_batch_run(hessian, V ,dim)
        elif self.coordinate_system == "E'-Optimized":
            H0 = V.T.dot(self.grid_of_hessians[0]).dot(V) / ((2 * np.pi * constants.c * 100) ** 2)
            for key in self.grid_of_hessians.keys():
                if key == 0: continue
                dim = self.n_vib
                hessian = self.grid_of_hessians[key]
                E += dE_Optimized_batch_run(hessian, V ,dim, H0)
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
        start = time.time()
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
            err = abs(E - old_E)
            err2 = abs(self.angles.sum())
            old_E = E
            
            if printing:
                logging.info('Normal mode localization: Cycle {:3d}    E: {:>25.7f}   change: {:>25.7f}  {:>10.5f}'\
                             .format(isweep, E, err, err2))
        end = time.time()
        logging.info('-------------------------------------------------------------------------------------------------------------------')
        logging.info('\nThe Jacobi sweeps have converged in {:.2f} s(wall)'.format(end - start))
    
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

#@jit(nopython=True)
def E_Optimized_batch_run(hessian, V, dim):
    E = 0
    H = V.T.dot(hessian).dot(V) / ((2 * np.pi * constants.c * 100) ** 2)
    for i in range(dim):
        for j in range(dim):
            if i < j:
                E += (H[i][j]) ** 2
    return E

#@jit(nopython=True)
def dE_Optimized_batch_run(hessian, V, dim, H0):
    E = 0
    H = V.T.dot(hessian).dot(V) / ((2 * np.pi * constants.c * 100) ** 2)
    for i in range(dim):
        for j in range(dim):
            if i < j:
                E += (H[i][j] - H0[i][j]) ** 2
    return E