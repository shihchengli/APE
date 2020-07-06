# -*- coding: utf-8 -*-
import math
import logging
import numpy as np
import pybel

import rmgpy.constants as constants

from ape.InternalCoordinates import get_RedundantCoords
from ape.exceptions import ConvergeError
from ape.common import SolvEig, mass_weighted_hessian

class HinderedRotor(object):

    def __init__(self, symbols, conformer, hessian, rotors_dict, linear, is_ts, n_vib):
        self.symbols = symbols
        self.conformer = conformer
        self.cart_coords = conformer.coordinates.value.reshape(-1,)
        self.mass = conformer.mass.value_si
        self.hessian = hessian
        self.rotors_dict = rotors_dict
        self.linear = linear
        self.is_ts = is_ts
        self.n_vib = n_vib

    def projectd_hessian(self):
        internal = get_RedundantCoords(self.symbols, self.cart_coords, self.rotors_dict)
        B = internal.B
        B_inv = internal.B_inv
        Bt_inv = internal.Bt_inv
        P = np.ones(B.shape[0],dtype=int)
        n_rotors = len(self.rotors_dict)
        if n_rotors != 0:
            P[-n_rotors:] = 0
        P = np.diag(P)
        projectd_hessian = B.T.dot(P).dot(Bt_inv).dot(self.hessian).dot(B_inv).dot(P).dot(B)
        return projectd_hessian

    def get_projected_out_freq(self, scan):
        """
        Calculate the vibrational frequency in the unit of cm^-1 of the internal rotation
        whose scan is provided by user.
        """
        internal = get_RedundantCoords(self.symbols, self.cart_coords)
        n_rotors = len(self.rotors_dict)
        B = internal.B[:-n_rotors]
        rotors_dict = self.rotors_dict
        # adding one dihedral internal related to torsion in B-matrix
        scan_indices = [scan[0]-1, scan[1]-1, scan[2]-1, scan[3]-1]
        val, grad = internal.calc_dihedral(internal.c3d, scan_indices, True)
        grad = grad.reshape(1,len(grad))
        B = np.concatenate((B, grad), axis=0)
        B_inv = np.linalg.pinv(B)
        Bt_inv = B_inv.T
        hessian = B.T.dot(Bt_inv).dot(self.hessian).dot(B_inv).dot(B)
        mwh = mass_weighted_hessian(self.conformer, hessian, self.linear, is_ts=self.is_ts)
        vib_freq, unweighted_v = SolvEig(mwh, self.mass, self.n_vib)

        ph = self.projectd_hessian()
        mwph = mass_weighted_hessian(self.conformer, ph, self.linear, is_ts=self.is_ts)
        projectd_vib_freq, projectd_unweighted_v = SolvEig(mwph, self.mass, self.n_vib)

        for i in vib_freq:
            match_freq = 0
            for j in projectd_vib_freq:
                if math.isclose(i,j,abs_tol=1) is True:
                    match_freq += 1
            if match_freq is 0:
                try:
                    index = vib_freq.tolist().index(i)
                    vector = unweighted_v[index]
                    magnitude = np.linalg.norm(vector)
                    reduced_mass = magnitude ** -2 / 1.660538921e-27 # in amu
                    projected_out_freq = i
                    break
                except ValueError:
                    pass
            if i == vib_freq[-1] and j == projectd_vib_freq[-1]:
                raise ConvergeError('Can\'t find the frequency of the hindered rotor whose scan is {scan}'.format(scan=scan))
        
        logging.info('The frequency of the hindered rotor whose scan is {scan} is {freq} cm-1'\
            .format(scan=scan, freq=projected_out_freq))
    
        return projected_out_freq, reduced_mass