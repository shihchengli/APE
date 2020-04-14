# -*- coding: utf-8 -*-
import math
import numpy as np
import pybel

from ape.InternalCoordinates import get_RedundantCoords
from ape.exceptions import ConvergeError

class HinderedRotor(object):

    def __init__(self, symbols, cart_coords, hessian, rotors_dict, mass=None, n_vib=None):
        self.symbols = symbols
        self.cart_coords = cart_coords
        self.hessian = hessian
        self.rotors_dict = rotors_dict
        self.mass = mass
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
        from ape import main
        internal = get_RedundantCoords(self.symbols, self.cart_coords)
        n_rotors = len(self.rotors_dict)
        B = internal.B[:-n_rotors]
        rotors_dict = self.rotors_dict
        #adding one dihedral internal related to torsion in B-matrix
        scan_indices = [scan[0]-1, scan[1]-1, scan[2]-1, scan[3]-1]
        val, grad = internal.calc_dihedral(internal.c3d, scan_indices, True)
        grad = grad.reshape(1,len(grad))
        B = np.concatenate((B, grad), axis=0)
        B_inv = np.linalg.pinv(B)
        Bt_inv = B_inv.T
        hessian = B.T.dot(Bt_inv).dot(self.hessian).dot(B_inv).dot(B)
        vib_freq, unweighted_v = main.SolvEig(hessian, self.mass, self.n_vib+1)

        projectd_hessian = self.projectd_hessian()
        projectd_vib_freq, projectd_unweighted_v = main.SolvEig(projectd_hessian, self.mass, self.n_vib+1)

        for i in vib_freq:
            match_freq = 0
            for j in projectd_vib_freq:
                if math.isclose(i,j,abs_tol=1) is True:
                    match_freq += 1
            if match_freq is 0:
                projected_out_freq = i
                break
            if i == vib_freq[-1] and j == projectd_vib_freq[-1]:
                raise ConvergeError('Can\'t find the frequency of the hindered rotor whose scan is {scan}'.format(scan=scan))
        
        print('The frequency of the hindered rotor whose scan is {scan} is {freq} cm-1'\
            .format(scan=scan, freq=projected_out_freq))
    
        return projected_out_freq
