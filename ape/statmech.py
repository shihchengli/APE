# -*- coding: utf-8 -*-
import os
import sys
import math
import logging
import numpy as np
from numpy import exp
from copy import deepcopy

import rmgpy.constants as constants
from rmgpy.statmech.vibration import HarmonicOscillator

from ape.sampling import SamplingJob
from ape.FitPES import from_sampling_result, cubic_spline_interpolations
from ape.schrodinger import SetAnharmonicH

class Statmech(object):
    """
    A class to solve shrodinger equation, evaluate partition function and related properties of 1-D PES by using statistical thermodynamics
    """
    def __init__(self, label, input_file, output_directory, Tlist=[298.15], P=100000, frequency_scale_factor=1, zpe_of_Hohf=None):
        self.label = label
        self.input_file = input_file
        self.output_directory = output_directory
        self.Tlist = Tlist
        self.P = P
        self.frequency_scale_factor = frequency_scale_factor
        self.result_info = list()
    
    def load_save(self):
        self.sampling = SamplingJob(self.label, self.input_file)
        self.sampling.parse()
        self.conformer = self.sampling.conformer
        self.csv_path = os.path.join(self.output_directory, '{}_samping_result.csv'.format(self.label))
        self.mode_dict, self.energy_dict, self.min_elect = from_sampling_result(self.csv_path)
        self.zpe_of_Hohf = self.sampling.zpe
        e0 = self.min_elect * constants.E_h * constants.Na + self.sampling.zpe
        self.conformer.E0 = (e0, "J/mol")
        for mode in self.conformer.modes:
            if isinstance(mode, HarmonicOscillator):
                frequencies = mode.frequencies.value_si
                mode.frequencies = (frequencies * self.frequency_scale_factor, "cm^-1")
        self.spin_multiplicity = self.conformer.spin_multiplicity
        self.optical_isomers =self.conformer.optical_isomers
        self.symbols = self.sampling.symbols

        # Solve SE of 1-D PES and calculate E S G Cp
        self.polynomial_dict = cubic_spline_interpolations(self.energy_dict, self.mode_dict)

        # Extract whether this system is QM/MM system or not
        self.is_QM_MM_INTERFACE = self.sampling.is_QM_MM_INTERFACE

    def calcThermoOfEachMode(self, eig, N, mode, T):
        beta = 1/(constants.kB*T) * constants.E_h
        Q = 0
        Q_vib = 0
        E = 0
        dQ = 0
        ddQ = 0
        for i in range(N):
            Ei = eig[i]
            Q += exp(-beta*Ei)
            dQ += Ei*exp(-beta*Ei)*beta/T
            ddQ += -2*Ei*exp(-beta*Ei)*beta/pow(T,2) + pow(Ei,2)*exp(-beta*Ei)*pow(beta,2)/pow(T,2)
            E += Ei*exp(-beta*Ei)
            if i == 0:
                zpve = Ei
            # Measuring energy relative to the zero point vibration frequency
            dE = Ei - zpve
            Q_vib += exp(-beta*dE)
        E /= Q
        is_tors = True if self.mode_dict[mode]['mode'] == 'tors' else False
        if is_tors:
            omega = self.mode_dict[mode]['symmetry_number']
            Q /= omega
            Q_vib /= omega
            dQ /= omega
            ddQ /= omega

        E0 = eig[0]
        v = (eig[1]-eig[0]) * constants.E_h / constants.h / (constants.c * 100) 
        #print(Q)

        F = -math.log(Q)/beta
        S = (E - F)/T
        Cv = (2/Q*dQ - T*pow(dQ/Q,2) + T/Q*ddQ)/beta

        return v, E0, E, S, F, Q, Q_vib, Cv

    def SolvEig(self, mode, T):
        Nbasis = 50
        Nbasis_prev = 0
        H_prev = None
        Qold =  np.log(sys.float_info[0])
        vold =  np.log(sys.float_info[0])
        converge = False
        while not converge:
            Nbasis += 1
            H = SetAnharmonicH(self.polynomial_dict, self.mode_dict, self.energy_dict, mode, Nbasis, N_prev=Nbasis_prev, H_prev=H_prev)
            Nbasis_prev = Nbasis
            H_prev = deepcopy(H)
            eig, v =np.linalg.eigh(H)
            v, E0, E, S, F, Q, Q_vib, Cv = self.calcThermoOfEachMode(eig, Nbasis, mode, T)

            if Qold == np.log(sys.float_info[0]):
                self.result_info.append("# \n# \t %d \t\t-\t\t-" % Nbasis) #first run
                logging.debug("# \t {} \t\t-\t\t-".format(Nbasis))
            else:
                self.result_info.append("# \n# \t %d \t\t %.10f \t\t %.10f" % (Nbasis, abs(Q-Qold), abs(v-vold)))
                logging.debug("# \t {:d} \t\t {:.10f} \t\t {:.10f}".format(Nbasis, abs(Q-Qold), abs(v-vold)))
            
            if ((abs(Q-Qold)<1e-4) and (abs(v-vold)<1e-2)):
                self.result_info.append("# Convergence criterion met")
                self.result_info.append("# ------------------------------------")
                converge = True
                self.result_info.append("# Frequency (cm-1): %.10f" % v)
                self.result_info.append("# Zero point vibrational energy (hartree): %.10f" % E0)
                self.result_info.append("# Energy (hartree): %.10f" % E )
                self.result_info.append("# Entropy (hartree/K): %.10f" % S)
                self.result_info.append("# Free energy (hartree): %.10f" % F)
                self.result_info.append("# Partition function: %.10f" % Q)
                hartree2kcalmol = constants.E_h * constants.Na / 4184
                E0 *= hartree2kcalmol
                E *= hartree2kcalmol
                S *= hartree2kcalmol * 1000
                F *= hartree2kcalmol
                Cv *= hartree2kcalmol * 1000
                '''
                print("Frequency (cm-1): ",v)
                print("Zero point vibrational energy (kcal/mol): ",E0)
                print("Energy (kcal/mol): ",E )
                print("Entropy (cal/mol/K): ",S)
                print("Free energy (kcal/mol): ",F)
                print("Partition function: ",Q)
                '''
            
            Qold = Q
            vold = v
        return v, E0, E, S, F, Q, Q_vib, Cv