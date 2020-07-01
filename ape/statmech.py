# -*- coding: utf-8 -*-
import sys
import math
import numpy as np
from numpy import exp
from copy import deepcopy

import rmgpy.constants as constants

from ape.FitPES import from_sampling_result, cubic_spline_interpolations
from ape.schrodinger import SetAnharmonicH

class Statmech(object):
    """
    The class to solve shrodinger equation, evaluate partition function and related properties of 1-D PES by using statistical thermodynamics
    """
    def __init__(self, conformer, polynomial_dict, mode_dict, energy_dict,T,P=100000):
        self.conformer = conformer
        self.polynomial_dict = polynomial_dict
        self.mode_dict = mode_dict
        self.energy_dict = energy_dict
        self.T = T
        self.P = P
    
    def calcThermoOfEachMode(self,eig,N,mode):
        T = self.T
        beta = 1/(constants.kB*T) * constants.E_h
        Q = 0
        E = 0
        dQ = 0
        ddQ = 0
        for i in range(N):
            Ei = eig[i]
            Q += exp(-beta*Ei)
            dQ += Ei*exp(-beta*Ei)*beta/T
            ddQ += -2*Ei*exp(-beta*Ei)*beta/pow(T,2) + pow(Ei,2)*exp(-beta*Ei)*pow(beta,2)/pow(T,2)
            E += Ei*exp(-beta*Ei)
        E /= Q
        is_tors = True if self.mode_dict[mode]['mode'] == 'tors' else False
        if is_tors:
            omega = self.mode_dict[mode]['symmetry_number']
            Q /= omega
            dQ /= omega
            ddQ /= omega

        E0 = eig[0]
        v = (eig[1]-eig[0]) * constants.E_h / constants.h / (constants.c * 100) 
        #print(Q)

        F = -math.log(Q)/beta
        S = (E - F)/T
        Cv = (2/Q*dQ - T*pow(dQ/Q,2) + T/Q*ddQ)/beta

        return v, E0, E, S, F, Q, Cv

    def SolvEig(self,mode):
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
            v, E0, E, S, F, Q, Cv = self.calcThermoOfEachMode(eig,Nbasis,mode)

            if Qold == np.log(sys.float_info[0]):
                print("\n\t %d \t\t-\t\t-" % Nbasis) #first run
            else:
                print("\n\t %d \t\t %.10f \t\t %.10f" % (Nbasis,abs(Q-Qold),abs(v-vold)))
            
            if ((abs(Q-Qold)<1e-4) and (abs(v-vold)<1e-2)):
                print(" Convergence criterion met\n","------------------------------------")
                converge = True
                print("Frequency (cm-1): %.10f" % v)
                print("Zero point vibrational energy (hartree): %.10f" % E0)
                print("Energy (hartree): %.10f" % E )
                print("Entropy (hartree/K): %.10f" % S)
                print("Free energy (hartree): %.10f" % F)
                print("Partition function: %.10f" % Q)
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
        return v, E0, E, S, F, Q, Cv