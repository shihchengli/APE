# -*- coding: utf-8 -*-
import sys
import math
import numpy as np
from numpy import exp
from copy import deepcopy

import rmgpy.constants as constants

from ape.FitPES import from_sampling_result, cubic_spline_interpolations
from ape.schrodinger import SetAnharmonicH

class ThermoJob(object):
    """
    The class to calculate thermodynamic properties, including E S G Cp
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

    def calcThermo(self, print_HOhf_result=False, zpe_of_Hohf=None):
        T = self.T
        P = self.P
        conformer = self.conformer
        print("Calculate global translation and rotation E, S")
        E_trans = 1.5 * constants.R * T / 4184
        S_trans = conformer.modes[0].get_entropy(T) / 4.184 - constants.R * math.log(P/101325) / 4.184
        Cv_trans = 1.5 * constants.R / 4184 * 1000  
        #F = conformer.modes[0].get_free_energy(T) / 4.184
        #Q = conformer.modes[0].get_partition_function

        E_rot = conformer.modes[1].get_enthalpy(T) / 4184
        S_rot = conformer.modes[1].get_entropy(T) / 4.184
        Cv_rot = conformer.modes[1].get_heat_capacity(T) / 4.184
        #F = conformer.modes[1].get_free_energy(T) / 4.184
        #Q = conformer.modes[1].get_partition_function

        print("Calculate internal E, S")
        E0 = 0
        E_int = 0
        S_int = 0
        #F_int = 0
        #Q_int = 1
        Cv_int = 0

        for mode in sorted(self.mode_dict.keys()):
            print("\n\t********** Mode ",mode," **********\n\n")
            v, e0, E, S, F, Q, Cv = self.SolvEig(mode)
            E0 += e0
            E_int += E
            S_int += S
            #F_int += F
            #Q_int *= Q
            Cv_int += Cv

        print("\n\t********** Final results **********\n\n")
        print("Temperature (K): ",T)
        print("Pressure (Pa): ",P)
        print("Zero point vibrational energy (kcal/mol): %.10f" % (E0))
        print("Translational energy (kcal/mol): %.10f" % (E_trans))
        print("Translational entropy (cal/mol/K): %.10f" % (S_trans))
        print("Translational Cv (cal/mol/K): %.10f" % (Cv_trans))
        print("Rotational energy (kcal/mol): %.10f" % (E_rot))
        print("Rotational entropy (cal/mol/K): %.10f" % (S_rot))
        print("Rotational Cv (cal/mol/K): %.10f" % (Cv_rot))
        print("Internal (rot+vib) energy (kcal/mol): %.10f" % (E_int))
        print("Internal (tor+vib) entropy (cal/mol/K): %.10f" % (S_int))
        print("Internal (tor+vib) Cv (cal/mol/K): %.10f" % (Cv_int))
        print("\n")
        print("Total energy (kcal/mol): %.10f" % (E_trans+E_rot+E_int))
        print("Total enthalpy (kcal/mol): %.10f" % (E_trans+E_rot+E_int+constants.kB*T*constants.Na/4184))
        print("Enthalpy H(%f K)-H(0 K) (kcal/mol):  %.10f" % (T, E_trans+E_rot+E_int+constants.kB*T*constants.Na/4184-E0))
        print("Total entropy (cal/mol/K): %.10f" % (S_trans+S_rot+S_int))
        print("Total Cv (cal/mol/K): %.10f" % (Cv_trans+Cv_rot+Cv_int))

        if print_HOhf_result:
            # compare to HOhf model
            E_vib = (conformer.modes[2].get_enthalpy(T) + zpe_of_Hohf) / 4184
            # E_vib should be calculated by freq...
            S_vib = conformer.modes[2].get_entropy(T) / 4.184
            print("\n")
            print("\n\t********** HOhf results **********\n\n")
            print("Translational energy (kcal/mol): %.10f" % (E_trans))
            print("Rotational energy (kcal/mol): %.10f" % (E_rot))
            print("Vibrational energy (kcal/mol): %.10f" % (E_vib))
            print("gas constant (RT): %.10f" % (constants.R * T / 4184))
            print("Translational entropy (cal/mol/K): %.10f" % (S_trans))
            print("Rotational entropy (cal/mol/K): %.10f" % (S_rot))
            print("Vibrational entropy (cal/mol/K): %.10f" % (S_vib))
            print("\n")
            print("Total energy (kcal/mol): %.10f" % (E_trans+E_rot+E_vib))
            print("Total enthalpy (kcal/mol): %.10f" % (E_trans+E_rot+E_vib+constants.R*T/4184))
            print("Enthalpy H(%f K)-H(0 K) (kcal/mol): %.10f" % (T, conformer.get_enthalpy(T)/4184))
            print("Total entropy (cal/mol/K): %.10f" % (S_trans+S_rot+S_vib))
            print("Total Cv (cal/mol/K): %.10f" % (conformer.get_heat_capacity(T) / 4.184))
    
    def calcQMMMThermo(self):
        T = self.T
        P = self.P

        print("Calculate internal E, S")
        E0 = 0
        E_int = 0
        S_int = 0
        #F_int = 0
        #Q_int = 1
        Cv_int = 0

        for mode in sorted(self.mode_dict.keys()):
            print("\n\t********** Mode ",mode," **********\n\n")
            v, e0, E, S, F, Q, Cv = self.SolvEig(mode)
            E0 += e0
            E_int += E
            S_int += S
            #F_int += F
            #Q_int *= Q
            Cv_int += Cv

        print("\n\t********** Final results **********\n\n")
        print("Temperature (K): ",T)
        print("Pressure (Pa): ",P)
        print("Zero point vibrational energy (kcal/mol): %.10f" % (E0))
        print("Internal (rot+vib) energy (kcal/mol): %.10f" % (E_int))
        print("Internal (tor+vib) entropy (cal/mol/K): %.10f" % (S_int))
        print("Internal (tor+vib) Cv (cal/mol/K): %.10f" % (Cv_int))

################################################################################

if __name__ == '__main__':
    csv_path = '../examples/propane_UMVT/samping_result.csv'
    freq_file = '../examples/propane_UMVT/propane.q.out'
    from main import APE
    ape = APE(freq_file)
    ape.parse()
    conformer = ape.conformer
    mode_dict, energy_dict = from_sampling_result(csv_path)
    polynomial_dict = cubic_spline_interpolations(energy_dict,mode_dict)
    thermo = ThermoJob(conformer, polynomial_dict, mode_dict, energy_dict,T=298.15,P=100000)
    thermo.calcQMMMThermo()
    #thermo.calcThermo(print_HOhf_result=True, zpe_of_Hohf=ape.zpe)