# -*- coding: utf-8 -*-
import os
import math
import logging

import rmgpy.constants as constants

from arkane.output import prettify

from ape.FitPES import from_sampling_result, cubic_spline_interpolations
from ape.schrodinger import SetAnharmonicH
from ape.statmech import Statmech

class ThermoJob(Statmech):
    """
    The class to calculate thermodynamic properties, including E S G Cp
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def calcThermo(self, T, print_HOhf_result=True):
        P = self.P
        conformer = self.conformer
        # Calculate global translation and rotation E, S
        E_trans = 1.5 * constants.R * T / 4184
        S_trans = conformer.modes[0].get_entropy(T) / 4.184 - constants.R * math.log(P / 101325) / 4.184
        Cv_trans = 1.5 * constants.R / 4184 * 1000
        Q_trans = conformer.modes[0].get_partition_function(T)

        E_rot = conformer.modes[1].get_enthalpy(T) / 4184
        S_rot = conformer.modes[1].get_entropy(T) / 4.184
        Cv_rot = conformer.modes[1].get_heat_capacity(T) / 4.184
        Q_rot = conformer.modes[1].get_partition_function(T)

        # logging.info("Calculate internal E, S")
        ZPE = 0
        E_int = 0
        S_int = 0
        # F_int = 0
        Q_int = 1
        Cv_int = 0

        for mode in sorted(self.mode_dict.keys()):
            self.result_info.append("\n# \t********** Mode {} **********".format(mode))
            v, e0, E, S, F, Q, Cv = self.SolvEig(mode, T)
            ZPE += e0
            E_int += E
            S_int += S
            # F_int += F
            Q_int *= Q
            Cv_int += Cv

        self.result_info.append("\n# \t********** Final results **********\n\n")
        self.result_info.append("# Temperature (K): %.2f" % (T))
        self.result_info.append("# Pressure (Pa): %.0f" % (P))
        self.result_info.append("# Zero point vibrational energy (kcal/mol): %.10f" % (ZPE))
        self.result_info.append("# Translational energy (kcal/mol): %.10f" % (E_trans))
        self.result_info.append("# Translational entropy (cal/mol/K): %.10f" % (S_trans))
        self.result_info.append("# Translational Cv (cal/mol/K): %.10f" % (Cv_trans))
        self.result_info.append("# Rotational energy (kcal/mol): %.10f" % (E_rot))
        self.result_info.append("# Rotational entropy (cal/mol/K): %.10f" % (S_rot))
        self.result_info.append("# Rotational Cv (cal/mol/K): %.10f" % (Cv_rot))
        self.result_info.append("# Internal (rot+vib) energy (kcal/mol): %.10f" % (E_int))
        self.result_info.append("# Internal (tor+vib) entropy (cal/mol/K): %.10f" % (S_int))
        self.result_info.append("# Internal (tor+vib) Cv (cal/mol/K): %.10f" % (Cv_int))
        self.result_info.append("\n")
        self.result_info.append("# Total energy (kcal/mol): %.10f" % (E_trans + E_rot + E_int))
        self.result_info.append("# Total enthalpy (kcal/mol): %.10f" % (E_trans + E_rot + E_int + constants.kB * T * constants.Na / 4184))
        self.result_info.append("# Enthalpy H(%f K)-H(0 K) (kcal/mol):  %.10f" % (T, E_trans + E_rot + E_int + constants.kB * T * constants.Na / 4184 - ZPE))
        self.result_info.append("# Total entropy (cal/mol/K): %.10f" % (S_trans + S_rot + S_int))
        self.result_info.append("# Total Cv (cal/mol/K): %.10f" % (Cv_trans + Cv_rot + Cv_int))

        if print_HOhf_result:
            # compare to HOhf model
            E_vib = (conformer.modes[2].get_enthalpy(T) + self.zpe_of_Hohf) / 4184
            # E_vib should be calculated by freq...
            S_vib = conformer.modes[2].get_entropy(T) / 4.184
            self.result_info.append("\n")
            self.result_info.append("\n# \t********** HOhf results **********\n\n")
            self.result_info.append("# Translational energy (kcal/mol): %.10f" % (E_trans))
            self.result_info.append("# Rotational energy (kcal/mol): %.10f" % (E_rot))
            self.result_info.append("# Vibrational energy (kcal/mol): %.10f" % (E_vib))
            self.result_info.append("# gas constant (RT): %.10f" % (constants.R * T / 4184))
            self.result_info.append("# Translational entropy (cal/mol/K): %.10f" % (S_trans))
            self.result_info.append("# Rotational entropy (cal/mol/K): %.10f" % (S_rot))
            self.result_info.append("# Vibrational entropy (cal/mol/K): %.10f" % (S_vib))
            self.result_info.append("\n")
            self.result_info.append("# Total energy (kcal/mol): %.10f" % (E_trans + E_rot + E_vib))
            self.result_info.append("# Total enthalpy (kcal/mol): %.10f" % (E_trans + E_rot + E_vib + constants.R * T / 4184))
            self.result_info.append("# Enthalpy H(%f K)-H(0 K) (kcal/mol): %.10f" % (T, conformer.get_enthalpy(T) / 4184))
            self.result_info.append("# Total entropy (cal/mol/K): %.10f" % (S_trans + S_rot + S_vib))
            self.result_info.append("# Total Cv (cal/mol/K): %.10f" % (conformer.get_heat_capacity(T) / 4.184))
        
        
        E0 = (self.conformer.E0.value_si - self.zpe_of_Hohf) * 0.001 / 4.184  + ZPE # in kcal/mol
        E = E_trans + E_rot + E_int # in kcal/mol
        S = S_trans + S_rot + S_int # in cal/mol/K
        F = (E + constants.R * T / 4184 - ZPE) - T * S * 0.001 + E0 # in kcal/mol
        Q = Q_trans * Q_rot * Q_int
        Cv = Cv_trans + Cv_rot + Cv_int # in cal/mol/K
        return  E0, E, S, F, Q, Cv
    
    def calcQMMMThermo(self, T, print_HOhf_result=True):
        P = self.P
        conformer = self.conformer
        logging.info("Calculate internal E, S")
        ZPE = 0
        E_int = 0
        S_int = 0
        # F_int = 0
        # Q_int = 1
        Cv_int = 0

        for mode in sorted(self.mode_dict.keys()):
            self.result_info.append("\n# \t********** Mode ",mode," **********\n\n")
            v, e0, E, S, F, Q, Cv = self.SolvEig(mode, T)
            ZPE += e0
            E_int += E
            S_int += S
            # F_int += F
            # Q_int *= Q
            Cv_int += Cv

        self.result_info.append("\n# \t********** Final results **********\n#\n")
        self.result_info.append("# Temperature (K): %.2f" % (T))
        self.result_info.append("# Pressure (Pa): %.0f" % (P))
        self.result_info.append("# Zero point vibrational energy (kcal/mol): %.10f" % (ZPE))
        self.result_info.append("# Internal (rot+vib) energy (kcal/mol): %.10f" % (E_int))
        self.result_info.append("# Internal (tor+vib) entropy (cal/mol/K): %.10f" % (S_int))
        self.result_info.append("# Internal (tor+vib) Cv (cal/mol/K): %.10f" % (Cv_int))

        if print_HOhf_result:
            # compare to HOhf model
            E_vib = (conformer.modes[2].get_enthalpy(T) + self.zpe_of_Hohf) / 4184
            # E_vib should be calculated by freq...
            S_vib = conformer.modes[2].get_entropy(T) / 4.184
            self.result_info.append("\n")
            self.result_info.append("\n# \t********** HOhf results **********\n\n")
            self.result_info.append("# Vibrational energy (kcal/mol): %.10f" % (E_vib))
            self.result_info.append("# Vibrational entropy (cal/mol/K): %.10f" % (S_vib))
        
        self.result_info.append('\n\n\n')

    
    def write_output(self):
        """
        Save the results of the ThermoJob to the `output.py` file located
        in `output_directory`.
        """
        output_file = os.path.join(self.output_directory, 'output.py')
        logging.info('Saving statistical mechanics parameters for {0}...'.format(self.label))
        f = open(output_file, 'a')

        conformer = self.conformer
        coordinates = conformer.coordinates.value_si * 1e10
        number = conformer.number.value_si

        f.write('# Coordinates for {0} in Input Orientation (angstroms):\n'.format(self.label))
        for i in range(coordinates.shape[0]):
            x = coordinates[i, 0]
            y = coordinates[i, 1]
            z = coordinates[i, 2]
            f.write('#   {0} {1:9.4f} {2:9.4f} {3:9.4f}\n'.format(self.symbols[i], x, y, z))
        f.write('\n')

        result = 'conformer(label={0!r}, E0={1!r}, modes={2!r}, spin_multiplicity={3:d}, optical_isomers={4:d}'.format(
            self.label,
            conformer.E0,
            conformer.modes,
            conformer.spin_multiplicity,
            conformer.optical_isomers,
        )
        try:
            result += ', frequency={0!r}'.format(self.sampling.imaginary_frequency)
        except AttributeError:
            pass
        result += ')'
        f.write('{0}\n\n'.format(prettify(result)))

        for line in self.result_info:
            line = line + '\n'
            f.write(line)
        f.write('\n')
        f.close()

    def execute(self):
        logging.info('Calculate thermodynamics for {0}'.format(self.label))
        for T in self.Tlist:
            self.result_info.append('\n\n# Thermodynamics for {0} at {1} K:\n'.format(self.label, T))
            if self.is_QM_MM_INTERFACE:
                self.calcQMMMThermo(T=T, print_HOhf_result=True)
            else:
                self.calcThermo(T=T, print_HOhf_result=True)
        self.write_output()