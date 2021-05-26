#!/usr/bin/env python3

"""
APE QChem module
Used to parse QChem output files
"""

import math
import logging
import os.path

import numpy as np

import rmgpy.constants as constants
from rmgpy.statmech import IdealGasTranslation, NonlinearRotor, LinearRotor, HarmonicOscillator, Conformer

from arkane.common import check_conformer_energy, get_element_mass
from arkane.exceptions import LogError
from arkane.ess.adapter import ESSAdapter
from arkane.ess.factory import register_ess_adapter

################################################################################


class QChemLog(ESSAdapter):
    """
    Represent an output file from QChem. The attribute `path` refers to the
    location on disk of the QChem output file of interest. Methods are provided
    to extract a variety of information into Arkane classes and/or NumPy
    arrays.
    """

    def __init__(self, path):
        self.path = path

    def job_is_finished(self):
        with open(self.path, 'r') as f:
            line = f.readline()
            while line != '':
                # This marks the end of the thermochemistry section
                if 'Thank you very much for using Q-Chem.' in line:
                    return True
                line = f.readline()
        return False

    def get_basis(self):
        """
        Return the basis of this job.
        """
        gen_basis_params = ''
        with open(self.path, 'r') as f:
            line = f.readline()
            while line != '':
                if line.startswith('$rem'):
                    while '$end' not in line:
                        if 'BASIS' in line.upper():
                            basis = line.split()[-1].upper()
                        line = f.readline()
                if line.startswith('$basis'): # basis defined by user
                    line = f.readline()
                    gen_basis_params = '\n$basis\n'
                    while '$end' not in line:
                        gen_basis_params += line
                        line = f.readline()
                    gen_basis_params += '$end\n'
                line = f.readline()
        return basis, gen_basis_params

    def is_unrestricted(self):
        is_unrestricted = False
        with open(self.path, 'r') as f:
            line = f.readline()
            while line != '':
                if line.startswith('$rem'):
                    while '$end' not in line:
                        if 'UNRESTRICTED' in line.upper():
                            OPTION = line.split()[1].upper()
                            if OPTION == 'TRUE':
                                is_unrestricted = True
                            else:
                                ISOTOPES_type = int(OPTION)
                        line = f.readline()
                    break
                line = f.readline()

        return is_unrestricted

    def is_QM_MM_INTERFACE(self):
        """
        Return the bool value.
        """
        is_QM_MM_INTERFACE = False

        with open(self.path, 'r') as f:
            line = f.readline()
            while line != '':
                if 'QM_MM_INTERFACE' in line.upper():
                    is_QM_MM_INTERFACE = True
                line = f.readline()

        return is_QM_MM_INTERFACE

    def get_QM_ATOMS(self):
        """
        Return the index of QM_atoms.
        """
        QM_atoms = []

        with open(self.path, 'r') as f:
            line = f.readline()
            while line != '':
                if '$QM_ATOMS' in line.upper():
                    line = f.readline()
                    while '$end' not in line and line != '\n':
                        QM_atoms.append(line.strip())
                        line = f.readline()
                    break
                line = f.readline()
            
        return QM_atoms

    def get_ISOTOPES(self):
        """
        Return the $ISOTOPES used in the frequency calculation.
        """
        ISOTOPES = {}
        ISOTOPES_type = False

        with open(self.path, 'r') as f:
            line = f.readline()
            while line != '':
                if line.startswith('$rem'):
                    while '$end' not in line:
                        if 'ISOTOPES' in line.upper():
                            OPTION = line.split()[1].upper()
                            if OPTION == 'TRUE':
                                ISOTOPES_type = True
                            else:
                                ISOTOPES_type = int(OPTION)
                        line = f.readline()
                if '$ISOTOPES' in line:
                    if ISOTOPES_type is True:
                        for i in range(3):
                            line = f.readline()
                        while '$end' not in line:
                            data = line.split()
                            atom_number = int(data[0])
                            mass = int(data[1])
                            ISOTOPES[atom_number] = mass
                            line = f.readline()
                        break
                line = f.readline()

        return ISOTOPES

    def get_force_field_params(self):
        """
        Return the sting of $force_field_params in the quantum calculation.
        """
        force_field_params = ''

        with open(self.path, 'r') as f:
            line = f.readline()
            while line != '':
                if '$force_field_params' in line:
                    line = f.readline()
                    while '$end' not in line:
                        force_field_params += line
                        line = f.readline()
                line = f.readline()
            
        return force_field_params
    
    def get_opt(self):
        """
        Return the sting of $opt in the quantum calculation.
        """
        opt = ''

        with open(self.path, 'r') as f:
            line = f.readline()
            while line != '':
                if '$opt' in line:
                    line = f.readline()
                    while '$end' not in line:
                        opt += line
                        line = f.readline()
                    break
                line = f.readline()
            
        return opt
    
    def get_fixed_molecule(self):
        """
        Return the string of fixed part in $moledule
        """
        # <Atom> <X> <Y> <Z> <MM atom type> <Bond 1> <Bond 2> <Bond 3> <Bond 4>
        # O 7.256000 1.298000 9.826000 -1  185  186  0 0
        # O 6.404000 1.114000 12.310000 -1  186  713  0 0
        # O 4.077000 1.069000 0.082000 -1  188  187  0 0
        # H 1.825000 1.405000 12.197000 -3  714  0  0 0
        # H 2.151000 1.129000 9.563000 -3  189  0  0 0
        # -----------------------------------

        fixed_molecule_string = ''

        n_atoms = len(self.get_QM_ATOMS())
        with open(self.path, 'r') as f:
            line = f.readline()
            while line != '':
                if '$MOLECULE' in line.upper():
                    for i in range(2):
                        line = f.readline()
                    for i in range(n_atoms):
                        line = f.readline()
                    while '$end' not in line and line != '\n':
                        fixed_molecule_string += line
                        line = f.readline()
                    break
                line = f.readline()

        return fixed_molecule_string
    
    def get_QM_USER_CONNECT(self):
        """
        Return a list of string with "<MM atom type> <Bond 1> <Bond 2> <Bond 3> <Bond 4>" in QM region.
        """
        QM_USER_CONNECT = []

        n_atoms = len(self.get_QM_ATOMS())
        with open(self.path, 'r') as f:
            line = f.readline()
            while line != '':
                if '$MOLECULE' in line.upper():

                    for i in range(2):
                        line = f.readline()
                    for i in range(n_atoms):
                        _, _, _, _, MM_atom_type, Bond_1, Bond_2, Bond_3, Bond_4 = line.split()
                        USER_CONNECT = '{}  {}  {}  {}  {}'.format(MM_atom_type, Bond_1, Bond_2, Bond_3, Bond_4)
                        QM_USER_CONNECT.append(USER_CONNECT)
                        line = f.readline()
                    break
                line = f.readline()
        
        return QM_USER_CONNECT

    def get_number_of_atoms(self):
        """
        Return the number of atoms in the molecular configuration used in
        the QChem output file.
        """
        n_atoms = 0

        with open(self.path, 'r') as f:
            line = f.readline()
            while line != '' and n_atoms == 0:
                # Automatically determine the number of atoms
                if 'Standard Nuclear Orientation' in line and n_atoms == 0:
                    for i in range(3):
                        line = f.readline()
                    while '----------------------------------------------------' not in line:
                        n_atoms += 1
                        line = f.readline()

                line = f.readline()

        return n_atoms

    def load_force_constant_matrix(self):
        """
        Return the force constant matrix (in Cartesian coordinates) from the
        QChem log file. If multiple such matrices are identified,
        only the last is returned. The units of the returned force constants
        are J/m^2. If no force constant matrix can be found in the log file,
        ``None`` is returned.
        """
        force = None

        if self.is_QM_MM_INTERFACE():
            n_atoms = len(self.get_QM_ATOMS()) + len(self.get_ISOTOPES())
        else:
            n_atoms = self.get_number_of_atoms()
        n_rows = n_atoms * 3
        with open(self.path, 'r') as f:
            line = f.readline()
            while line != '':
                # Read force constant matrix
                if 'Final Hessian.' in line or 'Hessian of the SCF Energy' in line:
                    force = np.zeros((n_rows, n_rows), np.float64)
                    for i in range(int(math.ceil(n_rows / 6.0))):
                        # Header row
                        line = f.readline()
                        # Matrix element rows
                        for j in range(n_rows):  # for j in range(i*6, Nrows):
                            data = f.readline().split()
                            for k in range(len(data) - 1):
                                force[j, i * 6 + k] = float(data[k + 1])
                                # F[i*5+k,j] = F[j,i*5+k]
                    # Convert from atomic units (Hartree/Bohr_radius^2) to J/m^2
                    force *= 4.35974417e-18 / 5.291772108e-11 ** 2
                if 'FINAL TENSOR RESULT' in line:
                    force = np.zeros((n_rows, n_rows), np.float64)
                    try:
                        for i in range(int(n_atoms)):
                            # Header row
                            line = f.readline()
                            # Matrix element rows
                            for j in range(3):
                                line = f.readline()
                                line = f.readline()
                                line = f.readline()
                                for k in range(int(n_atoms)):
                                    data = f.readline().split()
                                    for l in range(3):
                                        force[i * 3 + j, k * 3 + l] = float(data[l + 1])
                    except (ValueError, IndexError):
                        continue
                    finally:
                        # Convert from atomic units (Hartree/Bohr_radius^2) to J/m^2
                        force *= 4.35974417e-18 / 5.291772108e-11 ** 2

                line = f.readline()

        return force

    def load_geometry(self):

        """
        Return the optimum geometry of the molecular configuration from the
        QChem log file. If multiple such geometries are identified, only the
        last is returned.
        """
        atom, coord, number, mass = [], [], [], []

        with open(self.path) as f:
            log = f.readlines()

        # First check that the QChem job file (not necessarily a geometry optimization)
        # has successfully completed, if not an error is thrown
        completed_job = False
        for line in reversed(log):
            if 'Total job time:' in line:
                logging.debug('Found a successfully completed QChem Job')
                completed_job = True
                break

        if not completed_job:
            raise LogError('Could not find a successfully completed QChem job '
                           'in QChem output file {0}'.format(self.path))

        # Now look for the geometry.
        # Will return the final geometry in the file under Standard Nuclear Orientation.
        geometry_flag = False
        for i in reversed(range(len(log))):
            line = log[i]
            if 'Standard Nuclear Orientation' in line:
                for line in log[(i + 3):]:
                    if '------------' not in line:
                        data = line.split()
                        atom.append(data[1])
                        coord.append([float(c) for c in data[2:]])
                        geometry_flag = True
                    else:
                        break
                if geometry_flag:
                    break

        # Assign appropriate mass to each atom in the molecule
        for atom1 in atom:
            mass1, num1 = get_element_mass(atom1)
            mass.append(mass1)
            number.append(num1)
        coord = np.array(coord, np.float64)
        number = np.array(number, np.int)
        mass = np.array(mass, np.float64)
        if len(number) == 0 or len(coord) == 0 or len(mass) == 0:
            raise LogError('Unable to read atoms from QChem geometry output file {0}.'.format(self.path))
        
        if self.is_QM_MM_INTERFACE():
            QM_mass = []
            for i in self.get_QM_ATOMS():
                QM_mass.append(mass[int(i)-1])
            ISOTOPES = self.get_ISOTOPES()
            for i in sorted(ISOTOPES.keys()):
                QM_mass.append(ISOTOPES[i])
                #QM_mass.append(get_element_mass('H')[0])
            self.QM_mass = np.array(QM_mass, np.float64)

            QM_atom = []
            QM_coord = []
            for i in reversed(range(len(log))):
                line = log[i]
                if 'In VibMan new if statement' in line:
                    break

            geometry_flag = False
            for j in reversed(range(i)):
                line = log[j]
                if 'Standard Nuclear Orientation' in line:
                    for line in log[(j + 3):]:
                        if '------------' not in line:
                            data = line.split()
                            QM_atom.append(data[1])
                            QM_coord.append([float(c) for c in data[2:]])
                            geometry_flag = True
                        else:
                            break
                    if geometry_flag:
                        break
            self.QM_atom = QM_atom
            self.QM_coord = np.array(QM_coord, np.float64)
        
        return coord, number, mass

    def load_conformer(self, symmetry=None, spin_multiplicity=0, optical_isomers=None, label=''):
        """
        Load the molecular degree of freedom data from an output file created as the result of a
        QChem "Freq" calculation. As QChem's guess of the external symmetry number is not always correct,
        you can use the `symmetry` parameter to substitute your own value;
        if not provided, the value in the QChem output file will be adopted.
        """
        modes = []
        freq = []
        mmass = []
        rot = []
        inertia = []
        unscaled_frequencies = []
        e0 = 0.0
        if optical_isomers is None or symmetry is None:
            _optical_isomers, _symmetry, _ = self.get_symmetry_properties()
            if optical_isomers is None:
                optical_isomers = _optical_isomers
            if symmetry is None:
                symmetry = _symmetry
        with open(self.path, 'r') as f:
            line = f.readline()
            while line != '':
                # Read spin multiplicity if not explicitly given
                if '$molecule' in line.lower() and spin_multiplicity == 0:
                    line = f.readline()
                    if len(line.split()) == 2:
                        spin_multiplicity = int(float(line.split()[1]))
                        self.charge = int(float(line.split()[0]))
                        self.spin_multiplicity = spin_multiplicity
                        logging.debug(
                            'Conformer {0} is assigned a spin multiplicity of {1}'.format(label, spin_multiplicity))
                # The rest of the data we want is in the Thermochemistry section of the output
                elif 'VIBRATIONAL ANALYSIS' in line:
                    modes = []
                    line = f.readline()
                    while line != '':

                        # This marks the end of the thermochemistry section
                        if 'Thank you very much for using Q-Chem.' in line:
                            break

                        # Read vibrational modes
                        elif 'VIBRATIONAL FREQUENCIES (CM**-1)' in line:
                            frequencies = []
                            while 'STANDARD THERMODYNAMIC QUANTITIES AT' not in line:
                                if ' Frequency:' in line:
                                    if len(line.split()) == 4:
                                        frequencies.extend([float(d) for d in line.split()[-3:]])
                                    elif len(line.split()) == 3:
                                        frequencies.extend([float(d) for d in line.split()[-2:]])
                                    elif len(line.split()) == 2:
                                        frequencies.extend([float(d) for d in line.split()[-1:]])
                                line = f.readline()
                            line = f.readline()
                            # If there is an imaginary frequency, remove it
                            if frequencies[0] < 0.0:
                                frequencies = frequencies[1:]

                            unscaled_frequencies = frequencies
                            vibration = HarmonicOscillator(frequencies=(frequencies, "cm^-1"))
                            # modes.append(vibration)
                            freq.append(vibration)
                        # Read molecular mass for external translational modes
                        elif 'Molecular Mass:' in line:
                            mass = line.split()[2]
                            if mass == '************':
                                pass
                            else:
                                mass = float(mass)
                                translation = IdealGasTranslation(mass=(mass, "amu"))
                                # modes.append(translation)
                                # mmass.append(translation)
                                mmass = [translation]

                        # Read moments of inertia for external rotational modes, given in atomic units
                        elif 'Eigenvalues --' in line:
                            if line.split()[-1] == '************************************':
                                pass
                            else:
                                inertia = [float(d) for d in line.split()[-3:]]

                        # Read the next line in the file
                        line = f.readline()

                # Read the next line in the file
                line = f.readline()

                if len(inertia):
                    if inertia[0] == 0.0:
                        # If the first eigenvalue is 0, the rotor is linear
                        inertia.remove(0.0)
                        logging.debug('inertia is {}'.format(str(inertia)))
                        for i in range(2):
                            inertia[i] *= (constants.a0 / 1e-10) ** 2
                        inertia = np.sqrt(inertia[0] * inertia[1])
                        rotation = LinearRotor(inertia=(inertia, "amu*angstrom^2"), symmetry=symmetry)
                        rot.append(rotation)
                    else:
                        for i in range(3):
                            inertia[i] *= (constants.a0 / 1e-10) ** 2
                            rotation = NonlinearRotor(inertia=(inertia, "amu*angstrom^2"), symmetry=symmetry)
                            # modes.append(rotation)
                        rot.append(rotation)

                    inertia = []

        freq = [freq[-1]] # chosse the last Vibman isotope loop
        modes = mmass + rot + freq
        return Conformer(E0=(e0 * 0.001, "kJ/mol"), modes=modes, spin_multiplicity=spin_multiplicity,
                         optical_isomers=optical_isomers), unscaled_frequencies

    def load_energy(self, zpe_scale_factor=1.):
        """
        Load the energy in J/mol from a QChem log file. Only the last energy
        in the file is returned. The zero-point energy is *not* included in
        the returned value.
        """
        e_elect = None
        with open(self.path, 'r') as f:
            a = b = c = d = 0
            for line in f:
                if 'CCSD(T) total energy' in line:
                    a = float(line.split()[4]) * constants.E_h * constants.Na                
                if 'MP2         total energy' in line:
                    b = float(line.split()[4]) * constants.E_h * constants.Na
                if 'Final energy is' in line:
                    c = float(line.split()[3]) * constants.E_h * constants.Na
                if 'Total energy in the final basis set' in line:
                    d = float(line.split()[8]) * constants.E_h * constants.Na
                e_elect = a or b or c or d
        if e_elect is None:
            raise LogError('Unable to find energy in QChem output file {0}.'.format(self.path))
        return e_elect

    def load_zero_point_energy(self):
        """
        Load the unscaled zero-point energy in J/mol from a QChem output file.
        """
        zpe = None
        with open(self.path, 'r') as f:
            for line in f:
                if 'Zero point vibrational energy:' in line:
                    zpe = float(line.split()[4]) * 4184  # QChem's ZPE is in kcal/mol, convert to J/mol
                    logging.debug('ZPE is {}'.format(str(zpe)))
        if zpe is not None:
            return zpe
        else:
            raise LogError('Unable to find zero-point energy in QChem output file {0}.'.format(self.path))

    def load_scan_energies(self):
        """
        Extract the optimized energies in J/mol from a QChem log file, e.g. the
        result of a QChem "PES Scan" quantum chemistry calculation.
        """
        v_list = []
        angle = []
        read = False
        with open(self.path, 'r') as f:
            for line in f:
                if '-----------------' in line:
                    read = False
                if read:
                    values = [float(item) for item in line.split()]
                    angle.append(values[0])
                    v_list.append(values[1])
                if 'Summary of potential scan:' in line:
                    logging.info('found a successfully completed QChem Job')
                    read = True
                elif 'SCF failed to converge' in line:
                    raise LogError('QChem Job did not successfully complete: '
                                   'SCF failed to converge in file {0}.'.format(self.path))
        logging.info('   Assuming {0} is the output from a QChem PES scan...'.format(os.path.basename(self.path)))

        v_list = np.array(v_list, np.float64)
        # check to see if the scanlog indicates that one of your reacting species may not be the lowest energy conformer
        check_conformer_energy(v_list, self.path)

        # Adjust energies to be relative to minimum energy conformer
        # Also convert units from Hartree/particle to J/mol
        v_list -= np.min(v_list)
        v_list *= constants.E_h * constants.Na
        angle = np.arange(0.0, 2 * math.pi + 0.00001, 2 * math.pi / (len(v_list) - 1), np.float64)
        return v_list, angle

    def load_negative_frequency(self):
        """
        Return the imaginary frequency from a transition state frequency
        calculation in cm^-1.
        """
        nloop = 1
        loop = 1
        frequency = 0
        with open(self.path, 'r') as f:
            for line in f:
                if 'Number of loops through sets of isotopes' in line:
                    nloop = int((line.split()[-1]))
                if 'Vibman isotope loop' in line:
                    loop = int((line.split()[-1]))
                # Read imaginary frequency
                if ' Frequency:' in line:
                    frequency = float((line.split()[1]))
                    if nloop == loop:
                        break
        # Make sure the frequency is imaginary:
        if frequency < 0:
            return frequency
        else:
            raise LogError('Unable to find imaginary frequency in QChem output file {0}.'.format(self.path))

    def load_scan_pivot_atoms(self):
        """Not implemented for QChem"""
        raise NotImplementedError('The load_scan_pivot_atoms method is not implemented for QChem Logs')

    def load_scan_frozen_atoms(self):
        """Not implemented for QChem"""
        raise NotImplementedError('The load_scan_frozen_atoms method is not implemented for QChem Logs')

register_ess_adapter("QChemLog", QChemLog)