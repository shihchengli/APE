# -*- coding: utf-8 -*-

# [1] https://doi.org/10.1021/acs.jctc.5b01177 Thermodynamics of Anharmonic Systems: Uncoupled Mode Approximations for Molecules

import os
import csv
import logging
import numpy as np
from time import gmtime, strftime

import rmgpy.constants as constants
from rmgpy.statmech import HarmonicOscillator

from arkane.common import symbol_by_number
from arkane.statmech import is_linear

from arc.species.species import ARCSpecies

from ape.qchem import QChemLog
from ape.common import diagonalize_projected_hessian, get_internal_rotation_freq, sampling_along_torsion, sampling_along_vibration
from ape.InternalCoordinates import get_RedundantCoords, getXYZ
from ape.OptimalVibrations import OptVib
from ape.exceptions import InputError

class SamplingJob(object):
    """
    The SamplingJob class.
    """

    def __init__(self, label=None, input_file=None, output_directory=None, protocol=None, spin_multiplicity=None, charge=None, 
                 rem_variables_dict={}, gen_basis="", ncpus=None, is_ts=None, imaginary_bonds=None, rotors=None, thresh=0.01,
                 step_size_factor=1, coordinate_system='Normal Mode', nnl=None):
        self.label = label
        self.input_file = input_file
        self.output_directory = output_directory
        self.protocol = protocol
        self.spin_multiplicity = spin_multiplicity
        self.charge = charge
        self.rem_variables_dict = rem_variables_dict
        self.gen_basis = gen_basis
        self.ncpus = ncpus
        self.is_ts = is_ts
        self.imaginary_bonds = imaginary_bonds
        self.rotors = rotors
        self.thresh = thresh
        self.step_size_factor = step_size_factor
        self.coordinate_system = coordinate_system
        self.nnl = nnl

    def parse(self, save_log=True):
        """
        Parse QChem output file and crate the variables the sampling job needed.
        """
        Log = QChemLog(self.input_file)

        # Load force constant matrix
        self.hessian = Log.load_force_constant_matrix()

        # Load cartesian coordinate
        coordinates, number, mass = Log.load_geometry()

        # Create conformer class
        self.conformer, unscaled_frequencies = Log.load_conformer()

        # Define the sampling protocol
        if self.protocol is None:
            self.protocol = 'UMN'
        
        # Extract spin miltiplicity and number of optical isomers
        if self.spin_multiplicity is None:
            self.spin_multiplicity = self.conformer.spin_multiplicity

        # Extract net charge from QChem output file
        if self.charge is None:
            self.charge = Log.charge
        
        # Determine wheteher `UNRESTRICTED` variable should be used or not in QChem calculation
        self.unrestricted = Log.is_unrestricted()

        # Log some information related to QM/MM system
        self.is_QM_MM_INTERFACE = Log.is_QM_MM_INTERFACE()
        if self.is_QM_MM_INTERFACE:
            if self.coordinate_system != 'Normal Mode':
                raise NotImplementedError('{} not yet implemented for QMMM system.'.format(self.coordinate_system))
            self.QM_ATOMS = Log.get_QM_ATOMS()
            self.number_of_fixed_atoms = Log.get_number_of_atoms() - len(Log.get_QM_ATOMS())
            self.ISOTOPES = Log.get_ISOTOPES()
            self.nHcap = len(self.ISOTOPES)
            self.force_field_params = Log.get_force_field_params()
            self.opt = Log.get_opt()
            self.fixed_molecule_string = Log.get_fixed_molecule()
            self.QM_USER_CONNECT = Log.get_QM_USER_CONNECT()
            self.QM_mass = Log.QM_mass
            self.QM_coord = Log.QM_coord
            self.natom =  len(self.QM_ATOMS) + len(self.ISOTOPES)
            self.symbols = Log.QM_atom
            self.cart_coords = self.QM_coord.reshape(-1,)
            self.conformer.coordinates = (self.QM_coord, "angstroms")
            self.conformer.mass = (self.QM_mass, "amu")
            xyz = ''
            self.natoms_adsorbate = 0
            for i in range(len(self.QM_ATOMS)):
                if self.QM_USER_CONNECT[i].endswith('0  0  0  0'):
                    xyz += '{}\t{}\t\t{}\t\t{}'.format(self.symbols[i], self.cart_coords[3 * i], self.cart_coords[3 * i + 1], self.cart_coords[3 * i + 2])
                    self.natoms_adsorbate += 1
                    if i != self.natom-1: xyz += '\n'
            self.xyz = xyz
            if self.xyz == '':
                self.ARCSpecies = None
            else:
                self.ARCSpecies = ARCSpecies(label=self.label, xyz=self.xyz)
            if self.ncpus is None:
                # Default ncpus for QM/MM calculation
                self.ncpus = 8
            self.zpe = Log.load_zero_point_energy()
        else:
            self.nHcap = 0
            self.natoms_adsorbate = 0
            self.natom = Log.get_number_of_atoms()
            self.symbols = [symbol_by_number[i] for i in number]
            self.cart_coords = coordinates.reshape(-1,)
            self.conformer.coordinates = (coordinates, "angstroms")
            self.conformer.number = number
            self.conformer.mass = (mass, "amu")            
            self.xyz = getXYZ(self.symbols, self.cart_coords)
            self.ARCSpecies = ARCSpecies(label=self.label, xyz=self.xyz)
            if self.ncpus is None:
                self.ncpus = self.ARCSpecies.number_of_heavy_atoms
                if self.ncpus > 8: self.ncpus = 8
            self.zpe = Log.load_zero_point_energy()

        # Determine whether or not the species is linear from its 3D coordinates
        self.linearity = is_linear(self.conformer.coordinates.value)

        # Determine hindered rotors information
        if self.protocol == 'UMVT':
            if self.rotors is not None:
                self.rotors_dict = self.rotors
            else:
                self.rotors_dict = self.get_rotors_dict()
            if self.rotors_dict == {}:
                logging.info('No internal rotations are found for {label}'.format(label=self.label))
            self.n_rotors = len(self.rotors_dict)
        elif self.protocol == 'UMN':
            self.rotors_dict = []
            self.n_rotors = 0
        else:
            raise InputError('The protocol of {protocol} is invalid. Please use UMVT or UMN.'.format(protocol=self.protocol))

        # Determine whether this system is QM/MM system
        if self.is_QM_MM_INTERFACE:
            self.nmode = 3 * len(Log.get_QM_ATOMS()) - (1 if self.is_ts else 0)
            self.n_vib = 3 * len(Log.get_QM_ATOMS()) - self.n_rotors - (1 if self.is_ts else 0)
        else:        
            self.nmode = 3 * self.natom - (5 if self.linearity else 6) - (1 if self.is_ts else 0)
            self.n_vib = 3 * self.natom - (5 if self.linearity else 6) - self.n_rotors - (1 if self.is_ts else 0)

        # Create RedundantCoords object
        self.internal = get_RedundantCoords(self.label, self.symbols, self.cart_coords, nHcap=self.nHcap, natoms_adsorbate=self.natoms_adsorbate, imaginary_bonds=self.imaginary_bonds, save_log=save_log)
        
        # Create RedundantCoords object for torsional mode
        if self.protocol == 'UMVT':
            self.torsion_internal = get_RedundantCoords(self.label, self.symbols, self.cart_coords, self.rotors_dict, self.nHcap, self.natoms_adsorbate, imaginary_bonds=self.imaginary_bonds, save_log=save_log)
        
        # Extract imaginary frequency from transition state
        if self.is_ts:
            self.imaginary_frequency = Log.load_negative_frequency()

        # Determine max_loop
        if self.nnl is not None:
            self.max_nloop = int(1 / self.step_size_factor) * self.nnl
        else:
            self.max_nloop = 200

    def get_rotors_dict(self):
        """
        Determine possible unique rotors in the species to be treated as hindered rotors,
        taking into account all localized structures.
        The resulting rotors are saved in {'pivots': [1, 3], 'top': [3, 7], 'scan': [2, 1, 3, 7]} format.
        """
        rotors_dict = {}
        species = self.ARCSpecies
        if species is None:
            return rotors_dict
        species.determine_rotors()
        for i in species.rotors_dict:
            rotors_dict[i + 1] = {}
            pivots = species.rotors_dict[i]['pivots']
            top = species.rotors_dict[i]['top']
            scan = species.rotors_dict[i]['scan']
            rotors_dict[i + 1]['pivots'] = pivots 
            rotors_dict[i + 1]['top'] = top
            rotors_dict[i + 1]['scan'] = scan
        return rotors_dict

    def sampling(self, thresh=0.01, save_result=True, scan_res=10):
        """
        Sample 1-D PES of each mode.
        The sampling of UM-N was carried out symmetrically for each mode to the classical turning points or
        the energy rises more than 0.01 hartree (i.e., about 26 kJ/mol) compared with the reference stationary point. 
        The sampling of UM-VT was terminated when the torsional angle has been displaced by 2π.
        The energy cutoff energy could be changed by defining the value of thresh.
        The dictionary of sampling geometries, calculated energies and mode information will be returned.
        """
        logging.info('Sampling {}...'.format(self.label))
        xyz_dict = {}
        energy_dict = {}
        mode_dict = {}
        if not os.path.exists(self.output_directory):
            os.makedirs(self.output_directory)
        path = os.path.join(self.output_directory, 'output_file', self.label)
        if not os.path.exists(path):
            os.makedirs(path)
        
        # Determine the vibrational frequency and directional vector of each vibrational normal mode
        if self.protocol == 'UMVT' and self.n_rotors != 0:
            logging.info(self.torsion_internal.get_intco_log())
            rotors = [[rotor['pivots'], rotor['top']] for rotor in self.rotors_dict.values()]
            vib_freq, unweighted_v = diagonalize_projected_hessian(self.conformer, self.hessian, self.linearity, self.n_vib, rotors, label=self.label)
            logging.debug('\nFrequencies(cm-1) from projected Hessian: {}'.format(vib_freq))
            
            # Sample points along the 1-D PES of each torsion motion
            for i in range(self.n_rotors):
                mode = i + 1
                target_rotor = rotors[i]
                int_freq = get_internal_rotation_freq(self.conformer, self.hessian, target_rotor, rotors, self.linearity, self.n_vib, is_QM_MM_INTERFACE=self.is_QM_MM_INTERFACE, label=self.label)
                if self.is_QM_MM_INTERFACE:
                    XyzDictOfEachMode, EnergyDictOfEachMode, ModeDictOfEachMode, min_elect = sampling_along_torsion(self.symbols, self.cart_coords, mode, self.torsion_internal, self.conformer,
                    int_freq, self.rotors_dict, scan_res, path, self.ncpus, self.charge, self.spin_multiplicity, self.rem_variables_dict, self.gen_basis, self.is_QM_MM_INTERFACE, 
                    self.QM_USER_CONNECT, self.QM_ATOMS, self.force_field_params, self.fixed_molecule_string, self.opt, self.number_of_fixed_atoms, self.label)
                else:
                    XyzDictOfEachMode, EnergyDictOfEachMode, ModeDictOfEachMode, min_elect = sampling_along_torsion(self.symbols, self.cart_coords, mode, self.torsion_internal, self.conformer,
                    int_freq, self.rotors_dict, scan_res, path, self.ncpus, self.charge, self.spin_multiplicity, self.rem_variables_dict, self.gen_basis, label=self.label)
                xyz_dict[mode] = XyzDictOfEachMode
                energy_dict[mode] = EnergyDictOfEachMode
                mode_dict[mode] = ModeDictOfEachMode
        
        elif self.protocol == 'UMN' or self.n_rotors == 0:
            logging.info(self.internal.get_intco_log())
            vib_freq, unweighted_v = diagonalize_projected_hessian(self.conformer, self.hessian, self.linearity, self.n_vib, label=self.label)
            logging.debug('Vibrational frequencies of normal modes')
            for vib in vib_freq:
                logging.debug(vib)
        
        # Optimizing vibrational coordinates to modulate intermode coupling
        if self.coordinate_system != 'Normal Mode':
            logging.debug('\nVibrational coordinates setting...')
            if self.protocol == 'UMN' or self.n_rotors == 0:
                rotors = []
            optvib_path = os.path.join(self.output_directory, 'output_file', self.label, 'tmp')
            optvib = OptVib(self.symbols, self.nmode, self.coordinate_system, self.cart_coords, self.conformer, self.hessian, self.linearity, self.n_vib, rotors, 
                            self.label, optvib_path, self.ncpus, self.imaginary_bonds, self.charge, self.spin_multiplicity, self.rem_variables_dict, self.gen_basis, self.nHcap)
            if not os.path.exists(optvib_path):
                os.makedirs(optvib_path)
            vib_freq, unweighted_v = optvib.get_optvib()

        # Sample points along the 1-D PES of each vibration motion
        for i in range(self.nmode):
            if i in range(self.n_rotors): continue
            mode = i + 1
            vector = unweighted_v[i - self.n_rotors]
            freq = vib_freq[i - self.n_rotors]
            magnitude = np.linalg.norm(vector)
            reduced_mass = magnitude ** -2 / constants.amu # in amu
            step_size = np.sqrt(constants.hbar / (reduced_mass * constants.amu) / (freq * 2 * np.pi * constants.c * 100)) * 10 ** 10 * self.step_size_factor # in angstrom
            normalizes_vector = vector / magnitude
            if self.internal.nHcap is not None:
                new_nHcap = self.internal.nHcap - self.nHcap
                normalizes_vector = np.concatenate((normalizes_vector, [0, 0, 0] * new_nHcap), axis=None)
            qj = np.matmul(self.internal.B, normalizes_vector)
            qj = qj.reshape(-1,)
            if self.is_QM_MM_INTERFACE:
                XyzDictOfEachMode, EnergyDictOfEachMode, ModeDictOfEachMode, min_elect = sampling_along_vibration(self.symbols, self.cart_coords, mode, self.internal, qj, freq, reduced_mass,
                step_size, path, thresh, self.ncpus, self.charge, self.spin_multiplicity, self.rem_variables_dict, self.gen_basis, self.is_QM_MM_INTERFACE,
                self.QM_USER_CONNECT, self.QM_ATOMS, self.force_field_params, self.fixed_molecule_string, self.opt, self.number_of_fixed_atoms, max_nloop=self.max_nloop)
            else:
                XyzDictOfEachMode, EnergyDictOfEachMode, ModeDictOfEachMode, min_elect = sampling_along_vibration(self.symbols, self.cart_coords, mode, self.internal, qj, freq, reduced_mass, step_size,
                path, thresh, self.ncpus, self.charge, self.spin_multiplicity, self.rem_variables_dict, self.gen_basis, max_nloop=self.max_nloop)
            xyz_dict[mode] = XyzDictOfEachMode
            energy_dict[mode] = EnergyDictOfEachMode
            mode_dict[mode] = ModeDictOfEachMode

        # Add the ground-state energy (including zero-point energy) of the conformer
        # Convert the unit from hartree/particle to J/mol
        self.e_elect = min_elect
        e0 = self.e_elect * constants.E_h * constants.Na + self.zpe
        self.conformer.E0 = (e0, "J/mol")

        if save_result:
            if os.path.exists(self.csv_path):
                os.remove(self.csv_path)
            self.write_samping_result_to_csv_file(self.csv_path, mode_dict, energy_dict)

            path = os.path.join(self.output_directory, 'sampling', self.label)
            if not os.path.exists(path):
                os.makedirs(path)
            self.write_sampling_displaced_geometries(path, energy_dict, xyz_dict)

        return xyz_dict, energy_dict, mode_dict

    def write_samping_result_to_csv_file(self, csv_path, mode_dict, energy_dict):
        write_min_elect = False
        if os.path.exists(csv_path) is False:
            write_min_elect = True

        with open(csv_path, 'a') as f:
            writer = csv.writer(f)
            if write_min_elect:
                writer.writerow(['min_elect', self.e_elect])
            for mode in mode_dict.keys():
                if mode_dict[mode]['mode'] == 'tors':
                    is_tors = True
                    name = 'mode_{}_tors'.format(mode)
                else:
                    is_tors = False
                    name = 'mode_{}_vib'.format(mode)
                writer.writerow([name])
                if is_tors:
                    writer.writerow(['symmetry_number', mode_dict[mode]['symmetry_number']])
                writer.writerow(['M', mode_dict[mode]['M']])
                writer.writerow(['K', mode_dict[mode]['K']])
                writer.writerow(['step_size', mode_dict[mode]['step_size']])
                writer.writerow(['sample', 'total energy(HARTREE)'])
                for sample in sorted(energy_dict[mode].keys()):
                    writer.writerow([sample, energy_dict[mode][sample]])
            f.close()
            # logging.debug('Have saved the sampling result in {path}'.format(path=csv_path))
    
    def write_sampling_displaced_geometries(self, path, energy_dict, xyz_dict):
        # creat a format can be read by VMD software
        for mode in energy_dict.keys():
            txt_path = os.path.join(path, 'mode_{}.txt'.format(mode))
            with open(txt_path, 'w') as f:
                for sample in sorted(energy_dict[mode].keys()): 
                    content = record_script.format(natom=self.natom, sample=sample, e_elect=energy_dict[mode][sample], xyz=xyz_dict[mode][sample])
                    f.write(content)
                current_time = strftime("%Y-%m-%d %H:%M:%S", gmtime())
                f.write('\n    This sampling was finished on:   {time}'.format(time=current_time))
                f.write("""\n=------------------------------------------------------------------------------=""")
                f.write("""\nSampling finished.""")
                f.write("""\n=------------------------------------------------------------------------------=""")
                f.close()

    def execute(self):
        """
        Execute APE.
        """
        self.csv_path = os.path.join(self.output_directory, '{}_samping_result.csv'.format(self.label))
        if os.path.exists(self.csv_path):
            os.remove(self.csv_path)
        self.parse()
        self.sampling(thresh=self.thresh)

# creat a format can be read by VMD software
record_script ='''{natom}
# Point {sample} Energy = {e_elect}
{xyz}
'''
