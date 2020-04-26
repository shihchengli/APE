# -*- coding: utf-8 -*-

# [1] https://doi.org/10.1021/acs.jctc.5b01177 Thermodynamics of Anharmonic Systems: Uncoupled Mode Approximations for Molecules

"""
APE's main module.
To run APE through its API, first, first get a freq output file from Qchem, then call the .execute. For example:

    ape = APE(input_file)
    ape.execute()

"""

import os
import csv
import math
import numpy as np
import subprocess
import copy

import rmgpy.constants as constants

from arkane.common import symbol_by_number
from arkane.statmech import StatMechJob, determine_rotor_symmetry, is_linear

from arc.species.species import ARCSpecies

from ape.qchem import QChemLog
from ape.job import Job, record_script
from ape.torsion import HinderedRotor
from ape.InternalCoordinates import get_RedundantCoords, getXYZ
from ape.FitPES import cubic_spline_interpolations
from ape.calcThermo import ThermoJob
from ape.exceptions import InputError


class APE(object):
    """
    The main APE class.
    """

    def __init__(self,input_file, name=None, project_directory=None, protocol=None, multiplicity=None, charge = None,\
     external_symmetry=None, level_of_theory=None, basis=None, ncpus=None, imaginary_bonds=None):
        self.input_file = input_file
        self.name = name
        self.project_directory = project_directory if project_directory is not None\
            else os.path.abspath(os.path.dirname(input_file))
        self.protocol = protocol
        self.multiplicity = multiplicity
        self.charge = charge
        self.external_symmetry = external_symmetry
        self.level_of_theory = level_of_theory
        self.basis = basis
        self.ncpus = ncpus
        self.imaginary_bonds = imaginary_bonds

    def parse(self):
        Log = QChemLog(self.input_file)
        self.hessian = Log.load_force_constant_matrix()
        coordinates, number, mass = Log.load_geometry()
        self.conformer, unscaled_frequencies = Log.load_conformer()
        if self.name is None:
            self.name = self.input_file.split('/')[-1].split('.')[0]
        if self.multiplicity is None:
            self.multiplicity = Log.multiplicity
            #self.multiplicity = self.ARCSpecies.multiplicity
        if self.charge is None:
            self.charge = Log.charge
            #self.charge = 0

        self.is_QM_MM_INTERFACE = Log.is_QM_MM_INTERFACE()
        if self.is_QM_MM_INTERFACE:
            if self.imaginary_bonds is None:
                raise InputError('Lack of specified imaginary bonds to describe frustrated translation and rotation in QMMM system.')
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
            for i in range(len(self.QM_ATOMS)):
                if self.QM_USER_CONNECT[i].endswith('0  0  0  0'):
                    xyz += '{}\t{}\t\t{}\t\t{}'.format(self.symbols[i],self.cart_coords[3*i],self.cart_coords[3*i+1],self.cart_coords[3*i+2])
                    if i != self.natom-1: xyz += '\n'
            self.xyz = xyz
            self.ARCSpecies = ARCSpecies(label=self.name,xyz=self.xyz)
            if self.ncpus is None:
                self.ncpus = 8 # default cpu for QM/MM calculation
        else:
            self.natom = Log.get_number_of_atoms()
            self.symbols = [symbol_by_number[i] for i in number]
            self.cart_coords = coordinates.reshape(-1,)
            self.conformer.coordinates = (coordinates, "angstroms")
            self.conformer.mass = (mass, "amu")            
            self.xyz = getXYZ(self.symbols, self.cart_coords)
            self.ARCSpecies = ARCSpecies(label=self.name,xyz=self.xyz)
            if self.ncpus is None:
                self.ncpus = self.ARCSpecies.number_of_heavy_atoms
                if self.ncpus > 8: self.ncpus = 8
            # Below is information for thermodynamic caculation
            self.linearity = is_linear(self.conformer.coordinates.value)
            self.e_elect = Log.load_energy()
            self.zpe = Log.load_zero_point_energy()
            e0 = self.e_elect + self.zpe
            self.conformer.E0 = (e0, "J/mol")

        # Determine hindered rotors information
        if self.protocol == 'UMVT':
            self.rotors_dict = self.get_rotors_dict()
            self.n_rotors = len(self.rotors_dict)
        else:
            self.rotors_dict = []
            self.n_rotors = 0

        if self.is_QM_MM_INTERFACE:
            self.nmode = 3 * len(Log.get_QM_ATOMS())
            self.n_vib = 3 * len(Log.get_QM_ATOMS()) - self.n_rotors
        else:        
            self.nmode = 3 * self.natom - (5 if self.linearity else 6) 
            self.n_vib = 3 * self.natom - (5 if self.linearity else 6) - self.n_rotors
    
    def get_rotors_dict(self):
        rotors_dict = {}
        species = self.ARCSpecies
        species.determine_rotors()
        for i in species.rotors_dict:
            rotors_dict[i+1] = {}
            pivots = species.rotors_dict[i]['pivots']
            top = species.rotors_dict[i]['top']
            scan = species.rotors_dict[i]['scan']
            rotors_dict[i+1]['pivots'] = pivots 
            rotors_dict[i+1]['top'] = top
            rotors_dict[i+1]['scan'] = scan
        return rotors_dict
    
    def get_displaced_geometries_dict(self, vector, step_size, limit=15, torsion_ind=None):
        if torsion_ind is None:
            internal = get_RedundantCoords(self.symbols, self.cart_coords, imaginary_bonds=self.imaginary_bonds)
            magnitude = np.linalg.norm(vector)
            normalizes_vector = vector/magnitude
            qj = np.matmul(internal.B, normalizes_vector)
            x_dict = {0:self.cart_coords}

            positive_samples = range(limit)
            internal = get_RedundantCoords(self.symbols, self.cart_coords, imaginary_bonds=self.imaginary_bonds)
            if self.is_QM_MM_INTERFACE:
                internal.nHcap = self.nHcap
            for sample in positive_samples:
                #print('direction = 1')
                #print('ngrid =',sample+1)
                x_dict[sample+1] = x_dict[sample] + internal.transform_int_step((qj*step_size).reshape(-1,))
            
            negative_samples = list(range(-limit+1, 1))
            negative_samples.reverse()
            internal = get_RedundantCoords(self.symbols, self.cart_coords, imaginary_bonds=self.imaginary_bonds)
            if self.is_QM_MM_INTERFACE:
                internal.nHcap = self.nHcap
            for sample in negative_samples:
                #print('direction = -1')
                #print('ngrid =',-(sample-1))            
                x_dict[sample-1] = x_dict[sample] + internal.transform_int_step((-qj*step_size).reshape(-1,))
        else:
            rotors_dict = self.rotors_dict
            internal = get_RedundantCoords(self.symbols, self.cart_coords, rotors_dict, imaginary_bonds=self.imaginary_bonds)
            B = internal.B
            Bt_inv = np.linalg.pinv(B.dot(B.T)).dot(B)
            nrow = B.shape[0]
            qk = np.zeros(nrow, dtype=int)
            qk[torsion_ind] = 1
            x_dict = {0:self.cart_coords}
            limit = int(2*np.pi/step_size)
            positive_samples = range(limit)
            if self.is_QM_MM_INTERFACE:
                internal.nHcap = self.nHcap
            for sample in positive_samples:
                #print('direction = 1')
                #print('ngrid =',sample+1)
                x_dict[sample+1] = x_dict[sample] + internal.transform_int_step((qk*step_size).reshape(-1,))

        displaced_geometries_dict = {}
        for xi in x_dict:
            xyz = getXYZ(self.symbols, x_dict[xi])
            displaced_geometries_dict[xi] = xyz
        return displaced_geometries_dict

    def get_e_elect(self, xyz, path, file_name):
        if self.is_QM_MM_INTERFACE:
            QMMM_xyz_string = ''
            for i, xyz in enumerate(xyz.split('\n')):
                QMMM_xyz_string += " ".join([xyz, self.QM_USER_CONNECT[i]]) + '\n'
                if i == len(self.QM_ATOMS)-1:
                    break
            QMMM_xyz_string += self.fixed_molecule_string
            job = Job(QMMM_xyz_string, path, file_name,jobtype='sp', cpus=self.ncpus, charge=self.charge, multiplicity=self.multiplicity, level_of_theory=self.level_of_theory, basis=self.basis, QM_atoms=self.QM_ATOMS, force_field_params=self.force_field_params, opt=self.opt, number_of_fixed_atoms=self.number_of_fixed_atoms)
        else:
            job = Job(xyz, path, file_name,jobtype='sp', cpus=self.ncpus, charge=self.charge, multiplicity=self.multiplicity, level_of_theory=self.level_of_theory, basis=self.basis)
        job.write_input_file()
        job.submit()
        output_file_path = os.path.join(path, '{}.q.out'.format(file_name))
        e_elect = QChemLog(output_file_path).load_energy() / (constants.E_h * constants.Na) # in Hartree/particle
        return e_elect

    def sampling(self, thresh=0.05, save_result=True):
        xyz_dict = {}
        energy_dict = {}
        mode_dict = {}
        if not os.path.exists(self.project_directory):
            os.makedirs(self.project_directory)
        path = os.path.join(self.project_directory,'output_file')
        if not os.path.exists(path):
            os.makedirs(path)
        if self.protocol == 'UMVT' and self.n_rotors != 0:
            n_vib = self.n_vib
            if self.is_QM_MM_INTERFACE:
                n_vib -= 6 # due to frustrated translation and rotation 
            rotor = HinderedRotor(self.symbols, self.cart_coords, self.hessian, self.rotors_dict, self.conformer.mass.value_si, n_vib, self.imaginary_bonds)
            projected_hessian = rotor.projectd_hessian()          
            vib_freq, unweighted_v = SolvEig(projected_hessian, self.conformer.mass.value_si, self.n_vib)
            print('Frequencies(cm-1) from projected Hessian:',vib_freq)
            
            for i in range(self.n_rotors):
                xyz_dict[i+1] = {}
                energy_dict[i+1] = {}
                mode_dict[i+1] = {}
                pivots = self.rotors_dict[i+1]['pivots']
                top = self.rotors_dict[i+1]['top']
                scan = self.rotors_dict[i+1]['scan']
                scan_res = 10
                step_size = math.pi / (180/scan_res)
                #print('Sampling Mode ', (i+1))
                rotors_dict = self.rotors_dict
                internal = get_RedundantCoords(self.symbols, self.cart_coords, rotors_dict, imaginary_bonds=self.imaginary_bonds)
                scan_indices = internal.B_indices[-self.n_rotors:]
                torsion_ind = len(internal.B_indices) - self.n_rotors + scan_indices.index([ind-1 for ind in scan])
                displaced_geometries_dict = self.get_displaced_geometries_dict(vector=unweighted_v[i-self.n_rotors], step_size=step_size, torsion_ind=torsion_ind)
                mode_dict[i+1]['mode'] = 'tors'
                mode_dict[i+1]['M'] = self.conformer.get_internal_reduced_moment_of_inertia(pivots,top) * constants.Na * 1e23 # in amu*angstrom^2
                mode_dict[i+1]['K'] = (rotor.get_projected_out_freq(scan) * (2 * math.pi * constants.c * 100)) ** 2 # in 1/s^2
                mode_dict[i+1]['step_size'] = step_size # in radian                
                limit = int(360/scan_res)
                positive_samples = range(limit+1)
                for sample in positive_samples:
                    xyz = displaced_geometries_dict[sample]
                    file_name = 'tors_{}_{}'.format(i+1,sample)
                    e_elec = self.get_e_elect(xyz, path, file_name)
                    xyz_dict[i+1][sample] = xyz
                    if sample == 0:
                        energy_dict[i+1][sample] = 0
                        min_elect = e_elec
                    else: energy_dict[i+1][sample] = e_elec - min_elect
                    if e_elec - min_elect > thresh:
                        break
                v_list = [i * (constants.E_h * constants.Na) for i in energy_dict[i+1].values()] # in J/mol
                symmetry_number = determine_rotor_symmetry(v_list, self.name, scan)
                mode_dict[i+1]['symmetry_number'] = symmetry_number
        
        elif self.protocol == 'UMN':
            hessian = self.hessian
            mass = self.conformer.mass.value_si
            vib_freq, unweighted_v = SolvEig(hessian, mass, self.n_vib)
            print('Vibrational frequencies of normal modes: ',vib_freq)

        for i in range(self.nmode):
            if i in range(self.n_rotors): continue
            xyz_dict[i+1] = {}
            energy_dict[i+1] = {}
            mode_dict[i+1] = {}
            magnitude = np.linalg.norm(unweighted_v[i-self.n_rotors])
            reduced_mass = magnitude ** -2 / constants.amu # in amu
            step_size = np.sqrt(constants.hbar / (reduced_mass * constants.amu) / (vib_freq[i-self.n_rotors] * 2 * math.pi * constants.c * 100)) * 10 ** 10 # in angstrom
            #print('Sampling Mode ', (i+1))
            displaced_geometries_dict = self.get_displaced_geometries_dict(vector=unweighted_v[i-self.n_rotors], step_size=step_size)
            mode_dict[i+1]['mode'] = 'vib'
            mode_dict[i+1]['M'] = reduced_mass # in amu
            mode_dict[i+1]['K'] = (vib_freq[i-self.n_rotors] * (2 * math.pi * constants.c * 100)) ** 2 # in 1/s^2
            mode_dict[i+1]['step_size'] = step_size # in angstrom
            limit = (len(displaced_geometries_dict)-1)//2
            positive_samples = range(limit+1)
            for sample in positive_samples:
                xyz = displaced_geometries_dict[sample]
                file_name = 'vib_{}_{}'.format(i+1,sample)
                e_elec = self.get_e_elect(xyz, path, file_name)
                xyz_dict[i+1][sample] = xyz
                if sample == 0:
                    energy_dict[i+1][sample] = 0
                    min_elect = e_elec
                else: energy_dict[i+1][sample] = e_elec - min_elect
                if e_elec - min_elect > thresh:
                    break
            negative_samples = list(range(-limit+1, 1))
            negative_samples.reverse()
            for sample in negative_samples:
                sample -= 1
                xyz = displaced_geometries_dict[sample]
                file_name = 'vib_{}_{}'.format(i+1,sample)
                e_elec = self.get_e_elect(xyz, path, file_name)
                xyz_dict[i+1][sample] = xyz
                energy_dict[i+1][sample] = e_elec - min_elect
                if e_elec - min_elect > thresh:
                    break
        proc = subprocess.Popen(['rm {path}'.format(path=os.path.join(path,'input.qcin'))],shell=True)

        if save_result:
            path = self.project_directory
            self.write_samping_result_to_csv_file(path, mode_dict, energy_dict)

            path = os.path.join(path, 'plot')
            if not os.path.exists(path):
                os.makedirs(path)
            self.write_sampling_displaced_geometries(path, energy_dict, xyz_dict)

        return xyz_dict, energy_dict, mode_dict

    def write_samping_result_to_csv_file(self, path, mode_dict, energy_dict):
        csv_path = os.path.join(path, 'samping_result.csv')
        with open(csv_path, 'w') as f:
            writer = csv.writer(f)
            for i in range(self.nmode):
                is_tors = False
                if i in range(self.n_rotors): 
                    name = 'mode_{}_tors'.format(i+1)
                    is_tors = True
                else: name = 'mode_{}_vib'.format(i+1)
                writer = csv.writer(f)
                writer.writerow([name])
                if is_tors:
                    writer.writerow(['symmetry_number', mode_dict[i+1]['symmetry_number']])
                writer.writerow(['M', mode_dict[i+1]['M']])
                writer.writerow(['K', mode_dict[i+1]['K']])
                writer.writerow(['step_size', mode_dict[i+1]['step_size']])
                writer.writerow(['sample', 'total energy(HARTREE)'])
                for sample in sorted(energy_dict[i+1].keys()):
                    writer.writerow([sample, energy_dict[i+1][sample]])
            f.close()
            #print('Have saved the sampling result in {path}'.format(path=csv_path))
    
    def write_sampling_displaced_geometries(self, path, energy_dict, xyz_dict):
        # creat a format can be read by VMD software
        for i in range(self.nmode):
            txt_path = os.path.join(path, 'mode_{}.txt'.format(i+1))
            with open(txt_path, 'w') as f:
                for sample in sorted(energy_dict[i+1].keys()): 
                    content = record_script.format(natom=self.natom, sample=sample, e_elect=energy_dict[i+1][sample], xyz=xyz_dict[i+1][sample])
                    f.write(content)
                f.close()

    def execute(self):
        """
        Execute APE.
        """
        self.parse()
        xyz_dict, energy_dict, mode_dict = self.sampling()
        # Solve SE of 1-D PES and calculate E S G Cp
        polynomial_dict = cubic_spline_interpolations(energy_dict,mode_dict)
        thermo = ThermoJob(self.conformer, polynomial_dict, mode_dict, energy_dict, T=298.15, P=100000)
        if self.is_QM_MM_INTERFACE:
            thermo.calcQMMMThermo()
        else:
            thermo.calcThermo(print_HOhf_result=True, zpe_of_Hohf=self.zpe)
        

###################################################################################
def SolvEig(hessian, mass, n_vib):
    # Generate mass-weighted force constant matrix
    mass_3N_array = np.array([i for i in mass for j in range(3)])
    mass_mat = np.diag(mass_3N_array)
    inv_sq_mass_mat = np.linalg.inv(mass_mat**0.5)
    mass_weighted_hessian = inv_sq_mass_mat.dot(hessian).dot(inv_sq_mass_mat)
    eig, v = np.linalg.eigh(mass_weighted_hessian)
    vib_freq = np.sqrt(eig[-n_vib:]) / (2 * math.pi * constants.c * 100) # in cm^-1
    unweighted_v = np.matmul(inv_sq_mass_mat,v).T[-n_vib:]
    return vib_freq, unweighted_v

