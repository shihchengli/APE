# -*- coding: utf-8 -*-

"""
The parallel running version of APE.
To run APE in parallel through its API, first get a freq output file from Qchem, then call the .execute. For example:

    ape = Parallel_APE(input_file)
    ape.execute()

"""

import os
import csv
import numpy as np
import subprocess

import rmgpy.constants as constants

from ape.main import APE, SolvEig
from ape.torsion import HinderedRotor
from ape.sampling import sampling_along_torsion, sampling_along_vibration
from ape.FitPES import from_sampling_result,cubic_spline_interpolations
from ape.calcThermo import ThermoJob
from ape.exceptions import InputError, JobError
from parallel_ape.job import ParallelJob
from parallel_ape.PBS import check_job_status, delete_job

class Parallel_APE(APE):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.csv_path = os.path.join(self.project_directory, 'samping_result.csv')
        if os.path.exists(self.csv_path):
            os.remove(self.csv_path)

    def sampling(self, thresh=0.05, save_result=True, scan_res=10, sampling_mode=None):
        xyz_dict = {}
        energy_dict = {}
        mode_dict = {}
        mode = sampling_mode
        if mode is None:
            raise InputError('No specified sampling mode in parallel APE. Please assign one specific sampling mode number')
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
            
            if mode-1 in range(self.n_rotors):
                if self.is_QM_MM_INTERFACE:
                    XyzDictOfEachMode, EnergyDictOfEachMode, ModeDictOfEachMode = sampling_along_torsion(self.symbols, self.cart_coords, mode, self.internal, self.conformer, rotor, self.rotors_dict, scan_res, path, thresh, self.ncpus, self.charge, self.multiplicity, self.level_of_theory, self.basis, \
                    self.is_QM_MM_INTERFACE, self.nHcap, self.QM_USER_CONNECT, self.QM_ATOMS, self.force_field_params, self.fixed_molecule_string, self.opt, self.number_of_fixed_atoms)
                else:
                    XyzDictOfEachMode, EnergyDictOfEachMode, ModeDictOfEachMode = sampling_along_torsion(self.symbols, self.cart_coords, mode, self.internal, self.conformer, rotor, self.rotors_dict, scan_res, path, thresh, self.ncpus, self.charge, self.multiplicity, self.level_of_theory, self.basis)
        
        elif self.protocol == 'UMN':
            vib_freq, unweighted_v = SolvEig(self.hessian, self.conformer.mass.value_si, self.n_vib)
            print('Vibrational frequencies of normal modes: ',vib_freq)
        
        if mode-1 in range(self.nmode)[self.n_rotors:]:
            vector=unweighted_v[sampling_mode-1-self.n_rotors]
            freq = vib_freq[sampling_mode-1-self.n_rotors]
            magnitude = np.linalg.norm(vector)
            reduced_mass = magnitude ** -2 / constants.amu # in amu
            step_size = np.sqrt(constants.hbar / (reduced_mass * constants.amu) / (freq * 2 * np.pi * constants.c * 100)) * 10 ** 10 # in angstrom
            normalizes_vector = vector/magnitude
            qj = np.matmul(self.internal.B, normalizes_vector)
            if self.is_QM_MM_INTERFACE:
                XyzDictOfEachMode, EnergyDictOfEachMode, ModeDictOfEachMode = sampling_along_vibration(self.symbols, self.cart_coords, mode, self.internal, qj, freq, reduced_mass, self.rotors_dict, step_size, path, thresh, self.ncpus, self.charge, self.multiplicity, self.level_of_theory, self.basis, \
                self.is_QM_MM_INTERFACE, self.nHcap, self.QM_USER_CONNECT, self.QM_ATOMS, self.force_field_params, self.fixed_molecule_string, self.opt, self.number_of_fixed_atoms, max_nloop=15)
            else:
                XyzDictOfEachMode, EnergyDictOfEachMode, ModeDictOfEachMode = sampling_along_vibration(self.symbols, self.cart_coords, mode, self.internal, qj, freq, reduced_mass, self.rotors_dict, step_size, path, thresh, self.ncpus, self.charge, self.multiplicity, self.level_of_theory, self.basis, max_nloop=15)
        xyz_dict[mode] = XyzDictOfEachMode
        energy_dict[mode] = EnergyDictOfEachMode
        mode_dict[mode] = ModeDictOfEachMode

        if save_result:
            self.write_samping_result_to_csv_file(self.csv_path, mode_dict, energy_dict)

            path = os.path.join(self.project_directory, 'plot')
            if not os.path.exists(path):
                os.makedirs(path)
            self.write_sampling_displaced_geometries(path, energy_dict, xyz_dict)
        
        return XyzDictOfEachMode, EnergyDictOfEachMode, ModeDictOfEachMode

    def run(self):
        if not os.path.exists(self.project_directory):
            os.makedirs(self.project_directory)
        job_path = os.path.join(self.project_directory,'job')
        if not os.path.exists(job_path):
            os.makedirs(job_path)
        job_status_dict = {}

        freq_output = self.input_file
        ncpus = self.ncpus
        protocol = self.protocol
        imaginary_bonds_string = ''
        if self.imaginary_bonds is not None:
            for i, bond in enumerate(self.imaginary_bonds):
                imaginary_bonds_string += str(bond[0]) + '-' + str(bond[1])
                if i != len(self.imaginary_bonds) - 1:
                    imaginary_bonds_string += ','

        job = ParallelJob(job_path=job_path,input_file=freq_output,ncpus=ncpus,protocol=protocol,imaginary_bonds=imaginary_bonds_string)               
        for i in range(self.nmode):
            mode = i +1
            submit_filename = 'mode_{}.q.job'.format(mode)
            job.write_submit_file(submit_filename, sampling_mode=mode)
            job_status, job_id = job.submit(submit_filename)
            job_status_dict[job_id] = job_status
            if job_status == 'errored':
                for ID in job_status_dict.keys():
                    delete_job(ID)
                raise JobError('The submitting job, whose id is {job_id}, is errored.'.format(job_id=job_id))
        
        # write job status csv to track job progress
        csv_job_status_path = os.path.join(self.project_directory,'job.csv')
        self.write_job_status_csv(csv_job_status_path, job_status_dict)

        # Check jobs' status till every submitted jobs are finished properly.
        # Possible status: `done`, `running`, `queue`, `errored`
        Job_finished = False
        while not Job_finished:
            for job_id in job_status_dict.keys():
                job_status = check_job_status(job_id)
                if job_status_dict[job_id] != job_status:
                    if job_status == 'done':
                        mode = sorted(job_status_dict.keys()).index(job_id) + 1
                        path = os.path.join(self.project_directory, 'plot', 'mode_{}.txt'.format(mode))
                        if not os.path.exists(path):
                            job_status == 'errored'
                    job_status_dict[job_id] = job_status
                    # update jobs-tracking csv file
                    self.write_job_status_csv(csv_job_status_path, job_status_dict)
                if job_status == 'errored':
                    for ID in job_status_dict.keys():
                        delete_job(ID)
                    raise JobError('An error appers in {job_id} job'.format(job_id=job_id))
            job_status_set = set(job_status_dict.values())
            if len(job_status_set) == 1 and  job_status_set[0] == 'done':
                Job_finished = True
        #subprocess.Popen(['rm -rf {job_path}'.format(job_path=job_path)],shell=True) # delete the job folder

    def write_job_status_csv(self, csv_path, job_status_dict):
        with open(csv_path, 'w') as f:
            writer = csv.writer(f)
            writer.writerow(['mdoe','Job_id', 'status'])
            for i, ID in enumerate(sorted(job_status_dict.keys())):
                mode = i+1
                status = job_status_dict[ID]
                writer.writerow([mode, ID, status])
            f.close()

    def execute(self):
        """
        Execute APE in parallel.
        """
        self.parse()
        self.run()
        mode_dict, energy_dict = from_sampling_result(csv_path=self.csv_path)
        # Solve SE of 1-D PES and calculate E S G Cp
        polynomial_dict = cubic_spline_interpolations(energy_dict,mode_dict)
        thermo = ThermoJob(self.conformer, polynomial_dict, mode_dict, energy_dict, T=298.15, P=100000)
        if self.is_QM_MM_INTERFACE:
            thermo.calcQMMMThermo()
        else:
            thermo.calcThermo(print_HOhf_result=True, zpe_of_Hohf=self.zpe)