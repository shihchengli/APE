# -*- coding: utf-8 -*-
import os
import logging
import subprocess

from ape.qchem import QChemLog
from ape.job.inputs import fine, fine_zeolite, input_script

class Job(object):
    def __init__(self, xyz, path, file_name, jobtype, ncpus, charge=None, multiplicity=None, level_of_theory=None, basis=None, unrestricted=False, QM_atoms=None, force_field_params=None, opt=None, number_of_fixed_atoms=None):
        self.xyz = xyz
        self.path = path
        self.file_name = file_name
        self.jobtype = jobtype
        self.ncpus = ncpus
        self.charge = charge
        self.multiplicity = multiplicity
        self.level_of_theory = level_of_theory
        self.basis = basis
        self.unrestricted = unrestricted
        self.is_QM_MM_INTERFACE = False
        self.number_of_fixed_atoms = number_of_fixed_atoms

        # QMMM parameter
        if QM_atoms is not None:
            self.is_QM_MM_INTERFACE = True
            self.QM_atoms = QM_atoms
            self.force_field_params = force_field_params
            self.opt = opt

        if self.charge is None:
            self.charge = 0
        if self.multiplicity is None:
            self.multiplicity = 1
        if self.level_of_theory is None:
            self.level_of_theory = 'omegaB97X-D'
        if self.basis is None:
            if self.is_QM_MM_INTERFACE:
                self.basis = '6-311++G(3df,3pd)\n   basis2   6-31G*'
            else:
                self.basis = '6-311+G(2df,2pd)'

        self.input_path = os.path.join(self.path, '{}.qcin'.format(self.file_name))
        self.output_path = os.path.join(self.path, '{}.q.out'.format(self.file_name))

    def write_input_file(self):
        """
        Write a software-specific, job-specific input file.
        Save the file locally and also upload it to the server.
        """
        if self.is_QM_MM_INTERFACE:
            fine_string = fine_zeolite.format(AIMD_FIXED_ATOMS=self.number_of_fixed_atoms)
            QM_atoms = '\n$QM_ATOMS\n' + '\n'.join(self.QM_atoms) + '\n$end\n'
            force_field_params = '\n$force_field_params\n' + self.force_field_params + '$end\n'
            opt = '\n$opt\n' + self.opt + '$end\n'
        else:
            fine_string = fine
            QM_atoms = ''
            force_field_params = ''
            opt = ''

        if self.jobtype in {'opt', 'ts', 'sp'}:
            script = input_script.format(jobtype=self.jobtype, level_of_theory=self.level_of_theory, basis=self.basis,\
            unrestricted=self.unrestricted, fine=fine_string, QM_atoms=QM_atoms, force_field_params=force_field_params, opt=opt,\
            charge=self.charge, multiplicity=self.multiplicity, xyz=self.xyz)
        f = open(self.input_path, 'w')
        # logging.debug('self.input_path :',self.input_path))
        f.write(script)
        f.close()
    
    def submit(self):
        success = False
        if os.path.exists(self.output_path):
            log = QChemLog(self.output_path)
            success = log.job_is_finished()
        if success:
            file_name = '{}.q.out'.format(self.file_name)
            logging.info('{} exists, so this calculation is passed !'.format(file_name))
        else:
            proc = subprocess.Popen(['qchem -nt {cpus} {input_path} {output_path}'.format(cpus=self.ncpus, input_path=self.input_path, output_path=self.output_path)], shell=True)
            proc.wait()
        subprocess.Popen(['rm {input_path}'.format(input_path=self.input_path)], shell=True)
