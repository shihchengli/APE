# -*- coding: utf-8 -*-
import os
import logging
import subprocess

from ape.qchem import QChemLog
from ape.job.inputs import input_script
from ape.exceptions import JobError

class Job(object):
    def __init__(self, xyz, path, file_name, jobtype, ncpus, charge=None, multiplicity=None, rem_variables_dict=None, 
                 gen_basis="", QM_atoms=None, ISOTOPE=None, force_field_params=None, opt=None):
        self.xyz = xyz
        self.path = path
        self.file_name = file_name
        self.jobtype = jobtype
        self.ncpus = ncpus
        self.charge = charge
        self.multiplicity = multiplicity
        self.rem_variables_dict = rem_variables_dict
        self.gen_basis = gen_basis

        # QMMM parameter
        self.QM_atoms = QM_atoms
        self.ISOTOPE = ISOTOPE
        self.force_field_params = force_field_params
        self.opt = opt

        if self.charge is None:
            self.charge = 0
        if self.multiplicity is None:
            self.multiplicity = 1

        self.input_path = os.path.join(self.path, '{}.qcin'.format(self.file_name))
        self.output_path = os.path.join(self.path, '{}.q.out'.format(self.file_name))

    def write_input_file(self):
        """
        Write a software-specific, job-specific input file.
        Save the file locally and also upload it to the server.
        """
        QM_atoms, isotope_params, force_field_params, opt = "", "" , "", ""
        if self.QM_atoms is not None:
            QM_atoms = '\n$QM_ATOMS\n' + '\n'.join(self.QM_atoms) + '\n$end\n'
            force_field_params = '\n$force_field_params\n' + self.force_field_params + '$end\n'
            opt = '\n$opt\n' + self.opt + '$end\n'

        if self.ISOTOPE is not None:
            isotope_params = '\n$isotopes\n1       0\n' + str(len(self.ISOTOPE))
            for key in self.ISOTOPE:
                isotope_params += '\n{}      {}'.format(key, self.ISOTOPE[key])
            isotope_params += '\n$end\n'

        fine = ""
        for key in self.rem_variables_dict.keys():
            value = self.rem_variables_dict[key]
            fine += "\n   {0}   {1}".format(key, value)

        if self.jobtype in {'opt', 'ts', 'sp', 'freq'}:
            script = input_script.format(jobtype=self.jobtype, fine=fine, gen_basis=self.gen_basis, QM_atoms=QM_atoms, isotope=isotope_params,
            force_field_params=force_field_params, opt=opt, charge=self.charge, multiplicity=self.multiplicity, xyz=self.xyz)
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
            logging.debug('{} exists, so this calculation is passed!'.format(file_name))
        else:
            proc = subprocess.Popen(['qchem -nt {cpus} {input_path} {output_path}'.format(cpus=self.ncpus, input_path=self.input_path, output_path=self.output_path)], shell=True)
            proc.wait()
            log = QChemLog(self.output_path)
            success = log.job_is_finished()
            if not success:
                raise JobError('QChem job fails !')
        subprocess.Popen(['rm {input_path}'.format(input_path=self.input_path)], shell=True)

        # Remove efld file from QMMM calculation
        efld_path = os.path.join(self.path, '{}.q.out.efld'.format(self.file_name))
        if os.path.exists(efld_path):
            subprocess.Popen(['rm {efld_path}'.format(efld_path=efld_path)], shell=True)
