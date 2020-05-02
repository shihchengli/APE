# -*- coding: utf-8 -*-
import os
from parallel_ape.submit import submit_scripts
from parallel_ape.PBS import submit_job

class ParallelJob(object):
    def __init__(self, job_path, input_file, ncpus, protocol, imaginary_bonds=''):
        self.job_path = job_path
        self.input_file = input_file
        self.ncpus = ncpus
        self.protocol = protocol
        self.imaginary_bonds = imaginary_bonds # Ex: '-i 3-12,25-18'
    
    def write_submit_file(self, submit_filename, sampling_mode):
        """
        Write the Job's submit script.
        """
        submit = submit_scripts['parallel_ape'].format(input_file=self.input_file, ncpus=self.ncpus, protocol=self.protocol, imaginary_bonds=self.imaginary_bonds, sampling_mode=sampling_mode)
        
        if not os.path.isdir(self.job_path):
            os.makedirs(self.job_path)
        with  open(os.path.join(self.job_path, submit_filename), 'w') as f:
            f.write(submit)
    
    def submit(self, submit_filename):
        job_status, job_id = submit_job(submit_filename, remote_path=self.job_path)
        return job_status, job_id
        
