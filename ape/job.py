# -*- coding: utf-8 -*-
import os
import subprocess

class Job(object):
    def __init__(self, xyz, path, file_name, jobtype, cpus, charge=None, multiplicity=None, level_of_theory=None, basis=None):
        self.xyz = xyz
        self.path = path
        self.jobtype = jobtype
        self.cpus = cpus
        self.charge = charge
        self.multiplicity = multiplicity
        self.level_of_theory = level_of_theory
        self.basis = basis

        if self.cpus > 8:
            self.cpus = 8
        if self.charge is None:
            self.charge = 0
        if self.multiplicity is None:
            self.multiplicity = 1
        if self.level_of_theory is None:
            self.level_of_theory = 'omegaB97X-D'
        if self.basis is None:
            self.basis = '6-311+G(2df,2pd)'
        self.input_path = os.path.join(self.path, 'input.qcin')
        self.output_path = os.path.join(self.path, '{}.q.out'.format(file_name))

    def write_input_file(self):
        """
        Write a software-specific, job-specific input file.
        Save the file locally and also upload it to the server.
        """
        is_QM_MM_INTERFACE = False
        if is_QM_MM_INTERFACE:
            fine_string = fine_zeolite 
        else: fine_string = fine
        if self.jobtype in {'opt', 'ts', 'sp'}:
            script = input_script.format(jobtype=self.jobtype, level_of_theory=self.level_of_theory, basis=self.basis,\
            fine=fine, charge=self.charge, multiplicity=self.multiplicity, xyz=self.xyz)
        f = open(self.input_path, 'w')
        f.write(script)
        f.close()
    
    def submit(self):
        if os.path.exists(self.output_path):
            print('{} exists, so this calculation is passed !'.format(self.output_path))
            pass
        else:
            proc = subprocess.Popen(['qchem -nt {cpus} {input_path} {output_path}'.format(cpus=self.cpus,input_path=self.input_path,output_path=self.output_path)],shell=True)
            proc.wait()


###################################################################################
fine = """\n   max_scf_cycles   250
   geom_opt_max_cycles   1500
   basis2 6-31G*
   scf_algorithm rca_diis
   SYM_IGNORE  TRUE
   print_input   true
   geom_opt_dmax   80
   pop_mulliken false
   XC_GRID 000075000302"""

fine_zeolite = """\n   geom_opt_coords   0
   max_scf_cycles   250
   geom_opt_max_cycles   1500
   QM_MM_INTERFACE   Zeolite
   force_field   charmm27
   user_connect   TRUE
   symmetry   off
   sym_ignore   TRUE
   print_input   true
   qmmm_print   false
   qm_mm   TRUE
   qmmm_full_hessian   FALSE
   AIMD_FIXED_ATOMS 1422
   geom_opt_dmax   80
   pop_mulliken false
   isotopes true"""

input_script = """$rem
   JOBTYPE  {jobtype}
   EXCHANGE   {level_of_theory}
   BASIS   {basis}{fine}
$end

$molecule
{charge} {multiplicity}
{xyz}
"""

#creat a format can be read by VMD software
record_script ='''{natom}
# Point {sample} Energy = {e_elect}
{xyz}
'''
