# -*- coding: utf-8 -*-
import os
import subprocess

class Job(object):
    def __init__(self, xyz, path, file_name, jobtype, cpus, charge=None, multiplicity=None, level_of_theory=None, basis=None, QM_atoms=None, force_field_params=None, opt=None):
        self.xyz = xyz
        self.path = path
        self.jobtype = jobtype
        self.cpus = cpus
        self.charge = charge
        self.multiplicity = multiplicity
        self.level_of_theory = level_of_theory
        self.basis = basis
        self.is_QM_MM_INTERFACE = False

        # QMMM parameter
        if QM_atoms is not None:
            self.is_QM_MM_INTERFACE = True
            self.QM_atoms = QM_atoms
            self.force_field_params = force_field_params
            self.opt = opt

        if self.cpus > 8:
            self.cpus = 8
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
        self.input_path = os.path.join(self.path, 'input.qcin')
        self.output_path = os.path.join(self.path, '{}.q.out'.format(file_name))

    def write_input_file(self):
        """
        Write a software-specific, job-specific input file.
        Save the file locally and also upload it to the server.
        """
        if self.is_QM_MM_INTERFACE:
            fine_string = fine_zeolite
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
            fine=fine_string, QM_atoms=QM_atoms, force_field_params=force_field_params, opt=opt, charge=self.charge, multiplicity=self.multiplicity, xyz=self.xyz)
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
   thresh 14
   scf_convergence 7   
   AIMD_FIXED_ATOMS 1422
   geom_opt_dmax   80
   pop_mulliken false"""

input_script = """$rem
   JOBTYPE  {jobtype}
   EXCHANGE   {level_of_theory}
   BASIS   {basis}{fine}
$end
{QM_atoms}{force_field_params}{opt}
$molecule
{charge} {multiplicity}
{xyz}
$end
"""

#creat a format can be read by VMD software
record_script ='''{natom}
# Point {sample} Energy = {e_elect}
{xyz}
'''
