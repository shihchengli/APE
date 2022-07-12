#!/usr/bin/env python
# -*- coding: utf-8 -*-
# APE variables

coordinate_system = "Normal Mode" # coordinate_system include "Normal Mode", "E-Optimized" and "E'-Optimized"
cut_off_energy = 0.05 # Use 0.05 hartrees as cut_off_energy to ensure each mode sample to 4 natural lengths in each direction
step_size_factor = 1 # step_size = step_size_factor * natural_length
number_of_natural_length = 4 # sample till number_of_natural_length * natural_length in each direction of each mode
coordinate_type = "TRIC" # coordinate type include "RIC", "HDLC" and "TRIC" 

# QChem rem variables
exchange = 'omegaB97X-D'
basis = '6-311++G(3df,3pd)'
max_scf_cycles = 250
scf_algorithm = 'rca_diis'
sym_ignore = True
symmetry = False
print_input = True
mem_total = 20000
mem_static = 450
pop_mulliken = False
qm_mm_INTERFACE = 'Zeolite'
force_field = 'charmm27'
user_connect = True
qmmm_print = False
qm_mm = True

# It is recommended that the rotors in QM/MM system should be defined mannualy   
species('H-MFI_propane', './H-MFI_propane.q.out', protocol='UMVT',rotors={1: {'pivots': [1, 2], 'top': [1, 3, 4, 5], 'scan': [3, 1, 2, 6]}, 2: {'pivots': [2, 6], 'top': [6, 9, 10, 11], 'scan': [1, 2, 6, 9]}})

thermo('H-MFI_propane', Tlist=[298.15])