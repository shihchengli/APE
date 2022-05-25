#!/usr/bin/env python
# -*- coding: utf-8 -*-
# APE variables

coordinate_system = "Normal Mode"
cut_off_energy = 0.05
step_size_factor = 1
number_of_natural_length = 4
coordinate_type = "DLC"

# QChem rem variables
exchange = 'omegaB97X-D'
basis = '6-311+G(2df,2pd)'
max_scf_cycles = 250
scf_algorithm = 'rca_diis'
sym_ignore = True
symmetry = False
print_input = True
XC_GRID = '000075000302'
mem_total = 20000
mem_static = 450
pop_mulliken = False

species('1-propanol', './1-propanol.q.out', protocol='UMVT')


thermo('1-propanol', Tlist=[298.15])
