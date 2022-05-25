#!/usr/bin/env python
# -*- coding: utf-8 -*

# APE variables
coordinate_system = "Normal Mode" # coordinate_system include "Normal Mode", "E-Optimized" and "E'-Optimized"
cut_off_energy = 0.05 # Use 0.05 hartrees as cut_off_energy to ensure each mode sample to 4 natural lengths in each direction
step_size_factor = 1 # step_size = step_size_factor * natural_length
number_of_natural_length = 4 # sample till number_of_natural_length * natural_length in each direction of each mode
coordinate_type = "TRIC" # coordinate type include "DLC", "HDLC" and "TRIC" 


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


species('C2H5Cl', './C2H5Cl.q.out', protocol='UMVT')
transitionState('TS', './TS.q.out', protocol='UMVT',
rotors = {}
)

species('C2H4', './C2H4.q.out', protocol='UMVT')
species('HCl', './HCl.q.out', protocol='UMVT')

reaction(
    label = 'C2H5Cl=C2H4+HCl',
    reactants = ['C2H5Cl'],
    products = ['C2H4','HCl'],
    transitionState = 'TS',
    #tunneling='Eckart', #Tunneling Correction
)

kinetics('C2H5Cl=C2H4+HCl',
    Tmin = (663.15, 'K'),
    Tmax = (683.15, 'K'),
    Tcount =10,
    three_params =False,
)

