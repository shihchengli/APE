#!/usr/bin/env python
# -*- coding: utf-8 -*-

title = 'Example calculation of HCO + CH3OH = H2CO + CH2OH'

description = \
"""
This example is for a rigid rotor harmonic oscillator TST calculation calclated at the 
M08SO/MG3S level of DFT with a tight grid.  Tunneling is not included in this calculation.
"""

species('ch3oh', './ch3oh.out', protocol='UMVT')
transitionState('ts', './ts1.out', protocol='UMVT')
species('hco', './hco.out', protocol='UMVT')
'''
reaction(
    label = 'hco + ch3oh <=> hco + ch3oh',
    reactants = ['hco', 'ch3oh'],
    products = ['hco', 'ch3oh'], #product channels are only important if tunneling is computed
    transitionState = 'ts',
#    tunneling='Eckart', #we dont want to comput tunneling for this calculation
)

statmech('ts')
kinetics('hco + ch3oh <=> hco + ch3oh')
'''