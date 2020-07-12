#!/usr/bin/env python
# -*- coding: utf-8 -*-

title = 'Example calculation of NO + NO3 = 2NO2'

species('NO', './NO.q.out', protocol='UMVT')
species('NO3', './NO3.q.out', protocol='UMVT')
transitionState('TS', './TS.q.out', protocol='UMVT')
species('NO2', './NO2.q.out', protocol='UMVT')
"""
reaction(
    label = 'NO + NO3 <=> 2NO2',
    reactants = ['NO', 'NO3'],
    products = ['NO2', 'NO2'],
    transitionState = 'TS',
#    tunneling='Eckart', #we dont want to comput tunneling for this calculation
)

kinetics('NO + NO3 <=> 2NO2')
"""