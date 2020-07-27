#!/usr/bin/env python
# -*- coding: utf-8 -*-

title = 'Example calculation of C2H5OH = C2H4 + H2O'

species('C2H5OH', './C2H5OH.q.out', protocol='UMVT')
transitionState('TS', './TS.q.out', protocol='UMVT',
    rotors = {1: {'pivots': [1, 2], 'top': [1, 3, 4, 5], 'scan': [3, 1, 2, 7]}, 
              2: {'pivots': [1, 5], 'top': [1, 2, 3, 4], 'scan': [2, 1, 5, 6]},
              3: {'pivots': [5, 6], 'top': [6, 9], 'scan': [1, 5, 6 ,9]}},
)
species('C2H4', './C2H4.q.out', protocol='UMVT')
species('H2O', './H2O.q.out', protocol='UMVT')

reaction(
    label = 'C2H5OH <=> C2H4 + H2O',
    reactants = ['C2H5OH'],
    products = ['C2H4', 'H2O'],
    transitionState = 'TS',
#    tunneling='Eckart', #we dont want to comput tunneling for this calculation
)

kinetics('C2H5OH <=> C2H4 + H2O',
    Tmin = (293, 'K'),
    Tmax = (750, 'K'),
    Tcount = 30,
    three_params = True,
)

kinetics('C2H5OH <=> C2H4 + H2O',
    Tmin = (1300, 'K'),
    Tmax = (1510, 'K'),
    Tcount = 10,
    three_params = False,
)

kinetics('C2H5OH <=> C2H4 + H2O',
    Tmin = (1449, 'K'),
    Tmax = (1700, 'K'),
    Tcount = 10,
    three_params = False,
)

kinetics('C2H5OH <=> C2H4 + H2O',
    Tmin = (800, 'K'),
    Tmax = (1800, 'K'),
    Tcount = 30,
    three_params = True,
)
