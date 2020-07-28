#!/usr/bin/env python
# -*- coding: utf-8 -*-

title = 'Example calculation of propane thermodynamics with UMVT sampling scheme'

level_of_theory = 'omegaB97X-D'
basis = '6-311+G(2df,2pd)'
thresh = 0.05

species('propane', 'propane.q.out',
        protocol='UMVT')

thermo('propane', Tlist=[100, 200, 298.15, 300, 400, 500, 600, 700, 800, 900, 1000])
