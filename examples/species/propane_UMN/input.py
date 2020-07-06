# -*- coding: utf-8 -*-

title = 'Example calculation of propane thermodynamics with UMN sampling scheme'

species('propane', 'propane.q.out',
        protocol='UMN')

thermo('propane', Tlist=[100, 200, 298.15, 300, 400, 500, 600, 700, 800, 900, 1000])
