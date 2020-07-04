# -*- coding: utf-8 -*-

title = 'Example calculation of propane thermodynamics with UMN sampling scheme'

species('propane', 'propane.q.out',
        protocol='UMN')

thermo('propane')