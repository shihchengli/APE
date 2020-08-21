# -*- coding: utf-8 -*-

"""
Format of input files for Qchem software.
"""

fine = """\n   max_scf_cycles   250
   geom_opt_max_cycles   1500
   basis2 6-31G*
   scf_algorithm rca_diis
   SYM_IGNORE  TRUE
   print_input   true
   geom_opt_dmax   80
   pop_mulliken false
   XC_GRID 000075000302"""

freq_fine = """\n   max_scf_cycles   250
   SYM_IGNORE  TRUE
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
   AIMD_FIXED_ATOMS {AIMD_FIXED_ATOMS}
   geom_opt_dmax   80
   pop_mulliken false"""

input_script = """$rem
   JOBTYPE  {jobtype}
   EXCHANGE   {level_of_theory}
   BASIS   {basis}
   UNRESTRICTED   {unrestricted}{fine}
$end
{QM_atoms}{force_field_params}{opt}
$molecule
{charge} {multiplicity}
{xyz}
$end
"""