# -*- coding: utf-8 -*-

"""
Format of input files for Qchem software.
"""

input_script = """$rem
   JOBTYPE  {jobtype}{fine}
$end
{gen_basis}{QM_atoms}{force_field_params}{opt}
$molecule
{charge} {multiplicity}
{xyz}
$end
"""

rem_variable_list = ['BASIS_LIN_DEP_THRESH', 'BASIS', 'EXCHANGE', 'CORRELATION', 'ECP', 'METHOD', 'PURECART', 'BASISPROJTYPE', 'BASIS2', 
                     'DIIS_PRINT', 'DIIS_SUBSPACE_SIZE', 'MAX_DIIS_CYCLES', 'MAX_SCF_CYCLES', 'SCF_ALGORITHM', 'PSEUDO_CANONICAL', 'SCF_CONVERGENCE', 
                     'SCF_FINAL_PRINT', 'SCF_GUESS', 'SCF_GUESS_MIX', 'SCF_GUESS_PRINT', 'SCF_PRINT', 'THRESH_DIIS_SWITCH', 'THRESH', 'UNRESTRICTED', 
                     'VARTHRESH', 'DFT Options', 'FAST_XC', 'INC_DFT', 'INCDFT_GRIDDIFF_THRESH', 'INCDFT_DENDIFF_THRESH', 'INCDFT_DENDIFF_VARTHRESH', 
                     'XC_GRID', 'XC_SMART_GRID', 'CFMM_ORDER', 'DIRECT_SCF', 'EPAO_WEIGHTS', 'EPAO_ITERATE', 'GRAIN', 'INCFOCK', 'INTEGRAL_2E_OPR', 
                     'INTEGRALS_BUFFER', 'LIN_K', 'METECO', 'OMEGA', 'PAO_ALGORITHM', 'PAO_METHOD', 'RI_J', 'ARI', 'RI_K', 'ARI_R0', 'ARI_R1', 
                     'AO2MO_DISK', 'CD_ALGORITHM', 'CORE_CHARACTER', 'MEM_TOTAL', 'MEM_STATIC', 'N_FROZEN_CORE', 'N_FROZEN_VIRTUAL', 
                     'PRINT_CORE_CHARACTER', 'QM_MM_INTERFACE', 'FORCE_FIELD', 'USER_CONNECT','SYMMETRY', 'SYM_IGNORE', 'QMMM_PRINT', 'QM_MM', 
                     'AIMD_FIXED_ATOMS', 'POP_MULLIKEN', 'PRINT_INPUT']

