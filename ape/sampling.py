# -*- coding: utf-8 -*-

"""
A module to sample the geometries along given direction
"""

import os
import copy
import numpy as np
import rmgpy.constants as constants
from ape.job import Job
from ape.qchem import QChemLog
from ape.InternalCoordinates import get_RedundantCoords, getXYZ

def sampling_along_torsion(symbols, cart_coords, mode, internal_object, conformer, rotor, rotors_dict, scan_res, path, thresh, ncpus, charge=None, multiplicity=None, level_of_theory=None, basis=None, \
is_QM_MM_INTERFACE=None, nHcap=None, QM_USER_CONNECT=None, QM_ATOMS=None, force_field_params=None, fixed_molecule_string=None, opt=None, number_of_fixed_atoms=None):
    XyzDictOfEachMode = {}
    EnergyDictOfEachMode = {}
    ModeDictOfEachMode = {}

    pivots = rotors_dict[mode]['pivots']
    top = rotors_dict[mode]['top']
    scan = rotors_dict[mode]['scan']
    step_size = np.pi / (180/scan_res)
    projected_freq, reduced_mass = rotor.get_projected_out_freq(scan)

    ModeDictOfEachMode['mode'] = 'tors'
    ModeDictOfEachMode['M'] = conformer.get_internal_reduced_moment_of_inertia(pivots,top) * constants.Na * 1e23 # in amu*angstrom^2
    ModeDictOfEachMode['K'] = (projected_freq * (2 * np.pi * constants.c * 100)) ** 2 # in 1/s^2
    ModeDictOfEachMode['step_size'] = step_size # in radian

    n_rotors = len(rotors_dict)
    internal = copy.deepcopy(internal_object)
    scan_indices = internal.B_indices[-n_rotors:]
    torsion_ind = len(internal.B_indices) - n_rotors + scan_indices.index([ind-1 for ind in scan])    

    B = internal.B
    Bt_inv = np.linalg.pinv(B.dot(B.T)).dot(B)
    nrow = B.shape[0]
    qk = np.zeros(nrow, dtype=int)
    qk[torsion_ind] = 1
    nsample = int(360/scan_res) + 1

    initial_geometry = cart_coords
    cart_coords = initial_geometry.copy()
    Fail_in_torsion_sampling = False
    for sample in range(nsample):
        xyz = getXYZ(symbols, cart_coords)
        file_name = 'tors_{}_{}'.format(mode,sample)
        if is_QM_MM_INTERFACE:
            e_elec = get_e_elect(xyz, path, file_name, ncpus, charge, multiplicity, level_of_theory, basis, is_QM_MM_INTERFACE, \
            QM_USER_CONNECT, QM_ATOMS, force_field_params, fixed_molecule_string, opt, number_of_fixed_atoms)
        else:
            e_elec = get_e_elect(xyz, path, file_name, ncpus)
        XyzDictOfEachMode[sample] = xyz
        if sample == 0:
            EnergyDictOfEachMode[sample] = 0
            min_elect = e_elec
        else: EnergyDictOfEachMode[sample] = e_elec - min_elect
        if e_elec - min_elect > thresh:
            Fail_in_torsion_sampling = True
            print('Since the torsional barrier of mode {} is higher than {} hartree. \
            This mode will use harmonic basis to construct its hamiltonian matrix.'.format(mode,thresh))
            step_size = np.sqrt(constants.hbar / (reduced_mass * constants.amu) / (projected_freq * 2 * math.pi * constants.c * 100)) * 10 ** 10 # in angstrom
            XyzDictOfEachMode, EnergyDictOfEachMode, ModeDictOfEachMode = sampling_along_vibration(symbols, cart_coords, mode, internal, qk, projected_freq, reduced_mass, rotors_dict, step_size, path, thresh, ncpus, charge, multiplicity, level_of_theory, basis, max_nloop=15)
            break
        cart_coords += internal.transform_int_step((qk*step_size).reshape(-1,))
    if Fail_in_torsion_sampling is False:
        v_list = [i * (constants.E_h * constants.Na) for i in EnergyDictOfEachMode.values()] # in J/mol
        name = 'tors_{}'.format(mode)
        #symmetry_number = determine_rotor_symmetry(v_list, name, scan)
        symmetry_number = 3
        ModeDictOfEachMode['symmetry_number'] = symmetry_number
    return XyzDictOfEachMode, EnergyDictOfEachMode, ModeDictOfEachMode

def sampling_along_vibration(symbols, cart_coords, mode, internal_object, internal_vector, freq, reduced_mass, rotors_dict, step_size, path, thresh, ncpus, charge=None, multiplicity=None, level_of_theory=None, basis=None, \
is_QM_MM_INTERFACE=None, nHcap=None, QM_USER_CONNECT=None, QM_ATOMS=None, force_field_params=None, fixed_molecule_string=None, opt=None, number_of_fixed_atoms=None, max_nloop=15):
    XyzDictOfEachMode = {}
    EnergyDictOfEachMode = {}
    ModeDictOfEachMode = {}

    ModeDictOfEachMode['mode'] = 'vib'
    ModeDictOfEachMode['M'] = reduced_mass # in amu
    ModeDictOfEachMode['K'] = (freq * (2 * np.pi * constants.c * 100)) ** 2 # in 1/s^2
    ModeDictOfEachMode['step_size'] = step_size # in angstrom
    
    initial_geometry = cart_coords
    cart_coords = initial_geometry.copy()
    internal = copy.deepcopy(internal_object)
    qj = internal_vector

    sample = 0
    while True:
        xyz = getXYZ(symbols, cart_coords)
        file_name = 'vib_{}_{}'.format(mode,sample)
        if is_QM_MM_INTERFACE:
            e_elec = get_e_elect(xyz, path, file_name, ncpus, charge, multiplicity, level_of_theory, basis, is_QM_MM_INTERFACE, \
            QM_USER_CONNECT, QM_ATOMS, force_field_params, fixed_molecule_string, opt, number_of_fixed_atoms)
        else:
            e_elec = get_e_elect(xyz, path, file_name, ncpus)
        XyzDictOfEachMode[sample] = xyz
        if sample == 0:
            EnergyDictOfEachMode[sample] = 0
            min_elect = e_elec
        else: EnergyDictOfEachMode[sample] = e_elec - min_elect
        if e_elec - min_elect > thresh:
            break
        sample += 1
        if sample > max_nloop:
            break
        cart_coords += internal.transform_int_step((qj*step_size).reshape(-1,))
    
    cart_coords = initial_geometry.copy()
    internal = copy.deepcopy(internal_object)  
    cart_coords += internal.transform_int_step((-qj*step_size).reshape(-1,))
    sample = -1
    while True:
        xyz = getXYZ(symbols, cart_coords)
        file_name = 'vib_{}_{}'.format(mode,sample)
        if is_QM_MM_INTERFACE:
            e_elec = get_e_elect(xyz, path, file_name, ncpus, charge, multiplicity, level_of_theory, basis, is_QM_MM_INTERFACE, \
            QM_USER_CONNECT, QM_ATOMS, force_field_params, fixed_molecule_string, opt, number_of_fixed_atoms)
        else:
            e_elec = get_e_elect(xyz, path, file_name, ncpus)
        XyzDictOfEachMode[sample] = xyz
        EnergyDictOfEachMode[sample] = e_elec - min_elect
        if e_elec - min_elect > thresh:
            break
        sample -= 1
        if sample < -max_nloop:
            break
        cart_coords += internal.transform_int_step((-qj*step_size).reshape(-1,))
    return XyzDictOfEachMode, EnergyDictOfEachMode, ModeDictOfEachMode

def get_e_elect(xyz, path, file_name, ncpus, charge=None, multiplicity=None, level_of_theory=None, basis=None, is_QM_MM_INTERFACE=None, \
QM_USER_CONNECT=None, QM_ATOMS=None, force_field_params=None, fixed_molecule_string=None, opt=None, number_of_fixed_atoms=None):
    file_name = 'output'
    if is_QM_MM_INTERFACE:
        QMMM_xyz_string = ''
        for i, xyz in enumerate(xyz.split('\n')):
            QMMM_xyz_string += " ".join([xyz, QM_USER_CONNECT[i]]) + '\n'
            if i == len(QM_ATOMS)-1:
                break
        QMMM_xyz_string += fixed_molecule_string
        job = Job(QMMM_xyz_string, path, file_name,jobtype='sp', ncpus=ncpus, charge=charge, multiplicity=multiplicity, \
        level_of_theory=level_of_theory, basis=basis, QM_atoms=QM_ATOMS, force_field_params=force_field_params, opt=opt, \
        number_of_fixed_atoms=number_of_fixed_atoms)
    else:
        job = Job(xyz, path, file_name,jobtype='sp', ncpus=ncpus, charge=charge, multiplicity=multiplicity, \
        level_of_theory=level_of_theory, basis=basis)
    #job.write_input_file()
    job.submit()
    output_file_path = os.path.join(path, '{}.q.out'.format(file_name))
    e_elect = QChemLog(output_file_path).load_energy() / (constants.E_h * constants.Na) # in Hartree/particle
    return e_elect