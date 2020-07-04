# -*- coding: utf-8 -*-

"""
APE common module
"""

import os
import copy
import logging
import numpy as np
import rmgpy.constants as constants
from arkane.statmech import determine_rotor_symmetry
from ape.job.job import Job
from ape.qchem import QChemLog
from ape.InternalCoordinates import get_RedundantCoords, getXYZ
from ape.exceptions import SamplingError

def SolvEig(hessian, mass, n_vib):
    # Generate mass-weighted force constant matrix
    mass_3N_array = np.array([i for i in mass for j in range(3)])
    mass_mat = np.diag(mass_3N_array)
    inv_sq_mass_mat = np.linalg.inv(mass_mat**0.5)
    mass_weighted_hessian = inv_sq_mass_mat.dot(hessian).dot(inv_sq_mass_mat)
    eig, v = np.linalg.eigh(mass_weighted_hessian)
    vib_freq = np.sqrt(eig[-n_vib:]) / (2 * np.pi * constants.c * 100) # in cm^-1
    unweighted_v = np.matmul(inv_sq_mass_mat,v).T[-n_vib:]
    return vib_freq, unweighted_v

def sampling_along_torsion(symbols, cart_coords, mode, internal_object, conformer, rotor, rotors_dict, scan_res, path, thresh, ncpus, charge=None, multiplicity=None, level_of_theory=None, basis=None, unrestricted=None,\
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
            e_elec = get_e_elect(xyz, path, file_name, ncpus, charge, multiplicity, level_of_theory, basis, unrestricted, \
            is_QM_MM_INTERFACE, QM_USER_CONNECT, QM_ATOMS, force_field_params, fixed_molecule_string, opt, number_of_fixed_atoms)
        else:
            e_elec = get_e_elect(xyz, path, file_name, ncpus, charge, multiplicity, level_of_theory, basis, unrestricted)
        XyzDictOfEachMode[sample] = xyz
        if sample == 0:
            EnergyDictOfEachMode[sample] = 0
            min_elect = e_elec
        else: EnergyDictOfEachMode[sample] = e_elec - min_elect
        if e_elec - min_elect > thresh:
            Fail_in_torsion_sampling = True
            logging.info('Since the torsional barrier of mode {} is higher than {} hartree. \
            This mode will use harmonic basis to construct its hamiltonian matrix.'.format(mode,thresh))
            step_size = np.sqrt(constants.hbar / (reduced_mass * constants.amu) / (projected_freq * 2 * np.pi * constants.c * 100)) * 10 ** 10 # in angstrom
            XyzDictOfEachMode, EnergyDictOfEachMode, ModeDictOfEachMode, min_elect = sampling_along_vibration(symbols, cart_coords, mode, internal, qk, projected_freq, \
            reduced_mass, step_size, path, thresh, ncpus, charge, multiplicity, level_of_theory, basis, unrestricted)
            break
        cart_coords += internal.transform_int_step((qk*step_size).reshape(-1,))
    if Fail_in_torsion_sampling is False:
        v_list = [i * (constants.E_h * constants.Na) for i in EnergyDictOfEachMode.values()] # in J/mol
        name = 'tors_{}'.format(mode)
        symmetry_number = determine_rotor_symmetry(v_list, name, scan)
        #symmetry_number = 3
        ModeDictOfEachMode['symmetry_number'] = symmetry_number
    return XyzDictOfEachMode, EnergyDictOfEachMode, ModeDictOfEachMode, min_elect

def sampling_along_vibration(symbols, cart_coords, mode, internal_object, internal_vector, freq, reduced_mass, step_size, path, thresh, ncpus, charge=None, multiplicity=None, level_of_theory=None, basis=None, unrestricted=None, \
is_QM_MM_INTERFACE=None, nHcap=None, QM_USER_CONNECT=None, QM_ATOMS=None, force_field_params=None, fixed_molecule_string=None, opt=None, number_of_fixed_atoms=None, max_nloop=100):
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
            e_elec = get_e_elect(xyz, path, file_name, ncpus, charge, multiplicity, level_of_theory, basis, unrestricted, is_QM_MM_INTERFACE, \
            QM_USER_CONNECT, QM_ATOMS, force_field_params, fixed_molecule_string, opt, number_of_fixed_atoms)
        else:
            e_elec = get_e_elect(xyz, path, file_name, ncpus, charge, multiplicity, level_of_theory, basis, unrestricted)
        # Take the potential energy of stationary point as reference energy, min_elect
        if sample == 0:
            XyzDictOfEachMode[0] = xyz
            EnergyDictOfEachMode[sample] = 0
            min_elect = e_elec
        # The sampling of UM-N was carried out symmetrically for each mode to the classical turning points
        elif e_elec - min_elect < EnergyDictOfEachMode[sample-1]:
            if sample == 1:
                raise SamplingError('Sampling of mode {} fails. Make sure the directional vector of this normal mode is correct.'.format(mode))
            else:
                logging.info('Sampling of mode {} in positive direction is terminated at the classical turning points.'.format(mode))
                break
        elif e_elec - min_elect > 10:
            # Not include the wrong sampling point in sampling 1D-PES
            logging.warning('The potential energy of this point is too large. Sampling of point {} in mode {} might fail.'.format(sample, mode))
            break
        else:
            XyzDictOfEachMode[sample] = xyz
            EnergyDictOfEachMode[sample] = e_elec - min_elect
        if e_elec - min_elect > thresh:
            break
        sample += 1
        if sample > max_nloop:
            logging.warning('The energy of the end point is not above the cutoff value {thresh} hartree. Please increase the max_nloop value or increase step_size'.format(thresh=thresh))
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
            e_elec = get_e_elect(xyz, path, file_name, ncpus, charge, multiplicity, level_of_theory, basis, unrestricted, is_QM_MM_INTERFACE, \
            QM_USER_CONNECT, QM_ATOMS, force_field_params, fixed_molecule_string, opt, number_of_fixed_atoms)
        else:
            e_elec = get_e_elect(xyz, path, file_name, ncpus, charge, multiplicity, level_of_theory, basis, unrestricted)
        # The sampling of UM-N was carried out symmetrically for each mode to the classical turning points
        if e_elec - min_elect < EnergyDictOfEachMode[sample+1]:
            if sample == -1:
                raise SamplingError('Sampling of mode {} fails. Make sure the directional vector of this normal mode is correct.'.format(mode))
            else:
                logging.info('Sampling of mode {} in negative direction is terminated at the classical turning points'.format(mode))
                break
        elif e_elec - min_elect > 10:
            # Not include the wrong sampling point in sampling 1D-PES
            logging.warning('The potential energy of this point is too large. Sampling of point {} in mode {} might fail.'.format(sample, mode))
            break
        else:
            XyzDictOfEachMode[sample] = xyz
            EnergyDictOfEachMode[sample] = e_elec - min_elect
        if e_elec - min_elect > thresh:
            break
        sample -= 1
        if sample < -max_nloop:
            logging.warning('The energy of the end point is not above the cutoff value {thresh} hartree. Please increase the max_nloop value or increase step_size'.format(thresh=thresh))
            break
        cart_coords += internal.transform_int_step((-qj*step_size).reshape(-1,))
    return XyzDictOfEachMode, EnergyDictOfEachMode, ModeDictOfEachMode, min_elect

def get_e_elect(xyz, path, file_name, ncpus, charge=None, multiplicity=None, level_of_theory=None, basis=None, unrestricted=None, \
is_QM_MM_INTERFACE=None, QM_USER_CONNECT=None, QM_ATOMS=None, force_field_params=None, fixed_molecule_string=None, opt=None, number_of_fixed_atoms=None):
    #file_name = 'output'
    if is_QM_MM_INTERFACE:
        QMMM_xyz_string = ''
        for i, xyz in enumerate(xyz.split('\n')):
            QMMM_xyz_string += " ".join([xyz, QM_USER_CONNECT[i]]) + '\n'
            if i == len(QM_ATOMS)-1:
                break
        QMMM_xyz_string += fixed_molecule_string
        job = Job(QMMM_xyz_string, path, file_name,jobtype='sp', ncpus=ncpus, charge=charge, multiplicity=multiplicity, \
        level_of_theory=level_of_theory, basis=basis, unrestricted=unrestricted, QM_atoms=QM_ATOMS, \
        force_field_params=force_field_params, opt=opt, number_of_fixed_atoms=number_of_fixed_atoms)
    else:
        job = Job(xyz, path, file_name,jobtype='sp', ncpus=ncpus, charge=charge, multiplicity=multiplicity, \
        level_of_theory=level_of_theory, basis=basis, unrestricted=unrestricted)
    job.write_input_file()
    job.submit()
    output_file_path = os.path.join(path, '{}.q.out'.format(file_name))
    e_elect = QChemLog(output_file_path).load_energy() / (constants.E_h * constants.Na) # in Hartree/particle
    return e_elect