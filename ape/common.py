# -*- coding: utf-8 -*-

"""
APE common module
"""

import os
import math
import copy
import logging
import numpy as np
import rmgpy.constants as constants
from arkane.statmech import determine_rotor_symmetry
from ape.job.job import Job
from ape.qchem import QChemLog
from ape.InternalCoordinates import getXYZ
from ape.exceptions import SamplingError

def SolvEig(mass_weighted_hessian, mass, n_vib):
    """
    The list of vibrational frequencies is returned in cm^-1.
    The directional vectors are not mormalized. 
    """
    mass_3N_array = np.array([i for i in mass for j in range(3)])
    mass_mat = np.diag(mass_3N_array)
    inv_sq_mass_mat = np.linalg.inv(mass_mat**0.5)
    eig, v = np.linalg.eigh(mass_weighted_hessian)
    # Convert eigenvalues to vibrational frequencies in cm^-1
    vib_freq = np.sqrt(eig[-n_vib:]) / (2 * np.pi * constants.c * 100)
    unweighted_v = np.matmul(inv_sq_mass_mat,v).T[-n_vib:]
    return vib_freq, unweighted_v

def mass_weighted_hessian(conformer, hessian, linear, is_ts):
    """
    For a given `conformer` with associated force constant matrix `hessian` and the 
    linearity of themolecule `linear`, project out the nonvibrational modes from the 
    force constant matrix and use this to determine the vibrational frequencies. The
    list of vibrational frequencies is returned in cm^-1.
    Refer to Gaussian whitepaper (http://gaussian.com/vib/) for procedure to calculate
    harmonic oscillator vibrational frequencies using the force constant matrix.
    """
    mass = conformer.mass.value_si
    coordinates = conformer.coordinates.value
    if linear is None:
        linear = is_linear(coordinates)
        if linear:
            logging.info('Determined species {0} to be linear.'.format(label))
    n_atoms = len(conformer.mass.value)

    # Put origin in center of mass
    xm = 0.0
    ym = 0.0
    zm = 0.0
    totmass = 0.0
    for i in range(n_atoms):
        xm += mass[i] * coordinates[i, 0]
        ym += mass[i] * coordinates[i, 1]
        zm += mass[i] * coordinates[i, 2]
        totmass += mass[i]

    xm /= totmass
    ym /= totmass
    zm /= totmass

    for i in range(n_atoms):
        coordinates[i, 0] -= xm
        coordinates[i, 1] -= ym
        coordinates[i, 2] -= zm
    # Make vector with the root of the mass in amu for each atom
    amass = np.sqrt(mass / constants.amu)

    # Rotation matrix
    inertia = conformer.get_moment_of_inertia_tensor()
    inertia_xyz = np.linalg.eigh(inertia)[1]

    external = 6
    if linear:
        external = 5

    d = np.zeros((n_atoms * 3, external), np.float64)

    # Transform the coordinates to the principal axes
    p = np.dot(coordinates, inertia_xyz)

    for i in range(n_atoms):
        # Projection vectors for translation
        d[3 * i + 0, 0] = amass[i]
        d[3 * i + 1, 1] = amass[i]
        d[3 * i + 2, 2] = amass[i]

    # Construction of the projection vectors for external rotation
    for i in range(n_atoms):
        d[3 * i, 3] = (p[i, 1] * inertia_xyz[0, 2] - p[i, 2] * inertia_xyz[0, 1]) * amass[i]
        d[3 * i + 1, 3] = (p[i, 1] * inertia_xyz[1, 2] - p[i, 2] * inertia_xyz[1, 1]) * amass[i]
        d[3 * i + 2, 3] = (p[i, 1] * inertia_xyz[2, 2] - p[i, 2] * inertia_xyz[2, 1]) * amass[i]
        d[3 * i, 4] = (p[i, 2] * inertia_xyz[0, 0] - p[i, 0] * inertia_xyz[0, 2]) * amass[i]
        d[3 * i + 1, 4] = (p[i, 2] * inertia_xyz[1, 0] - p[i, 0] * inertia_xyz[1, 2]) * amass[i]
        d[3 * i + 2, 4] = (p[i, 2] * inertia_xyz[2, 0] - p[i, 0] * inertia_xyz[2, 2]) * amass[i]
        if not linear:
            d[3 * i, 5] = (p[i, 0] * inertia_xyz[0, 1] - p[i, 1] * inertia_xyz[0, 0]) * amass[i]
            d[3 * i + 1, 5] = (p[i, 0] * inertia_xyz[1, 1] - p[i, 1] * inertia_xyz[1, 0]) * amass[i]
            d[3 * i + 2, 5] = (p[i, 0] * inertia_xyz[2, 1] - p[i, 1] * inertia_xyz[2, 0]) * amass[i]

    # Make sure projection matrix is orthonormal

    inertia = np.identity(n_atoms * 3, np.float64)

    p = np.zeros((n_atoms * 3, 3 * n_atoms + external), np.float64)

    p[:, 0:external] = d[:, 0:external]
    p[:, external:external + 3 * n_atoms] = inertia[:, 0:3 * n_atoms]

    for i in range(3 * n_atoms + external):
        norm = 0.0
        for j in range(3 * n_atoms):
            norm += p[j, i] * p[j, i]
        for j in range(3 * n_atoms):
            if norm > 1E-15:
                p[j, i] /= np.sqrt(norm)
            else:
                p[j, i] = 0.0
        for j in range(i + 1, 3 * n_atoms + external):
            proj = 0.0
            for k in range(3 * n_atoms):
                proj += p[k, i] * p[k, j]
            for k in range(3 * n_atoms):
                p[k, j] -= proj * p[k, i]

    # Order p, there will be vectors that are 0.0
    i = 0
    while i < 3 * n_atoms:
        norm = 0.0
        for j in range(3 * n_atoms):
            norm += p[j, i] * p[j, i]
        if norm < 0.5:
            p[:, i:3 * n_atoms + external - 1] = p[:, i + 1:3 * n_atoms + external]
        else:
            i += 1

    # T is the transformation vector from cartesian to internal coordinates
    T = np.zeros((n_atoms * 3, 3 * n_atoms - external), np.float64)

    T[:, 0:3 * n_atoms - external] = p[:, external:3 * n_atoms]

    # Generate mass-weighted force constant matrix
    # This converts the axes to mass-weighted Cartesian axes
    # Units of Fm are J/m^2*kg = 1/s^2
    weighted_hessian = hessian.copy()
    for i in range(n_atoms):
        for j in range(n_atoms):
            for u in range(3):
                for v in range(3):
                    weighted_hessian[3 * i + u, 3 * j + v] /= math.sqrt(mass[i] * mass[j])

    hessian_int = np.dot(T.T, np.dot(weighted_hessian, T))

    # Get eigenvalues of internal force constant matrix, V = 3N-6 * 3N-6
    eig, v = np.linalg.eigh(hessian_int)

    logging.debug('Frequencies from internal Hessian')
    for i in range(3 * n_atoms - external):
        with np.warnings.catch_warnings():
            np.warnings.filterwarnings('ignore', r'invalid value encountered in sqrt')
            logging.debug(np.sqrt(eig[i]) / (2 * math.pi * constants.c * 100))

    # Normal modes in mass weighted cartesian coordinates
    vmw = np.dot(T, v)
    eigm = np.zeros((3 * n_atoms - external, 3 * n_atoms - external), np.float64)

    for i in range(3 * n_atoms - external):
        eigm[i, i] = eig[i]

    fm = np.dot(vmw, np.dot(eigm, vmw.T))
    return fm

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

    initial_geometry = cart_coords.copy()
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
            XyzDictOfEachMode, EnergyDictOfEachMode, ModeDictOfEachMode, min_elect = sampling_along_vibration(symbols, initial_geometry, mode, internal, qk, projected_freq, \
            reduced_mass, step_size, path, thresh, ncpus, charge, multiplicity, level_of_theory, basis, unrestricted)
            break
        cart_coords += internal.transform_int_step((qk*step_size).reshape(-1,))
    if Fail_in_torsion_sampling is False:
        v_list = [i * (constants.E_h * constants.Na) for i in EnergyDictOfEachMode.values()] # in J/mol
        name = 'tors_{}'.format(mode)
        symmetry_number = determine_rotor_symmetry(v_list, name, scan)
        # symmetry_number = 3
        ModeDictOfEachMode['symmetry_number'] = symmetry_number
    return XyzDictOfEachMode, EnergyDictOfEachMode, ModeDictOfEachMode, min_elect

def sampling_along_vibration(symbols, cart_coords, mode, internal_object, internal_vector, freq, reduced_mass, step_size, path, thresh, ncpus, charge=None, multiplicity=None, level_of_theory=None, basis=None, unrestricted=None, \
is_QM_MM_INTERFACE=None, nHcap=None, QM_USER_CONNECT=None, QM_ATOMS=None, force_field_params=None, fixed_molecule_string=None, opt=None, number_of_fixed_atoms=None, max_nloop=20):
    XyzDictOfEachMode = {}
    EnergyDictOfEachMode = {}
    ModeDictOfEachMode = {}

    ModeDictOfEachMode['mode'] = 'vib'
    ModeDictOfEachMode['M'] = reduced_mass # in amu
    ModeDictOfEachMode['K'] = (freq * (2 * np.pi * constants.c * 100)) ** 2 # in 1/s^2
    ModeDictOfEachMode['step_size'] = step_size # in angstrom
    
    initial_geometry = cart_coords.copy()
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
            logging.warning('The energy of the end point is not above the cutoff value {thresh} hartree. Please increase the max_nloop value or increase step_size.'.format(thresh=thresh))
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
                logging.info('Sampling of mode {} in negative direction is terminated at the classical turning points.'.format(mode))
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
            logging.warning('The energy of the end point is not above the cutoff value {thresh} hartree. Please increase the max_nloop value or increase step_size.'.format(thresh=thresh))
            break
        cart_coords += internal.transform_int_step((-qj*step_size).reshape(-1,))
    return XyzDictOfEachMode, EnergyDictOfEachMode, ModeDictOfEachMode, min_elect

def get_e_elect(xyz, path, file_name, ncpus, charge=None, multiplicity=None, level_of_theory=None, basis=None, unrestricted=None, \
is_QM_MM_INTERFACE=None, QM_USER_CONNECT=None, QM_ATOMS=None, force_field_params=None, fixed_molecule_string=None, opt=None, number_of_fixed_atoms=None):
    # file_name = 'output'
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