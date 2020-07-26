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

from arkane.statmech import is_linear, determine_rotor_symmetry

from ape.job.job import Job
from ape.qchem import QChemLog
from ape.InternalCoordinates import getXYZ
from ape.exceptions import SamplingError

def diagonalize_projected_hessian(conformer, hessian, linear, n_vib, rotors=[], get_projected_out_freqs=False, label=None):
    """
    For a given `conformer` with associated force constant matrix `hessian`, lists of
    rotor information `rotors`, `pivots`, and `top1`, and the linearity of the
    molecule `linear`, project out the nonvibrational modes from the force
    constant matrix and use this to determine the vibrational frequencies. The
    list of vibrational frequencies is returned in cm^-1. The list of directional vectors,
    in cartesian coordinates, is returned.
    Refer to Gaussian whitepaper (http://gaussian.com/vib/) for procedure to calculate
    harmonic oscillator vibrational frequencies using the force constant matrix.
    """
    n_rotors = 0
    for rotor in rotors:
        if len(rotor) == 8:
            n_rotors += 2
        elif len(rotor) == 5 and isinstance(rotor[1][0], list):
            n_rotors += len(rotor[1])
        else:
            n_rotors += 1

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

    # Now we can start thinking about projecting out the internal rotations
    d_int = np.zeros((3 * n_atoms, n_rotors), np.float64)

    counter = 0
    for i, rotor in enumerate(rotors):
        if len(rotor) == 5 and isinstance(rotor[1][0], list):
            scan_dir, pivots_list, tops, sigmas, semiclassical = rotor
        elif len(rotor) == 5:
            scanLog, pivots, top, symmetry, fit = rotor
            pivots_list = [pivots]
            tops = [top]
        elif len(rotor) == 2:
            pivots, top = rotor
            pivots_list = [pivots]
            tops = [top]
        elif len(rotor) == 3:
            pivots, top, symmetry = rotor
            pivots_list = [pivots]
            tops = [top]
        elif len(rotor) == 8:
            scan_dir, pivots1, top1, symmetry1, pivots2, top2, symmetry2, symmetry = rotor
            pivots_list = [pivots1, pivots2]
            tops = [top1, top2]
        else:
            raise ValueError("{} not a proper rotor format".format(rotor))
        for k in range(len(tops)):
            top = tops[k]
            pivots = pivots_list[k]
            # Determine pivot atom
            if pivots[0] in top:
                pivot1 = pivots[0]
                pivot2 = pivots[1]
            elif pivots[1] in top:
                pivot1 = pivots[1]
                pivot2 = pivots[0]
            else:
                raise ValueError('Could not determine pivot atom for rotor {}.'.format(label))
            # Projection vectors for internal rotation
            e12 = coordinates[pivot1 - 1, :] - coordinates[pivot2 - 1, :]
            for j in range(n_atoms):
                atom = j + 1
                if atom in top:
                    e31 = coordinates[atom - 1, :] - coordinates[pivot1 - 1, :]
                    d_int[3 * (atom - 1):3 * (atom - 1) + 3, counter] = np.cross(e31, e12) * amass[atom - 1]
                else:
                    e31 = coordinates[atom - 1, :] - coordinates[pivot2 - 1, :]
                    d_int[3 * (atom - 1):3 * (atom - 1) + 3, counter] = np.cross(e31, -e12) * amass[atom - 1]
            counter += 1

    # Normal modes in mass weighted cartesian coordinates
    vmw = np.dot(T, v)
    eigm = np.zeros((3 * n_atoms - external, 3 * n_atoms - external), np.float64)

    for i in range(3 * n_atoms - external):
        eigm[i, i] = eig[i]

    fm = np.dot(vmw, np.dot(eigm, vmw.T))

    # Internal rotations are not normal modes => project them on the normal modes and orthogonalize
    # d_int_proj =  (3N-6) x (3N) x (3N) x (Nrotors)
    d_int_proj = np.dot(vmw.T, d_int)

    # Reconstruct d_int
    for i in range(n_rotors):
        for j in range(3 * n_atoms):
            d_int[j, i] = 0
            for k in range(3 * n_atoms - external):
                d_int[j, i] += d_int_proj[k, i] * vmw[j, k]

    # Ortho normalize
    for i in range(n_rotors):
        norm = 0.0
        for j in range(3 * n_atoms):
            norm += d_int[j, i] * d_int[j, i]
        for j in range(3 * n_atoms):
            d_int[j, i] /= np.sqrt(norm)
        for j in range(i + 1, n_rotors):
            proj = 0.0
            for k in range(3 * n_atoms):
                proj += d_int[k, i] * d_int[k, j]
            for k in range(3 * n_atoms):
                d_int[k, j] -= proj * d_int[k, i]

    # Calculate the frequencies corresponding to the internal rotors
    int_proj = np.dot(fm, d_int)
    kmus = np.array([np.linalg.norm(int_proj[:, i]) for i in range(int_proj.shape[1])])
    int_rotor_freqs = np.sqrt(kmus) / (2.0 * math.pi * constants.c * 100.0)

    if get_projected_out_freqs:
        return int_rotor_freqs

    # Do the projection
    d_int_proj = np.dot(vmw.T, d_int)
    proj = np.dot(d_int, d_int.T)
    inertia = np.identity(n_atoms * 3, np.float64)
    proj = inertia - proj
    fm = np.dot(proj, np.dot(fm, proj))
    # Get eigenvalues of mass-weighted force constant matrix
    eig, v = np.linalg.eigh(fm)
    eig.sort()

    # Convert eigenvalues to vibrational frequencies in cm^-1
    # Only keep the modes that don't correspond to translation, rotation, or internal rotation

    logging.debug('Frequencies from projected Hessian')
    for i in range(3 * n_atoms):
        with np.warnings.catch_warnings():
            np.warnings.filterwarnings('ignore', r'invalid value encountered in sqrt')
            logging.debug(np.sqrt(eig[i]) / (2 * math.pi * constants.c * 100))

    # Convert eigenvalues to vibrational frequencies in cm^-1
    vib_freq = np.sqrt(eig[-n_vib:]) / (2 * np.pi * constants.c * 100)

    # Transforme directional vectors of normal modes from mass-weighted coordinates into cartesian coordinates
    mass_3N_array = np.array([i for i in mass for j in range(3)])
    mass_mat = np.diag(mass_3N_array)
    inv_sq_mass_mat = np.linalg.inv(mass_mat ** 0.5)
    unweighted_v = np.matmul(inv_sq_mass_mat, v).T[-n_vib:]
    return vib_freq, unweighted_v

def get_internal_rotation_freq(conformer, hessian, target_rotor, rotors, linear, n_vib, is_QM_MM_INTERFACE=False, get_reduced_mass=False):
    # If is QM/MM sysyem, delete six normal modes related to frustrated motion
    if is_QM_MM_INTERFACE:
        n_vib -= 6
    
    # Calculate eigenvalues and directional vectors of vibrational normal modes
    # Only keep the modes that don't correspond to translation, rotation, or internal rotation
    projected_vib_freq, projected_v = diagonalize_projected_hessian(conformer, hessian, linear, n_vib, rotors)

    # Create a list of rotors, which exclude the target rotor
    no_tar_rotors = list()
    for rotor in rotors:
        if rotor != target_rotor:
            no_tar_rotors.append(rotor)
    
    # Calculate eigenvalues and directional vectors of normal modes which include one internal rotation
    vib_freq, v = diagonalize_projected_hessian(conformer, hessian, linear, n_vib+1, no_tar_rotors)

    # Compare this two lists to get the vibrational frequency and directional vector of this internal rotation
    internal_rotation_freq = None
    for i, freqa in enumerate(vib_freq):
        match_freq = 0
        for freqb in projected_vib_freq[i:]:
            if math.isclose(freqa, freqb, abs_tol=5) is True:
                match_freq += 1
        if match_freq == 0:
            try:
                vector = v[i]
                magnitude = np.linalg.norm(vector)
                reduced_mass = magnitude ** -2 / 1.660538921e-27 # in amu
                internal_rotation_freq = freqa
                logging.info('The vibrational frequency of internal rotation whose pivot is {pivot} is {freq:.2f} cm^-1'.format(pivot=target_rotor[0], freq=freqa))
                break
            except ValueError:
                pass

    if internal_rotation_freq is None:
        raise ConvergeError('Can\'t find the frequency of the hindered rotor whose pivot is {pivot}'.format(pivot=target_rotor[0]))
    
    if get_reduced_mass:
        return internal_rotation_freq, reduced_mass

    return internal_rotation_freq

def sampling_along_torsion(symbols, cart_coords, mode, internal_object, conformer, int_freq, reduced_mass, rotors_dict, scan_res, \
        path, thresh, ncpus, charge=None, multiplicity=None, level_of_theory=None, basis=None, unrestricted=None, \
        is_QM_MM_INTERFACE=None, nHcap=None, QM_USER_CONNECT=None, QM_ATOMS=None, force_field_params=None, \
        fixed_molecule_string=None, opt=None, number_of_fixed_atoms=None):
    XyzDictOfEachMode = {}
    EnergyDictOfEachMode = {}
    ModeDictOfEachMode = {}

    # Extract rotor information
    pivots = rotors_dict[mode]['pivots']
    top = rotors_dict[mode]['top']
    scan = rotors_dict[mode]['scan']

    # Determine sampling step size
    step_size = np.pi / (180 / scan_res)

    # Save information of this mode
    ModeDictOfEachMode['mode'] = 'tors'
    ModeDictOfEachMode['M'] = conformer.get_internal_reduced_moment_of_inertia(pivots, top) * constants.Na * 1e23 # in amu*angstrom^2
    ModeDictOfEachMode['K'] = (int_freq * (2 * np.pi * constants.c * 100)) ** 2 # in 1/s^2
    ModeDictOfEachMode['step_size'] = step_size # in radian

    # Create the unit vector, qk, for displacement along the kth bond torsion
    n_rotors = len(rotors_dict)
    internal = copy.deepcopy(internal_object)
    scan_indices = internal.B_indices[-n_rotors:]
    torsion_ind = len(internal.B_indices) - n_rotors + scan_indices.index([ind-1 for ind in scan])    
    B = internal.B
    Bt_inv = np.linalg.pinv(B.dot(B.T)).dot(B)
    nrow = B.shape[0]
    qk = np.zeros(nrow, dtype=int)
    qk[torsion_ind] = 1

    # Start to sample 1-D PES
    nsample = int(360 / scan_res) + 1
    initial_geometry = cart_coords.copy()
    cart_coords = initial_geometry.copy()
    fail_in_torsion_sampling = False
    for sample in range(nsample):
        xyz = getXYZ(symbols, cart_coords)
        file_name = 'tors_{}_{}'.format(mode, sample)

        # Calculate electronic energy of each sampling point, and save the result
        if is_QM_MM_INTERFACE:
            e_elec = get_electronic_energy(xyz, path, file_name, ncpus, charge, multiplicity, level_of_theory, basis, unrestricted, \
                is_QM_MM_INTERFACE, QM_USER_CONNECT, QM_ATOMS, force_field_params, fixed_molecule_string, opt, number_of_fixed_atoms)
        else:
            e_elec = get_electronic_energy(xyz, path, file_name, ncpus, charge, multiplicity, level_of_theory, basis, unrestricted)
        XyzDictOfEachMode[sample] = xyz

        # Take the electronic energy of stationary point as reference energy, min_elect
        if sample == 0:
            EnergyDictOfEachMode[sample] = 0
            min_elect = e_elec
        else: EnergyDictOfEachMode[sample] = e_elec - min_elect

        # Check if this mode is a high-barrier hindered rotation
        if e_elec - min_elect > thresh:
            # Treat this mode as a vibrational normal mode
            fail_in_torsion_sampling = True
            logging.info('Since the torsional barrier of mode {} is higher than {} hartree. \
                This mode will use harmonic basis to construct its hamiltonian matrix.'.format(mode, thresh))
            step_size = np.sqrt(constants.hbar / (reduced_mass * constants.amu) / (int_freq * 2 * np.pi * constants.c * 100)) * 10 ** 10 # in angstrom
            XyzDictOfEachMode, EnergyDictOfEachMode, ModeDictOfEachMode, min_elect = sampling_along_vibration(symbols, initial_geometry, mode, internal, qk, \
                int_freq, reduced_mass, step_size, path, thresh, ncpus, charge, multiplicity, level_of_theory, basis, unrestricted)
            break
        
        # Update cartesian coordinate of each sampling point
        cart_coords += internal.transform_int_step((qk * step_size).reshape(-1,))
    
    # Determine the symmetry number of this internal rotation, and save the result
    if fail_in_torsion_sampling is False:
        v_list = [i * (constants.E_h * constants.Na) for i in EnergyDictOfEachMode.values()] # in J/mol
        label = 'tors_{}'.format(mode)
        #symmetry_number = determine_rotor_symmetry(v_list, label, pivots)
        symmetry_number = 3
        ModeDictOfEachMode['symmetry_number'] = symmetry_number

    return XyzDictOfEachMode, EnergyDictOfEachMode, ModeDictOfEachMode, min_elect

def sampling_along_vibration(symbols, cart_coords, mode, internal_object, internal_vector, freq, reduced_mass, step_size, path, \
        thresh, ncpus, charge=None, multiplicity=None, level_of_theory=None, basis=None, unrestricted=None, is_QM_MM_INTERFACE=None, \
        nHcap=None, QM_USER_CONNECT=None, QM_ATOMS=None, force_field_params=None, fixed_molecule_string=None, opt=None, number_of_fixed_atoms=None, max_nloop=50):
    XyzDictOfEachMode = {}
    EnergyDictOfEachMode = {}
    ModeDictOfEachMode = {}

    # Save information of this mode
    ModeDictOfEachMode['mode'] = 'vib'
    ModeDictOfEachMode['M'] = reduced_mass # in amu
    ModeDictOfEachMode['K'] = (freq * (2 * np.pi * constants.c * 100)) ** 2 # in 1/s^2
    ModeDictOfEachMode['step_size'] = step_size # in angstrom
    
    # The direction of mode j in internal coordinates
    qj = internal_vector

    # Start to sample 1-D PES
    initial_geometry = cart_coords.copy()
    cart_coords = initial_geometry.copy()
    internal = copy.deepcopy(internal_object)

    # Sample points in positive direction
    sample = 0
    while True:
        xyz = getXYZ(symbols, cart_coords)
        file_name = 'vib_{}_{}'.format(mode, sample)
        # Calculate electronic energy of each sampling point, and save the result
        if is_QM_MM_INTERFACE:
            e_elec = get_electronic_energy(xyz, path, file_name, ncpus, charge, multiplicity, level_of_theory, basis, unrestricted, is_QM_MM_INTERFACE, \
                QM_USER_CONNECT, QM_ATOMS, force_field_params, fixed_molecule_string, opt, number_of_fixed_atoms)
        else:
            e_elec = get_electronic_energy(xyz, path, file_name, ncpus, charge, multiplicity, level_of_theory, basis, unrestricted)

        # Take the electronic energy of stationary point as reference energy, min_elect
        if sample == 0:
            XyzDictOfEachMode[0] = xyz
            EnergyDictOfEachMode[sample] = 0
            min_elect = e_elec
        # The sampling of UM-N was carried out symmetrically for each mode to the classical turning points
        elif e_elec - min_elect < EnergyDictOfEachMode[sample - 1]:
            if sample == 1:
                raise SamplingError('Sampling of mode {} fails. Make sure the directional vector of this normal mode is correct.'.format(mode))
            else:
                logging.info('Sampling of mode {} in positive direction is terminated at the classical turning points.'.format(mode))
                break
        # Not include the wrong sampling point in sampling 1D-PES
        elif e_elec - min_elect > 10:
            logging.warning('The potential energy of this point is too large. Sampling of point {} in mode {} might fail.'.format(sample, mode))
            break
        else:
            XyzDictOfEachMode[sample] = xyz
            EnergyDictOfEachMode[sample] = e_elec - min_elect

        # Check if energy rises above a cutoff energy
        if e_elec - min_elect > thresh:
            break
        
        # If the number of sampling point is over max load, please check if this job is normal
        sample += 1
        if sample > max_nloop:
            logging.warning('The energy of the end point is not above the cutoff value {thresh} hartree. Please increase the max_nloop value or increase step_size.'.format(thresh=thresh))
            break
        
        # Update cartesian coordinate of each sampling point
        cart_coords += internal.transform_int_step((qj * step_size).reshape(-1,))
    
    # Sample points in negative direction
    cart_coords = initial_geometry.copy()
    internal = copy.deepcopy(internal_object)  
    cart_coords += internal.transform_int_step((-qj * step_size).reshape(-1,))
    sample = -1
    while True:
        xyz = getXYZ(symbols, cart_coords)
        file_name = 'vib_{}_{}'.format(mode, sample)
        # Calculate electronic energy of each sampling point, and save the result
        if is_QM_MM_INTERFACE:
            e_elec = get_electronic_energy(xyz, path, file_name, ncpus, charge, multiplicity, level_of_theory, basis, unrestricted, is_QM_MM_INTERFACE, \
                QM_USER_CONNECT, QM_ATOMS, force_field_params, fixed_molecule_string, opt, number_of_fixed_atoms)
        else:
            e_elec = get_electronic_energy(xyz, path, file_name, ncpus, charge, multiplicity, level_of_theory, basis, unrestricted)

        # The sampling of UM-N was carried out symmetrically for each mode to the classical turning points
        if e_elec - min_elect < EnergyDictOfEachMode[sample + 1]:
            if sample == -1:
                raise SamplingError('Sampling of mode {} fails. Make sure the directional vector of this normal mode is correct.'.format(mode))
            else:
                logging.info('Sampling of mode {} in negative direction is terminated at the classical turning points.'.format(mode))
                break
        # Not include the wrong sampling point in sampling 1D-PES
        elif e_elec - min_elect > 10:
            logging.warning('The potential energy of this point is too large. Sampling of point {} in mode {} might fail.'.format(sample, mode))
            break
        else:
            XyzDictOfEachMode[sample] = xyz
            EnergyDictOfEachMode[sample] = e_elec - min_elect

        # Check if energy rises above a cutoff energy
        if e_elec - min_elect > thresh:
            break
        
        # If the number of sampling point is over max load, please check if this job is normal
        sample -= 1
        if sample < -max_nloop:
            logging.warning('The energy of the end point is not above the cutoff value {thresh} hartree. Please increase the max_nloop value or increase step_size.'.format(thresh=thresh))
            break
        
        # Update cartesian coordinate of each sampling point
        cart_coords += internal.transform_int_step((-qj * step_size).reshape(-1,))

    return XyzDictOfEachMode, EnergyDictOfEachMode, ModeDictOfEachMode, min_elect

def get_electronic_energy(xyz, path, file_name, ncpus, charge=None, multiplicity=None, level_of_theory=None, basis=None, unrestricted=None, \
        is_QM_MM_INTERFACE=None, QM_USER_CONNECT=None, QM_ATOMS=None, force_field_params=None, fixed_molecule_string=None, opt=None, number_of_fixed_atoms=None):
    file_name = 'output'
    if is_QM_MM_INTERFACE:
        # Create geometry format of QM/MM system 
        # <Atom> <X> <Y> <Z> <MM atom type> <Bond 1> <Bond 2> <Bond 3> <Bond 4>
        # For example:
        # O 7.256000 1.298000 9.826000 -1  185  186  0 0
        # O 6.404000 1.114000 12.310000 -1  186  713  0 0
        # O 4.077000 1.069000 0.082000 -1  188  187  0 0
        # H 1.825000 1.405000 12.197000 -3  714  0  0 0
        # H 2.151000 1.129000 9.563000 -3  189  0  0 0
        # -----------------------------------
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
    
    # Write Q-Chem input file
    #job.write_input_file()

    # Job submission
    #job.submit()

    # Parse output file to get the calculated electronic energy
    output_file_path = os.path.join(path, '{}.q.out'.format(file_name))
    e_elect = QChemLog(output_file_path).load_energy() / (constants.E_h * constants.Na) # in Hartree/particle

    return e_elect