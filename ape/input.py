# -*- coding: utf-8 -*-

"""
This module contains functionality for parsing APE input files.
"""


import logging
import os.path

import numpy as np

from rmgpy.kinetics.model import TunnelingModel
from rmgpy.kinetics.tunneling import Wigner, Eckart
from rmgpy.statmech.conformer import Conformer
from rmgpy.statmech.rotation import LinearRotor, NonlinearRotor, KRotor, SphericalTopRotor
from rmgpy.statmech.torsion import HinderedRotor, FreeRotor
from rmgpy.statmech.translation import IdealGasTranslation
from rmgpy.statmech.vibration import HarmonicOscillator

from ape.species import Species, TransitionState
from ape.sampling import SamplingJob
from ape.qchem import QChemLog
from ape.thermo import ThermoJob
from ape.reaction import Reaction
from ape.kinetics import KineticsJob
from ape.exceptions import InputError
from ape.job.inputs import rem_variable_list

################################################################################

species_dict, transition_state_dict, reaction_dict = dict(), dict(), dict()
job_list = list()
directory, output_directory = str(), str()


def species(label, *args, **kwargs):
    """Load a species from an input file"""
    global species_dict, job_list, directory
    if label in species_dict:
        raise ValueError('Multiple occurrences of species with label {0!r}.'.format(label))
    logging.info('Loading species {0}...'.format(label))

    spec = Species(label=label)
    species_dict[label] = spec

    path = None
    if len(args) == 1:
        # The argument is a path to a conformer input file
        path = os.path.join(directory, args[0])
        spec.path = path
        job = SamplingJob(label=label, input_file=path, output_directory=output_directory)
        spec.conformer, unscaled_frequencies = QChemLog(path).load_conformer()
        logging.debug('Added species {0} to a sampling job.'.format(label))
        job_list.append(job)
    elif len(args) > 1:
        raise InputError('species {0} can only have two non-keyword argument '
                         'which should be the species label and the '
                         'path to a quantum file.'.format(spec.label))

    if len(kwargs) > 0:
        # The species parameters are given explicitly
        protocol = 'UMVT'
        E0 = None
        multiplicity = None
        charge = None
        rotors = None
        for key, value in kwargs.items():
            if key == 'protocol':
                protocol = value.upper()
            elif key == 'E0':
                E0 = value
            elif key == 'multiplicity':
                multiplicity = value
            elif key == 'charge':
                charge = value
            elif key == 'rotors':
                rotors = value
            else:
                raise TypeError('species() got an unexpected keyword argument {0!r}.'.format(key))

        spec.conformer.E0 = E0

        job.protocol = protocol
        job.multiplicity = multiplicity
        job.charge = charge
        job.rotors = rotors
    
    return spec

def transitionState(label, *args, **kwargs):
    """Load a transition state from an input file"""
    global transition_state_dict, job_list, directory
    if label in transition_state_dict:
        raise ValueError('Multiple occurrences of transition state with label {0!r}.'.format(label))
    logging.info('Loading transition state {0}...'.format(label))
    ts = TransitionState(label=label)
    transition_state_dict[label] = ts

    if len(args) == 1:
        # The argument is a path to a conformer input file
        path = os.path.join(directory, args[0])
        ts.path = path
        job = SamplingJob(label=label, input_file=path, output_directory=output_directory, is_ts=True)
        Log = QChemLog(path)
        ts.conformer, unscaled_frequencies = Log.load_conformer()
        ts.frequency = (Log.load_negative_frequency(), "cm^-1")
        job_list.append(job)

    elif len(args) == 0:
        # The species parameters are given explicitly
        E0 = None
        modes = []
        spin_multiplicity = 1
        optical_isomers = 1
        frequency = None
        for key, value in kwargs.items():
            if key == 'E0':
                E0 = value
            elif key == 'modes':
                modes = value
            elif key == 'spinMultiplicity':
                spin_multiplicity = value
            elif key == 'opticalIsomers':
                optical_isomers = value
            elif key == 'frequency':
                frequency = value
            else:
                raise TypeError('transition_state() got an unexpected keyword argument {0!r}.'.format(key))

        ts.conformer = Conformer(E0=E0, modes=modes, spin_multiplicity=spin_multiplicity,
                                 optical_isomers=optical_isomers)
        ts.frequency = frequency
    else:
        if len(args) == 0 and len(kwargs) == 0:
            raise InputError(
                'The transition_state needs to reference a quantum job file or contain kinetic information.')
        raise InputError('The transition_state can only link a quantum job or directly input information, not both.')

    if len(kwargs) > 0:
        # The species parameters are given explicitly
        protocol = 'UMVT'
        E0 = None
        rotors = None
        for key, value in kwargs.items():
            if key == 'protocol':
                protocol = value.upper()
            elif key == 'E0':
                E0 = value
            elif key == 'rotors':
                rotors = value
            else:
                raise TypeError('species() got an unexpected keyword argument {0!r}.'.format(key))
        
        if protocol == 'UMVT' and rotors is None:
            raise InputError('If the transition state is sampled by using UMVT algorithm, the rotors are needed to be specified.')
               
        job.protocol = protocol
        ts.conformer.E0 = E0
        job.rotors = rotors

    return ts

def reaction(label, reactants, products, transitionState=None, tunneling=''):
    """Load a reaction from an input file"""
    global reaction_dict, species_dict, transition_state_dict
    if label in reaction_dict:
        label = label + transitionState
        if label in reaction_dict:
            raise ValueError('Multiple occurrences of reaction with label {0!r}.'.format(label))
    logging.info('Loading reaction {0}...'.format(label))
    reactants = sorted([species_dict[spec] for spec in reactants])
    products = sorted([species_dict[spec] for spec in products])
    if transitionState:
        transitionState = transition_state_dict[transitionState]
    if transitionState and (tunneling == '' or tunneling is None):
        transitionState.tunneling = None
    elif tunneling.lower() == 'wigner':
        transitionState.tunneling = Wigner(frequency=None)
    elif tunneling.lower() == 'eckart':
        transitionState.tunneling = Eckart(frequency=None, E0_reac=None, E0_TS=None, E0_prod=None)

    elif transitionState and not isinstance(tunneling, TunnelingModel):
        raise ValueError('Unknown tunneling model {0!r}.'.format(tunneling))
    rxn = Reaction(label=label, reactants=reactants, products=products, transition_state=transitionState, output_directory=output_directory)

    if isinstance(rxn, Reaction):
        reaction_dict[label] = rxn

    return rxn

def thermo(label, Tlist=[298.15]):
    """Generate a thermo job"""
    global job_list, species_dict
    try:
        spec = species_dict[label]
    except KeyError:
        raise ValueError('Unknown species label {0!r} for thermo() job.'.format(label))
    for job in job_list:
        if job.label == label:
            input_file = job.input_file
    job = ThermoJob(label=label, input_file= input_file, output_directory=output_directory, Tlist=Tlist)
    job_list.append(job)

def kinetics(label, Tmin=None, Tmax=None, Tlist=None, Tcount=0, three_params=True):
    """Generate a kinetics job"""
    global job_list, reaction_dict
    try:
        rxn = reaction_dict[label]
    except KeyError:
        raise ValueError('Unknown reaction label {0!r} for kinetics() job.'.format(label))
    job = KineticsJob(reaction=rxn, Tmin=Tmin, Tmax=Tmax, Tcount=Tcount, Tlist=Tlist, three_params=three_params)
    job_list.append(job)


################################################################################


def load_input_file(path, output_path=None):
    """
    Load the APE input file located at `path` on disk, and return a list of
    the jobs defined in that file.
    """
    global species_dict, transition_state_dict, reaction_dict, job_list, directory, output_directory
    directory = os.path.dirname(path)
    output_directory = output_path
    # Clear module-level variables
    species_dict, transition_state_dict, reaction_dict = dict(), dict(), dict()
    job_list = []

    global_context = {'__builtins__': None}
    local_context = {
        '__builtins__': None,
        'True': True,
        'False': False,
        'range': range,
        # Statistical mechanics
        'IdealGasTranslation': IdealGasTranslation,
        'LinearRotor': LinearRotor,
        'NonlinearRotor': NonlinearRotor,
        'KRotor': KRotor,
        'SphericalTopRotor': SphericalTopRotor,
        'HarmonicOscillator': HarmonicOscillator,
        'HinderedRotor': HinderedRotor,
        'FreeRotor': FreeRotor,
        # Functions
        'reaction': reaction,
        'species': species,
        'transitionState': transitionState,
        # Jobs
        'kinetics': kinetics,
        'thermo': thermo,
    }

    with open(path, 'r') as f:
        try:
            exec(f.read(), global_context, local_context)
        except (NameError, TypeError, SyntaxError):
            logging.error('The input file {0!r} was invalid:'.format(path))
            raise
    
    gen_basis = local_context.get('gen_basis', "")
    thresh = local_context.get('cut_off_energy', 0.01)
    step_size_factor = local_context.get('step_size_factor', 1)
    nnl = local_context.get('number_of_natural_length', None)
    # coordinate_system include "Normal Mode", "E-Optimized" and "E'-Optimized"
    coordinate_system = local_context.get('coordinate_system', 'Normal Mode')
    rem_variables_dict = {}
    for key in local_context.keys():
        if key.upper() in rem_variable_list:
            rem_variables_dict[key.upper()] = local_context.get(key)
    if coordinate_system not in ["Normal Mode", "E-Optimized", "E'-Optimized"]:
        raise InputError("The value of coordinate_system should be Normal Mode, E-Optimized or E'-Optimized.")
    frequency_scale_factor = local_context.get('frequency_scale_factor', 1)

    for job in job_list:
        if isinstance(job, SamplingJob):
            job.rem_variables_dict = rem_variables_dict
            job.gen_basis = gen_basis
            job.thresh = thresh
            job.step_size_factor = step_size_factor
            job.coordinate_system = coordinate_system
            job.nnl = nnl
        if isinstance(job, ThermoJob):
            job.coordinate_system = coordinate_system
            job.frequency_scale_factor = frequency_scale_factor
        if isinstance(job, KineticsJob):
            job.reaction.frequency_scale_factor = frequency_scale_factor

    return job_list, reaction_dict, species_dict, transition_state_dict



