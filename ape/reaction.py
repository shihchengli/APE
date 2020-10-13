#!/usr/bin/env python3

import logging
import numpy as np
import math

import rmgpy.constants as constants
from rmgpy.exceptions import ReactionError
from rmgpy.reaction import Reaction as rmg_Reaction
from rmgpy.statmech.vibration import HarmonicOscillator

from ape.thermo import ThermoJob

class Reaction(object):
    def __init__(self, label='', reactants=None, products=None, transition_state=None, output_directory=None, frequency_scale_factor=1):
        self.label = label
        self.reactants = reactants
        self.products = products
        self.transition_state = transition_state
        self.output_directory = output_directory
        self.frequency_scale_factor = frequency_scale_factor
        self.kinetics = None
        self.thermo_dict = {}
    
    def rmg_Reaction(self):
        specs = self.reactants + self.products + [self.transition_state]
        for spec in specs:
            thermo = ThermoJob(spec.label, spec.path, output_directory=self.output_directory, frequency_scale_factor=self.frequency_scale_factor)
            thermo.load_save()
            spec.conformer.E0 = thermo.conformer.E0
            for mode in spec.conformer.modes:
                if isinstance(mode, HarmonicOscillator):
                    frequencies = mode.frequencies.value_si
                    mode.frequencies = (frequencies * self.frequency_scale_factor, "cm^-1")
        rxn = rmg_Reaction(label=self.label, reactants=self.reactants, products=self.products, transition_state=self.transition_state)
        return rxn
    
    def calcThermo(self, T, P=101325):
        self.thermo_dict[T] = {}
        specs = set(self.reactants + self.products + [self.transition_state])
        for spec in specs:
            logging.debug('    Calculating thermodynamics properties for {0} at {1} K'.format(spec.label, T))
            self.thermo_dict[T][spec.label] = {}
            thermo = ThermoJob(spec.label, spec.path, output_directory=self.output_directory, P=P, frequency_scale_factor=self.frequency_scale_factor)
            thermo.load_save()
            E0, E, S, F, Q, Cv = thermo.calcThermo(T, print_HOhf_result=False)
            self.thermo_dict[T][spec.label]['E0'] = E0
            self.thermo_dict[T][spec.label]['E'] = E
            self.thermo_dict[T][spec.label]['S'] = S
            self.thermo_dict[T][spec.label]['F'] = F
            self.thermo_dict[T][spec.label]['Q'] = Q
            self.thermo_dict[T][spec.label]['Cv'] = Cv
    
    def calculate_tst_rate_coefficient(self, T, P=101325):
        """
        Evaluate the forward rate coefficient for the reaction with
        corresponding transition state `TS` at temperature `T` in K using
        (canonical) transition state theory. The TST equation is
        .. math:: k(T) = \\kappa(T) \\frac{k_\\mathrm{B} T}{h} \\frac{Q^\\ddagger(T)}{Q^\\mathrm{A}(T) Q^\\mathrm{B}(T)} \\exp \\left( -\\frac{E_0}{k_\\mathrm{B} T} \\right)
        where :math:`Q^\\ddagger` is the partition function of the transition state,
        :math:`Q^\\mathrm{A}` and :math:`Q^\\mathrm{B}` are the partition function
        of the reactants, :math:`E_0` is the ground-state energy difference from
        the transition state to the reactants, :math:`T` is the absolute
        temperature, :math:`k_\\mathrm{B}` is the Boltzmann constant, and :math:`h`
        is the Planck constant. :math:`\\kappa(T)` is an optional tunneling
        correction.
        """
        # Calculate partition function of reactants and transition state
        if T not in self.thermo_dict.keys():
            self.calcThermo(T, P)
        # Determine TST rate constant at each temperature
        Qreac = 1.0
        E0 = 0.0
        for spec in self.reactants:
            Qreac *= self.thermo_dict[T][spec.label]['Q'] / (constants.R * T / P)
            E0 -= self.thermo_dict[T][spec.label]['E0'] * 4184
        Qts = self.thermo_dict[T][self.transition_state.label]['Q'] / (constants.R * T / P)
        E0 += self.thermo_dict[T][self.transition_state.label]['E0'] * 4184
        k = (constants.kB * T / constants.h * Qts / Qreac) * math.exp(-E0 / constants.R / T)

        # Apply tunneling correction
        k *= self.transition_state.calculate_tunneling_factor(T)

        return k

    def get_free_energy_of_reaction(self, T, P=101325):
        """
        Return the Gibbs free energy of reaction in J/mol evaluated at
        temperature `T` in K.
        """
        # Calculate Gibbs free energy of reactants and products
        if T not in self.thermo_dict.keys():
            self.calcThermo(T, P)
        dGrxn = 0.0
        for reactant in self.reactants:
            dGrxn -= self.thermo_dict[T][reactant.label]['F'] * 4184
        for product in self.products:
            dGrxn += self.thermo_dict[T][product.label]['F'] * 4184
        return dGrxn

    def get_equilibrium_constant(self, T, type='Kc'):
        """
        Return the equilibrium constant for the reaction at the specified
        temperature `T` in K. The `type` parameter lets	you specify the
        quantities used in the equilibrium constant: ``Ka`` for	activities,
        ``Kc`` for concentrations (default), or ``Kp`` for pressures. Note that
        this function currently assumes an ideal gas mixture.
        """
        # Use free energy of reaction to calculate Ka
        dGrxn = self.get_free_energy_of_reaction(T)
        K = np.exp(-dGrxn / constants.R / T)
        # Convert Ka to Kc or Kp if specified
        P0 = 1e5
        if type == 'Kc':
            # Convert from Ka to Kc; C0 is the reference concentration
            C0 = P0 / constants.R / T
            K *= C0 ** (len(self.products) - len(self.reactants))
        elif type == 'Kp':
            # Convert from Ka to Kp; P0 is the reference pressure
            K *= P0 ** (len(self.products) - len(self.reactants))
        elif type != 'Ka' and type != '':
            raise ReactionError('Invalid type "{0}" passed to Reaction.get_equilibrium_constant(); '
                                'should be "Ka", "Kc", or "Kp".'.format(type))
        if K == 0:
            raise ReactionError('Got equilibrium constant of 0')
        return K

    
