#!/usr/bin/env python3

"""
APE kinetics module
"""

import logging
import os.path
import string

import numpy as np

import rmgpy.quantity as quantity
from rmgpy.kinetics.arrhenius import Arrhenius
from rmgpy.kinetics.tunneling import Wigner, Eckart

from arkane.output import prettify
from arkane.kinetics import KineticsJob as RMG_KineticsJob
from arkane.kinetics import KineticsDrawer

class KineticsJob(object):
    """
    A representation of an APE kinetics job. This job is used to compute
    and save the high-pressure-limit kinetics information for a single reaction.

    `usedTST` - a boolean representing if TST was used to calculate the kinetics
                if kinetics is already given in the input, then it is False.
    `three_params` - a boolean representing if the modified three-parameter Arrhenius equation is used to calculate
                     high pressure kinetic rate coefficients. If it is False, the classical two-parameter Arrhenius
                     equation is used.
    """

    def __init__(self, reaction, Tmin=None, Tmax=None, Tlist=None, Tcount=0, three_params=True, ncpus=None):
        self.usedTST = False
        self.Tmin = Tmin if Tmin is not None else (298, 'K')
        self.Tmax = Tmax if Tmax is not None else (2500, 'K')
        self.Tcount = Tcount if Tcount > 3 else 50
        self.three_params = three_params
        self.ncpus = ncpus

        if Tlist is not None:
            self.Tlist = Tlist
            self.Tmin = (min(self.Tlist.value_si), 'K')
            self.Tmax = (max(self.Tlist.value_si), 'K')
            self.Tcount = len(self.Tlist.value_si)
        else:
            self.Tlist = (1 / np.linspace(1 / self.Tmax.value_si, 1 / self.Tmin.value_si, self.Tcount), 'K')

        self.reaction = reaction
        self.k_units = None

    @property
    def Tmin(self):
        """The minimum temperature at which the computed k(T) values are valid, or ``None`` if not defined."""
        return self._Tmin

    @Tmin.setter
    def Tmin(self, value):
        self._Tmin = quantity.Temperature(value)

    @property
    def Tmax(self):
        """The maximum temperature at which the computed k(T) values are valid, or ``None`` if not defined."""
        return self._Tmax

    @Tmax.setter
    def Tmax(self, value):
        self._Tmax = quantity.Temperature(value)

    @property
    def Tlist(self):
        """The temperatures at which the k(T) values are computed."""
        return self._Tlist

    @Tlist.setter
    def Tlist(self, value):
        self._Tlist = quantity.Temperature(value)
    
    def execute(self, output_directory=None, plot=True, print_HOhf_result=True):
        """
        Execute the kinetics job, saving the results within
        the `output_directory`.

        If `plot` is True, then plots of the raw and fitted values for the kinetics
        will be saved.
        """
        # Assign the ncpus for each Reaction object
        for rxn in self.reaction:
            rxn.ncpus = self.ncpus
        
        self.rmg_reaction = self.reaction.rmg_Reaction()
        self.generate_kinetics()
        if output_directory is not None:
            try:
                self.write_output(output_directory)
            except Exception as e:
                logging.warning("Could not write kinetics output file due to error: "
                                "{0} in reaction {1}".format(e, self.reaction.label))
        
            if plot:
                try:
                    self.plot(output_directory)
                except Exception as e:
                    logging.warning("Could not plot kinetics due to error: "
                                    "{0} in reaction {1}".format(e, self.reaction.label))
                try:
                    self.draw(output_directory)
                except Exception as e:
                    logging.warning("Could not draw reaction {1} due to error: {0}".format(e, self.reaction.label))
        logging.debug('Finished kinetics job for reaction {0}.'.format(self.reaction.label))

        if print_HOhf_result:
            f = open(os.path.join(output_directory, 'output.py'), 'a')
            f.write('\n#   By using harmonic approximation method\n')
            f.close()
            rmg_KineticsJob = RMG_KineticsJob(self.rmg_reaction, Tlist=self.Tlist, three_params=self.three_params)
            rmg_KineticsJob.generate_kinetics()
            rmg_KineticsJob.execute(output_directory, plot=False)

    def generate_kinetics(self):
        """
        Generate the kinetics data for the reaction and fit it to a modified Arrhenius model.
        """
        self.usedTST = True
        kinetics_class = 'Arrhenius'

        tunneling = self.reaction.transition_state.tunneling
        if isinstance(tunneling, Wigner) and tunneling.frequency is None:
            tunneling.frequency = (self.reaction.transition_state.frequency.value_si, "cm^-1")
        elif isinstance(tunneling, Eckart) and tunneling.frequency is None:
            tunneling.frequency = (self.reaction.transition_state.frequency.value_si, "cm^-1")
            tunneling.E0_reac = (sum([reactant.conformer.E0.value_si
                                      for reactant in self.reaction.reactants]) * 0.001, "kJ/mol")
            tunneling.E0_TS = (self.reaction.transition_state.conformer.E0.value_si * 0.001, "kJ/mol")
            tunneling.E0_prod = (sum([product.conformer.E0.value_si
                                      for product in self.reaction.products]) * 0.001, "kJ/mol")
        elif tunneling is not None:
            if tunneling.frequency is not None:
                # Frequency was given by the user
                pass
            else:
                raise ValueError('Unknown tunneling model {0!r} for reaction {1}.'.format(tunneling, self.reaction.label))
        logging.debug('Generating {0} kinetics model for {1}...'.format(kinetics_class, self.reaction.label))
        order = len(self.reaction.reactants)
        self.k_units = {1: 's^-1', 2: 'cm^3/(mol*s)', 3: 'cm^6/(mol^2*s)'}[order]
        self.K_eq_units = {2: 'mol^2/cm^6', 1: 'mol/cm^3', 0: '       ', -1: 'cm^3/mol', -2: 'cm^6/mol^2'}[
            len(self.reaction.products) - len(self.reaction.reactants)]
        self.k_r_units = {1: 's^-1', 2: 'cm^3/(mol*s)', 3: 'cm^6/(mol^2*s)'}[len(self.reaction.products)]

        # Initialize Object for Converting Units
        if self.K_eq_units != '       ':
            keq_unit_converter = quantity.Units(self.K_eq_units).get_conversion_factor_from_si()
        else:
            keq_unit_converter = 1
        
        factor = 1e6 ** (order - 1)

        self.k_list = np.zeros_like(self.Tlist.value_si)
        self.k0_list = np.zeros_like(self.Tlist.value_si)
        self.kappa_list = np.zeros_like(self.Tlist.value_si)
        self.Keq_list = np.zeros_like(self.Tlist.value_si)
        for i, T in enumerate(self.Tlist.value_si):
            logging.info('Generate kinetics for {0} at {1} K...'.format(self.reaction.label, T))
            tunneling = self.reaction.transition_state.tunneling
            self.reaction.transition_state.tunneling = None
            k0 = self.reaction.calculate_tst_rate_coefficient(T) * factor
            self.k0_list[i] = k0
            self.reaction.transition_state.tunneling = tunneling
            k = self.reaction.calculate_tst_rate_coefficient(T) * factor
            self.k_list[i] = k
            self.kappa_list[i] = k / k0
            tunneling = self.reaction.transition_state.tunneling
            self.Keq_list[i] = keq_unit_converter * self.reaction.get_equilibrium_constant(T) # returns SI units

        self.reaction.kinetics = Arrhenius().fit_to_data(self.Tlist.value_si, self.k_list, kunits=self.k_units,
                                                         three_params=self.three_params)

    def write_output(self, output_directory):
        """
        Save the results of the kinetics job to the `output.py` file located
        in `output_directory`.
        """
        reaction = self.reaction

        k0_revs, k_revs = [], []
        
        logging.info('Saving kinetics for {0}...'.format(reaction.label))

        f = open(os.path.join(output_directory, 'output.py'), 'a')

        if self.usedTST:
            # If TST is not used, eg. it was given in 'reaction', then this will throw an error.
            f.write('#   By using APE algorithm\n')
            f.write('#   ======= =========== =========== =========== ===============\n')
            f.write('#   Temp.   k (TST)     Tunneling   k (TST+T)   Units\n')
            f.write('#   ======= =========== =========== =========== ===============\n')

            t_list = self.Tlist.value_si
            
            for i, T in enumerate(t_list):
                k0 = self.k0_list[i]
                kappa = self.kappa_list[i]
                k = self.k_list[i]

                f.write('#    {0:4g} K {1:11.3e} {2:11g} {3:11.3e} {4}\n'.format(T, k0, kappa, k, self.k_units))
            f.write('#   ======= =========== =========== =========== ===============\n')
            f.write('\n\n')

            f.write('#   ======= ============ =========== ============ ============= =========\n')
            f.write('#   Temp.    Kc (eq)        Units     k_rev (TST) k_rev (TST+T)   Units\n')
            f.write('#   ======= ============ =========== ============ ============= =========\n')

            for i, T in enumerate(t_list):
                k = self.k_list[i]
                k0 = self.k0_list[i]
                K_eq = self.Keq_list[i]
                k0_rev = k0 / K_eq
                k_rev = k / K_eq
                k0_revs.append(k0_rev)
                k_revs.append(k_rev)
                f.write('#    {0:4g} K {1:11.3e}   {2}  {3:11.3e}   {4:11.3e}      {5}\n'.format(
                    T, K_eq, self.K_eq_units, k0_rev, k_rev, self.k_r_units))

            f.write('#   ======= ============ =========== ============ ============= =========\n')
            f.write('\n\n')

            kinetics_0_rev = Arrhenius().fit_to_data(t_list, np.array(k0_revs), kunits=self.k_r_units,
                                                     three_params=self.three_params)
            kinetics_rev = Arrhenius().fit_to_data(t_list, np.array(k_revs), kunits=self.k_r_units,
                                                   three_params=self.three_params)

            f.write('# k_rev (TST) = {0} \n'.format(kinetics_0_rev))
            f.write('# k_rev (TST+T) = {0} \n\n'.format(kinetics_rev))

        if self.three_params:
            f.write('# kinetics fitted using the modified three-parameter Arrhenius equation '
                    'k = A * (T/T0)^n * exp(-Ea/RT) \n')
        else:
            f.write('# kinetics fitted using the two-parameter Arrhenius equation k = A * exp(-Ea/RT) \n')

        # Reaction path degeneracy is INCLUDED in the kinetics itself!
        rxn_str = 'kinetics(label={0!r}, kinetics={1!r})'.format(reaction.label, reaction.kinetics)
        f.write('{0}\n\n'.format(prettify(rxn_str)))

        f.close()
    
    def plot(self, output_directory):
        """
        Plot both the raw kinetics data and the Arrhenius fit versus 
        temperature. The plot is saved to the file ``kinetics.pdf`` in the
        output directory. The plot is not generated if ``matplotlib`` is not
        installed.
        """
        import matplotlib.pyplot as plt

        f, ax = plt.subplots()
        if self.Tlist is not None:
            t_list = [t for t in self.Tlist.value_si]
        else:
            t_list = 1000.0 / np.arange(0.4, 3.35, 0.05)
        klist = np.zeros_like(t_list)
        klist2 = np.zeros_like(t_list)
        for i in range(len(t_list)):
            klist[i] = self.reaction.calculate_tst_rate_coefficient(t_list[i])
            klist2[i] = self.reaction.kinetics.get_rate_coefficient(t_list[i])

        order = len(self.reaction.reactants)
        klist *= 1e6 ** (order - 1)
        klist2 *= 1e6 ** (order - 1)
        t_list = [1000.0 / t for t in t_list]
        plt.semilogy(t_list, klist, 'ob', label='TST calculation')
        plt.semilogy(t_list, klist2, '-k', label='Fitted rate')
        plt.legend()
        reaction_str = '{0} {1} {2}'.format(
            ' + '.join([reactant.label for reactant in self.reaction.reactants]),
            '<=>', ' + '.join([product.label for product in self.reaction.products]))
        plt.title(reaction_str)
        plt.xlabel('1000 / Temperature (K^-1)')
        plt.ylabel('Rate coefficient ({0})'.format(self.k_units))

        plot_path = os.path.join(output_directory, 'plots')

        if not os.path.exists(plot_path):
            os.mkdir(plot_path)
        valid_chars = "-_.()<=> %s%s" % (string.ascii_letters, string.digits)
        filename = ''.join(c for c in reaction_str if c in valid_chars) + '.pdf'
        plt.savefig(os.path.join(plot_path, filename))
        plt.close()

    def draw(self, output_directory, file_format='pdf'):
        """
        Generate a PDF drawing of the reaction.
        This requires that Cairo and its Python wrapper be available; if not,
        the drawing is not generated.
        You may also generate different formats of drawings, by changing format to
        one of the following: `pdf`, `svg`, `png`.
        """

        drawing_path = os.path.join(output_directory, 'paths')

        if not os.path.exists(drawing_path):
            os.mkdir(drawing_path)
        valid_chars = "-_.()<=> %s%s" % (string.ascii_letters, string.digits)
        reaction_str = '{0} {1} {2}'.format(
            ' + '.join([reactant.label for reactant in self.reaction.reactants]),
            '<=>', ' + '.join([product.label for product in self.reaction.products]))
        filename = ''.join(c for c in reaction_str if c in valid_chars) + '.pdf'
        path = os.path.join(drawing_path, filename)

        KineticsDrawer().draw(self.reaction, file_format=file_format, path=path)