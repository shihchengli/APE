import numpy as np
import rmgpy.constants as constants
from rmgpy.statmech.torsion import FreeRotor, HinderedRotor

#Compare to the result of tables in https://aip.scitation.org/doi/pdf/10.1063/1.1723744

symmetry = 3
T = 298.15 
inertia = (1/0.8/(np.sqrt(8 * np.pi ** 3 * constants.kB * T)/(symmetry * constants.h)))**2
inertia /= constants.amu*10**(-20)

mode = FreeRotor(
    inertia=(inertia, "amu*angstrom^2"),
    symmetry=symmetry,
    )
Q = mode.get_partition_function(T)
print('FreeRotor')
print('1/Q: %.2f' % (1/Q))
print('V/RT: %.1f' % (0))
entropy = mode.get_entropy(T)
print('S: %.3f' % (entropy/4.184))
print('')
############################################################################

barrier = 16 * (constants.R*T) / 1000
mode = HinderedRotor(
    inertia=(inertia, "amu*angstrom^2"),
    symmetry=symmetry,
    barrier=(barrier, "kJ/mol"),
    quantum=True,
    )
print('HinderedRotor')
phi = np.arange(0.0, 2 * constants.pi + 0.0001, constants.pi / 18.)
potential = np.zeros_like(phi)
for i in range(phi.shape[0]):
    potential[i] = mode.get_potential(phi[i]) / (constants.Na * constants.E_h)
    #print(potential[i])
freq = mode.get_frequency() # in cm-1
#print(freq)
print('1/Q: %.2f' % (1/Q))
print('V/RT: %.1f' % (barrier/((constants.R*T) / 1000)))
k = (freq * (2 * np.pi * constants.c * 100)) ** 2 # in 1/s^2
#print('M: ',inertia)
#print('K: ',k)
entropy = mode.get_entropy(T)
print('S: %.3f' % (entropy/4.184))