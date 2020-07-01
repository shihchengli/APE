# -*- coding: utf-8 -*-
import numpy as np
from math import sqrt
from math import factorial as fact

import rmgpy.constants as constants

from ape.HarmonicBasis import IntXHmHnexp
from ape.FourierBasis import IntXPhimPhin

hbar1 = constants.hbar / constants.E_h # in hartree*s
hbar2 = constants.hbar * 10 ** 20 / constants.amu # in amu*angstrom^2/s
def Hmn(m, n, polynomial_dict, mode_dict, energy_dict, mode, is_tors):
    result = 0
    if is_tors:
        # use fourier basis function
        I = mode_dict[mode]['M'] # in amu*angstrom^2
        k = mode_dict[mode]['K'] # 1/s^2
        step_size = mode_dict[mode]['step_size']
        delta_q = sqrt(I) * step_size # in sqrt(amu)*angstrom
        L = np.pi * sqrt(I) # in sqrt(amu)*angstrom
        x2 = 0
        for i in sorted(polynomial_dict[mode].keys()):
            x1 = x2
            x2 += delta_q
            a = [polynomial_dict[mode][i][ind] for ind in ['ai','bi','ci','di']]
            result += IntXPhimPhin(m,n,x1,x2,L,a)
        if (m==n and m!=0):
            if (m%2==0): m /= 2
            else: m = (m+1)/2
            result += pow(m*np.pi/L,2)*(hbar1*hbar2)/2 # in hartree
        
    else:
        # use harmonic basis functions
        M = mode_dict[mode]['M'] # in amu
        k = mode_dict[mode]['K'] # 1/s^2
        step_size = mode_dict[mode]['step_size'] # in angstrom
        delta_q = sqrt(M) * step_size # in sqrt(amu)*angstrom
        wj = sqrt(k) # in 1/s
        P = -hbar1*wj/2 # in hartree
        R = sqrt(hbar2/wj) # in sqrt(amu)*angstrom
        samples = sorted(energy_dict[mode].keys())
        for i in sorted(polynomial_dict[mode].keys()):
            a = [polynomial_dict[mode][i][ind] for ind in ['ai','bi','ci','di']]
            a[0] = a[0]
            a[1] = a[1]*pow(R,1)
            a[2] = a[2]*pow(R,2)
            a[3] = a[3]*pow(R,3)

            if i == samples[0]:
                x1 = -np.inf
                x2 = delta_q*i/R
            elif i == samples[-1] + 1:
                x1 = delta_q*(i-1)/R
                x2 = np.inf
            else:
                x2 = delta_q*i/R
                x1 = delta_q*(i-1)/R

            result += IntXHmHnexp(m,n,x1,x2,a)
        result /= pow(2,m/2.0)*pow(2,n/2.0)*sqrt(fact(m))*sqrt(fact(n))*sqrt(np.pi)
        if m==n: result += -(1/2)*P*(2*m+1)
        elif m == (n+2): result += sqrt(m)*sqrt(m-1)*(1/2)*P
    return result

def SetAnharmonicH(polynomial_dict, mode_dict, energy_dict, mode, size, N_prev, H_prev):
    H = np.zeros((size, size), np.float64)
    is_tors = True if mode_dict[mode]['mode'] == 'tors' else False
    for m in range(size):
        for n in range(m+1):
            if m < N_prev and n < N_prev:
                Hmn_val = H_prev[m][n]
            else:
                Hmn_val = Hmn(m, n, polynomial_dict, mode_dict, energy_dict, mode, is_tors)
            H[m][n] = Hmn_val
            H[n][m] = Hmn_val
    return H