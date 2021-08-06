# -*- coding: utf-8 -*-
import logging
import numpy as np
from math import sqrt, sin, cos, acos
from math import factorial as fact

import rmgpy.constants as constants

from ape.HarmonicBasis import IntXHmHnexp
from ape.FourierBasis import IntXPhimPhin

# fix overflow problem
from decimal import Decimal as D
from decimal import getcontext
getcontext().prec = 15

import multiprocessing as mp

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

            negative_energy = check_negative_energy(a,x1,x2)
            if negative_energy:
                root1, root2 = negative_energy
                result += IntXPhimPhin(m,n,x1,root1,L,a)
                result += IntXPhimPhin(m,n,root2,x2,L,a)
            else:
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
                x2 = i
            elif i == samples[-1] + 1:
                x1 = (i-1)
                x2 = np.inf
            else:
                x2 = i
                x1 = (i-1)
            
            f1 = D.sqrt(D(fact(m)))
            f2 = D(pow(2,m/2.0)*pow(2,n/2.0)*sqrt(np.pi))
            f3 = D.sqrt(D(fact(n)))
            f = [f1,f2,f3]

            negative_energy = check_negative_energy(a,x1,x2)
            if negative_energy:
                root1, root2 = negative_energy
                result += IntXHmHnexp(m,n,x1,root1,a,f)
                result += IntXHmHnexp(m,n,root2,x2,a,f)
            else:
                result += IntXHmHnexp(m,n,x1,x2,a,f)

        if m==n: result += D(-(1/2)*P*(2*m+1))
        elif m == (n+2): result += D(sqrt(m)*sqrt(m-1)*(1/2)*P)

    return result

def SetAnharmonicH(polynomial_dict, mode_dict, energy_dict, mode, size, N_prev, H_prev):
    H_ind = []
    is_tors = True if mode_dict[mode]['mode'] == 'tors' else False
    pool = mp.Pool(mp.cpu_count())
    H = np.zeros((size, size), np.float64)
    if N_prev == 0:
        for m in range(size):
            for n in range(m+1):
                H_ind.append((m,n))
        H_temp = pool.starmap(Hmn, [(m, n, polynomial_dict, mode_dict, energy_dict, mode, is_tors) for m,n in H_ind])
        i = 0
        for m in range(size):
            for n in range(m+1):
                H[m][n] = H_temp[i]
                H[n][m] = H_temp[i]
                i += 1
        pool.close()
    elif size > N_prev:
        for m in range(N_prev):
            for n in range(m+1):
                H[m][n] = H_prev[m][n]
                H[n][m] = H_prev[n][m]
        m = size - 1
        for n in range(size):
            H_ind.append((m,n))
        H_temp = pool.starmap(Hmn, [(m, n, polynomial_dict, mode_dict, energy_dict, mode, is_tors) for m,n in H_ind])
        for n in range(size):
            H[m][n] = H_temp[n]
            H[n][m] = H_temp[n]
    pool.close()
    return H

def check_negative_energy(coeff, x1, x2):
    a = coeff[3]
    b = coeff[2]
    c = coeff[1]
    d = coeff[0]
    if a == 0:
        return False
    f = (3 * c / a - pow(b / a, 2)) / 3
    g = (2 * pow(b / a, 3) - (9 * b * c / pow(a, 2)) + (27 * d / a)) / 27
    h = pow(g, 2) / 4 + pow(f, 3) / 27

    if h < 0:
        i = sqrt(pow(g, 2) / 4 - h)
        j = pow(i, 1.0 / 3)
        k = acos(-g / (2 * i))
        l = -j
        m = cos(k / 3)
        n = sqrt(3) * sin(k / 3)
        p = -b / (3 * a)
        r1 = 2 * j * cos(k / 3) - b / (3 * a)
        r2 = l * (m + n) + p
        r3 = l * (m - n) + p
        if (r1 > x1 and r1 < x2) and (r2 > x1 and r2 < x2):
            if r1 > r2:
                root2 = r1
                root1 = r2
            else:
                root2 = r2
                root1 = r1
            return [root1, root2]
        elif (r1 > x1 and r1 < x2) and (r3 > x1 and r3 < x2):
            if r1>r3:
                root2 = r1
                root1 = r3
            else:
                root2 = r3
                root1 = r1
            return [root1, root2]
        elif (r2 > x1 and r2 < x2) and (r3 > x1 and r3 < x2):
            if r2 > r3:
                root2 = r2
                root1 = r3
            else:
                root2 = r3
                root1 = r2
            return [root1, root2]
    return False
