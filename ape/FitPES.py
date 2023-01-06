# -*- coding: utf-8 -*-
import csv
import numpy as np
import ast
import rmgpy.constants as constants

def from_sampling_result(csv_path):
    """
    Rebuild mode_dict and energy_dict by importing sampling_result.csv.
    """
    mode_dict = {}
    energy_dict = {}
    rotors = {}
    with open(csv_path, 'r') as f:
        reader = csv.reader(f, dialect='excel')
        for row in reader:
            if row[0].split('_')[0] == 'mode':
                nmode, mode=row[0].split('_')[1:]
                nmode = int(nmode)
                mode_dict[nmode] = {}
                mode_dict[nmode]['mode'] = mode
                energy_dict[nmode] = {}
            elif row[0] == 'job_keys': job_keys = ast.literal_eval(row[1])
            elif row[0] == 'min_elect': min_elect = float(row[1])
            elif row[0] == 'symmetry_number': mode_dict[nmode]['symmetry_number'] = int(row[1])
            elif row[0] == 'rotor': mode_dict[nmode]['rotor'] = ast.literal_eval(row[1]); rotors[nmode] = ast.literal_eval(row[1])
            elif row[0] == 'M': mode_dict[nmode]['M'] = float(row[1])
            elif row[0] == 'K': mode_dict[nmode]['K'] = float(row[1])
            elif row[0] == 'step_size': mode_dict[nmode]['step_size'] = float(row[1])
            elif row[0] == 'sample': continue
            else:
                ind, e_elect = int(row[0]), float(row[1])
                energy_dict[nmode][ind] = e_elect
    return mode_dict, energy_dict, min_elect, rotors, job_keys

def get_is_tors_list(mode_dict):
    """
    A list contains the boolean values of each mode to determine 
    if the normal mode is internal rotation or not.
    """
    is_tors_list = []
    for mode in sorted(mode_dict.keys()):
        is_tors = True if mode_dict[mode]['mode'] == 'tors' else False
        is_tors_list.append(is_tors)
    return is_tors_list

def get_k_list(mode_dict):
    """
    A list contains the second derivatives of each mode at the origin.
    """
    k_list = []
    for mode in sorted(mode_dict.keys()):
        k_list.append(mode_dict[mode]['K'])
    return k_list

def get_delta_q_list(mode_dict):
    delta_q_list = []
    for mode in sorted(mode_dict.keys()):
        M = mode_dict[mode]['M']
        step_size = mode_dict[mode]['step_size']
        delta_q = np.sqrt(M) * step_size
        delta_q_list.append(delta_q)
    return delta_q_list

unitConverter = constants.E_h / constants.amu * 10 ** 20 # from hartree/angstrom^2/amu to 1/s^2
def cubic_spline_interpolations(energy_dict,mode_dict):
    k_list = get_k_list(mode_dict)
    delta_q_list = get_delta_q_list(mode_dict)
    is_tors_list = get_is_tors_list(mode_dict)
    def delta_E(mode,i):
        delta_E = energy_dict[mode][i]-energy_dict[mode][i-1]
        return delta_E
    polynomial_dict = {}
    nmode = len(energy_dict)
    for mode in range(nmode):
        delta_q = delta_q_list[mode]
        k = k_list[mode]
        is_tors = is_tors_list[mode]
        mode += 1
        R = []
        samples = sorted(energy_dict[mode].keys())
        n_samples = len(energy_dict[mode])
        for i in range(n_samples):
            if i == 0:
                R.append([3*delta_E(mode,samples[i+1])/(delta_q**2)-0.5*k/unitConverter])   
            elif i != 0 and i != n_samples-1:
                R.append([3*(delta_E(mode,samples[i])/(delta_q**2)+delta_E(mode,samples[i+1])/(delta_q**2))])
            else:
                R.append([3*delta_E(mode,samples[i])/(delta_q**2)+0.5*k/unitConverter])
        R = np.matrix(R)
        L = []
        for i in range(n_samples):
            L_list = [0] * len(energy_dict[mode])
            if i == 0:
                L_list[i] = 2/delta_q
                L_list[i+1] = 1/delta_q
                L.append(L_list)
            elif i != 0 and i != len(energy_dict[mode])-1:
                L_list[i-1] = 1/delta_q
                L_list[i] = 2*(1/delta_q+1/delta_q)
                L_list[i+1] = 1/delta_q
                L.append(L_list)
            else:
                L_list[i-1] = 1/delta_q
                L_list[i] = 2/delta_q
                L.append(L_list)
        L = np.matrix(L)
        D = np.linalg.inv(L).dot(R).T.tolist()[0]

        each_polynomial_dict = {}
        for i in range(n_samples):
            n_step = samples[i]
            each_spline_polynomial_dict = {}
            if i == 0:
                pass
            else:
                a_i = energy_dict[mode][samples[i-1]]
                b_i = D[i-1]*delta_q
                c_i = 3*delta_E(mode,samples[i])-(2*D[i-1]+D[i])*delta_q
                d_i = -2*delta_E(mode,samples[i])+(D[i-1]+D[i])*delta_q

                ai = a_i-b_i*(n_step-1)+c_i*(n_step-1)**2-d_i*(n_step-1)**3
                bi = b_i*(1/delta_q)-2*c_i*(1/delta_q)*(n_step-1)+3*d_i*(1/delta_q)*(n_step-1)**2
                ci = c_i*((1/delta_q)**2)-3*d_i*(1/delta_q)**2*(n_step-1)
                di = d_i*(1/delta_q)**3

                each_spline_polynomial_dict['ai'] = ai
                each_spline_polynomial_dict['bi'] = bi
                each_spline_polynomial_dict['ci'] = ci
                each_spline_polynomial_dict['di'] = di
                each_polynomial_dict[n_step] = each_spline_polynomial_dict
        polynomial_dict[mode] = each_polynomial_dict
        if is_tors == False:
            step = delta_q
            #-inf to x_0
            ind = samples[0]
            delta_q = samples[0]*step
            ai = each_polynomial_dict[ind+1]['ai']
            bi = each_polynomial_dict[ind+1]['bi']
            ci = each_polynomial_dict[ind+1]['ci']
            di = each_polynomial_dict[ind+1]['di']
            ddy = 2*ci+6*di*delta_q
            dy = bi+2*ci*delta_q+3*di*delta_q**2
            each_spline_polynomial_dict = {}
            each_spline_polynomial_dict['ai'] = 0.5*ddy*delta_q**2-dy*delta_q+energy_dict[mode][samples[0]]
            each_spline_polynomial_dict['bi'] = -ddy*delta_q+dy
            each_spline_polynomial_dict['ci'] = 0.5*ddy
            each_spline_polynomial_dict['di'] = 0
            polynomial_dict[mode][ind] = each_spline_polynomial_dict
            
            #x_ndata to inf
            ind = samples[-1] + 1
            delta_q = samples[-1]*step
            ai = each_polynomial_dict[ind-1]['ai']
            bi = each_polynomial_dict[ind-1]['bi']
            ci = each_polynomial_dict[ind-1]['ci']
            di = each_polynomial_dict[ind-1]['di']
            ddy = 2*ci+6*di*delta_q
            dy = bi+2*ci*delta_q+3*di*delta_q**2
            each_spline_polynomial_dict = {}
            each_spline_polynomial_dict['ai'] = 0.5*ddy*delta_q**2-dy*delta_q+energy_dict[mode][samples[-1]]
            each_spline_polynomial_dict['bi'] = -ddy*delta_q+dy
            each_spline_polynomial_dict['ci'] = 0.5*ddy
            each_spline_polynomial_dict['di'] = 0
            polynomial_dict[mode][ind] = each_spline_polynomial_dict
    #print('Cubic spline interpolations completed.')
    return polynomial_dict

def plot(polynomial_dict, energy_dict, mode_dict, mode):
    """
    Can offer a MATLAB code to plot the 1-D potential energy surface to check if the cubic spline interpolations works well.
    The automatically generating ``.png``, ``.svg``, ``.pdf``, and ``.ps`` will be developed in the future.
    The purpose of this function is to debug and create a 1-D potential energy surface plot.
    """
    delta_q_list = get_delta_q_list(mode_dict)
    delta_q = delta_q_list[mode-1]
    samples = sorted(energy_dict[mode].keys())
    for i in sorted(polynomial_dict[mode].keys()):
        if i == samples[0]:
            x1 = -np.inf
            x2 = delta_q*i
        elif i == samples[-1] + 1:
            x1 = delta_q*(i-1)
            x2 = np.inf
        else:
            x2 = delta_q*i
            x1 = delta_q*(i-1)
        a = [polynomial_dict[mode][i][ind] for ind in ['ai','bi','ci','di']]
        print('ezplot(\'{}+{}*Q+{}*Q^2+{}*Q^3\', [{} {}])'.format(a[0],a[1],a[2],a[3],x1,x2))         

################################################################################

if __name__ == '__main__':
    csv_path = '../examples/propane_UMVT/propane_samping_result.csv'
    mode_dict, energy_dict, _ = from_sampling_result(csv_path)
    polynomial_dict = cubic_spline_interpolations(energy_dict,mode_dict)
    #A MATLAB code to plot the 1-D PES
    plot = plot(polynomial_dict,energy_dict,mode_dict,mode=2)