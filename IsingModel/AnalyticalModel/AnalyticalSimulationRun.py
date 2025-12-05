# -*- coding: utf-8 -*-

import pickle

import os
import sys
sys.path.append('../Parameters/')

from Parameters import RunParameters

from SolveMatrixEquation import create_and_solve_steady_state
from GetExpression import get_expression_steady_state, get_expression_mean_and_instantaneous_var_steady_state, get_expression_mean_and_var_steady_state

#%% Steady state from analytcial model

def get_steady_states(htf_list, J_list, eb_list, grid_prms, file=None):
    steady_states = {}
    for J, eb in zip(J_list, eb_list):
        for htf in htf_list:
            print(f'J={J}, eb={eb}, htf={htf}')
            run_prms = RunParameters(htf, eb, J, None, grid_prms)
            
            M, probabilities = create_and_solve_steady_state(run_prms, grid_prms)
            steady_states[f'J{J}eb{eb}htf{htf}'] = (M.states.str_states, probabilities)
            
    if not file==None:
        f = open(file+'_steady_states', 'wb')
        pickle.dump(steady_states, f)
        f.close()
            
    return steady_states

#%% Mean from analytical model

def get_means(htf_list, J_list, eb_list, grid_prms, T, N=None, file=None):
    means = {}
    for J, eb in zip(J_list, eb_list):
        for htf in htf_list:
            print(f'J={J}, eb={eb}, htf={htf}')
            run_prms = RunParameters(htf, eb, J, None, grid_prms)
            
            means[f'J{J}eb{eb}htf{htf}'] = get_expression_steady_state(run_prms, grid_prms, N=N)
            
    if not file==None:
        f = open(file+'_means', 'wb')
        pickle.dump(means, f)
        f.close()
            
    return means

#%% Instantaneous standard deviation from analytical model

def get_means_and_instantaneous_vars(htf_list, J_list, eb_list, grid_prms, T, N=None, file=None):
    means_and_vars = {}
    for J, eb in zip(J_list, eb_list):
        for htf in htf_list:
            print(f'J={J}, eb={eb}, htf={htf}')
            run_prms = RunParameters(htf, eb, J, None, grid_prms)
            
            mean, var = get_expression_mean_and_instantaneous_var_steady_state(run_prms, grid_prms, N=N)
            means_and_vars[f'J{J}eb{eb}htf{htf}'] = (mean, var)
            
    if not file==None:
        f = open(file+'_means_and_instantaneous_vars', 'wb')
        pickle.dump(means_and_vars, f)
        f.close()
            
    return means_and_vars

#%% Mean and standard deviation from analytical model

def get_means_and_vars(htf_list, J_list, eb_list, grid_prms, T, N=None, limit=False, file=None):
    means_and_vars = {}
    for J, eb in zip(J_list, eb_list):
        for htf in htf_list:
            print(f'J={J}, eb={eb}, htf={htf}')
            run_prms = RunParameters(htf, eb, J, None, grid_prms)
            
            mean, var = get_expression_mean_and_var_steady_state(run_prms, grid_prms, T, N=N, limit=limit)
            means_and_vars[f'J{J}eb{eb}htf{htf}'] = (mean, var)
            
    if not file==None:
        if limit:
            f = open(file+'_means_and_vars_limit', 'wb')
        else:
            f = open(file+'_means_and_vars', 'wb')
        pickle.dump(means_and_vars, f)
        f.close()
            
    return means_and_vars

#%% Correlation time from analytical model

def get_correlation_times(htf_list, J_list, eb_list, grid_prms, T, N=None, file=None):
    correlation_times = {}
    for J, eb in zip(J_list, eb_list):
        for htf in htf_list:
            print(f'J={J}, eb={eb}, htf={htf}')
            run_prms = RunParameters(htf, eb, J, None, grid_prms)
            
            _, _, _, tau = get_expression_mean_and_var_steady_state(run_prms, grid_prms, T, N=N, limit=True, return_all=True)
            correlation_times[f'J{J}eb{eb}htf{htf}'] = tau
            
    if not file==None:
        f = open(file+'_correlation_times', 'wb')
        pickle.dump(correlation_times, f)
        f.close()
            
    return correlation_times

#%%
if __name__ == '__main__':
    
    from Metadata import read_metadata
    
    if len(sys.argv)<2:
        print("Please supply a metadatafile as commandline argument")
        sys.exit()

    file_name    = sys.argv[1]
    file         = os.path.join(os.path.realpath('__file__'), f'../../Data/{file_name}') 

    grid_prms, sim_prms, varied_prms = read_metadata(file)
    
    f = open(file+'_shift_compensation', 'rb')
    shift_compensation = pickle.load(f)
    f.close()

    get_means_and_vars(varied_prms.htf_list, varied_prms.J_list, varied_prms.eb_list, grid_prms, sim_prms.T, N=shift_compensation.N_expr, file=file)
