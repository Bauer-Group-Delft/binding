# -*- coding: utf-8 -*-

import numpy as np
import scipy as sp

from FullMatrix import FullRatesAndEnergies, FullStates, FullMatrix

#%% Solving full matrix

def create(run_prms, grid_prms):
    R = FullRatesAndEnergies(run_prms.htf, run_prms.eb, run_prms.J, beta=grid_prms.beta, adapt_conc=grid_prms.do_htf_update)
    S = FullStates(grid_prms.L)
    M = FullMatrix(S, R, grid_prms.sites.flatten(), grid_prms.periodic)
    return M

def create_and_solve_steady_state(run_prms, grid_prms, compact=False):
    M = create(run_prms, grid_prms)
    P_steady = M.calculate_steady_state()
    if compact:
        P_steady = M.compact_states(P_steady)
    return M, P_steady

def solve_master_equation(mat, T, v0, steps=1):
    ts = np.linspace(0,T,steps+1)
    
    # Solution
    v = []
    for t in ts:
        v.append(sp.sparse.linalg.expm(mat*t)@v0)
    
    return ts, np.array(v)

def create_and_solve(run_prms, grid_prms, T, steps=1, compact=False):
    M = create(run_prms, grid_prms)
    mat = M.create_matrix(sparse=True)
    v0  = M.states.create_initial_condition(grid_prms.IC_mode, run_prms.IC_prms)
    
    t, v = solve_master_equation(mat, T, v0, steps)
        
    if compact:
        v_compact = M.compact_states(v)
        return M, t, v, v_compact
        
    return M, t, v