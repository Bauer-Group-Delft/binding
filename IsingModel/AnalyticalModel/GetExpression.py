# -*- coding: utf-8 -*-

import sys
sys.path.append('../Parameters/')

from Parameters import RunParameters

from SolveMatrixEquation import create, create_and_solve_steady_state
from MatrixExponentialIntegrals import matrix_exponential_integral, mean_variance_and_correlation_time_occupancy

#%% Expression estimation methods

'''
Expression can be calculated in different ways:
    - mean occupancy: This assumes that the rate of transcriptions is linearly 
      related to the number of TFs bound. Each TF contributes independently to 
      the transcription.
    - fraction of time with at least N TFs bound: This assumes that there is only 
      transcription when at least N TFs are bound. The rate of transcription 
      is constant, i.e. independent of how many more than N TFs are bound.
The average expression can be calculated in steady state or over the runtime.

'''

def get_expression_per_state(M, N):
    if N==None:
        # Use the occupancy
        occupancies = M.calculate_state_occupancies()
        expression  = occupancies
    else:
        # Use the probability of having at least N TFs bound
        occupancies = M.calculate_state_occupancies(normalize=False)
        expression  = (occupancies>=N)
    return expression

def get_expression_steady_state(run_prms, grid_prms, N=None):
    # In steady state the probabilities of the states are constant, so the time mean 
    # over the runtime is equal to the probability at any time.

    M, v = create_and_solve_steady_state(run_prms, grid_prms)
    expression = get_expression_per_state(M, N)

    return expression@v

def get_expression_mean_and_var_steady_state(run_prms, grid_prms, T, N=None, limit=False, return_all=False, diagonalize=False):
    M, v = create_and_solve_steady_state(run_prms, grid_prms)
    mat = M.create_matrix(sparse=False)
    expression = get_expression_per_state(M, N)
    
    return mean_variance_and_correlation_time_occupancy(mat, T, v, expression, limit=limit, return_all=return_all, diagonalize=diagonalize)

def get_expression_runtime(run_prms, grid_prms, T, N=None):
    # If the simulation does not start in steady state, the probabilities of the states
    # will change over time. To get the time mean over the runtime, we need to integrate
    # over the matrix exponential.
    
    M = create(run_prms, grid_prms)
    mat = M.create_matrix(sparse=False)
    integral = matrix_exponential_integral(mat, T)
    expression = get_expression_per_state(M, N)
    
    v0 = M.states.create_initial_condition(grid_prms.IC_mode, run_prms.IC_prms)
    
    return (expression@integral@v0)/T

def get_expression(eb, args):
    run_prms = RunParameters(args['htf'], eb, args['J'], args['IC_prms'], args['grid_prms'])
    
    if args['T']==None:
        return get_expression_steady_state(run_prms, args['grid_prms'], N=args['N'])
    else:
        return get_expression_runtime(run_prms, args['grid_prms'], args['T'], N=args['N'])
    