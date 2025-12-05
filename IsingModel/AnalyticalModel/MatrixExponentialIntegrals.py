# -*- coding: utf-8 -*-

import numpy as np
import scipy as sp
import scipy.sparse
import scipy.sparse.linalg

def matrix_exponential_integral(mat, T, mode=1):
    '''
    This function calculates the integral:
        
        int_0^T e^{M\tau} d\tau.
        
    This can be used to calculate the mean occupancy over a time T:
        
        (expression@integral@v0)/T
    '''
    
    if mode == 0: # Invertible matrix
        I = np.identity(len(mat))
        integral = np.linalg.inv(mat)@(sp.linalg.expm(mat*T)-I)
    
    elif mode == 1: # Matrix with one eigenvalue equal to 0
        eigenvalues, eigenvectors = np.linalg.eig(mat)
        
        order = np.argsort(np.abs(eigenvalues))
        D = np.diag(eigenvalues[order])
        P = eigenvectors[:,order]
        
        ## Test: accuracy of diagonalization, PDP^-1 should be equal to the original matrix.
        # print(mat)
        # print(eigenvectors@np.diag(eigenvalues)@np.linalg.inv(eigenvectors))
        # print(P@D@np.linalg.inv(P))
        ## --------------------------------------------------------------------------------
        
        D_int = np.zeros_like(D)
        D_int[0,0] = T
        D_int[1:,1:] = matrix_exponential_integral(D[1:,1:], T, mode=0)

        integral = P@D_int@np.linalg.inv(P)
    
    return integral

def covariance_matrix_exponential_integral_using_diagonalization(mat, T, limit=False):
    '''
    This function calculates the integral:
        
        (1/T) * int_0^T (T-\tau)e^{M\tau} d\tau.
    
    In the limit case T >> tau_n (limit = True), this simplifies to 
    
        int_0^infinty e^{M\tau} d\tau.
        
    This integral can be used to calculate the standard deviation for an averaging
    time T (limit = False):
        
        sigma_{nT}^2 = (2/T)*(occupancies@integral@(occupancies*P_steady))
        
    or to calculate the correlation time tau_n:
        
        tau_n = (occupancies@integral@(occupancies*P_steady))/sigma_n^2
    '''
    
    eigenvalues, eigenvectors = np.linalg.eig(mat)   
    order = np.argsort(np.abs(eigenvalues))
    eigs = eigenvalues[order]
    P = eigenvectors[:,order]
    # print('sum of eigenvectors elements:\n', np.sum(P, axis=0), np.sum(P))
    # print('first row of inverse eigenvector matrix:\n', np.linalg.inv(P)[0], '\n')
    
    ## Test: accuracy of diagonalization, PDP^-1 should be equal to the original matrix.
    # print(mat, '\n')
    # D = np.diag(eigs)
    # print(P@D@np.linalg.inv(P), '\n')
    # print(np.real(P@D@np.linalg.inv(P)), '\n')
    # print(np.abs(P@D@np.linalg.inv(P)), '\n')
    ## --------------------------------------------------------------------------------
    
    if eigs.dtype=='complex128':
        print('Warning: eigenvalues are complex!')
        # print(eigs)
        
        # eigs = np.real(eigs)#+np.round(np.imag(eigs),10)*1j
        # P = np.real(P)#+np.round(np.imag(P),10)*1j
        # print(eigs)

    D_int = np.zeros((len(eigs),len(eigs))).astype(eigs.dtype)
    if limit:
        D_int[1:,1:] = np.diag(-1/eigs[1:])
    else:
        D_int[1:,1:] = np.diag((-1/eigs[1:])*(1+(1-np.exp(eigs[1:]*T))/(eigs[1:]*T)))
        # So in the limit of large T this indeed becomes np.diag(-1/eigs[1:])
        
    integral = P@D_int@np.linalg.inv(P)
    return np.real(integral)
    

def covariance_matrix_exponential_integral(mat, P_steady, T, limit=False):
    '''
    This function calculates the integral:
        
        (1/T) * int_0^T (T-\tau)e^{M\tau} d\tau.
    
    In the limit case T >> tau_n (limit = True), this simplifies to 
    
        int_0^infinty e^{M\tau} d\tau.
        
    This integral can be used to calculate the standard deviation for an averaging
    time T (limit = False):
        
        sigma_{nT}^2 = (2/T)*(occupancies@integral@(occupancies*P_steady))
        
    or to calculate the correlation time tau_n:
        
        tau_n = (occupancies@integral@(occupancies*P_steady))/sigma_n^2
    '''
    
    correction_mat = np.repeat([P_steady], len(mat), axis=0).T
    new_mat = mat + correction_mat
    new_mat_inv = np.linalg.inv(new_mat)-correction_mat
    if limit:
        integral = -new_mat_inv
    else:
        integral = -(1/T)*(new_mat_inv@new_mat_inv)@(np.eye(len(mat))-sp.sparse.linalg.expm(mat*T))-(new_mat_inv)

    return integral

def mean_variance_and_correlation_time_occupancy(mat, T, P_steady, occupancies, limit=True, return_all=True, diagonalize=False):
    if diagonalize:
        integral = covariance_matrix_exponential_integral_using_diagonalization(mat, T, limit=limit)
    else:
        integral = covariance_matrix_exponential_integral(mat, P_steady, T, limit=limit)
    occupancy_integral = (occupancies@integral@(occupancies*P_steady))
    
    mu = occupancies@P_steady
    variance = (2/T)*occupancy_integral
    
    if variance<0:
        print(f'variance: {variance}')
        print(f'integral: {integral}')
        # raise ValueError(f'Variance should be nonnegative! {variance}')
    
    if return_all:
        sigma_squared = (occupancies-mu)**2@P_steady
        tau = occupancy_integral/sigma_squared
        # Note that the calculation for tau only makes sense in the limit case.
        
        return mu, sigma_squared, variance, tau
    
    return mu, variance

if __name__ == "__main__":
    from FullMatrix import FullRatesAndEnergies, FullStates, FullMatrix

    T = 1050000000.0000001
    
    # Variance higher than 1
    J   = 14.0 
    eb  = 3.0 
    htf = -17.16
    L   = 4
    sites = np.arange(L)
    periodic = True

    
    
    R = FullRatesAndEnergies(htf, eb, J, beta=1)
    S = FullStates(L)
    M = FullMatrix(S, R, sites, periodic)
    mat = M.create_matrix(sparse=False)
    P_steady = M.calculate_steady_state()
    occupancies = M.calculate_state_occupancies()
    
    integral = covariance_matrix_exponential_integral_using_diagonalization(mat, T, limit=False)
    old_variance = (2/T)*(occupancies@integral@(occupancies*P_steady))
    
    mean, variance = mean_variance_and_correlation_time_occupancy(mat, T, P_steady, occupancies, limit=False, return_all=False)
    
    print('mean:', mean)
    print('old variance:', old_variance)
    print('new variance:', variance)