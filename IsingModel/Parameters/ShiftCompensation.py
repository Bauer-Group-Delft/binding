# -*- coding: utf-8 -*-

import math as m
import numpy as np

import sys
sys.path.append('../AnalyticalModel/')

import ParametersFromLiterature as prms
from Bisection import bisection

from GetExpression import get_expression

class ShiftCompensation():
    def __init__(self, grid_parameters, J_list, htf_05, N_expr=None, T=None, IC_prms=None, tol=0.002):

        self.grid_prms = grid_parameters
        self.J_list    = J_list
        self.eb_list   = None
        self.htf_05    = htf_05
        self.N_expr    = N_expr
        self.T         = T 
        self.IC_prms   = IC_prms
        self.tol       = tol
        
        if self.N_expr==None:
            if grid_parameters.do_htf_update:
                fill_frac   = prms.htf_to_c(self.grid_prms.beta, self.htf_05)
                N_molecules = int(fill_frac*self.grid_prms.N_sites_tot)
                # Maximum number of molecules present on the field at 0.5 occupancy:
                N_present   = m.ceil((self.grid_prms.L**self.grid_prms.dim)/2)
                # Minimum (most negative) htf at 0.5 occupancy:
                adapted_htf = prms.c_to_htf(self.grid_prms.beta, (N_molecules - N_present)/(self.grid_prms.N_sites_tot - N_present))
                # To balance out the htf, we need eb = -htf, so if this is the minimum htf, then the maximum eb is given by:
                self.eb_max = -adapted_htf
            else:
                self.eb_max = -self.htf_05 
        else:
            self.eb_max = -self.htf_05-(1/self.grid_prms.beta)*np.log(2**(1/self.N_expr)-1)
        
    def find_eb(self):
        if np.any(self.eb_list == None):
            eb_list = []
            for j, J in enumerate(self.J_list):
                func_args = {'htf' : self.htf_05, 'J' : J, 'IC_prms' : self.IC_prms, \
                             'grid_prms' : self.grid_prms, 'T' : self.T, 'N' : self.N_expr}
                eb_list.append(bisection([0,self.eb_max], get_expression, func_args, tol=self.tol))  
            self.eb_list = np.array(eb_list)
        return self.eb_list
