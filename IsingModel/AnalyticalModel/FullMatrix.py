# -*- coding: utf-8 -*-

import numpy as np
import scipy as sp

#%% Rates and Energies

'''
The FullRatesAndEnergies class contains methods of calculating the rates and energies for the full matrix.
The mode supplied to this class will determine which of these methods is used. 
The parameters htf, eb and J must be supplied. Using the energies, the steady state probabilities
for the different grid states can be calculated.
'''

class FullRatesAndEnergies:
    def __init__(self, htf, eb, J, beta, adapt_conc=False):
        # If adapt concentration is True, htf should be a list of htf
        # values for every valid number of molecules present on the grid.
        
        self.beta = beta
        self.htf  = htf
        self.eb   = eb
        self.J    = J
        self.adapt_conc = adapt_conc
        
        self.kon  = np.exp(beta*htf)
        self.koff = lambda eb, N: np.exp(-beta*(eb+N*J))
        
            
    def count_neighbors(self, s, i, periodic):
        # Retruns the number of neightbors of site i in state s.
        # periodic is a boolean variable which indicates whether periodic BCs 
        # are used or not.

        L = len(s)
        if periodic or (i > 0 and i < (L-1)):
            n = s[i-1] + s[(i+1)%L] # Number of neighbors, using periodic BCs
        else: # non-periodic BCs
            if L>1:
                if i==0:
                    n = s[1]
                elif i==L-1:
                    n = s[L-2]
            else: 
                # If there is only one cell in the grid and the BCs are not periodic, 
                # there are no neighbors.
                n=0
        return n

    def determine_rates(self, s, i, binding_sites, periodic):
        '''
        The function determine_rates needs to be supplied with the current state s, the 
        index i of the grid site that changes, a list of indices of the binding sites, 
        and and a boolean indicating whether the grid is periodic. Whether k_on or k_off is 
        retruned depends on the state of grid site i (empty --> k_on, filled --> k_off).
        '''

        if type(s[0])==str:
            s = [int(site) for site in s]
        
        if s[i]: # s[i] == 1
            e = int(i in binding_sites)*self.eb  
            n = self.count_neighbors(s, i, periodic)
            return self.koff(e, n)
        else:  # s[i] == 0
            if self.adapt_conc:
                N_present = np.sum(s)
                if N_present >= len(self.kon):
                    return 0
                else:
                    return self.kon[N_present]
            else:
                return self.kon
            
            
    def determine_energies(self, s, binding_sites, periodic):
        '''
        The function determine_energies needs to be supplied with the current state s, 
        a list of indices of the binding sites, and and a boolean indicating whether 
        the grid is periodic.
        '''
        
        if type(s[0])==str:
            s = [int(site) for site in s]
        
        N_present = np.sum(s)
        if N_present == 0:
            return 0
        elif self.adapt_conc and N_present >= len(self.htf):
            # If the number of bound molecules is too large, the state has infinite energy.
            return np.inf

        # Add -htf for each filled site
        if self.adapt_conc:
            E = -np.sum(self.htf[:N_present])
        else:
            E = -self.htf*N_present

        # Add -eb for each of the filled binding sites
        E += -self.eb*np.sum(s[binding_sites])

        # Add -J for each pair of neighbors       
        if periodic:
            E += -np.sum(s*np.concatenate([s[1:],[s[0]]]))*self.J
        else:
            E += -np.sum(s[:-1]*s[1:])*self.J
            
        return E
    
    def calculate_steady_state_probabilities(self, states, binding_sites, periodic):  
        E = np.zeros(len(states))
        for i,s in enumerate(states):
            E[i] = self.determine_energies(s, binding_sites, periodic)
        p_steady = np.exp(-self.beta*E)
        return p_steady/np.sum(p_steady)
    
    
#%% States
    
'''

The FullStates class contains all methods relating to generating the states of the full
matrix and calculating their properties. 

'''

class FullStates:
    def __init__(self, L):
        self.L = L
        self.Nstates = 2**self.L
        self.str_states, self.int_states = self.create_states()
    
    def create_states(self):
        # Creates an array containing all possible states of a grid. 
        # The mode determines in what way the states are represented (str or int).

        str_states = np.array([f'{i:0{self.L}b}' for i in range(self.Nstates)])
        int_states = np.array([np.array(list(s)).astype(int) for s in str_states])
        return str_states, int_states # return states in both formats
    
            
    def calculate_state_occupancies(self, sites, normalize=True):
        # Calculate the mean occupancy of sites for each state
        
        sites = np.array(sites)
        state_occupancies = np.sum(self.int_states[:, sites], axis=1) # only works if sites is 1D numpy array
        if normalize:
            state_occupancies = state_occupancies/len(sites)
        return state_occupancies
    
    def compact_states(self, v, sites):
        # Combine states with the same number of TFs bound on the given sites.
        
        n = self.calculate_state_occupancies(sites, normalize=False)
        
        if len(v.shape)==1:
            v_compact = np.zeros(len(sites)+1)
            for i in range(len(n)):
                v_compact[n[i]] += v[i]
        else:
            v_compact = np.zeros((len(v), len(sites)+1))
            for i in range(len(n)):
                v_compact[:,n[i]] += v[:,i]
            
        return v_compact   
         
    
    def create_initial_condition(self, IC_mode, IC_prms):
        if IC_mode == 'fill fraction':
            # Create intial condition based on fill fraction, 
            # each state with the correct number of TFs has equal probability.
            v0 = self.IC_fill_frac(IC_prms[0])
            
        elif IC_mode == 'probability':
            # Use given probabilities.
            states, probabilities = IC_prms
            v0 = self.IC_probability(states, probabilities)
        
        elif IC_mode == 'grid':
            # Use a given initial grid.
            grid = IC_prms[0]
            v0 = np.zeros(self.Nstates)
            v0[self.str_states==''.join(grid)] = 1
            
        else:
            raise ValueError('Unknown initial condition mode: {self.IC_mode}.')
        
        return v0
    
    def IC_fill_frac(self, fill_frac):
        v0 = np.zeros(self.Nstates)
        v0[np.sum(self.int_states, axis=1)==round(fill_frac*self.L)] = 1
        return v0/np.sum(v0) # Normalization
        
    def IC_probability(self, states, probabilities):
        if np.all(states==self.str_states):
            v0 = probabilities
        else:
            v0 = np.zeros(self.Nstates)
            for i,state in enumerate(states):
                if state in self.str_states:
                    v0[state==self.str_states] += probabilities[i]
                else:
                    raise ValueError('Invalid state in initial condition: {state}.')
                    
        return v0
    
    
#%% Full Matrix

'''

The FullMatrix class contains all methods relating to generating the full matrix in 
sparse or dense format.

'''

class FullMatrix:
    
    def __init__(self, states, rates, binding_sites, periodic_BCs):
        self.states = states
        self.rates  = rates

        self.binding_sites = binding_sites
        self.periodic_BCs  = periodic_BCs

        self.sparse_mat    = None
        self.dense_mat     = None
        self.steady_state  = None
        
    def create_matrix(self, sparse):
        if sparse:
            if np.any(self.sparse_mat == None):
                self.sparse_mat = self.create_sparse_matrix()
            return self.sparse_mat
            
        else:
            if np.any(self.dense_mat == None):
                self.dense_mat = self.create_dense_matrix()
            return self.dense_mat
           
    def create_dense_matrix(self):
        # Method for generating full matrix
        mat = [[0 for i in range(self.states.Nstates)] for j in range(self.states.Nstates)]
        for s_idx, s in enumerate(self.states.str_states):
            for i in range(self.states.L):
                new_state = s[:i] + str(int(not int(s[i]))) + s[i+1:]
                mat[int(new_state, base=2)][s_idx] = self.rates.determine_rates(s, i, self.binding_sites, self.periodic_BCs)
        
        # diagonal elements
        mat -= np.diag(np.sum(mat, axis=0))
        
        return mat
    
    def create_sparse_matrix(self):
        data = np.array([])
        rows = np.array([])
        cols = np.array([])
        for s_idx, s in enumerate(self.states.str_states):
            cols = np.concatenate((cols, [s_idx]*self.states.L)) # next L states are in the same column
            for i in range(self.states.L):
                new_state = s[:i] + str(int(not int(s[i]))) + s[i+1:]
                rows = np.concatenate((rows, [np.where(self.states.str_states==new_state)[0][0]]))
                data = np.concatenate((data, [self.rates.determine_rates(s, i, self.binding_sites, self.periodic_BCs)]))
    
            # diagonal element
            data = np.concatenate((data, [-np.sum(data[np.where(cols==s_idx)])]))
            cols = np.concatenate((cols,[s_idx]))
            rows = np.concatenate((rows,[s_idx]))
    
        mat = sp.sparse.csc_array((data,(rows,cols)), shape=(self.states.Nstates,self.states.Nstates))
        #print(mat.toarray())
        
        return mat
    

    def calculate_state_occupancies(self, normalize=True):
        return self.states.calculate_state_occupancies(self.binding_sites, normalize)
    
    def compact_states(self, v):
        return self.states.compact_states(v, self.binding_sites)

    def calculate_steady_state(self):
        if not np.any(self.steady_state == None):
            return self.steady_state
        
        self.steady_state = self.rates.calculate_steady_state_probabilities(self.states.int_states, self.binding_sites, self.periodic_BCs)
        return self.steady_state