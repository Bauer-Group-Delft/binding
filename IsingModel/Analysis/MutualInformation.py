# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from scipy import special
import pickle

import os
import sys
sys.path.append('../AnalyticalModel')
sys.path.append('../Parameters')
sys.path.append('../Parameters/Bcd_data')

from FullMatrix import FullStates 
from Metadata import read_metadata
from FileReader import read_data
from bcd import read_bcd_data, discretize_concentrations

class MutualInformation:
    def __init__(self, file_name, analytical=True, stochastic=True, N_bins_gaussian=1000, N_bins_histogram=400, conc_samples=60, write=True):
        file = os.path.join(os.path.realpath('__file__/../..'), f'Data/{file_name}') 
        
        self.analytical = analytical
        self.stochastic = stochastic
        self.N_bins_gaussian = N_bins_gaussian
        self.N_bins_histogram = N_bins_histogram
        self.conc_samples = conc_samples
        self.read_required_info(file)
        
        self.MI_results = self.calculate_all_mutual_information()
        
        if write:
            f = open(file+f'_mutual_information_manuscript_{N_bins_gaussian}_{N_bins_histogram}', 'wb')
            pickle.dump(self.MI_results, f)
            f.close()
    
    def read_required_info(self, file):
        # Metadata
        self.Grid_prms, self.Sim_prms, self.Varied_prms = read_metadata(file)
    
        f = open(file+'_shift_compensation', 'rb')
        self.shift_compensation = pickle.load(f)
        f.close()
        
        self.J_len = len(self.Varied_prms.J_list)
    
        # Analytical results
        if self.analytical:
            f = open(file+'_steady_states', 'rb')
            all_steady_states = pickle.load(f)
            f.close()
        
            f = open(file+'_means_and_vars', 'rb')
            all_means_and_vars = pickle.load(f)
            f.close()

            self.steady_states  = np.array([[all_steady_states[f'J{J}eb{eb}htf{htf}'][1] for htf in self.Varied_prms.htf_list] for J, eb in zip(self.Varied_prms.J_list, self.Varied_prms.eb_list)])
            self.means_and_vars = np.array([[all_means_and_vars[f'J{J}eb{eb}htf{htf}']   for htf in self.Varied_prms.htf_list] for J, eb in zip(self.Varied_prms.J_list, self.Varied_prms.eb_list)])
        
        # Stochastic results
        if self.stochastic:
            self.data = [np.zeros((self.J_len, len(self.Varied_prms.htf_list), self.Sim_prms.runs))]
            for j, (J, eb) in enumerate(zip(self.Varied_prms.J_list, self.Varied_prms.eb_list)):
                for h, htf in enumerate(self.Varied_prms.htf_list):
                    dataset = read_data(file+f'_J{J}eb{eb}htf{htf}', self.Sim_prms.runs)
        
                    while len(self.data)<len(dataset):
                        self.data.append(np.zeros((self.J_len, len(self.Varied_prms.htf_list), self.Sim_prms.runs)))
            
                    for i in range(len(dataset)):
                        if self.shift_compensation.N_expr==None:
                            expression = np.array(dataset[i]['mean binding site occupancy'], dtype='float')
                        else:
                            expression = np.array([np.sum(np.array(x.split(' '), dtype='float')[self.shift_compensation.N_expr:]) for x in dataset[i]['fraction of time with N TFs bound']])
                            
                        self.data[i][j,h,:] = expression
        
        # Bicoid data
        TGB, meanTGB, varTGB, xs_full, mask = read_bcd_data('../../Parameters/Bcd_data/TG_normBcd', normalization_mode=1)
        self.rel_conc_bins, self.threshold_pos, _ = discretize_concentrations((TGB, meanTGB, varTGB, xs_full, mask), number_of_samples=self.conc_samples)
        
        self.meanTGB = meanTGB[mask]
        self.varTGB  = varTGB[mask]
        self.xs      = xs_full[mask]

    def mutual_information_dependence_on_J(self, distribution_mode, args):
        # C = Hunchback gene expression, S = Bicoid concentration, X = position.
        I_dict = {'P_S' : [], 'I_CS' : np.zeros(self.J_len),
                  'P_X' : [], 'I_CX' : np.zeros(self.J_len)}

        # P(S|X) is obtained from literature.
        # For every position (X), the possible concentrations (S) are assumed to be gaussian distributed. 
        # The mean and variance of the Gaussian distributions are obtained from expreiments in literature.
        P_S_given_X = gaussian_conditional_probability_distributions(self.rel_conc_bins, list(zip(self.meanTGB, self.varTGB)))

        # P(X) is a unifrom distribution. Each position has the same probability
        P_X         = 1/len(P_S_given_X)
        # P(S)
        P_S         = np.sum(P_S_given_X*P_X, axis=0)
        
        plt.plot((self.rel_conc_bins[:-1]+self.rel_conc_bins[1:])/2, P_S)
        plt.xlabel('Bicoid concentration, S')
        plt.ylabel('Total probability, P(S)')

        for j in range(self.J_len):

            # P(C|S) is obtained from analytical calculations or simulations
            if distribution_mode == 'steady state':
                steady_states, Grid_prms, N_expr = args
                P_C_given_S = steady_state_conditional_probability_distributions(steady_states[j], Grid_prms, N_expr)
            elif distribution_mode == 'gaussian':
                means_and_vars, N_bins = args
                I_dict['N bins'] = N_bins
                P_C_given_S = gaussian_conditional_probability_distributions(np.linspace(0,1,N_bins+1), means_and_vars[j])
            elif distribution_mode == 'histogram':
                data, N_bins = args
                I_dict['N bins'] = N_bins
                P_C_given_S = histogram_conditional_probability_distributions(np.linspace(0,1,N_bins+1), data[j])
            else:
                raise ValueError(f'Invalid distribution mode: {distribution_mode}')

            # P(C|X)
            P_C_given_X = P_S_given_X@P_C_given_S

            # Calculation for I(C,S)
            P_C_and_S   = P_C_given_S*np.tile(P_S, (np.size(P_C_given_S[0]),1)).T
            P_C         = np.sum(P_C_and_S, axis=0)
            I_dict['P_S'].append([P_C_given_S, P_C_and_S, P_C])
            I_dict['I_CS'][j] = calculate_mutual_information([P_C_given_S, P_C_and_S, P_C])

            # Calculation for I(C,X)
            P_C_and_X   = P_C_given_X*P_X
            P_C         = np.sum(P_C_and_X, axis=0)
            I_dict['P_X'].append([P_C_given_X, P_C_and_X, P_C])
            I_dict['I_CX'][j] = calculate_mutual_information([P_C_given_X, P_C_and_X, P_C])

        return I_dict

    def calculate_all_mutual_information(self):        
        ## Calculate mutual information
        I_dict = {'variables' : self.Varied_prms, 'steady state' : None, 'gaussian' : None, 'histogram' : None}
        
        if self.analytical:
            I_dict['steady state'] = self.mutual_information_dependence_on_J('steady state', (self.steady_states, self.Grid_prms, self.shift_compensation.N_expr))
            I_dict['gaussian']     = self.mutual_information_dependence_on_J('gaussian', (self.means_and_vars, self.N_bins_gaussian))

        if self.stochastic:
            I_dict['histogram'] = {}
            for i in range(len(self.data)):
                I_dict_h = self.mutual_information_dependence_on_J('histogram', (self.data[i], self.N_bins_histogram))
                if i == 0:
                    for key, val in I_dict_h.items():
                        I_dict['histogram'][key] = [val]
                else:
                    for key, val in I_dict_h.items():
                        I_dict['histogram'][key].append(val)
        
        return I_dict
    
    
## Probability distributions
def steady_state_conditional_probability_distributions(steady_states, grid_prms, N=None): 
    sites = grid_prms.sites.flatten()
    P_expression_given_htf = np.zeros((len(steady_states),len(sites)+1))

    for i in range(len(steady_states)):
        S = FullStates(grid_prms.L)
        P_expression_given_htf[i,:] = S.compact_states(steady_states[i], sites)
   
    if N == None:
        return P_expression_given_htf
    else:
        return np.array([[np.sum(P_expression_given_htf[i,:N]), np.sum(P_expression_given_htf[i,N:])] for i in range(len(steady_states))])

def gaussian_conditional_probability_distributions(y_bins, means_and_vars):
    N_bins = len(y_bins)-1
    P_y_given_x = np.zeros((len(means_and_vars),N_bins))

    for i in range(len(means_and_vars)):
        mu, var = means_and_vars[i]
        # Integrate over a segment dx to get the total probability of being within the segment. 
        # Note that since we renormalize later, we can leave out the factor 1/2.
        segment_integrals = special.erf((y_bins[1:]-mu)/np.sqrt(2*var))-special.erf((y_bins[:-1]-mu)/np.sqrt(2*var)) 
        # Renormalize to make sure the integral remains 1, since we only consider occupancies in the range [0,1], while
        # gaussian distributions have as domain the entire real number line.
        P_y_given_x[i,:] = segment_integrals/np.sum(segment_integrals)
        
    return P_y_given_x

def histogram_conditional_probability_distributions(y_bins, data):
    N_bins = len(y_bins)-1
    P_y_given_x = np.zeros((len(data),N_bins))

    for i in range(len(data)):
        hist_values, bin_edges = np.histogram(data[i], bins=y_bins, density=True)
        P_y_given_x[i,:] = hist_values*(1/N_bins) # Discretized distribution, adds up to 1
        
    return P_y_given_x


## Mutual Information
def calculate_mutual_information(Ps):
    P_y_given_x, P_y_and_x, P_y = Ps

    # If P_y_and_x is nonzero, so will P_y_given_x and P_y.
    # We should only take the logarithms of those terms, since all others are zero.
    # (We assume 0*log(0)=0, as is done in 'Element of information theory, second edition' 
    # by Thomas M. Cover and Joy A. Thomas.)

    nonzero_terms = np.where(P_y_and_x!=0)
    I = np.sum(P_y_and_x[nonzero_terms]*np.log2(P_y_given_x[nonzero_terms]/P_y[nonzero_terms[1]]))
    return I