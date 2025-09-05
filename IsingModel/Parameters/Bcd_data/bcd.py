import math as m
import numpy as np
import scipy as sp
import scipy.io

import os
import sys
sys.path.append('../.')

import ParametersFromLiterature as prms

def read_bcd_data(rel_path, normalization_mode=1):    
    # This is bcd data from Thomas Gregor, David Tank, Eric Wieschaus and Bill Bialek, kindly obtained from Thomas Gregor.
    # Reference to publication: Gregor T, Tank DW, Wieschaus EF, Bialek W. Probing the limits to positional information. 
    # Cell. 2007 Jul 13;130(1):153-64. doi: 10.1016/j.cell.2007.05.025. PMID: 17632062; PMCID: PMC2253670. 
    # The data can be used to replot figure 6 in this paper.
    
    file = os.path.join(os.path.realpath('__file__'), rel_path) 

    if normalization_mode == 1:
        TGB = sp.io.loadmat(file)['TGB_norm1']
    elif normalization_mode == 2:
        TGB = sp.io.loadmat(file)['TGB_norm2']
    else:
        raise ValueError(f'Invalid normalization mode: {normalization_mode}.')
        
    meanTGB = np.mean(TGB, axis=0)
    varTGB  = np.var(TGB, axis=0)
    
    xs_full = np.linspace(0.,1,200)
    mask    = np.where(np.logical_and(xs_full>=0.1, xs_full<=0.9))
    
    return TGB, meanTGB, varTGB, xs_full, mask

def discretize_concentrations(bcd_data, number_of_samples=60):  
    TGB, meanTGB, _, xs_full, mask = bcd_data
    
    rel_conc_bins = np.arange(m.ceil(np.max(TGB[:,mask])*number_of_samples)+1)/number_of_samples
    threshold_idx = np.argmin(np.abs(xs_full-prms.x_rel_threshold))
    threshold_pos = xs_full[threshold_idx]
    rel_threshold_conc = meanTGB[threshold_idx]
    
    return rel_conc_bins, threshold_pos, rel_threshold_conc