# -*- coding: utf-8 -*-

import numpy as np
from copy import deepcopy

class DataFile:
    
    def __init__(self, grid, file=None):
        self.grid = grid  
        
        if not file == None:
            self.file = file
                
            # Prepare datafile
            f_data = open(f'{self.file}_data.txt', 'a')
            f_data.write('run, initial grid state, final grid state, mean grid cell occupancy, mean binding site occupancy, fraction of time with N TFs bound, warning\n')
            f_data.close()
        
        self.setup()
        
    def setup(self):
        self.times  = []
        self.grids  = []
        self.N_bound = []
        self.N_bound_per_site = []
        
    def store_data(self, t):
        self.times.append(t)
        self.N_bound.append(self.grid.get_NBound())
        self.N_bound_per_site.append(self.grid.get_NBoundPerSite())
        self.grids.append(deepcopy(self.grid.get_Grid())) 
        # Use deepcopy to store the current state of the grid
        
    def write_run_data(self, run, transient=None, warning=None):
        if self.file == None:
            raise ValueError('No file assigned for writing data.')
        elif len(self.grids)==0:
            print('Warning: No grid data available.')
            return
            
        init_grid_str  = ''.join(np.array(self.grids[0], dtype='int').astype('str'))
        final_grid_str = ''.join(np.array(self.grids[-1], dtype='int').astype('str'))

        mean_grid_occupancy = time_mean_after_transient(self.get_N_present_per_site, transient)
        mean_site_occupancy = time_mean_after_transient(self.get_N_bound_per_site, transient)
        
        time_fraction_N_bound_str = ' '.join(self.get_N_bound_time_fraction(transient).astype('str'))
    
        line = f'{run+1}, {init_grid_str}, {final_grid_str}, {mean_grid_occupancy}, {mean_site_occupancy}, {time_fraction_N_bound_str}, '    

        if warning != None:
            line += warning
            
        f_data = open(f'{self.file}_data.txt', 'a')
        f_data.write(line + '\n')
        f_data.close()
        
    
    ## Getters
    def get_grids(self):
        return np.array(self.times), np.array(self.grids, dtype='int')
    
    def get_N_bound(self):
        return np.array(self.times), np.array(self.N_bound, dtype='int')
    
    def get_N_bound_per_site(self):
        return np.array(self.times), np.array(self.N_bound_per_site)
    
    def get_N_present(self):
        return np.array(self.times), np.sum(self.grids, axis=tuple(range(1,len(np.shape(self.grids)))))
    
    def get_N_present_per_site(self):
        return np.array(self.times), np.mean(self.grids, axis=tuple(range(1,len(np.shape(self.grids)))))
    
    def get_N_bound_time_fraction(self, transient=None):
        t, x = self.get_N_bound()
        if not transient == None:
            t, x = remove_transient(t, transient, x)
        return np.array([time_mean(t, x == N) for N in range(int(np.max(x))+1)])
        # Note that there will be less elements in time_fraction_N_bound than the number of binding sites 
        # if there is no moment where all bindign sites are bound simulateously. The missing entries should be 0.
    

## Functions for data processing

# Convert Gillespie output to output with regular timesteps.
def sample_output(t, output, T, samples):
    times = np.linspace(0,T,samples+1)
    idxs = np.searchsorted(t, times, 'right')
    return times, output[idxs-1]

def remove_transient(t, T_trans, x, start_at_0=False):
    idx = np.searchsorted(t, T_trans, side='right')
    if start_at_0:
        t = np.concatenate(([0], t[idx:]-T_trans))
    else:
        t = np.concatenate(([T_trans], t[idx:]))
    return t, x[idx-1:]

def time_mean(t, x): # weighted mean over time
    dts = t[1:]-t[:-1]
    # repeat dts to the match the shape of x
    dts_rep = np.array([np.full(np.shape(x[0]),dt) for dt in dts])
    return np.sum(dts_rep*x[:-1], axis=0)/(t[-1]-t[0])
    
def time_mean_after_transient(func, transient=None):
    t, x = func()
    if not transient == None:
        t, x = remove_transient(t, transient, x)
    return time_mean(t, x)