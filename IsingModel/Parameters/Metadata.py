# -*- coding: utf-8 -*-

#%% Read and write metadata

import pickle

def write_metadata(file, Grid_prms, Sim_prms, Varied_prms):
    f = open(file+'_metadata', 'wb')
    ## Note that previous contents of the file will be removed.
    pickle.dump(Grid_prms, f)
    pickle.dump(Sim_prms, f)
    pickle.dump(Varied_prms, f)
    f.close()
    
def read_metadata(file):
    f = open(file+'_metadata', 'rb')
    Grid_prms   = pickle.load(f)
    Sim_prms    = pickle.load(f)
    Varied_prms = pickle.load(f)
    f.close()
    
    return Grid_prms, Sim_prms, Varied_prms