# -*- coding: utf-8 -*-

import numpy as np
import scipy as sp
import scipy.ndimage

from Grid import TF_Grid

class TF_Grid_1D(TF_Grid):
    def __init__(self, run_prms, grid_prms):
        self.shape = [grid_prms.L] # TODO: would it be easier to make this [1,L], more methods can be general?
        super().__init__(run_prms, grid_prms)

    def create_EnvironmentMatrix(self):
        kernel = [1,0,1]
        if self.periodic:
            return sp.ndimage.convolve1d(self.Grid, kernel, mode='wrap')
        else:
            return sp.ndimage.convolve1d(self.Grid, kernel, mode='constant', cval=0)
    
    def create_BindingSiteMatrix(self):
        BindingSiteMatrix = np.zeros(self.shape)
        BindingSiteMatrix[self.sites[:,0]]=1
        return BindingSiteMatrix  
    
    
    def choose_Pos(self, randomNumber, rates):
        randomNumber, col_idx = self.choose_Idx(randomNumber, rates)
        return [col_idx]
    
    def update_GridPos(self, SiteIndex, increment):
        self.Grid[SiteIndex[0]] += increment
    
    def update_EnvironmentMatrix(self, neighbors, increment):
        # Note that a for loop is required here to prevent issues with small grids,
        # where a site might have the same neighbor on both sides due to the periodic BCs.
        for neighbor in neighbors:
            self.PS[neighbor[0]] += increment

    