# -*- coding: utf-8 -*-

import numpy as np
import scipy as sp
import scipy.signal

from Grid import TF_Grid

class TF_Grid_2D(TF_Grid):
    def __init__(self, run_prms, grid_prms):
        self.shape = [grid_prms.L, grid_prms.L]
        super().__init__(run_prms, grid_prms)
        
    def create_EnvironmentMatrix(self):
        kernel = [[0,1,0],[1,0,1],[0,1,0]]
        if self.periodic:
            return sp.signal.convolve2d(self.Grid, kernel, mode='same', boundary='wrap')
        else:
            return sp.signal.convolve2d(self.Grid, kernel, mode='same', boundary='fill', fillvalue=0)

    def create_BindingSiteMatrix(self):
        BindingSiteMatrix = np.zeros(self.shape)
        BindingSiteMatrix[self.sites[:,0],self.sites[:,1]]=1
        return BindingSiteMatrix

    
    def choose_Pos(self, randomNumber, rates):
        randomNumber, row_idx = self.choose_Idx(randomNumber, np.sum(rates,1))
        randomNumber, col_idx = self.choose_Idx(randomNumber, rates[row_idx,:])
        return [row_idx,col_idx]
    
    def update_GridPos(self, SiteIndex, increment):
        self.Grid[SiteIndex[0],SiteIndex[1]] += increment
    
    def update_EnvironmentMatrix(self, neighbors, increment):
        # Note that a for loop is required here to prevent issues with small grids,
        # where a site might have the same neighbor on both sides due to the periodic BCs.
        for neighbor in neighbors:            
            self.PS[neighbor[0],neighbor[1]] += increment
    