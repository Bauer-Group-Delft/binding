# -*- coding: utf-8 -*-

import numpy as np
import itertools as it

import ParametersFromLiterature as prms

class GridParameters():
    def __init__(self, dim, L, sites, IC_mode, periodic=False, do_htf_update=False, beta=prms.beta, grid_cell_size=prms.grid_cell_size, R=prms.R):
        '''
        Parameters
        ----------
        dim : integer
            Dimension of the grid. Supported dimensions: 1D (dim=1) and 2D (dim=2).
            
        L : integer
            The length of the grid in grid cells in each dimension.
            
        sites : 1D or 2D numpy array containing integers
            An array containing the positions of the binding sites. The positions
            can be given using linear (1D) indices or indices with the same dimension (nD)
            as the grid.
            In the case of 1D indices, sites will be a 1D array with length equal to the number
            of binding sites. In the case of nD indices, sites will be a 2D array where the
            length along the first dimension is equal to the number of binding sites and the
            length along the second dimension is equal to the dimension of the grid.

        IC_mode : string
            Mode for generating initial conditions. Options: 'fill fraction', 'probability', 'grid'.
            When the mode is 'fill fraction', the grid is filled randomly to obtain the given fill fraction (fraction of sites filled).
            When the mode is 'probability', a grid is randoly chosen from a list of grids according to the given probabilities.
            When the mode is 'grid', the given grid is used as initial condition.
            
        periodic : boolean
            Indicates whether periodic boundary conditions should be used (True) or not (False). The default is False.
            
        do_htf_update : boolean
            Indicates whether to change the htf value when molecules appear on or 
            disappear from the grid. The default is False.
            
        beta : float, optional
            Multiplication factor to convert energies to units of $k_BT$. The default is prms.beta.
            
        grid_cell_size : float, optional
            The size of a single grid cell in micrometer. The default is prms.grid_cell_size.
            
        R : float, optional
            The size of the nucleus in micrometer. The default is prms.R.

        Raises
        ------
        ValueError
            A ValueError is raised when the binding site indices do not have a valid dimension 
            (1D or dimension equal to that of the grid).

        Returns
        -------
        None.

        '''
        
        self.dim = dim        
        self.L = L
        
        if len(np.shape(sites)) == 1:
            self.sites = self.convert_indices(sites)
        elif np.shape(sites)[1]==dim:
            self.sites = sites
        else:
            raise ValueError('Invalid binding site indices.')
        
        self.beta            = beta
        self.periodic        = periodic
        self.do_htf_update   = do_htf_update
        self.IC_mode         = IC_mode

        V_site               = (grid_cell_size)**3          # um^3
        V_nucleus            = 4/3*np.pi*R**3               # um^3
        self.N_sites_tot     = V_nucleus/V_site
    
        self.time_conversion = prms.time_conversion(self.dim, grid_cell_size)
            

    def convert_indices(self, site_indices):
        '''
        The convert_indices method converts linear indices into indices with the
        same dimension as the grid.

        Parameters
        ----------
        site_indices : 1D numpy array contaning integers
            A 1D array containing the linear indices of the sites.

        Returns
        -------
        converted_site_indices : 2D numpy array containing integers
            An array containing the converted indices of each site.

        '''
        
        if self.dim == 1:
            return np.expand_dims(site_indices, axis=1)
        elif self.dim == 2:
            rows = (site_indices/self.L).astype(int)
            cols = (site_indices%self.L).astype(int)
            return np.array([rows,cols]).T
        else:
            ValueError('Grid dimension not supported: {dim}.')


class RunParameters():
    def __init__(self, htf, eb, J, IC_prms, grid_prms):
        '''
        Parameters
        ----------
        htf : float
            The energy that determines the on-rate. If grid_prms.do_htf_update is True,
            this value will be used to calculated the htf values for each different number
            of TFs bound. The supplied htf will in that case be taken to indicate the htf 
            when no TFs are bound, though the final htf value for 0 TFs bound might deviate
            slightly from the supplied value (see comment below).
            
        eb : float
            The binding site energy.
            
        J : float
            The binding energy between two neighboring TFs.
            
        IC_prms : list
            Parameters required for the chosen IC_mode.
            When the mode is 'fill fraction', IC_prms contains a float indicating the fill fraction.
            When the mode is 'probability', IC_prms contains a list containing the grid configurations (strings of 0s and 1s) and a list with the corresponding probabilities.
            When the mode is 'grid', IC_prms contains string of 0s and 1s indicating the grid configuration.

        Returns
        -------
        None.

        '''
        
        self.eb         = eb
        self.J          = J
        self.IC_prms    = IC_prms
        
        if grid_prms.do_htf_update:
            fill_frac   = prms.htf_to_c(grid_prms.beta, htf)
            N_molecules = round(fill_frac*grid_prms.N_sites_tot)
            self.htf    = [prms.c_to_htf(grid_prms.beta, (N_molecules - N_present)/(grid_prms.N_sites_tot - N_present))\
                                 for N_present in range(min(N_molecules, grid_prms.L**grid_prms.dim)+1)]
                
            # Note that here we cast N_molecules to an integer to prevent N_bound from becoming
            # larger than N_molecules (which would results in negative concentrations and errors
            # with the logarithm in the htf calculation). This means that the actual htf value 
            # in the simulation will be slightly lower than the set htf value and the range of
            # possible htf values is now discrete rather than continuous. These effects should
            # only be noticable for very low concentrations.
            
        else:
            self.htf    = htf


class SimulationParameters():
    def __init__(self, runs, T, T_trans=None, max_timesteps=int(1e6)):
        '''
        Parameters
        ----------
        runs : integer
            Number of runs of the simulation.
            
        T : float
            Runtime in time units.
            
        T_trans : float, optional
            Time of transient in time units, which should be removed when calculating model results. The default is None.
            
        max_timesteps : integer, optional
            Maximum number of timesteps for which to run the simulation. The default is int(1e6).

        Returns
        -------
        None.

        '''
        
        self.runs             = runs
        self.T                = T
        self.T_trans          = T_trans 
        self.max_timesteps    = max_timesteps


class VariedParameters():
    def __init__(self, htf_list, eb_list, J_list, IC_prms_list):
        '''             
        Parameters
        ----------
        htf_list : list of floats.
            List of the different htf values used for the simulations. 
            For each htf in htf_list, the simulations will be repeated for all 
            sets of eb and J. 
        eb_list : list of floats. 
            List of the different binding energies used for the simulations.
        J_list : list of floats.
            List of the different clustering strengths used for the simulations.
            J_list must have the same length as eb_list, where values with the same 
            index form a set. 
        IC_prms_list : dictionary of lists.
            Dictionary containing the different IC_prms used for the different simulations. 
            IC_prms_list contains IC_prms for every combination of htf, eb and J. 
            The keys are given by f'J{J}eb{eb}htf{htf}'.

        Returns
        -------
        None.

        '''
        
        self.htf_list         = htf_list
        self.J_list           = J_list
        self.eb_list          = eb_list
        self.IC_prms_list     = IC_prms_list
        
        self.zipped           = list(it.product(zip(J_list, eb_list), htf_list))