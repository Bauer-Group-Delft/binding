# -*- coding: utf-8 -*-

import random
import numpy as np
from abc import ABC, abstractmethod
from copy import deepcopy

# This is a grid template for grids where interactions between neighbors are important

'''
Inheritance stucture:

1.                  TF_Grid (Grid.py)

                      |          |
                      v          v
                      
2.  TF_Grid_1D (Grid1D.py)    TF_Grid_2D (Grid2D.py)

Methods defined on level 1: Init, setup, create_Grid, IC_fill_frac, IC_probability, 
update, update_Grid, choose_Idx, find_Neighbors, update_htf, update_Rates,
get_Grid, get_AppearRates, get_DisappearRates, get_NBound, get_NBoundPerSite, check_Grid.

Methods defined on level 2: Init, create_EnvironmentMatrix, create_BindingSiteMatrix, 
choose_Pos, update_GridPos, update_EnvironmentMatrix.

'''

class TF_Grid(ABC):
    def __init__(self, run_prms, grid_prms):
        '''
        The init method is called when a TF_Grid object is created.

        Parameters
        ----------
        htf : float or list of floats
            The energy that determines the on-rate. If do_htf_update is True, a list
            of floats must be supplied indicating the htf value when 0 to L TFs are
            bound.
            
        eb : float
            The binding site energy.
            
        J : float
            The binding energy between two neighboring TFs.
            
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
            
        periodic : boolean
            Indicates whether periodic boundary conditions should be used (True) or not (False).

        do_htf_update : boolean
            Indicates whether to change the htf value when molecules appear on or 
            disappear from the grid.
            
        IC_mode : string
            Mode for generating initial conditions. Options: 'fill fraction', 'probability', 'grid'.
            When the mode is 'fill fraction', the grid is filled randomly to obtain the given fill fraction (fraction of sites filled).
            When the mode is 'probability', a grid is randoly chosen from a list of grids according to the given probabilities.
            When the mode is 'grid', the given grid is used as initial condition.

        IC_prms : list
            Parameters required for the chosen IC_mode.
            When the mode is 'fill fraction', IC_prms contains a float indicating the fill fraction.
            When the mode is 'probability', IC_prms contains a list containing the grid configurations (strings of 0s and 1s) and a list with the corresponding probabilities.
            When the mode is 'grid', IC_prms contains string of 0s and 1s indicating the grid configuration.

        Raises
        ------
        ValueError
            If the dimensions of sites do not match the dimensions of the grid, 
            a ValueError is raised.

        Returns
        -------
        None.

        '''
        
        self.L                = grid_prms.L
        self.sites            = grid_prms.sites
        self.beta             = grid_prms.beta
        self.periodic         = grid_prms.periodic
        self.do_htf_update    = grid_prms.do_htf_update
        self.IC_mode          = grid_prms.IC_mode
        
        self.eb               = run_prms.eb
        self.J                = run_prms.J
        self.IC_prms          = run_prms.IC_prms
        
        if self.do_htf_update:
            self.adapted_htfs = run_prms.htf
            self.htf          = None
        else:
            self.htf          = run_prms.htf
        
        self.setup()    
    
    ## Functions for setting up the grid
    def setup(self):
        '''
        The setup method sets up the grid by creating a grid, an environment matrix,
        a binding site matrix and a matrix indicating where TFs are bound to binding sites.
        It then calculates the appearance and dissapeareance rates for each grid cell.
            
        This method is also used to reset the grid when reinitializing the simulation.
        Even if htf is changed, the grid will be created based on the same initial Fill fraction.
        In update rates, the htf is then updated based on the new number of TFs on the grid, 
        which should revert it back to its initial value.

        Returns
        -------
        None.

        '''
        
        self.Grid = self.create_Grid()        
        self.PS   = self.create_EnvironmentMatrix() # PS is important, as it stores the environment
        # print('PS:\n',self.PS)
        
        if len(self.sites):
            self.BindingSites  = self.create_BindingSiteMatrix()
            self.Bound         = np.where(self.BindingSites==1,self.Grid,0)
            # creates a matrix with the value of the grid at the binding sites and 0 at all other sites.
        else:
            self.BindingSites  = np.zeros(self.shape)
            self.Bound         = np.zeros(self.shape)
        
        # Calculate the rates for appearing and disappearing given the current configuration.
        self.update_Rates()

    def create_Grid(self):
        '''
        The create_Grid method creates an initial grid based on the the chosen initial
        condition mode.

        Returns
        -------
        grid : numpy array containing boolean values
            A matrix with the shape of the grid, that contains zeros or ones indicating 
            presence or absence of a TF, respectively.

        '''
        
        if self.IC_mode == 'fill fraction':
            # Fill a grid based on a given fill fraction
            grid = self.IC_fill_frac(self.IC_prms[0])
            
        elif self.IC_mode == 'probability':
            # Randomly choose a grid from a list of grids
            states, probabilities = self.IC_prms
            grid = self.IC_probability(states, probabilities)
        
        elif self.IC_mode == 'grid':
            # Use a given initial grid
            grid = self.IC_prms[0]
            
        else:
            raise ValueError('Unknown initial condition mode: {self.IC_mode}.')
        
        return grid
    
    def IC_fill_frac(self, fill_frac):
        '''
        The IC_fill_frac method randomly fills a grid matrix with transcription factors 
        to generate the initial condition. The number of TFs in the grid is determined
        based on the filling fraction, which is calculated as:
            P = exp(beta*htf)/(exp(beta*htf)+1).
            
        Parameters
        ----------
        fill_frac : float
            The fraction of frid cells that is filled. The fill fractions is used to 
            determine the number of TFs on the grid. The number of TFs is rounded to the nearest integer.
            
        Returns
        -------
        grid : numpy array containing boolean values
            A matrix with the shape of the grid, that contains zeros or ones indicating 
            presence or absence of a TF, respectively.
        '''
        
        N_sites = np.prod(self.shape)
        # Determine the number of TFs based on the fill fraction.
        N_TF = round(fill_frac*N_sites)
        # Determine where the TFs will be in the grid.
        TF_indices = random.sample(range(N_sites), N_TF) # randomly sample NTF indices.
        # Create a linear version of the grid and fill the chosen indices with ones.
        TF_List = np.zeros(N_sites)
        TF_List[TF_indices] = 1
        # Reshape the grid to a square.
        return np.reshape(TF_List, self.shape)
        
    def IC_probability(self, states, probabilities):
        '''
        The IC_probability method randomly chooses a grid state from the list of states
        based on the given probabilities.
        
        Parameters
        ----------
        states : list of strings
            A list of strings representing the possible grid states, one of which will
            be randomly chosen.
        probabilities : list of floats
            The probabilities of each of the grid states in states.
        
        Returns
        -------
        grid : numpy array containing boolean values
            A matrix with the shape of the grid, that contains zeros or ones indicating 
            presence or absence of a TF, respectively.
        '''
        
        str_grid = np.random.choice(states, p=probabilities)
        TF_List = np.array([int(s) for s in str_grid])
        return np.reshape(TF_List, self.shape)
    
    @abstractmethod
    def create_EnvironmentMatrix(self):
        '''
        The create_EnvironmentMatrix method calculates the number of TFs present in
        the environment of each site.
        This method should be overwritten in each subclass.

        Returns
        -------
        EnvironmentMatrix : numpy array containing integers
            A matrix with the shape of the grid, that contains integers indicating the
            number of neigboring sites that are filled for each site.
        
        '''
        pass
    
    @abstractmethod
    def create_BindingSiteMatrix(self):
        '''
        The create_BindingSiteMatrix method creates a matrix that indicates the position
        of the binding site.
        This method should be overwritten in each subclass.

        Returns
        -------
        BindingSiteMatrix : numpy array containing boolean values
            A matrix with the shape of the grid, that contains zeros or ones, 
            where ones indicate that the site is a binding site.
        '''
        pass
    
    
    ## Functions for updating the grid
    def update(self, *args, **kwargs):
        '''
        The update method updates the grid and rates after an event occurs.

        Returns
        -------
        None.

        '''
        
        # Update grid.
        self.update_Grid(*args, **kwargs)
        # Update appearance and disappearance rates.
        self.update_Rates()
        
    def update_Grid(self, randomNumber, rates, increment):
        '''
        The update method updates the grid, environment, rates and bound TFs after an
        event occurs. The way in which these should be adapted is indicated by increment.

        Parameters
        ----------
        randomNumber : float
            A random number that is unifromly distributed between 0 and the sum of rates.
        rates : numpy array
            An array with the same shape as the grid, containing for every position
            in the grid the appearance or disappearance rate of a TF.
        increment : integer (+1 or -1)
            Increment should be +1 for appearance and -1 for disappearance. 

        Returns
        -------
        None.

        '''
        
        # Update grid, randomly choose position to be updated.
        pos = self.choose_Pos(randomNumber,rates)
        self.update_GridPos(pos, increment)
        
        # Update neighbors for surrounding positions.
        neighbors = self.find_Neighbors(pos)
        self.update_EnvironmentMatrix(neighbors, increment)
        
        # Update Bound matrix
        self.Bound = np.where(self.BindingSites==1, self.Grid, 0)
    
    @abstractmethod
    def choose_Pos(self, randomNumber, rates):
        '''
        Choose a random position in the grid based on the random number.
        This method should be overwritten in each subclass.

        Parameters
        ----------
        randomNumber : float
            A random number that is unifromly distributed between 0 and the sum of rates.
        rates : numpy array
            An array with the same shape as the grid, containing for every position
            in the grid the appearance or disappearance rate of a TF.

        Returns
        -------
        Pos : 1D list containing integers
            The chosen position in the grid.

        '''
        pass
    
    
    def choose_Idx(self, randomNumber, rates):
        '''
        Choose a random row or column in the grid based on the random number.
        This method should be overwritten in each subclass.

        Parameters
        ----------
        randomNumber : float
            A random number that is unifromly distributed between 0 and the sum of rates.
        rates : numpy array
            An array with the same shape as the grid, containing for every row or column
            in the grid the appearance or disappearance rate of a TF.

        Raises
        ------
        ValueError
            If the chosen index is larger than the number of rows or columns, 
            a ValueError is raised.
            
            If the random number becomes negative, a ValueError is raised.

        Returns
        -------
        randomNumber : float
            The updated random number to be used for the next comparison.
        idx : integer
            The index of the chosen row or column.

        '''
        
        cumRates = np.cumsum(rates) + 1e-50 # TODO: Why do we add 1e-50?
        idx = np.searchsorted(cumRates, randomNumber, side='right') 
        # np.searchsorted uses bisection algorithm to efficiently search the list for 
        # the index of the first element that is larger or equal to a certain value.
        if idx > 0: 
            # Only subtract anything if there were rates smaller than the random number.
            randomNumber -= cumRates[idx-1]
        
        if idx >= len(cumRates):
            print('Error: Index chosen is invalid, make sure that random number is correct:', \
                  randomNumber < cumRates[-1])
            raise ValueError
            
        if randomNumber < 0:
            print('Error: Random number < 0:', randomNumber)
            raise ValueError
        
        return randomNumber, idx
    
    def find_Neighbors(self, site):
        '''
        The find_Neigbors method finds the neighbors of a given site.

        Parameters
        ----------
        site : list or 1D numpy array containing integers
            The site which should be updated.

        Returns
        -------
        neighbors : 2D numpy array containing integers
            The sites adjacent to SiteIndex

        '''
    
        neighbors = []
        for i in range(len(site)):
            for incr in [1, -1]:
                new_coords = np.array(deepcopy(site))
                new_coords[i] += incr # adapt only coordinate i.
                
                # Check whether new coordinates are valid.
                if self.periodic or (new_coords[i]>=0 and new_coords[i]<self.L):
                    neighbors.append(new_coords)
        
        if self.periodic:
            # Apply periodicity
            return (np.array(neighbors)+self.L)%self.L
            
        return np.array(neighbors)
    
    @abstractmethod
    def update_GridPos(self, SiteIndex, increment):
        '''
        The update_Grid method updates the grid at the site indicated by SiteIndex,
        by adding increment to the value at that site.
        This method should be overwritten in each subclass.
        
        TODO: These update functions can likely be made general for all subclasses by indexing
        using bollean arrays, though that might take up too much memory and not be
        beneficial.

        Parameters
        ----------
        SiteIndex : list or 1D numpy array containing integers
            The site which should be updated.
        increment : integer (+1 or -1)
            Increment should be +1 for appearance and -1 for disappearance. 

        Returns
        -------
        None.

        '''
        pass
        
    @abstractmethod
    def update_EnvironmentMatrix(self, neighbors, increment):
        '''
        The update_EnvironmentMatrix method updates the environment matrix at the
        sites indicated by neighbors, by adding increment to the value at those sites.
        This method should be overwritten in each subclass.

        Parameters
        ----------
        neighbors : 2D numpy array containing integers
            The sites of which the environment should be updated.
        increment : integer (+1 or -1)
            Increment should be +1 for appearance and -1 for disappearance. 

        Returns
        -------
        None.

        '''
        pass
    
    # required for modelling bath
    def update_htf(self):
        N_present = int(np.sum(self.Grid))
        self.htf = self.adapted_htfs[N_present]
        #print(self.htf)
    
    def update_Rates(self):
        '''
        The update_Rates method recalculates the appearance and disappearance rates.
        
        Parameters
        ----------
        SiteIndex : list or 1D numpy array containing integers
            The site which should be updated.
            
        neighbors : 2D numpy array containing integers
            The sites of which the environment should be updated.

        Returns
        -------
        None.

        '''
        
        # Recalculate htf for implementing a bath of TFs
        if self.do_htf_update:
            self.update_htf()
            # htf is only used for updating the rates (in particular the on-rate)
            # so we only need to update it here.
        
        self.AppearRates = (1-(self.Grid))*np.exp(self.beta*self.htf)
        self.DisappearRates = self.Grid*np.exp(-self.beta*(self.J*self.PS)) 
        self.DisappearRates = np.where(self.BindingSites==1, self.DisappearRates*np.exp(-self.beta*self.eb), self.DisappearRates)
    
    
        # TODO: Adapt this method such that rates are only recalulated for the sites 
        # indicated by neighbors. The method will then likely become subclass dependent.
        # sites = np.concatenate(SiteIndex, neighbors)
        # self.AppearRates[sites[:,0],sites[:,1]] = (1-(self.Grid[sites[:,0],sites[:,1]]))*np.exp(self.beta*self.htf)
        # self.DisappearRates[sites[:,0],sites[:,1]]  = self.Grid[sites[:,0],sites[:,1]]*np.exp(-self.beta*(self.J*self.PS[sites[:,0],sites[:,1]])) 
        # self.DisappearRates[sites[:,0],sites[:,1]]  = np.where(self.BindingSites[sites[:,0],sites[:,1]]==1, self.DisappearRates[sites[:,0],sites[:,1]]*np.exp(-self.beta*self.eb), self.DisappearRates[sites[:,0],sites[:,1]])

    
    ## Getters 
    def get_Grid(self):
        '''
        Returns
        -------
        Grid : numpy array containing boolean values
            A matrix with the shape of the grid, that contains zeros or ones indicating 
            presence or absence of a TF, respectively.

        '''
        return self.Grid
    
    def get_AppearRates(self):
        '''
        Returns
        -------
        AppearRates : numpy array of floats
            An array with the same shape as the grid, containing for every position
            in the grid the appearance rate of a TF.

        '''
        return self.AppearRates
    
    def get_DisappearRates(self):
        '''
        Returns
        -------
        DisappearRates : numpy array of floats
            An array with the same shape as the grid, containing for every position
            in the grid the disappearance rate of a TF.

        '''
        return self.DisappearRates    
    
    def get_NBound(self):
        '''
        Returns
        -------
        NBound : integer
            The number of transcription factors bound to a binding site.

        '''
        return np.sum(self.Bound)
    
    def get_NBoundPerSite(self):
        '''
        Returns
        -------
        NBoundPerSite : float
            The average number of transcription factors bound per binding site.

        '''
        return self.get_NBound()/len(self.sites)


    ## Other
    def check_Grid(self):
        '''
        The check_Grid method executes a check that make sure the grid values are 
        valid.

        Raises
        ------
        ValueError
            When the grid values exceed 1 or are negative, a ValueError is raised.

        Returns
        -------
        GO : boolean
            A boolean variable that indicates whether the grid is okay and the 
            simulation can continue.

        '''

        # Check that here is at most 1 TF per site.
        if np.any(self.Grid>1): 
            raise ValueError(f'Grid values are too high: {self.Grid}')
        
        # Check that grid values are positive
        if np.any(self.Grid<0): 
            raise ValueError(f'Grid values are negative: {self.Grid}')
        
        return True