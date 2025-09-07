# -*- coding: utf-8 -*-

import random
import numpy as np
from IPython.display import clear_output

'''
The gillespie simulation should be initialized with a grid object representing 
the grid to be simulated and an optional datafile object used to store data. The
simulation is initialized with time = 0. If a datafile is available, the 
bath information and initial state are stored.

The run method is used to run the simulation until a set end time or maximum 
number of timesteps is reached. Every iteration the simulation obtains the 
relevant rates from the grid, calculates the timestep using the update_time 
method, chooses and perfomes a reaction using the perform_reaction method and, 
if a datafile is available, stores the new data.

The simulation can be reset using the reinitialize method.
'''

class GillespieSim:
    
    def __init__(self, grid, datafile=None):
        '''
        The init method is called when a GillespieSim object is created.
        
        Parameters
        ----------
        grid : TF_Grid object
            A grid object, which represents the simulation grid.
            
            The grid obect must have the following methods: 
                setup, get_AppearRates, get_DisappearRates, update, check_Grid.

        datafile : DataFile object, optional
            A datafile object used to store the results of the simulation.
            
            The datafile object, if supplied, must have the following methods:
                set_grid, store_info, store_data.
                
            The default is None.
        
        Returns
        -------
        None.
        
        '''
        
        self.time = 0
        self.grid = grid
        self.datafile = datafile            
    
    def run(self, T, max_timesteps=int(1e6)):
        '''
        The run method runs a simulation of the simulation grid for a set time T.
        Note: How much the time advances in each step can be different. The
        simulation ends when the number of timesteps required to reach T is larger
        than max_timesteps.
        
        T : float
            The time to be simulated.
        max_timesteps : integer, optional
            The maximum number of timesteps for which the simulation should run. 
            The default is int(1e6).
        
        Returns
        -------
        None.
        '''
        
        if self.datafile:
            # Store initial state.
            self.datafile.store_data(self.time)
        
        for i in range(max_timesteps):
            #print(i) 
            
            # Get updated rates
            AppearRates    = self.grid.get_AppearRates()
            DisappearRates = self.grid.get_DisappearRates()
            
            # Calculate time of next event
            self.time = self.update_time(self.time, np.sum(DisappearRates) + np.sum(AppearRates))
            
            '''
            # TODO: make code able to output at multiple times
            # Check whether it is time to output results
            if self.time > T_output:
                # Write results to file
                if self.datafile: # TODO: remove this if statement from if self.time > T.
                    self.datafile.store_data(T_output)
                    self.datafile.write_run_data(run, transient, warning)
                    # TODO: make sure write_run_data also writes the time and this info can be read in correctly
                    # If we output here, we can remove the output from repeat_runs.
                i_out += 1 # TODO: initialize i_out to 0
                T_output = output_times[i_out] # next output time
                # TODO: make sure output times is given as argument
                
            '''
            
            # Check whether next event occurs within the simulated time.
            if self.time > T:
                if self.datafile:
                    self.datafile.store_data(T)
                return
            
            # Perform one step of the Gillespie algorithm
            self.perform_reaction(AppearRates, DisappearRates)

            # Check whether everything is okay.
            Grid_OK = self.grid.check_Grid()
            if not Grid_OK:
                # Something is wrong with the grid or a stopping condition is reached.
                return
            
            # Store data
            if self.datafile:
                self.datafile.store_data(self.time)
        
        warning = f'Warning: Maximum number of timesteps ({max_timesteps}) reached at time {self.time}/{T}.'
        print(warning)
        return warning

    # Gillespie time update
    def update_time(self, time, rate):
        
        '''
        The update_time method updates the simulation time according to the
        Gillespie algorithm:
            t = -(1/rate)*ln(randomNumber),
        where rate is the total reaction rate and randomNumber is a random
        number uniformly distributed between 0 and 1.

        Parameters
        ----------
        time : float
            The current simulation time.
            
        rate : float
            The sum of the rates of all possible events.

        Returns
        -------
        newTime : float
            The updated simulation time.
            
        '''
        
        randomNumber = random.random()
        if randomNumber == 0: 
            # Take logarithm of a very small number, since log(0)=-infty.
            # This essentially sets an upper limit to the time step size.
            # TODO: it might be better to use a random number 0 < x <= 1.
            # TODO: then also change random number generation below and use <= for comparison again.
            minRandomNumber = 1e-100
            dt = - np.log(minRandomNumber)/rate
            print('Warning: Random number equals 0:', randomNumber)
        else: 
            dt = - np.log(randomNumber)/rate
        newTime = time + dt
        return newTime
    
    
    # This is the core Gillespie function: given the number of available processes and rates, it picks an event
    def perform_reaction(self, AppearRates, DisappearRates):
        '''
        The perform_reaction method determines which event occurs according to 
        the gillespie algorithm:
            
        The grid is then updated according to the chosen event.

        Parameters
        ----------
        AppearRates : numpy array
            An array with the same shape as the grid, containing for every position
            in the grid the appearance rate of a TF.
            
        DisappearRates : numpy array
            An array with the same shape as the grid, containing for every position
            in the grid the disappearance rate of a TF.

        Raises
        ------
        ValueError
            When the random number is larger than sum of all rates, a ValueError is raised.

        Returns
        -------
        None.

        '''
        
        sumAppearRates    = np.sum(AppearRates)
        sumDisappearRates = np.sum(DisappearRates)
        
        # Pick a random number between zero and the sum of all rates.
        # The value of this random number determines the event that has occured.
        randomNumber = random.random()*(sumAppearRates + sumDisappearRates)
        # Note: rand.random() generates a uniformly distributed radom number 0 <= x < 1.
        # So to make the comparisons fair, we need to use x < rate instead of x <= rate.
        
        # Either a particle appeared ...
        if randomNumber < sumAppearRates:
            rates = AppearRates
            increment = +1
        # ... or it disappeared.
        elif randomNumber <  sumAppearRates + sumDisappearRates:
            randomNumber -= sumAppearRates
            rates = DisappearRates
            increment = -1
        else:
            raise ValueError("Random number is larger than sum of all rates: %3.2f %3.2f" \
                   %(randomNumber, sumAppearRates + sumDisappearRates))

        # Update a randomly chosen position in the grid.
        self.grid.update(randomNumber, rates, increment)
        
        
    def repeat_runs(self, runs, T, transient=None, max_timesteps=int(1e6)):
        '''
        Parameters
        ----------
        runs : integer
            The nuber of times to repeat the simulation.
        T : float
            The time to be simulated.
        max_timesteps : integer, optional
            The maximum number of timesteps for which the simulation should run. 
            The default is int(1e6).

        Returns
        -------
        None.

        '''
        
        for run in range(runs):
            print('run',run)
            
            # Reset simulation
            self.reinitialize()
            
            # Run simulation
            warning = self.run(T, max_timesteps)
            
            # Write results to file
            if self.datafile:
                self.datafile.write_run_data(run, transient, warning)
            
            clear_output(wait=True)
        
        
    def reinitialize(self):
        '''
        This function resets the simulation by using new initial conditions and
        deleting the stored data.

        Returns
        -------
        None.

        '''
        
        self.time = 0
        self.grid.setup()
        if self.datafile:
            self.datafile.setup()
