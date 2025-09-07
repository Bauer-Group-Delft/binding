
from mpi4py import MPI

import os
import sys
sys.path.append('../Parameters/')

from Metadata import read_metadata
from Parameters import RunParameters

from Grid1D import TF_Grid_1D
from Grid2D import TF_Grid_2D
from Data import DataFile
from Gillespie import GillespieSim

#%% Stochastic simulation

def simulation(file, htf, eb, J):
    grid_prms, sim_prms, variables = read_metadata(file)
    run_prms = RunParameters(htf, eb, J, variables.IC_prms_list[f'J{J}eb{eb}htf{htf}'], grid_prms)
    
    ## Generate grid
    if grid_prms.dim == 1:
        SimulationGrid = TF_Grid_1D(run_prms, grid_prms)
    elif grid_prms.dim == 2:
        SimulationGrid = TF_Grid_2D(run_prms, grid_prms)
    else:
        raise ValueError('Unsuported grid dimension: {grid_prms.dim}.')
    
    ## Generate datafile
    Datafile = DataFile(SimulationGrid, file+f'_J{J}eb{eb}htf{htf}')
    
    ## Generate simulation
    Simulation = GillespieSim(SimulationGrid, Datafile)
    
    ## Run repeated simulations
    Simulation.repeat_runs(sim_prms.runs, sim_prms.T, sim_prms.T_trans, sim_prms.max_timesteps)
    

def parallel_simulation(file): 
    
    comm = MPI.COMM_WORLD
    nprocs = comm.Get_size()
    myrank = comm.Get_rank()            
    
    print(f'Process {myrank+1}/{nprocs} started.')

    if myrank == 0:
        grid_prms, sim_prms, variables = read_metadata(file)
        print(f'Process {myrank+1}/{nprocs} loaded parameters')
    else:
        grid_prms = None
        sim_prms  = None
        variables = None
        print(f'Process {myrank+1}/{nprocs} ready to recieve parameters.')
     
    grid_prms = comm.bcast(grid_prms, root=0)
    sim_prms  = comm.bcast(sim_prms,  root=0)
    variables = comm.bcast(variables, root=0)    
    
    print(f'Process {myrank+1}/{nprocs} recieved parameters.')
            
    for (J, eb), htf in variables.zipped[myrank::nprocs]:
        print(f'J={J}, eb={eb}, htf={htf}')
        run_prms = RunParameters(htf, eb, J, variables.IC_prms_list[f'J{J}eb{eb}htf{htf}'], grid_prms)
        
        ## Generate grid
        if grid_prms.dim == 1:
            SimulationGrid = TF_Grid_1D(run_prms, grid_prms)
        elif grid_prms.dim == 2:
            SimulationGrid = TF_Grid_2D(run_prms, grid_prms)
        else:
            raise ValueError('Unsuported grid dimension: {grid_prms.dim}.')
        
        ## Generate datafile
        Datafile = DataFile(SimulationGrid, file+f'_J{J}eb{eb}htf{htf}')
        
        ## Generate simulation
        Simulation = GillespieSim(SimulationGrid, Datafile)
        
        ## Run repeated simulations
        Simulation.repeat_runs(sim_prms.runs, sim_prms.T, sim_prms.T_trans, sim_prms.max_timesteps)
        
        ## Attempt to clean memory
        del run_prms
        del SimulationGrid
        del Datafile
        del Simulation
    
    print(f'Process {myrank+1}/{nprocs} finshed.')
    
if __name__ == '__main__':
    '''
    ### Single metadata file
    file_name    = sys.argv[1]
    file         = os.path.join(os.path.realpath('__file__'), f'../../Data/{file_name}') 

    

    parallel_simulation(file)
    
    ## To run in conda terminal use:
    # mpiexec -n 6 python -m mpi4py 'C:/../StochasticModel/StochasticSimulationRun.py' <filename>
    ## Make sure to run the code from the folder where the file is located to prevent errors with loading modules.
    '''

    ### Time sweep metadata files
    file_name = sys.argv[1]
    
    for i in [100]: #[75,85,95]:
        print(f'Simulation {i} starting...')
        file = os.path.join(os.path.realpath('__file__'), f'../../Data/{file_name}_{i}') 
        parallel_simulation(file)