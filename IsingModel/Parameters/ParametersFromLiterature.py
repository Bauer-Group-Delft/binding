# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

#%% Bicoid concentration gradient

'''
The maximum concentrations found in literatrue vary widely (55 nM - 210 nM) [(Gregor, 2007a), 
(Abu-Arish, 2010), (Xu, 2015), (Fernandes, 2022)]. This wide range is likely due to different
methods of measurement and correction. The exact time of measurement might have played a big
role in this, since the concentration has been shown to decrease during interphase of the nuclear 
cycles due to an increase in the size of the nuclei.  

The length constants found in literature range from 0.15 - 0.25 [(Houchmandzadeh, 2002), 
(Gregor, 2007b), (Abu-Arish, 2010), (Liu, 2013), (Durrieu, 2018), (Fernandes, 2022)], 
therefore we will use the interediate value of 0.2 [(Houchmandzadeh, 2002), (Gregor, 2007b)].

(Linda used lambda = 0.244.) 

As a check it might be interesting to repeat the simulations with the parameter values which 
would give the largest concentration range (largest maximum concentration and smallest length
constant) to see how this affects our conculsions.

'''

c0_nM             = 140                            # (Abu-Arish et al., 2010)
min_c0_nM         = 55                             # +/- 3 nM (Gregor et al., 2007a)
max_c0_nM         = 210                            # (Fernandes et al., 2022)

N_Avogadro        = 6.02214076*10**23
c0                = c0_nM*N_Avogadro/(10**24)      # molecules/um^3 =~ 0.6*c0_nM
min_c0            = min_c0_nM*N_Avogadro/(10**24)
max_c0            = max_c0_nM*N_Avogadro/(10**24)

length_const      = 0.2                            # (Houchmandzadeh et al., 2002)
min_length_const  = 0.15                           # (Fernandes et al., 2022)
max_length_const  = 0.25                           # (Abu-Arish et al., 2010)

conc_gradient = lambda x, A, l : A*np.exp(-x/l)

c     = lambda x : conc_gradient(x, c0, length_const)           # in molecules/um^3
c_min = lambda x : conc_gradient(x, min_c0, min_length_const)
c_max = lambda x : conc_gradient(x, max_c0, max_length_const)

small_c_range = conc_gradient(np.array([0,1]), min_c0, max_length_const)
large_c_range = conc_gradient(np.array([0,1]), max_c0, min_length_const)

#%% Hunchback concentration profile

x_rel_threshold   = 0.48                           # (Gregor et al., 2007a)
c_05              = c(x_rel_threshold)
min_c_05          = c_min(x_rel_threshold)
max_c_05          = c_max(x_rel_threshold)
n                 = 5                              # (Gregor et al., 2007a)

Hill_eqn = lambda c, c_05, n : 1/(1+(c_05/c)**n)   # https://en.wikipedia.org/wiki/Hill_equation_(biochemistry)


#%% Length scale

'''
We chose the grid cell size to be 4 nm, which corresponds to approximately 12 bp. This is the 
smallest distance between the centers of two Bicoid binding sites in the hunchback P2 promoter
[(Fernandes, 2022) supplementary file 1].
'''

grid_cell_size_nm = 4                              # nm, 4 nm is around 12 bp
grid_cell_size    = grid_cell_size_nm*10**-3       # um
V_site            = (grid_cell_size)**3            # um^3


#%% htf values

beta       = 1
c_to_htf   = lambda beta, c   : 1/beta*np.log(c/(1-c))            # Note that c must be in molecules/grid cell
htf_to_c   = lambda beta, htf : 1/(1+np.exp(-beta*htf))           # Note that c will be in molecules/grid cell
# c_to_htf = lambda beta, c   : 1/beta*(np.log(c)-np.log(1-c))
# htf_to_c = lambda beta, htf : np.exp(beta*htf)/(np.exp(beta*htf)+1)

htf_range     = c_to_htf(beta, c(np.array([0,1]))*V_site)
min_htf_range = c_to_htf(beta, c_min(np.array([0,1]))*V_site)
max_htf_range = c_to_htf(beta, c_max(np.array([0,1]))*V_site)
htf_05        = c_to_htf(beta, c_05*V_site)
min_htf_05    = c_to_htf(beta, min_c_05*V_site)
max_htf_05    = c_to_htf(beta, max_c_05*V_site)

small_htf_range = c_to_htf(beta, small_c_range*V_site)
large_htf_range = c_to_htf(beta, large_c_range*V_site)


#%% Volume (only relevant when including depletion effects)

'''
The Diameter of the nuclei varies between approximately 5 and 10 um over nc10-14. 
In nuclear cycle 14 it is approximately 6 um [(Gregor, 2007b), figure 4c].
'''

R               = 3                            # radius nucleus in um (interphase nc 14, (Gregor et al., 2007b))
V_nucleus       = 4/3*np.pi*R**3               # um^3
N_sites_tot     = V_nucleus/V_site


#%% Time scale

'''
(Gregor, 2007a) found a lower bound for the diffusion constant of ~0.3 um^2/s 
[(Abu-Arish, 2010)]. The actual diffusion constant is likely closer to ~7 um^2/s
[(Abu-Arish, 2010), (Porcher, 2010)], although it has also been suggested that there
is a slow and a fast moving population of Bicoid, with diffusion constants 0.22 um^2/s
and 7.7 um^2/s, respectively [(Porcher, 2010)]. Since in our model we consider free 
moving Bicoid, which can bind to DNA, we will use a diffusion constant of 7 um^2/s.

The nuclear cycles 10-13 take about 10-20 minutes [(Porcher, 2010)], meaning that the
concentration measurements must occur at a timescale of minutes. Therefore, we will
use a simulation time of 10 minutes.
'''

dims            = 1                            # dimensions of the grid
D               = 7                            # um^2/s Diffusion constant of Bicoid (Abu-Arish et al., 2010)
runtime         = 600                          # seconds = 10 min (as used by Linda)

# MSD = 2mDt, where m is the number of dimensions in which diffusion occurs
time_conversion = lambda d, dx : dx**2/(2*(3-d)*D) 


#%%

def get_htf_from_number_of_molecules(N_molecules):
    c = np.array(N_molecules)/V_nucleus  # Concentration in molecules/um^3
    c_nM = c/(N_Avogadro/(10**24))       # Concentration in nM
    l = -1/np.log(c[1]/c[0])
    htf_range = c_to_htf(beta, c*V_site)
    print(f'For a gradient ranging from {N_molecules[0]} to {N_molecules[1]} molecules,\nthe initial concentration is {c_nM[0]} nM and the decay length is {l}.')
    print(f'The corresponding htf range is {htf_range}, with htf_05 equal to {c_to_htf(beta, conc_gradient(x_rel_threshold, c[0], l)*V_site)}.')

def plot_concentration_profiles():
    positions  = np.linspace(0,1,101)
    conc_range = c(positions)

    plt.figure(figsize=(6,4))
    plt.plot(positions, conc_range/c0, c='darkgreen', label='Bcd', linewidth=3)
    plt.fill_between(positions, c_min(positions)/min_c0, c_max(positions)/max_c0, color='darkgreen', alpha=0.3)
    plt.plot(positions, Hill_eqn(conc_range, c_05, n), c='deeppink', label='Hb', linewidth=3)
    plt.fill_between(positions, Hill_eqn(c_min(positions), min_c_05, n), Hill_eqn(c_max(positions), max_c_05, n), color='deeppink', alpha=0.3)

    plt.xlabel('relative position, x/L')
    plt.ylabel('relative concentration')
    plt.xlim([positions[0],positions[-1]])
    plt.xlim([0,1])
    plt.xticks(np.arange(6)*0.2)
    plt.ylim([0,1])
    plt.yticks(np.arange(6)*0.2)
    plt.legend()

    print(f'Minimum Bicoid concentration: ~{int(c(1)/0.6)} nM')
    
def plot_number_of_molecules():
    x_rel              = np.linspace(0,1,101)     # relative position along the embryo
    molecules_per_nucl = c(x_rel)*V_nucleus
    molecules_per_site = c(x_rel)*V_site 
    htf_along_embryo   = c_to_htf(beta, molecules_per_site)

    fig, ax = plt.subplots(1,3, figsize=(10,3))
    fig.tight_layout(pad=3.0)
    
    ax[0].plot(x_rel, molecules_per_nucl, c='darkgreen', linewidth=3)
    ax[0].fill_between(x_rel, c_min(x_rel)*V_nucleus, c_max(x_rel)*V_nucleus, color='darkgreen', alpha=0.3)
    ax[0].set_xlabel('relative position, x')
    ax[0].set_ylabel('Number of molecules per nucleus, N')
    ax[0].set_xlim([0,1])
    #ax[0].set_title('Number of molecules per nucleus')
    
    ax[1].plot(x_rel, molecules_per_site, c='darkgreen', linewidth=3)
    ax[1].fill_between(x_rel, c_min(x_rel)*V_site, c_max(x_rel)*V_site, color='darkgreen', alpha=0.3)
    ax[1].set_xlabel('relative position, x')
    ax[1].set_ylabel('Number of molecules per site, N')
    ax[1].set_xlim([0,1])
    #ax[1].set_title('Number of molecules per site (filling fraction)')
    
    ax[2].plot(x_rel, htf_along_embryo, c='darkgreen', linewidth=3, label='$h_{tf}$')
    ax[2].plot(x_rel, c_to_htf(beta, conc_gradient(x_rel, min_c0, max_length_const)*V_site), c='darkgreen', linewidth=1, label='small $h_{tf}$ range')
    ax[2].plot(x_rel, c_to_htf(beta, conc_gradient(x_rel, max_c0, min_length_const)*V_site), c='darkgreen', linewidth=1, label='large $h_{tf}$ range')
    ax[2].fill_between(x_rel, c_to_htf(beta, c_min(x_rel)*V_site), c_to_htf(beta, c_max(x_rel)*V_site), color='darkgreen', alpha=0.3)
    ax[2].scatter(x_rel_threshold, htf_05, c='darkgreen', label='switching point')
    ax[2].scatter(x_rel_threshold, min_htf_05, c='darkgreen', alpha=0.3)
    ax[2].scatter(x_rel_threshold, max_htf_05, c='darkgreen', alpha=0.3)
    ax[2].set_xlabel('relative position, x')
    ax[2].set_ylabel('$h_{tf}$')
    ax[2].legend()
    ax[2].set_xlim([0,1])
    #ax[2].set_title('$h_{tf}$')

    print(f'Total number of sites in a nucleus: {int(N_sites_tot)}')
    print(f'Number of molecules per nucleus: {int(c_min(0)*V_nucleus)}-{int(c_max(0)*V_nucleus)} ({int(c(0)*V_nucleus)}) (anterior) to {int(c_min(1)*V_nucleus)}-{int(c_max(1)*V_nucleus)} ({int(c(1)*V_nucleus)}) (posterior).')
    print(f'Relevant htf range (min): {min_htf_range[0]} to {min_htf_range[1]} ({[int(c_min(0)*V_nucleus), int(c_min(1)*V_nucleus)]}).')
    print(f'Relevant htf range (max): {max_htf_range[0]} to {max_htf_range[1]} ({[int(c_max(0)*V_nucleus), int(c_max(1)*V_nucleus)]}).')
    print(f'Relevant htf range (smallest): {small_htf_range[0]} to {small_htf_range[1]} ({(small_c_range*V_nucleus).astype(int)}).')
    print(f'Relevant htf range (largest): {large_htf_range[0]} to {large_htf_range[1]} ({(large_c_range*V_nucleus).astype(int)}).')
    print(f'Relevant htf range (intermediate): {htf_range[0]} to {htf_range[1]}.')
    print(f'The boundary of Hunchback expression occurs around relative position {x_rel_threshold}.\nAt this position, the htf is {min_htf_05} to {max_htf_05} ({htf_05}).')
    
def print_timescale_info():
    print(f'For a {dims}D grid with a grid cell size of {grid_cell_size} nm:')
    print(f'A time unit in the model corresponds to {time_conversion(dims, grid_cell_size)} seconds.')
    print(f'A runtime of {runtime} seconds corresponds to {runtime/time_conversion(dims, grid_cell_size)} time units.\n')