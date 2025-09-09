# %% INITIALIZATION
# Note: Assumes one is running from the ../HillMode/Code/ directory. 
####################################################################################################
### LOAD LIBRARIES
################################################
import numpy as np
import os
import scipy.io
import pandas as pd 
import json 
from scipy import special
from copy import deepcopy

####################################################################################################
### FUNCTIONS
################################################
### EXPERIMENT DATA
def read_bcd_data(rel_path, normalization_mode=1):
    file = os.path.join(os.getcwd(), rel_path)

    if normalization_mode == 1:
        bcd = scipy.io.loadmat(file)['TGB_norm1']
    elif normalization_mode == 2:
        bcd = scipy.io.loadmat(file)['TGB_norm2']
    else:
        raise ValueError(f'Invalid normalization mode: {normalization_mode}.')

    bcd_mean = np.mean(bcd, axis=0)
    bcd_var  = np.var(bcd, axis=0)

    xs_full = np.linspace(0.,1,len(bcd_mean))
    mask    = np.where(np.logical_and(xs_full>=0.1, xs_full<=0.9))

    return bcd, bcd_mean, bcd_var, xs_full, mask

def discretize_concentrations(bcd_data, number_of_samples=60): 
    bcd, _, _, _, mask = bcd_data
    rel_conc_bins = np.arange(np.ceil(np.max(bcd[:,mask])*number_of_samples)+1)/number_of_samples
    return rel_conc_bins


## PROBABILITY DISTRIBUTIONS
def gaussian_conditional_probability_distributions(y_bins, means_and_vars):
    C_bins = len(y_bins)-1
    P_Y_given_X = np.zeros((len(means_and_vars),C_bins))

    for i in range(len(means_and_vars)):
        mu, var = means_and_vars[i]
        # Integrate over a segment dx to get the total probability of being within the segment. 
        # Note that since we renormalize later, we can leave out the factor 1/2.
        segment_integrals = special.erf((y_bins[1:]-mu)/np.sqrt(2*var))-special.erf((y_bins[:-1]-mu)/np.sqrt(2*var)) 
        # Renormalize to make sure the integral remains 1, since we only consider occupancies in the range [0,1], 
        # while Gaussian distributions have as domain the entire real number line.
        P_Y_given_X[i,:] = segment_integrals/np.sum(segment_integrals)
        
    return P_Y_given_X + 1e-50

def get_all_probability_distributions(P_Y_given_X): # P(Y|X)

    # P(Y,X)
    P_X = 1/len(P_Y_given_X)
    P_Y_and_X = P_Y_given_X*P_X

    # P(Y)
    P_Y = np.sum(P_Y_and_X, axis=0) # axis 0 is htf values, axis 1 is probabilities

    return P_Y_given_X, P_Y_and_X, P_Y


### ENTROPY AND MUTUAL INFORMATION
def calculate_entropy(P):
    # (We assume 0*log(0)=0, as is done in 'Element of information theory, second edition' 
    # by Thomas M. Cover and Joy A. Thomas.)
    
    nonzero_terms = np.where(P!=0)
    H = -np.sum(P[nonzero_terms]*np.log2(P[nonzero_terms]))
    return H  

def calculate_conditional_entropy(P_Y_given_X, P_Y_and_X):
    # (We assume 0*log(0)=0, as is done in 'Element of information theory, second edition' 
    # by Thomas M. Cover and Joy A. Thomas.)
    
    nonzero_terms = np.where(P_Y_and_X!=0)
    H = -np.sum(P_Y_and_X[nonzero_terms]*np.log2(P_Y_given_X[nonzero_terms]))
    return H          

def calculate_mutual_information(Ps):
    P_Y_given_X, P_Y_and_X, P_Y = Ps

    # If P_Y_and_X is nonzero, so will P_Y_given_X and P_Y.
    # We should only take the logarithms of those terms, since all others are zero.
    # (We assume 0*log(0)=0, as is done in 'Element of information theory, second edition' 
    # by Thomas M. Cover and Joy A. Thomas.)

    nonzero_terms = np.where(P_Y_and_X!=0)
    I = np.sum(P_Y_and_X[nonzero_terms]*np.log2(P_Y_given_X[nonzero_terms]/P_Y[nonzero_terms[1]]))
    return I

## BINDING SITE DEVICES
def hill_function(S, h, k):
    return S**h / (S**h + k + 1e-50)

def calculate_mean_and_variance(h, k, Ss, tau, mode='koff', koff=1, kon=1, n=1):
    # Calculate mean and variance for the compressed gene expression from our Hill function model.
    # See Appendix A of the paper for the derivation of the variance.
    # 'mode' sets the type of calculation for the Hill model we want to use
    #   - 'koff': sets koff as constant, with a variable kon defined by k
    #   - 'kon': sets kon as constant, with a variable koff defined by k 
    #   - 'non-coop': a check for a non-cooperative binding devices
    #                 here, h sets the number of binding sites

    # Version with koff constant
    if mode == 'koff':
        Cmeans = hill_function(Ss, h, k)
        Cvars = (2./(tau*koff)) * ((Ss**h)/(Ss**h + k)) * (k /(Ss**h +k))**2

    # Version with kon constant
    elif mode == 'kon':
        Cmeans = hill_function(Ss, h, k)
        Cvars = (2./(tau*kon)) * (k*(Ss**h)) / ((k+Ss**h)**3)

    # Version with non-cooperative binding devices (for comparison)
    # Assumes constant koff as well 
    elif mode == 'non-coop':
        Cmeans = hill_function(Ss, n, k)
        Cvars = (2./(tau*h*koff)) * ((Ss**n)/(Ss**n + k)) * (k /(Ss**n +k))**2
        
    else:
        print("No valid mode for calculating the mean and variance of C was chosen.")
    return Cmeans, Cvars

## OTHER
def to_list_if_array(variable):
    if isinstance(variable, np.ndarray):
        return variable.tolist()
    return variable

###################################################################################################
# %% SET PARAMETERS
########################################################################################
# ---------------- GLOSSARY ----------------
# X:    normalized position along embryo (output)
# S:    normalized Bicoid concentration (input)
# C:    normalized hunchback expression (compressed input)

# koff: off rate for Bicoid [s-1]
# kon:  on rate for Bicoid [s-1]
# k:    dissociation equilibrium constant of the binding devices (koff/kon) [-]
# h:    cooperativity parameter / Hill coefficient of the binding sensors [-]
# tau:  τ, averaging (or measurement) time [s]

# S_bins:  number of bins for discretization of the concentrations
# C_bins:  number of bins discretization of the occupation levels
# C_cov:   minimal <C> coverage for the optimal h,k combinations

# hillmode: sets the type of calculation for the Hill model we want to use
#   - 'koff': sets koff as constant, with a variable kon defined by k
#   - 'kon': sets kon as constant, with a variable koff defined by k 
#   - 'non-coop': a check for a non-cooperative binding devices
#                 here, h sets the number of non-cooperative binding sites
#                 rather than the cooperativity. 

# Note: P_Ab denotes marginal distribution of A from P(A,B)
############################################
### Set variable parameters
ks      = np.array(list(np.logspace(-4.5, -1, 60)) +  list(np.linspace(0.11, 0.99, 65)))
hs      = np.linspace(1, 20, 200)
taus    = [10, 600]

hillmode = 'koff'   # see glossary; options: 'koff', 'kon', 'non-coop' 
kon     = 1         # value for kon if kon is set as constant 
koff    = 1         # value for koff if koff is set as constant 
n       = 1         # ONLY FOR 'non-coop' MODE; should be kept at 1 for full non-cooperativity

S_bins  = 60    # discretization of the concentrations
C_bins  = 100   # discretization of the occupation levels
C_cov   = 0.6   # minimal <C> coverage 

############################################
### Create savepath and save parameters
savedir  = f"hillmode-{hillmode}_tau{np.min(taus):.0f}-{np.max(taus):.0f}_h{np.min(hs):.0f}-{np.max(hs):.0f}_k{np.min(ks):.0E}-{np.max(ks):.0E}_kon{kon:.0f}_koff{koff:.0f}_n{n:d}_Sbins{S_bins:d}_Cbins{C_bins:d}_Ccov{C_cov:.2f}"
savepath = f"../Data/ModelData/{savedir}"
if not os.path.exists(savepath):
    os.makedirs(savepath)

parameter_dict = {
    "hillmode" : hillmode,
    "taus" : to_list_if_array(taus),
    "hs" : to_list_if_array(hs),
    "ks" : to_list_if_array(ks),
    "kon" : kon,
    "koff" : koff,
    "n" : n,
    "S_bins" : S_bins,
    "C_bins" : C_bins,
    "C_cov" : C_cov
}
with open(f"{savepath}/parameters.json", "w") as f:
    f.write(json.dumps(parameter_dict, indent=4))

########################################################################################
# %% CALCULATE MUTUAL INFORMATIONS FOR HILL MODEL
########################################################################################
##### 
### Load experiment data
file_path = '../Data/ExperimentData/TG_normBcd.mat'
bcd, bcd_mean, bcd_var, Xs_full, mask = read_bcd_data(file_path, normalization_mode=1)

############################################
### Get marginals P(X), P(S), and conditional P(S|X)
# Extract discretized position data X, P(X) :
Xs  = Xs_full[mask] # Positions from normalized embryo position of 0.1 to 0.9 (160 bins)
P_X = np.ones(len(Xs))*1/len(Xs) # P(X) marginal distribution for positions

# Extract S, discretized concentration data based on bcd data:
# this is between 0 and 1 -> ultimately results in more as S goes to 1.095 due to std fluctuations
rel_conc_bins = discretize_concentrations((bcd, bcd_mean, bcd_var, Xs_full, mask), number_of_samples=S_bins)
Ss = rel_conc_bins[:-1] + 1/S_bins/2 # mid points of the relative concentration bins

# Get P(S|X), P(S,X), P(S)
P_S_given_X = gaussian_conditional_probability_distributions(rel_conc_bins, list(zip(bcd_mean[mask], bcd_var[mask])))
P_S_and_X = P_S_given_X*P_X[:, np.newaxis]
P_Sx = np.sum(P_S_and_X, axis=0)


########################################################################################
##### Calculate I(C;s) and I(C;x) for different h,k combinations for each τ
# Note: entropy is also calculated but not used; can be used to check things. 

# Initialize arrays to store results
I_CX_test = np.zeros((len(hs), len(ks), len(taus)))
I_CS_test = np.zeros((len(hs), len(ks), len(taus)))
H_C_test = np.zeros((len(hs), len(ks), len(taus)))
H_C_given_X_test = np.zeros((len(hs), len(ks), len(taus)))
H_C_given_S_test = np.zeros((len(hs), len(ks), len(taus)))

# Calculate mutual informations  
for tau_idx, tau in enumerate(taus):
    for h_idx, h in enumerate(hs):
        for k_idx, k in enumerate(ks):
            # Calculate mean and variance for binding device from Hill function estimate
            Cmeans, Cvars = calculate_mean_and_variance(h, k, Ss, tau, mode=hillmode, koff=koff, kon=kon, n=n)
            
            # Calculate probability distributions related to C
            P_C_given_S = gaussian_conditional_probability_distributions(np.linspace(0, 1, C_bins + 1), list(zip(Cmeans, Cvars)))
            P_C_given_X = P_S_given_X @ P_C_given_S
            
            Ps_CS = get_all_probability_distributions(P_C_given_S)
            P_C_given_S, P_C_and_S, P_Cs = Ps_CS
            P_C_and_S = P_C_given_S * P_Sx[:, np.newaxis]

            P_Cs = np.sum(P_C_and_S, axis=0)
            Ps_CS = (P_C_given_S, P_C_and_S, P_Cs)
            Ps_CX = get_all_probability_distributions(P_C_given_X)
            P_C_given_X, P_C_and_X, P_Cx = Ps_CX

            # Check whether the normalization is OK
            if not np.allclose(np.sum(P_C_given_S, axis=1), 1) or not np.allclose(np.sum(P_C_given_X, axis=1), 1) or not np.isclose(np.sum(P_C_and_S), 1) or not np.isclose(np.sum(P_C_and_X), 1):
                print('Incorrect normalization, check:')
                print(f'sum(P_C_given_S) = {np.sum(P_C_given_S)}'); print(f'sum(P_C_given_X) = {np.sum(P_C_given_X)}')
                print(f'sum(P_C_and_S) = {np.sum(P_C_and_S)}'); print(f'sum(P_C_and_X) = {np.sum(P_C_and_X)}')

            # Calculate entropies
            H_C = calculate_entropy(Ps_CS[2])
            H_C_given_X = calculate_conditional_entropy(Ps_CX[0], Ps_CX[1])
            H_C_given_S = calculate_conditional_entropy(Ps_CS[0], Ps_CS[1])
            
            # Calculate mutual informations
            I_CX = calculate_mutual_information(Ps_CX)
            I_CS = calculate_mutual_information(Ps_CS)
            
            # Store results 
            H_C_test[h_idx, k_idx, tau_idx] = H_C
            I_CX_test[h_idx, k_idx, tau_idx] = I_CX
            I_CS_test[h_idx, k_idx, tau_idx] = I_CS
            H_C_given_X_test[h_idx, k_idx, tau_idx] = H_C_given_X
            H_C_given_S_test[h_idx, k_idx, tau_idx] = H_C_given_S

            # Print progress
            print(f'tau = {tau:.0f}, h = {h:.2f}, k = {k:.2E}, I(C;s) = {I_CS:.2f}, I(C;x) = {I_CX:.2f}')

# Save the data
result = np.column_stack(
    [np.array(np.meshgrid(hs,ks,taus)).T.reshape(-1, 3), 
     I_CX_test.reshape(-1), 
     I_CS_test.reshape(-1)])
MI_df = pd.DataFrame(result, columns=['h', 'k', 'tau', 'I_CX', 'I_CS'])
MI_df.to_csv(f'{savepath}/MIs_all.csv', index=False)

########################################################################################
# %% CALCULATE OPTIMAL MUTUAL INFORMATION FOR HILL MODEL
########################################################################################
##### 
# Below we want the find the combination of h,k that give the highest I(C;x) for a given I(C;s)
# We do this by picking a range of I(C;s) values, find all h,k combinations that give this value 
# (approximately) from those combinations, we pick the one with the best I(C;x)

# Here we find the indices of the h and k values which give a set of predefined I_CS_values
# Note: since different combinations of h,k can give the same I_CS value we use Nvals 
#       to limit the number of values we store.
Nvals = 100 
I_CS_values = np.arange(0.1, 5.0, 0.02)
I_CS_contours_values = np.full((Nvals, 2, len(I_CS_values), len(taus)), -999.)
for tau_idx, tau in enumerate(taus):
    for I_CS_idx, I_CS_value in enumerate(I_CS_values):
        valid_indices = np.where(np.abs(I_CS_test[:, :, tau_idx] - I_CS_value) < 0.01)
        if valid_indices[0].size > 0:
            maxlen = min(len(valid_indices[0]), Nvals)
            I_CS_contours_values[:maxlen, 0, I_CS_idx, tau_idx] = valid_indices[0][:maxlen]
            I_CS_contours_values[:maxlen, 1, I_CS_idx, tau_idx] = valid_indices[1][:maxlen]


# Now we loop over the valid_indices and find the one with the best I_CX - we also store h and k for this value
I_CX_at_contours = np.zeros((Nvals, len(I_CS_values), len(taus)))
optimal_I_CXs = np.zeros((len(I_CS_values), len(taus)))
optimal_hs = np.zeros((len(I_CS_values), len(taus)))
optimal_ks = np.zeros((len(I_CS_values), len(taus)))

for tau_idx, tau in enumerate(taus):
    for I_CS_idx, I_CS_value in enumerate(I_CS_values):
        maxlen = Nvals - len(np.where(I_CS_contours_values[:, 0, I_CS_idx, tau_idx] < -900)[0])
        if maxlen > 0:
            # Get the valid h and k indices for the current I_CS value
            indices_h = I_CS_contours_values[:maxlen, 0, I_CS_idx, tau_idx]
            indices_k = I_CS_contours_values[:maxlen, 1, I_CS_idx, tau_idx]
            
            # Retrieve I_CX values at these I(C;s) contours
            I_CX_at_contours[:maxlen, I_CS_idx, tau_idx] = I_CX_test[indices_h.astype('int'), indices_k.astype('int'), tau_idx]
            
            # Only consider the I_CX values which have at least a coverage of C_cov
            cur_hs = hs[indices_h.astype('int')]
            cur_ks = ks[indices_k.astype('int')]
            Cmean_coverage = np.zeros((maxlen))
            for idx in range(maxlen):
                Cmean = hill_function(Ss, cur_hs[idx], cur_ks[idx])
                Cmean_coverage[idx] = np.abs((np.max(Cmean) - np.min(Cmean)))
            valid_indices = np.where(Cmean_coverage > C_cov)[0]

            if maxlen != len(valid_indices):
                print(f'Did not find valid indices for tau = {taus[tau_idx]:.0f} at I(C;s) = {I_CS_value:.2f} with C mean coverage over {C_cov:.2f}')

            # Find the maximum I_CX value for the current I_CS value
            if len(valid_indices) > 0:
                valid_I_CX_at_contours = I_CX_at_contours[valid_indices, I_CS_idx, tau_idx]
                best_I_CX = np.max(valid_I_CX_at_contours)
                optimal_I_CXs[I_CS_idx, tau_idx] = best_I_CX

                if best_I_CX > I_CS_value:
                     print(f'For tau = {tau}: I(C|x) > I(C|s) (by more than 10e-3), should not be possible!')
                
                # Find the index of the maximum I_CX value
                best_index = np.where(np.abs(I_CX_at_contours[:maxlen, I_CS_idx, tau_idx] - best_I_CX) < 1e-5)[0][0]
                
                # Store the optimal h and k values
                optimal_hs[I_CS_idx, tau_idx] = hs[int(indices_h[best_index])]
                optimal_ks[I_CS_idx, tau_idx] = ks[int(indices_k[best_index])]

# Save the optimal MI data 
optimal_MI_df = pd.DataFrame(np.column_stack(
    [np.repeat(taus, len(I_CS_values)), 
     optimal_hs.reshape(-1), 
     optimal_ks.reshape(-1), 
     np.tile(I_CS_values, len(taus)), 
     optimal_I_CXs.reshape(-1)]),
     columns = ['tau', 'optimal_h', 'optimal_k', 'I_CS_value', 'optimal_I_CX'])
optimal_MI_df.to_csv(f'{savepath}/MIs_optimal.csv', index=False)

########################################################################################
##### Calculate the optimal I(C;x) and I(C;s) per h 
# Aside from the optimal I(C;x) per I(C;s), we are also interested in the optimal I(C;x) per h
# Calculate the optimal I(C;x) and I(C;s) per h for each tau
k_idx_ICX = np.argmax(I_CX_test, axis=1)
h_indices, tau_indices = np.indices(k_idx_ICX.shape)
optimalh_ks_ICX = ks[k_idx_ICX]
optimalh_I_CXs = I_CX_test[h_indices, k_idx_ICX, tau_indices]

k_idx_ICS = np.argmax(I_CS_test, axis=1)
h_indices, tau_indices = np.indices(k_idx_ICS.shape)
optimalh_ks_ICS = ks[k_idx_ICS]
optimalh_I_CSs = I_CS_test[h_indices, k_idx_ICS, tau_indices]

# Save the data
optimalh_MI_df = pd.DataFrame(np.column_stack(
    [np.repeat(taus, len(hs)), 
     np.tile(hs, len(taus)), 
     optimalh_ks_ICS.reshape(-1),  
     optimalh_I_CSs.reshape(-1),
     optimalh_ks_ICX.reshape(-1),  
     optimalh_I_CXs.reshape(-1)]),     
     columns = ['tau', 'h', 'optimalh_k_ICS', 'optimalh_I_CS', 'optimalh_k_ICX', 'optimalh_I_CX'])
optimalh_MI_df.to_csv(f'{savepath}/MIs_optimal_h.csv', index=False)

########################################################################################
# %% CALCULATE INFORMATION BOTTLENECK FOR GIVEN PARAMETERS
########################################################################################
##### 
# What this does is that it tries to find the optimal binding site for each lambda in
# the optimization goal: I(C,x) - lambda I(C,s). It uses the algorithm described in
# Tishby, N., Pereira, F. C., & Bialek, W. (2000). The information bottleneck method arXiv preprint physics/0004057.

# Set parameters
nruns       = 1 # number of runs, can be used to average for a more exact value
eps         = 2.2204e-16 # small constant to avoid division by zero
lambdas     = 1./(np.array([0.005, 0.01]+list(np.linspace(0.02,1,50)))[::-1]) # lambda_val
niterations = 100 # 

# Initialize storage arrays
P_C_given_S_evolve = {}  # store the evolution of P(C|s)s, the optimal binding site devices
I_CSs_IB = np.zeros((len(lambdas), nruns)) 
I_CXs_IB = np.zeros((len(lambdas), nruns)) 

# Calculate P(X|S) for later use
P_S_and_XT = P_S_and_X.T+eps
P_X_given_S = np.transpose(P_S_and_XT)/(P_Sx[None,:]+eps)

for run_idx in range(nruns):
    # initialize P(C|s) with random values
    P_C_given_S_0 = np.zeros((C_bins,len(Ss)))
    for i in range(len(Ss)): 
        sample = np.random.uniform(low=np.nextafter(0, 1), size=C_bins)
        P_C_given_S_0[:,i]=sample/np.sum(sample)
    
    # Initialize dictionary for each 
    P_C_given_S_evolve[C_bins] = np.zeros((len(lambdas), C_bins, len(Ss)))     
    for lambda_idx, lambda_val in enumerate(lambdas):  
        P_C_given_S = P_C_given_S_0
        P_C_and_S   = P_C_given_S_0*P_Sx 
        P_Sc        = np.sum(P_C_and_S,0) 

        # Iteratively update the distributions  
        for n in range(niterations): 
            # Calculate relevant probability distributions 
            P_C_and_X   = (np.dot(P_C_given_S,P_S_and_XT))
            P_Cs        = np.sum(P_C_and_S,1)
            P_X_given_C = np.transpose(P_C_and_X)/np.sum(P_C_and_X,1)

            # Calculate mutual information I(C;s) for current binding device P_C_given_S
            I_CSs_IB[lambda_idx,run_idx] = np.sum(np.sum(P_C_and_S * np.log2(P_C_and_S/(np.sum(P_C_and_S, 1)[:,np.newaxis] * np.sum(P_C_and_S,0)[np.newaxis,:]+eps)+eps)))
            
            # Calculate mutual information I(C;x)
            I_CXs_IB[lambda_idx,run_idx] = np.sum(np.sum(P_C_and_X * np.log2(P_C_and_X/(np.sum(P_C_and_X, 0)[np.newaxis,:] * np.sum(P_C_and_X,1)[:, np.newaxis]+eps)+eps)))

            # Calculate Kullback-Leibler divergence 
            DKL = np.sum(P_X_given_S[:,np.newaxis,:]*np.log(P_X_given_S[:,np.newaxis,:]/(P_X_given_C[:,:,np.newaxis]+eps)+eps),0)

            # Update distribution
            P_C_given_S = P_Cs[:,np.newaxis]*np.exp(-DKL*lambda_val)
            Zs          = np.sum(P_C_given_S,0) # normalization constant
            P_C_given_S = P_C_given_S/Zs[np.newaxis,:]
            P_C_and_S   = P_C_given_S*P_Sx[np.newaxis,:]

        # Save the P(C|s) distribution for current lambda value 
        P_C_given_S_evolve[C_bins][lambda_idx,:,:] = deepcopy(P_C_given_S)

# Check if the I(C;x) and I(C;s) values make sense
if not np.all(I_CSs_IB > I_CXs_IB):
    raise ValueError("Error in the information bottleneck algorithm, I(C;s) should always be larger than I(C;x).")

# Save the information bottleneck data 
np.savetxt(f'{savepath}/I_CX_IB.csv', I_CXs_IB)
np.savetxt(f'{savepath}/I_CS_IB.csv', I_CSs_IB)