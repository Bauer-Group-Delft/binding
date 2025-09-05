# -*- coding: utf-8 -*-

#%% Bisection

def bisection(var_init, func, func_args=None, tol=0.02):
 
    occ0 = func(var_init[0], func_args)
    occ1 = func(var_init[1], func_args)
    if ((occ0 >= 0.5 - tol/2) and (occ0 <= 0.5 + tol/2)):
        # If lower bound is close enough to 0.5 occupancy, retrun lower bound
        return var_init[0]
    elif ((occ1 >= 0.5 - tol/2) and (occ1 <= 0.5 + tol/2)):
        # If upper bound is close enough to 0.5 occupancy, retrun upper bound
        return var_init[1]
    else:
        # Check that range is valid for doing bisection
        if not (occ0<=0.5 and occ1>=0.5):
            raise ValueError(f'Invalid variable range: var_init = {var_init}, occupancies = [{occ0},{occ1}]')
    # Note that these checks make sure that the difference in occupancy to the halfway point is less than half 
    # the tolerance, while the bisection method below makes sure the distance between the obtained halfwaypoint 
    # and the actual halfway point is less than half the tolerance. As long as the slope around the halfwaypoint
    # is larger or equal to 1, this should not lead to problems, because the curves are increasing.
    
    # Bisection method
    while (var_init[1]-var_init[0]) > tol:

        var = (var_init[0]+var_init[1])/2            
        occupancy = func(var, func_args)

        if occupancy > 0.5:
            var_init[1] = var
        else:
            var_init[0] = var

    return (var_init[0]+var_init[1])/2 # Error = tol/2