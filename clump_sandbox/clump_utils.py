import sys
sys.dont_write_bytecode = True
import os
import numpy as np
from scipy.stats import poisson, geom
from scipy.optimize import minimize
import random
#====================================================================
def evaluate_clump_distribution(PARAM_DICT, genomes_well, x=None, mean=None):
    if (x == None):
        min_clump_size, max_clump_size = get_min_max_clump_size(PARAM_DICT, genomes_well, mean)
        x = np.arange(min_clump_size, max_clump_size+1)

    if (PARAM_DICT['distribution'] == 'geometric'):
        if (mean > .99999):
            print(" >> Error: mean of geometric distribution must be between 0 and 1. Got:" +str(mean))
            mean = .99999
        return np.asarray(geom.pmf(x, mean)).ravel()
    
    elif (PARAM_DICT['distribution'] == '1inf-geo'):
        if (mean > .99999):
            print(" >> Error: mean of geometric distribution must be between 0 and 1. Got:" +str(mean))
            mean = .99999
        
        y = geom.pmf(x, mean)
        f_1 = np.float32(PARAM_DICT['f_1'])
        y_rest = y[1:]
        y_rest_scaled = y_rest * ((1-f_1) / np.sum(y_rest))
        return np.concatenate(([f_1], y_rest_scaled)).ravel()
    
    elif (PARAM_DICT['distribution'] == 'poisson'):
        x = np.arange(0, max_clump_size) # Poisson needs to start at 0
        return np.asarray(poisson.pmf(x, mean)).ravel()

#====================================================================
def FitClumpsToDistribution(genomes_well,
                            PARAM_DICT,
                            mean):
    '''
    Parameters:
     - genomes_well <int>: How many virions there are
     - PARAM_DICT <dict>: The distribution of clump sizes
     - mean <float, None>: the parameter for the distribution
    
    Returns
     - NUM_CLUMPS <ndarray>: The fitted clump size frequencies. MAY CONTAIN FLOATS
     - best_nll <float>: The lowest negative log likelihood value obtained during fitting
    '''
    #----------------------------------------------------------------
    ''' Get max/min clump size, 1st/99th percentile '''

    min_clump_size, max_clump_size = get_min_max_clump_size(PARAM_DICT, genomes_well, mean)
    #----------------------------------------------------------------
    ''' Arrays to hold clump stuff '''
    
    NUM_CLUMPS = np.empty((int(PARAM_DICT['num_fits']), max_clump_size), float)
    NLL        = np.empty(int(PARAM_DICT['num_fits']), float)
    #----------------------------------------------------------------
    ''' Set up constraint '''

    def constraint(NUM_CLUMPS):
        total = 0
        for i in range(len(NUM_CLUMPS)):
            total += (i+1) * NUM_CLUMPS[i]
        return total - genomes_well # Must equal 0

    constraint_dict = {'type':'eq', 'fun': constraint}
    #----------------------------------------------------------------
    ''' Get y-values to fit to: normalized outputs of clump distribution '''
    
    # if (distribution_name == 'poisson'):
    #     # CLUMP_SIZES = np.arange(0, max_clump_size) # poisson has to start at 0
        
        
    #     # RATIOS = NormalizeDiscrete(clump_distribution.pmf(CLUMP_SIZES))
    #     CLUMP_RATIOS      = evaluate_clump_distribution(PARAM_DICT)
    #     NORMALIZED_RATIOS = NormalizeDiscrete(CLUMP_RATIOS) 
    #     CLUMP_SIZES = np.arange(1, max_clump_size+1)
    
    # elif (distribution_name == 'geometric'):
    #     CLUMP_SIZES = np.arange(1, max_clump_size+1) # geometric has to start at 1
    #     RATIOS = NormalizeDiscrete(clump_distribution.pmf(CLUMP_SIZES))

    CLUMP_RATIOS      = evaluate_clump_distribution(PARAM_DICT, genomes_well, mean=mean)
    NORMALIZED_RATIOS = NormalizeDiscrete(CLUMP_RATIOS)

    print("~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~")
    print("RATIOS =", NORMALIZED_RATIOS, np.shape(NORMALIZED_RATIOS))
    print("~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~")
    #----------------------------------------------------------------
    ''' Get the bounds: try to scale upper bound according to ratio - add wiggle room'''
    
    BOUNDS = []
    for clump_size in range(max_clump_size):
        if (clump_size < min_clump_size-1): # clump size is too small
            BOUNDS.append([0, .99])
       
        elif (clump_size == min_clump_size): # clump size is the minimum possible. Usually this is 1
            BOUNDS.append([0, genomes_well / (clump_size+1) * 10 * NORMALIZED_RATIOS[clump_size]])
        
        elif (clump_size == max_clump_size - 1): # clump size is the maximum possible
            BOUNDS.append([.01, genomes_well])
        
        else: # clump size in between minimum and maximum. Scale upper bound according to the distribution. x2 gives wiggle room
            bound = genomes_well / (clump_size+1) * 10 * NORMALIZED_RATIOS[clump_size]
            if (bound < 2):
                bound = 2
            BOUNDS.append([0, bound])

    # print("~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~")
    # print("BOUNDS =", BOUNDS, np.shape(BOUNDS))
    # print("~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~")
    #----------------------------------------------------------------
    ''' Repeatedly try fitting with perturbed x-values '''
    
    for trial in range(int(PARAM_DICT['num_fits'])):
        # Fit clump frequencies to the clump distribution
        initial_guess = []
        for j in range(NORMALIZED_RATIOS.shape[0]):
            val = genomes_well / (j+1) * NORMALIZED_RATIOS[j] + random.randint(-50, 50)
            if (val > 0):
                initial_guess.append(val)
            else:
                initial_guess.append(genomes_well / (j+1) * NORMALIZED_RATIOS[j])
        
        # print("~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~")
        # print("initial_guess =", initial_guess, np.shape(initial_guess))
        # print("~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~")

        result = minimize(NegLogLikeClump,
                          initial_guess,
                          constraints=constraint_dict,
                          bounds=BOUNDS,
                          args=(NORMALIZED_RATIOS),)

        NUM_CLUMPS[trial, :] = result.x

        # Evaluate negative log likelihood after fitting
        total = 0
        for j in range(len(NUM_CLUMPS[trial, :])):
            total += (j+1) * NUM_CLUMPS[trial, j]
        
        NLL[trial] = NegLogLikeClump(NUM_CLUMPS[trial, :], NORMALIZED_RATIOS)

        # print(" >> sum=", total, ", gen/well", genomes_well, ", nll=", NegLogLikeClump(NUM_CLUMPS[trial, :], RATIOS), ", NUM_CLUMPS", NUM_CLUMPS.shape)
    
    min_idx = np.where(NLL == np.min(NLL))
    
    # print("~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~")
    # print("FITTED NUM_CLUMPS =", NUM_CLUMPS[min_idx, :][0][0], np.shape(NUM_CLUMPS[min_idx, :][0][0]))
    # print("~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~")
    
    return NUM_CLUMPS[min_idx, :][0][0], NLL[min_idx]
#====================================================================
def CreateClumps(genomes_well, PARAM_DICT):
    '''
    Parameters:
     - genomes_well <int>: How many virions there are
     - PARAM_DICT <dict>: Dictionary of parameters, obtained from sandbox_config.yml
    
    Returns
     - CLUMP_IDS <list>: List of identifiers for each virion which determines which clump they're in
     - NUM_CLUMPS <list>: The frequencies of clumps of size 1+index in list
     - CLUMP_SIZES <list>: 1, 2, 3, ..., max_clump_size
     - max_clump_size <int>: The maximum number of virions in a clump
    '''
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ''' Prepare the parameters for the clump distribution '''

    mean = None
    if (PARAM_DICT['fixed_mean'] == True):
        mean = float(PARAM_DICT['mean'])
        #clump_size_distribution = get_clump_distribution(mean, PARAM_DICT['distribution'])
    
    elif (PARAM_DICT['distribution'] in ['geometric', '1inf-geo', 'poisson']):
        if (PARAM_DICT['distribution'] in ['geometric', '1inf-geo']):
            mean = float(PARAM_DICT['vMax_c']) / genomes_well
        
        elif (PARAM_DICT['distribution'] == 'poisson'):
            mean = genomes_well / float(PARAM_DICT['vMax_c'])
        
        #clump_size_distribution = get_clump_distribution(mean, PARAM_DICT['distribution'])
    
    print("GENOMES/WELL =", genomes_well, ", dist =", PARAM_DICT['distribution'], ", mean =", round(mean, 5), ", vMax =", PARAM_DICT['vMax'])
    print(" ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~")
    #----------------------------------------------------------------
    if (mean <= 10**-6): # Mean is too small, return only single virions
        print("AAAAAAHHHHH, mean=", mean)
        print(" ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~")
        print("NUM_CLUMPS =", np.array([genomes_well]), ", len =", 1, ", genomes_well =", genomes_well, " | len(CLUMP_IDS) =", np.arange(0, genomes_well, 1).shape[0])
        print(" ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~")
        return np.arange(0, genomes_well, 1), np.array([genomes_well]), np.array([1]), 1
    
    
    else:
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ''' Fit '''
        
        NUM_CLUMPS, best_nll = FitClumpsToDistribution(genomes_well,
                                                       PARAM_DICT,
                                                       mean) # Not rounded yet
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ''' Turn to integers '''

        total = 0
        for i in range(NUM_CLUMPS.shape[0]):
            total         += (i+1) * int(NUM_CLUMPS[i])
            NUM_CLUMPS[i]  = int(NUM_CLUMPS[i])
        total2 = 0
        if not (total == genomes_well):
            for i in range(1, len(NUM_CLUMPS)):
                total2 += (i+1) * NUM_CLUMPS[i]
            NUM_CLUMPS[0] = genomes_well - total2
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ''' Create clump ids '''

        max_clump_size = NUM_CLUMPS.shape[0]+1
        CLUMP_SIZES = np.arange(1, max_clump_size)
        CLUMP_IDS = [] # Integer assigned to each viron denoting which clump it is a part of

        clump_num = 0 # Assign each virion an integer 'ID', virions in the same clump will have the same ID
        index = 0
        for i in range(len(NUM_CLUMPS)): 
            for j in range(int(NUM_CLUMPS[i])):
                for k in range(CLUMP_SIZES[i]):
                    #print("clump size ", CLUMP_SIZES[i]+1, ", clump number = ", j, ", clump_num = ", clump_num, ", index = ", index)
                    index = index + 1
                    CLUMP_IDS.append(clump_num)
                clump_num = clump_num + 1
        print(" ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~")
        print("NUM_CLUMPS =", NUM_CLUMPS, ", len =", len(NUM_CLUMPS), ", genomes_well =", genomes_well, " | len(CLUMP_IDS) =", len(CLUMP_IDS))
        print(" ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~")
        return CLUMP_IDS, NUM_CLUMPS, CLUMP_SIZES, max_clump_size
#====================================================================
def NegLogLikeClump(NUM_CLUMPS, RATIOS):
    '''
    Parameters
     - NUM_CLUMPS <ndarray>: The frequencies of clump sizes
     - RATIOS <ndarray>: The normalized desired clump size frequencies obtained from "clump distribution"
    '''
    nll = 0
    epsilon = 10**-10
      
    for i in range(np.shape(NUM_CLUMPS)[0]): # Normalize
        model_result = NUM_CLUMPS[i] / sum(NUM_CLUMPS) # Estimated from clump distribution
        
        if (model_result < epsilon):
            model_result = epsilon
        if (model_result > 1 - epsilon):
            model_result = 1 - epsilon
            
        #print("num_clumps[i] * (i+1) / num_objects: ", num_clumps[i], i+1, num_objects)
        nll += RATIOS[i] * np.log(model_result) + (1 - RATIOS[i]) * np.log(1 - model_result)
    
    return -nll
#====================================================================
def NormalizeDiscrete(y):
    y = np.asarray(y)

    RATIOS = np.ones_like(y) * -999
    for i in range(y.shape[0]):
        RATIOS[i] = y[i] / np.sum(y)
    return RATIOS
#====================================================================
def get_min_max_clump_size(PARAM_DICT, genomes_well, mean=None):
    min_clump_size, max_clump_size = 1, 1
    
    if (PARAM_DICT['distribution'] == 'custom'):
        custom_values = np.asarray([np.float(x) for x in PARAM_DICT['custom_values']])
        max_clump_size = custom_values.shape[0]+1 # count how many elements are in the list of 'custom_values'
    
    else:
        if (PARAM_DICT['distribution'] in ['geometric', '1inf-geo', 'poisson']):
            if (PARAM_DICT['distribution'] in ['geometric', '1inf-geo']):
                f = geom(mean)
            elif (PARAM_DICT['distribution'] == 'poisson'):
                f = poisson(mean)
            
            while(True):
                if (f.cdf(min_clump_size) <= .01):
                    min_clump_size +=1
                if (f.cdf(max_clump_size) >= .99):
                    break
                else:
                    max_clump_size += 1

    print(" >> found max clump size =", max_clump_size, ", min clump size =", min_clump_size)

    if (max_clump_size > genomes_well):
        print(" >> >> Error. max clump size", max_clump_size, "> genomes/well", genomes_well, ". Set max clump size =", genomes_well)
        max_clump_size = genomes_well
    
    return min_clump_size, max_clump_size