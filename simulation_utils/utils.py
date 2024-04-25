import numpy as np
from scipy.stats import skewnorm
import yaml
#====================================================================
''' Parameter list (Initialized with already predicted values) '''
def PrepareParameters(config, simul_name, sheet, scale=None):
    assert simul_name in ['clump', 'comp', 'acc_dam', 'clump_comp', 'clump_acc_dam', 'null'], simul_name+" must be 'clump', 'comp', 'acc_dam', 'clump_comp', 'clump_acc_dam', or 'null'."
    #----------------------------------------------------------------
    PARAM_DICT = {}
    #----------------------------------------------------------------
    PARAM_DICT['sheet']           = config['SIMULATION_PARAMETERS']['sheet']
    PARAM_DICT['scale']           = config['SIMULATION_PARAMETERS']['scale']
    PARAM_DICT['num_simulations'] = config['SIMULATION_PARAMETERS']['num_simulations']
    PARAM_DICT['i_mean']  = config['GFP_VIRUS_PARAMETERS']['mu_inf'] #0
    PARAM_DICT['i_stdev'] = config['GFP_VIRUS_PARAMETERS']['std_inf'] #1
    PARAM_DICT['r_mean']  = config['CELL_PARAMETERS']['mu_res']
    PARAM_DICT['r_stdev'] = config['CELL_PARAMETERS']['std_res']
    PARAM_DICT['gamma']   = config['SIMULATION_PARAMETERS']['gamma']

    if (scale == None):
        PARAM_DICT['scale'] = config['SIMULATION_PARAMETERS']['scale'] #100
    else:
        PARAM_DICT['scale'] = scale
    #----------------------------------------------------------------
    if ('clump' in simul_name):
        PARAM_DICT['mean'] = config['CLUMP_PARAMETERS']['mean']
        PARAM_DICT['lb'] = config['CLUMP_PARAMETERS']['lb']
        PARAM_DICT['ub'] = config['CLUMP_PARAMETERS']['ub']
        PARAM_DICT['scheme'] = config['CLUMP_PARAMETERS']['scheme']
        PARAM_DICT['distribution'] = config['CLUMP_PARAMETERS']['distribution']
        vMAX = config['CLUMP_PARAMETERS']['vMAX']
        PARAM_DICT['vMax'] = float(vMAX[sheet]) / PARAM_DICT['scale']

    if ('comp' in simul_name):
        KAPPA = config['COMPENSATION_PARAMETERS']['kappa']
        vMAX  = config['COMPENSATION_PARAMETERS']['vMax']
        PARAM_DICT['kappa'] = float(KAPPA[sheet])
        if not ('clump' in simul_name):
            PARAM_DICT['vMax']  = float(vMAX[sheet]) / PARAM_DICT['scale']

    if ('acc_dam' in simul_name):
        BETA = config['ACCRUED_DAMAGE_PARAMETERS']['beta']
        vMAX = config['ACCRUED_DAMAGE_PARAMETERS']['vMax']
        PARAM_DICT['beta'] = float(BETA[sheet])
        if not ('clump' in simul_name):
            PARAM_DICT['vMax'] = float(vMAX[sheet]) / PARAM_DICT['scale']
    
    if (simul_name == 'null'):
        B     = config['NULL_PARAMETERS']['b']
        vMAX  = config['NULL_PARAMETERS']['vMax']
        PARAM_DICT['b']    = float(B[sheet])
        PARAM_DICT['vMax'] = float(vMAX[sheet]) / PARAM_DICT['scale']
    #----------------------------------------------------------------
    return PARAM_DICT
#====================================================================
def PrepareData(dose_inf_df, scale=None):
    num_zeros = 0
    for i in range(dose_inf_df.shape[0]):
        if ((dose_inf_df.loc[i, 'Genomes/well'] / scale) <= 0):
            num_zeros += 1
    print('---num_zeros', num_zeros)

    GEN_WELL_DATA = np.empty(dose_inf_df.shape[0] - num_zeros, int)  # Initial amount of GFP genomes in simulation
    GEN_CELL_DATA = np.empty(dose_inf_df.shape[0] - num_zeros, float)
    INF_CELL_DATA = np.empty(dose_inf_df.shape[0] - num_zeros, float)

    i = 0
    j = 0
    while (i < dose_inf_df.shape[0]):
        if (dose_inf_df.loc[dose_inf_df.shape[0] - i - 1, 'Genomes/well'] > 0):
            #print("i=", dose_inf_df.shape[0] - i - 1, ", dose_inf_df.loc[i, 'Genomes/well']=", dose_inf_df.loc[dose_inf_df.shape[0] - i - 1, 'Genomes/well'])
            if (scale==None):
                GEN_WELL_DATA[j] = int(dose_inf_df.loc[dose_inf_df.shape[0] - i - 1, 'Genomes/well'])
            else:
                GEN_WELL_DATA[j] = int(dose_inf_df.loc[dose_inf_df.shape[0] - i - 1, 'Genomes/well'] / scale)

            GEN_CELL_DATA[j] = dose_inf_df.loc[dose_inf_df.shape[0] - i - 1, 'Genomes/cell']
            INF_CELL_DATA[j] = dose_inf_df.loc[dose_inf_df.shape[0] - i - 1, 'IU/cell']
            j+=1
        i+=1
    
    return GEN_WELL_DATA, GEN_CELL_DATA, INF_CELL_DATA, num_zeros
#====================================================================
def PrepareParamList(GEN_CELL_DATA, simul_name, cell_count, PARAM_DICT):
    assert simul_name in ['clump', 'comp', 'acc_dam', 'clump_comp', 'clump_acc_dam', 'null'], simul_name+" must be 'clump', 'comp', 'acc_dam', 'clump_comp', 'clump_acc_dam' or 'null'."
    #----------------------------------------------------------------
    PARAM_LIST_STR = ['simul_name='+str(simul_name), 'num_simulations='+str(PARAM_DICT['num_simulations']), 'sheet='+str(PARAM_DICT['sheet']),
                      'scale='+str(PARAM_DICT['scale']), 'cell_count='+str(cell_count), 
                      'muG='+str(PARAM_DICT['i_mean']), 'sigmaG='+str(PARAM_DICT['i_stdev']),
                      'muR='+str(PARAM_DICT['r_mean']), 'sigmaR='+str(PARAM_DICT['r_stdev']), 
                      'gamma='+str(PARAM_DICT['gamma']), 'vMax='+str(PARAM_DICT['vMax'])]
    #----------------------------------------------------------------
    if ('clump' in simul_name):
        PARAM_LIST_STR.append('scheme='+str(PARAM_DICT['scheme']))
        PARAM_LIST_STR.append('distribution='+str(PARAM_DICT['distribution']))
        PARAM_LIST_STR.append('ub='+str(PARAM_DICT['ub']))
        PARAM_LIST_STR.append('lb='+str(PARAM_DICT['lb']))
        PARAM_LIST_STR.append('mean='+str(PARAM_DICT['mean']))

    if ('comp' in simul_name):
        PARAM_LIST_STR.append('kappa='+str(PARAM_DICT['kappa']))

    if ('acc_dam' in simul_name):
        PARAM_LIST_STR.append('beta='+str(PARAM_DICT['beta']))

    if (simul_name == 'null'):
        PARAM_LIST_STR.append('b='+str(PARAM_DICT['b']))
    #----------------------------------------------------------------
    for i in range(GEN_CELL_DATA.shape[0] - len(PARAM_LIST_STR)):
        PARAM_LIST_STR.append('FILLER')
    #----------------------------------------------------------------
    return PARAM_LIST_STR
#====================================================================
def Innoculation(Viron, Cell, gamma, acc_dam=False, beta=0):
    is_successful = False # Boolean value to see if the infection was successful

    resistivity = Cell.r
    if acc_dam:
        resistivity = Cell.r + beta * Cell.inter
        
    if (Viron.i > resistivity):
        is_successful = True
        #print(" - Successful infection | i = ", Viron.i, ", r = ", Cell.r, " | Cell.inter = ", Cell.inter)
        Cell.infG = True
        Cell.numG += 1

    else: # Second chance
        P_i = np.exp((Viron.i - resistivity) * gamma)
        if (np.random.uniform() <= P_i):
            is_successful = True
            Cell.infG = True
            Cell.numG += 1
            #print("Second chance | i = ", Viron.i, ", r = ", Cell.r, ", P_i = ", P_i)
        #else:
            #print("Failed infection | i = ", Viron.i, ", r = ", Cell.r)

    if acc_dam:
        #print("ACCC DAAMMM", Cell.inter)
        Cell.inter += 1
        #print("ACCC DAAMMM res=", round(resistivity, 3), ", inf=", round(Viron.i, 3), ", beta=", beta, ", gamma=", gamma, ", Cell.inter=", Cell.inter)
    return is_successful
#====================================================================
def model(x, params):
    try: # Get parameters
        g = params['gamma'].value
        n = params['n'].value
    except KeyError:
        g, n = params

    y = 1 - np.exp(-g * x**n)

    return y
#====================================================================
def ClumpMetrics2(NUM_CLUMPS):
    ''' Initialize lists '''
    CLUMP_IDS = [] # Integer assigned to each viron denoting which clump it is a part of
    CLUMP_SIZES = np.arange(1, np.shape(NUM_CLUMPS)[0]+1)
    #----------------------------------------------------------------
    clump_num = 0 # Assign each virion an integer 'ID', virions in the same clump will have the same ID
    index = 0
    for i in range(len(NUM_CLUMPS)): 
        for j in range(int(NUM_CLUMPS[i])):
            for k in range(CLUMP_SIZES[i]):
                #print("clump size ", CLUMP_SIZES[i]+1, ", clump number = ", j, ", clump_num = ", clump_num, ", index = ", index)
                index = index + 1
                CLUMP_IDS.append(clump_num)
            clump_num = clump_num + 1

    total = 0
    for i in range(len(NUM_CLUMPS)):
        total += (i+1) * NUM_CLUMPS[i]
    print("NUM_CLUMPS =", NUM_CLUMPS, ", largest clump =", len(NUM_CLUMPS)+1, ", virions in NUM_CLUMPS =", total, " | len(CLUMP_IDS) =", len(CLUMP_IDS))
    print("~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-")
    return CLUMP_IDS
#====================================================================
def CalcNumClumps(total_virions, max_virions_in_clump, diameter_nz,
                  mean_clump_diam=173.884, std_clump_diam=200, skew_clump_diam=3.569,
                  mean_virion_diam=230, lb_virion_diam=150, ub_virion_diam=300, std_virion_diam=35,
                  scheme='sphere_packing', distribution='fixed', tolerance=0.1):
    assert (scheme in ['linear', 'regular_polygon', 'sphere_packing']), "scheme must be 'linear', 'regular_polygon', or 'sphere_packing'."
    assert (distribution in ['uniform', 'normal', 'fixed']), "dist must be 'uniform', 'normal' or 'fixed'."
    #----------------------------------------------------------------
    num_clumps_of_size_i = np.zeros(max_virions_in_clump, int) # i refers to position in list, pos. 0 means clump size 1

    virions_used = 0

    while (virions_used < total_virions):
        virion_diam = GetVirionDiam(distribution, mean_virion_diam, lb_virion_diam, ub_virion_diam, std_virion_diam)
        broken = False
        err = 9999
        while (err >= tolerance):
            virions_in_clump = -999 # Put in scope

            diameter = skewnorm.rvs(skew_clump_diam, mean_clump_diam, std_clump_diam, size=1)[0] # Best parameters so far
            if (diameter < lb_virion_diam or (distribution=='fixed' and diameter <= 1.5*mean_virion_diam)):
                virions_in_clump = 1
                broken = True
                break
            elif (diameter > max(diameter_nz)):
                diameter = max(diameter_nz)
        
            virions_in_clump = GetVirionsInClumpFromDiameter(scheme, diameter, virion_diam)

            virions_in_clump_low = int(virions_in_clump)
            virions_in_clump_hi = virions_in_clump_low + 1

            err = min([abs(virions_in_clump - virions_in_clump_low), abs(virions_in_clump - virions_in_clump_hi)])

        #print("broken=", broken, ", clump_diam=", diameter, ", vir=", virions_in_clump)
        virions_in_clump = int(round(virions_in_clump))
        num_clumps_of_size_i[virions_in_clump-1] += 1

        virions_used += virions_in_clump
    # print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    # print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    #print("num_clumps_of_size_i=", num_clumps_of_size_i)
    # print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    # print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    #----------------------------------------------------------------
    ''' Clean up and correct any errors '''
    for valid_end in range(num_clumps_of_size_i.shape[0]-1, 0, -1):
        if (num_clumps_of_size_i[valid_end] > 0):
            break
    num_clumps_of_size_i = num_clumps_of_size_i[:valid_end+1]
    #print("num_clumps_of_size_i=", np.shape(num_clumps_of_size_i), num_clumps_of_size_i)

    total = 0
    for i in range(num_clumps_of_size_i.shape[0]):
        total += (i+1) * num_clumps_of_size_i[i]

    #print("pre  :", num_clumps_of_size_i, ", len:", num_clumps_of_size_i.shape[0], ',', total, ',', total_virions)
    num_clumps_of_size_i[0] += (total_virions-total)
    #print("post :", num_clumps_of_size_i, ", len:", num_clumps_of_size_i.shape[0], ',', total, ',', total_virions)

    if (num_clumps_of_size_i[0] < 0):
        num_clumps_of_size_i[0] += (num_clumps_of_size_i.shape[0]+1) * num_clumps_of_size_i[-1]
        num_clumps_of_size_i[-1] -= 1
    for valid_end in range(num_clumps_of_size_i.shape[0]-1, 0, -1):
        if (num_clumps_of_size_i[valid_end] > 0):
            break
    num_clumps_of_size_i = num_clumps_of_size_i[:valid_end+1]
    #print("final:", num_clumps_of_size_i, ", len:", num_clumps_of_size_i.shape[0], ',', total, ',', total_virions)








    # ''' WORKS FOR LIN AND POLY Clean up and correct any errors '''
    # for valid_end in range(num_clumps_of_size_i.shape[0]-1, 0, -1):
    #     if (num_clumps_of_size_i[valid_end] > 0):
    #         break
    # num_clumps_of_size_i = num_clumps_of_size_i[:valid_end+1]
    # print("num_clumps_of_size_i=", num_clumps_of_size_i)

    # total = 0
    # for i in range(num_clumps_of_size_i.shape[0]):
    #     total += (i+1) * num_clumps_of_size_i[i]

    # #print("pre  :", num_clumps_of_size_i, ", len:", num_clumps_of_size_i.shape[0], ',', total, ',', total_virions)
    # num_clumps_of_size_i[0] += (total_virions-total)
    # #print("post :", num_clumps_of_size_i, ", len:", num_clumps_of_size_i.shape[0], ',', total, ',', total_virions)

    # if (num_clumps_of_size_i[0] < 0):
    #     num_clumps_of_size_i[0] += (num_clumps_of_size_i.shape[0]+1) * num_clumps_of_size_i[-1]
    #     num_clumps_of_size_i[-1] -= 1
    # for valid_end in range(num_clumps_of_size_i.shape[0]-1, 0, -1):
    #     if (num_clumps_of_size_i[valid_end] > 0):
    #         break
    # num_clumps_of_size_i = num_clumps_of_size_i[:valid_end+1]
    # #print("final:", num_clumps_of_size_i, ", len:", num_clumps_of_size_i.shape[0], ',', total, ',', total_virions)

    return num_clumps_of_size_i
#====================================================================
def GenClumpFromDiameter(diameter, mean, lb, ub, scheme='linear', dist='uniform', stdev=35, target_x=None, target_prob=0.004,
                         tolerance=0.01):
    assert (scheme in ['linear', 'regular_polygon', 'sphere_packing']), "scheme must be 'linear', 'regular_polygon', or 'sphere_packing'."
    assert (dist in ['uniform', 'normal', 'fixed']), "dist must be 'uniform', 'normal' or 'fixed'."
    #----------------------------------------------------------------
    if ((dist=='uniform') or (dist=='normal')):
        err = 9999
        while (err >= tolerance):
            if (dist == 'uniform'):
                virion_diam = np.random.uniform(lb, ub, 1)[0]
                if (virion_diam >= diameter):
                    virions_in_clump = 1
                    break

            elif (dist == 'normal'):
                virion_diam = np.random.normal(mean, stdev, 1)[0]
                if (virion_diam >= diameter):
                    virions_in_clump = 1
                    break
                else:
                    if (virion_diam < lb):
                        virion_diam = lb
                    elif (virion_diam > ub):
                        virion_diam = ub
            
            if (scheme=='linear'):
                virions_in_clump = diameter / virion_diam
            elif (scheme=='regular_polygon'):
                virions_in_clump = np.pi / np.arcsin(virion_diam / diameter)
                #print("POLY:", virion_diam, diameter, np.arcsin(virion_diam / diameter), virions_in_clump)
            elif (scheme=='sphere_packing'):
                virions_in_clump = 0.64 * (diameter / virion_diam)**3
            
            virions_in_clump_low = int(virions_in_clump)
            virions_in_clump_hi  = virions_in_clump_low + 1

            err = min([abs(virions_in_clump - virions_in_clump_low), abs(virions_in_clump - virions_in_clump_hi)])

            print("diam=", diameter, ", virion_diam=", round(virion_diam, 5), ", virions_in_clump=", round(virions_in_clump, 5), ", err=", round(err, 5))
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    elif (dist=='fixed'):
        virion_diam = mean
        if (scheme=='linear'):
            virions_in_clump = diameter / virion_diam
        elif (scheme=='regular_polygon'):
            virions_in_clump = np.pi / np.arcsin(virion_diam / diameter)
            #print("POLY:", virion_diam, diameter, np.arcsin(virion_diam / diameter), virions_in_clump)
        elif (scheme=='sphere_packing'):
            virions_in_clump = 0.64 * (diameter / virion_diam)**3
    #------------------------------------------------------------
    num_virions = round(virions_in_clump)
    print("diam=", diameter, ", virion_diam=", round(virion_diam, 5), ", virions_in_clump=", round(virions_in_clump, 5), ", err=", round(err, 5))
    #----------------------------------------------------------------
    return num_virions, virion_diam

#====================================================================
def GetVirionDiam(distribution, mean_virion_diam=230, lb_virion_diam=150, ub_virion_diam=300, std_virion_diam=35):
    if (distribution == 'uniform'):
        virion_diam = np.random.uniform(lb_virion_diam, ub_virion_diam, 1)[0]

    elif (distribution == 'normal'):
        virion_diam = np.random.normal(mean_virion_diam, std_virion_diam, 1)[0]
        if (virion_diam < lb_virion_diam):
            virion_diam = lb_virion_diam
        elif (virion_diam > ub_virion_diam):
            virion_diam = ub_virion_diam
    
    elif (distribution=='fixed'):
        virion_diam = mean_virion_diam
    
    return virion_diam
#====================================================================
def GetVirionsInClumpFromDiameter(scheme, diameter, virion_diam):
    if (scheme=='linear'):
        virions_in_clump = diameter / virion_diam
    elif (scheme=='regular_polygon'):
        virions_in_clump = np.pi / np.arcsin(virion_diam / diameter)
        #print("POLY:", virion_diam, diameter, np.arcsin(virion_diam / diameter), virions_in_clump)
    elif (scheme=='sphere_packing'):
        virions_in_clump = 0.64 * (diameter / virion_diam)**3
    
    return virions_in_clump
#====================================================================
def Compensate(best_inf, Virion, kappa):
    #print("best_inf:", best_inf, type(best_inf), ",inf:", Virion.i, type(Virion.i), "k:", kappa, type(kappa))
    new_inf = 0
    if (Virion.i >= 0):
        new_inf = best_inf - kappa * Virion.i
    else:
        new_inf = best_inf + kappa * Virion.i
    return new_inf
#====================================================================
class Cell:
    def __init__(self, is_inf_G, num_infG, resistivity, interations):
        self.infG = is_inf_G     # Boolean value: 0=unifected, 1=infected with GFP
        self.numG = num_infG     # Integer: Number of GFP virons in cell
        self.r = resistivity     # Float: Represents the resistance of the cell to infection
        self.inter = interations # Integer: Keeps track of how many virions have interacted with the cell
#====================================================================
class Virion:
    def __init__(self, infectivity, marker, clump_num=None):
        self.i = infectivity # Float: Represents ability of viron to infect a cell
        self.strain = marker # Bool: GFP or Venus
        self.num = clump_num # int: Clump number to which virus is a part of