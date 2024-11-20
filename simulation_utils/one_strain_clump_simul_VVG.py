import sys
sys.dont_write_bytecode = True
import numpy as np
import os
import itertools

sys.path.append(os.getcwd())
from simulation_utils.utils import ClumpMetrics2, CalcNumClumps, Innoculation, Cell, Virion, Compensate, DiamFunc
from misc.misc_utils import FlattenMeans, Trapezoid
#====================================================================
def SimulateClump(GEN_WELL_DATA, PARAM_DICT, cell_count, clump_size_df, simul_name, sheet_name, save_clump_info=False, CLUMP_DIST_PARAMS_DICT=None):
    assert simul_name in ['clump', 'clump_comp', 'clump_acc_dam', 'var_clump_diam'], " >> Check config.yml\n >> 'simul_name' must be: 'clump', 'clump_acc_dam', 'clump_comp', 'var_clump_diam'. Got: "+simul_name
    assert PARAM_DICT['scheme'] in ['linear', 'regular_polygon', 'sphere_packing'], " >> Check config.yml\n >> 'scheme' must be: 'linear', 'regular_polygon', 'sphere_packing'. Got:"+PARAM_DICT['scheme']
    #----------------------------------------------------------------
    INF_WELL_SIMUL     = np.empty(GEN_WELL_DATA.shape[0], int)
    TOTAL_INTERACTIONS = np.empty(GEN_WELL_DATA.shape[0], int)
    #----------------------------------------------------------------
    ''' Prepare information for clumps '''
    if save_clump_info:
        CLUMP_DICT = {str(x): [] for x in GEN_WELL_DATA}
    #----------------------------------------------------------------
    ''' Import data '''
    diameter = clump_size_df.pop('d.nm')
    diameter = np.asarray(diameter)
    #----------------------------------------------------------------
    ''' Average '''
    cutoff = 91.28

    means_clump_size_df = clump_size_df.filter(like='Mean')
    stds_clump_size_df  = clump_size_df.filter(like='SD')

    mean_of_means = means_clump_size_df.mean(axis=1)
    mean_of_means = np.asarray(mean_of_means)

    mean_of_means = FlattenMeans(diameter, mean_of_means, cutoff)

    mean_of_means /= Trapezoid(diameter, mean_of_means) # Normalize by dividing by area  
    #----------------------------------------------------------------
    ''' Filter to include only non-zero portions '''
    diameter_nz = []
    mean_of_means_nz = []

    for i in range(mean_of_means.shape[0]):
        if (mean_of_means[i] > 0):
            diameter_nz.append(diameter[i])
            mean_of_means_nz.append(mean_of_means[i])
    #----------------------------------------------------------------
    ''' Calculate the maximum number of virions which could be in a clump '''
    max_virions_in_clump = -999
    if (PARAM_DICT['scheme'] == 'linear'):
        max_virions_in_clump = round(max(diameter_nz) / PARAM_DICT['lb'])
        print("linear | diam=", max(diameter_nz), ", max=", max_virions_in_clump)

    elif (PARAM_DICT['scheme']=='regular_polygon'):
        max_virions_in_clump = round(np.pi / np.arcsin(PARAM_DICT['lb'] / max(diameter_nz)))
        print("poly | diam=", max(diameter_nz), ", max=", max_virions_in_clump)

    elif (PARAM_DICT['scheme']=='sphere_packing'):
        max_virions_in_clump = round(0.64 * (max(diameter_nz) / PARAM_DICT['lb'])**3)
        print("sphere | diam=", max(diameter_nz), ", max=", max_virions_in_clump)
    #----------------------------------------------------------------
    ''' Simulate infections '''
    for init in range(len(GEN_WELL_DATA)): # Iterate thorugh each experiment
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        mean_clump_diam = -999 # Put in scope
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        print("=======================================================================")
        
        if (simul_name == 'clump'):
            mean_clump_diam = 173.844 # Derived from '2022_10_27_TB_size_distribution'
            print(simul_name, ", GEN_WELL_DATA[", init, "] =", GEN_WELL_DATA[init], "| vMax =", PARAM_DICT['vMax'], "| scheme =", PARAM_DICT['scheme'], "| scale =", PARAM_DICT['scale'])
        elif (simul_name == 'clump_comp'):
            mean_clump_diam = 173.844 # Derived from '2022_10_27_TB_size_distribution'
            print(simul_name, ", GEN_WELL_DATA[", init, "] =", GEN_WELL_DATA[init], "| vMax =", PARAM_DICT['vMax'], "| scheme =", PARAM_DICT['scheme'], "| kappa =", PARAM_DICT['kappa'], "| scale =", PARAM_DICT['scale'])
        elif (simul_name == 'clump_acc_dam'):
            mean_clump_diam = 173.844 # Derived from '2022_10_27_TB_size_distribution'
            print(simul_name, ", GEN_WELL_DATA[", init, "] =", GEN_WELL_DATA[init], "| vMax =", PARAM_DICT['vMax'], "| scheme =", PARAM_DICT['scheme'], "| beta =", PARAM_DICT['beta'], "| scale =", PARAM_DICT['scale'])
        
        elif (simul_name == 'var_clump_diam'):
            mean_clump_diam = DiamFunc(GEN_WELL_DATA[init], PARAM_DICT)
            print(simul_name, ", GEN_WELL_DATA[", init, "] =", GEN_WELL_DATA[init], "| vMax =", PARAM_DICT['vMax'], "| mean diam =", round(mean_clump_diam, 3), "(173.884) |", PARAM_DICT['scheme'], "| scale =", PARAM_DICT['scale'])
        
        print(" - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ")
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ''' Set up virus and cell populations '''
        GFP_POOL = []
        CELL_POOL = []
        num_infG = 0
        total = 0
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # CLUMP_IDS: 0,1,2,3,4,4,4,5,6,7,7,8, ...
        # CLUMP_NUMS: # size 1 clumps, # size 2 clumps, ...
        # CLUMP_SIZES: 1,2,3,...,max clump size
        # radius: max clump size
        CLUMP_NUMS = -99999999
        if ((sheet_name == 'use_with_size_distribution') and not (simul_name == 'var_clump_diam')):
            keys_   = list(CLUMP_DIST_PARAMS_DICT.keys())
            keys__  = [int(float(key) / PARAM_DICT['scale']) for key in keys_]
            key_idx = keys__.index(GEN_WELL_DATA[init])
            PARAMS  = CLUMP_DIST_PARAMS_DICT[keys_[key_idx]]
            print(" - - - - - skew =", PARAMS[0], ", mean =", PARAMS[1], ", stdev =", PARAMS[2], "- - - - - ")

            CLUMP_NUMS = CalcNumClumps(GEN_WELL_DATA[init], 
                                        max_virions_in_clump,
                                        diameter_nz,
                                        skew_clump_diam=PARAMS[0],
                                        mean_clump_diam=PARAMS[1],
                                        std_clump_diam=PARAMS[2],
                                        mean_virion_diam=PARAM_DICT['mean'], 
                                        lb_virion_diam=PARAM_DICT['lb'], 
                                        ub_virion_diam=PARAM_DICT['ub'],
                                        scheme=PARAM_DICT['scheme'], 
                                        distribution=PARAM_DICT['distribution'])
        else:
            CLUMP_NUMS = CalcNumClumps(GEN_WELL_DATA[init], 
                                        max_virions_in_clump,
                                        diameter_nz,
                                        mean_clump_diam=mean_clump_diam,
                                        mean_virion_diam=PARAM_DICT['mean'], 
                                        lb_virion_diam=PARAM_DICT['lb'], 
                                        ub_virion_diam=PARAM_DICT['ub'],
                                        scheme=PARAM_DICT['scheme'], 
                                        distribution=PARAM_DICT['distribution'])
            
        CLUMP_SIZES = np.arange(1, len(CLUMP_NUMS)+1)
        radius = len(CLUMP_NUMS)

        CLUMP_IDS = ClumpMetrics2(CLUMP_NUMS)
        
        if save_clump_info:
            CLUMP_DICT[str(GEN_WELL_DATA[init])] = CLUMP_NUMS.tolist()
        #================================================================
        ''' Create virus population '''

        for j in range(GEN_WELL_DATA[init]):
            i = np.random.normal(PARAM_DICT['i_mean'], PARAM_DICT['i_stdev'])
            #print("j = ", j, "CLUMP_IDS[j] = ", CLUMP_IDS[j])
            GFP_POOL.append(Virion(i, True, CLUMP_IDS[j])) # Assign each viron with an infectivity, marker, and clump ID
        #================================================================
        ''' Create a cell and then try to infect it '''

        for cell_num in range(cell_count): # Orig. range(cell_count)
            #print(" == == == ==", cell_num, "== == == ==")
            #print(" + + + Cell number ", cell_num, "+ + +")
            r = np.random.normal(PARAM_DICT['r_mean'], PARAM_DICT['r_stdev'])
            CELL_POOL.append(Cell(False, 0, r, 0))
            #============================================================
            if (len(GFP_POOL) > 0):
                num_interactions = np.random.poisson((GEN_WELL_DATA[init] / PARAM_DICT['vMax'])) # Calculate the number of clumps that will visit the cell
                
                VIRION_ATTACKERS = []
                VIRION_ATTACKERS_IDX = []
                CLUMPS_USED = []
                TO_REMOVE_IDX = []
                
                if (num_interactions > 0):
                    for _ in range(int(num_interactions)):
                        # - - - - - - - - - - - - - - - - - - - - - - - -
                        # Select the virions that will try to infect the cell
                        index = np.random.randint(0, len(GFP_POOL), size=1)[0] # Select a virion from the pool
                        target_clump_id = GFP_POOL[index].num
                        
                        if not (target_clump_id in CLUMPS_USED):  # Avoid creating another clump with the same ID
                            #print(" >> target_clump_id", target_clump_id, ", radius", radius)
                            for aa in range(index - radius, index + radius + 1):
                                if ((aa >= 0) and (aa < len(GFP_POOL))): 
                                    #print("aa:", aa, ", curr ID:", GFP_POOL[aa].num,",targ. ID:", target_clump_id,",idx:", index, ", int.:", num_interactions, ", rad:", radius, ", a:",a,",cell_num:",cell_num, "vir_left:", len(GFP_POOL))
                                    if (GFP_POOL[aa].num == target_clump_id): # Find virions with the same clump ID number
                                        VIRION_ATTACKERS.append(GFP_POOL[aa])
                                        VIRION_ATTACKERS_IDX.append(aa)
                        
                        CLUMPS_USED.append(target_clump_id)

                    #print(" >>", 'num_interactions =', num_interactions, ', ', _, len(VIRION_ATTACKERS), type(VIRION_ATTACKERS))
                    
                    try: # Flatten in case a list of lists was created
                        VIRION_ATTACKERS     = list(itertools.chain.from_iterable(VIRION_ATTACKERS)) 
                        VIRION_ATTACKERS_IDX = list(itertools.chain.from_iterable(VIRION_ATTACKERS_IDX))
                    except TypeError: # Don't flatten
                        pass
                    # - - - - - - - - - - - - - - - - - - - - - - - -
                    # Vanilla clump simulation
                    if (simul_name == 'clump' or simul_name == 'var_clump_diam'):
                        for b in range(len(VIRION_ATTACKERS)):
                            is_successful = Innoculation(VIRION_ATTACKERS[b], 
                                                        CELL_POOL[cell_num], 
                                                        PARAM_DICT['gamma'])
                            total += 1
                            if is_successful:
                                TO_REMOVE_IDX.append(VIRION_ATTACKERS_IDX[b])
                        #print("TO_REMOVE_IDX=", TO_REMOVE_IDX, ", len(GFP_POOL)=", len(GFP_POOL))
                    # - - - - - - - - - - - - - - - - - - - - - - - -
                    # Clumps + Compensation
                    elif (simul_name == 'clump_comp'):
                        best_inf = -10000
                        for idx in range(len(VIRION_ATTACKERS)):
                            if (VIRION_ATTACKERS[idx].i > best_inf):
                                best_inf = VIRION_ATTACKERS[idx].i # Identify the highest infectivity value

                        for b in range(len(VIRION_ATTACKERS)):
                            VIRION_ATTACKERS[b].i = Compensate(best_inf, VIRION_ATTACKERS[b], PARAM_DICT['kappa'])
                            is_successful = Innoculation(VIRION_ATTACKERS[b], 
                                                        CELL_POOL[cell_num], 
                                                        PARAM_DICT['gamma'])
                            total += 1
                            #print("aaa=", aaa, TO_REMOVE_IDX[aaa], "len(GFP_POOL)=", len(GFP_POOL))
                            if is_successful:
                                TO_REMOVE_IDX.append(VIRION_ATTACKERS_IDX[b])
                    # - - - - - - - - - - - - - - - - - - - - - - - -
                    # Clumps + Accrued damage
                    elif (simul_name == 'clump_acc_dam'):
                        for b in range(len(VIRION_ATTACKERS)):
                            is_successful = Innoculation(VIRION_ATTACKERS[b], 
                                                        CELL_POOL[cell_num], 
                                                        gamma=PARAM_DICT['gamma'], 
                                                        acc_dam=True, 
                                                        beta=PARAM_DICT['beta'])

                            total += 1
                            if is_successful:
                                TO_REMOVE_IDX.append(VIRION_ATTACKERS_IDX[b])    
                                #print("----- match ----- inf:", is_successful, ",len:", len(GFP_POOL))
                    # - - - - - - - - - - - - - - - - - - - - - - - -
                    if (PARAM_DICT['remove'] == 1):
                        # Remove virions which successfully infected from pool
                        TO_REMOVE_IDX = list(sorted(set(TO_REMOVE_IDX), reverse=True))

                        for idx in range(len(TO_REMOVE_IDX)):
                            if (len(TO_REMOVE_IDX) >= len(GFP_POOL)):
                                print('BREAKING idx_len=', len(TO_REMOVE_IDX), ", virions left=",len(GFP_POOL))
                                break
                            else:
                                #print("idx=",idx,", idx=", TO_REMOVE_IDX[idx], ', len=', len(TO_REMOVE_IDX), ", virions left=", len(GFP_POOL))
                                #GFP_POOL.remove(GFP_POOL[TO_REMOVE_IDX[idx]]) works but slow
                                #print("cell_num=", cell_num, ", len(GFP)=", len(GFP_POOL), ", len(REMOVE)=", len(TO_REMOVE_IDX), ", idx=", idx, ", r_idx=", TO_REMOVE_IDX[idx])
                                del GFP_POOL[TO_REMOVE_IDX[idx]]

                    #if (cell_num == 0 or cell_num == 500 or cell_num == 1500 or cell_num == 2299):
                    #    print("cell=", cell_num, ", poisson:", num_clump_interactions, ",len(CURR_CLUMP):", len(VIRIONS_IN_CURR_CLUMP), ", GFP_GENOMES[",init,"]:", GEN_WELL_DATA[init], ", len(GFP):", len(GFP_POOL), ", total/cell:", round(total/cell_count, 5))

                if (cell_num == 0 or cell_num == int(.33 * cell_count) or cell_num == int(.66 * cell_count) or cell_num == cell_count-1):
                    print("cell:", cell_num, ", poisson:", num_interactions, ", GFP_GENOMES[",init,"]:", GEN_WELL_DATA[init], ", len(GFP):", len(GFP_POOL), ", total/cell:", round(total/cell_count, 5), ", total/genomes:", round(total/GEN_WELL_DATA[init], 5))
            #--------------------------------------------------------
            else:
                print("[][][][] OUT OF VIRIONS [][][][]")
                break
        #============================================================
        for j in range(cell_count):
            if (CELL_POOL[j].infG == True):
                num_infG += 1
            
        #print(" + + + + + + + + + + + + + + ")
        #print(str(GEN_WELL_DATA[init]), ":", CLUMP_DICT[str(GEN_WELL_DATA[init])])
        #print(" + + + + + + + + + + + + + + ")
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        print("~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-")
        print("num infG =", num_infG, "| cell_count =", cell_count, "| avg. inter/cell =", round(total / cell_count, 4), "| avg. inter/virion =", round(total / GEN_WELL_DATA[init], 4))
        
        INF_WELL_SIMUL[init]     = num_infG
        TOTAL_INTERACTIONS[init] = total
    #----------------------------------------------------------------
    if save_clump_info:
        return INF_WELL_SIMUL, TOTAL_INTERACTIONS, CLUMP_DICT
    else:
        return INF_WELL_SIMUL, TOTAL_INTERACTIONS