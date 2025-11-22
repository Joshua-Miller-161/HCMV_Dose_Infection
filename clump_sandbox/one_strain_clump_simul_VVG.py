import sys
sys.dont_write_bytecode = True
import numpy as np
import os
import itertools
import yaml

sys.path.append(os.getcwd())
from simulation_utils.utils import Innoculation, Cell, Virion, Compensate
from clump_utils import CreateClumps
#====================================================================
def SimulateClump(GENOMES_WELL,
                  PARAM_DICT,
                  simul_name,
                  save_clump_info=False):
    '''
    Arguments
     - GENOMES_WELL <ndarray>: An array containing the amounts of virions you want to simulate. [10 virions, 20 virions, ..., 1000 virions]
     - PARAM_DICT <dict>: Dictionary of simulation parameters obtained from sandbox_config.yml
     - simul_name <str>: Name of the simulation type
     - save_clump_info [Bool]: Whether or not to save the clump size frequency distribution

    Returns
     - INF_WELL_SIMUL <ndarray>: An array with shape GENOMES_WELL which stores the numbers of infected cells at each genomes/well values in GENOMES_WELL
     - TOTAL_INTERACTIONS <ndarray>: An array with shape GENOMES_WELL which stores how many times all cells were visited at each genomes/well value in GENOMES_WELL
     - CLUMP_DICT [dict]: A dictionary with form {genomes/well: clump freq. dist., genomes/well: clump_freq. dist.,...}, keeps track of the clump distribution for genomes/well value in GENOMES_WELL
    '''
    assert simul_name in ['clump', 'clump_comp', 'clump_acc_dam'], " >> Check config.yml\n >> 'simul_name' must be: 'clump', 'clump_acc_dam', 'clump_comp', 'var_clump_diam'. Got: "+simul_name
    #----------------------------------------------------------------
    cell_count = int(PARAM_DICT['cell_count'])
    #----------------------------------------------------------------
    INF_WELL_SIMUL     = np.empty(GENOMES_WELL.shape[0], int)
    TOTAL_INTERACTIONS = np.empty(GENOMES_WELL.shape[0], int)
    #----------------------------------------------------------------
    ''' Prepare information for clumps '''
    if save_clump_info:
        CLUMP_DICT = {str(int(x)): [] for x in GENOMES_WELL}
    #----------------------------------------------------------------
    ''' Simulate infections '''
    for iter in range(len(GENOMES_WELL)): # Iterate thorugh each experiment
        print("===============================================================================")
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

        CLUMP_IDS, CLUMP_NUMS, CLUMP_SIZES, radius = CreateClumps(int(GENOMES_WELL[iter]),
                                                                  PARAM_DICT)

        if save_clump_info:
            CLUMP_DICT[str(int(GENOMES_WELL[iter]))] = CLUMP_NUMS.tolist()

        #================================================================
        ''' Create virus population '''

        for j in range(int(GENOMES_WELL[iter])):
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
                num_interactions = np.random.poisson((GENOMES_WELL[iter] / PARAM_DICT['vMax'])) # Calculate the number of clumps that will visit the cell
                
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

                    total += len(VIRION_ATTACKERS) # TEST
                    # - - - - - - - - - - - - - - - - - - - - - - - -
                    # Vanilla clump simulation
                    if (simul_name == 'clump' or simul_name == 'var_clump_diam'):
                        for b in range(len(VIRION_ATTACKERS)):
                            is_successful = Innoculation(VIRION_ATTACKERS[b], 
                                                         CELL_POOL[cell_num], 
                                                         PARAM_DICT['gamma'])
                            #total += 1 ORIGINAL
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
                        if (PARAM_DICT['remove_clump'] == 1):
                            for i in range(len(TO_REMOVE_IDX)):
                                for j in range(len(GFP_POOL)):
                                    if (GFP_POOL[TO_REMOVE_IDX[i]].num == GFP_POOL[j].num): # Find virions in same clump as successfully infecting virion
                                        TO_REMOVE_IDX.append(j)
                            
                            #print("TO_REMOVE_IDX", TO_REMOVE_IDX, type(TO_REMOVE_IDX))
                            TO_REMOVE_IDX = list(sorted(set(TO_REMOVE_IDX), reverse=True))

                        else:# Remove virions which successfully infected from pool
                            TO_REMOVE_IDX = list(sorted(set(TO_REMOVE_IDX), reverse=True))

                        for idx in range(len(TO_REMOVE_IDX)):
                            if (len(TO_REMOVE_IDX) >= len(GFP_POOL)):
                                print('BREAKING idx_len=', len(TO_REMOVE_IDX), ", virions left=",len(GFP_POOL))
                                break
                            else:
                                #print("idx=",idx,", idx=", TO_REMOVE_IDX[idx], ', len=', len(TO_REMOVE_IDX), ", virions left=", len(GFP_POOL))
                                #GFP_POOL.remove(GFP_POOL[TO_REMOVE_IDX[idx]]) works but slow
                                print("cell_num=", cell_num, ", len(GFP)=", len(GFP_POOL), ", len(REMOVE)=", len(TO_REMOVE_IDX), ", idx=", idx, ", r_idx=", TO_REMOVE_IDX[idx], ", clump_num=", GFP_POOL[TO_REMOVE_IDX[idx]].num)
                                del GFP_POOL[TO_REMOVE_IDX[idx]]

                    #if (cell_num == 0 or cell_num == 500 or cell_num == 1500 or cell_num == 2299):
                    #    print("cell=", cell_num, ", poisson:", num_clump_interactions, ",len(CURR_CLUMP):", len(VIRIONS_IN_CURR_CLUMP), ", GFP_GENOMES[",init,"]:", GENOMES_WELL[init], ", len(GFP):", len(GFP_POOL), ", total/cell:", round(total/cell_count, 5))

                if (cell_num == 0 or cell_num == int(.33 * cell_count) or cell_num == int(.66 * cell_count) or cell_num == cell_count-1):
                    print("cell:", cell_num, ", poisson:", num_interactions, ", GFP_GENOMES[",iter,"]:", GENOMES_WELL[iter], ", len(GFP):", len(GFP_POOL), ", total/cell:", round(total/cell_count, 5), ", total/genomes:", round(total/GENOMES_WELL[iter], 5))
            #--------------------------------------------------------
            else:
                print("[][][][] OUT OF VIRIONS [][][][]")
                break
        #============================================================
        for j in range(cell_count):
            if (CELL_POOL[j].infG == True):
                num_infG += 1
            
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        print("-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-")
        print("num infG =", num_infG, "| cell_count =", cell_count, "| avg. inter/cell =", round(total / cell_count, 4), "| avg. inter/virion =", round(total / GENOMES_WELL[iter], 4))
        
        INF_WELL_SIMUL[iter]     = num_infG
        TOTAL_INTERACTIONS[iter] = total
    #----------------------------------------------------------------
    if save_clump_info:
        return INF_WELL_SIMUL, TOTAL_INTERACTIONS, CLUMP_DICT
    else:
        return INF_WELL_SIMUL, TOTAL_INTERACTIONS