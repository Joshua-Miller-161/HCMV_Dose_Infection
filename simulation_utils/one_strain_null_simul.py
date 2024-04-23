import numpy as np
import os
import sys

sys.path.append(os.getcwd())
from simulation_utils.utils import Innoculation, Cell, Virion
#====================================================================
def SimulateNull(GEN_WELL_DATA, PARAM_DICT, cell_count):
    INF_WELL_SIMUL = np.empty(GEN_WELL_DATA.shape[0], int)
    #----------------------------------------------------------------
    for init in range(len(GEN_WELL_DATA)): # Iterate thorugh each experiment
        print("=======================================================================")
        print("GEN_WELL_DATA = ", GEN_WELL_DATA[init], "| Experiment ", init)
        print(" - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ")
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ''' Set up virus and cell populations '''
        GFP_POOL = []
        CELL_POOL = []
        num_infG = 0
        total = 0
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        for j in range(GEN_WELL_DATA[init]):
            i = np.random.normal(PARAM_DICT['i_mean'], PARAM_DICT['i_stdev'])
            GFP_POOL.append(Virion(i, True))
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        for k in range(cell_count):
            r = np.random.normal(PARAM_DICT['r_mean'], PARAM_DICT['r_stdev'])
            CELL_POOL.append(Cell(False, 0, r, 0))
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ''' ACTUAL SIMULATION PART '''
            if (len(GFP_POOL) > 0):
                num_interactions = np.random.poisson((len(GFP_POOL) / PARAM_DICT['vMax']) + PARAM_DICT['b']) # Calculate the number of virions that will visit the cell
                if (num_interactions >= len(GFP_POOL)):
                    num_interactions = len(GFP_POOL)
                
                VIRION_ATTACKERS = []
                VIRION_ATTACKERS_IDX = []
                TO_REMOVE_IDX = []
                #----------------------------------------------------
                ''' Calculate the maximum infectivity of the virions coming to the cell '''

                if (num_interactions < len(GFP_POOL)):
                    idx_created = 0
                    while(idx_created <= num_interactions):
                        index_g = np.random.randint(0, len(GFP_POOL))
                        if not (index_g in VIRION_ATTACKERS_IDX):
                            VIRION_ATTACKERS_IDX.append(index_g)
                            idx_created += 1

                else:
                    VIRION_ATTACKERS_IDX = np.arange(0, num_interactions)
                
                for idx in VIRION_ATTACKERS_IDX:
                    VIRION_ATTACKERS.append(GFP_POOL[idx])
                #----------------------------------------------------                                  
                ''' Now have virions try to infect the cell with the buffed infectivity '''
                is_successful = False
                for b in range(num_interactions):
                    is_successful = Innoculation(VIRION_ATTACKERS[b], CELL_POOL[k], gamma=PARAM_DICT['gamma'])
                    total += 1
                    if is_successful:
                        TO_REMOVE_IDX.append(VIRION_ATTACKERS_IDX[b])
                #----------------------------------------------------
                ''' Remove virions which have successfully infected cells '''
                TO_REMOVE_IDX = list(sorted(set(TO_REMOVE_IDX), reverse=True))

                for c in range(len(TO_REMOVE_IDX)):
                    if (len(TO_REMOVE_IDX) >= len(GFP_POOL)):
                        #print('BREAKING idx_len=', len(TO_REMOVE_IDX), ", virions left=",len(GFP_POOL))
                        break
                    else:
                        #print("c=",c,", idx=", TO_REMOVE_IDX[c], ', len=', len(TO_REMOVE_IDX), ", virions left=", len(GFP_POOL))
                        GFP_POOL.remove(GFP_POOL[TO_REMOVE_IDX[c]])
                #print('DONE + + + + + + ')
            else:
                print("[][][][] OUT OF VIRIONS [][][][]")
                break
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if (CELL_POOL[k].infG == True):
                num_infG = num_infG + 1
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        print("num infG =", num_infG, "| cell_count =", cell_count, "| interactions / cell =", round(total / cell_count, 5))
        
        INF_WELL_SIMUL[init] = num_infG
    #----------------------------------------------------------------
    return INF_WELL_SIMUL