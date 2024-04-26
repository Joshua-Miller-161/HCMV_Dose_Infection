import numpy as np
import os
import sys

sys.path.append(os.getcwd())
from simulation_utils.utils import Innoculation, Cell, Virion
#================================================================
def SimulateAccDam(GEN_WELL_DATA, PARAM_DICT, cell_count):
    INF_WELL_SIMUL = np.empty(GEN_WELL_DATA.shape[0], int)
    #----------------------------------------------------------------
    for init in range(len(GEN_WELL_DATA)): # Orig len(GEN_WELL_DATA) Iterate thorugh each experiment
        print("=======================================================================")
        print(PARAM_DICT['simul_name'], "| GEN_WELL_DATA[", init, "] =", GEN_WELL_DATA[init], "| gamma =", PARAM_DICT['gamma'], "| beta =", PARAM_DICT['beta'], "| vMax =", PARAM_DICT['vMax'], "| scale =", PARAM_DICT['scale'])
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
                #----------------------------------------------------
                ''' Determine which virions interact with the cell '''
                num_interactions = np.random.poisson(GEN_WELL_DATA[init] / PARAM_DICT['vMax']) # Works with GEN_WELL_DATA[init]. Calculate the number of virions that will visit the cell
                
                VIRION_ATTACKERS = []
                VIRION_ATTACKERS_IDX = []
                TO_REMOVE_IDX = []
                
                if (num_interactions > 0):

                    if (num_interactions < len(GFP_POOL)):
                        idx_created = 0
                        while(idx_created < num_interactions):
                            index_g = np.random.randint(0, len(GFP_POOL))
                            if not (index_g in VIRION_ATTACKERS_IDX):
                                VIRION_ATTACKERS_IDX.append(index_g)
                                idx_created += 1

                    else:
                        #print("num_interactions >= len(GFP_POOL)", "num_interactions:", num_interactions, ", len(GFP):", len(GFP_POOL))
                        VIRION_ATTACKERS_IDX = np.arange(0, len(GFP_POOL))
                    
                    for idx in VIRION_ATTACKERS_IDX:
                        VIRION_ATTACKERS.append(GFP_POOL[idx])
                    #----------------------------------------------------
                    ''' Have virions interact with cells '''
                    for a in range(len(VIRION_ATTACKERS)):
                        is_successful = Innoculation(VIRION_ATTACKERS[a], 
                                                     CELL_POOL[k], 
                                                     gamma=PARAM_DICT['gamma'], 
                                                     acc_dam=True, 
                                                     beta=PARAM_DICT['beta'])
                        total += 1
                        if is_successful:
                            TO_REMOVE_IDX.append(VIRION_ATTACKERS_IDX[a])
                    #----------------------------------------------------
                    if (PARAM_DICT['remove'] == True):
                        # Remove virions which have successfully the cell
                        TO_REMOVE_IDX = list(sorted(set(TO_REMOVE_IDX), reverse=True))

                        for c in range(len(TO_REMOVE_IDX)):
                            if (len(TO_REMOVE_IDX) >= len(GFP_POOL)):
                                #print('BREAKING idx_len=', len(TO_REMOVE_IDX), ", virions left=",len(GFP_POOL))
                                break
                            else:
                                #print("c=",c,", idx=", TO_REMOVE_IDX[c], ', len=', len(TO_REMOVE_IDX), ", virions left=", len(GFP_POOL))
                                #GFP_POOL.remove(GFP_POOL[TO_REMOVE_IDX[c]]) Works but slow
                                del GFP_POOL[TO_REMOVE_IDX[c]]

                if (k == 0 or k == 500 or k == 1500 or k == 2299):
                    print("cell=", k, ", poisson:", num_interactions, ",len(VIR_ATT):", len(VIRION_ATTACKERS), ", GFP_GENOMES[",init,"]:", GEN_WELL_DATA[init], ", len(GFP):", len(GFP_POOL), ", total:", total, ", Cell.inter:", CELL_POOL[k].inter)

            else:
                print("[][][][] OUT OF VIRIONS [][][][]")
                break
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if (CELL_POOL[k].infG == True):
                num_infG = num_infG + 1
        print("num infG = ", num_infG, " | cell_count = ", cell_count, " | GEN_WELL_DATA = ", GEN_WELL_DATA[init], " | avg. inter/cell = ", round(total / cell_count, 5))
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        INF_WELL_SIMUL[init] = num_infG
    #----------------------------------------------------------------
    return INF_WELL_SIMUL