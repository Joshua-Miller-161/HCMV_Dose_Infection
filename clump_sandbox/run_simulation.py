import sys
sys.dont_write_bytecode = True
import pandas as pd
import yaml
import os
import json
import numpy as np

sys.path.append(os.getcwd())
from misc_utils import PrepareParameters, PrepareParamList, MakeFilename
from one_strain_clump_simul_VVG import SimulateClump
from misc.misc_utils import ExtractParams, AverageDicts
#====================================================================
''' Get some of the parameters '''

with open('clump_sandbox/sandbox_config.yml', 'r') as c:
    config = yaml.load(c, Loader=yaml.FullLoader)

scale           = config['SIMULATION_PARAMETERS']['scale']
simul_name      = config['SIMULATION_PARAMETERS']['simul_name']
num_simulations = config['SIMULATION_PARAMETERS']['num_simulations']
#====================================================================
''' Simulate infections '''
def RunMultiple(num_simulations, simul_name, config, save_path=None, save_clump_info=False):
    assert simul_name in ['clump', 'clump_comp', 'clump_acc_dam'], " >> Check sandbox_config.yml\n >> 'simul_name' must be 'clump', 'clump_acc_dam', 'clump_comp'. Got: " + simul_name
    assert config['SIMULATION_PARAMETERS']['remove'] in range(2), " >> Check sandbox_config.yml\n >> 'remove' must be either 0 (False) or 1 (True). Got: "+ str(config['SIMULATION_PARAMETERS']['remove'])
    #----------------------------------------------------------------
    ''' Get parameters for the simulation '''
    PARAM_DICT = PrepareParameters(config, simul_name)

    #print(PARAM_DICT)
    #----------------------------------------------------------------
    ''' Genomes/well '''
    
    # GENOMES_WELL = np.asarray([10,20,30,40,50,60,70,80,90,
    #                            100,200,300,400,500,600,700,800,900,
    #                            1000,2000,3000,4000,5000,6000,7000,8000,9000,
    #                            10000,20000,30000,40000,50000,60000,70000,80000,90000]) / PARAM_DICT['scale']

    GENOMES_WELL = np.asarray([10,20,30,40,50,60,70,80,90,
                               100,200,300,400,500,600,700,800,900]) / PARAM_DICT['scale']
    #----------------------------------------------------------------
    ''' Cell count '''
    
    cell_count = int(PARAM_DICT['cell_count'])
    #----------------------------------------------------------------
    ''' Data to save '''
    dict_ = {}
    dict_['GFP genomes (scaled)'] = GENOMES_WELL
    #----------------------------------------------------------------
    if ('clump' in simul_name):
        if (save_path==None and simul_name == 'clump'):
            save_path = os.path.join(os.getcwd(), 'simulation_results/clump')
        elif (save_path==None and simul_name == 'clump_comp'):
            save_path = os.path.join(os.getcwd(), 'simulation_results/clump_comp')
        elif (save_path==None and simul_name == 'clump_acc_dam'):
            save_path = os.path.join(os.getcwd(), 'simulation_results/clump_acc_dam')

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        
        filename = MakeFilename(PARAM_DICT)
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        if not 'Parameters' in dict_.keys():
            PARAM_LIST_STR = PrepareParamList(GENOMES_WELL, simul_name, PARAM_DICT)
            dict_['Parameters'] = PARAM_LIST_STR
        
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        for simulation in range(num_simulations):
            print("<>><<>><<>><<>><<>><<>><<>>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>")
            print("<>><<>><<>><<>><<>><<>><<>>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>")
            print("<>><<>><<>><<>><<>><<>><<>>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>")
            print("<>><<>><<>><<>><<>><<>><<>>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>")
            print(" - - - - - - - - - - - - - - Simulation", simulation, " - - - - - - - - - - - - - - ")
            print("<>><<>><<>><<>><<>><<>><<>>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>")
            print("<>><<>><<>><<>><<>><<>><<>>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>")
            print("<>><<>><<>><<>><<>><<>><<>>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>")
            print("<>><<>><<>><<>><<>><<>><<>>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>")

            if save_clump_info:
                dict_['GFP IU run='+str(simulation)], dict_['total_interactions run='+str(simulation)], CLUMP_DICT = SimulateClump(GENOMES_WELL, 
                                                                                                                                   PARAM_DICT, 
                                                                                                                                   simul_name,
                                                                                                                                   save_clump_info=True)
                # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
                if (simulation == 0):
                    with open(os.path.join(os.path.join(save_path, 'clump_information'), filename+'_CLUMP'+'.json'), "w") as f:
                        json.dump(CLUMP_DICT, f)
                
                elif (simulation > 0):
                    CLUMP_DICT_ORIG = {}
                    with open(os.path.join(os.path.join(save_path, 'clump_information'), filename+'_CLUMP'+'.json'), "r") as f:
                        CLUMP_DICT_ORIG = json.load(f)
                    
                    CLUMP_DICT_AVG = AverageDicts(CLUMP_DICT, CLUMP_DICT_ORIG, simulation)
                    del(CLUMP_DICT)
                    del(CLUMP_DICT_ORIG)
                    with open(os.path.join(os.path.join(save_path, 'clump_information'), filename+'_CLUMP'+'.json'), "w") as f:
                        json.dump(CLUMP_DICT_AVG, f)
                # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            else:
                dict_['GFP IU run='+str(simulation)], dict_['total_interactions run='+str(simulation)] = SimulateClump(GENOMES_WELL, 
                                                                                                                       PARAM_DICT,
                                                                                                                       simul_name,
                                                                                                                       save_clump_info=False)
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            PARAM_LIST_STR = PrepareParamList(GENOMES_WELL, simul_name, PARAM_DICT)
            dict_['Parameters'] = PARAM_LIST_STR
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            SimulResults = pd.DataFrame.from_dict(dict_)
            SimulResults.to_csv(os.path.join(save_path, filename+'.csv'), index=False)
            print(" - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ")
            print(" >> Saved at:", os.path.join(save_path, filename+'.csv'))
    print("=======================================================================")
#====================================================================
RunMultiple(num_simulations, simul_name, config, save_clump_info=True)