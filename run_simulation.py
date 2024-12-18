import sys
sys.dont_write_bytecode = True
import pandas as pd
import yaml
import os
import json

sys.path.append(os.getcwd())
from simulation_utils.utils import PrepareData, PrepareParameters, PrepareParamList
#from simulation_utils.one_strain_clump_simul import SimulateClump
from simulation_utils.one_strain_clump_simul_VVG import SimulateClump
from simulation_utils.one_strain_accdam_simul import SimulateAccDam
from simulation_utils.one_strain_comp_simul import SimulateComp
from simulation_utils.one_strain_null_simul import SimulateNull
from misc.misc_utils import ExtractParams, AverageDicts
#====================================================================
''' Select simulation type & dataset '''

with open('config.yml', 'r') as c:
    config = yaml.load(c, Loader=yaml.FullLoader)

sheet           = config['SIMULATION_PARAMETERS']['sheet']
scale           = config['SIMULATION_PARAMETERS']['scale']
simul_name      = config['SIMULATION_PARAMETERS']['simul_name']
num_simulations = config['SIMULATION_PARAMETERS']['num_simulations']

SHEET_NAMES = ['2021_10_05 TB_GFP_epithelial', '2020_07_02 ME_GFP_fibroblast', 
               '2020_05_29 TR_GFP_fibroblast', '2021_07_13 GFP_TB_fibroblast', 
               '2020_08_12 TB_GFP_fibroblast', '2020_09_14 TR_GFP_epithelial',
               '2021_08_13 ME_mC_epithelial', '2022_11_02_TB_GFP_fib', 
               'use_with_size_distribution', 'use_with_size_dist_interp', 
               '2022_10_27_TB_size_distribution', '2022_10_27_TB_size_dist_interp']

dose_inf_df = pd.read_excel('data/Experimental_data_Ed_Josh.xlsx', sheet_name=SHEET_NAMES[sheet])

#print(dose_inf_df)
#====================================================================
''' Simulate infections '''
def RunMultiple(num_simulations, simul_name, config, dose_inf_df, sheet, save_path=None, scale=None, save_clump_info=False):
    assert simul_name in ['clump', 'comp', 'acc_dam', 'clump_comp', 'clump_acc_dam', 'var_clump_diam', 'null'], " >> Check config.yml\n >> 'simul_name' must be 'clump', 'comp', 'acc_dam', 'clump_acc_dam', 'clump_comp', 'var_clump_diam', or 'null'. Got: " + simul_name
    assert sheet in range(10), " >> Check config.yml\n >> 'sheet' must be an integer between 0 and 9. Got: "+ str(sheet)
    assert config['SIMULATION_PARAMETERS']['remove'] in range(2), " >> Check config.yml\n >> 'remove' must be either 0 (False) or 1 (True). Got: "+ str(config['SIMULATION_PARAMETERS']['remove'])
    #----------------------------------------------------------------
    ''' Get parameters for the simulation '''
    PARAM_DICT = PrepareParameters(config, simul_name, sheet, scale)

    print(PARAM_DICT)
    #----------------------------------------------------------------
    ''' Prepare the experimental data '''
    GEN_WELL_DATA, GEN_CELL_DATA, INF_CELL_DATA, num_zeros = PrepareData(dose_inf_df, PARAM_DICT['scale'])
    
    
    print(" ?     ?      ?     ? GEN_WELL_DATA ?      ?       ?     ?")
    print("")
    print(GEN_WELL_DATA)
    print("")
    print(" ?     ?      ?     ? GEN_WELL_DATA ?      ?       ?     ?")

    #----------------------------------------------------------------
    ''' Get cell count'''
    DATA_DICT = ExtractParams(dose_inf_df)
    cell_count = int(DATA_DICT['cell_count'] / PARAM_DICT['scale'])
    #----------------------------------------------------------------
    ''' Data to save '''
    dict_ = {}
    dict_['GFP genomes (scaled)'] = GEN_WELL_DATA
    #----------------------------------------------------------------
    if ('clump' in simul_name):
        if (save_path==None and simul_name == 'clump'):
            save_path = os.path.join(os.getcwd(), 'simulation_results/clump')
        elif (save_path==None and simul_name == 'clump_comp'):
            save_path = os.path.join(os.getcwd(), 'simulation_results/clump_comp')
        elif (save_path==None and simul_name == 'clump_acc_dam'):
            save_path = os.path.join(os.getcwd(), 'simulation_results/clump_acc_dam')
        
        elif (save_path==None and simul_name == 'var_clump_diam'):
            save_path = os.path.join(os.getcwd(), 'simulation_results/var_clump_diam')
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        scheme_short = ''
        if (PARAM_DICT['scheme']=='linear'):
            scheme_short='lin'
        elif (PARAM_DICT['scheme']=='regular_polygon'):
            scheme_short='poly'
        elif (PARAM_DICT['scheme']=='sphere_packing'):
            scheme_short='sp'

        dist_short = ''
        if (PARAM_DICT['distribution']=='normal'):
            dist_short = 'norm'
        elif (PARAM_DICT['distribution']=='uniform'):
            dist_short = 'uni'
        elif (PARAM_DICT['distribution']=='fixed'):
            dist_short = 'fix'
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        filename = ''
        if (simul_name == 'clump'):
            filename = "ClumpSimulVVG_"+SHEET_NAMES[sheet]+"_s="+str(PARAM_DICT['scale'])+"_vMax="+str(PARAM_DICT['vMax'])+"_"+scheme_short+"_"+dist_short+"_r="+str(PARAM_DICT['remove']) # Specify filename
        elif (simul_name == 'clump_comp'):
            filename = "ClumpCompSimulVVG_"+SHEET_NAMES[sheet]+"_s="+str(PARAM_DICT['scale'])+"_vMax="+str(PARAM_DICT['vMax'])+"_k="+str(PARAM_DICT['kappa'])+"_"+scheme_short+"_"+dist_short+"_r="+str(PARAM_DICT['remove']) # Specify filename
        elif (simul_name == 'clump_acc_dam'):
            filename = "ClumpAccDamSimulVVG_"+SHEET_NAMES[sheet]+"_s="+str(PARAM_DICT['scale'])+"_vMax="+str(PARAM_DICT['vMax'])+"_b="+str(PARAM_DICT['beta'])+"_"+scheme_short+"_"+dist_short+"_r="+str(PARAM_DICT['remove']) # Specify filename

        elif (simul_name == 'var_clump_diam'):
            func_short = ''
            if (PARAM_DICT['diameter_func']=='constant'):
                func_short = 'c'
            elif (PARAM_DICT['diameter_func']=='linear'):
                func_short = 'l'
            elif (PARAM_DICT['diameter_func']=='exponential'):
                func_short = 'e'
            filename = "VarClumpDiamSimul_"+SHEET_NAMES[sheet]+"_s="+str(PARAM_DICT['scale'])+"_vMax="+str(PARAM_DICT['vMax'])+"_f="+func_short+"_"+scheme_short+"_"+dist_short+"_r="+str(PARAM_DICT['remove']) # Specify filename
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        if not 'Parameters' in dict_.keys():
                PARAM_LIST_STR = PrepareParamList(GEN_CELL_DATA, simul_name, cell_count, PARAM_DICT)
                dict_['Parameters'] = PARAM_LIST_STR
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ''' Get clump_size_params_dict if necessary '''
        CLUMP_SIZE_PARAMS_DICT = None # Put in scope
        if (((sheet == 8) or (sheet == 9)) and not (simul_name == 'var_clump_diam')):
            CLUMP_SIZE_PARAMS_DICT = config['CLUMP_PARAMETERS']['clump_dist_params_dict']
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
            clump_size_df = 69
            if (SHEET_NAMES[sheet] == 'use_with_size_distribution'):
                clump_size_df = pd.read_excel('data/Experimental_data_Ed_Josh.xlsx', sheet_name='2022_10_27_TB_size_distribution')
            elif (SHEET_NAMES[sheet] == 'use_with_size_dist_interp'):
                clump_size_df = pd.read_excel('data/Experimental_data_Ed_Josh.xlsx', sheet_name='2022_10_27_TB_size_dist_interp')
            #print(clump_size_df)
            if save_clump_info:
                dict_['GFP IU run='+str(simulation)], dict_['total_interactions run='+str(simulation)], CLUMP_DICT = SimulateClump(GEN_WELL_DATA, 
                                                                                                                                   PARAM_DICT, 
                                                                                                                                   cell_count, 
                                                                                                                                   clump_size_df,
                                                                                                                                   simul_name,
                                                                                                                                   SHEET_NAMES[sheet],
                                                                                                                                   save_clump_info=True,
                                                                                                                                   CLUMP_DIST_PARAMS_DICT=CLUMP_SIZE_PARAMS_DICT)
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
                dict_['GFP IU run='+str(simulation)], dict_['total_interactions run='+str(simulation)] = SimulateClump(GEN_WELL_DATA, 
                                                                                                                       PARAM_DICT, 
                                                                                                                       cell_count,
                                                                                                                       clump_size_df,
                                                                                                                       simul_name,
                                                                                                                       SHEET_NAMES[sheet],
                                                                                                                       save_clump_info=False,
                                                                                                                       CLUMP_DIST_PARAMS_DICT=CLUMP_SIZE_PARAMS_DICT)
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            PARAM_LIST_STR = PrepareParamList(GEN_CELL_DATA, simul_name, cell_count, PARAM_DICT)
            dict_['Parameters'] = PARAM_LIST_STR
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            SimulResults = pd.DataFrame.from_dict(dict_)
            SimulResults.to_csv(os.path.join(save_path, filename+'_n='+str(num_simulations)+'.csv'), index=False)
            print(" - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ")
            print(" >> Saved at:", os.path.join(save_path, filename+'_n='+str(num_simulations)+'.csv'))
    #----------------------------------------------------------------
    elif (simul_name == 'comp'):
        if (save_path==None):
            save_path = os.path.join(os.getcwd(), 'simulation_results/comp')
        filename = "CompSimul_"+SHEET_NAMES[sheet]+"_s="+str(PARAM_DICT['scale'])+"_vMax="+str(PARAM_DICT['vMax'])+"_k="+str(PARAM_DICT['kappa'])+"_r="+str(PARAM_DICT['remove']) # Specify filename
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        if not 'Parameters' in dict_.keys():
                PARAM_LIST_STR = PrepareParamList(GEN_CELL_DATA, simul_name, cell_count, PARAM_DICT)
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
            dict_['GFP IU run='+str(simulation)], dict_['total_interactions run='+str(simulation)] = SimulateComp(GEN_WELL_DATA, 
                                                                                                                  PARAM_DICT, 
                                                                                                                  cell_count)
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            SimulResults = pd.DataFrame.from_dict(dict_)
            SimulResults.to_csv(os.path.join(save_path, filename+'_n='+str(num_simulations)+'.csv'), index=False)
            print(" - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ")
            print(" >> Saved at:", os.path.join(save_path, filename+'_n='+str(num_simulations)+'.csv'))
    #----------------------------------------------------------------
    elif (simul_name == 'acc_dam'):
        print("AHHHHHHHHH     ACC DAM")
        if (save_path==None):
            save_path = os.path.join(os.getcwd(), 'simulation_results/acc_dam')
        filename = "AccDamSimul_"+SHEET_NAMES[sheet]+"_s="+str(PARAM_DICT['scale'])+"_vMax="+str(PARAM_DICT['vMax'])+"_b="+str(PARAM_DICT['beta'])+"_r="+str(PARAM_DICT['remove']) # Specify filename
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        if not 'Parameters' in dict_.keys():
                PARAM_LIST_STR = PrepareParamList(GEN_CELL_DATA, simul_name, cell_count, PARAM_DICT)
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
            dict_['GFP IU run='+str(simulation)], dict_['total_interactions run='+str(simulation)] = SimulateAccDam(GEN_WELL_DATA, 
                                                                                                                    PARAM_DICT, 
                                                                                                                    cell_count)
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            SimulResults = pd.DataFrame.from_dict(dict_)
            SimulResults.to_csv(os.path.join(save_path, filename+'_n='+str(num_simulations)+'.csv'), index=False)
            print(" - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ")
            print(" >> Saved at:", os.path.join(save_path, filename+'_n='+str(num_simulations)+'.csv'))
    #----------------------------------------------------------------    
    elif (simul_name == 'null'):
        if (save_path==None):
            save_path = os.path.join(os.getcwd(), 'simulation_results/null')
        filename = "NullSimul_"+SHEET_NAMES[sheet]+"_s="+str(PARAM_DICT['scale'])+"_vMax="+str(PARAM_DICT['vMax'])+"_b="+str(PARAM_DICT['b'])+"_r="+str(PARAM_DICT['remove']) # Specify filename
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        if not 'Parameters' in dict_.keys():
                PARAM_LIST_STR = PrepareParamList(GEN_CELL_DATA, simul_name, cell_count, PARAM_DICT)
                dict_['Parameters'] = PARAM_LIST_STR
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        for simulation in range(num_simulations):
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            print("<>><<>><<>><<>><<>><<>><<>>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>")
            print("<>><<>><<>><<>><<>><<>><<>>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>")
            print("<>><<>><<>><<>><<>><<>><<>>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>")
            print("<>><<>><<>><<>><<>><<>><<>>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>")
            print(" - - - - - - - - - - - - - - Simulation", simulation, " - - - - - - - - - - - - - - ")
            print("<>><<>><<>><<>><<>><<>><<>>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>")
            print("<>><<>><<>><<>><<>><<>><<>>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>")
            print("<>><<>><<>><<>><<>><<>><<>>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>")
            print("<>><<>><<>><<>><<>><<>><<>>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>")
            dict_['GFP IU run='+str(simulation)], dict_['total_interactions run='+str(simulation)] = SimulateNull(GEN_WELL_DATA, 
                                                                                                                  PARAM_DICT, 
                                                                                                                  cell_count)
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            SimulResults = pd.DataFrame.from_dict(dict_)
            SimulResults.to_csv(os.path.join(save_path, filename+'_n='+str(num_simulations)+'.csv'), index=False)
            print(" - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ")
            print(" >> Saved at:", os.path.join(save_path, filename+'_n='+str(num_simulations)+'.csv'))
            
    print("=======================================================================")
#====================================================================
RunMultiple(num_simulations, simul_name, config, dose_inf_df, sheet, save_clump_info=True)