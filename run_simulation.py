import pandas as pd
import yaml
import os
import sys
import json

sys.path.append(os.getcwd())
from simulation_utils.utils import PrepareData, PrepareParameters, PrepareParamList
from simulation_utils.one_strain_clump_simul import SimulateClump
from simulation_utils.one_strain_accdam_simul import SimulateAccDam
from simulation_utils.one_strain_comp_simul import SimulateComp
from simulation_utils.one_strain_null_simul import SimulateNull
from misc.misc_utils import ExtractParams
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
               '2022_10_27_TB_size_distribution']
dose_inf_df = pd.read_excel('data/Experimental_data_Ed_Josh.xlsx', sheet_name=SHEET_NAMES[sheet])

#print(dose_inf_df)
#====================================================================
''' Simulate infections '''
def RunMultiple(num_simulations, simul_name, config, dose_inf_df, sheet, save_path=None, scale=None, save_clump_info=False):
    assert simul_name in ['clump', 'comp', 'acc_dam', 'clump_comp', 'clump_acc_dam', 'null'], simul_name+" must be 'clump', 'comp', 'acc_dam', 'clump_acc_dam', 'clump_comp', or 'null'."
    #----------------------------------------------------------------
    ''' Get parameters for the simulation '''
    PARAM_DICT = PrepareParameters(config, simul_name, sheet, scale)
    #----------------------------------------------------------------
    ''' Prepare the experimental data '''
    GEN_WELL_DATA, GEN_CELL_DATA, INF_CELL_DATA, num_zeros = PrepareData(dose_inf_df, PARAM_DICT['scale'])
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
            filename = "ClumpSimul_"+SHEET_NAMES[sheet]+"_s="+str(PARAM_DICT['scale'])+"_vMax="+str(PARAM_DICT['vMax'])+"_"+scheme_short+"_"+dist_short+"_r="+str(PARAM_DICT['remove']) # Specify filename
        elif (simul_name == 'clump_comp'):
            filename = "ClumpCompSimul_"+SHEET_NAMES[sheet]+"_s="+str(PARAM_DICT['scale'])+"_vMax="+str(PARAM_DICT['vMax'])+"_k="+str(PARAM_DICT['kappa'])+"_"+scheme_short+"_"+dist_short+"_r="+str(PARAM_DICT['remove']) # Specify filename
        elif (simul_name == 'clump_acc_dam'):
            filename = "ClumpAccDamSimul_"+SHEET_NAMES[sheet]+"_s="+str(PARAM_DICT['scale'])+"_vMax="+str(PARAM_DICT['vMax'])+"_b="+str(PARAM_DICT['beta'])+"_"+scheme_short+"_"+dist_short+"_r="+str(PARAM_DICT['remove']) # Specify filename
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
            clump_size_df = pd.read_excel('data/Experimental_data_Ed_Josh.xlsx', sheet_name='2022_10_27_TB_size_distribution')
            #print(clump_size_df)
            if save_clump_info:
                dict_['GFP IU run='+str(simulation)], dict_['total_interactions run='+str(simulation)], CLUMP_DICT = SimulateClump(GEN_WELL_DATA, 
                                                                                                                                   PARAM_DICT, 
                                                                                                                                   cell_count, 
                                                                                                                                   clump_size_df,
                                                                                                                                   simul_name,
                                                                                                                                   save_clump_info=True)
                with open(os.path.join(os.path.join(save_path, 'clump_information'), filename+'_run='+str(simulation)+'_CLUMP'+'.json'), "w") as f:
                    json.dump(CLUMP_DICT, f)
            else:
                dict_['GFP IU run='+str(simulation)], dict_['total_interactions run='+str(simulation)] = SimulateClump(GEN_WELL_DATA, 
                                                                                                                       PARAM_DICT, 
                                                                                                                       cell_count,
                                                                                                                       clump_size_df,
                                                                                                                       simul_name,
                                                                                                                       save_clump_info=False)
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
            SimulResults.to_csv(os.path.join(save_path, filename+'_n='+str(num_simulations)+'VAR_REMOVAL.csv'), index=False)
            print(" - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ")
            print(" >> Saved at:", os.path.join(save_path, filename+'_n='+str(num_simulations)+'VAL_REMOVAL.csv'))
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