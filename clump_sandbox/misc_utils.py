import sys
sys.dont_write_bytecode = True
import os
#====================================================================
def PrepareParameters(config, simul_name):
    assert simul_name in ['clump', 'clump_comp', 'clump_acc_dam'], simul_name+" must be 'clump', 'clump_comp', 'clump_acc_dam'."
    #----------------------------------------------------------------
    PARAM_DICT = {}
    #----------------------------------------------------------------
    PARAM_DICT['simul_name']      = config['SIMULATION_PARAMETERS']['simul_name']
    PARAM_DICT['scale']           = float(config['SIMULATION_PARAMETERS']['scale'])
    PARAM_DICT['cell_count']      = float(config['SIMULATION_PARAMETERS']['cell_count']) / float(config['SIMULATION_PARAMETERS']['scale'])
    PARAM_DICT['num_simulations'] = int(config['SIMULATION_PARAMETERS']['num_simulations'])
    PARAM_DICT['num_fits']        = int(config['SIMULATION_PARAMETERS']['num_fits'])
    PARAM_DICT['remove']          = int(config['SIMULATION_PARAMETERS']['remove'])
    PARAM_DICT['i_mean']  = float(config['GFP_VIRUS_PARAMETERS']['mu_inf']) #0
    PARAM_DICT['i_stdev'] = float(config['GFP_VIRUS_PARAMETERS']['std_inf']) #1
    PARAM_DICT['r_mean']  = float(config['CELL_PARAMETERS']['mu_res'])
    PARAM_DICT['r_stdev'] = float(config['CELL_PARAMETERS']['std_res'])
    PARAM_DICT['gamma']   = float(config['SIMULATION_PARAMETERS']['gamma'])
    #----------------------------------------------------------------
    if ('clump' in simul_name):
        PARAM_DICT['distribution'] = str(config['CLUMP_PARAMETERS']['distribution'])
        
        vMax   = config['CLUMP_PARAMETERS']['vMax']
        vMax_c = config['CLUMP_PARAMETERS']['vMax_c']

        PARAM_DICT['vMax']       = float(vMax) / float(PARAM_DICT['scale'])
        PARAM_DICT['vMax_c']     = float(vMax_c) / float(PARAM_DICT['scale'])
        PARAM_DICT['fixed_mean'] = bool(config['CLUMP_PARAMETERS']['fixed_mean'])
        PARAM_DICT['mean']       = float(config['CLUMP_PARAMETERS']['mean'])
        PARAM_DICT['f_1']        = float(config['CLUMP_PARAMETERS']['f_1'])

    if ('comp' in simul_name):
        KAPPA = config['COMPENSATION_PARAMETERS']['kappa']
        PARAM_DICT['kappa'] = float(KAPPA)
        if not ('clump' in simul_name):
            PARAM_DICT['vMax']  = float(vMax) / PARAM_DICT['scale']

    if ('acc_dam' in simul_name):
        beta = config['ACCRUED_DAMAGE_PARAMETERS']['beta']
        PARAM_DICT['beta'] = float(beta)
        if not ('clump' in simul_name):
            PARAM_DICT['vMax'] = float(vMax) / PARAM_DICT['scale']
  
    return PARAM_DICT

#====================================================================
def PrepareParamList(GENOMES_WELL, simul_name, PARAM_DICT):
    assert simul_name in ['clump', 'comp', 'acc_dam', 'clump_comp', 'clump_acc_dam', 'var_clump_diam', 'null'], simul_name+" must be 'clump', 'comp', 'acc_dam', 'clump_comp', 'clump_acc_dam', 'var_clump_diam', or 'null'."
    #----------------------------------------------------------------
    PARAM_LIST_STR = ['simul_name='+str(simul_name), 'num_simulations='+str(PARAM_DICT['num_simulations']),
                      'scale='+str(PARAM_DICT['scale']), 'cell_count='+str(PARAM_DICT['cell_count']), 
                      'muG='+str(PARAM_DICT['i_mean']), 'sigmaG='+str(PARAM_DICT['i_stdev']),
                      'muR='+str(PARAM_DICT['r_mean']), 'sigmaR='+str(PARAM_DICT['r_stdev']), 
                      'gamma='+str(PARAM_DICT['gamma']), 'vMax='+str(PARAM_DICT['vMax']), 'remove='+str(PARAM_DICT['remove'])]
    #----------------------------------------------------------------
    if ('clump' in simul_name):
        PARAM_LIST_STR.append('distribution='+str(PARAM_DICT['distribution']))
        PARAM_LIST_STR.append('fixed_mean='+str(PARAM_DICT['fixed_mean']))

        if (PARAM_DICT['fixed_mean'] == True):
            PARAM_LIST_STR.append('mean='+str(PARAM_DICT['mean']))
        else:
            PARAM_LIST_STR.append('vMax_c='+str(PARAM_DICT['vMax_c']))

        if (PARAM_DICT['distribution'] == '1inf-geo'):
            PARAM_LIST_STR.append('f_1='+str(PARAM_DICT['f_1']))

    if ('comp' in simul_name):
        PARAM_LIST_STR.append('kappa='+str(PARAM_DICT['kappa']))

    if ('acc_dam' in simul_name):
        PARAM_LIST_STR.append('beta='+str(PARAM_DICT['beta']))

    #----------------------------------------------------------------
    for i in range(GENOMES_WELL.shape[0] - len(PARAM_LIST_STR)):
        PARAM_LIST_STR.append('FILLER')
    #----------------------------------------------------------------
    return PARAM_LIST_STR
#====================================================================
def MakeFilename(PARAM_DICT):
    filename = ''

    if (PARAM_DICT['distribution'] == 'geometric'):
        dist_short = 'geo'
    elif (PARAM_DICT['distribution'] == '1inf-geo'):
        dist_short = '1geo'
    elif (PARAM_DICT['distribution'] == 'poisson'):
        dist_short = 'poi'

    if ('clump' in PARAM_DICT['simul_name']):
        
        param_text = 's='+str(PARAM_DICT['scale'])+'_vMax='+str(PARAM_DICT['vMax'])+'_f='+dist_short
        
        if (PARAM_DICT['fixed_mean'] == True or PARAM_DICT['fixed_mean'] == 'True' or PARAM_DICT['fixed_mean'] == 1 or PARAM_DICT['fixed_mean'] == '1'):
            param_text += '_fix=1_mean='+str(PARAM_DICT['mean'])
        else:
            param_text += '_vMax_c='+str(PARAM_DICT['vMax_c'])
        
        if (PARAM_DICT['distribution'] == '1inf-geo'):
            param_text += '_f_1='+str(PARAM_DICT['f_1'])

        if (PARAM_DICT['simul_name'] == 'clump'):
            filename = "ClumpSimulVVG_"+param_text+"_r="+str(PARAM_DICT['remove'])+"_n="+str(PARAM_DICT['num_simulations']) # Specify filename
        
        elif (PARAM_DICT['simul_name'] == 'clump_comp'):
            filename = "ClumpCompSimulVVG_"+param_text+"_k="+PARAM_DICT['kappa']+"_r="+str(PARAM_DICT['remove'])+"_n="+str(PARAM_DICT['num_simulations']) # Specify filename
        
        elif (PARAM_DICT['simul_name'] == 'clump_acc_dam'):
            filename = "ClumpAccDamSimulVVG_"+param_text+"_b="+str(PARAM_DICT['beta'])+"_r="+str(PARAM_DICT['remove'])+"_n="+str(PARAM_DICT['num_simulations']) # Specify filename

    return filename
#====================================================================