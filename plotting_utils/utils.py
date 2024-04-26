import numpy as np
import pandas as pd
from itertools import chain
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
def negLogLikeModel(params, xData, yData):
    new_xData = []
    new_yData = []

    for i in range(len(yData)):    
        if ((yData[i] != 0) and (xData[i] != 0)):
            new_xData.append(xData[i])
            new_yData.append(yData[i])

    model_result = model(new_xData, params)

    nll = 0
    epsilon = 10**-10
    for i in range(len(new_yData)):

        if (model_result[i] < epsilon):
            model_result[i] = epsilon

        if (model_result[i] > 1 - epsilon):
            model_result[i] = 1 - epsilon

        nll += new_yData[i] * np.log(model_result[i]) + (1 - new_yData[i]) * np.log(1 - model_result[i])

    return -nll
#====================================================================
def CombineSameGenWell(GEN_WELL_SIMUL, INF_WELL_SIMULS):
    KEYS = list(set(GEN_WELL_SIMUL))
    INF_WELL_SIMULS_DICT = {float(key):[] for key in KEYS}

    for i in range(np.shape(GEN_WELL_SIMUL)[0]):
        INF_WELL_SIMULS_DICT[GEN_WELL_SIMUL[i]].append(list(INF_WELL_SIMULS[:, i]))

    for key in INF_WELL_SIMULS_DICT.keys():
        INF_WELL_SIMULS_DICT[key] = np.asarray(list(chain.from_iterable(INF_WELL_SIMULS_DICT[key])))

    GEN_WELL_SIMUL_U  = np.asarray(list(INF_WELL_SIMULS_DICT.keys()))
    INF_WELL_SIMULS_U = np.asarray(list(INF_WELL_SIMULS_DICT.values()))

    IDX = np.argsort(GEN_WELL_SIMUL_U)
    GEN_WELL_SIMUL_U  = GEN_WELL_SIMUL_U[IDX]
    INF_WELL_SIMULS_U = INF_WELL_SIMULS_U[IDX, :]

    del(INF_WELL_SIMULS_DICT)

    print(type(GEN_WELL_SIMUL_U), GEN_WELL_SIMUL_U)
    return GEN_WELL_SIMUL_U, INF_WELL_SIMULS_U
#====================================================================
def MakeDataPretty(df, x_col_name, y_col_name, num_simulations):
    GEN_WELL  = np.asarray(df.loc[:, x_col_name])
    DATA = np.empty((num_simulations, df.shape[0]), float)

    for i in range(num_simulations):
        DATA[i, :] = df.loc[:, y_col_name+ ' run='+str(i)]
    
    GEN_WELL_U, DATA_U = CombineSameGenWell(GEN_WELL, DATA)
    
    STDEVS_U = np.std(DATA_U, axis=1)
    CIS_U = (1.96 * np.ones(np.shape(DATA_U)[0])) * STDEVS_U / (np.sqrt(np.shape(DATA_U)[1]) * np.ones(np.shape(DATA_U)[0]))
    MAXS_U = np.amax(DATA_U, axis=1)
    MINS_U = np.amin(DATA_U, axis=1)

    return GEN_WELL_U, DATA_U, CIS_U, MINS_U, MAXS_U
#====================================================================
def MakeFilename(PARAM_DICT_SIMUL, sheet_name):
    filename = ''
    if ('clump' in PARAM_DICT_SIMUL['simul_name']):
        scheme_short = ''
        if (PARAM_DICT_SIMUL['scheme']=='linear'):
            scheme_short='lin'
        elif (PARAM_DICT_SIMUL['scheme']=='regular_polygon'):
            scheme_short='poly'
        elif (PARAM_DICT_SIMUL['scheme']=='sphere_packing'):
            scheme_short='sp'

        dist_short = ''
        if (PARAM_DICT_SIMUL['distribution']=='normal'):
            dist_short = 'norm'
        elif (PARAM_DICT_SIMUL['distribution']=='uniform'):
            dist_short = 'uni'
        elif (PARAM_DICT_SIMUL['distribution']=='fixed'):
            dist_short = 'fix'
        
        if (PARAM_DICT_SIMUL['simul_name'] == 'clump'):
            filename = "ClumpSimul_"+sheet_name+"_s="+str(PARAM_DICT_SIMUL['scale'])+"_vMax="+str(PARAM_DICT_SIMUL['vMax'])+"_"+scheme_short+"_"+dist_short+"_r="+str(PARAM_DICT_SIMUL['remove'])+"_n="+str(PARAM_DICT_SIMUL['num_simulations']) # Specify filename
        elif (PARAM_DICT_SIMUL['simul_name'] == 'clump_comp'):
            filename = "ClumpCompSimul_"+sheet_name+"_s="+str(PARAM_DICT_SIMUL['scale'])+"_vMax="+str(PARAM_DICT_SIMUL['vMax'])+"_k="+str(PARAM_DICT_SIMUL['kappa'])+"_"+scheme_short+"_"+dist_short+"_r="+str(PARAM_DICT_SIMUL['remove'])+"_n="+str(PARAM_DICT_SIMUL['num_simulations']) # Specify filename
        elif (PARAM_DICT_SIMUL['simul_name'] == 'clump_acc_dam'):
            filename = "ClumpAccDamSimul_"+sheet_name+"_s="+str(PARAM_DICT_SIMUL['scale'])+"_vMax="+str(PARAM_DICT_SIMUL['vMax'])+"_b="+str(PARAM_DICT_SIMUL['beta'])+"_"+scheme_short+"_"+dist_short+"_r="+str(PARAM_DICT_SIMUL['remove'])+"_n="+str(PARAM_DICT_SIMUL['num_simulations']) # Specify filename

    elif (PARAM_DICT_SIMUL['simul_name'] == 'acc_dam'):
        filename = "AccDamSimul_"+sheet_name+"_s="+str(PARAM_DICT_SIMUL['scale'])+"_vMax="+str(PARAM_DICT_SIMUL['vMax'])+"_b="+str(PARAM_DICT_SIMUL['beta'])+"_r="+str(PARAM_DICT_SIMUL['remove'])+"_n="+str(PARAM_DICT_SIMUL['num_simulations'])

    elif (PARAM_DICT_SIMUL['simul_name'] == 'comp'):
        filename = "CompSimul_"+sheet_name+"_s="+str(PARAM_DICT_SIMUL['scale'])+"_vMax="+str(PARAM_DICT_SIMUL['vMax'])+"_k="+str(PARAM_DICT_SIMUL['kappa'])+"_r="+str(PARAM_DICT_SIMUL['remove'])+"_n="+str(PARAM_DICT_SIMUL['num_simulations'])

    elif (PARAM_DICT_SIMUL['simul_name'] == 'null'):
        filename = "NullSimul_"+sheet_name+"_s="+str(PARAM_DICT_SIMUL['scale'])+"_vMax="+str(PARAM_DICT_SIMUL['vMax'])+"_b="+str(PARAM_DICT_SIMUL['b'])+"_r="+str(PARAM_DICT_SIMUL['remove'])+"_n="+str(PARAM_DICT_SIMUL['num_simulations']) 
    return filename