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