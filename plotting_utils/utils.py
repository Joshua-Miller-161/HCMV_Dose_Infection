import numpy as np
import pandas as pd
from itertools import chain
import sys
import os
from lmfit import Parameters, minimize, report_fit
import matplotlib.pyplot as plt

sys.path.append(os.getcwd())
from misc.misc_utils import ExtractParams
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
#====================================================================
def PlotText(ax, PARAM_DICT_SIMUL, xMin, xMax, yMin, yMax):
    remove_str = ''
    if (PARAM_DICT_SIMUL['remove'] == 1):
        remove_str = 'Successful virions\nremoved'
    elif (PARAM_DICT_SIMUL['remove'] == 0):
        remove_str = 'Successful virions\nleft in'

    if ('clump' in PARAM_DICT_SIMUL['simul_name']):
        if (PARAM_DICT_SIMUL['simul_name'] == 'clump'):
            ax.text(1.1 * xMin, .1 * yMax, "Scale: 1/" + str(PARAM_DICT_SIMUL['scale']) + ", " + r'$ \gamma = $' + str(PARAM_DICT_SIMUL['gamma']) + "\nvMax = " + str(PARAM_DICT_SIMUL['vMax'])+"\n"+remove_str)
        elif (PARAM_DICT_SIMUL['simul_name'] == 'clump_comp'):
            ax.text(1.1 * xMin, .1 * yMax, "Scale: 1/" + str(PARAM_DICT_SIMUL['scale']) + ", " + r'$ \gamma = $' + str(PARAM_DICT_SIMUL['gamma']) + "\nvMax = " + str(PARAM_DICT_SIMUL['vMax'])+' '+r'$\kappa=$'+str(PARAM_DICT_SIMUL['kappa'])+"\n"+remove_str)
        elif (PARAM_DICT_SIMUL['simul_name'] == 'clump_acc_dam'):
            ax.text(1.1 * xMin, .1 * yMax, "Scale: 1/" + str(PARAM_DICT_SIMUL['scale']) + ", " + r'$ \gamma = $' + str(PARAM_DICT_SIMUL['gamma']) + "\nvMax = " + str(PARAM_DICT_SIMUL['vMax']) + ' ' + r'$\beta=$'+str(PARAM_DICT_SIMUL['beta'])+"\n"+remove_str)
        
        elif (PARAM_DICT_SIMUL['simul_name'] == 'var_clump_diam'):
            if (PARAM_DICT_SIMUL['diameter_func'] == 'constant'):
                ax.text(1.1 * xMin, .01 * yMax, "Scale: 1/" + str(PARAM_DICT_SIMUL['scale']) + ", " + r'$ \gamma = $' + str(PARAM_DICT_SIMUL['gamma']) + "\nvMax = " + str(PARAM_DICT_SIMUL['vMax']) + '\n' + r'$\mu_{clump\_diam}=$' + str(PARAM_DICT_SIMUL['mean_clump_diam'])+"\n"+remove_str)
            elif (PARAM_DICT_SIMUL['diameter_func'] == 'linear'):
                ax.text(1.1 * xMin, .01 * yMax, "Scale: 1/" + str(PARAM_DICT_SIMUL['scale']) + ", " + r'$ \gamma = $' + str(PARAM_DICT_SIMUL['gamma']) + "\nvMax = " + str(PARAM_DICT_SIMUL['vMax']) + "\nvMaxD = " + str(PARAM_DICT_SIMUL['vMaxD'])+", bD = " + str(PARAM_DICT_SIMUL['bD']) + "\n"+remove_str)
            elif (PARAM_DICT_SIMUL['diameter_func'] == 'exponential'):
                ax.text(1.1 * xMin, .01 * yMax, "Scale: 1/" + str(PARAM_DICT_SIMUL['scale']) + ", " + r'$ \gamma = $' + str(PARAM_DICT_SIMUL['gamma']) + "\nvMax = " + str(PARAM_DICT_SIMUL['vMax']) + "\nAASDAJSD" + "\n"+remove_str)

        ax.text(1.1 * xMin, .005 * yMax, PARAM_DICT_SIMUL['scheme'])
        ax.text(1.1 * xMin, .0025 * yMax, PARAM_DICT_SIMUL['distribution'])

    elif (PARAM_DICT_SIMUL['simul_name'] == 'acc_dam'):
        ax.text(1.1 * xMin, .01 * yMax, "Scale: 1/" + str(PARAM_DICT_SIMUL['scale']) + ", " + r'$ \gamma = $' + str(PARAM_DICT_SIMUL['gamma']) + '\nvMax = ' + str(PARAM_DICT_SIMUL['vMax']) + ", " + r'$\beta$ = ' + str(PARAM_DICT_SIMUL['beta'])+"\n"+remove_str)
        ax.text(1.1 * xMin, .1 * yMax, r'$\lambda_{num\_interactions}=\frac{genomes/well}{vMax}$')

    elif (PARAM_DICT_SIMUL['simul_name'] == 'comp'):
        ax.text(1.1 * xMin, .01 * yMax, "Scale: 1/" + str(PARAM_DICT_SIMUL['scale']) + ", " + r'$ \gamma = $' + str(PARAM_DICT_SIMUL['gamma']) + '\nvMax = ' + str(PARAM_DICT_SIMUL['vMax']) + ", " + r'$\kappa$ = ' + str(PARAM_DICT_SIMUL['kappa'])+"\n"+remove_str)
        ax.text(1.1 * xMin, .1 * yMax, r'$\lambda_{num\_interactions}=\frac{genomes/well}{vMax}$')

    elif (PARAM_DICT_SIMUL['simul_name'] == 'null'):
        ax.text(1.1 * xMin, .01 * yMax, "Scale: 1/" + str(PARAM_DICT_SIMUL['scale']) + ", " + r'$ \gamma = $' + str(PARAM_DICT_SIMUL['gamma']) + '\nvMax = ' + str(PARAM_DICT_SIMUL['vMax']) + ", b = "+str(PARAM_DICT_SIMUL['b'])+"\n"+remove_str)
        ax.text(1.1 * xMin, .1 * yMax, r'$\lambda_{num\_interactions}=\frac{genomes/well}{vMax}+b$')
#====================================================================
def BasicFormat(ax, xMin=10**-3, xMax=10**3, yMin=10**-6, yMax=1, xlabel='Genomes/cell', ylabel='Infections/cell'):
    ax.plot(np. linspace(xMin, xMax), np.linspace(yMin, yMax), 'k--', linewidth = 1)
    ax.set_xlabel(xlabel)
    ax.set_ylabel('Infections/cell')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(xMin, xMax)
    ax.set_ylim(yMin, yMax)
#====================================================================
def PlotSimul(ax, file_path, band_type, replacement_val, x_name='GFP genomes (scaled)', y_name='GFP IU', color='green', marker='s'):
    df_simul = pd.read_csv(file_path)
    PARAM_DICT_SIMUL = ExtractParams(df_simul)

    num_simulations = int(PARAM_DICT_SIMUL['num_simulations'])
    cell_count      = PARAM_DICT_SIMUL['cell_count']

    GEN_WELL_SIMUL = np.asarray(df_simul.loc[:, x_name])
    GEN_CELL_SIMUL = GEN_WELL_SIMUL / cell_count
    INF_WELL_SIMUL = np.empty((num_simulations, df_simul.shape[0]), float)

    for i in range(num_simulations):
        INF_WELL_SIMUL[i, :] = df_simul.loc[:, y_name+' run='+str(i)]

    INF_CELL_SIMUL_MEAN = np.mean(INF_WELL_SIMUL, axis=0) / cell_count
    print("||||||", np.shape(INF_WELL_SIMUL), np.shape(INF_CELL_SIMUL_MEAN), "|||||||")

    #--------------------------------------------------------------------
    GEN_WELL_SIMUL_U, INF_WELL_SIMUL_U, CIS_U, MINS_U, MAXS_U = MakeDataPretty(df_simul, 'GFP genomes (scaled)', 'GFP IU', num_simulations)

    GEN_CELL_SIMUL_U      = GEN_WELL_SIMUL_U / cell_count
    INF_WELL_SIMUL_MEAN_U = np.mean(INF_WELL_SIMUL_U, axis=1)
    INF_CELL_SIMUL_MEAN_U = INF_WELL_SIMUL_MEAN_U.ravel() / cell_count
    #--------------------------------------------------------------------
    for i in range(np.shape(INF_WELL_SIMUL_MEAN_U)[0]):
        if (band_type == 'minmax'):
            if (MINS_U[i] <= 0):
                MINS_U[i] = replacement_val * cell_count

        elif (band_type == 'CIs'):
            if (INF_WELL_SIMUL_MEAN_U[i] <= 0):
                INF_WELL_SIMUL_MEAN_U[i] = replacement_val * cell_count
    #--------------------------------------------------------------------
    if (num_simulations == 1):
        ax.scatter(GEN_CELL_SIMUL, INF_WELL_SIMUL / cell_count, s=80, facecolors='none', edgecolors=color, marker=marker)

    elif (num_simulations > 1):
        if (band_type == 'minmax'):
            ax.fill_between(GEN_CELL_SIMUL_U, MINS_U / cell_count, MAXS_U / cell_count, color='black', alpha=.3)
        elif (band_type == 'CIs'):
            ax.fill_between(GEN_CELL_SIMUL_U, (INF_WELL_SIMUL_MEAN_U - CIS_U) / cell_count, (INF_WELL_SIMUL_MEAN_U + CIS_U) / cell_count, color='black', alpha=.3)

        ax.plot(GEN_CELL_SIMUL_U, INF_CELL_SIMUL_MEAN_U, color='black', linestyle='-')
    #--------------------------------------------------------------------
    return GEN_CELL_SIMUL, INF_CELL_SIMUL_MEAN
#====================================================================
def PlotFit(ax, x_data, y_data, lower_idx=0, upper_idx=-1, yMin=10**-6, yMax=1, color='red', linestyle='-'):
    params = Parameters()
    params.add('gamma', value=.45, min=0, max=1, vary=True)
    params.add('n', value=1, min=0, max=3, vary=True)
    print("=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+")
    result  = minimize(negLogLikeModel, params, method = 'differential_evolution', args=(x_data[lower_idx:upper_idx+1], y_data[lower_idx:upper_idx+1]),)
    report_fit(result)
    print("=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+")
    g = result.params['gamma'].value
    n = result.params['n'].value

    y_model = model(x_data, result.params)

    ax.plot(x_data[lower_idx:upper_idx+1], y_model[lower_idx:upper_idx+1], color=color, linestyle=linestyle, linewidth = 1)
    ax.axvline(x=x_data[lower_idx], ymin=yMin, ymax=yMax, color='black', alpha=0.3, zorder=0)
    ax.axvline(x=x_data[upper_idx], ymin=yMin, ymax=yMax, color='black', alpha=0.3, zorder=1)

    return g, n
#====================================================================
# Shoutout Pablo on StackExchange
def add_subplot_axes(ax, rect, axisbg='w'):
    fig = plt.gcf()
    box = ax.get_position()
    width = box.width
    height = box.height
    inax_position  = ax.transAxes.transform(rect[0:2])
    transFigure = fig.transFigure.inverted()
    infig_position = transFigure.transform(inax_position)    
    x = infig_position[0]
    y = infig_position[1]
    width *= rect[2]
    height *= rect[3]  # <= Typo was here
    subax = fig.add_axes([x,y,width,height],facecolor='none')  # matplotlib 2.0+
    #subax = fig.add_axes([x,y,width,height])
    x_labelsize = subax.get_xticklabels()[0].get_size()
    y_labelsize = subax.get_yticklabels()[0].get_size()
    x_labelsize *= rect[2]**0.5
    y_labelsize *= rect[3]**0.5
    subax.xaxis.set_tick_params(labelsize=x_labelsize)
    subax.yaxis.set_tick_params(labelsize=y_labelsize)
    return subax