import sys
sys.dont_write_bytecode = True
import numpy as np
import pandas as pd
from itertools import chain
from scipy.stats import chisquare
import os
from lmfit import Parameters, minimize, report_fit
import matplotlib.pyplot as plt
import random

sys.path.append(os.getcwd())
from misc.misc_utils import ExtractParams
#====================================================================
def model(x, params):
    #print('MODEL, MODEL, MODEL, MODEL, MODEL')
    try: # Get parameters
        g = params['gamma'].value
        n = params['n'].value
    except KeyError:
        g, n = params

    y = 1 - np.exp(-g * x**n)

    return y
#====================================================================
def model_n1(x, params):
    #print('MODEL_N1, MODEL_N1, MODEL_N1, MODEL_N1, MODEL_N1')
    try: # Get parameters
        g = params['gamma'].value
    except KeyError:
        g = params

    y = 1 - np.exp(-g * x)

    return y
#====================================================================
def negLogLikeModel(params, xData, yData, fix_n):
    new_xData = []
    new_yData = []

    for i in range(len(yData)):    
        if ((yData[i] != 0) and (xData[i] != 0)):
            new_xData.append(xData[i])
            new_yData.append(yData[i])

    new_xData = np.asarray(new_xData)
    new_yData = np.asarray(new_yData)

    if fix_n:
        #print("ASFLJSHF:ASJFHAS:FJASF: fix_n=", fix_n, "SD:LKASDAS")
        model_result = model_n1(new_xData, params)
    else:
        #print("ASFLJSHF:ASJFHAS:FJASF: fix_n=", fix_n, "SD:LKASDAS")
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
            ax.text(1.1 * xMin, .01 * yMax, "Scale: 1/" + str(PARAM_DICT_SIMUL['scale']) + ", " + r'$ \gamma = $' + str(PARAM_DICT_SIMUL['gamma']) + "\nvMax = " + str(PARAM_DICT_SIMUL['vMax'])+"\n"+remove_str)
        elif (PARAM_DICT_SIMUL['simul_name'] == 'clump_comp'):
            ax.text(1.1 * xMin, .01 * yMax, "Scale: 1/" + str(PARAM_DICT_SIMUL['scale']) + ", " + r'$ \gamma = $' + str(PARAM_DICT_SIMUL['gamma']) + "\nvMax = " + str(PARAM_DICT_SIMUL['vMax'])+' '+r'$\kappa=$'+str(PARAM_DICT_SIMUL['kappa'])+"\n"+remove_str)
        elif (PARAM_DICT_SIMUL['simul_name'] == 'clump_acc_dam'):
            ax.text(1.1 * xMin, .01 * yMax, "Scale: 1/" + str(PARAM_DICT_SIMUL['scale']) + ", " + r'$ \gamma = $' + str(PARAM_DICT_SIMUL['gamma']) + "\nvMax = " + str(PARAM_DICT_SIMUL['vMax']) + ' ' + r'$\beta=$'+str(PARAM_DICT_SIMUL['beta'])+"\n"+remove_str)
        
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
def PlotSimul(ax, file_path, band_type, replacement_val, x_name='GFP genomes (scaled)', y_name='GFP IU', 
              color='green', marker='^', s=20, scatter=True, return_gen_well=False):
    df_simul = pd.read_csv(file_path)
    PARAM_DICT_SIMUL = ExtractParams(df_simul)

    num_simulations = 0
    for key in df_simul.keys():
        if ("IU run" in key):
            num_simulations += 1
    print("ASDOASJD:ASKDJASDKASDAS:LDKNAS:DKLASD:KLASD:LASD:LASD:ASDLK", num_simulations)
    cell_count     = PARAM_DICT_SIMUL['cell_count']

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
        ax.scatter(GEN_CELL_SIMUL, INF_WELL_SIMUL / cell_count, s=s, facecolors='none', edgecolors=color, marker=marker)

    elif (num_simulations > 1):
        if (band_type == 'minmax'):
            ax.fill_between(GEN_CELL_SIMUL_U, MINS_U / cell_count, MAXS_U / cell_count, color='black', alpha=.3)
        elif (band_type == 'CIs'):
            ax.fill_between(GEN_CELL_SIMUL_U, (INF_WELL_SIMUL_MEAN_U - CIS_U) / cell_count, (INF_WELL_SIMUL_MEAN_U + CIS_U) / cell_count, color='black', alpha=.3)

        ax.plot(GEN_CELL_SIMUL_U, INF_CELL_SIMUL_MEAN_U, color='black', linestyle='-')
    
    if (scatter == True):
        ax.scatter(GEN_CELL_SIMUL_U, INF_CELL_SIMUL_MEAN_U, s=s, color=color, marker=marker)
    #--------------------------------------------------------------------
    if not return_gen_well:
        return GEN_CELL_SIMUL, INF_CELL_SIMUL_MEAN
    else:
        return GEN_WELL_SIMUL, GEN_CELL_SIMUL, INF_WELL_SIMUL, INF_CELL_SIMUL_MEAN
#====================================================================
def PlotFit(ax, x_data, y_data, lower_idx=0, upper_idx=-1, 
            plot_vlines=False, yMin_1=0, yMax_1=1, yMin_2=0, yMax_2=1, 
            plot_fit_markers=False, s=5, marker_color='green', marker='^',
            g_min=0, g_max=1, n_min=0, n_max=3,
            color='red', linestyle='-', n_fits=1):
    
    G_VALS   = np.ones(n_fits, float)
    N_NALS   = np.ones(n_fits, float)
    NLL_VALS_GOOD = np.ones(n_fits, float)
    NLL_VALS_BAD  = np.ones(n_fits, float)
    for i in range(n_fits):
        params = Parameters()
        params.add('gamma', value=random.uniform(g_min, g_max), min=g_min, max=g_max)
        params.add('n', value=random.uniform(n_min, n_max), min=n_min, max=n_max)
        #print("=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+")
        result  = minimize(negLogLikeModel, params, method='differential_evolution', args=(x_data[lower_idx:upper_idx+1], y_data[lower_idx:upper_idx+1], False),)
        #report_fit(result)
        #print("=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+")
        G_VALS[i]        = result.params['gamma'].value
        N_NALS[i]        = result.params['n'].value
        NLL_VALS_GOOD[i] = negLogLikeModel(result.params, x_data[lower_idx:upper_idx+1], y_data[lower_idx:upper_idx+1], False)
        #----------------------------------------------------------------
        # Get bad fit by fixing n=1
        params_1 = Parameters()
        params_1.add('gamma', value=random.uniform(g_min, g_max), min=g_min, max=g_max)

        result_1  = minimize(negLogLikeModel, params_1, method='differential_evolution', args=(x_data[lower_idx:upper_idx+1], y_data[lower_idx:upper_idx+1], True),)

        NLL_VALS_BAD[i] = negLogLikeModel(result_1.params, x_data[lower_idx:upper_idx+1], y_data[lower_idx:upper_idx+1], True)

        #print('n_vary:', NLL_VALS_GOOD[i], ", n_fixed:", NLL_VALS_BAD[i])
    #----------------------------------------------------------------
    # Get p-value
    res = chisquare([NLL_VALS_GOOD, NLL_VALS_BAD], axis=None)
    p = res.pvalue
    #----------------------------------------------------------------

    #y_model = model(x_data, result.params)
    best_idx = np.argmin(NLL_VALS_GOOD)
    g_best, n_best = G_VALS[best_idx], N_NALS[best_idx]
    print("g_best =", round(g_best, 5), ", n_best =", round(n_best, 5), "best_idx =", best_idx, ", p =", p)
    y_model = 1 - np.exp(-g_best * x_data**n_best)

    ax.plot(x_data[lower_idx:upper_idx+1], y_model[lower_idx:upper_idx+1], color=color, linestyle=linestyle, linewidth = 1)
    
    if plot_vlines:
        ax.axvline(x=x_data[lower_idx], ymin=yMin_1, ymax=yMax_1, color=color, alpha=0.3, zorder=0)
        ax.axvline(x=x_data[upper_idx], ymin=yMin_2, ymax=yMax_2, color=color, alpha=0.3, zorder=1)
    
    if plot_fit_markers:
        ax.scatter(x_data[lower_idx:upper_idx+1], y_data[lower_idx:upper_idx+1], 
                   s=s, facecolors=marker_color, marker=marker)

    return g_best, n_best, p
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
#====================================================================
def PlotVirionsInClumpFreq(ax, CLUMP_DICT, GEN_WELL_SIMUL, legend_size=4, 
                           scatter=True, markersize=50, cell_count=None, round_digits=5):
    CLUMP_MARKERS = ['s', '^', 'o']

    NUM_CLUMPS_2  = CLUMP_DICT[str(GEN_WELL_SIMUL[len(GEN_WELL_SIMUL) - 1])]
    CLUMP_SIZES_2 = np.arange(1, len(NUM_CLUMPS_2)+1)
    if not scatter:
        if not (cell_count == None):
            ax.vlines(CLUMP_SIZES_2, 0, NUM_CLUMPS_2, colors='goldenrod', lw=3, alpha=0.7, zorder=0,
                      label = 'Gen./cell = ' + str(round(GEN_WELL_SIMUL[len(GEN_WELL_SIMUL) - 1] / cell_count, round_digits)))
        else:
            ax.vlines(CLUMP_SIZES_2, 0, NUM_CLUMPS_2, colors='goldenrod', lw=3, alpha=0.7, zorder=0,
                      label = 'Gen./well = ' + str(GEN_WELL_SIMUL[len(GEN_WELL_SIMUL) - 1]))
    else:
        ax.vlines(CLUMP_SIZES_2, 0, NUM_CLUMPS_2, colors='goldenrod', lw=3, alpha=0.7, zorder=0)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    NUM_CLUMPS_1  = CLUMP_DICT[str(GEN_WELL_SIMUL[int(len(GEN_WELL_SIMUL) / 2)])]
    CLUMP_SIZES_1 = np.arange(1, np.shape(NUM_CLUMPS_1)[0]+1)

    if not scatter:
        if not (cell_count == None):
            ax.vlines(CLUMP_SIZES_1, 0, NUM_CLUMPS_1, colors='r', lw=3, alpha=0.7, zorder=1,
                      label = 'Gen./cell = ' + str(round(GEN_WELL_SIMUL[int(len(GEN_WELL_SIMUL) / 2)] / cell_count, round_digits)))
        else:
            ax.vlines(CLUMP_SIZES_1, 0, NUM_CLUMPS_1, colors='r', lw=3, alpha=0.7, zorder=1,
                      label = 'Gen./well = ' + str(GEN_WELL_SIMUL[int(len(GEN_WELL_SIMUL) / 2)]))
    else:
        ax.vlines(CLUMP_SIZES_1, 0, NUM_CLUMPS_1, colors='r', lw=3, alpha=0.5, zorder=1)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    NUM_CLUMPS_0 = CLUMP_DICT[str(GEN_WELL_SIMUL[0])]
    CLUMP_SIZES_0 = np.arange(1, np.shape(NUM_CLUMPS_0)[0]+1)
    if not scatter:
        if not (cell_count == None):
            ax.vlines(CLUMP_SIZES_0, 0, NUM_CLUMPS_0, colors='b', lw=3, alpha=0.7, zorder=2,
                      label = 'Gen./cell = ' + str(round(GEN_WELL_SIMUL[0] / cell_count, round_digits)))
        else:
            ax.vlines(CLUMP_SIZES_0, 0, NUM_CLUMPS_0, colors='b', lw=3, alpha=0.7, zorder=2,
                      label = 'Gen./well = ' + str(GEN_WELL_SIMUL[0]))
    else:
        ax.vlines(CLUMP_SIZES_0, 0, NUM_CLUMPS_0, colors='b', lw=3, alpha=0.5, zorder=2)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if scatter:
        if (cell_count == None):
            ax.scatter(CLUMP_SIZES_2, NUM_CLUMPS_2, s = markersize, color='goldenrod', edgecolors='black', marker=CLUMP_MARKERS[2], label = 'Gen./well = ' + str(GEN_WELL_SIMUL[len(GEN_WELL_SIMUL) - 1]))
            ax.scatter(CLUMP_SIZES_1, NUM_CLUMPS_1, s = markersize, color='red', edgecolors ='black', marker=CLUMP_MARKERS[1], label='Gen./well = ' + str(GEN_WELL_SIMUL[int(len(GEN_WELL_SIMUL) / 2)]))
            ax.scatter(CLUMP_SIZES_0, NUM_CLUMPS_0, s = markersize, color='blue', edgecolors='black', marker=CLUMP_MARKERS[0], label='Gen./well = ' + str(GEN_WELL_SIMUL[0]))
        else:
            ax.scatter(CLUMP_SIZES_2, NUM_CLUMPS_2, s = markersize, color='goldenrod', edgecolors='black', marker=CLUMP_MARKERS[2], label = 'Gen./cell = ' + str(round(GEN_WELL_SIMUL[len(GEN_WELL_SIMUL) - 1] / cell_count, round_digits)))
            ax.scatter(CLUMP_SIZES_1, NUM_CLUMPS_1, s = markersize, color='red', edgecolors ='black', marker=CLUMP_MARKERS[1], label='Gen./cell = ' + str(round(GEN_WELL_SIMUL[int(len(GEN_WELL_SIMUL) / 2)] / cell_count, round_digits)))
            ax.scatter(CLUMP_SIZES_0, NUM_CLUMPS_0, s = markersize, color='blue', edgecolors='black', marker=CLUMP_MARKERS[0], label='Gen./cell = ' + str(round(GEN_WELL_SIMUL[0] / cell_count, round_digits))) 
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ax.legend(loc='upper right', prop={'size':legend_size})