import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import pandas as pd
from lmfit import Parameters, minimize, report_fit
import yaml
import sys
import os

sys.path.append(os.getcwd())

from misc.misc_utils import ExtractParams
from simulation_utils.utils import PrepareData
from plotting_utils.utils import negLogLikeModel, model, CombineSameGenWell
#====================================================================
''' Get simulation data '''

file = 'simulation_results/comp/CompSimul_2021_07_13 GFP_TB_fibroblast_s=100_vMax=800.0_k=-0.03_r=1_n=1.csv'
df_simul = pd.read_csv(file)
#--------------------------------------------------------------------
PARAM_DICT_SIMUL = ExtractParams(df_simul)

simul_name      = PARAM_DICT_SIMUL['simul_name']
num_simulations = int(PARAM_DICT_SIMUL['num_simulations'])
sheet           = int(PARAM_DICT_SIMUL['sheet'])
scale           = PARAM_DICT_SIMUL['scale']
i_mean          = PARAM_DICT_SIMUL['muG']
i_stdev         = PARAM_DICT_SIMUL['sigmaG']
r_mean          = PARAM_DICT_SIMUL['muR']
r_stdev         = PARAM_DICT_SIMUL['sigmaR']
gamma           = PARAM_DICT_SIMUL['gamma']
vMax            = PARAM_DICT_SIMUL['vMax']

assert simul_name in ['clump', 'comp', 'acc_dam', 'null', 'clump_comp', 'clump_acc_dam'], "Got: "+simul_name+". 'simul_name' must be one of 'clump', 'comp', 'acc_dam', 'null', 'clump_comp', or 'clump_acc_dam'."

if ('clump' in simul_name):
    scheme = PARAM_DICT_SIMUL['scheme']
    distribution = PARAM_DICT_SIMUL['distribution']
elif (simul_name == 'acc_dam'):
    beta = PARAM_DICT_SIMUL['beta']
elif (simul_name == 'comp'):
    kappa = PARAM_DICT_SIMUL['kappa']
elif (simul_name == 'null'):
    b = PARAM_DICT_SIMUL['b']
#====================================================================
''' Get corresponding experimental data '''

SHEET_NAMES = ['2021_10_05 TB_GFP_epithelial', '2020_07_02 ME_GFP_fibroblast', 
               '2020_05_29 TR_GFP_fibroblast', '2021_07_13 GFP_TB_fibroblast', 
               '2020_08_12 TB_GFP_fibroblast', '2020_09_14 TR_GFP_epithelial',
               '2022_11_02_TB_GFP_fib']
df_data = pd.read_excel('data/Experimental_data_Ed_Josh.xlsx', sheet_name=SHEET_NAMES[sheet])
#--------------------------------------------------------------------
PARAM_DICT_DATA = ExtractParams(df_data)
cell_count = int(PARAM_DICT_DATA['cell_count'] / scale)
#====================================================================
''' Create figure '''

fig = plt.figure(figsize=(6, 6), dpi=100, constrained_layout=False)
#====================================================================
with open('config.yml', 'r') as c:
    config = yaml.load(c, Loader=yaml.FullLoader)

LOWERS = config['PLOTTING']['LOWERS'] # [2, 20, 0, 6, 0, 0]
UPPERS = config['PLOTTING']['UPPERS'] # [36, -15, -1, 36, -1, -1]

print(LOWERS)
print(UPPERS)

MARKERS = config['PLOTTING']['markers_dict']#['o', '^', 's', 'D']
COLORS = config['PLOTTING']['colors_dict']
LETTERS = ['A', 'B', 'C']
#====================================================================
''' Prepare experimental and simulation data for plotting '''

GEN_WELL_SIMUL  = np.asarray(df_simul.loc[:, 'GFP genomes (scaled)'])
GEN_CELL_SIMUL  = GEN_WELL_SIMUL / cell_count
INF_WELL_SIMULS = np.empty((num_simulations, df_simul.shape[0]), float)

for i in range(num_simulations):
    INF_WELL_SIMULS[i, :] = df_simul.loc[:, 'GFP IU run='+str(i)]

print("||||||", np.shape(INF_WELL_SIMULS), np.shape(np.mean(INF_WELL_SIMULS, axis=0)), "|||||||")

INF_WELL_SIMULS_MEAN = np.mean(INF_WELL_SIMULS, axis=0)
INF_CELL_SIMULS_MEAN = INF_WELL_SIMULS.ravel() / cell_count

INF_WELL_SIMULS_STDEVS = np.std(INF_WELL_SIMULS, axis=0)
CIS = (1.96 * np.ones(np.shape(INF_WELL_SIMULS)[1])) * INF_WELL_SIMULS_STDEVS / (np.sqrt(np.shape(INF_WELL_SIMULS)[0]) * np.ones(np.shape(INF_WELL_SIMULS)[1]))

MAXS = np.amax(INF_WELL_SIMULS, axis=0)
MINS = np.amin(INF_WELL_SIMULS, axis=0)

#--------------------------------------------------------------------
GEN_WELL_SIMUL_U, INF_WELL_SIMULS_U = CombineSameGenWell(GEN_WELL_SIMUL, INF_WELL_SIMULS)

INF_WELL_SIMULS_MEAN_U = np.mean(INF_WELL_SIMULS_U, axis=1)
INF_CELL_SIMULS_MEAN_U = INF_WELL_SIMULS_MEAN_U.ravel() / cell_count

INF_WELL_SIMULS_STDEVS_U = np.std(INF_WELL_SIMULS_U, axis=1)
CIS_U = (1.96 * np.ones(np.shape(INF_WELL_SIMULS_U)[0])) * INF_WELL_SIMULS_STDEVS_U / (np.sqrt(np.shape(INF_WELL_SIMULS_U)[1]) * np.ones(np.shape(INF_WELL_SIMULS_U)[0]))

MAXS_U = np.amax(INF_WELL_SIMULS_U, axis=1)
MINS_U = np.amin(INF_WELL_SIMULS_U, axis=1)

print(np.shape(GEN_WELL_SIMUL_U), np.shape(INF_WELL_SIMULS_U), np.shape(MAXS_U))
#--------------------------------------------------------------------
GEN_WELL_DATA, GEN_CELL_DATA, INF_CELL_DATA, num_zeros = PrepareData(df_data, scale)
#====================================================================
''' Fit model to experimental and simulated data '''
params = Parameters()
params.add('gamma', value=.45, min=0, max=1, vary=True)
params.add('n', value=1, min=0, max=3, vary=True)
print("=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+")
result_data  = minimize(negLogLikeModel, params, method = 'differential_evolution', args=(GEN_CELL_DATA[LOWERS[sheet]:UPPERS[sheet]], INF_CELL_DATA[LOWERS[sheet]:UPPERS[sheet]]),)
report_fit(result_data)
print("=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+")
print("+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=")
result_simul = minimize(negLogLikeModel, params, method = 'differential_evolution', args=(GEN_CELL_SIMUL[LOWERS[sheet]:UPPERS[sheet]], INF_CELL_SIMULS_MEAN[LOWERS[sheet]:UPPERS[sheet]]),)
report_fit(result_simul)
print("=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+")
g_data = result_data.params['gamma'].value
n_data = result_data.params['n'].value
g_simul = result_simul.params['gamma'].value
n_simul = result_simul.params['n'].value
#g2 = result_simul_low.params['gamma'].value
#n2 = result_simul_low.params['n'].value

print("g_data = ", round(g_data,5), ", n_data = ", round(n_data,5), ", g_simul = ", round(g_simul,5), ", n_simul = ", round(n_simul,5), 'cell_count =', cell_count)#, ", g2 = ", g2, ", n2 = ", n2)

x_data = GEN_CELL_DATA
y_data = model(x_data, result_data.params)
x_simul = GEN_CELL_SIMUL
y_simul = model(x_simul, result_simul.params)
#x2 = GFP_GEN_CELL
#y2 = model(x2, result_simul_low.params)
print("=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+")
print("+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=")
#====================================================================
''' Plot simulation results and best fit curves '''

plt.scatter(GEN_CELL_DATA, INF_CELL_DATA, s=80, facecolors=COLORS['GFP'], edgecolors='none', marker=MARKERS['GFP'], alpha=.3)
plt.plot(x_data[LOWERS[sheet]:UPPERS[sheet]], y_data[LOWERS[sheet]:UPPERS[sheet]], 'b-.', linewidth = 2)

if (num_simulations == 1):
    plt.scatter(GEN_CELL_SIMUL, INF_CELL_SIMULS_MEAN.ravel(), s=80, facecolors='none', edgecolors=COLORS['GFP'], marker=MARKERS['GFP'])
elif (num_simulations > 1):
    #plt.fill_between(GEN_CELL_SIMUL, (INF_WELL_SIMULS_MEAN-CIS) / cell_count, (INF_WELL_SIMULS_MEAN+CIS) / cell_count, color='black', alpha=.3)
    
    #plt.fill_between(GEN_CELL_SIMUL, MINS / cell_count, MAXS / cell_count, color='black', alpha=.3)
    
    plt.fill_between(GEN_WELL_SIMUL_U / cell_count, MINS_U / cell_count, MAXS_U / cell_count, color='black', alpha=.3)
    plt.plot(GEN_WELL_SIMUL_U / cell_count, INF_CELL_SIMULS_MEAN_U, color='black', linestyle='-')
    #plt.scatter(GEN_WELL_SIMUL_U / cell_count, MINS_U / cell_count, color='black', alpha=.3)
    #plt.scatter(GEN_WELL_SIMUL_U / cell_count, MAXS_U / cell_count, color='black', alpha=.3)
    
    # for gen_well, inf_well_simul in INF_WELL_SIMULS_DICT.items():
    #     plt.scatter(np.asarray([float(gen_well)] * len(inf_well_simul)) / cell_count, 
    #                 inf_well_simul / cell_count,
    #                 s=20, facecolors='greenyellow', edgecolors='none', marker='o', alpha=.3)

plt.plot(x_simul[LOWERS[sheet]:UPPERS[sheet]], y_simul[LOWERS[sheet]:UPPERS[sheet]], 'r--', linewidth = 2)

#plt.scatter(GEN_CELL_DATA, INF_CELL_DATA, s=80, facecolors='green', edgecolors='none', marker='s')
#plt.plot(x_data[LOWERS[sheet]:UPPERS[sheet]], y_data[LOWERS[sheet]:UPPERS[sheet]], color='pink', linestyle='-.', linewidth = 2)
#====================================================================
''' Formatting '''
legendD = mlines.Line2D([], [], color='b', linestyle='-.', markerfacecolor='none', markeredgecolor=COLORS['GFP'], markerfacecoloralt='none', marker=MARKERS['GFP'],
                          markersize=10, label = "HCMV-TB (GFP) (data): n = " + str(round(n_data, 3)))
legendS = ''
if (num_simulations == 1):
    legendS = mlines.Line2D([], [], color='y', linestyle='--', markerfacecolor='none', markeredgecolor=COLORS['GFP'], markerfacecoloralt='none', marker=MARKERS['GFP'],
                            markersize=10, label= "Simulation: n = " + str(round(n_simul, 3)))
elif (num_simulations > 1):
    legendS = mlines.Line2D([], [], color='k', linestyle='-', label=r'Simul. mean, $\overline{n} = $' + str(round(n_simul, 3)) + " ("+str(num_simulations)+" runs)")
#--------------------------------------------------------------------
xMin = 10**-3
xMax = 10**3
yMin = 10**-6
yMax = 1

plt.legend(handles = [legendD, legendS], loc='upper left', prop={'size': 8})
plt.plot(np. linspace(xMin,xMax), np.linspace(yMin, yMax), 'k--', linewidth = 1)
plt.xlabel('Genomes/cell')
plt.ylabel('Infections/cell')
plt.xscale('log')
plt.yscale('log')
plt.xlim(xMin, xMax)
plt.ylim(yMin, yMax)
#--------------------------------------------------------------------
remove_str = ''
if (PARAM_DICT_SIMUL['remove'] == 1):
    remove_str = 'Successful virions\nremoved'
elif (PARAM_DICT_SIMUL['remove'] == 0):
    remove_str = 'Successful virions\nleft in'

if ('clump' in simul_name):
    if (simul_name == 'clump'):
        plt.title(SHEET_NAMES[sheet] + ' | Clumping')
        plt.text(1.1 * xMin, .1 * yMax, "Scale: 1/" + str(scale) + ", " + r'$ \gamma = $' + str(gamma) + "\n vMax = " + str(vMax)+"\n"+remove_str)
    elif (simul_name == 'clump_comp'):
        plt.title(SHEET_NAMES[sheet] + ' | Clumping + Compensation')
        plt.text(1.1 * xMin, .1 * yMax, "Scale: 1/" + str(scale) + ", " + r'$ \gamma = $' + str(gamma) + "\n vMax = " + str(vMax)+' '+r'$\kappa=$'+str(PARAM_DICT_SIMUL[kappa])+"\n"+remove_str)
    elif (simul_name == 'clump_acc_dam'):
        plt.title(SHEET_NAMES[sheet] + ' | Clumping + Acc. Damage')
        plt.text(1.1 * xMin, .1 * yMax, "Scale: 1/" + str(scale) + ", " + r'$ \gamma = $' + str(gamma) + "\n vMax = " + str(vMax) + ' ' + r'$\beta=$'+str(PARAM_DICT_SIMUL['beta'])+"\n"+remove_str)
    plt.text(1.1 * xMin, .05 * yMax, scheme)
    plt.text(1.1 * xMin, .025 * yMax, distribution)

elif (simul_name == 'acc_dam'):
    plt.text(1.1 * xMin, .01 * yMax, "Scale: 1/" + str(scale) + ", " + r'$ \gamma = $' + str(gamma) + '\n vMax = ' + str(vMax) + ", " + r'$\beta$ = ' + str(beta)+"\n"+remove_str)
    plt.text(1.1 * xMin, .1 * yMax, r'$\lambda_{num\_interactions}=\frac{genomes/well}{vMax}$')
    plt.title(SHEET_NAMES[sheet] + ' | Acc. Damage')

elif (simul_name == 'comp'):
    plt.text(1.1 * xMin, .01 * yMax, "Scale: 1/" + str(scale) + ", " + r'$ \gamma = $' + str(gamma) + '\n vMax = ' + str(vMax) + ", " + r'$\kappa$ = ' + str(kappa)+"\n"+remove_str)
    plt.text(1.1 * xMin, .1 * yMax, r'$\lambda_{num\_interactions}=\frac{genomes/well}{vMax}$')
    plt.title(SHEET_NAMES[sheet] + ' | Compensation')

elif (simul_name == 'null'):
    plt.text(1.1 * xMin, .01 * yMax, "Scale: 1/" + str(scale) + ", " + r'$ \gamma = $' + str(gamma) + '\n vMax = ' + str(vMax) + ", b = "+str(PARAM_DICT_SIMUL['b'])+"\n"+remove_str)
    plt.text(1.1 * xMin, .1 * yMax, r'$\lambda_{num\_interactions}=\frac{genomes/well}{vMax}+b$')
    plt.title(SHEET_NAMES[sheet] + ' | Null')

plt.text(.6 * xMin, 2.2 * yMax, '', fontsize=16, fontweight='bold', va='top', ha='right')
#====================================================================
plt.show()
#====================================================================
''' Save figure '''
filename = ''
if ('clump' in simul_name):
    scheme_short = ''
    if (PARAM_DICT_SIMUL['scheme']=='linear'):
        scheme_short='lin'
    elif (PARAM_DICT_SIMUL['scheme']=='regular_polygon'):
        scheme_short='poly'
    elif (PARAM_DICT_SIMUL['scheme']=='sphere_packing'):
        scheme_short='sp'

    dist_short = ''
    if (distribution=='normal'):
        dist_short = 'norm'
    elif (distribution=='uniform'):
        dist_short = 'uni'
    elif (distribution=='fixed'):
        dist_short = 'fix'
    
    if (simul_name == 'clump'):
        filename = "ClumpSimul_"+SHEET_NAMES[sheet]+"_s="+str(scale)+"_vMax="+str(vMax)+"_"+scheme_short+"_"+dist_short+"_r="+str(PARAM_DICT_SIMUL['remove'])+"_n="+str(num_simulations) # Specify filename
    elif (simul_name == 'clump_comp'):
        filename = "ClumpCompSimul_"+SHEET_NAMES[sheet]+"_s="+str(scale)+"_vMax="+str(vMax)+"_k="+str(PARAM_DICT_SIMUL['kappa'])+"_"+scheme_short+"_"+dist_short+"_r="+str(PARAM_DICT_SIMUL['remove'])+"_n="+str(num_simulations) # Specify filename
    elif (simul_name == 'clump_acc_dam'):
        filename = "ClumpAccDamSimul_"+SHEET_NAMES[sheet]+"_s="+str(scale)+"_vMax="+str(vMax)+"_b="+str(PARAM_DICT_SIMUL['beta'])+"_"+scheme_short+"_"+dist_short+"_r="+str(PARAM_DICT_SIMUL['remove'])+"_n="+str(num_simulations) # Specify filename

elif (simul_name == 'acc_dam'):
    filename = "AccDamSimul_"+SHEET_NAMES[sheet]+"_s="+str(scale)+"_vMax="+str(vMax)+"_b="+str(beta)+"_r="+str(PARAM_DICT_SIMUL['remove'])+"_n="+str(num_simulations)

elif (simul_name == 'comp'):
    filename = "CompSimul_"+SHEET_NAMES[sheet]+"_s="+str(scale)+"_vMax="+str(vMax)+"_k="+str(kappa)+"_r="+str(PARAM_DICT_SIMUL['remove'])+"_n="+str(num_simulations)

elif (simul_name == 'null'):
    filename = "NullSimul_"+SHEET_NAMES[sheet]+"_s="+str(scale)+"_vMax="+str(vMax)+"_b="+str(b)+"_r="+str(PARAM_DICT_SIMUL['remove'])+"_n="+str(num_simulations) 

fig.savefig(os.path.join(os.path.join(os.getcwd(), 'figs'), filename+".pdf"), bbox_inches = 'tight', pad_inches = 0) # Save figure in the new directory