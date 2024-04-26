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
from plotting_utils.utils import negLogLikeModel, model, CombineSameGenWell, MakeDataPretty, MakeFilename
#====================================================================
''' Get simulation data '''

file = 'simulation_results/comp/CompSimul_2022_11_02_TB_GFP_fib_s=50_vMax=800.0_k=-3.0_r=1_n=3.csv'
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
               '2021_08_13 ME_mC_epithelial', '2022_11_02_TB_GFP_fib']
df_data = pd.read_excel('data/Experimental_data_Ed_Josh.xlsx', sheet_name=SHEET_NAMES[sheet])

print(SHEET_NAMES[6])
#--------------------------------------------------------------------
PARAM_DICT_DATA = ExtractParams(df_data)
cell_count = int(PARAM_DICT_DATA['cell_count'] / scale)
#====================================================================
''' Create figure '''

fig, ax = plt.subplots(1, 1, figsize=(6, 6), dpi=100, constrained_layout=False)
#====================================================================
with open('config.yml', 'r') as c:
    config = yaml.load(c, Loader=yaml.FullLoader)

LOWERS = config['PLOTTING']['LOWERS'] # [2, 20, 0, 6, 0, 0]
UPPERS = config['PLOTTING']['UPPERS'] # [36, -15, -1, 36, -1, -1]

print(LOWERS)
print(UPPERS)

MARKERS = config['PLOTTING']['markers_dict']#['o', '^', 's', 'D']
COLORS = config['PLOTTING']['colors_dict']

if ('GFP' in SHEET_NAMES[sheet]):
    color  = COLORS['GFP']
    marker = MARKERS['GFP']
elif (('cherry' in SHEET_NAMES[sheet]) or ('mCherry' in SHEET_NAMES[sheet]) or ('mC' in SHEET_NAMES[sheet])):
    color = COLORS['cherry']
    marker = MARKERS['cherry']

LETTERS = ['A', 'B', 'C']

replacement_val = config['PLOTTING']['replacement_val']
band_type = config['PLOTTING']['band_type']
#====================================================================
''' Prepare experimental and simulation data for plotting '''

GEN_WELL_SIMUL = np.asarray(df_simul.loc[:, 'GFP genomes (scaled)'])
GEN_CELL_SIMUL = GEN_WELL_SIMUL / cell_count
INF_WELL_SIMUL = np.empty((num_simulations, df_simul.shape[0]), float)

for i in range(num_simulations):
    INF_WELL_SIMUL[i, :] = df_simul.loc[:, 'GFP IU run='+str(i)]

print("||||||", np.shape(INF_WELL_SIMUL), np.shape(np.mean(INF_WELL_SIMUL, axis=0)), "|||||||")

INF_WELL_SIMUL_MEAN = np.mean(INF_WELL_SIMUL, axis=0)
INF_CELL_SIMUL_MEAN = INF_WELL_SIMUL.ravel() / cell_count
#--------------------------------------------------------------------
GEN_WELL_SIMUL_U, INF_WELL_SIMUL_U, CIS_U, MINS_U, MAXS_U = MakeDataPretty(df_simul, 'GFP genomes (scaled)', 'GFP IU', num_simulations)

GEN_CELL_SIMUL_U      = GEN_WELL_SIMUL_U / cell_count
INF_WELL_SIMUL_MEAN_U = np.mean(INF_WELL_SIMUL_U, axis=1)
INF_CELL_SIMUL_MEAN_U = INF_WELL_SIMUL_MEAN_U.ravel() / cell_count
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
result_simul = minimize(negLogLikeModel, params, method = 'differential_evolution', args=(GEN_CELL_SIMUL[LOWERS[sheet]:UPPERS[sheet]], INF_CELL_SIMUL_MEAN[LOWERS[sheet]:UPPERS[sheet]]),)
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
''' Replace values of 0 with a replacement value '''
print("AHHHHHHHH ", np.shape(INF_WELL_SIMUL_U), np.shape(INF_WELL_SIMUL_MEAN_U))
for i in range(np.shape(INF_WELL_SIMUL_MEAN_U)[0]):
    if (band_type == 'minmax'):
        if (MINS_U[i] <= 0):
            MINS_U[i] = replacement_val * cell_count

    elif (band_type == 'CIs'):
        if (INF_WELL_SIMUL_MEAN_U[i] <= 0):
            INF_WELL_SIMUL_MEAN_U[i] = replacement_val * cell_count
#====================================================================
''' Plot data and simulation results and best fit curves '''

ax.scatter(GEN_CELL_DATA, INF_CELL_DATA, s=80, facecolors=color, edgecolors='none', marker=marker, alpha=.3)
ax.plot(x_data[LOWERS[sheet]:UPPERS[sheet]], y_data[LOWERS[sheet]:UPPERS[sheet]], 'b-.', linewidth = 2)
#--------------------------------------------------------------------
if (num_simulations == 1):
    ax.scatter(GEN_CELL_SIMUL, INF_WELL_SIMUL_MEAN.ravel() / cell_count, s=80, facecolors='none', edgecolors=color, marker=marker)

elif (num_simulations > 1):
    if (band_type == 'minmax'):
        ax.fill_between(GEN_WELL_SIMUL_U / cell_count, MINS_U / cell_count, MAXS_U / cell_count, color='black', alpha=.3)
    elif (band_type == 'CIs'):
        ax.fill_between(GEN_WELL_SIMUL_U / cell_count, (INF_WELL_SIMUL_MEAN_U - CIS_U) / cell_count, (INF_WELL_SIMUL_MEAN_U + CIS_U) / cell_count, color='black', alpha=.3)

    ax.plot(GEN_WELL_SIMUL_U / cell_count, INF_CELL_SIMUL_MEAN_U, color='black', linestyle='-')

ax.plot(x_simul[LOWERS[sheet]:UPPERS[sheet]], y_simul[LOWERS[sheet]:UPPERS[sheet]], 'r--', linewidth = 2)

#ax.scatter(GEN_CELL_DATA, INF_CELL_DATA, s=80, facecolors='green', edgecolors='none', marker='s')
#ax.plot(x_data[LOWERS[sheet]:UPPERS[sheet]], y_data[LOWERS[sheet]:UPPERS[sheet]], color='pink', linestyle='-.', linewidth = 2)
#====================================================================
''' Formatting '''
legendD = mlines.Line2D([], [], color='b', linestyle='-.', markerfacecolor=color, markeredgecolor='none', markerfacecoloralt='none', marker=marker, alpha=0.6,
                        markersize=10, label = "HCMV-TB (GFP) (data): n = " + str(round(n_data, 3)))
legendS = ''
if (num_simulations == 1):
    legendS = mlines.Line2D([], [], color='y', linestyle='--', markerfacecolor='none', markeredgecolor=color, markerfacecoloralt='none', marker=marker,
                            markersize=10, label= "Simulation: n = " + str(round(n_simul, 3)))
elif (num_simulations > 1):
    legendS = mlines.Line2D([], [], color='k', linestyle='-', label=r'Simul. mean, $\overline{n} = $' + str(round(n_simul, 3)) + " ("+str(num_simulations)+" runs)")
#--------------------------------------------------------------------
xMin = 10**-3
xMax = 10**3
yMin = 10**-6
yMax = 1

ax.legend(handles = [legendD, legendS], loc='upper left', prop={'size': 8})
ax.plot(np. linspace(xMin,xMax), np.linspace(yMin, yMax), 'k--', linewidth = 1)
ax.set_xlabel('Genomes/cell')
ax.set_ylabel('Infections/cell')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(xMin, xMax)
ax.set_ylim(yMin, yMax)
#--------------------------------------------------------------------
remove_str = ''
if (PARAM_DICT_SIMUL['remove'] == 1):
    remove_str = 'Successful virions\nremoved'
elif (PARAM_DICT_SIMUL['remove'] == 0):
    remove_str = 'Successful virions\nleft in'

if ('clump' in simul_name):
    if (simul_name == 'clump'):
        ax.set_title(SHEET_NAMES[sheet] + ' | Clumping')
        ax.text(1.1 * xMin, .1 * yMax, "Scale: 1/" + str(scale) + ", " + r'$ \gamma = $' + str(gamma) + "\nvMax = " + str(vMax)+"\n"+remove_str)
    elif (simul_name == 'clump_comp'):
        ax.set_title(SHEET_NAMES[sheet] + ' | Clumping + Compensation')
        ax.text(1.1 * xMin, .1 * yMax, "Scale: 1/" + str(scale) + ", " + r'$ \gamma = $' + str(gamma) + "\nvMax = " + str(vMax)+' '+r'$\kappa=$'+str(PARAM_DICT_SIMUL[kappa])+"\n"+remove_str)
    elif (simul_name == 'clump_acc_dam'):
        ax.set_title(SHEET_NAMES[sheet] + ' | Clumping + Acc. Damage')
        ax.text(1.1 * xMin, .1 * yMax, "Scale: 1/" + str(scale) + ", " + r'$ \gamma = $' + str(gamma) + "\nvMax = " + str(vMax) + ' ' + r'$\beta=$'+str(PARAM_DICT_SIMUL['beta'])+"\n"+remove_str)
    ax.text(1.1 * xMin, .05 * yMax, scheme)
    ax.text(1.1 * xMin, .025 * yMax, distribution)

elif (simul_name == 'acc_dam'):
    ax.text(1.1 * xMin, .01 * yMax, "Scale: 1/" + str(scale) + ", " + r'$ \gamma = $' + str(gamma) + '\nvMax = ' + str(vMax) + ", " + r'$\beta$ = ' + str(beta)+"\n"+remove_str)
    ax.text(1.1 * xMin, .1 * yMax, r'$\lambda_{num\_interactions}=\frac{genomes/well}{vMax}$')
    ax.set_title(SHEET_NAMES[sheet] + ' | Acc. Damage')

elif (simul_name == 'comp'):
    ax.text(1.1 * xMin, .01 * yMax, "Scale: 1/" + str(scale) + ", " + r'$ \gamma = $' + str(gamma) + '\nvMax = ' + str(vMax) + ", " + r'$\kappa$ = ' + str(kappa)+"\n"+remove_str)
    ax.text(1.1 * xMin, .1 * yMax, r'$\lambda_{num\_interactions}=\frac{genomes/well}{vMax}$')
    ax.set_title(SHEET_NAMES[sheet] + ' | Compensation')

elif (simul_name == 'null'):
    ax.text(1.1 * xMin, .01 * yMax, "Scale: 1/" + str(scale) + ", " + r'$ \gamma = $' + str(gamma) + '\nvMax = ' + str(vMax) + ", b = "+str(PARAM_DICT_SIMUL['b'])+"\n"+remove_str)
    ax.text(1.1 * xMin, .1 * yMax, r'$\lambda_{num\_interactions}=\frac{genomes/well}{vMax}+b$')
    ax.set_title(SHEET_NAMES[sheet] + ' | Null')

ax.text(.6 * xMin, 2.2 * yMax, '', fontsize=16, fontweight='bold', va='top', ha='right')
#====================================================================
''' Save figure '''
filename = MakeFilename(PARAM_DICT_SIMUL, SHEET_NAMES[sheet])
fig.savefig(os.path.join(os.path.join(os.getcwd(), 'figs'), filename+".pdf"), bbox_inches = 'tight', pad_inches = 0) # Save figure in the new directory
#====================================================================
''' Create figure '''
fig2, ax2 = plt.subplots(1,1,figsize=(6, 6), dpi=100, constrained_layout=False)
#--------------------------------------------------------------------
GEN_WELL_SIMUL_U, INTER_U, CIS_U, MINS_U, MAXS_U = MakeDataPretty(df_simul, 'GFP genomes (scaled)', 'total_interactions', num_simulations)

INTER_CELL = np.mean(INTER_U, axis=1) / cell_count
INTER_WELL = np.mean(INTER_U, axis=1)  / GEN_WELL_SIMUL_U

ax2.plot(GEN_WELL_SIMUL_U / cell_count, INTER_CELL,label='interactions / cell_count')
ax2.plot(GEN_WELL_SIMUL_U / cell_count, INTER_WELL, label='interactions / GENOMES/WELL')
#--------------------------------------------------------------------
xMin = 10**-3
xMax = 10**3
yMin = 10**(int(np.log10(min([min(INTER_CELL), min(INTER_WELL)])))-1)
yMax = 10**(int(np.log10(max([max(INTER_CELL), max(INTER_WELL)]))) + 1)

ax2.set_xlabel('Genomes/cell')
ax2.set_ylabel('Avg. interactions')
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlim(xMin, xMax)
ax2.set_ylim(yMin, yMax)
ax2.legend(loc='upper left')

remove_str = ''
if (PARAM_DICT_SIMUL['remove'] == 1):
    remove_str = 'Successful virions\nremoved'
elif (PARAM_DICT_SIMUL['remove'] == 0):
    remove_str = 'Successful virions\nleft in'

if ('clump' in simul_name):
    if (simul_name == 'clump'):
        ax2.set_title(SHEET_NAMES[sheet] + ' | Clumping')
        ax2.text(1.1 * xMin, .1 * yMax, "Scale: 1/" + str(scale) + ", " + r'$ \gamma = $' + str(gamma) + "\nvMax = " + str(vMax)+"\n"+remove_str)
    elif (simul_name == 'clump_comp'):
        ax2.set_title(SHEET_NAMES[sheet] + ' | Clumping + Compensation')
        ax2.text(1.1 * xMin, .1 * yMax, "Scale: 1/" + str(scale) + ", " + r'$ \gamma = $' + str(gamma) + "\nvMax = " + str(vMax)+' '+r'$\kappa=$'+str(PARAM_DICT_SIMUL[kappa])+"\n"+remove_str)
    elif (simul_name == 'clump_acc_dam'):
        ax2.set_title(SHEET_NAMES[sheet] + ' | Clumping + Acc. Damage')
        ax2.text(1.1 * xMin, .1 * yMax, "Scale: 1/" + str(scale) + ", " + r'$ \gamma = $' + str(gamma) + "\nvMax = " + str(vMax) + ' ' + r'$\beta=$'+str(PARAM_DICT_SIMUL['beta'])+"\n"+remove_str)
    ax2.text(1.1 * xMin, .05 * yMax, scheme)
    ax2.text(1.1 * xMin, .025 * yMax, distribution)

elif (simul_name == 'acc_dam'):
    ax2.text(1.1 * xMin, .01 * yMax, "Scale: 1/" + str(scale) + ", " + r'$ \gamma = $' + str(gamma) + '\nvMax = ' + str(vMax) + ", " + r'$\beta$ = ' + str(beta)+"\n"+remove_str)
    ax2.text(1.1 * xMin, .1 * yMax, r'$\lambda_{num\_interactions}=\frac{genomes/well}{vMax}$')
    ax2.set_title(SHEET_NAMES[sheet] + ' | Acc. Damage')

elif (simul_name == 'comp'):
    ax2.text(1.1 * xMin, .01 * yMax, "Scale: 1/" + str(scale) + ", " + r'$ \gamma = $' + str(gamma) + '\nvMax = ' + str(vMax) + ", " + r'$\kappa$ = ' + str(kappa)+"\n"+remove_str)
    ax2.text(1.1 * xMin, .1 * yMax, r'$\lambda_{num\_interactions}=\frac{genomes/well}{vMax}$')
    ax2.set_title(SHEET_NAMES[sheet] + ' | Compensation')

elif (simul_name == 'null'):
    ax2.text(1.1 * xMin, .01 * yMax, "Scale: 1/" + str(scale) + ", " + r'$ \gamma = $' + str(gamma) + '\nvMax = ' + str(vMax) + ", b = "+str(PARAM_DICT_SIMUL['b'])+"\n"+remove_str)
    ax2.text(1.1 * xMin, .1 * yMax, r'$\lambda_{num\_interactions}=\frac{genomes/well}{vMax}+b$')
    ax2.set_title(SHEET_NAMES[sheet] + ' | Null')
#====================================================================
''' Save figure '''
filename = MakeFilename(PARAM_DICT_SIMUL, SHEET_NAMES[sheet])
fig2.savefig(os.path.join(os.path.join(os.getcwd(), 'figs/interaction_figs'), "INTER_"+filename+".pdf"), bbox_inches = 'tight', pad_inches = 0) # Save figure in the new directory
plt.show()