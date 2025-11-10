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
from plotting_utils.utils import negLogLikeModel, model, CombineSameGenWell, MakeDataPretty, MakeFilename, PlotText, PlotFit, PlotSimul, BasicFormat
#====================================================================
''' Get simulation data '''

#file = "simulation_results/clump/ClumpSimulVVG_use_with_size_dist_interp_s=10_vMax=10000.0_sp_norm_r=1_n=10.csv"
#file = "simulation_results/clump/ClumpSimulVVG_use_with_size_dist_interp_s=10_vMax=10000.0_poly_norm_r=1_n=10.csv"
#file = "simulation_results/clump/ClumpSimulVVG_use_with_size_dist_interp_s=10_vMax=10000.0_lin_norm_r=1_n=10.csv"

#file = "simulation_results/clump/ClumpSimulVVG_use_with_size_dist_interp_s=10_vMax=10000.0_sp_fix_r=1_n=10.csv"
#file = "simulation_results/clump/ClumpSimulVVG_use_with_size_dist_interp_s=10_vMax=10000.0_poly_fix_r=1_n=10.csv"
#file = "simulation_results/clump/ClumpSimulVVG_use_with_size_dist_interp_s=10_vMax=10000.0_lin_fix_r=1_n=10.csv"

file = "/Users/joshuamiller/Python Files/HCMV_Dose_Infection/simulation_results/clump/ClumpSimulVVG_s=1.0_vMax=1000.0_vMax_c=1000.0_r=True_n=3.csv"
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

try:
    sheet = int(PARAM_DICT_SIMUL['sheet'])
except KeyError:
    sheet = ''

assert simul_name in ['clump', 'comp', 'acc_dam', 'null', 'clump_comp', 'clump_acc_dam', 'var_clump_diam'], "Got: "+simul_name+". 'simul_name' must be one of 'clump', 'comp', 'acc_dam', 'null', 'clump_comp', 'clump_acc_dam', 'var_clump_diam'."

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
               '2021_08_13 ME_mC_epithelial', '2022_11_02_TB_GFP_fib', 
               'use_with_size_distribution', 'use_with_size_dist_interp', 
               '2022_10_27_TB_size_distribution', '2022_10_27_TB_size_dist_interp']

df_data = pd.read_excel('data/Experimental_data_Ed_Josh.xlsx', sheet_name=SHEET_NAMES[sheet])
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

if (('GFP' in SHEET_NAMES[sheet]) or ('use' in SHEET_NAMES[sheet])):
    color  = COLORS['GFP']
    marker = MARKERS['GFP']
elif (('cherry' in SHEET_NAMES[sheet]) or ('mCherry' in SHEET_NAMES[sheet]) or ('mC' in SHEET_NAMES[sheet])):
    color = COLORS['cherry']
    marker = MARKERS['cherry']

LETTERS = ['A', 'B', 'C']

replacement_val = config['PLOTTING']['replacement_val']
band_type = config['PLOTTING']['band_type']
#====================================================================
''' Plot experimental data '''
GEN_WELL_DATA, GEN_CELL_DATA, INF_CELL_DATA, num_zeros = PrepareData(df_data, scale)

ax.scatter(GEN_CELL_DATA, INF_CELL_DATA, facecolors='none', edgecolors=color, marker=marker, alpha=1)

print("AHASDHALJDHAS:JDA:SLDKJA:SLKJDNA:LSKDA:", np.shape(GEN_CELL_DATA), np.shape(INF_CELL_DATA))

g_data, n_data, p_data = PlotFit(ax, GEN_CELL_DATA, INF_CELL_DATA, LOWERS[sheet], UPPERS[sheet], color='b', linestyle='-.')
#====================================================================
''' Plot simulation results '''
GEN_CELL_SIMUL, INF_CELL_SIMUL = PlotSimul(ax, file, band_type, replacement_val)

g_simul, n_simul, p_simul = PlotFit(ax, GEN_CELL_SIMUL, INF_CELL_SIMUL, LOWERS[sheet], UPPERS[sheet])

#====================================================================
''' Formatting '''
legendD = mlines.Line2D([], [], color='b', linestyle='-.', markerfacecolor=color, markeredgecolor='none', markerfacecoloralt='none', marker=marker, alpha=0.6,
                        markersize=10, label = "HCMV-TB (GFP) (data): n = " + str(round(n_data, 3)))
legendS = ''
if (num_simulations == 1):
    legendS = mlines.Line2D([], [], color='y', linestyle='--', markerfacecolor='none', markeredgecolor=color, markerfacecoloralt='none', marker=marker,
                            markersize=10, label= "Simulation: n = " + str(round(n_simul, 3)))
elif (num_simulations > 1):
    legendS = mlines.Line2D([], [], color='r', linestyle='-', label=r'Simul. mean, $\overline{n} = $' + str(round(n_simul, 3)) + " ("+str(num_simulations)+" runs)")
#--------------------------------------------------------------------
BasicFormat(ax)

xMin = 10**-3
xMax = 10**3
yMin = 10**-6
yMax = 1

ax.legend(handles = [legendD, legendS], loc='upper left', prop={'size': 8})

PlotText(ax, PARAM_DICT_SIMUL, xMin, xMax, yMin, yMax)

if ('clump' in simul_name):
    if (simul_name == 'clump'):
        ax.set_title(SHEET_NAMES[sheet] + ' | Clumping')
    elif (simul_name == 'clump_comp'):
        ax.set_title(SHEET_NAMES[sheet] + ' | Clumping + Compensation')
    elif (simul_name == 'clump_acc_dam'):
        ax.set_title(SHEET_NAMES[sheet] + ' | Clumping + Acc. Damage')
    elif (simul_name == 'var_clump_diam'):
        ax.set_title(SHEET_NAMES[sheet] + ' | Clumping + Var. Clump Diam.')
elif (simul_name == 'acc_dam'):
    ax.set_title(SHEET_NAMES[sheet] + ' | Acc. Damage')
elif (simul_name == 'comp'):
    ax.set_title(SHEET_NAMES[sheet] + ' | Compensation')
elif (simul_name == 'null'):
    ax.set_title(SHEET_NAMES[sheet] + ' | Null')

print("____________________________________________________________")
print(" >> cells", cell_count)
print("____________________________________________________________")

ax.text(.6 * xMin, 2.2 * yMax, '', fontsize=16, fontweight='bold', va='top', ha='right')
#====================================================================
''' Save figure '''
filename = MakeFilename(PARAM_DICT_SIMUL, SHEET_NAMES[sheet])
fig.savefig(os.path.join(os.path.join(os.getcwd(), 'figs'), filename+".pdf"), bbox_inches = 'tight', pad_inches = 0) # Save figure in the new directory

#====================================================================
#====================================================================
#====================================================================
#====================================================================
#====================================================================
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
xMin = 10**-4
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

PlotText(ax2, PARAM_DICT_SIMUL, xMin, xMax, yMin, yMax)

if ('clump' in simul_name):
    if (simul_name == 'clump'):
        ax2.set_title(SHEET_NAMES[sheet] + ' | Clumping')
    elif (simul_name == 'clump_comp'):
        ax2.set_title(SHEET_NAMES[sheet] + ' | Clumping + Compensation')
    elif (simul_name == 'clump_acc_dam'):
        ax2.set_title(SHEET_NAMES[sheet] + ' | Clumping + Acc. Damage')
elif (simul_name == 'acc_dam'):
    ax2.set_title(SHEET_NAMES[sheet] + ' | Acc. Damage')
elif (simul_name == 'comp'):
    ax2.set_title(SHEET_NAMES[sheet] + ' | Compensation')
elif (simul_name == 'null'):
    ax2.set_title(SHEET_NAMES[sheet] + ' | Null')
#====================================================================
''' Save figure '''
filename = MakeFilename(PARAM_DICT_SIMUL, SHEET_NAMES[sheet])
fig2.savefig(os.path.join(os.path.join(os.getcwd(), 'figs/interaction_figs'), "INTER_"+filename+".pdf"), bbox_inches = 'tight', pad_inches = 0) # Save figure in the new directory
plt.show()