import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import pandas as pd
import yaml
import sys
import os
import json

sys.path.append(os.getcwd())

from misc.misc_utils import Trapezoid, FlattenMeans, ExtractParams
from plotting_utils.utils import PlotSimul, PlotFit, PlotText
from simulation_utils.utils import PrepareData
#====================================================================
''' Files '''
# file_simul      = "simulation_results/clump/ClumpSimulVVG_use_with_size_dist_interp_s=10_vMax=10000.0_sp_fix_r=1_n=10.csv"
# file_simul_size = "simulation_results/clump/clump_information/ClumpSimulVVG_use_with_size_dist_interp_s=10_vMax=10000.0_sp_fix_r=1_CLUMP.json"

file_simul      = "simulation_results/clump/ClumpSimulVVG_use_with_size_distribution_s=10_vMax=10000.0_sp_fix_r=1_n=10.csv"
file_simul_size = "simulation_results/clump/clump_information/ClumpSimulVVG_use_with_size_distribution_s=10_vMax=10000.0_sp_fix_r=1_CLUMP.json"

file_data = "data/Experimental_data_Ed_Josh.xlsx"

SHEET_NAMES = ['2021_10_05 TB_GFP_epithelial', '2020_07_02 ME_GFP_fibroblast', 
               '2020_05_29 TR_GFP_fibroblast', '2021_07_13 GFP_TB_fibroblast', 
               '2020_08_12 TB_GFP_fibroblast', '2020_09_14 TR_GFP_epithelial',
               '2021_08_13 ME_mC_epithelial', '2022_11_02_TB_GFP_fib', 
               'use_with_size_distribution', 'use_with_size_dist_interp', 
               '2022_10_27_TB_size_distribution', '2022_10_27_TB_size_dist_interp']

assert "clump" or "Clump" in file_simul, "Simulation type must be 'clump', 'clump_comp', or 'clump_acc_dam'."
#====================================================================
''' Key parameters '''
flatten_means = True
cutoff = 91.28

yMin_1 = .35
yMax_1 = .75
yMin_2 = .65
yMax_2 = .95
#====================================================================
with open('config.yml', 'r') as c:
    config = yaml.load(c, Loader=yaml.FullLoader)

sheet = config['SIMULATION_PARAMETERS']['sheet']

mean = config['CLUMP_PARAMETERS']['mean']
lb   = config['CLUMP_PARAMETERS']['lb']
ub   = config['CLUMP_PARAMETERS']['ub']

LOWERS = config['PLOTTING']['LOWERS'] # [2, 20, 0, 6, 0, 0]
UPPERS = config['PLOTTING']['UPPERS'] # [36, -15, -1, 36, -1, -1]

print(LOWERS)
print(UPPERS)

MARKERS = config['PLOTTING']['markers_dict']#['o', '^', 's', 'D']
COLORS  = config['PLOTTING']['colors_dict']

if (('GFP' in SHEET_NAMES[sheet]) or (SHEET_NAMES[sheet] == 'use_with_size_distribution')):
    color  = COLORS['GFP']
    marker = MARKERS['GFP']
elif (('cherry' in SHEET_NAMES[sheet]) or ('mCherry' in SHEET_NAMES[sheet]) or ('mC' in SHEET_NAMES[sheet])):
    color = COLORS['cherry']
    marker = MARKERS['cherry']
else:
    color = 'green'
    marker = 's'

LETTERS = ['A', 'B', 'C']
replacement_val = config['PLOTTING']['replacement_val']
band_type = config['PLOTTING']['band_type']
#====================================================================
''' Get simulation data '''
df_simul = pd.read_csv(file_simul)
#====================================================================
''' Get simulated clump size data '''
CLUMP_DICT = {}
with open(file_simul_size, "r") as f:
    CLUMP_DICT = json.load(f)

print('==============================')
#print(CLUMP_DICT)
print('==============================')
#====================================================================
''' Get experimental data '''

df_data = pd.read_excel(file_data, sheet_name=SHEET_NAMES[sheet])

df_sizes_data = pd.read_excel(file_data, sheet_name='2022_10_27_TB_size_distribution')
diameter = df_sizes_data.pop('d.nm')
diameter = np.asarray(diameter)
#====================================================================
''' Average '''
means_df = df_sizes_data.filter(like='Mean')
stds_df  = df_sizes_data.filter(like='SD')

mean_of_means = means_df.mean(axis=1)
mean_of_means = np.asarray(mean_of_means)

if flatten_means:
    mean_of_means = FlattenMeans(diameter, mean_of_means, cutoff)

mean_of_means /= Trapezoid(diameter, mean_of_means) # Normalize by dividing by area 

#====================================================================
''' Create subplot array '''
fig = plt.figure(figsize=(10, 6), dpi=100, constrained_layout=False)
gs  = fig.add_gridspec(9, 5)
ax0 = fig.add_subplot(gs[:, :-2])
ax1 = fig.add_subplot(gs[0:3, -2:])
ax2 = fig.add_subplot(gs[3:, -2:])
fig.subplots_adjust(wspace=1, hspace=20)
#====================================================================
''' Get parameters '''
PARAM_DICT_SIMUL = ExtractParams(df_simul)

simul_name      = PARAM_DICT_SIMUL['simul_name']
scale = PARAM_DICT_SIMUL['scale']
i_mean = PARAM_DICT_SIMUL['muG'] #0
i_stdev = PARAM_DICT_SIMUL['sigmaG'] #1
r_mean = PARAM_DICT_SIMUL['muR']
r_stdev = PARAM_DICT_SIMUL['sigmaR']
gamma = PARAM_DICT_SIMUL['gamma']
vMax = PARAM_DICT_SIMUL['vMax']
scheme = PARAM_DICT_SIMUL['scheme']
distribution = PARAM_DICT_SIMUL['distribution']

num_simulations = 0
for key in df_simul.keys():
    if ("IU run" in key):
        num_simulations += 1
#--------------------------------------------------------------------
PARAM_DICT_DATA = ExtractParams(df_data)
cell_count = int(PARAM_DICT_DATA['cell_count'] / scale)
#====================================================================
''' Plot experimental data '''
GEN_WELL_DATA, GEN_CELL_DATA, INF_CELL_DATA, num_zeros = PrepareData(df_data, scale)

ax0.scatter(GEN_CELL_DATA, INF_CELL_DATA, facecolors='none', edgecolors=color, marker=marker, alpha=1)

print("AHASDHALJDHAS:JDA:SLDKJA:SLKJDNA:LSKDA:", np.shape(GEN_CELL_DATA), np.shape(INF_CELL_DATA))

g_data, n_data, p_data = PlotFit(ax0, GEN_CELL_DATA, INF_CELL_DATA, LOWERS[sheet], UPPERS[sheet], color='b', linestyle='-.',
                                 yMin_1=yMin_1, yMax_1=yMax_1, yMin_2=yMin_2, yMax_2=yMax_2,
                                 n_fits=100)

#--------------------------------------------------------------------
''' Plot simulation results '''
GEN_WELL_SIMUL, GEN_CELL_SIMUL, INF_WELL_SIMUL, INF_CELL_SIMUL = PlotSimul(ax0, file_simul, band_type, replacement_val, return_gen_well=True)

g_simul, n_simul, p_simul = PlotFit(ax0, GEN_CELL_SIMUL, INF_CELL_SIMUL, LOWERS[sheet], UPPERS[sheet],
                                    yMin_1=yMin_1, yMax_1=yMax_1, yMin_2=yMin_2, yMax_2=yMax_2)
#====================================================================
ax1.scatter(diameter, mean_of_means, facecolors='None', marker='o', edgecolor='k', label='Mean of all samples')
#====================================================================
''' Plot pdf of clump distribution '''
CLUMP_MARKERS = ['s', '^', 'o']

NUM_CLUMPS_2  = CLUMP_DICT[str(GEN_WELL_SIMUL[len(GEN_WELL_SIMUL) - 1])]
CLUMP_SIZES_2 = np.arange(1, len(NUM_CLUMPS_2)+1)

ax2.vlines(CLUMP_SIZES_2, 0, NUM_CLUMPS_2, colors='y', lw=3, alpha=0.7, zorder=0)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NUM_CLUMPS_1  = CLUMP_DICT[str(GEN_WELL_SIMUL[int(len(GEN_WELL_SIMUL) / 2)])]
CLUMP_SIZES_1 = np.arange(1, np.shape(NUM_CLUMPS_1)[0]+1)

ax2.vlines(CLUMP_SIZES_1, 0, NUM_CLUMPS_1, colors='r', lw=3, alpha=0.5, zorder=1)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NUM_CLUMPS_0 = CLUMP_DICT[str(GEN_WELL_SIMUL[0])]
CLUMP_SIZES_0 = np.arange(1, np.shape(NUM_CLUMPS_0)[0]+1)

ax2.vlines(CLUMP_SIZES_0, 0, NUM_CLUMPS_0, colors='b', lw=3, alpha=0.5, zorder=2)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ax2.scatter(CLUMP_SIZES_2, NUM_CLUMPS_2, s = 50, color='yellow', edgecolors='black', marker=CLUMP_MARKERS[2], label = 'Gen./well = ' + str(GEN_WELL_SIMUL[len(GEN_WELL_SIMUL) - 1]))
ax2.scatter(CLUMP_SIZES_1, NUM_CLUMPS_1, s = 50, color='red', edgecolors ='black', marker=CLUMP_MARKERS[1], label='Gen./well = ' + str(GEN_WELL_SIMUL[int(len(GEN_WELL_SIMUL) / 2)]))
ax2.scatter(CLUMP_SIZES_0, NUM_CLUMPS_0, s = 50, color='blue', edgecolors='black', marker=CLUMP_MARKERS[0], label='Gen./well = ' + str(GEN_WELL_SIMUL[0]))
#====================================================================
''' Formatting '''
#--------------------------------------------------------------------
xMin = 10**-3
xMax = 10**3
yMin = 10**-6
yMax = 1

legendD = mlines.Line2D([], [], color='b', linestyle='-.', markerfacecolor=color, markeredgecolor='none', markerfacecoloralt='none', marker=marker,
                          markersize=10, label = "HCMV-TB (GFP) (data): n = " + str(round(n_data, 3)) + ' (p='+str(round(p_data, 5)) + ')')
legendS = ''
if (num_simulations == 1):
    legendS = mlines.Line2D([], [], color='y', linestyle='--', markerfacecolor='none', markeredgecolor=color, markerfacecoloralt='none', marker=marker,
                            markersize=10, label= "Simulation: n = " + str(round(n_simul, 3)) + ' (p='+str(round(p_simul, 5)) + ')')
elif (num_simulations > 1):
    legendS = mlines.Line2D([], [], color='k', linestyle='-', label=r'Simul. mean, $\overline{n} = $' + str(round(n_simul, 3)) + ' (p='+str(round(p_simul, 5)) + ')' + " ("+str(num_simulations)+" runs)")

ax0.legend(handles = [legendD, legendS], loc='upper left', prop={'size': 8})
ax0.plot(np.linspace(xMin,xMax), np.linspace(yMin, yMax), 'k--', linewidth = 1)
ax0.set_xlabel('Genomes/cell')
ax0.set_ylabel('Infections/cell')
ax0.set_xscale('log')
ax0.set_yscale('log')
ax0.set_xlim(xMin, xMax)
ax0.set_ylim(yMin, yMax)

PlotText(ax0, PARAM_DICT_SIMUL, xMin, xMax, yMin, yMax)

ax0.text(.6 * xMin, 2.2 * yMax, 'A', fontsize=16, fontweight='bold', va='top', ha='right')

ax0.set_title(SHEET_NAMES[sheet] + ' | '+PARAM_DICT_SIMUL['simul_name'])
#--------------------------------------------------------------------
ax1.legend(loc='upper left')
ax1.set_ylim(0, 1.5 * max(mean_of_means))
ax1.set_xscale('log')
ax1.set_xlabel("Diameter (nm)")
ax1.set_ylabel("Norm. frequency")
ax1.text(0.2, 1.7 * max(mean_of_means), 'B', fontsize=16, fontweight='bold', va='top', ha='right')
#--------------------------------------------------------------------
yMin2 = 7 * 10**-1
yMax2 = 1000 * max(NUM_CLUMPS_2)
ax2.set_ylim(yMin2, yMax2)
ax2.set_yscale('log')
ax2.set_xlabel('Clump size (num. virions in clump)')
ax2.set_ylabel('Frequency')
ax2.legend(loc='upper right', prop = {'size' : 9})
ax2.grid()
ax2.minorticks_on()
#ax2.grid(which='minor', color='#999999', linestyle='-', alpha=0.2)
ax2.text(-1 * (np.shape(CLUMP_SIZES_2)[0] / 100), 3 * yMax2, 'C', fontsize=16, fontweight='bold', va='top', ha='right')
#====================================================================
plt.show()
#====================================================================
''' Save figure '''
scheme_short = ''
if (scheme=='linear'):
    scheme_short='lin'
elif (scheme=='regular_polygon'):
    scheme_short='poly'
elif (scheme=='sphere_packing'):
    scheme_short='sp'

dist_short = ''
if (distribution=='normal'):
    dist_short = 'norm'
elif (distribution=='uniform'):
    dist_short = 'uni'
elif (distribution=='fixed'):
    dist_short = 'fix'

filename = ''
if (simul_name == 'clump'):
    filename = "ClumpSimul_"+SHEET_NAMES[sheet]+"_s="+str(scale)+"_vMax="+str(vMax)+"_"+scheme_short+"_"+dist_short # Specify filename
elif (simul_name == 'clump_comp'):
    filename = "ClumpCompSimul_"+SHEET_NAMES[sheet]+"_s="+str(scale)+"_vMax="+str(vMax)+"_k="+str(PARAM_DICT_SIMUL['kappa'])+"_"+scheme_short+"_"+dist_short # Specify filename
elif (simul_name == 'clump_acc_dam'):
    filename = "ClumpAccDamSimul_"+SHEET_NAMES[sheet]+"_s="+str(scale)+"_vMax="+str(vMax)+"_b="+str(PARAM_DICT_SIMUL['beta'])+"_"+scheme_short+"_"+dist_short # Specify filename
elif (simul_name == 'var_clump_diam'):
    filename = "VarClumpDiam_"+SHEET_NAMES[sheet]+"_s="+str(scale)+"_vMax="+str(vMax)+scheme_short+"_"+dist_short # Specify filename

fig.savefig(os.path.join(os.path.join(os.getcwd(), 'figs'), filename+".pdf"), bbox_inches = 'tight', pad_inches = 0) # Save figure in the new directory