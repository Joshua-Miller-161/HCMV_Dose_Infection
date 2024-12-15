import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import pandas as pd
from scipy.stats import norm
import yaml
import sys
import os
import json

sys.path.append(os.getcwd())

from misc.misc_utils import ExtractParams
from simulation_utils.utils import PrepareData, PrepareParameters
from plotting_utils.utils import BasicFormat, PlotSimul, PlotFit, PlotText, add_subplot_axes, PlotVirionsInClumpFreq
#====================================================================
''' Make subplot array '''
fig, ax = plt.subplots(2, 2, figsize=(8, 8), dpi=100)
fig.subplots_adjust(wspace=.3, hspace=.2)

ax_data          = ax[0][0] 
ax_clump_acc_dam = ax[1][0]
ax_clump         = ax[0][1]
ax_clump_comp    = ax[1][1]

x_text  = 1.2*10**-3
y_text  = 5*10**-1
x_label = 1.2*10**-4
y_label = 1.5
yMin_1 = .45
yMax_1 = .70
yMin_2 = .65
yMax_2 = .95
x_legend = .01
y_legend = .95
xMax_clump = 70
yMax_clump = 10000000
#====================================================================
''' Files '''

SHEET_NAMES = ['2021_10_05 TB_GFP_epithelial', '2020_07_02 ME_GFP_fibroblast', 
               '2020_05_29 TR_GFP_fibroblast', '2021_07_13 GFP_TB_fibroblast', 
               '2020_08_12 TB_GFP_fibroblast', '2020_09_14 TR_GFP_epithelial',
               '2021_08_13 ME_mC_epithelial', '2022_11_02_TB_GFP_fib', 
               'use_with_size_distribution', 'use_with_size_dist_interp', 
               '2022_10_27_TB_size_distribution', '2022_10_27_TB_size_dist_interp']

data_file               = "data/Experimental_data_Ed_Josh.xlsx"
clump_file              = "simulation_results/clump/ClumpSimulVVG_use_with_size_dist_interp_s=10_vMax=10000.0_sp_fix_r=1_n=10.csv"
clump_info_file         = "simulation_results/clump/clump_information/ClumpSimulVVG_use_with_size_dist_interp_s=10_vMax=10000.0_sp_fix_r=1_CLUMP.json"
clump_comp_file         = "simulation_results/clump_comp/ClumpCompSimulVVG_use_with_size_dist_interp_s=10_vMax=10000.0_k=0.0_sp_fix_r=1_n=10.csv"
clump_comp_info_file    = "simulation_results/clump_comp/clump_information/ClumpCompSimulVVG_use_with_size_dist_interp_s=10_vMax=10000.0_k=0.0_sp_fix_r=1_CLUMP.json"
clump_acc_dam_file      = "simulation_results/clump_acc_dam/ClumpAccDamSimulVVG_use_with_size_dist_interp_s=10_vMax=10000.0_b=-2.0_sp_fix_r=1_n=10.csv"
clump_acc_dam_info_file = "simulation_results/clump_acc_dam/clump_information/ClumpAccDamSimulVVG_use_with_size_dist_interp_s=10_vMax=10000.0_b=-2.0_sp_fix_r=1_CLUMP.json"
#====================================================================
with open('config.yml', 'r') as c:
    config = yaml.load(c, Loader=yaml.FullLoader)
sheet  = config['SIMULATION_PARAMETERS']['sheet']

mu_res  = config['CELL_PARAMETERS']['mu_res']  # Mean cell resistivity
std_res = config['CELL_PARAMETERS']['std_res'] # Standard deviation of resisitivty

mu_inf  = config['GFP_VIRUS_PARAMETERS']['mu_inf']  # Mean cell resistivity
std_inf = config['GFP_VIRUS_PARAMETERS']['std_inf'] # Standard deviation of resisitivty

LOWERS = config['PLOTTING']['LOWERS'] # [2, 20, 0, 6, 0, 0]
UPPERS = config['PLOTTING']['UPPERS'] # [36, -15, -1, 36, -1, -1]

print(LOWERS)
print(UPPERS)

MARKERS = config['PLOTTING']['markers_dict']#['o', '^', 's', 'D']
COLORS = config['PLOTTING']['colors_dict']

if ('GFP' in SHEET_NAMES[sheet] or SHEET_NAMES[sheet] == 'use_with_size_distribution'):
    color  = COLORS['GFP']
    marker = MARKERS['GFP']
elif (('cherry' in SHEET_NAMES[sheet]) or ('mCherry' in SHEET_NAMES[sheet]) or ('mC' in SHEET_NAMES[sheet])):
    color = COLORS['cherry']
    marker = MARKERS['cherry']

LETTERS = ['A', 'B', 'C']

replacement_val = config['PLOTTING']['replacement_val']
band_type = config['PLOTTING']['band_type']

marker = '^'
marker_color = 'gray'
#====================================================================
#====================================================================
#====================================================================
''' Plot experimental data '''
df_data = pd.read_excel(data_file, sheet_name='use_with_size_dist_interp')

PARAM_DICT_DATA = ExtractParams(df_data)
cell_count = int(PARAM_DICT_DATA['cell_count'])

GEN_WELL_DATA, GEN_CELL_DATA, INF_CELL_DATA, num_zeros = PrepareData(df_data, 1)

#ax_data.scatter(GEN_CELL_DATA, INF_CELL_DATA, facecolors='none', edgecolors=color, marker=marker, alpha=1)
ax_data.scatter(GEN_CELL_DATA, INF_CELL_DATA, facecolors='grey', edgecolors='none', marker='^', alpha=1)

g_data, n_data, p_data = PlotFit(ax_data, GEN_CELL_DATA, INF_CELL_DATA, LOWERS[sheet], UPPERS[sheet], color='r', linestyle='-',
                                 plot_vlines=False, plot_fit_markers=True, s=40, n_fits=100)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
BasicFormat(ax_data)

legendData = mlines.Line2D([], [], color='none', linestyle='-',
                           markerfacecolor='gray', markeredgecolor='none', markerfacecoloralt='none', marker='^', markersize=10,
                           label = "Observed data")
legendN = mlines.Line2D([], [], color='r', linestyle='-',
                        label='n=' + str(round(n_data, 3)))

ax_data.legend(handles = [legendData, legendN], loc='upper left', prop={'size': 7}, bbox_to_anchor=(x_legend, y_legend))

ax_data.text(x_text, y_text, 'TB GFP fibroblast', fontsize=10, fontweight='bold')
ax_data.text(x_label, y_label, 'A', fontsize=16, fontweight='bold')
#====================================================================
''' Plot particle size data '''
size_df = pd.read_excel("data/Experimental_data_Ed_Josh.xlsx", sheet_name='2022_10_27_TB_size_distribution')

DIAMETERS    = size_df.pop('d.nm')
DIAM_MEANS_DICT = {}
for key in size_df.keys():
    if ('Mean' in key):
        strs = key.split('_')
        DIAM_MEANS_DICT[strs[1]] = np.asarray(size_df.loc[:, key].values)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
subpos    = [.6, 0.13, 0.35, 0.35] # [x,y,width,height]
ax_data_inset  = add_subplot_axes(ax_data, subpos)

#shapes = ['o', '*', 'D', '8', '^', 's', 'H', 'P']
#colors = ['r', 'g', 'b', 'm', 'c', 'y', 'k', 'pink']
colors = ['goldenrod', 'red', 'blue']
linestyles = ['-', '-.', ':']

keys_ = list(DIAM_MEANS_DICT.keys())
idx   = [0, 3, len(keys_)-1]
for i in range(len(idx)):
    print(keys_[i])
    #ax_data_inset.scatter(DIAMETERS, DIAM_MEANS_DICT[keys_[i]], facecolors='None', marker=shapes[i], edgecolor=colors[i], label='GEN./WELL: '+str(keys_[i]))
    ax_data_inset.plot(DIAMETERS, DIAM_MEANS_DICT[keys_[idx[i]]], color=colors[i], linestyle=linestyles[i], label='GEN./cell: '+str(round(float(keys_[idx[i]]) / cell_count, 3)))

ax_data_inset.axvline(230, 0, 1, color='black', linestyle='--', label=r'$230\mu m$')

ax_data_inset.set_xscale('log')
ax_data_inset.legend(loc='upper left', prop={'size': 4})
ax_data_inset.set_xlabel('Diameter (nm)', fontsize=6)
ax_data_inset.set_ylabel('Percent in the sample', fontsize=6)
#====================================================================
#====================================================================
#====================================================================
''' Plot clump '''
PARAM_DICT_SIMUL = ExtractParams(pd.read_csv(clump_file))
cell_count = float(PARAM_DICT_SIMUL['cell_count'])
scale      = float(PARAM_DICT_SIMUL['scale'])
num_simulations = PARAM_DICT_SIMUL['num_simulations']

GEN_WELL_SIMUL, GEN_CELL_SIMUL, INF_WELL_SIMUL, INF_CELL_SIMUL = PlotSimul(ax_clump, clump_file, band_type, replacement_val, 
                                                                           scatter=True, marker=marker, color=marker_color, return_gen_well=True)

g_clump, n_clump, p_clump = PlotFit(ax_clump, GEN_CELL_SIMUL, INF_CELL_SIMUL, LOWERS[sheet], UPPERS[sheet],
                                    plot_fit_markers=True, s=50, marker_color='green', marker=marker, 
                                    n_fits=100)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
BasicFormat(ax_clump)

ax_clump.text(x_text, y_text, 'Clump', fontsize=10, fontweight='bold')
ax_clump.text(x_label, y_label, 'B', fontsize=16, fontweight='bold')

legendSimul = mlines.Line2D([], [], color='k', linestyle='-', 
                            markerfacecolor=marker_color, markeredgecolor='none', markerfacecoloralt='none', marker=marker, markersize=10,
                            label = "Simulation ("+str(int(num_simulations))+" runs)")
legendN = mlines.Line2D([], [], color='r', linestyle='-',
                        label='n=' + str(round(n_clump, 3)))

ax_clump.legend(handles = [legendSimul, legendN], loc='upper left', prop={'size': 7}, bbox_to_anchor=(x_legend, y_legend))
#====================================================================
''' Plot virions in clump '''
subpos    = [.6, 0.13, 0.35, 0.35] # [x,y,width,height]
ax_clump_inset  = add_subplot_axes(ax_clump, subpos)

CLUMP_DICT = {}
with open(clump_info_file, "r") as f:
    CLUMP_DICT = json.load(f)
PlotVirionsInClumpFreq(ax_clump_inset, CLUMP_DICT, GEN_WELL_SIMUL,
                       scatter=False, cell_count=cell_count, round_digits=3)

ax_clump_inset.set_xlim(0, xMax_clump)
ax_clump_inset.set_ylim(.5, yMax_clump)
ax_clump_inset.set_xlabel('Virions in clump', fontsize=6)
ax_clump_inset.set_ylabel('Frequency', fontsize=6)
ax_clump_inset.set_yscale('log')
#====================================================================
#====================================================================
#====================================================================
''' Plot Accrued damage '''
PARAM_DICT_SIMUL = ExtractParams(pd.read_csv(clump_acc_dam_file))
cell_count = float(PARAM_DICT_SIMUL['cell_count'])
scale      = float(PARAM_DICT_SIMUL['scale'])
num_simulations = PARAM_DICT_SIMUL['num_simulations']

GEN_WELL_SIMUL, GEN_CELL_SIMUL, INF_WELL_SIMUL, INF_CELL_SIMUL = PlotSimul(ax_clump_acc_dam, clump_acc_dam_file, band_type, replacement_val, 
                                                                           scatter=True, color=marker_color, marker=marker, return_gen_well=True)

g_clump_acc_dam, n_clump_acc_dam, p_clump_acc_dam = PlotFit(ax_clump_acc_dam, GEN_CELL_SIMUL, INF_CELL_SIMUL, LOWERS[sheet], UPPERS[sheet],
                                                            plot_fit_markers=True, s=50, marker_color='green', marker=marker, 
                                                            n_fits=100)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
BasicFormat(ax_clump_acc_dam)

ax_clump_acc_dam.text(x_text, y_text, 'Clump + Acc. Damage', fontsize=10, fontweight='bold')
ax_clump_acc_dam.text(x_label, y_label, 'C', fontsize=16, fontweight='bold')

legendSimul = mlines.Line2D([], [], 
                            markerfacecolor=marker_color, markeredgecolor='none', markerfacecoloralt='none', marker=marker, markersize=10,
                            color='k', linestyle='-', 
                            label = "Simulation ("+str(int(num_simulations))+" runs)")
legendN = mlines.Line2D([], [], color='r', linestyle='-',
                        label='n=' + str(round(n_clump_acc_dam, 3)))

ax_clump_acc_dam.legend(handles = [legendSimul, legendN], loc='upper left', prop={'size': 7}, bbox_to_anchor=(x_legend, y_legend))
#====================================================================
''' Plot virions in clump '''
subpos    = [.6, 0.13, 0.35, 0.35] # [x,y,width,height]
ax_clump_acc_dam_inset  = add_subplot_axes(ax_clump_acc_dam, subpos)

CLUMP_DICT = {}
with open(clump_acc_dam_info_file, "r") as f:
    CLUMP_DICT = json.load(f)
PlotVirionsInClumpFreq(ax_clump_acc_dam_inset, CLUMP_DICT, GEN_WELL_SIMUL,
                       scatter=False, cell_count=cell_count, round_digits=3)

ax_clump_acc_dam_inset.set_xlim(0, xMax_clump)
ax_clump_acc_dam_inset.set_ylim(.5, yMax_clump)
ax_clump_acc_dam_inset.set_xlabel('Virions in clump', fontsize=6)
ax_clump_acc_dam_inset.set_ylabel('Frequency', fontsize=6)
ax_clump_acc_dam_inset.set_yscale('log')
#====================================================================
#====================================================================
#====================================================================
''' Plot compensation '''
PARAM_DICT_SIMUL = ExtractParams(pd.read_csv(clump_comp_file))
cell_count = float(PARAM_DICT_SIMUL['cell_count'])
scale      = float(PARAM_DICT_SIMUL['scale'])
num_simulations = PARAM_DICT_SIMUL['num_simulations']

GEN_WELL_SIMUL, GEN_CELL_SIMUL, INF_WELL_SIMUL, INF_CELL_SIMUL = PlotSimul(ax_clump_comp, clump_comp_file, band_type, replacement_val, 
                                                                           scatter=True, color=marker_color, marker=marker, return_gen_well=True)

g_clump_comp, n_clump_comp, p_clump_comp = PlotFit(ax_clump_comp, GEN_CELL_SIMUL, INF_CELL_SIMUL, LOWERS[sheet], UPPERS[sheet],
                                                   plot_fit_markers=True, s=50, marker_color='green', marker=marker, 
                                                   n_fits=100)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
BasicFormat(ax_clump_comp)

ax_clump_comp.text(x_text, y_text, 'Clump + compensation', fontsize=10, fontweight='bold')
ax_clump_comp.text(x_label, y_label, 'D', fontsize=16, fontweight='bold')

legendSimul = mlines.Line2D([], [], color='k', linestyle='-', 
                            markerfacecolor=marker_color, markeredgecolor='none', markerfacecoloralt='none', marker=marker, markersize=10,
                            label = "Simulation ("+str(int(num_simulations))+" runs)")
legendN = mlines.Line2D([], [], color='r', linestyle='-',
                        label='n=' + str(round(n_clump_comp, 3)))

ax_clump_comp.legend(handles = [legendSimul, legendN], loc='upper left', prop={'size': 7}, bbox_to_anchor=(x_legend, y_legend))
#====================================================================
''' Plot virions in clump '''
subpos    = [.6, 0.13, 0.35, 0.35] # [x,y,width,height]
ax_clump_comp_inset  = add_subplot_axes(ax_clump_comp, subpos)

CLUMP_DICT = {}
with open(clump_comp_info_file, "r") as f:
    CLUMP_DICT = json.load(f)
PlotVirionsInClumpFreq(ax_clump_comp_inset, CLUMP_DICT, GEN_WELL_SIMUL,
                       scatter=False, cell_count=cell_count, round_digits=3)

ax_clump_comp_inset.set_xlim(0, xMax_clump)
ax_clump_comp_inset.set_ylim(.5, yMax_clump)
ax_clump_comp_inset.set_xlabel('Virions in clump', fontsize=6)
ax_clump_comp_inset.set_ylabel('Frequency', fontsize=6)
ax_clump_comp_inset.set_yscale('log')
#====================================================================
plt.show()
#====================================================================
#fig.savefig(os.path.join(os.path.join(os.getcwd(), 'figs'), 'VVG_four_panel_clump_same_params.pdf'), bbox_inches='tight', pad_inches=0)
#====================================================================