import sys
sys.dont_write_bytecode = True
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib.ticker import MaxNLocator
import pandas as pd
import yaml
import os
import json

sys.path.append(os.getcwd())

from misc_utils import MakeFilename
from plot_utils import PlotText, PlotSimul
from misc.misc_utils import ExtractParams
from plotting_utils.utils import PlotFit
#====================================================================
''' Files '''

file_simul      = "/Users/joshuamiller/Python Files/HCMV_Dose_Infection/simulation_results/clump/ClumpSimulVVG_s=1.0_vMax=1000.0_f=poi_vMax_c=1000.0_r=1_n=100.csv"
file_simul_size = "/Users/joshuamiller/Python Files/HCMV_Dose_Infection/simulation_results/clump/clump_information/ClumpSimulVVG_s=1.0_vMax=1000.0_f=poi_vMax_c=1000.0_r=1_n=100_CLUMP.json"

assert "clump" or "Clump" in file_simul, "Simulation type must be 'clump', 'clump_comp', or 'clump_acc_dam'."

#====================================================================
with open('clump_sandbox/sandbox_config.yml', 'r') as c:
    config = yaml.load(c, Loader=yaml.FullLoader)

#lower = config['PLOTTING']['LOWERS'] # [2, 20, 0, 6, 0, 0]
#upper = config['PLOTTING']['UPPERS'] # [36, -15, -1, 36, -1, -1]
lb_lo = 0
ub_lo = 9
lb_hi = 10
ub_hi = 17

marker = '^'
color  = 'green'

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
#====================================================================
''' Get parameters '''

PARAM_DICT_SIMUL = ExtractParams(df_simul)

print(PARAM_DICT_SIMUL)

simul_name   = PARAM_DICT_SIMUL['simul_name']
scale        = PARAM_DICT_SIMUL['scale']
cell_count   = PARAM_DICT_SIMUL['cell_count']
i_mean       = PARAM_DICT_SIMUL['muG'] #0
i_stdev      = PARAM_DICT_SIMUL['sigmaG'] #1
r_mean       = PARAM_DICT_SIMUL['muR']
r_stdev      = PARAM_DICT_SIMUL['sigmaR']
gamma        = PARAM_DICT_SIMUL['gamma']
vMax         = PARAM_DICT_SIMUL['vMax']
distribution = PARAM_DICT_SIMUL['distribution']

try:
    vMax_c = PARAM_DICT_SIMUL['vMax_c']
    print(" >> Got vMax_c", vMax_c)
except KeyError:
    vMax_c = ''

try: 
    fixed = bool(PARAM_DICT_SIMUL['fixed_mean'])
except KeyError:
    fixed = True

try:
    mean = float(PARAM_DICT_SIMUL['mean'])
except KeyError:
    mean = ''

num_simulations = 0
for key in df_simul.keys():
    if ("IU run" in key):
        num_simulations += 1
#====================================================================












''' Plot simulation results '''

fig = plt.figure(figsize=(10, 6), dpi=100, constrained_layout=False)
gs  = fig.add_gridspec(9, 5)
ax0 = fig.add_subplot(gs[:, :-2])
ax1 = fig.add_subplot(gs[0:3, -2:])
ax2 = fig.add_subplot(gs[3:, -2:])
fig.subplots_adjust(wspace=1, hspace=20)
#--------------------------------------------------------------------
''' Plot simulations '''

GEN_WELL_SIMUL, GEN_CELL_SIMUL, INF_WELL_SIMUL, INF_CELL_SIMUL = PlotSimul(ax0, file_simul, band_type, replacement_val, return_gen_well=True, s=100)

g_simul, n_simul, p_simul = PlotFit(ax0, GEN_CELL_SIMUL, INF_CELL_SIMUL, lb_lo, ub_lo, plot_fit_markers=True, marker_color='pink', s=50, linewidth=2, linestyle='-.')

g_simul_hi, n_simul_hi, p_simul_hi = PlotFit(ax0, GEN_CELL_SIMUL, INF_CELL_SIMUL, lb_hi, ub_hi, plot_fit_markers=True, marker_color='skyblue', s=50, linewidth=2, linestyle='--', color='blue')
#--------------------------------------------------------------------
''' Plot the mean of the distribution '''

MEANS = np.empty_like(GEN_WELL_SIMUL, float)

if (fixed == True and not mean==''):
    MEANS = np.ones_like(MEANS) * float(mean)

else:
    if (distribution in ['geometric', '1inf-geo']):
        for i in range(GEN_WELL_SIMUL.shape[0]):
            MEANS[i] = vMax_c / GEN_WELL_SIMUL[i]
            if (MEANS[i] > 0.99999):
                MEANS[i] = 0.99999
    
    elif (distribution == 'poisson'):
        for i in range(GEN_WELL_SIMUL.shape[0]):
            MEANS[i] = GEN_WELL_SIMUL[i] / vMax_c

print(MEANS)
ax1.scatter(GEN_CELL_SIMUL, MEANS, facecolors='None', marker='o', edgecolor='k')
ax1.scatter(GEN_CELL_SIMUL[lb_lo:ub_lo+1], MEANS[lb_lo:ub_lo+1], facecolors='red', marker='o', s=5, label='n region')
ax1.scatter(GEN_CELL_SIMUL[lb_hi:ub_hi], MEANS[lb_hi:ub_hi], facecolors='blue', marker='o', s=5, label='n region (hi gen./cell)')
#====================================================================
''' Plot pdf of clump distribution '''

CLUMP_MARKERS = ['s', '^', 'o']

NUM_CLUMPS_2  = CLUMP_DICT[str(int(GEN_WELL_SIMUL[len(GEN_WELL_SIMUL) - 1]))]
CLUMP_SIZES_2 = np.arange(1, len(NUM_CLUMPS_2)+1)

ax2.vlines(CLUMP_SIZES_2, 0, NUM_CLUMPS_2 / np.sum(NUM_CLUMPS_2), colors='y', lw=3, alpha=0.7, zorder=0)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NUM_CLUMPS_1  = CLUMP_DICT[str(int(GEN_WELL_SIMUL[int(len(GEN_WELL_SIMUL) / 2)]))]
CLUMP_SIZES_1 = np.arange(1, np.shape(NUM_CLUMPS_1)[0]+1)

ax2.vlines(CLUMP_SIZES_1, 0, NUM_CLUMPS_1 / np.sum(NUM_CLUMPS_1), colors='r', lw=3, alpha=0.5, zorder=1)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NUM_CLUMPS_0 = CLUMP_DICT[str(int(GEN_WELL_SIMUL[0]))]
CLUMP_SIZES_0 = np.arange(1, np.shape(NUM_CLUMPS_0)[0]+1)

ax2.vlines(CLUMP_SIZES_0, 0, NUM_CLUMPS_0 / np.sum(NUM_CLUMPS_0), colors='b', lw=3, alpha=0.5, zorder=2)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ax2.scatter(CLUMP_SIZES_2, NUM_CLUMPS_2 / np.sum(NUM_CLUMPS_2), s = 50, color='yellow', edgecolors='black', marker=CLUMP_MARKERS[2], label = 'Gen./well = ' + str(GEN_WELL_SIMUL[len(GEN_WELL_SIMUL) - 1]))
ax2.scatter(CLUMP_SIZES_1, NUM_CLUMPS_1 / np.sum(NUM_CLUMPS_1), s = 50, color='red', edgecolors ='black', marker=CLUMP_MARKERS[1], label='Gen./well = ' + str(GEN_WELL_SIMUL[int(len(GEN_WELL_SIMUL) / 2)]))
ax2.scatter(CLUMP_SIZES_0, NUM_CLUMPS_0 / np.sum(NUM_CLUMPS_0), s = 50, color='blue', edgecolors='black', marker=CLUMP_MARKERS[0], label='Gen./well = ' + str(GEN_WELL_SIMUL[0]))
#====================================================================









#====================================================================
''' Formatting '''

xMin = 10**-3
xMax = 10**3
yMin = 10**-6
yMax = 1

legendS = ''
if (num_simulations == 1):
    legendS = mlines.Line2D([], [], color='y', linestyle='--', markerfacecolor='none', markeredgecolor=color, markerfacecoloralt='none', marker=marker,
                            markersize=10, label= "Simulation: n = " + str(round(n_simul, 3)))
elif (num_simulations > 1):
    # legendS = mlines.Line2D([], [], color='k', linestyle='-', label=r'Simul. mean, $\overline{n} = $' + str(round(n_simul, 3)) + " ("+str(num_simulations)+" runs)")
    legendS = mlines.Line2D([], [], color='k', linestyle='-', markeredgecolor='green', markerfacecolor='green', marker=marker, label="Simul. mean ("+str(num_simulations)+" runs)")

legend_n = mlines.Line2D([], [],  linestyle='None', markeredgecolor='pink', markerfacecolor='pink', marker=marker, label=r'$\overline{n} = $'+str(round(n_simul, 3)))

legend_n_hi = mlines.Line2D([], [],  linestyle='None', markeredgecolor='skyblue', markerfacecolor='skyblue', marker=marker, label=r'$\overline{n} = $'+str(round(n_simul_hi, 3)))

ax0.legend(handles = [legendS, legend_n, legend_n_hi], loc='upper left', prop={'size': 8})
ax0.plot(np.linspace(xMin, xMax), np.linspace(yMin, yMax), 'k--', linewidth = 1)
ax0.set_xlabel('Genomes/cell')
ax0.set_ylabel('Infections/cell')
ax0.set_xscale('log')
ax0.set_yscale('log')
ax0.set_xlim(xMin, xMax)
ax0.set_ylim(yMin, yMax)

PlotText(ax0, PARAM_DICT_SIMUL, xMin, xMax, yMin, yMax)

ax0.text(.6 * xMin, 2.2 * yMax, 'A', fontsize=16, fontweight='bold', va='top', ha='right')

ax0.set_title(PARAM_DICT_SIMUL['simul_name'])

ax0.grid()
#--------------------------------------------------------------------
ax1.set_ylim(0, 1.25 * max(MEANS))
ax1.set_xscale('log')
ax1.set_xlabel("Genomes/cell")

title_str =''
if (PARAM_DICT_SIMUL['distribution'] == 'geometric'):
    title_str = distribution+ ' '+r'$f(clump\_size) = (1-p)^{1-clump\_size}p$'
    if (fixed == True):
        title_str += '\n'+r'$p=$'+str(mean)
    else:
        title_str += '\n'+r'$p=\frac{vMax\_c}{genomes\_well}$'
    ax1.set_ylabel('p', rotation=0)

if (PARAM_DICT_SIMUL['distribution'] == '1inf-geo'):
    title_str = distribution+ ' '+r'$f(clump\_size) = f_1 - (1-p)^{1-clump\_size}p$'
    if (fixed == True):
        title_str += '\n'+r'$p=$'+str(mean)
    else:
        title_str += '\n'+r'$p=\frac{vMax\_c}{genomes\_well}$'
    ax1.set_ylabel('p', rotation=0)

elif (PARAM_DICT_SIMUL['distribution'] == 'poisson'):
    title_str = distribution+ '\n'+r'$f(clump\_size) = \frac{e^{-\lambda} \lambda^{clump\_size}}{clump\_size!}$'
    if (fixed == True and not mean==''):
        title_str += '\n'+r'$\lambda=$'+str(mean)
    else:
        title_str += '\n'+r'$\lambda=\frac{genomes\_well}{vMax\_c}$'
    ax1.set_ylabel(r'$\lambda$', rotation=0)

ax1.set_title(title_str)

ax1.legend(loc='best')
ax1.text(.85 * np.min(GEN_CELL_SIMUL), 1.5 * max(MEANS), 'B', fontsize=16, fontweight='bold', va='top', ha='right')
#--------------------------------------------------------------------
yMin2 = 0
yMax2 = 1.1
ax2.set_ylim(yMin2, yMax2)

def get_end(arr):
    end = 0
    for i in range(np.shape(arr)[0]):
        if (int(arr[i]) > 0):
            end += 1
    return end

ax2.set_xlim(.25, np.max([get_end(NUM_CLUMPS_2), get_end(NUM_CLUMPS_1), get_end(NUM_CLUMPS_0)])+5)
ax2.xaxis.set_major_locator(MaxNLocator(integer=True))
ax2.set_xlabel('Clump size (num. virions in clump)')
ax2.set_ylabel('Relative Frequency')
ax2.legend(loc='upper right', prop = {'size' : 9})
ax2.set_yticks([0,.25,.5,.75,1])
ax2.grid()
ax2.minorticks_on()
#ax2.grid(which='minor', color='#999999', linestyle='-', alpha=0.2)
ax2.text(0, 1.05 * yMax2, 'C', fontsize=16, fontweight='bold', va='top', ha='right')
#====================================================================

''' Save figure '''
filename = MakeFilename(PARAM_DICT_SIMUL)
fig.savefig(os.path.join(os.path.join(os.getcwd(), 'figs'), filename+".pdf"), bbox_inches = 'tight', pad_inches = 0) # Save figure in the new directory


















''' Create figure '''
fig_inter, ax_inter = plt.subplots(1,1,figsize=(6, 6), dpi=100, constrained_layout=False)
#--------------------------------------------------------------------

INTERACTIONS = np.empty((num_simulations, df_simul.shape[0]), float)

for i in range(num_simulations):
    INTERACTIONS[i, :] = df_simul.loc[:, 'total_interactions run='+str(i)]

INTER_CELL = np.mean(INTERACTIONS, axis=0) / cell_count
INTER_WELL = np.mean(INTERACTIONS, axis=0)  / GEN_WELL_SIMUL

ax_inter.scatter(GEN_CELL_SIMUL, INTER_CELL, label='interactions / cell_count', marker='s')
ax_inter.scatter(GEN_CELL_SIMUL, INTER_WELL, label='interactions / Genomes/well', marker='*')
#--------------------------------------------------------------------
xMin = 10**-3
xMax = 10**3
yMin = 10**(int(np.log10(min([min(INTER_CELL), min(INTER_WELL)])))-2)
yMax = 10**(int(np.log10(max([max(INTER_CELL), max(INTER_WELL)]))) + 1)

ax_inter.vlines([GEN_CELL_SIMUL[lb_lo], GEN_CELL_SIMUL[ub_lo]], yMin, yMax, color='red', linewidth=1, linestyle='--', label='n region')
ax_inter.vlines([GEN_CELL_SIMUL[lb_hi], GEN_CELL_SIMUL[ub_hi]], yMin, yMax, color='blue', linewidth=1, linestyle='--', label='n region (hi gen./cell)')

ax_inter.set_xlabel('Genomes/cell')
ax_inter.set_ylabel('Avg. interactions/cell')
ax_inter.set_xscale('log')
ax_inter.set_yscale('log')
ax_inter.set_xlim(xMin, xMax)
ax_inter.set_ylim(yMin, yMax)
ax_inter.grid()

ax_inter.set_title('Interactions of virions with cells')

PlotText(ax_inter, PARAM_DICT_SIMUL, xMin=1, xMax=xMax, yMin=yMin, yMax=1)

ax_inter.legend(loc='best')

fig_inter.savefig(os.path.join(os.path.join(os.getcwd(), 'figs/interaction_figs'), "INTER_"+filename+".pdf"), bbox_inches = 'tight', pad_inches = 0) # Save figure in the new directory

plt.show()