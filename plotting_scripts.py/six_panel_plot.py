import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import pandas as pd
from scipy.stats import norm
from lmfit import Parameters, minimize, report_fit
import yaml
import sys
import os

sys.path.append(os.getcwd())

from misc.misc_utils import ExtractParams
from simulation_utils.utils import PrepareData
from plotting_utils.utils import BasicFormat, PlotSimul, PlotFit, PlotText
#====================================================================
''' Make subplot array '''
fig = plt.figure(figsize=(8, 8), dpi=100, constrained_layout=False)
gs  = fig.add_gridspec(6, 2)
ax_data       = fig.add_subplot(gs[:2, 0])
ax_res_inf    = fig.add_subplot(gs[0, 1])
ax_clump_dist = fig.add_subplot(gs[1, 1])
ax_null       = fig.add_subplot(gs[2:4, 0])
ax_clump      = fig.add_subplot(gs[2:4, 1])
ax_comp       = fig.add_subplot(gs[4:6, 0])
ax_acc_dam    = fig.add_subplot(gs[4:6, 1])

fig.subplots_adjust(wspace=.4, hspace=.8)
#====================================================================
''' Files '''
SHEET_NAMES = ['2021_10_05 TB_GFP_epithelial', '2020_07_02 ME_GFP_fibroblast', 
               '2020_05_29 TR_GFP_fibroblast', '2021_07_13 GFP_TB_fibroblast', 
               '2020_08_12 TB_GFP_fibroblast', '2020_09_14 TR_GFP_epithelial',
               '2021_08_13 ME_mC_epithelial', '2022_11_02_TB_GFP_fib']

data_file    = "data/Experimental_data_Ed_Josh.xlsx"
null_file    = "simulation_results/null/NullSimul_2021_07_13 GFP_TB_fibroblast_s=50_vMax=1600.0_b=0.1_r=1_n=10.csv"
clump_file   = ""
comp_file    = "simulation_results/comp/CompSimul_2021_07_13 GFP_TB_fibroblast_s=50_vMax=1600.0_k=-0.03_r=1_n=10.csv"
acc_dam_file = "simulation_results/acc_dam/AccDamSimul_2021_07_13 GFP_TB_fibroblast_s=100_vMax=700.0_b=-0.8_r=1_n=10.csv"
#====================================================================
with open('config.yml', 'r') as c:
    config = yaml.load(c, Loader=yaml.FullLoader)
sheet  = config['SIMULATION_PARAMETERS']['sheet']
scale  = config['SIMULATION_PARAMETERS']['scale']

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
''' Get corresponding experimental data '''

SHEET_NAMES = ['2021_10_05 TB_GFP_epithelial', '2020_07_02 ME_GFP_fibroblast', 
               '2020_05_29 TR_GFP_fibroblast', '2021_07_13 GFP_TB_fibroblast', 
               '2020_08_12 TB_GFP_fibroblast', '2020_09_14 TR_GFP_epithelial',
               '2021_08_13 ME_mC_epithelial', '2022_11_02_TB_GFP_fib']
df_data = pd.read_excel('data/Experimental_data_Ed_Josh.xlsx', sheet_name=SHEET_NAMES[sheet])

#--------------------------------------------------------------------
PARAM_DICT_DATA = ExtractParams(df_data)
cell_count = int(PARAM_DICT_DATA['cell_count'] / scale)
#--------------------------------------------------------------------

GEN_WELL_DATA, GEN_CELL_DATA, INF_CELL_DATA, num_zeros = PrepareData(df_data, scale)
#====================================================================
''' Plot data '''

ax_data.scatter(GEN_CELL_DATA, INF_CELL_DATA, facecolors='none', edgecolors=color, marker=marker, alpha=1)

print("AHASDHALJDHAS:JDA:SLDKJA:SLKJDNA:LSKDA:", np.shape(GEN_CELL_DATA), np.shape(INF_CELL_DATA))

g_data, n_data = PlotFit(ax_data, GEN_CELL_DATA, INF_CELL_DATA, LOWERS[sheet], UPPERS[sheet])

BasicFormat(ax_data)

ax_data.text(10**-2, 2*10**-6, 'n='+str(round(n_data, 5))+'\n'+SHEET_NAMES[sheet])

#====================================================================
''' Plot resistivity and infectivity histograms '''
x = np.linspace(-15, 15, 100000)
ax_res_inf.plot(x, norm.pdf(x, mu_res, std_res), 'k-', label='Cell')
ax_res_inf.plot(x, norm.pdf(x, mu_inf, std_inf), color=color, linestyle='-', label='Virus')

ax_res_inf.set_ylabel('Freq.')
ax_res_inf.legend()
#====================================================================
''' Plot changing clump diameter '''



#====================================================================
''' Plot null '''

GEN_CELL_SIMUL, INF_CELL_SIMUL = PlotSimul(ax_null, null_file, band_type, replacement_val)

g_null, n_null = PlotFit(ax_null, GEN_CELL_SIMUL, INF_CELL_SIMUL, LOWERS[sheet], UPPERS[sheet])

BasicFormat(ax_null)

ax_null.text(10**-2, 2*10**-6, 'n='+str(round(n_null, 5))+'\nnull model')
#====================================================================




#====================================================================
''' Plot compensation '''

GEN_CELL_SIMUL, INF_CELL_SIMUL = PlotSimul(ax_comp, comp_file, band_type, replacement_val)

g_comp, n_comp = PlotFit(ax_comp, GEN_CELL_SIMUL, INF_CELL_SIMUL, LOWERS[sheet], UPPERS[sheet])

BasicFormat(ax_comp)

ax_comp.text(10**-2, 2*10**-6, 'n='+str(round(n_comp, 5))+'\nCompensation model')


#====================================================================


#====================================================================

''' Plot compensation '''

GEN_CELL_SIMUL, INF_CELL_SIMUL = PlotSimul(ax_acc_dam, acc_dam_file, band_type, replacement_val)

g_acc_dam, n_acc_dam = PlotFit(ax_acc_dam, GEN_CELL_SIMUL, INF_CELL_SIMUL, LOWERS[sheet], UPPERS[sheet])

BasicFormat(ax_acc_dam)

ax_acc_dam.text(10**-2, 2*10**-6, 'n='+str(round(n_acc_dam, 5))+'\nAcc. Damage model')


plt.show()