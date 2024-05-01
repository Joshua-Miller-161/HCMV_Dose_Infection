import numpy as np
from scipy.optimize import minimize
import scipy.stats as stats
import os
import sys
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import openpyxl
from lmfit import minimize, Parameters, Parameter, report_fit

from misc_utils import Trapezoid, CreateMuSigAmp, Model, Residuals, NegLogLike, FlattenMeans
#====================================================================
''' Initialize plot '''

fig, ax = plt.subplots(1, 2, figsize=(12, 6))
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
shapes = ['o', '*', 'D', '8', '^', 's', 'H', 'P']
colors = ['r', 'g', 'b', 'm', 'c', 'y', 'k', 'pink']
#====================================================================
''' Key parameters '''
num_gauss = 1
flatten_means = True
cutoff = 150
gen_well = 7
#====================================================================
''' Import data '''
gen_well_inf_df = pd.read_excel("data/Experimental_data_Ed_Josh.xlsx", sheet_name='2022_11_02_TB_GFP_fib')
size_df         = pd.read_excel("data/Experimental_data_Ed_Josh.xlsx", sheet_name='2022_10_27_TB_size_distribution')

DIAMETERS    = size_df.pop('d.nm')
DIAM_MEANS_DICT = {}
for key in size_df.keys():
    if ('Mean' in key):
        strs = key.split('_')
        DIAM_MEANS_DICT[strs[1]] = np.asarray(size_df.loc[:, key].values)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
keys_ = list(DIAM_MEANS_DICT.keys())
for i in range(len(keys_)):
    print(keys_[i])
    ax[0].scatter(DIAMETERS, DIAM_MEANS_DICT[keys_[i]], facecolors='None', marker=shapes[i], edgecolor=colors[i], label='GENOMES/WELL: '+str(keys_[i]))

ax[0].set_xscale('log')
ax[0].legend(loc='upper left')
ax[0].set_xlabel('Diameter (nm)')
ax[0].set_ylabel('Percent in the sample')
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cutoff = 150
for key in list(DIAM_MEANS_DICT.keys()):
    print(key)
    for i in range(np.shape(DIAM_MEANS_DICT[key])[0]):
        DIAM_MEANS_DICT[key] = FlattenMeans(DIAMETERS, DIAM_MEANS_DICT[key], cutoff)

        DIAM_MEANS_DICT[key] /= Trapezoid(DIAMETERS, DIAM_MEANS_DICT[key]) # Normalize by dividing by area 

print("ADASDASDASD", list(DIAM_MEANS_DICT.keys()))
GENOMES_WELL = list(DIAM_MEANS_DICT.keys())
#====================================================================
''' Pick which thing '''
PERC_IN_SAMPLE = DIAM_MEANS_DICT[GENOMES_WELL[gen_well]]
#====================================================================
params_norm = CreateMuSigAmp(num_gauss=num_gauss, model_type='normal')
result_norm = minimize(NegLogLike, params_norm, 
                       method = 'differential_evolution', 
                       args=(DIAMETERS, PERC_IN_SAMPLE, 'normal'),)

params_log = CreateMuSigAmp(num_gauss=num_gauss, model_type='lognormal')
result_log = minimize(NegLogLike, params_log, 
                       method = 'differential_evolution', 
                       args=(DIAMETERS, PERC_IN_SAMPLE, 'lognormal'),)

params_skew = CreateMuSigAmp(num_gauss=num_gauss, model_type='skewnormal', 
                             min_mu=60, max_mu=250, 
                             min_skew=20, max_skew=50000,
                             min_std=1, max_std=250,
                             min_amp=0.1, max_amp=2)
result_skew = minimize(NegLogLike, params_skew, 
                       method = 'differential_evolution', 
                       args=(DIAMETERS, PERC_IN_SAMPLE, 'skewnormal'),)

print(" - - - - Normal - - - - ")
report_fit(result_norm)
print("=================================")
print(" - - - - Log normal - - - - ")
report_fit(result_log)
print("=================================")
print(" - - - - Skew normal - - - - ")
report_fit(result_skew)
#====================================================================
''' Make predicted model curves '''
y_norm = Model(DIAMETERS, result_norm.params, 'normal')
y_log  = Model(DIAMETERS, result_log.params, 'lognormal')
y_skew = Model(DIAMETERS, result_skew.params, 'skewnormal')
ssr_norm = np.sum(np.power(y_norm - PERC_IN_SAMPLE, 2))
ssr_log  = np.sum(np.power(y_log - PERC_IN_SAMPLE, 2))
ssr_skew = np.sum(np.power(y_skew - PERC_IN_SAMPLE, 2))

x      = np.linspace(min(DIAMETERS), max(DIAMETERS), 1000000)
y_norm = Model(x, result_norm.params, 'normal')
y_log  = Model(x, result_log.params, 'lognormal')
y_skew = Model(x, result_skew.params, 'skewnormal')
#====================================================================

ax[1].scatter(DIAMETERS, PERC_IN_SAMPLE, facecolors='None', marker=shapes[gen_well], edgecolor=colors[gen_well], label='GENOMES/WELL: '+GENOMES_WELL[gen_well])

ax[1].plot(x, y_norm, 'r--', label=str(num_gauss) + r'$\mu$' + ' normal, SSR:'+str(round(ssr_norm, 8)))
ax[1].plot(x, y_log, 'b.-', label=str(num_gauss) + r'$\mu$' + ' lognormal, SSR:'+str(round(ssr_log, 8)))
ax[1].plot(x, y_skew, 'g-', label=str(num_gauss) + r'$\mu$' + ' skewnormal, SSR:'+str(round(ssr_skew, 8)))
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

ax[1].set_xscale('log')
ax[1].legend(loc='upper left')
ax[1].set_xlabel('Diameter (nm)')
ax[1].set_ylabel('Percent in the sample (normalized)')
plt.show()
#====================================================================