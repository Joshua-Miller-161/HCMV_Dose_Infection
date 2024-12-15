import sys
sys.dont_write_bytecode = True
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
''' Key parameters '''

num_gauss = 1
flatten_means = True
cutoff = 100
gen_well = 20
#====================================================================
''' Import data '''

gen_well_inf_df = pd.read_excel("data/Experimental_data_Ed_Josh.xlsx", sheet_name='2022_11_02_TB_GFP_fib')
size_df         = pd.read_excel("data/Experimental_data_Ed_Josh.xlsx", sheet_name='2022_10_27_TB_size_dist_interp')

orig_gen_well   = [15240000, 6131312, 2450310, 972623, 383425, 150101, 58345, 16367]

diameters = size_df.pop('d.nm')
DIAM_MEANS_DICT = {}
for key in size_df.keys():
    if ('Mean' in key):
        strs = key.split('_')
        DIAM_MEANS_DICT[strs[1]] = np.asarray(size_df.loc[:, key].values)

interp_gen_well = np.asarray([int(x) for x in list(DIAM_MEANS_DICT.keys())])


#====================================================================
''' Normalize '''

for key in list(DIAM_MEANS_DICT.keys()):
    print(key)
    for i in range(np.shape(DIAM_MEANS_DICT[key])[0]):
        DIAM_MEANS_DICT[key] = FlattenMeans(diameters, DIAM_MEANS_DICT[key], cutoff)

        DIAM_MEANS_DICT[key] /= Trapezoid(diameters, DIAM_MEANS_DICT[key]) # Normalize by dividing by area 

print(interp_gen_well)
#====================================================================
''' Pick which thing gen/well to analyze '''

perc_in_sample = DIAM_MEANS_DICT[str(interp_gen_well[gen_well])]
#====================================================================
''' Fit curves to data '''

params_norm = CreateMuSigAmp(num_gauss=num_gauss, model_type='normal')
result_norm = minimize(NegLogLike, params_norm, 
                       method = 'differential_evolution', 
                       args=(diameters, perc_in_sample, 'normal'),)

params_log = CreateMuSigAmp(num_gauss=num_gauss, model_type='lognormal')
result_log = minimize(NegLogLike, params_log, 
                       method = 'differential_evolution', 
                       args=(diameters, perc_in_sample, 'lognormal'),)

params_skew = CreateMuSigAmp(num_gauss=num_gauss, model_type='skewnormal', 
                             min_mu=30, max_mu=350, 
                             min_skew=1, max_skew=200000,
                             min_std=1, max_std=140,
                             min_amp=0.1, max_amp=2)
result_skew = minimize(NegLogLike, params_skew, 
                       method = 'differential_evolution', 
                       args=(diameters, perc_in_sample, 'skewnormal'),)

print(" - - - - Normal", interp_gen_well[gen_well], "- - - - ")
report_fit(result_norm)
print("=================================")
print(" - - - - Log normal", interp_gen_well[gen_well], "- - - - ")
report_fit(result_log)
print("=================================")
print(" - - - - Skew normal", interp_gen_well[gen_well], " - - - - ")
report_fit(result_skew)
#====================================================================
''' Make predicted model curves '''

y_norm = Model(diameters, result_norm.params, 'normal')
y_log  = Model(diameters, result_log.params, 'lognormal')
y_skew = Model(diameters, result_skew.params, 'skewnormal')
ssr_norm = np.sum(np.power(y_norm - perc_in_sample, 2))
ssr_log  = np.sum(np.power(y_log - perc_in_sample, 2))
ssr_skew = np.sum(np.power(y_skew - perc_in_sample, 2))

x      = np.linspace(min(diameters), max(diameters), 1000000)
y_norm = Model(x, result_norm.params, 'normal')
y_log  = Model(x, result_log.params, 'lognormal')
y_skew = Model(x, result_skew.params, 'skewnormal')
#====================================================================
''' Plot '''

fig, ax = plt.subplots(1, 2, figsize=(12, 6))

shapes = ['1', '*', '2', '8', '3', 's', '4', '+']

cmap = plt.get_cmap('turbo')
colors_interp = cmap(np.linspace(0, 1, interp_gen_well.shape[0]))
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

_ = 0
for i in range(interp_gen_well.shape[0]):
    if (interp_gen_well[i] in orig_gen_well):
        if (i == gen_well):
            ax[0].scatter(diameters, DIAM_MEANS_DICT[str(interp_gen_well[i])], facecolors=colors_interp[i], marker='D', edgecolor='black')
        else:
            ax[0].scatter(diameters, DIAM_MEANS_DICT[str(interp_gen_well[i])], facecolors=colors_interp[i], marker=shapes[_], edgecolor=colors_interp[i])
        _ += 1
    
    if (i == gen_well):
        ax[0].plot(diameters, DIAM_MEANS_DICT[str(interp_gen_well[i])], linewidth=1, color=colors_interp[i], label='GENOMES/WELL: '+r'$\mathbf{'+str(interp_gen_well[i])+'}$')
    else:
        ax[0].plot(diameters, DIAM_MEANS_DICT[str(interp_gen_well[i])], linewidth=1, color=colors_interp[i], label='GENOMES/WELL: '+str(interp_gen_well[i]))
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

ax[1].scatter(diameters, perc_in_sample, facecolors=colors_interp[gen_well], marker='D', edgecolor='black', label='GENOMES/WELL: '+r'$\mathbf{'+str(interp_gen_well[gen_well])+'}$')

ax[1].plot(x, y_norm, 'r--', label=str(num_gauss) + r'$\mu$' + ' normal, SSR:'+str(round(ssr_norm, 8)))
ax[1].plot(x, y_log, 'b.-', label=str(num_gauss) + r'$\mu$' + ' lognormal, SSR:'+str(round(ssr_log, 8)))
ax[1].plot(x, y_skew, 'g-', label=str(num_gauss) + r'$\mu$' + ' skewnormal, SSR:'+str(round(ssr_skew, 8)))
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

ax[0].set_xscale('log')
ax[0].legend(loc='upper left', prop={'size': 8})
ax[0].set_xlabel('Diameter (nm)')
ax[0].set_ylabel('Percent in the sample')

ax[1].set_xscale('log')
ax[1].legend(loc='upper left')
ax[1].set_xlabel('Diameter (nm)')
ax[1].set_ylabel('Percent in the sample (normalized)')
plt.show()
#====================================================================