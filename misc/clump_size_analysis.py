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
cutoff = 91.28
#====================================================================
''' Import data '''
df = pd.read_excel("/Users/joshuamiller/Documents/Montana Research/CLUMP_AND_INFECTION_SAME_SAMPLE/2022_10_27_TB_dose_response_size_distribution_Josh.xlsx",
                    sheet_name='Means and sd')
diameter = df.pop('d.nm')
diameter = np.asarray(diameter)
print(type(diameter))
#====================================================================
''' Initialize plot '''
num_cols = 2

fig, ax = plt.subplots(1, num_cols, figsize=(10, 5))
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
means  = ['Mean_1', 'Mean_4', 'Mean_7', 'Mean_16', 'Mean_19', 'Mean_23']
shapes = ['o', '*', 'D', '8', '^', 's']
colors = ['r', 'g', 'b', 'm', 'c', 'y']
#====================================================================
''' Average '''
means_df = df.filter(like='Mean')
stds_df  = df.filter(like='SD')

mean_of_means = means_df.mean(axis=1)
mean_of_means = np.asarray(mean_of_means)

if flatten_means:
    mean_of_means = FlattenMeans(diameter, mean_of_means, cutoff)

mean_of_means /= Trapezoid(diameter, mean_of_means) # Normalize by dividing by area  
#====================================================================
params_norm = CreateMuSigAmp(num_gauss=num_gauss, model_type='normal')
result_norm = minimize(NegLogLike, params_norm, 
                       method = 'differential_evolution', 
                       args=(diameter, mean_of_means, 'normal'),)

params_log = CreateMuSigAmp(num_gauss=num_gauss, model_type='lognormal')
result_log = minimize(NegLogLike, params_log, 
                       method = 'differential_evolution', 
                       args=(diameter, mean_of_means, 'lognormal'),)

params_skew = CreateMuSigAmp(num_gauss=num_gauss, model_type='skewnormal')
result_skew = minimize(NegLogLike, params_skew, 
                       method = 'differential_evolution', 
                       args=(diameter, mean_of_means, 'skewnormal'),)

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
y_norm = Model(diameter, result_norm.params, 'normal')
y_log  = Model(diameter, result_log.params, 'lognormal')
y_skew = Model(diameter, result_skew.params, 'skewnormal')
ssr_norm = np.sum(np.power(y_norm - mean_of_means, 2))
ssr_log  = np.sum(np.power(y_log - mean_of_means, 2))
ssr_skew = np.sum(np.power(y_skew - mean_of_means, 2))

x      = np.linspace(min(diameter), max(diameter), 1000)
y_norm = Model(x, result_norm.params, 'normal')
y_log  = Model(x, result_log.params, 'lognormal')
y_skew = Model(x, result_skew.params, 'skewnormal')
#====================================================================
''' Plot data '''
for i in range(len(means)):
    ax[0].scatter(diameter, means_df[means[i]], 
                  facecolors='None', marker=shapes[i], edgecolor=colors[i], 
                  label=means[i])

ax[1].scatter(diameter, mean_of_means, facecolors='None', marker='o', edgecolor='k', label='Mean of all samples')

x = np.linspace(min(diameter), max(diameter), 1000)
ax[1].plot(x, y_norm, 'r--', label=str(num_gauss) + r'$\mu$' + ' normal, SSR:'+str(round(ssr_norm, 8)))
ax[1].plot(x, y_log, 'b.-', label=str(num_gauss) + r'$\mu$' + ' lognormal, SSR:'+str(round(ssr_log, 8)))
ax[1].plot(x, y_skew, 'g-', label=str(num_gauss) + r'$\mu$' + ' skewnormal, SSR:'+str(round(ssr_skew, 8)))
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
for i in range(num_cols):
    ax[i].set_xscale('log')
    ax[i].legend(loc='upper left')
    ax[i].set_xlabel('Diameter (nm)')
ax[0].set_ylabel('Percent in the sample')
ax[1].set_ylabel('Percent in the sample (normalized)')

plt.show()
#====================================================================