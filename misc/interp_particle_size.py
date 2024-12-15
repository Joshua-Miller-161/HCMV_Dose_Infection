import sys
sys.dont_write_bytecode = True
import numpy as np
from scipy.interpolate import interp1d
import scipy.stats as stats
import os
import pandas as pd
import matplotlib.pyplot as plt
from lmfit import minimize, Parameters, Parameter, report_fit
#====================================================================
''' Set up data '''

interp_genomes_well = np.array([11258853, 8311591, 4519604, 3329074, 1802153, 1324443, 713714, 
                                523323, 280707, 205347, 109631, 80010, 42513, 30952, 22516])

size_df = pd.read_excel("data/Experimental_data_Ed_Josh.xlsx", sheet_name='2022_10_27_TB_size_distribution')

diameters = size_df.loc[:, 'd.nm'].values

particle_size_means_df = size_df.filter(like='Mean')

print(particle_size_means_df.head())

genomes_well_str = []
for key in particle_size_means_df.keys():
    if ('Mean' in key):
        genomes_well_str.append(str(key.split('_')[1]))
genomes_well = np.asarray([int(x) for x in genomes_well_str])
print(genomes_well)
print("========================================================")
#====================================================================
''' Convert to dictionary '''

count = 1
for gen_well in interp_genomes_well:
    particle_size_means_df.insert(count, 'Mean_'+str(gen_well), -9)
    count += 1
    if (count % 3 == 0):
        count += 1

print(particle_size_means_df.head(5))
#====================================================================
''' Iterate through dataframe and interpolate consecutive values '''

for row_idx in range(particle_size_means_df.shape[0]):
    print(" >> ------ >> size", size_df['d.nm'].iloc[row_idx], "<< ------ <<")

    for i in range(genomes_well.shape[0]-1):
        endpoints = genomes_well[i:i+2]

        print(" >> endpoints", endpoints, ", freq.", particle_size_means_df.loc[row_idx, 'Mean_'+genomes_well_str[i]], particle_size_means_df.loc[row_idx, 'Mean_'+genomes_well_str[i+1]])

        interp_x = interp_genomes_well[(endpoints[1] < interp_genomes_well) & (interp_genomes_well < endpoints[0])]

        f = interp1d(endpoints, 
                     [particle_size_means_df.loc[row_idx, 'Mean_'+genomes_well_str[i]], 
                      particle_size_means_df.loc[row_idx, 'Mean_'+genomes_well_str[i+1]]], 
                     kind='linear',
                     fill_value="extrapolate")

        interp_y = f(interp_x)

        for j in range(interp_y.shape[0]):
            #print(" >> j=", j, interp_x[j], interp_y[j])
            particle_size_means_df.loc[row_idx, 'Mean_'+str(interp_x[j])] = interp_y[j]

print("+++++++++++++++++++++++++++++++++++++++")
print("+++++++++++++++++++++++++++++++++++++++")
print("+++++++++++++++++++++++++++++++++++++++")
print(particle_size_means_df.head(20))

particle_size_means_df.to_excel("data/2022_10_27_size_dist_interp.xlsx")
#====================================================================
''' Plot '''

fig, ax = plt.subplots(1, 2, figsize=(14, 7))
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
shapes = ['o', '*', 'D', '8', '^', 's', 'H', 'P']

cmap = plt.get_cmap('turbo')
colors_orig   = cmap(np.linspace(0, 1, len(genomes_well_str)))
colors_interp = cmap(np.linspace(0, 1, particle_size_means_df.shape[1]))
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Original data
for i in range(len(genomes_well_str)):
    ax[0].scatter(diameters, particle_size_means_df.loc[:, 'Mean_'+genomes_well_str[i]], facecolors='None', marker=shapes[i], edgecolor=colors_orig[i], label='GENOMES/WELL: '+str(genomes_well_str[i]))
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Original + interpolated data
for i in range(len(genomes_well_str)):
    ax[1].scatter(diameters, particle_size_means_df.loc[:, 'Mean_'+genomes_well_str[i]], facecolors='None', marker=shapes[i], edgecolor=colors_orig[i])

for j in range(particle_size_means_df.shape[1]):
    ax[1].plot(diameters, particle_size_means_df.iloc[:, j], linewidth=1, color=colors_interp[j], label='GENOMES/WELL: '+str(list(particle_size_means_df.keys())[j])[5:])
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ax[0].set_xscale('log')
ax[0].legend(loc='upper left')
ax[0].set_xlabel('Diameter (nm)')
ax[0].set_ylabel('Percent in the sample')
ax[0].set_title('Original')

ax[1].set_xscale('log')
ax[1].legend(loc='upper left', prop={'size': 8})
ax[1].set_xlabel('Diameter (nm)')
ax[1].set_ylabel('Percent in the sample')
ax[1].set_title('Interpolated')

#fig.savefig('figs/size_dist_interp.pdf', bbox_inches='tight', pad_inches=0)
#====================================================================
# n = 10

# # Generate n evenly spaced values between 0 and 1
# values = np.linspace(0, 1, n)

# # Get the rainbow colormap
# cmap = plt.get_cmap('turbo')

# # Get n colors from the colormap
# colors = cmap(values)

# # Example data
# x = np.linspace(0, 10, 100)
# y = [np.sin(x + i) for i in range(n)]

# # Plotting
# fig2, ax2 = plt.subplots()
# for i in range(n):
#     ax2.plot(x, y[i], color=colors[i], label=f'Line {i+1}')

# ax2.legend()

plt.show()