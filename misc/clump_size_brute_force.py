import numpy as np
from scipy.optimize import minimize
from scipy.stats import norm, skewnorm
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import openpyxl
import yaml

sys.path.append(os.getcwd())
from misc.misc_utils import Trapezoid, FlattenMeans, FindStdev, FindStdev2, GenerateClumpDiameter, SkewNormalPDF
from simulation_utils.utils import GenClumpFromDiameter
#====================================================================
''' Get parameters '''
with open('config.yml', 'r') as c:
    config = yaml.load(c, Loader=yaml.FullLoader)

mean = config['CLUMP_PARAMETERS']['mean']
lb = config['CLUMP_PARAMETERS']['lb']
ub = config['CLUMP_PARAMETERS']['ub']
#====================================================================
''' Key parameters '''
num_gauss = 1
flatten_means = True
cutoff = 91.28

total_virions = 1000000
scheme = 'linear'
distribution = 'normal'
#====================================================================
''' Import data '''
df = pd.read_excel("data/Experimental_data_Ed_Josh.xlsx", sheet_name='2022_10_27_TB_size_distribution')
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
''' Filter to include only non-zero portions '''
diameter_nz = []
mean_of_means_nz = []

for i in range(mean_of_means.shape[0]):
    if (mean_of_means[i] > 0):
        diameter_nz.append(diameter[i])
        mean_of_means_nz.append(mean_of_means[i])
#====================================================================
''' Generate clumps using a pre-defined number of virions which most 
    closely matches the target distribution '''

max_virions_in_clump = -999
if (scheme == 'linear'):
    max_virions_in_clump = round(max(diameter_nz) / lb)
    print("diam=", max(diameter_nz), ", max=", max_virions_in_clump)

elif (scheme=='regular_polygon'):
    max_virions_in_clump = round(np.pi / np.arcsin(lb / max(diameter_nz)))
    print("diam=", max(diameter_nz), ", max=", max_virions_in_clump)
#--------------------------------------------------------------------
''' Calculate standard dev. necesary to put lb or ub at target_prob. '''
if (distribution == 'normal'):
    std_lb = FindStdev(mean, lb, 0.004)
    std_ub = FindStdev(mean, ub, 0.004)
#std2 = FindStdev2(diameter, mean_of_means, mean, lb, ub) # Also gives 70
#--------------------------------------------------------------------
def CalcNumClumps(total_virions, mean, lb, ub, max_virions_in_clump, diameter_nz, scheme='linear', dist='normal'):
    num_clumps_of_size_i = np.zeros(max_virions_in_clump) # i refers to position in list, pos. 0 means clump size 1

    virions_used = 0

    while (virions_used < total_virions):
        # 173.884216 (init = 2450)
        # std_1:   199.692402 (init = 100)
        # amp_1:   0.92133552 (init = 0.995)
        # skew_1:  3.56933583 (init = 0.995)
        diam = skewnorm.rvs(3.569, 173.884, 200, size=1)[0] # Best parameters so far
        if (diam < lb):
            diam = lb
        elif (diam > max(diameter_nz)):
            diam = diam > max(diameter_nz)
        
        num_virions, virion_diam = GenClumpFromDiameter(diam, mean, lb, ub, scheme='linear', dist='normal')
        
        num_clumps_of_size_i[num_virions-1] += 1

        virions_used += num_virions
    # print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    # print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    # print(num_clumps_of_size_i)
    # print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    # print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    return num_clumps_of_size_i
#====================================================================
num_clumps_of_size_i = CalcNumClumps(total_virions, mean, lb, ub, max_virions_in_clump, diameter_nz, scheme='linear', dist='normal')

for extra_lim in range(num_clumps_of_size_i.shape[0]-1, 0, -1):
    if not (num_clumps_of_size_i[extra_lim] <= 0):
        break

meta = [''] * (extra_lim+1)

meta[0] = scheme
meta[1] = distribution
meta[2] = 'mean='+str(mean)
meta[3] = 'lb='+str(lb)
meta[4] = 'ub='+str(ub)

print(num_clumps_of_size_i[:extra_lim+1])
df = pd.DataFrame(zip(*[np.arange(1, extra_lim+2), num_clumps_of_size_i[:extra_lim+1], meta]),
                  columns=['Clump_size', 'Number of clumps', 'Metadata'])

#df.to_csv('/Users/joshuamiller/Documents/Montana Research/CLUMP_AND_INFECTION_SAME_SAMPLE/EstimatedClumps_mu=230.csv', index=False)
#====================================================================
''' Plot data '''
for i in range(len(means)):
    ax[0].scatter(diameter, means_df[means[i]], 
                  facecolors='None', marker=shapes[i], edgecolor=colors[i], 
                  label=means[i])

ax[1].scatter(diameter, mean_of_means, facecolors='None', marker='o', edgecolor='k', label='Mean of all samples')
ax[1].scatter(diameter_nz, mean_of_means_nz, facecolors='None', marker='o', edgecolor='b')
ax[1].plot(diameter, norm.pdf(diameter, mean, std_lb) * (max(mean_of_means) / max(norm.pdf(diameter, mean, std_lb))), 'r-')
ax[1].plot(diameter, norm.pdf(diameter, mean, std_ub) * (max(mean_of_means) / max(norm.pdf(diameter, mean, std_ub))), 'b--')
ax[1].plot(diameter, norm.pdf(diameter, mean, .5 * (std_lb+std_ub)) * (max(mean_of_means) / max(norm.pdf(diameter, mean, .5 * (std_lb+std_ub)))), 'g-.')
ax[1].vlines(mean, 0, 1.5 * max(mean_of_means), label=r'$\mu$='+str(mean))
ax[1].vlines(lb, 0, 1.5 * max(mean_of_means), label='lb='+str(lb))
ax[1].vlines(ub, 0, 1.5 * max(mean_of_means), label='ub='+str(ub))
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
for i in range(num_cols):
    ax[i].set_xscale('log')
    ax[i].legend(loc='upper left')
    ax[i].set_xlabel('Diameter (nm)')
ax[0].set_ylabel('Percent in the sample')
ax[1].set_ylabel('Percent in the sample (normalized)')

plt.show()
#====================================================================