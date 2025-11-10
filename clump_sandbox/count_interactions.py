import sys
sys.dont_write_bytecode = True
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import geom, poisson
import numpy as np
from collections import Counter
import random
import itertools
#====================================================================
''' Parameters '''

num_simulations = 100
cell_count   = 100 # number of cells
genomes_well = np.arange(10, 110, 10) # number of virions that could attack the cells, 10, 20, 30, ..., 100

p = .5 # Geometric dist. parameter
#lam = 1 # Poisson dist. parameter

num_visits_per_cell_param = 2 # The parameter that controls how many of virions and/or clumps that will visit a cell
#====================================================================
''' Simulate interactions '''

interactions       = np.empty((num_simulations, genomes_well.shape[0]), np.float32) # Array to store number of interactions
visitor_sizes_mean = np.empty((num_simulations, genomes_well.shape[0]), np.float32) # Array to store the average clump size of the visiting virions/clumps
visitor_sizes_std  = np.empty((num_simulations, genomes_well.shape[0]), np.float32)

for simul in range(num_simulations):
    print("##############################################################################")
    print("##############################################################################")
    print("##############################################################################")
    print("                            simulation", simul)
    print("##############################################################################")
    print("##############################################################################")
    print("##############################################################################")
    for i in range(genomes_well.shape[0]): # Go through each genomes/well value
        
        print("_______________________________________________________________________________")
        print(" << >> << >> << >> genomes/well", genomes_well[i], ", p=", p, ", cell_count=", cell_count, '<< >> << >> << >>')
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # We must now sample single virions and/or clumps from the distribution
        # but, ENSURE THE TOTAL NUMBER OF VIRIONS IN ALL CLUMPS = genomes/well
        # to do this, take a clump size one at a time from the distribution until
        # the sum of the clump sizes equals genomes/well

        counter = 0
        sampled_clump_sizes = [] # empty list to hold clump sizes. The sum of this list must be genomes/well
        while True:
            sampled_clump_size = geom.rvs(p, size=1)[0] # Get 1 clump size

            counter += sampled_clump_size
            
            if (counter == genomes_well[i]): # Exit loop when correct sum is reached
                sampled_clump_sizes.append(sampled_clump_size)
                break

            elif (counter < genomes_well[i]): # Check if cumulative sum of clump sizes is bigger than/genomes well
                sampled_clump_sizes.append(sampled_clump_size) # If not, then we can safely put the clump in
            
            elif (counter > genomes_well[i]):
                counter -= sampled_clump_size # Backtrack if genomes/well was overshot
            
        print(" >> sampled_clump_sizes:", sampled_clump_sizes)
        print(" _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _")

        # Count the number times each clump size occurred
        count_dict = Counter(list(sorted(sampled_clump_sizes))) # Sort clump sizes and count them
        result     = {str(k): v for k, v in count_dict.items()} # Organize sampled clump sizes into a dictionary {'1': frequency, '2': frequency, ...}
        prod       = sum([int(key) * int(value) for key, value in count_dict.items()]) # Multiply frequency times clump size

        print(" >> Frequencies:", result)
        print(" >> sum(freq. * clump size) =", prod)
        print(" _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _")
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Iterate through cells and count how many total interactions they get

        total_interactions = 0
        all_visitors = [] # List to hold all visitors that interacted with the cells

        for cell in range(cell_count):
            
            # Use a poisson distribution to see how many virions and/or clumps will visit the cell
            num_visitors = poisson.rvs(num_visits_per_cell_param, size=1)[0] 
            
            if (num_visitors > 0 and (num_visitors <= len(sampled_clump_sizes))): # second condition b/c
                clump_sizes_of_visitors = list(random.sample(sampled_clump_sizes, num_visitors)) # Randomly pick num_visitor visitors (individual virions of clumps)

                all_visitors.append(clump_sizes_of_visitors)

                #print(" >> >> cell:", cell, ", num_visitors:", num_visitors, ", clump size of visitors:", clump_sizes_of_visitors)
                #print(" -> -> -> -> -> -> interactions=sum(clump sizes of visitors) =", sum(clump_sizes_of_visitors), "<- <- <- <- <- <-")

                total_interactions += sum(clump_sizes_of_visitors)
        
        interactions[simul, i] = total_interactions

        flat_visitors                = np.asarray(list(itertools.chain(*all_visitors))) # Flatten to a single list
        visitor_sizes_mean[simul, i] = np.mean(flat_visitors) # Compute mean
        visitor_sizes_std[simul, i]  = np.std(flat_visitors)  # Compute standard deviation

        print(" >> total interactions:", total_interactions, ", avg. visitor size:", visitor_sizes_mean[simul, i], ", stdev:", visitor_sizes_std[simul, i])












#====================================================================
''' Plot '''

fig = plt.figure(figsize=(12, 7), dpi=100, layout='tight')
gs  = fig.add_gridspec(nrows=3, ncols=4)
ax0 = fig.add_subplot(gs[0, 0:2])
ax1 = fig.add_subplot(gs[1:, 0:2])
ax2 = fig.add_subplot(gs[:2, 2:])
ax3 = fig.add_subplot(gs[2, 2])
ax4 = fig.add_subplot(gs[2, 3])
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot distribution

clump_sizes = np.arange(1,11) # 1, 2, 3, ..., 10
clump_freqs = geom.pmf(clump_sizes, p) # PDF

ax0.vlines(clump_sizes, 0, clump_freqs, colors='b', lw=3, alpha=0.5, zorder=2)
ax0.scatter(clump_sizes, clump_freqs, s=50, color='b', edgecolors='black', marker='o', label = 'p='+str(p))
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot Poisson distribution

num_visitors  = np.arange(0, 11) # 0, 1, 2, 3, ..., 10
visitor_freqs = poisson.pmf(num_visitors, num_visits_per_cell_param) # PDF

ax1.vlines(num_visitors, 0, visitor_freqs, colors='r', lw=3, alpha=0.5, zorder=2)
ax1.scatter(num_visitors, visitor_freqs, s=50, color='r', edgecolors='black', marker='*', label = r'$\lambda=$'+str(num_visits_per_cell_param))
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot interactions

mean_interactions_cell = np.mean(interactions, axis=0) / cell_count
cis = 1.96 * np.std(interactions/cell_count, axis=0) / np.sqrt(num_simulations)

ax2.plot(genomes_well/cell_count, mean_interactions_cell, color='black', linestyle='--', marker='s')
ax2.fill_between(genomes_well/cell_count, (mean_interactions_cell-cis), (mean_interactions_cell+cis), color='k', alpha=.1)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot mean clump size of visitors

mean_visitor_sizes_mean = np.mean(visitor_sizes_mean, axis=0)
cis = 1.96 * np.std(visitor_sizes_mean/cell_count, axis=0) / np.sqrt(num_simulations)

ax3.plot(genomes_well/cell_count, mean_visitor_sizes_mean, color='green', linestyle='--', marker='^')
ax3.fill_between(genomes_well/cell_count, (mean_visitor_sizes_mean - cis), (mean_visitor_sizes_mean+cis), color='green', alpha=.1)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

mean_visitor_sizes_std = np.mean(visitor_sizes_std, axis=0)
cis = 1.96 * np.std(visitor_sizes_std, axis=0) / np.sqrt(num_simulations)

ax4.plot(genomes_well/cell_count, mean_visitor_sizes_std, color='magenta', linestyle='--', marker='v')
ax4.fill_between(genomes_well/cell_count, (mean_visitor_sizes_std - cis), (mean_visitor_sizes_std+cis), color='magenta', alpha=.1)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Format plots

ax0.set_title('geometric')
ax0.legend(loc='best')
ax0.set_xlabel('clump size (virions in clump)')
ax0.set_ylabel('f(clump size)')
ax0.grid()


ax1.set_title('Poisson')
ax1.legend(loc='best')
ax1.set_xlabel('num. visitors')
ax1.set_ylabel('f(number of visitors/cell)')
ax1.grid()


ax2.set_title('Interactions/cell, cell_count='+str(cell_count))
ax2.set_xlabel('genomes/cell')
ax2.set_ylabel('interactions/cell')
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlim(10**-2, 10**2)


ax3.set_xlabel('genomes/cell')
ax3.set_ylabel('mean clump size')
ax3.set_xscale('log')
ax3.set_yscale('log')
ax3.set_xlim(10**-2, 10**2)


ax4.set_xlabel('genomes/cell')
ax4.set_ylabel('stdev. clump size')
ax4.set_xscale('log')
ax4.set_yscale('log')
ax4.set_xlim(10**-2, 10**2)

plt.show()




# 1-inflated geometric distribution
# y = np.asarray(geom.pmf(x, p[i]))

# y[0] = .7
# y_rest = y[1:]
# y_rest_scaled = y_rest * ((1-y[0]) / np.sum(y_rest))

# geo_1inf[i, :] = np.concatenate(([.7], y_rest_scaled))

# #some example data
# x = np.linspace(0.1, 9.9, 20)
# y = 3.0 * x
# #some confidence interval
# ci = 1.96 * np.std(y)/np.sqrt(len(x))

# fig, ax = plt.subplots()
# ax.plot(x,y)
# ax.fill_between(x, (y-ci), (y+ci), color='b', alpha=.1)