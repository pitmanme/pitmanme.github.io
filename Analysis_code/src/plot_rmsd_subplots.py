import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from numpy import array
from matplotlib.pyplot import *
import pandas as pd
#-------------------------------------------------------------------------------------#
# Functions
#-------------------------------------------------------------------------------------#
def sliding_mean(data_array, window=12):
    '''
    This function takes an array of numbers and smoothes them out
    '''
    data_array = array(data_array)  
    new_list = []  
    for i in range(len(data_array)):  
        indices = range(max(i - window + 1, 0),  
                        min(i + window + 1, len(data_array)))  
        avg = 0  
        for j in indices:  
            avg += data_array[j]  
        avg /= float(len(indices))  
        new_list.append(avg)  
          
    return array(new_list)
    
def get_cmap(n, name='jet'):
    '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct
    RGB color; the keyword argument name must be a standard mpl colormap name.'''
    return plt.cm.get_cmap(name, n)

#-------------------------------------------------------------------------------------#
# Inputs
#-------------------------------------------------------------------------------------#

#File to load
#ranks = list(range(1, 5+1))
ranks = [1, 2, 3, 4, 5]
trials = [1]

DIR = '/Users/mpitman/work/protein_sim/hpc'
file = 'rmsd_peptide.xvg'
ids = [f'trial{i}' for i in trials]
cols = [f'rank{i}' for i in ranks]

figsize = (10, len(ranks)*2)
fname = 'peptide_rmsd_subplots.pdf'

# Make df with ranks as columns and rows as trials
df = pd.DataFrame(columns = cols, index = ids)
for id in ids:
    data = [f'{DIR}/rank{i}/{id}/{file}' for i in ranks]
    print(data)
    df.loc[id] = data


# Load all values into a dictionary
objects = {}
for j in df.index:
    for i in df.columns:
        inputs = np.loadtxt(df.at[j, i], comments=['#', '@'])
        x = inputs[:, 0]
        y = inputs[:, 1]
        y_avg = sliding_mean(y)
        objects[f'{i}_{j}'] = {'x': x, 'y' : y, 'y_avg': y_avg}

#-------------------------------------------------------------------------------------#
# Driver code.
#-------------------------------------------------------------------------------------#

# Make as many subplots as ranks
fig, axs = plt.subplots(len(ranks), sharex=True, figsize=figsize)
plt.gcf().subplots_adjust(bottom=0.15)

x_max = [0]
y_max = [0]

cmap = get_cmap(len(objects.keys()))
iter = 0

for rank in ranks:
    for i in objects.keys():
        if i.startswith(f"rank{rank}"):
            #key_list.append(i)
            # Generate random color samples from a cmap
            c = cmap(iter)
            
            # Plot RMSD values
            axs[rank-1].plot(objects[i]['x'], objects[i]['y'], color=c, lw=1, alpha=0.5)
            axs[rank-1].plot(objects[i]['x'], objects[i]['y_avg'], color=c, lw=2, label = i)
            
            if np.max(objects[i]['y']) > y_max:
                y_max = np.max(objects[i]['y'])
                
            if np.max(objects[i]['x']) > x_max:
                x_max = np.max(objects[i]['x'])
            iter = iter + 1

            # Make Legend.
            l = axs[rank-1].legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0., fontsize = 16)

            # change the line width for the legend
            for line in l.get_lines():
                line.set_linewidth(4.0)



# Set graph controls
for ax in axs.reshape(-1):
    ax.set_ylim([0, y_max + 0.05 * y_max])
    ax.set_xlim([0, x_max])
    ax.tick_params(axis='both', which='major', labelsize=18)

    ax.set_ylabel('RMSD (nm)',fontsize=18)
    
plt.xlabel('Time (ns)', fontsize=18)
plt.tight_layout()

plt.savefig(fname)

#plt.show()

