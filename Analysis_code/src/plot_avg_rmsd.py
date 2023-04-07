import numpy as np
from numpy import array

import matplotlib.pyplot as plt
from matplotlib.pyplot import *
import matplotlib.cm as cm

import pandas as pd
#-------------------------------------------------------------------------------------#
# Fuctions
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
ranks = [1, 2]
trials = [1, 2]

DIR = '/Users/mpitman/work/protein_sim/hpc'
file = 'rmsd_peptide.xvg'

# What are you calculation the RMSD of?
group_calculated = 'Peptide'

fname = f'test_{group_calculated}_avg_rmsd.pdf'

ids = [f'trial{i}' for i in trials]
cols = [f'rank{i}' for i in ranks]

# Make df with ranks as columns and rows as trials
df = pd.DataFrame(columns = cols, index = ids)
for id in ids:
    # You may need to edit the directory location.
    data = [f'{DIR}/rank{i}/{id}/{file}' for i in ranks]
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


# control plot size
plt.figure(figsize=(20,5))


#-------------------------------------------------------------------------------------#
#     Edits not needed below here (unless you change variable names)
#-------------------------------------------------------------------------------------#

plt.gcf().subplots_adjust(bottom=0.15)

x_max = [0]
y_max = [0]

cmap = get_cmap(len(objects.keys()))
iter = 0


# Group the data (trials) by rank
for rank in ranks:
    key_list = []
    for key in objects:
        if key.startswith(f"rank{rank}"):
            key_list.append(key)

    # rank_obj is a list of dicts
    rank_obj = [objects[str(x)] for x in key_list]
    
    # Calculate the maximum time step
    timestep_max = [0]
    for i in trials:
        if np.max(len(rank_obj[i-1]['x'])) > timestep_max:
            timestep_max = np.max(len(rank_obj[i-1]['x']))
            longest_sim = i
    
    # Make empty numpy array with trials cols and ts rows
    y_vals = np.empty((timestep_max, len(trials)))
    y_vals[:] = np.nan
    
    # Fill numpy array with vals
    for i in trials:
        y_vals[0:len(rank_obj[i-1]['x']), i-1] = rank_obj[i-1]['y']
        
    # Calculate average values along rows, ignore NaNs
    y_avg = np.nanmean(y_vals, axis = 1)
    y_std = np.nanstd(y_vals, axis = 1)
    #print(f'y_avg = {y_avg}')
    #print(f'y_std = {y_std}')
    
    # Smooth avg and std
    y_avg_smooth = sliding_mean(y_avg)
    y_std_smooth = sliding_mean(y_std)
    
    
    # Generate random color samples from a cmap
    c = cmap(iter)
    
    # Plot uncertainty fill
    plt.fill_between(rank_obj[longest_sim-1]['x'], y_avg_smooth - y_std_smooth,
            y_avg_smooth + y_std_smooth, color=c, alpha=0.2)
    
    # Plot average smoothed line
    plt.plot(rank_obj[longest_sim-1]['x'], y_avg_smooth, color=c, lw=2, label = f'Rank{rank}')
    
    # Plot light raw data
    #plt.plot(objects[i]['x'], objects[i]['y'], color=c, lw=1, alpha=0.5)
    
    if np.max(y_avg_smooth + y_std_smooth) > y_max:
        y_max = np.max(y_avg_smooth + y_std_smooth)
        
    if np.max(rank_obj[longest_sim-1]['x']) > x_max:
        x_max = np.max(rank_obj[longest_sim-1]['x'])
    iter = iter + 1
    

# Make Legend.
l = legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0., fontsize = 20)
# change the line width for the legend
for line in l.get_lines():
    line.set_linewidth(4.0)

plt.tick_params(axis='both', which='major', labelsize=20)

# Set graph controls
plt.ylim([0, y_max + 0.05 * y_max])
plt.xlim([0, x_max])

plt.xlabel('Time (ns)', fontsize=20)
plt.ylabel(f'{group_calculated} RMSD (nm)',fontsize=20)
plt.tight_layout()

plt.savefig(fname)
print(f"Figure saved to file {fname}.")
plt.show()

