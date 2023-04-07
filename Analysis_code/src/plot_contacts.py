import matplotlib.pyplot as plt
from matplotlib.pyplot import *
from mpl_toolkits.axes_grid1 import make_axes_locatable

import pandas as pd
from contact_analysis import *
#-------------------------------------------------------------------------------------#
# Contact_Analysis: A tool to plot and analyze average contacts across simulations
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, see http://www.gnu.org/licenses/
#
# Authors: Dr Mary Pitman, 2023
#
#-------------------------------------------------------------------------------------#
"""
Contact_Analysis:

Notes on use:

  To run, edit the input variabls in 'Define Analysis Inputs' and either place your
  simulation and structure files in the directories as describes in
  'Collect Simulation and structure files' or edit 'Collect Simulation and structure
  files' to reflect where your files are stored. Unique simulations are called ranks
  here with multiple trials performed per rank. This current version supports more
  than 1 rank (simulation) and as few as 1 trial.
  
  In the some earlier versions of Matplotlib the colorbar function used here does not
  work. Set plt_cbar to False, update the code, or upgrade your Matplotlib.
  
  Verified package versions that work for this software:
  
  pandas 1.5.3
  mdanalysis  2.4.2
  matplotlib 3.7.1
  python 3.11.0

Usage:
  python plot_contacts.py
  
"""
#-------------------------------------------------------------------------------------#
#   Define Analysis Inputs:
#-------------------------------------------------------------------------------------#

#Define Simulations to analyse
ranks = [1, 2, 3, 4, 5]   # No. of types of simulations
trials = [1]  # No. of trials per simulation type.

# Base directory
DIR = '/Users/mpitman/work/protein_sim/hpc'
# Common name of simulation files
sim_file = 'protein_compressed.xtc'
# Common name of structurefiles
struc_file = 'protein_t0.pdb'

# Define which groups to calculate distances between
group1 = "chainID A and name CA"                      # y-axes
group2 = "chainID B and name CA and resid 129:243"    # x-axes

# Define titles for groups
title_group1 = "Pose"
title_group2 = "F-box only protein 44"

# Radius for contact cutoff:
radius = 8.0  # in Angstroms

# Define what timesteps to analyze
# 'step' and 'end' can also be defined
start = 500

# Plot color bar?
plt_cbar = True
# Show plot?
show_fig = True
# Save fig?
save_fig = True

# Name to save figure as:
fname = 'test_avg_contacts.pdf'

# Save of figure:
figsize=(6, len(ranks)*2.3)

# Font sizes:
titles_size = 18
labels_size = 14

#-------------------------------------------------------------------------------------#
#   Collect Simulation and structure files:
'''
   The default directory format is DIR/rank*/trial*/
   The simulation and structure files are located in the default directory.
   Edit below for your purposes!
'''
#-------------------------------------------------------------------------------------#

# Collect simulation files into dataframe
ids = [f'trial{i}' for i in trials]
cols = [f'rank{i}' for i in ranks]

df_sims = pd.DataFrame(columns = cols, index = ids)
for id in ids:
    data = [f'{DIR}/rank{i}/{id}/{sim_file}' for i in ranks]
    df_sims.loc[id] = data
    
# Collect structure files
df_strucs = pd.DataFrame(columns = cols, index = ids)
for id in ids:
    data = [f'{DIR}/rank{i}/{id}/{struc_file}' for i in ranks]
    df_strucs.loc[id] = data

# Organize inputs.
dfs, df_labels, groups = package_inputs(df_sims, df_strucs, ids, cols, group1, group2)


#-------------------------------------------------------------------------------------#
#   Driver code:
#   Edit below this point if you know what you are doing.
#-------------------------------------------------------------------------------------#

# Make as many subplots as ranks
fig, axs = plt.subplots(len(ranks), figsize=figsize)
plt.gcf().subplots_adjust(bottom=0.15)

# Calculate average contact maps
selection1, selection2, axs, iml = plot_avg_contact_maps(
                                                        dfs,
                                                        df_labels,
                                                        groups,
                                                        axs,
                                                        radius=radius,
                                                        start = start
                                                        )
# Edit ticks and labels for plots.
axs = graph_controls(selection1, selection2, axs)

# Put colorbar above subplots
if plt_cbar == True:
    divider = make_axes_locatable(axs[0])
    cax = divider.append_axes("top", size="7%", pad="8%")
    fig.add_axes(cax)
    clb = plt.colorbar(
                      iml, cax=cax, orientation='horizontal',
                      shrink=0.5, location='top', pad = 0.08
                      )
    clb.ax.tick_params(labelsize = labels_size)
    clb.ax.set_title(f'P(Cα contact < {radius} $\AA$)', fontsize = labels_size)

#Control plot titles
iter = 0
for ax in axs.reshape(-1):
    ax.set_ylabel(f'{title_group1} {ranks[iter]}', fontsize = titles_size)
    iter = iter + 1
    
axs[len(ranks)-1].set_xlabel(title_group2, fontsize = titles_size)

plt.xticks(rotation = 90)

if save_fig == True:
    plt.savefig(fname)
    
if show_fig == True:
    plt.show()
