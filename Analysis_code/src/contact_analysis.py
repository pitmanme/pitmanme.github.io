import MDAnalysis as mda
from MDAnalysis.analysis import distances
import sys

import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.pyplot import *
from mpl_toolkits.axes_grid1 import make_axes_locatable
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
"""
  Verified package versions for this script:
  
  pandas 1.5.3
  mdanalysis  2.4.2
  matplotlib 3.7.1
  python 3.11.0
"""
#-------------------------------------------------------------------------------------#

def contacts_within_cutoff(u, group_a, group_b, radius=8.0, start = 0, end=None, step=None):
    avg_contact = np.zeros((len(group_a.positions), len(group_b.positions)))
    div = 1.0
    
    for ts in u.trajectory[start: end: step]:
        # calculate distances between group_a and group_b
        dist = distances.distance_array(group_a.positions, group_b.positions, box=u.dimensions)
        # Make binary based on radius cutoff
        dist[dist < radius] = 1.0
        dist[dist >= radius] = 0.0
        avg_contact = np.add(dist, avg_contact)
        div = div + 1.0
    avg_contact = avg_contact / div

    return avg_contact


def read_traj(structure_file, trajectory_file, groups):
    group1 = groups[0]
    group2 = groups[1]
    u = mda.Universe(structure_file, trajectory_file)
    # Chain A in the peptide, chain B is the multimer
    chA = u.select_atoms(group1)
    chB = u.select_atoms(group2)
    return u, chA, chB


def get_shape(selection1, selection2):
    return [len(selection1), len(selection2)]
    

def package_inputs(df_sims, df_strucs, ids, cols, group1, group2):
    dfs = [ df_sims, df_strucs]
    df_labels = [ids, cols]
    groups = [group1, group2]
    return dfs, df_labels, groups
    

def plot_avg_contact_maps(dfs, df_labels, groups, axs, radius=8.0, cmap='inferno', start = 0, end=None, step=None):
    # Unpack inputs
    df_sims = dfs[0]
    df_strucs = dfs[1]
    ids = df_labels[0]
    cols = df_labels[1]

    # Get shape of data for zeros array to initialize data collection
    u_test, selection1, selection2 = read_traj(
                                              df_strucs.loc[ids[0]].at[cols[0]],
                                              df_sims.loc[ids[0]].at[cols[0]],
                                              groups
                                              )
    data_shape = get_shape(selection1, selection2)

    counter = 0
    # Calculate the average contact per column, over trials
    for col_sim, col_struc in zip(df_sims, df_strucs):
        avg_arr = np.zeros(data_shape)
        iter = 0.0
        for sim, struc in zip(df_sims[col_sim], df_strucs[col_struc]):
            #print(f'sim: {sim}')
            #print(f'struc: {struc}')
            u, chA, chB = read_traj(struc, sim, groups)
            
            # Test that the assumption of atom number holds here:
            try:
                #print(data_shape)
                #print([len(chA), len(chB)])
                data_shape == [len(chA), len(chB)]
            except Exception as e:
                print(
                      '''
                      Original assumption of algorithm does not hold here.
                      Use of 'df_strucs.loc[ids[0]].at[cols[0]]' assumes that
                      the simulations have the same number of atoms for each selection
                      type. If this is not the case, you need to rewrite to vary the
                      index(here, 0) of each for loop below.
                      '''
                      )
                sys.exit()
            
            avg_contact = contacts_within_cutoff(u, chA, chB, radius = radius, start = start, end = None, step = None)
            avg_arr = np.add(avg_contact, avg_arr)

            iter = iter + 1.0
            
        # Normalize based on number of trials.
        avg_arr = avg_arr / iter
        #print(f'The max is {np.max(avg_arr)}')
        
        # Make subplots for number of ranks
        # Contact map
        iml = axs[counter].imshow(avg_arr, cmap=cmap, vmax = 1.0, vmin = 0.0)
        # Get the maximum contact value overall
        avg_max = [0]
        if np.max(avg_arr) > avg_max:
            avg_max = np.max(avg_arr)
            
        counter = counter + 1
        
    return chA, chB, axs, iml


def graph_controls(chA, chB, axs, label_frequency = 5, minor_l_alpha = 0.3, fontsize = 14):
    # Set graph controls for each subplot
    chA_names =  list(zip(chA.residues.resnames, chA.residues.resids))
    chB_names =  list(zip(chB.residues.resnames, chB.residues.resids))
    label_frequency = 5 # Plot a label every five indices

    for ax in axs.reshape(-1):
        # Major ticks
        ax.set_xticks(np.arange(0, chB.residues.resids[-1:][0]-chB.residues.resids[0], label_frequency))
        ax.set_yticks(np.arange(0, chA.residues.resids[-1:][0]-chA.residues.resids[0], label_frequency))

        # Labels for major ticks
        ax.set_xticklabels(
                          np.arange(chB.residues.resids[0], chB.residues.resids[-1:][0], label_frequency),
                          rotation=90, size = fontsize
                          )
        ax.set_yticklabels(
                          np.arange(chA.residues.resids[0], chA.residues.resids[-1:][0], label_frequency),
                          size = fontsize
                          )

        ax.set_xticks(np.arange(-.5, len(chB_names), 1), minor=True)
        ax.set_yticks(np.arange(-.5, len(chA_names), 1), minor=True)
        ax.grid(which='minor', color='w', linestyle='-', linewidth=0.4, alpha = minor_l_alpha)
        
        # Remove minor ticks
        ax.tick_params(which='minor', bottom=False, left=False)
    
    return axs

