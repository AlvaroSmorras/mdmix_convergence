#!/usr/bin/env python

# explanation of opendx here:
# https://www.ics.uci.edu/~dock/manuals/apbs/html/user-guide/x2674.html
from gridData import Grid
# to install
# conda install -c conda-forge griddataformats
# pip install gridDataFormats
import numpy as np
import sys
import os
import pandas as pd

def parse_hotspots_from_pdb(file):
    # Cutre-function to parse the pdb with highest densities
    # Will crash if coordinates are very big and columns merge
    d = []
    with open(file) as pdbfile:
        for line in pdbfile:
            if line.startswith('ATOM'):
                line = line.strip().split()
                atomi, x, y, z, density = line[1], line[5], line[6], line[7], line[8]
                d.append({'coords':np.array([float(x), float(y), float(z)]), 'density':density})
    # sorting them by density so the cluster representative is chosen based on that
    d.sort(key=lambda x:x['density'], reverse=True)
    return d

def cluster_to_pseudoatoms(cluster_d):
    # Auxiliary function to print the pymol command to create pseudo atoms from the clusters (just to visulize)
    for name, v in cluster_d.items():
        print('pseudoatom %s, pos=[%s, %s, %s]' %(name, v['coords'][0], v['coords'][1], v['coords'][2]))

def cluster_data_KDTree(a, thr=0.1):
    # function from internet to cluster using a kdetree
    # https://stackoverflow.com/questions/49953817/cluster-data-based-on-distance-threshold
    from scipy.spatial import cKDTree as KDTree
    t = KDTree(a)
    mask = np.ones(a.shape[:1], bool)
    idx = 0
    nxt = 1
    while nxt:
        mask[t.query_ball_point(a[idx], thr)] = False
        nxt = mask[idx:].argmax()
        mask[idx] = True
        idx += nxt
    return a[mask]

def cluster_hotspots(hotspot_dict):
    # function to cluster the hotspots based on distance and to recover them formated
    # This is necessary because cpptraj gives many points close to each other and makes for redundant assessments
    # With this we might want to increase the threshold of density to something lower than 80% of the max
    cmatrix = np.array([x['coords'] for x in hotspot_dict])
    #dmatrix = np.array([[np.linalg.norm(i-j) for i in cmatrix] for j in cmatrix])
    cluster = cluster_data_KDTree(cmatrix, thr=2)
    cluster_indexes = [np.where((cmatrix == vect).all(axis=1))[0][0] for vect in cluster]
    cluster_d = {}
    for i, ci in enumerate(cluster_indexes):
        cluster_d['Cluster '+str(i)] = {'coords':cluster[i], 'density':hotspot_dict[ci]['density']}
    return cluster_d

def find_grid_value_at_coordinates(grid, coordinates):
    # Function to obtain the value in the grid of certain coordinates
    # If grid is the path to a grid, try to load it. If not, supose its a grid
    if type(grid) == str and os.path.isfile(grid):
        grid = Grid(grid)
    # transform coordinates to index
    xi, yi, zi = ( int((coordinates[i]-grid.origin[i])/grid.delta[i]//1) for i in range(3))
    return g.grid[xi][yi][zi]
def density_in_grid_on_cluster_points(grid, cluster_d ):
    # If grid is the path to a grid, try to load it. If not, supose its a grid
    if type(grid) == str and os.path.isfile(grid):
        grid = Grid(grid)
    densities = {}
    for clust_name, clust_data in cluster_d.items():
        density = find_grid_value_at_coordinates(grid, clust_data['coords'])
        densities[clust_name] = density
    return densities

def table_initialization(grid_paths):
    # now I need this grid to be repeated n times, one for each cluster
    names = [x.split('/')[-1] for x in grid_paths]
    replicas = sorted(list(set(['Replica '+x.split('_')[0] for x in names])))
    replicas.remove('Replica full')
    nanoseconds = sorted(list(set([int(x.split('_')[-1].split('.')[0]) for x in names])))
    base_df = pd.DataFrame(index=replicas, columns=nanoseconds)
    
    for name in names:
        if name.startswith('full_sampling'):
            full, lluf, solvent, probe, nanoseconds = name.split('_')
            replica = replicas
        else:
            replica, solvent, probe, sampling, nanoseconds = name.split('_')
            replica = 'Replica '+replica
        nanoseconds = nanoseconds.split('.')[0]
        base_df.loc[replica, int(nanoseconds)] = 0
    print(base_df)
    

if __name__=='__main__':
    #sys.argv[1]
    
    ### params
    dgrids_folder = 'dgrids/SOCS1_AF/ETA/'
    ###
    
    hotspots = parse_hotspots_from_pdb(dgrids_folder+'top_density_ETA_WAT.pdb')
    cluster_d = cluster_hotspots(hotspots)
    g = Grid(dgrids_folder+'full_sampling_ETA_WAT_300.dx')
    #cluster_to_pseudoatoms(cluster_d)
    import glob
    grid_paths = glob.glob(dgrids_folder+'/*ETA_WAT*dx')
    table_initialization(grid_paths)    