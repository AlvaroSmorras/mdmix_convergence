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
import glob


def mkdir_if_missing(dir):
    # Auxiliar function to create directories
    if not os.path.isdir(dir): os.mkdir(dir)

def parse_hotspots_from_pdb(file, density_threshold):
    # Cutre-function to parse the pdb with highest densities
    # Will crash if coordinates are very big and columns merge
    d = []
    with open(file) as pdbfile:
        for line in pdbfile:
            if line.startswith('ATOM'):
                line = line.strip().split()
                atomi, x, y, z, density = line[1], line[5], line[6], line[7], line[8]
                if float(density) >= density_threshold: # filter-out densities lower than threshold
                    d.append({'coords':np.array([float(x), float(y), float(z)]), 'density':density})
    # sorting them by density so the cluster representative is chosen based on that
    d.sort(key=lambda x:x['density'], reverse=True)
    return d

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

def cluster_hotspots(hotspot_dict, distance_threshold):
    # function to cluster the hotspots based on distance and to recover them formated
    # This is necessary because cpptraj gives many points close to each other and makes for redundant assessments
    # With this we might want to increase the threshold of density to something lower than 80% of the max
    cmatrix = np.array([x['coords'] for x in hotspot_dict])
    #dmatrix = np.array([[np.linalg.norm(i-j) for i in cmatrix] for j in cmatrix])
    cluster = cluster_data_KDTree(cmatrix, thr=distance_threshold)
    cluster_indexes = [np.where((cmatrix == vect).all(axis=1))[0][0] for vect in cluster]
    cluster_d = {}
    for i, ci in enumerate(cluster_indexes):
        cluster_d['Cluster_'+str(i)] = {'coords':cluster[i], 'density':hotspot_dict[ci]['density']}
    return cluster_d

def find_grid_value_at_coordinates(grid, coordinates):
    # Function to obtain the value in the grid of certain coordinates
    # If grid is the path to a grid, try to load it. If not, supose its a grid
    if type(grid) == str and os.path.isfile(grid):
        grid = Grid(grid)
    # transform coordinates to index
    xi, yi, zi = ( int((coordinates[i]-grid.origin[i])/grid.delta[i]//1) for i in range(3))
    return grid.grid[xi][yi][zi]
def density_in_grid_on_cluster_points(grid, cluster_d ):
    # If grid is the path to a grid, try to load it. If not, supose its a grid
    if type(grid) == str and os.path.isfile(grid):
        grid = Grid(grid)
    densities = {}
    for clust_name, clust_data in cluster_d.items():
        density = find_grid_value_at_coordinates(grid, clust_data['coords'])
        densities[clust_name] = density
    return densities

def iterate_grids_and_clusters(grid_paths, cluster_d):
    # Main function to look for grid values corresponding to clusters

    # Parsing variables and initializing base dataframe
    names = [x.split('/')[-1] for x in grid_paths]
    dgrid_dir = '/'.join(grid_paths[0].split('/')[:-1])
    replicas = sorted(list(set(['Replica_'+x.split('_')[0] for x in names])))
    replicas.remove('Replica_full-sampling')
    nanoseconds = sorted(list(set([str(x.split('_')[-1].split('.')[0]) for x in names])))+['full-sampling']
    base_df = pd.DataFrame(columns=nanoseconds)

    # External loop will be grids, as they take the most time to read
    for name in names:
        # full-sampling is a special case as it has no replicas and will fill a whole column
        grid = Grid('/'.join([dgrid_dir, name]))
        if name.startswith('full-sampling'): 
            replica, solvent, probe, nanoseconds = name.split('_')
            nanoseconds = 'full-sampling'
        else:
            replica, solvent, probe, sampling, nanoseconds = name.split('_')
        nanoseconds = nanoseconds.split('.')[0]
        cluster_densities = density_in_grid_on_cluster_points(grid, cluster_d)
        for cluster_n, density in cluster_densities.items():
            # again full-sampling being special
            if name.startswith('full-sampling'): 
                row_name = ['-'.join([cluster_n,sampling,replica]) for replica in replicas]
            else:
                row_name = '-'.join([cluster_n,sampling,'Replica_%s'%replica])
            # fill dataframe with data            
            base_df.loc[row_name, nanoseconds] = density
    # Comodity to have things sorted by cluster
    base_df.sort_index(inplace=True)
    return base_df

def parse_yaml(yaml_input):
    import yaml
    # Function to parse input yaml
    # It also merges the different dictionaries from each block together
    with open(yaml_input) as file:
        parameters = yaml.load(file, Loader=yaml.FullLoader)
    
    return parameters

def iterate_solvents_and_probes(parameters):
    mkdir_if_missing(parameters['output directory'])
    for solvent in parameters['solvents']:
        grids_folder = parameters['input folder'] + solvent +'/'
        for probe in parameters[solvent]:
            print('Working on %s_%s'%(solvent, probe))
            hotspots_file = grids_folder+f'top_density_{solvent}_{probe}.pdb'
            hotspots = parse_hotspots_from_pdb(hotspots_file, parameters['minimum hotspot density'])
            cluster_d = cluster_hotspots(hotspots, distance_threshold=parameters['hotspot clustering distance threshold'])
            #print(cluster_d)
            grid_paths = sorted(glob.glob(grids_folder+f'/*{solvent}_{probe}*dx'))
            densities_dataframe = iterate_grids_and_clusters(grid_paths, cluster_d)
            densities_dataframe.to_csv(parameters['output directory']+parameters['output prefix']+f'_{solvent}_{probe}.csv')


            

if __name__=='__main__':
        #Check usage
    if len(sys.argv) !=2:
        print(f'Usage: {sys.argv[0]} input.yaml')
        exit(1)
    #Read input
    input_yaml_file = sys.argv[1]
    parameters = parse_yaml(input_yaml_file)
    iterate_solvents_and_probes(parameters)
