#!/usr/bin/env python
import sys
import numpy as np
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

def cluster_to_pseudoatoms(cluster_d, outfile):
    if type(outfile) == str:
        fh = open(outfile, 'w')
    else:
        fh = outfile
    # Auxiliary function to print the pymol command to create pseudo atoms from the clusters (just to visulize)
    for name, v in cluster_d.items():
        fh.write(f"pseudoatom {name}, pos=[{v['coords'][0]}, {v['coords'][1]}, {v['coords'][2]}], b={v['density']}\n")
    fh.close()

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

def parse_hotspots_from_pdb(file, density_threshold = 0.8):
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
def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Take cpptraj density pdb and return clustered points for pymol visualization')
    parser.add_argument('-H', type=str, help='cpptraj pdb output from density calculation', required=True)
    parser.add_argument('-dist', type=float, help='Distance threshold to cluster the hotspots. Default: 2.0 A', nargs='?', default=2)
    parser.add_argument('-dens', type=float,help='Minumum density threshold to select hotspots. Default: 0.8', nargs='?', default=0.8)
    parser.add_argument('-out', help='Output pymol script for visualization', nargs='?', default=sys.stdout)
    args = parser.parse_args()
    return args
if __name__=='__main__':
    args = parse_args()

    hotspots_dict = parse_hotspots_from_pdb(args.H, args.dens)
    clusters = cluster_hotspots(hotspots_dict, args.dist)
    cluster_to_pseudoatoms(clusters, args.out)
