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
    with open(outfile, 'w') as fh:
        # Auxiliary function to print the pymol command to create pseudo atoms from the clusters (just to visulize)
        for name, v in cluster_d.items():
            fh.write(f"pseudoatom {name}, pos=[{v['coords'][0]}, {v['coords'][1]}, {v['coords'][2]}]\n")

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

if __name__=='__main__':
    if len(sys.argv) != 4:
        print(f'usage: {sys.argv[0]} hotspots.pdb distance_threshold pymol_script.pml')
        exit(1)
    hotspots_pdb = sys.argv[1]
    distance_threshold = sys.argv[2]
    output_pml = sys.argv[3]

    hotspots_dict = parse_hotspots_from_pdb(hotspots_pdb)
    clusters = cluster_hotspots(hotspots_dict, distance_threshold)
    cluster_to_pseudoatoms(clusters, output_pml)
