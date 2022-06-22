#!/usr/bin/env python

# explanation of opendx here:
# https://www.ics.uci.edu/~dock/manuals/apbs/html/user-guide/x2674.html
from gridData import Grid
# to install
# conda install -c conda-forge griddataformats
# pip install gridDataFormats
import numpy as np
from math import log
import sys

def parse_hotspots_from_pdb(file):
    d = {}
    with open(file) as pdbfile:
        for line in pdbfile:
            if line.startswith('ATOM'):
                line = line.strip().split()
                atomi, x, y, z, density = line[1], line[5], line[6], line[7], line[8]
                d['Hotspot '+atomi] = {'coords':np.array([float(x), float(y), float(z)]), 'density':density}
    return d

def cluster_data_KDTree(a, thr=0.1):
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
    cmatrix = np.array([hotspot_dict[x]['coords'] for x in hotspot_dict])
    dmatrix = np.array([[np.linalg.norm(i-j) for i in cmatrix] for j in cmatrix])
    cluster = cluster_data_KDTree(cmatrix, thr=2)
    print(cluster)
    print(np.shape(cluster))

if __name__=='__main__':
    #sys.argv[1]
    d = parse_hotspots_from_pdb('dgrids/SOCS1_AF/ETA/top_density_ETA_CT.pdb')
    cluster_hotspots(d)