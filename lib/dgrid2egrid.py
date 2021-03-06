#!/usr/bin/env python

# explanation of opendx here:
# https://www.ics.uci.edu/~dock/manuals/apbs/html/user-guide/x2674.html
from gridData import Grid
import numpy as np
from math import log
import sys

def info_from_dx(dx_files, outfile = False):
    g = load_and_sum_dgrid_files(dx_files)
    if outfile: outfh = open(outfile, 'w')
    xyz = []
    values = []
    non_solvent_accesible = 0
    for xi in range(g.grid.shape[0]):
        for yi in range(g.grid.shape[1]):
            for zi in range(g.grid.shape[2]):

                x = (g.origin[0]+g.delta[0]*xi)
                y = (g.origin[1]+g.delta[1]*yi)
                z = (g.origin[2]+g.delta[2]*zi)
                value = g.grid[xi][yi][zi]

                if outfile: outfh.write('%s\t%s\t%s\t%s\n' %(value,x,y,z))
                xyz.append((x,y,z))
                values.append(value)
                if value == 0: non_solvent_accesible += 1
    solvent_accessible = (g.grid.shape[0]*g.grid.shape[1]*g.grid.shape[2]) - non_solvent_accesible
    if outfile: outfh.close()
    return xyz, values, solvent_accessible, g

def dgrid2egrid(g, out_grid, expected_counts):
    T = 300
    kb = 0.0019872041
    for xi in range(g.grid.shape[0]):
            for yi in range(g.grid.shape[1]):
                for zi in range(g.grid.shape[2]):
                    count = g.grid[xi][yi][zi]
                    if count == 0:
                        fenergy = '999'
                    else:
                        fenergy = round(-kb*T*log(count/expected_counts),6)
                    g.grid[xi][yi][zi] = fenergy
    g.export(out_grid, type='float')

def load_and_sum_dgrid_files(dx_files):
    grids = []
    for dx_file in dx_files:
        grids.append(Grid(dx_file))
    return sum(grids)


if __name__ == '__main__':
    dgrid_files = [sys.argv[x] for x in range(1, len(sys.argv)-1)]
    egrid_file = sys.argv[-1]
    xyz_list, counts_list, solvent_accesible_gridpoints, sum_grid = info_from_dx(dgrid_files)
    expected_counts = sum(counts_list)/solvent_accesible_gridpoints
    dgrid2egrid(sum_grid, expected_counts = expected_counts, out_grid = egrid_file)
