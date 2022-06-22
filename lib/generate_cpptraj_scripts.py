#!/usr/bin/env python

import numpy as np
import sys
import yaml
import os

def parse_yaml(yaml_input):
    # Function to parse input yaml
    # It also merges the different dictionaries from each block together
    with open(yaml_input) as file:
        parameters = yaml.load(file, Loader=yaml.FullLoader)
    
    #squash yaml, information from each block together (at the moment 'Data', 'Grid' and 'Sampling')
    new_dict = {}
    for key, data in parameters.items():
        new_dict[key] = {}
        for x in data:
            for k,d in x.items():
                new_dict[key][k] = d
    return new_dict

def mkdir_if_missing(dir):
    # Auxiliar function to create directories
    if not os.path.isdir(dir): os.mkdir(dir)

def flatten_list_of_lists(input_list):
    # Auxiliar function to flatten a list of list to a single list
    return [item for sublist in input_list for item in sublist]

def initialite_traj_pool(data_parameters, check_files = False):
    # Auxiliar function to obtain the list of available trajectories. Usefull to check if they exist too
    files = {}
    for solvent in data_parameters['solvents']:
        files[solvent] = []
        for replica in range(1, data_parameters['replicas']+1):
            r = []
            for step in range(1, data_parameters['nanoseconds']+1):
                nc_file = '%s/%s_%s/md%s.nc'%(data_parameters['data directory'], solvent, replica, step)
                if check_files:
                    if not os.path.isfile(nc_file):
                        print('%s does not exist, please check the input' %nc_file)
                        exit(1)
                r.append(nc_file)
            files[solvent].append(r)
    return files

def get_grid_center(grid_parameters):
    # Auxiliar function to calculate the grid center from the grid origin, the steps to each side and the delta
    xo, yo, zo = grid_parameters['coordinates origin'].split()
    xc = float(xo)+(int(grid_parameters['dx'])*float(grid_parameters['delta'])/2)
    yc = float(yo)+(int(grid_parameters['dy'])*float(grid_parameters['delta'])/2)
    zc = float(zo)+(int(grid_parameters['dz'])*float(grid_parameters['delta'])/2)
    r = '%s %s %s' %(xc, yc, zc)
    return r

def create_sampling(sampling_parameters, trajectory_files):
    # Create the partially sampled lists of trajectories depending on the input
    r = {}
    partial_steps = sorted([int(x) for x in sampling_parameters['Sampling steps']])
    for solvent, replicas in trajectory_files.items():
            if sampling_parameters['Cross replica']: # merge the pool of trajectories from the replicas
                merged_trajectory_pool = flatten_list_of_lists(replicas)
                if partial_steps[-1] > len(merged_trajectory_pool):
                    print('Sampling steps are higher than the acumulated steps with all replicas.')
                    exit(1)
                for partial_sample in partial_steps:
                    # Is no replace good? or i'm increasing the sampling artificially
                    sample_files = [np.random.choice(merged_trajectory_pool, size=partial_sample, replace=False)
                                                                    for mr in range(sampling_parameters['Meta-replicas'])]
                    r.setdefault(solvent, dict())['CrossR_'+str(partial_sample)] = sample_files
            if sampling_parameters['Intra replica']:  # Create separate pools for each replica, NOT YET TESTED
                if partial_steps[-1] > len(replicas[0]):
                    print('Sampling steps are higher than the trajectories in each replica.')
                    exit(1)
                for partial_sample in partial_steps:
                    sample_files = [frames[:partial_sample] for frames in replicas for mr in range(sampling_parameters['Meta-replicas'])]
                    r.setdefault(solvent, dict())['IntraR_'+str(partial_sample)] = sample_files
    return r

def write_cpptraj_files(parameters, files_to_sample):
    # Main loop function to write the cpptraj input files 
    # template to format, do we want to normalize? Check!
    grid_cmd_template = 'grid {out_dxname} {dx} {delta} {dy} {delta} {dz} {delta} gridcenter {center_coords} {mask} normframe\n' 
    mkdir_if_missing(parameters['Sampling']['Output directory'])
    # cpptraj needs the grid center, so if input has only origin, we can calculate the center
    if not parameters['Grid']['coordinates center'] and parameters['Grid']['coordinates origin']:
        parameters['Grid']['coordinates center'] = get_grid_center(parameters['Grid'])
    for s_i, solvent in enumerate(parameters['Data']['solvents']):
            if not solvent in parameters['Data']:
                # Is it better to specify it on the input or have a dictionary?
                print('The probe masks for %s solvent have not been specified in the input' %solvent)
                print('You need to add the mask values for the probes in your solvent')
                exit(1)
            mkdir_if_missing('/'.join([parameters['Sampling']['Output directory'],solvent]))
            for sampling_n, sampling_f in files_to_sample[solvent].items():
                for r_i, meta_replica in enumerate(sampling_f):
                    with open('/'.join([parameters['Sampling']['Output directory'],solvent,'%s_'%(str(r_i+1))+sampling_n+ '.ptraj']), 'w') as out_file:
                        out_file.write('parm %s/%s\n'%(parameters['Data']['data directory'],parameters['Data']['topologies'][s_i]))
                        [out_file.write('trajin %s\n'%x) for x in meta_replica]
                        for probe, mask in parameters['Data'][solvent].items():
                            # A single readin of the trajectories will create all the probe grids, shoud I resample acroos probes too?
                            out_file.write(grid_cmd_template.format(out_dxname='%s/%s/%s_%s_%s_%s.dx'%(parameters['Sampling']['Output grids directory'], solvent,
                                                                                                    str(r_i+1),solvent, probe, sampling_n),
                                                                                dx=parameters['Grid']['dx'],
                                                                                dy=parameters['Grid']['dy'],
                                                                                dz=parameters['Grid']['dz'],
                                                                                delta=parameters['Grid']['delta'],
                                                                                center_coords=parameters['Grid']['coordinates center'],
                                                                                mask=mask))

if __name__ == '__main__':
    #Check usage
    if len(sys.argv) !=2:
        print('Usage: generate_cpptraj_scripts.py input.yaml')
        exit(1)
    #Read input
    yaml_input = sys.argv[1]
    parameters = parse_yaml(yaml_input) # Parse input
    trajectory_files = initialite_traj_pool(parameters['Data']) # Where are the trajectories?
    files_to_sample = create_sampling(parameters['Sampling'], trajectory_files) # Lets create the sampling 
    write_cpptraj_files(parameters, files_to_sample) # write the inputs for cpptraj

