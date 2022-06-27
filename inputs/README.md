# Inputs

### 1. YAML inputs
Most of the python scripts in this project use yaml inputs. I have done this to give me flexibility on how i want to organize the data and be adaptable to people inputing their own data. As I said, this repo aims to be a tool for other people to check the convergence of their mdmix data.

#### Example of _generate_cpptraj_scripts.py_ input
This input is taken from SOCS1_AF system, and explained in further detail

```{yaml}
# Generic input to generate cpptraj scripts sampling different trajectories
Grid: 
  # If you don't know the grid center, put the origin.
  # cpptraj uses the grid center, but then the dx specify the origin.
  - coordinates center: 9.5 11.0 -8.0
  - coordinates origin: False
  - delta: 0.5 # suposed to be equal across 3 dimensions
  - dx: 156 # I separated the steps for each dimension as some boxes can be prismatic
  - dy: 156
  - dz: 156

Data:
  - replicas: 3
  - nanoseconds: 100
  - data directory: raw_data/SOCS1_AF
  - topologies: [SOCS1_AF_ETA.prmtop, SOCS1_AF_MAM.prmtop]
   # list of the topologies to pair with the solvents. Please maintain the same order.
  - solvents: [ETA, MAM] # list of the solvents to lookfor
  # specify here the probes you want to test and their mask
    ETA: 
      CT: :ETA@C1
      OH: :ETA@O1
      WAT: :WAT@O
    MAM:
      CT: :MAM@C2
      N: :MAM@N1
      O: :MAM@O1
      WAT: :WAT@O

Sampling:
  # Pool all trajectories from replicas together
  - Cross replica: True 
  # Separate the replicas. To see how the densities evolve in time.
  # Ideally from Uniform distribution to concentrated hotspots.
  - Intra replica: False 
  - Meta-replicas: 10
  - Sampling steps: [5, 10, 25, 50, 100, 150, 200, 250, 300]
  - Replacement: True 
  - Output directory: inputs/SOCS1_AF
  - Output grids directory: dgrids/SOCS1_AF
```

#### Example of _density_across_grids_pointwise.py_ input

```{yaml}
# Input for calculating density across grids in hotspots
---
input folder: dgrids/SOCS1_AF/
solvents: [ETA, MAM]
ETA: [CT, OH, WAT]
MAM: [CT, O, N, WAT]
hotspot clustering distance threshold: 2
output directory: outputs/SOCS1_AF/
output prefix: SOCS1_AF
```