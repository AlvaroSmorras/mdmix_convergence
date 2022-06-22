# MDMix convergence anlysis:
## How many nanoseconds should I run?
The increase in computational power has made us launch longer and longer MDMix simulations, but is it necessary?
Maybe we are still short of sampling when using some solvents. How does the solvent concentration affect convergence? and the protein flexibility?
This small project aims to answer some of this questions and to provide real data to chose how long you want your MDMix simulations to run.


### Step 0. Test MDMix trajectories
The project will use a sample of mdmix trajectories from different proteins, solvents and lengths with the aim to cover as much as variability as possible.
Ideally we would have absurdily long MDMix simulations (i.e. >500 ns per replica) for each solvent. As I don't have this data at the moment, i'll leave all the code here so the analysis can be extended in the future.

Data available at the moment:

| System | Solvents | Nanoseconds | Replicas |
|--------|----------|-------------|----------|
| [SOCS1 ~ AlphaFold2](https://www.uniprot.org/uniprot/O15524) | ETA, MAM | 100 | 3 |
| [PFKFB3 ~ 6HVI](https://www.uniprot.org/uniprot/Q16875) | ETA, CLE, PYR5 | 100 | 3 |
| [p62 ZZ-Domain ~ 6KHZ](https://www.uniprot.org/uniprot/Q13501) | ETA, ION5, WAT | 200 | 3 |


### Step 1. Creating the density and energy grids
First of all we need to create the density/energy grids using increasingly big number of trajectories.

```{bash}
python lib/generate_cpptraj_scripts.py lib/inputs/cpptraj_sampling.yml
```

The input yaml file has the following format (i'll show the one for SOCS1_AF as an example)

```{yaml}
# Generic input to generate cpptraj scripts sampling different trajectories
Grid: # If you don't know the grid center, put the origin.
  - coordinates center: 9.5 11.0 -8.0
  - coordinates origin: False
  - delta: 0.5
  - dx: 156
  - dy: 156
  - dz: 156

Data:
  - replicas: 3
  - nanoseconds: 100
  - data directory: raw_data/SOCS1_AF
  - topologies: [SOCS1_AF_ETA.prmtop, SOCS1_AF_MAM.prmtop]
  - solvents: [ETA, MAM]
    ETA: # specify here the probes you want to test and their mask
      CT: :ETA@C1
      OH: :ETA@O1
      WAT: :WAT@O
    MAM:
      CT: :MAM@C2
      N: :MAM@N1
      O: :MAM@O1
      WAT: :WAT@O

Sampling:
  - Cross replica: True
  - Intra replica: False
  - Meta-replicas: 5
  - Sampling steps: 5 10 25 50 100 150 200 250 300
  - Output directory: lib/inputs/SOCS1_AF
  - Output grids directory: dgrids/SOCS1_AF
```
