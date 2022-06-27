# MDMix convergence anlysis:
## How many nanoseconds should I run?
The increase in computational power has made us launch longer and longer MDMix simulations, but is it necessary?
Maybe we are still short of sampling when using some solvents. How does the solvent concentration affect convergence? and the protein flexibility?
This small project aims to answer some of this questions and to provide real data to chose how long you want your MDMix simulations to run.


### Step 0. Test MDMix trajectories
The project will use a sample of mdmix trajectories from different proteins, solvents and lengths with the aim to cover as much as variability as possible.
Ideally we would have absurdily long MDMix simulations (i.e. >500 ns per replica) for each solvent. As I don't have this data at the moment, i'll leave all the code here so the analysis can be extended to other systems and solvents in the future.

Data available at the moment:

| System | Solvents | Nanoseconds | Replicas | Prot Size (aa) |
|--------|----------|-------------|----------|----------------|
| [SOCS1 ~ AlphaFold2](https://www.uniprot.org/uniprot/O15524) | ETA, MAM | 100 | 3 | 211 |
| [PFKFB3 ~ 6HVI](https://www.uniprot.org/uniprot/Q16875) | ETA, CLE, PYR5 | 100 | 3 | 446 |
| [p62 ZZ-Domain ~ 6KHZ](https://www.uniprot.org/uniprot/Q13501) | ETA, ION5, WAT | 200 | 3 | 56 |


### Step 1. Creating the density and energy grids
First of all we need to create the density/energy grids using increasingly big number of trajectories.

We'll calculate the density grids using cpptraj, so first we will generate the input files.
```{bash}
python lib/generate_cpptraj_scripts.py inputs/{system}\_cpptraj_sampling.yml
```

The input yaml file has the following format (i'll show the one for SOCS1_AF with detailed explanation as an example)

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

Once we have the input files, its time to launch cpptraj to obtain the files. The filesystem will have the following structure:
* cpptraj inputs in inputs/{system}/{solvent}/{metareplica}\_{sampling_method}\_{nanoseconds}.ptraj
* density grids in dgrids/{system}/{solvent}/{metareplica}\_{sampling_method}\_{nanoseconds}.dx

##### Cpptraj inputs to density grids
```{bash}
for input in inputs/*/*/*ptraj; do cpptraj -i $input; done
```

##### Density grids to energy grids (WIP)
```{bash}
```

### Step 2. Obtaining energies in certain points to check convergence
There are many options on how to check the convergence.
IDEAS:
- Check average discrepancy of dGB in each point across the grid.
- Check energy in certain points (hotspots).
      - Cpptraj have the option to return highest densities (>80\% (?)) we can use those.
      - Too many points together in hotspots, use clustering 