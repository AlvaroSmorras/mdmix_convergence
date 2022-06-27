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
python lib/generate_cpptraj_scripts.py inputs/{system}_cpptraj_sampling.yml
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