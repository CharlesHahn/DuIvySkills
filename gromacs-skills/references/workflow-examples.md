# GROMACS Workflow Examples

This document provides common simulation workflows with step-by-step examples.

## Workflow 1: Basic Protein Simulation

### Overview
Set up and run a basic protein simulation in water with ions.

### Prerequisites
- Protein structure file (PDB format)
- GROMACS installed
- Force field selected

### Step 1: Generate Topology

```bash
# Generate topology from PDB
gmx pdb2gmx -f protein.pdb -o protein_processed.gro -p topol.top -water tip3p

# Select force field (e.g., amber99sb-ildn)
# Select water model (e.g., TIP3P)
```

**Output files**:
- `protein_processed.gro`: Processed coordinates
- `topol.top`: Topology file
- `posre.itp`: Position restraints (optional)

### Step 2: Define Simulation Box

```bash
# Define cubic box with 1.0 nm padding
gmx editconf -f protein_processed.gro -o protein_box.gro -c -d 1.0 -bt cubic
```

**Output files**:
- `protein_box.gro`: Protein in box

### Step 3: Add Solvent

```bash
# Add water molecules
gmx solvate -cp protein_box.gro -cs spc216.gro -o protein_solv.gro -p topol.top
```

**Output files**:
- `protein_solv.gro`: Solvated system
- `topol.top`: Updated topology with water

### Step 4: Add Ions

```bash
# Generate .tpr for genion
gmx grompp -f mdp/ions.mdp -c protein_solv.gro -p topol.top -o ions.tpr

# Replace water with ions (neutralize + add salt)
echo "SOL" | gmx genion -s ions.tpr -p topol.top -o protein_ions.gro -pname NA -nname CL -neutral -conc 0.15
```

**Output files**:
- `protein_ions.gro`: System with ions
- `topol.top`: Updated topology with ions

### Step 5: Energy Minimization

```bash
# Generate .tpr for energy minimization
gmx grompp -f mdp/minim.mdp -c protein_ions.gro -p topol.top -o em.tpr -maxwarn 2

# Run energy minimization
gmx mdrun -s em.tpr -deffnm em
```

**Output files**:
- `em.gro`: Minimized structure
- `em.trr`: Trajectory
- `em.edr`: Energy file
- `em.log`: Log file

**Check energy minimization**:
```bash
# Extract potential energy
echo "Potential" | gmx energy -f em.edr -o em_potential.xvg

# Check if potential energy has converged (should reach minimum)
```

### Step 6: NVT Equilibration

```bash
# Generate .tpr for NVT equilibration
gmx grompp -f mdp/nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr

# Run NVT equilibration
gmx mdrun -s nvt.tpr -deffnm nvt
```

**Output files**:
- `nvt.gro`: Equilibrated structure
- `nvt.trr`: Trajectory
- `nvt.edr`: Energy file
- `nvt.log`: Log file

**Check temperature**:
```bash
# Extract temperature
echo "Temperature" | gmx energy -f nvt.edr -o temperature.xvg

# Temperature should stabilize around target value
```

### Step 7: NPT Equilibration

```bash
# Generate .tpr for NPT equilibration
gmx grompp -f mdp/npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr

# Run NPT equilibration
gmx mdrun -s npt.tpr -deffnm npt
```

**Output files**:
- `npt.gro`: Equilibrated structure
- `npt.trr`: Trajectory
- `npt.edr`: Energy file
- `npt.log`: Log file

**Check pressure and density**:
```bash
# Extract pressure
echo "Pressure" | gmx energy -f npt.edr -o pressure.xvg

# Extract density
echo "Density" | gmx energy -f npt.edr -o density.xvg

# Pressure and density should stabilize
```

### Step 8: Production MD

```bash
# Generate .tpr for production MD
gmx grompp -f mdp/md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr

# Run production MD
gmx mdrun -s md.tpr -deffnm md
```

**Output files**:
- `md.gro`: Final structure
- `md.trr`: Trajectory (full precision)
- `md.xtc`: Trajectory (compressed)
- `md.edr`: Energy file
- `md.log`: Log file
- `md.cpt`: Checkpoint file

**Check energy conservation**:
```bash
# Extract total energy
echo "Total-Energy" | gmx energy -f md.edr -o total_energy.xvg

# Energy should be conserved (small fluctuations)
```

## Workflow 2: Protein-Ligand Complex

### Overview
Set up simulation of protein-ligand complex.

### Step 1: Prepare Ligand

```bash
# Generate ligand topology
gmx pdb2gmx -f ligand.pdb -o ligand.gro -p ligand.top -ff amber99sb-ildn

# Or use external tools (e.g., antechamber, CGenFF)
```

### Step 2: Merge Protein and Ligand

```bash
# Combine protein and ligand coordinates
gmx editconf -f protein.gro -o protein_ligand.gro

# Manually edit file or use cat to add ligand
# Or use PyMOL/VMD to merge structures

# Merge topologies manually or using specialized tools
```

### Step 3: Continue from Workflow 1

After merging, follow steps 2-8 from Workflow 1.

### Step 4: Create Index File

```bash
# Create index file with protein and ligand groups
gmx make_ndx -f md.tpr -o index.ndx

# Create protein-ligand group
echo "Protein | Ligand\nq" | gmx make_ndx -f md.tpr -o index.ndx
```

## Workflow 3: Trajectory Analysis

### Overview
Analyze simulation trajectory.

### Step 1: Check Trajectory

```bash
# Check trajectory files
gmx check -s md.tpr -f md.xtc -e md.edr
```

### Step 2: PBC Correction

```bash
# Create index file
gmx make_ndx -f md.tpr -o index.ndx

# Add center group (find center atom)
echo "Protein" | dit find_center -f npt.gro index.ndx

# Add center group to index file
echo "\n[ center ]\n<atom_number>\n" >> index.ndx

# Apply PBC correction
echo "center\nProtein" | gmx trjconv -s md.tpr -f md.xtc -o center.xtc -n index.ndx -pbc atom -center

# Keep molecules intact
echo "Protein" | gmx trjconv -s md.tpr -f center.xtc -o mol.xtc -n index.ndx -pbc mol -ur compact

# Remove translation and rotation
echo "Backbone\nProtein" | gmx trjconv -s md.tpr -f mol.xtc -o fit.xtc -n index.ndx -fit rot+trans
```

### Step 3: RMSD Analysis

```bash
# Calculate RMSD relative to initial structure
echo "Backbone\nProtein" | gmx rms -s md.tpr -f fit.xtc -o rmsd.xvg

# Calculate RMSD relative to average structure
echo "Backbone\nBackbone" | gmx rms -s md.tpr -f fit.xtc -o rmsd_avg.xvg -ref average.pdb
```

### Step 4: RMSF Analysis

```bash
# Calculate per-atom RMSF
echo "Backbone" | gmx rmsf -s md.tpr -f fit.xtc -o rmsf.xvg

# Calculate per-residue RMSF
echo "Backbone" | gmx rmsf -s md.tpr -f fit.xtc -o rmsf_res.xvg -res
```

### Step 5: Radius of Gyration

```bash
# Calculate radius of gyration
echo "Protein" | gmx gyrate -s md.tpr -f fit.xtc -o rg.xvg
```

### Step 6: Hydrogen Bond Analysis

```bash
# Analyze intra-protein hydrogen bonds
echo "Protein\nProtein" | gmx hbond -s md.tpr -f fit.xtc -num hb.xvg

# Analyze protein-water hydrogen bonds
echo "Protein\nWater" | gmx hbond -s md.tpr -f fit.xtc -num hb_water.xvg
```

### Step 7: Secondary Structure

```bash
# Assign secondary structure
echo "Protein" | gmx do_dssp -s md.tpr -f fit.xtc -sc ss.xpm -ss ss.xvg
```

### Step 8: Energy Analysis

```bash
# Extract potential energy
echo "Potential" | gmx energy -f md.edr -o potential.xvg

# Extract temperature
echo "Temperature" | gmx energy -f npt.edr -o temperature.xvg

# Extract pressure
echo "Pressure" | gmx energy -f npt.edr -o pressure.xvg

# Extract density
echo "Density" | gmx energy -f npt.edr -o density.xvg
```

### Step 9: Visualization with DuIvyTools

After analysis, use DuIvyTools to visualize GROMACS output files:

```bash
# Visualize RMSD
dit xvg_show -f rmsd.xvg -x "Time (ns)" -y "RMSD (nm)" -t "RMSD Analysis"

# Visualize RMSF
dit xvg_show -f rmsf_res.xvg -x "Residue Number" -y "RMSF (nm)" -t "Per-Residue RMSF"

# Visualize radius of gyration
dit xvg_show -f rg.xvg -x "Time (ns)" -y "Rg (nm)" -t "Radius of Gyration"

# Visualize hydrogen bonds
dit xvg_show -f hb.xvg -x "Time (ns)" -y "Number of H-bonds" -t "Hydrogen Bonds"

# Visualize secondary structure matrix
dit xpm_show -f ss.xpm -t "Secondary Structure" -cmap nipy_spectral

# Visualize energy components
dit xvg_show -f potential.xvg -x "Time (ps)" -y "Potential Energy (kJ/mol)" -t "Potential Energy"
dit xvg_show -f temperature.xvg -x "Time (ps)" -y "Temperature (K)" -t "Temperature"

# Compare multiple analyses
dit xvg_compare -f rmsd.xvg rg.xvg -c 1 1 -l "RMSD" "Rg" -x "Time (ns)"
```

**Note**: Use the `duivytools-helper` skill for comprehensive visualization options and advanced plotting features.

## Workflow 4: PCA Analysis

### Overview
Perform principal component analysis.

### Step 1: Calculate Covariance Matrix

```bash
# Calculate covariance matrix
echo "C-alpha\nC-alpha" | gmx covar -s md.tpr -f fit.xtc -o eigenvalues.xvg -v eigenvectors.trr -xpma covar.xpm
```

### Step 2: Analyze Eigenvalues

```bash
# Plot eigenvalues
dit xvg_show -f eigenvalues.xvg -x "PC Number" -y "Eigenvalue" -t "Eigenvalue Spectrum"

# Calculate contribution percentages
# Total variance = sum of all eigenvalues
# PC contribution = (individual eigenvalue / total variance) × 100%
```

### Step 3: Project Trajectory

```bash
# Project onto PC1
echo "C-alpha\nC-alpha" | gmx anaeig -s md.tpr -f fit.xtc -v eigenvectors.trr -first 1 -last 1 -proj pc1.xvg

# Project onto PC2
echo "C-alpha\nC-alpha" | gmx anaeig -s md.tpr -f fit.xtc -v eigenvectors.trr -first 2 -last 2 -proj pc2.xvg

# Generate 2D projection
echo "C-alpha\nC-alpha" | gmx anaeig -s md.tpr -f fit.xtc -v eigenvectors.trr -first 1 -last 2 -2d 2dproj.xvg
```

### Step 4: Extract Extreme Conformations

```bash
# Extract extreme conformations along PC1
echo "C-alpha\nC-alpha" | gmx anaeig -s md.tpr -f fit.xtc -v eigenvectors.trr -first 1 -last 1 -extreme pc1_extreme.pdb

# Extract average structure
echo "C-alpha\nC-alpha" | gmx anaeig -s md.tpr -f fit.xtc -v eigenvectors.trr -first 1 -last 1 -average average.pdb
```

### Step 5: Visualization with DuIvyTools

```bash
# Visualize eigenvalues
dit xvg_show -f eigenvalues.xvg -x "PC Number" -y "Eigenvalue" -t "Eigenvalue Spectrum"

# Visualize PC projections
dit xvg_show -f pc1.xvg -x "Time (ps)" -y "PC1 Coordinate" -t "PC1 Projection"
dit xvg_show -f pc2.xvg -x "Time (ps)" -y "PC2 Coordinate" -t "PC2 Projection"

# Visualize PC1 vs PC2 scatter plot
dit xvg_show_scatter -f pc1.xvg pc2.xvg -c 1 1 -x "PC1" -y "PC2" -t "PC1 vs PC2"

# Visualize 2D projection
dit xvg_show -f 2dproj.xvg -x "Time (ps)" -y "PC1" -y2 "PC2" -t "2D Projection"
```

**Note**: Use the `duivytools-helper` skill for comprehensive visualization options and advanced plotting features.

## Workflow 5: Free Energy Landscape

### Overview
Generate free energy landscape.

### Method 1: RMSD + Gyrate

```bash
# Calculate RMSD
echo "Backbone\nProtein" | gmx rms -s md.tpr -f fit.xtc -o rmsd.xvg

# Calculate radius of gyration
echo "Protein" | gmx gyrate -s md.tpr -f fit.xtc -o gyrate.xvg

# Combine data
dit xvg_combine -f rmsd.xvg gyrate.xvg -c 0,1 1 -o sham.xvg

# Generate FEL
gmx sham -tsham 310 -nlevels 100 -f sham.xvg -ls gibbs.xpm -g gibbs.log

# Visualize
dit xpm_show -f gibbs.xpm -m contour -cmap jet -o fel.png -x "RMSD (nm)" -y "Rg (nm)"
```

### Method 2: Principal Components

```bash
# Calculate PCA (see Workflow 4)

# Combine PC projections
dit xvg_combine -f pc1.xvg pc2.xvg -c 0,1 1 -o sham.xvg

# Generate FEL
gmx sham -tsham 310 -nlevels 100 -f sham.xvg -ls gibbs.xpm -g gibbs.log

# Visualize
dit xpm_show -f gibbs.xpm -m contour -cmap jet -o fel.png -x "PC1" -y "PC2" -t "Free Energy Landscape"
```

### Extract Lowest Energy Conformations

```bash
# View energy minima
cat gibbs.log

# Check frames corresponding to minima
cat bindex.ndx

# Extract conformation at specific time
echo "Protein" | gmx trjconv -f md.xtc -s md.tpr -b <time> -e <time> -o min1.pdb
```

## Workflow 6: Membrane Protein Simulation

### Overview
Set up simulation of membrane protein.

### Step 1: Prepare Protein

```bash
# Generate topology
gmx pdb2gmx -f protein.pdb -o protein.gro -p topol.top -water tip3p
```

### Step 2: Insert into Membrane

```bash
# Use membrane insertion tool (e.g., CHARMM-GUI, insane.py, MemGen)
# Or manually place protein in membrane

# Example using insane.py:
insane.py -f protein.gro -o membrane.gro -p topol.top -l POPC -u POPC -x 10 -y 10 -z 10 -sol TIP3P -salt 0.15
```

### Step 3: Energy Minimization

```bash
# Generate .tpr
gmx grompp -f mdp/minim_membrane.mdp -c membrane.gro -p topol.top -o em.tpr -maxwarn 2

# Run minimization
gmx mdrun -s em.tpr -deffnm em
```

### Step 4: Equilibration

```bash
# NVT equilibration
gmx grompp -f mdp/nvt_membrane.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -s nvt.tpr -deffnm nvt

# NPT equilibration (semisotropic pressure)
gmx grompp -f mdp/npt_membrane.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -s npt.tpr -deffnm npt
```

### Step 5: Production MD

```bash
# Generate .tpr
gmx grompp -f mdp/md_membrane.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr

# Run production
gmx mdrun -s md.tpr -deffnm md
```

## Workflow 7: Umbrella Sampling

### Overview
Calculate free energy profile using umbrella sampling.

### Step 1: Steered MD

```bash
# Generate .tpr for pulling
gmx grompp -f mdp/pull.mdp -c complex.gro -p topol.top -o pull.tpr

# Run steered MD
gmx mdrun -s pull.tpr -deffnm pull
```

### Step 2: Extract Windows

```bash
# Extract frames at regular intervals
gmx trjconv -s pull.tpr -f pull.xtc -o window_0.gro -dump 100
gmx trjconv -s pull.tpr -f pull.xtc -o window_1.gro -dump 200
# ... etc
```

### Step 3: Run Umbrella Sampling

```bash
# For each window
for i in {0..20}; do
    # Generate .tpr with specific pull coordinate
    gmx grompp -f mdp/umbrella.mdp -c window_${i}.gro -p topol.top -o umbrella_${i}.tpr

    # Run umbrella sampling
    gmx mdrun -s umbrella_${i}.tpr -deffnm umbrella_${i}
done
```

### Step 4: WHAM Analysis

```bash
# Create tpr list
ls umbrella_*.tpr > tpr_files.dat

# Create pullf file list
ls umbrella_*/pullf.xvg > pullf_files.dat

# Run WHAM
gmx wham -it tpr_files.dat -if pullf_files.dat -o profile.xvg
```

## Workflow 8: Replica Exchange

### Overview
Run replica exchange molecular dynamics.

### Step 1: Generate .tpr Files for Different Temperatures

```bash
# For each temperature
for temp in 300 310 320 330 340; do
    # Generate .mdp with specific temperature
    sed "s/ref_t.*/ref_t = $temp/" mdp/remd.mdp > remd_${temp}.mdp

    # Generate .tpr
    gmx grompp -f remd_${temp}.mdp -c npt.gro -p topol.top -o remd_${temp}.tpr
done
```

### Step 2: Run REMD

```bash
# Run replica exchange
gmx mdrun -multidir remd_300 remd_310 remd_320 remd_330 remd_340 -deffnm remd -replex 1000
```

## Workflow 9: Clustering Analysis

### Overview
Cluster trajectory conformations.

### Step 1: Preprocess Trajectory

```bash
# Fit trajectory
echo "Backbone\nProtein" | gmx trjconv -s md.tpr -f md.xtc -o fit.xtc -fit rot+trans

# Remove rotation and translation
echo "Backbone\nProtein" | gmx trjconv -s md.tpr -f fit.xtc -o nofit.xtc
```

### Step 2: Calculate RMSD Matrix

```bash
# Calculate distance matrix
echo "Backbone" | gmx cluster -s md.tpr -f nofit.xtc -dm rmsd.xpm -method linkage -cutoff 0.2
```

### Step 3: Cluster Analysis

```bash
# Perform clustering
echo "Backbone" | gmx cluster -s md.tpr -f nofit.xtc -dm rmsd.xpm -method gromos -cutoff 0.2 -cl clusters.xpm -o cluster_size.xvg
```

### Step 4: Extract Cluster Representatives

```bash
# Extract representative structures from each cluster
echo "Backbone" | gmx cluster -s md.tpr -f nofit.xtc -method gromos -cutoff 0.2 -clid cluster_id -cl cluster_rep.pdb
```

## Workflow 10: Restart Simulation

### Overview
Restart simulation from checkpoint.

### Step 1: Checkpoint File

```bash
# Checkpoint file (.cpt) is created automatically
# Default: every 15 minutes
```

### Step 2: Continue Simulation

```bash
# Use .cpt file as input
gmx grompp -f mdp/md.mdp -c md.cpt -t md.cpt -e md.edr -p topol.top -o md_restart.tpr

# Run continued simulation
gmx mdrun -s md_restart.tpr -deffnm md_restart
```

### Step 3: Merge Trajectories

```bash
# Merge original and continued trajectories
gmx trjcat -f md.xtc md_restart.xtc -o md_combined.xtc -s md.tpr

# Merge energy files
gmx eneconv -f md.edr md_restart.edr -o md_combined.edr
```

## Common MDP File Parameters

### Energy Minimization (minim.mdp)
```
integrator  = steep
emtol       = 1000.0
emstep      = 0.01
nsteps      = 50000
nstlist     = 10
cutoff-scheme = Verlet
ns_type     = grid
rlist       = 1.0
rcoulomb    = 1.0
rvdw        = 1.0
pbc         = xyz
constraints = none
```

### NVT Equilibration (nvt.mdp)
```
integrator  = md
dt          = 0.002
nsteps      = 50000
nstlist     = 10
cutoff-scheme = Verlet
ns_type     = grid
rlist       = 1.0
rcoulomb    = 1.0
rvdw        = 1.0
pbc         = xyz
constraints = h-bonds
continuation = no
constraint_algorithm = lincs
temperature coupling = V-rescale
tc-grps     = Protein Non-Protein
tau_t       = 0.1     0.1
ref_t       = 300     300
```

### NPT Equilibration (npt.mdp)
```
integrator  = md
dt          = 0.002
nsteps      = 50000
nstlist     = 10
cutoff-scheme = Verlet
ns_type     = grid
rlist       = 1.0
rcoulomb    = 1.0
rvdw        = 1.0
pbc         = xyz
constraints = h-bonds
continuation = no
constraint_algorithm = lincs
temperature coupling = V-rescale
tc-grps     = Protein Non-Protein
tau_t       = 0.1     0.1
ref_t       = 300     300
pressure coupling = Parrinello-Rahman
pcoupl      = Parrinello-Rahman
pcoupltype  = isotropic
tau_p       = 2.0
ref_p       = 1.0
compressibility = 4.5e-5
```

### Production MD (md.mdp)
```
integrator  = md
dt          = 0.002
nsteps      = 50000000
nstlist     = 10
cutoff-scheme = Verlet
ns_type     = grid
rlist       = 1.0
rcoulomb    = 1.0
rvdw        = 1.0
pbc         = xyz
constraints = h-bonds
continuation = yes
constraint_algorithm = lincs
temperature coupling = V-rescale
tc-grps     = Protein Non-Protein
tau_t       = 0.1     0.1
ref_t       = 300     300
pressure coupling = Parrinello-Rahman
pcoupl      = Parrinello-Rahman
pcoupltype  = isotropic
tau_p       = 2.0
ref_p       = 1.0
compressibility = 4.5e-5
```

## Tips and Best Practices

1. **Always backup** files before running simulations
2. **Check energy minimization** before proceeding
3. **Monitor equilibration** (temperature, pressure, density)
4. **Verify energy conservation** in production runs
5. **Use appropriate time steps** (0.002 ps with H-bonds)
6. **Adjust cutoff distances** for your system
7. **Document all parameters** for reproducibility
8. **Test on small systems** before large simulations
9. **Use checkpoint files** for long simulations
10. **Monitor performance** and optimize parameters
11. **Visualize results** using DuIvyTools for publication-quality figures
12. **Check output files** (.xvg, .xpm) after each analysis step

## Visualization with DuIvyTools

After GROMACS analysis, use DuIvyTools to visualize output files:

**For .xvg files** (time-series data):
```bash
# Basic plotting
dit xvg_show -f rmsd.xvg -x "Time (ns)" -y "RMSD (nm)"

# Multiple file comparison
dit xvg_compare -f rmsd1.xvg rmsd2.xvg -c 1 1 -l "Run1" "Run2"

# Statistical analysis
dit xvg_ave -f rmsd.xvg -o rmsd_stats.xvg
```

**For .xpm files** (matrices, heatmaps):
```bash
# Basic heatmap
dit xpm_show -f dccm.xpm -m contour -cmap bwr

# 3D visualization
dit xpm_show -f fel.xpm -m 3d -eg plotly

# Difference between matrices
dit xpm_diff -f matrix1.xpm matrix2.xpm -o diff.xpm
```

**Call the `duivytools-helper` skill** for comprehensive visualization options and advanced plotting features.

## Troubleshooting

### Energy Not Converging in Minimization
- Increase nsteps
- Check for clashes
- Use different minimization algorithm (steepest descent vs conjugate gradient)

### Temperature Not Stabilizing
- Increase equilibration time
- Check thermostat parameters
- Verify system is properly solvated

### Pressure Not Stabilizing
- Increase equilibration time
- Check barostat parameters
- Verify box size is appropriate

### Energy Drift in Production
- Check energy conservation
- Verify dt is not too large
- Check constraints
- Reduce time step if needed

### Domain Decomposition Fails
- Reduce number of MPI ranks
- Adjust domain decomposition parameters
- Use fewer threads

## Getting Help

For detailed parameters, use:
```bash
gmx <command> -h
```

For version-specific documentation:
```bash
web_search: "site:manual.gromacs.org gmx <command> <version>"
```