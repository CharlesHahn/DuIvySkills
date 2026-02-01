---
name: gromacs
description: "Comprehensive guide for GROMACS molecular dynamics simulation package. Use when Claude needs to execute GROMACS commands or understand GROMACS functionality for: (1) Setting up molecular dynamics simulations (topology generation, energy minimization, equilibration, production), (2) Analyzing simulation trajectories (RMSD, RMSF, hydrogen bonds, distances, clustering), (3) Processing trajectory files (format conversion, PBC correction, fitting), (4) Extracting and analyzing energy and thermodynamic data, (5) Creating index groups and selecting atoms, (6) Any GROMACS-related simulation or analysis tasks"
---

# GROMACS

GROMACS is a versatile package to perform molecular dynamics, i.e. simulate the Newtonian equations of motion for systems with hundreds to millions of particles. This skill provides command reference and workflow guidance for GROMACS simulations and analysis.

## Quick Start

### Check GROMACS Version

```bash
gmx --version
```

Note the version number for compatibility and documentation lookup.

### Get Command Help

**Priority 1: Local help (fastest, version-matched)**
```bash
gmx <command> -h
```

**Priority 2: Online documentation (detailed, official)**
```bash
# Search with version number
web_search: "site:manual.gromacs.org gmx <command> <version>"
# Example: "site:manual.gromacs.org gmx rms 2024.3"
```

**Always use local help first** for detailed parameters and options.

## Command Categories

### 1. Topology and Structure

| Command | Description |
|---------|-------------|
| `pdb2gmx` | Generate GROMACS topology from PDB structure |
| `editconf` | Edit structure files (box definition, dimensions) |
| `solvate` | Add solvent molecules to simulation box |
| `insert-molecules` | Insert molecules into simulation box |
| `genrestr` | Generate position restraint files |

### 2. Simulation Setup

| Command | Description |
|---------|-------------|
| `grompp` | Generate simulation run input file (.tpr) from topology and parameters |
| `mdrun` | Run molecular dynamics simulation |
| `mdrun-mpi` | Run simulation with MPI parallelization |

### 3. Energy Analysis

| Command | Description |
|---------|-------------|
| `energy` | Extract energy components from energy file |
| `eneconv` | Convert energy file formats |
| `bar` | Free energy calculations using Bennett acceptance ratio |

### 4. Trajectory Analysis

| Command | Description |
|---------|-------------|
| `rms` | Calculate root mean square deviation (RMSD) |
| `rmsf` | Calculate root mean square fluctuation (RMSF) |
| `gyrate` | Calculate radius of gyration |
| `hbond` | Analyze hydrogen bonds |
| `distance` | Calculate distances between atom groups |
| `angle` | Calculate angles between atom groups |
| `dihedral` | Calculate dihedral angles |
| `cluster` | Cluster analysis of conformations |
| `mindist` | Calculate minimum distance between groups |
| `sasa` | Calculate solvent accessible surface area |
| `principal` | Principal component analysis of molecules |
| `do_dssp` | Secondary structure analysis (DSSP) |

### 5. Structural Analysis

| Command | Description |
|---------|-------------|
| `covar` | Covariance matrix analysis (for PCA) |
| `anaeig` | Analyze eigenvectors from covariance analysis |
| `mdmat` | Calculate distance matrix between residues |
| `sham` | Generate free energy landscapes |
| `order` | Calculate order parameters |
| `rotacf` | Calculate rotational autocorrelation functions |
| `dielectric` | Calculate dielectric properties |

### 6. Trajectory Processing

| Command | Description |
|---------|-------------|
| `trjconv` | Convert trajectory formats, apply PBC correction, fit structures |
| `trjcat` | Concatenate multiple trajectory files |
| `trjorder` | Reorder atoms in trajectory |
| `dump` | Dump trajectory to text format |
| `eneconv` | Convert energy file formats |

### 7. Index and Selection

| Command | Description |
|---------|-------------|
| `make_ndx` | Create or edit index files |
| `select` | Select atoms using selection language |
| `genrestr` | Generate position restraints |
| `genion` | Replace water molecules with ions |

### 8. Tools and Utilities

| Command | Description |
|---------|-------------|
| `xpm2ps` | Convert XPM matrix files to PostScript |
| `x2top` | Generate topology from coordinates |
| `check` | Check simulation files for errors |
| `wham` | Weighted histogram analysis method |
| `tune_pme` | Tune PME parameters for performance |

## Common Parameters

### Input/Output Files

- `-f INPUT`: Input trajectory/coordinate file
- `-s TOPOLOGY`: Input topology file (.tpr)
- `-n INDEX`: Input index file (.ndx)
- `-c INPUT`: Input structure file
- `-o OUTPUT`: Output file
- `-p TOPOL`: Output topology file

### Time Selection

- `-b TIME`: Start time (ps)
- `-e TIME`: End time (ps)
- `-dt TIME`: Time step for output (ps)

### Trajectory Processing

- `-pbc TYPE`: Periodic boundary condition handling (none, mol, atom, com, nojump, cluster)
- `-center TYPE`: Center coordinates
- `-fit TYPE`: Fit trajectory (none, xyz, rotxy, xy, progression, transxy, rot+trans)
- `-ur TYPE`: Unit cell representation (compact, triclinic, rect)

### Analysis Options

- `-type TYPE`: Calculation type
- `-ng NUMBER`: Number of groups
- `-frame FRAME`: Specific frame number

### Output Control

- `-nice LEVEL`: Nice level for CPU usage
- `-nt NUMBER`: Number of threads
| `-ntmpi NUMBER`: Number of MPI threads
| `-mpi`: Enable MPI
| `-nb TYPE`: Neighbor searching type (cpu, gpu)

### Visualization

- `-xvg TYPE`: Plot format (none, xmgrace, xmgr, matplotlib)
| `-xvg FORMAT`: XVG format (none, xmgrace, xmgr, matplotlib)

## Version Detection

### Get Version Information

```bash
gmx --version
```

Example output:
```
GROMACS version:    2024.3
GROMACS modification: 2024.3-1
Precision: mixed
MPI: enabled
OpenMP: enabled
GPU support: enabled
FFT library: fftw-MPI
```

### Check for Specific Features

```bash
# Check for GPU support
gmx mdrun -h | grep -i gpu

# Check for MPI support
gmx --version | grep -i mpi

# Check for OpenMP support
gmx --version | grep -i openmp
```

## Getting Help

### Method 1: Local Help (Recommended)

Get help for any command:

```bash
gmx <command> -h
```

Example:
```bash
gmx rms -h
gmx mdrun -h
gmx pdb2gmx -h
```

**When to use local help**:
- User has GROMACS installed
- Need version-specific parameters
- Quick parameter lookup
- Offline work

### Method 2: Online Documentation

Search GROMACS official manual:

```bash
# General search
web_search: "site:manual.gromacs.org gmx <command>"

# Version-specific search
web_search: "site:manual.gromacs.org gmx <command> <version>"

# Examples
web_search: "site:manual.gromacs.org gmx rms 2024.3"
web_search: "site:manual.gromacs.org gmx mdrun 2023"
web_search: "site:manual.gromacs.org gmx pdb2gmx current"
```

**When to use online search**:
- User doesn't have GROMACS installed
- Need detailed explanations and examples
- Need to compare parameters across versions
- Troubleshooting specific issues

### Help Workflow

1. **User requests GROMACS task**
2. **Check GROMACS version**: `gmx --version`
3. **Get basic info** from this skill's command categories
4. **Get detailed parameters**: `gmx <command> -h`
5. **If more detail needed**: Search online documentation with version

## Common Workflows

### 1. Basic Simulation Setup

```
pdb2gmx → editconf → solvate → genion → grompp → mdrun
```

### 2. Energy Minimization

```
grompp (emin.mdp) → mdrun → energy
```

### 3. Equilibration (NVT/NPT)

```
grompp (nvt.mdp) → mdrun → energy
grompp (npt.mdp) → mdrun → energy
```

### 4. Production Run

```
grompp (md.mdp) → mdrun → analysis
```

### 5. Trajectory Analysis

```
trjconv (PBC correction) → rms/rmsf/gyrate/hbond → visualization with DuIvyTools
```

### 6. PCA Analysis

```
trjconv (fit) → covar → anaeig → visualization with DuIvyTools
```

**Note**: Use `duivytools-skills` skill for visualizing all GROMACS output files (.xvg, .xpm)

## File Formats

### Input Formats

- **.pdb**: Protein Data Bank format
- **.gro**: GROMACS coordinate format
- **.tpr**: GROMACS portable run input (topology + parameters)
- **.xtc**: Compressed trajectory (lossy, good for long simulations)
- **.trr**: Full precision trajectory
- **.ndx**: Index file (atom groups)
- **.mdp**: Molecular dynamics parameter file
- **.top**: Topology file

### Output Formats

- **.xvg**: Grace/XVG plot format (ASCII) - Time-series data, energies, analysis results
- **.xpm**: Pixel map format (matrices, heatmaps) - 2D matrices for correlation, contact maps, free energy landscapes
- **.edr**: Energy file (binary)
- **.log**: Log file
- **.cpt**: Checkpoint file (for continuation)

### Visualization of GROMACS Output Files

**For visualization of GROMACS output files (.xvg, .xpm, etc.), use DuIvyTools:**

Call the `duivytools-skills` skill when you need to:
- **Visualize .xvg files**: Plot RMSD, RMSF, energies, hydrogen bonds, distances, angles, PC projections, etc.
- **Visualize .xpm files**: Display matrices, heatmaps (DCCM, distance contact maps, free energy landscapes, secondary structure)
- **Analyze data**: Calculate averages, distributions, correlations from XVG files
- **Generate publication-quality plots**: Create professional figures with customizable styles

**Common visualization examples**:
```bash
# Visualize XVG data (RMSD, RMSF, energy, etc.)
dit xvg_show -f rmsd.xvg -x "Time (ns)" -y "RMSD (nm)"

# Visualize XPM matrix (DCCM, FEL, etc.)
dit xpm_show -f dccm.xpm -m contour -cmap bwr -zmin -1 -zmax 1

# Compare multiple XVG files
dit xvg_compare -f rmsd1.xvg rmsd2.xvg -c 1 1 -l "Run1" "Run2"
```

**When to use DuIvyTools**:
- After any GROMACS analysis that produces .xvg or .xpm files
- When user requests visualization or plots
- For publication-quality figures with customizable styling
- For statistical analysis of GROMACS output data

## Performance Optimization

### Parallel Execution

```bash
# OpenMP (multithreading)
gmx mdrun -nt 4 -s topol.tpr -deffnm md

# MPI (multinode)
mpirun -np 4 gmx_mpi mdrun -s topol.tpr -deffnm md

# Hybrid MPI+OpenMP
mpirun -np 2 gmx_mpi mdrun -nt 2 -s topol.tpr -deffnm md
```

### GPU Acceleration

```bash
# GPU for non-bonded interactions
gmx mdrun -ntomp 4 -nb gpu -s topol.tpr -deffnm md

# GPU for update
gmx mdrun -ntomp 4 -nb gpu -update gpu -s topol.tpr -deffnm md
```

### Memory Optimization

```bash
# Verlet buffer tuning
gmx mdrun -ntomp 4 -nb gpu -tunepme -s topol.tpr -deffnm md
```

## Troubleshooting

### Common Issues

**"Fatal error: No such group"**: Index group not found
- Solution: Create proper index file with `make_ndx`

**"Back Off! I just backed up md.log"**: Simulation crashed
- Solution: Check energy conservation, timestep, constraints

**"Fatal error: Domain decomposition error"**: Parallel setup issue
- Solution: Adjust MPI ranks or domain decomposition parameters

**"Fatal error: Number of coordinates in coordinate file does not match topology"**: Mismatched files
- Solution: Regenerate topology or coordinate file

## Best Practices

- **Always check version** before using commands
- **Use local help first** for accurate parameters
- **Test on small systems** before large simulations
- **Monitor energy conservation** during production runs
- **Back up files** before modifications
- **Document all parameters** and settings
- **Validate results** with multiple analysis methods
- **Use appropriate output frequency** to balance performance and detail
- **Visualize results** using DuIvyTools for publication-quality figures
- **Check output files** (.xvg, .xpm) after each analysis step

## Reference Documentation

For detailed information, consult these references:

- **[Command Categories](references/command-categories.md)** - Complete list of GROMACS commands by category
- **[Common Parameters](references/common-parameters.md)** - Detailed explanation of frequently used parameters
- **[Version Compatibility](references/version-compatibility.md)** - Version differences and migration guide
- **[Workflow Examples](references/workflow-examples.md)** - Common simulation workflows with examples

## Additional Resources

- **Official Documentation**: https://manual.gromacs.org/
- **GitHub Repository**: https://github.com/gromacs/gromacs
- **User Forum**: https://gromacs.bioexcel.eu/
- **Tutorials**: http://www.mdtutorials.com/gmx/
- **Citation**: Abraham et al., "GROMACS: High performance molecular simulations through multi-level parallelism from laptops to supercomputers", SoftwareX 1-2, 19-25 (2015)

## Important Notes

- **Version matters**: Commands and parameters can vary between GROMACS versions
- **Always check help**: Run `gmx <command> -h` for accurate parameter information
- **Context awareness**: This skill provides overview and workflow guidance; detailed parameters should be obtained via local help or online documentation
- **Performance tuning**: Optimize parameters for your hardware and system size
- **Validation**: Always validate simulation results with appropriate analysis