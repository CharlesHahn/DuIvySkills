# Common GROMACS Parameters

This document explains frequently used GROMACS parameters across different commands.

## Input/Output File Parameters

### `-f INPUT`
**Input trajectory or structure file**

- **Purpose**: Specify primary input file
- **Common formats**: .xtc, .trr, .gro, .pdb, .edr
- **Usage examples**:
  ```bash
  gmx rms -s md.tpr -f md.xtc -o rmsd.xvg
  gmx energy -f md.edr -o energy.xvg
  ```

### `-s INPUT`
**Input topology file (.tpr)**

- **Purpose**: Specify portable run input file containing topology and parameters
- **Format**: .tpr only
- **Usage examples**:
  ```bash
  gmx rms -s md.tpr -f md.xtc -o rmsd.xvg
  gmx trjconv -s md.tpr -f md.xtc -o fit.xtc
  ```

### `-n INPUT`
**Input index file (.ndx)**

- **Purpose**: Specify index file with atom groups
- **Format**: .ndx only
- **Usage examples**:
  ```bash
  echo "Protein" | gmx rms -s md.tpr -f md.xtc -n index.ndx -o rmsd.xvg
  gmx make_ndx -f md.tpr -o index.ndx
  ```

### `-c INPUT`
**Input structure file**

- **Purpose**: Specify coordinate file (used in grompp, editconf)
- **Common formats**: .gro, .pdb
- **Usage examples**:
  ```bash
  gmx editconf -f protein.pdb -o protein.gro
  gmx grompp -f mdp/minim.mdp -c ions.gro -p topol.top -o em.tpr
  ```

### `-p INPUT`
**Input topology file (.top)**

- **Purpose**: Specify topology file (used in grompp, genion)
- **Format**: .top only
- **Usage examples**:
  ```bash
  gmx grompp -f mdp/minim.mdp -c ions.gro -p topol.top -o em.tpr
  echo "SOL" | gmx genion -s ions.tpr -p topol.top -o ions.gro -pname NA -nname CL -neutral
  ```

### `-o OUTPUT`
**Output file**

- **Purpose**: Specify primary output file
- **Format depends on command**: .xvg, .xtc, .trr, .gro, .pdb, .edr, etc.
- **Usage examples**:
  ```bash
  gmx rms -s md.tpr -f md.xtc -o rmsd.xvg
  gmx trjconv -s md.tpr -f md.xtc -o fit.xtc
  ```

### `-deffnm BASENAME`
**Default filename base**

- **Purpose**: Set default filename base for output files
- **Common usage**: In mdrun to simplify output naming
- **Usage examples**:
  ```bash
  gmx mdrun -s em.tpr -deffnm em  # Creates em.log, em.edr, em.trr, etc.
  gmx mdrun -s md.tpr -deffnm md  # Creates md.log, md.edr, md.xtc, etc.
  ```

## Time Selection Parameters

### `-b TIME`
**Start time**

- **Purpose**: Specify start time for analysis
- **Unit**: picoseconds (ps)
- **Usage examples**:
  ```bash
  # Analyze from 10 ns onwards
  gmx rms -s md.tpr -f md.xtc -b 10000 -o rmsd.xvg
  ```

### `-e TIME`
**End time**

- **Purpose**: Specify end time for analysis
- **Unit**: picoseconds (ps)
- **Usage examples**:
  ```bash
  # Analyze first 50 ns
  gmx rms -s md.tpr -f md.xtc -e 50000 -o rmsd.xvg
  ```

### `-dt TIME`
**Time step**

- **Purpose**: Specify time step for output (frame selection)
- **Unit**: picoseconds (ps)
- **Usage examples**:
  ```bash
  # Output every 100 ps
  gmx trjconv -s md.tpr -f md.xtc -o reduced.xtc -dt 100
  ```

## Trajectory Processing Parameters

### `-pbc TYPE`
**Periodic boundary condition handling**

- **Purpose**: Handle periodic boundary conditions in trajectories
- **Options**:
  - `none`: No PBC correction
  - `mol`: Keep molecules whole
  - `atom`: Keep atoms whole (default)
  - `com`: Center on center of mass
  - `nojump`: Prevent atoms from jumping across PBC
  - `cluster`: Cluster molecules together
- **Usage examples**:
  ```bash
  # Keep molecules intact
  echo "Protein" | gmx trjconv -s md.tpr -f md.xtc -o mol.xtc -pbc mol

  # Prevent jumps
  echo "Protein" | gmx trjconv -s md.tpr -f md.xtc -o nojump.xtc -pbc nojump

  # Center on COM
  echo "Protein" | gmx trjconv -s md.tpr -f md.xtc -o com.xtc -pbc com
  ```

### `-center TYPE`
**Center coordinates**

- **Purpose**: Center simulation box on specified group
- **Usage examples**:
  ```bash
  # Center on protein
  echo "Protein\nProtein" | gmx trjconv -s md.tpr -f md.xtc -o centered.xtc -center

  # Center on center atom
  echo "center\nProtein" | gmx trjconv -s md.tpr -f md.xtc -o centered.xtc -center
  ```

### `-fit TYPE`
**Fit trajectory to reference**

- **Purpose**: Fit trajectory to reference structure (remove translation/rotation)
- **Options**:
  - `none`: No fitting
  - `xyz`: Fit on x, y, z coordinates
  - `rotxy`: Fit on rotation in x-y plane
  - `xy`: Fit on x-y plane only
  - `progressive`: Progressive fitting
  - `transxy`: Translation in x-y plane only
  - `rot+trans`: Rotation + translation (most common)
- **Usage examples**:
  ```bash
  # Remove translation and rotation
  echo "Backbone\nProtein" | gmx trjconv -s md.tpr -f md.xtc -o fit.xtc -fit rot+trans

  # Fit on backbone
  echo "Backbone\nProtein" | gmx trjconv -s md.tpr -f md.xtc -o fit.xtc -fit rot+trans
  ```

### `-ur TYPE`
**Unit cell representation**

- **Purpose**: Control how unit cell is represented
- **Options**:
  - `compact`: Compact representation (default)
  - `triclinic`: Full triclinic box
  - `rect`: Rectangular box
- **Usage examples**:
  ```bash
  echo "Protein" | gmx trjconv -s md.tpr -f md.xtc -o rect.xtc -ur rect
  ```

## Performance Parameters

### `-nt NUMBER`
**Number of threads**

- **Purpose**: Specify number of CPU threads
- **Usage examples**:
  ```bash
  # Use 4 threads
  gmx mdrun -s md.tpr -deffnm md -nt 4

  # Use 8 threads
  gmx mdrun -s md.tpr -deffnm md -nt 8
  ```

### `-ntmpi NUMBER`
**Number of MPI threads**

- **Purpose**: Specify number of MPI ranks
- **Usage examples**:
  ```bash
  # Use 4 MPI ranks
  mpirun -np 4 gmx_mpi mdrun -s md.tpr -deffnm md

  # Use 8 MPI ranks with 2 OpenMP threads each
  mpirun -np 8 gmx_mpi mdrun -ntomp 2 -s md.tpr -deffnm md
  ```

### `-ntomp NUMBER`
**OpenMP threads per MPI rank**

- **Purpose**: Specify number of OpenMP threads per MPI rank
- **Usage examples**:
  ```bash
  # 2 MPI ranks, 4 OpenMP threads each
  mpirun -np 2 gmx_mpi mdrun -ntomp 4 -s md.tpr -deffnm md

  # 4 MPI ranks, 2 OpenMP threads each
  mpirun -np 4 gmx_mpi mdrun -ntomp 2 -s md.tpr -deffnm md
  ```

### `-nb TYPE`
**Neighbor searching type**

- **Purpose**: Specify hardware for neighbor searching
- **Options**:
  - `cpu`: Use CPU
  - `gpu`: Use GPU (default when GPU available)
- **Usage examples**:
  ```bash
  # Use GPU for neighbor searching
  gmx mdrun -s md.tpr -deffnm md -ntomp 4 -nb gpu

  # Use CPU for neighbor searching
  gmx mdrun -s md.tpr -deffnm md -ntomp 4 -nb cpu
  ```

### `-update TYPE`
**Update method**

- **Purpose**: Specify hardware for constraint updates
- **Options**:
  - `cpu`: Use CPU (default)
  - `gpu`: Use GPU
- **Usage examples**:
  ```bash
  # Use GPU for updates
  gmx mdrun -s md.tpr -deffnm md -ntomp 4 -nb gpu -update gpu
  ```

### `-pme TYPE`
**PME method**

- **Purpose**: Specify hardware for PME electrostatics
- **Options**:
  - `cpu`: Use CPU
  - `gpu`: Use GPU
- **Usage examples**:
  ```bash
  # Use GPU for PME
  gmx mdrun -s md.tpr -deffnm md -ntomp 4 -nb gpu -pme gpu
  ```

### `-bonded TYPE`
**Bonded interactions**

- **Purpose**: Specify hardware for bonded interactions
- **Options**:
  - `cpu`: Use CPU (default)
  - `gpu`: Use GPU
- **Usage examples**:
  ```bash
  # Use GPU for bonded interactions
  gmx mdrun -s md.tpr -deffnm md -ntomp 4 -bonded gpu
  ```

### `-tunepme`
**Tune PME parameters**

- **Purpose**: Automatically tune PME parameters for optimal performance
- **Usage examples**:
  ```bash
  # Tune PME parameters
  gmx mdrun -s md.tpr -deffnm md -ntomp 4 -tunepme
  ```

### `-gpu_id ID`
**GPU device ID**

- **Purpose**: Specify which GPU to use
- **Usage examples**:
  ```bash
  # Use GPU 0
  gmx mdrun -s md.tpr -deffnm md -ntomp 4 -gpu_id 0

  # Use GPU 1
  gmx mdrun -s md.tpr -deffnm md -ntomp 4 -gpu_id 1
  ```

## Analysis Parameters

### `-type TYPE`
**Calculation type**

- **Purpose**: Specify type of calculation
- **Usage examples**:
  ```bash
  # Different types depending on command
  gmx rms -s md.tpr -f md.xtc -o rmsd.xvg -type fit
  ```

### `-ng NUMBER`
**Number of groups**

- **Purpose**: Specify number of groups for analysis
- **Usage examples**:
  ```bash
  # Analyze hydrogen bonds with 2 groups
  echo "Protein\nProtein" | gmx hbond -s md.tpr -f md.xtc -num hb.xvg
  ```

### `-frame FRAME`
**Specific frame number**

- **Purpose**: Specify specific frame for analysis
- **Usage examples**:
  ```bash
  # Extract frame 1000
  gmx trjconv -s md.tpr -f md.xtc -o frame1000.pdb -dump 1000
  ```

## Output Control Parameters

### `-nice LEVEL`
**Nice level for CPU usage**

- **Purpose**: Set process priority (nice level)
- **Range**: -20 (highest priority) to 19 (lowest priority)
- **Usage examples**:
  ```bash
  # Run with low priority
  gmx mdrun -s md.tpr -deffnm md -nice 10

  # Run with high priority
  gmx mdrun -s md.tpr -deffnm md -nice -5
  ```

### `-mpi`
**Enable MPI**

- **Purpose**: Enable MPI parallelization
- **Usage examples**:
  ```bash
  # Run with MPI
  mpirun -np 4 gmx_mpi mdrun -s md.tpr -deffnm md
  ```

## Visualization Parameters

### `-xvg TYPE`
**Plot format**

- **Purpose**: Specify output plot format for XVG files
- **Options**:
  - `none`: No plot format
  - `xmgrace`: Grace/Xmgr format (default)
  - `xmgr`: Xmgr format
  - `matplotlib`: Matplotlib format
- **Usage examples**:
  ```bash
  # Output in Grace format
  gmx rms -s md.tpr -f md.xtc -o rmsd.xvg -xvg xmgrace

  # Output in Matplotlib format
  gmx rms -s md.tpr -f md.xtc -o rmsd.xvg -xvg matplotlib
  ```

## Specialized Parameters

### `-maxwarn NUMBER`
**Maximum number of warnings to ignore**

- **Purpose**: Ignore warnings in grompp (use with caution)
- **Usage examples**:
  ```bash
  # Ignore up to 3 warnings
  gmx grompp -f mdp/md.mdp -c npt.gro -p topol.top -o md.tpr -maxwarn 3
  ```

### `-cpt NUMBER`
**Checkpoint interval**

- **Purpose**: Set checkpoint interval in minutes
- **Usage examples**:
  ```bash
  # Write checkpoint every 15 minutes
  gmx mdrun -s md.tpr -deffnm md -cpt 15

  # Write checkpoint every 30 minutes
  gmx mdrun -s md.tpr -deffnm md -cpt 30
  ```

### `-nsteps NUMBER`
**Override number of steps**

- **Purpose**: Override number of steps from .mdp file
- **Usage examples**:
  ```bash
  # Run 100000 steps instead of .mdp setting
  gmx mdrun -s md.tpr -deffnm md -nsteps 100000
  ```

### `-replex NUMBER`
**Number of replica exchange steps**

- **Purpose**: Set number of steps between replica exchange attempts
- **Usage examples**:
  ```bash
  # Exchange every 1000 steps
  gmx mdrun -s md.tpr -deffnm md -replex 1000
  ```

## Restart Parameters

### `-t INPUT`
**Input trajectory for restart**

- **Purpose**: Specify input trajectory for restarting simulation
- **Usage examples**:
  ```bash
  # Restart from checkpoint
  gmx grompp -f mdp/md.mdp -c md.cpt -t md.cpt -p topol.top -o md.tpr
  ```

### `-e INPUT`
**Input energy file for restart**

- **Purpose**: Specify input energy file for restarting simulation
- **Usage examples**:
  ```bash
  # Restart from checkpoint
  gmx grompp -f mdp/md.mdp -c md.cpt -t md.cpt -e md.edr -p topol.top -o md.tpr
  ```

## Index and Selection Parameters

### `-select STRING`
**Selection string**

- **Purpose**: Select atoms using selection language
- **Usage examples**:
  ```bash
  # Select residues 1-100
  gmx select -s md.tpr -select "resid 1 to 100" -on res1-100.ndx

  # Select backbone atoms
  gmx select -s md.tpr -select "backbone" -on backbone.ndx

  # Select protein and exclude hydrogens
  gmx select -s md.tpr -select "protein and not hydrogen" -on protein_noH.ndx
  ```

## Reference Structure Parameters

### `-r INPUT`
**Reference structure for fitting**

- **Purpose**: Specify reference structure for fitting
- **Usage examples**:
  ```bash
  # Use average structure as reference
  gmx grompp -f mdp/md.mdp -c npt.gro -r average.pdb -p topol.top -o md.tpr
  ```

## Output Processing Parameters

### `-po OUTPUT`
**Processed topology**

- **Purpose**: Output processed topology file
- **Usage examples**:
  ```bash
  # Output processed topology
  gmx grompp -f mdp/md.mdp -c npt.gro -p topol.top -o md.tpr -po processed.top
  ```

## Energy Parameters

### `-dp`
**Remove drift from energy**

- **Purpose**: Remove linear drift from energy data
- **Usage examples**:
  ```bash
  # Remove drift from potential energy
  echo "Potential" | gmx energy -f md.edr -o potential.xvg -dp
  ```

### `-ee`
**Estimate error**

- **Purpose**: Estimate error in energy (requires multiple simulations)
- **Usage examples**:
  ```bash
  # Estimate error
  echo "Potential" | gmx energy -f md.edr -o potential.xvg -ee
  ```

## Parameter Groups

### Input File Group
- `-f`, `-s`, `-n`, `-c`, `-p`, `-t`, `-e`, `-r`

### Output File Group
- `-o`, `-deffnm`, `-po`

### Time Selection Group
- `-b`, `-e`, `-dt`

### Trajectory Processing Group
- `-pbc`, `-center`, `-fit`, `-ur`

### Performance Group
- `-nt`, `-ntmpi`, `-ntomp`, `-nb`, `-update`, `-pme`, `-bonded`, `-gpu_id`

### Control Group
- `-nice`, `-maxwarn`, `-cpt`, `-nsteps`

### Visualization Group
- `-xvg`

## Tips and Best Practices

1. **Always use absolute paths** when running on clusters or complex directory structures
2. **Check file compatibility** before running (use `gmx check`)
3. **Use appropriate time units** (picoseconds for -b, -e, -dt)
4. **Monitor performance** with different thread/GPU settings
5. **Save checkpoint files** regularly for long simulations
6. **Document all parameters** for reproducibility
7. **Test on small systems** before large simulations
8. **Use appropriate file formats** (.xtc for storage, .trr for analysis)