# GROMACS Command Categories

This document provides a comprehensive list of GROMACS commands organized by functional category.

## 1. Topology and Structure

### `pdb2gmx`
Generate GROMACS topology from PDB structure.

**Purpose**: Convert PDB files to GROMACS topology and coordinate files.

**Key parameters**:
- `-f INPUT`: Input PDB file
- `-o OUTPUT`: Output GROMACS coordinate file
- `-p OUTPUT`: Output topology file
- `-i OUTPUT`: Output position restraint file
- `-ff FORCEFIELD`: Force field selection
- `-water WATERMODEL`: Water model selection
- `-ignh`: Ignore hydrogen atoms in PDB
- `-merge`: Merge non-interactive atoms
- `-chainsepitp`: Generate separate chain identifiers

**Common workflows**:
- Basic: `pdb2gmx -f protein.pdb -o protein.gro -p topol.top`
- With water: `pdb2gmx -f protein.pdb -o protein.gro -p topol.top -water tip3p`

### `editconf`
Edit structure files (box definition, dimensions).

**Purpose**: Define simulation box, modify box dimensions, center molecules.

**Key parameters**:
- `-f INPUT`: Input structure file
- `-o OUTPUT`: Output structure file
- `-bt TYPE`: Box type (cubic, dodecahedron, octahedron, triclinic)
- `-d DISTANCE`: Distance between molecule and box edges
- `-c CENTER`: Center coordinates
- `-center`: Center molecule in box
- `-box SIZE`: Box dimensions (x y z)
- `-angles ANGLES`: Box angles (for triclinic)

**Common workflows**:
- Define cubic box: `editconf -f protein.gro -o protein_box.gro -bt cubic -d 1.0`
- Center molecule: `editconf -f protein.gro -o protein_centered.gro -center`

### `solvate`
Add solvent molecules to simulation box.

**Purpose**: Solvate simulation system with water molecules.

**Key parameters**:
- `-cp INPUT`: Input structure (protein in box)
- `-cs INPUT`: Solvent structure (e.g., spc216.gro)
- `-o OUTPUT`: Output solvated structure
- `-p TOPOLOGY`: Output topology file

**Common workflows**:
- Basic solvation: `solvate -cp protein_box.gro -cs spc216.gro -o protein_solv.gro -p topol.top`

### `insert-molecules`
Insert molecules into simulation box.

**Purpose**: Insert multiple copies of molecules into simulation box.

**Key parameters**:
- `-f INPUT`: Input structure
- `-ci INPUT`: Insert structure
| `-nmol NUMBER`: Number of molecules to insert
| `-o OUTPUT`: Output structure
| `-radius DISTANCE`: Minimum distance between molecules

**Common workflows**:
- Insert ions: `insert-molecules -f protein_solv.gro -ci ions.gro -nmol 10 -o protein_ions.gro`

### `genrestr`
Generate position restraint files.

**Purpose**: Create position restraint files for equilibration.

**Key parameters**:
- `-f INPUT`: Input structure
- `-o OUTPUT`: Output position restraint file
- `-i INPUT`: Index file with restrained atoms
| `-fc NUMBER`: Force constant

**Common workflows**:
- Restrain protein: `genrestr -f protein.gro -o posre.itp -f index.ndx -fc 1000`

## 2. Simulation Setup

### `grompp`
Generate simulation run input file (.tpr) from topology and parameters.

**Purpose**: Compile topology, parameters, and simulation conditions into .tpr file.

**Key parameters**:
- `-f INPUT`: Input parameter file (.mdp)
- `-c INPUT`: Input structure file
- `-p INPUT`: Input topology file
| `-o OUTPUT`: Output .tpr file
| `-t INPUT`: Input trajectory (for restart)
| `-e INPUT`: Energy file (for restart)
| `-r INPUT`: Reference structure for fitting
| `-n INPUT`: Index file
| `-po OUTPUT`: Processed topology
| `-maxwarn NUMBER`: Maximum number of warnings to ignore

**Common workflows**:
- Basic: `grompp -f mdp/minim.mdp -c ions.gro -p topol.top -o em.tpr`
- With index: `grompp -f mdp/nvt.mdp -c npt.gro -p topol.top -n index.ndx -o nvt.tpr`
- Continue: `grompp -f mdp/md.mdp -c md.cpt -t md.cpt -p topol.top -o md.tpr -n index.ndx`

### `mdrun`
Run molecular dynamics simulation.

**Purpose**: Execute molecular dynamics simulation.

**Key parameters**:
- `-s INPUT`: Input .tpr file
| `-deffnm BASENAME`: Default filename base
| `-o OUTPUT`: Output trajectory
| `-x OUTPUT`: Compressed trajectory
| `-c OUTPUT`: Final coordinates
| `-e OUTPUT`: Output energy file
| `-g OUTPUT`: Output log file
| `-cpt NUMBER`: Checkpoint interval (minutes)
| `-nt NUMBER`: Number of threads
| `-ntmpi NUMBER`: Number of MPI threads
| `-nb TYPE`: Neighbor searching (cpu, gpu)
| `-update TYPE`: Update method (cpu, gpu)
| `-pme TYPE`: PME method (cpu, gpu)
| `-bonded TYPE`: Bonded interactions (cpu, gpu)
| `-ntomp NUMBER`: OpenMP threads per MPI rank
| `-tunepme`: Tune PME parameters
| `-resethway`: Reset hardware interactions
| `-gpu_id ID`: GPU device ID
| `-nsteps NUMBER**: Override number of steps
| `-replex NUMBER**: Number of replica exchange steps

**Common workflows**:
- Basic: `mdrun -s em.tpr -deffnm em`
- Multi-threaded: `mdrun -s md.tpr -deffnm md -nt 4`
- GPU: `mdrun -s md.tpr -deffnm md -ntomp 4 -nb gpu`
- MPI: `mpirun -np 4 gmx_mpi mdrun -s md.tpr -deffnm md`

### `mdrun-mpi`
Run simulation with MPI parallelization.

**Purpose**: MPI-optimized version of mdrun for multi-node execution.

**Key parameters**: Same as `mdrun`, with additional MPI-specific options.

**Common workflows**:
- Basic MPI: `mpirun -np 4 gmx_mpi mdrun -s md.tpr -deffnm md`
- Hybrid MPI+OpenMP: `mpirun -np 2 gmx_mpi mdrun -ntomp 2 -s md.tpr -deffnm md`

## 3. Energy Analysis

### `energy`
Extract energy components from energy file.

**Purpose**: Extract and analyze energy terms from .edr file.

**Key parameters**:
- `-f INPUT`: Input energy file
| `-o OUTPUT`: Output XVG file
| `-dp`: Remove drift from energy
| `-ee`: Estimate error (requires multiple energies)
| `-b TIME`: Start time
| `-e TIME`: End time

**Common workflows**:
- Extract potential energy: `echo "Potential" | gmx energy -f md.edr -o potential.xvg`
- Extract temperature: `echo "Temperature" | gmx energy -f npt.edr -o temp.xvg`
- Extract pressure: `echo "Pressure" | gmx energy -f npt.edr -o press.xvg`

### `eneconv`
Convert energy file formats.

**Purpose**: Convert between different energy file formats.

**Key parameters**:
- `-f INPUT`: Input energy file
| `-o OUTPUT`: Output energy file

**Common workflows**:
- Basic: `eneconv -f md.edr -o md_new.edr`

### `bar`
Free energy calculations using Bennett acceptance ratio.

**Purpose**: Calculate free energy differences from thermodynamic integration.

**Key parameters**:
| `-f INPUT`: Input energy files
| `-o OUTPUT`: Output XVG file
| `-i INPUT`: Input lambda values

**Common workflows**:
- Basic BAR: `gmx bar -f lambda0.edr lambda1.edr -o dg.xvg`

## 4. Trajectory Analysis

### `rms`
Calculate root mean square deviation (RMSD).

**Purpose**: Calculate RMSD of structures relative to reference.

**Key parameters**:
- `-s INPUT`: Input .tpr file
| `-f INPUT`: Input trajectory
| `-o OUTPUT`: Output XVG file
| `-n INPUT`: Index file
| `-type TYPE`: Calculation type (fit, nofit)
| `-pbc TYPE`: PBC handling

**Common workflows**:
- RMSD to initial structure: `echo "Backbone\nProtein" | gmx rms -s md.tpr -f md.xtc -o rmsd.xvg`
- RMSD to average: `echo "Backbone\nProtein" | gmx rms -s md.tpr -f md.xtc -o rmsd_avg.xvg -ref average.pdb`

### `rmsf`
Calculate root mean square fluctuation (RMSF).

**Purpose**: Calculate per-atom RMSF.

**Key parameters**:
- `-s INPUT`: Input .tpr file
| `-f INPUT`: Input trajectory
| `-o OUTPUT`: Output XVG file
| `-n INPUT`: Index file
| `-res`: Output per-residue

**Common workflows**:
- Per-atom RMSF: `echo "Backbone" | gmx rmsf -s md.tpr -f md.xtc -o rmsf.xvg`
- Per-residue RMSF: `echo "Backbone" | gmx rmsf -s md.tpr -f md.xtc -o rmsf_res.xvg -res`

### `gyrate`
Calculate radius of gyration.

**Purpose**: Calculate radius of gyration and components.

**Key parameters**:
- `-s INPUT`: Input .tpr file
| `-f INPUT`: Input trajectory
| `-o OUTPUT`: Output XVG file
| `-n INPUT`: Index file

**Common workflows**:
- Basic: `echo "Protein" | gmx gyrate -s md.tpr -f md.xtc -o rg.xvg`

### `hbond`
Analyze hydrogen bonds.

**Purpose**: Calculate and analyze hydrogen bonds.

**Key parameters**:
- `-s INPUT`: Input .tpr file
| `-f INPUT`: Input trajectory
| `-num OUTPUT`: Output hydrogen bond count
| `-dist OUTPUT`: Output hydrogen bond distances
| `-ang OUTPUT`: Output hydrogen bond angles
| `-hbn OUTPUT`: Output hydrogen bond list
| `-hbm OUTPUT`: Output hydrogen bond matrix
| `-n INPUT`: Index file

**Common workflows**:
- Count hydrogen bonds: `echo "Protein\nProtein" | gmx hbond -s md.tpr -f md.xtc -num hb.xvg`
- Distance analysis: `echo "Protein\nProtein" | gmx hbond -s md.tpr -f md.xtc -dist hb_dist.xvg`

### `distance`
Calculate distances between atom groups.

**Purpose**: Calculate distance between two atom groups over time.

**Key parameters**:
- `-s INPUT`: Input .tpr file
| `-f INPUT`: Input trajectory
| `-n INPUT`: Index file
| `-oall OUTPUT`: Output all distances
| `-oallxyz OUTPUT`: Output x,y,z components
| `-select STRING`: Selection string

**Common workflows**:
- Basic: `echo "Group1\nGroup2" | gmx distance -s md.tpr -f md.xtc -n index.ndx -oall dist.xvg`

### `angle`
Calculate angles between atom groups.

**Purpose**: Calculate angle between three atom groups.

**Key parameters**:
- `-s INPUT`: Input .tpr file
| `-f INPUT`: Input trajectory
| `-n INPUT`: Index file
| `-oall OUTPUT`: Output all angles

**Common workflows**:
- Basic: `echo "Group1\nGroup2\nGroup3" | gmx angle -s md.tpr -f md.xtc -n index.ndx -oall angle.xvg`

### `dihedral`
Calculate dihedral angles.

**Purpose**: Calculate dihedral angles between four atom groups.

**Key parameters**:
- `-s INPUT`: Input .tpr file
| `-f INPUT`: Input trajectory
| `-n INPUT`: Index file
| `-oall OUTPUT`: Output all dihedrals

**Common workflows**:
- Basic: `echo "Group1\nGroup2\nGroup3\nGroup4" | gmx dihedral -s md.tpr -f md.xtc -n index.ndx -oall dih.xvg`

### `cluster`
Cluster analysis of conformations.

**Purpose**: Cluster trajectory frames based on RMSD.

**Key parameters**:
- `-s INPUT`: Input .tpr file
| `-f INPUT`: Input trajectory
| `-dm INPUT`: Distance matrix
| `-method METHOD`: Clustering method (linkage, gromos, etc.)
| `-cutoff DISTANCE`: Cutoff distance
| `-cl OUTPUT`: Output cluster file
| `-o OUTPUT`: Output XVG file

**Common workflows**:
- Basic: `echo "Backbone" | gmx cluster -s md.tpr -f md.xtc -method gromos -cutoff 0.2 -cl clusters.xpm -o cluster_size.xvg`

### `mindist`
Calculate minimum distance between groups.

**Purpose**: Calculate minimum distance between atom groups.

**Key parameters**:
- `-s INPUT`: Input .tpr file
| `-f INPUT`: Input trajectory
| `-n INPUT`: Index file
| `-od OUTPUT`: Output minimum distance
| `-pi OUTPUT`: Output periodic images

**Common workflows**:
- Basic: `echo "Group1\nGroup2" | gmx mindist -s md.tpr -f md.xtc -n index.ndx -od mindist.xvg`

### `sasa`
Calculate solvent accessible surface area.

**Purpose**: Calculate SASA for selected atoms.

**Key parameters**:
- `-s INPUT`: Input .tpr file
| `-f INPUT`: Input trajectory
| `-n INPUT`: Index file
| `-o OUTPUT`: Output XVG file
| `-surface OUTPUT`: Output SASA per atom
| `-probe RADIUS`: Probe radius

**Common workflows**:
- Basic: `echo "Protein" | gmx sasa -s md.tpr -f md.xtc -o sasa.xvg`

### `principal`
Principal component analysis of molecules.

**Purpose**: Calculate principal axes and moments of inertia.

**Key parameters**:
- `-s INPUT`: Input .tpr file
| `-f INPUT`: Input trajectory
| `-n INPUT`: Index file
| `-a OUTPUT`: Output XVG file with axes
| `-i OUTPUT`: Output moments of inertia

**Common workflows**:
- Basic: `echo "Protein" | gmx principal -s md.tpr -f md.xtc -a axes.xvg -i inertia.xvg`

### `do_dssp`
Secondary structure analysis (DSSP).

**Purpose**: Assign secondary structure using DSSP algorithm.

**Key parameters**:
- `-s INPUT`: Input .tpr file
| `-f INPUT`: Input trajectory
| `-sc OUTPUT`: Output secondary structure per residue
| `-ss OUTPUT`: Output secondary structure fractions
| `-a OUTPUT`: Output secondary structure area

**Common workflows**:
- Basic: `echo "Protein" | gmx do_dssp -s md.tpr -f md.xtc -sc ss.xpm -ss ss.xvg`

## 5. Structural Analysis

### `covar`
Covariance matrix analysis (for PCA).

**Purpose**: Calculate covariance matrix of atomic fluctuations.

**Key parameters**:
- `-s INPUT`: Input .tpr file
| `-f INPUT`: Input trajectory
| `-n INPUT`: Index file
| `-o OUTPUT`: Output eigenvalues
| `-v OUTPUT`: Output eigenvectors
| `-xpma OUTPUT`: Output covariance matrix
| `-ascii OUTPUT`: Output ASCII covariance matrix

**Common workflows**:
- Basic: `echo "C-alpha\nC-alpha" | gmx covar -s md.tpr -f md.xtc -o eigen.xvg -v eigenvec.trr -xpma covar.xpm`

### `anaeig`
Analyze eigenvectors from covariance analysis.

**Purpose**: Project trajectory onto eigenvectors and analyze PCA.

**Key parameters**:
- `-s INPUT`: Input .tpr file
| `-f INPUT`: Input trajectory
| `-v INPUT`: Input eigenvectors
| `-first NUMBER`: First eigenvector to analyze
| `-last NUMBER`: Last eigenvector to analyze
| `-proj OUTPUT`: Output projection
| `-2d OUTPUT`: Output 2D projection
| `-extreme OUTPUT`: Output extreme structures
| `-overlap OUTPUT`: Output overlap matrix

**Common workflows**:
- Project on PC1: `echo "C-alpha\nC-alpha" | gmx anaeig -s md.tpr -f md.xtc -v eigenvec.trr -first 1 -last 1 -proj pc1.xvg`
- 2D projection: `echo "C-alpha\nC-alpha" | gmx anaeig -s md.tpr -f md.xtc -v eigenvec.trr -first 1 -last 2 -2d 2dproj.xvg`

### `mdmat`
Calculate distance matrix between residues.

**Purpose**: Calculate average distance matrix between residues.

**Key parameters**:
- `-s INPUT`: Input .tpr file
| `-f INPUT`: Input trajectory
| `-n INPUT`: Index file
| `-mean OUTPUT`: Output mean distance matrix

**Common workflows**:
- Basic: `echo "Protein" | gmx mdmat -s md.tpr -f md.xtc -mean rdcm.xpm`

### `sham`
Generate free energy landscapes.

**Purpose**: Generate free energy landscapes from reaction coordinates.

**Key parameters**:
- `-f INPUT`: Input XVG file with reaction coordinates
| `-tsham TEMP`: Temperature
| `-ls OUTPUT`: Output Gibbs free energy
| `-g OUTPUT`: Output log file
| `-nlevels NUMBER`: Number of levels

**Common workflows**:
- Basic: `gmx sham -f sham.xvg -tsham 310 -nlevels 100 -ls gibbs.xpm`

### `order`
Calculate order parameters.

**Purpose**: Calculate lipid order parameters.

**Key parameters**:
- `-s INPUT`: Input .tpr file
| `-f INPUT`: Input trajectory
| `-n INPUT`: Index file
| `-d OUTPUT`: Output deuterium order parameters

**Common workflows**:
- Basic: `echo "Lipid" | gmx order -s md.tpr -f md.xtc -o order.xvg`

### `rotacf`
Calculate rotational autocorrelation functions.

**Purpose**: Calculate rotational correlation times.

**Key parameters**:
- `-s INPUT`: Input .tpr file
| `-f INPUT`: Input trajectory
| `-n INPUT`: Index file
| `-o OUTPUT`: Output correlation functions

**Common workflows**:
- Basic: `echo "Protein" | gmx rotacf -s md.tpr -f md.xtc -o rotacf.xvg`

### `dielectric`
Calculate dielectric properties.

**Purpose**: Calculate dielectric constant and properties.

**Key parameters**:
- `-f INPUT`: Input dipole moment trajectory
| `-d OUTPUT`: Output dielectric constant
| `-o OUTPUT`: Output correlation functions

**Common workflows**:
- Basic: `gmx dielectric -f dipole.xvg -d dielectric.xvg`

## 6. Trajectory Processing

### `trjconv`
Convert trajectory formats, apply PBC correction, fit structures.

**Purpose**: Process trajectory files (format conversion, PBC correction, fitting).

**Key parameters**:
- `-s INPUT`: Input .tpr file
| `-f INPUT`: Input trajectory
| `-o OUTPUT`: Output trajectory
| `-n INPUT`: Index file
| `-pbc TYPE`: PBC handling
| `-center TYPE`: Centering
| `-fit TYPE`: Fitting
| `-ur TYPE`: Unit cell representation
| `-b TIME`: Start time
| `-e TIME`: End time
| `-dt TIME`: Time step
| `-pbc mol`: Keep molecules intact

**Common workflows**:
- PBC correction: `echo "Protein\nProtein" | gmx trjconv -s md.tpr -f md.xtc -o fit.xtc -pbc mol -center`
- Format conversion: `echo "Protein" | gmx trjconv -s md.tpr -f md.xtc -o md.pdb`
- Remove translation/rotation: `echo "Backbone\nProtein" | gmx trjconv -s md.tpr -f md.xtc -o fit.xtc -fit rot+trans`

### `trjcat`
Concatenate multiple trajectory files.

**Purpose**: Combine multiple trajectory files into one.

**Key parameters**:
- `-f INPUT`: Input trajectories
| `-o OUTPUT`: Output trajectory
| `-s INPUT`: Input .tpr file
| `-cat`: Concatenate without time correction

**Common workflows**:
- Basic: `gmx trjcat -f part1.xtc part2.xtc -o combined.xtc -s md.tpr`

### `trjorder`
Reorder atoms in trajectory.

**Purpose**: Reorder atoms in trajectory file.

**Key parameters**:
- `-s INPUT`: Input .tpr file
| `-f INPUT`: Input trajectory
| `-o OUTPUT`: Output trajectory
| `-n INPUT`: Index file

**Common workflows**:
- Basic: `echo "Protein" | gmx trjorder -s md.tpr -f md.xtc -o ordered.xtc`

### `dump`
Dump trajectory to text format.

**Purpose**: Convert trajectory to ASCII text format.

**Key parameters**:
- `-s INPUT`: Input .tpr file
| `-f INPUT`: Input trajectory
| `-o OUTPUT`: Output text file

**Common workflows**:
- Basic: `gmx dump -s md.tpr -f md.xtc -o md.txt`

### `eneconv`
Convert energy file formats.

**Purpose**: Convert between different energy file formats.

**Key parameters**:
- `-f INPUT`: Input energy file
| `-o OUTPUT`: Output energy file

**Common workflows**:
- Basic: `eneconv -f md.edr -o md_new.edr`

## 7. Index and Selection

### `make_ndx`
Create or edit index files.

**Purpose**: Create and modify index files with atom groups.

**Key parameters**:
- `-f INPUT`: Input structure or .tpr file
| `-o OUTPUT`: Output index file
| `-n INPUT`: Input index file
| `-t`: Use selection language

**Common workflows**:
- Create basic index: `gmx make_ndx -f md.tpr -o index.ndx`
- Add protein-ligand group: `echo "Protein | Ligand\nq" | gmx make_ndx -f md.tpr -o index.ndx`
- Add residue range: `echo "r 1-100\nname 10 Residues1-100\nq" | gmx make_ndx -f md.tpr -o index.ndx`

### `select`
Select atoms using selection language.

**Purpose**: Select atoms using powerful selection syntax.

**Key parameters**:
- `-s INPUT`: Input .tpr file
| `-select STRING`: Selection string
| `-on OUTPUT`: Output index file
| `-os OUTPUT`: Output structure
| `-om OUTPUT`: Output PDB with selections

**Common workflows**:
- Basic selection: `gmx select -s md.tpr -select "resid 1 to 100 and name CA" -on ca_res1-100.ndx`
- Output structure: `gmx select -s md.tpr -select "protein" -os protein.pdb`

### `genrestr`
Generate position restraints.

**Purpose**: Create position restraint files for equilibration.

**Key parameters**:
- `-f INPUT`: Input structure
| `-o OUTPUT`: Output position restraint file
| `-i INPUT`: Index file
| `-fc NUMBER`: Force constant

**Common workflows**:
- Restrain protein backbone: `echo "Backbone" | gmx genrestr -f protein.gro -o posre.itp -fc 1000`

### `genion`
Replace water molecules with ions.

**Purpose**: Replace solvent molecules with ions to neutralize system.

**Key parameters**:
- `-s INPUT`: Input .tpr file
| `-p INPUT`: Input topology
| `-o OUTPUT`: Output structure
| `-pname NAME`: Positive ion name
| `-nname NAME`: Negative ion name
| `-pname NAME`: Positive ion name
| `-nname NAME`: Negative ion name
| `-conc CONC`: Ion concentration
| `-neutral`: Neutralize system
| `-p NUMBER`: Number of positive ions
| `-n NUMBER`: Number of negative ions

**Common workflows**:
- Neutralize: `echo "SOL" | gmx genion -s ions.tpr -p topol.top -o ions.gro -pname NA -nname CL -neutral`
- Add specific ions: `echo "SOL" | gmx genion -s ions.tpr -p topol.top -o ions.gro -pname NA -nname CL -p 10 -n 10`

## 8. Tools and Utilities

### `xpm2ps`
Convert XPM matrix files to PostScript.

**Purpose**: Convert XPM matrix files to PostScript for visualization.

**Key parameters**:
- `-f INPUT`: Input XPM file
| `-o OUTPUT`: Output PostScript file
| `-rainbow`: Use rainbow colors
| `-title STRING`: Plot title

**Common workflows**:
- Basic: `gmx xpm2ps -f matrix.xpm -o matrix.eps`

### `x2top`
Generate topology from coordinates.

**Purpose**: Generate topology from coordinates using force field.

**Key parameters**:
- `-f INPUT`: Input coordinates
| `-o OUTPUT`: Output topology
| `-ff FORCEFIELD`: Force field

**Common workflows**:
- Basic: `gmx x2top -f molecule.gro -o molecule.top -ff oplsaa`

### `check`
Check simulation files for errors.

**Purpose**: Validate simulation files and check for errors.

**Key parameters**:
- `-s INPUT`: Input .tpr file
| `-f INPUT`: Input trajectory
| `-e INPUT`: Input energy file

**Common workflows**:
- Basic: `gmx check -s md.tpr -f md.xtc -e md.edr`

### `wham`
Weighted histogram analysis method.

**Purpose**: Calculate free energy profiles from umbrella sampling.

**Key parameters**:
- `-it INPUT`: Input tpr files
| `-if INPUT`: Input pullf files
| `-o OUTPUT`: Output profile

**Common workflows**:
- Basic: `gmx wham -it tprs.dat -if pullf-files.dat -o profile.xvg`

### `tune_pme`
Tune PME parameters for performance.

**Purpose**: Optimize PME parameters for specific hardware.

**Key parameters**:
- `-s INPUT`: Input .tpr file
| `-np NUMBER`: Number of processors

**Common workflows**:
- Basic: `gmx tune_pme -s md.tpr -np 4`

## Command Summary by Category

- **Topology/Structure**: pdb2gmx, editconf, solvate, insert-molecules, genrestr
- **Simulation Setup**: grompp, mdrun, mdrun-mpi
- **Energy Analysis**: energy, eneconv, bar
- **Trajectory Analysis**: rms, rmsf, gyrate, hbond, distance, angle, dihedral, cluster, mindist, sasa, principal, do_dssp
- **Structural Analysis**: covar, anaeig, mdmat, sham, order, rotacf, dielectric
- **Trajectory Processing**: trjconv, trjcat, trjorder, dump, eneconv
- **Index/Selection**: make_ndx, select, genrestr, genion
- **Tools/Utilities**: xpm2ps, x2top, check, wham, tune_pme