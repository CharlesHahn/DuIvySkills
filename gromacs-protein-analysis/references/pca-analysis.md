# Principal Component Analysis (PCA) of Protein Dynamics

## Overview

Principal Component Analysis (PCA) identifies collective motions and dominant conformational changes in proteins by decomposing protein motion into principal components. It reduces dimensionality while preserving the most important dynamic information.

## When to Use PCA

- Identify collective motions and conformational transitions
- Extract major modes of protein flexibility
- Reduce dimensionality for further analysis
- Analyze protein domain movements
- Compare conformational changes between simulations
- Understand functional motions

## Prerequisites

- PBC-corrected trajectory (see [PBC Correction Guide](pbc-correction.md))
- Topology file (.tpr)
- Index file (.ndx) with appropriate atom groups
- Trajectory should be fitted to remove overall translation/rotation

## Workflow

### Step 1: Calculate Covariance Matrix

Compute covariance matrix and extract eigenvectors/eigenvalues:

```bash
echo -e "C-alpha\nC-alpha\n" | gmx covar -s md.tpr -f md.xtc -o eigenvalues.xvg -v eigenvectors.trr -xpma covar.xpm
```

- First input: Select reference group for alignment (C-alpha)
- Second input: Select group for covariance calculation (C-alpha)

**Output files**:
- `eigenvalues.xvg`: Eigenvalues from covariance matrix diagonalization
- `eigenvectors.trr`: Eigenvectors stored as trajectory
- `covar.xpm`: Covariance matrix (can be used for DCCM)

### Step 2: Analyze Eigenvalue Distribution

Examine eigenvalue distribution to understand collective motions:

```bash
# Plot eigenvalues
dit xvg_show -f eigenvalues.xvg -x "PC Number" -y "Eigenvalue" -t "Eigenvalue Spectrum"
```

Calculate contribution percentages:

```bash
# Read eigenvalues.xvg
# First column: PC number
# Second column: Eigenvalue

# Total variance = sum of all eigenvalues
# PC contribution = (individual eigenvalue / total variance) × 100%

# Example calculation for first 3 PCs:
# PC1% = eigenvalue1 / total × 100%
# PC2% = eigenvalue2 / total × 100%
# PC3% = eigenvalue3 / total × 100%
# Cumulative% = (PC1 + PC2 + PC3) / total × 100%
```

**Interpretation**:
- First few PCs often capture 50-80% of total motion
- Sharp drop in eigenvalues indicates dominant collective motions
- Uniform eigenvalues suggest no dominant collective motion

### Step 3: Project Trajectory onto Principal Components

Project trajectory onto first few principal components:

```bash
# Project onto PC1
echo -e "C-alpha\nC-alpha\n" | gmx anaeig -s md.tpr -f md.xtc -v eigenvectors.trr -first 1 -last 1 -proj pc1.xvg

# Project onto PC2
echo -e "C-alpha\nC-alpha\n" | gmx anaeig -s md.tpr -f md.xtc -v eigenvectors.trr -first 2 -last 2 -proj pc2.xvg

# Project onto PC3
echo -e "C-alpha\nC-alpha\n" | gmx anaeig -s md.tpr -f md.xtc -v eigenvectors.trr -first 3 -last 3 -proj pc3.xvg
```

**Output files**:
- `pc1.xvg`, `pc2.xvg`, `pc3.xvg`: Projections onto respective PCs

### Step 4: Visualize PC Projections

Generate plots of PC projections over time:

```bash
# PC1 projection
dit xvg_show -f pc1.xvg -x "Time (ps)" -y "PC1 Coordinate" -t "PC1 Projection"

# PC2 projection
dit xvg_show -f pc2.xvg -x "Time (ps)" -y "PC2 Coordinate" -t "PC2 Projection"

# PC1 vs PC2 scatter plot
dit xvg_show_scatter -f pc1.xvg pc2.xvg -c 1 1 -x "PC1" -y "PC2" -t "PC1 vs PC2"
```

### Step 5: Extract Extreme Conformations

Extract structures corresponding to extreme values along principal components:

```bash
# Extract extreme conformations along PC1
echo -e "C-alpha\nC-alpha\n" | gmx anaeig -s md.tpr -f md.xtc -v eigenvectors.trr -first 1 -last 1 -extreme pc1_extreme.pdb

# Extract average structure
echo -e "C-alpha\nC-alpha\n" | gmx anaeig -s md.tpr -f md.xtc -v eigenvectors.trr -first 1 -last 1 -average average.pdb
```

**Output files**:
- `pc1_extreme.pdb`: Structures with maximum deviation along PC1
- `average.pdb`: Average structure over the trajectory

### Step 6: 2D Projection (Optional)

Generate 2D projection of trajectory onto PC1-PC2 plane:

```bash
echo -e "C-alpha\nC-alpha\n" | gmx anaeig -s md.tpr -f md.xtc -v eigenvectors.trr -first 1 -last 2 -2d 2dproj.xvg
```

**Output file**:
- `2dproj.xvg`: 2D projection (PC1 vs PC2) over time

## Output Files

- **eigenvalues.xvg**: Eigenvalues from covariance matrix
- **eigenvectors.trr**: Eigenvectors stored as trajectory
- **covar.xpm**: Covariance matrix
- **pc1.xvg, pc2.xvg, pc3.xvg**: Projections onto principal components
- **2dproj.xvg**: 2D projection (PC1 vs PC2)
- **pc1_extreme.pdb**: Extreme conformations along PC1
- **average.pdb**: Average structure

## Interpretation Guidelines

### Eigenvalue Analysis

- **First PC**: Largest eigenvalue, captures most dominant motion
- **Cumulative variance**: Sum of first N eigenvalues / total variance
- **Variance fraction**: Individual PC contribution (eigenvalue / total)
- **Dimensionality**: Number of PCs needed to explain 80-90% of variance

### PC Projections

- **Amplitude**: Range of PC projection indicates magnitude of motion
- **Timescales**: Oscillation frequency relates to motion timescale
- **Transitions**: Jumps or shifts indicate conformational changes
- **Stationarity**: Constant mean indicates equilibrium sampling

### Extreme Conformations

- **Compare structures**: Examine differences between extreme conformations
- **Identify moving regions**: Locate atoms contributing most to PC motion
- **Functional relevance**: Relate motions to protein function

### Biological Interpretation

- **Domain motions**: Large-scale movements between domains
- **Loop flexibility**: Local motions in loops and termini
- **Allosteric pathways**: Correlated motions indicating allosteric communication
- **Functional motions**: Motions related to protein activity (e.g., hinge bending, channel opening)

## Common Issues and Solutions

### Issue: Eigenvalues show uniform distribution

**Possible causes**:
- Insufficient simulation time
- Excessive noise from high-frequency motions
- No dominant collective motion (protein may be rigid)

**Solutions**:
- Extend simulation time
- Apply time averaging or filtering
- Use larger atom groups (backbone instead of C-alpha)

### Issue: PC projections show no clear pattern

**Solutions**:
- Check trajectory fitting (remove translation/rotation)
- Verify atom selection consistency
- Consider using different atom groups

### Issue: Extreme conformations look unrealistic

**Solutions**:
- Verify simulation quality
- Check for numerical instabilities
- Consider using multiple PCs for extraction

### Issue: First few PCs don't capture significant variance

**Solutions**:
- May indicate protein is rigid or motion is evenly distributed
- Consider using more PCs for analysis
- Check simulation time and sampling

## Tips and Best Practices

- **Atom selection**: C-alpha is common for efficiency. Use backbone for detail.
- **Time selection**: Exclude equilibration. Use production phase only.
- **Trajectory fitting**: Always remove translation/rotation before PCA.
- **Number of PCs**: Analyze first 3-5 PCs, which typically capture most motion.
- **Validation**: Compare with experimental data if available.
- **Visualization**: Use multiple visualization methods for comprehensive analysis.

## Advanced Analysis

### Time-dependent PCA

Calculate PCA for different time windows:

```bash
# PCA for first half
echo -e "C-alpha\nC-alpha\n" | gmx covar -s md.tpr -f md.xtc -b 0 -e 50000 -o eigenvalues_early.xvg -v eigenvectors_early.trr

# PCA for second half
echo -e "C-alpha\nC-alpha\n" | gmx covar -s md.tpr -f md.xtc -b 50000 -e 100000 -o eigenvalues_late.xvg -v eigenvectors_late.trr
```

Compare eigenvalues and eigenvectors to identify changes in collective motions.

### Overlap Analysis

Calculate overlap between eigenvectors from different simulations:

```bash
echo -e "C-alpha\nC-alpha\n" | gmx anaeig -s md.tpr -f md.xtc -v eigenvectors1.trr -v2 eigenvectors2.trr -overlap overlap.xvg
```

Quantifies similarity between collective motions.

### Quasiharmonic Analysis

Calculate quasiharmonic entropy from PCA:

```bash
# Use eigenvalues to calculate entropy
# S_qh = (k_B/2) * sum(ln(k_B*T/λ_i))
```

Requires additional processing of eigenvalues.

### Essential Dynamics

Focus on essential subspace (first few PCs):

```bash
# Project trajectory onto essential subspace
# Analyze sampling within subspace
# Compare with experimental data
```

## Related Analyses

- **DCCM**: Analyzes correlated motions, complements PCA
- **FEL**: Maps conformational landscape using PCs as reaction coordinates
- **RMSF**: Provides per-residue flexibility information
- **Normal Mode Analysis**: Compares theoretical modes with PCA-derived modes

## Visualization Enhancements

### 3D Visualization

Project eigenvectors onto protein structure using PyMOL or VMD:

```bash
# Use gmx anaeig to generate eigenvector PDB
echo -e "C-alpha\nC-alpha\n" | gmx anaeig -s md.tpr -f md.xtc -v eigenvectors.trr -first 1 -last 1 -filt eigenvector1.pdb
```

### Porcupine Plots

Generate porcupine plots showing direction and magnitude of atomic motions:

```bash
# Requires external tools (e.g., NMWiz in VMD)
# Visualize eigenvectors as arrows on protein structure
```

### Cross-correlation of PCs

Analyze correlations between PCs:

```bash
# Calculate cross-correlation matrix of PC projections
# Use DuIvyTools or custom scripts
```

## Comparison with Experimental Data

### Compare with NMR Order Parameters

Project NMR order parameters onto PCA modes to validate motions.

### Compare with Crystal Structures

Compare PC1 extreme conformations with different crystal structures (e.g., open/closed states).

### Compare with SAXS

Calculate theoretical SAXS curves from PC-derived conformations and compare with experimental data.

## References

For theoretical background, see:
- Amadei et al. (1993) "Essential dynamics of proteins"
- García (1992) "Large-amplitude collective motions in proteins"
- Berendsen and Hayward (2000) "Collective variables and molecular dynamics"