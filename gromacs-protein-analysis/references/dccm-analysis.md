# Dynamics Cross-Correlation Matrix (DCCM) Analysis

## Overview

Dynamics Cross-Correlation Matrix (DCCM) analyzes correlated motions between atomic pairs in a protein. Positive correlations (red) indicate atoms moving together, while negative correlations (blue) indicate opposite movements.

## When to Use DCCM

- Identify allosteric communication pathways
- Study correlated motions in protein dynamics
- Analyze domain movements
- Understand collective behavior in protein complexes
- Validate experimental observations of correlated motions

## Prerequisites

- PBC-corrected trajectory (see [PBC Correction Guide](pbc-correction.md))
- Topology file (.tpr)
- Index file (.ndx) with appropriate atom groups
- Trajectory should be fitted to remove overall translation/rotation

## Workflow

### Step 1: Calculate Covariance Matrix

Compute the covariance matrix from the trajectory:

```bash
echo -e "C-alpha\nC-alpha\n" | gmx covar -s md.tpr -f md.xtc -o eigenvalues.xvg -v eigenvectors.trr -xpma covar.xpm -ascii covar.dat
```

- First input: Select reference group for alignment (C-alpha)
- Second input: Select group for covariance calculation (C-alpha)

**Output files**:
- `eigenvalues.xvg`: Eigenvalues from covariance matrix diagonalization
- `eigenvectors.trr`: Eigenvectors stored as trajectory
- `covar.xpm`: Covariance matrix in XPM format
- `covar.dat`: Covariance matrix in ASCII format

### Step 2: Convert to DCCM

Convert covariance matrix to dynamic cross-correlation matrix:

```bash
dit dccm_ascii -f covar.dat -o dccm.xpm
```

This step requires DuIvyTools. The DCCM values range from -1 (perfectly anti-correlated) to +1 (perfectly correlated).

### Step 3: Visualize DCCM

Generate heatmap visualization of the DCCM:

```bash
dit xpm_show -f dccm.xpm -o dccm.png -zmin -1 -zmax 1 -cmap bwr -m contour
```

**Parameters**:
- `-zmin -1 -zmax 1`: Set color range for correlation values
- `-cmap bwr`: Use blue-white-red colormap (blue=negative, white=zero, red=positive)
- `-m contour`: Use contour plot mode

### Step 4: Analyze Correlation Patterns

Examine the DCCM visualization for:

- **Positive correlations (red)**: Atoms moving together
- **Negative correlations (blue)**: Atoms moving in opposite directions
- **Diagonal elements**: Always +1 (self-correlation)
- **Symmetry**: Matrix should be symmetric across diagonal
- **Patterns**: Blocks or stripes indicating correlated domains or pathways

## Output Files

- **covar.dat**: Raw covariance matrix (ASCII format)
- **covar.xpm**: Covariance matrix (XPM format)
- **dccm.xpm**: Dynamic cross-correlation matrix
- **dccm.png**: Visualization of DCCM
- **eigenvalues.xvg**: Eigenvalues (useful for PCA)
- **eigenvectors.trr**: Eigenvectors (useful for PCA)

## Interpretation Guidelines

### Correlation Strength

- **> 0.5**: Strong positive correlation
- **0.3 to 0.5**: Moderate positive correlation
- **0.1 to 0.3**: Weak positive correlation
- **-0.1 to 0.1**: No significant correlation
- **-0.3 to -0.1**: Weak negative correlation
- **-0.5 to -0.3**: Moderate negative correlation
- **< -0.5**: Strong negative correlation

### Biological Interpretation

- **Strong positive correlations**: Rigid domains moving together, allosteric coupling
- **Strong negative correlations**: Antagonistic motions, hinge movements, breathing modes
- **No correlations**: Independent motions, flexible regions
- **Long-range correlations**: Potential allosteric communication pathways
- **Local correlations**: Secondary structure stability

## Common Issues and Solutions

### Issue: DCCM shows uniform values near zero

**Possible causes**:
- Insufficient simulation time or poor sampling
- Too few atoms in selection
- Excessive noise from high-frequency motions

**Solutions**:
- Extend simulation time
- Use larger atom groups (e.g., full backbone instead of C-alpha only)
- Apply time averaging or filter high-frequency motions

### Issue: DCCM appears asymmetric

**Solution**: Ensure proper fitting and alignment before analysis. The DCCM matrix should be symmetric.

### Issue: DCCM values exceed [-1, 1] range

**Solution**: This should not happen with properly calculated DCCM. Verify the `dit dccm_ascii` command and input data.

### Issue: Diagonal elements are not 1

**Solution**: This indicates an error in calculation. Check the covariance matrix calculation step.

## Tips and Best Practices

- **Atom selection**: C-alpha atoms are commonly used for efficiency. For detailed analysis, use backbone or all protein atoms.
- **Time selection**: Exclude equilibration phase. Use production phase only.
- **Trajectory fitting**: Always remove overall translation/rotation before analysis.
- **Averaging**: For noisy data, consider time-averaged DCCM over simulation windows.
- **Validation**: Compare DCCM patterns with known protein behavior or experimental data.

## Advanced Analysis

### Time-dependent DCCM

Calculate DCCM for different time windows to study evolution of correlations:

```bash
# Divide trajectory into windows
echo -e "C-alpha\nC-alpha\n" | gmx covar -s md.tpr -f md.xtc -b 0 -e 5000 -o eigenvalues_0-5ns.xvg -xpma covar_0-5ns.xpm -ascii covar_0-5ns.dat
echo -e "C-alpha\nC-alpha\n" | gmx covar -s md.tpr -f md.xtc -b 5000 -e 10000 -o eigenvalues_5-10ns.xvg -xpma covar_5-10ns.xpm -ascii covar_5-10ns.dat
```

Convert each window to DCCM and compare.

### Domain-specific DCCM

Calculate DCCM for specific domains or regions:

```bash
# Create index file with domain groups
echo -e "r 1-100\nname 10 Domain1\nr 101-200\nname 11 Domain2\nq\n" | gmx make_ndx -f md.tpr -o domains.ndx
```

Then use specific groups for covariance calculation.

### Network Analysis

Convert DCCM to correlation network for advanced analysis:

1. Apply threshold to identify significant correlations
2. Construct network with atoms as nodes and correlations as edges
3. Analyze network properties (centrality, communities, pathways)

## Related Analyses

- **PCA**: Identifies collective motions, complements DCCM analysis
- **RMSF**: Provides per-residue flexibility information
- **Contact Analysis**: Identifies persistent contacts alongside correlated motions
- **FEL**: Maps conformational states identified by correlated motions

## Visualization Enhancements

### Custom Color Schemes

Use different colormaps to highlight features:

```bash
# Coolwarm (diverging)
dit xpm_show -f dccm.xpm -o dccm_coolwarm.png -zmin -1 -zmax 1 -cmap coolwarm

# Viridis (sequential, alternative view)
dit xpm_show -f dccm.xpm -o dccm_viridis.png -zmin -1 -zmax 1 -cmap viridis
```

### Zooming on Regions

Focus on specific regions of the DCCM:

```bash
# Zoom on residues 50-100
dit xpm_show -f dccm.xpm -o dccm_zoom.png -xmin 50 -xmax 100 -ymin 50 -ymax 100 -zmin -1 -zmax 1 -cmap bwr
```

### 3D Visualization

Project correlations onto protein structure using PyMOL or VMD plugins.

## References

For theoretical background, see:
- Lange and Grubmüller (2006) "Full correlation analysis of protein dynamics"
- Hünenberger et al. (1995) "Dynamical properties of the solvent and the solute in proteins"