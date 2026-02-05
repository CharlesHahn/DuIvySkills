# Radius of Gyration (Gyrate) Analysis

## Overview

Radius of gyration (Rg) measures the mass-weighted root mean square distance of atoms from the center of mass. Gyrate analysis assesses protein compactness and provides insights into folding/unfolding transitions, global conformational changes, and overall protein size.

## When to Use Gyrate

- Monitor protein compactness during simulation
- Detect folding/unfolding transitions
- Assess global conformational changes
- Compare compactness between different states or mutants
- Validate protein stability
- Study domain movements or expansion/contraction
- Monitor protein aggregation propensity

## Prerequisites

- Trajectory file (.xtc/.trr) - PBC correction optional (see [PBC Correction Guide](pbc-correction.md) if needed)
- Topology file (.tpr)
- Index file (.ndx) with appropriate atom groups

## Workflow

### Step 1: Calculate Radius of Gyration

Calculate Rg for different atom groups:

```bash
# Rg for protein backbone
echo -e "Backbone\n" | gmx gyrate -s md.tpr -f md.xtc -o gyrate_backbone.xvg

# Rg for C-alpha atoms
echo -e "C-alpha\n" | gmx gyrate -s md.tpr -f md.xtc -o gyrate_calpha.xvg

# Rg for all protein atoms
echo -e "Protein\n" | gmx gyrate -s md.tpr -f md.xtc -o gyrate_protein.xvg

# Rg for protein-ligand complex
echo -e "Protein_Lig\n" | gmx gyrate -s md.tpr -f md.xtc -o gyrate_complex.xvg
```

**Parameters**:
- `-s md.tpr`: Structure file
- `-f md.xtc`: Trajectory file
- `-o gyrate_*.xvg`: Output file for Rg data

**Output files**:
- `gyrate_backbone.xvg`: Backbone Rg over time
- `gyrate_calpha.xvg`: C-alpha Rg over time
- `gyrate_protein.xvg`: Full protein Rg over time
- `gyrate_complex.xvg`: Complex Rg over time

### Step 2: Calculate Per-Axis Rg (Optional)

Calculate Rg for each principal axis separately:

```bash
echo -e "Protein\n" | gmx gyrate -s md.tpr -f md.xtc -o gyrate_axes.xvg -nm
```

**Parameters**:
- `-nm`: Output Rg for each axis (x, y, z)

**Output file**:
- `gyrate_axes.xvg`: Rg for each axis over time

**Columns in output**:
- Column 1: Time
- Column 2: Total Rg
- Column 3: Rg along x-axis
- Column 4: Rg along y-axis
- Column 5: Rg along z-axis

### Step 3: Calculate Rg for Specific Domains (Optional)

Calculate Rg for protein domains:

```bash
# Create index file with domain groups
echo -e "r 1-100\nname 10 Domain1\nr 101-200\nname 11 Domain2\nq\n" | gmx make_ndx -f md.tpr -o domains.ndx

# Calculate Rg for Domain1
echo -e "Domain1\n" | gmx gyrate -s md.tpr -f md.xtc -n domains.ndx -o gyrate_domain1.xvg
```

### Step 4: Visualize Rg

Generate Rg plots over time:

```bash
# Basic Rg plot
dit xvg_show -f gyrate_protein.xvg -x "Time (ns)" -y "Rg (nm)" -t "Radius of Gyration"

# Rg with smoothed line
dit xvg_show -f gyrate_protein.xvg -x "Time (ns)" -y "Rg (nm)" -t "Radius of Gyration" -m line --smooth

# Compare multiple Rg curves
dit xvg_show -f gyrate_backbone.xvg gyrate_calpha.xvg gyrate_protein.xvg -x "Time (ns)" -y "Rg (nm)" -t "Rg Comparison"

# Plot per-axis Rg
dit xvg_show -f gyrate_axes.xvg -c 1:3 -x "Time (ns)" -y "Rg (nm)" -t "Rg per Axis" --multi
```

### Step 5: Calculate Rg Distribution (Optional)

Analyze Rg distribution to understand compactness sampling:

```bash
# Extract Rg values
# Calculate histogram
# Fit to distribution model
```

## Output Files

- **gyrate_backbone.xvg**: Backbone Rg over time
- **gyrate_calpha.xvg**: C-alpha Rg over time
- **gyrate_protein.xvg**: Full protein Rg over time
- **gyrate_complex.xvg**: Complex Rg over time
- **gyrate_axes.xvg**: Rg for each axis over time
- **gyrate_domain*.xvg**: Domain-specific Rg

## Interpretation Guidelines

### Rg Values

Typical Rg ranges for different protein sizes:

- **< 1.0 nm**: Small proteins (< 50 residues)
- **1.0-1.5 nm**: Small to medium proteins (50-100 residues)
- **1.5-2.0 nm**: Medium proteins (100-200 residues)
- **2.0-2.5 nm**: Medium to large proteins (200-300 residues)
- **2.5-3.0 nm**: Large proteins (300-500 residues)
- **> 3.0 nm**: Very large proteins (> 500 residues)

Empirical relationship: Rg ≈ 0.395 × N^0.588 (for globular proteins)

### Rg Pattern Analysis

- **Stable Rg**: Fluctuates around constant value (~2-5% variation)
- **Decreasing Rg**: Protein becoming more compact (folding or collapse)
- **Increasing Rg**: Protein expanding (unfolding or denaturation)
- **Abrupt jumps**: Large conformational change or PBC artifacts
- **Drift**: Continuous increase/decrease indicates non-equilibrium

### Atom Selection Impact

- **Backbone**: Reflects overall protein shape, commonly used
- **C-alpha**: Similar to backbone, computationally efficient
- **All atoms**: Slightly larger values due to side chain volume
- **Heavy atoms**: Similar to all atoms, excludes hydrogens

### Per-Axis Analysis

- **Anisotropic Rg**: Different Rg values for x, y, z axes
- **Isotropic Rg**: Similar Rg values for all axes (spherical protein)
- **Shape information**: Ratio of axes indicates protein shape (elongated vs. spherical)

### Compactness Assessment

- **Compact globular**: Low Rg relative to size, isotropic
- **Extended conformation**: High Rg relative to size, anisotropic
- **Partially unfolded**: Intermediate Rg, may show drift
- **Fully unfolded**: High Rg, isotropic (random coil)

## Common Issues and Solutions

### Issue: Rg shows abrupt jumps

**Possible causes**:
- PBC artifacts (molecules crossing box boundaries)
- Protein centering issues
- Numerical instabilities

**Solutions**:
- Apply PBC correction: `gmx trjconv -pbc mol -center`
- Use `-pbc nojump` flag in `trjconv`
- Check trajectory for artifacts

### Issue: Rg continuously increases

**Possible causes**:
- Protein unfolding
- Simulation not converged
- Inadequate force field

**Solutions**:
- Check protein structure for unfolding
- Extend simulation time
- Verify force field compatibility

### Issue: Rg values are unexpectedly low

**Possible causes**:
- Incorrect atom selection
- Protein collapsed (over-compact)
- Trajectory centering issues

**Solutions**:
- Verify atom selection groups
- Check protein structure for unnatural collapse
- Verify trajectory processing

### Issue: Per-axis Rg shows large anisotropy

**Possible causes**:
- Protein naturally elongated (e.g., fibrous proteins)
- Trajectory not properly rotated
- Simulation box shape effects

**Solutions**:
- Check if anisotropy is expected for protein type
- Verify trajectory fitting
- Check simulation box dimensions

## Tips and Best Practices

- **Atom selection**: Backbone Rg is standard for overall compactness. Use C-alpha for efficiency.
- **Trajectory fitting**: Not strictly necessary for Rg, but recommended for consistency.
- **Time selection**: Use production phase only (exclude equilibration).
- **Multiple selections**: Calculate Rg for different atom groups to understand different aspects.
- **Per-axis analysis**: Use to identify anisotropic shape changes.
- **Statistical analysis**: Calculate average Rg, standard deviation, and confidence intervals.
- **Correlation with RMSD**: Combine Rg and RMSD for comprehensive stability assessment.
- **Visual inspection**: Examine trajectory to interpret large Rg changes.

## Advanced Analysis

### Time-dependent Rg Distribution

Calculate Rg distribution for different time windows:

```bash
# Extract Rg for first half and second half
# Compare distributions using histograms or KDE
# Identify changes in compactness over time
```

### Rg Difference Between Simulations

Compare Rg between different conditions:

```bash
# Calculate Rg for each simulation
# Plot side-by-side for comparison
# Use statistical tests to assess significance
```

### Correlation with RMSD

Analyze relationship between Rg and RMSD:

```bash
# Plot Rg vs. RMSD
# Calculate correlation coefficient
# Identify conformational transitions
```

### Principal Axes Analysis

Analyze protein shape using principal axes:

```bash
# Extract principal axes from covariance matrix
# Calculate axial ratios
# Relate to protein shape and function
```

### Folding/Unfolding Kinetics

Monitor Rg during folding/unfolding:

```bash
# Use Rg as reaction coordinate
- Identify intermediate states
- Calculate transition rates
- Map folding pathways
```

## Related Analyses

- **RMSD**: Provides structural stability, complement with Rg
- **RMSF**: Provides local flexibility, Rg provides global compactness
- **SASA**: Related to compactness and solvent exposure
- **PCA**: Identifies collective motions affecting Rg
- **FEL**: Use Rg as reaction coordinate for conformational landscape

## Visualization Enhancements

### Add Compactness Threshold

Mark compactness threshold on Rg plot:

```bash
dit xvg_show -f gyrate_protein.xvg -x "Time (ns)" -y "Rg (nm)" -t "Radius of Gyration" --hline 2.0 --hlabel "Compactness Threshold"
```

### Multiple Comparison Plots

Compare Rg from multiple simulations:

```bash
dit xvg_show -f gyrate_wt.xvg gyrate_mut1.xvg gyrate_mut2.xvg -x "Time (ns)" -y "Rg (nm)" -t "Rg Comparison"
```

### Rg vs. RMSD Scatter Plot

Analyze relationship between Rg and RMSD:

```bash
# Use custom script to plot Rg vs. RMSD
# Identify conformational states
```

### Shape Analysis

Visualize protein shape changes:

```bash
# Extract extreme conformations (min/max Rg)
# Superimpose and compare
# Identify regions responsible for expansion/contraction
```

## Interpretation Examples

### Example 1: Stable Globular Protein

- Rg stable at ~1.8 nm
- Small fluctuations (~0.05 nm)
- Isotropic per-axis Rg
- **Interpretation**: Protein is well-folded and stable

### Example 2: Unfolding Protein

- Rg increases from 1.8 nm to 3.5 nm
- Continuous drift throughout simulation
- Becomes more isotropic as it unfolds
- **Interpretation**: Protein is unfolding, may need different conditions

### Example 3: Conformational Change

- Rg stable at 1.8 nm for first 50 ns
- Sudden jump to 2.2 nm at 50 ns
- Stable at new level for remaining 50 ns
- **Interpretation**: Protein undergoes conformational transition to more extended state

### Example 4: Elongated Protein

- Rg stable at 2.5 nm
- Strong anisotropy (x: 2.0, y: 2.8, z: 3.0 nm)
- **Interpretation**: Protein is naturally elongated (e.g., fibrous protein)

## Comparison with Experimental Data

### Small-Angle X-ray Scattering (SAXS)

Compare Rg with SAXS data:

- Extract Rg from SAXS scattering curve
- Compare with simulation Rg
- Good agreement validates simulation

### Hydrodynamic Radius

Relate Rg to hydrodynamic radius (Rh):

```
Rh ≈ Rg / 0.775 (for globular proteins)
```

Compare with experimental Rh from DLS or viscometry.

### Crystallographic Packing

Compare simulation Rg with crystal structure:

- Crystal structure may be more compact due to packing
- Simulation should show slightly larger Rg in solution

## Theoretical Background

### Rg Calculation

Rg is calculated as:

```
Rg² = Σmi × ri² / Σmi
```

where mi is mass of atom i and ri is distance from center of mass.

### Relationship with Protein Size

For globular proteins:

```
Rg ≈ R0 × N^ν
```

where:
- N = number of residues
- R0 ≈ 0.395 nm
- ν ≈ 0.588 (for unfolded: ν ≈ 0.588, for folded: ν ≈ 0.33)

### Anisotropy Factor

Calculate anisotropy factor:

```
A = 1 - 3 × (Rg_x × Rg_y × Rg_z) / (Rg_x² + Rg_y² + Rg_z²)^1.5
```

- A ≈ 0: Isotropic (spherical)
- A → 1: Highly anisotropic (elongated)

## References

For theoretical background, see:
- Flory (1953) "Principles of Polymer Chemistry"
- Miller et al. (2002) "Radius of gyration of proteins"
- Kohn et al. (2004) "Dynamics of protein folding"