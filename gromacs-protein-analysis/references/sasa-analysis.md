# Solvent Accessible Surface Area (SASA) Analysis

## Overview

Solvent Accessible Surface Area (SASA) measures the surface area of a protein accessible to solvent molecules. SASA analysis provides insights into protein-solvent interactions, hydrophobic/hydrophilic surface distribution, binding site accessibility, and protein unfolding.

## When to Use SASA

- Analyze solvent exposure of protein surfaces
- Identify hydrophobic and hydrophilic regions
- Study ligand binding site accessibility
- Monitor protein unfolding or aggregation
- Assess protein solubility
- Compare surface properties between states or mutants
- Validate protein structure stability

## Prerequisites

- Trajectory file (.xtc/.trr) - PBC correction optional (see [PBC Correction Guide](pbc-correction.md) if needed)
- Topology file (.tpr)
- Index file (.ndx) with appropriate atom groups
- GROMACS version with SASA support

## Workflow

### Step 1: Calculate Total SASA

Calculate SASA for different atom groups:

```bash
# SASA for protein backbone
echo -e "Backbone\n" | gmx sasa -s md.tpr -f md.xtc -o sas_backbone.xvg -tu ns

# SASA for C-alpha atoms
echo -e "C-alpha\n" | gmx sasa -s md.tpr -f md.xtc -o sas_calpha.xvg -tu ns

# SASA for all protein atoms
echo -e "Protein\n" | gmx sasa -s md.tpr -f md.xtc -o sas_protein.xvg -tu ns

# SASA for protein-ligand complex
echo -e "Protein_Lig\n" | gmx sasa -s md.tpr -f md.xtc -o sas_complex.xvg -tu ns
```

**Parameters**:
- `-s md.tpr`: Structure file
- `-f md.xtc`: Trajectory file
- `-o sas_*.xvg`: Output file for SASA data
- `-tu ns`: Time unit in nanoseconds

**Output files**:
- `sas_backbone.xvg`: Backbone SASA over time
- `sas_calpha.xvg`: C-alpha SASA over time
- `sas_protein.xvg`: Full protein SASA over time
- `sas_complex.xvg`: Complex SASA over time

### Step 2: Calculate Hydrophobic and Hydrophilic SASA

Separate SASA into hydrophobic and hydrophilic components:

```bash
# Total SASA with hydrophobic/hydrophilic breakdown
echo -e "Protein\n" | gmx sasa -s md.tpr -f md.xtc -o sas_protein.xvg -output sas_protein_components.xvg -surface Protein -tu ns
```

**Parameters**:
- `-output`: Output file with component breakdown
- `-surface`: Surface group for SASA calculation

**Output file**:
- `sas_protein_components.xvg`: Contains multiple columns:
  - Column 1: Time
  - Column 2: Total SASA
  - Column 3: Hydrophobic SASA
  - Column 4: Hydrophilic SASA

### Step 3: Calculate Per-Residue SASA (Optional)

Calculate SASA for each residue:

```bash
# Per-residue SASA
echo -e "Protein\n" | gmx sasa -s md.tpr -f md.xtc -o sas_residue.xvg -or sas_per_residue.xvg -res -tu ns
```

**Parameters**:
- `-or`: Output per-residue SASA
- `-res`: Calculate per-residue SASA

**Output files**:
- `sas_residue.xvg`: Total SASA over time
- `sas_per_residue.xvg`: Per-residue SASA (averaged over time)

### Step 4: Calculate SASA for Specific Domains (Optional)

Calculate SASA for protein domains:

```bash
# Create index file with domain groups
echo -e "r 1-100\nname 10 Domain1\nr 101-200\nname 11 Domain2\nq\n" | gmx make_ndx -f md.tpr -o domains.ndx

# Calculate SASA for Domain1
echo -e "Domain1\n" | gmx sasa -s md.tpr -f md.xtc -n domains.ndx -o sas_domain1.xvg
```

### Step 5: Adjust Probe Radius (Optional)

Use different probe radii for SASA calculation:

```bash
# Default probe radius: 0.14 nm (water)
# Use larger probe for analysis with larger solvent molecules
echo -e "Protein\n" | gmx sasa -s md.tpr -f md.xtc -o sas_large_probe.xvg -probe 0.2
```

**Parameters**:
- `-probe`: Probe radius in nm (default: 0.14 nm)

### Step 6: Visualize SASA

Generate SASA plots over time:

```bash
# Basic SASA plot
dit xvg_show -f sas_protein.xvg -x "Time (ns)" -y "SASA (nm²)" -t "Total SASA"

# Hydrophobic/hydrophilic breakdown
dit xvg_show -f sas_protein_components.xvg -c 1:4 -x "Time (ns)" -y "SASA (nm²)" -t "SASA Components" --multi

# Per-residue SASA
dit xvg_show -f sas_per_residue.xvg -x "Residue Number" -y "SASA (nm²)" -t "Per-Residue SASA"

# Compare multiple SASA curves
dit xvg_show -f sas_wt.xvg sas_mut1.xvg -x "Time (ns)" -y "SASA (nm²)" -t "SASA Comparison"
```

## Output Files

- **sas_protein.xvg**: Total SASA over time
- **sas_backbone.xvg**: Backbone SASA over time
- **sas_calpha.xvg**: C-alpha SASA over time
- **sas_complex.xvg**: Complex SASA over time
- **sas_protein_components.xvg**: Total, hydrophobic, and hydrophilic SASA
- **sas_per_residue.xvg**: Per-residue SASA (averaged)
- **sas_domain*.xvg**: Domain-specific SASA

## Interpretation Guidelines

### SASA Values

Typical SASA ranges for different protein sizes:

- **< 50 nm²**: Small proteins (< 50 residues)
- **50-100 nm²**: Small to medium proteins (50-100 residues)
- **100-200 nm²**: Medium proteins (100-200 residues)
- **200-300 nm²**: Medium to large proteins (200-300 residues)
- **300-500 nm²**: Large proteins (300-500 residues)
- **> 500 nm²**: Very large proteins (> 500 residues)

Empirical relationship: SASA ≈ 6.3 × N^0.73 (for globular proteins)

### SASA Pattern Analysis

- **Stable SASA**: Fluctuates around constant value (~5-10% variation)
- **Decreasing SASA**: Protein becoming more buried (folding or collapse)
- **Increasing SASA**: Protein becoming more exposed (unfolding or denaturation)
- **Abrupt changes**: Large conformational change or ligand binding/unbinding
- **Drift**: Continuous increase indicates unfolding

### Hydrophobic vs. Hydrophilic SASA

- **Hydrophobic SASA**: Non-polar surface area (typically 40-60% of total)
- **Hydrophilic SASA**: Polar surface area (typically 40-60% of total)
- **Hydrophobic/Hydrophilic ratio**: Indicates surface character
  - > 1.0: Predominantly hydrophobic (likely to aggregate)
  - ~1.0: Balanced surface (typical for soluble proteins)
  - < 1.0: Predominantly hydrophilic (highly soluble)

### Per-Residue SASA

- **Low SASA (< 0.1 nm²)**: Buried residues (protein core)
- **Moderate SASA (0.1-0.5 nm²)**: Partially exposed residues
- **High SASA (> 0.5 nm²)**: Surface-exposed residues

### Binding Site Analysis

- **Active site SASA**: Often moderate (0.2-0.4 nm²) for accessibility
- **Ligand binding**: Decrease in SASA upon ligand binding
- **Allosteric sites**: May show dynamic SASA changes

## Common Issues and Solutions

### Issue: SASA shows abrupt jumps

**Possible causes**:
- PBC artifacts (molecules crossing box boundaries)
- Protein centering issues
- Numerical instabilities

**Solutions**:
- Apply PBC correction: `gmx trjconv -pbc mol -center`
- Use `-pbc nojump` flag in `trjconv`
- Check trajectory for artifacts

### Issue: SASA continuously increases

**Possible causes**:
- Protein unfolding
- Simulation not converged
- Inadequate force field

**Solutions**:
- Check protein structure for unfolding
- Extend simulation time
- Verify force field compatibility

### Issue: SASA values are unexpectedly low

**Possible causes**:
- Incorrect atom selection
- Protein over-compact
- Probe radius too small

**Solutions**:
- Verify atom selection groups
- Check protein structure
- Adjust probe radius if needed

### Issue: Hydrophobic/hydrophilic breakdown not working

**Possible causes**:
- GROMACS version doesn't support component output
- Missing surface group specification

**Solutions**:
- Check GROMACS version (requires 2020 or later)
- Use `-surface` parameter to specify surface group

## Tips and Best Practices

- **Atom selection**: All-atom SASA is standard. Use backbone for efficiency.
- **Probe radius**: Default 0.14 nm (water) is appropriate for most analyses.
- **Time selection**: Use production phase only (exclude equilibration).
- **Multiple selections**: Calculate SASA for different atom groups.
- **Component analysis**: Separate hydrophobic and hydrophilic SASA for deeper insights.
- **Per-residue analysis**: Useful for identifying binding sites or aggregation-prone regions.
- **Statistical analysis**: Calculate average SASA, standard deviation, and confidence intervals.
- **Correlation with other metrics**: Combine with Rg, RMSD for comprehensive analysis.

## Advanced Analysis

### Time-dependent SASA Distribution

Calculate SASA distribution for different time windows:

```bash
# Extract SASA for first half and second half
# Compare distributions using histograms or KDE
# Identify changes in solvent exposure over time
```

### SASA Difference Between Simulations

Compare SASA between different conditions:

```bash
# Calculate SASA for each simulation
# Plot side-by-side for comparison
# Use statistical tests to assess significance
```

### SASA vs. Rg Correlation

Analyze relationship between SASA and Rg:

```bash
# Plot SASA vs. Rg
# Calculate correlation coefficient
# Identify unfolding intermediates
```

### Binding Site Accessibility

Monitor SASA of binding site residues:

```bash
# Extract binding site residues
# Calculate SASA for these residues over time
# Assess accessibility changes upon ligand binding
```

### Aggregation Propensity

Identify aggregation-prone regions:

```bash
# Calculate per-residue SASA
# Identify hydrophobic residues with high SASA
- Use aggregation prediction tools (e.g., TANGO, AGGRESCAN)
```

## Related Analyses

- **RMSD**: Provides structural stability, complement with SASA
- **Rg**: Related to compactness and SASA
- **RMSF**: Local flexibility may correlate with SASA changes
- **Hydrogen bonds**: Solvent-protein hydrogen bonds relate to hydrophilic SASA
- **FEL**: Use SASA as reaction coordinate for conformational landscape

## Visualization Enhancements

### Add Solvent Accessibility Threshold

Mark SASA threshold on plot:

```bash
dit xvg_show -f sas_protein.xvg -x "Time (ns)" -y "SASA (nm²)" -t "Total SASA" --hline 150 --hlabel "Native State"
```

### Hydrophobic/Hydrophilic Ratio

Calculate and plot hydrophobic/hydrophilic ratio:

```bash
# Extract hydrophobic and hydrophilic SASA from components file
# Calculate ratio: Hydrophobic / Hydrophilic
# Plot ratio over time
```

### Surface Mapping

Map SASA onto protein structure:

```bash
# Use PyMOL or VMD to color structure by SASA
# Red = high SASA (exposed)
# Blue = low SASA (buried)
```

### SASA Heatmap

Create heatmap of per-residue SASA for multiple simulations:

```bash
# Combine per-residue SASA data from multiple simulations
# Generate heatmap showing surface exposure patterns
```

## Interpretation Examples

### Example 1: Stable Soluble Protein

- SASA stable at ~150 nm²
- Small fluctuations (~10 nm²)
- Balanced hydrophobic/hydrophilic ratio (~1.0)
- **Interpretation**: Protein is well-folded and soluble

### Example 2: Unfolding Protein

- SASA increases from 150 nm² to 300 nm²
- Continuous drift throughout simulation
- Hydrophobic SASA increases significantly
- **Interpretation**: Protein is unfolding, exposing hydrophobic core

### Example 3: Ligand Binding

- SASA decreases by ~20 nm² upon ligand binding
- Decrease occurs at binding site residues
- **Interpretation**: Ligand binds to protein, burying surface area

### Example 4: Aggregation-Prone Protein

- High hydrophobic SASA (> 60% of total)
- Large SASA for hydrophobic surface residues
- **Interpretation**: Protein has aggregation propensity due to exposed hydrophobic surface

## Comparison with Experimental Data

### Hydrogen-Deuterium Exchange (HDX)

Compare SASA with HDX protection factors:

- Low SASA regions should be protected (slow exchange)
- High SASA regions should be accessible (fast exchange)
- Validates simulation surface exposure

### Small-Angle X-ray Scattering (SAXS)

Relate SASA to SAXS scattering:

- SASA affects scattering profile
- Compare simulation SASA with SAXS-derived values
- Good agreement validates simulation

### NMR Chemical Shifts

Compare SASA with NMR chemical shift perturbations:

- Surface-exposed residues show larger chemical shift changes
- Buried residues show minimal changes
- Validates simulation solvent accessibility

## Theoretical Background

### SASA Calculation

SASA is calculated using rolling probe method:

1. Place probe sphere (radius 0.14 nm for water) around each atom
2. Calculate surface area accessible to probe
3. Sum contributions from all atoms

### Relationship with Protein Size

For globular proteins:

```
SASA ≈ A0 × N^α
```

where:
- N = number of residues
- A0 ≈ 6.3 nm²
- α ≈ 0.73

### Hydrophobicity Scales

Residues classified as hydrophobic or hydrophilic based on scales:

- **Hydrophobic**: Ala, Val, Leu, Ile, Met, Phe, Trp, Pro, Tyr
- **Hydrophilic**: Asp, Glu, Asn, Gln, Lys, Arg, His, Ser, Thr, Cys

## Solvent Accessibility and Protein Function

### Active Sites

- Often have moderate SASA (balance of accessibility and specificity)
- May show dynamic SASA changes during catalysis

### Membrane Proteins

- Have distinct SASA patterns (hydrophobic transmembrane regions)
- May require analysis with membrane-mimicking solvents

### Protein-Protein Interactions

- Interface residues often have reduced SASA in complex
- Identify binding sites by SASA decrease upon complex formation

## References

For theoretical background, see:
- Lee and Richards (1971) "Interpretation of protein structures: estimation of static accessibility"
- Shrake and Rupley (1973) "Environment and exposure to solvent of protein atoms"
- Richmond (1984) "Solvent accessible surface area and excluded volume in proteins"