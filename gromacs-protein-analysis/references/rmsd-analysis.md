# Root Mean Square Deviation (RMSD) Analysis

## Overview

Root Mean Square Deviation (RMSD) measures the average distance between atoms in a protein structure relative to a reference structure. RMSD is the most common metric for assessing structural stability, identifying equilibration phases, and evaluating simulation convergence.

## When to Use RMSD

- Assess overall structural stability during simulation
- Identify equilibration phase (when RMSD stabilizes)
- Evaluate simulation convergence
- Compare different simulation conditions or mutants
- Monitor protein unfolding or large conformational changes
- Validate simulation against experimental structures

## Prerequisites

- Trajectory file (.xtc/.trr) - PBC correction optional (see [PBC Correction Guide](pbc-correction.md) if RMSD shows abrupt jumps)
- Topology file (.tpr)
- Reference structure (usually the first frame or a known structure)
- Index file (.ndx) with appropriate atom groups

## Workflow

### Step 1: Fit Trajectory to Remove Overall Motion

Fit trajectory to remove overall translation and rotation:

```bash
echo -e "Protein\nProtein\n" | gmx trjconv -s md.tpr -f md.xtc -o fit.xtc -pbc nojump -fit rot+trans
```

- First input: Select group for centering (Protein)
- Second input: Select group for fitting (Protein)

**Output file**:
- `fit.xtc`: Fitted trajectory with overall motion removed

### Step 2: Calculate RMSD

Calculate RMSD using the fitted trajectory:

```bash
# RMSD for protein backbone
echo -e "Backbone\nBackbone\n" | gmx rms -s md.tpr -f fit.xtc -o rmsd_backbone.xvg -tu ns

# RMSD for C-alpha atoms
echo -e "C-alpha\nC-alpha\n" | gmx rms -s md.tpr -f fit.xtc -o rmsd_calpha.xvg -tu ns

# RMSD for all protein atoms
echo -e "Protein\nProtein\n" | gmx rms -s md.tpr -f fit.xtc -o rmsd_protein.xvg -tu ns
```

**Parameters**:
- `-s md.tpr`: Structure file for reference (uses first frame by default)
- `-f fit.xtc`: Fitted trajectory file
- `-o rmsd_*.xvg`: Output file for RMSD data
- `-tu ns`: Time unit in nanoseconds

**Output files**:
- `rmsd_backbone.xvg`: RMSD for backbone atoms
- `rmsd_calpha.xvg`: RMSD for C-alpha atoms
- `rmsd_protein.xvg`: RMSD for all protein atoms

### Step 3: Use Custom Reference Structure (Optional)

Calculate RMSD relative to a specific reference structure:

```bash
echo -e "Backbone\nBackbone\n" | gmx rms -s reference.pdb -f fit.xtc -o rmsd_to_ref.xvg -tu ns
```

- `-s reference.pdb`: Use custom reference structure instead of first frame

### Step 4: Calculate RMSD for Specific Domains (Optional)

Calculate RMSD for protein domains or subunits:

```bash
# Create index file with domain groups
echo -e "r 1-100\nname 10 Domain1\nr 101-200\nname 11 Domain2\nq\n" | gmx make_ndx -f md.tpr -o domains.ndx

# Calculate RMSD for Domain1
echo -e "Domain1\nDomain1\n" | gmx rms -s md.tpr -f fit.xtc -n domains.ndx -o rmsd_domain1.xvg
```

### Step 5: Visualize RMSD

Generate RMSD plots over time:

```bash
# Basic RMSD plot
dit xvg_show -f rmsd_backbone.xvg -x "Time (ns)" -y "RMSD (nm)" -t "Backbone RMSD"

# RMSD with smoothed line
dit xvg_show -f rmsd_backbone.xvg -x "Time (ns)" -y "RMSD (nm)" -t "Backbone RMSD" -m line --smooth

# Compare multiple RMSD curves
dit xvg_show -f rmsd_backbone.xvg rmsd_calpha.xvg rmsd_protein.xvg -x "Time (ns)" -y "RMSD (nm)" -t "RMSD Comparison"
```

## Output Files

- **rmsd_backbone.xvg**: Backbone RMSD over time
- **rmsd_calpha.xvg**: C-alpha RMSD over time
- **rmsd_protein.xvg**: Full protein RMSD over time
- **rmsd_to_ref.xvg**: RMSD relative to custom reference
- **rmsd_domain*.xvg**: Domain-specific RMSD

## Interpretation Guidelines

### RMSD Values

Typical RMSD ranges for stable protein simulations:

- **< 0.1 nm (1 Å)**: Very stable, minimal conformational change
- **0.1-0.2 nm (1-2 Å)**: Good stability, typical for well-folded proteins
- **0.2-0.3 nm (2-3 Å)**: Moderate flexibility, acceptable for many proteins
- **0.3-0.5 nm (3-5 Å)**: Significant flexibility or large conformational change
- **> 0.5 nm (5 Å)**: Large structural change, possible unfolding or rearrangement

### RMSD Pattern Analysis

- **Equilibration phase**: RMSD increases rapidly in first few ns, then stabilizes
- **Stable phase**: RMSD fluctuates around a constant value
- **Convergence**: RMSD reaches plateau and shows stable fluctuations
- **Drift**: Continuous increase in RMSD indicates poor convergence
- **Jumps**: Sudden RMSD jumps may indicate PBC artifacts or conformational changes

### Atom Selection Impact

- **Backbone**: Reflects overall protein stability, commonly used
- **C-alpha**: Similar to backbone, computationally efficient
- **All atoms**: More sensitive to side chain motions, larger values
- **Heavy atoms**: Excludes hydrogens, reduces noise

### Equilibration Assessment

Determine equilibration time from RMSD:

1. Plot RMSD over time
2. Identify when RMSD reaches plateau
3. Use statistical tests (e.g., block averaging) to confirm convergence
4. Exclude equilibration phase from production analysis

## Common Issues and Solutions

### Issue: RMSD shows abrupt jumps

**Possible causes**:
- PBC artifacts (molecules crossing box boundaries)
- Protein centering issues
- Numerical instabilities

**Solutions**:
- Apply PBC correction: `gmx trjconv -pbc mol -center`
- Use `-pbc nojump` flag in `trjconv`
- Fit trajectory before RMSD calculation

### Issue: RMSD continuously increases without plateau

**Possible causes**:
- Simulation not converged
- Protein unfolding
- Inadequate force field parameters

**Solutions**:
- Extend simulation time
- Check protein structure for unfolding
- Verify force field compatibility

### Issue: RMSD values are unexpectedly high

**Possible causes**:
- Wrong reference structure
- Incorrect atom selection
- Large conformational change

**Solutions**:
- Verify reference structure matches simulation system
- Check atom selection groups
- Examine trajectory for large-scale motions

### Issue: RMSD shows very low values (< 0.05 nm)

**Possible causes**:
- Reference structure is from same trajectory (first frame)
- Protein is very rigid
- Fitting removes too much motion

**Solutions**:
- Use different reference structure (e.g., crystal structure)
- Check if RMSD reflects real structural changes
- Consider alternative metrics (e.g., RMSF)

## Tips and Best Practices

- **Atom selection**: Backbone RMSD is standard for overall stability. Use C-alpha for efficiency.
- **Reference structure**: Use first frame for self-convergence, use crystal structure for validation.
- **Trajectory fitting**: Always fit trajectory before RMSD to remove overall motion.
- **Time selection**: Plot full trajectory to identify equilibration phase. Use production phase for analysis.
- **Multiple selections**: Calculate RMSD for different atom groups to understand different aspects of stability.
- **Statistical analysis**: Calculate average RMSD, standard deviation, and confidence intervals for production phase.
- **Visual inspection**: Combine RMSD analysis with trajectory visualization to interpret large deviations.

## Advanced Analysis

### Time-averaged RMSD

Calculate rolling average RMSD to identify trends:

```bash
# Use custom script or tools to calculate moving average
# Example: 100 ps window moving average
```

### RMSD Distribution

Analyze RMSD distribution to understand sampling:

```bash
# Extract RMSD values
# Calculate histogram
# Fit to distribution model (e.g., Gaussian)
```

### RMSD Difference Between Simulations

Compare RMSD between different simulations:

```bash
# Calculate RMSD for each simulation
# Plot side-by-side for comparison
# Use statistical tests (t-test, KS test) to assess significance
```

### Per-residue RMSD

Calculate RMSD contribution per residue (similar to RMSF):

```bash
# This is typically done using RMSF analysis
# See RMSF Analysis Guide for details
```

## Related Analyses

- **RMSF**: Provides per-residue flexibility, complements RMSD
- **Gyrate**: Assesses protein compactness and unfolding
- **PCA**: Identifies collective motions beyond overall RMSD
- **SASA**: Monitors solvent accessibility and surface changes
- **PBC Correction**: Required if RMSD shows jumps due to PBC artifacts

## Visualization Enhancements

### Add Equilibration Line

Mark equilibration phase on RMSD plot:

```bash
dit xvg_show -f rmsd_backbone.xvg -x "Time (ns)" -y "RMSD (nm)" -t "Backbone RMSD" --vline 5 --vlabel "Equilibration"
```

### Multiple Comparison Plots

Compare RMSD from multiple simulations:

```bash
dit xvg_show -f rmsd_sim1.xvg rmsd_sim2.xvg rmsd_sim3.xvg -x "Time (ns)" -y "RMSD (nm)" -t "RMSD Comparison"
```

### RMSD vs Temperature

Monitor RMSD at different simulation conditions:

```bash
# Plot RMSD for simulations at different temperatures
# Analyze temperature dependence of stability
```

## Interpretation Examples

### Example 1: Well-Equilibrated Simulation

- RMSD increases from 0 to 0.15 nm in first 2 ns
- Stabilizes at ~0.15 nm for remaining 98 ns
- Small fluctuations (~0.02 nm) around mean
- **Interpretation**: Simulation equilibrated after 2 ns, good stability

### Example 2: Unfolding Protein

- RMSD increases gradually throughout simulation
- Reaches > 0.8 nm by end of simulation
- No clear plateau
- **Interpretation**: Protein is unfolding, may need longer simulation or different conditions

### Example 3: Conformational Change

- RMSD stable at 0.2 nm for first 50 ns
- Sudden jump to 0.4 nm at 50 ns
- Stable at new level for remaining 50 ns
- **Interpretation**: Protein undergoes conformational transition, both states stable

## References

For theoretical background, see:
- Marti-Renom et al. (2002) "Comparative protein structure modeling of genes and genomes"
- Kuhlman and Baker (2000) "Native protein sequences are close to optimal for their structures"