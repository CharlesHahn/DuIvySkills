# Free Energy Landscape (FEL) Analysis

## Overview

Free Energy Landscape (FEL) maps the energy surface of protein conformations, identifying stable states, energy barriers, and transition pathways. Two common approaches use either RMSD + Gyrate or Principal Components (PCs) as reaction coordinates.

## When to Use FEL

- Identify stable conformational states
- Quantify energy barriers between states
- Understand conformational transitions
- Analyze protein folding/unfolding
- Study domain movements and allostery
- Compare free energy surfaces between conditions

## Two Approaches

### Method 1: RMSD + Gyrate

Use structural deviation and compactness as reaction coordinates.

**Advantages**:
- Intuitive interpretation
- Directly relates to structural changes
- Easy to calculate

**Disadvantages**:
- May not capture subtle conformational changes
- Sensitive to reference structure selection

### Method 2: Principal Components

Use PC projections as reaction coordinates.

**Advantages**:
- Captures collective motions
- Based on natural coordinates from dynamics
- Reveals dominant conformational changes

**Disadvantages**:
- Requires prior PCA analysis
- Less intuitive interpretation
- May miss non-dominant but important motions

## Prerequisites

### For Method 1 (RMSD + Gyrate)
- Trajectory file (.xtc/.trr) - PBC correction optional
- RMSD data (rmsd.xvg)
- Gyrate data (gyrate.xvg)
- Topology file (.tpr)

### For Method 2 (PCA)
- Trajectory file (.xtc/.trr) - PBC correction optional
- PCA eigenvectors (eigenvectors.trr)
- PC projections (pc1.xvg, pc2.xvg)
- Topology file (.tpr)

## Workflow: Method 1 (RMSD + Gyrate)

### Step 1: Calculate RMSD

Calculate RMSD relative to reference structure:

```bash
echo -e "Backbone\nProtein\n" | gmx rms -s md.tpr -f md.xtc -o rmsd.xvg
```

- First input: Select reference for alignment (Backbone)
- Second input: Select group for RMSD calculation (Protein)

**Output**: `rmsd.xvg` (time, RMSD)

### Step 2: Calculate Gyrate

Calculate radius of gyration:

```bash
echo -e "Protein\n" | gmx gyrate -s md.tpr -f md.xtc -o gyrate.xvg
```

- Input: Select Protein group

**Output**: `gyrate.xvg` (time, Rg, Rg_x, Rg_y, Rg_z)
- Column 1: Time
- Column 2: Total Rg
- Columns 3-5: X, Y, Z components

### Step 3: Combine RMSD and Gyrate Data

Combine RMSD and total Rg into single file for sham analysis:

```bash
dit xvg_combine -f rmsd.xvg gyrate.xvg -c 0,1 1 -o sham.xvg
```

**Parameters**:
- `-c 0,1 1`: Extract columns (time, RMSD) from first file and column 1 (Rg) from second file
- Creates `sham.xvg` with columns: time, RMSD, Rg

### Step 4: Generate Free Energy Landscape

Use gmx sham to calculate FEL:

```bash
gmx sham -tsham 310 -nlevels 100 -f sham.xvg -ls gibbs.xpm -g gibbs.log -lsh enthalpy.xpm -lss entropy.xpm
```

**Parameters**:
- `-tsham 310`: Temperature in Kelvin
- `-nlevels 100`: Number of energy levels (resolution)
- `-f sham.xvg`: Input file with reaction coordinates
- `-ls gibbs.xpm`: Output Gibbs free energy landscape
- `-g gibbs.log`: Log file with energy minima
- `-lsh enthalpy.xpm`: Enthalpy landscape
- `-lss entropy.xpm`: Entropy landscape

**Output files**:
- `gibbs.xpm`: Gibbs free energy landscape (XPM format)
- `gibbs.log`: Log file with energy minima and indices
- `enthalpy.xpm`: Enthalpy landscape
- `entropy.xpm`: Entropy landscape
- `bindex.ndx`: Index file with frame assignments
- `ener.xvg`: Energy values

### Step 5: Visualize FEL

Generate visualization of free energy landscape:

```bash
dit xpm_show -f gibbs.xpm -m contour -cmap jet -o fel.png
```

**Parameters**:
- `-m contour`: Contour plot mode
- `-cmap jet`: Jet colormap (blue=low energy, red=high energy)
- Set axis labels: `-x "RMSD (nm)" -y "Rg (nm)"`

### Step 6: Extract Lowest Energy Conformations

Identify and extract structures corresponding to energy minima:

```bash
# View energy minima
cat gibbs.log

# Example output:
# Minimum 0 at index 26 energy 7.589
# Minimum 1 at index 46 energy 7.589
# Minimum 2 at index 50 energy 7.589
# Minimum 3 at index 56 energy 7.589
# Minimum 4 at index 141 energy 5.803

# Check frames corresponding to minima
cat bindex.ndx

# Example output:
# [ 26 ]
# 1274
# [ 46 ]
# 2
# [ 141 ]
# 4
# 1282
```

Extract conformation at specific time:

```bash
# Extract frame 1274
echo -e "Protein\n" | gmx trjconv -f md.xtc -s md.tpr -b <time_at_frame_1274> -e <time_at_frame_1274> -o min1.pdb
```

To find time at specific frame, check sham.xvg:

```bash
# Frame number corresponds to line number in sham.xvg
# Time at frame N = value in column 0 of line N in sham.xvg
```

## Workflow: Method 2 (Principal Components)

### Step 1: Perform PCA Analysis

Complete PCA analysis (see [PCA Analysis Guide](pca-analysis.md)) to obtain:
- `eigenvectors.trr`: Eigenvectors
- `pc1.xvg`: Projection onto PC1
- `pc2.xvg`: Projection onto PC2

### Step 2: Combine PC Projections

Combine PC1 and PC2 projections:

```bash
dit xvg_combine -f pc1.xvg pc2.xvg -c 0,1 1 -o sham.xvg
```

**Parameters**:
- `-c 0,1 1`: Extract columns (time, PC1) from first file and column 1 (PC2) from second file
- Creates `sham.xvg` with columns: time, PC1, PC2

### Step 3: Generate Free Energy Landscape

Same as Method 1, Step 4:

```bash
gmx sham -tsham 310 -nlevels 100 -f sham.xvg -ls gibbs.xpm -g gibbs.log -lsh enthalpy.xpm -lss entropy.xpm
```

### Step 4: Visualize FEL

Generate visualization:

```bash
dit xpm_show -f gibbs.xpm -m contour -cmap jet -o fel.png -x "PC1" -y "PC2" -t "Free Energy Landscape"
```

### Step 5: Extract Lowest Energy Conformations

Same as Method 1, Steps 5-6.

## Output Files

- **sham.xvg**: Combined reaction coordinate data
- **gibbs.xpm**: Gibbs free energy landscape (XPM format)
- **gibbs.log**: Energy minima and indices
- **enthalpy.xpm**: Enthalpy landscape
- **entropy.xpm**: Entropy landscape
- **bindex.ndx**: Frame assignments for energy minima
- **ener.xvg**: Energy values
- **fel.png**: Visualization of free energy landscape
- **min*.pdb**: Structures at energy minima

## Interpretation Guidelines

### Energy Minima

- **Deep minima**: Stable conformational states
- **Shallow minima**: Metastable states or transient conformations
- **Multiple minima**: Existence of multiple stable states
- **Single minimum**: Single dominant conformation

### Energy Barriers

- **High barriers**: Slow transitions, rare events
- **Low barriers**: Rapid transitions, frequent interconversions
- **Barrier height**: Related to transition rates (k ∝ exp(-ΔG/kT))

### Landscape Topology

- **Funnel shape**: Protein folding landscape
- **Multiple funnels**: Multiple folding pathways
- **Rugged landscape**: Many local minima, glassy behavior
- **Smooth landscape**: Few barriers, rapid dynamics

### Reaction Coordinates

**RMSD + Gyrate**:
- **High RMSD, High Rg**: Unfolded or extended conformations
- **Low RMSD, Low Rg**: Compact, native-like conformations
- **Low RMSD, High Rg**: Unfolded but structured
- **High RMSD, Low Rg**: Misfolded or collapsed structures

**Principal Components**:
- **PC1 axis**: Dominant motion
- **PC2 axis**: Second dominant motion
- **Regions**: Different conformational states along collective motions

## Common Issues and Solutions

### Issue: FEL shows unrealistic energy barriers (> 20 kT)

**Possible causes**:
- Insufficient sampling
- Incorrect time range selection
- Temperature mismatch

**Solutions**:
- Extend simulation time
- Use consistent time range
- Verify temperature setting in sham command

### Issue: FEL appears noisy or fragmented

**Solutions**:
- Increase number of levels (`-nlevels`)
- Extend simulation time for better sampling
- Smooth data before analysis

### Issue: Multiple minima merge into single basin

**Solutions**:
- Increase resolution (more levels)
- Use different reaction coordinates
- Check if minima are physically distinct

### Issue: FEL shows no clear minima

**Solutions**:
- Protein may be rigid (single conformation)
- Check simulation quality
- Consider using different reaction coordinates

## Tips and Best Practices

- **Time selection**: Use production phase only, exclude equilibration
- **Temperature**: Match simulation temperature in sham command
- **Resolution**: Adjust `-nlevels` for appropriate detail
- **Sampling**: Ensure sufficient sampling of conformational space
- **Comparison**: Compare FELs from different conditions or mutants
- **Validation**: Relate energy minima to experimental structures

## Advanced Analysis

### Time-dependent FEL

Calculate FEL for different time windows:

```bash
# FEL for first half
dit xvg_combine -f rmsd_early.xvg gyrate_early.xvg -c 0,1 1 -o sham_early.xvg
gmx sham -tsham 310 -nlevels 100 -f sham_early.xvg -ls gibbs_early.xpm

# FEL for second half
dit xvg_combine -f rmsd_late.xvg gyrate_late.xvg -c 0,1 1 -o sham_late.xvg
gmx sham -tsham 310 -nlevels 100 -f sham_late.xvg -ls gibbs_late.xvg
```

Compare FELs to identify changes in conformational landscape.

### Difference FEL

Calculate difference between two FELs:

```bash
dit xpm_diff -f gibbs_early.xpm gibbs_late.xpm -o gibbs_diff.xpm
```

Identifies regions where free energy changed significantly.

### Pathway Analysis

Identify minimum energy paths between states:

```bash
# Use specialized tools (e.g., string method, nudged elastic band)
# Requires additional processing
```

### Multi-dimensional FEL

Extend to 3D using additional reaction coordinates:

```bash
# Combine three coordinates
dit xvg_combine -f pc1.xvg pc2.xvg pc3.xvg -c 0,1 1 1 -o sham_3d.xvg
gmx sham -tsham 310 -nlevels 100 -f sham_3d.xvg -ls gibbs_3d.xpm
```

Visualize 3D FEL using plotly:

```bash
dit xpm_show -f gibbs_3d.xpm -m 3d -eg plotly -o fel_3d.html
```

## Related Analyses

- **PCA**: Provides collective motions for FEL
- **DCCM**: Identifies correlated motions in energy minima
- **RMSF**: Analyzes flexibility at energy minima
- **Contact Analysis**: Identifies contacts in different states

## Visualization Enhancements

### Custom Color Schemes

Use different colormaps:

```bash
# Jet (default)
dit xpm_show -f gibbs.xpm -m contour -cmap jet -o fel_jet.png

# Viridis (perceptually uniform)
dit xpm_show -f gibbs.xpm -m contour -cmap viridis -o fel_viridis.png

# Plasma (alternative)
dit xpm_show -f gibbs.xpm -m contour -cmap plasma -o fel_plasma.png
```

### 3D Visualization

Generate 3D surface plot:

```bash
dit xpm_show -f gibbs.xpm -m 3d -eg plotly -o fel_3d.html
```

### Overlay Conformations

Extract and overlay structures from different minima using PyMOL or VMD.

### Energy Profiles

Extract energy profiles along specific paths:

```bash
# Extract energy values along specific PC1 value
# Requires custom processing of gibbs.xpm
```

## Comparison with Experimental Data

### Compare with NMR Ensemble

Overlay NMR structures on FEL to validate conformational sampling.

### Compare with Crystal Structures

Locate crystal structures on FEL to understand their energetic context.

### Compare with SAXS

Calculate theoretical SAXS curves from FEL-derived conformations.

## References

For theoretical background, see:
- Protein folding theory and energy landscapes
- Transition state theory and barrier crossing
- Conformational dynamics and free energy calculations