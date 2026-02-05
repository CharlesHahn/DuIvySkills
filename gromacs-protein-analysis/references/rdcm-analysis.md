# Residue Distance Contact Matrix (RDCM) Analysis

## Overview

Residue Distance Contact Matrix (RDCM) calculates average distances between residue pairs throughout the simulation, providing insights into inter-residue spatial relationships, long-range contacts, and protein folding patterns.

## When to Use RDCM

- Analyze inter-residue contacts and distances
- Identify long-range interactions
- Study protein folding/unfolding transitions
- Analyze domain-domain contacts
- Validate contact maps from experimental data
- Study conformational changes

## Prerequisites

- Trajectory file (.xtc/.trr) - PBC correction optional
- Topology file (.tpr)
- Index file (.ndx) with appropriate atom groups

## Workflow

### Step 1: Calculate Distance Matrix

Compute average distance matrix between residues:

```bash
echo -e "Protein\n" | gmx mdmat -f md.xtc -s md.tpr -mean rdcm.xpm
```

- Input: Select Protein group for analysis

**Parameters**:
- `-f`: Input trajectory file
- `-s`: Input topology file
- `-mean`: Calculate average distance matrix
- Optional: `-b` and `-e` to specify time range

**Output file**:
- `rdcm.xpm`: Average distance matrix in XPM format

### Step 2: Visualize Distance Matrix

Generate heatmap visualization of the distance matrix:

```bash
dit xpm_show -f rdcm.xpm -o rdcm.png
```

**Optional parameters**:
- `-cmap viridis`: Use viridis colormap (default works well)
- `-xmin`, `-xmax`, `-ymin`, `-ymax`: Zoom on specific regions
- `-m contour`: Use contour plot mode

### Step 3: Analyze Distance Patterns

Examine the distance matrix for:

- **Diagonal**: Short distances (residues close in sequence)
- **Off-diagonal patterns**: Long-range contacts between distant residues
- **Blue regions**: Short distances (close contacts)
- **Red/yellow regions**: Long distances (no contact)
- **Changes over time**: Compare matrices from different time windows

## Output Files

- **rdcm.xpm**: Average distance matrix (XPM format)
- **rdcm.png**: Visualization of distance matrix

## Interpretation Guidelines

### Distance Scale

The color scale typically represents distances in nanometers:
- **Dark blue**: < 0.5 nm (close contact)
- **Blue**: 0.5-1.0 nm (medium distance)
- **Green/Yellow**: 1.0-2.0 nm (no direct contact)
- **Red**: > 2.0 nm (distant)

### Contact Definition

Common contact definitions:
- **< 0.45 nm**: Direct atomic contact
- **< 0.6 nm**: Short-range interaction
- **< 1.0 nm**: Medium-range interaction

### Structural Features

- **Secondary structure**: Shows patterns for α-helices and β-sheets
- **Domain contacts**: Off-diagonal blocks indicate domain-domain interactions
- **Long-range contacts**: Contacts between distant sequence positions indicate tertiary structure
- **Flexible regions**: Variable distances indicate flexible loops or unstructured regions

## Common Issues and Solutions

### Issue: Distance matrix shows uniform high distances

**Possible causes**:
- Protein is unfolded or disordered
- Incorrect atom selection
- Simulation artifacts

**Solutions**:
- Verify protein structure in simulation
- Check index file selections
- Inspect trajectory for artifacts

### Issue: Distance matrix appears noisy

**Solutions**:
- Increase simulation time for better averaging
- Apply time averaging over windows
- Use larger time steps for analysis

### Issue: Unexpected distance patterns

**Solutions**:
- Verify protein sequence and structure
- Check for simulation instability
- Compare with known structures

## Tips and Best Practices

- **Time selection**: Use consistent time ranges for reproducibility
- **Averaging**: Longer simulations give better average distance matrices
- **Comparison**: Compare with experimental contact maps (NMR, cryo-EM)
- **Domain analysis**: Create domain-specific index groups for detailed analysis
- **Time evolution**: Calculate matrices for different time windows to study dynamics

## Advanced Analysis

### Time-dependent RDCM

Calculate distance matrices for different time windows:

```bash
# Early simulation (0-50 ns)
echo -e "Protein\n" | gmx mdmat -f md.xtc -s md.tpr -b 0 -e 50000 -mean rdcm_early.xpm

# Late simulation (50-100 ns)
echo -e "Protein\n" | gmx mdmat -f md.xtc -s md.tpr -b 50000 -e 100000 -mean rdcm_late.xpm
```

Compare matrices to identify conformational changes.

### Distance Distributions

For specific residue pairs, analyze distance distributions:

```bash
# Create index file with specific residues
echo -e "r 50\nname 10 Res50\nr 100\nname 11 Res100\nq\n" | gmx make_ndx -f md.tpr -o pair.ndx

# Calculate distance over time
echo -e "Res50\nRes100\n" | gmx distance -f md.xtc -s md.tpr -n pair.ndx -oall dist_50_100.xvg
```

### Contact Maps

Convert distance matrix to binary contact map:

```bash
# Use DuIvyTools to threshold distances
# (Requires custom script or manual processing)
```

Define contacts as residue pairs with average distance < threshold (typically 0.45-0.6 nm).

### Domain-specific Analysis

Calculate distance matrices for specific domains:

```bash
# Create domain index file
echo -e "r 1-100\nname 10 Domain1\nr 101-200\nname 11 Domain2\nq\n" | gmx make_ndx -f md.tpr -o domains.ndx

# Calculate inter-domain distances
echo -e "Domain1\nDomain2\n" | gmx mdmat -f md.xtc -s md.tpr -mean domain_distance.xpm
```

## Related Analyses

- **DCCM**: Analyzes correlated motions alongside distances
- **RMSF**: Provides per-residue flexibility information
- **Hydrogen bond analysis**: Identifies specific interactions
- **FEL**: Maps conformational states based on contact patterns

## Visualization Enhancements

### Custom Color Schemes

Use different colormaps to highlight features:

```bash
# Viridis (default, perceptually uniform)
dit xpm_show -f rdcm.xpm -o rdcm_viridis.png -cmap viridis

# Plasma (alternative perceptually uniform)
dit xpm_show -f rdcm.xpm -o rdcm_plasma.png -cmap plasma

# Coolwarm (diverging, good for differences)
dit xpm_show -f rdcm.xpm -o rdcm_coolwarm.png -cmap coolwarm
```

### Zooming on Regions

Focus on specific regions:

```bash
# Zoom on residues 50-150
dit xpm_show -f rdcm.xpm -o rdcm_zoom.png -xmin 50 -xmax 150 -ymin 50 -ymax 150
```

### Difference Maps

Compare distance matrices from different conditions:

```bash
# Calculate difference between two matrices
dit xpm_diff -f rdcm_early.xpm rdcm_late.xpm -o rdcm_diff.xpm
```

### Thresholding

Apply distance threshold to create contact maps:

```bash
# Convert to CSV and threshold
dit xpm2csv -f rdcm.xpm -o rdcm.csv
# Then process CSV to create binary contact map
```

## Interpretation Examples

### Alpha-Helical Proteins

Distance matrix shows characteristic pattern:
- Strong diagonal signal (residues i, i+3, i+4 in helix)
- Periodic pattern of short distances

### Beta-Sheet Proteins

Distance matrix shows:
- Strong off-diagonal signals for beta-strand pairs
- Anti-parallel or parallel strand patterns

### Multi-domain Proteins

Distance matrix shows:
- Clear separation between domains
- Inter-domain contacts in off-diagonal regions
- Variable distances indicating domain flexibility

## Validation and Comparison

### Compare with Experimental Data

- **NMR**: Compare with NOE distance constraints
- **Cryo-EM**: Compare with inter-residue distances in structure
- **Cross-linking**: Compare with cross-link distance constraints

### Compare with Static Structures

Calculate distance matrix from crystal structure:

```bash
# Convert PDB to gro
gmx editconf -f structure.pdb -o structure.gro

# Calculate distance matrix
echo -e "Protein\n" | gmx mdmat -f structure.gro -s structure.gro -mean rdcr_static.xpm
```

Compare with simulation-derived matrix to validate simulation.

## References

For theoretical background, see:
- Contact analysis in molecular dynamics simulations
- Protein folding and contact map analysis
- Domain organization and long-range contacts