---
name: gromacs-protein-analysis
description: "Comprehensive guide for analyzing and visualizing protein molecular dynamics simulation results from GROMACS. Use when Claude needs to perform trajectory analysis including: (1) Periodic boundary condition correction for protein or protein-ligand complexes, (2) Dynamics Cross-Correlation Matrix (DCCM) analysis to study correlated atomic motions, (3) Residue distance contact matrix (RDCM) to analyze inter-residue contacts, (4) Principal Component Analysis (PCA) to identify collective motions, (5) Free Energy Landscape (FEL) mapping using RMSD/Gyrate or PCA to understand conformational states"
---

# GROMACS Protein Analysis

This skill provides comprehensive workflows for analyzing protein molecular dynamics simulation results from GROMACS. It covers five major analysis types commonly used in protein dynamics studies.

## Prerequisites

- GROMACS simulation completed with trajectory files (.xtc/.trr) and topology file (.tpr)
- Basic understanding of GROMACS commands
- For visualization: DuIvyTools skill (call `duivytools-skills` when needed)

## Analysis Types

### 1. Periodic Boundary Condition (PBC) Correction

Correct trajectory for PBC artifacts to prevent molecules from crossing box boundaries and ensure proper alignment for downstream analysis.

**Purpose**: Remove PBC artifacts, center protein/ligand, eliminate overall translation/rotation

**Input**: Trajectory file (.xtc), topology file (.tpr), index file (.ndx)

**Output**: Corrected trajectory (fit.xtc), corrected topology (fit.tpr)

**When to use**: Before any analysis when molecules cross box boundaries or when RMSD shows abrupt jumps

**Detailed workflow**: See [PBC Correction Guide](references/pbc-correction.md)

### 2. Dynamics Cross-Correlation Matrix (DCCM)

Analyze correlated motions between atomic pairs to identify coordinated movements in the protein.

**Purpose**: Identify pairs of atoms that move together (positive correlation) or opposite (negative correlation)

**Input**: Trajectory file (.xtc), topology file (.tpr), index file (.ndx)

**Output**: Covariance matrix (covar.dat), DCCM matrix (dccm.xpm), visualization (dccm.png)

**Visualization**: Use `duivytools-skills` skill to generate heatmaps

**Detailed workflow**: See [DCCM Analysis Guide](references/dccm-analysis.md)

### 3. Residue Distance Contact Matrix (RDCM)

Calculate average distances between residue pairs to analyze inter-residue contacts and spatial relationships.

**Purpose**: Map residue-residue distances, identify long-range contacts, analyze protein structure

**Input**: Trajectory file (.xtc), topology file (.tpr), index file (.ndx)

**Output**: Distance contact matrix (rdcm.xpm), visualization (rdcm.png)

**Visualization**: Use `duivytools-skills` skill to generate heatmaps

**Detailed workflow**: See [RDCM Analysis Guide](references/rdcm-analysis.md)

### 4. Principal Component Analysis (PCA)

Identify collective motions and dominant conformational changes by decomposing protein motion into principal components.

**Purpose**: Extract major collective motions, analyze conformational flexibility, reduce dimensionality

**Input**: Trajectory file (.xtc), topology file (.tpr), index file (.ndx)

**Output**: Eigenvalues (eigenvalues.xvg), eigenvectors (eigenvectors.trr), projections (pc1.xvg, pc2.xvg)

**Key metrics**: First few PCs' contribution percentages

**Visualization**: Use `duivytools-skills` skill to plot eigenvalues and projections

**Detailed workflow**: See [PCA Analysis Guide](references/pca-analysis.md)

### 5. Free Energy Landscape (FEL)

Map free energy surfaces to understand conformational states and transitions using either RMSD/Gyrate or PCA.

**Purpose**: Identify stable conformations, quantify energy barriers, understand conformational landscape

**Method 1 - RMSD + Gyrate**: Use structural deviation and compactness as reaction coordinates

**Method 2 - PCA**: Use principal components as reaction coordinates

**Input**: RMSD data (rmsd.xvg), Gyrate data (gyrate.xvg) OR PC projections (pc1.xvg, pc2.xvg)

**Output**: Free energy landscape (gibbs.xpm), energy minima (gibbs.log), frame indices (bindex.ndx), visualization (fel.png)

**Visualization**: Use `duivytools-skills` skill to generate 2D/3D FEL plots

**Detailed workflow**: See [FEL Analysis Guide](references/fel-analysis.md)

## General Workflow

### Before Any Analysis

1. **Check trajectory quality**: Visual inspection of trajectory for artifacts
2. **PBC correction** (if needed): Use PBC correction workflow
3. **Ensure consistent atom selection**: Use same index groups for related analyses

### Analysis Sequence

Typical analysis sequence:

```
1. PBC correction (if needed)
   ↓
2. DCCM analysis
   ↓
3. RDCM analysis
   ↓
4. PCA analysis
   ↓
5. FEL analysis (using RMSD/Gyrate OR PCA)
```

### After Each Analysis

- **Verify output files**: Check that all expected files are generated
- **Visual inspection**: Use appropriate visualization to assess results
- **Documentation**: Record analysis parameters and observations

## Key Considerations

### Atom Selection

- **Backbone**: For overall protein motion and RMSD analysis
- **C-alpha**: For PCA and DCCM (reduces computational cost)
- **Protein**: For full protein analysis
- **Protein_Lig**: For protein-ligand complexes

### Time Selection

- **Equilibration phase**: Exclude initial equilibration period (typically first 10-20% of simulation)
- **Production phase**: Use production phase for analysis
- **Consistency**: Use same time range for related analyses

### Index Groups

- Create appropriate index groups using `gmx make_ndx`
- Ensure index groups match analysis requirements
- Document index group compositions

## When to Call DuIvyTools-Skills

Call the `duivytools-skills` skill for visualization tasks:

- **XVG files**: Plot RMSD, RMSF, energies, hydrogen bonds, gyrate
- **XPM files**: Visualize DCCM, RDCM, FEL matrices
- **Projections**: Plot PC1, PC2 projections
- **Statistical analysis**: Calculate averages, distributions

## Troubleshooting

### Common Issues

**RMSD shows abrupt jumps**: PBC artifacts - apply PBC correction

**DCCM values all near zero**: Check atom selection, ensure sufficient dynamics

**PCA shows uniform eigenvalues**: May indicate no dominant collective motion or excessive noise

**FEL shows unrealistic barriers**: Check time range selection, ensure sufficient sampling

**Files don't match**: Verify tpr and xtc have same atom count, use `gmx convert-tpr` if needed

## Reference Documentation

For detailed step-by-step workflows, consult these references:

- **[PBC Correction Guide](references/pbc-correction.md)** - Complete PBC correction workflow
- **[DCCM Analysis Guide](references/dccm-analysis.md)** - DCCM calculation and interpretation
- **[RDCM Analysis Guide](references/rdcm-analysis.md)** - Distance contact matrix analysis
- **[PCA Analysis Guide](references/pca-analysis.md)** - Principal component analysis workflow
- **[FEL Analysis Guide](references/fel-analysis.md)** - Free energy landscape mapping

## Quick Reference

### Required Input Files

- **Trajectory**: .xtc or .trr file from GROMACS simulation
- **Topology**: .tpr file (must match trajectory atom count)
- **Index**: .ndx file with custom atom groups

### Common Output Files

- **XVG**: Time-series data (RMSD, eigenvalues, projections)
- **XPM**: Matrix data (DCCM, RDCM, FEL)
- **TRR**: Vector data (eigenvectors)
- **PDB**: Structure files (average, extreme conformations)

### Analysis Order Recommendation

1. Correct PBC if needed
2. Perform DCCM and RDCM for contact analysis
3. Conduct PCA for collective motion analysis
4. Generate FEL for conformational landscape

## Best Practices

- **Always backup** original trajectory files before modifications
- **Use consistent time ranges** for all related analyses
- **Document all parameters** and index group selections
- **Visualize intermediate results** to catch issues early
- **Verify atom count consistency** between tpr and xtc files
- **Check statistical convergence** before interpreting results
