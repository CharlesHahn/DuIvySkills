---
name: gromacs-protein-analysis
description: "Comprehensive guide for analyzing and visualizing protein molecular dynamics simulation results from GROMACS. Use when Claude needs to perform trajectory analysis including: (1) Periodic boundary condition correction for protein or protein-ligand complexes, (2) RMSD analysis to measure structural stability and convergence, (3) RMSF analysis to evaluate per-residue flexibility, (4) Gyrate analysis to assess protein compactness and folding, (5) SASA analysis to study solvent accessibility and surface properties, (6) Dynamics Cross-Correlation Matrix (DCCM) analysis to study correlated atomic motions, (7) Residue distance contact matrix (RDCM) to analyze inter-residue contacts, (8) Principal Component Analysis (PCA) to identify collective motions, (9) Free Energy Landscape (FEL) mapping using RMSD/Gyrate or PCA to understand conformational states"
---

# GROMACS Protein Analysis

This skill provides comprehensive workflows for analyzing protein molecular dynamics simulation results from GROMACS. It covers five major analysis types commonly used in protein dynamics studies.

## Prerequisites

- GROMACS simulation completed with trajectory files (.xtc/.trr) and topology file (.tpr)
- Understanding of GROMACS commands (call `gromacs-skills` when needed)
- For visualization: DuIvyTools skill (call `duivytools-skills` when needed)

## Analysis Types

### 1. Periodic Boundary Condition (PBC) Correction

Correct trajectory for PBC artifacts to prevent molecules from crossing box boundaries and ensure proper alignment for downstream analysis.

**Purpose**: Remove PBC artifacts, center protein/ligand, eliminate overall translation/rotation

**Input**: Trajectory file (.xtc), topology file (.tpr), index file (.ndx)

**Output**: Corrected trajectory (fit.xtc), corrected topology (fit.tpr)

**When to use**: Before any analysis when molecules cross box boundaries or when RMSD shows abrupt jumps

**Detailed workflow**: See [PBC Correction Guide](references/pbc-correction.md)

### 2. Root Mean Square Deviation (RMSD)

Calculate RMSD to measure structural stability and assess simulation convergence.

**Purpose**: Monitor structural stability, identify equilibration phase, assess simulation convergence, compare structures

**Input**: Trajectory file (.xtc), topology file (.tpr), reference structure

**Output**: RMSD time series (rmsd.xvg), visualization (rmsd.png)

**Visualization**: Use `duivytools-skills` skill to plot RMSD over time

**Detailed workflow**: See [RMSD Analysis Guide](references/rmsd-analysis.md)

### 3. Root Mean Square Fluctuation (RMSF)

Calculate RMSF to evaluate per-residue flexibility and identify flexible/rigid regions.

**Purpose**: Identify flexible regions, assess local stability, analyze loop dynamics, compare residue flexibility

**Input**: Trajectory file (.xtc), topology file (.tpr), index file (.ndx)

**Output**: RMSF per residue (rmsf.xvg), B-factors (bfactor.pdb), visualization (rmsf.png)

**Visualization**: Use `duivytools-skills` skill to plot RMSF per residue

**Detailed workflow**: See [RMSF Analysis Guide](references/rmsf-analysis.md)

### 4. Radius of Gyration (Gyrate)

Calculate radius of gyration to assess protein compactness and folding state.

**Purpose**: Monitor protein compactness, detect unfolding/folding transitions, assess overall size changes

**Input**: Trajectory file (.xtc), topology file (.tpr), index file (.ndx)

**Output**: Gyrate time series (gyrate.xvg), per-axis data (gyrate_axes.xvg), visualization (gyrate.png)

**Visualization**: Use `duivytools-skills` skill to plot gyrate over time

**Detailed workflow**: See [Gyrate Analysis Guide](references/gyrate-analysis.md)

### 5. Solvent Accessible Surface Area (SASA)

Calculate SASA to study solvent accessibility and surface properties of the protein.

**Purpose**: Analyze solvent exposure, identify hydrophobic/hydrophilic surfaces, study ligand binding sites, monitor protein unfolding

**Input**: Trajectory file (.xtc), topology file (.tpr), index file (.ndx)

**Output**: Total SASA time series (sas.xvg), per-residue SASA (sas_per_residue.xvg), visualization (sas.png)

**Visualization**: Use `duivytools-skills` skill to plot SASA over time

**Detailed workflow**: See [SASA Analysis Guide](references/sasa-analysis.md)

### 6. Dynamics Cross-Correlation Matrix (DCCM)

Analyze correlated motions between atomic pairs to identify coordinated movements in the protein.

**Purpose**: Identify pairs of atoms that move together (positive correlation) or opposite (negative correlation)

**Input**: Trajectory file (.xtc), topology file (.tpr), index file (.ndx)

**Output**: Covariance matrix (covar.dat), DCCM matrix (dccm.xpm), visualization (dccm.png)

**Visualization**: Use `duivytools-skills` skill to generate heatmaps

**Detailed workflow**: See [DCCM Analysis Guide](references/dccm-analysis.md)

### 7. Residue Distance Contact Matrix (RDCM)

Calculate average distances between residue pairs to analyze inter-residue contacts and spatial relationships.

**Purpose**: Map residue-residue distances, identify long-range contacts, analyze protein structure

**Input**: Trajectory file (.xtc), topology file (.tpr), index file (.ndx)

**Output**: Distance contact matrix (rdcm.xpm), visualization (rdcm.png)

**Visualization**: Use `duivytools-skills` skill to generate heatmaps

**Detailed workflow**: See [RDCM Analysis Guide](references/rdcm-analysis.md)

### 8. Principal Component Analysis (PCA)

Identify collective motions and dominant conformational changes by decomposing protein motion into principal components.

**Purpose**: Extract major collective motions, analyze conformational flexibility, reduce dimensionality

**Input**: Trajectory file (.xtc), topology file (.tpr), index file (.ndx)

**Output**: Eigenvalues (eigenvalues.xvg), eigenvectors (eigenvectors.trr), projections (pc1.xvg, pc2.xvg)

**Key metrics**: First few PCs' contribution percentages

**Visualization**: Use `duivytools-skills` skill to plot eigenvalues and projections

**Detailed workflow**: See [PCA Analysis Guide](references/pca-analysis.md)

### 9. Free Energy Landscape (FEL)

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

### Analysis Independence

Most analyses are independent and can be performed in any order:

**Independent analyses** (no dependencies):
- RMSD, RMSF, Gyrate, SASA - Basic stability and property analysis
- DCCM, RDCM - Contact and correlation analysis
- PCA - Collective motion analysis

**Dependent analyses** (require other analyses):
- FEL - Requires RMSD/Gyrate OR PCA results as reaction coordinates

### Recommended Pre-Analysis Step

**PBC Correction** (optional but recommended):
- Consider PBC correction if RMSD shows abrupt jumps or molecules cross box boundaries
- PBC correction is not always necessary - only apply when indicated by trajectory quality
- Analysis issues may have other causes beyond PBC artifacts

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

### Basic Analysis

- **[PBC Correction Guide](references/pbc-correction.md)** - Complete PBC correction workflow
- **[RMSD Analysis Guide](references/rmsd-analysis.md)** - Root mean square deviation for stability assessment
- **[RMSF Analysis Guide](references/rmsf-analysis.md)** - Root mean square fluctuation for flexibility analysis
- **[Gyrate Analysis Guide](references/gyrate-analysis.md)** - Radius of gyration for compactness analysis
- **[SASA Analysis Guide](references/sasa-analysis.md)** - Solvent accessible surface area for surface properties

### Advanced Analysis

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

Analyses can be performed flexibly based on your research questions:

**For basic stability assessment**: Start with RMSD, RMSF, Gyrate, SASA (any order)

**For contact and correlation analysis**: Perform DCCM and RDCM (independent)

**For collective motion analysis**: Perform PCA (independent)

**For conformational landscape**: Generate FEL after obtaining RMSD/Gyrate or PCA data

**Note**: Consider PBC correction only if RMSD shows abrupt jumps or trajectory quality issues or User asks. Many analyses work fine without PBC correction.

## Best Practices

- **Never overwrite** existing files - use unique output filenames
- **Use consistent time ranges** for all related analyses
- **Document all parameters** and index group selections
- **Visualize intermediate results** to catch issues early
- **Verify atom count consistency** between tpr and xtc files
- **Check statistical convergence** before interpreting results
