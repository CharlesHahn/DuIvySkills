---
name: duivytools-skills
description: "CLI tool for visualizing and analyzing GROMACS MD simulation results. Use when Claude needs to: (1) Visualize XVG data (RMSD, RMSF, energy, gyrate, hydrogen bonds), (2) Display XPM matrices (DCCM, FEL, DSSP), (3) Manipulate NDX index files, (4) Calculate statistics or export data, (5) Generate Ramachandran plots"
---

# DuIvyTools (DIT)

> **CRITICAL: Always check command help before using**
>
> DuIvyTools is actively developed. **ALWAYS** run `dit <command> -h` to get the most accurate and up-to-date information for any command before using it.

DuIvyTools is a command-line tool for analyzing and visualizing GROMACS molecular dynamics simulation results. It provides ~30 commands for processing XVG, XPM, NDX, and other common GROMACS output files.

## Quick Start

### Installation

```bash
pip install DuIvyTools
```

Or use Tsinghua mirror (faster in China):
```bash
pip install DuIvyTools -i https://pypi.tuna.tsinghua.edu.cn/simple
```

### Verify Installation

```bash
dit
```

### Get Help

```bash
dit                    # All available commands
dit -h                 # Global parameters
dit xvg_show -h        # Specific command help (ALWAYS do this first!)
```

## Core Concepts

### Units and Data Conversion

- **GROMACS default units**: Time = ps, Distance = nm, Energy = kJ/mol
- **DO NOT convert units by default**: Preserve original units unless user explicitly requests conversion
- **When converting**: Update both data AND axis labels consistently

### Column and Row Indexing

- **Columns**: 0-based indexing (0 = first column)
- **Rows**: 0-based indexing (0 = first row)
- **Atom indices (NDX)**: 1-based (GROMACS convention)
- **Ranges**: Left-closed, right-open (e.g., `1-7` means columns 1,2,3,4,5,6)

### XPM Data Types

- **Discrete** (e.g., DSSP): Each pixel represents a discrete state
- **Continuous** (e.g., DCCM, FEL): Each pixel represents a continuous value
- **xpm_diff**: Only for Continuous types, NOT for Discrete types like DSSP

## Commands

### XVG Commands (Data Visualization)

| Command | Description |
|---------|-------------|
| `xvg_show` | Display XVG data with automatic legend/axis parsing |
| `xvg_compare` | Compare multiple XVG files with flexible column selection |
| `xvg_ave` | Calculate average, std dev, std err for each column |
| `xvg_show_distribution` | Show data distribution as histogram/PDF/CDF |
| `xvg_show_stack` | Display stacked area plots |
| `xvg_show_scatter` | Create 2D or 3D scatter plots |
| `xvg_box_compare` | Compare distributions using violin and scatter plots |
| `xvg_combine` | Combine data from multiple XVG files into one |
| `xvg_ave_bar` | Calculate and display averages across parallel experiments |
| `xvg_rama` | Draw Ramachandran plots from phi/psi angle data |
| `xvg_energy_compute` | Calculate protein-ligand interaction energies |

### XPM Commands (Matrix Visualization)

| Command | Description |
|---------|-------------|
| `xpm_show` | Visualize XPM matrices with 4 plotting engines |
| `xpm2csv` | Convert XPM to CSV format (x, y, value) |
| `xpm2dat` | Convert XPM to M×N matrix format |
| `xpm_diff` | Calculate difference between two XPM files (Continuous only) |
| `xpm_merge` | Merge two XPM files diagonally |

### NDX Commands (Index Files)

| Command | Description |
|---------|-------------|
| `ndx_add` | Add new index groups to NDX file |
| `ndx_split` | Split one index group into multiple groups |
| `ndx_show` | Display all index group names |
| `ndx_rm_dup` | Remove duplicate index groups |

### Utility Commands

| Command | Description |
|---------|-------------|
| `mdp_gen` | Generate GROMACS MDP file templates |
| `show_style` | Display/generate plotting style configuration files |
| `find_center` | Find geometric center of atom groups |
| `dccm_ascii` | Convert ASCII covariance matrix to DCCM XPM |
| `dssp` | Convert GROMACS 2023 DSSP format to XPM/XVG |

For complete parameter documentation, see [commands-reference.md](references/commands-reference.md).

## Plotting Engines

| Engine | Use Case |
|--------|----------|
| `matplotlib` (default) | Most comprehensive, all modes supported |
| `plotly` | Interactive plots, web integration |
| `gnuplot` | High-quality output, requires installation |
| `plotext` | Terminal-based, quick preview |

Use `-eg` parameter to specify engine: `dit xvg_show -f data.xvg -eg plotly`

## Common Workflows

### Visualize XVG Data

```bash
# Basic display
dit xvg_show -f rmsd.xvg

# Compare multiple files
dit xvg_compare -f rmsd.xvg gyrate.xvg -c 1 1,2 -l RMSD Gyrate

# Calculate average (from row 2000 to end)
dit xvg_ave -f rmsd.xvg -b 2000

# Show distribution
dit xvg_show_distribution -f gyrate.xvg -c 1
```

### Display XPM Matrix

```bash
# Basic display
dit xpm_show -f dccm.xpm -cmap coolwarm -zmin -1 -zmax 1

# 3D mode
dit xpm_show -f fel.xpm -m 3d -x PC1 -y PC2 -z Energy

# Contour plot
dit xpm_show -f dccm.xpm -m contour -cmap coolwarm

# Interactive plotly
dit xpm_show -f fel.xpm -eg plotly -m 3d
```

For more real-world scenarios, see [examples.md](references/examples.md).

## Reference Documentation

- **[commands-reference.md](references/commands-reference.md)** - Complete command parameters
- **[examples.md](references/examples.md)** - Real-world usage scenarios

## Resources

- GitHub: https://github.com/CharlesHahn/DuIvyTools
- Documentation: https://duivytools.readthedocs.io/
- Citation: https://doi.org/10.5281/zenodo.6339993
