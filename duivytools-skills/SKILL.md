---
name: duivytools-skills
description: "Comprehensive tool for GROMACS molecular dynamics simulation analysis and visualization. Use when Claude needs to process, visualize, or analyze GROMACS result files including XVG data visualization (RMSD, RMSF, energy, hydrogen bonds, gyrate), XPM matrix visualization (DCCM, FEL, DSSP secondary structure), NDX index file manipulation, statistical analysis (averages, distributions, correlations), Ramachandran plots, and batch processing for multiple simulation results"
---

# DuIvyTools (DIT)

> **⚠️ CRITICAL: Always check command help before using**
>
> DuIvyTools is actively developed. **ALWAYS** run `dit <command> -h` to get the most accurate and up-to-date information for any command before using it. This ensures you have the correct parameters and options for your version.

DuIvyTools is a command-line tool for analyzing and visualizing GROMACS molecular dynamics simulation results. It provides ~30 commands for processing XVG, XPM, NDX, and other common GROMACS output files.

## Quick Start

### Installation

```bash
pip install DuIvyTools -i https://pypi.tuna.tsinghua.edu.cn/simple
```

Verify installation:
```bash
dit
```

### Basic Usage

**IMPORTANT: Always check command help before using**

DuIvyTools is actively developed and parameters may vary between versions. **ALWAYS** use `dit <command> -h` to get the most accurate and up-to-date information for any command before using it.

```bash
# Get all available commands
dit

# Get global parameters
dit -h

# Get help for specific command (ALWAYS do this first!)
dit <command> -h
```

For example:
```bash
dit xvg_show -h
dit xpm_show -h
dit ndx_add -h
```

## Common Commands

### XVG Files (Data Plots)

- `xvg_show` - Display XVG data (simplest, shows all columns)
- `xvg_compare` - Compare multiple XVG files with flexible column selection
- `xvg_ave` - Calculate average, std dev, std err for each column
- `xvg_show_distribution` - Show data distribution (histogram/PDF/CDF)
- `xvg_show_stack` - Stacked area plots (e.g., secondary structure content)
- `xvg_show_scatter` - Scatter plots (2D or 3D with color mapping)
- `xvg_box_compare` - Violin + scatter plots for distribution comparison
- `xvg_combine` - Combine data from multiple XVG files into one
- `xvg_ave_bar` - Calculate averages across parallel experiments
- `xvg_rama` - Ramachandran plots from phi/psi angles
- `xvg_energy_compute` - Calculate protein-ligand interaction energies

### XPM Files (Matrix/Heatmaps)

- `xpm_show` - Visualize XPM matrices (4 engines: matplotlib/plotly/gnuplot/plotext)
- `xpm2csv` - Convert XPM to CSV format (x, y, value)
- `xpm2dat` - Convert XPM to M×N matrix format
- `xpm_diff` - Calculate difference between two XPM files
- `xpm_merge` - Merge two XPM files diagonally

### NDX Files (Index Groups)

- `ndx_add` - Add new index groups to NDX file
- `ndx_split` - Split one index group into multiple groups
- `ndx_show` - Display all index group names
- `ndx_rm_dup` - Remove duplicate index groups

### Other Utilities

- `mdp_gen` - Generate GROMACS MDP file templates
- `show_style` - Display/generate plotting style configuration files
- `find_center` - Find geometric center of atom groups
- `dccm_ascii` - Convert ASCII covariance matrix to DCCM XPM
- `dssp` - Convert GROMACS 2023 DSSP format to XPM/XVG

## Key Concepts

### Column and Row Indexing

- **Columns**: 0-based indexing (0 = first column)
- **Rows**: 0-based indexing (0 = first row)
- **Atom indices (NDX)**: 1-based (GROMACS convention)
- **Ranges**: Left-closed, right-open (e.g., `1-7` means columns 1,2,3,4,5,6)

### Common Parameters

**Input/Output**:
- `-f INPUT [INPUT ...]` - Input files (space-separated groups, comma-separated within groups)
- `-o OUTPUT` - Output file name
- `-ns` - Don't display figure (useful for batch processing)

**Column Selection**:
- `-c COLUMNS [COLUMNS ...]` - Select columns (0-based)
  - Example: `-c 1-7,10 0,1,4` = first file: cols 1-6 + col 10; second file: cols 0,1,4

**Data Range**:
- `-b BEGIN` - Start row (inclusive)
- `-e END` - End row (exclusive)
- `-dt DT` - Step size (default 1)

**Labels**:
- `-l LEGENDS [LEGENDS ...]` - Legend labels (match number of selected columns)
- `-x XLABEL`, `-y YLABEL`, `-z ZLABEL` - Axis labels
- `-t TITLE` - Figure title
- Supports LaTeX: `"$\Delta G_{energy}$"`

**Data Transformation**:
- `-xs`, `-ys`, `-zs` - Multiply data by factor (default 1.0)
- `-xp`, `-yp`, `-zp` - Add offset (default 0.0)

**Plotting**:
- `-eg {matplotlib,plotly,gnuplot,plotext}` - Plotting engine (default: matplotlib)
- `-cmap COLORMAP` - Color map (matplotlib/plotly)
- `-smv [{,CI,origin}]` - Show moving average (with CI or origin background)
- `-ws WINDOWSIZE` - Moving average window size (default 50)
- `--alpha ALPHA` - Transparency
- `-csv CSV` - Export data to CSV

### Plotting Engines

- **matplotlib** (default): Most comprehensive, supports all modes and parameters
- **plotly**: Interactive plots, good for 3D and contour plots
- **gnuplot**: High-performance, requires gnuplot installation, generates scripts
- **plotext**: Terminal-based, for simple quick plots

### XPM Visualization Modes

- **imshow** (matplotlib default): Fast image display
- **pcolormesh**: Grid-based visualization
- **3d**: 3D surface plots
- **contour**: Contour plots

## Usage Patterns

### Analyze Single XVG File

```bash
# Simple display
dit xvg_show -f rmsd.xvg

# With labels
dit xvg_show -f rmsd.xvg -x "Time (ns)" -y "RMSD (nm)" -t "RMSD Analysis"

# Convert ps to ns
dit xvg_show -f rmsd.xvg -xs 0.001 -x "Time (ns)"
```

### Compare Multiple Files

```bash
# Compare two files, different columns
dit xvg_compare -f rmsd.xvg gyrate.xvg -c 1 1,2 \
  -l RMSD Gyrate Gyrate_X -x "Time (ns)"

# With moving average
dit xvg_compare -f energy.xvg -c 1,3 -l "LJ(SR)" "Coulomb(SR)" -smv

# Export to CSV
dit xvg_compare -f rmsd.xvg -c 1 -ns -csv rmsd_data.csv
```

### Visualize XPM Matrix

```bash
# Basic display
dit xpm_show -f dssp.xpm

# 3D mode
dit xpm_show -f fel.xpm -m 3d -x PC1 -y PC2 -z Energy -t FEL

# Contour plot
dit xpm_show -f dccm.xpm -m contour -cmap coolwarm -zmin -1 -zmax 1

# Crop region (pixels)
dit xpm_show -f dssp.xpm -xmin 100 -xmax 200 -ymin 500 -ymax 1000

# Interactive plotly
dit xpm_show -f fel.xpm -eg plotly -m 3d
```

### Statistical Analysis

```bash
# Calculate averages (rows 2000 to end)
dit xvg_ave -f rmsd.xvg -b 2000

# Show distribution
dit xvg_show_distribution -f gyrate.xvg -c 1,2

# PDF/CDF
dit xvg_show_distribution -f gyrate.xvg -c 1,2 -m pdf
```

### Manipulate NDX Files

```bash
# Show all groups
dit ndx_show -f index.ndx

# Add new group (residues 1-10)
dit ndx_add -f index.ndx -o test.ndx -al lig -c 1-10

# Split group into 2
dit ndx_split -f index.ndx -al Protein 2
```

### Generate MDP Template

```bash
dit mdp_gen -o nvt.mdp
```

## Advanced Topics

### Custom Plotting Styles

Create `dit_mplstyle.mplstyle` in working directory:

```python
axes.labelsize: 12
lines.linewidth: 2
figure.dpi: 100
savefig.dpi: 300
image.cmap: coolwarm
```

DIT automatically loads this file.

### Batch Processing

```bash
# Generate multiple plots
for file in md*_rmsd.xvg; do
  dit xvg_show -f "$file" -ns -o "${file%.xvg}.png"
done
```

### Interpolation

```bash
# Linear interpolation, 10x
dit xpm_show -f fel.xpm -ip linear -ipf 10

# Bicubic (matplotlib imshow)
dit xpm_show -f fel.xpm -m imshow -ip bicubic
```

## Important Notes

### Installation Requirements

- Python 3.6+
- matplotlib (default)
- numpy, pandas

### Optional Dependencies

- `pip install plotly` - For interactive plots
- Install gnuplot - For gnuplot engine
- `pip install plotext` - For terminal plots

### Column Indexing

- **Always 0-based** for XVG/XPM columns and rows
- **1-based** for NDX atom indices (GROMACS convention)
- `-e` parameter is exclusive (not included in calculation)

### Data Selection

```bash
# Select columns 1-6 and column 10
dit xvg_show -f rmsd.xvg -c 1-7,10

# Rows 100-199 (every 2nd row)
dit xvg_show -f rmsd.xvg -b 100 -e 200 -dt 2
```

### LaTeX Support

```bash
dit xvg_show -f energy.xvg -l "$\Delta G_{binding}$" -x "Time ($ns$)"
```

## When to Use This Skill

Use DuIvyTools when you need to:

1. **Visualize XVG data**: RMSD, RMSF, energy, hydrogen bonds, gyrate, etc.
2. **Analyze XPM matrices**: DCCM, FEL, DSSP secondary structure, contact maps
3. **Process NDX files**: Add/split/show index groups
4. **Calculate statistics**: Averages, distributions, correlations
5. **Generate plots**: Line, scatter, violin, stacked, 3D, contour, heatmaps
6. **Batch process**: Multiple simulation results with consistent styling
7. **Export data**: Convert to CSV/DAT for other software

## Reference Documentation

For detailed information, consult these references:

- **[Installation & Configuration](references/installation.md)** - Installation methods, plotting engine setup, custom styles
- **[Commands Reference](references/commands-reference.md)** - Complete documentation for all commands, parameters, and options
- **[Usage Examples](references/examples.md)** - Real-world scenarios and complete workflows

## Get Help

```bash
# All commands
dit

# Global parameters
dit -h

# Specific command
dit xvg_show -h
dit xpm_show -h
dit ndx_add -h
```

## Citation

If you use DuIvyTools in your research, please cite:

https://doi.org/10.5281/zenodo.6339993

## Resources

- GitHub: https://github.com/CharlesHahn/DuIvyTools
- Documentation: https://duivytools.readthedocs.io/
- PyPI: https://pypi.org/project/DuIvyTools/
