# Periodic Boundary Condition (PBC) Correction

## Overview

Periodic boundary condition (PBC) correction removes artifacts from molecules crossing simulation box boundaries. This correction may be necessary when trajectory quality issues are observed (e.g., abrupt RMSD jumps). Many analyses work fine without PBC correction.

## When to Apply PBC Correction

- Molecules (protein or ligand) cross box boundaries
- RMSD/RMSF plots show abrupt jumps or discontinuities
- Visualization shows fragmented molecules

Note: PBC correction is optional. Only apply when you observe the above issues. Analysis problems may have other causes beyond PBC artifacts.

## Prerequisites

- Input trajectory file (.xtc)
- Topology file (.tpr)
- Index file (.ndx) with appropriate atom groups

## Workflow

### Step 1: Prepare Index File

Create or modify index file to include necessary groups:

```bash
# Create basic index file
gmx make_ndx -f md.tpr -o index.ndx

# Add protein-ligand group if needed
echo -e "Protein | Ligand\nq\n" | gmx make_ndx -f md.tpr -o index.ndx
```

### Step 2: Find Center Atom

Find atom closest to geometric center for proper centering:

```bash
# Use DuIvyTools to find center atom
echo -e "Protein\n" | dit find_center -f npt.gro index.ndx
```

The output shows the atom number closest to the center. Note this atom number.

### Step 3: Add Center Group to Index File

Manually add a center group to the index file with the identified atom:

```bash
# Add center group to index.ndx
echo -e "\n[ center ]\n<atom_number>\n" >> index.ndx
```

Replace `<atom_number>` with the actual atom number from Step 2.

### Step 4: Center the Trajectory

Center the system using the center atom:

```bash
echo -e "center\nProtein_Lig\n" | gmx trjconv -s md.tpr -f md.xtc -o center.xtc -n index.ndx -pbc atom -center
```

- First input: Select center group
- Second input: Select output group (Protein_Lig)

### Step 5: Ensure Molecular Integrity

Keep molecules intact across periodic boundaries:

```bash
echo -e "Protein_Lig\n" | gmx trjconv -s md.tpr -f center.xtc -o mol.xtc -n index.ndx -pbc mol -ur compact
```

- Input: Select group to keep intact (Protein_Lig)

### Step 6: Remove Overall Translation and Rotation

Fit trajectory to remove overall motion (optional but recommended):

```bash
echo -e "Backbone\nProtein_Lig\n" | gmx trjconv -s md.tpr -f mol.xtc -o fit.xtc -n index.ndx -fit rot+trans
```

- First input: Select reference structure for fitting (Backbone)
- Second input: Select output group (Protein_Lig)

### Step 7: Update Topology File

Update topology file to match corrected trajectory atom count:

```bash
echo -e "Protein_Lig\n" | gmx convert-tpr -s md.tpr -o fit.tpr -n index.ndx
```

- Input: Select same group used in Step 6

### Step 8: Verify Correction

Visual inspection of corrected trajectory:

```bash
# Convert to PDB for visual inspection
echo -e "Protein_Lig\n" | gmx trjconv -s md.tpr -f fit.xtc -o fit.pdb -dt 1000 -n index.ndx
```

Use PyMOL or VMD to inspect the PDB file. Check that:
- Molecules are intact and not fragmented
- Protein/ligand remains centered
- No molecules cross box boundaries

## Output Files

- **center.xtc**: Centered trajectory
- **mol.xtc**: Trajectory with intact molecules
- **fit.xtc**: Final corrected trajectory (no translation/rotation)
- **fit.tpr**: Corrected topology file
- **fit.pdb**: Sample structure for visual inspection

## Common Issues and Solutions

### Issue: Molecules still fragmented after correction

**Solution**: Check that the correct atom group is selected in Step 5. Ensure the group includes all atoms that should remain intact.

### Issue: Trajectory still shows PBC artifacts

**Solution**: Verify the center atom selection in Step 2. The center atom should be close to the geometric center of the system.

### Issue: tpr and xtc atom count mismatch

**Solution**: Always run Step 7 to update topology file after any group selection changes.

### Issue: RMSD still shows jumps after correction

**Solution**: Check if the reference structure for RMSD calculation is appropriate. Consider using the average structure as reference.

## Tips and Best Practices

- **Use appropriate time intervals**: When converting to PDB for inspection, use `-dt` parameter to limit output size
- **Document atom selections**: Keep track of which atoms are in each index group
- **Consistent group selection**: Use the same groups for all related analyses
- **Visual inspection**: Always visually inspect corrected trajectory before downstream analysis
- **Never overwrite existing files**: Use unique output filenames to preserve original data

## Related Analyses

PBC correction may benefit analyses where RMSD shows abrupt jumps:
- RMSD/RMSF analysis
- DCCM analysis
- PCA analysis
- FEL analysis
- Any structural analysis requiring aligned trajectories

Note: PBC correction is optional and not always necessary. Apply only when trajectory quality issues are observed.

## Alternative Approaches

### Alternative 1: Single-step correction

Combine centering and molecular integrity in one step:

```bash
echo -e "center\nProtein_Lig\n" | gmx trjconv -s md.tpr -f md.xtc -o fit.xtc -n index.ndx -pbc mol -center
```

### Alternative 2: Use protein center

Use protein mass center instead of specific atom:

```bash
echo -e "Protein\nProtein_Lig\n" | gmx trjconv -s md.tpr -f md.xtc -o fit.xtc -n index.ndx -pbc mol -center
```

## Visualization

Use DuIvyTools skill to visualize corrected trajectory:

```bash
# Plot RMSD to verify correction
dit xvg_show -f rmsd_corrected.xvg -x "Time (ns)" -y "RMSD (nm)"
```

The corrected RMSD should show smooth progression without abrupt jumps.