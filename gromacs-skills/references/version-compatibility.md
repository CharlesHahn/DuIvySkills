# GROMACS Version Compatibility

This document covers version differences, migration guides, and compatibility issues across GROMACS versions.

## Version Detection

### Check GROMACS Version

```bash
gmx --version
```

**Example output**:
```
GROMACS version:    2024.3
GROMACS modification: 2024.3-1
Precision: mixed
MPI: enabled
OpenMP: enabled
GPU support: enabled
FFT library: fftw-MPI
```

### Version Numbering

GROMACS uses year-based versioning:
- **YYYY.R**: Year.Release (e.g., 2024.3)
- **Old format**: Major.Minor (e.g., 5.1.4, 4.6.7)

## Major Version Changes

### GROMACS 2024.x (Current)

**Key changes**:
- Enhanced GPU support and performance
- Improved PME and bonded interactions on GPU
- New selection language features
- Enhanced trajectory compression
- Improved analysis tools

**Compatibility**:
- Topology files compatible with 2023.x
- Trajectory files compatible with 2023.x
- Some .mdp parameters added/changed

**Migration from 2023.x**:
- Generally compatible, no major changes needed
- Check new GPU parameters for optimization

### GROMACS 2023.x

**Key changes**:
- New selection language (more powerful)
- Enhanced DSSP integration
- Improved energy output formats
- New analysis commands

**Compatibility**:
- Topology files mostly compatible with 2022.x
- Trajectory files compatible
- Selection syntax changes (old syntax still supported)

**Migration from 2022.x**:
- Selection language updates (check complex selections)
- New analysis features available

### GROMACS 2022.x

**Key changes**:
- GPU-accelerated PME
- Improved domain decomposition
- New force field support
- Enhanced hydrogen bond analysis

**Compatibility**:
- Topology files compatible with 2021.x
- Trajectory files compatible
- Some .mdp parameter defaults changed

**Migration from 2021.x**:
- Check .mdp parameter defaults
- GPU parameters may need adjustment

### GROMACS 2021.x

**Key changes**:
- Major GPU acceleration improvements
- New neighbor searching algorithm
- Enhanced replica exchange
- Improved energy conservation

**Compatibility**:
- Topology files compatible with 2020.x
- Trajectory files compatible
- Performance parameters changed

**Migration from 2020.x**:
- Update GPU parameters for better performance
- Check energy conservation

### GROMACS 2020.x

**Key changes**:
- New Verlet cutoff scheme (default)
- Improved PME performance
- New pull code features
- Enhanced cluster analysis

**Compatibility**:
- Topology files compatible with 2019.x
- Trajectory files compatible
- Verlet scheme requires .mdp changes

**Migration from 2019.x**:
- Switch to Verlet cutoff scheme (recommended)
- Update .mdp parameters for Verlet scheme

### GROMACS 2019.x

**Key changes**:
- Improved GPU support
- New energy minimization algorithms
- Enhanced trajectory compression
- New analysis features

**Compatibility**:
- Topology files compatible with 2018.x
- Trajectory files compatible
- Some analysis parameters changed

**Migration from 2018.x**:
- Check analysis command parameters
- Update GPU settings

### GROMACS 2018.x

**Key changes**:
- Major GPU acceleration
- New PME GPU implementation
- Improved domain decomposition
- Enhanced selection syntax

**Compatibility**:
- Topology files compatible with 2016.x
- Trajectory files compatible
- GPU parameters changed significantly

**Migration from 2016.x**:
- GPU support requires new parameters
- Check domain decomposition settings

### GROMACS 2016.x

**Key changes**:
- New neighbor searching (Verlet)
- Improved PME performance
- New energy minimization
- Enhanced analysis tools

**Compatibility**:
- Topology files compatible with 5.x
- Trajectory files compatible
- Some .mdp parameters deprecated

**Migration from 5.x**:
- Update deprecated .mdp parameters
- Use new neighbor searching

### GROMACS 5.x

**Key changes from 4.x**:
- New topology format
- New .mdp file format
- Improved energy conservation
- New analysis commands

**Compatibility**:
- **Not compatible** with 4.x topology files
- Trajectory files compatible
- Requires topology regeneration

**Migration from 4.x**:
- Regenerate topology with `pdb2gmx`
- Update .mdp files
- Check force field compatibility

## Command Changes by Version

### 2024.x New Commands
- Enhanced selection features
- New GPU options
- Improved analysis tools

### 2023.x New Commands
- `gmx select` with enhanced features
- New `do_dssp` integration
- Enhanced cluster analysis

### 2022.x New Commands
- New GPU acceleration options
- Enhanced pull code
- New energy analysis features

### 2021.x New Commands
- Enhanced replica exchange
- New energy minimization options
- Improved trajectory processing

### 2020.x New Commands
- New Verlet scheme options
- Enhanced PME features
- New cluster methods

### 2019.x New Commands
- Enhanced hydrogen bond analysis
- New energy extraction options
- Improved trajectory compression

### 2018.x New Commands
- New GPU parameters
- Enhanced domain decomposition
- New selection features

### 2016.x New Commands
- New neighbor searching options
- Enhanced energy minimization
- New analysis features

### 5.x New Commands (compared to 4.x)
- `gmx energy` (enhanced)
- `gmx rms` (new options)
- `gmx rmsf` (new options)
- Many analysis commands enhanced

## Parameter Changes

### .mdp File Parameters

#### Verlet Cutoff Scheme (2020.x+)

**Old (group scheme)**:
```
cutoff-scheme = Group
ns_type = grid
nstlist = 10
rlist = 1.0
```

**New (Verlet scheme)**:
```
cutoff-scheme = Verlet
ns_type = grid
nstlist = 10
rlist = 1.0
rcoulomb = 1.0
rvdw = 1.0
```

#### GPU Parameters (2018.x+)

**New parameters**:
```
nb.gpu = 0
pme.gpu = 0
update.gpu = 0
bonded.gpu = 0
```

#### Pull Code Parameters

**Newer versions**:
```
pull-coord1-type = umbrella
pull-coord1-geometry = distance
pull-coord1-dim = Y Y Y
pull-coord1-k = 1000
pull-coord1-groups = 1 2
```

### Trajectory Parameters

#### Compression Options

**Newer versions support**:
```
xtc-precision = 1000
compressed-x-grps = Protein
```

#### Output Frequency

**Consistent across versions**:
```
nstxout = 0
nstvout = 0
nstenergy = 1000
nstlog = 1000
```

## Topology File Compatibility

### Format Changes

#### GROMACS 5.x vs 4.x

**5.x topology format**:
```
[ defaults ]
; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ
1               2               no              0.5     0.8333

[ atomtypes ]
;name  at.num mass charge ptype sigma epsilon
```

**4.x topology format** (deprecated):
```
[ defaults ]
; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ
1               2               yes             0.5     0.8333
```

### Force Field Compatibility

**Compatible force fields**:
- AMBER: Compatible across versions
- CHARMM: Compatible across versions
- OPLS-AA: Compatible across versions
- GROMOS: Compatible across versions
- MARTINI: Compatible across versions

**Note**: Some force field parameters may update with new GROMACS versions

## Trajectory File Compatibility

### Format Compatibility

**.xtc files**:
- Generally compatible across versions
- Compression algorithms may differ
- Use `trjconv` for format conversion if needed

**.trr files**:
- Compatible across versions
- Full precision format

**.edr files**:
- Generally compatible
- Some energy terms may differ
- Use `eneconv` for conversion if needed

**.tpr files**:
- **NOT compatible** across major versions
- Must regenerate with `grompp` for each version

## Common Compatibility Issues

### Issue 1: "Fatal error: Cannot read .tpr file"

**Cause**: .tpr file from different GROMACS version

**Solution**:
```bash
# Regenerate .tpr file
gmx grompp -f mdp/md.mdp -c structure.gro -p topol.top -o md.tpr
```

### Issue 2: "Fatal error: Unknown mdp option"

**Cause**: Using .mdp file from different version

**Solution**:
- Check .mdp parameters for current version
- Use `gmx grompp -h` to see valid options
- Update deprecated parameters

### Issue 3: "Fatal error: Domain decomposition failed"

**Cause**: Domain decomposition parameters incompatible

**Solution**:
```bash
# Try fewer MPI ranks
mpirun -np 2 gmx_mpi mdrun -s md.tpr -deffnm md

# Or use different domain decomposition
mpirun -np 4 gmx_mpi mdrun -s md.tpr -deffnm md -dd 2 2 2
```

### Issue 4: "Fatal error: GPU initialization failed"

**Cause**: GPU parameters incompatible

**Solution**:
```bash
# Check GPU availability
gmx mdrun -h | grep -i gpu

# Try CPU instead
gmx mdrun -s md.tpr -deffnm md -nb cpu -pme cpu

# Update GPU parameters for current version
```

### Issue 5: Energy drift in new version

**Cause**: Default parameter changes

**Solution**:
- Check .mdp parameters against old version
- Use same parameters as old simulation
- Verify energy conservation

## Version-Specific Features

### 2024.x Features
- Enhanced GPU performance
- New selection language features
- Improved trajectory compression
- Enhanced analysis tools

### 2023.x Features
- New selection language
- Enhanced DSSP integration
- Improved energy output
- New cluster methods

### 2022.x Features
- GPU-accelerated PME
- Improved domain decomposition
- New pull code features
- Enhanced hydrogen bond analysis

### 2021.x Features
- Major GPU improvements
- New neighbor searching
- Enhanced replica exchange
- Improved energy conservation

### 2020.x Features
- Verlet cutoff scheme (default)
- Improved PME performance
- New pull code
- Enhanced cluster analysis

### 2019.x Features
- Improved GPU support
- New energy minimization
- Enhanced trajectory compression
- New analysis features

### 2018.x Features
- Major GPU acceleration
- New PME GPU implementation
- Improved domain decomposition
- Enhanced selection syntax

## Migration Checklist

### When Upgrading GROMACS

1. **Check version compatibility**
   ```bash
   gmx --version
   ```

2. **Backup all simulation files**
   - .tpr files
   - .mdp files
   - Topology files

3. **Regenerate .tpr files**
   ```bash
   gmx grompp -f mdp/md.mdp -c structure.gro -p topol.top -o md.tpr
   ```

4. **Update .mdp parameters**
   - Check for deprecated parameters
   - Update GPU parameters
   - Verify cutoff scheme

5. **Test on small system**
   - Run short test simulation
   - Verify energy conservation
   - Check output files

6. **Update analysis scripts**
   - Check command parameters
   - Update selection syntax if needed
   - Verify output formats

7. **Document changes**
   - Keep record of version changes
   - Document parameter updates
   - Save old configurations

## Best Practices for Version Management

1. **Document GROMACS version** in simulation directories
2. **Keep .mdp files** with version notes
3. **Use version control** for simulation scripts
4. **Test upgrades** on small systems first
5. **Keep backups** of old version files
6. **Use consistent versions** within a project
7. **Check compatibility** before major upgrades

## Getting Help for Version-Specific Issues

### Local Help
```bash
# Get help for specific version
gmx <command> -h

# Check available options
gmx grompp -h
gmx mdrun -h
```

### Online Documentation
```bash
# Search for version-specific documentation
web_search: "site:manual.gromacs.org gmx <command> <version>"

# Examples
web_search: "site:manual.gromacs.org gmx mdrun 2024.3"
web_search: "site:manual.gromacs.org gmx grompp 2023"
web_search: "site:manual.gromacs.org mdp parameters 2022"
```

### Version Migration Guides
- Check GROMACS release notes
- Read migration guides for major version changes
- Consult GROMACS user forums

## Version History Summary

| Version | Year | Key Features | Compatibility |
|---------|------|--------------|---------------|
| 2024.x | 2024 | Enhanced GPU, new selection | Compatible with 2023.x |
| 2023.x | 2023 | New selection language, DSSP | Compatible with 2022.x |
| 2022.x | 2022 | GPU PME, improved DD | Compatible with 2021.x |
| 2021.x | 2021 | Major GPU improvements | Compatible with 2020.x |
| 2020.x | 2020 | Verlet scheme, new PME | Compatible with 2019.x |
| 2019.x | 2019 | Improved GPU, new analysis | Compatible with 2018.x |
| 2018.x | 2018 | Major GPU acceleration | Compatible with 2016.x |
| 2016.x | 2016 | Verlet neighbor search | Compatible with 5.x |
| 5.x | 2012-2016 | New topology format | Not compatible with 4.x |
| 4.x | 2008-2012 | Legacy format | Deprecated |

## Important Notes

- **Always check version** before running simulations
- **Regenerate .tpr files** when changing versions
- **Update .mdp parameters** for new features
- **Test compatibility** on small systems
- **Document version** in simulation directories
- **Use local help** (`gmx <command> -h`) for accurate parameters
- **Search online** with version number for specific issues