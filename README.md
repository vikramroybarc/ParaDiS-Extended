# ParaDiS

ParaDiS Public Release Version 4.0 — Local Enhanced Version

ParaDiS (Parallel Dislocation Simulator) is a simulation tool that performs direct numerical simulation of dislocation ensembles, the carriers of plasticity, to predict the strength of crystalline materials from the fundamental physics of defect motion, evolution, and interaction.

This repository is based on the LLNL ParaDiS public release and contains additional local developments extending the original code base. In particular, the current version includes:

1. **`MobilityLaw_BccnlEshlby`** — an additional BCC non-linear mobility-law implementation with Eshelby-related treatment.
2. **Drift-mode subcycling integrator** — an enhanced time-integration capability that introduces subcycling for drift-mode evolution.

These additions are local enhancements and are not part of the original LLNL ParaDiS 4.0 public release.

The ParaDiS code has been successfully deployed on high-performance computing architectures and used to study the origins of strength and strain hardening for cubic crystals, the strength of micro-pillars, and irradiated materials at LLNL. The original ParaDiS code has been deployed on more than one hundred thousand CPUs with over ten million active degrees of freedom.


## Installation

### Quick Start

The public distribution of ParaDiS is built using conventional makefiles.

1. Select a target system (for example, `SYS=gcc`) and compilation mode (`MODE=SERIAL` or `MODE=PARALLEL`) in `makefile.setup`.
2. Compile the code from the root of the distribution using the `make` command.
3. Successful compilation will produce binary executables in the `./bin` directory.
4. Test your installation by running an example:

```bash
./bin/paradis tests/frank_read_src.ctrl
```

### Detailed Instructions

The build system is comprised of a number of makefiles:

* `makefile`         : overall makefile to build ParaDiS and supporting utilities
* `makefile.setup`   : used to enable/disable various application features
* `makefile.sys`     : system-specific build settings for supported systems
* `makefile.srcs`    : complete list of source files used to build the system
* `src/makefile`     : makefile that builds ParaDiS
* `src/makefile.dep` : makefile dependencies, auto-generated via `make depend`
* `tests/makefile`   : makefile for cleaning and managing the test directory
* `utils/makefile`   : makefile for building the various utility applications
* `ext/makefile`     : makefile for building external dependencies

If you want or need to customize the build for your machine, you may need to adjust `makefile.sys` and `makefile.setup` before executing `make`.

Executing `make` from the main directory will build the main ParaDiS executable and all supporting tools. If you only want to build the main ParaDiS application, execute `make` from the main source directory (`src/`).

#### `makefile.sys`

Prior to building ParaDiS, you will need to identify the target system and environment for the application. Several target systems are preset and are listed below:

```text
linux.intel     Linux systems using native Intel compilers
linux           Generic Linux system
gcc             Generic system build using GNU compilers
aix             IBM AIX systems using native compilers
mac             macOS / MacBook Pro
bgp             LC BlueGene/P systems
bgq             LC BlueGene/Q systems
mc-cc           Stanford ME Linux system using Intel compilers
wcr             Stanford ME Linux system using Intel compilers
cygwin          Linux-compatible environment for Windows
xt4             Cray XT4 systems
```

If you are attempting to build ParaDiS on a system that is not listed above, you may need to copy and/or adjust an existing system configuration in `makefile.sys` and set the `SYS=` parameter in `makefile.setup`.

#### `makefile.setup`

All user-specific features and settings are enabled or disabled through `makefile.setup`. Here you can specify whether the application runs serially or in parallel via MPI, set the optimization level, enable X-Windows display support, and configure other build options.

The main parameters you need to select are:

```text
# Required parameter: identify the target host machine
SYS=[linux.intel | linux | gcc | aix | mac | bgp | bgq | mc-cc | wcr | cygwin | xt4]

MODE=SERIAL     # sets execution mode to serial
MODE=PARALLEL   # sets execution mode to parallel; requires MPI

XLIB_MODE=ON    # enables X-Window visualization; not recommended for production runs
XLIB_MODE=OFF   # disables X-Window visualization
```

All other parameters and settings are detailed in `makefile.setup`.

## Directory Structure

Brief description of the directories within this distribution:

* `./bin`      : executable applications, created during build
* `./src`      : C/C++ source files (`*.cc`)
* `./include`  : C/C++ include files (`*.h`)
* `./obj`      : object-file directories for parallel and serial builds
* `./docs`     : supporting documentation
* `./inputs`   : generic input files, including Rijm tables, FMM tables, gnuplot files, and X-Windows defaults
* `./tests`    : example tests (`*.ctrl`, `*.data`, `*.sh`)
* `./utils`    : support utilities
* `./tools`    : support tools

## Applications and Tools

The following applications and support utilities are created when building ParaDiS from the root directory:

* `./bin/paradis`        : main ParaDiS simulator
* `./bin/paradisconvert` : conversion utility for older, unsupported control and data files
* `./bin/calcdensity`    : dislocation density calculator
* `./bin/ctablegen`      : utility for creating PBC image-correction tables
* `./bin/ctablegenp`     : parallel/MPI utility for creating PBC image-correction tables
* `./bin/paradisgen`     : utility for creating initial random dislocation networks
* `./bin/paradisrepart`  : utility for repartitioning existing domain decompositions
* `./bin/stresstablegen` : utility for creating far-field stress tables

## Simulation Examples

There are numerous examples of simulation input files in the `tests/` directory.

ParaDiS is mainly controlled by two input files:

* a `.data` file, which specifies the simulation domain and initial dislocation network;
* a `.ctrl` file, which specifies material parameters, loading conditions, numerical parameters, mobility-law settings, and integration options.

Local simulations using `MobilityLaw_BccnlEshlby` or the drift-mode subcycling integrator may require control parameters that are not present in control files from the original LLNL public distribution.

Most of the control files in the `./tests` directory are configured for parallel execution using eight processors. This enables the tests to be run on a contemporary multicore workstation. The total number of MPI processes is controlled by the domain parameters:

```text
numXdoms = 2  # number of X-axis domains
numYdoms = 2  # number of Y-axis domains
numZdoms = 2  # number of Z-axis domains
```

The total domain count is:

```text
numXdoms * numYdoms * numZdoms
```

The domain count specified in the control file must be consistent with the number of MPI tasks launched for the simulation.

For the domain configuration above, an MPI-enabled ParaDiS run can be launched with:

```bash
mpirun -n 8 ./bin/paradis ./tests/mg-cAxis.ctrl
```

### Input-File Paths

Several test files use precomputed FMM tables located in the `./inputs` directory.

The paths in the ParaDiS control files generally assume that ParaDiS is launched from the main repository directory, for example:

```bash
./bin/paradis tests/example.ctrl
```

If ParaDiS is launched from another working directory, paths to the control file, data file, and any auxiliary files must be adjusted accordingly.

### FMM Correction Tables

Some examples use Fast Multipole Method (FMM) expansion tables. FMM image-correction tables can be generated using the `ctablegen` utility.

For an isotropic simulation:

```bash
./bin/ctablegen \
    -nu 3.327533e-01 \
    -mu 6.488424e+10 \
    -mporder 2 \
    -torder 5 \
    -outfile inputs/fmm-ctab.data
```

See the ParaDiS user's guide for additional information on constructing FMM correction tables.

### Initial Dislocation Networks

Initial dislocation networks can be generated using:

```bash
./bin/paradisgen
```

## Compatibility with the LLNL Base Version

This repository is derived from the LLNL ParaDiS Public Release Version 4.0.

The core ParaDiS functionality remains based on that release, while the local source tree contains additional developments, including:

```text
MobilityLaw_BccnlEshlby
Drift-mode subcycling integrator
```

As a result:

* simulations relying only on standard ParaDiS 4.0 functionality are intended to remain compatible with the LLNL base implementation;
* simulations using the local mobility law require the enhanced source tree;
* simulations using drift-mode subcycling require the enhanced integration implementation and its associated control parameters;
* input files containing locally introduced parameters may not run with an unmodified LLNL ParaDiS 4.0 executable.

When reproducing simulation results, users should therefore identify both the ParaDiS base version and whether these local enhancements were enabled.

## Citation

When publishing results obtained using this ParaDiS-derived code, please cite the original ParaDiS publication:

```bibtex
@article{arsenlis2007enabling,
  title={Enabling strain hardening simulations with dislocation dynamics},
  author={Arsenlis, Athanasios and Cai, Wei and Tang, Meijie and Rhee, Moono and Oppelstrup, Tomas and Hommes, Gregg and Pierce, Tom G and Bulatov, Vasily V},
  journal={Modelling and Simulation in Materials Science and Engineering},
  volume={15},
  number={6},
  pages={553},
  year={2007},
  publisher={IOP Publishing}
}
```

If results depend specifically on `MobilityLaw_BccnlEshlby` or the drift-mode subcycling integrator, the local implementation/version used for the simulations should also be identified in the associated publication or reproducibility documentation.

## License

ParaDiS is released under the BSD-3 license. See [LICENSE](LICENSE) for details.

The local enhancements in this repository should be used and distributed in accordance with the applicable repository license and any notices associated with the modified source files.

LLNL-CODE-853453

# Paradis_Local

Local ParaDiS development based on LLNL ParaDiS Public Release Version 4.0, including the `MobilityLaw_BccnlEshlby` mobility law and drift-mode subcycling integration enhancements.
