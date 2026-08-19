# ParaDiS-Extended

This repository is a modified version of **LLNL ParaDiS Public Release Version 4.0**.

The base ParaDiS source code is available from the LLNL ParaDiS repository. This README documents only the major changes introduced in **ParaDiS-Extended**.

## Changes from LLNL ParaDiS 4.0

### 1. Subcycling Integrator

A force-based **subcycling integration scheme** has been added, following the drift-mode subcycling approach used in ExaDiS.

The implementation includes:

- Separation of segment–segment interactions into force groups according to interaction distance
- Different update frequencies for different force groups
- RKF-based time integration within the subcycling framework
- Support for both linear and nonlinear mobility laws
- Automatic determination of force-group interaction radii
- Migration of interactions between groups when required by timestep/error constraints

The main implementation is contained in the subcycling-related source files, including `SubcyclingIntegrator.cc`.

### 2. Nonlinear BCC Mobility with Eshelby Resistance

A nonlinear BCC mobility law with resistance from **Eshelby inclusions** has been added.

The implementation:

- Extends the nonlinear BCC mobility formulation to dislocations interacting with Eshelby inclusions
- Determines the portions of dislocation segments intersecting inclusions
- Applies inclusion resistance using nodal shape-function weighting
- Handles edge, screw, and junction contributions within the nonlinear mobility calculation

The primary implementation is:

```

MobilityLaw_BCC_nl_Eshelby_SegmentResist.cc

````

### 3. Eshelby-Compatible Subcycling

The subcycling implementation is designed to operate with simulations containing Eshelby inclusions, allowing inclusion forces and modified mobility behavior to be used together with force subcycling.

## Building

The build procedure is the same as the base LLNL ParaDiS distribution. Configure `makefile.setup` and compile from the repository root:

```bash
make
````

## Base Code

This repository is derived from **ParaDiS Public Release Version 4.0** developed at Lawrence Livermore National Laboratory.

For the original ParaDiS documentation, build instructions, examples, and citation information, refer to the upstream LLNL ParaDiS distribution.

## License

The original ParaDiS code is released under the BSD-3 license. See [LICENSE](LICENSE) for details.

