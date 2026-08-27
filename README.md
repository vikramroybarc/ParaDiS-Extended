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

```text
MobilityLaw_BCC_nl_Eshelby_SegmentResist.cc
```

### 3. Eshelby-Compatible Subcycling

The subcycling implementation is designed to operate with simulations containing Eshelby inclusions, allowing inclusion forces and modified mobility behavior to be used together with force subcycling.

### 4. BCC Binary Junction Node Splitting

BCC three-arm binary junctions can now be split by both the serial and parallel multi-node topology paths, following the pinned [ExaDiS/OpenDiS topology implementation](https://github.com/llnl/exadis/tree/351607fa0f04aa34d29ee253438edde994473fd6/src/topology_types). The implementation identifies eligible non-planar BCC junctions, evaluates the glide-arm split energetics, marks the resulting degree-two corner node so remeshing preserves the physical junction geometry, and propagates that constraint across MPI domains.

The feature is controlled by `split3node` (enabled by default). Set `useParallelSplitMultiNode = 0` for the serial topology path or `1` for the parallel path in a `PARALLEL` build. Regression inputs are [`tests/bcc_binary_junction_node.ctrl`](tests/bcc_binary_junction_node.ctrl) and [`tests/bcc_binary_junction_node.data`](tests/bcc_binary_junction_node.data). See [`docs/BCC_Binary_Junction_Node_Splitting.tex`](docs/BCC_Binary_Junction_Node_Splitting.tex) for implementation and validation details.

## Building

The build procedure is the same as the base LLNL ParaDiS distribution. Configure `makefile.setup` and compile from the repository root:

```bash
make
```

## BCC Binary-Junction Regression Test

Run the one-step fixture from the `tests` directory:

```bash
cd tests
../bin/paradis -d bcc_binary_junction_node.data bcc_binary_junction_node.ctrl
```

The supplied control file uses the serial topology path. Change `useParallelSplitMultiNode` to `1` to exercise the parallel topology path in a `PARALLEL` build. To validate cross-domain constraint propagation, run the same fixture on two ranks:

```bash
mpirun -n 2 ../bin/paradis -doms 2 1 1 \
  -d bcc_binary_junction_node.data bcc_binary_junction_node.ctrl
```

A successful one-step run writes a restart containing five nodes: the three original pinned degree-one endpoints, one unconstrained degree-three node, and one degree-two node carrying constraint `16` (`CORNER_NODE`).

## Detailed implementation guides

- [Data structures and memory ownership (source)](docs/Data_Structures_and_Memory_Ownership.tex) documents the per-rank ownership tree, pooled nodes, arm arrays, segment views, cells, MPI ghosts and migration, allocation/free paths, and pointer invalidation rules.
- [Control-file to simulation code guide (source)](docs/Control_File_to_ParaDiS_Simulation_code.tex) documents parameter registration, parsing, precedence, validation, MPI broadcast, restart headers, and first-cycle consumers.
- [Topology and remeshing code guide (source)](docs/Topology_and_Remesh_code.tex) traces node splitting/merging, collision and MPI synchronization, mesh coarsening, and remesh rules 2 and 3.
- [BCC binary junction node splitting (overview)](docs/BCC_Binary_Junction_Node_Splitting.tex) describes the binary-junction test case and its topology behavior.
- [BCC binary junction node splitting (source walkthrough)](docs/BCC_Binary_Junction_Node_Splitting_code.tex) provides the detailed code-level explanation for that test.
- [Eshelby implementation notes (source)](docs/Eshelby_Implementation.tex) documents the Eshelby inclusion implementation.
- [Eshelby/FMM implementation notes (source)](docs/eshelby_fmm_paradis.tex) explains the ParaDiS Eshelby and FMM integration.
- [Nonlinear BCC mobility with Eshelby resistance](docs/MobilityLaw_BCC_nl_Eshelby_Resist.tex) documents the nonlinear BCC mobility law and resistance terms.
- [Subcycling integrator notes](docs/SubcyclingIntegrator.md) records the subcycling integrator design and usage notes.

### Documentation PDFs

- [ParaDiS 4.0 User Guide](docs/ParaDiS_4.0_UserGuide.pdf)
- [Data structures and memory ownership](docs/Data_Structures_and_Memory_Ownership.pdf)
- [BCC binary junction node splitting](docs/BCC_Binary_Junction_Node_Splitting.pdf)
- [BCC binary junction node splitting code guide](docs/BCC_Binary_Junction_Node_Splitting_code.pdf)
- [Eshelby implementation](docs/Eshelby_Implementation.pdf)
- [Eshelby/FMM notes](docs/Eshlby_FMM.pdf)
- [Nonlinear BCC mobility with Eshelby resistance](docs/MobilityLaw_BCC_nl_Eshelby_Resist.pdf)

### Documentation patch

- [BCC nonlinear Eshelby segment-weighting patch](docs/BCC_nl_Eshelby_segment_weighting.patch)

## Base Code

This repository is derived from **ParaDiS Public Release Version 4.0** developed at Lawrence Livermore National Laboratory.

For the original ParaDiS documentation, build instructions, examples, and citation information, refer to the upstream LLNL ParaDiS distribution.

## License

The original ParaDiS code is released under the BSD-3 license. See [LICENSE](LICENSE) for details.
