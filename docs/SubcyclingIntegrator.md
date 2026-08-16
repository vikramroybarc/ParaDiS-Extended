# Drift-mode subcycling integrator

This ParaDiS tree provides a force-subcycling timestep integrator based on
the drift algorithm in the pinned
[ExaDiS implementation](https://github.com/LLNL/exadis/blob/351607fa0f04aa34d29ee253438edde994473fd6/src/integrator_types/integrator_subcycling.h).
Drift mode is always enabled: only group 0 advances the physical network.
For every group-0 RKF stage, ParaDiS recomputes the live short-range,
self, applied, inclusion, and long-range contributions, adds the cached
segment-pair contributions from groups 1 through N-1, and then evaluates
mobility using that total nodal force.  Higher-group RKF steps only refresh
their force cache and adaptive timestep; their trial positions are restored.

This total-force mobility evaluation is required for nonlinear mobility
laws.  The implementation has been exercised with these selectors:

- `FCC_0`
- `BCC_0`
- `BCC_0b_Eshelby`
- `BCC_nl`
- `BCC_nl_Eshelby`
- `BCC_nl_Eshelby_Resist`
- `BCC_nl_Eshelby_SegmentResist` (alias for
  `BCC_nl_Eshelby_Resist`, whose source file carries the longer name)

The three Eshelby selectors require a build with `ESHELBY` enabled.

## Control-file settings

Select the integrator with:

```text
timestepIntegrator = subcycling
```

The optional controls and their defaults are:

```text
subcyclingNumGroups     = 5
subcyclingRadii         = [ 0.0 0.0 0.0 0.0 ]
subcyclingRtolThreshold = 1.0
subcyclingRtolRelative  = 0.1
subcyclingNextDT        = [ 0.0 0.0 0.0 0.0 0.0 ]
```

There may be two through five groups.  A run with all active radii set to
zero generates ExaDiS-style thresholds from `minSeg`, `rann`, and `maxSeg`.
Otherwise, provide `subcyclingNumGroups - 1` finite, nonnegative,
nondecreasing radii.  Shared-endpoint segment pairs and pairs separated by
less than one Burgers-vector unit are always placed in group 0.

`subcyclingNextDT` holds the adaptive RKF suggestion for each group.  Zero
entries start from ParaDiS `deltaTT`; accepted values are written to restart
control files so subgroup timestep history survives a restart.  The usual
`rTol`, `maxDT`, `dtIncrementFact`, `dtDecrementFact`, `dtExponent`, and
`dtVariableAdjustment` controls also apply.

## Force and state handling

Interaction groups are classified once at the beginning of each outer
cycle using stable endpoint tags and periodic segment distance.  The lists
and force caches are call-local, so remeshing and topology changes between
cycles cannot leave stale pointers.  Long-range and non-pair forces remain
in group 0.  ParaDiS performs a final full `NodeForce(FULL)` and mobility
evaluation after the accepted outer step, restoring coherent nodal forces,
per-arm forces, Eshelby intersection data, and velocities for collision and
topology code.  Outer-cycle old positions and velocities are preserved for
plastic-strain bookkeeping.

## Current restrictions

- `elasticinteraction` must be enabled.
- `includeInertia` is rejected because virtual subgroup trials require a
  separate inertial velocity history.
- `FULL_N2_FORCES` and `_ARLFEM` builds are rejected for this integrator.
- In parallel runs, subgroup mobility failures, error estimates, pair
  counts, and the resulting schedule are synchronized across MPI ranks.

