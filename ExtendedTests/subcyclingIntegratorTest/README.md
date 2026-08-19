# ParaDiS Integrator Comparison

This folder contains a comparison between the **subcycling integrator** implemented in ParaDiS-Extended and the standard **trapezoid integrator**.

Both simulations were run from the same initial configuration for approximately **18 minutes of wall-clock time**.

## Files

- `density_comparison.png` — dislocation density vs. simulation time
- `density_relative_error.png` — relative density difference between the two integrators
- `plastic_strain_comparison.png` — plastic strain vs. simulation time
- `plastic_strain_relative_error.png` — relative plastic-strain difference
- `density_comparison.csv` — interpolated density comparison data
- `plastic_strain_comparison.csv` — interpolated plastic-strain comparison data
- `compare_integrators.py` — Python script used to generate the comparison

## Performance

For approximately the same wall-clock time:

- **Subcycling** reached a simulation time of about `1.7e-5`
- **Trapezoid** reached a simulation time of about `9.5e-6`

This corresponds to an effective speedup of approximately:

\[
\frac{1.7\times10^{-5}}{9.5\times10^{-6}} \approx 1.8
\]

Thus, for this test case, the subcycling integrator advances the physical simulation approximately **1.8× faster** than the trapezoid integrator.

## Accuracy Comparison

### Plastic strain

The plastic-strain histories from the two integrators remain close over their common simulation-time interval.

- Maximum relative difference is approximately **6%**
- Typical difference is of the order of a few percent
- No continuously growing divergence is observed

### Dislocation density

The density agreement is very good during the early part of the simulation.

- Up to approximately `5.4e-6`, the density difference remains below about **0.5%**
- A transient divergence appears after approximately `5.5e-6`
- The maximum density difference is approximately **17%**
- The difference later decreases again to only a few percent

The temporary density divergence is likely associated with differences in the timing of discrete dislocation-network events such as junction formation, annihilation, collisions, or remeshing.

## Conclusion

For this test case, the subcycling integrator provides a significant computational speedup while retaining good agreement with the trapezoid solution in plastic strain and early-time dislocation-density evolution.

Further validation should compare the detailed network topology near the onset of the transient density divergence around `5.5e-6`.

## Generating the Plots

Run:

```bash
python compare_integrators.py /path/to/subcycling_results /path/to/trapezoid_results
```

Each results directory should contain:

```text
time_Plastic_strain
density
```

The comparison is performed at common physical simulation times using interpolation where necessary.
