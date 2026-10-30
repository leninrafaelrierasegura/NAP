# Rough-forcing matrix experiment

On branch `rough/operator-based-coefficients`, the default coefficient
source is `computed`: directly generated Cheb–Padé root/factor rows from
`matlab/computed_fractional_roots.csv`, with **no beta interpolation**.
The default experiments use order eight. Order four and the preserved
BRASIL baselines are compared in `rough_operator_comparison.Rmd`.

`conv_in_h_rough.Rmd` displays the completed experiment. Its results are stored
in `data_files/rough_h_results.rds`, and PNG/PDF figures are in
`data_files/rough/`. Rendering the notebook does not rerun the experiment.

From the project root, reproduce the numerical results with:

```r
library(Matrix)
source("fun/rough/load.R")
result <- rough_run_experiment(workers = 3L)
rough_save_figures(result)
```

The runtime dependencies are Matrix, MetricGraph, sp (for the older graph
constructor), statmod, Rcpp,
RcppEigen, and ggplot2. Rendering additionally uses rmarkdown and knitr.
Use an R installation containing those packages; the verified environment
uses `/usr/bin/Rscript --vanilla` (R 4.5.2). A C++ compiler is needed when
recomputing, but not when rendering saved results.

## Files

- `problem.R`: uniform graph mesh, analytic sine–hat loads, mass projection,
  stable spatial projection residual, exact response integrals, and certified
  infinite-tail summation.
- `matrix_solver.R`: positive inverse approximation, bracketed partial
  fractions, shifted FEM factorizations, R matrix step, and compiled-solver
  interface. The constant partial-fraction term is retained.
- `matrix_evolution.cpp`: every backward Euler step of the full FEM vector,
  with reused sparse factorizations, exact time-moment error accumulation,
  and compensated summation. Exact signed nodal reflections reduce the
  shifted FEM matrices to their invariant quarter-loop subspace. No discrete
  eigendecomposition or modal evolution is used.
- `experiment.R`: final settings, complete cases, source-checked checkpoints,
  result provenance, and convergence-rate regressions.
- `plotting.R`: raw and logarithmically corrected error figures.
- `load.R`: loads definitions without starting computation.

The infinite forcing is fixed across all meshes and alpha values. Modes whose
loads vanish on a particular mesh still contribute to the exact error. A
fixed finite truncation of the continuum forcing would be a different
experiment.

Development scripts, independent checks, and their reports are kept separately
in `/home/rierasl/Desktop/DOCS/NAP/rough_experiment/all_alphas/`.
The final alpha sequence is `seq(0.5, 1.9, 0.2)`, with diffusion fixed at 1.
This branch's independent comparison work is saved in
`/home/rierasl/Desktop/DOCS/NAP/rough_operator_comparison/`.
The earlier directory preserves the BRASIL development record.

The selector `coefficient_source="computed"` requires a direct row at the
requested beta and rejects unavailable rows. `coefficient_source="brasil"`
is retained for baseline reproduction and requires rSPDE. The explicit
legacy selector `coefficient_source="operator"` reads the installed rSPDE
tables with the original interpolation; it is **not** used for this
branch's reported results. Orders above eight are unavailable in the direct
coefficient artifact.

The temporal experiments at small alpha exhibit rational error floors with
the direct coefficients, even at order eight. Passing matrix correctness
checks does not mean that those temporal slopes match the BRASIL baseline.
