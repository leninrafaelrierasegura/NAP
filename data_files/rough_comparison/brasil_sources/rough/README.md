# Rough-forcing matrix experiment

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
constructor), rSPDE with the positive order-14 BRASIL tables, statmod, Rcpp,
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
Order 14 supplies positive BRASIL coefficients across that entire sequence;
the saved validation compares it with order 12 and the exact fractional rates.
