# Rough-forcing temporal experiment

`conv_in_tau_rough.Rmd` displays completed matrix-solver results for
`alpha = seq(0.5, 1.9, 0.2)`. Results are in
`data_files/rough_tau_results.rds`; figures are in `data_files/rough_tau/`.

From the project root:

```r
library(Matrix)
source("fun/rough/tau/load.R")
result <- rough_tau_run_experiment(workers = 3L)
rough_tau_save_figures(result)
```

The same fixed infinite L2 forcing and exact solution used in the spatial
rough experiment are retained. Each curve fixes its mesh, diffusion and
rational order. Diffusion is normalized by
`nu_alpha = (kappa^2 + 4*pi^2)^(-alpha/2)` so the slowest forced mode has rate
one for every alpha. Actual steps are `T/N` for `N = 5,10,20,40,80,160`;
regressions use the three finest steps.

`experiment.R` assembles consistent P1 mass and stiffness matrices directly
on the invariant quarter-loop and evaluates every backward Euler step with
the existing shifted-matrix backend. The factor four in matrices, loads and
error integrals accounts for the full reflected loop. Matrix equivalence to
MetricGraph assembly and full-graph iterates is tested independently.

No modal recurrence evolves the numerical solution. Continuum modal formulas
are used for the exact solution and error integration, including a certified
infinite tail. The notebook evaluates the full continuum error, rather than
an error against a fixed finite-dimensional reference solution.

`plotting.R` gives one panel for every alpha so the nearly overlapping time
curves remain visible. Rendering the notebook only loads saved results.

The tested runtime uses `/usr/bin/Rscript --vanilla` with the same R packages
as the spatial experiment. A C++ compiler is required for recomputing, but
not for displaying saved results. The alpha=0.5 curve needs a very fine mesh;
recomputation uses about 2 GB of working memory for that curve.

Independent validation, refinement checks and reproduction scripts are kept
in `/home/rierasl/Desktop/DOCS/NAP/rough_tau_experiment/`.
