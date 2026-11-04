# BRASIL matrix experiments

Four standalone R Markdown experiments use only installed `rSPDE` BRASIL
coefficient tables. Reusable coefficient, assembly, matrix-solver and error-diagnostic functions
are in `brasilfunc.Rmd`, including the sparse matrix time-stepping kernel.
Each notebook contains its own computation loops, regressions, plots and saves.

Open `NAP.Rproj`, then click Knit on any experiment. Each notebook follows
this order: purl and source the functions, chunk options, libraries, parameters,
heavy computation (`eval=FALSE`), load the errors, and plot. Knitting skips the
heavy computation and redraws the figures from saved results. No credentials
are needed. `docs/index.html` is the rendered entry point.

`purl()` regenerates `brasilfunc.R` from `brasilfunc.Rmd` before `source()` loads
it. Both files are kept in the project; edit the R Markdown source.

| Data | Notebook | Current settings |
| --- | --- | --- |
| Smooth, space | `cov_in_h_mat.Rmd` | m=8; h=2^(-7:-3); tau=2e-5; diffusion=0.001 |
| Smooth, time | `conv_in_tau_mat.Rmd` | m=8; h=0.001; N=200,63,20,6,2; diffusion=0.001 |
| Smooth, rational order | `conv_in_m_isolated.Rmd` | m=2:8; h=0.02; tau=0.005; diffusion=1 |
| Rough, space with sampled reference | `cov_in_h_rough_as_in_mat.Rmd` | m=14; J=5:9; reference cutoff=16; reference h=2^-18 |

All studies use kappa=4 and T=1. Each notebook defines its own alpha grid
and experiment settings in its parameter chunk.

## Rendering

Rendering needs `rmarkdown`, `knitr`, `ggplot2`, and Pandoc. The isolated-order
scalar plots also use `rSPDE` and `scales`. The combined rational/spatial figure
is created manually in the `eval=FALSE` chunk of `about.Rmd` from both studies'
saved results. Creation uses `gridExtra`, `tikzDevice`, a LaTeX installation
with TikZ and `mathpazo`, and `pdftoppm` (Poppler), and saves a 10 × 4 inch
PDF, PNG, and standalone TeX source. Knitting `about.Rmd` only displays the
saved PNG and does not require those creation tools. From the project root:

```r
rmarkdown::render("conv_in_m_isolated.Rmd", output_dir="docs")
# Or render the full site:
rmarkdown::render_site()
```

## Explicit recomputation

Edit the **Parameters** chunk, then run the notebook chunks from top to bottom.
Execute `recompute` manually when you want to run the experiment. It contains
the mesh assembly, refinement loops, error evaluation, and data saving.
The following chunks load the saved errors, calculate rates, draw plots, and
save the figures. Knitting runs those plotting chunks without rerunning the PDE.

Plots and tables use the variables in the parameter chunk directly. Keep those
parameters consistent with the saved errors when knitting. No settings are
extracted from saved metadata. No experiment parameter changes are needed in
`brasilfunc.Rmd`.

Recomputation needs Matrix, MetricGraph, rSPDE and ggplot2;
the smooth and rough matrix kernels additionally need Rcpp, RcppEigen and a C++
compiler.
The smooth notebooks compile their time-stepping/error kernel once in the
parent session before launching cases. Both smooth notebooks use nine workers,
one per alpha. In the spatial study each worker handles the current mesh; in
the temporal study each worker handles all five time grids. Each worker advances
every time step in C++ and uses the same centered reference-error formula. The R implementation remains
available for numerical regression checks. Changing `workers` affects the next
run; knitting only loads and plots saved results.

The installed rSPDE must provide the requested `m_brasil<m>t` tables.
Invalid signed coefficient rows are rejected, not altered. The shared graph
helper supports both current and older MetricGraph APIs.

The rational evolution always uses shifted sparse matrix solves and retains
the BRASIL constant term. The isolated rational-order study uses a complete numerical
FEM eigensystem only for its reference and independent diagnostic. The rough
spatial study uses invariant quarter-loop matrices with full-graph mass
normalization. The project contains only the four selected studies and their
saved data, figures, and rendered reports, together with shared functions.

Run the focused integration checks from the project root:

```sh
/usr/bin/Rscript --vanilla test_brasil.R
```

The sampled-reference rough spatial study uses a fixed finite-series reference
on a fine space-time mesh. Its cutoff, reference resolution, time steps, and
rational order were checked in DOCS/NAP; only the selected run is included here.
