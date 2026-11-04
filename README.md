# BRASIL matrix experiments

Four standalone R Markdown experiments use only installed `rSPDE` BRASIL
coefficient tables. Reusable coefficient, assembly, matrix-solver and error-diagnostic functions
are in `brasilfunc.Rmd`, including the sparse matrix time-stepping kernel.
Each notebook contains its own computation loops, regressions, plots and saves.

Open `NAP.Rproj`, then click Knit on any experiment. Knitting reads the
committed results and images; it never starts an experiment, compiles C++,
changes saved results, or contacts a messaging service. No credentials are needed.
`docs/index.html` is the rendered entry point.

| Data | Notebook | Current settings |
| --- | --- | --- |
| Smooth, space | `cov_in_h_mat.Rmd` | m=8; h=2^(-7:-3); tau=2e-5; diffusion=0.001 |
| Smooth, time | `conv_in_tau_mat.Rmd` | m=8; h=0.001; N=200,80,32,13,6,2; diffusion=0.001 |
| Smooth, rational order | `conv_in_m_isolated.Rmd` | m=1:8; h=0.02; tau=0.005; diffusion=1 |
| Rough, space with sampled reference | `cov_in_h_rough_as_in_mat.Rmd` | m=14; J=5:9; reference cutoff=16; reference h=2^-18 |

All studies use alpha=0.5,0.7,...,1.9, kappa=4 and T=1. Each notebook lists
the settings used for its saved results.

## Rendering

R packages `rmarkdown` and `knitr`, plus Pandoc, are sufficient to display
saved results. From the project root:

```r
rmarkdown::render("conv_in_m_isolated.Rmd", output_dir="docs")
# Or render all six pages:
rmarkdown::render_site()
```

## Explicit recomputation

In each notebook, edit **Parameters for the next run**, then manually run the
`setup`, `parameters`, and `recompute` chunks in order. The `recompute` chunk contains the full experiment computation directly: mesh
assembly, refinement loops, error evaluation, regression, plotting, and saving.
It uses the notebook variables and overwrites its saved data and figures.
Knit afterward to refresh the report.

The parameter and recompute chunks have `eval=FALSE`: Knit skips both and
shows the existing saved results, even after you edit settings for a future
run. Displayed parameter tables are read from the saved result metadata.
No experiment parameter changes are needed in `brasilfunc.Rmd`.

Recomputation needs Matrix, MetricGraph, rSPDE, ggplot2 and sp;
the rough matrix kernel additionally needs Rcpp, RcppEigen and a C++ compiler.
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
