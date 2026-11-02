# BRASIL matrix experiments

Eight standalone R Markdown experiments use only installed `rSPDE` BRASIL
coefficient tables. Reusable coefficient, assembly, matrix-solver and error-diagnostic functions
are in `brasilfunc.Rmd`, including the sparse matrix time-stepping kernel.
Each notebook contains its own computation loops, regressions, plots and saves.

Open `NAP.Rproj`, then click Knit on any experiment. Knitting reads the
committed results and images; it never starts an experiment, compiles C++,
changes saved results, or contacts a messaging service. No credentials are needed.
`docs/index.html` is the rendered entry point.

| Data | Notebook | Current settings |
| --- | --- | --- |
| Smooth, space | `cov_in_h_mat.Rmd` | m=8; h=1/N; tau=1e-5; diffusion=0.001 |
| Smooth, time | `conv_in_tau_mat.Rmd` | m=8; h=0.02; tau=1/N; diffusion=0.001 |
| Smooth, rational order | `conv_in_m_isolated.Rmd` | m=1:8; h=0.02; tau=1e-5; diffusion=0.001 |
| Rough, space | `conv_in_h_rough.Rmd` | m=14; J=5:9; alpha-specific time steps; diffusion=1 |
| Rough, space with sampled reference | `cov_in_h_rough_as_in_mat.Rmd` | m=14; J=10:3; tau=1e-7; diffusion=0.25; reference cutoff=20 |
| Rough, time with sampled reference | `conv_in_tau_rough_as_in_mat.Rmd` | m=14; fixed h=2^-20; N=1024,512,256,128,64,32,16,8; diffusion=0.25 |
| Rough, time | `conv_in_tau_rough.Rmd` | m=20 for alpha=0.5, otherwise 14; N=5,10,20,40,80,160 |
| Rough, rational order | `conv_in_m_isolated_rough.Rmd` | m=1:9; h=2^-20; tau=1e-7; diffusion=0.25 |

The smooth h/tau studies use the same eight counts `N=200,104,54,28,14,7,4,2`, giving h/tau values from 0.005 to 0.5. The isolated m study uses their fixed h and tau.

The two sampled-reference rough studies and isolated rough m study use
alpha=0.7,0.9,...,1.7 and diffusion=0.25. Their eight shared h/tau values are
2^(-10:-3), spanning 128 as required by dyadic refinement. The isolated
study matches the fixed h=2^-20 and tau=1e-7 of those two studies. Both
sampled references use cutoff 20, reference h=2^-22 and reference tau=1e-5.
All sweep points and m=1:9 are retained; the common exclusion of alpha=0.5
and 1.9 and the full diagnostic screen are recorded in DOCS/NAP.

All studies use kappa=4 and T=1. The smooth studies and original rough
h/tau studies retain alpha=0.5,0.7,...,1.9. The original
`conv_in_tau_rough.Rmd` keeps its alpha-specific mesh and diffusion.

## Rendering

R packages `rmarkdown` and `knitr`, plus Pandoc, are sufficient to display
saved results. From the project root:

```r
rmarkdown::render("conv_in_m_isolated.Rmd", output_dir="docs")
# Or render all ten pages:
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

Recomputation needs Matrix, MetricGraph, rSPDE, ggplot2, statmod and sp;
the rough matrix kernel additionally needs Rcpp, RcppEigen and a C++ compiler.
The installed rSPDE must provide the requested `m_brasil<m>t` tables.
Invalid signed coefficient rows are rejected, not altered. The shared graph
helper supports both current and older MetricGraph APIs.

The rational evolution always uses shifted sparse matrix solves and retains
the BRASIL constant term. The smooth isolated study uses a complete numerical
FEM eigensystem for its reference and independent diagnostic. The three
consistent rough studies use an exact restriction to every forced P1 sine
mode, with full-graph mass normalization. Unforced FEM coordinates remain
exactly zero. The rough isolated reference is streamed; its rational solver
still uses shifted sparse matrix solves, checked by an independent product
recurrence. Physical FEM and evolved dimensions are recorded separately. Its parameter diagnostics remain in DOCS/NAP; the project
contains only the selected run and installed BRASIL coefficient rows.

Run the focused integration checks from the project root:

```sh
/usr/bin/Rscript --vanilla test_brasil.R
```

The sampled-reference rough spatial study uses a fixed finite-series reference
on a fine space-time mesh. Its cutoff, reference resolution, time steps, and
rational order were checked in DOCS/NAP; only the selected run is included here.

The consistent sampled-reference rough studies retain the original infinite
forcing and sampled P1 reference-mass/right-endpoint-time norm. The exact
invariant reduction was checked against the original matrices, loads,
full-eigensystem reference, and sampled norm on small meshes. It avoids large
mesh arrays and stored reference histories at the shared fine settings.
Reference sensitivity, all candidate parameters, and reproducible diagnostics
are retained in
`/home/rierasl/Desktop/DOCS/NAP/consistent_rough_parameters_2026-11-02`.
