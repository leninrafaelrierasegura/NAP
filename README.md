# BRASIL matrix experiments

Five standalone R Markdown experiments use only installed `rSPDE` BRASIL
coefficient tables. All theory and implementation functions are in
`brasilfunc.Rmd`, including the sparse matrix time-stepping kernel.

Open `NAP.Rproj`, then click Knit on any experiment. Knitting reads the
committed results and images; it never starts an experiment, compiles C++,
changes saved results, or contacts a messaging service. No credentials are needed.
`docs/index.html` is the rendered entry point.

| Data | Notebook | Current settings |
| --- | --- | --- |
| Smooth, space | `cov_in_h_mat.Rmd` | m=8; h=2^(-7:-3); tau=2e-5; diffusion=0.001 |
| Smooth, time | `conv_in_tau_mat.Rmd` | m=8; h=0.001; N=200,80,32,13,6,2; diffusion=0.001 |
| Smooth, rational order | `conv_in_m_isolated.Rmd` | m=1:14; h=0.02; tau=0.005; diffusion=1 |
| Rough, space | `conv_in_h_rough.Rmd` | m=14; J=5:9; alpha-specific time steps; diffusion=1 |
| Rough, time | `conv_in_tau_rough.Rmd` | m=20 for alpha=0.5, otherwise 14; N=5,10,20,40,80,160 |

All studies use alpha=0.5,0.7,...,1.9, kappa=4 and T=1. The rough temporal
study keeps its alpha-specific mesh and diffusion. Each notebook lists them.

## Rendering

R packages `rmarkdown` and `knitr`, plus Pandoc, are sufficient to display
saved results. From the project root:

```r
rmarkdown::render("conv_in_m_isolated.Rmd", output_dir="docs")
# Or render all seven pages:
rmarkdown::render_site()
```

## Explicit recomputation

The final `eval=FALSE` chunk in each experiment calls its runner from
`brasilfunc.Rmd`. Run that chunk manually to replace the saved results.
Recomputation needs Matrix, MetricGraph, rSPDE, ggplot2, statmod and sp;
the rough matrix kernel additionally needs Rcpp, RcppEigen and a C++ compiler.
The installed rSPDE must provide the requested `m_brasil<m>t` tables.
Invalid signed coefficient rows are rejected, not altered. The shared graph
helper supports both current and older MetricGraph APIs.

The rational evolution always uses shifted sparse matrix solves and retains
the BRASIL constant term. The isolated study uses a complete numerical
FEM eigensystem only for its reference and independent diagnostic.

Run the focused integration checks from the project root:

```sh
/usr/bin/Rscript --vanilla test_brasil.R
```
