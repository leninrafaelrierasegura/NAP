# Fixed-interval BRASIL study

Open `conv_in_m_isolated_fixed_int.Rmd`, edit its settings, then run all chunks
or Knit. Use the project R environment with Matrix, MetricGraph, here, ggplot2,
jsonlite, knitr and rmarkdown. On this computer the working system installation
is `/usr/bin/Rscript --vanilla`; the old R-4.2 renv library is empty.

The default interval is `[min(kappa^2/lambda_max, 1/diagnostic_y_max), 1]`.
It is selected once after assembling the complete FEM spectrum and stays fixed
for every order and exponent. An explicit `fixed_x_min` can widen the interval;
settings that miss the full spectrum or scalar diagnostic are rejected.
Coefficients are computed directly at every requested beta; no table lookup,
beta interpolation, audit artifacts, or saved result files are used.

## Python setup

The generator requires numpy, scipy, gmpy2, baryrat and flamp. Create a local
environment (ignored by Git), then install its dependencies:

```bash
python3 -m venv --system-site-packages .venv-fixed-int
.venv-fixed-int/bin/python -m pip install -r fun/fixed_int/requirements.txt
.venv-fixed-int/bin/python fun/fixed_int/generate_coefficients.py --check
```

The R helper automatically selects `.venv-fixed-int/bin/python` (or the Windows
equivalent). To use a different Python environment, set `python` in the notebook
or set `NAP_BRASIL_PYTHON`. It must have all of the above dependencies.
There is no automatic package installation during execution.

The default fit uses 70 decimal digits and relative equioscillation tolerance
`1e-8`. It checks convergence and signed alternation of the extrema, converts
the barycentric representation at the same precision, and then rounds only
the final source poles, residues and constant to double precision. The bank
retains high-precision coefficients and the barycentric representation in
memory for inspection. Fitting jobs and their JSON exchange files are temporary
and removed after loading. Every run fits the current settings afresh.

## Helpers and mathematics

- `coefficients.R`: interval coverage, Python bridge, coefficient validation and
  direct lookup of freshly generated coefficients.
- `generate_coefficients.py`: parallel high-precision BRASIL generation and
  conversion/reconstruction checks. It can also run without R via its JSON CLI.
- `experiment.R`: the same tadpole data and exact discrete reference as the
  original isolated study, plus an explicit coefficient-based sparse solver
  and independent full-eigenbasis diagnostic.

Shared mesh/load/reference definitions are read from
`fun/isolated_rational_convergence.R` into a private R environment. They and the
original study are not modified. With MetricGraph 1.1.2, which predates the
notebook's `edges` constructor, a compatibility path calls MetricGraph's own
P1 FEM assembler on the prescribed two-edge tadpole. It uses the same physical
lengths, uniform mesh subdivisions, consistent mass matrix and analytic loads.

For `S(y)=k+sum(r_j/(y-p_j))`, the step is
`H(y)=S(y)/(S(y)+tau*diffusion*kappa^(2*beta))`. Real scalar bracketing gives
the time-step poles and residues. The update includes the mass solve from
`k/(k+tau*diffusion*kappa^(2*beta))`; it does not expand polynomial products.
All shifted Cholesky factors are reused throughout the time history.

The reference and rational solution share all spatial and temporal data.
Only the coefficients vary. The optional installed-table comparison uses that
same solver and skips unsupported/invalid rows. It never limits the orders
available to fresh coefficients. Open points are below the empirical
resolution indicator; only resolved points enter the fitted slopes.

Primary BRASIL implementation/API:
<https://baryrat.readthedocs.io/en/latest/index.html#baryrat.brasil>.

Render from NAP without running the site build:

```bash
/usr/bin/Rscript --vanilla -e \
  'rmarkdown::render("conv_in_m_isolated_fixed_int.Rmd", output_dir="docs", envir=new.env())'
```

Run the integration checks with:

```bash
/usr/bin/Rscript --vanilla test/test_fixed_int.R
```

They compare the shifted sparse solves with the direct dense inverse-power
scheme for noncommuting mass/stiffness matrices, verify complete time
histories, and exercise non-tabulated exponents, stale-interval rejection,
the integer-power bypass and an order above the package-table limit.
