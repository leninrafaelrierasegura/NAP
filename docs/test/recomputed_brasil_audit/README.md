# Recomputed BRASIL comparison

This experiment reads `conv_in_m_isolated.Rmd`'s settings chunk and holds its
FEM operator, initial data, forcing, backward Euler grid and error norm fixed.
It does not modify or execute the notebook, its setup, or any Slack calls.

The coefficient families are:

- `table`: the installed rSPDE tables, including the constant term.
- `same_interval`: freshly computed BRASIL approximations to `x^beta` on
  `[10^(-(5+m)/2), 1]`, degree `(m,m)`. This separates regeneration/precision
  from interval selection.
- `spectrum_interval`: freshly computed BRASIL approximations on the **fixed**
  interval `[min(0.001, kappa^2/lambda_max), 1]`. This includes the full discrete
  spectrum and the notebook's fixed scalar diagnostic. These coefficients are
  specific to the chosen mesh and scaling.
- `repair_example`: only the reported `(m,beta)=(15,0.95),(20,0.85)` cases,
  generated on the order-dependent interval. These do not enter the current
  `m<=14` convergence comparison.

Generation uses baryrat 2.1.2's BRASIL implementation with 70 decimal digits,
equioscillation tolerance `1e-8`, at most 2000 adjustment iterations, and 40
golden-section iterations for extrema. The conversion from barycentric form
to inverse-power poles/residues is performed at the same precision. The sparse
solver receives the resulting coefficients rounded to double precision.

The installed MetricGraph 1.1.2 predates the notebook's graph-construction API.
The audit uses its own `assemble_fem` routine on the prescribed tadpole
topology and uniform P1 subdivisions, with the notebook's analytic load
functions. It retains the full numerical eigensystem; no continuum-mode
assumption is introduced. This compatibility path is recorded in `context.rds`.

## Reproduce

From the NAP root, use the system R installation:

```bash
/usr/bin/Rscript --vanilla test/audit_recomputed_brasil.R --prepare
```

Use a Python environment with numpy, scipy and gmpy2. For this run,
`/home/rierasl/miniconda3/envs/fenicsenv/bin/python` already had those packages.
The two additional packages were installed in a temporary directory:

```bash
/home/rierasl/miniconda3/envs/fenicsenv/bin/python -m pip install \
  --target /tmp/nap-brasil-deps --no-deps baryrat==2.1.2 flamp==1.0.1
/home/rierasl/miniconda3/envs/fenicsenv/bin/python test/recompute_brasil.py \
  --jobs test/recomputed_brasil_audit/jobs.csv \
  --output test/recomputed_brasil_audit \
  --deps-dir /tmp/nap-brasil-deps --workers 2
/usr/bin/Rscript --vanilla test/audit_recomputed_brasil.R --compare
```

`fits.jsonl` retains the high-precision barycentric representation and
converted coefficients. `coefficients.csv` supplies double-precision values
to R. `fit_diagnostics.csv` records convergence, signs and conversion errors.
`comparison.csv`, `slopes.csv`, `matched_slopes.csv`, the PNG and `report.txt`
contain the matched scalar and space-time comparisons. Open plot markers are
below the notebook's empirical resolution indicator; slopes over those points
must not be interpreted as resolved convergence rates.

The fixed-spectrum fit changes the approximation domain and does not establish
a mesh-uniform improvement of the root-exponential asymptotic rate.

Primary method/API reference:
<https://baryrat.readthedocs.io/en/latest/index.html#baryrat.brasil>.
