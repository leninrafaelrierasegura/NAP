# Audit of the supplied rSPDE/BRASIL generation scripts

The supplied originals (`ratapprox.py`, `brasil.m`, `rational_coeffs.m`) were
preserved. SHA-256 hashes are recorded in `generator_audit/original_sha256.txt`.

## Provenance

The Python script uses the same algorithm and order-dependent interval as our
successful table-reproduction pilot: degree `(m,m)` BRASIL applied to `x^beta`
on `[10^(-(5+m)/2),1]`. The MATLAB script builds rows with the same layout as
rSPDE: beta, residues, poles, constant.

Historical provenance cannot be established from these files alone. No source
commit or generation metadata is included. The filenames were absent from both
the public rSPDE tree at the installed commit
`088624e45edf5c1525f237a26c3ea07d30a3cb4e` and the inspected current `devel` tree
`6cd7cc549e6ea7b8b9e8d90c1ab771bb2648f8a3`.
They could be precursor scripts or one stage of a larger pipeline.
**As written, they cannot generate the correctly normalized installed tables
without an additional normalization step.**

## Confirmed normalization error

If the fitted rational function is

    B(x) = G * product(x-z_j) / product(x-t_j),

`ratapprox.py` exports only `t_j` and `z_j`, discarding `G`. MATLAB's `poly()`
reconstructs monic polynomials. Reversing their coefficient arrays implements
the reciprocal-variable change correctly, but does not recover `G`.
Consequently `residue(c,b)` represents `B(1/y)/G`, rather than `B(1/y)`.

Actual results from running the supplied Python script unchanged and the same
MATLAB `poly`/`residue` expressions:

| m | beta | Missing gain G | Supplied conversion at y=1 | Installed table at y=1 |
|---|---|---|---|---|
| 8 | 0.25 | 1.8026915 | 0.5547076424 | 0.9999667752 |
| 8 | 0.65 | 7.9534015 | 0.1257320897 | 0.9999977954 |
| 8 | 0.95 | 116.02348 | 0.0086189447 | 0.9999999539 |

With the gain restored, all six tested MATLAB conversions reconstruct their
Python barycentric functions within approximately `8e-13` on sampled grids.
The other cases are `(14,.95)`, `(15,.95)` and `(20,.85)`.
These executions use baryrat 2.1.2; the version used for the historical package
tables is unknown. The two negative table entries were **not reproduced** by
these executions. Their exact historical cause remains unverified.

## Improvements

1. Preserve normalization. A minimal change would export `B.gain()` and multiply
   the MATLAB numerator by it. The preferred conversion bypasses polynomial
   expansion entirely. If `B.polres()` returns x-variable poles `t` and residues
   `v`, the inverse-power representation is

       p = 1/t
       r = -v/t^2
       k = B(0)

   These operations must be done before rounding to double precision.
   [`../recompute_brasil.py`](../recompute_brasil.py) already implements this
   conversion at 70 decimal digits and retains the original barycentric data.
2. Request convergence information, check equioscillation, record versions,
   interval, order, exponent and precision, and validate signs and reconstruction.
   The supplied script uses default tolerances and does not inspect a returned
   convergence status. Monomial expansion is an avoidable conditioning risk;
   this audit does not establish it as the cause of the old bad table entries.
3. Parameterize the Python executable. The hard-coded Mac path in `brasil.m`
   cannot run on the current Linux host. Check the subprocess exit status before
   reading outputs and use distinct paths per beta/order/run to prevent stale
   files or concurrent-write collisions.
4. Fix or replace the final MATLAB plot. It hard-codes degree four despite
   `m=8`, and assigns the zero polynomial to the denominator and pole polynomial
   to the numerator. A valid diagnostic evaluates the actual fitted function
   and samples the actual fitting interval. The plot runs after the table save.
5. For the isolated fixed-mesh study, fitting the finite spectral interval is
   the principal source of the large observed improvement. Correcting these
   generator defects alone does not imply a large improvement over already
   correctly normalized table rows. See the previous matched comparison in
   [`../recomputed_brasil_audit/README.md`](../recomputed_brasil_audit/README.md).

## Reproduce the supplied-generator audit

From the NAP root, after installing the temporary dependencies documented in the
previous recomputation audit:

```bash
/usr/bin/Rscript --vanilla -e 'print(packageVersion("rSPDE"))'
/home/rierasl/miniconda3/envs/fenicsenv/bin/python \
  test/brasil/audit_original_generators.py --deps-dir /tmp/nap-brasil-deps
/home/rierasl/bin/matlab -batch \
  "addpath('test/brasil'); audit_matlab_conversion"
```

The exported installed rows are retained in
`generator_audit/installed_coefficients.csv`. The Python audit uses them when
present and can still run without them. Numeric results are in
`generator_audit/normalization.csv` and `generator_audit/matlab_conversion.csv`.
Case folders contain original stdout, root CSVs and raw MATLAB conversion data.

Primary API reference:
<https://baryrat.readthedocs.io/en/latest/index.html>.
