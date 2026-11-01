# NAFDEMG matrix-study helpers

These Rmd files originate from the helpers in the NAFDEMG checkout used by
`cov_in_h_mat.Rmd` and `conv_in_tau_mat.Rmd`. They are kept in this
subdirectory so NAP's existing modal functions are unchanged. The fractional
evolution helper now reads rSPDE BRASIL coefficients through
`fun/brasil_coefficients.R`, retaining the shifted sparse Cholesky solves and
adding the constant remainder's mass-matrix solve. Original settings and
execution flags are preserved, including manually evaluated computation
chunks, reference grids, parallel workers, and output locations.

The original output root is `~/Desktop/DOCS/NAFDEMG/diffusion_experiments`.
`NAFDEMG_PRIVATE_DIR` can override it, as in the original implementation.
Generated R/hash files support NAP's existing `source_if_updated()` loader.
