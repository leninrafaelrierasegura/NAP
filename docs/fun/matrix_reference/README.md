# NAFDEMG matrix-study helpers

These Rmd files are verbatim copies of the helpers in the NAFDEMG checkout
used by `cov_in_h_mat.Rmd` and `conv_in_tau_mat.Rmd`. They are kept in this
subdirectory so NAP's existing modal functions are unchanged. The notebooks
only adapt their helper import paths; all original settings and execution
flags are preserved, including manually evaluated computation chunks,
package coefficients, reference grids, parallel workers, and private outputs.

The original output root is `~/Desktop/DOCS/NAFDEMG/diffusion_experiments`.
`NAFDEMG_PRIVATE_DIR` can override it, as in the original implementation.
Generated R/hash files support NAP's existing `source_if_updated()` loader.
