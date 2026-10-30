# Matched rough-forcing coefficient comparison

`brasil_h_results.rds`, `brasil_tau_results.rds`, their figures, and
`brasil_sources/rough/` preserve the inputs before this branch changed the
coefficient family. The baseline commit was `56b4ea5` (`tau rough added`).

`operator_m4_h.rds`, `operator_m8_h.rds`, `operator_m4_tau.rds`, and
`operator_m8_tau.rds` contain actual matrix-solver runs with the same
physical and refinement settings. The active source is
`matlab/computed_fractional_roots.csv`; rows are selected directly at beta,
without interpolation. This is the local operator-based Cheb–Padé artifact,
not the installed rSPDE table and not a BRASIL table.

`comparison.rds` combines the matched results and independent error
decompositions. CSV tables and PNG/PDF overlays are generated from that
object. `rough_operator_comparison.Rmd` explains the results.

The branch's main h/tau notebooks display the order-eight runs. The
preserved BRASIL results remain available independently of those outputs.
Development checks and execution records are stored in
`/home/rierasl/Desktop/DOCS/NAP/rough_operator_comparison/`.
