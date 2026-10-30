# Operator root/factor coefficient comparison

Completed 30 October 2026 on branch `rough/operator-based-coefficients`.

## Coefficient decision

The original request named the rSPDE operator tables. After checking their exponent grid, the user explicitly selected the local directly generated beta rows instead. The active source is `matlab/computed_fractional_roots.csv`, with its generation record in the adjacent JSON file. No beta interpolation or nearest-row substitution is used.

Orders four and eight were both run. The main spatial and temporal notebooks display order eight. These are operator-based Cheb–Padé roots/factor coefficients, not BRASIL residue tables and not the interpolated installed rSPDE table rows. All requested alpha values 0.5,0.7,...,1.9 are included.

## Preserved baseline and matched physics

Original BRASIL result objects, notebook sources, figures and helper sources are preserved in `data_files/rough_comparison/`. The baseline commit is 56b4ea5. Archived helper fingerprints match the provenance in both BRASIL result files.
The comparison code asserts equality of meshes, alphas, actual refinement targets/counts, diffusion, final time, shift, and fit ranges. Only coefficient family/order differs. The original infinite forcing, full continuum norm, tail bound and matrix evolution are retained.

## Results

At order eight, the finest spatial error ratios (operator / BRASIL) range from 0.996460 to 1.062502. The largest absolute change in corrected spatial slope is 0.035935.
Order four produces substantial rational contamination in the spatial curves. Order eight gives qualitatively similar spatial curves but its rational solution error is not negligible for every alpha.

| alpha| order| brasil_corrected| operator_corrected| finest_error_ratio|
|-----:|-----:|----------------:|------------------:|------------------:|
|   0.5|     4|        0.4763304|          0.4208213|          1.0804220|
|   0.7|     4|        0.6889651|          0.5391860|          1.2326316|
|   0.9|     4|        0.8963402|          0.5760259|          1.5638231|
|   1.1|     4|        1.1034983|          0.5483542|          2.1693530|
|   1.3|     4|        1.3145852|          0.5281124|          2.9944819|
|   1.5|     4|        1.5275283|          0.6168184|          3.5517512|
|   1.7|     4|        1.7037366|          0.9482364|          2.8486754|
|   1.9|     4|        1.7887415|          1.6167618|          1.2676571|
|   0.5|     8|        0.4763304|          0.4743278|          1.0041161|
|   0.7|     8|        0.6889651|          0.6814739|          1.0144378|
|   0.9|     8|        0.8963402|          0.8748292|          1.0378711|
|   1.1|     8|        1.1034983|          1.0686092|          1.0605857|
|   1.3|     8|        1.3145852|          1.2786502|          1.0625022|
|   1.5|     8|        1.5275283|          1.5102508|          1.0320207|
|   1.7|     8|        1.7037366|          1.7076706|          0.9973914|
|   1.9|     8|        1.7887415|          1.7915645|          0.9964603|

The temporal comparison is mixed: order eight is close to first order for alpha>=1.1, while alpha=0.5,0.7,0.9 show rational error floors on the original fine-step grid. Order four shows a larger floor. These findings are not labeled a matrix-solver failure.

| alpha| order| brasil_slope| operator_slope| finest_error_ratio|
|-----:|-----:|------------:|--------------:|------------------:|
|   0.5|     4|    0.9938899|      0.2687347|          3.6458095|
|   0.7|     4|    0.9942516|      0.6547694|          1.6933726|
|   0.9|     4|    0.9941478|      0.8814647|          1.1840974|
|   1.1|     4|    0.9940303|      0.9591847|          1.0531778|
|   1.3|     4|    0.9941227|      0.9841036|          1.0149640|
|   1.5|     4|    0.9940315|      0.9916843|          1.0034927|
|   1.7|     4|    0.9935611|      0.9932532|          1.0004696|
|   1.9|     4|    0.9939443|      0.9939217|          1.0000389|
|   0.5|     8|    0.9938899|      0.6588258|          1.6808887|
|   0.7|     8|    0.9942516|      0.9178137|          1.1208431|
|   0.9|     8|    0.9941478|      0.9750575|          1.0287326|
|   1.1|     8|    0.9940303|      0.9909429|          1.0045859|
|   1.3|     8|    0.9941227|      0.9937315|          1.0005791|
|   1.5|     8|    0.9940315|      0.9940174|          1.0000204|
|   1.7|     8|    0.9935611|      0.9935631|          0.9999969|
|   1.9|     8|    0.9939443|      0.9939447|          0.9999993|

## Numerical correctness checks

- Direct lookup, positivity and partial-fraction reconstruction were checked for all 64 order/alpha combinations (orders 1 through 8), with three representative step sizes each.
- A beta absent from the direct file is rejected; interpolation is not used as a fallback.
- Sixteen order/alpha fixtures compare the new positive partial fractions against the original polynomial matrix update, supplied with the same direct roots.
- Full-graph and lifted symmetry-reduced iterates agree on the same fixtures.
- All 176 convergence-case errors were checked against a 55-digit oracle evaluating the original root/factor PRODUCT directly, independently of the production partial-fraction conversion.
- Maximum matrix/oracle relative norm discrepancy: 3.2700085e-06.
- Independent spatial, temporal and rational solution-error components diagnose the floors; no false assumption that they must be negligible is imposed.
- Source helpers and notebook chunks parse, and the two results notebooks plus the comparison report render successfully.
- Baseline source/data and the pre-existing unrelated working-tree changes are preserved. No network publishing or messages were performed.

## Extra floor confirmation

At alpha=0.5 and order eight, additional actual matrix runs at N=320 and N=640 gave errors 0.0030923462 and 0.0029496738. The independently calculated continuous-time spatial/rational floor is 0.00290052384. Both extra runs agree with the independent oracle within 3e-6 relative error. They are excluded from the original regression range.

## Why direct rows do not guarantee sufficient accuracy

The coefficient artifact was generated for z^(beta-1) on [10^(-(5+m)/2),1]. The transformed inverse function is R(X)=X^(-1)*P(X^(-1))/Q(X^(-1)). Very fine small-alpha meshes have relevant inverse eigenvalues far below the generation interval, and finite-degree approximation errors also remain inside it. Direct evaluation at beta removes interpolation error but does not remove this rational approximation error.

## Artifacts and reproduction

- Report source: `rough_operator_comparison.Rmd`.
- Rendered report: `docs/rough_operator_comparison.html`.
- Main notebook outputs: `docs/conv_in_h_rough.html` and `docs/conv_in_tau_rough.html`.
- Raw candidates, baselines, tables, overlays, oracle diagnostics and extra floor probes: `data_files/rough_comparison/`.
- Automated test: `test/test_rough_operator_coefficients.R`.
- `run_h.R` and `run_tau.R` reproduce actual matrix cases, with source-checked checkpoints.
- `certify.py` and `oracle.py` perform independent high-precision checks.
- `floor_probe.R` performs the extra time refinements.
- The system R `/usr/bin/Rscript --vanilla` is used; the global R environment and project renv configuration are unchanged.

Conclusion: the available operator coefficients recover similar spatial curves at order eight, but do not reproduce all BRASIL temporal rates under the same original physical and refinement settings.
