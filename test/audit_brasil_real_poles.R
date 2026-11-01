# Run from NAP:
# /usr/bin/Rscript --vanilla test/audit_brasil_real_poles.R
# Reuses the coefficient audit and writes additional results under
# test/brasil_coefficient_audit/. Does not run the notebook or send messages.
#
# Mathematical guarantee, for ANY order m:
# p1 < ... < pm < 0, rj > 0, k >= 0, a > 0.
# D(y) = S(y)+a has D'(y) = -sum(rj/(y-pj)^2) < 0.
# On (-Inf,p1), D decreases from k+a>0 to -Inf.
# On (p[j-1],p[j]), it decreases from +Inf to -Inf.
# Hence exactly m simple real roots:
# q1 < p1 < q2 < p2 < ... < qm < pm < 0.
# These exhaust the degree-m denominator. Original poles pj are removable
# in H=S/(S+a), because H tends to 1 at each pj.
# Residue wj = S(qj)/S'(qj) = -a/S'(qj)
#            = a/sum(ri/(qj-pi)^2) > 0.
# Thus complex qj or wj are impossible under these assumptions in exact
# arithmetic. Floating-point root reconstruction can still fail at extremes.
# All physical parameters enter through a=tau*diffusion*kappa^(2*beta).

audit_brasil_real_poles <- function(project = ".") {
  base <- new.env(parent = globalenv())
  baseline_output <- capture.output(sys.source(
    file.path(project, "test/audit_brasil_coefficients.R"), base))
  baseline <- base$brasil_audit
  settings <- baseline$metadata$settings
  helper <- new.env(parent = globalenv())
  sys.source(file.path(project, "fun/brasil_coefficients.R"), helper)
  m <- settings$m_fixed
  tab <- get(paste0("m_brasil", m, "t"), asNamespace("rSPDE"))
  p <- as.matrix(tab[paste0("p", seq_len(m))])
  r <- as.matrix(tab[paste0("r", seq_len(m))])
  # Check a common column ordering across ALL rows. Linear interpolation then
  # preserves every strict pole gap, negative pole and positive residue.
  index <- order(p[1L, ])
  ordered <- p[, index, drop = FALSE]
  gaps <- ordered[, -1L, drop = FALSE] - ordered[, -m, drop = FALSE]
  stopifnot(all(is.finite(c(p, r, tab$k))), all(p < 0), all(r > 0),
            all(tab$k >= 0), all(gaps > 0), all(diff(tab$alpha) > 0))

  rows <- list()
  check_case <- function(beta, a, kind) {
    cf <- helper$brasil_coefficients(m, beta)
    result <- tryCatch(helper$brasil_partial_fractions(cf, time_step = a, scaling = 1),
                       error = function(e) e)
    if (inherits(result, "error")) {
      rows[[length(rows)+1L]] <<- data.frame(kind = kind, beta = beta, a = a,
        status = "numerical failure", complex_q = NA, complex_w = NA,
        max_imaginary_q = NA_real_, max_imaginary_w = NA_real_,
        max_reconstruction_error = NA_real_, detail = conditionMessage(result))
      return(invisible(NULL))
    }
    stopifnot(length(result$p) == m, length(result$r) == m,
      all(is.finite(c(result$p, result$r, result$k))),
      !is.complex(result$p), !is.complex(result$r),
      all(result$p < cf$poles),
      all(result$p[-1L] > head(cf$poles, -1L)), all(result$r > 0))
    rows[[length(rows)+1L]] <<- data.frame(kind = kind, beta = beta, a = a,
      status = "pass", complex_q = is.complex(result$p), complex_w = is.complex(result$r),
      max_imaginary_q = max(abs(Im(result$p))),
      max_imaginary_w = max(abs(Im(result$r))),
      max_reconstruction_error = result$reconstruction_error, detail = "")
  }
  # Notebook combinations: all alpha values and all actual time steps.
  for (beta in settings$alpha_vector/2) for (tau in baseline$metadata$tau_vector)
    check_case(beta, tau*settings$diffusion*settings$kappa^(2*beta), "notebook")
  # 999 stored beta rows, with a spanning ten orders of magnitude.
  for (beta in tab$alpha) for (a in c(1e-8, 1e-5, 1e-2, 1, 100))
    check_case(beta, a, "table row")
  # Every adjacent pair also gets a representative off-grid root calculation.
  midpoints <- head(tab$alpha, -1L) + diff(tab$alpha)/2
  for (beta in midpoints) check_case(beta, 1e-3, "interpolated midpoint")
  cases <- do.call(rbind, rows)
  report <- capture.output({
    cat("BRASIL real-pole/residue audit\n")
    cat("rSPDE", baseline$metadata$rSPDE_version, "; m =", m, "\n")
    cat(nrow(tab), "table rows have distinct ordered negative poles, positive residues, k>=0.\n")
    cat("Linear interpolation preserves those conditions throughout [",
        min(tab$alpha), ",", max(tab$alpha), "].\n")
    cat("Therefore ALL finite a>0 have real negative q_j and real positive w_j mathematically.\n")
    cat("Numerical case counts:\n")
    print(table(cases$kind, cases$status))
    cat("Complex q cases:", sum(cases$complex_q, na.rm = TRUE), "\n")
    cat("Complex w cases:", sum(cases$complex_w, na.rm = TRUE), "\n")
    cat("Maximum reconstruction discrepancy:", max(cases$max_reconstruction_error, na.rm = TRUE), "\n")
    failures <- cases[cases$status != "pass", , drop = FALSE]
    if (nrow(failures)) {
      cat("Numerical failures (distinct from complex roots):\n")
      print(failures, row.names = FALSE)
    }
    cat("No guarantee is made for invalid input coefficients, nonpositive a, or floating-point extremes.\n")
  })
  directory <- file.path(project, "test/brasil_coefficient_audit")
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(cases, file.path(directory, "real_poles_cases.csv"), row.names = FALSE)
  writeLines(report, file.path(directory, "real_poles_report.txt"))
  cat(paste(report, collapse = "\n"), "\n")
  invisible(list(cases = cases, ordered_pole_gaps_min = min(gaps),
                 table_beta_range = range(tab$alpha)))
}

brasil_real_poles_audit <- audit_brasil_real_poles()
