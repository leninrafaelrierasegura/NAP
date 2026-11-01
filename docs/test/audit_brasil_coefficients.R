# Audit the coefficients used by conv_in_tau_mat.Rmd.
# Run from NAP (use system R: the project renv library may be empty):
#   /usr/bin/Rscript --vanilla test/audit_brasil_coefficients.R
#   /usr/bin/Rscript --vanilla test/audit_brasil_coefficients.R --write-results
# Optional, with a newer rSPDE exposing WL2:
#   /usr/bin/Rscript --vanilla test/audit_brasil_coefficients.R --wl2 --write-results
# Or source this file from NAP; results remain in `brasil_audit`.
# No notebook execution, package installation, saved-study updates, or messages.
#
# BRASIL tables are obtained with rSPDE:::get_rational_coefficients() or directly
# as rSPDE:::m_brasil8t; interp_rational_coefficients() evaluates their columns.
# rational.coefficients.wl2() FITS a different, weighted-L2 approximation. It
# does not retrieve BRASIL tables. k_term=TRUE adds a fitted constant; it does
# not change the fitting objective into BRASIL's minimax objective.
# Primary reference (also documents s, k_term, and the exponent conventions):
# https://davidbolin.github.io/rSPDE/reference/rational.coefficients.wl2.html
#
# With beta = notebook alpha/2, the table represents
#   S(y) = k + sum_j r_j/(y-p_j) ~ y^(-beta), y >= 1,
#   B(x) = S(1/x) = k + sum_j r_j*x/(1-p_j*x) ~ x^beta.
# B(0)=k is essential: dropping it changes the approximation and time step.
# For normalized eigenvalues y=lambda/kappa^2 the backward Euler multiplier is
#   H(y)=S(y)/(S(y)+tau*diffusion*kappa^(2*beta)).
# All approximation-error maxima below are sampled, not certified suprema.

audit_brasil_coefficients <- function(project = ".", run_wl2 = FALSE,
                                      write_results = FALSE) {
  project <- normalizePath(project, mustWork = TRUE)
  notebook <- file.path(project, "conv_in_tau_mat.Rmd")
  if (!file.exists(notebook)) stop("Run from the NAP root or set project explicitly.")
  if (!requireNamespace("rSPDE", quietly = TRUE))
    stop("rSPDE is unavailable in this R library. Try /usr/bin/Rscript --vanilla.")
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Matrix is required.")
  ns <- asNamespace("rSPDE")
  helper <- new.env(parent = globalenv())
  sys.source(file.path(project, "fun/brasil_coefficients.R"), helper)
  # The production module attaches Matrix. Bind its S4 solve generic locally
  # so this standalone audit gets the same dispatch without attaching packages.
  helper$solve <- Matrix::solve

  # Read only the relevant numeric assignments from R chunks. Never run the
  # notebook setup: it contains Slack calls and expensive experiment code.
  lines <- readLines(notebook, warn = FALSE)
  code <- character()
  inside <- FALSE
  for (line in lines) {
    if (grepl("^```\\{r", line)) { inside <- TRUE; next }
    if (grepl("^```\\s*$", line)) { inside <- FALSE; next }
    if (inside) code <- c(code, line)
  }
  settings <- new.env(parent = baseenv())
  wanted <- c("T_final", "diffusion", "kappa", "tau_targets", "h_fixed",
              "m_fixed", "alpha_vector")
  for (expr in parse(text = code)) {
    if (is.call(expr) && as.character(expr[[1]])[1] %in% c("<-", "=") &&
        is.symbol(expr[[2]]) && as.character(expr[[2]]) %in% wanted)
      eval(expr, settings)
  }
  stopifnot(all(wanted %in% ls(settings)))
  settings <- as.list(settings)
  stopifnot(settings$m_fixed %in% helper$brasil_supported_orders(),
            settings$kappa > 0, settings$diffusion > 0,
            all(settings$alpha_vector > 0 & settings$alpha_vector < 2))
  betas <- settings$alpha_vector/2
  m <- settings$m_fixed
  # Match the notebook's T/N steps; these current targets are integer ratios.
  ratios <- settings$T_final/settings$tau_targets
  stopifnot(max(abs(ratios-round(ratios))) < 1e-10)
  taus <- settings$T_final/round(ratios)

  # Load just the two production solver definitions, without their setup code.
  solver_file <- file.path(project, "fun/matrix_reference/fractional_evolution.R")
  solver_names <- c("my.fractional.operators.frac", "my.solver.frac")
  for (expr in parse(solver_file)) {
    if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
        is.symbol(expr[[2]]) && as.character(expr[[2]]) %in% solver_names &&
        is.call(expr[[3]]) && identical(expr[[3]][[1]], as.name("function")))
      eval(expr, helper)
  }
  stopifnot(all(solver_names %in% ls(helper)))
  stopifnot(grepl("brasil_coefficients", paste(deparse(
    body(helper$my.fractional.operators.frac)), collapse = " "), fixed = TRUE))

  evaluate <- function(y, cf) {
    cf$k + as.numeric((1/outer(y, cf$p, "-")) %*% cf$r)
  }
  table_at <- function(order, beta) {
    tab <- get(paste0("m_brasil", order, "t"), ns)
    index <- which(abs(tab$alpha-beta) <= 32*.Machine$double.eps)
    at <- function(column) if (length(index) == 1L) tab[[column]][index] else
      stats::approx(tab$alpha, tab[[column]], xout = beta)$y
    list(k = unname(at("k")),
         p = unname(vapply(paste0("p", seq_len(order)), at, numeric(1))),
         r = unname(vapply(paste0("r", seq_len(order)), at, numeric(1))),
         exact_row = length(index) == 1L)
  }
  sorted_values <- function(cf) c(cf$k, cf$p[order(cf$p)], cf$r[order(cf$p)])

  # Inventory every installed order, including signed/duplicate/nonfinite rows.
  inventory <- do.call(rbind, lapply(helper$brasil_supported_orders(), function(order) {
    tab <- get(paste0("m_brasil", order, "t"), ns)
    getter <- get("get_rational_coefficients", ns)(order, "brasil")
    stopifnot(identical(tab, getter))
    p <- as.matrix(tab[paste0("p", seq_len(order))])
    r <- as.matrix(tab[paste0("r", seq_len(order))])
    valid <- rowSums(!is.finite(p) | p >= 0 | !is.finite(r) | r <= 0) == 0 &
      is.finite(tab$k) & tab$k >= 0 & apply(p, 1, anyDuplicated) == 0L
    data.frame(m = order, rows = nrow(tab), beta_min = min(tab$alpha),
      beta_max = max(tab$alpha), invalid_rows = sum(!valid))
  }))

  # Small linear-FEM interval: independent generalized eigensystem verifies the
  # actual sparse solver, all alpha/tau settings, normalization and mass term.
  # This is an algebraic solver check, not a rerun of the tadpole convergence study.
  n <- 9L
  h <- 1/(n-1L)
  C <- G <- matrix(0, n, n)
  for (j in seq_len(n-1L)) {
    ij <- j:(j+1L)
    C[ij, ij] <- C[ij, ij] + h/6*matrix(c(2, 1, 1, 2), 2)
    G[ij, ij] <- G[ij, ij] + 1/h*matrix(c(1, -1, -1, 1), 2)
  }
  L <- settings$kappa^2*C + G
  R <- chol(C)
  Rinv <- solve(R)
  eig <- eigen(t(Rinv) %*% L %*% Rinv, symmetric = TRUE)
  V <- Rinv %*% eig$vectors
  rhs <- sin(seq_len(n)) + seq_len(n)/n
  sparse_C <- Matrix::Matrix(C, sparse = TRUE)
  sparse_L <- Matrix::Matrix(L, sparse = TRUE)
  # Includes the range [1, 1+12/(kappa*h_fixed)^2] as well as a longer tail.
  # The FE bound is a uniform-interval benchmark, not an assembled tadpole spectrum.
  y_bound <- 1 + 12/(settings$kappa*settings$h_fixed)^2
  y_mesh <- exp(seq(0, log(y_bound), length.out = 4001L))
  y <- sort(unique(c(y_mesh,
                     10^seq(0, 12, length.out = 4001L))))
  coeff_rows <- case_rows <- transfer_rows <- list()
  for (beta in betas) {
    raw <- table_at(m, beta)
    cf <- helper$brasil_coefficients(m, beta)
    adapted <- list(k = cf$constant, p = cf$poles, r = cf$residues)
    # Sorting must preserve residue/pole pairing and every stored double.
    stopifnot(identical(sorted_values(raw), sorted_values(adapted)))
    package_cf <- get("interp_rational_coefficients", ns)(m, "brasil", "linear", beta)
    package_gap <- max(abs(sorted_values(raw)-sorted_values(package_cf))/
                         pmax(1, abs(sorted_values(raw))))
    stopifnot(package_gap < 1e-11)
    S <- evaluate(y, raw)
    S_mesh <- evaluate(y_mesh, raw)
    stopifnot(max(abs(S-helper$brasil_inverse(y, cf))) < 1e-12,
              max(abs(S-helper$brasil_power(1/y, cf))) < 1e-12,
              identical(helper$brasil_power(0, cf), cf$constant))
    coeff_rows[[length(coeff_rows)+1L]] <- data.frame(
      alpha = 2*beta, beta = beta, m = m, j = seq_len(m),
      pole = cf$poles, residue = cf$residues, constant = cf$constant,
      exact_table_row = raw$exact_row)
    case_rows[[length(case_rows)+1L]] <- data.frame(
      alpha = 2*beta, beta = beta, m = m, exact_table_row = raw$exact_row,
      package_linear_scaled_gap = package_gap, constant = cf$constant,
      sampled_max_abs_error_mesh_benchmark = max(abs(S_mesh-y_mesh^(-beta))),
      sampled_max_relative_error_mesh_benchmark = max(abs(S_mesh/y_mesh^(-beta)-1)),
      sampled_max_abs_error_y_1_to_1e12 = max(abs(S-y^(-beta))),
      sampled_max_relative_error_y_1_to_1e12 = max(abs(S/y^(-beta)-1)),
      error_at_y_1 = abs(evaluate(1, raw)-1),
      error_at_x_0 = cf$constant)
    for (tau in taus) {
      strength <- tau*settings$diffusion*settings$kappa^(2*beta)
      transfer <- helper$brasil_partial_fractions(cf, tau,
        settings$diffusion*settings$kappa^(2*beta))
      transfer_gap <- max(abs(evaluate(y, transfer)-S/(S+strength)))
      H_mesh <- S_mesh/(S_mesh+strength)
      H_exact <- 1/(1+strength*y_mesh^beta)
      H_without_k <- (S_mesh-raw$k)/(S_mesh-raw$k+strength)
      op <- helper$my.fractional.operators.frac(sparse_L, beta, sparse_C,
        scale.factor = settings$kappa^2, m = m, time_step = tau,
        diffusion = settings$diffusion)
      Seig <- evaluate(eig$values/settings$kappa^2, raw)
      expected <- as.numeric(V %*% (Seig/(Seig+strength)*as.numeric(crossprod(V, rhs))))
      actual <- as.numeric(helper$my.solver.frac(op, rhs))
      solver_gap <- max(abs(actual-expected))/max(1, max(abs(expected)))
      stopifnot(transfer_gap < 2e-12, solver_gap < 2e-10,
                length(op$residues) == m, length(op$factorizations) == m,
                op$remainder > 0, !is.null(op$mass_factor))
      transfer_rows[[length(transfer_rows)+1L]] <- data.frame(
        alpha = 2*beta, beta = beta, m = m, tau = tau,
        strength = strength, transfer_constant = transfer$k,
        reconstruction_max_abs_gap = transfer_gap, solver_scaled_gap = solver_gap,
        sampled_BE_multiplier_abs_error_mesh_benchmark = max(abs(H_mesh-H_exact)),
        sampled_BE_change_if_k_dropped_mesh_benchmark = max(abs(H_mesh-H_without_k)))
    }
  }

  # Off-grid beta distinguishes linear interpolation from spline defaults.
  off_beta <- 1/3
  off <- helper$brasil_coefficients(m, off_beta)
  raw <- table_at(m, off_beta)
  linear <- get("interp_rational_coefficients", ns)(m, "brasil", "linear", off_beta)
  spline <- get("interp_rational_coefficients", ns)(m, "brasil", "spline", off_beta)
  off_gap <- max(abs(sorted_values(raw)-sorted_values(linear))/
                  pmax(1, abs(sorted_values(raw))))
  stopifnot(!raw$exact_row, off_gap < 1e-11,
    identical(sorted_values(raw), sorted_values(list(k = off$constant,
              p = off$poles, r = off$residues))))
  interpolation <- data.frame(beta = off_beta, m = m,
    package_linear_scaled_gap = off_gap,
    linear_vs_spline_sampled_function_gap = max(abs(evaluate(y, linear)-evaluate(y, spline))))

  wl2_available <- "rational.coefficients.wl2" %in% getNamespaceExports("rSPDE")
  wl2 <- list(status = if (!wl2_available) "unavailable in installed rSPDE" else
                if (!run_wl2) "available; use --wl2 to fit a comparison" else "requested")
  if (run_wl2 && wl2_available) {
    # A finite interval permits a nonzero constant. d=1 for a metric graph.
    # Covariance alpha=beta (q=0) targets x^beta. Operator alpha=notebook alpha
    # also targets x^(alpha/2), but has m+1 poles, so is not order matched.
    # Keep WL2 to one representative case because it performs an optimization.
    beta <- betas[which.min(abs(betas-.65))]
    fun <- getExportedValue("rSPDE", "rational.coefficients.wl2")
    args <- list(alpha = beta, d = 1, m = m, x_min = 1/y_bound,
                 type = "covariance", seed = 1L)
    has_constant_option <- "k_term" %in% names(formals(fun))
    if (has_constant_option) args$k_term <- TRUE
    wl2 <- tryCatch({
      fit <- do.call(fun, args)
      if (is.null(fit$k)) fit$k <- 0
      stopifnot(length(fit$p) == m, length(fit$r) == m,
                all(is.finite(c(fit$p, fit$r, fit$k))), all(fit$p < 1),
                all(fit$r >= 0), fit$k >= 0, is.null(fit$q) || fit$q == 0)
      yy <- exp(seq(0, log(y_bound), length.out = 4001L))
      tabulated <- table_at(m, beta)
      list(status = "fitted separately; these are not tabulated BRASIL coefficients",
        arguments = args, k_term_supported = has_constant_option, fit = fit,
        comparison = data.frame(beta = beta, m = m,
          wl2_reported_relative_weighted_L2_error = fit$rel_err,
          wl2_sampled_max_abs_error = max(abs(evaluate(yy, fit)-yy^(-beta))),
          brasil_sampled_max_abs_error = max(abs(evaluate(yy, tabulated)-yy^(-beta))),
          sampled_function_gap = max(abs(evaluate(yy, fit)-evaluate(yy, tabulated)))))
    }, error = function(e) list(status = "WL2 comparison failed", error = conditionMessage(e)))
  }
  result <- list(
    metadata = list(target = notebook, settings = settings, tau_vector = taus,
      R = R.version.string, rSPDE_version = as.character(utils::packageVersion("rSPDE")),
      rSPDE_path = find.package("rSPDE"),
      rSPDE_description = utils::packageDescription("rSPDE"),
      source_md5 = tools::md5sum(c(notebook, solver_file,
        file.path(project, "fun/brasil_coefficients.R"))),
      mesh_benchmark_y_upper = y_bound,
      session = utils::sessionInfo()),
    inventory = inventory, coefficients = do.call(rbind, coeff_rows),
    approximation = do.call(rbind, case_rows), transfer = do.call(rbind, transfer_rows),
    interpolation = interpolation, wl2 = wl2)
  report <- capture.output({
    cat("BRASIL coefficient audit: conv_in_tau_mat.Rmd\n")
    cat("rSPDE", result$metadata$rSPDE_version, "at", result$metadata$rSPDE_path, "\n")
    cat("m =", m, "; beta = notebook alpha/2; diffusion =", settings$diffusion, "\n")
    cat("Direct extraction: get('m_brasil8t', asNamespace('rSPDE')) [for m=8].\n")
    cat("Package getter: rSPDE:::get_rational_coefficients(8, 'brasil').\n")
    cat("Interpolation: rSPDE:::interp_rational_coefficients(8, 'brasil', 'linear', beta).\n")
    cat("Keep pole/residue pairs together and retain k.\n\n")
    cat("PASS: helper coefficients exactly equal raw table rows/interpolants.\n")
    cat("PASS: reciprocal evaluation and B(0)=k.\n")
    cat("PASS:", nrow(result$transfer), "production sparse-solver/time-step cases.\n")
    cat("Worst transfer reconstruction gap:", max(result$transfer$reconstruction_max_abs_gap), "\n")
    cat("Worst sparse-solver scaled gap:", max(result$transfer$solver_scaled_gap), "\n\n")
    cat("Mesh benchmark: y in [1,", format(y_bound), "] from 1+12/(kappa*h_fixed)^2.\n")
    cat("This is a uniform-interval FE bound, not the assembled tadpole spectrum.\n")
    cat("Largest sampled BE multiplier error on this benchmark:",
        max(result$transfer$sampled_BE_multiplier_abs_error_mesh_benchmark), "\n")
    cat("Largest sampled BE multiplier change if k is omitted:",
        max(result$transfer$sampled_BE_change_if_k_dropped_mesh_benchmark), "\n\n")
    print(result$approximation, row.names = FALSE, digits = 5)
    cat("\nOff-grid interpolation check (beta=1/3):\n")
    print(interpolation, row.names = FALSE)
    cat("\nFull package table inventory (invalid means signs, duplicates or nonfinite values):\n")
    print(inventory, row.names = FALSE)
    cat("\nWL2:", result$wl2$status, "\n")
    if (!is.null(result$wl2$error)) cat(result$wl2$error, "\n")
    if (!is.null(result$wl2$comparison)) print(result$wl2$comparison, row.names = FALSE)
    cat("WL2 fits coefficients; it cannot retrieve tabulated BRASIL coefficients.\n")
    cat("WL2 covariance alpha=beta targets the same power; operator alpha=2*beta has m+1 poles.\n")
    cat("k_term=TRUE fits k on a finite interval; it does not select BRASIL.\n")
    cat("Approximation errors are sampled; relative tail errors can grow because k>0.\n")
    cat("This audit does not certify minimax optimality or rerun the full convergence study.\n")
  })
  cat(paste(report, collapse = "\n"), "\n")
  if (write_results) {
    directory <- file.path(project, "test/brasil_coefficient_audit")
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    for (name in c("inventory", "coefficients", "approximation", "transfer", "interpolation"))
      utils::write.csv(result[[name]], file.path(directory, paste0(name, ".csv")), row.names = FALSE)
    saveRDS(result, file.path(directory, "audit.rds"))
    writeLines(report, file.path(directory, "report.txt"))
    cat("Results written to", directory, "\n")
  }
  invisible(result)
}

brasil_audit <- audit_brasil_coefficients(
  run_wl2 = "--wl2" %in% commandArgs(trailingOnly = TRUE),
  write_results = "--write-results" %in% commandArgs(trailingOnly = TRUE))
