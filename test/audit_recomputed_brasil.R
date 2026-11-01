# Recompute/audit coefficients without editing the study or installed packages.
# From NAP:
# /usr/bin/Rscript --vanilla test/audit_recomputed_brasil.R --prepare
# <Python BRASIL generation; see test/recomputed_brasil_audit/README.md>
# /usr/bin/Rscript --vanilla test/audit_recomputed_brasil.R --compare
# Outputs stay in test/recomputed_brasil_audit/.

library(Matrix)
audit_dir <- "test/recomputed_brasil_audit"
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
helper <- new.env(parent = globalenv())
sys.source("fun/isolated_rational_convergence.R", helper)

read_settings <- function() {
  text <- readLines("conv_in_m_isolated.Rmd")
  first <- which(grepl("^```\\{r settings\\}", text))
  last <- first + which(grepl("^```$", text[(first+1L):length(text)]))[1L]
  settings <- new.env(parent = baseenv())
  eval(parse(text = text[(first+1L):(last-1L)]), settings)
  as.list(settings)
}

# The installed MetricGraph 1.1.2 lacks the notebook's newer `edges` constructor
# and set_manual_edge_lengths(). In that case use its SAME P1 FEM assembler on
# the prescribed tadpole topology and uniform subdivisions, with a thin graph
# interface for the notebook's analytic loads. No continuum modal replacement.
adapt_graph_api <- function() {
  if ("edges" %in% names(formals(MetricGraph::metric_graph$public_methods$initialize)))
    return("native notebook MetricGraph API")
  helper$gets.graph.tadpole <- function(h) {
    counts <- ceiling(c(1, 2)/h)
    n <- as.integer(2L + sum(counts-1L))
    paths <- list(c(1L, seq.int(3L, counts[1]+1L), 2L),
                  c(2L, seq.int(counts[1]+2L, n), 2L))
    E <- do.call(rbind, lapply(paths, function(nodes) cbind(head(nodes,-1L),tail(nodes,-1L))))
    lengths <- rep(c(1, 2)/counts, counts)
    fem <- get("assemble_fem", asNamespace("MetricGraph"))(E, lengths, n)
    PtE <- do.call(rbind, lapply(1:2, function(edge)
      cbind(edge, seq_len(counts[edge]-1L)/counts[edge])))
    list(mesh = c(fem, list(PtE = PtE)), compute_fem = function() invisible(NULL),
         get_edge_lengths = function() c(1, 2))
  }
  "MetricGraph::assemble_fem compatibility path; same topology, lengths, P1 mesh and loads"
}

scalar_diagnostics <- function(cf, beta, m, eigenvalues, kappa) {
  grids <- list(common = exp(seq(0, log(1e3), length.out = 20001L)),
    spectrum = eigenvalues/kappa^2,
    moving = 10^seq(0, (5+m)/2, length.out = 20001L))
  vapply(grids, function(y) max(abs(helper$brasil_inverse(y, cf)-y^(-beta))), numeric(1))
}

if ("--prepare" %in% commandArgs(TRUE)) {
  settings <- read_settings()
  compatibility <- adapt_graph_api()
  data <- helper$isolated_tadpole_data(settings$h, settings$kappa)
  spectrum <- helper$isolated_discrete_spectrum(data$C, data$L)
  N <- round(settings$T_final/settings$tau)
  stopifnot(abs(N*settings$tau-settings$T_final) < 1e-12)
  times <- (0:N)*settings$tau
  references <- lapply(settings$alpha_vector, function(alpha)
    helper$isolated_exact_be_reference(data$C, spectrum, data$initial,
      data$spatial_load, times, alpha, settings$diffusion))
  control_reference <- helper$isolated_exact_be_reference(data$C, spectrum,
    data$initial, data$spatial_load, times, 2, settings$diffusion)
  control <- helper$isolated_rational_error(data, spectrum, control_reference, 1L)
  rows <- list()
  for (reference in references) for (m in seq_len(settings$m_max)) {
    row <- helper$isolated_rational_error(data, spectrum, reference, m)
    cf <- helper$brasil_coefficients(m, reference$alpha/2)
    scalar <- scalar_diagnostics(cf, reference$alpha/2, m, spectrum$values, settings$kappa)
    row$family <- "table"
    row$scalar_common <- scalar["common"]
    row$scalar_spectrum <- scalar["spectrum"]
    row$scalar_moving <- scalar["moving"]
    rows[[length(rows)+1L]] <- row
  }
  baseline <- do.call(rbind, rows)
  utils::write.csv(baseline, file.path(audit_dir, "baseline.csv"), row.names = FALSE)
  # The finite interval covers BOTH the actual entire spectrum and the scalar
  # [1,1000] diagnostic. Its domain stays fixed across m.
  xmin <- min(1e-3, settings$kappa^2/max(spectrum$values))
  jobs <- expand.grid(family = c("same_interval", "spectrum_interval"),
    beta = settings$alpha_vector/2, m = seq_len(settings$m_max), stringsAsFactors = FALSE)
  jobs$x_min <- ifelse(jobs$family == "same_interval", 10^(-(5+jobs$m)/2), xmin)
  jobs <- rbind(jobs, data.frame(family = "repair_example", beta = c(.95, .85),
    m = c(15L, 20L), x_min = 10^(-(5+c(15L,20L))/2)))
  jobs$beta <- sprintf("%.12g", jobs$beta)
  jobs$x_min <- sprintf("%.17g", jobs$x_min)
  utils::write.csv(jobs, file.path(audit_dir, "jobs.csv"), row.names = FALSE)
  context <- list(settings = settings, compatibility = compatibility, data = data,
    spectrum = spectrum, references = references, baseline = baseline, control = control,
    xmin = xmin, R = R.version.string, session = sessionInfo(),
    source_md5 = tools::md5sum(c("conv_in_m_isolated.Rmd", "fun/isolated_rational_convergence.R",
      "fun/brasil_coefficients.R")))
  saveRDS(context, file.path(audit_dir, "context.rds"))
  cat("Prepared", nrow(jobs), "BRASIL fits; FEM dimension", nrow(data$C), "\n")
  cat("Normalized eigenvalue range:", range(spectrum$values/settings$kappa^2), "\n")
  cat("Fixed fitting interval: [", xmin, ",1]\n")
  cat("Classical BE control error:", control$error, "\n")
  cat(compatibility, "\n")
}

if ("--compare" %in% commandArgs(TRUE)) {
  context <- readRDS(file.path(audit_dir, "context.rds"))
  stopifnot(identical(unname(context$source_md5),
    unname(tools::md5sum(names(context$source_md5)))))
  coefficients <- read.csv(file.path(audit_dir, "coefficients.csv"))
  fits <- read.csv(file.path(audit_dir, "fit_diagnostics.csv"))
  stopifnot(all(fits$converged %in% c(TRUE,"True")),
    all(fits$valid_positive_structure %in% c(TRUE,"True")),
    max(fits$mp_conversion_gap) < 1e-40,
    max(fits$double_conversion_gap) < 1e-12)
  original <- helper$brasil_coefficients
  data <- context$data
  spectrum <- context$spectrum
  settings <- context$settings
  rows <- list(context$baseline)
  for (family in c("same_interval", "spectrum_interval")) {
    helper$brasil_coefficients <- local({
      chosen <- family
      function(m, beta) {
        part <- coefficients[coefficients$family == chosen & coefficients$m == m &
          abs(coefficients$beta-beta) < 1e-12, ]
        part <- part[order(part$pole), ]
        stopifnot(nrow(part) == m, all(part$pole < 0), all(part$residue > 0),
                  all(part$constant >= 0), !anyDuplicated(part$pole))
        list(poles = part$pole, residues = part$residue, constant = part$constant[1L],
             m = m, beta = beta)
      }
    })
    for (reference in context$references) for (m in seq_len(settings$m_max)) {
      row <- helper$isolated_rational_error(data, spectrum, reference, m)
      cf <- helper$brasil_coefficients(m, reference$alpha/2)
      scalar <- scalar_diagnostics(cf, reference$alpha/2, m, spectrum$values, settings$kappa)
      row$family <- family
      row$scalar_common <- scalar["common"]
      row$scalar_spectrum <- scalar["spectrum"]
      row$scalar_moving <- scalar["moving"]
      rows[[length(rows)+1L]] <- row
    }
  }
  helper$brasil_coefficients <- original
  results <- do.call(rbind, rows)
  results$noise_indicator <- pmax(results$absolute_error_disagreement, context$control$error)
  results$resolved <- results$product_error > 100*results$noise_indicator
  fit_slope <- function(d) if (nrow(d) >= 2L)
    unname(coef(lm(log(error) ~ sqrt(m), d))[2L]) else NA_real_
  slopes <- do.call(rbind, lapply(split(results, list(results$family,results$alpha)), function(d) {
    resolved <- d[d$resolved, ]
    data.frame(family = d$family[1L], alpha = d$alpha[1L],
      all_orders_slope = fit_slope(d), resolved_only_slope = fit_slope(resolved),
      resolved_orders = nrow(resolved), error_m4 = d$error[d$m==4L],
      error_m8 = d$error[d$m==8L], error_m14 = d$error[d$m==14L])
  }))
  # Matched resolved comparisons prevent floor-contaminated slopes from being
  # mistaken for an improvement (or deterioration) in convergence.
  matched <- do.call(rbind, lapply(settings$alpha_vector, function(alpha) {
    old <- results[results$alpha == alpha & results$family == "table", ]
    fresh <- results[results$alpha == alpha & results$family == "spectrum_interval", ]
    use <- old$m %in% fresh$m[fresh$resolved] & old$resolved
    orders <- old$m[use]
    data.frame(alpha = alpha, matched_orders = paste(orders, collapse = ","),
      table_slope = fit_slope(old[use, ]),
      spectrum_fit_slope = fit_slope(fresh[fresh$m %in% orders, ]),
      improvement_m4 = old$error[old$m==4L]/fresh$error[fresh$m==4L],
      improvement_m8 = old$error[old$m==8L]/fresh$error[fresh$m==8L])
  }))
  repairs <- do.call(rbind, lapply(list(c(15,.95), c(20,.85)), function(case) {
    m <- case[1L]
    beta <- case[2L]
    table <- get(paste0("m_brasil",m,"t"), asNamespace("rSPDE"))
    at <- which.min(abs(table$alpha-beta))
    old <- list(constant = table$k[at],
      poles = as.numeric(table[at,paste0("p",seq_len(m))]),
      residues = as.numeric(table[at,paste0("r",seq_len(m))]))
    part <- coefficients[coefficients$family == "repair_example" &
      coefficients$m == m & abs(coefficients$beta-beta)<1e-12, ]
    fresh <- list(constant = part$constant[1L], poles = part$pole, residues = part$residue)
    old_error <- scalar_diagnostics(old,beta,m,spectrum$values,settings$kappa)
    fresh_error <- scalar_diagnostics(fresh,beta,m,spectrum$values,settings$kappa)
    data.frame(m=m,beta=beta,table_constant=old$constant,fresh_constant=fresh$constant,
      table_min_residue=min(old$residues),fresh_min_residue=min(fresh$residues),
      table_scalar_moving=old_error["moving"],fresh_scalar_moving=fresh_error["moving"],
      table_scalar_spectrum=old_error["spectrum"],fresh_scalar_spectrum=fresh_error["spectrum"])
  }))
  utils::write.csv(results, file.path(audit_dir, "comparison.csv"), row.names = FALSE)
  utils::write.csv(slopes, file.path(audit_dir, "slopes.csv"), row.names = FALSE)
  utils::write.csv(matched, file.path(audit_dir, "matched_slopes.csv"), row.names = FALSE)
  utils::write.csv(repairs, file.path(audit_dir, "repaired_examples.csv"), row.names = FALSE)
  library(ggplot2)
  plot <- ggplot(results, aes(sqrt(m), error, colour = family, group = family)) +
    geom_line() + geom_point(aes(shape = resolved), size = 1.4) +
    facet_wrap(~alpha, scales = "free_y", ncol = 4) +
    scale_y_log10() + scale_shape_manual(values = c(`FALSE`=1,`TRUE`=16)) +
    labs(x = "sqrt(m)", y = "Isolated space-time error", colour = "Coefficients",
      title = "Fresh BRASIL coefficients on the same FEM/backward Euler problem",
      subtitle = "Open points: below the notebook's empirical resolution indicator") +
    theme_minimal() + theme(legend.position = "bottom")
  ggsave(file.path(audit_dir, "convergence_comparison.png"), plot, width = 13, height = 7, dpi = 160)
  report <- capture.output({
    cat("Recomputed BRASIL audit of conv_in_m_isolated.Rmd\n")
    cat("Settings: h=",settings$h," tau=",settings$tau," m_max=",settings$m_max,"\n")
    cat(context$compatibility, "\n")
    print(spectrum$diagnostics, row.names = FALSE)
    cat("Fixed fit interval [",context$xmin,",1] covers all discrete eigenvalues and [0.001,1].\n")
    cat("All",nrow(fits),"fits converged and have distinct negative poles, positive residues, k>=0.\n")
    cat("Maximum MP barycentric/partial-fraction conversion gap:", max(fits$mp_conversion_gap), "\n")
    cat("Maximum gap after conversion to double:", max(fits$double_conversion_gap), "\n")
    cat("Maximum sparse/eigenbasis error disagreement:", max(results$absolute_error_disagreement), "\n")
    cat("Classical BE control:",context$control$error,"\n\n")
    print(slopes, row.names = FALSE, digits = 6)
    cat("\nSlopes over identical resolved orders, table versus spectrum-specific fitting:\n")
    print(matched, row.names = FALSE, digits = 6)
    cat("\nRepair-example fits (outside the notebook's m<=14 range):\n")
    print(repairs, row.names = FALSE, digits = 8)
    cat("\nSame-interval regeneration tests precision alone. Spectrum-interval fitting changes the approximation domain.\n")
    cat("Spectrum-specific coefficients must be checked/refitted if h or kappa changes.\n")
    cat("No claim about a mesh-uniform asymptotic rate follows from fixed-spectrum results.\n")
    cat("High-order points at the numerical floor must not be used to estimate a convergence exponent.\n")
  })
  writeLines(report, file.path(audit_dir, "report.txt"))
  cat(paste(report, collapse="\n"), "\n")
  saveRDS(list(context=context, results=results, slopes=slopes, matched=matched, fits=fits),
          file.path(audit_dir, "comparison.rds"))
}
