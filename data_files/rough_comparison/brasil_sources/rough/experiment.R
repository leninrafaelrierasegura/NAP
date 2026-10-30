rough_parameters <- function() {
  list(J = 5:9, alpha = seq(0.5, 1.9, 0.2),
       time_steps = c(5e-4, 1.25e-4, 1/30000, 1e-5, 1/320000, 1e-6, 4e-7, 2e-7),
       rational_order = 14L, symmetry_reduction = TRUE,
       T_final = 1, diffusion = 1, kappa = 4,
       projection_quadrature = 16L, tail_tolerance = 1e-13,
       fit_J = 7:9, forcing = "sum_{j>=1} sin(2^j*pi*x)/j on loop; zero on stem",
       reconstruction = "right-constant backward Euler",
       error_integration = "analytic time moments, spatial projection residual, bounded infinite tail")
}

rough_run_experiment <- function(parameters = rough_parameters(),
    output_file = file.path("data_files", "rough_h_results.rds"),
    checkpoint_dir = tempfile("rough-checkpoints-"), workers = 1L,
    root = getwd(), cache_dir = tempdir()) {
  stopifnot(length(parameters$alpha) == length(parameters$time_steps),
            all(parameters$fit_J %in% parameters$J), workers >= 1L,
            workers == as.integer(workers))
  dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  rough_compile_solver(root, cache_dir)
  cases <- expand.grid(J = parameters$J, alpha_index = seq_along(parameters$alpha))
  identity <- list(parameters = parameters,
    sources = tools::md5sum(list.files(file.path(root, "fun/rough"),
                                      pattern = "\\.(R|cpp)$", full.names = TRUE)))
  run_case <- function(i) {
    J <- cases$J[i]
    a <- cases$alpha_index[i]
    alpha <- parameters$alpha[a]
    checkpoint <- file.path(checkpoint_dir, sprintf("J%02d-alpha%03d.rds", J, round(100*alpha)))
    if (file.exists(checkpoint)) {
      saved <- readRDS(checkpoint)
      if (!identical(saved$identity, identity)) stop("Checkpoint source/settings mismatch: ", checkpoint)
      return(saved$value)
    }
    message(sprintf("Matrix experiment: J=%s, h=%g, alpha=%g", J, 2^-J, alpha))
    mesh <- rough_mesh(J, parameters$kappa)
    projection <- rough_projection(mesh, parameters$projection_quadrature)
    value <- rough_solve(mesh, alpha, parameters$time_steps[a],
      order = parameters$rational_order, diffusion = parameters$diffusion,
      T_final = parameters$T_final, projection = projection,
      tail_tolerance = parameters$tail_tolerance,
      symmetry_reduction = parameters$symmetry_reduction)
    value$final_state <- value$states <- NULL
    value$J <- J
    value$alpha <- alpha
    value$n_nodes <- nrow(mesh$C)
    saveRDS(list(identity = identity, value = value), checkpoint)
    message(sprintf("Finished J=%s alpha=%g: E=%.9g (%.2f s)", J, alpha, value$error, value$elapsed))
    value
  }
  safe_case <- function(i) tryCatch(run_case(i), error = function(e)
    list(failed = TRUE, message = conditionMessage(e), case = i))
  values <- if (workers > 1L && .Platform$OS.type == "unix")
    parallel::mclapply(seq_len(nrow(cases)), safe_case,
                       mc.cores = min(workers, nrow(cases)), mc.preschedule = FALSE) else
    lapply(seq_len(nrow(cases)), safe_case)
  failed <- vapply(values, function(x) !is.list(x) || isTRUE(x$failed), logical(1))
  if (any(failed)) stop("Rough experiment failed: ", paste(vapply(values[failed],
    function(x) if (is.list(x)) x$message else as.character(x), character(1)), collapse = "; "))
  errors <- matrix(vapply(values, `[[`, numeric(1), "error"),
                    nrow = length(parameters$J), ncol = length(parameters$alpha))
  corrected <- errors * parameters$J
  fit <- which(parameters$J %in% parameters$fit_J)
  slope <- function(y) unname(coef(stats::lm(log10(y[fit]) ~ log10(2^-parameters$J[fit])))[2L])
  rates <- data.frame(alpha = parameters$alpha,
    theoretical = parameters$alpha, raw = apply(errors, 2L, slope),
    corrected = apply(corrected, 2L, slope))
  rates$relative_error_percent <- 100*abs(rates$corrected-rates$alpha)/rates$alpha
  result <- list(parameters = parameters, errors = errors, corrected_errors = corrected,
    rates = rates, cases = values, provenance = identity,
    generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z", tz = "Asia/Riyadh"),
    session_info = utils::sessionInfo())
  saveRDS(result, output_file)
  result
}
