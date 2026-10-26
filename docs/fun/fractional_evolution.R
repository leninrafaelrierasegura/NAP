# remotes::install_github("davidbolin/rspde", ref = "devel")
# remotes::install_github("davidbolin/metricgraph", ref = "devel")
library(rSPDE)
library(MetricGraph)
library(grateful)

library(ggplot2)
library(reshape2)
library(plotly)

# Function to compute the roots and factor for the rational approximation
rational_supported_orders <- function(coefficient_source = "package") {
  coefficient_source <- match.arg(coefficient_source, c("package", "computed"))
  if (coefficient_source == "package") 1:4 else 1:8
}

read_computed_fractional_roots <- function() {
  path <- here::here("matlab", "computed_fractional_roots.csv")
  if (!file.exists(path)) stop("Missing computed coefficients: ", path)
  tab <- read.csv(path, stringsAsFactors = FALSE)
  if (!identical(names(tab), c("m", "beta", "kind", "index", "value")) ||
      !all(vapply(tab[c("m", "beta", "index", "value")], is.numeric, logical(1))) ||
      !all(is.finite(as.matrix(tab[c("m", "beta", "index", "value")]))) ||
      !all(tab$m %in% 1:8) || !all(tab$beta >= .25 & tab$beta < 1) ||
      !all(tab$kind %in% c("factor", "rc", "rb")) ||
      anyDuplicated(tab[c("m", "beta", "kind", "index")]))
    stop("Malformed computed coefficient artifact.")
  tab
}

validate_fractional_roots <- function(roots, m) {
  all_roots <- c(roots$pr_roots, roots$pl_roots)
  if (length(roots$pr_roots) != m || length(roots$pl_roots) != m+1L ||
      length(roots$factor) != 1L || !is.finite(roots$factor) || roots$factor <= 0 ||
      !is.numeric(all_roots) || is.complex(all_roots) ||
      any(!is.finite(all_roots)) || any(all_roots >= 0))
    stop("Fractional coefficients need a positive factor and real negative roots.")
  numerator <- sort(roots$pr_roots)
  denominator <- sort(roots$pl_roots)
  if (any(denominator[-length(denominator)] >= numerator) ||
      any(numerator >= denominator[-1L]))
    stop("Fractional coefficient roots do not interlace; this approximation is unsafe.")
  roots
}

my.get.roots <- function(m,
                         beta, # smoothness parameter, beta = alpha/2
                         interpolation = c("linear", "cubic"),
                         coefficient_source = "package"
                         ) {
  coefficient_source <- match.arg(coefficient_source, c("package", "computed"))
  interpolation <- match.arg(interpolation)
  orders <- rational_supported_orders(coefficient_source)
  if (length(m) != 1L || !is.finite(m) || !(m %in% orders))
    stop("Coefficient source '", coefficient_source, "' supports m = ",
         paste(orders, collapse = ", "), " only.")
  if (length(beta) != 1L || !is.finite(beta) || beta < 0.25 || beta >= 1)
    stop("Fractional evolution tables require 0.25 <= beta < 1; beta = 1 uses the exact solver.")
  if (coefficient_source == "computed") {
    tab <- read_computed_fractional_roots()
    rows <- tab[tab$m == m & abs(tab$beta-beta) <= 32*.Machine$double.eps*max(1, abs(beta)), ]
    if (!nrow(rows)) stop("Computed coefficients are unavailable at beta = ",
      format(beta, digits = 16), ". Available beta: ",
      paste(sort(unique(tab$beta)), collapse = ", "),
      ". Generate and validate new coefficients; beta interpolation is disabled.")
    values <- function(kind, indices) {
      part <- rows[rows$kind == kind, ]
      if (nrow(part) != length(indices) || !setequal(part$index, indices))
        stop("Computed coefficient entries are incomplete for m = ", m, ".")
      part$value[match(indices, part$index)]
    }
    return(validate_fractional_roots(list(factor = values("factor", 0),
      pr_roots = values("rc", seq_len(m)), pl_roots = values("rb", seq_len(m+1L))), m))
  }
  mt <- get(paste0("m", m, "table"), envir = asNamespace("rSPDE"))
  interpolate <- function(column) {
    if (interpolation == "linear") approx(mt$beta, mt[[column]], beta)$y else
      splinefun(mt$beta, mt[[column]], method = "fmm")(beta)
  }
  rb <- rep(0, m + 1)
  rc <- rep(0, m)
  if(m == 1) {
    rc = interpolate("rc")
  } else {
    rc = sapply(1:m, function(i) {
      interpolate(paste0("rc.", i))
    })
  }
  rb = sapply(1:(m+1), function(i) {
    interpolate(paste0("rb.", i))
  })
  factor = interpolate("factor")
  return(validate_fractional_roots(list(pl_roots = rb, # roots \{r_{2j}\}_{j=1}^{m+1}
              pr_roots = rc, # roots \{r_{1i}\}_{i=1}^m
              factor = factor # this is c_m/b_{m+1}
              ), m))
}

# Function to compute polynomial coefficients from roots
poly.from.roots <- function(roots) {
  coef <- 1
  for (r in roots) coef <- c(-r*coef, 0)+c(0, coef)
  return(coef) # Descending powers; the constant coefficient is one.
}

# Function to compute the parameters for the partial fraction decomposition
compute.partial.fraction.param <- function(factor, # c_m/b_{m+1}
                                           pr_roots, # roots \{r_{1i}\}_{i=1}^m
                                           pl_roots, # roots \{r_{2j}\}_{j=1}^{m+1}
                                           time_step, # \tau
                                           scaling # diffusion * kappa^(2*beta)
                                           ) {
  strength <- (scaling*time_step)/factor
  if (!is.finite(strength) || strength <= 0)
    stop("Invalid or underflowed scaled time step in the rational update.")
  pr_coef <- c(0, poly.from.roots(pr_roots))
  pl_coef <- poly.from.roots(pl_roots)
  factor_pr_coef <- pr_coef
  pr_plus_pl_coef <- factor_pr_coef + strength * pl_coef
  #res <- gsignal::residue(factor_pr_coef, pr_plus_pl_coef)
  # Coefficients passed here are in descending powers.
  # polyroot() expects ascending powers.
  p <- polyroot(rev(pr_plus_pl_coef))
  
  # Verify that these poles are numerically real and negative.
  stopifnot(
    all(abs(Im(p)) <= 1e-10 * pmax(1, abs(Re(p)))),
    all(Re(p) < 0)
  )
  p <- Re(p)
  
  derivative <- function(roots, x) sum(vapply(seq_along(roots), function(j)
    -roots[j]*prod(1-roots[-j]*x), numeric(1)))
  # At a pole P=-strength*Q: use that identity instead of evaluating the
  # nearly zero numerator. Products/derivatives retain the root representation.
  r <- vapply(p, function(x) {
    -strength*prod(1-pl_roots*x)/
      (derivative(pr_roots, x)+strength*derivative(pl_roots, x))
  }, numeric(1))
  if (length(p) != length(pl_roots) || any(!is.finite(r)) || any(r <= 0))
    stop("Unstable rational update: residues must be finite and positive.")
  # Guard against a misleading real-pole result with inaccurate reconstruction.
  z <- 10^seq(0, 12, length.out = 65L)
  P <- vapply(z, function(x) prod(1-pr_roots*x), numeric(1))
  Q <- vapply(z, function(x) prod(1-pl_roots*x), numeric(1))
  expected <- P/(P+strength*Q)
  reconstructed <- rowSums(vapply(seq_along(p), function(j) r[j]/(z-p[j]), numeric(length(z))))
  if (any(!is.finite(reconstructed)) || max(abs(reconstructed-expected)) > 5e-9)
    stop("Unstable rational update: partial-fraction reconstruction failed.")
  
  res <- list(r = r, p = p, k = numeric(0))
  return(list(r = res$r, # residues \{a_k\}_{k=1}^{m+1}
              p = res$p, # poles \{p_k\}_{k=1}^{m+1}
              k = res$k # remainder r
              ))
}

# Function to compute the fractional operator
my.fractional.operators.frac <- function(L, # Laplacian matrix
                                         beta, # smoothness parameter beta
                                         C, # mass matrix (not lumped)
                                         scale.factor, # scaling parameter = kappa^2
                                         m = 1, # order supported by the selected source
                                         time_step, # physical time step = tau
                                         diffusion = 1, # coefficient outside the fractional power
                                         root_interpolation = "linear",
                                         coefficient_source = "package"
                                         ) {
  coefficient_source <- match.arg(coefficient_source, c("package", "computed"))
  stopifnot(length(diffusion) == 1L, is.finite(diffusion), diffusion > 0,
            length(time_step) == 1L, is.finite(time_step), time_step > 0,
            length(beta) == 1L, is.finite(beta), beta > 0, beta <= 1,
            length(scale.factor) == 1L, is.finite(scale.factor), scale.factor > 0)
  L <- L / scale.factor
  if(beta == 1){
    L <- L * scale.factor^beta
    return(list(C = C, # mass matrix
                L = L, # Laplacian matrix scaled
                m = m, # rational order
                beta = beta, # smoothness parameter
                diffusion = diffusion,
                LHS = C + time_step * diffusion * L,
                factorization = Matrix::Cholesky(
                  Matrix::forceSymmetric(C + time_step * diffusion * L), LDL = FALSE)
                ))
  } else {
    scaling <- diffusion * scale.factor^beta
    roots <- my.get.roots(m, beta, interpolation = root_interpolation,
                          coefficient_source = coefficient_source)
    poles_rs_k <- compute.partial.fraction.param(factor = roots$factor,
                                                 pr_roots = roots$pr_roots,
                                                 pl_roots = roots$pl_roots,
                                                 time_step = time_step, 
                                                 scaling = scaling)

    partial_fraction_terms <- list()
    for (i in 1:(m+1)) {
      # Here is where the terms in the sum in eq 12 are computed
      partial_fraction_terms[[i]] <- (L - poles_rs_k$p[i] * C)#/poles_rs_k$r[i]
      }
    return(list(C = C, # mass matrix
                L = L, # Laplacian matrix scaled
                m = m, # rational order
                beta = beta, # smoothness parameter
                partial_fraction_terms = partial_fraction_terms, # partial fraction terms
                residues = poles_rs_k$r, # residues \{a_k\}_{k=1}^{m+1}
                poles = poles_rs_k$p,
                coefficient_source = coefficient_source,
                diffusion = diffusion,
                factorizations = lapply(partial_fraction_terms, function(A) {
                  Matrix::Cholesky(Matrix::forceSymmetric(A), LDL = FALSE)
                })
                ))
  }
}

# Function to solve the iteration
my.solver.frac <- function(obj, # object returned by my.fractional.operators.frac()
                           v # vector to be solved for
                           ){
  beta <- obj$beta
  m <- obj$m
  if (beta == 1){
    return(if (!is.null(obj$factorization))
      solve(obj$factorization, v, system = "A") else solve(obj$LHS, v)
           )
  } else {
    partial_fraction_terms <- obj$partial_fraction_terms
    residues <- obj$residues
    output <- v*0
    for (i in seq_len(m + 1L)) {
      solved <- if (!is.null(obj$factorizations)) {
        solve(obj$factorizations[[i]], v, system = "A")
      } else solve(partial_fraction_terms[[i]], v)
      output <- output + residues[i] * solved
    }
    return(output # solve the linear system using the partial fraction decomposition
           )
  }
}

# my.solver.frac <- function(obj, v) {
#   beta <- obj$beta
#   m <- obj$m
#   
#   if (beta == 1) {
#     return(solve(obj$LHS, v))
#   } 
#   
#   partial_fraction_terms <- obj$partial_fraction_terms
#   residues <- obj$residues
#   
#   # --- PARALLEL PART ---
#   results <- parallel::mclapply(
#     1:(m+1),
#     function(i) residues[i] * solve(partial_fraction_terms[[i]], v),
#     mc.cores = parallel::detectCores()  
#   )
#   
#   return(Reduce(`+`, results))
# }

solve_fractional_evolution <- function(my_op_frac, time_step, time_seq, val_at_0, RHST) {
  CC <- my_op_frac$C
  SOL <- matrix(NA, nrow = nrow(CC), ncol = length(time_seq))
  SOL[, 1] <- val_at_0
  for (k in 1:(length(time_seq) - 1)) {
    rhs <- CC %*% SOL[, k] + time_step * RHST[, k + 1]
    SOL[, k + 1] <- as.matrix(my.solver.frac(my_op_frac, rhs))
  }
  return(SOL)
}

estimate_error_memory <- function(
    n_reference_nodes, n_reference_samples, n_blocks,
    n_solver_steps = numeric(0), matrix_copies = 6,
    baseline_gib = 0) {

  positive_counts <- function(x) {
    is.numeric(x) && length(x) > 0L &&
      all(is.finite(x)) && all(x >= 1) && all(x == floor(x))
  }
  stopifnot(
    positive_counts(n_reference_nodes), length(n_reference_nodes) == 1L,
    positive_counts(n_reference_samples), length(n_reference_samples) == 1L,
    positive_counts(n_blocks),
    is.numeric(n_solver_steps),
    length(n_solver_steps) == 0L || positive_counts(n_solver_steps),
    is.numeric(matrix_copies), length(matrix_copies) == 1L,
    is.finite(matrix_copies), matrix_copies >= 1,
    is.numeric(baseline_gib), length(baseline_gib) == 1L,
    is.finite(baseline_gib), baseline_gib >= 0
  )

  # Promote counts before multiplication to avoid integer overflow.
  samples <- as.double(n_reference_samples)
  blocks <- pmin(as.double(n_blocks), samples)
  per_block <- ceiling(samples / blocks)
  gib <- 1024^3
  matrix_bytes <- 8 * as.double(n_reference_nodes) * per_block
  block_bytes <- matrix_copies * matrix_bytes + 128 * per_block
  time_bytes <- 8 * (samples + 1 +
    2 * sum(as.double(n_solver_steps) + 1))

  data.frame(
    requested_blocks = as.double(n_blocks),
    n_blocks = blocks,
    samples_per_block = per_block,
    matrix_gib = matrix_bytes / gib,
    block_workspace_gib = block_bytes / gib,
    time_vectors_gib = time_bytes / gib,
    baseline_gib = baseline_gib,
    estimated_gib = (block_bytes + time_bytes) / gib + baseline_gib
  )
}

solve_fractional_evolution_error <- function(
    my_op_frac, time_step, time_seq, val_at_0,
    spatial_load, g_average, Psi, reference_C,
    reference_times, exact_at_time, n_blocks = 1000L,
    verbose = FALSE) {

  stopifnot(
    length(time_step) == 1L,
    is.finite(time_step), time_step > 0,
    length(time_seq) >= 2L,
    length(reference_times) >= 2L,
    all(is.finite(time_seq)),
    all(is.finite(reference_times)),
    all(diff(time_seq) > 0),
    all(diff(reference_times) > 0),
    reference_times[1] == time_seq[1],
    tail(reference_times, 1) == tail(time_seq, 1),
    length(g_average) == length(time_seq),
    length(val_at_0) == nrow(my_op_frac$C),
    length(spatial_load) == length(val_at_0),
    ncol(Psi) == length(val_at_0),
    nrow(Psi) == nrow(reference_C),
    nrow(reference_C) == ncol(reference_C),
    is.function(exact_at_time),
    length(n_blocks) == 1L,
    is.finite(n_blocks), n_blocks >= 1,
    n_blocks == floor(n_blocks),
    is.logical(verbose), length(verbose) == 1L, !is.na(verbose)
  )

  n_samples <- length(reference_times) - 1L
  n_blocks <- min(n_blocks, n_samples)
  samples_per_block <- n_samples / n_blocks
  n_steps <- length(time_seq) - 1L
  u <- as.numeric(val_at_0)
  spatial_load <- as.numeric(spatial_load)
  solver_step <- 0L
  projected_step <- -1L
  u_ref <- NULL
  next_ref <- 2L
  error_sum <- 0
  started <- proc.time()[["elapsed"]]

  if (verbose) {
    message(sprintf(
      "Error accumulation: %s reference samples in %s blocks; up to %s samples/block (%.2f GiB per dense matrix).",
      format(n_samples, big.mark = ",", scientific = FALSE, trim = TRUE),
      format(n_blocks, big.mark = ",", scientific = FALSE, trim = TRUE),
      format(ceiling(samples_per_block), big.mark = ",", scientific = FALSE, trim = TRUE),
      8 * nrow(reference_C) * ceiling(samples_per_block) / 1024^3
    ))
  }

  for (block in seq_len(n_blocks)) {
    last_ref <- if (block == n_blocks) {
      length(reference_times)
    } else {
      1L + floor(block * samples_per_block)
    }
    ids <- seq.int(next_ref, last_ref)
    times <- reference_times[ids]
    target_steps <- findInterval(times, time_seq, left.open = TRUE)
    runs <- rle(target_steps)

    error <- as.matrix(exact_at_time(times))
    stopifnot(
      nrow(error) == nrow(reference_C),
      ncol(error) == length(ids)
    )

    first_column <- 1L
    for (run in seq_along(runs$values)) {
      target <- runs$values[run]
      while (solver_step < target) {
        solver_step <- solver_step + 1L
        rhs <- my_op_frac$C %*% u +
          time_step * spatial_load * g_average[solver_step + 1L]
        u <- as.numeric(my.solver.frac(my_op_frac, rhs))
      }

      if (projected_step != target) {
        u_ref <- as.numeric(Psi %*% u)
        projected_step <- target
      }

      last_column <- first_column + runs$lengths[run] - 1L
      columns <- seq.int(first_column, last_column)
      error[, columns] <- error[, columns, drop = FALSE] - u_ref
      first_column <- last_column + 1L
    }

    squared <- colSums(error * (reference_C %*% error))
    weights <- times - reference_times[ids - 1L]
    error_sum <- error_sum + sum(weights * squared)
    rm(error)
    next_ref <- last_ref + 1L

    if (verbose) {
      processed <- last_ref - 1L
      elapsed <- proc.time()[["elapsed"]] - started
      remaining <- elapsed * (n_samples / processed - 1)
      message(sprintf(
        "Block %s/%s | %.1f%% | %s/%s samples | solver step %s/%s | elapsed %.1f min | ETA %.1f min",
        block, n_blocks, 100 * processed / n_samples,
        format(processed, big.mark = ",", scientific = FALSE, trim = TRUE),
        format(n_samples, big.mark = ",", scientific = FALSE, trim = TRUE),
        solver_step, n_steps, elapsed / 60, remaining / 60
      ))
    }
  }

  stopifnot(next_ref == length(reference_times) + 1L,
            solver_step == n_steps)
  sqrt(error_sum)
}

# Function to construct a piecewise constant projection of approximated values
# construct_piecewise_projection <- function(projected_U_approx, time_seq, overkill_time_seq) {
#   projected_U_piecewise <- matrix(NA, nrow = nrow(projected_U_approx), ncol = length(overkill_time_seq))
  
#   # Assign value at t = 0
#   projected_U_piecewise[, which(overkill_time_seq == 0)] <- projected_U_approx[, 1]
  
#   # Assign values for intervals (t_{k-1}, t_k]
#   for (k in 2:length(time_seq)) {
#     idxs <- which(overkill_time_seq > time_seq[k - 1] & overkill_time_seq <= time_seq[k])
#     projected_U_piecewise[, idxs] <- projected_U_approx[, k]
#   }
  
#   return(projected_U_piecewise)
# }
construct_piecewise_projection <- function(
    projected_U_approx, time_seq, overkill_time_seq) {

  stopifnot(
    length(time_seq) >= 2L,
    ncol(projected_U_approx) == length(time_seq),
    all(is.finite(time_seq)),
    all(diff(time_seq) > 0),
    all(is.finite(overkill_time_seq)),
    all(overkill_time_seq >= time_seq[1]),
    all(overkill_time_seq <= tail(time_seq, 1))
  )

  columns <- findInterval(
    overkill_time_seq,
    time_seq,
    left.open = TRUE
  ) + 1L

  projected_U_approx[, columns, drop = FALSE]
}

make_time_grid <- function(T_final, target_step) {
  stopifnot(
    length(T_final) == 1L, is.finite(T_final), T_final > 0,
    length(target_step) == 1L,
    is.finite(target_step), target_step > 0
  )
  N <- ceiling(T_final / target_step)
  time_step <- T_final / N
  time_seq <- seq(0, T_final, length.out = N + 1)
  
  return(list(N = N, time_step = time_step, time_seq = time_seq))
}

average_sine_forcing <- function(time_seq, A = 1, omega = pi) {
  dt <- diff(time_seq)

  stopifnot(
    length(time_seq) >= 2L,
    all(is.finite(time_seq)),
    all(dt > 0)
  )

  midpoint <- head(time_seq, -1) + dt/2
  z <- omega * dt/2

  # sin(z)/z, with its limiting value at z = 0.
  sinc_z <- rep(1, length(z))
  nonzero <- z != 0
  sinc_z[nonzero] <- sin(z[nonzero]) / z[nonzero]

  averages <- A * sin(omega * midpoint) * sinc_z

  # First entry is unused by the evolution solver.
  c(0, averages)
}
