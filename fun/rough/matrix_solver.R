# This is the same backward Euler update as matrix_reference:
# U[k+1] = H(C^-1 K / kappa^2) C^-1 (C U[k] + tau b).
# The operator tables use the same root/factor interpolation as my.get.roots()
# in fun/matrix_reference. BRASIL remains available for baseline checks.
rough_rational_coefficients <- function(beta, order = 8L,
    coefficient_source = c("computed", "operator", "brasil")) {
  coefficient_source <- match.arg(coefficient_source)
  if (coefficient_source %in% c("computed","operator")) {
    stopifnot(length(beta)==1L,is.finite(beta),beta>0,beta<1,length(order)==1L)
    if (coefficient_source=="computed") {
      stopifnot(order %in% 1:8)
      table_name <- file.path("matlab","computed_fractional_roots.csv")
      table <- utils::read.csv(table_name,stringsAsFactors=FALSE)
      rows <- table[table$m==order & abs(table$beta-beta)<32*.Machine$double.eps,]
      if (!nrow(rows)) stop("No directly computed coefficient row at beta = ",beta)
      values <- function(kind,indices) {
        row <- rows[rows$kind==kind,]
        if (nrow(row)!=length(indices) || !setequal(row$index,indices))
          stop("Missing or duplicate computed coefficient entries.")
        row$value[match(indices,row$index)]
      }
      rc <- values("rc",seq_len(order))
      rb <- values("rb",seq_len(order+1L))
      factor <- values("factor",0)
      interpolation <- "none"
    } else {
      stopifnot(order %in% 1:4)
      table_name <- paste0("m",order,"table")
      table <- get(table_name,envir=asNamespace("rSPDE"))
      if (beta < min(table$beta) || beta > max(table$beta))
        stop("Requested beta lies outside the operator coefficient table.")
      value <- function(name) stats::approx(table$beta,table[[name]],xout=beta)$y
      rc <- if (order==1L && !"rc.1" %in% names(table)) value("rc") else
        vapply(paste0("rc.",seq_len(order)),value,numeric(1))
      rb <- vapply(paste0("rb.",seq_len(order+1L)),value,numeric(1))
      factor <- value("factor")
      interpolation <- "linear"
    }
    stopifnot(all(is.finite(c(rc,rb,factor))),all(rc<0),all(rb<0),factor>0)
    # Algebraic partial fractions of the SAME root/factor inverse function:
    # factor*prod(1-rc*x)/prod(1-rb*x) = sum residues/(x-poles).
    # This is a conversion of representation, not a different approximation.
    poles <- 1/rb
    residues <- vapply(seq_along(poles),function(i) {
      x <- poles[i]
      factor*prod(1-rc*x)/(-rb[i]*prod(1-rb[-i]*x))
    },numeric(1))
    stopifnot(all(is.finite(residues)),all(residues>0))
    index <- order(poles)
    coefficients <- list(poles=poles[index],residues=residues[index],constant=0,
      order=as.integer(order),beta=beta,coefficient_source=coefficient_source,table=table_name,
      roots=list(pr_roots=rc,pl_roots=rb,factor=factor),interpolation=interpolation)
    if (coefficient_source=="computed") coefficients$table_md5 <- unname(tools::md5sum(table_name))
    x <- 10^seq(0,12,length.out=257L)
    expected <- vapply(x,function(z)factor*prod(1-rc*z)/prod(1-rb*z),numeric(1))
    actual <- vapply(x,function(z)sum(coefficients$residues/(z-coefficients$poles)),numeric(1))
    discrepancy <- max(abs(actual/expected-1))
    if (!is.finite(discrepancy) || discrepancy>2e-12)
      stop("Root/factor to partial-fraction conversion failed.")
    coefficients$inverse_reconstruction_error <- discrepancy
    return(coefficients)
  }
  stopifnot(beta > 0, beta < 1, order %in% 1:20)
  table <- get(paste0("m_brasil", order, "t"), envir = asNamespace("rSPDE"))
  row <- which(abs(table$alpha-beta) < 32*.Machine$double.eps)
  if (length(row) != 1L) stop("No exact BRASIL coefficient row at beta = ", beta)
  poles <- as.numeric(table[row, paste0("p", seq_len(order))])
  residues <- as.numeric(table[row, paste0("r", seq_len(order))])
  constant <- table$k[row]
  stopifnot(all(is.finite(c(poles, residues, constant))), all(poles < 0),
            all(residues > 0), constant >= 0)
  index <- order(poles)
  list(poles = poles[index], residues = residues[index], constant = constant,
       order = as.integer(order), beta = beta,coefficient_source="brasil",
       table=paste0("m_brasil",order,"t"))
}

rough_transfer <- function(coefficients, tau, diffusion = 1, kappa = 4) {
  p <- coefficients$poles
  r <- coefficients$residues
  strength <- tau * diffusion * kappa^(2*coefficients$beta)
  stopifnot(is.finite(strength), strength > 0)
  constant <- coefficients$constant
  # R(x) = constant + sum r_i/(x-p_i); H(x) = R(x)/(R(x)+strength).
  # R is strictly decreasing between poles. Exactly one H pole lies in
  # each interval (-Inf,p_1), (p_1,p_2), ... . Root coordinates normalized
  # to each interval keep uniroot accurate across many orders of magnitude.
  poles <- numeric(length(p))
  for (i in seq_along(p)) {
    upper <- p[i]
    lower <- if (i == 1L) 2*upper else p[i-1L]
    f <- function(x) constant + strength + sum(r/(x-p))
    if (i == 1L) while (f(lower) <= 0) lower <- 2*lower
    width <- upper-lower
    normalized <- function(z) {
      if (z == 0 && i > 1L) return(.Machine$double.xmax)
      if (z == 1) return(-.Machine$double.xmax)
      f(lower + width*z)
    }
    root <- stats::uniroot(normalized, c(0, 1), tol = 2*.Machine$double.eps,
                            maxiter = 1000L)$root
    poles[i] <- lower + width*root
  }
  residues <- vapply(poles, function(z) strength/sum(r/(z-p)^2), numeric(1))
  remainder <- constant/(constant+strength)
  stopifnot(all(is.finite(c(poles,residues,remainder))), all(poles < 0),
            all(residues > 0), remainder >= 0, remainder < 1)
  x <- 10^seq(0, 10, length.out = 257L)
  inverse <- vapply(x, function(z) constant + sum(r/(z-p)), numeric(1))
  expected <- inverse/(inverse+strength)
  actual <- remainder + vapply(x, function(z) sum(residues/(z-poles)), numeric(1))
  reconstruction_error <- max(abs(actual-expected))
  if (reconstruction_error > 2e-12)
    stop("Inaccurate partial fractions: ", format(reconstruction_error, digits = 5))
  list(poles = poles, residues = residues, remainder = remainder,
       reconstruction_error = reconstruction_error, coefficients = coefficients)
}

rough_operator <- function(mesh, alpha, tau, order = 8L, diffusion = 1,
                           coefficient_source = "computed") {
  transfer <- rough_transfer(rough_rational_coefficients(alpha/2, order,coefficient_source),
                             tau, diffusion, mesh$kappa)
  A <- lapply(transfer$poles, function(p) mesh$K/mesh$kappa^2 - p*mesh$C)
  # Compatible with my.solver.frac() in fun/matrix_reference. The constant
  # term is additional because these inverse approximants have nonzero limit.
  list(C = mesh$C, partial_fraction_terms = A,
       factorizations = lapply(A, function(A)
         Matrix::Cholesky(Matrix::forceSymmetric(A), LDL = FALSE)),
       mass_factor = Matrix::Cholesky(Matrix::forceSymmetric(mesh$C), LDL = FALSE),
       beta = alpha/2, m = length(A)-1L, residues = transfer$residues,
       remainder = transfer$remainder, transfer = transfer, order = order)
}

rough_matrix_step <- function(operator, rhs) {
  output <- operator$remainder * as.numeric(solve(operator$mass_factor, rhs, system = "A"))
  for (i in seq_along(operator$residues)) output <- output + operator$residues[i] *
    as.numeric(solve(operator$factorizations[[i]], rhs, system = "A"))
  output
}

rough_compile_solver <- function(root = getwd(), cache_dir = tempdir()) {
  if (!exists("rough_matrix_evolve_cpp", mode = "function", inherits = TRUE))
    Rcpp::sourceCpp(file.path(root, "fun/rough/matrix_evolution.cpp"),
                    cacheDir = cache_dir, showOutput = FALSE, verbose = FALSE)
  invisible(TRUE)
}

# Four signed reflections of the quarter-loop (0,1/2), with zero stem.
# This is an exact invariant FEM subspace for every forcing mode sin(2^j*pi*x).
# It groups nodal coefficients; it does not use discrete eigenfunctions.
rough_symmetry_embedding <- function(mesh) {
  n <- as.integer(1/mesh$h)
  quarter <- n/2L
  r <- seq_len(quarter-1L)
  rows <- mesh$loop_nodes[1L+c(r,n-r,n+r,2L*n-r)]
  Matrix::sparseMatrix(i = rows, j = rep(r,4L),
    x = rep(c(1,-1,1,-1),each=length(r)),
    dims = c(nrow(mesh$C),length(r)))
}

rough_solve <- function(mesh, alpha, time_step, order = 8L, diffusion = 1,
                        T_final = 1, projection = rough_projection(mesh),
                        tail_tolerance = 1e-13, return_states = FALSE,
                        symmetry_reduction = TRUE,coefficient_source = "computed") {
  stopifnot(is.finite(time_step), time_step > 0, T_final > 0)
  N <- ceiling(T_final/time_step)
  tau <- T_final/N
  solver_mesh <- mesh
  solver_projection <- projection$P
  embedding <- NULL
  if (symmetry_reduction) {
    embedding <- rough_symmetry_embedding(mesh)
    solver_mesh$C <- Matrix::crossprod(embedding, mesh$C %*% embedding)
    solver_mesh$K <- Matrix::crossprod(embedding, mesh$K %*% embedding)
    solver_mesh$load <- as.numeric(Matrix::crossprod(embedding,mesh$load))
    factor <- Matrix::Cholesky(Matrix::forceSymmetric(solver_mesh$C),LDL=FALSE)
    solver_projection <- as.matrix(solve(factor,
      Matrix::crossprod(embedding,mesh$B),system="A"))
    if (max(abs(embedding %*% solver_projection-projection$P)) > 1e-11)
      stop("The continuum projections do not lie in the selected symmetry subspace.")
  }
  operator <- rough_operator(solver_mesh, alpha, tau, order, diffusion,coefficient_source)
  tail <- rough_tail(mesh$J, alpha, diffusion, mesh$kappa, T_final, tail_tolerance)
  general <- function(A) as(as(A, "generalMatrix"), "CsparseMatrix")
  started <- proc.time()[["elapsed"]]
  result <- rough_matrix_evolve_cpp(general(solver_mesh$C),
    lapply(operator$partial_fraction_terms, general), operator$residues,
    operator$remainder, solver_mesh$load, solver_projection, projection$residual_gram,
    rough_rates(mesh$j, alpha, diffusion, mesh$kappa), mesh$j,
    tau, N, tail$squared, return_states)
  result$elapsed <- proc.time()[["elapsed"]] - started
  result$n_solver_nodes <- nrow(solver_mesh$C)
  if (!is.null(embedding)) {
    result$final_state <- as.numeric(embedding %*% result$final_state)
    if (return_states) result$states <- as.matrix(embedding %*% result$states)
  }
  result$tau <- tau
  result$N <- N
  result$order <- order
  result$coefficient_source <- coefficient_source
  result$coefficient_table <- operator$transfer$coefficients$table
  result$tail <- tail
  result$transfer_reconstruction_error <- operator$transfer$reconstruction_error
  result
}
