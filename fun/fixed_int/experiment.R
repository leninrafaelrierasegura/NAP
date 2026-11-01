# New helpers stay local to this study. Shared mesh/load/reference definitions
# are read into a private environment; no original coefficient functions change.
source(here::here("fun", "fixed_int", "coefficients.R"), local = TRUE)
.fixed_int_base <- new.env(parent = environment())
sys.source(here::here("fun", "isolated_rational_convergence.R"), .fixed_int_base)
.fixed_int_base$solve <- Matrix::solve
.fixed_int_base$isSymmetric <- Matrix::isSymmetric
.fixed_int_base$crossprod <- Matrix::crossprod

fixed_int_discrete_spectrum <- function(...) .fixed_int_base$isolated_discrete_spectrum(...)
fixed_int_reference <- function(...) .fixed_int_base$isolated_exact_be_reference(...)

fixed_int_tadpole_data <- function(h, kappa = 4) {
  stopifnot(length(h) == 1L, is.finite(h), h > 0, length(kappa) == 1L,
            is.finite(kappa), kappa > 0)
  if (3/h > 1200) stop("Choose h >= 0.0025 for this complete-eigensystem study.")
  if ("edges" %in% names(formals(MetricGraph::metric_graph$public_methods$initialize))) {
    result <- .fixed_int_base$isolated_tadpole_data(h, kappa)
    result$mesh_backend <- "MetricGraph graph/mesh API"
    return(result)
  }
  # MetricGraph 1.1.2 predates `edges` and set_manual_edge_lengths(). Use its
  # own P1 assembler on exactly the prescribed tadpole topology/edge lengths.
  # A thin graph interface lets the original analytic load code run unchanged.
  compatibility <- new.env(parent = .fixed_int_base)
  compatibility$gets.graph.tadpole <- function(h) {
    counts <- as.integer(ceiling(c(1,2)/h))
    n <- as.integer(2+sum(counts-1L))
    tail_nodes <- if (counts[1L] > 1L) seq.int(3L, counts[1L]+1L) else integer()
    loop_nodes <- if (counts[2L] > 1L) seq.int(counts[1L]+2L,n) else integer()
    paths <- list(c(1L,tail_nodes,2L), c(2L,loop_nodes,2L))
    E <- do.call(rbind,lapply(paths,function(nodes) cbind(head(nodes,-1L),tail(nodes,-1L))))
    lengths <- rep(c(1,2)/counts,counts)
    fem <- get("assemble_fem",asNamespace("MetricGraph"))(E,lengths,n)
    PtE <- do.call(rbind,lapply(1:2,function(edge)
      cbind(edge=rep(edge,counts[edge]-1L), position=seq_len(counts[edge]-1L)/counts[edge])))
    list(mesh=c(fem,list(PtE=PtE)), compute_fem=function() invisible(NULL),
         get_edge_lengths=function() c(1,2))
  }
  build <- .fixed_int_base$isolated_tadpole_data
  environment(build) <- compatibility
  result <- build(h,kappa)
  result$mesh_backend <- "MetricGraph P1 assembler (compatibility with version 1.1.2)"
  result
}

fixed_int_operator <- function(L, C, beta, kappa, tau, diffusion = 1, coefficients = NULL) {
  stopifnot(length(beta)==1L,is.finite(beta),beta>0,beta<=1,
    length(kappa)==1L,is.finite(kappa),kappa>0,length(tau)==1L,is.finite(tau),tau>0,
    length(diffusion)==1L,is.finite(diffusion),diffusion>0,
    identical(dim(L),dim(C)),Matrix::isSymmetric(L,tol=1e-12),Matrix::isSymmetric(C,tol=1e-12))
  if (beta==1) return(list(C=C,beta=beta,
    factorization=Matrix::Cholesky(Matrix::forceSymmetric(C+tau*diffusion*L),LDL=FALSE)))
  fixed_int_validate_coefficients(coefficients)
  if (abs(coefficients$beta-beta)>32*.Machine$double.eps)
    stop("The coefficient exponent does not match the physical exponent.")
  transfer <- .fixed_int_base$brasil_partial_fractions(coefficients,tau,diffusion*kappa^(2*beta))
  shifted <- lapply(transfer$p,function(q) Matrix::forceSymmetric(L/kappa^2-q*C))
  list(C=C,beta=beta,m=coefficients$m,residues=transfer$r,poles=transfer$p,
    remainder=transfer$k,partial_fraction_terms=shifted,
    mass_factor=if(transfer$k!=0) Matrix::Cholesky(Matrix::forceSymmetric(C),LDL=FALSE) else NULL,
    factorizations=lapply(shifted,Matrix::Cholesky,LDL=FALSE),
    reconstruction_error=transfer$reconstruction_error)
}

fixed_int_solve <- function(operator, rhs) .fixed_int_base$my.solver.frac(operator,rhs)
fixed_int_inverse <- function(y, coefficients) .fixed_int_base$brasil_inverse(y,coefficients)
fixed_int_power <- function(x, coefficients) .fixed_int_base$brasil_power(x,coefficients)

fixed_int_error <- function(data, spectrum, reference, m, bank = NULL, coefficients = NULL,
                            source = "Fixed interval", excitation_tolerance = 1e-8) {
  stopifnot(length(m)==1L,is.finite(m),m>=1,m==floor(m),
    nrow(reference$states)==nrow(data$C),ncol(reference$states)==length(reference$time_seq))
  beta <- reference$alpha/2
  if (beta!=1 && is.null(coefficients)) coefficients <- fixed_int_coefficients(bank,m,beta)
  if (beta!=1 && coefficients$m!=m) stop("The coefficient order does not match m.")
  normalized <- spectrum$values/data$kappa^2
  if (beta!=1 && !is.null(coefficients$x_min) &&
      max(normalized)>1/coefficients$x_min*(1+64*.Machine$double.eps))
    stop("The coefficients do not cover the entire discrete spectrum; regenerate the bank.")
  operator <- fixed_int_operator(data$L,data$C,beta,data$kappa,reference$dt,
                                 reference$diffusion,coefficients)
  if (beta==1) {
    rates <- reference$diffusion*spectrum$values
  } else {
    S <- fixed_int_inverse(normalized,coefficients)
    if (any(!is.finite(S)) || any(S<=0)) stop("Invalid inverse-power approximation on the spectrum.")
    rates <- reference$diffusion*data$kappa^(2*beta)/S
  }
  mu <- 1/(1+reference$dt*rates)
  u <- as.numeric(data$initial)
  y <- reference$initial_coefficients
  total <- spectral_total <- compensation <- spectral_compensation <- 0
  for (step in seq_len(length(reference$time_seq)-1L)) {
    forcing <- reference$g_average[step+1L]
    rhs <- data$C%*%u+reference$dt*data$spatial_load*forcing
    u <- as.numeric(fixed_int_solve(operator,rhs))
    delta <- reference$states[,step+1L]-u
    value <- reference$dt*sum(delta*(data$C%*%delta))
    adjusted <- value-compensation
    updated <- total+adjusted
    compensation <- (updated-total)-adjusted
    total <- updated
    # Independent recurrence evaluates S directly and never uses the poles or
    # residues of H. All numerical eigenvectors, including weak modes, remain.
    y <- mu*(y+reference$dt*reference$forcing_coefficients*forcing)
    dy <- reference$coefficients[,step+1L]-y
    value <- reference$dt*sum(dy^2)
    adjusted <- value-spectral_compensation
    updated <- spectral_total+adjusted
    spectral_compensation <- (updated-spectral_total)-adjusted
    spectral_total <- updated
  }
  sparse_error <- sqrt(max(0,total))
  product_error <- sqrt(max(0,spectral_total))
  gap <- u-as.numeric(spectrum$vectors%*%y)
  relative_rate_error <- rates/reference$exact_rates-1
  excited <- abs(reference$initial_coefficients)>excitation_tolerance*max(abs(reference$initial_coefficients)) |
    abs(reference$forcing_coefficients)>excitation_tolerance*max(abs(reference$forcing_coefficients))
  data.frame(alpha=reference$alpha,m=m,source=source,error=sparse_error,product_error=product_error,
    absolute_error_disagreement=abs(sparse_error-product_error),
    final_state_gap_C=sqrt(max(0,sum(gap*(data$C%*%gap)))),
    max_relative_rate_error_full=max(abs(relative_rate_error)),
    max_relative_rate_error_excited=if(any(excited)) max(abs(relative_rate_error[excited])) else 0,
    reconstruction_error=if(beta==1) 0 else operator$reconstruction_error)
}

fixed_int_table_coefficients <- function(m,beta) .fixed_int_base$brasil_coefficients(m,beta)

fixed_int_mark_resolution <- function(results, control_error) {
  results$noise_indicator <- pmax(results$absolute_error_disagreement,control_error,.Machine$double.xmin)
  results$resolved <- results$product_error>100*results$noise_indicator
  results
}

fixed_int_slopes <- function(results) {
  fit <- function(d, predictor) {
    if(nrow(d)<2L) return(NA_real_)
    unname(coef(lm(log(error)~predictor,d))[2L])
  }
  do.call(rbind,lapply(split(results,results$alpha),function(d) {
    alpha <- d$alpha[1L]
    d <- d[d$resolved & is.finite(d$error) & d$error>0,]
    data.frame(alpha=alpha,resolved_orders=nrow(d),
      first_m=if(nrow(d)) min(d$m) else NA_integer_,last_m=if(nrow(d)) max(d$m) else NA_integer_,
      slope_log_error_vs_sqrt_m=fit(d,sqrt(d$m)),slope_log_error_vs_m=fit(d,d$m))
  }))
}
