experiment_grid <- function(T_final, target_step) {
  stopifnot(length(T_final) == 1L, is.finite(T_final), T_final > 0,
            length(target_step) == 1L, is.finite(target_step), target_step > 0)
  ratio <- T_final/target_step
  # Snap only a quotient already very close to an integer. Subtracting a
  # relative tolerance can remove whole steps on grids with ~1e15 intervals.
  nearest <- round(ratio)
  if (abs(ratio-nearest) <= min(.01, 32*.Machine$double.eps*max(1, ratio)))
    ratio <- nearest
  N <- ceiling(ratio)
  if (N > 2^52) stop("Time grid exceeds the supported exact step-index range.")
  list(N = max(1, N), time_step = T_final/max(1, N), T_final = T_final)
}

uniform_tadpole_h <- function(h) {
  stopifnot(all(is.finite(h)), all(h > 0))
  inverse <- 1/h
  1/pmax(1, ceiling(inverse-32*.Machine$double.eps*pmax(1, inverse)))
}

uniform_tadpole_basis <- function(h, kappa, modes, N_finite) {
  stopifnot(length(kappa) == 1L, is.finite(kappa), kappa > 0,
            length(N_finite) == 1L, N_finite >= 0,
            N_finite == as.integer(N_finite), length(modes) > 0,
            all(modes == as.integer(modes)), !anyDuplicated(modes))
  catalogue <- unlist(lapply(0:N_finite, function(j)
    rep(j, if (j > 0 && j %% 2L == 0L) 2L else 1L)))
  stopifnot(all(modes >= 1L), all(modes <= length(catalogue)))
  h <- uniform_tadpole_h(h)
  w <- catalogue[modes]*pi/2
  z <- w*h
  if (any(z >= pi/2)) stop("Uniform modal backend needs at least four elements per wavelength.")
  residual <- numeric(length(z))
  small <- abs(z) < .1
  x <- z[small]^2
  residual[small] <- x^2*(1/720+x*(1/3024+x*(7/259200+x*(
    8.350702795147226e-7+x*(-9.986163466652883e-8+x*(
    -1.826716236438457e-8+x*(-1.465654963502789e-9)))))))
  zz <- z[!small]
  residual[!small] <- 1-3*(sin(zz/2)/(zz/2))^4/(2+cos(zz))
  list(h = h, eigenvalues = kappa^2+w^2,
       discrete_eigenvalues = kappa^2+12*sin(z/2)^2/(h^2*(2+cos(z))),
       gram = diag(1-residual, length(z)),
       residual_gram = diag(residual, length(z)))
}

discrete_uniform_rule <- function(N, order = 64L) {
  stopifnot(length(N) == 1L, is.finite(N), N >= 1, N == round(N),
            N <= 2^52, length(order) == 1L, is.finite(order), order >= 2L,
            order == as.integer(order), order <= 512L)
  order <- min(N, order)
  if (order == 1L) return(list(nodes = 1, weights = 1))
  j <- seq_len(order-1L)
  off <- .5*j*sqrt((1-(j/N)^2)/(4*j^2-1))
  J <- diag(.5, order)
  J[cbind(j, j+1L)] <- J[cbind(j+1L, j)] <- off
  e <- eigen(J, symmetric = TRUE)
  list(nodes = pmin(N, pmax(1, N*e$values+.5)), weights = e$vectors[1, ]^2)
}

rational_modal_rates <- function(eigenvalues, alpha, m, diffusion, kappa,
                                 root_interpolation = "linear", coefficient_source = "computed") {
  beta <- alpha/2
  stopifnot(length(alpha) == 1L, is.finite(alpha), alpha > 0, alpha <= 2,
            length(diffusion) == 1L, is.finite(diffusion), diffusion > 0,
            all(is.finite(eigenvalues)), all(eigenvalues >= kappa^2))
  # m = 0 is reserved for error-floor diagnostics, never a study order.
  if (beta == 1 || m == 0) return(diffusion*eigenvalues^beta)
  roots <- my.get.roots(m, beta, interpolation = root_interpolation,
                        coefficient_source = coefficient_source)
  z <- eigenvalues/kappa^2
  inverse <- roots$factor*vapply(z, function(x)
    prod(1-roots$pr_roots*x)/prod(1-roots$pl_roots*x), numeric(1))
  if (any(!is.finite(inverse)) || any(inverse <= 0)) stop("Invalid rational modal rates.")
  diffusion*kappa^(2*beta)/inverse
}

uniform_modal_coefficients <- function(k, time_step, rates, initial, forcing,
                                       A = 1, omega = pi) {
  stopifnot(all(is.finite(k)), all(k >= 0), time_step > 0,
            length(rates) == length(initial), length(initial) == length(forcing),
            all(is.finite(c(rates, initial, forcing, A, omega))), all(rates > 0))
  logq <- -log1p(time_step*rates)
  theta <- omega*time_step
  sinc <- if (theta == 0) 1 else sin(theta/2)/(theta/2)
  result <- exp(outer(logq, k))*initial
  # Stable exp(i*x)-1: R does not implement expm1 for complex numbers.
  expm1_i <- function(x) -2*sin(x/2)^2+1i*sin(x)
  for (i in seq_along(rates)) {
    response <- exp(logq[i])*time_step*A*sinc*Im(exp(1i*theta/2)*
      (expm1_i(theta*k)-expm1(logq[i]*k))/(expm1_i(theta)-expm1(logq[i])))
    result[i, ] <- result[i, ]+forcing[i]*response
  }
  result
}

uniform_modal_difference <- function(k, fraction, time_step, exact_rates, numerical_rates,
                                      initial, forcing, A = 1, omega = pi) {
  stopifnot(length(fraction) == 1L, is.finite(fraction), fraction >= 0, fraction <= 1,
            all(is.finite(k)), all(k >= 1), length(time_step) == 1L,
            is.finite(time_step), time_step > 0,
            length(exact_rates) == length(numerical_rates),
            length(initial) == length(exact_rates), length(forcing) == length(initial),
            all(is.finite(c(exact_rates,numerical_rates,initial,forcing,A,omega))),
            all(exact_rates > 0), all(numerical_rates > 0))
  dt <- time_step
  offset <- (1-fraction)*dt
  t <- k*dt-offset
  theta <- omega*dt
  # The cancellation-free expression below is needed for the high-order,
  # tiny-step regime. Coarse intervals have a well-resolved direct difference.
  if (abs(theta) >= .1) return(
    exact_modal_coefficients(t, exact_rates, initial, forcing, A, omega)-
      uniform_modal_coefficients(k, dt, numerical_rates, initial, forcing, A, omega))
  expm1_i <- function(x) -2*sin(x/2)^2+1i*sin(x)
  x <- dt*numerical_rates
  remainder <- x-log1p(x)
  small <- x < .001
  if (any(small)) {
    xx <- x[small]
    # x-log(1+x), retaining the O(x^2) term when direct subtraction rounds to 0.
    value <- numeric(length(xx))
    for (j in 2:12) value <- value+(-1)^j*xx^j/j
    remainder[small] <- value
  }
  logq <- -log1p(x)
  dr <- exact_rates-numerical_rates
  u <- theta/2
  u2 <- u*u
  # exp(i*u)/sinc(u)-1 = u*cot(u)-1 + i*u.
  phi_minus_one <- -u2*(1/3+u2*(1/45+u2*(2/945+u2*(1/4725+u2*2/93555))))+1i*u
  out <- matrix(0, length(initial), length(k))
  for (i in seq_along(initial)) {
    exponent_difference <- k*(dt*dr[i]+remainder[i])-exact_rates[i]*offset
    # Factor out the larger exponential so large decay never gives 0*Inf.
    positive <- exponent_difference > 0
    initial_difference <- numeric(length(k))
    initial_difference[!positive] <- -exp(-exact_rates[i]*t[!positive])*
      expm1(exponent_difference[!positive])
    initial_difference[positive] <- exp(logq[i]*k[positive])*
      expm1(-exponent_difference[positive])
    out[i, ] <- initial[i]*initial_difference
    if (omega != 0 && A != 0 && forcing[i] != 0) {
      den_exact <- exact_rates[i]+1i*omega
      den_numerical <- numerical_rates[i]*(1+phi_minus_one)+1i*omega
      inverse_difference <- (-dr[i]+numerical_rates[i]*phi_minus_one)/
        (den_exact*den_numerical)
      numerator <- expm1_i(theta*k)-expm1(logq[i]*k)
      numerator_difference <- -exp(1i*omega*t)*expm1_i(omega*offset)-initial_difference
      out[i, ] <- out[i, ]+forcing[i]*A*Im(
        inverse_difference*numerator+numerator_difference/den_exact)
    }
  }
  out
}

uniform_tadpole_error <- function(h, alpha, m, diffusion, grid, kappa, modes,
                                  N_finite, initial, forcing, A = 1, omega = pi,
                                  quadrature_order = 4L, step_sum_order = 64L,
                                  root_interpolation = "linear", coefficient_source = "computed") {
  basis <- uniform_tadpole_basis(h, kappa, modes, N_finite)
  N <- grid$N
  T <- grid$T_final
  stopifnot(abs(N*grid$time_step-T) <= 100*.Machine$double.eps*max(1, T))
  rule <- discrete_uniform_rule(N, step_sum_order)
  time_rule <- gauss_legendre_rule(quadrature_order)
  rates <- rational_modal_rates(basis$discrete_eigenvalues, alpha, m,
                                diffusion, kappa, root_interpolation, coefficient_source)
  if (N > step_sum_order &&
      (max(c(rates, diffusion*basis$eigenvalues^(alpha/2)))*T > step_sum_order^2/8 ||
       abs(omega)*T > step_sum_order/2))
    stop("Step-sum quadrature may underresolve the dynamics: increase step_sum_order or use direct stepping.")
  total <- 0
  for (j in seq_along(time_rule$nodes)) {
    t <- (rule$nodes-1+time_rule$nodes[j])*grid$time_step
    exact <- exact_modal_coefficients(t, diffusion*basis$eigenvalues^(alpha/2),
                                      initial, forcing, A, omega)
    delta <- uniform_modal_difference(
      rule$nodes, time_rule$nodes[j], grid$time_step,
      diffusion*basis$eigenvalues^(alpha/2), rates, initial, forcing, A, omega)
    squared <- colSums(delta*(basis$gram %*% delta))+
      colSums(exact*(basis$residual_gram %*% exact))
    total <- total+time_rule$weights[j]*sum(rule$weights*squared)
  }
  if (!is.finite(total) || total < 0) stop("Invalid uniform modal error.")
  sqrt(T*total)
}

tadpole_convergence_error <- function(h, time_step, m, alpha, parameters,
                                     quadrature_order = parameters$quadrature_order,
                                     step_sum_order = parameters$step_sum_order) {
  uniform_tadpole_error(
    h = h, alpha = alpha, m = m, diffusion = parameters$diffusion,
    grid = experiment_grid(parameters$T_final, time_step),
    kappa = parameters$kappa, modes = parameters$modes,
    N_finite = parameters$N_finite, initial = parameters$initial,
    forcing = parameters$forcing, A = parameters$A, omega = parameters$omega,
    quadrature_order = quadrature_order, step_sum_order = step_sum_order,
    coefficient_source = "computed"
  )
}
