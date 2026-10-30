# Fixed L2 forcing on the length-(1,2) tadpole. All j >= J loads vanish
# exactly on the uniform mesh h = 2^-J; they remain in the continuum error.
rough_mesh <- function(J, kappa = 4) {
  stopifnot(length(J) == 1L, is.finite(J), J == as.integer(J), J >= 2L,
            J <= 20L, length(kappa) == 1L, is.finite(kappa), kappa > 0)
  theta <- seq(-pi, pi, length.out = 129L)
  edges <- list(rbind(c(0, 0), c(1, 0)),
                cbind(1 + 1/pi + cos(theta)/pi, sin(theta)/pi))
  constructor <- MetricGraph::metric_graph$public_methods$initialize
  if ("edges" %in% names(formals(constructor))) {
    graph <- MetricGraph::metric_graph$new(edges = edges, verbose = 0)
    graph$set_manual_edge_lengths(edge_lengths = c(1, 2))
  } else {
    lines <- sp::SpatialLines(lapply(seq_along(edges), function(i)
      sp::Lines(list(sp::Line(edges[[i]])), ID = as.character(i))))
    graph <- MetricGraph::metric_graph$new(lines = lines, remove_circles = FALSE,
                                            verbose = FALSE)
    graph$edge_lengths <- c(1, 2)
  }
  h <- 2^-J
  graph$build_mesh(h = h)
  graph$compute_fem()
  C <- graph$mesh$C
  K <- kappa^2 * C + graph$mesh$G
  PtE <- graph$mesh$PtE
  stopifnot(nrow(C) == 3L * 2^J, max(abs(graph$mesh$h_e - h)) < 1e-13)
  loop <- which(PtE[, 1L] == 2L)
  loop <- loop[order(PtE[loop, 2L])]
  loop_nodes <- c(2L, 2L + loop, 2L)
  x <- seq(0, 2, by = h)
  stopifnot(length(loop_nodes) == length(x))
  j <- seq_len(J - 1L)
  theta_j <- pi * 2^j * h
  # Analytic integral of sin(w*x) against each assembled periodic hat.
  B <- matrix(0, nrow(C), length(j))
  B[2L + loop, ] <- outer(2 * PtE[loop, 2L], 2^j, function(x, k) sinpi(x*k)) *
    rep(h * (sin(theta_j/2)/(theta_j/2))^2, each = length(loop))
  load <- as.numeric(B %*% (1/j))
  list(J = J, h = h, C = C, K = K, graph = graph, B = B,
       load = load, j = j, loop_nodes = loop_nodes, loop_x = x,
       kappa = kappa)
}

rough_projection <- function(mesh, quadrature_order = 16L) {
  stopifnot(quadrature_order >= 8L, quadrature_order == as.integer(quadrature_order))
  Cfactor <- Matrix::Cholesky(Matrix::forceSymmetric(mesh$C), LDL = FALSE)
  P <- as.matrix(solve(Cfactor, mesh$B, system = "A"))
  # Integrate the projection residual itself, rather than subtracting two
  # almost equal Gram matrices. Resolved modes have at least four cells/cycle.
  rule <- statmod::gauss.quad(quadrature_order, kind = "legendre")
  z <- (rule$nodes + 1)/2
  weights <- rule$weights/2
  nodes <- mesh$loop_nodes
  left <- head(nodes, -1L)
  right <- tail(nodes, -1L)
  R <- matrix(0, length(mesh$j), length(mesh$j))
  for (q in seq_along(z)) {
    x <- head(mesh$loop_x, -1L) + mesh$h*z[q]
    exact <- outer(x, 2^mesh$j, function(x, k) sinpi(x*k))
    residual <- exact - ((1-z[q])*P[left, , drop = FALSE] + z[q]*P[right, , drop = FALSE])
    R <- R + mesh$h * weights[q] * crossprod(residual)
  }
  # Include any roundoff leakage of the mass projection onto the stem.
  stem <- c(1L, 2L + which(mesh$graph$mesh$PtE[, 1L] == 1L), 2L)
  for (q in seq_along(z)) {
    residual <- (1-z[q])*P[head(stem, -1L), , drop = FALSE] +
      z[q]*P[tail(stem, -1L), , drop = FALSE]
    R <- R + mesh$h * weights[q] * crossprod(residual)
  }
  list(P = P, residual_gram = R)
}

rough_rates <- function(j, alpha, diffusion = 1, kappa = 4) {
  stopifnot(all(is.finite(j)), all(j >= 1), all(j == floor(j)),
            length(alpha) == 1L, is.finite(alpha), alpha > 0, alpha < 2,
            length(diffusion) == 1L, is.finite(diffusion), diffusion > 0,
            length(kappa) == 1L, is.finite(kappa), kappa > 0)
  diffusion * (kappa^2 + pi^2 * 2^(2*j))^(alpha/2)
}

# Integral_0^T (1-exp(-a*t))^2 dt, stable also for small a*T.
rough_response_integral <- function(a, T_final = 1) {
  stopifnot(all(is.finite(a)), all(a > 0), T_final > 0, is.finite(T_final))
  z <- a * T_final
  answer <- T_final * (1 + 2*expm1(-z)/z - expm1(-2*z)/(2*z))
  small <- z < .1
  if (any(small)) {
    zz <- z[small]
    value <- rep(0, length(zz))
    for (n in 2:18) value <- value + (-1)^n * (2^n-2) * zz^n/factorial(n+1)
    answer[small] <- T_final * value
  }
  answer
}

rough_tail <- function(J, alpha, diffusion = 1, kappa = 4, T_final = 1,
                       relative_tolerance = 1e-13) {
  stopifnot(J >= 1L, J == as.integer(J), relative_tolerance > 0,
            relative_tolerance < 1)
  total <- 0
  for (j in J:450L) {
    a <- rough_rates(j, alpha, diffusion, kappa)
    term <- rough_response_integral(a, T_final)/(j*a)^2
    total <- total + term
    # Every omitted term is at most T/(nu^2*pi^(2*alpha))
    # * 2^(-2*alpha*j)/j^2. Bound the infinite remainder geometrically.
    next_j <- j + 1L
    bound <- T_final/(diffusion^2*pi^(2*alpha)) *
      2^(-2*alpha*next_j)/next_j^2/(1-2^(-2*alpha))
    if (bound <= relative_tolerance * total)
      return(list(squared = total, omitted_squared_bound = bound, last_j = j))
  }
  stop("Infinite forcing tail did not meet its requested tolerance.")
}
