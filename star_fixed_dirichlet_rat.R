# Spectral fractional diffusion with mixed Dirichlet/Kirchhoff conditions. Run: source("star_fixed_dirichlet_rat.R")
# Compatible stationary initial profile, f(s,t) = sin(t), and fixed Dirichlet values.
# The animation is kept in star_fixed_dirichlet_rat_animation and displayed in interactive R.
# This script reads brasilfunc.Rmd; it does not regenerate brasilfunc.R.

.star_dir <- local({
  source_file <- sys.frame(1)$ofile
  if (is.null(source_file)) {
    argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    source_file <- if (length(argument)) sub("^--file=", "", argument[1]) else "star_fixed_dirichlet_rat.R"
  }
  dirname(normalizePath(source_file, mustWork = TRUE))
})

library(MetricGraph)
library(Matrix)

# Extract the notebook functions to a temporary file, leaving project files alone.
.star_solver <- local({
  extracted <- tempfile(fileext = ".R")
  on.exit(unlink(extracted))
  knitr::purl(file.path(.star_dir, "brasilfunc.Rmd"), output = extracted,
              documentation = 0L, quiet = TRUE)
  env <- new.env(parent = globalenv())
  sys.source(extracted, envir = env)
  env
})

edges_roots_of_unity <- function(n, length = 10) {
  x_coords <- cos(2 * pi * (0:(n-1)) / n)
  y_coords <- sin(2 * pi * (0:(n-1)) / n)
  return(lapply(1:n, function(i) {
    rbind(c(0, 0), c(length*x_coords[i], length*y_coords[i]))
  }))
}

# Fractional extension of (5.2), defined by harmonic lifting:
# y = w + H*g,  w|D = 0,  y_t + A_D^(alpha/2)*w = f.
# A_D is -Delta_Gamma+c0 with homogeneous Dirichlet conditions at D and
# Kirchhoff conditions elsewhere. H*g is its local harmonic extension:
# (H*g)|D=g, K_FF*(H*g)_F + K_FD*g=0.
# Hence C_FF*w_t + C_FF*A_D^(alpha/2)*w = b_F - (C*H)_F*g_t.
# This specifies the nonhomogeneous spectral model; it is not obtained by
# overwriting entries of an all-Kirchhoff fractional solution. At alpha=2
# it reduces to the original mixed-boundary heat equation.
# Background on harmonic lifting for nonhomogeneous spectral operators:
# https://arxiv.org/abs/1703.05256 (adapted here to the metric-graph FEM).
n <- 6L
my_length <- 10
h <- 0.1
T_final <- 2 * pi
time_step <- 0.04
alpha <- 1
m <- 12L  # BRASIL rational approximation order; beta = alpha/2.
c0 <- 0
forcing_amplitude <- 1  # f(s,t) = forcing_amplitude * sin(t).

# Leaf numbers follow the supplied edge order: angles 2*pi*(j-1)/n.
dirichlet_leaves <- c(1L, 3L)
fixed_dirichlet_values <- c(1, 0.5)
dirichlet_values <- function(t) fixed_dirichlet_values
stopifnot(alpha > 0, alpha <= 2, c0 >= 0, n >= 2,
          length(dirichlet_leaves) > 0L,
          !anyDuplicated(dirichlet_leaves),
          all(dirichlet_leaves %in% seq_len(n)),
          is.finite(forcing_amplitude),
          length(fixed_dirichlet_values) == length(dirichlet_leaves),
          all(is.finite(fixed_dirichlet_values)))

edges <- edges_roots_of_unity(n, length = my_length)
graph <- metric_graph$new(edges = edges)
graph$build_mesh(h = h)
graph$compute_fem()
C <- graph$mesh$C
K <- graph$mesh$G + c0 * C

# Identify prescribed mesh nodes geometrically rather than assuming a node order.
D <- vapply(dirichlet_leaves, function(j) {
  endpoint <- edges[[j]][2L, ]
  distance <- rowSums(sweep(graph$mesh$V[, 1:2, drop = FALSE], 2, endpoint, "-")^2)
  ids <- which(distance < 1e-16)
  stopifnot(length(ids) == 1L)
  as.integer(ids)
}, integer(1))
F <- setdiff(seq_len(nrow(C)), D)
C_FF <- C[F, F, drop = FALSE]
C_FD <- C[F, D, drop = FALSE]
K_FF <- K[F, F, drop = FALSE]
K_FD <- K[F, D, drop = FALSE]

grid <- .star_solver$make_time_grid(T_final, time_step)
time <- grid$time_seq
tau <- grid$time_step
boundary_values <- vapply(time, dirichlet_values, numeric(length(D)))
boundary_values <- matrix(boundary_values, nrow = length(D))
stopifnot(all(is.finite(boundary_values)))
# Build one lifting basis vector per prescribed leaf.
H <- matrix(0, nrow(C), length(D))
H[D, ] <- diag(length(D))
H[F, ] <- as.matrix(Matrix::solve(K_FF, -K_FD))
lifting <- H %*% boundary_values
lifting_load <- (C %*% H)[F, , drop = FALSE]

# With at least one Dirichlet vertex on a connected graph of total length L,
# ||v||^2 <= L^2*||v'||^2 for v|D=0. Thus this positive scale is <= lambda_min
# even when c0=0, keeping the normalized spectrum in BRASIL's [1, infinity).
# No artificial reaction shift is added to the PDE.
operator_scale <- c0 + 1 / sum(graph$get_edge_lengths())^2
operator <- .star_solver$my.fractional.operators.frac(
  L = K_FF, beta = alpha / 2, C = C_FF, m = m,
  scale.factor = operator_scale, time_step = tau, diffusion = 1)
spatial_load <- as.numeric(C %*% rep(1, nrow(C)))
g_average <- .star_solver$average_sine_forcing(
  time, A = forcing_amplitude, omega = 1)

# All three experiments start with w=0: either y(0)=0 and g(0)=0,
# or the fixed-boundary experiment starts from y(0)=H*g.
remainder <- matrix(0, length(F), length(time))
solution <- lifting
for (step in seq_len(grid$N)) {
  delta_g <- boundary_values[, step + 1L] - boundary_values[, step]
  rhs <- C_FF %*% remainder[, step] +
    tau * spatial_load[F] * g_average[step + 1L] - lifting_load %*% delta_g
  remainder[, step + 1L] <- as.numeric(.star_solver$my.solver.frac(operator, rhs))
  solution[F, step + 1L] <- remainder[, step + 1L] + lifting[F, step + 1L]
}
stopifnot(all(is.finite(solution)))
message(sprintf("Fractional star: %d Dirichlet leaves, %d free unknowns, %d steps (alpha = %g, m = %d).",
                length(D), length(F), grid$N, alpha, m))

# Simple 3D animation: the star lies at z=0 and solution height is z=u.
# NA separates mesh segments, avoiding lines between different spokes.
columns <- unique(as.integer(round(seq(1, ncol(solution), length.out = 81))))
segment_nodes <- head(as.vector(t(cbind(graph$mesh$E, NA_integer_))), -1L)
x <- graph$mesh$V[segment_nodes, 1L]
y <- graph$mesh$V[segment_nodes, 2L]
z_range <- range(c(0, solution))
z_range <- z_range + c(-0.05, 0.05) * max(diff(z_range), 1e-6)
frame_names <- as.character(seq_along(columns))
star_fixed_dirichlet_rat_animation <- plotly::plot_ly() |>
  plotly::add_trace(x = x, y = y, z = rep(0, length(x)),
                    type = "scatter3d", mode = "lines", name = "Star graph",
                    line = list(color = "grey65", width = 3), hoverinfo = "skip") |>
  plotly::add_trace(x = x, y = y, z = solution[segment_nodes, 1L],
                    type = "scatter3d", mode = "lines", name = "u(s,t)",
                    line = list(color = "#2166ac", width = 6), connectgaps = FALSE) |>
  plotly::add_trace(x = graph$mesh$V[D, 1L], y = graph$mesh$V[D, 2L],
                    z = boundary_values[, 1L], type = "scatter3d", mode = "markers",
                    name = "Dirichlet leaves", marker = list(color = "#d73027", size = 5)) |>
  plotly::layout(
    title = "Fractional star (alpha = 1): f(s,t) = sin(t), fixed Dirichlet leaves",
    scene = list(xaxis = list(title = "x"), yaxis = list(title = "y"),
                 zaxis = list(title = "u", range = z_range, autorange = FALSE),
                 aspectmode = "manual", aspectratio = list(x = 1, y = 1, z = 0.6)),
    updatemenus = list(list(type = "buttons", showactive = FALSE,
      buttons = list(
        list(label = "Play", method = "animate", args = list(NULL,
          list(fromcurrent = TRUE, frame = list(duration = 80, redraw = TRUE),
               transition = list(duration = 0)))),
        list(label = "Pause", method = "animate", args = list(list(NULL),
          list(mode = "immediate", frame = list(duration = 0, redraw = FALSE))))))),
    sliders = list(list(currentvalue = list(prefix = "t = "),
      steps = lapply(seq_along(columns), function(i) list(
        label = sprintf("%.2f", time[columns[i]]), method = "animate",
        args = list(list(frame_names[i]), list(mode = "immediate",
          frame = list(duration = 0, redraw = TRUE), transition = list(duration = 0))))))))
star_fixed_dirichlet_rat_animation <- plotly::plotly_build(star_fixed_dirichlet_rat_animation)
star_fixed_dirichlet_rat_animation$x$frames <- lapply(seq_along(columns), function(i) list(
  # Plotly trace indices are zero-based; update the solution and boundary markers.
  name = frame_names[i], traces = I(c(1L, 2L)),
  data = list(list(type = "scatter3d", z = solution[segment_nodes, columns[i]]),
              list(type = "scatter3d", z = I(boundary_values[, columns[i]])))))
if (interactive()) print(star_fixed_dirichlet_rat_animation)
