# Mixed Dirichlet/Kirchhoff heat equation. Run: source("star_dirichlet.R")
# Zero initial data; selected leaves have time-dependent prescribed values.
# The animation is kept in star_dirichlet_animation and displayed in interactive R.
# This script reads brasilfunc.Rmd; it does not regenerate brasilfunc.R.

.star_dir <- local({
  source_file <- sys.frame(1)$ofile
  if (is.null(source_file)) {
    argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    source_file <- if (length(argument)) sub("^--file=", "", argument[1]) else "star_dirichlet.R"
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

# Equation (5.2): y_t - y'' + c0*y = f, y(0) = 0.
# Dirichlet data are prescribed on selected leaves; all other vertices have
# Kirchhoff conditions (at a free leaf this is a zero Neumann condition).
n <- 6L
my_length <- 10
h <- 0.1
T_final <- 2 * pi
time_step <- 0.04
alpha <- 2
c0 <- 0
forcing_amplitude <- 0  # f(s,t) = forcing_amplitude * sin(t); 0 isolates boundary forcing.

# Leaf numbers follow the supplied edge order: angles 2*pi*(j-1)/n.
dirichlet_leaves <- c(1L, 3L)
boundary_amplitudes <- c(1, 0.5)
dirichlet_values <- function(t) boundary_amplitudes * sin(t)
stopifnot(alpha == 2, c0 >= 0, n >= 2,
          length(dirichlet_leaves) > 0L,
          !anyDuplicated(dirichlet_leaves),
          all(dirichlet_leaves %in% seq_len(n)),
          length(boundary_amplitudes) == length(dirichlet_leaves))

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
stopifnot(all(is.finite(boundary_values)), max(abs(boundary_values[, 1L])) < 1e-12)

# The ordinary (beta=1) branch factors C_FF + tau*K_FF once.
operator <- .star_solver$my.fractional.operators.frac(
  L = K_FF, beta = alpha / 2, C = C_FF,
  scale.factor = 1, time_step = tau, diffusion = 1)
spatial_load <- as.numeric(C %*% rep(1, nrow(C)))
g_average <- .star_solver$average_sine_forcing(
  time, A = forcing_amplitude, omega = 1)
solution <- matrix(0, nrow(C), length(time))
for (step in seq_len(grid$N)) {
  g_old <- boundary_values[, step]
  g_new <- boundary_values[, step + 1L]
  # Eliminate the prescribed coefficients BEFORE solving the free equations.
  # The mass correction is essential when the boundary values vary in time.
  rhs <- C_FF %*% solution[F, step] +
    tau * spatial_load[F] * g_average[step + 1L] -
    C_FD %*% (g_new - g_old) - tau * K_FD %*% g_new
  solution[F, step + 1L] <- as.numeric(.star_solver$my.solver.frac(operator, rhs))
  solution[D, step + 1L] <- g_new
}
stopifnot(all(is.finite(solution)))
message(sprintf("Star: %d Dirichlet leaves, %d free unknowns, %d steps (alpha = 2).",
                length(D), length(F), grid$N))

# Simple 3D animation: the star lies at z=0 and solution height is z=u.
# NA separates mesh segments, avoiding lines between different spokes.
columns <- unique(as.integer(round(seq(1, ncol(solution), length.out = 81))))
segment_nodes <- head(as.vector(t(cbind(graph$mesh$E, NA_integer_))), -1L)
x <- graph$mesh$V[segment_nodes, 1L]
y <- graph$mesh$V[segment_nodes, 2L]
z_range <- range(c(0, solution))
z_range <- z_range + c(-0.05, 0.05) * max(diff(z_range), 1e-6)
frame_names <- as.character(seq_along(columns))
star_dirichlet_animation <- plotly::plot_ly() |>
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
    title = "Star graph: zero initial data, prescribed Dirichlet leaves",
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
star_dirichlet_animation <- plotly::plotly_build(star_dirichlet_animation)
star_dirichlet_animation$x$frames <- lapply(seq_along(columns), function(i) list(
  # Plotly trace indices are zero-based; update the solution and boundary markers.
  name = frame_names[i], traces = I(c(1L, 2L)),
  data = list(list(type = "scatter3d", z = solution[segment_nodes, columns[i]]),
              list(type = "scatter3d", z = I(boundary_values[, columns[i]])))))
if (interactive()) print(star_dirichlet_animation)
