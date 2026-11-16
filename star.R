# Star-graph experiment. Run: source("star.R")
# u_t + diffusion * (kappa^2 - Delta_Gamma)^(alpha/2) u = sin(t),
# u(s, 0) = amplitude * exp(-d_Gamma(s, center)^2 / (2 * sigma^2)).
# The animation is kept in star_animation and displayed in interactive R.
# This script reads brasilfunc.Rmd; it does not regenerate brasilfunc.R.

.star_dir <- local({
  source_file <- sys.frame(1)$ofile
  if (is.null(source_file)) {
    argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    source_file <- if (length(argument)) sub("^--file=", "", argument[1]) else "star.R"
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

# Editable experiment parameters. alpha = 2 selects ordinary backward Euler;
# 0 < alpha < 2 selects the BRASIL fractional matrix solver (kappa > 0).
n <- 6L
my_length <- 10
h <- 0.1
T_final <- 2 * pi
time_step <- 0.04
alpha <- 1
kappa <- 1
diffusion <- 1
m <- 6L
sigma <- 0.7
amplitude <- 1

edges <- edges_roots_of_unity(n, length = my_length)
graph <- metric_graph$new(edges = edges)
center_vertex <- which(rowSums(graph$V[, 1:2, drop = FALSE]^2) < 1e-20)
stopifnot(length(center_vertex) == 1L, nrow(graph$V) == n + 1L,
          sigma > 0, alpha > 0, alpha <= 2, kappa >= 0,
          alpha == 2 || kappa > 0)
graph$build_mesh(h = h)
graph$compute_fem()
C <- graph$mesh$C
L <- kappa^2 * C + graph$mesh$G

# L2-project the Gaussian using a finer quadrature mesh. On these straight
# spokes, distance from (0,0) is precisely intrinsic distance from the hub.
quadrature_graph <- graph$clone(deep = TRUE)
quadrature_graph$build_mesh(h = h / 4)
quadrature_graph$compute_fem()
radius <- sqrt(rowSums(quadrature_graph$mesh$V[, 1:2, drop = FALSE]^2))
gaussian <- amplitude * exp(-radius^2 / (2 * sigma^2))
Psi <- graph$fem_basis(quadrature_graph$get_mesh_locations())
initial_load <- Matrix::crossprod(Psi, quadrature_graph$mesh$C %*% gaussian)
u0 <- as.numeric(Matrix::solve(C, initial_load))

grid <- .star_solver$make_time_grid(T_final, time_step)
time <- grid$time_seq
tau <- grid$time_step
operator <- .star_solver$my.fractional.operators.frac(
  L = L, beta = alpha / 2, C = C,
  scale.factor = if (alpha == 2) 1 else kappa^2,
  m = m, time_step = tau, diffusion = diffusion)

# f(s,t) = sin(t) on every edge. Its spatial FEM load is C * 1.
# Use the exact time average of sin(t) over each backward Euler step.
spatial_load <- as.numeric(C %*% rep(1, nrow(C)))
g_average <- .star_solver$average_sine_forcing(time, A = 1, omega = 1)
solution <- matrix(0, nrow(C), length(time))
solution[, 1L] <- u0
for (step in seq_len(grid$N)) {
  rhs <- C %*% solution[, step] + tau * spatial_load * g_average[step + 1L]
  solution[, step + 1L] <- as.numeric(.star_solver$my.solver.frac(operator, rhs))
}
stopifnot(all(is.finite(solution)))
message(sprintf("Star: %d spokes, %d unknowns, %d time steps (alpha = %g).",
                n, nrow(C), grid$N, alpha))

# Simple 3D animation: the star lies at z=0 and solution height is z=u.
# NA separates mesh segments, avoiding lines between different spokes.
columns <- unique(as.integer(round(seq(1, ncol(solution), length.out = 81))))
segment_nodes <- head(as.vector(t(cbind(graph$mesh$E, NA_integer_))), -1L)
x <- graph$mesh$V[segment_nodes, 1L]
y <- graph$mesh$V[segment_nodes, 2L]
z_range <- range(c(0, solution))
z_range <- z_range + c(-0.05, 0.05) * max(diff(z_range), 1e-6)
frame_names <- as.character(seq_along(columns))
star_animation <- plotly::plot_ly() |>
  plotly::add_trace(x = x, y = y, z = rep(0, length(x)),
                    type = "scatter3d", mode = "lines", name = "Star graph",
                    line = list(color = "grey65", width = 3), hoverinfo = "skip") |>
  plotly::add_trace(x = x, y = y, z = solution[segment_nodes, 1L],
                    type = "scatter3d", mode = "lines", name = "u(s,t)",
                    line = list(color = "#2166ac", width = 6), connectgaps = FALSE) |>
  plotly::layout(
    title = "Star graph: Gaussian initial condition, f(s,t) = sin(t)",
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
star_animation <- plotly::plotly_build(star_animation)
star_animation$x$frames <- lapply(seq_along(columns), function(i) list(
  # Preserve a JSON array for the zero-based solution trace index.
  name = frame_names[i], traces = I(1L),
  data = list(list(type = "scatter3d", z = solution[segment_nodes, columns[i]]))))
if (interactive()) print(star_animation)
