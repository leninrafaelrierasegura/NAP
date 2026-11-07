# Gaussian diffusion on MetricGraph's logo, using the notebooks' matrix solver.
# Run from RStudio with source("diff_anim.R"), or from a terminal with:
#   /usr/bin/Rscript --vanilla diff_anim.R
# The last block writes introductory ordinary diffusion and a separate
# fractional BRASIL animation as self-contained HTML files. Set options(NAP.diff_anim.run = FALSE) to load only functions.

.diff_anim_dir <- local({
  source_file <- sys.frame(1)$ofile
  if (is.null(source_file)) {
    argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    source_file <- if (length(argument)) sub("^--file=", "", argument[1]) else "diff_anim.R"
  }
  dirname(normalizePath(source_file, mustWork = TRUE))
})

# Keep the shared notebook functions together in their own environment; no
# notebook computation or regeneration of brasilfunc.R is needed here.
.diff_anim_solver <- new.env(parent = globalenv())
sys.source(file.path(.diff_anim_dir, "brasilfunc.R"), envir = .diff_anim_solver)

solve_gaussian_diffusion <- function(graph, h = 0.1, T_final = 10,
                                     time_step = 0.02, alpha = 2, m = 6L,
                                     kappa = 0, diffusion = 1,
                                     center_vertex = 11L, sigma = 0.7,
                                     amplitude = 1, projection_refinement = 4L) {
  stopifnot(length(h) == 1L, is.finite(h), h > 0,
            length(alpha) == 1L, is.finite(alpha), alpha > 0, alpha <= 2,
            length(kappa) == 1L, is.finite(kappa), kappa >= 0,
            length(diffusion) == 1L, is.finite(diffusion), diffusion > 0,
            length(center_vertex) == 1L, is.finite(center_vertex),
            center_vertex == as.integer(center_vertex),
            center_vertex >= 1L, center_vertex <= nrow(graph$V),
            length(sigma) == 1L, is.finite(sigma), sigma > 0,
            length(amplitude) == 1L, is.finite(amplitude), amplitude > 0,
            length(projection_refinement) == 1L,
            is.finite(projection_refinement), projection_refinement >= 1,
            projection_refinement == as.integer(projection_refinement))
  if (alpha < 2 && kappa == 0)
    stop("The fractional BRASIL solver requires kappa > 0; use alpha = 2 for unshifted heat diffusion.")
  graph$build_mesh(h = h)
  graph$compute_fem()
  C <- graph$mesh$C
  G <- graph$mesh$G
  L <- kappa^2 * C + G
  grid <- .diff_anim_solver$make_time_grid(T_final, time_step)

  # Intrinsic Gaussian: u0(s) = amplitude * exp(-d_Gamma(s,v0)^2/(2*sigma^2)),
  # where d_Gamma is shortest-path arclength on the metric graph. sigma is in
  # graph-length units. Coordinates do not enter this initial condition.
  # Approximate its L2 projection using MetricGraph's finer mesh, consistent
  # mass matrix and coarse FEM basis, as in the reference-mesh notebooks.
  quadrature_graph <- graph$clone(deep = TRUE)
  quadrature_graph$build_mesh(h = h / projection_refinement)
  quadrature_graph$compute_fem()
  graph$compute_geodist(obs = FALSE)
  vertex_distance <- graph$geo_dist[[".vertices"]][center_vertex, ]
  locations <- quadrature_graph$get_mesh_locations()
  edge <- locations[, 1L]
  arclength <- locations[, 2L] * graph$edge_lengths[edge]
  # A path from v0 to an edge interior enters through one of its endpoints.
  distance <- pmin(vertex_distance[graph$E[edge, 1L]] + arclength,
                   vertex_distance[graph$E[edge, 2L]] + graph$edge_lengths[edge] - arclength)
  # Unreachable components have distance Inf and hence exactly zero initial data.
  gaussian <- amplitude * exp(-distance^2 / (2 * sigma^2))
  Psi <- graph$fem_basis(locations)
  initial_load <- Matrix::crossprod(Psi, quadrature_graph$mesh$C %*% gaussian)
  initial <- as.numeric(Matrix::solve(C, initial_load))

  # This is the same factorized sparse BRASIL/backward Euler operator used by
  # conv_in_m_isolated, conv_in_tau_mat, and both cov_in_h notebooks. In
  # particular it retains the BRASIL constant term; alpha=2 uses ordinary BE.
  # For alpha=2 the scale factor cancels algebraically inside the shared
  # operator. Use 1 so kappa=0 is allowed without dividing by zero; then
  # L=G and (C + tau*diffusion*G) U^{k+1} = C U^k conserves graph mass.
  started <- proc.time()[["elapsed"]]
  operator <- .diff_anim_solver$my.fractional.operators.frac(
    L = L, beta = alpha / 2, C = C, scale.factor = if (alpha == 2) 1 else kappa^2,
    m = m, time_step = grid$time_step, diffusion = diffusion)
  setup_seconds <- proc.time()[["elapsed"]] - started
  solution <- matrix(0, nrow(C), length(grid$time_seq))
  solution[, 1L] <- initial
  started <- proc.time()[["elapsed"]]
  for (step in seq_len(grid$N)) {
    # f = 0, so the entire right-hand side is C U^k.
    solution[, step + 1L] <- as.numeric(.diff_anim_solver$my.solver.frac(
      operator, C %*% solution[, step]))
  }
  solve_seconds <- proc.time()[["elapsed"]] - started
  if (any(!is.finite(solution))) stop("The solver produced non-finite values.")
  mass <- Matrix::colSums(C %*% solution)
  l2_norm <- sqrt(Matrix::colSums(solution * (C %*% solution)))
  list(graph = graph, time = grid$time_seq, solution = solution,
       diagnostics = data.frame(time = grid$time_seq, mass = mass,
                                mass_fraction = mass / mass[1L], l2_norm = l2_norm),
       timings = c(factorization = setup_seconds, time_stepping = solve_seconds),
       parameters = list(h = h, time_step = grid$time_step, T_final = T_final,
                         alpha = alpha, m = m, kappa = kappa, diffusion = diffusion,
                         center_vertex = center_vertex, sigma = sigma, amplitude = amplitude,
                         projection_refinement = projection_refinement))
}

animate_graph_diffusion <- function(result, max_frames = 151L, frame_duration = 80) {
  stopifnot(length(max_frames) == 1L, is.finite(max_frames), max_frames >= 2,
            length(frame_duration) == 1L, is.finite(frame_duration), frame_duration > 0)
  graph <- result$graph
  U <- result$solution
  stopifnot(nrow(U) == nrow(graph$mesh$V), ncol(U) == length(result$time))
  columns <- unique(as.integer(round(seq(1L, ncol(U),
                                        length.out = min(max_frames, ncol(U))))))
  times <- result$time[columns]
  U <- U[, columns, drop = FALSE]

  # Draw each actual mesh segment separately: an NA separates segments, so
  # neither branches nor disconnected components acquire artificial lines.
  edge_nodes <- as.vector(t(cbind(graph$mesh$E, NA_integer_)))
  xy <- graph$mesh$V[, 1:2, drop = FALSE]
  x <- xy[edge_nodes, 1L]
  y <- xy[edge_nodes, 2L]
  support_nodes <- unique(c(graph$mesh$ind,
                           round(seq(1L, nrow(xy), length.out = min(120L, nrow(xy))))))
  support_x <- as.vector(t(cbind(xy[support_nodes, 1L], xy[support_nodes, 1L], NA_real_)))
  support_y <- as.vector(t(cbind(xy[support_nodes, 2L], xy[support_nodes, 2L], NA_real_)))
  support_z <- function(u) as.vector(t(cbind(0, u[support_nodes], NA_real_)))
  ordinary <- result$parameters$alpha == 2
  conserving <- ordinary && result$parameters$kappa == 0
  heading <- if (conserving) "Ordinary diffusion on a graph" else if (ordinary)
    "Ordinary diffusion with decay" else "Fractional diffusion with decay"
  explanation <- if (conserving)
    "A pulse spreads along edges and through junctions. Its total amount stays constant." else
    "A pulse spreads along the graph while its total amount decreases."
  title_at <- function(i) sprintf(
    "%s | t = %.3f<br><sup>%s<br>Total amount: %.4f (%.1f%% of initial)</sup>",
    heading, times[i], explanation, result$diagnostics$mass[columns[i]],
    100 * result$diagnostics$mass_fraction[columns[i]])
  z_range <- range(c(0, U))
  z_range <- z_range + c(-0.03, 0.08) * max(diff(z_range), 1e-6)
  frame_names <- sprintf("time_%04d", seq_along(times))
  animation_options <- list(mode = "immediate", fromcurrent = TRUE,
                            frame = list(duration = frame_duration, redraw = TRUE),
                            transition = list(duration = 0))
  steps <- lapply(seq_along(times), function(i) list(
    label = sprintf("%.3f", times[i]), method = "animate",
    args = list(list(frame_names[i]), list(mode = "immediate",
                frame = list(duration = 0, redraw = TRUE), transition = list(duration = 0)))))

  p <- plotly::plot_ly()
  p <- plotly::add_trace(p, x = x, y = y, z = rep(0, length(x)),
                         type = "scatter3d", mode = "lines", name = "Graph",
                         connectgaps = FALSE, hoverinfo = "skip",
                         line = list(color = "#30343b", width = 3))
  p <- plotly::add_trace(p, x = x, y = y, z = U[edge_nodes, 1L],
                         type = "scatter3d", mode = "lines", name = "Concentration",
                         connectgaps = FALSE, line = list(color = "#cb432f", width = 5),
                         hovertemplate = "x: %{x:.2f}<br>y: %{y:.2f}<br>u: %{z:.4f}<extra></extra>")
  p <- plotly::add_trace(p, x = support_x, y = support_y, z = support_z(U[, 1L]),
                         type = "scatter3d", mode = "lines", showlegend = FALSE,
                         connectgaps = FALSE, hoverinfo = "skip",
                         line = list(color = "#b0b7c3", width = 1))
  p <- plotly::layout(p,
    title = list(text = title_at(1L), x = 0.03, font = list(size = 17)),
    margin = list(t = 115, b = 90),
    scene = list(xaxis = list(title = "x", range = range(xy[, 1L]), autorange = FALSE),
                 yaxis = list(title = "y", range = range(xy[, 2L]), autorange = FALSE),
                 zaxis = list(title = "Concentration", range = z_range, autorange = FALSE),
                 aspectmode = "manual", aspectratio = list(x = 1.25, y = 1, z = 0.6),
                 camera = list(eye = list(x = 1.3, y = -1.7, z = 1.2)),
                 uirevision = "keep-camera"),
    sliders = list(list(active = 0, currentvalue = list(prefix = "Time: "),
                        pad = list(t = 25), steps = steps)),
    updatemenus = list(list(type = "buttons", direction = "left", showactive = FALSE,
                            x = 0, y = 0, yanchor = "top", pad = list(t = 65),
                            buttons = list(
                              list(label = "Play", method = "animate",
                                   args = list(NULL, animation_options)),
                              list(label = "Pause", method = "animate",
                                   args = list(list(NULL), list(mode = "immediate",
                                     frame = list(duration = 0, redraw = FALSE),
                                     transition = list(duration = 0))))))))
  p <- plotly::plotly_build(p)
  # Plotly uses zero-based trace indices. Only the solution and its vertical
  # supports change; the graph stays in the plane z=0. WebGL needs redraw=TRUE.
  p$x$frames <- lapply(seq_along(times), function(i) list(
    name = frame_names[i], traces = c(1L, 2L),
    data = list(list(type = "scatter3d", z = U[edge_nodes, i]),
                list(type = "scatter3d", z = support_z(U[, i]))),
    layout = list(title = list(text = title_at(i)))))
  plotly::config(p, displaylogo = FALSE)
}

# Save only the HTML; temporary dependency copies are removed after embedding.
save_diffusion_animation <- function(animation, output_file, title) {
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  widget_libdir <- tempfile("logo_diffusion_libs_", tmpdir = dirname(output_file))
  on.exit(unlink(widget_libdir, recursive = TRUE), add = TRUE)
  htmlwidgets::saveWidget(animation, output_file, selfcontained = TRUE,
                          libdir = widget_libdir, title = title)
  message("Animation saved to: ", output_file)
}

if (isTRUE(getOption("NAP.diff_anim.run", TRUE))) {
  library(MetricGraph)
  # First demonstration: u_t = diffusion * Delta_Gamma u, with Kirchhoff
  # conditions and zero forcing. alpha=2, kappa=0 means no decay of total mass.
  # The longer interval makes spreading into neighboring branches visible.
  h <- 0.1
  time_step <- 0.02
  T_final <- 10
  alpha <- 2
  kappa <- 0
  diffusion <- 1
  gaussian_sigma <- 0.7

  edges <- logo_lines()
  logo_graph <- metric_graph$new(edges = edges, perform_merges = TRUE)
  # Top of the A in GRAPH: the vertex at (5, 4), currently vertex 11.
  gaussian_center_vertex <- which(logo_graph$V[, 1L] == 5 & logo_graph$V[, 2L] == 4)
  stopifnot(length(gaussian_center_vertex) == 1L)
  logo_diffusion <- solve_gaussian_diffusion(
    logo_graph, h = h, T_final = T_final, time_step = time_step,
    alpha = alpha, kappa = kappa, diffusion = diffusion,
    center_vertex = gaussian_center_vertex, sigma = gaussian_sigma)
  logo_diffusion_animation <- animate_graph_diffusion(logo_diffusion)
  output_file <- file.path(.diff_anim_dir, "data_files", "logo_diffusion.html")
  save_diffusion_animation(logo_diffusion_animation, output_file,
                           "Ordinary diffusion on the MetricGraph logo")
  message(sprintf("Ordinary diffusion: %d unknowns, %d time steps; relative mass drift %.3g.",
                  nrow(logo_diffusion$solution), ncol(logo_diffusion$solution) - 1L,
                  max(abs(logo_diffusion$diagnostics$mass_fraction - 1))))

  # Second demonstration: preserve the previous fractional BRASIL experiment.
  # It starts from the same intrinsic Gaussian; alpha=1, kappa=1 adds decay.
  fractional_alpha <- 1
  fractional_kappa <- 1
  fractional_m <- 6L
  fractional_T_final <- 2
  logo_fractional_diffusion <- solve_gaussian_diffusion(
    logo_graph$clone(deep = TRUE), h = h, T_final = fractional_T_final,
    time_step = time_step, alpha = fractional_alpha, m = fractional_m,
    kappa = fractional_kappa, diffusion = diffusion,
    center_vertex = gaussian_center_vertex, sigma = gaussian_sigma)
  logo_fractional_diffusion_animation <- animate_graph_diffusion(logo_fractional_diffusion)
  fractional_output_file <- file.path(.diff_anim_dir, "data_files", "logo_diffusion_fractional.html")
  save_diffusion_animation(logo_fractional_diffusion_animation, fractional_output_file,
                           "Fractional BRASIL diffusion on the MetricGraph logo")
  if (interactive()) print(logo_diffusion_animation)
}
