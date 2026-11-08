# Gaussian diffusion on MetricGraph's logo, using the notebooks' matrix solver.
# Run from RStudio with source("diff_anim.R"), or from a terminal with:
#   /usr/bin/Rscript --vanilla diff_anim.R
# The last block writes ordinary and fractional BRASIL animations, each with
# either one or five initial Gaussian pulses, as self-contained HTML files.
# Named Plotly objects are also saved as .RData here and in the presentation
# project's data_files/ directory, ready for load() in an R Markdown slide.
# Set options(NAP.diff_anim.run = FALSE) to load only functions.

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
            length(center_vertex) >= 1L, all(is.finite(center_vertex)),
            all(center_vertex == as.integer(center_vertex)),
            all(center_vertex >= 1L), all(center_vertex <= nrow(graph$V)),
            !anyDuplicated(center_vertex),
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

  # Sum of intrinsic Gaussians, one per entry of center_vertex:
  # u0(s) = amplitude * sum_v exp(-d_Gamma(s,v)^2/(2*sigma^2)).
  # d_Gamma is shortest-path arclength; sigma is in graph-length units and
  # amplitude is the height of each Gaussian before they are summed.
  # Coordinates do not enter this initial condition.
  # Approximate its L2 projection using MetricGraph's finer mesh, consistent
  # mass matrix and coarse FEM basis, as in the reference-mesh notebooks.
  quadrature_graph <- graph$clone(deep = TRUE)
  quadrature_graph$build_mesh(h = h / projection_refinement)
  quadrature_graph$compute_fem()
  graph$compute_geodist(obs = FALSE)
  locations <- quadrature_graph$get_mesh_locations()
  edge <- locations[, 1L]
  arclength <- locations[, 2L] * graph$edge_lengths[edge]
  gaussian <- numeric(nrow(locations))
  for (vertex in center_vertex) {
    vertex_distance <- graph$geo_dist[[".vertices"]][vertex, ]
    # A path to an edge interior enters through one of its endpoints.
    distance <- pmin(vertex_distance[graph$E[edge, 1L]] + arclength,
                     vertex_distance[graph$E[edge, 2L]] + graph$edge_lengths[edge] - arclength)
    # Unreachable components have distance Inf and contribute zero.
    gaussian <- gaussian + amplitude * exp(-distance^2 / (2 * sigma^2))
  }
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

# Reuse the same physical parameters and an equally fine initial-projection
# mesh when changing the solver mesh. This isolates the effect of spatial h.
solve_coarse_diffusion <- function(result, h = 1) {
  stopifnot(length(h) == 1L, is.finite(h), h > 0)
  parameters <- result$parameters
  projection_h <- parameters$h / parameters$projection_refinement
  parameters$projection_refinement <- as.integer(max(1, ceiling(h / projection_h)))
  parameters$h <- h
  do.call(solve_gaussian_diffusion,
          c(list(graph = result$graph$clone(deep = TRUE)), parameters))
}

animate_graph_diffusion <- function(result, max_frames = 151L, frame_duration = 80,
                                    comparison = NULL) {
  stopifnot(length(max_frames) == 1L, is.finite(max_frames), max_frames >= 2,
            length(frame_duration) == 1L, is.finite(frame_duration), frame_duration > 0)
  graph <- result$graph
  U <- result$solution
  stopifnot(nrow(U) == nrow(graph$mesh$V), ncol(U) == length(result$time))
  columns <- unique(as.integer(round(seq(1L, ncol(U),
                                        length.out = min(max_frames, ncol(U))))))
  times <- result$time[columns]
  U <- U[, columns, drop = FALSE]
  comparison_U <- NULL
  if (!is.null(comparison)) {
    stopifnot(isTRUE(all.equal(result$time, comparison$time)),
              identical(graph$E, comparison$graph$E),
              isTRUE(all.equal(graph$V, comparison$graph$V)))
    # Reconstruct the independently computed coarse P1 solution at the fine
    # plotting locations. It is still solved only on its own h=1 mesh; using
    # common locations follows the curved graph geometry for both traces.
    comparison_U <- as.matrix(comparison$graph$fem_basis(graph$get_mesh_locations()) %*%
                                comparison$solution[, columns, drop = FALSE])
  }

  # Draw each actual mesh segment separately: an NA separates segments, so
  # neither branches nor disconnected components acquire artificial lines.
  edge_nodes <- head(as.vector(t(cbind(graph$mesh$E, NA_integer_))), -1L)
  xy <- graph$mesh$V[, 1:2, drop = FALSE]
  x <- xy[edge_nodes, 1L]
  y <- xy[edge_nodes, 2L]
  support_nodes <- unique(c(graph$mesh$ind,
                           round(seq(1L, nrow(xy), length.out = min(120L, nrow(xy))))))
  support_x <- as.vector(t(cbind(xy[support_nodes, 1L], xy[support_nodes, 1L], NA_real_)))
  support_y <- as.vector(t(cbind(xy[support_nodes, 2L], xy[support_nodes, 2L], NA_real_)))
  # Omit trailing separators, which plotly_build otherwise removes from the
  # initial traces but not from our manually supplied animation frames.
  support_x <- head(support_x, -1L)
  support_y <- head(support_y, -1L)
  support_z <- function(u) head(as.vector(t(cbind(0, u[support_nodes], NA_real_))), -1L)
  z_range <- range(c(0, U, comparison_U))
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
                         type = "scatter3d", mode = "lines",
                         name = if (is.null(comparison)) "Concentration" else
                           sprintf("h = %g", result$parameters$h),
                         connectgaps = FALSE, line = list(color = "#cb432f", width = 5),
                         hovertemplate = "x: %{x:.2f}<br>y: %{y:.2f}<br>u: %{z:.4f}<extra></extra>")
  p <- plotly::add_trace(p, x = support_x, y = support_y, z = support_z(U[, 1L]),
                         type = "scatter3d", mode = "lines", showlegend = FALSE,
                         connectgaps = FALSE, hoverinfo = "skip",
                         line = list(color = "#b0b7c3", width = 1))
  if (!is.null(comparison_U)) {
    p <- plotly::add_trace(p, x = x, y = y, z = comparison_U[edge_nodes, 1L],
                           type = "scatter3d", mode = "lines",
                           name = sprintf("h = %g", comparison$parameters$h),
                           connectgaps = FALSE,
                           line = list(color = "#1769aa", width = 5),
                           hovertemplate = "x: %{x:.2f}<br>y: %{y:.2f}<br>u: %{z:.4f}<extra></extra>")
  }
  p <- plotly::layout(p,
    # Fill the available slide area: no heading/caption space, a closer camera,
    # and a horizontal legend that does not reserve a wide side margin.
    margin = list(l = 0, r = 0, t = 0, b = 80),
    legend = list(orientation = "h", x = 0.02, y = 0.99,
                  xanchor = "left", yanchor = "top"),
    scene = list(xaxis = list(title = "x", range = range(xy[, 1L]), autorange = FALSE),
                 yaxis = list(title = "y", range = range(xy[, 2L]), autorange = FALSE),
                 zaxis = list(title = "Concentration", range = z_range, autorange = FALSE),
                 aspectmode = "manual", aspectratio = list(x = 1.25, y = 1, z = 0.6),
                 camera = list(eye = list(x = 1.0, y = -1.3, z = 0.92),
                               center = list(x = 0, y = 0, z = -0.05)),
                 uirevision = "keep-camera"),
    sliders = list(list(active = 0, x = 0.1, len = 0.87, xanchor = "left",
                        currentvalue = list(prefix = "Time: "),
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
  # Plotly uses zero-based trace indices. Both solutions and the fine-solution
  # supports change; the graph stays in the plane z=0. WebGL needs redraw=TRUE.
  p$x$frames <- lapply(seq_along(times), function(i) {
    traces <- c(1L, 2L)
    data <- list(list(type = "scatter3d", z = U[edge_nodes, i]),
                 list(type = "scatter3d", z = support_z(U[, i])))
    if (!is.null(comparison_U)) {
      traces <- c(traces, 3L)
      data[[3L]] <- list(type = "scatter3d", z = comparison_U[edge_nodes, i])
    }
    list(name = frame_names[i], traces = traces, data = data)
  })
  plotly::config(p, displaylogo = FALSE)
}

# Same named-object workflow as NAFDEMG's save_dual_for_presentation(),
# without sourcing its unrelated plotting/Python setup. The project path is
# anchored to this script so it is independent of the current working directory.
save_dual_for_presentation <- function(obj,
                                       project_dir = file.path(.diff_anim_dir, "data_files"),
                                       presentation_dir = file.path(dirname(.diff_anim_dir),
                                                                    "leninPresentations", "data_files")) {
  object_name <- deparse(substitute(obj))
  stopifnot(inherits(obj, "plotly"), length(object_name) == 1L)
  dir.create(project_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(presentation_dir, recursive = TRUE, showWarnings = FALSE)
  filename <- paste0(object_name, ".RData")
  project_file <- file.path(project_dir, filename)
  presentation_file <- file.path(presentation_dir, filename)
  objects <- new.env(parent = emptyenv())
  assign(object_name, obj, envir = objects)
  save(list = object_name, file = project_file, envir = objects)
  if (!file.copy(project_file, presentation_file, overwrite = TRUE))
    stop("Could not copy the Plotly object to: ", presentation_file)
  message("Saved Plotly object in:\n - ", project_file, "\n - ", presentation_file)
  invisible(c(project = project_file, presentation = presentation_file))
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
  save_dual_for_presentation(logo_diffusion_animation)
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
  save_dual_for_presentation(logo_fractional_diffusion_animation)
  # Two additional experiments: equal-width, unit-amplitude Gaussians on
  # A, M, C, R and H. M uses its downward-pointing notch, and R uses the
  # straight-segment junction in GRAPH. All centers are vertices; coordinates
  # below only identify them, and distances remain intrinsic.
  gaussian_centers <- data.frame(
    letter = c("A", "M", "C", "R", "H"),
    x = c(5, 1.5, 8, 2, 8),
    y = c(4, 5, 6, 2, 2))
  gaussian_center_vertices <- vapply(seq_len(nrow(gaussian_centers)), function(i) {
    vertex <- which(logo_graph$V[, 1L] == gaussian_centers$x[i] &
                      logo_graph$V[, 2L] == gaussian_centers$y[i])
    stopifnot(length(vertex) == 1L)
    as.integer(vertex)
  }, integer(1))
  names(gaussian_center_vertices) <- gaussian_centers$letter

  logo_five_diffusion <- solve_gaussian_diffusion(
    logo_graph$clone(deep = TRUE), h = h, T_final = T_final, time_step = time_step,
    alpha = alpha, kappa = kappa, diffusion = diffusion,
    center_vertex = gaussian_center_vertices, sigma = gaussian_sigma)
  # Overlay an independent h=1 solve, with the same time grid and Gaussian
  # quadrature resolution, on each red h=0.1 solution.
  coarse_h <- 1
  logo_five_coarse_diffusion <- solve_coarse_diffusion(logo_five_diffusion, h = coarse_h)
  logo_five_diffusion_animation <- animate_graph_diffusion(
    logo_five_diffusion, comparison = logo_five_coarse_diffusion)
  five_output_file <- file.path(.diff_anim_dir, "data_files", "logo_diffusion_five_gaussians.html")
  save_diffusion_animation(logo_five_diffusion_animation, five_output_file,
                           "Ordinary graph diffusion: five Gaussian pulses")
  save_dual_for_presentation(logo_five_diffusion_animation)
  message(sprintf("Five-pulse diffusion: relative mass drift %.3g.",
                  max(abs(logo_five_diffusion$diagnostics$mass_fraction - 1))))

  logo_five_fractional_diffusion <- solve_gaussian_diffusion(
    logo_graph$clone(deep = TRUE), h = h, T_final = fractional_T_final,
    time_step = time_step, alpha = fractional_alpha, m = fractional_m,
    kappa = fractional_kappa, diffusion = diffusion,
    center_vertex = gaussian_center_vertices, sigma = gaussian_sigma)
  logo_five_coarse_fractional_diffusion <- solve_coarse_diffusion(
    logo_five_fractional_diffusion, h = coarse_h)
  logo_five_fractional_diffusion_animation <- animate_graph_diffusion(
    logo_five_fractional_diffusion, comparison = logo_five_coarse_fractional_diffusion)
  five_fractional_output_file <- file.path(
    .diff_anim_dir, "data_files", "logo_diffusion_fractional_five_gaussians.html")
  save_diffusion_animation(logo_five_fractional_diffusion_animation, five_fractional_output_file,
                           "Fractional BRASIL graph diffusion: five Gaussian pulses")
  save_dual_for_presentation(logo_five_fractional_diffusion_animation)
  if (interactive()) print(logo_five_diffusion_animation)
}
