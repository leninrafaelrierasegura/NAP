# Pure (zero-shift) fractional diffusion of five intrinsic Gaussians.
# Run with the R installation containing MetricGraph and rSPDE:
#   /usr/bin/Rscript --vanilla true_frac_diff.R
# or source("true_frac_diff.R") in that RStudio session.
# Set options(NAP.true_frac_diff.run = FALSE) to define functions only.
#
# Equation: u_t + nu*(-Delta_Gamma)^(alpha/2) u = 0, kappa = 0.
# Continuity and Kirchhoff conditions; beta = alpha/2. No forcing or decay
# of total mass. Spatial variations smooth toward each component's mean.
# A graph may be disconnected, with one constant null mode per component.
# Apply kappa0func.Rmd separately on each connected block; the installed
# MetricGraph logo assembles to one connected component.
#
# As in diff_anim.R: geodesic Gaussians, refined-mesh L2 projection,
# consistent P1 FEM, backward Euler, cached BRASIL shifted sparse solves,
# fixed temperature colors, a time slider, and h=0.1 / h=1 overlay.
# New outputs are confined to this project's data_files directory;
# no original animation, notebook, or presentation file is overwritten.

.true_frac_dir <- local({
  source_file <- sys.frame(1)$ofile
  if (is.null(source_file)) {
    argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    source_file <- if (length(argument)) sub("^--file=", "", argument[1]) else "true_frac_diff.R"
  }
  dirname(normalizePath(source_file, mustWork = TRUE))
})

for (.true_frac_package in c("Matrix", "MetricGraph", "rSPDE", "knitr", "plotly", "htmlwidgets")) {
  if (!requireNamespace(.true_frac_package, quietly = TRUE))
    stop("Missing package ", .true_frac_package,
         ". Use /usr/bin/Rscript or an R session with the animation dependencies installed.")
}
rm(.true_frac_package)

# Extract function definitions only, into a temporary file. The notebook's
# demonstration chunks have purl=FALSE. Never regenerate brasilfunc.R or
# leave a stale kappa0func.R next to the notebook. Use a private environment.
.true_frac_solver <- local({
  extracted <- tempfile("kappa0_functions_", fileext = ".R")
  on.exit(unlink(extracted), add = TRUE)
  knitr::purl(file.path(.true_frac_dir, "kappa0func.Rmd"), output = extracted,
              quiet = TRUE, documentation = 0L)
  env <- new.env(parent = globalenv())
  sys.source(extracted, envir = env)
  env
})

solve_true_fractional_diffusion <- function(graph, h = 0.1, T_final = 10,
                                     time_step = 0.02, alpha = 1, m = 12L,
                                     diffusion = 1,
                                     center_vertex = 11L, sigma = 0.7,
                                     amplitude = 1, projection_refinement = 4L) {
  stopifnot(length(h) == 1L, is.finite(h), h > 0,
            length(alpha) == 1L, is.finite(alpha), alpha > 0, alpha <= 2,
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
  graph$build_mesh(h = h)
  graph$compute_fem()
  C <- graph$mesh$C
  G <- graph$mesh$G
  grid <- .true_frac_solver$kappa0_time_grid(T_final, time_step)

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

  # The logo need not be connected. Each component has its own null mode,
  # certified gap bound, and cached copy of the kappa0func.Rmd operator.
  # D_v + max(edge length) bounds the metric diameter: any interior point
  # is at most half its edge length from a vertex. Total component length
  # also bounds the diameter, so take the smaller of these two bounds.
  started <- proc.time()[["elapsed"]]
  vertex_distance <- graph$geo_dist[[".vertices"]]
  vertex_component <- integer(nrow(graph$V))
  component <- 0L
  for (v in seq_len(nrow(graph$V))) if (vertex_component[v] == 0L) {
    component <- component + 1L
    vertex_component[is.finite(vertex_distance[v, ])] <- component
  }
  node_edge <- as.integer(graph$get_mesh_locations()[, 1L])
  node_component <- vertex_component[graph$E[node_edge, 1L]]
  operators <- lapply(seq_len(component), function(id) {
    vertices <- which(vertex_component == id)
    edges <- which(vertex_component[graph$E[, 1L]] == id)
    nodes <- which(node_component == id)
    volume <- sum(graph$edge_lengths[edges])
    diameter_bound <- min(volume, max(vertex_distance[vertices, vertices]) +
                                  max(graph$edge_lengths[edges]))
    gamma_star <- 2/(volume * diameter_bound)
    op <- .true_frac_solver$my.fractional.operators.kappa0(
      G = G[nodes, nodes, drop = FALSE], beta = alpha/2,
      C = C[nodes, nodes, drop = FALSE], scale.factor = gamma_star,
      m = m, time_step = grid$time_step, diffusion = diffusion)
    if (abs(op$volume - volume) > 1e-10 * max(1, volume))
      stop("Component mesh volume does not match its physical edge lengths.")
    list(nodes = nodes, operator = op, volume = volume,
         diameter_bound = diameter_bound, gamma_star = gamma_star)
  })
  setup_seconds <- proc.time()[["elapsed"]] - started
  solution <- matrix(0, nrow(C), length(grid$times))
  solution[, 1L] <- initial
  started <- proc.time()[["elapsed"]]
  for (step in seq_len(grid$N)) {
    # f = 0. Each block conserves its component mass; there is no kappa
    # parameter and no positive shift hidden in the normalization scale.
    for (block in operators) {
      nodes <- block$nodes
      solution[nodes, step + 1L] <- .true_frac_solver$my.solver.kappa0(
        block$operator, block$operator$C %*% solution[nodes, step])
    }
  }
  solve_seconds <- proc.time()[["elapsed"]] - started
  if (any(!is.finite(solution))) stop("The solver produced non-finite values.")
  mass <- Matrix::colSums(C %*% solution)
  l2_norm <- sqrt(Matrix::colSums(solution * (C %*% solution)))
  component_mass <- vapply(operators, function(block)
    as.numeric(Matrix::colSums(block$operator$C %*% solution[block$nodes, , drop = FALSE])),
    numeric(ncol(solution)))
  colnames(component_mass) <- paste0("component_", seq_along(operators))
  component_drift <- sweep(component_mass, 2L, component_mass[1L, ], "-")
  relative_component_drift <- abs(component_drift) /
    rep(pmax(1, abs(component_mass[1L, ])), each = nrow(component_mass))
  if (max(relative_component_drift) > 1e-9)
    stop("Component mass balance failed: drift = ", max(relative_component_drift))
  component_info <- do.call(rbind, lapply(seq_along(operators), function(i) {
    block <- operators[[i]]
    data.frame(component = i, nodes = length(block$nodes), volume = block$volume,
               diameter_bound = block$diameter_bound, gamma_star = block$gamma_star,
               initial_mass = component_mass[1L, i],
               equilibrium_mean = component_mass[1L, i]/block$volume,
               max_absolute_mass_drift = max(abs(component_drift[, i])))
  }))
  # The equilibrium is constant on each connected component.
  centered <- solution
  for (i in seq_along(operators)) {
    nodes <- operators[[i]]$nodes
    centered[nodes, ] <- centered[nodes, , drop = FALSE] - component_info$equilibrium_mean[i]
  }
  centered_l2 <- sqrt(Matrix::colSums(centered * (C %*% centered)))
  if (any(diff(centered_l2) > 1e-10 * max(1, centered_l2[1L])))
    stop("The centered L2 norm unexpectedly increased under unforced diffusion.")
  list(graph = graph, time = grid$times, solution = solution,
       component_info = component_info, component_mass = component_mass,
       diagnostics = data.frame(time = grid$times, mass = mass,
                                mass_fraction = mass / mass[1L], l2_norm = l2_norm,
                                centered_l2_norm = centered_l2),
       timings = c(factorization = setup_seconds, time_stepping = solve_seconds),
       parameters = list(h = h, time_step = grid$time_step, T_final = T_final,
                         alpha = alpha, m = m, kappa = 0, diffusion = diffusion,
                         center_vertex = center_vertex, sigma = sigma, amplitude = amplitude,
                         projection_refinement = projection_refinement))
}

animate_true_fractional_diffusion <- function(result, max_frames = 151L, frame_duration = 80,
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
  # One fixed scale for both meshes and every frame: blue is cold, yellow
  # intermediate, and red hot. Units are relative to the initial amplitude.
  temperature_colors <- list(list(0, "#2166ac"), list(0.5, "#ffe34d"),
                             list(1, "#d7191c"))
  temperature_max <- max(c(0, U, comparison_U))
  temperature_line <- function(values, width, showscale = FALSE) list(
    color = values, colorscale = temperature_colors,
    cmin = 0, cmax = temperature_max, cauto = FALSE,
    width = width, showscale = showscale,
    colorbar = list(title = list(text = "Temperature<br>(relative units)", side = "right"),
                    x = 1.01, y = 0.55, len = 0.75, thickness = 12,
                    outlinewidth = 0, tickfont = list(size = 11)))
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
                         name = if (is.null(comparison)) "Temperature" else
                           sprintf("h = %g", result$parameters$h),
                         legendgroup = "fine", connectgaps = FALSE,
                         line = temperature_line(U[edge_nodes, 1L], width = 7, showscale = TRUE),
                         hovertemplate = "x: %{x:.2f}<br>y: %{y:.2f}<br>Temperature: %{z:.4f}<extra></extra>")
  p <- plotly::add_trace(p, x = support_x, y = support_y, z = support_z(U[, 1L]),
                         type = "scatter3d", mode = "lines", showlegend = FALSE,
                         connectgaps = FALSE, hoverinfo = "skip",
                         line = list(color = "#b0b7c3", width = 1))
  if (!is.null(comparison_U)) {
    p <- plotly::add_trace(p, x = x, y = y, z = comparison_U[edge_nodes, 1L],
                           type = "scatter3d", mode = "lines",
                           name = sprintf("h = %g (nodes marked)", comparison$parameters$h),
                           legendgroup = "coarse", connectgaps = FALSE,
                           line = temperature_line(comparison_U[edge_nodes, 1L], width = 3),
                           hovertemplate = "x: %{x:.2f}<br>y: %{y:.2f}<br>Temperature: %{z:.4f}<extra></extra>")
    # Temperature determines color for both curves. Mark the actual coarse
    # FEM nodes so the h=1 solution remains identifiable without using dashes.
    p <- plotly::add_trace(p, x = comparison$graph$mesh$V[, 1L],
                           y = comparison$graph$mesh$V[, 2L],
                           z = comparison$solution[, columns[1L]],
                           type = "scatter3d", mode = "markers",
                           legendgroup = "coarse", showlegend = FALSE,
                           marker = list(size = 3, color = comparison$solution[, columns[1L]],
                             colorscale = temperature_colors, cmin = 0, cmax = temperature_max,
                             cauto = FALSE, showscale = FALSE,
                             line = list(color = "#202124", width = 1)),
                           hovertemplate = "x: %{x:.2f}<br>y: %{y:.2f}<br>Temperature: %{z:.4f}<extra>coarse node</extra>")
  }
  p <- plotly::layout(p,
    # Keep the close view and compact legend; reserve space for the colorbar.
    margin = list(l = 0, r = 95, t = 0, b = 80),
    annotations = list(list(
      x = 0.98, y = 0.98, xref = "paper", yref = "paper", showarrow = FALSE,
      xanchor = "right", yanchor = "top",
      text = sprintf(paste0("Pure fractional diffusion &nbsp; \u03ba = 0 &nbsp; \u03b1 = %g",
                            "<br>Mass conserved on each component (relative drift %.1e)"),
                     result$parameters$alpha,
                     max(abs(result$diagnostics$mass_fraction - 1))),
      font = list(size = 17, color = "#202124"),
      bgcolor = "rgba(255,255,255,0.85)", borderpad = 4)),
    legend = list(orientation = "h", x = 0.02, y = 0.99,
                  xanchor = "left", yanchor = "top"),
    scene = list(xaxis = list(title = "x", range = range(xy[, 1L]), autorange = FALSE),
                 yaxis = list(title = "y", range = range(xy[, 2L]), autorange = FALSE),
                 zaxis = list(title = "Temperature", range = z_range, autorange = FALSE),
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
    data <- list(list(type = "scatter3d", z = U[edge_nodes, i],
                      line = list(color = U[edge_nodes, i])),
                 list(type = "scatter3d", z = support_z(U[, i])))
    if (!is.null(comparison_U)) {
      traces <- c(traces, 3L, 4L)
      data[[3L]] <- list(type = "scatter3d", z = comparison_U[edge_nodes, i],
                        line = list(color = comparison_U[edge_nodes, i]))
      data[[4L]] <- list(type = "scatter3d", z = comparison$solution[, columns[i]],
                        marker = list(color = comparison$solution[, columns[i]]))
    }
    list(name = frame_names[i], traces = traces, data = data)
  })
  plotly::config(p, displaylogo = FALSE)
}

# Save only the HTML; temporary dependency copies are removed after embedding.
save_true_fractional_animation <- function(animation, output_file, title) {
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  widget_libdir <- tempfile("true_frac_diff_libs_", tmpdir = dirname(output_file))
  on.exit(unlink(widget_libdir, recursive = TRUE), add = TRUE)
  htmlwidgets::saveWidget(animation, output_file, selfcontained = TRUE,
                          libdir = widget_libdir, title = title)
  message("Animation saved to: ", output_file)
}

# Same five centers as logo_diffusion_fractional_five_gaussians.html.
# These coordinates only select vertices; Gaussian distance is arc length.
run_true_fractional_demo <- function(h = 0.1, coarse_h = 1, T_final = 10,
                                     time_step = 0.02, alpha = 1, m = 12L,
                                     diffusion = 1, sigma = 0.7,
                                     output_dir = file.path(.true_frac_dir, "data_files")) {
  graph <- MetricGraph::metric_graph$new(edges = MetricGraph::logo_lines(),
                                         perform_merges = TRUE)
  centers <- data.frame(letter = c("A", "M", "C", "R", "H"),
                        x = c(5, 1.5, 8, 2, 8), y = c(4, 5, 6, 2, 2))
  vertices <- vapply(seq_len(nrow(centers)), function(i) {
    v <- which(abs(graph$V[, 1L] - centers$x[i]) < 1e-10 &
               abs(graph$V[, 2L] - centers$y[i]) < 1e-10)
    stopifnot(length(v) == 1L)
    as.integer(v)
  }, integer(1))
  names(vertices) <- centers$letter
  # m=12 resolves the rational response over the conservative component gap
  # scales. It remains editable independently of alpha, h, tau, and nu.
  fine <- solve_true_fractional_diffusion(graph$clone(deep = TRUE), h = h,
      T_final = T_final, time_step = time_step, alpha = alpha, m = m,
      diffusion = diffusion, center_vertex = vertices, sigma = sigma)
  coarse <- solve_true_fractional_diffusion(graph$clone(deep = TRUE), h = coarse_h,
      T_final = T_final, time_step = time_step, alpha = alpha, m = m,
      diffusion = diffusion, center_vertex = vertices, sigma = sigma,
      projection_refinement = as.integer(max(1, ceiling(coarse_h/(h/4)))))
  animation <- animate_true_fractional_diffusion(fine, comparison = coarse)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  html_file <- file.path(output_dir, "true_frac_diff_five_gaussians.html")
  save_true_fractional_animation(animation, html_file,
      "Pure fractional diffusion: five Gaussians, zero shift, conserved mass")
  # A named Plotly object, ready for load() and inclusion in an Rmd slide.
  true_frac_diff_animation <- animation
  save(true_frac_diff_animation,
       file = file.path(output_dir, "true_frac_diff_animation.RData"))
  diagnostics <- data.frame(fine$diagnostics, fine$component_mass)
  utils::write.csv(diagnostics, file.path(output_dir, "true_frac_diff_mass.csv"), row.names = FALSE)
  utils::write.csv(fine$component_info, file.path(output_dir, "true_frac_diff_components.csv"), row.names = FALSE)
  message(sprintf(paste0("Pure fractional diffusion: alpha=%g, kappa=0; %d components; ",
                         "%d fine unknowns, %d coarse unknowns; %d time steps; relative mass drift %.3g."),
                  alpha, nrow(fine$component_info), nrow(fine$solution), nrow(coarse$solution),
                  length(fine$time)-1L, max(abs(fine$diagnostics$mass_fraction - 1))))
  print(fine$component_info)
  if (interactive()) print(animation)
  invisible(list(fine = fine, coarse = coarse, animation = animation, html_file = html_file))
}

if (isTRUE(getOption("NAP.true_frac_diff.run", TRUE))) {
  true_frac_diff <- run_true_fractional_demo()
  true_frac_diff_animation <- true_frac_diff$animation
}
