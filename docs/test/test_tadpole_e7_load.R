# Run from the project root: Rscript --vanilla test/test_tadpole_e7_load.R
# Load only function definitions, without the plotting/package setup.
helpers <- new.env(parent = baseenv())
for (expr in parse("fun/matrix_reference/graphs_and_meshes.R")) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      is.call(expr[[3]]) && identical(expr[[3]][[1]], as.name("function"))) {
    eval(expr, helpers)
  }
}

# Each fixture specifies global nodes independently of the assembly routine.
# PtE stores edge numbers and relative positions; vertices 1 and 2 are the
# stem's free endpoint and the common junction, respectively.
fixtures <- list(
  no_interior_nodes = list(
    PtE = matrix(numeric(0), ncol = 2),
    edges = list(list(x = c(0, 1), nodes = c(1L, 2L)),
                 list(x = c(0, 2), nodes = c(2L, 2L)))),
  nonuniform_shuffled_nodes = list(
    PtE = rbind(c(2, .65), c(1, .15), c(2, .1), c(1, .8), c(2, .4)),
    edges = list(list(x = c(0, .15, .8, 1), nodes = c(1L, 4L, 6L, 2L)),
                 list(x = c(0, .2, .8, 1.3, 2), nodes = c(2L, 5L, 7L, 3L, 2L)))),
  small_intervals = list(
    PtE = rbind(c(1, 1e-6), c(2, .5), c(1, .3), c(2, 1 - 1e-6)),
    edges = list(list(x = c(0, 1e-6, .3, 1), nodes = c(1L, 3L, 5L, 2L)),
                 list(x = c(0, 1, 2 - 2e-6, 2), nodes = c(2L, 4L, 6L, 2L))))
)

for (name in names(fixtures)) {
  fixture <- fixtures[[name]]
  n <- 2L + nrow(fixture$PtE)
  graph <- list(mesh = list(C = diag(n), PtE = fixture$PtE),
                get_edge_lengths = function() c(1, 2))
  actual <- helpers$assemble_e7_load_tadpole(graph)

  # Independently interpolate each global hat on each edge, then integrate
  # e7 times that hat numerically. Split at knots to resolve every corner.
  expected <- vapply(seq_len(n), function(i) {
    sum(vapply(fixture$edges, function(edge) {
      hat <- stats::approxfun(edge$x, as.numeric(edge$nodes == i))
      sum(vapply(seq_len(length(edge$x) - 1L), function(j) {
        stats::integrate(function(x) sqrt(2/3) * cos(2*pi*x) * hat(x),
                         edge$x[j], edge$x[j + 1L],
                         abs.tol = 1e-13, rel.tol = 1e-11)$value
      }, numeric(1)))
    }, numeric(1)))
  }, numeric(1))

  max_error <- max(abs(actual - expected))
  stopifnot(length(actual) == n, all(is.finite(actual)), max_error < 1e-12)
  # Partition of unity: the sum equals the integral of e7 over the graph (0).
  stopifnot(abs(sum(actual)) < 1e-12)
  cat(sprintf("%s: maximum load error %.3g\n", name, max_error))
}

# Check normalization against the eigenfunction enumeration used by the notebooks.
graph <- list(mesh = list(C = diag(5),
                          PtE = rbind(c(1, .25), c(2, .2), c(2, .7))),
              get_edge_lengths = function() c(1, 2))
e7 <- helpers$gets.eigen.params(N_finite = 4, kappa = 4, graph = graph)$EIGENFUN[, 7]
stopifnot(max(abs(e7 - sqrt(2/3) * cos(2*pi*c(0, 0, .25, .4, 1.4)))) < 1e-14)
cat("Seventh-mode normalization and analytic load assembly passed.\n")
