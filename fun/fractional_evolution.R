## ------------------------------------------------------------------------------------------------------------------------------------------------
rational_supported_orders <- function(coefficient_source = "package") {
  coefficient_source <- match.arg(coefficient_source, c("package", "computed"))
  if (coefficient_source == "package") 1:4 else 1:8
}


## ------------------------------------------------------------------------------------------------------------------------------------------------
read_computed_fractional_roots <- function() {
  path <- here::here("matlab", "computed_fractional_roots.csv")
  if (!file.exists(path)) stop("Missing computed coefficients: ", path)
  tab <- read.csv(path, stringsAsFactors = FALSE)
  if (!identical(names(tab), c("m", "beta", "kind", "index", "value")) ||
      !all(vapply(tab[c("m", "beta", "index", "value")], is.numeric, logical(1))) ||
      !all(is.finite(as.matrix(tab[c("m", "beta", "index", "value")]))) ||
      !all(tab$m %in% 1:8) || !all(tab$beta >= .25 & tab$beta < 1) ||
      !all(tab$kind %in% c("factor", "rc", "rb")) ||
      anyDuplicated(tab[c("m", "beta", "kind", "index")]))
    stop("Malformed computed coefficient artifact.")
  tab
}


## ------------------------------------------------------------------------------------------------------------------------------------------------
validate_fractional_roots <- function(roots, m) {
  all_roots <- c(roots$pr_roots, roots$pl_roots)
  if (length(roots$pr_roots) != m || length(roots$pl_roots) != m+1L ||
      length(roots$factor) != 1L || !is.finite(roots$factor) || roots$factor <= 0 ||
      !is.numeric(all_roots) || is.complex(all_roots) ||
      any(!is.finite(all_roots)) || any(all_roots >= 0))
    stop("Fractional coefficients need a positive factor and real negative roots.")
  numerator <- sort(roots$pr_roots)
  denominator <- sort(roots$pl_roots)
  if (any(denominator[-length(denominator)] >= numerator) ||
      any(numerator >= denominator[-1L]))
    stop("Fractional coefficient roots do not interlace; this approximation is unsafe.")
  roots
}


## ------------------------------------------------------------------------------------------------------------------------------------------------
my.get.roots <- function(m,
                         beta, # smoothness parameter, beta = alpha/2
                         interpolation = c("linear", "cubic"),
                         coefficient_source = "package"
                         ) {
  coefficient_source <- match.arg(coefficient_source, c("package", "computed"))
  interpolation <- match.arg(interpolation)
  orders <- rational_supported_orders(coefficient_source)
  if (length(m) != 1L || !is.finite(m) || !(m %in% orders))
    stop("Coefficient source '", coefficient_source, "' supports m = ",
         paste(orders, collapse = ", "), " only.")
  if (length(beta) != 1L || !is.finite(beta) || beta < 0.25 || beta >= 1)
    stop("Fractional evolution tables require 0.25 <= beta < 1; beta = 1 uses the exact solver.")
  if (coefficient_source == "computed") {
    tab <- read_computed_fractional_roots()
    rows <- tab[tab$m == m & abs(tab$beta-beta) <= 32*.Machine$double.eps*max(1, abs(beta)), ]
    if (!nrow(rows)) stop("Computed coefficients are unavailable at beta = ",
      format(beta, digits = 16), ". Available beta: ",
      paste(sort(unique(tab$beta)), collapse = ", "),
      ". Generate and validate new coefficients; beta interpolation is disabled.")
    values <- function(kind, indices) {
      part <- rows[rows$kind == kind, ]
      if (nrow(part) != length(indices) || !setequal(part$index, indices))
        stop("Computed coefficient entries are incomplete for m = ", m, ".")
      part$value[match(indices, part$index)]
    }
    return(validate_fractional_roots(list(factor = values("factor", 0),
      pr_roots = values("rc", seq_len(m)), pl_roots = values("rb", seq_len(m+1L))), m))
  }
  mt <- get(paste0("m", m, "table"), envir = asNamespace("rSPDE"))
  interpolate <- function(column) {
    if (interpolation == "linear") approx(mt$beta, mt[[column]], beta)$y else
      splinefun(mt$beta, mt[[column]], method = "fmm")(beta)
  }
  rb <- rep(0, m + 1)
  rc <- rep(0, m)
  if(m == 1) {
    rc = interpolate("rc")
  } else {
    rc = sapply(1:m, function(i) {
      interpolate(paste0("rc.", i))
    })
  }
  rb = sapply(1:(m+1), function(i) {
    interpolate(paste0("rb.", i))
  })
  factor = interpolate("factor")
  return(validate_fractional_roots(list(pl_roots = rb, # roots \{r_{2j}\}_{j=1}^{m+1}
              pr_roots = rc, # roots \{r_{1i}\}_{i=1}^m
              factor = factor # this is c_m/b_{m+1}
              ), m))
}

