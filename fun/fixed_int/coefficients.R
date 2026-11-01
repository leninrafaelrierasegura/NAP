# Fresh BRASIL fits on one common interval; no package tables or interpolation.
fixed_int_validate_coefficients <- function(coefficients) {
  cf <- coefficients
  values <- c(cf$poles, cf$residues, cf$constant)
  if (length(cf$m) != 1L || !is.finite(cf$m) || cf$m < 1L || cf$m != floor(cf$m) ||
      length(cf$beta) != 1L || !is.finite(cf$beta) || cf$beta <= 0 || cf$beta >= 1 ||
      length(cf$poles) != cf$m || length(cf$residues) != cf$m || length(cf$constant) != 1L ||
      !is.numeric(values) || is.complex(values) || any(!is.finite(values)) ||
      any(cf$poles >= 0) || any(cf$residues <= 0) || cf$constant < 0 ||
      any(diff(cf$poles) <= 0))
    stop("BRASIL coefficients require ordered distinct negative poles, positive residues and k >= 0.")
  invisible(cf)
}

fixed_int_interval <- function(spectrum, kappa, x_min = NULL, diagnostic_y_max = 1e3) {
  stopifnot(length(kappa) == 1L, is.finite(kappa), kappa > 0,
    length(diagnostic_y_max) == 1L, is.finite(diagnostic_y_max), diagnostic_y_max > 1,
    length(spectrum$values) > 0L, all(is.finite(spectrum$values)), all(spectrum$values > 0))
  normalized <- spectrum$values/kappa^2
  if (min(normalized) < 1-1e-8)
    stop("The normalized FEM spectrum must have lower endpoint 1 (up to eigensolver roundoff).")
  required_min <- min(1/max(normalized), 1/diagnostic_y_max)
  if (is.null(x_min)) x_min <- required_min
  if (length(x_min) != 1L || !is.finite(x_min) || x_min <= 0 || x_min >= 1)
    stop("fixed_x_min must be NULL or a number strictly between 0 and 1.")
  if (x_min > required_min*(1+32*.Machine$double.eps))
    stop("The requested interval misses the FEM spectrum or scalar diagnostic. Use fixed_x_min <= ",
         format(required_min, digits = 17), ".")
  list(x_min = x_min, x_max = 1, y_min = 1, y_max = 1/x_min,
       spectral_y_min = min(normalized), spectral_y_max = max(normalized),
       diagnostic_y_max = diagnostic_y_max)
}

fixed_int_generate <- function(betas, orders, interval, digits = 70L, tolerance = 1e-8,
                               max_iterations = 2000L, workers = 2L,
                               python = Sys.getenv("NAP_BRASIL_PYTHON", unset = "")) {
  stopifnot(is.numeric(betas), all(is.finite(betas)), all(betas > 0 & betas < 1),
    !anyDuplicated(betas), length(orders) > 0L, all(is.finite(orders)),
    all(orders >= 1 & orders == floor(orders)), !anyDuplicated(orders),
    length(digits) == 1L, is.finite(digits), digits >= 30, digits == floor(digits),
    length(tolerance) == 1L, is.finite(tolerance), tolerance > 0, tolerance < 1,
    length(max_iterations) == 1L, is.finite(max_iterations), max_iterations >= 1,
    max_iterations == floor(max_iterations), length(workers) == 1L, is.finite(workers),
    workers >= 1, workers == floor(workers), is.list(interval),
    is.finite(interval$x_min), interval$x_min > 0, interval$x_min < 1)
  bank <- list(interval = interval, betas = betas, orders = as.integer(orders),
    coefficients = list(), diagnostics = data.frame(), toolchain = NULL)
  # alpha=2 needs no fractional approximation or Python environment.
  if (!length(betas)) return(bank)
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("R package jsonlite is required.")
  if (!nzchar(python)) {
    local_python <- here::here(".venv-fixed-int", if (.Platform$OS.type == "windows")
      "Scripts/python.exe" else "bin/python")
    python <- if (file.exists(local_python)) local_python else unname(Sys.which("python3"))
  }
  if (!file.exists(python)) python <- unname(Sys.which(python))
  if (!nzchar(python) || !file.exists(python))
    stop("Python was not found. See fun/fixed_int/README.md or set NAP_BRASIL_PYTHON.")
  generator <- here::here("fun", "fixed_int", "generate_coefficients.py")
  files <- c(request = tempfile(fileext = ".json"), output = tempfile(fileext = ".json"),
             log = tempfile(fileext = ".log"))
  on.exit(unlink(files), add = TRUE)
  status <- system2(python, c(shQuote(generator), "--check"), stdout = files["log"], stderr = files["log"])
  if (status != 0L)
    stop("BRASIL Python dependencies are unavailable. See fun/fixed_int/README.md.\n",
         paste(tail(readLines(files["log"], warn = FALSE), 10L), collapse = "\n"))
  request <- list(schema_version = 1L,
    x_min = format(interval$x_min, digits = 17, scientific = TRUE, trim = TRUE),
    betas = I(vapply(betas, format, character(1), digits = 17, scientific = TRUE, trim = TRUE)),
    orders = I(as.integer(orders)), digits = as.integer(digits), tolerance = tolerance,
    max_iterations = as.integer(max_iterations))
  jsonlite::write_json(request, files["request"], auto_unbox = TRUE, digits = 17)
  message("Computing ", length(betas)*length(orders), " BRASIL fits on [",
          format(interval$x_min, digits = 8), ", 1] at ", digits, " decimal digits.")
  status <- system2(python, c(shQuote(generator), "--request", shQuote(files["request"]),
    "--output", shQuote(files["output"]), "--workers", as.integer(workers)),
    stdout = files["log"], stderr = files["log"])
  if (status != 0L || !file.exists(files["output"]))
    stop("Fixed-interval BRASIL fitting failed:\n",
         paste(tail(readLines(files["log"], warn = FALSE), 20L), collapse = "\n"))
  generated <- jsonlite::read_json(files["output"], simplifyVector = FALSE)
  if (!identical(generated$schema_version, 1L) ||
      length(generated$fits) != length(betas)*length(orders)) stop("Malformed coefficient-generation result.")
  bank$toolchain <- generated$toolchain
  bank$request <- generated$request
  diagnostic_rows <- list()
  for (fit in generated$fits) {
    cf <- list(m = as.integer(fit$m), beta = as.numeric(fit$beta),
      poles = as.numeric(unlist(fit$poles)), residues = as.numeric(unlist(fit$residues)),
      constant = fit$constant, x_min = as.numeric(fit$x_min),
      high_precision = fit$high_precision)
    fixed_int_validate_coefficients(cf)
    if (!cf$m %in% orders || !any(abs(betas-cf$beta) <= 32*.Machine$double.eps) ||
        cf$x_min != interval$x_min || !isTRUE(fit$diagnostics$converged))
      stop("Generated coefficients do not match the requested orders, exponents or fixed interval.")
    bank$coefficients[[length(bank$coefficients)+1L]] <- cf
    diagnostic_rows[[length(diagnostic_rows)+1L]] <- data.frame(
      beta = cf$beta, m = cf$m, x_min = cf$x_min, constant = cf$constant,
      as.data.frame(fit$diagnostics))
  }
  bank$diagnostics <- do.call(rbind, diagnostic_rows)
  if (anyDuplicated(bank$diagnostics[c("beta", "m")])) stop("Duplicate generated coefficients.")
  if (any(!is.finite(as.matrix(bank$diagnostics))) ||
      max(bank$diagnostics$deviation) > tolerance*2 ||
      max(bank$diagnostics$equioscillation_gap) > tolerance*2 ||
      max(bank$diagnostics$mp_conversion_gap) > 10^(-digits/2) ||
      max(bank$diagnostics$double_conversion_gap) > 1e-12)
    stop("Generated coefficient diagnostics failed validation.")
  bank
}

fixed_int_coefficients <- function(bank, m, beta) {
  index <- which(vapply(bank$coefficients, function(cf)
    cf$m == m && abs(cf$beta-beta) <= 32*.Machine$double.eps*max(1,abs(beta)), logical(1)))
  if (length(index) != 1L)
    stop("No unique freshly computed coefficients for m=", m, ", beta=", beta,
         ". Regenerate the coefficient bank; interpolation is not used.")
  cf <- bank$coefficients[[index]]
  fixed_int_validate_coefficients(cf)
  if (cf$x_min != bank$interval$x_min) stop("The coefficient interval has changed; regenerate the bank.")
  cf
}
