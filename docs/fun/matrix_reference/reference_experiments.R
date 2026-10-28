## -----------------------------------------------------------------------------
run_reference_cases <- function(cases, FUN, workers = 1L) {
  stopifnot(is.function(FUN), length(workers) == 1L, is.finite(workers),
            workers >= 1L, workers == as.integer(workers))
  if (!length(cases)) return(list())
  workers <- min(as.integer(workers), length(cases))
  task <- function(x) tryCatch(list(value = FUN(x)),
                              error = function(e) list(error = conditionMessage(e)))
  if (workers > 1L) {
    if (.Platform$OS.type != "unix")
      stop("Parallel reference cases require Unix fork workers; set workers = 1 on this platform.")
    results <- parallel::mclapply(cases, task, mc.cores = workers, mc.preschedule = FALSE)
  } else results <- lapply(cases, task)
  failed <- which(vapply(results, function(x)
    inherits(x, "try-error") || !is.list(x) || !is.null(x$error) ||
      !("value" %in% names(x)), logical(1)))
  if (length(failed)) {
    i <- failed[1L]
    detail <- if (is.list(results[[i]])) results[[i]]$error else as.character(results[[i]])
    stop("Reference case ", i, " failed: ", paste(detail, collapse = " "), call. = FALSE)
  }
  lapply(results, `[[`, "value")
}

reference_experiment_run <- function(study, parameters, run_id = NULL) {
  stopifnot(length(study) == 1L, study %in% c("h", "m", "tau"), is.list(parameters))
  project <- Sys.getenv("NAFDEMG_PRIVATE_DIR", unset = "")
  if (!nzchar(project)) project <- path.expand("~/Desktop/DOCS/NAFDEMG")
  parent <- file.path(project, "diffusion_experiments", study)
  if (!dir.exists(parent) && !dir.create(parent, recursive = TRUE))
    stop("Cannot create private experiment directory: ", parent)
  if (is.null(run_id)) {
    run_dir <- tempfile(pattern = paste0(format(Sys.time(), "%Y%m%dT%H%M%S", tz = "UTC"), "-"),
                        tmpdir = parent)
  } else {
    stopifnot(length(run_id) == 1L, !is.na(run_id),
              grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", run_id))
    run_dir <- file.path(parent, run_id)
  }
  metadata <- file.path(run_dir, "parameters.rds")
  if (dir.exists(run_dir)) {
    if (!file.exists(metadata) || !identical(readRDS(metadata), parameters))
      stop("Existing run parameters do not match: ", run_dir)
  } else {
    if (!dir.create(run_dir)) stop("Cannot create run directory: ", run_dir)
    saveRDS(parameters, metadata)
  }
  normalizePath(run_dir, mustWork = TRUE)
}
