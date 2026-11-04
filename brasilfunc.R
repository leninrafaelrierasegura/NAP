## ----brasil-functions-0, eval=TRUE--------------------------------------------
brasil_supported_orders <- function() {
    tables <- grep("^m_brasil[0-9]+t$", ls(asNamespace("rSPDE")), value = TRUE)
    sort(as.integer(sub("^m_brasil([0-9]+)t$", "\\1", tables)))
}

brasil_coefficients <- function(m, beta) {
    stopifnot(length(m) == 1L, is.finite(m), m %in% brasil_supported_orders(), length(beta) == 1L, is.finite(beta),
        beta > 0, beta < 1)
    table <- get(paste0("m_brasil", m, "t"), envir = asNamespace("rSPDE"))
    tolerance <- 32 * .Machine$double.eps * max(1, abs(beta))
    index <- which(abs(table$alpha - beta) <= tolerance)
    value <- function(column) {
        if (length(index) == 1L)
            table[[column]][index]
        else stats::approx(table$alpha, table[[column]], xout = beta)$y
    }
    poles <- vapply(paste0("p", seq_len(m)), value, numeric(1))
    residues <- vapply(paste0("r", seq_len(m)), value, numeric(1))
    constant <- value("k")
    if (any(!is.finite(c(poles, residues, constant))) || any(poles >= 0) || any(residues <= 0) || constant < 0 ||
        anyDuplicated(poles))
        stop("BRASIL table at m = ", m, ", beta = ", beta, " requires distinct negative poles, positive residues, and a nonnegative constant.")
    index <- order(poles)
    list(poles = unname(poles[index]), residues = unname(residues[index]), constant = unname(constant), m = as.integer(m),
        beta = beta)
}

brasil_inverse <- function(y, coefficients) {
    stopifnot(all(is.finite(y)), all(y > 0))
    coefficients$constant + as.numeric((1/outer(y, coefficients$poles, "-")) %*% coefficients$residues)
}

brasil_power <- function(x, coefficients) {
    stopifnot(all(is.finite(x)), all(x >= 0))
    coefficients$constant + as.numeric((x/(1 - outer(x, coefficients$poles))) %*% coefficients$residues)
}

brasil_partial_fractions <- function(coefficients, time_step, scaling) {
    strength <- time_step * scaling
    stopifnot(length(strength) == 1L, is.finite(strength), strength > 0)
    p <- coefficients$poles
    r <- coefficients$residues
    constant <- coefficients$constant
    poles <- numeric(length(p))
    denominator <- function(y) constant + strength + sum(r/(y - p))
    for (i in seq_along(p)) {
        upper <- p[i]
        lower <- if (i == 1L)
            2 * upper
        else p[i - 1L]
        if (i == 1L)
            while (denominator(lower) <= 0) lower <- 2 * lower
        width <- upper - lower
        normalized <- function(t) {
            if (t == 0 && i > 1L)
                return(.Machine$double.xmax)
            if (t == 1)
                return(-.Machine$double.xmax)
            denominator(lower + width * t)
        }
        root <- stats::uniroot(normalized, c(0, 1), tol = 2 * .Machine$double.eps, maxiter = 1000L)$root
        poles[i] <- lower + width * root
    }
    residues <- vapply(poles, function(z) strength/sum(r/(z - p)^2), numeric(1))
    remainder <- constant/(constant + strength)
    stopifnot(all(is.finite(c(poles, residues, remainder))), all(poles < 0), all(residues > 0), remainder >= 0,
        remainder < 1)
    y <- 10^seq(0, 12, length.out = 257L)
    inverse <- brasil_inverse(y, coefficients)
    expected <- inverse/(inverse + strength)
    actual <- remainder + as.numeric((1/outer(y, poles, "-")) %*% residues)
    reconstruction_error <- max(abs(actual - expected))
    if (reconstruction_error > 2e-12)
        stop("Inaccurate BRASIL time-step partial fractions: ", reconstruction_error)
    list(p = poles, r = residues, k = remainder, reconstruction_error = reconstruction_error)
}


## ----brasil-functions-1, eval=TRUE--------------------------------------------
my.fractional.operators.frac <- function(L, beta, C, scale.factor, m = 1, time_step, diffusion = 1) {
    stopifnot(length(diffusion) == 1L, is.finite(diffusion), diffusion > 0, length(time_step) == 1L, is.finite(time_step),
        time_step > 0, length(beta) == 1L, is.finite(beta), beta > 0, beta <= 1, length(scale.factor) == 1L, is.finite(scale.factor),
        scale.factor > 0)
    L <- L/scale.factor
    if (beta == 1) {
        L <- L * scale.factor^beta
        lhs <- C + time_step * diffusion * L
        return(list(C = C, L = L, m = m, beta = beta, diffusion = diffusion, LHS = lhs, factorization = Matrix::Cholesky(Matrix::forceSymmetric(lhs),
            LDL = FALSE)))
    }
    scaling <- diffusion * scale.factor^beta
    coefficients <- brasil_coefficients(m, beta)
    poles_rs_k <- brasil_partial_fractions(coefficients, time_step, scaling)
    partial_fraction_terms <- lapply(poles_rs_k$p, function(p) L - p * C)
    remainder <- if (length(poles_rs_k$k))
        poles_rs_k$k
    else 0
    list(C = C, L = L, m = m, beta = beta, diffusion = diffusion, partial_fraction_terms = partial_fraction_terms,
        residues = poles_rs_k$r, poles = poles_rs_k$p, remainder = remainder, mass_factor = if (remainder != 0) Matrix::Cholesky(Matrix::forceSymmetric(C),
            LDL = FALSE) else NULL, factorizations = lapply(partial_fraction_terms, function(A) Matrix::Cholesky(Matrix::forceSymmetric(A),
            LDL = FALSE)))
}

my.solver.frac <- function(obj, v) {
    if (obj$beta == 1) {
        return(if (!is.null(obj$factorization)) Matrix::solve(obj$factorization, v, system = "A") else Matrix::solve(obj$LHS, v))
    }
    partial_fraction_terms <- obj$partial_fraction_terms
    residues <- obj$residues
    output <- v * 0
    if (!is.null(obj$remainder) && obj$remainder != 0) {
        mass_solved <- if (!is.null(obj$mass_factor))
            Matrix::solve(obj$mass_factor, v, system = "A")
        else Matrix::solve(obj$C, v)
        output <- obj$remainder * mass_solved
    }
    for (i in seq_along(residues)) {
        solved <- if (!is.null(obj$factorizations)) {
            Matrix::solve(obj$factorizations[[i]], v, system = "A")
        }
        else Matrix::solve(partial_fraction_terms[[i]], v)
        output <- output + residues[i] * solved
    }
    return(output)
}

prepare_reference_error <- function(Psi, reference_C, reference_basis) {
    E <- as.matrix(reference_basis)
    stopifnot(nrow(Psi) == nrow(reference_C), nrow(reference_C) == ncol(reference_C), nrow(E) == nrow(Psi), ncol(E) >=
        1L, all(is.finite(E)))
    metric <- Matrix::forceSymmetric(as(Matrix::Matrix(Matrix::crossprod(Psi, reference_C %*% Psi), sparse = TRUE),
        "CsparseMatrix"))
    cross <- Matrix::crossprod(Psi, reference_C %*% E)
    factor <- Matrix::Cholesky(metric, LDL = FALSE)
    projection <- as.matrix(Matrix::solve(factor, cross, system = "A"))
    residual <- E - as.matrix(Psi %*% projection)
    list(metric = metric, projection = projection, residual_gram = as.matrix(Matrix::crossprod(residual, reference_C %*%
        residual)), exact_gram = as.matrix(Matrix::crossprod(E, reference_C %*% E)))
}

reference_coefficient_moments <- function(coefficients, weights) {
    coefficients <- as.matrix(coefficients)
    stopifnot(nrow(coefficients) >= 1L, ncol(coefficients) == length(weights), length(weights) >= 1L, all(is.finite(coefficients)),
        all(is.finite(weights)), all(weights > 0))
    weight <- sum(weights)
    centered <- sweep(coefficients, 1L, coefficients[, 1L], `-`)
    shift <- as.numeric(centered %*% (weights/weight))
    deviations <- sweep(centered, 1L, shift, `-`)
    list(weight = weight, mean = coefficients[, 1L] + shift, scatter = tcrossprod(sweep(deviations, 2L, sqrt(weights),
        `*`)))
}

reference_interval_error <- function(u, moments, projection) {
    delta <- as.numeric(u - projection$projection %*% moments$mean)
    metric_product <- as.numeric(projection$metric %*% delta)
    residual_product <- as.numeric(projection$residual_gram %*% moments$mean)
    variance_terms <- projection$exact_gram * moments$scatter
    value <- moments$weight * (sum(delta * metric_product) + sum(moments$mean * residual_product)) + sum(variance_terms)
    scale <- moments$weight * (sum(abs(delta * metric_product)) + sum(abs(moments$mean * residual_product))) +
        sum(abs(variance_terms))
    if (!is.finite(value) || value < -128 * .Machine$double.eps * scale)
        stop("Invalid reference-grid squared error.")
    max(0, value)
}

# R reference retained for numerical regression checks.
solve_fractional_evolution_error_R <- function(my_op_frac, time_step, time_seq, val_at_0, spatial_load, g_average,
    reference_times, coefficients_at_time, error_projection, verbose = FALSE) {
    stopifnot(length(time_step) == 1L, is.finite(time_step), time_step > 0,
        length(time_seq) >= 2L, length(reference_times) >= 2L,
        all(is.finite(time_seq)), all(is.finite(reference_times)),
        all(diff(time_seq) > 0), all(diff(reference_times) > 0),
        reference_times[1] == time_seq[1], tail(reference_times, 1) == tail(time_seq, 1),
        length(g_average) == length(time_seq), length(val_at_0) == nrow(my_op_frac$C),
        length(spatial_load) == length(val_at_0), is.function(coefficients_at_time),
        is.logical(verbose), length(verbose) == 1L, !is.na(verbose),
        nrow(error_projection$metric) == length(val_at_0),
        ncol(error_projection$metric) == length(val_at_0),
        nrow(error_projection$projection) == length(val_at_0))
    n_modes <- ncol(error_projection$projection)
    stopifnot(identical(dim(error_projection$exact_gram), c(n_modes, n_modes)),
        identical(dim(error_projection$residual_gram), c(n_modes, n_modes)))
    times <- reference_times[-1L]
    weights <- diff(reference_times)
    coefficients <- as.matrix(coefficients_at_time(times))
    stopifnot(nrow(coefficients) == n_modes, ncol(coefficients) == length(times),
        all(is.finite(coefficients)))
    runs <- rle(findInterval(times, time_seq, left.open = TRUE))
    u <- as.numeric(val_at_0)
    spatial_load <- as.numeric(spatial_load)
    solver_step <- 0L
    first_column <- 1L
    total <- compensation <- 0
    started <- proc.time()[["elapsed"]]
    for (run in seq_along(runs$values)) {
        target <- runs$values[run]
        while (solver_step < target) {
            solver_step <- solver_step + 1L
            rhs <- my_op_frac$C %*% u + time_step * spatial_load * g_average[solver_step + 1L]
            u <- as.numeric(my.solver.frac(my_op_frac, rhs))
        }
        last_column <- first_column + runs$lengths[run] - 1L
        columns <- seq.int(first_column, last_column)
        moments <- reference_coefficient_moments(coefficients[, columns, drop = FALSE], weights[columns])
        value <- reference_interval_error(u, moments, error_projection)
        adjusted <- value - compensation
        updated <- total + adjusted
        compensation <- (updated - total) - adjusted
        total <- updated
        first_column <- last_column + 1L
    }
    stopifnot(first_column == length(times) + 1L, solver_step == length(time_seq) - 1L)
    if (verbose)
        message(sprintf("Reference error: %s time samples, %s solver steps, %.2f s elapsed",
            length(times), solver_step, proc.time()[["elapsed"]] - started))
    sqrt(total)
}

# Compile once in the parent session; forked workers inherit the loaded kernel.
smooth_compile_solver <- function(cache_dir=tempdir()) {
    target <- environment(smooth_compile_solver)
    if (!exists("smooth_reference_evolve_cpp", envir=target, inherits=FALSE))
        Rcpp::sourceCpp(code=r"----(
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppEigen.h>
#include <memory>
#include <vector>
#include <cmath>
#include <limits>
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::SparseMatrix;
typedef Eigen::SimplicialLLT<SparseMatrix<double> > SmoothFactor;

// Backward Euler/BRASIL updates and the same centered sampled-error identity
// as the R reference implementation. Every solver step is advanced.
// [[Rcpp::export]]
Rcpp::List smooth_reference_evolve_cpp(
    const Eigen::MappedSparseMatrix<double>& C_in,
    Rcpp::List matrices, VectorXd residues, double remainder,
    VectorXd initial, VectorXd load, VectorXd forcing,
    MatrixXd projection,
    const Eigen::MappedSparseMatrix<double>& metric_in,
    MatrixXd residual_gram, MatrixXd exact_gram, MatrixXd coefficients,
    Rcpp::IntegerVector target_steps, VectorXd weights,
    double tau, int n_steps) {
    Eigen::setNbThreads(1);
    SparseMatrix<double> C=C_in, metric=metric_in;
    const int n=C.rows(), modes=projection.cols();
    const int samples=target_steps.size(), m=residues.size();
    if (n_steps<1 || !std::isfinite(tau) || tau<=0 || samples<1 ||
        C.cols()!=n || matrices.size()!=m || initial.size()!=n || load.size()!=n ||
        forcing.size()!=n_steps+1 || projection.rows()!=n || modes<1 ||
        coefficients.rows()!=modes || coefficients.cols()!=samples ||
        weights.size()!=samples || metric.rows()!=n || metric.cols()!=n ||
        residual_gram.rows()!=modes || residual_gram.cols()!=modes ||
        exact_gram.rows()!=modes || exact_gram.cols()!=modes)
        Rcpp::stop("Invalid smooth-reference dimensions or time grid.");
    SmoothFactor mass_factor;
    if (remainder!=0) {
        mass_factor.compute(C);
        if (mass_factor.info()!=Eigen::Success)
            Rcpp::stop("Mass factorization failed.");
    }
    std::vector<std::unique_ptr<SmoothFactor> > factors;
    for (int i=0;i<m;++i) {
        SparseMatrix<double> A=Rcpp::as<SparseMatrix<double> >(matrices[i]);
        if (A.rows()!=n || A.cols()!=n)
            Rcpp::stop("Invalid shifted matrix dimensions.");
        factors.emplace_back(new SmoothFactor());
        factors.back()->compute(A);
        if (factors.back()->info()!=Eigen::Success)
            Rcpp::stop("Shifted factorization failed.");
    }
    VectorXd u=initial, rhs(n), next(n), delta(n), metric_product(n);
    VectorXd anchor(modes), shift(modes), mean(modes), centered(modes);
    VectorXd deviation(modes), residual_product(modes);
    MatrixXd scatter(modes,modes);
    int step=0, first=0;
    double total=0, compensation=0;
    while (first<samples) {
        Rcpp::checkUserInterrupt();
        const int goal=target_steps[first];
        if (goal<1 || goal<step || goal>n_steps)
            Rcpp::stop("Invalid reference time-to-step mapping.");
        while (step<goal) {
            if ((step & 4095)==0) Rcpp::checkUserInterrupt();
            ++step;
            rhs.noalias()=C*u;
            rhs.noalias()+=tau*load*forcing[step];
            if (remainder!=0) next=remainder*mass_factor.solve(rhs);
            else next.setZero();
            for (int i=0;i<m;++i)
                next.noalias()+=residues[i]*factors[i]->solve(rhs);
            u.swap(next);
        }
        int last=first+1;
        while (last<samples && target_steps[last]==goal) ++last;
        long double weight_sum=0;
        for (int q=first;q<last;++q) {
            if (!std::isfinite(weights[q]) || weights[q]<=0)
                Rcpp::stop("Reference weights must be positive and finite.");
            weight_sum+=weights[q];
        }
        const double weight=static_cast<double>(weight_sum);
        anchor=coefficients.col(first);
        shift.setZero();
        for (int q=first;q<last;++q) {
            if ((q & 4095)==0) Rcpp::checkUserInterrupt();
            centered=coefficients.col(q)-anchor;
            shift.noalias()+=(weights[q]/weight)*centered;
        }
        mean=anchor+shift;
        scatter.setZero();
        for (int q=first;q<last;++q) {
            centered=coefficients.col(q)-anchor;
            deviation=std::sqrt(weights[q])*(centered-shift);
            scatter.noalias()+=deviation*deviation.transpose();
        }
        delta.noalias()=u-projection*mean;
        metric_product.noalias()=metric*delta;
        residual_product.noalias()=residual_gram*mean;
        const double variance=exact_gram.cwiseProduct(scatter).sum();
        double value=weight*(delta.dot(metric_product)+mean.dot(residual_product))+variance;
        const double scale=weight*(delta.cwiseProduct(metric_product).cwiseAbs().sum()+
            mean.cwiseProduct(residual_product).cwiseAbs().sum())+
            exact_gram.cwiseProduct(scatter).cwiseAbs().sum();
        if (!std::isfinite(value) || value < -128*std::numeric_limits<double>::epsilon()*scale)
            Rcpp::stop("Invalid reference-grid squared error.");
        value=std::max(0.0,value);
        const double adjusted=value-compensation;
        const double updated=total+adjusted;
        compensation=(updated-total)-adjusted;
        total=updated;
        first=last;
    }
    if (step!=n_steps) Rcpp::stop("Reference grid did not reach the final solver step.");
    return Rcpp::List::create(Rcpp::Named("error")=std::sqrt(total),
        Rcpp::Named("final_state")=u, Rcpp::Named("solver_steps")=step);
}
)----", env=target, cacheDir=cache_dir, showOutput=FALSE, verbose=FALSE)
    invisible(TRUE)
}

solve_fractional_evolution_error <- function(my_op_frac, time_step, time_seq, val_at_0, spatial_load, g_average,
    reference_times, coefficients_at_time, error_projection, verbose = FALSE) {
    stopifnot(length(time_step) == 1L, is.finite(time_step), time_step > 0,
        length(time_seq) >= 2L, length(reference_times) >= 2L,
        all(is.finite(time_seq)), all(is.finite(reference_times)),
        all(diff(time_seq) > 0), all(diff(reference_times) > 0),
        reference_times[1] == time_seq[1], tail(reference_times, 1) == tail(time_seq, 1),
        length(g_average) == length(time_seq), length(val_at_0) == nrow(my_op_frac$C),
        length(spatial_load) == length(val_at_0), is.function(coefficients_at_time),
        is.logical(verbose), length(verbose) == 1L, !is.na(verbose),
        nrow(error_projection$metric) == length(val_at_0),
        ncol(error_projection$metric) == length(val_at_0),
        nrow(error_projection$projection) == length(val_at_0))
    n_modes <- ncol(error_projection$projection)
    stopifnot(identical(dim(error_projection$exact_gram), c(n_modes, n_modes)),
        identical(dim(error_projection$residual_gram), c(n_modes, n_modes)))
    times <- reference_times[-1L]
    coefficients <- as.matrix(coefficients_at_time(times))
    stopifnot(nrow(coefficients) == n_modes, ncol(coefficients) == length(times),
        all(is.finite(coefficients)))
    smooth_compile_solver()
    if (my_op_frac$beta == 1) {
        matrices <- list(my_op_frac$LHS)
        residues <- 1
        remainder <- 0
    } else {
        matrices <- my_op_frac$partial_fraction_terms
        residues <- my_op_frac$residues
        remainder <- my_op_frac$remainder
    }
    general <- function(A) as(as(A, "generalMatrix"), "CsparseMatrix")
    started <- proc.time()[["elapsed"]]
    result <- smooth_reference_evolve_cpp(general(my_op_frac$C), lapply(matrices, general),
        residues, remainder, as.numeric(val_at_0), as.numeric(spatial_load), as.numeric(g_average),
        error_projection$projection, general(error_projection$metric),
        error_projection$residual_gram, error_projection$exact_gram, coefficients,
        as.integer(findInterval(times, time_seq, left.open=TRUE)), diff(reference_times),
        time_step, length(time_seq)-1L)
    if (verbose)
        message(sprintf("Reference error (C++): %s time samples, %s solver steps, %.2f s elapsed",
            length(times), result$solver_steps, proc.time()[["elapsed"]]-started))
    result$error
}

make_time_grid <- function(T_final, target_step) {
    stopifnot(length(T_final) == 1L, is.finite(T_final), T_final > 0, length(target_step) == 1L, is.finite(target_step),
        target_step > 0)
    N <- ceiling(T_final/target_step)
    time_step <- T_final/N
    time_seq <- seq(0, T_final, length.out = N + 1)
    return(list(N = N, time_step = time_step, time_seq = time_seq))
}

average_sine_forcing <- function(time_seq, A = 1, omega = pi) {
    dt <- diff(time_seq)
    stopifnot(length(time_seq) >= 2L, all(is.finite(time_seq)), all(dt > 0))
    midpoint <- head(time_seq, -1) + dt/2
    z <- omega * dt/2
    sinc_z <- rep(1, length(z))
    nonzero <- z != 0
    sinc_z[nonzero] <- sin(z[nonzero])/z[nonzero]
    averages <- A * sin(omega * midpoint) * sinc_z
    c(0, averages)
}


## ----brasil-functions-2, eval=TRUE--------------------------------------------
gets.graph.tadpole <- function(h) {
  stopifnot(length(h)==1L, is.finite(h), h>0)
  if ("edges" %in% names(formals(MetricGraph::metric_graph$public_methods$initialize))) {
    theta <- seq(-pi,pi,length.out=10000L)
    edges <- list(rbind(c(0,0),c(1,0)), cbind(1+1/pi+cos(theta)/pi,sin(theta)/pi))
    graph <- MetricGraph::metric_graph$new(edges=edges,verbose=0)
    graph$set_manual_edge_lengths(edge_lengths=c(1,2))
    graph$build_mesh(h=h)
    return(graph)
  }
  counts <- as.integer(ceiling(c(1,2)/h))
  n <- as.integer(2+sum(counts-1L))
  inside <- list(if(counts[1]>1L) seq.int(3L,counts[1]+1L) else integer(),
                 if(counts[2]>1L) seq.int(counts[1]+2L,n) else integer())
  paths <- list(c(1L,inside[[1]],2L), c(2L,inside[[2]],2L))
  E <- do.call(rbind,lapply(paths,function(v)cbind(head(v,-1L),tail(v,-1L))))
  lengths <- rep(c(1,2)/counts,counts)
  fem <- get("assemble_fem",asNamespace("MetricGraph"))(E,lengths,n)
  PtE <- do.call(rbind,lapply(1:2,function(e)cbind(e,seq_len(counts[e]-1L)/counts[e])))
  locations <- rbind(c(1,0),c(2,0),PtE)
  basis <- function(points) {
    edge <- as.integer(points[,1]); t <- points[,2]*counts[edge]
    cell <- pmin(floor(t),counts[edge]-1L); weight <- t-cell
    left <- right <- integer(nrow(points))
    for(e in 1:2) {
      ids <- which(edge==e)
      left[ids] <- paths[[e]][cell[ids]+1L]
      right[ids] <- paths[[e]][cell[ids]+2L]
    }
    Matrix::sparseMatrix(i=rep(seq_len(nrow(points)),2L),j=c(left,right),
      x=c(1-weight,weight),dims=c(nrow(points),n))
  }
  list(mesh=c(fem,list(PtE=PtE,h_e=lengths)),compute_fem=function()invisible(NULL),
       get_edge_lengths=function()c(1,2),get_mesh_locations=function()locations,fem_basis=basis)
}

tadpole.eig <- function(k, graph) {
    x1 <- c(0, graph$get_edge_lengths()[1] * graph$mesh$PtE[graph$mesh$PtE[, 1] == 1, 2])
    x2 <- c(0, graph$get_edge_lengths()[2] * graph$mesh$PtE[graph$mesh$PtE[, 1] == 2, 2])
    if (k == 0) {
        f.e1 <- rep(1, length(x1))
        f.e2 <- rep(1, length(x2))
        f1 = c(f.e1[1], f.e2[1], f.e1[-1], f.e2[-1])
        f = list(phi = f1/sqrt(3))
    }
    else {
        f.e1 <- -2 * sin(pi * k * 1/2) * cos(pi * k * x1/2)
        f.e2 <- sin(pi * k * x2/2)
        f1 = c(f.e1[1], f.e2[1], f.e1[-1], f.e2[-1])
        if ((k%%2) == 1) {
            f = list(phi = f1/sqrt(3))
        }
        else {
            f.e1 <- (-1)^{
                k/2
            } * cos(pi * k * x1/2)
            f.e2 <- cos(pi * k * x2/2)
            f2 = c(f.e1[1], f.e2[1], f.e1[-1], f.e2[-1])
            f <- list(phi = f1, psi = f2/sqrt(3/2))
        }
    }
    return(f)
}

gets.eigen.params <- function(N_finite = 4, kappa = 1, alpha = 0.5, graph) {
    EIGENVAL <- NULL
    EIGENVAL_ALPHA <- NULL
    EIGENVAL_MINUS_ALPHA <- NULL
    EIGENFUN <- NULL
    INDEX <- NULL
    for (j in 0:N_finite) {
        lambda_j <- kappa^2 + (j * pi/2)^2
        lambda_j_alpha_half <- lambda_j^(alpha/2)
        lambda_j_minus_alpha_half <- lambda_j^(-alpha/2)
        e_j <- tadpole.eig(j, graph)$phi
        EIGENVAL <- c(EIGENVAL, lambda_j)
        EIGENVAL_ALPHA <- c(EIGENVAL_ALPHA, lambda_j_alpha_half)
        EIGENVAL_MINUS_ALPHA <- c(EIGENVAL_MINUS_ALPHA, lambda_j_minus_alpha_half)
        EIGENFUN <- cbind(EIGENFUN, e_j)
        INDEX <- c(INDEX, j)
        if (j > 0 && (j%%2 == 0)) {
            lambda_j <- kappa^2 + (j * pi/2)^2
            lambda_j_alpha_half <- lambda_j^(alpha/2)
            lambda_j_minus_alpha_half <- lambda_j^(-alpha/2)
            e_j <- tadpole.eig(j, graph)$psi
            EIGENVAL <- c(EIGENVAL, lambda_j)
            EIGENVAL_ALPHA <- c(EIGENVAL_ALPHA, lambda_j_alpha_half)
            EIGENVAL_MINUS_ALPHA <- c(EIGENVAL_MINUS_ALPHA, lambda_j_minus_alpha_half)
            EIGENFUN <- cbind(EIGENFUN, e_j)
            INDEX <- c(INDEX, j + 0.1)
        }
    }
    return(list(EIGENVAL = EIGENVAL, EIGENVAL_ALPHA = EIGENVAL_ALPHA, EIGENVAL_MINUS_ALPHA = EIGENVAL_MINUS_ALPHA,
        EIGENFUN = EIGENFUN, INDEX = INDEX))
}

trig_hat_integrals <- function(a, b, w, A = 0, B = 0) {
    h <- b - a
    z <- w * h/2
    theta <- w * (a + b)/2
    if (abs(z) < 0.001) {
        S <- 1 - z^2/6 + z^4/120 - z^6/5040
        Q <- z/3 - z^3/30 + z^5/840 - z^7/45360
    }
    else {
        S <- sin(z)/z
        Q <- (sin(z) - z * cos(z))/z^2
    }
    v <- A * cos(theta) + B * sin(theta)
    q <- A * sin(theta) - B * cos(theta)
    h/2 * c(left = S * v + Q * q, right = S * v - Q * q)
}

project_u0_tadpole <- function(graph, amplitude = 10) {
    C <- graph$mesh$C
    PtE <- graph$mesh$PtE
    lengths <- graph$get_edge_lengths()
    stopifnot(!is.null(C), length(lengths) == 2L, nrow(C) == 2L + nrow(PtE), length(amplitude) == 1L, is.finite(amplitude))
    local_load <- function(a, b, w, A, B) {
        h <- b - a
        z <- w * h/2
        theta <- w * (a + b)/2
        if (abs(z) < 0.001) {
            S <- 1 - z^2/6 + z^4/120 - z^6/5040
            Q <- z/3 - z^3/30 + z^5/840 - z^7/45360
        }
        else {
            S <- sin(z)/z
            Q <- (sin(z) - z * cos(z))/z^2
        }
        v <- A * cos(theta) + B * sin(theta)
        q <- A * sin(theta) - B * cos(theta)
        h/2 * c(S * v + Q * q, S * v - Q * q)
    }
    b0 <- numeric(nrow(C))
    w <- 3 * pi/2
    for (edge in 1:2) {
        idx <- which(PtE[, 1] == edge)
        idx <- idx[order(PtE[idx, 2])]
        x <- c(0, lengths[edge] * PtE[idx, 2], lengths[edge])
        start_node <- if (edge == 1L)
            1L
        else 2L
        nodes <- c(start_node, 2L + idx, 2L)
        A <- if (edge == 1L)
            2 * amplitude/sqrt(3)
        else 0
        B <- if (edge == 2L)
            amplitude/sqrt(3)
        else 0
        for (j in seq_len(length(x) - 1L)) {
            load <- local_load(x[j], x[j + 1L], w, A, B)
            b0[nodes[j]] <- b0[nodes[j]] + load[1]
            b0[nodes[j + 1L]] <- b0[nodes[j + 1L]] + load[2]
        }
    }
    as.numeric(Matrix::solve(C, b0))
}

assemble_e7_load_tadpole <- function(graph) {
    C <- graph$mesh$C
    PtE <- graph$mesh$PtE
    lengths <- graph$get_edge_lengths()
    stopifnot(!is.null(C), length(lengths) == 2L, nrow(C) == 2L + nrow(PtE))
    b7 <- numeric(nrow(C))
    for (edge in 1:2) {
        idx <- which(PtE[, 1] == edge)
        idx <- idx[order(PtE[idx, 2])]
        x <- c(0, lengths[edge] * PtE[idx, 2], lengths[edge])
        start_node <- if (edge == 1L)
            1L
        else 2L
        nodes <- c(start_node, 2L + idx, 2L)
        for (j in seq_len(length(x) - 1L)) {
            load <- trig_hat_integrals(x[j], x[j + 1L], w = 2 * pi, A = sqrt(2/3), B = 0)
            b7[nodes[j]] <- b7[nodes[j]] + load[1]
            b7[nodes[j + 1L]] <- b7[nodes[j + 1L]] + load[2]
        }
    }
    b7
}

decayed_sine_convolution <- function(t, rate, A = 1, omega = pi) {
    stopifnot(all(is.finite(t)), all(t >= 0), length(rate) == 1L, is.finite(rate), rate >= 0, is.finite(A), is.finite(omega))
    z <- (rate + (0+1i) * omega) * t
    small <- Mod(z) < 0.5
    answer <- numeric(length(t))
    if (any(small)) {
        zz <- z[small]
        term <- rep(1 + (0+0i), length(zz))
        series <- rep(0 + (0+0i), length(zz))
        for (k in 1:24) {
            term <- term * zz/(k + 1)
            series <- series + term
        }
        answer[small] <- A * exp(-rate * t[small]) * t[small] * Im(series)
    }
    if (any(!small)) {
        tt <- t[!small]
        answer[!small] <- A * Im((exp((0+1i) * omega * tt) - exp(-rate * tt))/(rate + (0+1i) * omega))
    }
    answer
}

exact_modal_coefficients <- function(t, rates, initial, forcing, A = 1, omega = pi) {
    stopifnot(length(rates) == length(initial), length(rates) == length(forcing), all(is.finite(rates)), all(rates >
        0))
    out <- matrix(0, nrow = length(rates), ncol = length(t))
    for (j in seq_along(rates)) {
        out[j, ] <- initial[j] * exp(-rates[j] * t) + forcing[j] * decayed_sine_convolution(t, rates[j], A = A,
            omega = omega)
    }
    out
}

run_reference_cases <- function(cases, FUN, workers = 1L) {
    stopifnot(is.function(FUN), length(workers) == 1L, is.finite(workers), workers >= 1L, workers == as.integer(workers))
    if (!length(cases))
        return(list())
    workers <- min(as.integer(workers), length(cases))
    task <- function(x) tryCatch(list(value = FUN(x)), error = function(e) list(error = conditionMessage(e)))
    if (workers > 1L) {
        if (.Platform$OS.type != "unix")
            stop("Parallel reference cases require Unix fork workers; set workers = 1 on this platform.")
        results <- parallel::mclapply(cases, task, mc.cores = workers, mc.preschedule = FALSE)
    }
    else results <- lapply(cases, task)
    failed <- which(vapply(results, function(x) inherits(x, "try-error") || !is.list(x) || !is.null(x$error) ||
        !("value" %in% names(x)), logical(1)))
    if (length(failed)) {
        i <- failed[1L]
        detail <- if (is.list(results[[i]]))
            results[[i]]$error
        else as.character(results[[i]])
        stop("Reference case ", i, " failed: ", paste(detail, collapse = " "), call. = FALSE)
    }
    lapply(results, `[[`, "value")
}


## ----brasil-functions-3, eval=TRUE--------------------------------------------


isolated_discrete_spectrum <- function(C, L, max_dimension = 1200L, max_dense_mib = 1024, residual_tolerance = 1e-08,
    orthogonality_tolerance = 1e-10) {
    n <- nrow(C)
    stopifnot(n >= 1L, ncol(C) == n, identical(dim(L), dim(C)), Matrix::isSymmetric(C, tol = 1e-12), Matrix::isSymmetric(L, tol = 1e-12),
        length(max_dimension) == 1L, is.finite(max_dimension), max_dimension >= 1, length(max_dense_mib) == 1L,
        is.finite(max_dense_mib), max_dense_mib > 0, is.finite(residual_tolerance), residual_tolerance > 0, is.finite(orthogonality_tolerance),
        orthogonality_tolerance > 0)
    estimated_mib <- 20 * 8 * as.double(n)^2/1024^2
    if (n > max_dimension || estimated_mib > max_dense_mib)
        stop("Complete eigendecomposition refused before dense allocation: dimension ", n, ", estimated temporary storage ",
            round(estimated_mib, 1), " MiB. Choose a moderate isolated-study mesh or explicitly raise the guards.")
    started <- proc.time()[["elapsed"]]
    dense_C <- as.matrix(C)
    dense_L <- as.matrix(L)
    stopifnot(all(is.finite(dense_C)), all(is.finite(dense_L)))
    R <- chol(dense_C)
    Rinv <- backsolve(R, diag(n))
    symmetric <- Matrix::crossprod(Rinv, dense_L %*% Rinv)
    ee <- eigen((symmetric + t(symmetric))/2, symmetric = TRUE)
    stopifnot(all(is.finite(ee$values)), all(ee$values > 0))
    V <- backsolve(R, ee$vectors)
    CV <- C %*% V
    LV <- L %*% V
    residual <- as.matrix(LV - sweep(as.matrix(CV), 2L, ee$values, `*`))
    denominator <- sqrt(colSums(as.matrix(LV)^2)) + abs(ee$values) * sqrt(colSums(as.matrix(CV)^2))
    residual_error <- max(sqrt(colSums(residual^2))/pmax(denominator, .Machine$double.xmin))
    orthogonality_error <- max(abs(Matrix::crossprod(V, CV) - diag(n)))
    if (residual_error > residual_tolerance || orthogonality_error > orthogonality_tolerance)
        stop("Numerical eigensystem failed validation: residual = ", residual_error, ", C-orthogonality = ", orthogonality_error)
    list(values = ee$values, vectors = V, diagnostics = data.frame(dimension = n, eigenvalue_min = min(ee$values),
        eigenvalue_max = max(ee$values), relative_residual = residual_error, C_orthogonality_error = orthogonality_error,
        estimated_dense_mib = estimated_mib, elapsed_seconds = proc.time()[["elapsed"]] - started))
}

isolated_exact_be_reference <- function(C, spectrum, initial, spatial_load, time_seq, alpha, diffusion = 1, A = 1,
    omega = pi, max_history_mib = 256) {
    n <- nrow(C)
    stopifnot(length(initial) == n, length(spatial_load) == n, all(is.finite(initial)), all(is.finite(spatial_load)),
        length(time_seq) >= 2L, all(is.finite(time_seq)), all(diff(time_seq) > 0), length(alpha) == 1L, is.finite(alpha),
        alpha > 0, alpha <= 2, length(diffusion) == 1L, is.finite(diffusion), diffusion > 0, length(A) == 1L, is.finite(A),
        length(omega) == 1L, is.finite(omega), is.finite(max_history_mib), max_history_mib > 0, identical(dim(spectrum$vectors),
            c(n, n)), length(spectrum$values) == n)
    dt <- diff(time_seq)[1L]
    stopifnot(max(abs(diff(time_seq) - dt)) <= 128 * .Machine$double.eps * max(1, abs(tail(time_seq, 1L))))
    history_mib <- 2 * 8 * as.double(n) * length(time_seq)/1024^2
    if (history_mib > max_history_mib)
        stop("Exact reference history exceeds max_history_mib; use a moderate time grid.")
    V <- spectrum$vectors
    x <- as.numeric(Matrix::crossprod(V, C %*% initial))
    f <- as.numeric(Matrix::crossprod(V, spatial_load))
    rates <- diffusion * spectrum$values^(alpha/2)
    mu <- 1/(1 + dt * rates)
    g <- average_sine_forcing(time_seq, A = A, omega = omega)
    states <- coefficients <- matrix(0, n, length(time_seq))
    states[, 1L] <- initial
    coefficients[, 1L] <- x
    y <- x
    for (k in seq_len(length(time_seq) - 1L)) {
        y <- mu * (y + dt * f * g[k + 1L])
        coefficients[, k + 1L] <- y
        states[, k + 1L] <- as.numeric(V %*% y)
    }
    list(states = states, coefficients = coefficients, initial_coefficients = x, forcing_coefficients = f, exact_rates = rates,
        g_average = g, dt = dt, time_seq = time_seq, alpha = alpha, diffusion = diffusion, history_mib = history_mib)
}

isolated_tadpole_data <- function(h, kappa = 4, initial_amplitude = 10, forcing_amplitude = 10) {
    stopifnot(length(h) == 1L, is.finite(h), h > 0, length(kappa) == 1L, is.finite(kappa), kappa > 0, is.finite(initial_amplitude),
        is.finite(forcing_amplitude))
    graph <- gets.graph.tadpole(h)
    graph$compute_fem()
    C <- graph$mesh$C
    L <- kappa^2 * C + graph$mesh$G
    PtE <- graph$mesh$PtE
    lengths <- graph$get_edge_lengths()
    stopifnot(length(lengths) == 2L, max(abs(lengths - c(1, 2))) < 1e-12, nrow(C) == 2L + nrow(PtE))
    load <- numeric(nrow(C))
    element_lengths <- numeric()
    for (edge in 1:2) {
        idx <- which(PtE[, 1L] == edge)
        idx <- idx[order(PtE[idx, 2L])]
        coordinates <- c(0, lengths[edge] * PtE[idx, 2L], lengths[edge])
        nodes <- c(if (edge == 1L) 1L else 2L, 2L + idx, 2L)
        element_lengths <- c(element_lengths, diff(coordinates))
        for (i in seq_len(length(coordinates) - 1L)) {
            pair <- trig_hat_integrals(coordinates[i], coordinates[i + 1L], 2 * pi, A = forcing_amplitude * sqrt(2/3))
            load[nodes[i]] <- load[nodes[i]] + pair[1L]
            load[nodes[i + 1L]] <- load[nodes[i + 1L]] + pair[2L]
        }
    }
    list(C = C, L = L, initial = project_u0_tadpole(graph, initial_amplitude), spatial_load = load, kappa = kappa,
        h_requested = h, h_min = min(element_lengths), h_max = max(element_lengths))
}

isolated_rational_rates <- function(eigenvalues, alpha, m, kappa, diffusion=1, coefficients=NULL) {
  if(alpha==2) return(diffusion*eigenvalues)
  if(is.null(coefficients)) coefficients <- brasil_coefficients(m,alpha/2)
  diffusion*kappa^alpha/brasil_inverse(eigenvalues/kappa^2,coefficients)
}

isolated_beta_bracket <- function(m,beta) {
  table <- get(paste0("m_brasil",m,"t"),asNamespace("rSPDE"))
  grid <- table$alpha; tolerance <- 32*.Machine$double.eps*max(1,abs(beta))
  list(exact=any(abs(grid-beta)<=tolerance),lower=max(grid[grid<=beta+tolerance]),
       upper=min(grid[grid>=beta-tolerance]))
}


isolated_rational_error <- function(data,spectrum,reference,m,excitation_tolerance=1e-8) {
    started <- proc.time()[["elapsed"]]
    stopifnot(nrow(reference$states) == nrow(data$C), ncol(reference$states) == length(reference$time_seq))
    coefficients <- if(reference$alpha==2) NULL else brasil_coefficients(m,reference$alpha/2)
    op <- my.fractional.operators.frac(data$L,reference$alpha/2,data$C,data$kappa^2,
      m,reference$dt,reference$diffusion)
    rates <- isolated_rational_rates(spectrum$values,reference$alpha,m,data$kappa,
      reference$diffusion,coefficients)
    mu <- 1/(1 + reference$dt * rates)
    y <- reference$initial_coefficients
    u <- as.numeric(data$initial)
    total <- spectral_total <- compensation <- spectral_compensation <- 0
    for (k in seq_len(length(reference$time_seq) - 1L)) {
        forcing <- reference$g_average[k + 1L]
        rhs <- data$C %*% u + reference$dt * data$spatial_load * forcing
        u <- as.numeric(my.solver.frac(op, rhs))
        delta <- reference$states[, k + 1L] - u
        value <- reference$dt * sum(delta * (data$C %*% delta))
        adjusted <- value - compensation
        updated <- total + adjusted
        compensation <- (updated - total) - adjusted
        total <- updated
        y <- mu * (y + reference$dt * reference$forcing_coefficients * forcing)
        dy <- reference$coefficients[, k + 1L] - y
        value <- reference$dt * sum(dy^2)
        adjusted <- value - spectral_compensation
        updated <- spectral_total + adjusted
        spectral_compensation <- (updated - spectral_total) - adjusted
        spectral_total <- updated
    }
    sparse_error <- sqrt(max(0, total))
    product_error <- sqrt(max(0, spectral_total))
    gap <- u - as.numeric(spectrum$vectors %*% y)
    final_gap <- sqrt(max(0, sum(gap * (data$C %*% gap))))
    excited <- abs(reference$initial_coefficients) > excitation_tolerance * max(abs(reference$initial_coefficients)) |
        abs(reference$forcing_coefficients) > excitation_tolerance * max(abs(reference$forcing_coefficients))
    relative_rate_error <- rates/reference$exact_rates - 1
    bracket <- if (reference$alpha == 2)
        list(exact = TRUE, lower = 1, upper = 1)
    else isolated_beta_bracket(m,reference$alpha/2)
    data.frame(alpha = reference$alpha, diffusion = reference$diffusion, m = m, error = sparse_error, product_error = product_error,
        absolute_error_disagreement = abs(sparse_error - product_error), relative_error_disagreement = abs(sparse_error -
            product_error)/max(product_error, .Machine$double.xmin), final_state_gap_C = final_gap, max_relative_rate_error_full = max(abs(relative_rate_error)),
        max_relative_rate_error_excited = if (any(excited))
            max(abs(relative_rate_error[excited]))
        else 0, excited_numerical_modes = sum(excited), coefficient_source = "brasil", coefficient_beta_exact = bracket$exact,
        beta_exact_in_package_table = bracket$exact, beta_table_lower = bracket$lower, beta_table_upper = bracket$upper, elapsed_seconds = proc.time()[["elapsed"]] -
            started)
}



## ----brasil-functions-4, eval=TRUE--------------------------------------------
rough_rates <- function(j, alpha, diffusion = 1, kappa = 4) {
    stopifnot(all(is.finite(j)), all(j >= 1), all(j == floor(j)), length(alpha) == 1L, is.finite(alpha), alpha >
        0, alpha < 2, length(diffusion) == 1L, is.finite(diffusion), diffusion > 0, length(kappa) == 1L, is.finite(kappa),
        kappa > 0)
    diffusion * (kappa^2 + pi^2 * 2^(2 * j))^(alpha/2)
}

rough_rational_coefficients <- function(beta, order = 16L) {
    stopifnot(beta > 0, beta < 1, order %in% 1:20)
    table <- get(paste0("m_brasil", order, "t"), envir = asNamespace("rSPDE"))
    row <- which(abs(table$alpha - beta) < 32 * .Machine$double.eps)
    if (length(row) != 1L)
        stop("No exact BRASIL coefficient row at beta = ", beta)
    poles <- as.numeric(table[row, paste0("p", seq_len(order))])
    residues <- as.numeric(table[row, paste0("r", seq_len(order))])
    constant <- table$k[row]
    stopifnot(all(is.finite(c(poles, residues, constant))), all(poles < 0), all(residues > 0), constant >= 0)
    index <- order(poles)
    list(poles = poles[index], residues = residues[index], constant = constant, order = as.integer(order), beta = beta)
}


rough_tau_mesh <- function(J, kappa = 4) {
    stopifnot(length(J) == 1L, is.finite(J), J == as.integer(J), J >= 3L, J <= 22L, length(kappa) == 1L, is.finite(kappa),
        kappa > 0)
    h <- 2^-J
    cells <- as.integer(2^(J - 1L))
    n <- cells - 1L
    C <- Matrix::bandSparse(n, k = c(-1L, 0L, 1L), diagonals = list(rep(2 * h/3, n - 1L), rep(8 * h/3, n), rep(2 *
        h/3, n - 1L)))
    G <- Matrix::bandSparse(n, k = c(-1L, 0L, 1L), diagonals = list(rep(-4/h, n - 1L), rep(8/h, n), rep(-4/h, n -
        1L)))
    j <- seq_len(J - 1L)
    x <- seq_len(n) * h
    theta <- pi * 2^j * h
    B <- outer(x, 2^j, function(x, k) sinpi(x * k)) * rep(4 * h * (sin(theta/2)/(theta/2))^2, each = n)
    list(J = J, h = h, C = C, K = kappa^2 * C + G, B = B, j = j, kappa = kappa, load = as.numeric(B %*% (1/j)),
        x = x)
}




## ----brasil-functions-rough-sampled, eval=TRUE--------------------------------
rough_sampled_reference <- function(J_reference, cutoff, max_nodes = 2^22,
                                    max_basis_mib = 2048) {
    stopifnot(length(J_reference)==1L,is.finite(J_reference),
      J_reference==as.integer(J_reference),J_reference>=4L,J_reference<=23L,
      length(cutoff)==1L,is.finite(cutoff),cutoff==as.integer(cutoff),cutoff>=1L,
      J_reference>=cutoff+2L,max_nodes>0,max_basis_mib>0)
    cells <- 2^(J_reference-1L)
    n <- cells-1L
    if (n>max_nodes || 8*n*cutoff/1024^2>max_basis_mib)
        stop("Reference mesh/basis exceeds its allocation guard.")
    h <- 2^-J_reference
    C <- Matrix::bandSparse(n,k=c(-1L,0L,1L),
      diagonals=list(rep(2*h/3,n-1L),rep(8*h/3,n),rep(2*h/3,n-1L)))
    x <- seq_len(n)*h
    E <- outer(x,2^seq_len(cutoff),function(x,k) sinpi(x*k))
    list(J=as.integer(J_reference),h=h,C=C,x=x,basis=E,cutoff=as.integer(cutoff))
}

rough_sampled_fem_basis <- function(J, reference_x) {
    stopifnot(length(J)==1L,is.finite(J),J==as.integer(J),J>=3L,J<=20L,
      all(is.finite(reference_x)),all(reference_x>0 & reference_x<.5))
    position <- reference_x*2^J
    left <- as.integer(floor(position))
    fraction <- position-left
    n <- 2^(J-1L)-1L
    rows <- seq_along(position)
    keep_left <- left>=1L & left<=n
    keep_right <- left+1L>=1L & left+1L<=n
    Matrix::drop0(Matrix::sparseMatrix(
      i=c(rows[keep_left],rows[keep_right]),
      j=c(left[keep_left],left[keep_right]+1L),
      x=c(1-fraction[keep_left],fraction[keep_right]),
      dims=c(length(position),n)))
}

rough_sampled_coefficients <- function(times,cutoff,alpha,diffusion=1,kappa=4) {
    stopifnot(all(is.finite(times)),all(times>=0),cutoff>=1L,cutoff==as.integer(cutoff))
    j <- seq_len(cutoff)
    rates <- rough_rates(j,alpha,diffusion,kappa)
    -expm1(-outer(rates,times))/(j*rates)
}

rough_compile_sampled_solver <- function(cache_dir=tempdir()) {
    target <- environment(rough_compile_sampled_solver)
    if (!exists("rough_sampled_evolve_cpp",envir=target,inherits=FALSE))
        Rcpp::sourceCpp(code=r"----(
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppEigen.h>
#include <memory>
#include <vector>
#include <cmath>
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::SparseMatrix;
typedef Eigen::SimplicialLDLT<SparseMatrix<double> > SampledFactor;

// Same shifted FEM update and sampled reference quadratic norm as the
// smooth notebook. No analytic time moments or continuum-tail sum enters.
// [[Rcpp::export]]
Rcpp::List rough_sampled_evolve_cpp(
    const Eigen::MappedSparseMatrix<double>& C_in, Rcpp::List matrices,
    Eigen::VectorXd residues, double remainder,
    Eigen::VectorXd load, Eigen::MatrixXd projection,
    const Eigen::MappedSparseMatrix<double>& metric_in,
    Eigen::MatrixXd residual_gram, Eigen::MatrixXd coefficients,
    Rcpp::IntegerVector target_steps, Eigen::VectorXd weights,
    double tau, int n_steps) {
  SparseMatrix<double> C=C_in, metric=metric_in;
  const int n=C.rows(), modes=projection.cols(), m=residues.size();
  const int samples=target_steps.size();
  if (n_steps<1 || tau<=0 || matrices.size()!=m || load.size()!=n ||
      projection.rows()!=n || coefficients.rows()!=modes ||
      coefficients.cols()!=samples || weights.size()!=samples ||
      residual_gram.rows()!=modes || residual_gram.cols()!=modes ||
      metric.rows()!=n || metric.cols()!=n)
    Rcpp::stop("Invalid sampled-reference dimensions.");
  SampledFactor mass_factor;
  mass_factor.compute(C);
  if (mass_factor.info()!=Eigen::Success) Rcpp::stop("Mass factorization failed.");
  std::vector<std::unique_ptr<SampledFactor> > factors;
  for (int i=0;i<m;++i) {
    SparseMatrix<double> A=Rcpp::as<SparseMatrix<double> >(matrices[i]);
    factors.emplace_back(new SampledFactor());
    factors.back()->compute(A);
    if (factors.back()->info()!=Eigen::Success) Rcpp::stop("Shifted factorization failed.");
  }
  VectorXd u=VectorXd::Zero(n),rhs(n),next(n),delta(n);
  int step=0;
  double total=0,compensation=0;
  for(int q=0;q<samples;++q) {
    if ((q & 4095)==0) Rcpp::checkUserInterrupt();
    int goal=target_steps[q];
    if(goal<step || goal>n_steps || goal<1 || !(weights[q]>0))
      Rcpp::stop("Invalid reference time-to-step mapping.");
    while(step<goal) {
      if ((step & 4095)==0) Rcpp::checkUserInterrupt();
      rhs.noalias()=C*u;
      rhs.noalias()+=tau*load;
      next=remainder*mass_factor.solve(rhs);
      for(int i=0;i<m;++i) next.noalias()+=residues[i]*factors[i]->solve(rhs);
      u.swap(next);
      ++step;
    }
    delta.noalias()=u-projection*coefficients.col(q);
    double value=weights[q]*(delta.dot(metric*delta)+
      coefficients.col(q).dot(residual_gram*coefficients.col(q)));
    if(!std::isfinite(value) || value<0) Rcpp::stop("Invalid sampled squared error.");
    double adjusted=value-compensation,updated=total+adjusted;
    compensation=(updated-total)-adjusted;
    total=updated;
  }
  if(step!=n_steps) Rcpp::stop("Reference grid did not reach the final solver step.");
  return Rcpp::List::create(Rcpp::Named("error")=std::sqrt(total),
    Rcpp::Named("sampled_squared")=total,Rcpp::Named("final_state")=u);
}
)----",env=target,cacheDir=cache_dir,showOutput=FALSE,verbose=FALSE)
    invisible(TRUE)
}

rough_sampled_error <- function(mesh,reference,error_projection,alpha,time_step,
    reference_time_step,m=14L,diffusion=1,T_final=1,
    max_solver_steps=1e7,max_reference_samples=2e6,max_coefficient_mib=512) {
    stopifnot(reference$cutoff>=mesh$J-1L,reference$J>mesh$J,
      length(T_final)==1L,is.finite(T_final),T_final>0,
      length(time_step)==1L,is.finite(time_step),time_step>0,
      length(reference_time_step)==1L,is.finite(reference_time_step),reference_time_step>0,
      length(max_solver_steps)==1L,is.finite(max_solver_steps),max_solver_steps>0,
      length(max_reference_samples)==1L,is.finite(max_reference_samples),max_reference_samples>0,
      length(max_coefficient_mib)==1L,is.finite(max_coefficient_mib),max_coefficient_mib>0)
    solver_N <- ceiling(T_final/time_step)
    reference_N <- ceiling(T_final/reference_time_step)
    if (!is.finite(solver_N) || !is.finite(reference_N) ||
        solver_N>max_solver_steps || reference_N>max_reference_samples ||
        8*reference_N*reference$cutoff/1024^2>max_coefficient_mib)
        stop("Time grid/reference coefficient matrix exceeds its allocation guard.")
    grid <- make_time_grid(T_final,time_step)
    reference_grid <- make_time_grid(T_final,reference_time_step)
    coefficients <- rough_sampled_coefficients(reference_grid$time_seq[-1L],
      reference$cutoff,alpha,diffusion,mesh$kappa)
    target_steps <- findInterval(reference_grid$time_seq[-1L],grid$time_seq,left.open=TRUE)
    op <- my.fractional.operators.frac(mesh$K,alpha/2,mesh$C,mesh$kappa^2,
      m,grid$time_step,diffusion)
    general <- function(A) as(as(A,"generalMatrix"),"CsparseMatrix")
    started <- proc.time()[["elapsed"]]
    result <- rough_sampled_evolve_cpp(general(mesh$C),
      lapply(op$partial_fraction_terms,general),op$residues,op$remainder,
      mesh$load,error_projection$projection,general(error_projection$metric),
      error_projection$residual_gram,coefficients,as.integer(target_steps),
      diff(reference_grid$time_seq),grid$time_step,grid$N)
    result$elapsed <- proc.time()[["elapsed"]]-started
    result$tau <- grid$time_step
    result$reference_tau <- reference_grid$time_step
    result$N <- grid$N
    result$reference_N <- reference_grid$N
    result
}


