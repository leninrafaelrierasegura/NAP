# Run from the project root. Checks use small matrices and temporary outputs.
f <- tempfile(fileext='.R')
knitr::purl('brasilfunc.Rmd',output=f,quiet=TRUE)
e <- new.env(parent=globalenv()); sys.source(f,e); unlink(f)
stopifnot(all(vapply(as.list(e),is.function,logical(1))))

# Independently evaluate the supplied source coefficients in a complete FEM
# eigenbasis, then compare with shifted sparse solves on the same matrices.
worst <- 0
for(h in c(.25,.1)) {
  data <- e$isolated_tadpole_data(h)
  spectrum <- e$isolated_discrete_spectrum(data$C,data$L)
  rhs <- as.numeric(data$C%*%data$initial+.03*data$spatial_load)
  for(beta in c(seq(.25,.95,.1),1/3)) for(m in c(1L,4L,8L,14L)) {
    cf <- e$brasil_coefficients(m,beta)
    raw <- get(paste0('m_brasil',m,'t'),asNamespace('rSPDE'))
    column <- function(name) stats::approx(raw$alpha,raw[[name]],beta)$y
    y <- spectrum$values/data$kappa^2
    inverse <- rep(column('k'),length(y))
    for(j in seq_len(m))inverse <- inverse+column(paste0('r',j))/(y-column(paste0('p',j)))
    stopifnot(max(abs(e$brasil_inverse(y,cf)-inverse))<1e-12,
      max(abs(e$brasil_power(1/y,cf)-inverse))<1e-12)
    for(tau in c(2e-5,.2)) for(nu in c(.001,1)) {
      op <- e$my.fractional.operators.frac(data$L,beta,data$C,data$kappa^2,m,tau,nu)
      expected <- as.numeric(spectrum$vectors%*%((inverse/(inverse+tau*nu*data$kappa^(2*beta)))*
        as.numeric(Matrix::crossprod(spectrum$vectors,rhs))))
      actual <- as.numeric(e$my.solver.frac(op,rhs))
      gap <- max(abs(expected-actual)); worst <- max(worst,gap)
      stopifnot(gap<2e-10*max(1,max(abs(expected))),length(op$residues)==m,op$remainder>0)
    }
  }
  op <- e$my.fractional.operators.frac(data$L,1,data$C,16,4L,.03,.7)
  stopifnot(max(abs(e$my.solver.frac(op,rhs)-Matrix::solve(data$C+.03*.7*data$L,rhs)))<1e-11)
}
cat('BRASIL table extraction and independent sparse/dense comparisons pass; max gap:',worst,'\n')

# Check complete time histories, including a non-tabulated exponent and alpha=2.
for(alpha in c(.5,2/3,1.9,2)) {
  ref <- e$isolated_exact_be_reference(data$C,spectrum,data$initial,data$spatial_load,
    seq(0,1,length.out=21),alpha,.3)
  for(m in c(1L,4L)) {
    value <- e$isolated_rational_error(data,spectrum,ref,m)
    stopifnot(value$absolute_error_disagreement<1e-9,value$final_state_gap_C<1e-9)
    if(alpha==2)stopifnot(value$error<1e-10)
  }
}

# Preserved saved settings and errors are checked before any reproduction.
h <- readRDS('data_files/rough_h_results.rds')
t <- readRDS('data_files/rough_tau_results.rds')
iso <- readRDS('data_files/isolated_m_results.rds')
stopifnot(identical(h$parameters,e$rough_parameters()),identical(t$parameters,e$rough_tau_parameters()),
  isTRUE(all.equal(iso$parameters,e$isolated_parameters())),nrow(iso$results)==112L)
# Reproduce a cheap original spatial case using the embedded matrix kernel.
e$rough_compile_solver()
J <- h$parameters$J[1]; alpha <- h$parameters$alpha[1]
mesh <- e$rough_mesh(J)
value <- e$rough_solve(mesh,alpha,h$parameters$time_steps[1],order=h$parameters$rational_order)
stopifnot(abs(value$error/h$errors[1,1]-1)<1e-9)
# Full and invariant-subspace matrix evolution must agree.
full <- e$rough_solve(mesh,alpha,.1,order=4L,symmetry_reduction=FALSE,return_states=TRUE)
reduced <- e$rough_solve(mesh,alpha,.1,order=4L,return_states=TRUE)
stopifnot(abs(full$error-reduced$error)<1e-11)
# Quarter-loop temporal solver uses the same shifted matrix update.
quarter <- e$rough_tau_mesh(J)
projection <- e$rough_tau_projection(quarter)
value <- e$rough_tau_solve(quarter,alpha,10L,4L,1,return_states=TRUE)
op <- e$rough_operator(quarter,alpha,.1,4L,1)
u <- rep(0,nrow(quarter$C))
for(k in 1:10) {
  u <- e$rough_matrix_step(op,as.numeric(quarter$C%*%u+.1*quarter$load))
  stopifnot(max(abs(u-value$states[,k+1L]))<1e-11)
}
# P1 interpolation maps reproduce nodal values and constants.
graph <- e$gets.graph.tadpole(.1);graph$compute_fem()
if(!'edges'%in%names(formals(MetricGraph::metric_graph$public_methods$initialize))) {
  B <- graph$fem_basis(graph$get_mesh_locations())
  stopifnot(max(abs(B-Matrix::Diagonal(nrow(B))))<1e-13)
}
cat('Time histories, preserved parameters, rough baseline and symmetry checks pass.\n')
