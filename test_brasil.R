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
stopifnot(nrow(iso$results)==112L,all(h$parameters$alpha==t$parameters$alpha))
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

# Exercise the actual notebook chunks with edited settings. All output paths are
# redirected to a temporary directory; committed results must stay byte-identical.
chunk_code <- function(file,label) {
  lines <- readLines(file,warn=FALSE)
  start <- grep(paste0("^```\\{r ",label,"([,}])"),lines)
  stopifnot(length(start)==1L)
  end <- start+which(lines[(start+1L):length(lines)]=="```")[1L]
  paste(lines[(start+1L):(end-1L)],collapse="\n")
}
manual_run <- function(file,changes,runner) {
  output <- tempfile("brasil-workflow-");dir.create(output)
  env <- new.env(parent=e)
  eval(parse(text=chunk_code(file,"parameters")),env)
  for(name in names(changes))assign(name,changes[[name]],env)
  env[[runner]] <- local({
    target <- runner;directory <- output
    function(...) {
      args <- list(...);args$output_file <- file.path(directory,"results.rds")
      if(target%in%c("run_smooth_h","run_smooth_tau"))args$figure_file <- file.path(directory,"plot.png")
      if(target=="run_isolated_brasil")args$figure_dir <- file.path(directory,"figures")
      do.call(e[[target]],args)
    }
  })
  env$rough_save_figures <- function(result)e$rough_save_figures(result,file.path(output,"figures"))
  env$rough_tau_save_figures <- function(result)e$rough_tau_save_figures(result,file.path(output,"figures"))
  eval(parse(text=chunk_code(file,"recompute")),env)
  stopifnot(file.exists(file.path(output,"results.rds")),
    length(list.files(output,pattern="\\.png$",recursive=TRUE))>=1L)
  if(grepl("smooth",runner)) {
    saved <- new.env();load(file.path(output,"results.rds"),saved)
    result <- list(parameters=attr(saved$errors_projected,"parameters"),errors=saved$errors_projected)
  } else result <- readRDS(file.path(output,"results.rds"))
  list(result=result,environment=env)
}
original_files <- list.files("data_files",recursive=TRUE,full.names=TRUE)
original_hashes <- tools::md5sum(original_files)
common <- list(T_final=.1,kappa=3,diffusion=.2,alpha_vector=c(.5,.9),
  m_fixed=2L,overkill_h=.0625,overkill_time_step=.025,
  initial_amplitude=2,forcing_amplitude=3,AAA=.7,OMEGA=2,workers=1L)
x <- manual_run("cov_in_h_mat.Rmd",c(common,list(h_vector=c(.25,.125),time_step=.05)),"run_smooth_h")
stopifnot(identical(x$result$parameters$axis,c(.25,.125)),x$result$parameters$fixed$tau==.05,
  x$result$parameters$T==.1,x$result$parameters$diffusion==.2,
  x$result$parameters$initial_amplitude==2,x$result$parameters$OMEGA==2,
  identical(dim(x$result$errors),c(2L,2L)))
x <- manual_run("conv_in_tau_mat.Rmd",c(common,list(tau_vector=c(.05,.025),h_fixed=.25)),"run_smooth_tau")
stopifnot(identical(x$result$parameters$axis,c(.05,.025)),x$result$parameters$fixed$h==.25,
  x$result$parameters$kappa==3,x$result$parameters$m==2L,
  x$result$parameters$forcing_amplitude==3,identical(dim(x$result$errors),c(2L,2L)))
# Printing the scalar diagnostic figures stays on a temporary graphics device.
grDevices::pdf(tempfile(fileext=".pdf"))
x <- manual_run("conv_in_m_isolated.Rmd",list(h=.25,tau=.05,m_max=2L,alpha_vector=c(.5,.9),
  T_final=.1,kappa=3,diffusion=.2,initial_amplitude=2,forcing_amplitude=3,AAA=.7,OMEGA=2),"run_isolated_brasil")
grDevices::dev.off()
stopifnot(x$result$parameters$h==.25,x$result$parameters$tau==.05,
  x$result$parameters$m_max==2L,x$result$parameters$T_final==.1,
  x$result$parameters$AAA==.7,x$result$parameters$OMEGA==2,
  identical(dim(x$result$errors),c(2L,2L)),max(x$result$results$absolute_error_disagreement)<1e-10)
x <- manual_run("conv_in_h_rough.Rmd",list(alpha_vector=c(.5,.9),J_vector=3:4,
  time_steps=c(.05,.025),m_fixed=2L,diffusion=.2,T_final=.1,kappa=3,fit_J=3:4,
  projection_quadrature=8L,tail_tolerance=1e-8),"rough_run_experiment")
stopifnot(identical(x$result$parameters$J,3:4),x$result$parameters$rational_order==2L,
  x$result$parameters$T_final==.1,x$result$parameters$diffusion==.2,
  identical(dim(x$result$errors),c(2L,2L)))
x <- manual_run("conv_in_tau_rough.Rmd",list(alpha_vector=c(.5,.9),J_vector=c(4L,5L),
  N_vector=c(2L,4L),m_fixed=c(2L,2L),diffusion=c(.2,.3),T_final=.1,kappa=3,
  fit_indices=1:2,projection_quadrature=8L,tail_tolerance=1e-8),"rough_tau_run_experiment")
stopifnot(identical(x$result$parameters$N,c(2L,4L)),identical(x$result$parameters$diffusion,c(.2,.3)),
  identical(x$result$tau,c(.05,.025)),identical(dim(x$result$errors),c(2L,2L)),
  identical(original_hashes,tools::md5sum(original_files)))
cat("All five actual notebook workflows honor edited parameters and save data/figures; committed results unchanged.\n")
