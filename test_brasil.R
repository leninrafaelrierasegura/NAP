# Run from the project root. Checks use small matrices and temporary outputs.
grDevices::pdf(file=NULL)
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
stopifnot(nrow(iso$results)==length(iso$parameters$alpha_vector)*iso$parameters$m_max,
  all(h$parameters$alpha==t$parameters$alpha))
sampled_h <- readRDS('data_files/rough_h_as_in_mat_results.rds')
sampled_tau <- readRDS('data_files/rough_tau_as_in_mat_results.rds')
rough_m <- readRDS('data_files/isolated_m_rough_results.rds')
stopifnot(identical(sampled_h$parameters$alpha,sampled_tau$parameters$alpha),
  identical(sampled_h$parameters$alpha,rough_m$parameters$alpha_vector),
  identical(sampled_h$parameters$diffusion,sampled_tau$parameters$diffusion),
  identical(sampled_h$parameters$diffusion,rough_m$parameters$diffusion),
  identical(sampled_h$parameters$h,sampled_tau$tau),length(sampled_tau$tau)==8L,
  max(sampled_tau$tau)/min(sampled_tau$tau)==128,
  all(sampled_tau$parameters$h==rough_m$parameters$h),
  all(sampled_h$parameters$time_steps==rough_m$parameters$tau),
  nrow(rough_m$results)==length(rough_m$parameters$alpha_vector)*rough_m$parameters$m_max,
  all(rough_m$results$resolved))
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

# Verify the new coordinate restriction against physical FEM matrices and an
# independently assembled reference-mass norm, without using the closed-form
# cross products on the physical side of the comparison.
physical <- e$rough_tau_mesh(5L,3)
invariant <- e$rough_invariant_mesh(5L,3)
V <- outer(physical$x,2^invariant$j,function(x,k)sinpi(x*k))/
  rep(sqrt(invariant$mode_mass),each=nrow(physical$C))
reference <- e$rough_sampled_reference(12L,10L)
Psi <- e$rough_sampled_fem_basis(5L,reference$x)
P <- e$rough_invariant_projection(invariant,e$rough_invariant_reference(12L,10L))
stopifnot(max(abs(Matrix::crossprod(V,physical$C%*%V)-diag(4L)))<1e-12,
  max(abs(Matrix::crossprod(V,physical$K%*%V)-as.matrix(invariant$K)))<1e-9,
  max(abs(as.numeric(Matrix::crossprod(V,physical$load))-invariant$load))<1e-12,
  max(abs(Matrix::crossprod(Psi%*%V,reference$C%*%reference$basis)-P$projection))<1e-12)
spectrum_rough <- e$isolated_discrete_spectrum(physical$C,physical$K)
data_rough <- list(C=physical$C,L=physical$K,initial=numeric(nrow(physical$C)),
  spatial_load=physical$load,kappa=3)
reference_rough <- e$isolated_rough_be_reference(data_rough,spectrum_rough,seq(0,.1,by=.005),.7,.2)
original_error <- e$isolated_rational_error(data_rough,spectrum_rough,reference_rough,3L)
e$rough_compile_invariant_solver()
invariant_error <- e$rough_invariant_isolated_error(invariant,.7,.005,3L,.2,.1)
stopifnot(abs(original_error$error-invariant_error$error)<1e-12)

# Exercise the actual notebook chunks with edited settings. All output paths are
# redirected to a temporary directory; committed results must stay byte-identical.
chunk_code <- function(file,label) {
  lines <- readLines(file,warn=FALSE)
  start <- grep(paste0("^```\\{r ",label,"([,}])"),lines)
  stopifnot(length(start)==1L)
  end <- start+which(lines[(start+1L):length(lines)]=="```")[1L]
  paste(lines[(start+1L):(end-1L)],collapse="\n")
}
manual_run <- function(file,changes,kind) {
  output <- tempfile("brasil-workflow-");dir.create(output)
  env <- new.env(parent=e)
  eval(parse(text=chunk_code(file,"parameters")),env)
  for(name in names(changes))assign(name,changes[[name]],env)
  # Evaluate the actual visible notebook code; only route its literal paths
  # into a temporary directory. No experiment runner is substituted.
  stopifnot(any(grepl("^```\\{r recompute, eval=FALSE}",readLines(file,warn=FALSE))))
  code <- chunk_code(file,"recompute")
  code <- gsub('"data_files/',paste0('"',output,'/data_files/'),code,fixed=TRUE)
  eval(parse(text=code),env)
  stopifnot(file.exists(env$output_file),
    length(list.files(output,pattern="\\.png$",recursive=TRUE))>=1L)
  if(grepl("smooth",kind)) {
    saved <- new.env();load(env$output_file,saved)
    result <- list(parameters=attr(saved$errors_projected,"parameters"),errors=saved$errors_projected)
  } else result <- readRDS(env$output_file)
  list(result=result,environment=env)
}
original_files <- list.files("data_files",recursive=TRUE,full.names=TRUE)
original_hashes <- tools::md5sum(original_files)
common <- list(T_final=.1,kappa=3,diffusion=.2,alpha_vector=c(.5,.9),
  m_fixed=2L,overkill_h=.0625,overkill_time_step=.025,
  initial_amplitude=2,forcing_amplitude=3,AAA=.7,OMEGA=2,workers=1L)
x <- manual_run("cov_in_h_mat.Rmd",c(common,list(h_vector=c(.25,.125),time_step=.05)),"smooth-h")
stopifnot(identical(x$result$parameters$axis,c(.25,.125)),x$result$parameters$fixed$tau==.05,
  x$result$parameters$T==.1,x$result$parameters$diffusion==.2,
  x$result$parameters$initial_amplitude==2,x$result$parameters$OMEGA==2,
  identical(dim(x$result$errors),c(2L,2L)))
x <- manual_run("conv_in_tau_mat.Rmd",c(common,list(tau_vector=c(.05,.025),h_fixed=.25)),"smooth-tau")
stopifnot(identical(x$result$parameters$axis,c(.05,.025)),x$result$parameters$fixed$h==.25,
  x$result$parameters$kappa==3,x$result$parameters$m==2L,
  x$result$parameters$forcing_amplitude==3,identical(dim(x$result$errors),c(2L,2L)))
# Printing the scalar diagnostic figures stays on a temporary graphics device.
grDevices::pdf(tempfile(fileext=".pdf"))
x <- manual_run("conv_in_m_isolated.Rmd",list(h=.25,tau=.05,m_max=2L,alpha_vector=c(.5,.9),
  T_final=.1,kappa=3,diffusion=.2,initial_amplitude=2,forcing_amplitude=3,AAA=.7,OMEGA=2),"isolated-m")
grDevices::dev.off()
stopifnot(x$result$parameters$h==.25,x$result$parameters$tau==.05,
  x$result$parameters$m_max==2L,x$result$parameters$T_final==.1,
  x$result$parameters$AAA==.7,x$result$parameters$OMEGA==2,
  identical(dim(x$result$errors),c(2L,2L)),max(x$result$results$absolute_error_disagreement)<1e-10)
x <- manual_run("conv_in_h_rough.Rmd",list(alpha_vector=c(.5,.9),J_vector=3:4,
  time_steps=c(.05,.025),m_fixed=2L,diffusion=.2,T_final=.1,kappa=3,fit_J=3:4,
  projection_quadrature=8L,tail_tolerance=1e-8),"rough-h")
stopifnot(identical(x$result$parameters$J,3:4),x$result$parameters$rational_order==2L,
  x$result$parameters$T_final==.1,x$result$parameters$diffusion==.2,
  identical(dim(x$result$errors),c(2L,2L)))
x <- manual_run("conv_in_tau_rough.Rmd",list(alpha_vector=c(.5,.9),J_vector=c(4L,5L),
  N_vector=c(2L,4L),m_fixed=c(2L,2L),diffusion=c(.2,.3),T_final=.1,kappa=3,
  fit_indices=1:2,projection_quadrature=8L,tail_tolerance=1e-8),"rough-tau")
stopifnot(identical(x$result$parameters$N,c(2L,4L)),identical(x$result$parameters$diffusion,c(.2,.3)),
  identical(x$result$tau,c(.05,.025)),identical(dim(x$result$errors),c(2L,2L)),
  identical(original_hashes,tools::md5sum(original_files)))
sampled_common <- list(alpha_vector=c(.7,1.7),T_final=.1,kappa=3,diffusion=.2,
  m_fixed=3L,reference_cutoff=16L,reference_J=18L,reference_time_step=.002,
  max_tail_squared_fraction=.01)
x <- manual_run("conv_in_tau_rough_as_in_mat.Rmd",c(sampled_common,list(J_fixed=5L,
  J_vector=c(5L,5L),N_vector=c(4L,2L),fit_indices=1:2)),"sampled-rough-tau")
stopifnot(identical(dim(x$result$errors),c(2L,2L)),x$result$parameters$diffusion==.2,
  all(x$result$parameters$J==5L),identical(x$result$tau,c(.025,.05)))
x <- manual_run("cov_in_h_rough_as_in_mat.Rmd",c(sampled_common,list(J_vector=4:3,
  time_steps=c(.005,.005),fit_J=4:3,workers=1L)),"sampled-rough-h")
stopifnot(identical(dim(x$result$errors),c(2L,2L)),all(x$result$parameters$time_steps==.005))
x <- manual_run("conv_in_m_isolated_rough.Rmd",list(alpha_vector=c(.7,1.7),
  T_final=.1,kappa=3,diffusion=.2,J=5L,tau=.005,m_max=3L),"isolated-rough-m")
stopifnot(identical(dim(x$result$errors),c(3L,2L)),x$result$parameters$J==5L,
  x$result$parameters$tau==.005,max(x$result$results$absolute_error_disagreement)<1e-10,
  identical(original_hashes,tools::md5sum(original_files)))
cat("All eight actual notebook workflows honor edited parameters and save data/figures; committed results unchanged.\n")
