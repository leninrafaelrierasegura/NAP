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

# Independently integrate the full reference-grid error, including reference
# times on solver boundaries, nonuniform weights, and skipped solver intervals.
coarse_C <- Matrix::Diagonal(2L)
coarse_L <- Matrix::Diagonal(x=c(4,9))
Psi <- Matrix::Matrix(rbind(c(1,0),c(.5,.5),c(0,1)),sparse=TRUE)
reference_C <- Matrix::Diagonal(x=c(1,2,1))
reference_basis <- rbind(c(1,0),c(0,1),c(0,.5))
projection <- e$prepare_reference_error(Psi,reference_C,reference_basis)
solver_times <- c(0,.1,.2,.3)
initial <- c(.3,-.1);load <- c(.2,-.1);forcing <- c(0,1,.7,.5)
op <- e$my.fractional.operators.frac(coarse_L,.35,coarse_C,4,4L,.1,.2)
coefficients <- function(t) rbind(exp(-t),sin(2*t))
history <- matrix(0,2L,length(solver_times));history[,1L] <- initial
for(k in seq_len(length(solver_times)-1L))
  history[,k+1L] <- as.numeric(e$my.solver.frac(op,
    coarse_C%*%history[,k]+.1*load*forcing[k+1L]))
for(reference_times in list(c(0,.02,.07,.1,.18,.2,.24,.3),c(0,.05,.3),c(0,.3))) {
  times <- reference_times[-1L]
  steps <- vapply(times,function(t)which(solver_times[-1L]>=t)[1L],integer(1))
  difference <- reference_basis%*%coefficients(times)-Psi%*%history[,steps+1L,drop=FALSE]
  expected <- sqrt(sum(diff(reference_times)*Matrix::colSums(difference*(reference_C%*%difference))))
  actual <- e$solve_fractional_evolution_error(op,.1,solver_times,initial,load,forcing,
    reference_times,coefficients,projection)
  stopifnot(abs(actual-expected)<1e-12*max(1,expected))
}
cat('Reference-grid errors agree with independent direct integration.\n')

# Check the saved results retained for publication.
h <- readRDS('data_files/rough_h_as_in_mat_results.rds')
iso <- readRDS('data_files/isolated_m_results.rds')
stopifnot(nrow(iso$results)==iso$parameters$m_max*length(iso$parameters$alpha_vector),
  identical(h$parameters$study,"rough_h_sampled_reference"),
  identical(dim(h$errors),c(length(h$parameters$J),length(h$parameters$alpha))),
  all(is.finite(h$errors)),all(h$errors>0),
  all(h$max_tail_squared_fraction<=h$parameters$max_tail_squared_fraction))
# P1 interpolation maps reproduce nodal values and constants.
graph <- e$gets.graph.tadpole(.1);graph$compute_fem()
if(!'edges'%in%names(formals(MetricGraph::metric_graph$public_methods$initialize))) {
  B <- graph$fem_basis(graph$get_mesh_locations())
  stopifnot(max(abs(B-Matrix::Diagonal(nrow(B))))<1e-13)
}
cat('Time histories, saved publication results and interpolation checks pass.\n')

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
  eval(parse(text=chunk_code(file,"libraries")),env)
  eval(parse(text=chunk_code(file,"parameters")),env)
  for(name in names(changes))assign(name,changes[[name]],env)
  # Run the actual computation, load, and plotting chunks. Redirect their output
  # paths so checks never overwrite committed experiment results or figures.
  stopifnot(any(grepl("^```\\{r recompute, eval=FALSE}",readLines(file,warn=FALSE))))
  run_chunk <- function(label) {
    code <- gsub('"data_files/',paste0('"',output,'/data_files/'),chunk_code(file,label),fixed=TRUE)
    eval(parse(text=code),env)
  }
  run_chunk("recompute")
  stopifnot(file.exists(env$output_file))
  if(grepl("smooth",kind)) {
    saved <- new.env();load(env$output_file,saved)
    result <- list(parameters=attr(saved$errors_projected,"parameters"),errors=saved$errors_projected)
  } else result <- readRDS(env$output_file)
  lines <- readLines(file,warn=FALSE)
  first <- grep("^```\\{r load-results[,}]",lines)
  headers <- grep("^```\\{r ",lines)
  labels <- sub("^```\\{r ([^,}]+).*","\\1",lines[headers[headers>=first]])
  for(label in labels)run_chunk(label)
  stopifnot(length(list.files(output,pattern="\\.png$",recursive=TRUE))>=1L)
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
x <- manual_run("cov_in_h_rough_as_in_mat.Rmd",list(alpha_vector=c(.5,.9),J_vector=3:4,
  time_steps=c(.05,.025),m_fixed=2L,diffusion=.2,T_final=.1,kappa=3,fit_J=3:4,
  reference_cutoff=10L,reference_J=12L,reference_time_step=.0125,
  max_tail_squared_fraction=1),"rough-h-sampled")
stopifnot(identical(x$result$parameters$J,3:4),x$result$parameters$m==2L,
  x$result$parameters$T_final==.1,x$result$parameters$diffusion==.2,
  x$result$parameters$reference_cutoff==10L,x$result$parameters$reference_J==12L,
  identical(dim(x$result$errors),c(2L,2L)),
  identical(original_hashes,tools::md5sum(original_files)))
cat("All four notebook workflows compute, load, and plot with edited parameters; committed results unchanged.\n")
