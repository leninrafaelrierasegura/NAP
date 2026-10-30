# Consistent P1 matrices on the invariant quarter-loop, with factor four
# accounting for the four signed reflections on the full length-two loop.
rough_tau_mesh <- function(J, kappa = 4) {
  stopifnot(length(J)==1L, is.finite(J), J==as.integer(J), J>=3L, J<=22L,
            length(kappa)==1L, is.finite(kappa), kappa>0)
  h <- 2^-J
  cells <- as.integer(2^(J-1L))
  n <- cells-1L
  C <- Matrix::bandSparse(n,k=c(-1L,0L,1L),
    diagonals=list(rep(2*h/3,n-1L),rep(8*h/3,n),rep(2*h/3,n-1L)))
  G <- Matrix::bandSparse(n,k=c(-1L,0L,1L),
    diagonals=list(rep(-4/h,n-1L),rep(8/h,n),rep(-4/h,n-1L)))
  j <- seq_len(J-1L)
  x <- seq_len(n)*h
  theta <- pi*2^j*h
  B <- outer(x,2^j,function(x,k)sinpi(x*k)) *
    rep(4*h*(sin(theta/2)/(theta/2))^2,each=n)
  list(J=J,h=h,C=C,K=kappa^2*C+G,B=B,j=j,kappa=kappa,
       load=as.numeric(B%*%(1/j)),x=x)
}

rough_tau_projection <- function(mesh, quadrature_order = 16L) {
  stopifnot(quadrature_order>=8L,quadrature_order==as.integer(quadrature_order))
  factor <- Matrix::Cholesky(Matrix::forceSymmetric(mesh$C),LDL=FALSE)
  P <- as.matrix(solve(factor,mesh$B,system="A"))
  rule <- statmod::gauss.quad(quadrature_order,"legendre")
  z <- (rule$nodes+1)/2
  n <- nrow(P)
  left <- rbind(rep(0,ncol(P)),P)
  right <- rbind(P,rep(0,ncol(P)))
  starts <- (0:n)*mesh$h
  residual_gram <- matrix(0,ncol(P),ncol(P))
  for(q in seq_along(z)) {
    exact <- outer(starts+mesh$h*z[q],2^mesh$j,function(x,k)sinpi(x*k))
    residual <- exact-((1-z[q])*left+z[q]*right)
    residual_gram <- residual_gram+2*mesh$h*rule$weights[q]*crossprod(residual)
  }
  list(P=P,residual_gram=residual_gram)
}

rough_tau_solve <- function(mesh,alpha,N,order,diffusion,T_final=1,
    projection=rough_tau_projection(mesh),tail_tolerance=1e-13,
    return_states=FALSE,coefficient_source="computed") {
  stopifnot(N>=1,N==as.integer(N),is.finite(T_final),T_final>0)
  tau <- T_final/N
  operator <- rough_operator(mesh,alpha,tau,order,diffusion,coefficient_source)
  tail <- rough_tail(mesh$J,alpha,diffusion,mesh$kappa,T_final,tail_tolerance)
  general <- function(A)as(as(A,"generalMatrix"),"CsparseMatrix")
  started <- proc.time()[["elapsed"]]
  result <- rough_matrix_evolve_cpp(general(mesh$C),
    lapply(operator$partial_fraction_terms,general),operator$residues,
    operator$remainder,mesh$load,projection$P,projection$residual_gram,
    rough_rates(mesh$j,alpha,diffusion,mesh$kappa),mesh$j,
    tau,as.integer(N),tail$squared,return_states)
  result$elapsed <- proc.time()[["elapsed"]]-started
  result$tau <- tau
  result$N <- N
  result$order <- order
  result$coefficient_source <- coefficient_source
  result$coefficient_table <- operator$transfer$coefficients$table
  result$alpha <- alpha
  result$diffusion <- diffusion
  result$J <- mesh$J
  result$n_solver_nodes <- nrow(mesh$C)
  result$tail <- tail
  result$transfer_reconstruction_error <- operator$transfer$reconstruction_error
  result
}

rough_tau_parameters <- function() {
  alpha <- seq(0.5,1.9,0.2)
  kappa <- 4
  list(alpha=alpha,J=c(20L,16L,13L,11L,10L,9L,8L,8L),
    diffusion=(kappa^2+4*pi^2)^(-alpha/2),
    rational_orders=rep(8L,8L),coefficient_source="computed",N=c(5L,10L,20L,40L,80L,160L),
    fit_indices=4:6,T_final=1,kappa=kappa,projection_quadrature=16L,
    tail_tolerance=1e-13,forcing="same fixed infinite L2 forcing as conv_in_h_rough",
    reconstruction="right-constant backward Euler",
    error_integration="exact time moments, stable projection residual, certified infinite tail")
}

rough_tau_run_experiment <- function(parameters=rough_tau_parameters(),
    output_file="data_files/rough_tau_results.rds",checkpoint_dir=tempfile("rough-tau-"),
    workers=3L,root=getwd(),cache_dir=tempdir()) {
  stopifnot(length(parameters$alpha)==length(parameters$J),
    length(parameters$alpha)==length(parameters$diffusion),
    length(parameters$alpha)==length(parameters$rational_orders),
    all(parameters$fit_indices%in%seq_along(parameters$N)),workers>=1)
  dir.create(checkpoint_dir,recursive=TRUE,showWarnings=FALSE)
  dir.create(dirname(output_file),recursive=TRUE,showWarnings=FALSE)
  rough_compile_solver(root,cache_dir)
  files <- c(list.files(file.path(root,"fun/rough"),pattern="\\.(R|cpp)$",full.names=TRUE),
    list.files(file.path(root,"fun/rough/tau"),pattern="\\.R$",full.names=TRUE))
  provenance <- list(parameters=parameters,sources=tools::md5sum(files),
    coefficient_file=if(parameters$coefficient_source=="computed")
      tools::md5sum(file.path(root,"matlab/computed_fractional_roots.csv")) else NULL)
  run_alpha <- function(a) {
    checkpoint <- file.path(checkpoint_dir,sprintf("alpha%03d.rds",round(100*parameters$alpha[a])))
    if(file.exists(checkpoint)) {
      saved <- readRDS(checkpoint)
      if(!identical(saved$provenance,provenance))stop("Temporal checkpoint source/settings mismatch.")
      return(saved$values)
    }
    message(sprintf("Preparing tau curve: alpha=%.1f, J=%d",parameters$alpha[a],parameters$J[a]))
    mesh <- rough_tau_mesh(parameters$J[a],parameters$kappa)
    projection <- rough_tau_projection(mesh,parameters$projection_quadrature)
    values <- lapply(parameters$N,function(N) {
      value <- rough_tau_solve(mesh,parameters$alpha[a],N,parameters$rational_orders[a],
        parameters$diffusion[a],parameters$T_final,projection,parameters$tail_tolerance,
        coefficient_source=parameters$coefficient_source)
      value$states <- value$final_state <- NULL
      message(sprintf("alpha=%.1f, N=%d, error=%.10g, %.2fs",
        parameters$alpha[a],N,value$error,value$elapsed))
      value
    })
    saveRDS(list(provenance=provenance,values=values),checkpoint)
    values
  }
  safe_alpha <- function(a)tryCatch(run_alpha(a),error=function(e)
    list(failed=TRUE,message=conditionMessage(e)))
  curves <- if(workers>1L && .Platform$OS.type=="unix")
    parallel::mclapply(seq_along(parameters$alpha),safe_alpha,
      mc.cores=min(workers,length(parameters$alpha)),mc.preschedule=FALSE) else
    lapply(seq_along(parameters$alpha),safe_alpha)
  failed <- vapply(curves,function(x)!is.list(x)||isTRUE(x$failed)||inherits(x,"try-error"),logical(1))
  if(any(failed))stop("Temporal matrix experiment failed: ",paste(vapply(curves[failed],
    function(x)if(is.list(x))x$message else as.character(x),character(1)),collapse="; "))
  errors <- vapply(curves,function(x)vapply(x,`[[`,numeric(1),"error"),numeric(length(parameters$N)))
  tau <- parameters$T_final/parameters$N
  fit <- parameters$fit_indices
  slope <- function(e)unname(coef(stats::lm(log10(e[fit])~log10(tau[fit])))[2])
  observed <- apply(errors,2L,slope)
  result <- list(parameters=parameters,tau=tau,errors=errors,
    rates=data.frame(alpha=parameters$alpha,theoretical=1,observed=observed,
      relative_error_percent=100*abs(observed-1)),curves=curves,provenance=provenance,
    generated_at=format(Sys.time(),"%Y-%m-%d %H:%M:%S %Z",tz="Asia/Riyadh"),
    session_info=utils::sessionInfo())
  saveRDS(result,output_file)
  result
}
