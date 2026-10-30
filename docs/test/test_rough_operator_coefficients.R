# Run from project root with /usr/bin/Rscript --vanilla.
library(Matrix)
source('fun/rough/tau/load.R')
rough_compile_solver(cache_dir='/tmp/rough-operator-rcpp')
coefficients <- read.csv('matlab/computed_fractional_roots.csv')
for(m in 1:8)for(alpha in rough_parameters()$alpha) {
  c <- rough_rational_coefficients(alpha/2,m,'computed')
  stopifnot(c$interpolation=='none',c$coefficient_source=='computed',length(c$poles)==m+1L)
  rows <- coefficients[coefficients$m==m&abs(coefficients$beta-alpha/2)<1e-12,]
  expected <- rows$value[rows$kind=='factor']
  stopifnot(identical(c$roots$factor,expected))
  for(tau in c(.2,.00625,2e-7)) {
    t <- rough_transfer(c,tau)
    stopifnot(t$remainder==0,all(t$residues>0),length(t$poles)==m+1L)
  }
}
stopifnot(inherits(try(rough_rational_coefficients(.333,8,'computed'),silent=TRUE),'try-error'))

# Independent comparison with the old polynomial transfer formula, supplied
# with the SAME direct roots, rather than with interpolated package roots.
reference <- new.env(parent=.GlobalEnv)
for(expr in parse('fun/matrix_reference/fractional_evolution.R'))
  if(is.call(expr)&&identical(expr[[1]],as.name('<-'))&&is.call(expr[[3]])&&
     identical(expr[[3]][[1]],as.name('function')))eval(expr,reference)
mesh <- rough_mesh(4L)
for(m in c(4L,8L))for(alpha in rough_parameters()$alpha) {
  c <- rough_rational_coefficients(alpha/2,m,'computed')
  tau <- .002
  old <- reference$my.fractional.operators.frac(mesh$K,alpha/2,mesh$C,16,
    m=m,time_step=tau,fractional_roots=c$roots)
  current <- rough_operator(mesh,alpha,tau,m)
  rhs <- as.numeric(mesh$C%*%seq_len(nrow(mesh$C))/nrow(mesh$C)+tau*mesh$load)
  a <- rough_matrix_step(current,rhs)
  b <- as.numeric(reference$my.solver.frac(old,rhs))
  stopifnot(max(abs(a-b))<1e-10)
  reduced <- rough_solve(mesh,alpha,tau,m,T_final=.02,return_states=TRUE)
  full <- rough_solve(mesh,alpha,tau,m,T_final=.02,return_states=TRUE,symmetry_reduction=FALSE)
  stopifnot(max(abs(reduced$states-full$states))<1e-12,
            abs(reduced$error-full$error)<1e-12)
  cat('PASS direct root/factor transfer and full matrix equivalence: m=',m,
      ', alpha=',alpha,'\n',sep='')
}
cat('All direct coefficient, positive-pole, original-formula, and symmetry checks passed.\n')
