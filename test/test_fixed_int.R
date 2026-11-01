# From NAP: /usr/bin/Rscript --vanilla test/test_fixed_int.R
# Requires the Python environment described in fun/fixed_int/README.md.
source("fun/fixed_int/experiment.R")
rejects <- function(expr) inherits(try(force(expr),silent=TRUE),"try-error")
data <- fixed_int_tadpole_data(.25,4)
spectrum <- fixed_int_discrete_spectrum(data$C,data$L)
interval <- fixed_int_interval(spectrum,4)
stopifnot(interval$x_min==1e-3,interval$y_max>=max(spectrum$values/16),
  rejects(fixed_int_interval(spectrum,4,x_min=.01)))
betas <- c(.25,1/3,.95)
bank <- fixed_int_generate(betas,c(1L,3L),interval,workers=2L)
stopifnot(length(bank$coefficients)==6L,
  all(bank$diagnostics$x_min==interval$x_min),
  all(bank$diagnostics$equioscillation_gap<2e-8),
  rejects(fixed_int_coefficients(bank,1L,1/3+1e-6)))
stale <- bank
stale$interval$x_min <- interval$x_min/2
stopifnot(rejects(fixed_int_coefficients(stale,1L,1/3)))

time_seq <- seq(0,.2,by=.01)
C <- as.matrix(data$C)
L <- as.matrix(data$L)
A <- solve(C,L)/16
rhs <- sin(seq_len(nrow(C)))
stopifnot(max(abs(C%*%L-L%*%C))>1e-4)
worst_dense_gap <- worst_history_gap <- 0
for(beta in betas) {
  reference <- fixed_int_reference(data$C,spectrum,data$initial,data$spatial_load,time_seq,2*beta,1)
  for(m in c(1L,3L)) {
    cf <- fixed_int_coefficients(bank,m,beta)
    op <- fixed_int_operator(data$L,data$C,beta,4,.01,1,cf)
    stopifnot(!is.complex(op$poles),!is.complex(op$residues),all(op$poles<0),
      all(op$residues>0),length(op$factorizations)==m,op$remainder>0,!is.null(op$mass_factor))
    # Direct dense inverse-power scheme: same algebra as equation (9), without
    # the time-step partial fractions or any assumption that C and L commute.
    S <- diag(cf$constant,nrow(C))
    for(j in seq_len(m)) S <- S+cf$residues[j]*solve(A-cf$poles[j]*diag(nrow(C)))
    expected <- solve(S+.01*4^(2*beta)*diag(nrow(C)),S%*%solve(C,rhs))
    actual <- as.numeric(fixed_int_solve(op,rhs))
    gap <- max(abs(actual-expected))
    worst_dense_gap <- max(worst_dense_gap,gap)
    stopifnot(gap<2e-10)
    result <- fixed_int_error(data,spectrum,reference,m,bank=bank)
    worst_history_gap <- max(worst_history_gap,result$absolute_error_disagreement)
    stopifnot(result$absolute_error_disagreement<1e-9,result$final_state_gap_C<1e-9)
    yy <- 10^seq(0,3,length.out=257)
    stopifnot(max(abs(fixed_int_inverse(yy,cf)-fixed_int_power(1/yy,cf)))<1e-12,
      identical(fixed_int_power(0,cf),cf$constant))
  }
}

# Stale spectral coverage must fail after a mesh refinement.
fine <- fixed_int_tadpole_data(.02,4)
fine_spectrum <- fixed_int_discrete_spectrum(fine$C,fine$L)
fine_reference <- fixed_int_reference(fine$C,fine_spectrum,fine$initial,fine$spatial_load,
                                     time_seq,2/3,1)
stopifnot(rejects(fixed_int_error(fine,fine_spectrum,fine_reference,1L,bank=bank)),
  fixed_int_interval(fine_spectrum,4)$x_min<interval$x_min)

# Integer powers do not invoke Python and do not require any coefficient rows.
empty <- fixed_int_generate(numeric(),1L,interval,python="/missing/python")
reference <- fixed_int_reference(data$C,spectrum,data$initial,data$spatial_load,time_seq,2,1)
integer <- fixed_int_error(data,spectrum,reference,1L,bank=empty)
stopifnot(integer$error<1e-10,!fixed_int_mark_resolution(integer,integer$error)$resolved,
  fixed_int_slopes(fixed_int_mark_resolution(integer,integer$error))$alpha==2)

# Fresh orders exceed the package table limit: there is no table lookup.
high <- fixed_int_generate(.65,21L,interval,workers=1L)
high_cf <- fixed_int_coefficients(high,21L,.65)
high_op <- fixed_int_operator(data$L,data$C,.65,4,.01,1,high_cf)
stopifnot(length(high_op$poles)==21L,all(high_op$poles<0),all(high_op$residues>0))
cat("PASS: fresh fitting, non-tabulated beta, fixed interval, stale-bank rejection, integer bypass, m=21.\n")
cat("PASS: dense matrix scheme and full time-history comparison; maximum gaps:",
    worst_dense_gap,worst_history_gap,"\n")
