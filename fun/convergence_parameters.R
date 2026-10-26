# Common problem and quadrature for the three tadpole convergence notebooks.
# Changing these values requires rerunning the computation chunks.
convergence_parameters <- list(
  T_final = 1, diffusion = 1, kappa = 4,
  modes = c(5L, 7L), N_finite = 4L,
  initial = c(10, 0), forcing = c(0, 10), A = 1, omega = pi,
  quadrature_order = 4L, step_sum_order = 64L
)
alpha_vector <- c(0.5, 0.9, 1.3)
study_range <- 10^seq(-4, -2, length.out = 6)
m_fixed <- 8L
m_vector <- 1:4
h_fixed <- 1e-5
time_step_fixed <- 1e-10
