## ------------------------------------------------------------------------------------------------------------------------------------------------
# decayed_sine_convolution <- function(t, rate, A = 1, omega = pi) {
#   stopifnot(all(is.finite(t)), all(t >= 0), length(rate) == 1L,
#             is.finite(rate), rate >= 0, is.finite(A), is.finite(omega))
#   z <- (rate + 1i * omega) * t
#   small <- Mod(z) < 0.5
#   answer <- numeric(length(t))
#   if (any(small)) {
#     # exp(-rate*t) * (exp(z)-1)/z; subtract the real constant before
#     # taking the imaginary part so the sine response is O(t^2) accurately.
#     zz <- z[small]
#     term <- rep(1 + 0i, length(zz))
#     series <- rep(0 + 0i, length(zz))
#     for (k in 1:24) {
#       term <- term * zz / (k + 1)
#       series <- series + term
#     }
#     answer[small] <- A * exp(-rate * t[small]) * t[small] * Im(series)
#   }
#   if (any(!small)) {
#     tt <- t[!small]
#     answer[!small] <- A * Im(
#       (exp(1i * omega * tt) - exp(-rate * tt)) / (rate + 1i * omega)
#     )
#   }
#   answer
# }


## ------------------------------------------------------------------------------------------------------------------------------------------------
# exact_modal_coefficients <- function(t, rates, initial, forcing, A = 1, omega = pi) {
#   stopifnot(length(rates) == length(initial), length(rates) == length(forcing),
#             all(is.finite(rates)), all(rates > 0))
#   out <- matrix(0, nrow = length(rates), ncol = length(t))
#   for (j in seq_along(rates)) {
#     out[j, ] <- initial[j] * exp(-rates[j] * t) + forcing[j] *
#       decayed_sine_convolution(t, rates[j], A = A, omega = omega)
#   }
#   out
# }


## ------------------------------------------------------------------------------------------------------------------------------------------------
# gauss_legendre_rule <- function(order = 4L) {
#   stopifnot(length(order) == 1L, is.finite(order), order == as.integer(order),
#             order >= 2L, order <= 32L)
#   k <- seq_len(order - 1L)
#   off <- k / sqrt(4 * k^2 - 1)
#   J <- matrix(0, order, order)
#   J[cbind(k, k + 1L)] <- off
#   J[cbind(k + 1L, k)] <- off
#   ev <- eigen(J, symmetric = TRUE)
#   ix <- order(ev$values)
#   # Nodes/weights on [0,1], with weights summing to one.
#   list(nodes = (ev$values[ix] + 1) / 2, weights = ev$vectors[1, ix]^2)
# }

