# Analytic L2 projection of an equally weighted, normalized sum of tadpole modes.
# Source after brasilfunc.R; mode numbering matches gets.eigen.params().
project_initial_modes_tadpole <- function(graph, initial_modes, amplitude = 10) {
    C <- graph$mesh$C
    PtE <- graph$mesh$PtE
    lengths <- graph$get_edge_lengths()
    stopifnot(length(initial_modes) > 0L, all(is.finite(initial_modes)),
        all(initial_modes >= 1L), all(initial_modes == as.integer(initial_modes)),
        !anyDuplicated(initial_modes), length(amplitude) == 1L, is.finite(amplitude),
        !is.null(C), length(lengths) == 2L, max(abs(lengths - c(1, 2))) < 1e-12,
        nrow(C) == 2L + nrow(PtE))

    # e1 is constant; each positive frequency has phi, and even frequencies
    # additionally have psi. Thus e5=phi_3 and e7=psi_4.
    mode_table <- data.frame(k = 0L, kind = "constant")
    for (k in seq_len(max(initial_modes))) {
        mode_table <- rbind(mode_table, data.frame(k = k, kind = "phi"))
        if (k %% 2L == 0L)
            mode_table <- rbind(mode_table, data.frame(k = k, kind = "psi"))
    }
    load_mode <- function(j) {
        mode <- mode_table[j, ]
        k <- mode$k
        w <- k * pi/2
        load <- numeric(nrow(C))
        for (edge in 1:2) {
            ids <- which(PtE[, 1] == edge)
            ids <- ids[order(PtE[ids, 2])]
            x <- c(0, lengths[edge] * PtE[ids, 2], lengths[edge])
            nodes <- c(if (edge == 1L) 1L else 2L, 2L + ids, 2L)
            A <- B <- 0
            if (k == 0L) {
                A <- 1/sqrt(3)
            } else if (mode$kind == "psi") {
                A <- (if (edge == 1L) (-1)^(k/2) else 1)/sqrt(3/2)
            } else {
                normalization <- if (k %% 2L == 1L) sqrt(3) else 1
                if (edge == 1L) {
                    A <- -2 * sin(pi * k/2)/normalization
                } else {
                    B <- 1/normalization
                }
            }
            for (i in seq_len(length(x) - 1L)) {
                pair <- trig_hat_integrals(x[i], x[i + 1L], w, A, B)
                load[nodes[i]] <- load[nodes[i]] + pair[1L]
                load[nodes[i + 1L]] <- load[nodes[i + 1L]] + pair[2L]
            }
        }
        load
    }
    mode_loads <- vapply(as.integer(initial_modes), load_mode, numeric(nrow(C)))
    load <- amplitude * rowSums(mode_loads)/sqrt(length(initial_modes))
    as.numeric(Matrix::solve(C, load))
}
