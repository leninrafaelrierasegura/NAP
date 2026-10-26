read_mat_with_labels <- function(m, folder = "matlab") {
  # File paths
  info_file <- file.path(folder, sprintf("m%dinfo.txt", m))
  bin_file  <- file.path(folder, sprintf("m%d.bin", m))
  
  # 1. Read dimensions
  dims <- scan(info_file, what = integer(), quiet = TRUE)
  n_rows <- dims[1]
  n_cols <- dims[2]
  
  # 2. Read binary data
  con <- file(bin_file, "rb")
  data <- readBin(con, what = "double", n = n_rows * n_cols, size = 8, endian = "little")
  close(con)
  
  # 3. Reshape into matrix
  mat <- matrix(data, nrow = n_rows, ncol = n_cols, byrow = FALSE)
  
  # 4. Build column names to match rSPDE exactly
  rc_names <- if (m == 1) {
    "rc"                     # special case: m=1 uses just "rc"
  } else {
    sprintf("rc.%d", 1:m)     # otherwise rc.1, rc.2, ...
  }
  rb_names <- sprintf("rb.%d", 1:(m+1))
  
  col_names <- c("beta", "factor", rc_names, rb_names)
  colnames(mat) <- col_names
  
  # Return as data.frame
  as.data.frame(mat)
}

# Example:
m1_data <- read_mat_with_labels(4, folder = "matlab")
head(m1_data)

head(rSPDE:::m4table)

all_data <- lapply(1:8, read_mat_with_labels, folder = "matlab")

saveRDS(all_data, file = "matlab/chebfun_tables.RDS")
