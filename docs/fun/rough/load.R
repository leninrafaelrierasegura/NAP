# Source this file from the project root; no computation or notification occurs.
for (.rough_file in c("problem.R", "matrix_solver.R", "experiment.R", "plotting.R", "comparison.R"))
  source(file.path("fun", "rough", .rough_file), local = TRUE)
rm(.rough_file)
