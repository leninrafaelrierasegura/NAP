## ------------------------------------------------------------------------------------------------------------------------------------------------
# myggsave <- function(plot,
#                      width = 9.22,
#                      height = 7.05,
#                      use_mathpazo = TRUE) {
# 
#   dir_to_save <- here::here("data_files/tikzpic")
#   obj_name <- deparse(substitute(plot))
# 
#   # --- Filename suffix ---
#   suffix <- if (use_mathpazo) "_with_mathpazo" else ""
# 
#   tex_name <- file.path(dir_to_save, paste0(obj_name, suffix, ".tex"))
#   pdf_name <- file.path(dir_to_save, paste0(obj_name, suffix, ".pdf"))
#   png_name <- file.path(dir_to_save, paste0(obj_name, suffix, ".png"))
# 
#   # Create directory if needed
#   if (!dir.exists(dir_to_save)) dir.create(dir_to_save, recursive = TRUE)
# 
#   # --- TikZ + LaTeX options ---
#   tikz_packages <- c(
#     "\\usepackage{tikz}",
#     "\\usetikzlibrary{calc}",
#     "\\usepackage[T1]{fontenc}"
#   )
# 
#   if (use_mathpazo) {
#     tikz_packages <- c(tikz_packages, "\\usepackage[sc]{mathpazo}")
#   }
# 
#   options(
#     tikzDocumentDeclaration = "\\documentclass[tikz,border=0pt]{standalone}",
#     tikzLatexPackages = tikz_packages
#   )
# 
#   # Save TikZ plot
#   tikzDevice::tikz(tex_name, standAlone = TRUE, width = width, height = height)
#   print(plot)
#   dev.off()
# 
#   # Compile to PDF
#   system(
#     paste0("pdflatex -output-directory=", dir_to_save, " ", tex_name),
#     ignore.stdout = TRUE,
#     ignore.stderr = TRUE
#   )
# 
#   # --- Convert PDF to PNG (ImageMagick, no trimming) ---
#   cmd <- paste(
#     "convert -density 400",
#     shQuote(pdf_name),
#     shQuote(png_name)
#   )
# 
#   system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
# 
#   if (!file.exists(png_name)) {
#     stop("PNG conversion failed. Likely ImageMagick PDF policy restriction.")
#   }
# 
#   # --- Clean auxiliary files ---
#   aux_ext <- c(".aux", ".log", ".tex")
#   for (ext in aux_ext) {
#     f <- file.path(dir_to_save, paste0(obj_name, suffix, ext))
#     if (file.exists(f)) file.remove(f)
#   }
# 
#   ras_files <- list.files(
#     dir_to_save,
#     pattern = paste0(obj_name, suffix, "_ras[0-9]+\\.png$"),
#     full.names = TRUE
#   )
#   if (length(ras_files) > 0) file.remove(ras_files)
# 
#   # --- Copy PNG to presentation folder ---
#   pres_dir <- path.expand("~/Desktop/leninPresentations/data_files")
#   if (!dir.exists(pres_dir)) dir.create(pres_dir, recursive = TRUE)
# 
#   pres_png <- file.path(pres_dir, basename(png_name))
#   file.copy(png_name, pres_png, overwrite = TRUE)
# 
#   # --- Output ---
#   cat("PDF saved at:", pdf_name, "\n",
#       "PNG saved at:", png_name, "\n",
#       "PNG copied to:", pres_png, "\n")
# }

