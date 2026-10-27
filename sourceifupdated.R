source_if_updated <- function(name) {
  rmd_file  <- here::here(paste0(name, ".Rmd"))
  hash_file <- here::here(paste0(name, ".hash"))
  r_file    <- here::here(paste0(name, ".R"))
  purl_file <- here::here("purl_log.txt")
  
  new_hash <- digest::digest(file = rmd_file, algo = "md5")
  
  old_hash <- if (file.exists(hash_file)) {
    readLines(hash_file, warn = FALSE)
  } else {
    ""
  }
  
  if (!identical(new_hash, old_hash)) {
    capture.output(
      knitr::purl(rmd_file, output = r_file),
      file = purl_file
    )
    writeLines(new_hash, hash_file)
  }
  
  source(r_file)
}