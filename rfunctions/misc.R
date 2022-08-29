#' Save objects
#'
#' As devtools::use_data but choosing dest_dir 
#'  
#' @export
mysave <- function (..., dir = ".", overwrite = FALSE,
    compress = "bzip2") {

  objs <- eval(substitute(alist(...)))
  objs <- vapply(objs, as.character, character(1))
  paths <- file.path(dir, paste0(objs, ".rda"))

  if (any(file.exists(paths)) & !overwrite) {
    existing_file <- objs[file.exists(paths)]
    stop(paste0(existing_file, " already exists.\n", "Use overwrite = TRUE."))
  }
  message("Saving ", paste(unlist(objs), collapse = ", "),
    " as ", paste(basename(paths), collapse = ", "), " to ",
    dir)
  envir <- parent.frame()
  mapply(save, list = objs, file = paths, MoreArgs = list(envir = envir,
      compress = compress))
  invisible()
}


myload <- function (..., dir = ".") {

  objs <- eval(substitute(alist(...)))
  objs <- vapply(objs, as.character, character(1))
  paths <- file.path(dir, paste0(objs, ".rda"))

  if (any(!file.exists(paths))) {
    stop(paste0(existing_file, " does not exist\n"))
  }
  lapply(paths, load, .GlobalEnv)
  invisible()
}


source_dir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "*.R")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}
