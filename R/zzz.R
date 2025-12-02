.onLoad <- function(libname, pkgname) {
  dir_cache <- .getCacheDir()
  .ensureDir(dir_cache)

  myops <- list(FALSE)
  name_pkg <- .getPkgName()
  names(myops) <- paste0(name_pkg, ".clean_cache")
  op <- options()
  toset <- !(names(myops) %in% names(op))
  if (any(toset)) options(myops[toset])

  env <- new.env(parent = emptyenv())
  env$gis <- NULL
  env$n_hic <- 0
  env$n_annotation <- 0
  env$n_track <- 0
  env$n_concatemer <- 0
  env$n_hic_under <- 0
  env$chroms_add <- NULL
  env$chroms_sub <- NULL
  env$grs_range <- NULL
  assign(".env", env, envir = asNamespace(pkgname))

  invisible()
}

.onUnload <- function(libpath) {
  dir_cache <- .getCacheDir()
  name_pkg <- .getPkgName()
  if (getOption(paste0(name_pkg, ".clean_cache"))) {
    unlink(dir_cache, recursive = TRUE)
  }
}
