.onLoad <- function(libname, pkgname) {
  dir_cache <- get_cache_dir()
  ensure_dir(dir_cache)

  myops <- list(FALSE)
  name_pkg <- get_pkg_name()
  names(myops) <- glue::glue("{name_pkg}.clean_cache")
  op <- options()
  toset <- !(names(myops) %in% names(op))
  if (any(toset)) options(myops[toset])

  env <- new.env(parent = emptyenv())
  env$gis <- NULL
  env$n_hic <- 0
  env$n_annotation <- 0
  env$n_track <- 0
  env$n_concatemer <- 0
  env$chrom_add <- NULL
  env$chrom_sub <- NULL
  env$grs_range <- NULL
  assign(".env", env, envir = asNamespace(pkgname))

  invisible()
}

.onUnload <- function(libpath) {
  dir_cache <- get_cache_dir()
  name_pkg <- get_pkg_name()
  if (getOption(glue::glue("{name_pkg}.clean_cache"))) {
    unlink(dir_cache, recursive = TRUE)
  }
}
