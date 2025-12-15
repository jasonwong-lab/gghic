#' @rdname resolution
methods::setGeneric("resolution", function(x) standardGeneric("resolution"))

#' @rdname focus
methods::setGeneric("focus", function(x) standardGeneric("focus"))

#' @rdname features
methods::setGeneric("features", function(x, name) standardGeneric("features"))

#' @rdname features-set
#' @param ... Additional arguments (not used).
methods::setGeneric(
  "features<-",
  function(x, name, ..., value) standardGeneric("features<-")
)

#' @rdname as_tibble
methods::setGeneric(
  "as_tibble", function(x, which = "interactions") standardGeneric("as_tibble")
)

#' @rdname gghic
#' @param ... Additional arguments (not used).
methods::setGeneric("gghic", function(x, ...) standardGeneric("gghic"))

#' @rdname build
#' @param ... Additional arguments (not used).
methods::setGeneric("build", function(x, ...) standardGeneric("build"))

#' @rdname tidy
#' @param ... Additional arguments (not used).
methods::setGeneric("tidy", function(x, ...) standardGeneric("tidy"))

#' @rdname select
#' @param ... Additional arguments (not used).
methods::setGeneric("select", function(x, ...) standardGeneric("select"))

#' @rdname gghypergraph
#' @param ... Additional arguments (not used).
methods::setGeneric(
  "gghypergraph", function(x, ...) standardGeneric("gghypergraph")
)

#' @rdname hypergraphData
methods::setGeneric(
  "hypergraphData",
  function(x, which = c("selected", "tidied")) standardGeneric("hypergraphData")
)
