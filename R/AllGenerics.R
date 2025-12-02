#' Generic functions for ChromatinContacts objects
#'
#' @description
#' Generic functions to access and manipulate ChromatinContacts objects.
#'
#' @name generics
#' @keywords internal
NULL

#' @rdname generics
methods::setGeneric("resolution", function(x) standardGeneric("resolution"))

#' @rdname generics
methods::setGeneric("focus", function(x) standardGeneric("focus"))

#' @rdname generics
methods::setGeneric("features", function(x, name) standardGeneric("features"))

#' @rdname generics
methods::setGeneric(
  "features<-",
  function(x, name, ..., value) standardGeneric("features<-")
)

#' @rdname generics
methods::setGeneric(
  "as_tibble", function(x, which = "interactions") standardGeneric("as_tibble")
)

methods::setGeneric("gghic", function(x, ...) standardGeneric("gghic"))
