#' Convert ChromatinContacts to tibble
#'
#' @name as_tibble
#' @aliases as_tibble,ChromatinContacts-method
#' @description
#' Converts interaction or feature data from ChromatinContacts to tibble.
#'
#' @param x ChromatinContacts object.
#' @param which Character. Data to convert: `"interactions"` (default), or
#'   feature name (`"TADs"`, `"loops"`, `"compartments"`, `"multi_contacts"`).
#'
#' @return Tibble with requested data.
#'
#' @examples
#' \dontrun{
#' cc <- ChromatinContacts("sample.cool") |> import()
#' df <- as_tibble(cc)
#' tads_df <- as_tibble(cc, which = "TADs")
#' }
#'
#' @export
methods::setMethod(
  "as_tibble", "ChromatinContacts",
  function(
    x,
    which = c(
      "interactions", "compartments", "TADs", "loops", "multi_contacts"
    )
  ) {
    which <- match.arg(which)

    obj <- methods::slot(x, which)

    if (is.null(obj)) {
      warning("Slot '", which, "' is NULL. Returning empty tibble.")
      return(tibble::tibble())
    }

    obj |>
      as.data.frame() |>
      tibble::as_tibble()
  }
)
