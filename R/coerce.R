#' Convert ChromatinContacts to tibble
#'
#' @description
#' Converts interaction data or feature data from a `ChromatinContacts` object
#' to a tibble (data frame).
#'
#' @param x A `ChromatinContacts` object.
#' @param which Character string. Which data to convert: `"interactions"`
#'   (default) for Hi-C interactions, or a feature name like `"TADs"`,
#'   `"loops"`, `"compartments"`, or `"multi_contacts"`.
#'
#' @return A tibble containing the requested data.
#'
#' @examples
#' \dontrun{
#' cc <- ChromatinContacts("sample.cool") |> import()
#'
#' # Convert interactions to tibble
#' df <- as_tibble(cc)
#'
#' # Convert specific feature
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
