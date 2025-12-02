# Code from https://github.com/clauswilke/relayer

#' Rename aesthetic mappings in a ggplot2 layer
#'
#' @description
#' Advanced utility function for remapping aesthetic names in a ggplot2 layer.
#' Used internally by gghic for creating complex multi-layer visualizations
#' with custom aesthetics.
#'
#' @param layer A ggplot2 layer object.
#' @param new_aes Named character vector. Maps new aesthetic names to original
#'   names (e.g., `c("fill2" = "fill")`).
#' @param mapping Optional aesthetic mapping. If NULL, created from `new_aes`.
#'
#' @return A modified ggplot2 layer with remapped aesthetics.
#'
#' @details
#' This function creates a new ggproto object that remaps aesthetics, allowing
#' multiple layers to use different aesthetic names (e.g., `fill` and `fill2`)
#' with independent scales.
#'
#' This is primarily an internal function used by `geom_hic_under()` to allow
#' overlaying two heatmaps with independent color scales.
#'
#' @examples
#' \dontrun{
#' # Typically used internally
#' layer <- geom_tile(aes(fill = score))
#' renamed_layer <- renameGeomAes(layer, c("fill2" = "fill"))
#' }
#'
#' @keywords internal
#' @export
renameGeomAes <- function(layer, new_aes, mapping = NULL) {
  mapping <- mapping %||% .makeAesMapping(new_aes)

  geom <- layer$geom

  # default_aes_new <- do.call(
  #   aes, plyr::defaults(.makeAesNullMapping(new_aes), geom$default_aes)
  # )

  default_aes_new <- do.call(
    aes, purrr::list_modify(geom$default_aes, !!!.makeAesNullMapping(new_aes))
  )

  draw_key_new <- function(data, params, size) {
    data <- .renameData(data, new_aes)
    geom$draw_key(data, params, size)
  }

  new_geom <- ggproto(
    "remapped", geom,
    default_aes = default_aes_new,
    draw_key = draw_key_new,
    draw_panel = function(self, data, panel_params, coord, ...) {
      data <- .remapData(data, mapping)
      ggplot2::ggproto_parent(geom, self)$draw_panel(data, panel_params, coord, ...)
    },
    aesthetics = function(self) {
      union(ggplot2::ggproto_parent(geom, self)$aesthetics(), new_aes)
    }
  )

  layer$geom <- new_geom
  layer
}

.remapData <- function(data, mapping) {
  evaled <- lapply(mapping, rlang::eval_tidy, data = data)
  evaled <- lapply(evaled, unname)

  for (col in names(evaled)) {
    data[[col]] <- evaled[[col]]
  }
  data
}

.renameData <- function(data, aes_names) {
  old_aes <- names(aes_names)
  if (is.null(old_aes)) {
    return(data)
  }

  new_aes <- unname(aes_names)
  new_names <- names(data)

  for (i in seq_along(new_aes)) {
    new_names[names(data) == new_aes[i]] <- old_aes[i]
  }
  names(data) <- new_names

  data
}

.makeAesMapping <- function(aes_vect) {
  old_aes <- names(aes_vect)
  if (is.null(old_aes)) {
    return(NULL)
  }

  x <- NULL
  for (i in seq_along(old_aes)) {
    x <- c(x, rlang::list2(!!old_aes[i] := sym(aes_vect[i])))
  }

  do.call(aes, x)
}

.makeAesNullMapping <- function(aes_vect) {
  x <- NULL
  for (name in aes_vect) {
    x <- c(x, rlang::list2(!!name := NULL))
  }
  if (is.null(x)) {
    NULL
  } else {
    do.call(aes, x)
  }
}
