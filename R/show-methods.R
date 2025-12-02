methods::setMethod("show", "ChromatinContacts", function(object) {
  cat("ChromatinContacts object\n")
  cat(strrep("-", 50), "\n", sep = "")

  cat("File: ", basename(object@cooler_path), "\n", sep = "")
  if (!is.null(object@resolution)) {
    cat(
      "Resolution: ", format(object@resolution, big.mark = ","), " bp\n",
      sep = ""
    )
  }

  if (!is.null(object@seqinfo)) {
    n_seq <- length(Seqinfo::seqnames(object@seqinfo))
    cat(
      "Sequences: ", n_seq, " (",
      paste(head(Seqinfo::seqnames(object@seqinfo), 3), collapse = ", "),
      if (n_seq > 3) ", ..." else "", ")\n",
      sep = ""
    )
  }

  if (!is.null(object@focus) && !is.character(object@focus)) {
    n_focus <- length(object@focus)
    cat(
      "Focus: ", n_focus, " region", ifelse(n_focus != 1, "s", ""), "\n",
      sep = ""
    )

    n_show <- min(3, n_focus)
    for (i in seq_len(n_show)) {
      anchor1 <- InteractionSet::anchors(object@focus[i], type = "first")
      anchor2 <- InteractionSet::anchors(object@focus[i], type = "second")

      str_1 <- paste(
        as.character(GenomicRanges::seqnames(anchor1)),
        format(GenomicRanges::start(anchor1), big.mark = ","),
        format(GenomicRanges::end(anchor1), big.mark = ","),
        sep = ":"
      )
      str_2 <- paste(
        as.character(GenomicRanges::seqnames(anchor2)),
        format(GenomicRanges::start(anchor2), big.mark = ","),
        format(GenomicRanges::end(anchor2), big.mark = ","),
        sep = ":"
      )

      cat("  [", i, "] ", str_1, " <-> ", str_2, "\n", sep = "")
    }
    if (n_focus > 3) {
      cat("  ... and ", n_focus - 3, " more\n", sep = "")
    }
  } else if (is.character(object@focus)) {
    cat("Focus: ", object@focus, "\n", sep = "")
  } else {
    cat("Focus: genome-wide\n")
  }

  if (!is.null(object@interactions)) {
    n_interactions <- length(object@interactions)
    cat(
      "Interactions: ", format(n_interactions, big.mark = ","),
      " interaction", ifelse(n_interactions != 1, "s", ""), "\n",
      sep = ""
    )

    mcols_names <- names(S4Vectors::mcols(object@interactions))
    if (length(mcols_names) > 0) {
      cat(
        "  Metadata columns: ", paste(mcols_names, collapse = ", "), "\n",
        sep = ""
      )
    }
  } else {
    cat("Interactions: not loaded (use import() to load)\n")
  }

  features_list <- list(
    Compartments = object@compartments,
    TADs = object@TADs,
    Loops = object@loops,
    `Multi-contacts` = object@multi_contacts,
    Tracks = object@tracks
  )

  features_present <- purrr::keep(
    features_list, ~ !is.null(.x) & length(.x) > 0
  )

  if (length(features_present) > 0) {
    cat("Features:\n")
  }

  for (name in names(features_present)) {
    obj <- features_present[[name]]
    n_features <- length(obj)
    cat(
      "  ", name, ": ", format(n_features, big.mark = ","),
      " region", ifelse(n_features != 1, "s", ""), "\n",
      sep = ""
    )
  }

  if (!is.null(object@tracks) && length(object@tracks) > 0) {
    tracks <- object@tracks
    n_tracks <- length(tracks)
    cat(
      "  Tracks: ", n_tracks, " track", ifelse(n_tracks != 1, "s", ""), "\n",
      sep = ""
    )
    tnames <- names(tracks)
    if (n_tracks > 0) {
      n_show <- min(5, n_tracks)
      for (i in seq_len(n_show)) {
        track <- tracks[[i]]
        n_ranges <- length(track)
        cat(
          "    [", i, "] ", if (!is.null(tnames)) tnames[i] else "", ": ",
          format(n_ranges, big.mark = ","), " range",
          ifelse(n_ranges != 1, "s", ""), "\n",
          sep = ""
        )
      }
      if (n_tracks > 5) {
        cat("    ... and ", n_tracks - 5, " more\n", sep = "")
      }
    }
  }


  cat(strrep("-", 50), "\n", sep = "")
})
