devtools::load_all()

data <- readRDS("~/Downloads/grs_hic.rds") |>
  scale_data("balanced", log10)

concatemers <- readRDS("~/Downloads/concatemers.rds")

ggplot(
  data,
  aes(
    seqnames1 = seqnames1, start1 = start1, end1 = end1,
    seqnames2 = seqnames2, start2 = start2, end2 = end2,
    fill = score
  )
) +
  geom_hic(rasterize = TRUE) +
  geom_concatemer(concatemer_granges = concatemers, group_identifier = "read_idx")
