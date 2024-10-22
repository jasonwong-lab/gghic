"""Utility functions."""

# %%
import tempfile
import subprocess

import cooler
import pandas as pd


# %%
def subset_chr_cooler(
    clr,
    chrs=None,
    filename=None,
):
    if chrs is None:
        chrs = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    bins = clr.bins()[:]
    pixels = clr.pixels()[:]

    bins_filtered = bins[bins["chrom"].isin(chrs)]
    pixels_filtered = pixels[
        pixels["bin1_id"].isin(bins_filtered.index)
        & pixels["bin2_id"].isin(bins_filtered.index)
    ]
    max_bin_id_pixels = max(
        pixels_filtered["bin1_id"].max(), pixels_filtered["bin2_id"].max()
    )
    max_idx_bins = bins_filtered.index[-1]

    if max_bin_id_pixels >= max_idx_bins:
        bins_filtered = pd.concat(
            [
                bins_filtered,
                bins[bins_filtered.index[-1] + 1 : bins_filtered.index[-1] + 2],
            ],
            ignore_index=False,
        )
    index_mapping = {old: new for new, old in enumerate(bins_filtered.index)}
    # Update the index of bins_filtered
    bins_filtered.index = range(len(bins_filtered))
    # Update bin1_id and bin2_id in pixels_filtered
    pixels_filtered.loc[:, "bin1_id"] = pixels_filtered["bin1_id"].map(index_mapping)
    pixels_filtered.loc[:, "bin2_id"] = pixels_filtered["bin2_id"].map(index_mapping)

    if filename is not None:
        f = filename
    else:
        f = tempfile.NamedTemporaryFile(suffix=".cool", delete=False).name

    cooler.create_cooler(
        f,
        bins=bins_filtered,
        pixels=pixels_filtered,
    )

    subp = subprocess.Popen(
        f"cooler balance --force {f}",
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        encoding="utf-8",
    )
    (out, err) = subp.communicate()
    print(f"stdout: {out}")
    print(f"stderr: {err}")

    clr_new = cooler.Cooler(f)

    return clr_new


# %%
