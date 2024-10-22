"""Subset the GYB cooler file to include only chromosomes 4 and 11 for the purpose of writing an R package."""

# %%
from pathlib import Path

import cooler

from utils import subset_chr_cooler


# %%
FILE_COOLER = Path("~/projects/ONT/analysis/3c/cooler/GYB/GYB.chrs.mcool").expanduser()
RESOLUTIONS = [100000, 5000]


# %%
for resolution in RESOLUTIONS:
    clr = cooler.Cooler(f"{FILE_COOLER}::/resolutions/{resolution}")
    clr_chr4_11 = subset_chr_cooler(clr, ["chr4", "chr11"], f"../inst/extdata/cooler/chr4_11-{int(resolution / 1000)}kb.cool")


# %%
