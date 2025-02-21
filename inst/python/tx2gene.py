#!/usr/bin/env python3

import argparse
import gzip
import re
import sys
from pathlib import Path
from typing import TextIO

PATTERNS = {
    "gene_id": re.compile(r'gene_id "(.*?)";'),
    "transcript_id": re.compile(r'transcript_id "(.*?)";'),
    "gene_name": re.compile(r'gene_name "(.*?)";'),
    "gene_type": re.compile(r'gene_type "(.*?)";'),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate a tx2gene file from GTF",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i", "--input", type=Path, required=True, help="Input GTF file path"
    )
    parser.add_argument(
        "-o", "--output", type=Path, required=True, help="Output file path"
    )
    parser.add_argument(
        "-g", "--gzip", action="store_true", help="Input file is gzip compressed"
    )
    return parser.parse_args()


def open_file(file_path: Path, is_gz: bool) -> TextIO:
    if is_gz or str(file_path).endswith(".gz"):
        return gzip.open(file_path, "rt")
    return open(file_path)


def parse_gtf_line(line: str, fields: dict[str, dict[str, str]]) -> None:
    if line.startswith("#"):
        return

    parts = line.strip().split("\t")
    if parts[2] == "gene":
        return

    chromosome = parts[0]
    attributes = parts[8]

    matches = {field: pattern.search(attributes) for field, pattern in PATTERNS.items()}

    if not all(matches.values()):
        missing = [f for f, m in matches.items() if not m]
        sys.stderr.write(f"Warning: line missing fields: {', '.join(missing)}\n")

    transcript_id = matches["transcript_id"].group(1)
    gene_id = matches["gene_id"].group(1)
    gene_name = matches["gene_name"].group(1)
    gene_type = matches["gene_type"].group(1)

    fields["tx_to_gene"][transcript_id] = gene_id
    fields["tx_to_name"][transcript_id] = gene_name
    fields["gene_to_chr"][gene_id] = chromosome
    fields["gene_to_type"][gene_id] = gene_type


def write_output(output_path: Path, fields: dict[str, dict[str, str]]) -> None:
    with open(output_path, "w") as f:
        for transcript_id, gene_id in fields["tx_to_gene"].items():
            chromosome = fields["gene_to_chr"][gene_id]
            gene_name = fields["tx_to_name"][transcript_id]
            gene_type = fields["gene_to_type"][gene_id]
            f.write(
                f"{chromosome}\t{gene_id}\t{gene_name}\t{transcript_id}\t{gene_type}\n"
            )


def main() -> None:
    args = parse_args()

    fields = {"tx_to_gene": {}, "tx_to_name": {}, "gene_to_chr": {}, "gene_to_type": {}}

    with open_file(args.input, args.gzip) as f:
        for line in f:
            parse_gtf_line(line, fields)

    write_output(args.output, fields)


if __name__ == "__main__":
    main()
