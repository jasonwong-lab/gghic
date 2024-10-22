# %%
import re
import os
import gzip
import logging


# %%
logging.basicConfig(level=logging.WARNING, format='%(levelname)s: %(message)s')

REGEX_PATTERNS = {
    'gene_id': re.compile(r'gene_id "(.*?)";'),
    'transcript_id': re.compile(r'transcript_id "(.*?)";'),
    'transcript_name': re.compile(r'transcript_name "(.*?)";'),
    'transcript_type': re.compile(r'transcript_type "(.*?)";'),
    'gene_name': re.compile(r'gene_name "(.*?)";'),
    'gene_type': re.compile(r'gene_type "(.*?)";')
}


# %%
def read_file(file_path, is_gzip=False):
    """Open a file safely, supporting gzip if necessary."""
    try:
        if is_gzip or file_path.endswith(".gz"):
            return gzip.open(file_path, "rt", encoding='utf-8')
        else:
            return open(file_path, "r", encoding='utf-8')
    except IOError as e:
        logging.error(f"Error opening file {file_path}: {e}")
        return None


# %%
def parse_attributes(attributes):
    """Extract values from attributes using predefined regex."""
    result = {}
    for key, pattern in REGEX_PATTERNS.items():
        match = pattern.search(attributes)
        if match:
            result[key] = match.group(1)
        else:
            result[key] = 'unknown'
    return result


# %%
def generate_tx2gene(input_file, output_file, is_gzip=False):
    tx_file_path = os.path.expanduser(input_file)
    tx_file = read_file(tx_file_path, is_gzip)
    if not tx_file:
        return None

    tx_to_gene = {}
    tx_to_gname = {}
    tx_to_tname = {}
    tx_to_type = {}
    tx_to_cse = {}
    gene_to_type = {}

    with tx_file:
        for line in tx_file:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9 or fields[2] == "gene":
                continue

            attr = parse_attributes(fields[8])
            if 'unknown' in attr.values():
                logging.warning(f"Warning: gtf line '{line.strip()}' missing required fields")
                continue

            chromosome = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            tx_to_gene[attr['transcript_id']] = attr['gene_id']
            tx_to_gname[attr['transcript_id']] = attr['gene_name']
            tx_to_tname[attr['transcript_id']] = attr['transcript_name']
            tx_to_type[attr['transcript_id']] = attr['transcript_type']
            tx_to_cse[attr['transcript_id']] = (chromosome, start, end)
            gene_to_type[attr['gene_id']] = attr['gene_type']

    output_path = os.path.expanduser(output_file)
    try:
        with open(output_path, "w", encoding='utf-8') as new_file:
            for t, g in tx_to_gene.items():
                line_new = "\t".join(
                    [
                        tx_to_cse[t][0],
                        str(tx_to_cse[t][1]),
                        str(tx_to_cse[t][2]),
                        g,
                        tx_to_gname[t],
                        t,
                        tx_to_tname[t],
                        gene_to_type[g],
                        tx_to_type[t],
                    ]
                ) + "\n"
                new_file.write(line_new)
    except IOError as e:
        logging.error(f"Error writing to file {output_file}: {e}")


# %%
