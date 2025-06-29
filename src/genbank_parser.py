from Bio import SeqIO
import os
import pandas as pd
from Bio.Seq import Seq


def parse_genbank(file_path, id):
    """
    Parses a GenBank file and extracts its first sequence record.
    """
    # Check if the file exists before attempting to open it
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Error: The file '{file_path}' was not found.")

    try:
        # Read the GenBank file
        with open(file_path, "r") as input_handle:

            gen = SeqIO.parse(input_handle, "genbank")
            record_gb = next(gen)  # Get the first record

            if record_gb is None:
                raise ValueError(f"Error: The file '{file_path}' does not contain valid GenBank records.")
            return record_gb

    except Exception as e:
        raise ValueError(f"Error while parsing '{file_path}': {e}")


def build_df_from_record(record_gb):
    data = list()
    for feature in record_gb.features:
        feature_data = {
                        # "Locus": record_gb.id,
                        "Type": feature.type,
                        "Location": feature.location,
                        "Qualifiers": feature.qualifiers,
                        "Strand": feature.location.strand,
                        "Gene": feature.qualifiers.get("gene", ["Unknown"])[0]
                    }
        data.append(feature_data)

    df = pd.DataFrame(data)
    return df


def get_seq(location, seq):
    if "join" in str(location):
        ranges = [(int(part.start), int(part.end)) for part in location.parts]
    else:
        ranges = [(location.start, location.end)]
    gene_seq = ''
    for range in ranges:
        gene_seq = gene_seq + seq[range[0]:range[1]+1]

    if location.strand == -1:
        gene_seq = Seq(gene_seq).reverse_complement()
    return gene_seq


def count_ca_rate(location, seq):
    gene_seq = get_seq(location, seq)
    counter_C = gene_seq.count('C')
    counter_A = gene_seq.count('A')
    ca_rate = (counter_A + counter_C) / len(gene_seq) * 100
    return ca_rate
