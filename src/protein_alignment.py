from Bio import Align
from Bio.Align import substitution_matrices
from Bio.codonalign.codonseq import CodonSeq, cal_dn_ds
from src import models


def align_proteins(seq1, seq2):
    """
    Aligns two protein sequences using the BLOSUM80 matrix.

    :param seq1: First protein sequence.
    :param seq2: Second protein sequence.
    :return: Aligned protein sequences.
    """
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM80")
    aligner.mode = "global"  # Global alignment

    alignments = aligner.align(seq1, seq2)
    best_alignment = alignments[0]  # Take the best scoring alignment

    return best_alignment


def back_translate(protein_alignment, seq):
    base = ''
    n = 0
    for aa in protein_alignment:
        if aa == "-":
            base = base + "---"
        else:
            base = base + seq[n:n + 3]
            n = n + 3
    return base


def compute_dn_ds(cds1_base, cds2_base):
    """
    Computes dN, dS, and dN/dS ratio for two codon-aligned sequences.

    :param cds1: Codon-aligned sequence 1.
    :param cds2: Codon-aligned sequence 2.
    :return: dN, dS, dN/dS ratio, and selection type.
    """
    seq1 = CodonSeq(cds1_base)
    seq2 = CodonSeq(cds2_base)

    dN, dS = cal_dn_ds(seq1, seq2)
    dn_ds_ratio = float(dN / dS) if dS > 0 else float("inf")  # Convert "Infinity" to float

    # Round to 3 decimal places
    rounded_ratio = round(dn_ds_ratio, 3) if dn_ds_ratio != float("inf") else dn_ds_ratio

    # Determine selection type
    if dn_ds_ratio == float("inf"):
        selection_type = "Undefined"
    elif rounded_ratio > 1:
        selection_type = "Positive Selection"
    elif rounded_ratio == 1:
        selection_type = "Neutral Evolution"
    else:
        selection_type = "Negative Selection"

    return models.DNDS(dN, dS, dn_ds_ratio, selection_type)
