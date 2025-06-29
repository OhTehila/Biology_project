from itertools import product
from Bio.Data.CodonTable import standard_dna_table


BASES = ['A', 'T', 'C', 'G']


def generate_codon_table():
    """
    Generates a dictionary of codons and their corresponding amino acids.
    """
    codon_table = {}

    for codon in map("".join, product(BASES, repeat=3)):  # Generate all possible codons
        amino_acid = standard_dna_table.forward_table.get(codon, '*')  # '*' represents stop codon - default
        codon_table[codon] = amino_acid
    return codon_table


def count_synonymous_positions(codon, codon_table):
    """
    Calculates the number of synonymous positions for a given codon.
    """
    synonymous_positions = 0

    for i in range(3):  # Iterate over each nucleotide position
        non_synonymous_base = 0
        original_base = codon[i]
        for new_base in BASES:
            if new_base != original_base:  # Change only one nucleotide
                mutated_codon = codon[:i] + new_base + codon[i+1:]
                if codon_table[mutated_codon] != codon_table[codon]:  # Check if amino acid remains unchanged
                    non_synonymous_base = non_synonymous_base + 1
                # else:
                #     synonymous_positions += 1
        if non_synonymous_base > 0:
            synonymous_positions = synonymous_positions + 1

    assert (0 <= synonymous_positions <= 3)

    return synonymous_positions


def create_synonymous_positions_dict():
    """
    Computes the synonymous positions for all codons and returns a dictionary.
    """
    codon_table = generate_codon_table()
    synonymous_dict = {codon: count_synonymous_positions(codon, codon_table) for codon in codon_table}

    return synonymous_dict
