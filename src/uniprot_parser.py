import os
import pandas as pd


def read_uniprot_to_dataframe(file_path="data/uniprot_bacillus_subtilis.xlsx"):
    """
    Reads a UniProt Excel (.xlsx) file and returns a DataFrame.

    :param file_path: Path to the Excel file.
    :return: DataFrame containing the UniProt data.
    :raises FileNotFoundError: If the specified file does not exist.
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"The file '{file_path}' was not found.")

    # Load the Excel file into a DataFrame
    df = pd.read_excel(file_path)
    print(f"File '{file_path}' loaded successfully!")
    return df


def extract_uniprot_genes(df, genbank_genes):
    """
    Extracts gene names from UniProt based on GenBank gene names.
    If a gene from a row exists in GenBank, only that name is added.
    If no match is found, random name from row is added.

    Args:
        df (pd.DataFrame): DataFrame containing UniProt data with 'Gene Names' column.
        genbank_genes (set): Set of gene names from GenBank.

    Returns:
        set: A set of unique gene names from UniProt.
    """

    uniprot_genes = set()
    for genes in df['Gene Names'].dropna():
        gene_list = genes.split(' ')  # Split multiple gene names into individual names
        found_in_genbank = False  # Check if any gene in this row is in GenBank
        for gene in gene_list:
            if gene in genbank_genes:
                uniprot_genes.add(gene)  # If one matches, add only it and stop checking this row
                found_in_genbank = True
                break
        if not found_in_genbank:  # If none of the names are in GenBank, add first names
            uniprot_genes.add(gene_list[0])
    return uniprot_genes
