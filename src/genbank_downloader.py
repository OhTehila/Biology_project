import os
from Bio import Entrez

# Set email for Entrez access (required by NCBI)
Entrez.email = "simcha8199@gmail.com"


def download_genbank_files(genome_ids, data_dir):
    """
    Downloads GenBank files for given accession numbers and saves them in the specified directory.
    """
    os.makedirs(data_dir, exist_ok=True)

    for name, accession in genome_ids.items():
        output_file = os.path.join(data_dir, f"{name}.gb")

        if os.path.exists(output_file):
            print(f"File {name}.gb already exists, skipping: {output_file}")
            continue  # Skip downloading if the file already exists

        print(f"Downloading {name} (Accession Number: {accession})...")

        with Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text") as in_handle:
            with open(output_file, "w") as out_handle:
                out_handle.write(in_handle.read())

        print(f"Saved: {output_file}")

    print("All downloads completed successfully!")
