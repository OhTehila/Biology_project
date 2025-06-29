import os
import requests


def download_uniprot_excel(file_path="data/uniprot_bacillus_subtilis.xlsx"):

    """
    Downloads the UniProt Excel (.xlsx) file for Bacillus subtilis with specific columns and saves it in the 'data' directory.

    :param file_path: Path to save the Excel file (default: 'data/uniprot_bacillus_subtilis.xlsx').
    """
    if os.path.exists(file_path):
        print(f"File uniprot_bacillus_subtilis.xlsx already exists, skipping: {file_path}")
        return  # Skip downloading if the file already exists

    # URL with correct column names
    url = (
        "https://rest.uniprot.org/uniprotkb/stream?"
        "query=organism_id:224308&format=xlsx&fields="
        "reviewed,id,protein_name,gene_names,organism_name,length,ft_transmem,sequence"
    )

    os.makedirs(os.path.dirname(file_path), exist_ok=True)

    response = requests.get(url)
    if response.status_code == 200:
        with open(file_path, "wb") as file:
            file.write(response.content)
        print(f"UniProt Excel file downloaded successfully at '{file_path}'!")
    else:
        print(f"Failed to download UniProt file: {response.status_code}")
