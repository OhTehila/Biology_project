# **Genomic Analysis Project**

## **Overview**

This project is designed for genomic data analysis. It processes and analyzes genetic sequences from GenBank and UniProt, calculates dN/dS ratios, and performs protein alignment. The structure is modular, with each script handling a specific aspect of the workflow.

## **Setup Instructions**

### **1. Install dependencies**

pip install -r requirements.txt

### **2. Run the project**

python main.py

## **Project Structure**

├── data
│   ├── Bacillus_subtilis.gb
├── src
│   ├── genbank_downloader.py
│   ├── genbank_parser.py
│   ├── models.py
│   ├── diagram_maker.py
│   ├── uniprot_downloader.py
│   ├── uniprot_parser.py
│   ├── codon_analysis.py
│   ├── protein_alignment.py
│   └── __init__.py
├── main.py
├── requirements.txt
├── README.md

## **Modules Description**

- **genbank_downloader.py**  
  Downloads GenBank files containing genome annotations for Bacillus subtilis.

- **genbank_parser.py**  
  Parses GenBank files, extracting gene and CDS information for structured analysis.

- **models.py**  
  Defines core data structures for genes and CDS, including methods for data transformation and dN/dS calculations.  
  - `Gene`: Represents a gene with attributes like name, type, location, and length.  
  - `CDS(Gene)`: Inherits from `Gene`, representing coding sequences with additional details such as protein ID and function.  
  - `DNDS`: Stores dN/dS ratio calculations and determines selection type.  
  - `DNDSdf`: Represents dN/dS comparisons between two coding sequences.

- **diagram_maker.py**  
  Generates visual representations of genomic data, including gene positions and sequence alignments.

- **uniprot_downloader.py**  
  Downloads UniProt protein sequence data for Bacillus subtilis, including functional annotations.

- **uniprot_parser.py**  
  Processes UniProt files and extracts protein-related information, matching genes with GenBank data.

- **codon_analysis.py**  
  Performs codon usage and dN/dS ratio analysis to study evolutionary selection pressure on genes.

- **protein_alignment.py**  
  Aligns protein sequences to identify similarities and differences between genes across different samples or species.

- **main.py**  
  The main entry point for running the project. Calls necessary modules to perform data retrieval, processing, and analysis.

- **requirements.txt**  
  Lists required Python dependencies for the project.

- **README.md**  
  Contains project documentation, setup instructions, and an overview of functionalities.

## **Key Dependencies**

- `biopython` - Used for handling biological sequence data, including parsing GenBank files and working with DNA/protein sequences.
- `pandas` - Used for handling tabular data, organizing genomic and protein data, and performing transformations.
- `matplotlib` - Generates visualizations of genomic data, such as gene structures and alignment plots.
- `seaborn` - Enhances data visualization, particularly for statistical analysis of genomic metrics.
- `openpyxl` - Reads and writes Excel files, used for exporting processed data into structured reports.
- `requests` - Handles HTTP requests to fetch GenBank and UniProt data programmatically.
- `matplotlib-venn` - Creates Venn diagrams to visualize overlaps between gene sets from different sources.

## **Additional Notes**
- The project follows an object-oriented approach, utilizing `models.py` to manage genomic data representation.
- The parsing modules (`genbank_parser.py`, `uniprot_parser.py`) work independently but integrate through shared data structures.
