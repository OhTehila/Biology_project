import os
from src import genbank_downloader, genbank_parser, models, diagram_maker, uniprot_downloader, uniprot_parser, codon_analysis, protein_alignment
from collections import namedtuple
from collections import defaultdict
from itertools import combinations
import re
import pandas as pd
from Bio import SeqRecord


def partA(data_dir, results_dir):
    """
    Main execution function for analyzing the Bacillus subtilis genome.
    """

    print("\n===== Starting Part 1: Collect & Analyze Genome Data =====\n")

    # Define genome name and ID
    genome_name, genome_id = "Bacillus_subtilis", "NC_000964.3"

    # Download GenBank file
    genbank_downloader.download_genbank_files({genome_name: genome_id}, data_dir)

    # Validate the file exists
    file_path = os.path.join(data_dir, f"{genome_name}.gb")
    assert os.path.exists(file_path), f"Error: File {file_path} does not exist!"

    # Parse GenBank file
    record_gb = genbank_parser.parse_genbank(file_path, genome_id)
    assert type(record_gb) is SeqRecord.SeqRecord, "Error: record_gb not found."

    # Convert GenBank record into a DataFrame
    df = genbank_parser.build_df_from_record(record_gb)

    # 1.1.1 - Count the number of features by type
    print("\n[1.1] Counting Genome Features")

    type_dict = dict(df['Type'].value_counts())
    print("Feature Type dict:", type_dict)

    # 1.1.2 - Count pseudo-genes
    print("\n[1.1.2] Counting Pseudo-genes")

    df['Pseudo'] = df['Qualifiers'].apply(lambda qual: 'pseudo' in qual)
    pseudo_genes = df[(df['Type'] == 'gene') & df['Pseudo']]
    print(f"Number of Pseudo-genes: {len(pseudo_genes)}")

    # 1.2.1 - Compute Gene Lengths
    print("\n[1.2.1] Computing Gene Lengths")

    newdf = pd.DataFrame()

    for _, row in df.iterrows():
        gene = models.Gene(genome_id, row)
        gene_df = gene.to_dataframe()
        gene_df['Location'] = [row['Location']]

        newdf = pd.concat([newdf, gene_df], ignore_index=True)

    cds_genes = newdf[newdf['Type'] == 'CDS']
    avg_cds_length = cds_genes['Length'].mean()
    print(f"Average CDS Length: {avg_cds_length:.4f} bp")

    other_genes = newdf[~newdf.index.isin(cds_genes.index) & ~newdf.index.isin(pseudo_genes.index)]
    assert len(newdf) == (len(other_genes) + len(cds_genes) + len(pseudo_genes)), \
        f"Mismatch in gene counts! Total: {len(newdf)}, Expected: {len(other_genes) + len(cds_genes) + len(pseudo_genes)}"

    results_dir_A = os.path.join(results_dir, "results_A")
    os.makedirs(results_dir_A, exist_ok=True)
    cds_genes.to_csv(os.path.join(results_dir_A, "cds_genes.csv"), index=False)

    # 1.2.2 - Plot histograms of gene lengths
    print("\n[1.2.2] Plotting Gene Length Histograms")
    df['Length'] = df['Location'].apply(lambda loc: len(loc))
    protein_lengths = df[df["Type"] == "CDS"]["Length"].tolist()
    pseudo_lengths = df[(df["Type"] == "gene") & (df["Pseudo"])]["Length"].tolist()
    other_gene_lengths = df[(df["Type"] == "gene") & (~df["Pseudo"])]["Length"].tolist()

    data_lists = [protein_lengths, pseudo_lengths, other_gene_lengths]
    titles = ["CDS", "pseudo", "other genes"]
    length_diagram = os.path.join(results_dir_A, "length_diagram.png")
    diagram_maker.multiple_histograms(data_lists, length_diagram, titles, "Length (bp)", "Count", bins=50)

    # Compute %CA for genes
    Location = namedtuple('Location', ['start', 'end', 'strand'])
    seq = record_gb.seq.upper()
    location = Location(start=0, end=len(seq) + 1, strand=1)

    # 1.3.1 - Compute CA% for the full genome
    print("\n[1.3.1] Computing %CA in the Full Genome")
    ca_rate = genbank_parser.count_ca_rate(location, seq)
    print(f"%CA in Full Genome: {ca_rate:.4f}%")

    # Computing %CA for Each Gene, add start and end to dataframe
    newdf['CA_rate'] = newdf["Location"].apply(lambda loc: genbank_parser.count_ca_rate(loc, seq))
    # newdf.drop('Location')
    newdf.drop('Location', axis=1, inplace=True)

    # 1.3.2(a) - Compute %CA in CDS Genes
    print("\n[1.3.2](a) Computing %CA in CDS Genes")
    cds_ca_rate = newdf[newdf['Type'] == 'CDS']
    avg_cds_ca = cds_ca_rate['CA_rate'].mean()
    print(f"Average %CA in CDS Genes: {avg_cds_ca:.4f}%")

    # 1.3.2(b) - Compute %CA in all Genes
    print("\n[1.3.2](b) Computing %CA in all Genes")
    avg_ca = newdf['CA_rate'].mean()
    print(f"Average %CA in all Genes: {avg_ca:.4f}%")

    # 1.3.3 - Compare %CA values
    print("\n[1.3.3] Comparing %CA Between Full Genome, CDS Genes and All Genes")
    print(f"Full Genome CA%: {ca_rate:.4f}%")
    print(f"CDS Genes CA%: {avg_cds_ca:.4f}%")
    print(f"All Genes CA%: {avg_ca:.4f}%")

    # 1.3.4 - Plot %CA Histogram for CDS Genes
    print("\n[1.3.4] Plotting %CA Distribution Histogram for CDS Genes")
    data = cds_ca_rate['CA_rate']
    title = "Distribution of %CA in CDS Genes"
    CDS_CA_dist_diagram = os.path.join(results_dir_A, "CDS_CA_distribution.png")

    diagram_maker.histogram(data, CDS_CA_dist_diagram, title, "CA Distribution", "Frequency", bins=20)

    data = newdf['CA_rate']
    title = "Distribution of %CA in all Genes"
    genes_CA_dist_diagram = os.path.join(results_dir_A, "genes_CA_distribution.png")

    diagram_maker.histogram(data, genes_CA_dist_diagram, title, "CA Distribution", "Frequency", bins=20, color='skyblue', edgecolor='black')

    # 1.3.5 - Report the Top & Bottom 5 Genes by %CA
    print("\n[1.3.5] Reporting Top & Bottom 5 Genes by %CA")

    sorted_df = cds_ca_rate.sort_values(by='CA_rate')
    sub_df = sorted_df[['Gene', 'Start', 'End', 'Strand', 'CA_rate']]
    print(f"\nTop 5 Genes with Lowest %CA:\n{sub_df.head(5)}")
    print(f"\nTop 5 Genes with Highest %CA:\n{sub_df.tail(5)}")

    # Save final sorted dataset
    sorted_final_df = newdf.sort_values(by='Start')
    sorted_final_df.to_csv(os.path.join(results_dir_A, "part_a.csv"), index='Start')
    print("\nFinal Processed Data Frame saved successfully in 'results/results_A/part_a.csv'.")


def partB(data_dir, results_dir):
    """
    Main execution function for analyzing UniProt data for Bacillus subtilis.
    """

    print("\n===== Starting Part 2: Analyze UniProt Data for Bacillus subtilis =====\n")
    # Download and Read UniProt Excel
    print("\nDownloading and Reading UniProt Excel File")
    uniprot_downloader.download_uniprot_excel()
    df = uniprot_parser.read_uniprot_to_dataframe()
    print("len: ", len(df))

    # [2.1] Comparing GenBank and UniProt Gene Names
    print("\n[2.1] Comparing GenBank and UniProt Gene Names")

    # Step 1: Extract gene names from GenBank
    cds_genes = pd.read_csv('results/results_A/cds_genes.csv')
    genbank_genes = set(cds_genes['Gene'].dropna().unique())
    print(f"Total unique genes in GenBank: {len(genbank_genes)}")

    # Step 2: Extract and clean gene names from UniProt
    uniprot_genes = uniprot_parser.extract_uniprot_genes(df, genbank_genes)
    print(f"Total unique genes in UniProt: {len(uniprot_genes)}")

    # Step 3: Find differences between the two sets
    only_in_genbank = genbank_genes - uniprot_genes
    only_in_uniprot = uniprot_genes - genbank_genes
    in_both = genbank_genes & uniprot_genes

    # Print results
    print(f"Genes only in GenBank: {len(only_in_genbank)}")
    print(f"Genes only in UniProt: {len(only_in_uniprot)}")
    print(f"Genes in Both: {len(in_both)}")

    results_dir_B = os.path.join(results_dir, "results_B")
    os.makedirs(results_dir_B, exist_ok=True)

    # Save Venn Diagram
    data_lists = [genbank_genes, uniprot_genes]
    labels = ["Length (bp)", "Count"]
    title = "'Gene Overlap Between GenBank and UniProt'"
    venn_diagram = os.path.join(results_dir_B, "venn_genbank_uniprot.png")
    diagram_maker.venn_diagram(data_lists, labels, venn_diagram, title)

    # Save Bar Plot
    titles = ['Only in GenBank', 'Only in UniProt', 'In Both']
    data = [len(only_in_genbank), len(only_in_uniprot), len(in_both)]
    title = 'Comparison of Gene Sets Between GenBank and UniProt'
    bar_plot = os.path.join(results_dir_B, "barplot_genbank_uniprot.png")
    diagram_maker.bar_plot(titles, data, bar_plot, title, 'Categories', 'Number of Genes')

    # Save Pie Chart
    title = 'Gene Distribution Between GenBank and UniProt'
    pie_chart = os.path.join(results_dir_B, "piechart_genbank_uniprot.png")
    diagram_maker.pie_chart(titles, data, pie_chart, title)

    # [2.2] Finding 5 Longest and Shortest Proteins in UniProt
    print("\n[2.2] Finding 5 Longest and Shortest Proteins in UniProt")

    # Sort DataFrame by protein length in descending order to find the 5 longest proteins
    longest_proteins = df[['Protein names', 'Length']].sort_values(by='Length', ascending=False).head(5)
    print("\nTop 5 Longest Proteins:")
    print(longest_proteins)

    # Sort DataFrame by protein length in ascending order to find the 5 shortest proteins
    shortest_proteins = df[['Protein names', 'Length']].sort_values(by='Length', ascending=True).head(5)
    print("\nTop 5 Shortest Proteins:")
    print(shortest_proteins)

    # [2.3] Extracting and Analyzing Transmembrane Regions
    print("\n[2.3] Extracting and Analyzing Transmembrane Regions")

    # [2.3.1] Filter proteins with transmembrane regions
    transmembrane_proteins = df.dropna(subset=["Transmembrane"])
    transmembrane_count = len(transmembrane_proteins)
    print(f"[2.3.1] Total proteins with transmembrane regions: {transmembrane_count}")

    # Extract transmembrane regions
    transmembrane_sequences = list()
    for row in transmembrane_proteins.itertuples():
        # Find all patterns like TRANSMEM 29..49 even if there are multiple in one cell
        regions = re.findall(r'TRANSMEM (\d+)\.\.(\d+)', row.Transmembrane)
        for start, end in regions:
            start, end = int(start), int(end)
            transmembrane_sequences.append(row.Sequence[start-1:end])

    # [2.3.2] Total number of transmembrane regions
    total_transmembrane_regions = len(transmembrane_sequences)
    print(f"[2.3.2] Total transmembrane regions found: {total_transmembrane_regions}")

    # Calculate lengths of transmembrane regions
    transmembrane_lengths = [len(seq) for seq in transmembrane_sequences]

    # Statistics for transmembrane lengths
    if transmembrane_lengths:
        min_length = min(transmembrane_lengths)
        max_length = max(transmembrane_lengths)
        avg_length = sum(transmembrane_lengths) / len(transmembrane_lengths)
    else:
        min_length = max_length = avg_length = 0

    print("[2.3.3] Statistics for transmembrane lengths:")
    print(f"Minimum length of transmembrane region: {min_length}")
    print(f"Maximum length of transmembrane region: {max_length}")
    print(f"Average length of transmembrane region: {avg_length:.2f}")

    # Save histogram
    data = transmembrane_lengths
    title = 'Distribution of Transmembrane Region Lengths'
    transmembrane_histogram = os.path.join(results_dir_B, "transmembrane_lengths_histogram.png")

    diagram_maker.no_line_histogram(data, transmembrane_histogram, title, 'Length of Transmembrane Regions', "Frequency", bins=30)

    #  [2.4] Analyzing Hydrophobic Amino Acids in Transmembrane Regions
    print("\n[2.4] Analyzing Hydrophobic Amino Acids in Transmembrane Regions")

    # Step 1: Define hydrophobic amino acids
    hydrophobic_aa = {'A', 'F', 'L', 'I', 'V', 'M', 'P', 'W'}

    total_hydrophobic_percentages = list()

    # Step 2: Calculate hydrophobic content for each transmembrane region
    for sequence in transmembrane_sequences:
        hydrophobic_count = sum(1 for aa in sequence if aa in hydrophobic_aa)
        hydrophobic_percentage = (hydrophobic_count / len(sequence)) * 100
        total_hydrophobic_percentages.append(hydrophobic_percentage)

    # Step 3: Calculate average hydrophobic percentage
    average_hydrophobic_percentage = sum(total_hydrophobic_percentages) / len(total_hydrophobic_percentages)
    print(f"Average hydrophobic content in transmembrane regions: {average_hydrophobic_percentage:.2f}%")

    # Print some examples
    print("Sample hydrophobic content in transmembrane regions:", total_hydrophobic_percentages[:10])

    # Step 4: Save histogram for hydrophobic content distribution
    data = total_hydrophobic_percentages
    title = 'Distribution of Hydrophobic Content in Transmembrane Regions'
    hydrophobic_histogram = os.path.join(results_dir_B, "hydrophobic_content_histogram.png")

    diagram_maker.no_line_histogram(data, hydrophobic_histogram, title, 'Hydrophobic Content (%)', "Frequency")


def partC(data_dir, results_dir):
    """
    Main execution function for analyzing coronavirus genomes.
    """

    print("\n===== Starting Part 3: Codon & Genome Analysis =====\n")

    # 3.1 - Compute synonymous positions for all codons
    print("\n[3.1] Computing synolennymous positions dictionary")

    synonymous_dict = codon_analysis.create_synonymous_positions_dict()
    print(synonymous_dict)

    # 3.2 - Compare the two coronavirus genomes
    print("\n[3.2] Comparing coronavirus genomes")

    genome_ids = {
        "virus_2020": "NC_045512.2",
        "virus_2025": "PV009232.1"
    }

    # Load GenBank files
    genbank_downloader.download_genbank_files(genome_ids, data_dir)

    record_gbs = []
    for genome_id in genome_ids.keys():
        file_path = os.path.join(data_dir, f"{genome_id}.gb")
        assert os.path.exists(file_path), f"Error: File {file_path} does not exist!"
        record_gb = genbank_parser.parse_genbank(file_path, genome_ids[genome_id])
        record_gbs.append(record_gb)

    assert len(record_gbs) >= len(genome_ids),  "Error: Number of records is less than the expected number of files!"

    # 3.2.1 - Print genome lengths
    print("\n[3.2.1] Print genome lengths")

    genome_names = dict()
    for record in record_gbs:
        virus_name = next((key for key, val in genome_ids.items() if val.startswith(record.name)), "Unknown")
        genome_names[record.id] = virus_name

    for record in record_gbs:
        print(f"Genome length of {genome_names[record.id]} ({record.id}): {len(record.seq)} bp")

    # 3.2.2 - Count total genes and protein-coding (CDS) genes
    print("\n[3.2.2] Counting genes and protein-coding (CDS) genes")

    record_dfs = {record_gb.id: genbank_parser.build_df_from_record(record_gb) for record_gb in record_gbs}

    for record_id, record_df in record_dfs.items():
        type_dict = dict(record_df['Type'].value_counts())
        print(f"{genome_names[record_id]} ({record_id}): Total Genes: {type_dict['gene']}, Protein-Coding Genes (CDS): {type_dict['CDS']}")

    # 3.2.3 - Compare shared and unique genes
    print("\n[3.2.3] Comparing shared and unique genes")

    # Step 1: Collect gene sets for all records dynamically
    record_genes = dict()
    for record_id, record_df in record_dfs.items():
        genes = set(record_df[record_df['Gene'] != "Unknown"]['Gene'])
        record_genes[record_id] = genes

    # Step 2: Find shared and unique genes
    shared_genes = set.intersection(*record_genes.values())  # Genes present in all records
    unique_genes = defaultdict(set)  # Stores unique genes for each record

    for record_id, genes in record_genes.items():
        unique_genes[record_id] = genes - shared_genes  # Unique to this record

    # Print results
    print(f"Shared genes across all records: {shared_genes}, amount: {len(shared_genes)}")

    for record_id, unique_set in unique_genes.items():
        print(f"Genes unique to {genome_names[record_id]} ({record_id}): {unique_set if len(unique_set) > 0 else '{}'}, amount: {len(unique_set)}")

    # Step 3.2.4: Compute dN/dS for shared genes
    print("\n[3.2.4] Computing dN/dS for shared genes")

    records_seq = {record_gb.id: record_gb.seq for record_gb in record_gbs}

    # Extract all CDS sequences
    record_cds_dfs = dict()

    for record_id, record_df in record_dfs.items():
        record_cds_df = record_df.copy()
        record_cds_df = record_cds_df[record_cds_df['Type'] == 'CDS']
        record_cds_df['Seq'] = record_cds_df['Location'].apply(lambda loc: genbank_parser.get_seq(loc, records_seq[record_id]))
        record_cds_df['Translation'] = record_cds_df['Qualifiers'].apply(lambda qual: qual.get("translation", [""])[0] if qual else "")

        record_cds_df['Length'] = record_cds_df["Location"].apply(lambda loc: len(loc))
        record_cds_df = record_cds_df.sort_values(by=['Gene', 'Length'], ascending=[True, True])
        record_cds_df = record_cds_df.drop_duplicates(subset='Gene', keep='first')
        record_cds_df = record_cds_df.set_index('Gene')

        record_cds_dfs[record_id] = record_cds_df

    dnds_df = pd.DataFrame()

    # Iterate over all pairs of dataframes
    for record_id_1, record_id_2 in combinations(record_cds_dfs.keys(), 2):
        df_1 = record_cds_dfs[record_id_1]
        df_2 = record_cds_dfs[record_id_2]

        # Find shared genes between the two dataframes
        shared_genes = df_1.index.intersection(df_2.index)
        for gene in shared_genes:
            cds1 = df_1.loc[gene]
            cds2 = df_2.loc[gene]

            # Align proteins and back-translate to DNA
            best_alignment = protein_alignment.align_proteins(cds1['Translation'], cds2['Translation'])
            target, query = best_alignment[0], best_alignment[1]
            cds1_base = protein_alignment.back_translate(target, cds1['Seq'])
            cds2_base = protein_alignment.back_translate(query, cds2['Seq'])
            DNDS_data = protein_alignment.compute_dn_ds(cds1_base, cds2_base)

            # Create a DNDSdf object
            CDS1, CDS2 = models.CDS(record_id_1, cds1), models.CDS(record_id_2, cds2)
            dnds_entry = models.DNDSdf(CDS1, CDS2, DNDS_data)

            # Convert to DataFrame and append to the main DataFrame
            dnds_df = pd.concat([dnds_df, dnds_entry.to_dataframe()], ignore_index=True)

    # Save dN/dS results to CSV
    results_dir_C = os.path.join(results_dir, "results_C")
    os.makedirs(results_dir_C, exist_ok=True)

    dnds_df.to_csv(os.path.join(results_dir_C, "part_c.csv"), index=False)
    print("\nFinal Processed Data Frame saved successfully in 'results/part_c.csv'.")


def main():

    # Define project directories
    project_root = os.path.abspath(os.path.dirname(__file__))
    os.path.exists(project_root)  # sanity check
    data_dir = os.path.join(project_root, "data")
    results_dir = os.path.join(project_root, "results")
    os.makedirs(results_dir, exist_ok=True)  # sanity check

    partA(data_dir, results_dir)
    partB(data_dir, results_dir)
    partC(data_dir, results_dir)


if __name__ == "__main__":
    main()
