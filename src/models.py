import pandas as pd


class Gene:
    """
    Represents a general Gene object.
    """

    def __init__(self, record_id, gene_data):
        """
        Constructor for Gene objects.

        :param record_id: Identifier of the virus (GenBank record).
        :param gene_data: Dictionary containing gene information.
        """
        self.virus_id = record_id
        self.name = gene_data['Qualifiers'].get('gene', ["Unknown"])[0]
        self.type = gene_data['Type']
        # self.location = gene_data['Location']
        # self.qualifiers = gene_data['Qualifiers']
        self.length = len(gene_data['Location'])
        self.Start = gene_data['Location'].start
        self.End = gene_data['Location'].end
        self.Strand = gene_data['Strand']

    def to_dict(self):
        """
        Converts the Gene object into a dictionary for easy CSV export.

        :return: Dictionary representation of the Gene object.
        """
        return {
            "Virus_ID": self.virus_id,
            "Gene": self.name,
            "Type": self.type,
            "Start": self.Start,
            "End": self.End,
            "Strand": self.Strand,
            "Length": self.length,
        }

    def to_dataframe(self):
        gene_dict = self.to_dict()
        return pd.DataFrame([gene_dict])


class CDS(Gene):
    """
    Represents a CDS (Coding Sequence) extracted from a GenBank file.
    """

    def __init__(self, record_id, cds):
        """
        Constructor for CDS objects.

        :param record_id: Identifier of the virus (GenBank record).
        :param cds: Dictionary containing CDS information.
        """
        assert cds['Type'] == 'CDS', "Error: Type is not CDS."

        # Initialize the base class (Gene)
        super().__init__(record_id, cds)

        # Specific CDS information
        self.protein_id = cds['Qualifiers'].get('protein_id', ["Unknown"])[0]
        self.function = cds['Qualifiers'].get('product', ["Unknown Product"])[0]
        self.details = cds['Qualifiers'].get('note', "No additional info")

    def to_dict(self):
        """
        Converts the CDS object into a dictionary for easy CSV export.

        :return: Dictionary representation of the CDS object.
        """
        base_dict = super().to_dict()
        base_dict.update({
            "Protein_ID": self.protein_id,
            "Function": self.function,
            "Details": self.details
        })
        return base_dict


class DNDS:
    """
    Represents a dN/dS calculation result for a pair of aligned CDS sequences.
    """

    def __init__(self, dN, dS, dn_ds_ratio, selection_type):
        """
        Constructor for DNDS objects.

        :param dN: Non-synonymous substitutions.
        :param dS: Synonymous substitutions.
        :param dn_ds_ratio: dN/dS ratio.
        :param selection_type: Type of selection (Positive, Neutral, Negative).
        """
        self.dN = dN
        self.dS = dS
        self.dn_ds_ratio = dn_ds_ratio
        self.selection_type = selection_type

    def __repr__(self):
        """Returns a string representation of the DNDS object for debugging and logging."""
        return (f"DNDS(dN: {self.dN:.3f}, dS: {self.dS:.3f}, dN/dS: {self.dn_ds_ratio}, "
                f"Selection: {self.selection_type})")

    def to_dict(self):
        """Converts the DNDS object to a dictionary for easy CSV export."""
        return {
            "dN": self.dN,
            "dS": self.dS,
            "dN/dS": self.dn_ds_ratio,
            "Selection Type": self.selection_type
        }


class DNDSdf:
    """
    A structured class to generate a dN/dS DataFrame row from two CDS records.

    Attributes:
    - gene: The shared gene name.
    - cds1: The CDS object from the first genome.
    - cds2: The CDS object from the second genome.
    - dnds: A DNDS object containing dN, dS, dN/dS ratio, and selection type.

    Methods:
    - to_dataframe(): Converts the DNDSdf object into a pandas DataFrame row,
                      merging attributes from both CDS objects and the dN/dS results.
    """

    def __init__(self, CDS1, CDS2, DNDS_data):
        """
        Initializes a DNDSdf object.

        :param CDS1: The CDS object from the first genome.
        :param CDS2: The CDS object from the second genome.
        :param DNDS_data: A DNDS object containing dN, dS, dN/dS ratio, and selection type.
        """
        self.gene = CDS1.name
        self.cds1 = CDS1
        self.cds2 = CDS2
        self.dnds = DNDS_data

    def to_dataframe(self):
        """
        Converts the DNDSdf object into a pandas DataFrame row.

        This method merges:
        - CDS information from both genomes (CDS1 and CDS2).
        - dN/dS computation results from the DNDS object.

        :return: A pandas DataFrame containing a single row.
        """
        cds1_dict = self.cds1.to_dict()
        cds2_dict = self.cds2.to_dict()
        dnds_dict = self.dnds.to_dict()
        row = {
            **{f"{key}_1": value for key, value in cds1_dict.items()},
            **{f"{key}_2": value for key, value in cds2_dict.items()},
            **{key: value for key, value in dnds_dict.items()}
        }

        return pd.DataFrame([row])
