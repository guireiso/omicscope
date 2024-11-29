"""   Import PatternLab Output

PatternLab V allows user to export Excel file with filtered identified
proteins, sample information, biological and technical replicates.

Here, we advise users to export data with at least 1 unique peptide/protein
and protein abundance normalized by XIC.

@author: Reis-de-Oliveira G <guioliveirareis@gmail.com>
"""
from copy import copy

import numpy as np
import pandas as pd
import re
from typing import Dict, List, Tuple


class Input:
    def __init__(self, Table):
        """  PatternLab V output for OmicScope input

        Args:
            Table (str): Path to PatternLab V output (.xlsx)
        """
        self.Table = Table
        patternlab = self.PatternLab()
        self.assay = patternlab[0]  # Expression data
        self.pdata = patternlab[1]  # Phenotype data
        self.rdata = patternlab[2]  # row/gene data
        # Extract analyzed conditions
        self.Conditions = list(self.pdata['Condition'].drop_duplicates())

    def PatternLab(self):
        """
        Parses the content of a .plp file and generates a structured DataFrame.
        
        Args:
        - plp_content (str): Content of the .plp file as a string.
        
        Returns:
        - pd.DataFrame: Generated DataFrame with parsed data.
        """
        with open(self.Table, "r") as file:
            plp_content = file.read()
        class_description_dict = {}
        sparse_matrix_data = []
        index_mapping = {}
        secondary_labels = {}

        current_section = None
        file_name = None
        normalization_factors = []

        for line in plp_content.strip().split("\n"):
            line = line.strip()

            # Handle section headers
            if line.startswith("###Description"):
                current_section = "Description"
                continue
            elif line.startswith("###SparseMatrix"):
                current_section = "SparseMatrix"
                continue
            elif line.startswith("###Index"):
                current_section = "Index"
                continue
            elif line.startswith("###SecondaryLabels"):
                current_section = "SecondaryLabels"
                continue
            elif line.startswith("#ClassDescription"):
                match = re.match(r"#ClassDescription\s+(\d+)\s+(.+)", line)
                if match:
                    class_description_dict[int(match.group(1))] = match.group(2).strip()
                continue

            # Parse content based on the current section
            if current_section == "SparseMatrix":
                if line.startswith("#"):
                    # Extract file name if it's a comment line in SparseMatrix
                    file_name = line.split("\\")[-1]
                    continue
                tokens = line.split()
                if len(tokens) < 2:  # Skip malformed lines
                    continue
                class_id = int(tokens[0])
                class_desc = class_description_dict.get(class_id, "Unknown")
                row = {"FileName": f"{class_desc}:: {file_name}"}
                for token in tokens[1:]:  # Skip the ClassID (first token)
                    index, value = map(float, token.split(":"))
                    row[int(index)] = value
                sparse_matrix_data.append(row)
            elif current_section == "Index":
                parts = line.split("\t", maxsplit=2)
                index_id = int(parts[0])
                locus = parts[1].strip()
                description = parts[2].strip() if len(parts) > 2 else ""
                index_mapping[index_id] = {"Locus": locus, "Description": description}
            elif current_section == "SecondaryLabels":
                parts = line.split("\t")
                file_name = parts[1].strip()
                normalization_factor = float(parts[2])
                secondary_labels[file_name] = normalization_factor
                normalization_factors.append(normalization_factor)
        min_norm_factor = min(normalization_factors)
        # Create DataFrame
        data = {}
        for entry in sparse_matrix_data:
            file_name = entry.pop("FileName")
            for index_id, value in entry.items():
                protein_info = index_mapping.get(index_id, {"Locus": f"Protein_{index_id}", "Description": ""})
                locus = protein_info["Locus"]
                description = protein_info["Description"]

                if locus not in data:
                    data[locus] = {"Locus": locus, "Description": description}

                data[locus][file_name] = value*min_norm_factor if value > 0 else np.nan

        # Convert to DataFrame
        df = pd.DataFrame.from_dict(data, orient="index")

        # Move "Description" to the last column
        description_col = df.pop("Description")
        df["Description"] = description_col

        # Reset index for cleaner presentation
        df.reset_index(drop=True, inplace=True)

        #OmicScope data:
        rdata = df.rename({'Locus':'Accession'})[['Locus', 'Description']]
        rdata.columns = ['Accession', 'Description']
        rdata['gene_name'] = rdata['Description'].str.split('GN=').str[-1].str.split(' ').str[0]
        # assay
        assay = df.iloc[:, ~df.columns.isin(['Locus', 'Description'])]
        # pdata
        pdata = pd.DataFrame(assay.columns, columns=['Sample'])
        pdata['Condition'] = pdata['Sample'].str.split(':: ').str[0]
        pdata.Sample = pdata['Sample'].str.split(':: ').str[1]
        pdata['Biological'] = pdata.index
        pdata.Sample = pdata.Biological.astype(str) + '_' + pdata.Sample
        #assay.columns
        assay.columns = list(pdata.Sample)
        return [assay, pdata, rdata]
