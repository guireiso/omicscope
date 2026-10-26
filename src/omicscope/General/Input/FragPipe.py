
from copy import copy

import numpy as np
import pandas as pd

from .utils import drop_decoys_contaminants, read_table, require_pdata


class Input:

    def __init__(self, Table, pdata=None):
        """ FragPipe output
        FragPipe output is composed by a combined_protein.tsv file and it contains
        several parameters that can be used to filter the data. For TMT/iTRAQ
        experiments, use the TMT-Integrator protein report
        (tmt-report/abundance_protein_MD.tsv, or the _None/_GN variants).

        Args:
            Table (str): Path to the combined_protein.tsv file (label-free) or
            the TMT-Integrator abundance_protein_*.tsv file (isobaric labelling).
            pdata (str): Path to the phenotype data, which is a .xlsx or .xls
            file that contains all the features associated to the samples
        """
        self.Table = Table
        self.rdata_assay()
        self.pdata = require_pdata(pdata, 'FragPipe output',
                                   'the experiment names of combined_protein.tsv (e.g. "IDHmut_1" for '
                                   '"IDHmut_1 MaxLFQ Intensity") or the sample names of the TMT-Integrator report')
        self.Conditions = list(self.pdata['Condition'].drop_duplicates())

    def rdata_assay(self):
        data = read_table(self.Table)
        if 'ReferenceIntensity' in data.columns:
            return self.tmt_integrator(data)
        data = data.rename(columns={'Protein ID':'Accession',
                      'Gene': 'gene_name'},)
        data = drop_decoys_contaminants(data, 'Protein')
        data = data.set_index(['Accession'])
        # '<sample> MaxLFQ Intensity'; fall back to '<sample> Intensity'
        for suffix in [' MaxLFQ Intensity', ' Intensity']:
            columns = [c for c in data.columns if c.endswith(suffix) and
                       not c[:-len(suffix)].endswith((' Unique', ' Razor', ' Total', ' MaxLFQ'))]
            if len(columns) > 0:
                break
        if len(columns) == 0:
            raise KeyError('OmicScope did not find "<sample> MaxLFQ Intensity" or "<sample> Intensity" columns '
                           'in combined_protein.tsv. Please, enable label-free quantification (IonQuant) in FragPipe.')
        assay = data[columns]
        assay.columns = [c[:-len(suffix)] for c in columns]

        rdata = data.iloc[:,~data.columns.str.contains('Intensity')]
        rdata = copy(rdata.reset_index())
        rdata['gene_name'] = rdata['gene_name'].fillna(rdata['Accession'])
        self.rdata = rdata
        self.assay = copy(assay)

    def tmt_integrator(self, data):
        """TMT-Integrator report: log2 abundances in the columns after ReferenceIntensity."""
        accession = 'ProteinID' if 'ProteinID' in data.columns else 'Index'
        data = drop_decoys_contaminants(data, accession)
        data = data.drop_duplicates(accession)
        samples = list(data.columns[list(data.columns).index('ReferenceIntensity') + 1:])
        assay = data.set_index(accession)[samples].apply(pd.to_numeric, errors='coerce')
        # OmicScope log-transforms abundances itself, so values are brought back to linear scale
        assay = np.power(2, assay)
        assay.index.name = 'Accession'
        rdata = data.drop(columns=samples).rename(columns={accession: 'Accession', 'Gene': 'gene_name'})
        if 'gene_name' not in rdata.columns:
            rdata['gene_name'] = rdata['Index'] if accession != 'Index' else rdata['Accession']
        rdata['gene_name'] = rdata['gene_name'].fillna(rdata['Accession'])
        self.rdata = rdata.reset_index(drop=True)
        self.assay = assay
