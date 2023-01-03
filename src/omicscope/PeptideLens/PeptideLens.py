"""  PeptideLens Module
OmicScope allows user to import data and perform statistical analysis
for differential gene/protein expression. This is the key module of
all OmicScope workflow, being the input for all other modules.

"""
from .Omicscope_peptides import OmicScope_Peptide
import pandas as pd
from copy import copy
from sklearn.neighbors import LocalOutlierFactor
import requests
import sys
import json
import numpy as np


class PeptideLens(OmicScope_Peptide):
    def __init__(self, Table, ControlGroup, Method, FoldChange_cutoff=0,
                 PValue_cutoff=1, logTransformed=False, ExcludeKeratins=True,
                 **kwargs):
        super().__init__(Table, ControlGroup, Method, FoldChange_cutoff=0,
                         PValue_cutoff=1, logTransformed=False, ExcludeKeratins=True,
                         **kwargs)
        self.PValue_cutoff = PValue_cutoff
        self.Uniprot = self.rest_uniprot()
        self.quant_Uniprot = self.Merge_Uniprot()
        self.quant_Uniprot = self.motif()
        self.find_outliers()

    def get_url(self, url, **kwargs):
        response = requests.get(url, **kwargs)
        if not response.ok:
            print(response.text)
            response.raise_for_status()
            sys.exit()
        return response

    def rest_uniprot(self):
        import numpy as np
        df = self.quant_data
        accessions = list(df.Accession.drop_duplicates())
        WEBSITE_API = "https://rest.uniprot.org/"
        joined = ",".join(accessions)
        r = self.get_url(url=(f"{WEBSITE_API}/uniprotkb/accessions?accessions={joined}&fields=id,"
                              "gene_primary,ft_mod_res,ft_carbohyd,sequence"))
        results = json.loads(r.text)['results']
        fasta = []
        accession = []
        gene_name = []
        modification = []
        position = []
        for i in results:
            fasta.append(i['sequence']['value'])
            accession.append(i['primaryAccession'])
            modification.append([i['features'][a]['description'] for a, x in enumerate(i['features'])])
            position.append([i['features'][a]['location']['start']['value']-1 for a, x in enumerate(i['features'])])
            try:
                gene_name.append(i['genes'][0]['geneName']['value'])
            except KeyError:
                gene_name.append(np.nan)

        Protein_Annotation = pd.DataFrame(zip(accession, gene_name, fasta, modification, position),
                                          columns=['Accession', 'gene_name', 'Protein_Sequence',
                                                   'Modification', 'Uniprot_MOD_POS'])
        Protein_Annotation = Protein_Annotation.explode(['Modification', 'Uniprot_MOD_POS'])
        Protein_Annotation = Protein_Annotation.dropna()
        Protein_Annotation['Uniprot_MOD_POS'] = Protein_Annotation.Uniprot_MOD_POS.astype(int)
        Protein_Annotation['Protein_Sequence'] = Protein_Annotation['Protein_Sequence'].astype(str)
        Protein_Annotation['Residue'] = [y[x] for x, y in zip(Protein_Annotation['Uniprot_MOD_POS'],
                                                              Protein_Annotation['Protein_Sequence'])]
        return Protein_Annotation

    def Merge_Uniprot(self):
        df = copy(self.quant_data)
        parser = self.Uniprot[['Accession', 'Protein_Sequence', 'Modification', 'Uniprot_MOD_POS', 'Residue', 'gene_name']]
        df = df.merge(parser, on='Accession', how='left')
        position = []
        for i, j in zip(df.Sequence, df.Protein_Sequence.astype(str)):
            position.append(j.find(i))
        df['Peptide_Start'] = position
        df['Peptide_Start'] = df['Peptide_Start']+1
        df['peptide_lenght'] = df['Sequence'].str.len()
        df['Peptide_End'] = df['Peptide_Start'] + df['peptide_lenght'] - 1
        df = df.drop_duplicates()
        df = df.dropna()
        df['Uniprot_MOD_POS'] = df['Uniprot_MOD_POS'] + 1
        df['Protein_Mod'] = df['gene_name_y']+'_'+df.Residue+df.Uniprot_MOD_POS.astype(int).astype(str)
        df = df[(df['Uniprot_MOD_POS'] > df['Peptide_Start']) & (df['Uniprot_MOD_POS'] < df['Peptide_End'])]
        df = df[df.columns.drop(list(df.filter(regex='gene_name_')))]
        return df

    def motif(self):
        df = self.quant_Uniprot
        motif = []
        for i, seq in zip(df.Uniprot_MOD_POS, df.Protein_Sequence):
            pos = int(i-1+10)
            seq1 = '__________'+seq+'__________'
            sequence = seq1[pos-10:pos+11]
            motif.append(sequence[:10]+sequence[10].lower()+sequence[11:])
        df['motif'] = motif
        return df

    def find_outliers(self):
        df = self.quant_data
        proteins = list(df.Accession.drop_duplicates())
        fc = []
        dfs = []
        for i in proteins:
            df1 = df
            df1 = df1[df1['Accession'] == i]
            df1 = df1[df1['-log10(p)'] >= -np.log10(self.PValue_cutoff)]
            fc.append(df1['log2(fc)'])
            dfs.append(df1)
            del df1
        outlier_detec = []
        for i in fc:
            try:
                i = np.array(i)
                i = i.reshape(-1, 1)
                n_nei = int(len(i)/2)
                clf = LocalOutlierFactor(n_neighbors=n_nei)
                df3 = clf.fit_predict(i)
            except:
                df3 = 'NA'
            outlier_detec.append(df3)
            del i, df3

        data = []
        for out, i, p in zip(outlier_detec, dfs, proteins):
            if out == 'NA':
                i['outlier'] = 1
                i['Accession'] = p
                data.append(i)
            else:
                i['outlier'] = out
                i['Accession'] = p
                data.append(i)
        data = pd.concat(data)
        data = df.merge(data, how='left')
        data = data.merge(self.quant_Uniprot, how='left')
        data = data.sort_values(['pvalue', 'outlier'], ignore_index=True)
        self.annotated_quant_data = data
        PTM = data
        PTM = PTM[(PTM['outlier'] == -1)]
        PTM = PTM.dropna()
        self.predicted_PTMs = PTM
        return data

    def savefile(self, Path: str, data='annotated_quant_data'):
        from copy import copy
        PepLens = copy(self)
        experimental = PepLens.experimental
        string = '-'.join(PepLens.Conditions)
        try:
            if data == 'annotated_quant_data':
                dfAsString = PepLens.annotated_quant_data[['gene_name', 'Accession', 'pvalue', 'log2(fc)', 'TotalMean']].to_csv(
                    sep='\t', index=False)
            elif data == 'quant_data':
                dfAsString = PepLens.quant_data[['gene_name', 'Accession', 'pvalue', 'log2(fc)', 'TotalMean']].to_csv(
                    sep='\t', index=False)
            elif data == 'predicted_PTMs':
                dfAsString = PepLens.predicted_PTMs[['gene_name', 'Accession', 'pvalue', 'log2(fc)', 'TotalMean']].to_csv(
                    sep='\t', index=False)

            with open(Path + '/' + string + '.omics', 'w') as f:
                f.write("OmicScope v1.0.0" + "\n" +
                        "This file is the output performed by OmicScope pipeline and can be used as input" +
                        " for group comparisons having the controling group used as used according to Hubble." +
                        "Please, cite: Reis-de-Oliveira G, Martins-de-Souza D. OmicScope: an Comprehensive Python " +
                        "library for Systems Biology Visualization" +
                        '\nControlGroup:' + '\t' + PepLens.ctrl + '\n' +
                        'Experimental:' + '\t' + '\t'.join(experimental) + '\n' +
                        'Expression:\n' + '-------\n' +
                        dfAsString)
        except UnboundLocalError:
            raise UnboundLocalError(
                'You only can export the file with quantification data (quant_data), or annotated_quant_data, or the predicted PTMs')
