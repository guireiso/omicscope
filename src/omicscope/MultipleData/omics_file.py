""" The .omics file: one study (differential results and, optionally, enrichment) for Nebula.

Layout (tab separated)::

    OmicScope
    <citation line>
    ControlGroup:   <control>
    Experimental:   <experimental group(s)>
    Statistics:     <p-value column used to call regulated entities>
    DataType:       proteome | transcriptome           (2.0)
    Replicates:     <condition>=<n> ...                 (2.0)
    Expression:
    -------
    <table: gene_name, Accession, <p column>, log2(fc), TotalMean[, pvalue, Mean difference, SE]>

    Enrichment Analysis:
    -------
    <table: Gene_set, Term, [NES,] Overlap, Adjusted P-value, Genes>

The lines marked 2.0 and the extra columns are optional: files written by earlier versions are read as before.
`Mean difference` (difference of log2 means) and its standard error `SE` let Nebula combine studies by effect
size (meta_analysis); `Replicates` documents the sample size of each study.
"""
import io
import os

import pandas as pd

CITATION = ('This file is the output performed by OmicScope pipeline and can be used as input for group comparisons '
            'having the controlling group used as used according to OmicScope. Please, cite: Reis-de-Oliveira G, et al. '
            'OmicScope unravels systems-level insights from quantitative proteomics data. Nat Commun 15, 6510 (2024).')
TRANSCRIPTOME_METHODS = ('RPKM',)


def data_type(omicscope):
    """'transcriptome' for RNA-seq inputs, 'proteome' otherwise."""
    method = str(getattr(omicscope, 'Method', '') or omicscope.Params.get('Params', {}).get('ImportMethod', ''))
    return 'transcriptome' if method.split(' ')[0] in TRANSCRIPTOME_METHODS else 'proteome'


def omics_text(omicscope, enrichment=None):
    """Text of the .omics file of an OmicScope analysis (and of its EnrichmentScope results, if given)."""
    qd = omicscope.quant_data
    columns = ['gene_name', 'Accession', omicscope.pvalue, 'log2(fc)', 'TotalMean']
    columns += [c for c in ('pvalue', 'Mean difference', 'SE') if c in qd.columns and c not in columns]
    expression = qd[columns].to_csv(sep='\t', index=False)
    lines = ['OmicScope', CITATION,
             f'ControlGroup:\t{omicscope.ctrl}',
             'Experimental:\t' + '\t'.join(omicscope.experimental),
             f'Statistics:\t{omicscope.pvalue}',
             f'DataType:\t{data_type(omicscope)}']
    pdata = getattr(omicscope, 'pdata', None)
    if pdata is not None and 'Condition' in pdata.columns:
        unit = 'Biological' if 'Biological' in pdata.columns else 'Sample'
        n = pdata.groupby('Condition')[unit].nunique()
        lines.append('Replicates:\t' + '\t'.join(f'{c}={int(n[c])}' for c in omicscope.Conditions if c in n.index))
    text = '\n'.join(lines) + '\nExpression:\n-------\n' + expression
    if enrichment is not None:
        cols = ['Gene_set', 'Term', 'NES', 'Overlap', 'Adjusted P-value', 'Genes'] if enrichment.Analysis == 'GSEA' \
            else ['Gene_set', 'Term', 'Overlap', 'Adjusted P-value', 'Genes']
        text += '\nEnrichment Analysis:\n-------\n' + enrichment.results[cols].to_csv(sep='\t', index=False)
    return text.replace('\r', '').replace(os.linesep, '\n')


def parse_omics(text):
    """Read the text of a .omics file.

    Returns:
        dict: group, label (experimental group), pvalue (column name), original (expression table),
        enrichment (table or None), metadata (dict with DataType and Replicates when present).
    """
    lines = text.replace('\r\n', '\n').replace('\r', '\n').split('\n')
    marks = [i for i, line in enumerate(lines) if line.strip() == '-------']
    if not marks:
        raise ValueError('This file is not a valid .omics file (no data section found).')
    meta, group, pvalue = {}, None, None
    for line in lines[:marks[0]]:
        key, _, value = line.partition(':')
        value = value.strip('\t').strip()
        if key == 'Experimental':
            group = value.split('\t')[-1]
        elif key == 'Statistics':
            pvalue = value
        elif key == 'DataType':
            meta['DataType'] = value
        elif key == 'Replicates':
            meta['Replicates'] = dict(item.split('=', 1) for item in value.split('\t') if '=' in item)
    end = len(lines)
    for i in range(marks[0] + 1, len(lines)):
        if lines[i].startswith('Enrichment Analysis'):
            end = i
            break
    block = '\n'.join(line for line in lines[marks[0] + 1:end] if line.strip())
    original = pd.read_csv(io.StringIO(block), sep='\t')
    enrichment = None
    later = [m for m in marks if m > end]
    if later:
        table = '\n'.join(line for line in lines[later[0] + 1:] if line.strip())
        if table:
            enrichment = pd.read_csv(io.StringIO(table), sep='\t')
            genes = enrichment['Genes'].astype(str).str.replace("'", '').str.replace('[', '').str.replace(']', '')
            enrichment['Genes'] = genes.str.split(', ').apply(lambda g: [x for x in g if x and x != 'nan'])
    meta.setdefault('DataType', 'proteome')
    return {'group': group, 'label': group, 'pvalue': pvalue, 'original': original, 'enrichment': enrichment,
            'metadata': meta}
