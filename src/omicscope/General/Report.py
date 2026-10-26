""" Reporting: a methods paragraph and a script that reproduces the analysis.

methods_text writes the statistical workflow actually applied (taken from the object's settings and
from Params, which records the decisions made at run time, such as the variance model chosen
automatically) with the matching references, ready to adapt for a manuscript.
reproducible_script returns Python code that reruns the same analysis with the omicscope package.
"""
import re

import omicscope

FDR_NAMES = {'fdr_tsbh': 'the two-stage Benjamini-Krieger-Yekutieli procedure',
             'fdr_bh': 'the Benjamini-Hochberg procedure', 'qvalue': "Storey's q-value"}
REFERENCES = {
    'omicscope': 'Reis-de-Oliveira G, et al. OmicScope unravels systems-level insights from quantitative proteomics '
                 'data. Nat Commun 15, 6510 (2024).',
    'limma': 'Smyth GK. Linear models and empirical Bayes methods for assessing differential expression in '
             'microarray experiments. Stat Appl Genet Mol Biol 3, 3 (2004).',
    'fdr_bh': 'Benjamini Y, Hochberg Y. Controlling the false discovery rate. J R Stat Soc B 57, 289-300 (1995).',
    'fdr_tsbh': 'Benjamini Y, Krieger AM, Yekutieli D. Adaptive linear step-up procedures that control the false '
                'discovery rate. Biometrika 93, 491-507 (2006).',
    'qvalue': 'Storey JD, Tibshirani R. Statistical significance for genomewide studies. PNAS 100, 9440-9445 (2003).',
    'storey2005': 'Storey JD, et al. Significance analysis of time course microarray experiments. PNAS 102, '
                  '12837-12842 (2005).',
    'welch': 'Welch BL. The generalization of Student\'s problem when several different population variances are '
             'involved. Biometrika 34, 28-35 (1947).',
    'gameshowell': 'Games PA, Howell JF. Pairwise multiple comparison procedures with unequal n\'s and/or variances. '
                   'J Educ Stat 1, 113-125 (1976).',
    'minprob': 'Lazar C, et al. Accounting for the multiple natures of missing values in label-free quantitative '
               'proteomics data sets. J Proteome Res 15, 1116-1125 (2016).',
    'gseapy': 'Fang Z, et al. GSEApy: a comprehensive package for performing gene set enrichment analysis in Python. '
              'Bioinformatics 39, btac757 (2023).',
    'gsea': 'Subramanian A, et al. Gene set enrichment analysis. PNAS 102, 15545-15550 (2005).',
    'contaminants': 'Frankenfield AM, et al. Protein contaminants matter: building universal protein contaminant '
                    'libraries for DDA and DIA proteomics. J Proteome Res 21, 2104-2113 (2022).',
}


def _params(self):
    return getattr(self, 'Params', {}).get('Params', {})


def methods_text(self, enrichment=None):
    """Methods paragraph describing the analysis that was run, with references.

    Args:
        enrichment (EnrichmentScope, optional): enrichment run on this object, described as well.

    Returns:
        str: the paragraph followed by a numbered reference list.
    """
    p, refs, s = _params(self), ['omicscope'], []
    method = p.get('ImportMethod', getattr(self, 'Method', 'the selected'))
    n_samples = len(self.pdata) if hasattr(self, 'pdata') else None
    source = ('Quantitative data (OmicScope General format)' if str(method).startswith('General')
              else f'Quantitative data exported from {method}')
    s.append(f"{source} were analysed with OmicScope v{omicscope.__version__} [1]"
             + (f" ({n_samples} biological samples, conditions: {', '.join(map(str, self.Conditions))})."
                if n_samples else '.'))
    if 'Stats_Workflow_1' in p and 'performed' in p.get('Stats_Workflow_1', ''):
        s.append('Technical replicates were averaged and proteins quantified in at least two biological replicates '
                 'of every condition were tested.')
        if getattr(self, 'ExcludeContaminants', False):
            s.append('Common contaminants were removed [%d].' % (len(refs) + 1))
            refs.append('contaminants')
        norm = getattr(self, 'normalization_method', None)
        if norm:
            s.append(f'Abundances were {norm}-normalized.')
        imp = getattr(self, 'imputation_method', None)
        if imp == 'minprob':
            s.append('Missing values were imputed by left-censored random draws (MinProb) [%d].' % (len(refs) + 1))
            refs.append('minprob')
        elif imp:
            s.append(f'Missing values were imputed ({imp}).')
        if getattr(self, 'logTransform', False):
            s.append('Abundances were log2-transformed.')
        if getattr(self, 'batch_correction', False):
            s.append(f"Batch effects ({getattr(self, 'batch_column', 'Batch')}) were removed with a linear model "
                     'that preserves the condition effect.')
    design = str(getattr(self, 'ExperimentalDesign', 'static')).lower()
    if design == 'longitudinal':
        s.append('Changes over time were tested with natural cubic splines '
                 f"({getattr(self, 'degrees_of_freedom', 2)} degrees of freedom) comparing models with and "
                 'without condition-specific trajectories [%d].' % (len(refs) + 1))
        refs.append('storey2005')
    elif 'Stats_Workflow_3' in p:
        test = p['Stats_Workflow_3'].replace('OmicScope performed ', '').replace('T-test', 't-test')
        test = test.replace('moderated (empirical Bayes) ', 'moderated ')
        sentence = f'Differential abundance was assessed with a {test}'
        if getattr(self, '_use_moderated', False):
            sentence += ', with protein variances shrunk towards a common prior (empirical Bayes) [%d]' % (len(refs) + 1)
            refs.append('limma')
        if getattr(self, '_use_welch', False):
            sentence += ' [%d]' % (len(refs) + 1)
            refs.append('welch')
        s.append(sentence + '.')
        if len(self.Conditions) > 2:
            post = 'Games-Howell' if getattr(self, '_use_welch', False) else 'Tukey'
            s.append(f'Pairwise differences of significant proteins were located with {post} post-hoc tests'
                     + (' [%d].' % (len(refs) + 1) if post == 'Games-Howell' else '.'))
            if post == 'Games-Howell':
                refs.append('gameshowell')
        ratio = re.search(r'= ([0-9.]+)', p.get('Stats_Variance_Model', ''))
        if ratio and getattr(self, 'welch', None) == 'auto':
            s.append(('Welch' if getattr(self, '_use_welch', False) else 'Pooled-variance') + ' tests were chosen '
                     f'automatically: the median protein variance differed {float(ratio.group(1)):.1f}-fold between '
                     'conditions (Welch is used above 2-fold when every condition has at least four replicates).')
    fdr = getattr(self, 'fdr', 'fdr_bh')
    s.append(f'P-values were adjusted for multiple testing with {FDR_NAMES.get(fdr, fdr)} [%d]' % (len(refs) + 1)
             + f"; proteins with {self.pvalue} <= {self.PValue_cutoff}"
             + (f" and |log2 fold change| >= {self.FoldChange_cutoff}" if getattr(self, 'FoldChange_cutoff', 0) else '')
             + ' were considered differentially abundant.')
    refs.append(fdr if fdr in REFERENCES else 'fdr_bh')
    on_off = re.match(r'(\d+)', p.get('DataHandling_BeforeStat_4', ''))
    if getattr(self, 'detection_layer', False) and on_off and int(on_off.group(1)) > 0:
        s.append(f'{on_off.group(1)} proteins detected in some conditions but absent from others were tested on '
                 'detection frequency '
                 '(Fisher-Freeman-Halton exact test, FDR controlled separately) and their fold changes bounded by '
                 'the detection limit.')
    if enrichment is not None:
        dbs = ', '.join(map(str, enrichment.dbs if isinstance(enrichment.dbs, (list, tuple)) else [enrichment.dbs]))
        if enrichment.Analysis == 'GSEA':
            rank = {'signed_p': 'sign(log2 fold change) x -log10(p)', 'log2fc': 'log2 fold change',
                    'pi_score': 'log2 fold change x -log10(p)'}.get(getattr(enrichment, 'rank_by', 'signed_p'))
            s.append(f'Pre-ranked gene set enrichment analysis [%d] of all quantified proteins ranked by {rank} was '
                     f'run with GSEApy [%d] against {dbs}; gene sets with FDR <= {enrichment.padjust_cutoff} '
                     'were reported.' % (len(refs) + 1, len(refs) + 2))
            refs += ['gsea', 'gseapy']
        else:
            s.append(f'Over-representation analysis of the differentially abundant proteins was run with GSEApy '
                     f'[%d] against {dbs}; terms with adjusted p <= {enrichment.padjust_cutoff} were reported.'
                     % (len(refs) + 1))
            refs.append('gseapy')
    text = ' '.join(s)
    return text + '\n\nReferences\n' + '\n'.join(f'{i}. {REFERENCES[r]}' for i, r in enumerate(refs, 1))


def reproducible_script(self, table='your_data_file', pdata=None, enrichment=None):
    """Python code that reruns this analysis with the omicscope package.

    Args:
        table (str, optional): path of the quantitative file to write in the script.
        pdata (str, optional): path of the pdata file, if one was used.
        enrichment (EnrichmentScope, optional): enrichment to rerun as well.

    Returns:
        str: the script.
    """
    p = _params(self)
    method = getattr(self, 'Method', None) or p.get('ImportMethod', 'General').split(' ')[0]
    args = {'Table': table, 'Method': method, 'ControlGroup': getattr(self, 'ControlGroup', None),
            'ExperimentalDesign': getattr(self, 'ExperimentalDesign', 'static'), 'pvalue': self.pvalue,
            'PValue_cutoff': self.PValue_cutoff, 'FoldChange_cutoff': getattr(self, 'FoldChange_cutoff', 0),
            'normalization_method': getattr(self, 'normalization_method', None),
            'imputation_method': getattr(self, 'imputation_method', None),
            'logTransform': getattr(self, 'logTransform', True),
            'ExcludeContaminants': getattr(self, 'ExcludeContaminants', True),
            'independent_ttest': getattr(self, 'ind_variables', True), 'welch': getattr(self, 'welch', 'auto'),
            'moderated': getattr(self, 'moderated', True), 'fdr': getattr(self, 'fdr', 'fdr_tsbh'),
            'detection_layer': getattr(self, 'detection_layer', True), 'min_detected': getattr(self, 'min_detected', 2)}
    if str(args['ExperimentalDesign']).lower() == 'longitudinal':
        args['degrees_of_freedom'] = getattr(self, 'degrees_of_freedom', 2)
    if getattr(self, 'batch_correction', False):
        args['batch_correction'], args['batch_column'] = True, getattr(self, 'batch_column', 'Batch')
    if pdata:
        args['pdata'] = pdata
    lines = [f'# OmicScope v{omicscope.__version__}', 'import omicscope as omics', '',
             'data = omics.OmicScope(']
    lines += [f'    {k}={v!r},' for k, v in args.items()]
    lines += [')', "data.quant_data.to_csv('quantitative_results.csv', index=False)"]
    if enrichment is not None:
        dbs = enrichment.dbs if isinstance(enrichment.dbs, (list, tuple)) else [enrichment.dbs]
        extra = ''
        if enrichment.Analysis == 'GSEA':
            extra = f", rank_by={getattr(enrichment, 'rank_by', 'signed_p')!r}"
        lines += ['', f'enrichment = omics.EnrichmentScope(data, Analysis={enrichment.Analysis!r}, dbs={list(dbs)!r}, '
                      f'padjust_cutoff={enrichment.padjust_cutoff!r}, organism={str(enrichment.organism)!r}{extra})',
                  "enrichment.results.to_csv('enrichment_results.csv', index=False)"]
    return '\n'.join(lines) + '\n'
