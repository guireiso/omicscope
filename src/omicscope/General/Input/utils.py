"""   Helpers shared by the import methods

Reading tables and pdata from paths or uploaded files, removing decoy and
contaminant entries and matching pdata sample names to the sample columns
found in the quantitative data.
"""

import re

import pandas as pd

# File extensions that software keep (or not) in run names.
RUN_EXTENSIONS = ('.raw', '.d', '.dia', '.wiff', '.wiff2', '.wiff.scan', '.mzml', '.mzxml',
                  '.mgf', '.htrms', '.ms2', '.tdf', '.mzb', '.quant')

# Prefixes used by search engines/FASTA databases to tag decoys and contaminants.
DECOY_CONTAMINANT = re.compile(r'(?:^|;)\s*(?:REV__|CON__|rev_|Reverse_|DECOY_|decoy_|'
                               r'contam_|contaminant_|Cont_|CONT_)')


def _rewind(file):
    """Uploaded files (e.g. from Streamlit) can be read more than once."""
    if hasattr(file, 'seek'):
        file.seek(0)


def read_table(Table):
    """Read a tab, comma or semicolon separated table (path or file-like)."""
    name = str(getattr(Table, 'name', Table)).lower()
    if name.endswith('.parquet'):
        _rewind(Table)
        try:
            return pd.read_parquet(Table)
        except ImportError:
            raise ImportError('Reading .parquet files requires the "pyarrow" package '
                              '(pip install pyarrow). Alternatively, use the .tsv output '
                              '(e.g. report.pg_matrix.tsv from DIA-NN).')
    if name.endswith(('.xlsx', '.xls')):
        _rewind(Table)
        return pd.read_excel(Table)
    _rewind(Table)
    data = pd.read_csv(Table, sep='\t', low_memory=False)
    if data.shape[1] == 1:
        _rewind(Table)
        data = pd.read_csv(Table, sep=None, engine='python')
        if data.shape[1] == 1:
            raise ValueError('OmicScope could not split the columns of the quantitative data. '
                             'Please, confirm the file is tab-separated (tsv) or comma-separated (csv).')
    return data


def read_pdata(pdata):
    """Read pdata from a DataFrame, an Excel file or a csv/tsv file."""
    if isinstance(pdata, pd.DataFrame):
        return pdata.copy()
    _rewind(pdata)
    try:
        data = pd.read_excel(pdata)
    except ValueError:
        _rewind(pdata)
        data = pd.read_csv(pdata, sep=None, engine='python')
    if data.shape[1] == 1:
        raise KeyError('OmicScope does not import pdata correctly. Please, confirm if data is xlsx, xls or csv.')
    if not {'Sample', 'Condition'}.issubset(data.columns):
        raise KeyError('pdata must contain at least the columns "Sample" and "Condition" '
                       f'(found: {", ".join(map(str, data.columns))}).')
    return data


def require_pdata(pdata, software, sample_hint):
    """Raise a clear error when a method that does not encode conditions gets no pdata."""
    if pdata is None:
        raise ValueError(f'{software} does not contain the sample conditions, so a pdata file is required '
                         '(columns "Sample", "Condition" and, optionally, "Biological"). '
                         f'Sample names must be {sample_hint}.')
    return read_pdata(pdata)


def drop_decoys_contaminants(df, column):
    """Remove rows whose identifiers carry decoy or contaminant prefixes."""
    return df[~df[column].astype(str).str.contains(DECOY_CONTAMINANT)]


def clean_run_name(name):
    """Drop folders and raw-file extensions from a run name.

    'D:\\data\\run_01.d.dia' -> 'run_01'
    """
    name = str(name).strip().replace('\\', '/').rstrip('/').split('/')[-1]
    changed = True
    while changed:
        changed = False
        for ext in RUN_EXTENSIONS:
            if name.lower().endswith(ext) and len(name) > len(ext):
                name = name[:-len(ext)]
                changed = True
    return name


def clean_run_columns(columns):
    """Clean run names, keeping the originals if cleaning would create duplicates."""
    cleaned = [clean_run_name(c) for c in columns]
    if len(set(cleaned)) < len(cleaned):
        return list(columns)
    return cleaned


def _sample_key(name):
    return clean_run_name(name).lower()


def match_samples(assay, pdata):
    """Rename assay columns so that they match the pdata sample names.

    Exact names are kept. Otherwise, names are compared without folders,
    raw-file extensions (.raw, .d, .mzML, ...), surrounding spaces and case.
    """
    pdata = pdata.copy()
    if 'Sample' not in pdata.columns:
        return assay, pdata
    pdata['Sample'] = pdata['Sample'].astype(str).str.strip()
    if 'Condition' in pdata.columns and pdata['Condition'].dtype == object:
        pdata['Condition'] = pdata['Condition'].str.strip()
    assay = assay.copy()
    assay.columns = [str(c) for c in assay.columns]
    missing = [s for s in pdata['Sample'] if s not in assay.columns]
    if missing:
        keys = {}
        for column in assay.columns:
            keys.setdefault(_sample_key(column), []).append(column)
        samples = set(pdata['Sample'])
        rename = {}
        for sample in missing:
            hits = [h for h in keys.get(_sample_key(sample), []) if h not in samples]
            if len(hits) == 1:
                rename[hits[0]] = sample
        assay = assay.rename(columns=rename)
    return assay, pdata


def missing_samples_message(missing, available):
    available = [str(a) for a in available]
    shown = ', '.join(available[:30]) + (', ...' if len(available) > 30 else '')
    return ('These pdata samples were not found in the quantitative data: '
            f'{", ".join(map(str, missing))}.\n'
            f'Sample names found in the data: {shown}.\n'
            'Edit the "Sample" column of pdata so that each name matches one of the names above. '
            'Folders, file extensions (.raw, .d, .mzML, ...) and letter case are ignored when matching.')
