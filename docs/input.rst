Input
=====

*OmicScope* offers ten methods for integrating data into its pipeline.
Seven of them read the output of proteomics software directly, one reads
transcriptomics tables, and two (General and Snapshot) accept data from
any source:

**Important:** Examples for each data format, together with a matching
pdata file, can be downloaded in the `OmicScope Web App
<https://omicscope.ib.unicamp.br/>`__ (sidebar, *Example files*).

.. list-table::
   :header-rows: 1
   :widths: 22 18 38 22

   * - Software
     - Method
     - File to export
     - pdata
   * - Progenesis QI for Proteomics
     - ``'Progenesis'``
     - Protein measurements (.csv; .xls/.xlsx also accepted)
     - optional
   * - PatternLab V
     - ``'PatternLab'``
     - .plp file from the XIC module
     - optional
   * - MaxQuant
     - ``'MaxQuant'``
     - proteinGroups.txt (LFQ, intensity, iBAQ or TMT reporter intensities)
     - **required**
   * - DIA-NN
     - ``'DIA-NN'``
     - report.pg_matrix.tsv (recommended) or report.tsv
     - **required**
   * - Spectronaut
     - ``'Spectronaut'``
     - Protein pivot report, or a long report with R.Condition
     - required for the pivot report
   * - FragPipe
     - ``'FragPipe'``
     - combined_protein.tsv (label-free) or TMT-Integrator
       abundance_protein_*.tsv (TMT)
     - **required**
   * - Proteome Discoverer
     - ``'ProteomeDiscoverer'``
     - Proteins table exported to Excel (label-free or TMT)
     - optional
   * - Transcriptomics
     - ``'RPKM'``
     - Tab-separated RPKM/FPKM/TPM table
     - **required**
   * - Any source
     - ``'General'``
     - Excel workbook with assay, rdata and pdata sheets
     - included
   * - Pre-analysed data
     - ``'Snapshot'``
     - Excel/csv with fold changes and p-values
     - not used

-  **Progenesis QI for Proteomics** (Waters) identifies and quantifies
   proteins in DIA experiments. Export the protein table with normalized
   abundances.
-  **PatternLab V** is an integrated environment for shotgun proteomics.
   Export the .plp file from the XIC section (filters can be applied
   before exporting). *This import method was developed in collaboration
   with the PatternLab developers.*
-  **MaxQuant**, **DIA-NN**, **Spectronaut**, **FragPipe** and
   **Proteome Discoverer** are covered in their own sections below.
-  **General** imports data from any platform (including other omics)
   organised in an Excel workbook with three sheets: assay, rdata and
   pdata.
-  **Snapshot** imports results that were already analysed elsewhere
   (fold change, p-value and groups).

When a method needs pdata, the ``Sample`` column must contain the sample
(run) names used in the quantitative file. OmicScope ignores folders,
file extensions (``.raw``, ``.d``, ``.mzML``, ``.wiff``, ...), letter
case and surrounding spaces when it matches the two, so
``D:\data\Sample_01.raw`` in the data file matches ``Sample_01`` in
pdata. If a pdata sample still cannot be found, the error message lists
the sample names available in the data file. Samples that are in the data
file but not in pdata are left out of the analysis, which is a simple way
to analyse a subset of runs.

For more information about the formatting of these import methods, see
the appropriate sections below.

Import OmicScope
----------------

First, OmicScope package must be imported in the Python programming
environment.

.. code:: python

    import omicscope as omics


.. parsed-literal::

    OmicScope v 1.4.0 For help: https://omicscope.readthedocs.io/en/latest/ or https://omicscope.ib.unicamp.brIf you use  in published research, please cite:
    'Reis-de-Oliveira, G., et al (2024). OmicScope unravels systems-level insights from quantitative proteomics data 
    
    

Progenesis QI for Proteomics
----------------------------

Progenesis exports protein quantitation data in a CSV file containing
information about samples, protein groups, and quantitative values.

*OmicScope* imports Progenesis output and extracts the abundance levels
of each protein (assay), the features of each protein (rdata), and
features of each sample (pdata). *OmicScope* can also accept Excel
spreadsheets (with extensions .xls or .xlsx) that contain a **single
sheet** for the Progenesis workflow, as many users may use Excel to
visualize and handle data.

.. code:: python

    progenesis = omics.OmicScope('../../tests/data/proteins/progenesis.xls', Method='Progenesis')


.. parsed-literal::

    User already performed statistical analysis
    OmicScope identifies: 697 deregulations
    

**Only for OmicScope Package (not available in OmicScope App)**

Since Progenesis exports certain information about sample groupings,
*OmicScope* allows the user to input an Excel file containing all this
information using the pdata argument (for more information about pdata
format, see below). Furthermore, users can filter identifications based
on the minimum number of unique peptides by specifying the parameter
``UniquePeptides`` (recommended: ``UniquePeptides = 1``).

.. code:: python

    progenesis_uniquepepfilt = omics.OmicScope('../../tests/data/proteins/progenesis.xls', Method='Progenesis', UniquePeptides=1)
    print('Original proteomics data: ' + str(len(progenesis.quant_data)) + '\n'+
          'Filtered proteomics data: ' + str(len(progenesis_uniquepepfilt.quant_data))
          )


.. parsed-literal::

    User already performed statistical analysis
    OmicScope identifies: 582 deregulations
    Original proteomics data: 2179
    Filtered proteomics data: 1797
    

**IMPORTANT**: Progenesis performs differential proteomics analyses
based on preset groups, and *OmicScope* takes these statistical analyses
into account. However, if the user has a specific experimental design,
*OmicScope* Statistical Workflow can be used by renaming two columns in
the original .csv file, as follows:

-  “Anova (p)” → “Original Anova (p)”
-  “q Value” → “Original q Value”

PatternLab
----------

PatternLab exports an plp file, which contains
the same type of information as Progenesis, including assay, pdata, and
rdata. However, this exported file does not include differential
proteomics statistics. Therefore, *OmicScope* automatically performs
statistical analyses for PatternLab data.

For users aiming to use Statistical Workflows from T-Fold workflow in PatternLab,
we advise using "Snapshot" import method.

.. code:: python

    plv = omics.OmicScope('../../tests/data/proteins/patternlab.plp', Method='PatternLab')

MaxQuant
--------

**File:** ``combined/txt/proteinGroups.txt``. **pdata:** required.

proteinGroups.txt has protein features and abundances but not the sample
conditions, so a pdata file is required. See the `pdata section
<#pdata>`__ below for its format. The ``Sample`` column holds the
experiment names typed in MaxQuant, i.e. the text after the quantification
prefix: ``CoV1`` for the column ``LFQ intensity CoV1``.

OmicScope removes reverse hits, potential contaminants, proteins only
identified by site and proteins with Andromeda score <= 0. The first gene
of each protein group is used as ``gene_name``; when ``Gene names`` is
empty, the gene is taken from the FASTA header (``GN=``) or, failing that,
the accession is used.

Choose the abundance columns with ``quant_strategy``:

-  ``'LFQ intensity'`` (default)
-  ``'Intensity'``
-  ``'iBAQ'``
-  ``'Reporter intensity corrected'`` (TMT/iTRAQ)
-  ``'Reporter intensity'`` (TMT/iTRAQ, not corrected for isotopic
   impurities)

If the requested columns are missing, OmicScope uses ``Reporter intensity
corrected`` (isobaric labelling) or ``Intensity`` instead. In TMT
experiments, sample names are ``<channel> <experiment>`` (for example
``1 Batch1`` for the column ``Reporter intensity corrected 1 Batch1``),
or only ``<channel>`` when there is a single experiment.

.. code:: python

    maxquant = omics.OmicScope('proteinGroups.txt', Method='MaxQuant',
                               pdata='MaxQuant_pdata.xlsx')
    # TMT
    tmt = omics.OmicScope('proteinGroups.txt', Method='MaxQuant',
                          pdata='TMT_pdata.xlsx',
                          quant_strategy='Reporter intensity corrected')

**Troubleshooting:**

-  Check that the quantitative columns contain valid values (MaxQuant
   writes zeros when LFQ could not be computed).
-  OmicScope expects the column labels of MaxQuant 1.6 and later
   (``Majority protein IDs``, ``Fasta headers``, ``Gene names``).

DIA-NN
------

**File:** ``report.pg_matrix.tsv`` (recommended) or ``report.tsv``.
**pdata:** required.

DIA-NN does not store the sample conditions, so a pdata file is required.
OmicScope accepts:

-  **report.pg_matrix.tsv** (DIA-NN 1.8, 1.9 and 2.x): one row per protein
   group and one column per run. It is small and loads fast. The values
   are the MaxLFQ protein-group quantities, already filtered by DIA-NN at
   1% FDR.
-  **report.tsv** (main report of DIA-NN 1.8 and 1.9): one row per
   precursor and run. OmicScope uses ``PG.MaxLFQ`` and keeps rows with
   ``Q.Value`` <= 0.01 and ``PG.Q.Value`` <= 0.01 (change the latter with
   ``pg_qvalue``). These files can be larger than 1 GB.
-  **report.parquet** (main report of DIA-NN 2.x) is read only when the
   ``pyarrow`` package is installed. Otherwise use
   ``report.pg_matrix.tsv``, which DIA-NN 2.x still writes.

Sample names in pdata are the run names: the raw file names, as in the
``Run`` column of report.tsv or the column headers of the matrix. Folders
and extensions are ignored, so ``D:\Raw\HeLa_01.raw`` (matrix header),
``HeLa_01`` (``Run`` column) and ``HeLa_01.raw`` all match the pdata
sample ``HeLa_01``. Decoys and protein groups tagged as contaminants
(``Cont_``, ``CON__``, ``contam_`` prefixes) are removed.

.. code:: python

    diann = omics.OmicScope('report.pg_matrix.tsv', Method='DIA-NN',
                            pdata='DIANN_pdata.xlsx')

Spectronaut
-----------

**File:** a protein-level report exported from the *Report* perspective.
**pdata:** required for pivot reports only.

OmicScope reads the protein-group quantity (``PG.Quantity``) from two
layouts:

-  **Pivot report** (a pivot schema with protein-group columns): one
   row per protein group and one ``[n] <run>.PG.Quantity`` column per run.
   It does not contain the conditions, so pdata is required. Sample names
   are the run names, i.e. the text between ``[n]`` and ``.PG.Quantity``
   (``HeLa_01`` for ``[1] HeLa_01.PG.Quantity``).
-  **Long report** (for example the *BGS Factory Report (Normal)*
   schema, or a custom schema): one row per protein group (or precursor)
   and run. Include at least ``R.Condition``, ``R.FileName``,
   ``R.Replicate``, ``PG.ProteinGroups``, ``PG.Genes`` and
   ``PG.Quantity``. Conditions and replicates are read from
   ``R.Condition`` and ``R.Replicate``, so pdata is optional; pass one to
   change the design.

Values written as ``Filtered`` or ``NaN`` are treated as missing, and
decimal commas are accepted. Decoys (``PG.IsDecoy``/``EG.IsDecoy``) and
protein groups tagged as contaminants are removed.

.. code:: python

    # Long report: conditions taken from R.Condition
    spectronaut = omics.OmicScope('Spectronaut_Report.tsv', Method='Spectronaut',
                                  ControlGroup='Control')
    # Pivot report
    spectronaut = omics.OmicScope('Spectronaut_Report_Pivot.tsv', Method='Spectronaut',
                                  pdata='Spectronaut_pdata.xlsx')

FragPipe
--------

**File:** ``combined_protein.tsv`` (label-free) or
``tmt-report/abundance_protein_MD.tsv`` (TMT). **pdata:** required.

-  **Label-free (IonQuant):** OmicScope uses the ``<sample> MaxLFQ
   Intensity`` columns, or ``<sample> Intensity`` when MaxLFQ was not
   computed. Sample names are the experiment names set in the FragPipe
   workflow (``IDHmut_1`` for ``IDHmut_1 MaxLFQ Intensity``).
-  **TMT (TMT-Integrator):** use the protein report
   ``abundance_protein_MD.tsv`` (``_None`` and ``_GN`` normalisation
   variants also work). Sample names are the names given in the
   annotation file, i.e. the columns after ``ReferenceIntensity``.
   TMT-Integrator reports log2 abundances; OmicScope converts them back to
   linear values so that the default ``logTransform=True`` gives the
   original log2 values.

Decoys (``rev_``) and contaminants (``contam_``) are removed.

.. code:: python

    fragpipe = omics.OmicScope('combined_protein.tsv', Method='FragPipe',
                               pdata='FragPipe_pdata.xlsx')

Proteome Discoverer
-------------------

**File:** the *Proteins* table exported to Excel (File > Export > To
Microsoft Excel). **pdata:** optional.

Proteome Discoverer (PD) writes samples, protein features and abundances
in a single sheet. Export at least ``Accession``, ``Description`` and the
``Abundances (Normalized)`` (or ``Abundance``) columns. To import the PD
statistics instead of running OmicScope's tests, also export ``Abundance
Ratio P-Value`` and ``Abundance Ratio Adj. P-Value``.

Conditions are read from the abundance headers. In ``Abundances
(Normalized): F17: Sample, SCZ_5, CBD (SCZ)`` the last study factor
(``CBD (SCZ)``) is the condition and the factor before it (``SCZ_5``) is
the biological replicate, so files sharing it are averaged as technical
replicates. With a single study factor (common in TMT, e.g. ``F1: 126,
Sample, CTRL``), each file/channel is one biological replicate. Only
master proteins are kept, and proteins flagged in the ``Contaminant``
column are removed.

.. code:: python

    pdisc = omics.OmicScope('ProteomeDiscoverer_proteins.xlsx', Method='ProteomeDiscoverer')

RPKM (transcriptomics)
----------------------

**File:** tab-separated table with ``gene_name``, ``Accession`` and one
column per sample holding RPKM, FPKM or TPM values. **pdata:** required.

Values are transformed with log2(value + 1) before the analysis. Sample
names in pdata are the column names of the table.

.. code:: python

    rpkm = omics.OmicScope('RPKM.tsv', Method='RPKM', pdata='RPKM_pdata.xlsx')

Not supported directly
----------------------

Tools without a dedicated method can be imported with ``General``: put
protein abundances in the first sheet (one column per sample), protein
information in the second (at least ``Accession`` and ``gene_name``) and
pdata in the third. This applies, for example, to PEAKS (protein table
with ``Area`` columns), AlphaDIA (``pg.matrix.tsv``), quantms/MSstats
protein-level results and Perseus matrices. Spectral-count exports such as
Scaffold samples reports are not recommended for intensity-based
statistics.

General
-------

The General workflow allows users to analyze data generated by other
platforms, including Genomics and Transcriptomics. To do this, users
need to organize an Excel file into three sheets: assay, rdata, and
pdata.

-  **Assay:** Contains the abundance of N proteins (rows) from M samples
   (columns).
-  **Rdata:** Includes N proteins (rows) with their respective features
   within each column.
-  **Pdata:** Contains M samples (rows) with their respective
   characteristics, such as conditions, as well as the organization of
   biological and technical replicates.

For more information about how to properly format and import each of
these sheets, see the respective sections below.

The pdata sheet can also be supplied as a separate file with
``pdata=``. The workbook then only needs the assay and rdata sheets. For
a time-course example in both layouts, see the *Longitudinal example*
section of the pdata page.

.. code:: python

    general = omics.OmicScope('../../tests/data/proteins/general.xlsx', Method='General')

Assay
~~~~~

The assay sheet should contain the abundance data for each
protein/feature/transcript. The first row contains the sample names for
each of the abundance values below.

.. code:: python

    import pandas as pd
    
    assay = pd.read_excel('../../tests/data/proteins/general.xlsx', sheet_name=0)
    # Slicing example to facilitate visualization
    assay.head().iloc[:,0:5]




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>VCC_HB_1_1_2020</th>
          <th>VCC_HB_1_2</th>
          <th>VCC_HB_2_1</th>
          <th>VCC_HB_2_1_2</th>
          <th>VCC_HB_3_1</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>2.938847e+04</td>
          <td>3.110927e+04</td>
          <td>2.521807e+04</td>
          <td>3.090703e+04</td>
          <td>2.383499e+04</td>
        </tr>
        <tr>
          <th>1</th>
          <td>7.081308e+04</td>
          <td>6.446946e+04</td>
          <td>5.825493e+04</td>
          <td>5.931610e+04</td>
          <td>6.309095e+04</td>
        </tr>
        <tr>
          <th>2</th>
          <td>1.007536e+05</td>
          <td>1.011999e+05</td>
          <td>7.301329e+04</td>
          <td>7.349391e+04</td>
          <td>9.766835e+04</td>
        </tr>
        <tr>
          <th>3</th>
          <td>2.588031e+04</td>
          <td>3.769105e+04</td>
          <td>2.992691e+04</td>
          <td>3.460095e+04</td>
          <td>2.596320e+04</td>
        </tr>
        <tr>
          <th>4</th>
          <td>1.019192e+06</td>
          <td>1.109406e+06</td>
          <td>1.060396e+06</td>
          <td>1.078239e+06</td>
          <td>1.003426e+06</td>
        </tr>
      </tbody>
    </table>
    </div>



rdata
~~~~~

The rdata sheet needs to have at least two columns: ‘Accession’ and
‘Description’.

1. **Accession:** An array of unique values that represent the proteins
   in the assay dataframe.
2. **Description:** The header from UniProt Fasta.

Optionally, user may add “gene_name” column for alternative names.

.. code:: python

    rdata = pd.read_excel('../../tests/data/proteins/general.xlsx', sheet_name=1)
    rdata.head(3)




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>Accession</th>
          <th>Peptide count</th>
          <th>Unique peptides</th>
          <th>Confidence score</th>
          <th>Anova (p)</th>
          <th>q Value</th>
          <th>Max fold change</th>
          <th>Power</th>
          <th>Highest mean condition</th>
          <th>Lowest mean condition</th>
          <th>Description</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>P0DJI8</td>
          <td>1</td>
          <td>1</td>
          <td>6.8809</td>
          <td>0.000000e+00</td>
          <td>0.000000</td>
          <td>2.192654</td>
          <td>1.000000</td>
          <td>COVID</td>
          <td>CTRL</td>
          <td>Serum amyloid A-1 protein OS=Homo sapiens OX=9...</td>
        </tr>
        <tr>
          <th>1</th>
          <td>P63313</td>
          <td>2</td>
          <td>0</td>
          <td>24.1939</td>
          <td>0.000000e+00</td>
          <td>0.000000</td>
          <td>3.823799</td>
          <td>1.000000</td>
          <td>COVID</td>
          <td>CTRL</td>
          <td>Thymosin beta-10 OS=Homo sapiens OX=9606 GN=TM...</td>
        </tr>
        <tr>
          <th>2</th>
          <td>P03886</td>
          <td>3</td>
          <td>0</td>
          <td>24.0213</td>
          <td>1.299387e-07</td>
          <td>0.000041</td>
          <td>1.386199</td>
          <td>0.999998</td>
          <td>CTRL</td>
          <td>COVID</td>
          <td>NADH-ubiquinone oxidoreductase chain 1 OS=Homo...</td>
        </tr>
      </tbody>
    </table>
    </div>



pdata
~~~~~

Pdata contains a description of each sample analyzed in the workflow.
Pdata must have at least the following 3 columns: ‘Sample’, ‘Condition’,
and ‘Biological’.

1. **Sample:** The name of each sample to be analyzed, matching those in
   the first row of the Assay sheet.
2. **Condition:** Respective group for each sample. All technical and
   biological replicates belonging to an experimental condition should
   have the same identifier here.
3. **Biological:** Respective biological replicate for each sample. If
   two or more technical replicates were used for a single biological
   replicate, those replicates should have the same identifier here.

When performing longitudinal analysis, users must also include a
``TimeCourse`` column containing the day/hour/time/etc. associated with
each sample.

See the example below for how to construct a pdata sheet. In this
example, there are two groups being compared: COVID *vs.* CTRL. COVID
contains 12 biological replicates, CTRL contains 7 biological
replicates. All replicates were injected twice for two instrumental
replicates. These replicates will be averaged and not considered
individual samples for T-Test purposes.

.. code:: python

    pdata = pd.read_excel('../../tests/data/proteins/general.xlsx', sheet_name=2)
    pdata




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>Sample</th>
          <th>Condition</th>
          <th>Biological</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>VCC_HB_1_1_2020</td>
          <td>COVID</td>
          <td>1</td>
        </tr>
        <tr>
          <th>1</th>
          <td>VCC_HB_1_2</td>
          <td>COVID</td>
          <td>1</td>
        </tr>
        <tr>
          <th>2</th>
          <td>VCC_HB_2_1</td>
          <td>COVID</td>
          <td>2</td>
        </tr>
        <tr>
          <th>3</th>
          <td>VCC_HB_2_1_2</td>
          <td>COVID</td>
          <td>2</td>
        </tr>
        <tr>
          <th>4</th>
          <td>VCC_HB_3_1</td>
          <td>COVID</td>
          <td>3</td>
        </tr>
        <tr>
          <th>5</th>
          <td>VCC_HB_3_1_2</td>
          <td>COVID</td>
          <td>3</td>
        </tr>
        <tr>
          <th>6</th>
          <td>VCC_HB_4_1</td>
          <td>COVID</td>
          <td>4</td>
        </tr>
        <tr>
          <th>7</th>
          <td>VCC_HB_4_1_2</td>
          <td>COVID</td>
          <td>4</td>
        </tr>
        <tr>
          <th>8</th>
          <td>VCC_HB_5_1</td>
          <td>COVID</td>
          <td>5</td>
        </tr>
        <tr>
          <th>9</th>
          <td>VCC_HB_5_1_2</td>
          <td>COVID</td>
          <td>5</td>
        </tr>
        <tr>
          <th>10</th>
          <td>VCC_HB_6_1</td>
          <td>COVID</td>
          <td>6</td>
        </tr>
        <tr>
          <th>11</th>
          <td>VCC_HB_6_1_2</td>
          <td>COVID</td>
          <td>6</td>
        </tr>
        <tr>
          <th>12</th>
          <td>VCC_HB_7_1</td>
          <td>COVID</td>
          <td>7</td>
        </tr>
        <tr>
          <th>13</th>
          <td>VCC_HB_7_1_2</td>
          <td>COVID</td>
          <td>7</td>
        </tr>
        <tr>
          <th>14</th>
          <td>VCC_HB_8_1</td>
          <td>COVID</td>
          <td>8</td>
        </tr>
        <tr>
          <th>15</th>
          <td>VCC_HB_8_1_2</td>
          <td>COVID</td>
          <td>8</td>
        </tr>
        <tr>
          <th>16</th>
          <td>VCC_HB_9_1</td>
          <td>COVID</td>
          <td>9</td>
        </tr>
        <tr>
          <th>17</th>
          <td>VCC_HB_9_1_2</td>
          <td>COVID</td>
          <td>9</td>
        </tr>
        <tr>
          <th>18</th>
          <td>VCC_HB_10_1</td>
          <td>COVID</td>
          <td>10</td>
        </tr>
        <tr>
          <th>19</th>
          <td>VCC_HB_10_1_2_</td>
          <td>COVID</td>
          <td>10</td>
        </tr>
        <tr>
          <th>20</th>
          <td>VCC_HB_11_1</td>
          <td>COVID</td>
          <td>11</td>
        </tr>
        <tr>
          <th>21</th>
          <td>VCC_HB_11_1_2_</td>
          <td>COVID</td>
          <td>11</td>
        </tr>
        <tr>
          <th>22</th>
          <td>VCC_HB_12_1</td>
          <td>COVID</td>
          <td>12</td>
        </tr>
        <tr>
          <th>23</th>
          <td>VCC_HB_12_1_2_</td>
          <td>COVID</td>
          <td>12</td>
        </tr>
        <tr>
          <th>24</th>
          <td>VCC_HB_A_1</td>
          <td>CTRL</td>
          <td>1</td>
        </tr>
        <tr>
          <th>25</th>
          <td>VCC_HB_A_1_2</td>
          <td>CTRL</td>
          <td>1</td>
        </tr>
        <tr>
          <th>26</th>
          <td>VCC_HB_B_1</td>
          <td>CTRL</td>
          <td>2</td>
        </tr>
        <tr>
          <th>27</th>
          <td>VCC_HB_B_1_2</td>
          <td>CTRL</td>
          <td>2</td>
        </tr>
        <tr>
          <th>28</th>
          <td>VCC_HB_C_1</td>
          <td>CTRL</td>
          <td>3</td>
        </tr>
        <tr>
          <th>29</th>
          <td>VCC_HB_C_1_2</td>
          <td>CTRL</td>
          <td>3</td>
        </tr>
        <tr>
          <th>30</th>
          <td>VCC_HB_D_1</td>
          <td>CTRL</td>
          <td>4</td>
        </tr>
        <tr>
          <th>31</th>
          <td>VCC_HB_D_1_2</td>
          <td>CTRL</td>
          <td>4</td>
        </tr>
        <tr>
          <th>32</th>
          <td>VCC_HB_E_1</td>
          <td>CTRL</td>
          <td>5</td>
        </tr>
        <tr>
          <th>33</th>
          <td>VCC_HB_E_1_2</td>
          <td>CTRL</td>
          <td>5</td>
        </tr>
        <tr>
          <th>34</th>
          <td>VCC_HB_F_1</td>
          <td>CTRL</td>
          <td>6</td>
        </tr>
        <tr>
          <th>35</th>
          <td>VCC_HB_F_1_2</td>
          <td>CTRL</td>
          <td>6</td>
        </tr>
        <tr>
          <th>36</th>
          <td>VCC_HB_G_1</td>
          <td>CTRL</td>
          <td>7</td>
        </tr>
        <tr>
          <th>37</th>
          <td>VCC_HB_G_1_2</td>
          <td>CTRL</td>
          <td>7</td>
        </tr>
      </tbody>
    </table>
    </div>



For detailed instructions on constructing pdata and integrating it into
your experimental design, please refer to the page titled `How to Make
Pdata <https://omicscope.readthedocs.io/en/latest/pdata.html>`__.

Snapshot
--------

The Snapshot method is an alternative option in OmicScope for analyzing
multiple ’omics studies by importing pre-analyzed data from other
platforms.

To use the Snapshot method, the user needs to upload a CSV or Excel file
organized as follows:

1. First row: **ControlGroup: LIST_YOUR_CONTROL_HERE**
2. Second row: **Experimental:
   LIST_YOUR_EXPERIMENTAL_GROUPS_SEPARATED_BY_COMMAS**
3. Third row: A table header containing the following values:
   ‘Accession’, ‘gene_name’, ‘log2(fc)’, and either ‘pvalue’ or
   ‘pAdjusted’.
4. Subsequent rows: The molecular data to fill the columns listed in the
   third row.

It is important to note that Snapshot contains a comparatively limited
amount of information, which means that not all plots and enrichment
analyses will be available. Nevertheless, once the data is imported into
OmicScope, it can still be exported as an .omics file and used in the
Nebula module.

Additional Informations
-----------------------

Users can also define any of the following additional parameters that
are in the OmicScope function to optimize their analysis.

1.  **ControlGroup** (default, ``ControlGroup = None``): Users can
    define a control group to perform comparisons against a specific
    group. The name of this group has to be explicitly defined in the
    ‘Conditions’ column on the pdata table.
2.  **ExperimentalDesign** (default, ``ExperimentalDesign = 'static'``)
    (options: ‘static’, ‘longitudinal’): Comparisons among independent
    groups are called static experimental designs. However, if the
    experiment takes into account several time points of related
    samples, then it is performing a longitudinal experimental design.
    **Note:** in this case, the pdata table must present a ‘TimeCourse’
    column.
3.  **pvalue** (default, ``pvalue = 'pAdjusted'``) (options: ‘pvalue’,
    ‘pAdjusted’, ‘pTukey’): Defines the type of statistics used to
    report differentially regulated proteins. The options are nominal
    p-value (‘pvalue’), Benjamini-Hochberg adjusted p-value
    (‘pAdjusted’), or Tukey post-hoc correction (‘pTukey’, only
    available for multiple group comparisons in static experiments).
4.  **PValue_cutoff** (default = ``PValue_cutoff = 0.05``): Statistical
    cutoff to consider proteins differentially regulated.
5.  **normalization_method** (default =
    ``normalization_method = None``): Certain data may require a
    normalization preprocessing step. OmicScope offers three methods of
    normalization: ‘average’, ‘median’, ‘quantile’. Defaults to None.
6.  **imputation_method** (default = ``imputation_method = None``): Some
    data may require data imputation to handle null values as a
    preprocessing step. OmicScope provides three methods of data
    imputation: ‘mean’, ‘median’, ‘knn’. Defaults to None.
7.  **FoldChange_cutoff** (default, ``FoldChange_cutoff = 0``): Cutoff
    of the absolute abundance ratio to consider a protein to be
    differentially regulated. 0 indicates that p-values alone are
    sufficient to determine dysregulation.
8.  **logTransform** (default, ``logTransform = True``): Usually,
    analysis software reports abundance in nominal values, requiring a
    log-transformation of the values to normalize abundance data. If
    users performed transformation before the OmicScope workflow, set
    logTransformed=True.
9.  **ExcludeContaminants** (default, ``ExcludeContaminants = True``):
    Recently, Frankenfield (2022) evaluated the most common contaminants
    found in proteomics workflows. By default, OmicScope removes them
    from analyses. If this is not desired, OmicScope can leave them in
    the final results with ExcludeContaminants=False.
10. **degrees_of_freedom** (default, ``degrees_of_freedom = 2``): For
    longitudinal experiments, users can optimize this parameter
    according to their study, choosing a greater degree of freedom to
    perform the subsequent statistical analyses. Note that
    ExperimentalDesign and pdata must still be appropriately configured.
11. **independent_ttest** (default, ``independent_ttest = True``): If
    running a t-test, the user can specify if data sampling was
    independent (True) or paired (False).


