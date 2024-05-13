
Input
=====

*OmicScope* offers five methods for integrating data into its pipeline, three of which rely on proteomic software for protein identification and quantitative analyses:


* 
  **Progenesis QI for Proteomics** *(\ ``Method = 'Progenesis'``\ )*\ : Progenesis QI for Proteomics (Waters Corporation) is a software that enables protein quantification and identification (via APEX3D and ProteinLynx Global Server) for experiments that use Data Independent Acquisition (DIA).

* 
  **PatternLab V** *(\ ``Method = 'PatternLab'``\ )*\ : PatternLab V is an integrated computational environment for analyzing shotgun proteomic data, and is considered one of the best options for quantitative proteomics using data-dependent acquisition, due to its high-confidence parameters for protein quantitation and identification.

* 
  **MaxQuant** *(\ ``Method = 'MaxQuant'``\ )*\ : MaxQuant is the most widely used software for quantitative proteomics, offering users a range of parameter options for quantitative analyses.

* 
  **DIA-NN** *(\ ``Method = 'DIA-NN'``\ )*\ : DIA-NN is a popular software that performs protein identification and quantification for DIA experiments, offering users a variety of parameter options for quantitative analyses.

* 
  **General** *(\ ``Method = 'General'``\ )*\ : To import data from other sources, this generic input method requires users to import an Excel workbook with 3 sheets: quantitative values (sheet1 = assay), protein features (sheet2 = rdata), and sample information (sheet3 = pdata).

* 
  **Snapshot** *(\ ``Method = 'Snapshot'``\ )*\ : The Snapshot method is the simplest approach to be used in the OmicScope workflow. It involves using a single, concise Excel sheet that contains essential information about proteins in a study, including fold change, p-value, and grouping.

For more information about the formatting of these import methods, see the appropriate sections below.

Import OmicScope
----------------

First, OmicScope package must be imported in the Python programming environment.

.. code-block:: python

   import omicscope as omics

.. code-block::

   OmicScope v 1.3.12 For help: Insert
   If you use  in published research, please cite: 'XXXXX'
   Reis-de-Oliveira G, Martins-de-Souza D. OmicScope: from quantitative proteomics to systems biology.




Progenesis QI for Proteomics
----------------------------

Progenesis exports protein quantitation data in a csv file containing information about samples, protein groups, and quantitative values. 

*OmicScope* imports Progenesis output and extracts the abundance levels of each protein (assay), the features of each protein (rdata), and features of each sample (pdata). *OmicScope*  can also accept Excel spreadsheets (with extensions .xls or .xlsx) that contain a **single sheet** for the Progenesis workflow, as many users may use Excel to visualize and handle data.

.. code-block:: python

   progenesis = omics.OmicScope('../../tests/data/proteins/progenesis.xls', Method='Progenesis')

.. code-block::

   User already performed statistical analysis
   OmicScope identifies: 697 deregulations



Since Progenesis exports certain information about sample groupings, *OmicScope* allows the user to input an Excel file containing all this information using pdata argument (for more information about pdata format, see below). Furthermore, users can filter identifications based on the minimum number of unique peptides by specifying the parameter ``UniquePeptides`` (recommended: ``UniquePeptides = 1``\ ).

.. code-block:: python

   progenesis_uniquepepfilt = omics.OmicScope('../../tests/data/proteins/progenesis.xls', Method='Progenesis', UniquePeptides=1)
   print('Original proteomics data: ' + str(len(progenesis.quant_data)) + '\n'+
         'Filtered proteomics data: ' + str(len(progenesis_uniquepepfilt.quant_data))
         )

.. code-block::

   User already performed statistical analysis
   OmicScope identifies: 582 deregulations
   Original proteomics data: 2179
   Filtered proteomics data: 1797



**Note**\ : Progenesis performs differential proteomics analyses based on preset groups, and *OmicScope* takes these statistical analyses into account. However, if the user has a specific experimental design, *OmicScope* Statistical Workflow can be used by renaming two columns in the original .csv file, as follows:


* "Anova (p)" → "Original Anova (p)"
* "q Value" → "Original q Value"

PatternLab
----------

PatternLab exports an Excel file with an .xlsx extension, which contains the same type of information as Progenesis, including assay, pdata, and rdata. However, this exported file does not include differential proteomics statistics. Therefore, *OmicScope* automatically performs statistical analyses for PatternLab data.

.. code-block:: python

   plv = omics.OmicScope('../../tests/data/proteins/patternlab.xlsx', Method='PatternLab')

.. code-block::

   Anova test was performed!
   OmicScope performed statistical analysis (Static workflow)



MaxQuant and DIA-NN
-------------------

While MaxQuant exports the **proteinGroups.txt** file, which contains a comprehensive description of the assay and rdata, DIA-NN exports the main output containing the same information. However, since in both cases the pdata is missing, these methods **requires** an additional Excel file for pdata. See the pdata section below for how to format this file.

.. code-block:: python

   maxquant = omics.OmicScope('../../tests/data/proteins/MQ.txt', Method='MaxQuant',
                              pdata='../../tests/data/proteins/MQ_pdata.xlsx')

.. code-block::

   Anova test was performed!
   OmicScope performed statistical analysis (Static workflow)



General
-------

The General workflow allows users to analyze data generated by other platforms, including Genomics and Transcriptomics. To do this, users need to organize an Excel file into three sheets: assay, rdata, and pdata. 


* **Assay:** Contains the abundance of N proteins (rows) from M samples (columns).
* **Rdata:** Includes N proteins (rows) with their respective features within each column.
* **Pdata:** Contains M samples (rows) with their respective characteristics, such as conditions as well as the organization of biological and technical replicates.

For more information about how to properly format and import each of these sheets, see the respective sections below.

.. code-block:: python

   general = omics.OmicScope('../../tests/data/proteins/general.xlsx', Method='General')

.. code-block::

   User already performed statistical analysis
   OmicScope identifies: 697 deregulations



Assay
^^^^^

The assay sheet should contain the abundance data for each protein/feature/transcript. The first row contains the sample names for each of the abundance values below.

.. code-block:: python

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
^^^^^

The rdata sheet needs to have at least two columns: 'Accession' and 'Description'.


#. **Accession:** An array of unique values that represent the proteins in the assay dataframe.
#. **Description:** The header from UniProt Fasta.

Optionally, user may add "gene_name" column for alternative names.

.. code-block:: python

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
^^^^^

Pdata contains a description of each sample analyzed in the workflow. Pdata must have at least the following 3 columns: 'Sample', 'Condition', and 'Biological'.


#. **Sample:** The name of each sample to be analysed, matching those in the first row of the Assay sheet.
#. **Condition:** Respective group for each sample. All technical and biological replicates belonging to an experimental condition should have the same identifier here.
#. **Biological:** Respective biological replicate for each sample. If two or more technical replicates were used for a single biological replicate, those replicates should have the same identifier here.

When performing longitudinal analysis, users must also include a ``TimeCourse`` column containing the day/hour/time/etc. associated with each sample.

See the example below for how to construct a pdata sheet. In this example, there are two groups being compared: COVID *vs.* CTRL. COVID contains 12 biological replicates, CTRL contains 7 biological replicates. All replicates were injected twice for two instrumental replicates. These replicates will be averaged and not considered individual samples for T-Test purposes.

.. code-block:: python

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


Snapshot
--------

The Snapshot method is an alternative option in OmicScope that allows for the analysis of multiple 'omics studies, importing pre-analyzed data from other platforms.

To use the Snapshot method, the user needs to upload a CSV or Excel file organized as follows:


#. First row: **ControlGroup: LIST_YOUR_CONTROL_HERE**
#. Second row: **Experimental: LIST_YOUR_EXPERIMENTAL_GROUPS_SEPARATED_BY_COMMAS**
#. Third row: A table header containing the following values: 'Accession', 'gene_name', 'log2(fc)', and either 'pvalue' or 'pAdjusted'. 
#. Subsequent rows: The molecular data to fill the columns listed in the third row.

It is important to note that Snapshot contains a comparatively limited amount of information, which means that not all plots and enrichment analyses will be available. Nevertheless, once the data is imported into OmicScope, it can still be exported as an .omics file and used in the Nebula module.

Additional Informations
-----------------------

Users can also define any of the following additional parameters that are in the OmicScope function to optimize their analysis.


#. 
   **ControlGroup** (default, ``ControlGroup = None``\ ): Users can define a control group to perform comparisons against a specific group. The name of this group has to be explicitly defined in the 'Conditions' column on the pdata table.

#. 
   **ExperimentalDesign** (default, ``ExperimentalDesign = 'static'``\ ) (options: 'static', 'longitudinal'): Comparisons among independent groups are called static experimental designs. However, if the experiment takes into account several time points of related samples, then it is performing a longitudinal experimental design. **Note:** in this case, the pdata table must present a 'TimeCourse' column.

#. 
   **pvalue** (default, ``pvalue = 'pAdjusted'``\ ) (options: 'pvalue', 'pAdjusted', 'pTukey'): Defines the type of statistics used to report differentially regulated proteins. The options are nominal p-value ('pvalue'), Benjamini-Hochberg adjusted p-value ('pAdjusted'), or Tukey post-hoc correction ('pTukey', only available for multiple group comparisons in static experiments).

#. 
   **PValue_cutoff** (default = ``PValue_cutoff = 0.05``\ ): Statistical cutoff to consider proteins differentially regulated.

#. 
   **normalization_method** (default, ``normalization_method = None``\ ): User can normalize data using OmicScope. Nowadays, OmicScope offers normalization using "average", "median", and "quantile" approaches.

#. 
   **imputation_method** (default, ``imputation_method = None``\ ): User can replace null values using three imputation methods: "mean", "knn", and "median".

#. 
   **FoldChange_cutoff** (default, ``FoldChange_cutoff = 0``\ ): Cutoff of the absolute abundance ratio to consider a protein to be differentially regulated. 0 indicates that p-values alone are sufficient to determine dysregulation.

#. 
   **logTransform** (default, ``logTransform = True``\ ): Usually, analysis software reports abundance in nominal values, requiring a log-transformation of the values to normalize abundance data.

#. 
   **ExcludeContaminants** (default, ``ExcludeContaminants = True``\ ): Recently, Frankenfield (2022) evaluated the most common contaminants found in proteomics workflows. By default, OmicScope removes them from analyses. If this is not desired, OmicScope can leave them in the final results with ExcludeContaminants=False.

#. 
   **degrees_of_freedom** (default, ``degrees_of_freedom = 2``\ ): For longitudinal experiments, users can optimize this parameter according to their study, choosing a greater degree of freedom to perform the subsequent statistical analyses. Note that ExperimentalDesign and pdata must still be appropriately configured.

#. 
   **independent_ttest** (default, ``independent_ttest = True``\ ): If running a t-test, the user can specify if data sampling was independent (True) or paired (False).
