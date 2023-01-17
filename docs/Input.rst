=====
Input
=====

*OmicScope* have four methods to input data into its pipeline, of which three are based on proteomic software that performs protein identification and quantitation. The methods are:


* 
  **Progenesis QI for proteomics** *(Method = 'Progenesis')*\ : while associated with ProteinLynx Global Server (PLGS), Progenesis is a software that allows both proteins quantitation and identification of experiments that performed Data Independente Acquisition (DIA).  

* 
  **PatternLab V** *(Method = 'PatternLab')*\ : PatternLab V is an integrated computational environment for analyzing shotgun proteomic data, being considered one of the best options to perform quantitative proteomics from data-dependent acquisition, since its high confidence parameters for protein quantitation and identification.

* 
  **MaxQuant** *(Method = 'MaxQuant')*\ : MaxQuant is the most used software for quantitative proteomics, offering users several parameter options to run quantitative analysis.

* 
  **General** *(Method = 'General')*\ : A generic input method, in which users must use an Excel file to specify the quantitative values (sheet1 = assay), protein features (sheet2 = rdata), and sample information (sheet3 = pdata).

Import OmicScope
----------------

.. code-block:: python

   import omicscope as omics

.. code-block::

   OmicScope v 1.0.0 For help: Insert
   If you use  in published research, please cite: 'lallalala'
   Reis-de-Oliveira G, Martins-de-Souza D. OmicScope: an Comprehensive Python library for Systems Biology Visualization.




Progenesis QI for Proteomics
----------------------------

Progenesis exports a .csv extension file, which contains information regarding samples, proteins and quantitative values. For the Progenesis workflow, *OmicScope* imports this files and extract abundance levels of each protein (assay), the features of each protein (rdata), and features of each sample (pdata). 

**OBS:** Since Progenesis performs differential proteomics, OmicScope take into account the statistical analysis performed by Progenesis; however, if user have an specific experimental design and needs OmicScope Statistical Workflow, user have to rewrite in the original .csv file:


* "Anova (p)" --> "Original Anova (p)" 
* "q Value" --> "Original q Value"

.. code-block:: python

   progenesis = omics.OmicScope('../tests/data/proteins/progenesis.csv', Method = 'Progenesis')

.. code-block::

   User already performed statistical analysis



Since Progenesis do not export all information regarding samples, OmicScope allows user to input an Excel file containing all information regarding samples. This file will be used to categorize samples, perform statistical analysis, and generate figures. Additionally, users can filter data based on the minimum number of unique peptides adding to parameter function 'UniquePeptides' (suggestion: UniquePeptides = 1)

.. code-block:: python

   progenesis = omics.OmicScope('../tests/data/proteins/progenesis.csv', Method = 'Progenesis',
    pdata = '../tests/data/proteins/progenesis_pdata.xls', UniquePeptides = 1)

.. code-block::

   User already performed statistical analysis



PatternLab
----------

PatternLab exports an Excel file (.xlsx extension), containing in general the same information as Progenesis: assay, pdata and rdata. However, in this export there is no differential proteomics statistics, then OmicScope will automatically perform statistical analysis. Despite users can filter proteins manually excel file, OmicScope **requires** 'filtering_method' parameter to filter and quantify the identified proteins. OmicScope filters proteins based on 1) percentage of times that an protein was identified among samples (insert an integer number); or 2) exclude proteins that was identified at least in 2 samples per group (insert 'minimum').

.. code-block:: python

   plv = omics.OmicScope('../tests/data/proteins/patternlab.xlsx', Method = 'PatternLab', filtering_method = 70)

.. code-block::

   C:\Users/Guilherme/omicscope/omicscope/src\omicscope\Omicscope.py:134: FutureWarning: In a future version, the Index constructor will not infer numeric dtypes when passed object-dtype sequences (matching Series behavior)
     expression = expression.set_index(pdata).T


   Anova test was performed!
   OmicScope performed statistical analysis (Static workflow)



MaxQuant
--------

MaxQuant exports a **proteinGroups** file (.txt extension), containing a comprehensive description of assay and rdata. Due to missing pdata, MaxQuant workflow is **requires** 'filtering_method' parameter and an Excel file for pdata. OmicScope filters proteins based on 1) percentage of times that an protein was identified among samples (insert an integer number); or 2) exclude proteins that was identified at least in 2 samples per group (insert 'minimum').

.. code-block:: python

   maxquant = omics.OmicScope('../tests/data/proteins/MQ.txt', Method='MaxQuant',
               pdata='../tests/data/proteins/MQ_pdata.xlsx', filtering_method=70)

.. code-block::

   Anova test was performed!
   OmicScope performed statistical analysis (Static workflow)



General
-------

General workflow allows user to analyse data generated in other plataforms, even for Transcriptomics and Metabolomics. For that, users have to organize an Excel file in three (3) sheets containing, respectivelly, assay, rdata and pdata. 

**Assay** contains abundance of the N proteins (rows) among M samples (columns); **Rdata** has N proteins (rows) with their respective features among columns; **Pdata** has M samples (rows) with their respective characteristics (such as conditions, biological and technical replicates).

The following sections shows examples of how each sheet must be described.

.. code-block:: python

   general = omics.OmicScope('../tests/data/proteins/general.xls', Method='General')

.. code-block::

   Independent T-test was carried out!
   OmicScope performed statistical analysis (Static workflow)



Assay
^^^^^

.. code-block:: python

   import pandas as pd

   assay = pd.read_excel('../tests/data/proteins/general.xls', sheet_name=0)
   assay


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
         <th>VCC_KO_1_VINO</th>
         <th>VCC_KO_1_VINO_2</th>
         <th>VCC_KO_1_VINO_29102021</th>
         <th>VCC_KO_1_VINO_29102021_3</th>
         <th>VCC_KO_2_VINO</th>
         <th>VCC_KO_2_VINO_2</th>
         <th>VCC_KO_2_VINO_29102021</th>
         <th>VCC_KO_2_VINO_29102021_3</th>
         <th>VCC_KO_3_VINO</th>
         <th>VCC_KO_3_VINO_2</th>
         <th>...</th>
         <th>VCC_WT_2_VIN_29102021</th>
         <th>VCC_WT_2_VIN_29102021_2</th>
         <th>VCC_WT_3_VIN</th>
         <th>VCC_WT_3_VIN_2</th>
         <th>VCC_WT_3_VIN_29102021</th>
         <th>VCC_WT_3_VIN_29102021_2</th>
         <th>VCC_WT_4_VIN</th>
         <th>VCC_WT_4_VIN_2</th>
         <th>VCC_WT_4_VIN_29102021</th>
         <th>VCC_WT_4_VIN_29102021_2</th>
       </tr>
     </thead>
     <tbody>
       <tr>
         <th>0</th>
         <td>61282.526104</td>
         <td>58475.057832</td>
         <td>66491.864803</td>
         <td>63965.456771</td>
         <td>58599.602771</td>
         <td>58349.651075</td>
         <td>61126.678243</td>
         <td>61396.041785</td>
         <td>55983.435295</td>
         <td>55382.566170</td>
         <td>...</td>
         <td>82171.713393</td>
         <td>86964.333856</td>
         <td>83896.220644</td>
         <td>85960.705463</td>
         <td>123508.762577</td>
         <td>77645.954774</td>
         <td>83303.856481</td>
         <td>87632.085234</td>
         <td>78080.558618</td>
         <td>81497.447186</td>
       </tr>
       <tr>
         <th>1</th>
         <td>48284.094432</td>
         <td>51659.072375</td>
         <td>48700.892150</td>
         <td>55211.947643</td>
         <td>51033.426081</td>
         <td>50100.916082</td>
         <td>54566.724267</td>
         <td>50468.832724</td>
         <td>54797.997214</td>
         <td>52039.446331</td>
         <td>...</td>
         <td>58684.503206</td>
         <td>71913.438722</td>
         <td>71047.636656</td>
         <td>71125.976724</td>
         <td>53174.444736</td>
         <td>79038.061177</td>
         <td>67214.986877</td>
         <td>68608.124964</td>
         <td>65715.209981</td>
         <td>75314.101558</td>
       </tr>
       <tr>
         <th>2</th>
         <td>8275.498103</td>
         <td>7672.835670</td>
         <td>7676.683705</td>
         <td>7388.702687</td>
         <td>8971.608574</td>
         <td>8993.363424</td>
         <td>8689.472709</td>
         <td>9342.557740</td>
         <td>8261.663352</td>
         <td>7056.970146</td>
         <td>...</td>
         <td>33309.128490</td>
         <td>17392.234792</td>
         <td>10650.392858</td>
         <td>10640.789093</td>
         <td>14516.837540</td>
         <td>12384.828169</td>
         <td>20016.681999</td>
         <td>18983.880260</td>
         <td>19210.197630</td>
         <td>16118.917424</td>
       </tr>
       <tr>
         <th>3</th>
         <td>283603.747996</td>
         <td>275358.163322</td>
         <td>264519.003841</td>
         <td>322882.142746</td>
         <td>200863.590415</td>
         <td>222174.322464</td>
         <td>200538.991041</td>
         <td>266430.806302</td>
         <td>201782.520396</td>
         <td>174000.923670</td>
         <td>...</td>
         <td>898091.979181</td>
         <td>686655.971644</td>
         <td>420550.143562</td>
         <td>401333.316279</td>
         <td>552460.098385</td>
         <td>529880.936082</td>
         <td>438354.668416</td>
         <td>419538.761093</td>
         <td>487150.346242</td>
         <td>328164.625834</td>
       </tr>
       <tr>
         <th>4</th>
         <td>87324.461931</td>
         <td>93193.890073</td>
         <td>87119.771902</td>
         <td>92960.354306</td>
         <td>63819.952903</td>
         <td>71969.767523</td>
         <td>66863.673529</td>
         <td>67127.229702</td>
         <td>84533.473807</td>
         <td>90097.134209</td>
         <td>...</td>
         <td>135941.135022</td>
         <td>88903.637078</td>
         <td>150594.063275</td>
         <td>128800.719644</td>
         <td>105642.253308</td>
         <td>127004.578001</td>
         <td>131777.765141</td>
         <td>129648.706863</td>
         <td>112132.161616</td>
         <td>91798.715011</td>
       </tr>
       <tr>
         <th>...</th>
         <td>...</td>
         <td>...</td>
         <td>...</td>
         <td>...</td>
         <td>...</td>
         <td>...</td>
         <td>...</td>
         <td>...</td>
         <td>...</td>
         <td>...</td>
         <td>...</td>
         <td>...</td>
         <td>...</td>
         <td>...</td>
         <td>...</td>
         <td>...</td>
         <td>...</td>
         <td>...</td>
         <td>...</td>
         <td>...</td>
         <td>...</td>
       </tr>
       <tr>
         <th>1625</th>
         <td>3041.008709</td>
         <td>3324.828994</td>
         <td>2395.469265</td>
         <td>2308.737050</td>
         <td>3024.876139</td>
         <td>3307.253531</td>
         <td>1482.271672</td>
         <td>2080.586651</td>
         <td>3404.841525</td>
         <td>2913.123049</td>
         <td>...</td>
         <td>1626.248130</td>
         <td>554.882338</td>
         <td>3246.652696</td>
         <td>3765.351514</td>
         <td>1786.243434</td>
         <td>590.597996</td>
         <td>3119.301412</td>
         <td>3526.340539</td>
         <td>3128.571684</td>
         <td>2227.247013</td>
       </tr>
       <tr>
         <th>1626</th>
         <td>356867.255801</td>
         <td>348689.935124</td>
         <td>346851.549311</td>
         <td>372927.779495</td>
         <td>345165.462002</td>
         <td>360979.669247</td>
         <td>364735.213928</td>
         <td>320466.392034</td>
         <td>379752.627090</td>
         <td>390026.201243</td>
         <td>...</td>
         <td>282947.076495</td>
         <td>503996.690429</td>
         <td>349771.334353</td>
         <td>345677.687551</td>
         <td>173491.701831</td>
         <td>399923.485429</td>
         <td>359197.687162</td>
         <td>376166.710301</td>
         <td>326953.732596</td>
         <td>401299.676304</td>
       </tr>
       <tr>
         <th>1627</th>
         <td>26291.382233</td>
         <td>27847.865002</td>
         <td>28356.816852</td>
         <td>28826.128188</td>
         <td>30888.249387</td>
         <td>29509.525712</td>
         <td>32908.353274</td>
         <td>28202.412855</td>
         <td>19738.878606</td>
         <td>23828.137321</td>
         <td>...</td>
         <td>9079.839066</td>
         <td>27535.650419</td>
         <td>26673.002539</td>
         <td>25324.427145</td>
         <td>9715.138527</td>
         <td>31024.734948</td>
         <td>28805.776472</td>
         <td>28705.478299</td>
         <td>19537.250425</td>
         <td>33589.138308</td>
       </tr>
       <tr>
         <th>1628</th>
         <td>373635.872897</td>
         <td>374435.718688</td>
         <td>425780.144847</td>
         <td>414410.635963</td>
         <td>321142.352638</td>
         <td>372596.419505</td>
         <td>425871.626524</td>
         <td>355517.091009</td>
         <td>314295.114249</td>
         <td>346018.826251</td>
         <td>...</td>
         <td>174652.041234</td>
         <td>457759.006886</td>
         <td>326655.080904</td>
         <td>383970.132213</td>
         <td>238890.714726</td>
         <td>456183.199148</td>
         <td>400575.244035</td>
         <td>388277.379826</td>
         <td>389082.294175</td>
         <td>456536.266350</td>
       </tr>
       <tr>
         <th>1629</th>
         <td>4364.240925</td>
         <td>3584.293089</td>
         <td>3645.068279</td>
         <td>3990.684871</td>
         <td>4012.145214</td>
         <td>3629.953428</td>
         <td>4386.101259</td>
         <td>4647.649644</td>
         <td>3009.701602</td>
         <td>2356.009793</td>
         <td>...</td>
         <td>6102.804264</td>
         <td>4960.095760</td>
         <td>2735.980209</td>
         <td>2392.566347</td>
         <td>4644.415659</td>
         <td>2274.800507</td>
         <td>2240.936668</td>
         <td>2566.819595</td>
         <td>3239.444465</td>
         <td>2344.484279</td>
       </tr>
     </tbody>
   </table>
   <p>1630 rows × 32 columns</p>
   </div>


rdata
^^^^^

Rdata need to have at least two columns: 'Accession' and 'Description'.


#. 'Accession': is an unique value that represents those proteins in dataframe.
#. 'Description': is the fasta header from Uniprot.

.. code-block:: python

   rdata = pd.read_excel('../tests/data/proteins/general.xls', sheet_name=1)
   rdata


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
         <td>Q61823</td>
         <td>8</td>
         <td>1</td>
         <td>44.7130</td>
         <td>1.439696</td>
         <td>1.000000</td>
         <td>WT</td>
         <td>KO</td>
         <td>Programmed cell death protein 4 OS=Mus musculu...</td>
       </tr>
       <tr>
         <th>1</th>
         <td>Q91V61</td>
         <td>6</td>
         <td>0</td>
         <td>30.6978</td>
         <td>1.309501</td>
         <td>1.000000</td>
         <td>WT</td>
         <td>KO</td>
         <td>Sideroflexin-3 OS=Mus musculus OX=10090 GN=Sfx...</td>
       </tr>
       <tr>
         <th>2</th>
         <td>Q3TMQ6</td>
         <td>1</td>
         <td>0</td>
         <td>12.8896</td>
         <td>2.049949</td>
         <td>1.000000</td>
         <td>WT</td>
         <td>KO</td>
         <td>Angiogenin-4 OS=Mus musculus OX=10090 GN=Ang4 ...</td>
       </tr>
       <tr>
         <th>3</th>
         <td>Q8JZQ2</td>
         <td>4</td>
         <td>1</td>
         <td>27.5190</td>
         <td>2.126119</td>
         <td>0.999997</td>
         <td>WT</td>
         <td>KO</td>
         <td>AFG3-like protein 2 OS=Mus musculus OX=10090 G...</td>
       </tr>
       <tr>
         <th>4</th>
         <td>O89053</td>
         <td>7</td>
         <td>3</td>
         <td>47.6594</td>
         <td>1.459878</td>
         <td>0.999993</td>
         <td>WT</td>
         <td>KO</td>
         <td>Coronin-1A OS=Mus musculus OX=10090 GN=Coro1a ...</td>
       </tr>
       <tr>
         <th>...</th>
         <td>...</td>
         <td>...</td>
         <td>...</td>
         <td>...</td>
         <td>...</td>
         <td>...</td>
         <td>...</td>
         <td>...</td>
         <td>...</td>
       </tr>
       <tr>
         <th>1625</th>
         <td>Q7TST0</td>
         <td>1</td>
         <td>0</td>
         <td>5.3525</td>
         <td>1.119898</td>
         <td>0.050005</td>
         <td>WT</td>
         <td>KO</td>
         <td>Butyrophilin-like protein 1 OS=Mus musculus OX...</td>
       </tr>
       <tr>
         <th>1626</th>
         <td>P27659</td>
         <td>22</td>
         <td>7</td>
         <td>194.1972</td>
         <td>1.025275</td>
         <td>0.050002</td>
         <td>WT</td>
         <td>KO</td>
         <td>60S ribosomal protein L3 OS=Mus musculus OX=10...</td>
       </tr>
       <tr>
         <th>1627</th>
         <td>Q62148</td>
         <td>4</td>
         <td>1</td>
         <td>33.2507</td>
         <td>1.039149</td>
         <td>0.050002</td>
         <td>WT</td>
         <td>KO</td>
         <td>Retinal dehydrogenase 2 OS=Mus musculus OX=100...</td>
       </tr>
       <tr>
         <th>1628</th>
         <td>J3QM76</td>
         <td>4</td>
         <td>0</td>
         <td>22.3837</td>
         <td>1.021277</td>
         <td>0.050001</td>
         <td>WT</td>
         <td>KO</td>
         <td>Coiled-coil domain-containing protein 179 OS=M...</td>
       </tr>
       <tr>
         <th>1629</th>
         <td>P63024;P63044</td>
         <td>2</td>
         <td>0</td>
         <td>14.0456</td>
         <td>1.053974</td>
         <td>0.050000</td>
         <td>WT</td>
         <td>KO</td>
         <td>Vesicle-associated membrane protein 3 OS=Mus m...</td>
       </tr>
     </tbody>
   </table>
   <p>1630 rows × 9 columns</p>
   </div>


pdata
^^^^^

Pdata presents a description of each sample analysed. Pdata must have at least 3 columns, 'Sample', 'Condition', 'Biological'.


#. 'Sample': identifier of each sample analysed
#. 'Condition': respective group for each sample.
#. 'Biological': respective biological replicate for each sample.

While performing longitudinal analysis, users must input 'TimeCourse' column showing day/hour/time associated with respective sample.

.. code-block:: python

   pdata = pd.read_excel('../tests/data/proteins/general.xls', sheet_name=2)
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
         <th>TechRep</th>
       </tr>
     </thead>
     <tbody>
       <tr>
         <th>0</th>
         <td>VCC_KO_1_VINO</td>
         <td>KO</td>
         <td>1</td>
         <td>1</td>
       </tr>
       <tr>
         <th>1</th>
         <td>VCC_KO_1_VINO_2</td>
         <td>KO</td>
         <td>2</td>
         <td>1</td>
       </tr>
       <tr>
         <th>2</th>
         <td>VCC_KO_1_VINO_29102021</td>
         <td>KO</td>
         <td>3</td>
         <td>1</td>
       </tr>
       <tr>
         <th>3</th>
         <td>VCC_KO_1_VINO_29102021_3</td>
         <td>KO</td>
         <td>4</td>
         <td>1</td>
       </tr>
       <tr>
         <th>4</th>
         <td>VCC_KO_2_VINO</td>
         <td>KO</td>
         <td>5</td>
         <td>1</td>
       </tr>
       <tr>
         <th>5</th>
         <td>VCC_KO_2_VINO_2</td>
         <td>KO</td>
         <td>6</td>
         <td>1</td>
       </tr>
       <tr>
         <th>6</th>
         <td>VCC_KO_2_VINO_29102021</td>
         <td>KO</td>
         <td>7</td>
         <td>1</td>
       </tr>
       <tr>
         <th>7</th>
         <td>VCC_KO_2_VINO_29102021_3</td>
         <td>KO</td>
         <td>8</td>
         <td>1</td>
       </tr>
       <tr>
         <th>8</th>
         <td>VCC_KO_3_VINO</td>
         <td>KO</td>
         <td>9</td>
         <td>1</td>
       </tr>
       <tr>
         <th>9</th>
         <td>VCC_KO_3_VINO_2</td>
         <td>KO</td>
         <td>10</td>
         <td>1</td>
       </tr>
       <tr>
         <th>10</th>
         <td>VCC_KO_3_VINO_29102021</td>
         <td>KO</td>
         <td>11</td>
         <td>1</td>
       </tr>
       <tr>
         <th>11</th>
         <td>VCC_KO_3_VINO_29102021_3</td>
         <td>KO</td>
         <td>12</td>
         <td>1</td>
       </tr>
       <tr>
         <th>12</th>
         <td>VCC_KO_4_VINO</td>
         <td>KO</td>
         <td>13</td>
         <td>1</td>
       </tr>
       <tr>
         <th>13</th>
         <td>VCC_KO_4_VINO_2</td>
         <td>WT</td>
         <td>14</td>
         <td>1</td>
       </tr>
       <tr>
         <th>14</th>
         <td>VCC_KO_4_VINO_29102021</td>
         <td>WT</td>
         <td>15</td>
         <td>1</td>
       </tr>
       <tr>
         <th>15</th>
         <td>VCC_KO_4_VINO_29102021_3</td>
         <td>WT</td>
         <td>16</td>
         <td>1</td>
       </tr>
       <tr>
         <th>16</th>
         <td>VCC_WT_1_VIN</td>
         <td>WT</td>
         <td>1</td>
         <td>1</td>
       </tr>
       <tr>
         <th>17</th>
         <td>VCC_WT_1_VIN_2</td>
         <td>WT</td>
         <td>2</td>
         <td>1</td>
       </tr>
       <tr>
         <th>18</th>
         <td>VCC_WT_1_VIN_29102021</td>
         <td>WT</td>
         <td>3</td>
         <td>1</td>
       </tr>
       <tr>
         <th>19</th>
         <td>VCC_WT_1_VIN_29102021_2</td>
         <td>WT</td>
         <td>4</td>
         <td>1</td>
       </tr>
       <tr>
         <th>20</th>
         <td>VCC_WT_2_VIN</td>
         <td>WT</td>
         <td>5</td>
         <td>1</td>
       </tr>
       <tr>
         <th>21</th>
         <td>VCC_WT_2_VIN_2</td>
         <td>WT</td>
         <td>6</td>
         <td>1</td>
       </tr>
       <tr>
         <th>22</th>
         <td>VCC_WT_2_VIN_29102021</td>
         <td>WT</td>
         <td>7</td>
         <td>1</td>
       </tr>
       <tr>
         <th>23</th>
         <td>VCC_WT_2_VIN_29102021_2</td>
         <td>WT</td>
         <td>8</td>
         <td>1</td>
       </tr>
       <tr>
         <th>24</th>
         <td>VCC_WT_3_VIN</td>
         <td>WT</td>
         <td>9</td>
         <td>1</td>
       </tr>
       <tr>
         <th>25</th>
         <td>VCC_WT_3_VIN_2</td>
         <td>WT</td>
         <td>10</td>
         <td>1</td>
       </tr>
       <tr>
         <th>26</th>
         <td>VCC_WT_3_VIN_29102021</td>
         <td>WT</td>
         <td>11</td>
         <td>1</td>
       </tr>
       <tr>
         <th>27</th>
         <td>VCC_WT_3_VIN_29102021_2</td>
         <td>WT</td>
         <td>12</td>
         <td>1</td>
       </tr>
       <tr>
         <th>28</th>
         <td>VCC_WT_4_VIN</td>
         <td>WT</td>
         <td>13</td>
         <td>1</td>
       </tr>
       <tr>
         <th>29</th>
         <td>VCC_WT_4_VIN_2</td>
         <td>WT</td>
         <td>14</td>
         <td>1</td>
       </tr>
       <tr>
         <th>30</th>
         <td>VCC_WT_4_VIN_29102021</td>
         <td>WT</td>
         <td>15</td>
         <td>1</td>
       </tr>
       <tr>
         <th>31</th>
         <td>VCC_WT_4_VIN_29102021_2</td>
         <td>WT</td>
         <td>16</td>
         <td>1</td>
       </tr>
     </tbody>
   </table>
   </div>


Additional Informations
-----------------------

User can also define and optimize some extra parameters that are in function OmicScope.


#. 
   **ControlGroup** (default = None): User can define control group ('ControlGroup=None', default) to perform comparisons against an specific group (this group have to be explicit in column Conditions on pdata table)

#. 
   **ExperimentalDesign** (default = 'static'): comparison among independent groups are called 'static' experimental design. On the other hand, if the experiment take into account several time points, than it is performing an 'longitudinal' experimental design (in this case, pdata table must present 'TimeCourse' column).

#. 
   **pvalue** (defaul = 'pAdjusted'): defines the kinds of statitics that will be used to report differentially regulated proteins, which the options are: nominal p-value ('pvalue'); Benjamini-Hochberg Adjusted p-value ('pAdjusted'); or Tukey post-hoc correction ('pTukey', just for multiple group comparisons in static experiments).

#. 
   **PValue_cutoff** (default = 0.05): Statistical cutoff to consider proteins differentially regulated. 

#. 
   **FoldChange_cutoff** (default = 0): cutoff of abundance ratio to consider proteins differentially regulated. 

#. 
   **logTransformed** (default = False): Usually softwares report abundance in their nominal values, requiring a log-transformation of the values. If user perform transformation before OmicScope workflow, logTransformed=True.

#. 
   **ExcludeKeratins** (default = True): Since keratins are considered sample contaminant in most studies, OmicScope can exclude than from final results.

#. 
   **degrees_of_freedom** (default = 2 ): For longitudinal analysis, users can optmize the parameters according to their studies choosing a greater degree of freedom to perform analysis.
