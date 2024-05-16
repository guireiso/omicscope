EnrichmentScope Object
======================

`omics.EnrichmentScope() <https://omicscope.readthedocs.io/en/latest/reference/omicscope.html#omicscope.EnrichmentScope>`__

The EnrichmentScope is a dedicated module within OmicScope designed for
conducting *in silico* enrichment analyses. It offers two primary types
of enrichment analyses: Over-Representation Analysis (ORA) and Gene-Set
Enrichment Analysis (GSEA). These analyses are powered by the
`Enrichr <https://maayanlab.cloud/Enrichr/>`__ libraries, providing
access to over 224 databases for thorough analysis.

To perform enrichment analysis, users must initially import and perform
statistical analysis using the OmicScope module.

The OmicScope object is used to perform Enrichment Analysis in the
EnrichmentScope module. It’s important to note that ORA uses
differentially regulated proteins found in OmicScope, while the GSEA
algorithm uses all quantified proteins to perform the statistical
analysis.

To perform enrichment analysis, users also need to select appropriate
databases, with ``KEGG_2021_Human`` used by default. Users can also
define alternative background sizes, target organisms, or pAdjusted
cutoffs to consider terms enriched.

.. code:: ipython3

    import omicscope as omics
    
    data = omics.OmicScope('../../tests/data/proteins/progenesis.xls', Method = 'Progenesis')
    
    ora = omics.EnrichmentScope(data, Analysis='ORA', dbs = ['KEGG_2021_Human'])


.. parsed-literal::

    OmicScope v 1.4.0 For help: https://omicscope.readthedocs.io/en/latest/ or https://omicscope.ib.unicamp.brIf you use  in published research, please cite:
    'Reis-de-Oliveira, G., et al (2024). OmicScope unravels systems-level insights from quantitative proteomics data 
    
    User already performed statistical analysis
    OmicScope identifies: 697 deregulations
    

Enrichment results
------------------

Enrichment results are stored in ``object.results`` as a table
(DataFrame), with the following columns:

+-----------------------------------+-----------------------------------+
| Column Name                       | Description                       |
+===================================+===================================+
| Gene_set                          | Gene set library used for         |
|                                   | enrichment analysis               |
+-----------------------------------+-----------------------------------+
| Term                              | Enriched term                     |
+-----------------------------------+-----------------------------------+
| Overlap                           | Ratio of proteins overlapped in   |
|                                   | the experimental gene list and    |
|                                   | the total number of genes in the  |
|                                   | library term                      |
+-----------------------------------+-----------------------------------+
| P-value                           | Nominal p-value from Fisher’s     |
|                                   | exact test                        |
+-----------------------------------+-----------------------------------+
| Adjusted P-value                  | Adjusted p-value according to     |
|                                   | Benjamini-Hochberg (pAdjusted)    |
+-----------------------------------+-----------------------------------+
| Combined Score *(ORA only)*       | Score for the enrichment analysis |
+-----------------------------------+-----------------------------------+
| Genes                             | Genes overlapped between          |
|                                   | experimental data and the         |
|                                   | database                          |
+-----------------------------------+-----------------------------------+
| -log10(pAdj)                      | Log-transformed pAdjusted value   |
+-----------------------------------+-----------------------------------+
| N_Proteins                        | Number of proteins overlapped     |
|                                   | between the experimental gene     |
|                                   | list and the target library term  |
+-----------------------------------+-----------------------------------+
| Regulation                        | Log2(foldchange) of each protein  |
|                                   | overlapped                        |
+-----------------------------------+-----------------------------------+
| Down-regulated                    | Number of down-regulated proteins |
+-----------------------------------+-----------------------------------+
| Up-regulated                      | Number of up-regulated proteins   |
+-----------------------------------+-----------------------------------+
| Normalized Enrichment Score (NES) | An attempt to predict the effect  |
| *(GSEA only)*                     | of proteins on pathways (specific |
|                                   | to GSEA analysis)                 |
+-----------------------------------+-----------------------------------+

.. code:: ipython3

    ora.results.head(4)




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
          <th>index</th>
          <th>Gene_set</th>
          <th>Term</th>
          <th>Overlap</th>
          <th>P-value</th>
          <th>Adjusted P-value</th>
          <th>Old P-value</th>
          <th>Old Adjusted P-value</th>
          <th>Odds Ratio</th>
          <th>Combined Score</th>
          <th>Genes</th>
          <th>-log10(pAdj)</th>
          <th>N_Proteins</th>
          <th>regulation</th>
          <th>down-regulated</th>
          <th>up-regulated</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>0</td>
          <td>KEGG_2021_Human</td>
          <td>Parkinson disease</td>
          <td>58/249</td>
          <td>1.704579e-31</td>
          <td>4.789868e-29</td>
          <td>0</td>
          <td>0</td>
          <td>9.082385</td>
          <td>643.458087</td>
          <td>[NDUFA11, CALML3, COX6A1, UBE2L3, TUBB8, UCHL1...</td>
          <td>28.319676</td>
          <td>58</td>
          <td>[0.2670808325175823, -0.10715415448907055, 0.7...</td>
          <td>33</td>
          <td>25</td>
        </tr>
        <tr>
          <th>1</th>
          <td>1</td>
          <td>KEGG_2021_Human</td>
          <td>Pathways of neurodegeneration</td>
          <td>78/475</td>
          <td>6.471702e-31</td>
          <td>9.092742e-29</td>
          <td>0</td>
          <td>0</td>
          <td>6.000855</td>
          <td>417.135594</td>
          <td>[NDUFA11, CALML3, ATP2A1, COX6A1, UBE2L3, TUBB...</td>
          <td>28.041305</td>
          <td>78</td>
          <td>[0.2670808325175823, -0.10715415448907055, -0....</td>
          <td>51</td>
          <td>27</td>
        </tr>
        <tr>
          <th>2</th>
          <td>2</td>
          <td>KEGG_2021_Human</td>
          <td>Prion disease</td>
          <td>54/273</td>
          <td>1.174929e-25</td>
          <td>1.100517e-23</td>
          <td>0</td>
          <td>0</td>
          <td>7.318264</td>
          <td>420.093386</td>
          <td>[NDUFA11, COX6A1, TUBB8, PPP3CB, TUBB6, PPP3CC...</td>
          <td>22.958403</td>
          <td>54</td>
          <td>[0.2670808325175823, 0.7932637717587971, -0.33...</td>
          <td>29</td>
          <td>25</td>
        </tr>
        <tr>
          <th>3</th>
          <td>3</td>
          <td>KEGG_2021_Human</td>
          <td>Amyotrophic lateral sclerosis</td>
          <td>61/364</td>
          <td>8.377698e-25</td>
          <td>5.885333e-23</td>
          <td>0</td>
          <td>0</td>
          <td>6.014281</td>
          <td>333.426032</td>
          <td>[NDUFA11, COX6A1, ACTG1, TUBB8, ACTR1A, PPP3CB...</td>
          <td>22.230229</td>
          <td>61</td>
          <td>[0.2670808325175823, 0.7932637717587971, -0.22...</td>
          <td>38</td>
          <td>23</td>
        </tr>
      </tbody>
    </table>
    </div>



Background - ORA only
~~~~~~~~~~~~~~~~~~~~~

When conducting Over-Representation Analysis (ORA), the background gene
list assumes a pivotal role in enrichment analysis by serving as the
reference set against which the experimental gene list is compared. To
put it simply, the background gene list encompasses all the genes or
proteins that could potentially be present in the experimental dataset.

By default, when ``background = None``, EnrichmentScope includes all
genes found in the database as part of the background. Alternatively,
users have the option to set ``background = True`` to encompass all
proteins identified in the experiment. They can also use
``background = int`` to specify the background size, which could be, for
instance, the reviewed human proteome in the case of human experiments
(although this is not recommended). Another option is to define
background = ``[ListOfGenes]`` to specify a particular gene set for
comparative analysis.

}Plots and Figures
------------------

EnrichmentScope introduces a variety of figures that aim to integrate
the enrichment outcomes with the differentially regulated proteins in
biological systems.

Users can choose between saving the generated plots in vector format
(using ``vector=True``) or in .png format (with ``vector=False``). They
have the flexibility to set the desired figure resolution (using
``dpi=300``) and specify a file path for saving the plots. Moreover,
users can adjust the color schemes of the plots using the “palettes”
command, selecting color palettes from Matplotlib. These customizable
options empower users to create informative and visually appealing
visualizations that cater to their specific requirements and preferences

Dotplot - `object.dotplot() <https://omicscope.readthedocs.io/en/latest/reference/enrichmentvis.html#omicscope.EnrichmentAnalysis.EnrichmentVisualization.dotplot>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``dotplot`` function ranks enriched terms on the y-axis based on
their adjusted p-values, while the x-axis represents the adjusted
p-values. Additionally, the size of each dot is proportional to
-log10(pAdjusted), providing an indication of the significance of the
enrichment. Furthermore, the color of each dot is coded based on the
number of proteins used in the enrichment analysis.

**How to interpret**: The positioning of each dot on the plot indicates
the statistical significance of the term, with more statistically
significant terms located towards the top-right side of the plot.
Additionally, the color of each dot corresponds to the number of
proteins associated with that term, with darker blue indicating a higher
number of associated proteins.

.. code:: ipython3

    ora.dotplot(dpi=90, palette='PuBu')



.. image:: enrichmentscope_files%5Cenrichmentscope_7_0.png


Heatmap - `object.Heatmap() <https://omicscope.readthedocs.io/en/latest/reference/enrichmentvis.html#omicscope.EnrichmentAnalysis.EnrichmentVisualization.heatmap>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The heatmap is a valuable tool within the EnrichmentScope workflow,
aiding in the visualization of proteins that are shared between enriched
terms, helping to reduce data redundancy. In this heatmap, proteins are
depicted on the y-axis, while terms are assigned to the x-axis.

By default, the heatmap colors are mapped according to the adjusted
p-value. However, users have the option to color each protein based on
its fold-change by setting ``foldchange=True``.

**How to interpret**: When looking for specific proteins, users can
identify the specific pathways (terms) associated with those proteins.
Conversely, when exploring several pathways, users can observe the group
of proteins that are shared between those pathways (terms). In the
examples provided below, we highlight the default parameters and color
coding based on fold change.

.. code:: ipython3

    ora.heatmap(linewidths=0.5)



.. image:: enrichmentscope_files%5Cenrichmentscope_9_0.png


.. code:: ipython3

    # color based on protein fold-change
    ora.heatmap(linewidths=0.5, foldchange=True)



.. image:: enrichmentscope_files%5Cenrichmentscope_10_0.png


Number of DEPs - `object.number_deps() <https://omicscope.readthedocs.io/en/latest/reference/enrichmentvis.html#omicscope.EnrichmentAnalysis.EnrichmentVisualization.number_deps>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``number_deps`` function counts the number of up- and down-regulated
entities (x-axis) and plots them according to each enriched term
(y-axis). In this plot, sizes indicate the number of proteins found in
each group.

**How to interpret**: For users performing ORA and GSEA analyses,
questions often arise about the number of up- and down-regulated
proteins associated with each term.

.. code:: ipython3

    ora.number_deps(palette=['firebrick','darkcyan'] ,dpi = 90)



.. image:: enrichmentscope_files%5Cenrichmentscope_12_0.png


Enrichment Network - `object.enrichment_network() <https://omicscope.readthedocs.io/en/latest/reference/enrichmentvis.html#omicscope.EnrichmentAnalysis.EnrichmentVisualization.enrichment_network>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In proteomics, major pathways frequently share several proteins, and
visualizing pathways and proteins together in a network can be highly
informative.

The Enrichment Network function visually connects terms to their
associated proteins. In this visualization, terms are depicted in gray,
and the node size is proportional to ``-log10(p-adjusted)``. Proteins
are represented uniformly in size and are color-coded based on their
fold-change. Labels can be added to the plot by using the
``labels``\ =True option (default: ``False``).

**Note**: Note: Visualizing graphs can be complex, particularly when
dealing with substantial amounts of information. To achieve the best
visualization possible, several software options, such as Cytoscape and
Gephi, have been specifically designed for this purpose. Users can
export the plot to these external tools by specifying
``save=PATH_TO_SAVE``.

.. code:: ipython3

    ora.enrichment_network(top = 10, dpi = 90)



.. image:: enrichmentscope_files%5Cenrichmentscope_14_0.png




.. parsed-literal::

    [<networkx.classes.graph.Graph at 0x182bed80d90>]



Enrichment Map - `object.enrichment_map() <https://omicscope.readthedocs.io/en/latest/reference/enrichmentvis.html#omicscope.EnrichmentAnalysis.EnrichmentVisualization.enrichment_map>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An advantageous aspect of employing graphical representations in
enrichment analysis is their ability to reduce data redundancy. The
``enrichment_map`` function takes advantage of this by rendering nodes
as terms and edges as similarity scores, typically calculated using
statistical metrics such as Jaccard similarity (default). If users opt
to enable ``modules=True``, the Louvain method is utilized to identify
communities within the network. Each community is assigned a unique
term, typically the one with the highest degree, to describe the
community when ``labels=True`` is specified.

Similar to the ``enrichment_network`` function, users can easily export
the generated enrichment map to external tools for further exploration
and visualization by adding ``save=PATH_TO_SAVE``.

**How to interpret**: While aiming to investigate pathways that share
proteins, users can look inside modules to identify pathways that
present high similarity regarding protein presence. On the other hand,
while avoiding redundancy, users can look for the node that presents a
higher degree (number of connections) inside each module and/or a lower
p-value and consider that node to represent the whole module.

.. code:: ipython3

    ora.enrichment_map(dpi=90, modules=True)



.. image:: enrichmentscope_files%5Cenrichmentscope_16_0.png




.. parsed-literal::

    [<networkx.classes.graph.Graph at 0x182bc0941d0>]



