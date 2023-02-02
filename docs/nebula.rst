
.. code-block:: python

   import sys
   sys.path.insert(1, 'C:/Users/Guilherme/omicscope/omicscope/src/')

Nebula
======

Nebula is an OmicScope module that provides multi-data integration, providing several functions for analyzing independent proteomics experiments. Due to the *General* input workflow performed by the OmicScope function, Nebula is also capable of analyzing multi-omics data.

To run the Nebula workflow, the input must be a folder with .omics files, which are exported from the OmicScope and EnrichmentScope workflows.

Exporting **.omics** files - *object.savefile()*
--------------------------------------------------------

OmicScope and EnrichmentScope have the **savefile function** to export *quantitative data* and *quantitative and enrichment data*\ , respectively. To export, the user only has to call the savefile function and add the folder path, from which data will be saved.

.. code-block:: python

   # OmicScope Example
   import omicscope as omics

   df = omics.OmicScope('../tests/data/proteins/progenesis.csv', Method='Progenesis', ControlGroup='WT')
   df.savefile(PATH_TO_SAVE)

   # EnrichmentScope Example
   # WARNING: EnrichmentScope also includes QUANTITATIVE DATA
   enr = omics.EnrichmentScope(df)
   enr.savefile(PATH_TO_SAVE)

By default, the file's name is the condition extracted from OmicScope/EnrichmentScope . In the case above, the exported file name is: 'WT-KO.omics' since those are the conditions reported by the user.

.. code-block:: python

   df.Conditions

Nebula Object
=============

Nebula processes quantitative data and (if applied) the enrichment data. From the .omics file, Nebula extracts the experimental group (which the user can modify in the .omics file) and uses that to name the respective experiment. Nebula then reports the number of groups/experiments imported, their names, and if there is enrichment data to be included.

.. code-block:: python

   import omicscope as omics

   nebula = omics.Nebula('../tests/data/MultipleGroups/omics_file/')

.. code-block::

   You imported your data successfully!
           Data description:
           1. N groups imported: 4
           2. Groups: COVID,Covid_HB,NEURONS,SH_DIF
           3. N groups with enchment data: 4




The names for each group also can be changed in command line, as shown below:

.. code-block:: python

   nebula.groups = ['Astrocytes', 'Human_Brain', 'Neurons', 'SHSY5Y']
   nebula.groups

.. code-block::

   ['Astrocytes', 'Human_Brain', 'Neurons', 'SHSY5Y']




Figures and plots
=================

Barplot - *object.barplot()*
--------------------------------

Nebula barplot shows quantified and differentially regulated proteins/genes across all studies. 

.. code-block:: python

   nebula.barplot(dpi=90)


.. image:: nebula_files/nebula_12_0.png
   :target: nebula_files/nebula_12_0.png
   :alt: png


Enrichment Dotplot - *object.dotplot_enrichment()*
------------------------------------------------------

If the *.omics file* contains enrichment results, it is possible to compare the enrichment for each group with the dotplot_enrichment() function. By default, according to p-value, a list of the top 5 Terms for each group is created and then used to filter each enrichment data to be compared. 

.. code-block:: python

   nebula.dotplot_enrichment(top = 5, dpi=90)


.. image:: nebula_files/nebula_14_0.png
   :target: nebula_files/nebula_14_0.png
   :alt: png


Differentially regulated - *object.diff_reg()*
--------------------------------------------------

The comparison between groups can be performed only at differentially regulated levels, showing the number of proteins that are up- and down-regulated.

.. code-block:: python

   nebula.diff_reg(dpi=90)


.. image:: nebula_files/nebula_16_0.png
   :target: nebula_files/nebula_16_0.png
   :alt: png


Protein Overlap - *object.protein_overlap()*
------------------------------------------------

The Venn Diagram is a classical plot used to visualize the overlap and uniqueness between groups. Despite several tools that quickly reproduce venn diagrams (such as Interactive Venn), these plots are limited in the number of groups that can be compared, since all geometric figures need to overlap one another.

Since it is not uncommon for proteomics studies to evaluate several groups, Nebula plots upset plot. In upset plot, several groups can be compared at once; in the low-left barplot is described the number of entities associated with each group; in the up-right barplot is shown intersection size for each comparison, which are highlighted in the colored and linked circles in the frame.

The protein overlap function performs comparisons between all groups at the protein level.

.. code-block:: python

   nebula.protein_overlap(dpi=90)


.. image:: nebula_files/nebula_18_0.png
   :target: nebula_files/nebula_18_0.png
   :alt: png


Enrichment Overlap - *object.protein_overlap()*
---------------------------------------------------

Working in the same way that protein_overlap, enrichment_overlap performs the same visualization for terms, that was assigned for enrichment analysis. 

.. code-block:: python

   nebula.enrichment_overlap(dpi=90)


.. image:: nebula_files/nebula_20_0.png
   :target: nebula_files/nebula_20_0.png
   :alt: png


Pearson Correlation - *object.correlation()*
------------------------------------------------

Pearson correlation can be used to evaluate how much the proteome (or differentially regulated proteins) shares similarities in protein levels. To make the data easier to visualize, we plot the pair-wise comparison in a heatmap with hierarchical clustering.

.. code-block:: python


   nebula.correlation(dpi=90)


.. image:: nebula_files/nebula_22_0.png
   :target: nebula_files/nebula_22_0.png
   :alt: png


Fisher's test - *object.fisher_heatmap()*
---------------------------------------------

To determine the statistical significance of the comparison between the groups, a pair-wise Fisher's exact test could be applied. A heatmap is plotted with the values, and labels are shown in the log10-scale transformation.

.. code-block:: python

   nebula.fisher_heatmap(pvalue = 0.05,dpi=90)


.. image:: nebula_files/nebula_24_0.png
   :target: nebula_files/nebula_24_0.png
   :alt: png



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
         <th>Astrocytes</th>
         <th>Human_Brain</th>
         <th>Neurons</th>
         <th>SHSY5Y</th>
       </tr>
     </thead>
     <tbody>
       <tr>
         <th>Astrocytes</th>
         <td>0.000000</td>
         <td>0.660228</td>
         <td>0.389940</td>
         <td>0.210782</td>
       </tr>
       <tr>
         <th>Human_Brain</th>
         <td>0.660228</td>
         <td>0.000000</td>
         <td>0.281677</td>
         <td>0.088847</td>
       </tr>
       <tr>
         <th>Neurons</th>
         <td>0.389940</td>
         <td>0.281677</td>
         <td>0.000000</td>
         <td>0.325381</td>
       </tr>
       <tr>
         <th>SHSY5Y</th>
         <td>0.210782</td>
         <td>0.088847</td>
         <td>0.325381</td>
         <td>0.000000</td>
       </tr>
     </tbody>
   </table>
   </div>


Protein Network - *object.network()*
----------------------------------------

Network visualization provides an overview of individual proteins shared among groups. Using a systems biology approach, network theory can help find communities/modules and group information based on similarities. Since there are several programs designed to plot graphs/networks, Nebula exports the information as a .graphml file, which can be imported into Cytoscape, Gephi, and so on.

.. code-block:: python

   nebula.network(dpi=90)


.. image:: nebula_files/nebula_26_0.png
   :target: nebula_files/nebula_26_0.png
   :alt: png


.. code-block::

   <networkx.classes.graph.Graph at 0x1d2dd558f90>




Group Network - *object.group_network()*
--------------------------------------------

Network function can be very slow, due to several proteins that must be plotted on the graph, so Nebula also has the group_network function. This function filters the proteins based on p-value (default: protein_pvalue=0.05), followed by a pair-wise Fisher test, which is used to link each group according to the p-value cutoff (default: graph_pvalue=0.05). The links are labeled on the log10-scale.

.. code-block:: python


   nebula.group_network(protein_pvalue=1, graph_pvalue=0.05, dpi=90)


.. image:: nebula_files/nebula_28_0.png
   :target: nebula_files/nebula_28_0.png
   :alt: png


.. code-block::

   <networkx.classes.graph.Graph at 0x1d2de6fce10>




Circular graphs - *object.circular_path()*
----------------------------------------------

The circular plot was designed to compare groups that were enriched for a determined term according to their respective differentially regulated proteins. Additionally, the proteins are plotted with their respective regulations, being up-(red) or down-regulated (blue).

**ATTENTION**\ : To use circular_path, the system must have R installed with the circlize package.

.. code-block:: python

   nebula.circular_path('Amyotrophic lateral sclerosis')


.. image:: nebula_files/nebula_30_0.png
   :target: nebula_files/nebula_30_0.png
   :alt: png


Circos plot - *object.circos_plot()*
----------------------------------------

Circos is a software designed to visualize complex data (e.g. data from multiple groups) in circular mode. In Nebula, circos was implemented to allow visualization of proteins differentially regulated at once, highlighting those ones that are shared among groups (darkcyan links) as well as the regulation of the proteins (on edge heatmap). If enrichment analysis is present in the .omics file, circos_plot incorporates the shared enrichment terms (black links) to give an idea of the number of pathways shared between groups.

**ATTENTION**\ : To use circos_plot, the system must have Perl installed and configured according to Circos software and the appropriate configuration for the system.

.. code-block:: python

   nebula.circos_plot()


.. image:: nebula_files/nebula_32_0.png
   :target: nebula_files/nebula_32_0.png
   :alt: png

