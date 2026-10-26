How to make pdata
=================

Statistical Analysis Overview
-----------------------------

OmicScope offers a toolkit for conducting differential proteomics
analyses, covering statistical approaches for both static and
longitudinal experimental designs (as illustrated in the figure below).
By default, OmicScope assumes a static workflow (designated as
‘ExperimentalDesign=static’). In this mode, it employs t-tests or
analysis of variance (ANOVA) for statistical analysis. For longitudinal
analyses (designated as ‘ExperimentalDesign=longitudinal’), OmicScope
assumes protein abundance varies over time according to a natural cubic
spline, as suggested by Storey’s 2005. This method assesses differences
within and between groups over time.

After calculating nominal p-values, OmicScope applies the
Benjamini-Hochberg correction to account for multiple hypothesis testing
and reports the adjusted p-value (pAdjusted). Alternatively, users have
the flexibility to import statistical analyses from other software tools
by including a “pvalue” or “pAdjusted” column in the rdata using the
General input method.

Pdata Role
----------

Pdata (also known as phenotype data or metadata) plays a crucial role in
allowing OmicScope to correctly conduct statistical analysis. To perform
this task, Pdata must contain as much information as possible to compare
2 or more groups, different time courses, or different classes
throughout time.

When importing data into OmicScope, users must select the appropriate
method for data handling. While performing static analysis for
Progenesis, Proteome Discoverer, PatternLab and Spectronaut (long
report with R.Condition), OmicScope can automatically identify and
select groups based on the software’s output. The General method keeps
pdata in the third sheet of its workbook. However, when performing
longitudinal analysis or using MaxQuant, DIA-NN, FragPipe, RPKM or a
Spectronaut pivot report, OmicScope requires a pdata file to select the
appropriate test.

The ``Sample`` column must list the sample names used in the
quantitative file (run or experiment names). OmicScope ignores folders,
file extensions (``.raw``, ``.d``, ``.mzML``, ...), letter case and
surrounding spaces when matching them, and reports the available names
when a sample cannot be found. pdata can be an Excel file (.xlsx/.xls),
a comma-, semicolon- or tab-separated text file, or a pandas DataFrame.

For all cases, OmicScope allows users to incorporate external pdata into
the workflow, which helps tailor the statistical analysis to the
specific experimental design (as described below).

Static Experimental Design
--------------------------

Static Workflow
~~~~~~~~~~~~~~~

Most proteomics experiments aim to compare proteomic signatures between
independent groups, which is why OmicScope defaults to a static
experimental design (designated as ``'ExperimentalDesign=static'``). The
static workflow involves two main statistical tests: the t-test and
ANOVA.

When comparing two groups, OmicScope conducts an independent t-test
(when ``independent_ttest=True``) if the groups are independent, or a
paired t-test (when ``independent_ttest=False``) if the groups are
related. For comparisons involving more than two groups, OmicScope
employs a one-way ANOVA. Additionally, for proteins with a pAdjusted
value less than 0.05, OmicScope performs a Tukey post-hoc correction.
This helps identify and highlight the groups with significant
differences.

Static Pdata
~~~~~~~~~~~~

To create a pdata file for running the static workflow, users should
include the following columns:

1. **Sample:** The name of each sample to be analyzed, matching those in
   the first row of the Assay sheet.
2. **Condition:** Respective group for each sample. All technical and
   biological replicates belonging to an experimental condition should
   have the same identifier here.
3. **Biological:** Respective biological replicate for each sample. If
   two or more technical replicates were used for a single biological
   replicate, those replicates should have the same identifier here.

Pdata Example 1
~~~~~~~~~~~~~~~

In the example below, each sample is assigned to a specific Condition,
and each biological replicate is reported. Here, two distinct conditions
were documented, and each biological replicate was acquired twice. Once
this ‘pdata’ is integrated into the OmicScope workflow, the process
involves calculating the mean of technical replicates and performing an
independent t-test as the statistical test.

.. code:: python

    # Pdata for static experimental designs
    import pandas as pd
    pdata = pd.read_excel('..\\..\\tests\\data\\proteins\\general.xlsx', sheet_name=2)
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



.. code:: python

    print('Number of Conditions: ' + str(len(pdata.Condition.drop_duplicates())))


.. parsed-literal::

    Number of Conditions: 2
    

Longitudinal Experimental Design
--------------------------------

Longitudinal Workflow
~~~~~~~~~~~~~~~~~~~~~

To accommodate the potential complexities of longitudinal experimental
designs, OmicScope categorizes these experiments into two primary types:

1. *Within-group experiments*: These designs aim to identify
   differentially regulated proteins over time within a single group.
2. *Between-group experiments*: These designs aim to detect differential
   protein regulation over time by comparing different groups.

Pdata workflow
~~~~~~~~~~~~~~

OmicScope manages these distinctions much like the static workflow,
examining the number of conditions (#conditions) in the ‘Condition’
column. It selects “Within-group” if the #conditions is equal to 1, and
“Between-group” if the #conditions exceed 1. Additionally, in the
longitudinal workflow, the user is **required to add a “TimeCourse”**
column to define the sampling frequency of the study.

Pdata Example 2
~~~~~~~~~~~~~~~

In the example below, the ‘pdata’ contains two distinct groups (12
Control and 12 Treatment) in the ‘Condition’ column, indicating a
Between-group analysis. Additionally, the ``TimeCourse`` column includes
4 time points, and each biological replicate was acquired twice.

.. code:: python

    pdata = pd.read_excel('..\\../tests/data/proteins/longitudinal_pdata.xlsx', sheet_name=0)
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
          <th>TimeCourse</th>
          <th>Biological</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>Sample1_Day1_Bio1_1</td>
          <td>Control</td>
          <td>1</td>
          <td>1</td>
        </tr>
        <tr>
          <th>1</th>
          <td>Sample1_Day1_Bio1_2</td>
          <td>Control</td>
          <td>1</td>
          <td>1</td>
        </tr>
        <tr>
          <th>2</th>
          <td>Sample2_Day1_Bio2_1</td>
          <td>Control</td>
          <td>1</td>
          <td>2</td>
        </tr>
        <tr>
          <th>3</th>
          <td>Sample2_Day1_Bio2_2</td>
          <td>Control</td>
          <td>1</td>
          <td>2</td>
        </tr>
        <tr>
          <th>4</th>
          <td>Sample3_Day1_Bio3_1</td>
          <td>Control</td>
          <td>1</td>
          <td>3</td>
        </tr>
        <tr>
          <th>5</th>
          <td>Sample3_Day1_Bio3_2</td>
          <td>Control</td>
          <td>1</td>
          <td>3</td>
        </tr>
        <tr>
          <th>6</th>
          <td>Sample4_Day2_Bio1_1</td>
          <td>Control</td>
          <td>3</td>
          <td>4</td>
        </tr>
        <tr>
          <th>7</th>
          <td>Sample4_Day2_Bio1_2</td>
          <td>Control</td>
          <td>3</td>
          <td>4</td>
        </tr>
        <tr>
          <th>8</th>
          <td>Sample5_Day2_Bio2_1</td>
          <td>Control</td>
          <td>3</td>
          <td>5</td>
        </tr>
        <tr>
          <th>9</th>
          <td>Sample5_Day2_Bio2_2</td>
          <td>Control</td>
          <td>3</td>
          <td>5</td>
        </tr>
        <tr>
          <th>10</th>
          <td>Sample6_Day2_Bio3_1</td>
          <td>Control</td>
          <td>3</td>
          <td>6</td>
        </tr>
        <tr>
          <th>11</th>
          <td>Sample6_Day2_Bio3_2</td>
          <td>Control</td>
          <td>3</td>
          <td>6</td>
        </tr>
        <tr>
          <th>12</th>
          <td>Sample7_Day3_Bio1_1</td>
          <td>Control</td>
          <td>5</td>
          <td>7</td>
        </tr>
        <tr>
          <th>13</th>
          <td>Sample7_Day3_Bio1_2</td>
          <td>Control</td>
          <td>5</td>
          <td>7</td>
        </tr>
        <tr>
          <th>14</th>
          <td>Sample8_Day3_Bio2_1</td>
          <td>Control</td>
          <td>5</td>
          <td>8</td>
        </tr>
        <tr>
          <th>15</th>
          <td>Sample8_Day3_Bio2_2</td>
          <td>Control</td>
          <td>5</td>
          <td>8</td>
        </tr>
        <tr>
          <th>16</th>
          <td>Sample9_Day3_Bio3_1</td>
          <td>Control</td>
          <td>5</td>
          <td>9</td>
        </tr>
        <tr>
          <th>17</th>
          <td>Sample9_Day3_Bio3_2</td>
          <td>Control</td>
          <td>5</td>
          <td>9</td>
        </tr>
        <tr>
          <th>18</th>
          <td>Sample10_Day4_Bio1_1</td>
          <td>Control</td>
          <td>7</td>
          <td>10</td>
        </tr>
        <tr>
          <th>19</th>
          <td>Sample10_Day4_Bio1_2</td>
          <td>Control</td>
          <td>7</td>
          <td>10</td>
        </tr>
        <tr>
          <th>20</th>
          <td>Sample11_Day4_Bio2_1</td>
          <td>Control</td>
          <td>7</td>
          <td>11</td>
        </tr>
        <tr>
          <th>21</th>
          <td>Sample11_Day4_Bio2_2</td>
          <td>Control</td>
          <td>7</td>
          <td>11</td>
        </tr>
        <tr>
          <th>22</th>
          <td>Sample12_Day5_Bio3_1</td>
          <td>Control</td>
          <td>7</td>
          <td>12</td>
        </tr>
        <tr>
          <th>23</th>
          <td>Sample12_Day5_Bio3_2</td>
          <td>Control</td>
          <td>7</td>
          <td>12</td>
        </tr>
        <tr>
          <th>24</th>
          <td>Sample13_Day1_Bio1_1</td>
          <td>Treatment</td>
          <td>1</td>
          <td>13</td>
        </tr>
        <tr>
          <th>25</th>
          <td>Sample13_Day1_Bio1_2</td>
          <td>Treatment</td>
          <td>1</td>
          <td>13</td>
        </tr>
        <tr>
          <th>26</th>
          <td>Sample14_Day1_Bio2_1</td>
          <td>Treatment</td>
          <td>1</td>
          <td>14</td>
        </tr>
        <tr>
          <th>27</th>
          <td>Sample14_Day1_Bio2_2</td>
          <td>Treatment</td>
          <td>1</td>
          <td>14</td>
        </tr>
        <tr>
          <th>28</th>
          <td>Sample15_Day1_Bio3_1</td>
          <td>Treatment</td>
          <td>1</td>
          <td>15</td>
        </tr>
        <tr>
          <th>29</th>
          <td>Sample15_Day1_Bio3_2</td>
          <td>Treatment</td>
          <td>1</td>
          <td>15</td>
        </tr>
        <tr>
          <th>30</th>
          <td>Sample16_Day2_Bio1_1</td>
          <td>Treatment</td>
          <td>3</td>
          <td>16</td>
        </tr>
        <tr>
          <th>31</th>
          <td>Sample16_Day2_Bio1_2</td>
          <td>Treatment</td>
          <td>3</td>
          <td>16</td>
        </tr>
        <tr>
          <th>32</th>
          <td>Sample17_Day2_Bio2_1</td>
          <td>Treatment</td>
          <td>3</td>
          <td>17</td>
        </tr>
        <tr>
          <th>33</th>
          <td>Sample17_Day2_Bio2_2</td>
          <td>Treatment</td>
          <td>3</td>
          <td>17</td>
        </tr>
        <tr>
          <th>34</th>
          <td>Sample18_Day2_Bio3_1</td>
          <td>Treatment</td>
          <td>3</td>
          <td>18</td>
        </tr>
        <tr>
          <th>35</th>
          <td>Sample18_Day2_Bio3_2</td>
          <td>Treatment</td>
          <td>3</td>
          <td>18</td>
        </tr>
        <tr>
          <th>36</th>
          <td>Sample19_Day3_Bio1_1</td>
          <td>Treatment</td>
          <td>5</td>
          <td>19</td>
        </tr>
        <tr>
          <th>37</th>
          <td>Sample19_Day3_Bio1_2</td>
          <td>Treatment</td>
          <td>5</td>
          <td>19</td>
        </tr>
        <tr>
          <th>38</th>
          <td>Sample20_Day3_Bio2_1</td>
          <td>Treatment</td>
          <td>5</td>
          <td>20</td>
        </tr>
        <tr>
          <th>39</th>
          <td>Sample20_Day3_Bio2_2</td>
          <td>Treatment</td>
          <td>5</td>
          <td>20</td>
        </tr>
        <tr>
          <th>40</th>
          <td>Sample21_Day3_Bio3_1</td>
          <td>Treatment</td>
          <td>5</td>
          <td>21</td>
        </tr>
        <tr>
          <th>41</th>
          <td>Sample21_Day3_Bio3_2</td>
          <td>Treatment</td>
          <td>5</td>
          <td>21</td>
        </tr>
        <tr>
          <th>42</th>
          <td>Sample22_Day4_Bio1_1</td>
          <td>Treatment</td>
          <td>7</td>
          <td>22</td>
        </tr>
        <tr>
          <th>43</th>
          <td>Sample22_Day4_Bio1_2</td>
          <td>Treatment</td>
          <td>7</td>
          <td>22</td>
        </tr>
        <tr>
          <th>44</th>
          <td>Sample23_Day4_Bio2_1</td>
          <td>Treatment</td>
          <td>7</td>
          <td>23</td>
        </tr>
        <tr>
          <th>45</th>
          <td>Sample23_Day4_Bio2_2</td>
          <td>Treatment</td>
          <td>7</td>
          <td>23</td>
        </tr>
        <tr>
          <th>46</th>
          <td>Sample24_Day5_Bio3_1</td>
          <td>Treatment</td>
          <td>7</td>
          <td>24</td>
        </tr>
        <tr>
          <th>47</th>
          <td>Sample24_Day5_Bio3_2</td>
          <td>Treatment</td>
          <td>7</td>
          <td>24</td>
        </tr>
      </tbody>
    </table>
    </div>



Pdata Example 3
~~~~~~~~~~~~~~~

It’s important to note that in some cases researchers may employ
independent or related sampling over time. Independent sampling involves
evaluating different individuals over time, while related sampling
entails assessing the same individuals repeatedly. As OmicScope assumes
independent sampling by default, it’s essential to add a fifth column
labeled “Individual” if the experimental design involves related
sampling. This column associates each sample with its respective
individual number.

Using the example provided, when conducting related sampling, the user
should add the ``Individual`` column to associate each biological sample
with the corresponding individual.

.. code:: python

    import pandas as pd
    pdata = pd.read_excel('../../tests/data/proteins/longitudinal_pdata.xlsx', sheet_name=1)
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
          <th>TimeCourse</th>
          <th>Individual</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>Sample1_Day1_Bio1_1</td>
          <td>Control</td>
          <td>1</td>
          <td>1</td>
          <td>1</td>
        </tr>
        <tr>
          <th>1</th>
          <td>Sample1_Day1_Bio1_2</td>
          <td>Control</td>
          <td>1</td>
          <td>1</td>
          <td>1</td>
        </tr>
        <tr>
          <th>2</th>
          <td>Sample2_Day1_Bio2_1</td>
          <td>Control</td>
          <td>2</td>
          <td>1</td>
          <td>2</td>
        </tr>
        <tr>
          <th>3</th>
          <td>Sample2_Day1_Bio2_2</td>
          <td>Control</td>
          <td>2</td>
          <td>1</td>
          <td>2</td>
        </tr>
        <tr>
          <th>4</th>
          <td>Sample3_Day1_Bio3_1</td>
          <td>Control</td>
          <td>3</td>
          <td>1</td>
          <td>3</td>
        </tr>
        <tr>
          <th>5</th>
          <td>Sample3_Day1_Bio3_2</td>
          <td>Control</td>
          <td>3</td>
          <td>1</td>
          <td>3</td>
        </tr>
        <tr>
          <th>6</th>
          <td>Sample4_Day2_Bio1_1</td>
          <td>Control</td>
          <td>4</td>
          <td>3</td>
          <td>1</td>
        </tr>
        <tr>
          <th>7</th>
          <td>Sample4_Day2_Bio1_2</td>
          <td>Control</td>
          <td>4</td>
          <td>3</td>
          <td>1</td>
        </tr>
        <tr>
          <th>8</th>
          <td>Sample5_Day2_Bio2_1</td>
          <td>Control</td>
          <td>5</td>
          <td>3</td>
          <td>2</td>
        </tr>
        <tr>
          <th>9</th>
          <td>Sample5_Day2_Bio2_2</td>
          <td>Control</td>
          <td>5</td>
          <td>3</td>
          <td>2</td>
        </tr>
        <tr>
          <th>10</th>
          <td>Sample6_Day2_Bio3_1</td>
          <td>Control</td>
          <td>6</td>
          <td>3</td>
          <td>3</td>
        </tr>
        <tr>
          <th>11</th>
          <td>Sample6_Day2_Bio3_2</td>
          <td>Control</td>
          <td>6</td>
          <td>3</td>
          <td>3</td>
        </tr>
        <tr>
          <th>12</th>
          <td>Sample7_Day3_Bio1_1</td>
          <td>Control</td>
          <td>7</td>
          <td>5</td>
          <td>1</td>
        </tr>
        <tr>
          <th>13</th>
          <td>Sample7_Day3_Bio1_2</td>
          <td>Control</td>
          <td>7</td>
          <td>5</td>
          <td>1</td>
        </tr>
        <tr>
          <th>14</th>
          <td>Sample8_Day3_Bio2_1</td>
          <td>Control</td>
          <td>8</td>
          <td>5</td>
          <td>2</td>
        </tr>
        <tr>
          <th>15</th>
          <td>Sample8_Day3_Bio2_2</td>
          <td>Control</td>
          <td>8</td>
          <td>5</td>
          <td>2</td>
        </tr>
        <tr>
          <th>16</th>
          <td>Sample9_Day3_Bio3_1</td>
          <td>Control</td>
          <td>9</td>
          <td>5</td>
          <td>3</td>
        </tr>
        <tr>
          <th>17</th>
          <td>Sample9_Day3_Bio3_2</td>
          <td>Control</td>
          <td>9</td>
          <td>5</td>
          <td>3</td>
        </tr>
        <tr>
          <th>18</th>
          <td>Sample10_Day4_Bio1_1</td>
          <td>Control</td>
          <td>10</td>
          <td>7</td>
          <td>1</td>
        </tr>
        <tr>
          <th>19</th>
          <td>Sample10_Day4_Bio1_2</td>
          <td>Control</td>
          <td>10</td>
          <td>7</td>
          <td>1</td>
        </tr>
        <tr>
          <th>20</th>
          <td>Sample11_Day4_Bio2_1</td>
          <td>Control</td>
          <td>11</td>
          <td>7</td>
          <td>2</td>
        </tr>
        <tr>
          <th>21</th>
          <td>Sample11_Day4_Bio2_2</td>
          <td>Control</td>
          <td>11</td>
          <td>7</td>
          <td>2</td>
        </tr>
        <tr>
          <th>22</th>
          <td>Sample12_Day5_Bio3_1</td>
          <td>Control</td>
          <td>12</td>
          <td>7</td>
          <td>3</td>
        </tr>
        <tr>
          <th>23</th>
          <td>Sample12_Day5_Bio3_2</td>
          <td>Control</td>
          <td>12</td>
          <td>7</td>
          <td>3</td>
        </tr>
        <tr>
          <th>24</th>
          <td>Sample13_Day1_Bio1_1</td>
          <td>Treatment</td>
          <td>13</td>
          <td>1</td>
          <td>4</td>
        </tr>
        <tr>
          <th>25</th>
          <td>Sample13_Day1_Bio1_2</td>
          <td>Treatment</td>
          <td>13</td>
          <td>1</td>
          <td>4</td>
        </tr>
        <tr>
          <th>26</th>
          <td>Sample14_Day1_Bio2_1</td>
          <td>Treatment</td>
          <td>14</td>
          <td>1</td>
          <td>5</td>
        </tr>
        <tr>
          <th>27</th>
          <td>Sample14_Day1_Bio2_2</td>
          <td>Treatment</td>
          <td>14</td>
          <td>1</td>
          <td>5</td>
        </tr>
        <tr>
          <th>28</th>
          <td>Sample15_Day1_Bio3_1</td>
          <td>Treatment</td>
          <td>15</td>
          <td>1</td>
          <td>6</td>
        </tr>
        <tr>
          <th>29</th>
          <td>Sample15_Day1_Bio3_2</td>
          <td>Treatment</td>
          <td>15</td>
          <td>1</td>
          <td>6</td>
        </tr>
        <tr>
          <th>30</th>
          <td>Sample16_Day2_Bio1_1</td>
          <td>Treatment</td>
          <td>16</td>
          <td>3</td>
          <td>4</td>
        </tr>
        <tr>
          <th>31</th>
          <td>Sample16_Day2_Bio1_2</td>
          <td>Treatment</td>
          <td>16</td>
          <td>3</td>
          <td>4</td>
        </tr>
        <tr>
          <th>32</th>
          <td>Sample17_Day2_Bio2_1</td>
          <td>Treatment</td>
          <td>17</td>
          <td>3</td>
          <td>5</td>
        </tr>
        <tr>
          <th>33</th>
          <td>Sample17_Day2_Bio2_2</td>
          <td>Treatment</td>
          <td>17</td>
          <td>3</td>
          <td>5</td>
        </tr>
        <tr>
          <th>34</th>
          <td>Sample18_Day2_Bio3_1</td>
          <td>Treatment</td>
          <td>18</td>
          <td>3</td>
          <td>6</td>
        </tr>
        <tr>
          <th>35</th>
          <td>Sample18_Day2_Bio3_2</td>
          <td>Treatment</td>
          <td>18</td>
          <td>3</td>
          <td>6</td>
        </tr>
        <tr>
          <th>36</th>
          <td>Sample19_Day3_Bio1_1</td>
          <td>Treatment</td>
          <td>19</td>
          <td>5</td>
          <td>4</td>
        </tr>
        <tr>
          <th>37</th>
          <td>Sample19_Day3_Bio1_2</td>
          <td>Treatment</td>
          <td>19</td>
          <td>5</td>
          <td>4</td>
        </tr>
        <tr>
          <th>38</th>
          <td>Sample20_Day3_Bio2_1</td>
          <td>Treatment</td>
          <td>20</td>
          <td>5</td>
          <td>5</td>
        </tr>
        <tr>
          <th>39</th>
          <td>Sample20_Day3_Bio2_2</td>
          <td>Treatment</td>
          <td>20</td>
          <td>5</td>
          <td>5</td>
        </tr>
        <tr>
          <th>40</th>
          <td>Sample21_Day3_Bio3_1</td>
          <td>Treatment</td>
          <td>21</td>
          <td>5</td>
          <td>6</td>
        </tr>
        <tr>
          <th>41</th>
          <td>Sample21_Day3_Bio3_2</td>
          <td>Treatment</td>
          <td>21</td>
          <td>5</td>
          <td>6</td>
        </tr>
        <tr>
          <th>42</th>
          <td>Sample22_Day4_Bio1_1</td>
          <td>Treatment</td>
          <td>22</td>
          <td>7</td>
          <td>4</td>
        </tr>
        <tr>
          <th>43</th>
          <td>Sample22_Day4_Bio1_2</td>
          <td>Treatment</td>
          <td>22</td>
          <td>7</td>
          <td>4</td>
        </tr>
        <tr>
          <th>44</th>
          <td>Sample23_Day4_Bio2_1</td>
          <td>Treatment</td>
          <td>23</td>
          <td>7</td>
          <td>5</td>
        </tr>
        <tr>
          <th>45</th>
          <td>Sample23_Day4_Bio2_2</td>
          <td>Treatment</td>
          <td>23</td>
          <td>7</td>
          <td>5</td>
        </tr>
        <tr>
          <th>46</th>
          <td>Sample24_Day5_Bio3_1</td>
          <td>Treatment</td>
          <td>24</td>
          <td>7</td>
          <td>6</td>
        </tr>
        <tr>
          <th>47</th>
          <td>Sample24_Day5_Bio3_2</td>
          <td>Treatment</td>
          <td>24</td>
          <td>7</td>
          <td>6</td>
        </tr>
      </tbody>
    </table>
    </div>


Longitudinal example
~~~~~~~~~~~~~~~~~~~~

A complete longitudinal dataset in the General format is available at
``tests/data/proteins/longitudinal_example.xlsx`` (also offered in the
OmicScope app, under *Pdata*, as "Download Longitudinal example"). It
mimics a time-course experiment in which cells from 8 donors were either
left untreated (``Control``, donors 1-4) or treated (``Treated``, donors
5-8) and sampled at 0, 6, 12, 24 and 48 hours, giving 40 samples and
2,000 human proteins.

The workbook has the three sheets read by the General method:

- **assay**: protein abundances, one column per sample. Values below
  the detection limit are left empty.
- **rdata**: ``Accession`` and ``Description`` (with ``GN=``) of each
  protein, plus peptide counts.
- **pdata**: ``Sample``, ``Condition``, ``Biological``, ``TimeCourse``
  (in hours) and ``Individual``. Each donor is measured at every time
  point, so the ``Individual`` column tells OmicScope that the samples
  are related.

A fourth sheet, ``truth``, is ignored by OmicScope. The data are
simulated, and this sheet lists the profile each protein was given, so
you can compare it with the results. Most proteins do not change. The
rest rise briefly in treated cells (peak at 6-12 h), fall and stay down,
rise only late (24-48 h), move in opposite directions in the two groups,
change over time in the same way in both groups, or are higher in treated
cells at every time point.

.. code:: python

    import omicscope as omics

    longitudinal = omics.OmicScope('../../tests/data/proteins/longitudinal_example.xlsx',
                                   Method='General',
                                   ExperimentalDesign='longitudinal',
                                   degrees_of_freedom=3)

The ``pAdjusted`` column answers the main question of a between-group
design: does the time profile of the protein differ between Control and
Treated? With ``degrees_of_freedom=3`` it finds about 99% of the proteins
whose profiles differ between the groups, and about 5% of the proteins it
reports are false positives, matching the 5% cutoff. With the default of
2 degrees of freedom, the spline is too stiff to follow the short peak at
6-12 h, so fewer than half of those proteins are found. Proteins that are
only higher in Treated, with no change over time, are not reported by
``pAdjusted``. Proteins that change in the same way in both groups are not
reported either. They appear in ``pTime`` / ``pAdjusted_time`` instead,
which test the time trend shared by the groups.

The same data are also available as two files: the quantitative data
(sheets assay and rdata only) and a separate pdata file. The pdata file
then provides the third sheet:

.. code:: python

    longitudinal = omics.OmicScope('../../tests/data/proteins/longitudinal_example_assay.xlsx',
                                   pdata='../../tests/data/proteins/longitudinal_example_pdata.xlsx',
                                   Method='General',
                                   ExperimentalDesign='longitudinal',
                                   degrees_of_freedom=3)

Things to keep in mind when you build your own longitudinal pdata:

- ``TimeCourse`` must be numeric (for example, hours or days), because
  OmicScope fits a curve over time.
- Give each individual a unique number across all conditions (here 1-8,
  not 1-4 in each group). If the numbers restart in each group,
  OmicScope assumes that Control donor 1 and Treated donor 1 are the same
  person. In this example that makes OmicScope miss about a third of the
  regulated proteins.
- Leave out the ``Individual`` column only if every sample comes from a
  different subject. If you remove it from this example, OmicScope
  ignores the differences between donors and finds only a third of the
  regulated proteins.
- Do not add extra columns to a longitudinal pdata. OmicScope adds every
  column other than ``Sample``, ``Biological`` and ``Individual`` to the
  statistical model.
- Remember to set ``ExperimentalDesign='longitudinal'``. With the default
  static design, OmicScope compares the groups as if all 40 samples came
  from different people, and many unchanged proteins are reported as
  regulated.

