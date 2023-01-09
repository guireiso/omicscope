========
Welcome to OmicScope!
========

OmicScope: from quantitative proteomics to Systems Biology!
---------------------------------------------------------------

.. start-badges

.. list-table::
    :stub-columns: 1

    * - docs
      - |docs|
    * - tests
      - | |github-actions| |requires|
        | |codecov|
    * - package
      - | |version| |wheel| |supported-versions| |supported-implementations|
        | |commits-since|
.. |docs| image:: https://readthedocs.org/projects/omicscope/badge/?style=flat
    :target: https://omicscope.readthedocs.io/
    :alt: Documentation Status

.. |github-actions| image:: https://github.com/guireiso/omicscope/actions/workflows/github-actions.yml/badge.svg
    :alt: GitHub Actions Build Status
    :target: https://github.com/guireiso/omicscope/actions

.. |requires| image:: https://requires.io/github/guireiso/omicscope/requirements.svg?branch=main
    :alt: Requirements Status
    :target: https://requires.io/github/guireiso/omicscope/requirements/?branch=main

.. |codecov| image:: https://codecov.io/gh/guireiso/omicscope/branch/main/graphs/badge.svg?branch=main
    :alt: Coverage Status
    :target: https://codecov.io/github/guireiso/omicscope

.. |version| image:: https://img.shields.io/pypi/v/omicscope.svg
    :alt: PyPI Package latest release
    :target: https://pypi.org/project/omicscope

.. |wheel| image:: https://img.shields.io/pypi/wheel/omicscope.svg
    :alt: PyPI Wheel
    :target: https://pypi.org/project/omicscope

.. |supported-versions| image:: https://img.shields.io/pypi/pyversions/omicscope.svg
    :alt: Supported versions
    :target: https://pypi.org/project/omicscope

.. |supported-implementations| image:: https://img.shields.io/pypi/implementation/omicscope.svg
    :alt: Supported implementations
    :target: https://pypi.org/project/omicscope

.. |commits-since| image:: https://img.shields.io/github/commits-since/guireiso/omicscope/v1.0.0.svg
    :alt: Commits since latest release
    :target: https://github.com/guireiso/omicscope/compare/v1.0.0...main



.. end-badges

Introduction
==============

*OmicScope* is a comprehensive workflow designed to analyse and provide insights regarding quantitative proteomics data. Todate, *OmicScope* works with data generated from [Progenesis QI for Proteomics](https://www.nonlinear.com/progenesis/qi-for-proteomics/), [MaxQuant](https://www.maxquant.org/), and [PatternLab V](http://www.patternlabforproteomics.org/). Additionally, a fourth generical input can be used for further analysis, allowing users to run software coming from different plataforms (including transcriptomics and metabolomics).

OmicScope can perform differential expression analysis in both static and longitudinal experimental designs. For static experiments, proteins differentially regulated are determined by t-tests (2 group comparison) or One-way ANOVA  (>2 group comparison); while for longitudinal analysis, OmicScope perform the pipeline suggested by [Storey, 2005](https://www.pnas.org/doi/10.1073/pnas.0504609102). 


Citation
============
::

    Reis-de-Oliveira G, Martins-de-Souza, G. OmicScope: from quantitative proteomics to Systems Biology!, 
    Journal, 2023;


Installation
============

::

    pip install omicscope

You can also install the in-development version with::

    pip install https://github.com/guireiso/omicscope/archive/main.zip


Documentation
=============

https://omicscope.readthedocs.io/