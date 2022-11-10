========
Overview
========

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

OmicScope: from DIA proteomics to Systems Biology visualisation

* Free software: MIT license

Installation
============

::

    pip install omicscope

You can also install the in-development version with::

    pip install https://github.com/guireiso/omicscope/archive/main.zip


Documentation
=============


https://omicscope.readthedocs.io/


Development
===========

To run all the tests run::

    tox

Note, to combine the coverage data from all the tox environments run:

.. list-table::
    :widths: 10 90
    :stub-columns: 1

    - - Windows
      - ::

            set PYTEST_ADDOPTS=--cov-append
            tox

    - - Other
      - ::

            PYTEST_ADDOPTS=--cov-append tox
