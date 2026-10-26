============
Installation
============

OmicScope needs Python 3.12 or newer. Install it from PyPI::

    pip install omicscope

or add it to a project managed with `uv <https://docs.astral.sh/uv/>`__::

    uv add omicscope

To try the development version::

    pip install https://github.com/guireiso/omicscope/archive/main.zip

All dependencies ship as binary wheels for Windows, macOS and Linux, so no compiler is needed.

Development setup
-----------------

Clone the repository and let uv create the environment from the lock file::

    git clone https://github.com/guireiso/omicscope
    cd omicscope
    uv sync
    uv run pytest

The documentation is built with::

    uv sync --group docs
    uv run sphinx-build docs docs/_build/html

Web application
---------------

The OmicScope web application (https://omicscope.ib.unicamp.br) runs the same package in the browser, with
interactive figures, and needs no installation.
