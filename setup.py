#!/usr/bin/env python
# -*- encoding: utf-8 -*-

import io
import re
from glob import glob
from os.path import basename
from os.path import dirname
from os.path import join
from os.path import splitext

from setuptools import find_packages
from setuptools import setup


def read(*names, **kwargs):
    with io.open(join(dirname(__file__), *names), encoding=kwargs.get('encoding', 'utf8')) as fh:
        return fh.read()


setup(
    name='omicscope',
    version='1.0.0',
    license='MIT',
    description='OmicScope: from DIA proteomics to Systems Biology visualisation',
    long_description='{}\n{}'.format(
        re.compile('^.. start-badges.*^.. end-badges', re.M | re.S).sub('', read('README.rst')),
        re.sub(':[a-z]+:`~?(.*?)`', r'``\1``', read('CHANGELOG.rst')),
    ),
    author='Guilherme Reis-de-Oliveira',
    author_email='guioliveirareis@gmail.com',
    url='https://github.com/guireiso/omicscope',
    packages=find_packages('src'),
    package_dir={'': 'src'},
    py_modules=[splitext(basename(path))[0] for path in glob('src/*.py')],
    include_package_data=True,
    zip_safe=False,
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Operating System :: Unix',
        'Operating System :: POSIX',
        'Operating System :: Microsoft :: Windows',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Topic :: Utilities',
    ],
    project_urls={
        'Documentation': 'https://omicscope.readthedocs.io/',
        'Changelog': 'https://omicscope.readthedocs.io/en/latest/changelog.html',
        'Issue Tracker': 'https://github.com/guireiso/omicscope/issues',
    },
    keywords=[
        'Proteomics',
        'Omics',
        'Differential Expression',
        'Systems Biology'
    ],
    python_requires='>=3.7',
    install_requires=[
        "setuptools>=30.3.0",
        "wheel",
        "adjustText==0.7.3",
        "altair==4.2.0",
        "arrow==1.2.3",
        "asttokens==2.1.0",
        "attrs==22.1.0",
        "backcall==0.2.0",
        "binaryornot==0.4.4",
        "certifi==2022.9.24",
        "chardet==5.0.0",
        "charset-normalizer==2.1.1",
        "click==8.1.3",
        "colorama==0.4.6",
        "contourpy==1.0.6",
        "cookiecutter==2.1.1",
        "cycler==0.11.0",
        "debugpy==1.6.3",
        "decorator==5.1.1",
        "distlib==0.3.6",
        "entrypoints==0.4",
        "et-xmlfile==1.1.0",
        "executing==1.2.0",
        "filelock==3.8.0",
        "fonttools==4.38.0",
        "gseapy==1.0.3",
        "idna==3.4",
        "iniconfig==1.1.1",
        "ipykernel==6.17.1",
        "ipykernel",
        "ipython==8.6.0",
        "jedi==0.18.1",
        "Jinja2==3.1.2",
        "jinja2-time==0.2.0",
        "joblib==1.2.0",
        "jsonschema==4.17.0",
        "jupyter_client==7.4.4",
        "jupyter_core==5.0.0",
        "kiwisolver==1.4.4",
        "kneed==0.8.1",
        "MarkupSafe==2.1.1",
        "matplotlib==3.6.2",
        "matplotlib-inline==0.1.6",
        "nest-asyncio==1.5.6",
        "networkx==2.8.8",
        "numpy==1.23.4",
        "openpyxl==3.0.10",
        "packaging==21.3",
        "pandas==1.5.1",
        "parso==0.8.3",
        "patsy==0.5.3",
        "pickleshare==0.7.5",
        "Pillow==9.3.0",
        "platformdirs==2.5.3",
        "pluggy==1.0.0",
        "prompt-toolkit==3.0.32",
        "psutil==5.9.4",
        "pure-eval==0.2.2",
        "py==1.11.0",
        "Pygments==2.13.0",
        "pyparsing==3.0.9",
        "pyrsistent==0.19.2",
        "pytest==7.2.0",
        "python-dateutil==2.8.2",
        "python-slugify==6.1.2",
        "pytz==2022.6",
        "pywin32==305",
        "pyvis",
        "PyYAML==6.0",
        "pyzmq==24.0.1",
        "requests==2.28.1",
        "scikit-learn==1.1.3",
        "scipy==1.9.3",
        "seaborn==0.12.1",
        "six==1.16.0",
        "stack-data==0.6.0",
        "statsmodels==0.13.5",
        "text-unidecode==1.3",
        "threadpoolctl==3.1.0",
        "toolz==0.12.0",
        "tornado==6.2",
        "tox==3.27.0",
        "traitlets==5.5.0",
        "UpSetPlot==0.6.1",
        "urllib3==1.26.12",
        "virtualenv==20.16.6",
        "wcwidth==0.2.5",
        "xlrd==2.0.1"
    ],
    extras_require={
        # eg:
        #   'rst': ['docutils>=0.11'],
        #   ':python_version=="2.6"': ['argparse'],
    },
    entry_points={
        'console_scripts': [
            'omicscope = omicscope.cli:main',
        ]
    },
)
