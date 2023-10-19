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
    version='1.3.11',
    license='MIT',
    description='OmicScope: from quantitative proteomics to systems biology.',
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
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
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
        "adjustText==0.8",
        "altair==4.2.0",
        "flake8==6.0.0",
        "gseapy==1.0.3",
        "kneed==0.8.1",
        "matplotlib==3.6.2",
        "matplotlib-inline==0.1.6",
        "networkx==2.8.8",
        "numpy==1.23.4",
        "openpyxl==3.0.10",
        "packaging==21.3",
        "pandas==1.5.1",
        "Pillow==9.3.0",
        "pkginfo==1.9.4",
        "pycirclize==0.4.0",
        "scikit-learn==1.1.3",
        "scipy==1.9.3",
        "seaborn==0.12.1",
        "statsmodels==0.13.5",
        "tox==3.27.0",
        "UpSetPlot==0.6.1",
        "xlrd==2.0.1",
        "zipp==3.11.0",
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
