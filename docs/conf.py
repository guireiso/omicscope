from __future__ import unicode_literals

import os

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.coverage',
    'sphinx.ext.doctest',
    'sphinx.ext.extlinks',
    'sphinx.ext.ifconfig',
    'sphinx.ext.napoleon',
    'sphinx.ext.todo',
    'sphinx.ext.viewcode',
    'sphinx-rtd-theme'
]
source_suffix = '.rst'
master_doc = 'index'
project = 'OmicScope'
year = '2022'
author = 'Guilherme Reis-de-Oliveira'
copyright = '{0}, {1}'.format(year, author)
version = release = '1.3.11'

pygments_style = 'trac'
templates_path = ['.']
extlinks = {
    'issue': ('https://github.com/guireiso/omicscope/issues/%s', '#'),
    'pr': ('https://github.com/guireiso/omicscope/pull/%s', 'PR #'),
}
# on_rtd is whether we are on readthedocs.org
on_rtd = os.environ.get('READTHEDOCS', None) == 'True'

if not on_rtd:  # only set the theme if we're building docs locally
    html_theme = 'sphinx_rtd_theme'

html_use_smartypants = True
html_last_updated_fmt = '%b %d, %Y'
html_split_index = False
html_theme = 'sphinx_rtd_theme'
html_sidebars = {
    '**': ['searchbox.html', 'globaltoc.html', 'sourcelink.html'],
}
html_theme_options = {
    'logo_only': True,
    'display_version': True,
}
html_theme_path = ["."]
html_logo = "logo.png"

html_short_title = '%s-%s' % (project, version)


napoleon_use_ivar = True
napoleon_use_rtype = False
napoleon_use_param = False
