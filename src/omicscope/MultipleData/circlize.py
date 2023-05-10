import io
from copy import copy

import pandas as pd


def deps(self, pvalue=0.05):
    data = copy(self)
    groups = data.groups
    regulation = data.original
    regulation = [x[x[data.pvalue] < pvalue] for x in regulation]
    regulation = [x[['gene_name', 'log2(fc)']] for x in regulation]
    regulation = [x.rename(columns={'log2(fc)': y}) for x, y in zip(regulation, groups)]
    regulation = pd.concat(regulation)
    regulation = regulation.groupby('gene_name').sum().reset_index()
    return regulation


def enrichment_filtering(self, Term_enriched):
    data = copy(self.enrichment)
    data = [x[x['Term'].str.contains(Term_enriched)] for x in data]
    data = pd.concat(data)
    deps = list(data['Genes'])
    deps = [x.replace("'", '') for x in deps]
    deps = [x.replace("[", '') for x in deps]
    deps = [x.replace("]", '') for x in deps]
    deps = [x.split(', ') for x in deps]
    deps = sum(deps, [])
    deps = list(set(deps))
    return deps


def deps_matrix(df):
    deps = copy(df)
    deps[deps > 0] = 1
    deps[deps < 0] = 0.5
    deps[deps.isna()] = 0
    return deps


def color_matrix(df, colors):
    colmat = copy(df)
    for i, z in zip(range(0, len(colmat.columns)), colors):
        colmat.iloc[:, i] = z
    return colmat


def circlize(matrix, colmat, colors, labels, width=3000, height=3000,
             save=None, vector=True):
    import os
    import inspect
    import subprocess
    import pickle
    import shutil
    import matplotlib.pyplot as plt

    circlize_path = os.path.abspath(inspect.getfile(circlize))
    circlize_path = circlize_path.removesuffix('circlize.py')
    wdir = os.getcwd()
    wdir = wdir + '\\'
    if save is None:
        vector = False

    dictionary = {"matrix": matrix,
                  "colmat": colmat,
                  "colors": list(colors),
                  "labels": labels,
                  "width": width,
                  "height": height,
                  "save": wdir,
                  "vector": vector
                  }
    pickle.dump(dictionary, open(wdir+'circlize.p', 'wb'))
    subprocess.run('Rscript '+circlize_path+'circlize.R', shell=True,
                   text=True, capture_output=True)

    if save is None:
        figure = plt.imread(wdir+'_my_plot.png')
        fig, ax = plt.subplots()
        im = ax.imshow(figure)
        plt.axis('off')
        plt.show()
        os.remove(wdir+'_my_plot.png')
    else:
        if vector is True:
            extension = '.svg'
        else:
            extension = '.png'
        shutil.move(wdir+'_my_plot'+extension, save+'_my_plot'+extension)
        print('Your file was saved in ' + save)
    os.remove(wdir+'circlize.p')
