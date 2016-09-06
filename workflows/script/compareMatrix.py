#!/usr/bin/env python

"""
using seaborn jointplot to plot kde on log probability
cmd:

./compareMatrix.py Ke_tadmodel30.hdf5 Nan_tadmodel0.hdf5 0.05 "\$Model_{old}\$" "\$Model_{new}\$" modelComp.pdf

"""
from __future__ import print_function
import matplotlib
matplotlib.use('Agg')
import numpy as np
import pandas as pd
import seaborn as sns
import alab.matrix
import sys
sns.set(style="white")

if len(sys.argv) >= 6:
    modelfile1 = sys.argv[1]
    modelfile2 = sys.argv[2]
    cutoff     = float(sys.argv[3])
    name1      = sys.argv[4]
    name2      = sys.argv[5]
    outputname    = sys.argv[6]
else:
    print(sys.argv[0]+"<modelfile1> <modelfile2> <cutoff> <x-axis name> <y-axis name> <output>")
    sys.exit(0)


#modelfile1 = "Ke_tadmodel30.hdf5"
#modelfile2 = "Nan_tadmodel0.hdf5"
#cutoff = 0.05
#name1  = "$Model_{old}$"
#name2  = "$Model_{new}$"
#outputname = "complog0.05.pdf"

mx = alab.matrix.contactmatrix(modelfile1)
my = alab.matrix.contactmatrix(modelfile2)

n = len(mx)
x = mx.matrix[np.triu_indices(n,1)]
y = my.matrix[np.triu_indices(n,1)]

x1 = pd.Series(np.log(x[(x > cutoff) & (y > cutoff)]),name = name1)
y1 = pd.Series(np.log(y[(x > cutoff) & (y > cutoff)]),name = name2)

g = sns.jointplot(x1,y1,kind="kde",size=7,space=0)

g.savefig(outputname)
