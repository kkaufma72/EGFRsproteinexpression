# Analysisexample
estimating the degree of association of genes with EGFRs protein expression

import itertools
import os
import warnings
import copy
                                                                                                                                
import numpy as np
import pandas as pd
import plotly as pl
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
#import yaml
import pickle
from scipy.cluster.hierarchy import fcluster, linkage, to_tree, leaders, leaves_list, ClusterNode, fcluster, ward
from scipy.spatial.distance import pdist, squareform
import sklearn
import seaborn as sns

import matplotlib as mpl
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

warnings.filterwarnings("ignore")

pl.offline.init_notebook_mode(connected=True)

from plotly.offline import init_notebook_mode
init_notebook_mode(connected = True)
np.random.seed(7678)

import matplotlib.image as mpi    
import datetime

import time

from IPython.display import display, HTML

display(HTML(data="""
<style>
    div#notebook-container    { width: 95%; }
    div#menubar-container     { width: 65%; }
    div#maintoolbar-container { width: 99%; }
</style>
""")) 

# pd.set_option('display.max_columns', None)
# pd.set_option("max_rows", None)
