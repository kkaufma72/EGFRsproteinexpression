from numpy import asarray, exp, finfo, isnan, log, log2, nan, sign, sqrt, unique, power, e, pi

#from rpy2.robjects.numpy2ri import numpy2ri
# import rpy2.robjects.numpy2ri as numpy2ri
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
from rpy2.robjects import r, pandas2ri

import numpy as np
from scipy.stats import pearsonr

eps = finfo(float).eps

# ro.conversion.py2ri = numpy2ri

pandas2ri.activate()

mass = importr("MASS")
    
def universal_metric(x, y, n_grid=24, delta = 1.0):
    
    pearson_correlation = pearsonr(x, y)[0]

    if isnan(pearson_correlation) or unique(x).size == 1 or unique(y).size == 1:

        return nan

    else:

        pearson_correlation_abs = abs(pearson_correlation)

       # xr = pandas2ri.py2ri(x)
       # yr = pandas2ri.py2ri(y)

        #print('IC: sum x:{} sum y:{}'.format(np.sum(x), np.sum(y)))
        
        #print('bandwidth x: {}'.format(mass.bcv(x)[0]))
        bandwidth_x = delta * mass.bcv(x)[0] * (1 - pearson_correlation_abs * 0.75)

        #print('bandwidth y: {}'.format(mass.bcv(y)[0]))
        bandwidth_y = delta * mass.bcv(y)[0] * (1 - pearson_correlation_abs * 0.75)

        #print('bandwidth... done')
        
        fxy = (
            asarray(
                mass.kde2d(
                    x, y, asarray((bandwidth_x, bandwidth_y)), n=asarray((n_grid,))
                )[2]
            )
            + eps
        )

        pxy = fxy/(fxy.sum())
        px = pxy.sum(axis=1)
        px = px/px.sum()
        py = pxy.sum(axis=0)
        py = py/py.sum()
        hxy = - (pxy * log2(pxy)).sum() 
        hx = - (px * log2(px)).sum() 
        hy = - (py * log2(py)).sum()
        mi = hx + hy - hxy
        # This is the universal metric (see e.g. Kraskov et al https://arxiv.org/pdf/q-bio/0311039.pdf
                                  # Li et al https://arxiv.org/abs/cs/0111054 and Bennett et al https://arxiv.org/abs/1006.3520
        nmi = mi/np.max([hx, hy])

        # The universal metric is normalized is a similar way to the original Linfoot's Information Coefficient
        # (https://www.sciencedirect.com/science/article/pii/S001999585790116X)
        # using the entropy and mutual information for a Gaussian distribution (see e.g. Example 8.5.1 in
        # Elements of Information Theory 2nd ed - T. Cover, J. Thomas Wiley, 2006)
        
        ent_gauss_x = 0.5 * log2(2 * pi * e * np.std(x)**2)
        ent_gauss_y = 0.5 * log2(2 * pi * e * np.std(y)**2)

        IC_nmi = sign(pearson_correlation) * sqrt(1 - power(2.0, -2.0 * nmi * np.max([ent_gauss_x, ent_gauss_y])))

        # The problem with this normalization is that in makes the IC dependent of the sigma of the input data and
        # this can produce nan's for vectors with low entropy.

        # print('sx: {} sy: {} mi: {}: nmi: {} IC_nmi: {}'.format(np.std(x), np.std(y), mi, nmi, IC_nmi))
        
        return IC_nmi

            

            
