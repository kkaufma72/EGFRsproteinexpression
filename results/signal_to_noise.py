import numpy as np

def signal_to_noise(x, y):
    vals = x.unique()
    c0 = x.index[x == vals[0]]
    c1 = x.index[x == vals[1]]
    m0 = y.loc[c0].mean()
    m1 = y.loc[c1].mean()
    s0 = y.loc[c0].std()
    s1 = y.loc[c1].std()
    if s0 == 0:
        if m0 == 0:
            s0 = 0.2
        else:
            s0 = 0.2 * np.abs(m0)
    if s1 == 0:
        if m1 == 0:
            s1 = 0.2
        else:
            s1 = 0.2 * np.abs(m1)
        
    return (m0 - m1)/(s0 + s1)
