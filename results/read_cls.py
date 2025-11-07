import pandas as pd

def read_cls(filename, sample_names = None):
    # read a phenotype file in CLS format and retuns the phenotype as a Pandas list

    with open(filename) as cls_file:
        line = cls_file.readlines()
        if line[0].strip() == '#numeric':
            name = line[1].strip().split(sep='#')[1]
            phen = line[2].strip().split(sep=' ')
            phen = pd.Series(phen, name = name)
            phen = phen.astype('float64')

        else: # categorical
            ser = filename.split("/")
            name = ser[-1].split(".")[0]
            phen = line[2].strip().split(sep=' ')
            phen = pd.Series(phen, name = name)

        if sample_names is not None:
            phen.index = sample_names

    return(phen)
