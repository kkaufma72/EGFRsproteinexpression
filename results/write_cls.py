import pandas as pd

def write_cls(series, filename, phenotype_type = 'categorical'):
    # write a series (phenotype) to a CLS file

    if phenotype_type == 'categorical' or phenotype_type == 'binary':
        
        size = len(series)
        classes = list(series.unique())
        classes_string = ' '.join(classes)
        phen_string = ' '.join(series)
        num_classes = len(classes)

        with open(filename, mode = 'w') as output_file:

            output_file.writelines('{} {} 1\n'.format(size, num_classes))
            output_file.writelines('#{}\n'.format(classes_string))
            output_file.writelines('{}'.format(phen_string))
          
    elif phenotype_type == 'continuous':

        phen_string = ' '.join(series.astype(str))
        
        with open(filename, mode = 'w') as output_file:

            output_file.writelines('#numeric\n')
            output_file.writelines('#{}\n'.format(series.name))
            output_file.writelines('{}'.format(phen_string))
          
    else:
        raise('ERROR: Unknown phenotype_type: {}'.format(phenotype_type))

    return