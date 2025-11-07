from scipy import stats
from scipy.stats import pearsonr 
from numpy.random import random_sample, seed, normal, choice, random_sample
import pandas as pd
from numpy import sqrt, round, exp

def simulated_annealing_projection(
    dist_matrix,
    x,
    y,
    arch_names,
    SA_max_iter = 20000, 
    SA_sigma = 0.40, 
    SA_init_temp = 0.0075,
    SA_random_seed = 2343):

    # Initialize
    seed(SA_random_seed)
    N = dist_matrix.shape[0]

    cor_vs_t = pd.Series(0, index = range(SA_max_iter))
    proj_dist_matrix = pd.DataFrame(0, index = dist_matrix.index, columns = dist_matrix.columns)
    
    for arch1 in arch_names:
        for arch2 in arch_names:
            proj_dist_matrix.loc[arch1, arch2] = sqrt((x[arch1] - x[arch2])**2 + (y[arch1] - y[arch2])**2)
        
    for arch1 in arch_names:
        for arch2 in arch_names:
            proj_dist_matrix.loc[arch1, arch2] = sqrt((x[arch1] - x[arch2])**2 + (y[arch1] - y[arch2])**2)
        
    # Compute initial correlation between projected and original distance matrices
    cor_before = pearsonr(dist_matrix.values.flatten(), proj_dist_matrix.values.flatten())[0]

    for t in range(SA_max_iter):
        SA_temp = SA_init_temp*(1 - t/(SA_max_iter + 1))    # lower temperature linearly at every step
        if t % 100 == 0:
            print('time: {}/{} correlation: {}'.format(t, SA_max_iter, round(cor_before, 5)))     
        
        # Generate proposed new a and y values for random node
        node = choice(arch_names)
        new_x = normal(x[node], SA_sigma, 1)
        new_y = normal(y[node], SA_sigma, 1)
    
        # Compute tentative new proj distance matrix with new values
        tent_proj_dist_matrix = proj_dist_matrix
        for arch2 in arch_names:
            tent_proj_dist_matrix.loc[node, arch2] = sqrt((new_x - x[arch2])**2 + (new_y - y[arch2])**2)
            tent_proj_dist_matrix.loc[arch2, node] = tent_proj_dist_matrix.loc[node, arch2]
            if arch2 == node:
                tent_proj_dist_matrix.loc[node, arch2] = 0.0
    
        cor_after = pearsonr(dist_matrix.values.flatten(), tent_proj_dist_matrix.values.flatten())[0]
    
        # Compute acceptance probability (notice increase of correlation produces prob > 1 gurantees acceptance)
        prob = exp((1/SA_temp)*(cor_after - cor_before))
    
        #print('before: {} after: {} prob: {}'.format(round(cor_before, 5), round(cor_after, 5), round(prob, 5)))  
    
        # Compare uniform random number vs. prob
        if random_sample((1,)) < prob:  # move accepted
            #print('move accepted')
            cor_vs_t.iloc[t] = cor_after
            cor_before = cor_after
            x[node] = new_x
            y[node] = new_y
            proj_dist_matrix = tent_proj_dist_matrix
    
        else: # move rejected
            #print('move rejected')
            cor_vs_t.iloc[t] = cor_before
    
    # Return time series, final theta configuration and projected distance matrix
    return {'cor_vs_t': cor_vs_t, 'x': x, 'y': y, 'proj_dist_matrix': proj_dist_matrix}
