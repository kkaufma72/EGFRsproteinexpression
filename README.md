# %% [markdown]
# <table style="border:2px solid white;" cellspacing="0" cellpadding="0" border-collapse: collapse; border-spacing: 0;>
#   <tr> 
#     <th style="background-color:white"> <img src="media/CCAL.png" width=225 height=225></th>
#     <th style="background-color:white"> <img src="media/logoMoores.jpg" width=175 height=175></th>
#     <th style="background-color:white"> <img src="media/UCSD_School_of_Medicine_logo.png" width=175 height=175></th> 
#   </tr>
# </table>

# %% [markdown]
# # Analysis Example: Estimating the degree of association of genes with EGFR’s protein expression using the MCC entropy. 
# 
# ### In this example we use the protein expression profile of [EGFR](https://en.wikipedia.org/wiki/Epidermal_growth_factor_receptor) (RPPA) as a reference phenotype and we match to it the gene expression (RNA-Seq) profiles of 16,184 genes in 893 cancer cell lines using the MCC entropy. The relevant datasets were downloaded from the publicly available Broad Institute’s DepMap [depmap.org](http://depmap.org). The computation of the MCC entropies takes only about 2.8 seconds. 
# 
# ### We also randomly permuted the EGFR RPPA profile 1,000 times and matched all the genes profiles against these randomized profiles as part of a permutation test [69] to estimate the statistical significance of the scores, namely, the p-values and False Discovery Rates (FDRs). 

# %%


# %%


# %% [markdown]
# 

# %%


# %% [markdown]
# ## <span style="color: blue;"> Set up notebook environment and import pcal library</span>
# 
# ### <span style="color: blue;"> The pcal library contains all of our functions and computational methods</span>
# 
# 

# %%
%load_ext autoreload
%autoreload 2

import pcal
from notebook_set_up import *

print(datetime.datetime.now())
start = time.process_time() 

# %% [markdown]
# ## <span style="color: blue;"> Create results directory to store the analysis results that will be generated</span>

# %%
results_dir = 'results'

if not os.path.exists('{}'.format(results_dir)):
    os.mkdir('{}'.format(results_dir))
    print('Results directory has been created: {}'.format(results_dir))
else:    
    print('Results directory already exists: {}'.format(results_dir))

# %% [markdown]
# ## <span style="color: blue;"> Read DepMap Reverse Phase Protein Array (RPPA) dataset </span>
# 
# ### <span style="color: blue;"> Reverse Phase Protein Array (RRPA) are an antibody-based high throughput proteomic approach that enables the simultaneous quantification of large numbers of proteins. The dataset has proteins as rows and samples as columns and each entry shows the protein expression levels of a given protein in a given sample  ([Barretina et al](https://www.nature.com/articles/nature11003))
#  </span>

# %%
RPPA = pcal.read_dataset('data/RPPA_2024.gct')
RPPA

# %% [markdown]
# ## <span style="color: blue;"> Read DepMap Gene Expression (mRNA) RNA-Seq dataset and perform basic preprocessing</span>
# 
# ### <span style="color: blue;"> A processed RNA-Seq dataset is a collection of data generated from RNA sequencing technology that captures a snapshot of which genes are being expressed (transcribed into RNA) in a biological sample at a given time. The dataset has genes as rows and samples as columns and each entry shows the gene expression levels of a given gene in a given sample  ([Barretina et al](https://www.nature.com/articles/nature11003))
#  </span>

# %%
mRNA = pcal.read_dataset('data/mRNA_2024.gct')
print('Initial dataset size before preprocessing {}'.format(mRNA.shape))

# replace NAs with zeroes

mRNA.fillna(value=0, inplace=True) 

# merge duplicate gene entries into one with the max value

mRNA = mRNA.groupby(mRNA.index, axis=0).max()

# exclude genes that do not have at leat 50% of entries different from zero

non_zero_genes = []    
pct_non_zero_thres = 0.50
for gene in mRNA.index:
    if np.sum(mRNA.loc[gene,:] > 0)/mRNA.shape[1] >= pct_non_zero_thres:
        non_zero_genes.append(gene)
        
mRNA = mRNA.loc[non_zero_genes, :]          

# exclude genes with standard deviation less than 0.01

stds = mRNA.std(axis = 1)
mRNA = mRNA.loc[stds >= 0.01, :]

print('Final dataset size after preprocessing {}'.format(mRNA.shape))

mRNA

# %% [markdown]
# ## <span style="color: blue;"> Match gene expression profiles against the protein expression profile of EGFR using the MCC entropy as a measure of association </span>
# 
# ### <span style="color: blue;"> The generated heatmap shows  the top/bottom 25 genes most positively/negatively associated with the EGFR protein profile where the reference and genes profiles have been standardized and clipped at 2.25 standard deviations. The color mapping shows red for feature values corresponding to values above the mean, and blue for values below the mean, while white represents the mean. The top hit, namely tho most associated gene (mRNA) with the EGFR protein profile is EGFR itself (MCCe = 0.90 ± 0.015). </span>
# 

# %%
EGFR_RPPA_vs_mRNA_scores =  pcal.match_target_vs_features(
                                    target = RPPA.loc['EGFR', :],
                                    target_type = 'continuous',
                                    target_ascending = False,
                                    truncate_feature_names_n_chars = 1000,
                                    features = mRNA,
                                    features_type = 'continuous',
                                    normalize_features = True,
                                    metric = pcal.MCC_entropy_multi,
                                    normalize_features_clims = [-2.25, 2.25],
                                    features_scores_ascending = False,
                                    CI_n_bootstraps = 100,
                                    p_val_n_permutations = 1,
                                    title = '',
                                    title_font_size = 12,
                                    n_features_plot = 25,
                                    font_scale = 0.65,
                                    row_compression = 1.0,
                                    figsize = 'auto',
                                    plot_dpi = 500,
                                    fill_na_with_zeroes = True,
                                    print_timings = True,
                                    cluster_within_category = False,
                                    filepath_prefix = '{}/EGFR_RPPA_vs_mRNA'.format(results_dir)) 

EGFR_RPPA_vs_mRNA_scores

# %%


# %% [markdown]
# ## <span style="color: blue;"> Plot the MCC entropies for all genes (sorted) </span>
# 

# %%
fig, axs = plt.subplots(1, 2, figsize=(16, 5), width_ratios = [1, 0])

axs[0].scatter(range(len(EGFR_RPPA_vs_mRNA_scores)), EGFR_RPPA_vs_mRNA_scores.loc[:, 'MCC Entropy'], marker = '.', s = 50, color = 'blue') #, linewidth = 3)
axs[0].set_xlabel('Gene Rank', fontsize = 16)
axs[0].set_ylabel('$MCC$ Entropy', fontsize = 16)


axs[0].set_xticks(ticks = [x for x in np.arange(0, 16184, 1000)],
                   labels = [x for x in np.arange(0, 16184, 1000)], fontsize=14) 

axs[0].set_yticks(ticks = [y for y in np.arange(-0.6, 1.1, 0.1)], 
                   labels = [np.round(y, 1) for y in np.arange(-0.6, 1.1, 0.1)], fontsize=16) 


axs[1].set_visible(False)

axs[0].plot((0, 16184), (0, 0), color = 'gray', linestyle = '--')

fig.savefig(fname = '{}/MCC_for_all_genes.png'.format(results_dir))


# %% [markdown]
# ## <span style="color: blue;"> Generate a Cumulative Entropy Diagram of the EGFR reference profile and the features (genes) </span>
# 
# ### <span style="color: blue;"> This diagram shows the angular distance proportional to the MCC entropy and the radial distance proportional to the cumulative entropy of the marginal variables. The Cumulative Entropy Diagram (CED) helps to visualize the relationships between a reference variable and a group of additional variables, and their intrinsic cumulative entropies. The reference profile (EGFR protein) is shinw in red and the genes (mRNA) in orange. A subset of EGFR-relevant genes is shown in blue. </span>
# 
# 
# 

# %%
EGFR_relevant_genes = ['EGFR', 'ANXA1', 'ANXA2', 'BCAR3', 'CAV1', 'CLDN1', 'CYR61', 'ITGA2', 
                       'EPHA2', 'LAMC2', 'MET', 'MYOF', 'PLAU', 'TGFBI', 'YAP1', 'AMOTL2', 'TINAGL1']


pcal.cumulative_entropy_diagram(reference =  RPPA.loc['EGFR', :],
                           features = mRNA,
                            features_to_highlight = EGFR_relevant_genes,
                            features_to_highlight_color = 'blue',
                            reference_name = 'EGFR (protein)', 
                            figsize = (8, 8),
                            assoc_metric = pcal.MCC_entropy,
                            assoc_metric_name = 'MCC Entropy',
                            marginal_metric = pcal.cumulative_entropy,
                            marginal_metric_name = 'Cumulative Entropy',
                            marker_size = 500,
                            proportional_marker_size = True,
                            proportional_marker_size_exp = 3,
                            reference_marker_size = 400,
                            reference_marker_color = 'red', 
                            reference_text_size = 6,
                            reference_text_color = 'red',
                            feature_text_size = 6,
                            feature_text_color = 'blue',
                            fea_label_dist_factor = 1.05,
                            title = '',
                            features_to_label = EGFR_relevant_genes,
                            marker_color = 'orange', 
                            max_angle = np.pi,
                            min_angle = 0.0,
                            max_assoc = 1.0,
                            min_assoc = 0,
                            max_marg = 0.41,
                            show_grid = True,
                            grid_color = 'grey',
                            facecolor = '#EBF0F7',
                            output_file = '{}/EGFR_cumulative_entropy_diagram.png'.format(results_dir),
                            plot_dpi = 400)

# %%


# %%


# %%


# %%


# %%


# %%


# %%


