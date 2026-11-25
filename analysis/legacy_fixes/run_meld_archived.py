# Archived MELD exploratory source, distinct from the damaged historical R plot.
# Set MELD_PARAMETER_SEARCH_FILE and MELD_OUTPUT_FILE before running.
# This source requires the original expression/metadata and parameter-search file.
import os

# %%
import pandas as pd
import numpy as np
import graphtools as gt
import phate
import magic
import scprep
import meld
import cmocean
import sklearn
import scipy
import seaborn as sns

# setting defaults for matplotlib font sizes
import matplotlib.pyplot as plt
plt.rc('font', size=14)

# making sure plots & clusters are reproducible
np.random.seed(42)

# Archived IPython-only line: %load_ext autoreload
# Archived IPython-only line: %autoreload 2

import diffxpy.api as de

# %%
data = pd.read_csv("./1_exp_matrix.csv")

# %%
data
#行 基因 列 细胞

# %%
data = data.T

# %%
data.columns=data.iloc[0]

# %%
data= data.drop(['Unnamed: 0'])

# %%
data
#行 细胞  列 基因

# %%
metadata = pd.read_csv("./1_metadata.csv",index_col=0)

# %%
metadata
#行 细胞

# %%
#data = scprep.filter.filter_rare_genes(data)

# %%
#data, metadata = scprep.filter.filter_library_size(
#    data, metadata, cutoff=15000, 
#    keep_cells='below')

# %%
metadata1=metadata[(metadata['orig.ident']=='Lin-con') | (metadata['orig.ident']=='IR')]

# %%
metadata1

# %%
a=metadata1.index.tolist()

# %%
data=data.loc[a,]

# %%
data

# %%
metadata=metadata1

# %%
data_libnorm, libsize = scprep.normalize.library_size_normalize(data, return_library_size=True)

# %%
#type(metadata1.index)

# %%
metadata['library_size'] = libsize
metadata.head()

# %%
data_sqrt = np.sqrt(data_libnorm)

# %%
fig, ax = plt.subplots(1)

groups, counts = np.unique(metadata['orig.ident'], return_counts=True)
for i, c in enumerate(counts):
    ax.bar(i, c)
    
ax.set_xticks(np.arange(i+1))
ax.set_xticklabels(groups)
ax.set_ylabel('# cells')

fig.tight_layout()

# %%
data_pca = scprep.reduce.pca(data_sqrt)

phate_op = phate.PHATE(knn=10, decay=10, n_jobs=-1)
data_phate = phate_op.fit_transform(data_pca)

# %%
scprep.plot.scatter2d(data_phate, c=metadata['orig.ident'], 
                      legend_anchor=(1,1), figsize=(7,5), s=10, label_prefix='PHATE', ticks=False)

# %%
scprep.plot.scatter2d(data_phate, c=metadata['celltype'], 
                      legend_anchor=(1,1), figsize=(12,5), s=10, label_prefix='PHATE', ticks=False)

# %%
metadata['IR'] = [1 if sl.startswith('IR') else 0 for sl in metadata['orig.ident']]

# %%
metadata['sample']=metadata['orig.ident'] 

# %%
metadata['sample2']=metadata['sample']+'A'

# %%
import tempfile
with tempfile.TemporaryDirectory() as tempdir:
        #filepath = os.path.join(tempdir, "parameter_search_results.pickle.gz")
        #scprep.io.download.download_url(URL, filepath)
        #results = pd.read_pickle(filepath)
        results = pd.read_pickle(os.environ["MELD_PARAMETER_SEARCH_FILE"])

# We want to take the average of each set of random seeds for each combination of beta and knn values
results_wide = results.groupby(['beta', 'knn']).mean().sort_values(by='MSE').reset_index()

# %%
results.groupby(['beta', 'knn']).mean()

# %%
results.groupby(['beta', 'knn']).mean().sort_values(by='MSE')

# %%
results

# %%
results_wide

# %%
ax = scprep.plot.scatter(results_wide['beta'], results_wide['knn'], 
                         s=50, c=results_wide['MSE'], vmax=0.006, cmap='inferno_r')

# Highlight the top performing combination with a large red dot
top_result = results_wide.sort_values('MSE').iloc[0]
ax.scatter(top_result['beta'], top_result['knn'], c='r', s=100, linewidth=1, edgecolor='k')

# %%
top_result

# %%
G = gt.Graph(data_pca, knn=int(top_result['knn']), use_pygsp=True)

# %%
meld_op = meld.MELD(beta=top_result['beta'])
sample_densities = meld_op.fit_transform(G, sample_labels=metadata['sample2'])

# %%
def replicate_normalize_densities(sample_densities, replicate):
    replicates = np.unique(replicate)
    sample_likelihoods = sample_densities.copy()
    for rep in replicates:
        curr_cols = sample_densities.columns[[col.endswith(rep) for col in sample_densities.columns]]
        sample_likelihoods[curr_cols] = sklearn.preprocessing.normalize(sample_densities[curr_cols], norm='l1')
    return sample_likelihoods

# %%
sample_likelihoods = sample_densities.copy()

# %%
sample_likelihoods

# %%
curr_cols = sample_densities.columns[[col.endswith('A') for col in sample_densities.columns]]
sample_likelihoods[curr_cols] = sklearn.preprocessing.normalize(sample_densities[curr_cols], norm='l1')

# %%
sample_likelihoods

# %%
#指定实验因素
experimental_samples = ['IRA']
curr_sample = experimental_samples[0]
scprep.plot.scatter2d(data_phate, c=sample_likelihoods[curr_sample],cmap=meld.get_meld_cmap(),
                          vmin=0, vmax=1,
                          title=curr_sample, ticks=False)

# %%
fig, axes = plt.subplots(1,2, figsize=(8.7,4))

scprep.plot.scatter2d(data_phate, c=sample_likelihoods[experimental_samples].mean(axis=1), 
                      cmap=meld.get_meld_cmap(), vmin=0, vmax=1,
                      title='Mean', ticks=False, ax=axes[0])
scprep.plot.scatter2d(data_phate, c=sample_likelihoods[experimental_samples].std(axis=1), vmin=0, 
                      cmap='inferno', title='St. Dev.', ticks=False, ax=axes[1])

fig.tight_layout()

# %%
experimental_samples

# %%
metadata['chd_likelihood'] = sample_likelihoods[experimental_samples].mean(axis=1).values

# %%
sample_cmap = {'chdA': '#fb6a4a',
                 'chdB': '#de2d26',
                 'chdC': '#a50f15',
                 'tyrA': '#6baed6',
                 'tyrB': '#3182bd',
                 'tyrC': '#08519c'}

# %%
metadata['clusterID'] = scprep.utils.sort_clusters_by_values(metadata['celltype'], metadata['chd_likelihood'])

# %%
#横坐标 分组 纵坐标 chd 颜色 样本

fig, ax = plt.subplots(1, figsize=(10,10))

# See example usage: https://scprep.readthedocs.io/en/stable/examples/jitter.html
scprep.plot.jitter(metadata['clusterID'], metadata['chd_likelihood'], c=metadata['IR'], 
                   legend=False, plot_means=False, xlabel=False, ylabel='Mean chd likelihood',
                   ax=ax)

### This code will plot the ratio of tyr:chd cells per cluster
means = metadata.groupby('clusterID')['IR'].mean()
ax.scatter(means.index, means - np.mean(metadata['IR']) + 0.5, edgecolor='k', s=100)

# Axis tick labels
ax.set_xticklabels(metadata.set_index('clusterID')['celltype'].drop_duplicates().sort_index(), rotation=90)
ax.set_ylim(0,1)

fig.tight_layout()

# %%
metadata

# %%
metadata.to_csv(os.environ["MELD_OUTPUT_FILE"])

# %%
# Get cluster indicides and number of cells per cluster
clusters, counts = np.unique(metadata['clusterID'], return_counts=True)

# Keep cluster labels with at least 1% of the data
clusters = clusters[counts > data.shape[0] * 0.01]

# %%
data_cluster_phate = {}

for cluster in clusters:
    curr_data = data_pca.loc[metadata['clusterID'] == cluster]
    data_cluster_phate[cluster] = phate.PHATE(verbose=0).fit_transform(curr_data)

# %%
fig,axes= plt.subplots(5,4, figsize=(4*3, 5*3))

for i , ax in enumerate(axes.flatten()):
    if not i < len(clusters):
        ax.axis('off')
        continue
    curr_cluster = clusters[i]
    curr_phate = data_cluster_phate[curr_cluster]
    
    scprep.plot.scatter2d(curr_phate, 
                          c=metadata['chd_likelihood'].loc[metadata['clusterID'] == curr_cluster], 
                          cmap=meld.get_meld_cmap(), vmin=0, vmax=1,
                         ax=ax, ticks=False, 
                          title='Cluster {} ({})'.format(curr_cluster, curr_phate.shape[0]), 
                          legend=False, fontsize=10)

# %%
np.random.seed(0)
vfc_op_per_cluster = {}

for cluster in np.unique(metadata['clusterID']):
    curr_G = gt.Graph(data_pca.loc[metadata['clusterID'] == cluster], use_pygsp=True)
    curr_G.compute_fourier_basis()
    curr_sample_labels = metadata['IR'].loc[metadata['clusterID'] == cluster]
    curr_likelihood = metadata['chd_likelihood'].loc[metadata['clusterID'] == cluster]
    curr_vfc = meld.VertexFrequencyCluster(n_clusters = 3)
    curr_vfc.fit_transform(curr_G, curr_sample_labels, curr_likelihood)
    vfc_op_per_cluster[cluster] = curr_vfc

# %%
subclustering_results = {}
for cluster in np.unique(metadata['clusterID']):
    curr_vfc = vfc_op_per_cluster[cluster]
    clusters_by_n = {}
    for n in [2,3,4,5]:
        clusters_by_n[n] = curr_vfc.predict(n)
    subclustering_results[cluster] = clusters_by_n

# %%



