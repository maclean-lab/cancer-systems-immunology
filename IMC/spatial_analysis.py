# Filename: spatial_analysis.py
# Description: perform co-occurrence analysis between select immune clusters
# Author: Xiaojun Wu
# Email: xiaojunw@usc.edu

# %%
# load libraries
import os
from typing import Literal

import numpy as np
import scipy
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import ListedColormap
import seaborn as sns

import scanpy as sc
import squidpy as sq

# %%
# load data
adata_full = sc.read_h5ad('data/imc_cells_subclustered.h5ad')
adata_full.obs['cell_id'] = adata_full.obs['cell_id'].astype('string')
adata_full.obs['ObjectNumber'] = adata_full.obs['ObjectNumber'].astype(int)
adata_full.obs['batch'] = adata_full.obs['batch'].astype(int)
adata_full.obs['patient_id'] = adata_full.obs['patient_id'].astype(int)
adata_full.obs['pub_id'] = adata_full.obs['pub_id'].astype(int)
adata_full.obs['roi'] = adata_full.obs['roi'].astype(int)
metacluster_labels = adata_full.obs['cluster'].cat.categories
subcluster_labels = adata_full.obs['subcluster'].cat.categories

# load subject information
subject_info = pd.read_csv('data/subject_info.csv', index_col=0)
subject_info = subject_info.drop_duplicates(subset=['pub_id', 'responder'])
subject_info['pub_id'] = subject_info['pub_id'].astype(int)
subject_info['patient_id'] = subject_info['patient_id'].astype('string')
subject_info['tissue'] = subject_info['tissue'].astype('string')
subject_info['HR_status'] = subject_info['HR_status'].astype('category')
subject_info['responder'] = subject_info['responder'].astype('boolean')
subject_info = subject_info.sort_values('pub_id')
subject_info.index = subject_info['pub_id']
subject_info.index.name = None

# update responder information
adata_full.obs = adata_full.obs.merge(subject_info[['pub_id', 'responder']],
                                      on='pub_id', how='left')
adata_full.obs.drop(columns=['responder_x'], inplace=True)
adata_full.obs.rename(columns={'responder_y': 'responder'}, inplace=True)

# move cell_id to obs_names
adata_full.obs_names = adata_full.obs['cell_id']
adata_full.obs.drop(columns=['cell_id'], inplace=True)

# update protein names
adata_full.var_names = adata_full.var_names.astype('string')

# add spatial coordinates to obsm
adata_full.obsm['spatial'] = adata_full.obs.loc[:, ['Pos_X', 'Pos_Y']].values

# split data by ROI
adata_all = {}

for pub_id in adata_full.obs['pub_id'].unique():
    adata_patient = adata_full[adata_full.obs['pub_id'] == pub_id].copy()

    for group in adata_patient.obs['group'].unique():
        adata_group = adata_patient[adata_patient.obs['group'] == group].copy()

        for roi in adata_group.obs['roi'].unique():
            adata = adata_group[adata_group.obs['roi'] == roi].copy()
            adata.obs['cluster'] = \
                adata.obs['cluster'].cat.set_categories(metacluster_labels)
            adata.obs['subcluster'] = \
                adata.obs['subcluster'].cat.set_categories(subcluster_labels)
            adata_all[(pub_id, group, roi)] = adata

adata_all = dict(sorted(adata_all.items()))

output_dir = os.path.join('.', 'co-occurrence')
if not os.path.exists(output_dir):
    os.makedirs(output_dir)


# %%
# plot spatial distribution of clusters
# plot for metaclusters
# can we get colors from the returned axes?
figure_path = os.path.join(output_dir, 'metaclusters_spatial.pdf')
with PdfPages(figure_path) as pdf:
    for (pub_id, group, roi), adata in adata_all.items():
        sq.pl.spatial_scatter(adata, shape=None, color='cluster',
                              title=f'R-{pub_id:02d}, {group}, {roi:03d}',
                              alpha=0.5, edgecolor='none')
        pdf.savefig()
        plt.close()

# plot for subclusters
figure_path = os.path.join(output_dir, 'subclusters_spatial.pdf')
with PdfPages(figure_path) as pdf:
    for (pub_id, group, roi), adata in adata_all.items():
        sq.pl.spatial_scatter(adata, shape=None, color='subcluster',
                              title=f'R-{pub_id:02d}, {group}, {roi:03d}',
                              alpha=0.5, edgecolor='none')
        pdf.savefig()
        plt.close()

# %%
# compute co-occurrence scores
co_occurrence_dists = [0, 10, 25] + list(range(50, 301, 50))
co_occurrence_columns = ['pub_id', 'responder', 'group', 'roi', 'dist',
                         'source', 'target', 'score']
sample_co_occurrence_scores = []

for (pub_id, group, roi), adata in adata_all.items():
    adata = adata.copy()

    sample_scores, _ = sq.gr.co_occurrence(
        adata, cluster_key='cluster', interval=co_occurrence_dists,
        copy=True, show_progress_bar=False)

    cluster_labels = adata.obs['cluster'].cat.categories

    # convert co-occurrence matrices to a DataFrame
    # get indices of computed values and map to cluster labels or distances
    scores_mask = ~np.isnan(sample_scores)
    scores_long = pd.DataFrame(np.column_stack(np.where(scores_mask)))
    scores_long.columns = ['source', 'target', 'dist']
    scores_long['source'] = scores_long['source'].map(
        lambda x: cluster_labels[x])
    scores_long['target'] = scores_long['target'].map(
        lambda x: cluster_labels[x])
    scores_long['dist'] = scores_long['dist'].map(
        lambda x: co_occurrence_dists[x + 1])

    # add the corresponding scores
    scores_long['score'] = sample_scores[scores_mask]

    # add pub_id, group, roi, response status
    scores_long['pub_id'] = pub_id
    scores_long['group'] = group
    scores_long['roi'] = roi
    scores_long['responder'] = subject_info.loc[pub_id, 'responder']

    # append co-occurrence scores of the current sample
    scores_long = scores_long[co_occurrence_columns]
    sample_co_occurrence_scores.append(scores_long)

# combine co-occurrence scores of all samples
co_occurrence_scores = pd.concat(sample_co_occurrence_scores,
                                 ignore_index=True)
co_occurrence_scores['pub_id'] = co_occurrence_scores['pub_id'].astype(
    'category')
co_occurrence_scores['responder'] = co_occurrence_scores['responder'].astype(
    'boolean')
co_occurrence_scores['group'] = pd.Categorical(
    co_occurrence_scores['group'], categories=['BL', 'C1D1', 'WK8'],
    ordered=True)
co_occurrence_scores['dist'] = co_occurrence_scores['dist'].astype('category')
co_occurrence_scores['source'] = co_occurrence_scores['source'].astype(
    'category')
co_occurrence_scores['target'] = co_occurrence_scores['target'].astype(
    'category')

# cluster pairs for downstream analyses
cluster_pairs = [('CD8+ T cell', 'macrophage'),
                 ('CD8+ T cell', 'MDSC'),
                 ('B cell', 'CD8+ T cell'),
                 ('B cell', 'B cell'),
                 ('B cell', 'dendritic cell')]


# %%
# define function to plot co-occurrence scores
def plot_co_occurrence_scores(dist_scores: pd.DataFrame, dist: int | float,
                              by: Literal['subject', 'time_point']):
    figure_path = os.path.join(
        output_dir, f'co_occurrence_by_{by}_dist_{dist}.pdf')

    with PdfPages(figure_path) as pdf:
        for source, target in cluster_pairs:
            # get scores for the current pair of clusters
            pair_scores = dist_scores[(dist_scores['source'] == source)
                                      & (dist_scores['target'] == target)]

            fig, axes = plt.subplots(
                1, 2, figsize=(10, 4), sharey=True,
                width_ratios=[responder_ids.size, non_responder_ids.size])

            # plot responder scores
            responder_scores = pair_scores[pair_scores['responder']].copy()
            responder_scores['pub_id'] = responder_scores[
                'pub_id'].cat.set_categories(responder_ids, ordered=True)
            axes[0].set_title('Responders')
            if by == 'subject':
                sns.stripplot(data=responder_scores, x='pub_id',
                              y='score', hue='group', jitter=False, dodge=True,
                              size=marker_size, log_scale=True,
                              palette=custom_palette, ax=axes[0])
            else:
                sns.stripplot(data=responder_scores, x='group',
                              y='score', hue='group', jitter=True,
                              size=marker_size, log_scale=True,
                              palette=custom_palette, ax=axes[0])

            # plot non-responder scores
            non_responder_scores = pair_scores[
                ~pair_scores['responder']].copy()
            non_responder_scores['pub_id'] = non_responder_scores[
                'pub_id'].cat.set_categories(non_responder_ids, ordered=True)
            if by == 'subject':
                sns.stripplot(data=non_responder_scores, x='pub_id',
                              y='score', hue='group', jitter=False, dodge=True,
                              size=marker_size, log_scale=True,
                              palette=custom_palette, ax=axes[1])
            else:
                sns.stripplot(data=non_responder_scores, x='group',
                              y='score', hue='group', jitter=True,
                              size=marker_size, log_scale=True,
                              palette=custom_palette, ax=axes[1])
            axes[1].set_title('Non-responders')

            # add horizontal line at y=1 and remove legend inside each axes
            for ax in axes:
                ax.axhline(y=1, color='#777777', linestyle='--')

                ax_legnd = ax.get_legend()
                if ax_legnd is not None:
                    ax.get_legend().remove()

            # add figure title
            fig.suptitle(f'{source} -- {target}, dist: {dist}')

            # add legend to the right of the figure
            if by == 'subject':
                handles, labels = axes[0].get_legend_handles_labels()
                fig.legend(handles, labels, title='Time point',
                           loc='center right', bbox_to_anchor=(1, 0.5),
                           frameon=True)

            # save current page
            plt.tight_layout(rect=[0, 0, 0.9, 1])
            pdf.savefig(fig)
            plt.close(fig)


# %%
# plot co-occurrence scores at each distance
# one PDF per distance, all cluster pairs in one PDF
responder_ids = subject_info[subject_info['responder']].index
non_responder_ids = subject_info[~subject_info['responder']].index
custom_palette = ['#a6cee3', '#1f78b4', '#08306b']
marker_size = 8

for dist in co_occurrence_dists[1:]:
    # get scores for the current distance and filter out subjects without
    # known response status
    dist_scores = co_occurrence_scores[
        (co_occurrence_scores['dist'] == dist)
        & (~co_occurrence_scores['responder'].isna())].copy()

    # plot scores of each subject across all time points
    plot_co_occurrence_scores(dist_scores, dist, 'subject')
    # plot scores of each time point across all subjects
    plot_co_occurrence_scores(dist_scores, dist, 'time_point')


# %%
# plot co-occurrence scores vs distances (R-02 only)
subject_scores = co_occurrence_scores[co_occurrence_scores['pub_id'] == 2]
figure_path = os.path.join(output_dir, 'co_occurrence_vs_distance.pdf')

with PdfPages(figure_path) as pdf:
    for source, target in cluster_pairs:
        pair_scores = subject_scores[
            (subject_scores['source'] == source)
            & (subject_scores['target'] == target)]

        num_group_samples = pair_scores.groupby('group')['roi'].nunique()
        num_group_samples = num_group_samples[num_group_samples > 0]
        fig, axes = plt.subplots(1, num_group_samples.size,
                                 figsize=(num_group_samples.size * 4, 4),
                                 sharey=True)
        for i, group in enumerate(num_group_samples.index):
            group_scores = pair_scores[pair_scores['group'] == group].copy()
            sns.lineplot(data=group_scores, x='dist', y='score', hue='roi',
                         palette=sns.color_palette('rainbow'), ax=axes[i])
            axes[i].set_title(group)

            # add horizontal line at y=1
            axes[i].axhline(y=1, color='#777777', linestyle='--')

        # add figure title
        fig.suptitle(f'{source} -- {target}, R-02')

        # save current page
        plt.tight_layout(rect=[0, 0, 1, 1])
        pdf.savefig(fig)
        plt.close(fig)

# %%
# run t-test for co-occurrence scores between time points for individual
# subjects
co_occurrence_t_tests = pd.DataFrame(
    columns=['pub_id', 'responder', 'source', 'target', 'groups',
             'num_samples_group_1', 'num_samples_group_2', 'hypothesis',
             'statistic', 'p_value'])
group_pairs = [('BL', 'C1D1'), ('C1D1', 'WK8')]
test_dist = 50
test_hypotheses = ['two-sided', 'less', 'greater']

for source, target in cluster_pairs:
    for group_1, group_2 in group_pairs:
        for pub_id in co_occurrence_scores['pub_id'].cat.categories:
            subject_scores = co_occurrence_scores[
                (co_occurrence_scores['pub_id'] == pub_id)
                & (co_occurrence_scores['dist'] == test_dist)
                & (co_occurrence_scores['source'] == source)
                & (co_occurrence_scores['target'] == target)]
            group_1_scores = subject_scores.loc[
                subject_scores['group'] == group_1, 'score'].to_numpy()
            group_1_scores = group_1_scores[~np.isnan(group_1_scores)]
            group_2_scores = subject_scores.loc[
                subject_scores['group'] == group_2, 'score'].to_numpy()
            group_2_scores = group_2_scores[~np.isnan(group_2_scores)]

            if group_1_scores.size < 2 or group_2_scores.size < 2:
                continue

            for hypothesis in test_hypotheses:
                test_result = scipy.stats.ttest_ind(
                    group_1_scores, group_2_scores, alternative=hypothesis)
                new_row = pd.DataFrame([{
                    'pub_id': pub_id,
                    'responder': subject_info.loc[pub_id, 'responder'],
                    'source': source,
                    'target': target,
                    'groups': f'{group_1} vs {group_2}',
                    'num_samples_group_1': group_1_scores.size,
                    'num_samples_group_2': group_2_scores.size,
                    'hypothesis': hypothesis,
                    'statistic': test_result.statistic,
                    'p_value': test_result.pvalue
                }])
                co_occurrence_t_tests = pd.concat(
                    [co_occurrence_t_tests, new_row], ignore_index=True)

output_path = os.path.join(output_dir, 'co_occurrence_t_tests.csv')
co_occurrence_t_tests.to_csv(output_path, index=False)


# %%
# run ANOVA and Kruskal-Wallis tests for co-occurrence scores from all time
# points, for either all responders or all non-responders
co_occurrence_tests = pd.DataFrame(
    columns=[
        'responder', 'source', 'target', 'bl_size', 'c1d1_size', 'wk8_size',
        'bl_mean', 'c1d1_mean', 'wk8_mean', 'f_stat', 'f_p_value',
        'bl_median', 'c1d1_median', 'wk8_median', 'kw_stat', 'kw_p_value'])
test_dist = 10

for is_responder in [True, False]:
    for source, target in cluster_pairs:
        response_scores = co_occurrence_scores[
            (co_occurrence_scores['responder'] == is_responder)
            & (co_occurrence_scores['dist'] == test_dist)
            & (co_occurrence_scores['source'] == source)
            & (co_occurrence_scores['target'] == target)]
        bl_scores = response_scores.loc[
            response_scores['group'] == 'BL', 'score'].to_numpy()
        bl_scores = bl_scores[~np.isnan(bl_scores)]
        c1d1_scores = response_scores.loc[
            response_scores['group'] == 'C1D1', 'score'].to_numpy()
        c1d1_scores = c1d1_scores[~np.isnan(c1d1_scores)]
        wk8_scores = response_scores.loc[
            response_scores['group'] == 'WK8', 'score'].to_numpy()
        wk8_scores = wk8_scores[~np.isnan(wk8_scores)]

        # scipy documentation recommends at least 5 samples per group for the
        # Kruskal-Wallis test
        if min(bl_scores.size, c1d1_scores.size, wk8_scores.size) < 5:
            continue

        anova_result = scipy.stats.f_oneway(bl_scores, c1d1_scores, wk8_scores)
        kw_result = scipy.stats.kruskal(bl_scores, c1d1_scores, wk8_scores)
        new_row = pd.DataFrame([{
            'responder': is_responder,
            'source': source,
            'target': target,
            'bl_size': bl_scores.size,
            'c1d1_size': c1d1_scores.size,
            'wk8_size': wk8_scores.size,
            'bl_mean': np.mean(bl_scores),
            'c1d1_mean': np.mean(c1d1_scores),
            'wk8_mean': np.mean(wk8_scores),
            'f_stat': anova_result.statistic,
            'f_p_value': anova_result.pvalue,
            'bl_median': np.median(bl_scores),
            'c1d1_median': np.median(c1d1_scores),
            'wk8_median': np.median(wk8_scores),
            'kw_stat': kw_result.statistic,
            'kw_p_value': kw_result.pvalue
        }])
        co_occurrence_tests = pd.concat(
            [co_occurrence_tests, new_row], ignore_index=True)

output_path = os.path.join(output_dir, 'co_occurrence_anova_kw.csv')
co_occurrence_tests.to_csv(output_path, index=False)

# %%
# plot CD8+ T cells and macrophages on spatial coordinates
# the following cells are shown in the same color:
# - CD8+ T cells
# - macrophages within 10 microns of any CD8+ T cell
# - all other macrophages
figure_path = os.path.join(output_dir, 'R-02_cd8t_vs_macrophage.pdf')
neighbor_dist = 10
cd8t_neighbor_labels = {True: 'macrophage (<=10)', False: 'macrophage (>10)'}
custom_palette = ListedColormap(['C2', 'C6', 'blueviolet'])

with PdfPages(figure_path) as pdf:
    for (pub_id, group, roi), adata in adata_all.items():
        # skip if not R-02
        if pub_id != 2:
            continue

        # keep only CD8+ T cells and macrophages for plotting
        adata = adata[
            adata.obs['cluster'].isin(['CD8+ T cell', 'macrophage']), :].copy()

        # determine which macrophages close neighbors of CD8+ T cells
        sq.gr.spatial_neighbors(adata, coord_type='generic',
                                n_neighs=adata.shape[0], radius=neighbor_dist)
        neighbors = adata.obsp['spatial_connectivities']

        # assign labels to all cells for plotting
        # - CD8+ T cell
        # - macrophage (<=10)
        # - macrophage (>10)
        adata.obs['plot_label'] = adata.obs['cluster'].astype('string')
        cd8t_indices = np.where(
            adata.obs['cluster'].values == 'CD8+ T cell')[0]

        for i, cell_id in enumerate(adata.obs_names):
            # no need to re-assign label for a CD8+ T cell
            if adata.obs.loc[cell_id, 'cluster'] == 'CD8+ T cell':
                continue

            # get neighbors of the current macrophage
            cell_neighbors = np.squeeze(neighbors[i, :].toarray())

            # assign label based on whether any CD8+ T cells are neighbors
            is_cd8t_neighbor = np.any(cell_neighbors[cd8t_indices])
            adata.obs.loc[cell_id, 'plot_label'] = \
                cd8t_neighbor_labels[is_cd8t_neighbor]

        # make the plot and save
        sq.pl.spatial_scatter(adata, shape=None, color='plot_label',
                              title=f'R-{pub_id:02d}, {group}, {roi:03d}',
                              size=20, edgecolor='none',
                              palette=custom_palette)
        pdf.savefig()
        plt.close()
