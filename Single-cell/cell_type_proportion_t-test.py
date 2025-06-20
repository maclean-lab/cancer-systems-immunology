import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad
from scipy import stats
from statsmodels.stats.multitest import fdrcorrection

#read files
B = sc.read('/lung_mets_data/finalized/B_subclusters.h5ad')
MM = sc.read('/yingtong/lung_mets_data/finalized/mature_myeloid_subclusters.h5ad')
MDSC = sc.read('/lung_mets_data/finalized/MDSCs_subclusters.h5ad')
NK = sc.read('/lung_mets_data/finalized/NK_subclusters.h5ad')
T = sc.read('/lung_mets_data/T_cells_subclusters.h5ad')

CD45 = B.concatenate(MM, MDSC, NK, T_rd, batch_key = 'major_cluster', batch_categories =
                     ['B cells', 'Mature myeloid', 'MDSCs', 'NK cells', 'T cells'])

cell_number = pd.crosstab(CD45.obs['clusters'],CD45.obs['treatment'])
cell_number_sample = pd.crosstab(CD45.obs['clusters'], CD45.obs['sample_name'])
cell_type = list(cell_number.index)

v_samples = ['V1', 'V2a', 'V2b', 'V2c', 'V4a']
e_samples = ['E2', 'E3a', 'E3b', 'E4a']
ep_samples = ['EP1', 'EP2a', 'EP2b', 'EP3b']
ec_samples = ['EC1', 'EC3a', 'EC4b', 'EC4c']
epc_samples = ['EPC1', 'EPC3b', 'EPC4a']
pc_samples = ['PC1', 'PC2a', 'PC2b', 'PC4a']

def make_perct_list(adata, c, t, sample_set): # c=cell type, t= treatment
    perct_list = []
    for b in sample_set:
        c_n_obs= adata[(adata.obs['clusters'] == c) & (adata.obs["treatment"] == t)& (adata.obs['sample_name'] == b),:].n_obs
        cd45_n_obs = adata[(adata.obs['treatment'] == t) & (adata.obs['sample_name'] == b)].n_obs
        if cd45_n_obs == 0:
            continue
        else:
            perct_list.append(c_n_obs/cd45_n_obs)
    return(perct_list)

# V vs E
VE_test = pd.DataFrame(index = cell_type)
ve_tstats = []
ve_pval = []

for c in cell_type:
    #print(c)
    v = make_perct_list(CD45, c, "V", v_samples)
    e = make_perct_list(CD45, c, "E", e_samples)
    results = stats.ttest_ind(v, e)
    ve_tstats.append(results.statistic)
    ve_pval.append(results.pvalue)
    #print(v)
    #print(e)
    #print(results)

VE_test['V_cell_number'] = cell_number['V']
VE_test['E_cell_number'] = cell_number['E']
VE_test['t-stats'] = ve_tstats
VE_test['pval'] = ve_pval
ve_FDR = fdrcorrection(ve_pval)
VE_test['FDR'] = ve_FDR[1]

VE_test.to_csv('/lung_mets_data/finalized/VvsE_simple_T_test.csv')

# V vs EPC
VEPC_test = pd.DataFrame(index = cell_type)
vepc_tstats = []
vepc_pval = []

for c in cell_type:
    #print(c)
    v = make_perct_list(CD45, c, "V", v_samples)
    epc = make_perct_list(CD45, c, "EPC", epc_samples)
    results = stats.ttest_ind(v, epc)
    vepc_tstats.append(results.statistic)
    vepc_pval.append(results.pvalue)
    #print(v)
    #print(e)
    #print(results)

VEPC_test['V_cell_number'] = cell_number['V']
VEPC_test['EPC_cell_number'] = cell_number['EPC']
VEPC_test['t-stats'] = vepc_tstats
VEPC_test['pval'] = vepc_pval
vepc_FDR = fdrcorrection(vepc_pval)
VEPC_test['FDR'] = vepc_FDR[1]

VEPC_test.to_csv('/lung_mets_data/finalized/VvsEPC_simple_T_test.csv')

# E vs EPC
EEPC_test = pd.DataFrame(index = cell_type)
eepc_tstats = []
eepc_pval = []

for c in cell_type:
    #print(c)
    e = make_perct_list(CD45, c, "E", e_samples)
    epc = make_perct_list(CD45, c, "EPC", epc_samples)
    results = stats.ttest_ind(e, epc)
    eepc_tstats.append(results.statistic)
    eepc_pval.append(results.pvalue)
    #print(v)
    #print(e)
    #print(results)

EEPC_test['E_cell_number'] = cell_number['E']
EEPC_test['EPC_cell_number'] = cell_number['EPC']
EEPC_test['t-stats'] = eepc_tstats
EEPC_test['pval'] = eepc_pval
eepc_FDR = fdrcorrection(eepc_pval)
EEPC_test['FDR'] = eepc_FDR[1]

EEPC_test.to_csv('/lung_mets_data/finalized/EvsEPC_simple_T_test.csv')
