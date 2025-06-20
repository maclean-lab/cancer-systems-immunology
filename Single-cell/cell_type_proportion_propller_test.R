library(Seurat)
library(SeuratData)
library(SeuratDisk)

library(speckle)
library(limma)
library(ggplot2)

MDSCs <- LoadH5Seurat("/lung_mets_data/finalized/MDSCs_subclusters.h5seurat")
Idents(object = MDSCs) <- 'clusters'


NKs <- LoadH5Seurat("/lung_mets_data/finalized/NK_subclusters.h5seurat")
Idents(object = NKs) <- 'clusters'


Ts <- LoadH5Seurat("/lung_mets_data/T_cells_subclusters.h5seurat")
Idents(object = Ts) <- 'clusters'

B <- LoadH5Seurat("/lung_mets_data/finalized/B_subclusters.h5seurat")
Idents(object = B) <- 'clusters'



MM <- LoadH5Seurat("/lung_mets_data/finalized/mature_myeloid_subclusters.h5seurat")
Idents(object = MM) <- 'clusters'

adata_all = merge(x = MDSCs, y = list(NKs, Ts, B, MM))
adata_all
table(Idents(adata_all))
write.table(table(Idents(adata_all)), file = '/lung_mets_data/CD45_cell_number.csv')

adata_all$sample <- adata_all$sample_name
adata_all$group <- adata_all$treatment
adata_all$clusters <- Idents(adata_all)

# V vs E
adata_VE <- subset(x = adata_all, subset = group %in% c("V", "E"))
table_VE <-propeller(x = adata_VE, sample = 'sample_name', group = 'treatment', transform="logit")


# V vs EPC
adata_VEPC <- subset(x = adata_all, subset = group %in% c("V", "EPC"))
table_VEPC <-propeller(x = adata_VEPC, sample = 'sample_name', group = 'treatment', transform="asin")



#V vs PC
adata_VPC <- subset(x = adata_all, subset = group %in% c("V", "PC"))
table_VPC <-propeller(x = adata_VPC, sample = 'sample_name', group = 'treatment', transform="logit")




# E vs EPC
adata_EEPC <- subset(x = adata_all, subset = group %in% c("E", "EPC"))
table_EEPC <-propeller(x = adata_EEPC, sample = 'sample_name', group = 'treatment', transform="logit")

write.table(table_VE, file = '/lung_mets_data/propeller_ttest_VE.csv', col.names = TRUE,
            row.names = FALSE, sep = ",")
write.table(table_VEPC, file = '/lung_mets_data/propeller_ttest_VEPC.csv', col.names = TRUE,
            row.names = FALSE, sep = ",")
write.table(table_VPC, file = '/lung_mets_data/propeller_ttest_VPC.csv', col.names = TRUE,
            row.names = FALSE, sep = ",")
write.table(table_EEPC, file = '/lung_mets_data/propeller_ttest_EEPC.csv', col.names = TRUE,
            row.names = FALSE, sep = ",")