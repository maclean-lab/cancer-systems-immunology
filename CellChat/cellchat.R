#load libraries
library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)

#set working directory to where scanpy/seurat data objects are stored
#if data is from scanpy (python), then it will need to be converted to a seurat (R) object (https://satijalab.org/seurat/archive/v2.4/conversion_vignette)
setwd("~/filepath") 

#load data based on which cell types to include
adata <- LoadH5Seurat("adata.h5seurat")
adata
adata$clusters

#define which cell types to use in analysis
Idents(object = adata) <- "clusters"
adata$clusters
adata$clusters <- factor(adata$clusters, levels=c("Cancer", "Endothelial", "Lipofibroblasts", "Mature myeloid", "MDSCs", "T cells", "B cells", "NK cells"))

#break into treatment groups
adata_V = subset(x = adata, subset = treatment == "V")
adata_V

adata_E = subset(x = adata, subset = treatment == "E")
adata_E

adata_EPC = subset(x = adata, subset = treatment == "EPC")
adata_EPC

adata_PC = subset(x = adata, subset = treatment == "PC")
adata_PC

#run CellChat
cellchat_V <- createCellChat(object = adata_V, group.by = "clusters")
cellchat_E <- createCellChat(object = adata_E, group.by = "clusters")
cellchat_EPC <- createCellChat(object = adata_EPC, group.by = "clusters")
cellchat_PC <- createCellChat(object = adata_PC, group.by = "clusters")

show(table(cellchat_V@idents))
show(table(cellchat_E@idents))
show(table(cellchat_EPC@idents))
show(table(cellchat_PC@idents))

CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)

#show the database
dplyr::glimpse(CellChatDB$interaction)

CellChatDB.use <- CellChatDB

#set the used database in the object
cellchat_V@DB <- CellChatDB.use
cellchat_E@DB <- CellChatDB.use
cellchat_EPC@DB <- CellChatDB.use
cellchat_PC@DB <- CellChatDB.use

cellchat_V <- subsetData(cellchat_V)
cellchat_E <- subsetData(cellchat_E)
cellchat_EPC <- subsetData(cellchat_EPC)
cellchat_PC <- subsetData(cellchat_PC)

future::plan("multisession", workers = 6) # do parallel

cellchat_V <- identifyOverExpressedGenes(cellchat_V)
cellchat_V <- identifyOverExpressedInteractions(cellchat_V)

cellchat_E <- identifyOverExpressedGenes(cellchat_E)
cellchat_E <- identifyOverExpressedInteractions(cellchat_E)

cellchat_EPC <- identifyOverExpressedGenes(cellchat_EPC)
cellchat_EPC <- identifyOverExpressedInteractions(cellchat_EPC)

cellchat_PC <- identifyOverExpressedGenes(cellchat_PC)
cellchat_PC <- identifyOverExpressedInteractions(cellchat_PC)

cellchat_V <- projectData(cellchat_V, PPI.mouse)
cellchat_E <- projectData(cellchat_E, PPI.mouse)
cellchat_EPC <- projectData(cellchat_EPC, PPI.mouse)
cellchat_PC <- projectData(cellchat_PC, PPI.mouse)

cellchat_V <- computeCommunProb(cellchat_V)
cellchat_E <- computeCommunProb(cellchat_E)
cellchat_EPC <- computeCommunProb(cellchat_EPC)
cellchat_PC <- computeCommunProb(cellchat_PC)

cellchat_V <- filterCommunication(cellchat_V, min.cells = 1)
cellchat_E <- filterCommunication(cellchat_E, min.cells = 1)
cellchat_EPC <- filterCommunication(cellchat_EPC, min.cells = 1)
cellchat_PC <- filterCommunication(cellchat_PC, min.cells = 1)

cellchat_V <- computeCommunProbPathway(cellchat_V)
cellchat_E <- computeCommunProbPathway(cellchat_E)
cellchat_EPC <- computeCommunProbPathway(cellchat_EPC)
cellchat_PC <- computeCommunProbPathway(cellchat_PC)

cellchat_V <- aggregateNet(cellchat_V)
cellchat_E <- aggregateNet(cellchat_E)
cellchat_EPC <- aggregateNet(cellchat_EPC)
cellchat_PC <- aggregateNet(cellchat_PC)

future.seed=TRUE

cellchat_V <- netAnalysis_computeCentrality(cellchat_V, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
cellchat_E <- netAnalysis_computeCentrality(cellchat_E, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
cellchat_EPC <- netAnalysis_computeCentrality(cellchat_EPC, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
cellchat_PC <- netAnalysis_computeCentrality(cellchat_PC, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways

object.list_E <- list(V = cellchat_V, E = cellchat_E)
object.list_EPC <- list(V = cellchat_V, EPC = cellchat_EPC)
object.list_PC <- list(V = cellchat_V, PC = cellchat_PC)
object.list_all <- list(V = cellchat_V, E = cellchat_E, EPC = cellchat_EPC, PC = cellchat_PC)
cellchat_comparison_E <- mergeCellChat(object.list_E, add.names = names(object.list_E))
cellchat_comparison_EPC <- mergeCellChat(object.list_EPC, add.names = names(object.list_EPC))
cellchat_comparison_PC <- mergeCellChat(object.list_PC, add.names = names(object.list_PC))
cellchat_comparison_all <- mergeCellChat(object.list_all, add.names = names(object.list_all))

#merge the following slots for treatment comparison: 'data.signaling','images','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'
cellchat_comparison_E
cellchat_comparison_EPC
cellchat_comparison_PC

#create signaling tables
df_V.net <- subsetCommunication(cellchat_V)
df_E.net <- subsetCommunication(cellchat_E)
df_EPC.net <- subsetCommunication(cellchat_EPC)
df_PC.net <- subsetCommunication(cellchat_PC)

write.table(df_V.net, file = 'V.csv', col.names = TRUE, row.names = FALSE, sep = ",")
write.table(df_E.net, file = 'E.csv', col.names = TRUE, row.names = FALSE, sep = ",")
write.table(df_EPC.net, file = 'EPC.csv', col.names = TRUE, row.names = FALSE, sep = ",")
write.table(df_PC.net, file = 'PC.csv', col.names = TRUE, row.names = FALSE, sep = ",")

#heatmap based on a merged object
gg_E <- netVisual_heatmap(cellchat_comparison_E, measure = "weight", title = "Differential interaction strength: V vs E")
gg_EPC <- netVisual_heatmap(cellchat_comparison_EPC, measure = "weight", title = "Differential interaction strength: V vs EPC")
gg_PC <- netVisual_heatmap(cellchat_comparison_PC, measure = "weight", title = "Differential interaction strength: V vs PC")

pathways.show <- c("SELPLG") 
LR.show <- "SELPLG_SELL" # show one ligand-receptor pair

weight.max <- getMaxWeight(object.list_all, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets

plotcolors = c('#e6762b', '#2b9746', '#c22b2c', '#7e5146', '#866099', '#b4b637', '#286ca0', '#c970a1')

#show circos plots
gg1 <- netVisual_individual(cellchat_V, signaling = pathways.show, pairLR.use = LR.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 15, color.use = plotcolors)
gg2 <-netVisual_individual(cellchat_E, signaling = pathways.show, pairLR.use = LR.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 15, color.use = plotcolors)
gg3 <-netVisual_individual(cellchat_EPC, signaling = pathways.show, pairLR.use = LR.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 15, color.use = plotcolors)
gg4 <-netVisual_individual(cellchat_PC, signaling = pathways.show, pairLR.use = LR.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 15, color.use = plotcolors)

cellchat_comparison_all@meta$datasets = factor(cellchat_comparison_all@meta$datasets, levels = c("V", "E", "EPC", "PC")) # set factor level
plotGeneExpression(cellchat_comparison_all, signaling = pathways.show, split.by = "datasets", colors.ggplot = T)

#show heatmaps
gg1 <- netVisual_heatmap(cellchat_V, measure = "weight", color.use = plotcolors, color.heatmap = c("#FFFFFF", "#b2182b"), font.size = 26, font.size.title = 26)
gg1
ggnew <- netVisual_circle(cellchat_V@net$weight, vertex.label.cex = 2, edge.weight.max = 5, edge.width.max = 5, color.use = plotcolors)
ggnew
