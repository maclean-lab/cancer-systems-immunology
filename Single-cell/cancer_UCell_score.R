library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(ggplot2)
library(UCell)

SaveFigure <- function(plots, name, type = "png", width, height, res){
  if(type == "png") {
    png(paste0(name, ".", type),
        width = width, height = height, units = "in", res = 200)
  } else {
    pdf(paste0(name, ".", type),
        width = width, height = height)
  }
  print(plots)
  dev.off()
}

Cancer <- LoadH5Seurat("lung_mets_data/lung_mets_cancer_rd.h5seurat")

Cancer$treatment <- factor(x = Cancer$treatment, levels = c('V', 'E', 'EP', 'EC', 'EPC', 'PC'))

markers <- list()

markers$stemness <- c('Cd44','Cd24a','Abcg2','Prom1','Itga6','Itgb3', 'Lgr5','B3galt5','Cd70','Procr',
                      'Thy1','Muc1','B4galnt1','Nectin4','Lin28a','Aldh1a1','Aldh1a3','Epcam')

Cancer <- AddModuleScore_UCell(Cancer, features = markers)
signature.names <- paste0(names(markers), "_UCell")

vp1 <- VlnPlot(Cancer, features = signature.names, group.by = "treatment", 
               cols = c('#b5b5b1', '#202d4f', '#7092bf', '#7ec6bb', '#682a7a', '#e0555f'))
SaveFigure(vp1, "stemness_Ucell", width = 6, height = 6)

#PAM50
markers <- list()
markers$pam50_basal = c('Pttg1', 'Cdc20', 'Orc6', 'Kif2c', 'Ube2c', 'Melk', 'Birc5', 'Nuf2', 'Cep55', 'Exo1', 'Cenpf', 'Ndc80', 'Tyms',
                        'Ube2t', 'Anln', 'Ccnb1', 'Rrm2', 'Mki67')

Cancer <- AddModuleScore_UCell(Cancer, features = markers)
signature.names <- paste0(names(markers), "_UCell")

vp2 <- VlnPlot(Cancer, features = signature.names, group.by = "treatment",
               cols = c('#b5b5b1', '#202d4f', '#7092bf', '#7ec6bb', '#682a7a', '#e0555f'))
SaveFigure(vp2, "PAM50_Basal_Ucell", width = 6, height = 6)



markers <- list()
markers$pam50_myo = c('Mia', 'Foxc1', 'Actr3b', 'Ccne1', 'Phgdh', 'Sfrp1', 'Myc', 'Cdh3', 'Krt5', 'Krt17', 'Krt14', 'Egfr')

Cancer <- AddModuleScore_UCell(Cancer, features = markers)
signature.names <- paste0(names(markers), "_UCell")

vp2 <- VlnPlot(Cancer, features = signature.names, group.by = "treatment",
               cols = c('#b5b5b1', '#202d4f', '#7092bf', '#7ec6bb', '#682a7a', '#e0555f'))
SaveFigure(vp2, "PAM50_Myo_Ucell", width = 6, height = 6)


markers <- list()
markers$pam50_luminal = c('Cdc6', 'Esr1', 'Bcl2', 'Foxa1', 'Cxxc5','Mlph','Mapt', 'Slc39a6', 'Nat1', 'Mdm2', 'Pgr', 'Mmp11', 'Blvra')

Cancer <- AddModuleScore_UCell(Cancer, features = markers)
signature.names <- paste0(names(markers), "_UCell")

vp2 <- VlnPlot(Cancer, features = signature.names, group.by = "treatment",
               cols = c('#b5b5b1', '#202d4f', '#7092bf', '#7ec6bb', '#682a7a', '#e0555f'))
SaveFigure(vp2, "PAM50_Luminal_Ucell", width = 6, height = 6)



markers <- list()
markers$pam50_her = c('Fgfr4','Grb7', 'Erbb2')

Cancer <- AddModuleScore_UCell(Cancer, features = markers)
signature.names <- paste0(names(markers), "_UCell")

vp2 <- VlnPlot(Cancer, features = signature.names, group.by = "treatment",
               cols = c('#b5b5b1', '#202d4f', '#7092bf', '#7ec6bb', '#682a7a', '#e0555f'))
SaveFigure(vp2, "PAM50_Her2_Ucell", width = 6, height = 6)



#EMT genes
markers <- list()
markers$mesenchymal = c('Vim', 'Cdh2', 'Foxc2', 'Snai1', 'Snai2', 'Twist1', 'Fn1', 'Itgb6', 'Mmp2', 'Mmp3', 'Mmp9', 'Sox10', 
                        'Gclc')
markers$epithelial = c('Cdh1', 'Dsp', 'Ocln')

Cancer <- AddModuleScore_UCell(Cancer, features = markers)
signature.names <- paste0(names(markers), "_UCell")

vp3<-VlnPlot(Cancer, features = signature.names, group.by = "treatment",
             cols = c('#b5b5b1', '#202d4f', '#7092bf', '#7ec6bb', '#682a7a', '#e0555f'))
SaveFigure(vp3, "EMT_Ucell", width = 12, height = 6)

## test significance
Cancer_V = subset(x = Cancer, subset = treatment == "V")
Cancer_E = subset(x = Cancer,  subset = treatment == "E")
Cancer_EPC = subset(x = Cancer,  subset = treatment == "EPC")
Cancer_PC = subset(x = Cancer,  subset = treatment == "PC")
Cancer_EC = subset(x = Cancer,  subset = treatment == "EC")
Cancer_EP = subset(x = Cancer,  subset = treatment == "EP")


### stemness
t.test(Cancer_V$stemness_UCell, Cancer_E$stemness_UCell)
t.test(Cancer_V$stemness_UCell, Cancer_EPC$stemness_UCell)
t.test(Cancer_V$stemness_UCell, Cancer_EP$stemness_UCell)
t.test(Cancer_V$stemness_UCell, Cancer_EC$stemness_UCell)
t.test(Cancer_V$stemness_UCell, Cancer_PC$stemness_UCell)

### pam50
t.test(Cancer_V$pam50_UCell, Cancer_E$pam50_UCell)
t.test(Cancer_V$pam50_UCell, Cancer_EPC$pam50_UCell)
t.test(Cancer_V$pam50_UCell, Cancer_EP$pam50_UCell)
t.test(Cancer_V$pam50_UCell, Cancer_EC$pam50_UCell)
t.test(Cancer_V$pam50_UCell, Cancer_PC$pam50_UCell)

### pam50 Basal
t.test(Cancer_V$pam50_basal_UCell, Cancer_E$pam50_basal_UCell)
t.test(Cancer_V$pam50_basal_UCell, Cancer_EPC$pam50_basal_UCell)
t.test(Cancer_V$pam50_basal_UCell, Cancer_EP$pam50_basal_UCell)
t.test(Cancer_V$pam50_basal_UCell, Cancer_EC$pam50_basal_UCell)
t.test(Cancer_V$pam50_basal_UCell, Cancer_PC$pam50_basal_UCell)

### pam50 Myo
t.test(Cancer_V$pam50_myo_UCell, Cancer_E$pam50_myo_UCell)
t.test(Cancer_V$pam50_myo_UCell, Cancer_EPC$pam50_myo_UCell)
t.test(Cancer_V$pam50_myo_UCell, Cancer_EP$pam50_myo_UCell)
t.test(Cancer_V$pam50_myo_UCell, Cancer_EC$pam50_myo_UCell)
t.test(Cancer_V$pam50_myo_UCell, Cancer_PC$pam50_myo_UCell)

### pam50 Luminal
t.test(Cancer_V$pam50_luminal_UCell, Cancer_E$pam50_luminal_UCell)
t.test(Cancer_V$pam50_luminal_UCell, Cancer_EPC$pam50_luminal_UCell)
t.test(Cancer_V$pam50_luminal_UCell, Cancer_EP$pam50_luminal_UCell)
t.test(Cancer_V$pam50_luminal_UCell, Cancer_EC$pam50_luminal_UCell)
t.test(Cancer_V$pam50_luminal_UCell, Cancer_PC$pam50_luminal_UCell)

### pam50 Her2
t.test(Cancer_V$pam50_her_UCell, Cancer_E$pam50_her_UCell)
t.test(Cancer_V$pam50_her_UCell, Cancer_EPC$pam50_her_UCell)
t.test(Cancer_V$pam50_her_UCell, Cancer_EP$pam50_her_UCell)
t.test(Cancer_V$pam50_her_UCell, Cancer_EC$pam50_her_UCell)
t.test(Cancer_V$pam50_her_UCell, Cancer_PC$pam50_her_UCell)


### Messenchymal marker
t.test(Cancer_V$mesenchymal_UCell, Cancer_E$mesenchymal_UCell)
t.test(Cancer_V$mesenchymal_UCell, Cancer_EPC$mesenchymal_UCell)
t.test(Cancer_V$mesenchymal_UCell, Cancer_EP$mesenchymal_UCell)
t.test(Cancer_V$mesenchymal_UCell, Cancer_EC$mesenchymal_UCell)
t.test(Cancer_V$mesenchymal_UCell, Cancer_PC$mesenchymal_UCell)