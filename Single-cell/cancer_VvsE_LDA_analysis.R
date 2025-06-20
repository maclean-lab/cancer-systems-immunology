library(scater)
library(tidyverse)
library(reticulate)
library(zellkonverter)
library(MASS)
library(ggplot2)
library(emdist)


# read data (turn h5ad file to sce)
sce <- readH5AD("/lung_mets_data/lung_mets_cancer_rd.h5ad")

cancer.ve = sce[, sce$treatment %in% c('V', 'E')]

# use the DE genes from scanpy and choose top genes for LDA (DE genes from seurat also work)
DEgenes <- read.csv('/lung_mets_cancer_DE_VvsE_remove_doublets_remove_Gm_mt.csv')
sorted_DEgenes = DEgenes[order(DEgenes$E_pvals_adj), ]
DElist <- sorted_DEgenes[1:250, 1]

# prepare for LDA
cancer.ve.de = cancer.ve[rownames(cancer.ve) %in% DElist,]
data <- cancer.ve.de@assays@data@listData[["log"]]
rownames(data) <- rownames(cancer.ve.de@assays@data@listData[["X"]])
data <- t(data)
treatment <- as.character(cancer.ve.de@colData@listData[["treatment"]])

## LDA 
#scale each predictor variable
X <- scale(data)

#make this reproducible
set.seed(1)

#Use 80% of dataset as training set and remaining 20% as testing set
sample <- sample(c(TRUE, FALSE), nrow(X), replace=TRUE, prob=c(0.8,0.2))
train <- as.data.frame(X[sample, ])
test <- as.data.frame(X[!sample, ])
train_y <- treatment[sample]
test_y <- treatment[!sample]

model <- lda(train_y~. , data = train)
predicted <- predict(model, test)
# test accuracy
mean(predicted$class==test_y) 
# training accuracy
predicted_train <- predict(model, train)
mean(predicted_train$class == train_y)

# visulization
lda_plot <- cbind(as.data.frame(predict(model)$x), as.data.frame(train_y))
#lda_plot$train_y <- sapply(lda_plot$train_y, as.numeric)

ggplot(lda_plot, aes(x = LD1, fill = train_y)) + 
  geom_histogram(alpha = 0.5, position = "identity")

# identify top genes that are used for separation
df<- as.data.frame(model$scaling)
gene_df = data.frame(matrix(nrow = nrow(df), ncol = 2))
gene_df['gene'] <- rownames(df)
gene_df['LD1'] <- df['LD1']
for (i in c(1:nrow(gene_df))){
  gene_df[i, 3] <- gsub('[`]', '', gene_df[i, 3])
}
sort_gene_df <- gene_df[order(abs(gene_df$LD1), decreasing = TRUE), ]
sort_gene_df$gene <- factor(sort_gene_df$gene, levels = sort_gene_df$gene[order(-abs(sort_gene_df$LD1))])


## Top 10 DE genes from LDA perform LDA
set.seed(10)
top10 <- sort_gene_df[1:10,]$gene
top10 <- as.character(top10)


cancer.top10.de = cancer.ve[rownames(cancer.ve) %in% top10,]
data <- cancer.top10.de@assays@data@listData[["log"]]
rownames(data) <- rownames(cancer.top10.de@assays@data@listData[["X"]])
data <- t(data)
treatment <- as.character(cancer.top10.de@colData@listData[["treatment"]])
X <- scale(data)

sample <- sample(c(TRUE, FALSE), nrow(X), replace=TRUE, prob=c(0.8,0.2))
train <- as.data.frame(X[sample, ])
test <- as.data.frame(X[!sample, ])
train_y <- treatment[sample]
test_y <- treatment[!sample]

model <- lda(train_y~. , data = train)
predicted <- predict(model, test)

# test accuracy
mean(predicted$class==test_y) 
# training accuracy
predicted_train <- predict(model, train)
mean(predicted_train$class == train_y)

# visulization
lda_plot_e <- cbind(as.data.frame(predict(model)$x), as.data.frame(train_y))

ve <- ggplot(lda_plot_e, aes(x = LD1, fill = train_y)) + 
  geom_histogram(alpha = 0.6, position = "identity") + theme_bw() + 
  theme(text = element_text(family = "sans", size = 10),
        axis.title = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  scale_fill_manual(values = c("#202d4f", "#b5b5b1"))

ggsave("/LDA_plot/LDA_cancer_hist_ve.jpg", width = 5, height = 4)

# Earth mover's distance
v_LD1 <- cbind(rep(1, times = length(predict(model)$x[train_y == 'V'])), predict(model)$x[train_y == 'V'])
e_LD1 <- cbind(rep(1, times = length(predict(model)$x[train_y == 'E'])), predict(model)$x[train_y == 'E'])

emd_dist <- emd(as.matrix(v_LD1), as.matrix(e_LD1), dist="euclidean")
emd_dist
