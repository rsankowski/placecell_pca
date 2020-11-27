#dir.create("bin")
#dir.create("data")
library(tidyverse)
library(readxl)
library(Rtsne)
library(umap)

ctrl_cells <- read_excel("data/20201124_Roman_PlaceCells.xlsx", sheet = 1, col_names = F)
mumt_cells <- read_excel("data/20201124_Roman_PlaceCells.xlsx", sheet = 2, col_names = F)
sham_cells <- read_excel("data/20201124_Roman_PlaceCells.xlsx", sheet = 3, col_names = F)
clp_cells <- read_excel("data/20201124_Roman_PlaceCells.xlsx", sheet = 4, col_names = F)

ctrl_list <- list()
for (i in 0:(ncol(ctrl_cells)/20-1)) {
  min <- 1 + 20*i
  max <- 20 + 20*i
  ctrl_list[[i+1]] <- as.matrix(ctrl_cells[,min:max])
}

mumt_list <- list()
for (i in 0:(ncol(mumt_cells)/20-1)) {
  min <- 1 + 20*i
  max <- 20 + 20*i
  mumt_list[[i+1]] <- as.matrix(mumt_cells[,min:max])
}

sham_list <- list()
for (i in 0:(ncol(sham_cells)/20-1)) {
  min <- 1 + 20*i
  max <- 20 + 20*i
  sham_list[[i+1]] <- as.matrix(sham_cells[,min:max])
}

clp_list <- list()
for (i in 0:(ncol(clp_cells)/20-1)) {
  min <- 1 + 20*i
  max <- 20 + 20*i
  clp_list[[i+1]] <- as.matrix(clp_cells[,min:max])
}

#flatten the matrices and connect the membranes
ctrl_list <- lapply(ctrl_list, as.vector) %>%
  as.data.frame()
colnames(ctrl_list) <- paste0("ctrl", 1:ncol(ctrl_list))

mumt_list <- lapply(mumt_list, as.vector) %>%
  as.data.frame()
colnames(mumt_list) <- paste0("mumt", 1:ncol(mumt_list))

sham_list <- lapply(sham_list, as.vector) %>%
  as.data.frame()
colnames(sham_list) <- paste0("sham", 1:ncol(sham_list))

clp_list <- lapply(clp_list, as.vector) %>%
  as.data.frame()
colnames(clp_list) <- paste0("clp", 1:ncol(clp_list))

all_list <- cbind(ctrl_list, mumt_list, sham_list, clp_list)

pca_lst <- prcomp(t(as.matrix(all_list)), scale. = T)

#scree plot
plot(pca_lst)
eigs <- pca_lst$sdev^2
plot(1:20, eigs[1:20]/sum(eigs))

#tsne
tsne_lst <- Rtsne(t(as.matrix(all_list)))
tsne_lst_10pc <- Rtsne(pca_lst$x[,1:10])

#umap
umap_lst <- umap(t(as.matrix(all_list)))
umap_lst_10pc <- umap(pca_lst$x[,1:10])

#build data frame
df <- data.frame(pca_lst$x[,1:2], condition=gsub("[0-9]+", "",rownames(pca_lst$x)))
df_tsne <- data.frame(tsne_lst$Y, condition=gsub("[0-9]+", "",rownames(pca_lst$x)))
df_tsne_10pc <- data.frame(tsne_lst_10pc$Y, condition=gsub("[0-9]+", "",rownames(pca_lst$x)))
df_umap <- data.frame(umap_lst$layout, condition=gsub("[0-9]+", "",rownames(pca_lst$x)))
df_umap_10pc <- data.frame(umap_lst_10pc$layout, condition=gsub("[0-9]+", "",rownames(pca_lst$x)))

#Exploratory plots
#pca
ggplot(df, aes(PC1, PC2, color=condition)) +
  geom_point()

#plot tsne
ggplot(df_tsne, aes(df_tsne[[1]], df_tsne[[2]], color=condition)) +
  geom_point()

ggplot(df_tsne_10, aes(df_tsne[[1]], df_tsne[[2]], color=condition)) +
  geom_point()

#plot umap
ggplot(df_umap, aes(df_umap[[1]], df_umap[[2]], color=condition)) +
  geom_point()

ggplot(df_umap_10pc, aes(df_umap_10pc[[1]], df_umap_10pc[[2]], color=condition)) +
  geom_point()

#nice plot
ggplot(df_umap, aes(df_umap[[1]], df_umap[[2]], fill=condition)) +
  geom_point(pch=21, stroke=0.25, size=6) +
  theme_void() +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "UMAP of individual place cells", x="UMAP1", y="UMAP2")

ggsave("plots/umap_full_dataset.pdf", useDingbats=F)

#tsne
ggplot(df_tsne, aes(df_tsne[[1]], df_tsne[[2]], fill=condition)) +
  geom_point(pch=21, stroke=0.25, size=6) +
  theme_void() +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "t-SNE of individual place cells")

ggsave("plots/tsne_full_dataset.pdf", useDingbats=F)

