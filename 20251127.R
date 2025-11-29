##### Use autoencoder to embed DepMap genomics data
library(tidyverse)
library(data.table)
library(ggpubr)
library(h2o)

##### Aim 1: Use RNA expression only for AE embedding

# load data: RNA expression
load("CCLE_24Q2_GE_match_sample_info.RData")
load("sample_info_match_biomarkers.RData")

# remove low variance genes
ge_var <- apply(ccle_ge_match_sam, 2, var, na.rm = TRUE) # calculate variance for each gene
plot(density(ge_var)) # plot variance density plot
abline(v = quantile(ge_var, 0.2)) # select a cutoff
sum(ge_var > quantile(ge_var, 0.2)) # how many genes passed variance cutoff: 15322
ge_dat_hv <- ccle_ge_match_sam[, which(ge_var > quantile(ge_var, 0.2))] # select high variance genes

# remove any cell lines (rows) that contains NAs
ge_dat_hv <- ge_dat_hv |> na.omit()

# standardize the data
ge_dat_hv <- scale(ge_dat_hv) |> as.data.frame()

# Hyperparameter search grid to find best AE model
h2o.init(nthreads = 20, max_mem_size = "60G")
rna_h2o <- as.h2o(ge_dat_hv) # samples as rows, genes as columns
hyper_grid <- list(hidden = list(
  c(1000, 500, 1000),
  c(1000, 500, 200, 500, 1000), 
  c(1000, 100, 50, 100, 1000),
  c(1000, 500, 100, 500, 1000), 
  c(2000, 1000, 200, 1000, 2000), 
  c(1000, 500, 300, 500, 1000)
))
ae_grid <- h2o.grid(
  algorithm = "deeplearning",
  x = colnames(rna_h2o), # genes as features
  training_frame = rna_h2o,
  grid_id = "RNA_ae1",
  autoencoder = TRUE,
  activation = "TanhWithDropout", # tried "Rectifier", but gave a lot of failed models maybe due to exploding gradients
  # hidden_dropout_ratios = c(0.5), # dropout to avoid overfitting
  # initial_weight_distribution = "Uniform",
  hyper_params = hyper_grid,
  # rate = 0.005, # don't set if using adaptive learning rate (default)
  nesterov_accelerated_gradient = TRUE, # NAG is faster than classical momentum
  epochs = 30,
  stopping_rounds = 5, # stop if no improvement for 5 epochs
  # sparse = TRUE,
  # ignore_const_cols = FALSE,
  seed = 524
)
ae_grid
# Print grid details
h2o.getGrid("RNA_ae1", sort_by = "mse", decreasing = FALSE)
# save as a text file
write.table(
  h2o.getGrid("RNA_ae1", sort_by = "mse", decreasing = FALSE)@summary_table,
  file = "./results/output_files/RNA_AE_hyperparameter_grid_search_results.txt",
  sep = "\t", row.names = FALSE, quote = FALSE
)
# hidden	model_ids	mse
# [1000, 100, 50, 100, 1000]	RNA_ae1_model_3	0.0192287161230024
# [1000, 500, 100, 500, 1000]	RNA_ae1_model_4	0.0192292069522191
# [1000, 500, 200, 500, 1000]	RNA_ae1_model_2	0.0192297674229018
# [1000, 500, 300, 500, 1000]	RNA_ae1_model_7	0.0192298687808162
# [2000, 1000, 200, 1000, 2000]	RNA_ae1_model_5	0.0192308907700851
# [1000, 500, 1000]	RNA_ae1_model_1	0.0192826357160933
# decided to use: hidden = c(1000, 500, 200, 500, 1000) for reasonable complexity and good performance

# AE learning with best parameters
ae_model <- h2o.deeplearning(
  x = colnames(rna_h2o), # genes as features
  training_frame = rna_h2o,
  autoencoder = TRUE,
  hidden = c(1000, 500, 200, 500, 1000),
  hidden_dropout_ratios = c(0.5, 0.5, 0.5, 0.5, 0.5),
  activation = "TanhWithDropout",
  epochs = 30,
  stopping_rounds = 5, # stop if no improvement for 5 epochs
  nesterov_accelerated_gradient = TRUE, # NAG is faster than classical momentum
  seed = 524
)
ae_model
# extract the compressed features
compressed_features <- h2o.deepfeatures(ae_model, rna_h2o, layer = 3) |> as.data.frame()
rownames(compressed_features) <- rownames(ge_dat_hv)
# save compressed features
setwd("~/Google Drive/My Drive/JCBC/Projects/20251120_DepMap_Autoencoder/")
fwrite(
  compressed_features,
  file = "./results/output_files/RNA_AE_compressed_features.csv",
  sep = ",", row.names = TRUE, quote = FALSE
)
h2o.shutdown()

##### Aim 2: Use RNA + mutation + CRISPR for AE embedding

# load data: RNA expression, mutation, CRISPR
load("CCLE_24Q2_GE_match_sample_info.RData")
load("DepMap_24Q2_Chronos_match_sample_info.RData")
load("CCLE_24Q2_DAMMUT_match_sample_info.RData")
load("CCLE_24Q2_HOTMUT_match_sample_info.RData")
load("sample_info_match_biomarkers.RData")

# RNA: remove low variance genes
ge_var <- apply(ccle_ge_match_sam, 2, var, na.rm = TRUE) # calculate variance for each gene
plot(density(ge_var)) # plot variance density plot
abline(v = quantile(ge_var, 0.2)) # select a cutoff
sum(ge_var > quantile(ge_var, 0.2)) # how many genes passed variance cutoff: 15322
ge_dat_hv <- ccle_ge_match_sam[, which(ge_var > quantile(ge_var, 0.2))] # select high variance genes

# RNA:standardize the data
ge_dat_hv <- scale(ge_dat_hv) |> as.data.frame()

# DAMMUT: remove low instance mutations
sum(apply(na.omit(dammut_dat_match_sam), 2, function(x) sum(x != 0)) > 20) # 635 mutations with more than 20 instances
dam_mut_dat_filt <- dammut_dat_match_sam[, which(apply(na.omit(dammut_dat_match_sam), 2, function(x) sum(x != 0)) > 20)]
colnames(dam_mut_dat_filt)

# HOTMUT: remove low instance mutations
sum(apply(na.omit(hotmut_dat_match_sam), 2, function(x) sum(x != 0)) > 20) # 18 mutations with more than 20 instances
hot_mut_dat_filt <- hotmut_dat_match_sam[, which(apply(na.omit(hotmut_dat_match_sam), 2, function(x) sum(x != 0)) > 20)]
colnames(hot_mut_dat_filt)

# CRISPR: remove low variance genes
crispr_var <- apply(chronos_dat_match_sam, 2, var, na.rm = TRUE) # calculate variance for each gene
plot(density(crispr_var)) # plot variance density plot
abline(v = quantile(crispr_var, 0.2)) # select a cutoff
sum(crispr_var > quantile(crispr_var, 0.2)) # how many genes passed variance cutoff: 14748
crispr_dat_hv <- chronos_dat_match_sam[, which(crispr_var > quantile(crispr_var, 0.2))] # select high variance genes

# CRISPR: standardize the data
crispr_dat_hv <- scale(crispr_dat_hv) |> as.data.frame()

# combine RNA + mutation + CRISPR data, adding suffix to column names to avoid duplicated names
colnames(ge_dat_hv) <- paste0(colnames(ge_dat_hv),"_RNA")
colnames(dam_mut_dat_filt) <- paste0(colnames(dam_mut_dat_filt),"_DAMMUT")
colnames(hot_mut_dat_filt) <- paste0(colnames(hot_mut_dat_filt),"_HOTMUT")
colnames(crispr_dat_hv) <- paste0(colnames(crispr_dat_hv),"_CRISPR")
multi_omics_dat <- cbind(
  ge_dat_hv,
  dam_mut_dat_filt,
  hot_mut_dat_filt,
  crispr_dat_hv
)
dim(multi_omics_dat) # 1959 samples, 30723 features

# remove any cell lines (rows) that contains > 80% NAs
multi_omics_dat <- multi_omics_dat[rowSums(is.na(multi_omics_dat)) <= ncol(multi_omics_dat) * 0.8, ]

# Hyperparameter search grid to find best AE model
h2o.init(nthreads = 20, max_mem_size = "60G")
rna_h2o <- as.h2o(multi_omics_dat) # samples as rows, genes as columns
hyper_grid <- list(hidden = list(
  c(1000, 500, 1000),
  c(1000, 500, 200, 500, 1000), 
  c(1000, 100, 50, 100, 1000), 
  c(1000, 500, 100, 500, 1000), 
  c(2000, 1000, 200, 1000, 2000), 
  c(1000, 500, 300, 500, 1000)
))
ae_grid <- h2o.grid(
  algorithm = "deeplearning",
  x = colnames(rna_h2o), # genes as features
  training_frame = rna_h2o,
  grid_id = "RNA_ae1",
  autoencoder = TRUE,
  activation = "TanhWithDropout", # tried "Rectifier", but gave a lot of failed models maybe due to exploding gradients
  # hidden_dropout_ratios = c(0.5), # dropout to avoid overfitting
  # initial_weight_distribution = "Uniform",
  hyper_params = hyper_grid,
  # rate = 0.005, 
  nesterov_accelerated_gradient = TRUE, # NAG is faster than classical momentum
  epochs = 30,
  stopping_rounds = 5, # stop if no improvement for 5 epochs
  # sparse = TRUE,
  # ignore_const_cols = FALSE,
  seed = 524
)
ae_grid
# Print grid details
h2o.getGrid("RNA_ae1", sort_by = "mse", decreasing = FALSE)
# save as a text file
write.table(
  h2o.getGrid("RNA_ae1", sort_by = "mse", decreasing = FALSE)@summary_table,
  file = "./results/output_files/MultiOmics_AE_hyperparameter_grid_search_results.txt",
  sep = "\t", row.names = FALSE, quote = FALSE
)
# hidden	model_ids	mse
# [1000, 100, 50, 100, 1000]	RNA_ae1_model_3	0.0163892185907279
# [1000, 500, 300, 500, 1000]	RNA_ae1_model_6	0.0163901573573774
# [1000, 500, 1000]	RNA_ae1_model_1	0.0164147952330418
# decided to use: hidden = c(1000, 500, 300, 500, 1000) for reasonable complexity and good performance

# AE learning with best parameters
ae_model <- h2o.deeplearning(
  x = colnames(rna_h2o), # genes as features
  training_frame = rna_h2o,
  autoencoder = TRUE,
  hidden = c(1000, 500, 300, 500, 1000),
  hidden_dropout_ratios = c(0.5, 0.5, 0.5, 0.5, 0.5),
  activation = "TanhWithDropout",
  epochs = 30,
  stopping_rounds = 5, # stop if no improvement for 5 epochs
  nesterov_accelerated_gradient = TRUE, # NAG is faster than classical momentum
  seed = 524
)
ae_model
# extract the compressed features
compressed_features <- h2o.deepfeatures(ae_model, rna_h2o, layer = 3) |> as.data.frame()
rownames(compressed_features) <- rownames(multi_omics_dat)
# save compressed features
fwrite(
  compressed_features,
  file = "./results/output_files/MultiOmics_AE_compressed_features.csv",
  sep = ",", row.names = TRUE, quote = FALSE
)
h2o.shutdown()

##### Aim 3: Compare RNA only or multi-omics AE embedding clustering with known lineages or common mutations

### RNA only AE UMAP
setwd("~/Library/CloudStorage/GoogleDrive-jaychung0524@gmail.com/My Drive/JCBC/Light_Horse/Bioinformatics/tools/DepMap_24Q2")
load("sample_info_match_biomarkers.RData")
setwd("~/Google Drive/My Drive/JCBC/Projects/20251120_DepMap_Autoencoder/")
RNA_AE_compressed_features <- fread(
  "./results/output_files/RNA_AE_compressed_features.csv",
  data.table = FALSE, 
) |> column_to_rownames(var = "V1")
# remove cell lines (rows) that are all zeros (NA data)
RNA_AE_compressed_features <- RNA_AE_compressed_features[rowSums(RNA_AE_compressed_features) != 0, ]
sample_info <- sample_info[match(rownames(RNA_AE_compressed_features), sample_info$StrippedCellLineName), ]
all(rownames(RNA_AE_compressed_features) == sample_info$StrippedCellLineName)
library(umap)
umap_dat <- umap(
  RNA_AE_compressed_features, 
  n_components = 2,
  n_neighbors = 50,
  random_state = 524,
  verbose = TRUE,
  n_threads = 15
)
umap_dat <- as.data.frame(umap_dat$layout)
umap_dat$cell_line <- sample_info$StrippedCellLineName
umap_dat$lineage <- sample_info$OncotreeLineage
umap_dat$subtype <- sample_info$OncotreeSubtype
library(ggrepel)
library(microViz)
ggplot(umap_dat, aes(V1, V2, color = lineage, fill = lineage)) +
  geom_point(size = 3, alpha = 0.5) + 
  scale_fill_manual(values = distinct_palette(n = length(unique(umap_dat$lineage)), pal = "brewerPlus", add = "lightgrey")) + 
  scale_color_manual(values = distinct_palette(n = length(unique(umap_dat$lineage)), pal = "brewerPlus", add = "lightgrey")) + 
  theme_classic(base_size = 20) + 
  # guides(fill = guide_legend(ncol = 1)) + 
  # add lineage labels at the center of each lineage cluster
  geom_label_repel(
    data = umap_dat %>%
      group_by(lineage) %>%
      summarize(V1 = median(V1), V2 = median(V2)),
    aes(label = lineage, color = lineage),
    size = 4,
    # color = "black",
    fill = "white",
    alpha = 0.7,
    show.legend = FALSE, 
    max.overlaps = 100
  ) +
  theme(legend.position = "none", 
        legend.title = element_text(size = 10), 
        legend.text = element_text(size = 10)
  ) +
  labs(x = "UMAP1", 
       y = "UMAP2", 
       title = "UMAP of AE compressed CCLE RNA expression")
ggsave("./results/plots/CCLE_RNA_AE_UMAP_lineage.png", width = 10, height = 8, dpi = 300)

### Full RNA without AE UMAP
setwd("~/Library/CloudStorage/GoogleDrive-jaychung0524@gmail.com/My Drive/JCBC/Light_Horse/Bioinformatics/tools/DepMap_24Q2")
load("CCLE_24Q2_GE_match_sample_info.RData")
load("sample_info_match_biomarkers.RData")
ge_var <- apply(ccle_ge_match_sam, 2, var, na.rm = TRUE) # calculate variance for each gene
plot(density(ge_var)) # plot variance density plot
abline(v = quantile(ge_var, 0.2)) # select a cutoff
sum(ge_var > quantile(ge_var, 0.2)) # how many genes passed variance cutoff: 15322
ge_dat_hv <- ccle_ge_match_sam[, which(ge_var > quantile(ge_var, 0.2))] # select high variance genes
# standardize the data
ge_dat_hv <- scale(ge_dat_hv) |> as.data.frame() |> na.omit()
sample_info <- sample_info[match(rownames(ge_dat_hv), sample_info$StrippedCellLineName), ]
all(rownames(ge_dat_hv) == sample_info$StrippedCellLineName)
library(umap)
umap_dat <- umap(
  ge_dat_hv, 
  n_components = 2,
  n_neighbors = 50,
  random_state = 524,
  verbose = TRUE,
  n_threads = 15
)
umap_dat <- as.data.frame(umap_dat$layout)
umap_dat$cell_line <- sample_info$StrippedCellLineName
umap_dat$lineage <- sample_info$OncotreeLineage
umap_dat$subtype <- sample_info$OncotreeSubtype
library(ggrepel)
library(microViz)
ggplot(umap_dat, aes(V1, -V2, color = lineage, fill = lineage)) +
  geom_point(size = 3, alpha = 0.5) + 
  scale_fill_manual(values = distinct_palette(n = length(unique(umap_dat$lineage)), pal = "brewerPlus", add = "lightgrey")) + 
  scale_color_manual(values = distinct_palette(n = length(unique(umap_dat$lineage)), pal = "brewerPlus", add = "lightgrey")) + 
  theme_classic(base_size = 20) + 
  # guides(fill = guide_legend(ncol = 1)) + 
  # add lineage labels at the center of each lineage cluster
  geom_label_repel(
    data = umap_dat %>%
      group_by(lineage) %>%
      summarize(V1 = median(V1), V2 = median(V2)),
    aes(label = lineage, color = lineage),
    size = 4,
    # color = "black",
    fill = "white",
    alpha = 0.7,
    show.legend = FALSE, 
    max.overlaps = 100
  ) +
  theme(legend.position = "none", 
        legend.title = element_text(size = 10), 
        legend.text = element_text(size = 10)
  ) +
  labs(x = "UMAP1", 
       y = "UMAP2", 
       title = "UMAP of full CCLE RNA expression (no AE)")
setwd("~/Google Drive/My Drive/JCBC/Projects/20251120_DepMap_Autoencoder/")
ggsave("./results/plots/CCLE_RNA_no_AE_UMAP_lineage.png", width = 10, height = 8, dpi = 300)

### Multi-omics AE UMAP
setwd("~/Library/CloudStorage/GoogleDrive-jaychung0524@gmail.com/My Drive/JCBC/Light_Horse/Bioinformatics/tools/DepMap_24Q2")
load("sample_info_match_biomarkers.RData")
setwd("~/Google Drive/My Drive/JCBC/Projects/20251120_DepMap_Autoencoder/")
MultiOmics_AE_compressed_features <- fread(
  "./results/output_files/MultiOmics_AE_compressed_features.csv",
  data.table = FALSE, 
) |> column_to_rownames(var = "V1")
# remove cell lines (rows) that are all zeros (NA data)
MultiOmics_AE_compressed_features <- MultiOmics_AE_compressed_features[rowSums(MultiOmics_AE_compressed_features) != 0, ]
sample_info <- sample_info[match(rownames(MultiOmics_AE_compressed_features), sample_info$StrippedCellLineName), ]
all(rownames(MultiOmics_AE_compressed_features) == sample_info$StrippedCellLineName)
library(umap)
umap_dat <- umap(
  MultiOmics_AE_compressed_features, 
  n_components = 2,
  n_neighbors = 50,
  random_state = 524,
  verbose = TRUE,
  n_threads = 15
)
umap_dat <- as.data.frame(umap_dat$layout)
umap_dat$cell_line <- sample_info$StrippedCellLineName
umap_dat$lineage <- sample_info$OncotreeLineage
umap_dat$subtype <- sample_info$OncotreeSubtype
library(ggrepel)
library(microViz)
ggplot(umap_dat, aes(-V1, -V2, color = lineage, fill = lineage)) +
  geom_point(size = 3, alpha = 0.5) + 
  scale_fill_manual(values = distinct_palette(n = length(unique(umap_dat$lineage)), pal = "brewerPlus", add = "lightgrey")) + 
  scale_color_manual(values = distinct_palette(n = length(unique(umap_dat$lineage)), pal = "brewerPlus", add = "lightgrey")) + 
  theme_classic(base_size = 20) + 
  # guides(fill = guide_legend(ncol = 1)) + 
  # add lineage labels at the center of each lineage cluster
  geom_label_repel(
    data = umap_dat %>%
      group_by(lineage) %>%
      summarize(V1 = median(V1), V2 = median(V2)),
    aes(label = lineage, color = lineage),
    size = 4,
    # color = "black",
    fill = "white",
    alpha = 0.7,
    show.legend = FALSE, 
    max.overlaps = 100
  ) +
  theme(legend.position = "none", 
        legend.title = element_text(size = 10), 
        legend.text = element_text(size = 10)
  ) +
  labs(x = "UMAP1", 
       y = "UMAP2", 
       title = "UMAP of AE compressed CCLE Multi-omics")
setwd("~/Google Drive/My Drive/JCBC/Projects/20251120_DepMap_Autoencoder/")
ggsave("./results/plots/CCLE_MultiOmics_AE_UMAP_lineage.png", width = 10, height = 8, dpi = 300)

### Conclusion 1: RNA deep features preserve original transcriptomic landscape across DepMap cell lines
### Conclusion 2: Multi-OMICS deep features reveal richer cell state structure than RNA alone
### Next steps: Prediction accuracy comparison on independent test data: deep features vs full OMICS
