##### Use autoencoder to embed DepMap genomics data
library(tidyverse)
library(data.table)
library(ggpubr)
library(h2o)

setwd("~/Google Drive/My Drive/JCBC/Projects/20251120_DepMap_Autoencoder/")
dir.create("./results/plots", recursive = TRUE, showWarnings = FALSE)
dir.create("./results/output_files", recursive = TRUE, showWarnings = FALSE)
dir.create("./data", recursive = TRUE, showWarnings = FALSE)
dir.create("./codes", recursive = TRUE, showWarnings = FALSE)

##### Aim 1: Use RNA expression only for AE embedding

# load data: RNA expression
setwd("~/Library/CloudStorage/GoogleDrive-jaychung0524@gmail.com/My Drive/JCBC/Light_Horse/Bioinformatics/tools/DepMap_24Q2")
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