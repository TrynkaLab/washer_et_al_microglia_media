# Initial QC of data - filter, normalisation 
# Input: Cellranger matrices
# Output: Seurat object(s)
.libPaths("/software/teamtrynka/ma23/R4.1/libs") # path to libraries

library(tidyverse)
library(Seurat)
library(patchwork)
library(future)
library(SCpubr)
# change the current plan to access parallelization
plan("multiprocess", workers = 14)
plan()
options(future.globals.maxSize= 4097152000) # 2Gb
options(stringsAsFactors = FALSE)

dir.create("../../data/results/1.3.QC", recursive = T)
cols = c("#0B00F6","#12130F","#3C71C6","#20F8F6","#AA00FF","#F40006")

#### functions --------

filter_seurat = function(seurat_object){
  
  message("Performing filter by number of genes and mitochondrial percentage.")
  seurat_object = subset(seurat_object, subset = nFeature_RNA > 1000  & percent.mt < 10 & nCount_RNA < 40000 & vireo_doublets == "singlet")
  message("Now the object has ", dim(seurat_object)[1], " genes and ", dim(seurat_object)[2], " cells.")
  return(seurat_object)
}


calculate_PCA_UMAP_neighbors_clusters_merged = function(seurat_object){
  seurat_object = RunPCA(seurat_object, verbose = FALSE)
  
  seurat_object = FindNeighbors(seurat_object, dims = 1:30, verbose = FALSE)
  seurat_object = FindClusters(seurat_object, 
                               resolution = 0.3, # Selected in 1.2.Determine_clustering_params - choosing lower loses CRM and pvM info
                               verbose = FALSE)
  
  seurat_object = RunUMAP(seurat_object, dims = 1:30, verbose = FALSE)
  
  
  
  return(seurat_object)
}

calculate_PCA_UMAP_neighbors_clusters_integrated = function(seurat_object){
  DefaultAssay(seurat_object) = "integrated"
  seurat_object = RunPCA(seurat_object, verbose = FALSE)
  
  seurat_object = FindNeighbors(seurat_object, dims = 1:30, verbose = FALSE)
  seurat_object = FindClusters(seurat_object, 
                               resolution = 0.2, # Selected in 1.1.Determine_clustering_params - choosing lower loses CRM and pvM info
                               verbose = FALSE)
  
  seurat_object = RunUMAP(seurat_object, dims = 1:30, verbose = FALSE)
  
  
  
  return(seurat_object)
}
UMAPS_all_libs = function(seurat_object) {
  p1 = DimPlot(seurat_object, label = TRUE) + NoLegend() + ggtitle("Clusters")
  p2 = DimPlot(seurat_object, label = F,  group.by = "orig.ident",
               cols = cols, shuffle = T)  + ggtitle("Library")
  p2[[1]]$layers[[1]]$aes_params$alpha = .4 
  p3 = FeaturePlot(seurat_object, label = F, features = "percent.mt", order = TRUE)  + ggtitle("Mit %") 
  
  p4 = FeaturePlot(seurat_object, label = F, features  = "nFeature_RNA", order = TRUE) + 
    ggtitle("Number of genes")
  p5 = FeaturePlot(seurat_object, label = F, features  = "nCount_RNA", order = TRUE) + 
    ggtitle("UMI counts")
  
  p6 = DimPlot(seurat_object, label = F,  group.by = "Phase")  + ggtitle("Cell cycle")
  p6[[1]]$layers[[1]]$aes_params$alpha = .4 
  
  
  p = (p1 | p2) / (p3 | p4) / (p5 | p6) 
  
  return(p)
}


### start of analysis ---------
# list of microglia input data from different media (all uncoated)
media_names = c("IGBN","IM","IMBN","ITGBN","ITM_ADMEM","ITMG")
# Loads a list of gene expression (feature) counts per cell
cellranger_data = list()
for(media in media_names){
  cellranger_data[[media]] =  Read10X(data.dir = paste0("../../data/",media,"/outs/filtered_feature_bc_matrix/"),
                                      strip.suffix = T)  ## remove trailing -1 from cell barcodes
  
}

# Initialize the Seurat object with the raw data. No filtering yet
seurat_list = list()
for(media in media_names){
  seurat_list[[media]] = CreateSeuratObject(counts = cellranger_data[[media]] , project = media)
  print(median(seurat_list[[media]]$nCount_RNA)) # Median UMI counts - as cellRanger summary
  print(median(seurat_list[[media]]$nFeature_RNA)) # Median genes per cell - as cellRanger summary
  seurat_list[[media]]$cell = paste(seurat_list[[media]]$orig.ident,rownames(seurat_list[[media]]@meta.data),sep="_")
  
}
rm(cellranger_data)

for(media in media_names){
  
  seurat_list[[media]][["percent.mt"]] = PercentageFeatureSet(seurat_list[[media]], pattern = "^MT-")
  
}

##### Check count metrics before filtering ####

input_merged = merge(seurat_list$IGBN, list(seurat_list$IM, seurat_list$IMBN, seurat_list$ITGBN,
                                            seurat_list$ITM_ADMEM, seurat_list$ITMG),
                     add.cell.ids=media_names)


# Visualize QC metrics as a violin plot - Features are genes and counts are UMIs
VlnPlot(input_merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)

# scatterplots
plot1 = FeatureScatter(input_merged, feature1 = "nCount_RNA", feature2 = "percent.mt") + NoLegend()  +
  ylab("% of mitochondrial genes") +
  xlab("UMI counts")
plot2 = FeatureScatter(input_merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  ylab("Number of genes") +
  xlab("UMI counts")

plot1 + plot2 

# density plots
p1 = input_merged@meta.data %>% 
  ggplot(aes(color=orig.ident, x=nCount_RNA, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 40000, color = "red") +
  ylab("Cell density") +
  xlab("UMI counts") + 
  NoLegend()
p2 = input_merged@meta.data %>% 
  ggplot(aes(color=orig.ident, x=nFeature_RNA, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 1000,  color = "red") +
  ylab("Cell density") +
  xlab("Number of genes") + 
  NoLegend()

p3 = input_merged@meta.data %>% 
  ggplot(aes(color=orig.ident, x=percent.mt, fill=orig.ident)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  geom_vline(xintercept = 10) +
  ylab("Cell density") +
  xlab("% of mitochondrial genes")


p1 + p2 + p3


######### End of count check
rm ( p1, p2, p3, plot1, plot2)

######### Filtering #########

### We filter out cells that have less than 1000 unique feature counts and over 40,000 UMI counts
### We filter out cells that have > 10% mitochondrial counts (the usual one for human cells [it depends on cell type and tissue], although many here seem above)

# doublet
vireo_doublets = read.csv("../../data/results/1.1.Determine_doublets/vireo_doublets_all.csv") # from 1.1.Determine_doublets.R

for(media in media_names){
  seurat_list[[media]]@meta.data$original_cell = rownames(seurat_list[[media]]@meta.data)
  seurat_list[[media]]@meta.data = merge(seurat_list[[media]]@meta.data,vireo_doublets, by = "cell" )
  seurat_list[[media]]@meta.data$vireo_doublets = "singlet"
  seurat_list[[media]]@meta.data[seurat_list[[media]]@meta.data$donor_id == "doublet","vireo_doublets"] = "doublet"
  seurat_list[[media]]@meta.data[seurat_list[[media]]@meta.data$donor_id == "unassigned","vireo_doublets"] = "unassigned"
  rownames(seurat_list[[media]]@meta.data) = seurat_list[[media]]@meta.data$original_cell
}

input_merged = merge(seurat_list$IGBN, list(seurat_list$IM, seurat_list$IMBN, seurat_list$ITGBN,
                                            seurat_list$ITM_ADMEM, seurat_list$ITMG),
                     add.cell.ids=media_names)
# saveRDS(input_merged, "../../data/results/1.3.QC/media_merged_unfiltered_unnormalised.rds")

for(media in media_names){
  seurat_list[[media]] = filter_seurat(seurat_list[[media]])
}

## Scoring genes by cell cycle gene expression - from Tirosh et al., 2016 ----------
# Pre-normalising and saving results because there are not enough bins if run on un-normalised data or SCTransformed

copy_list = seurat_list
for(media in media_names){
  copy_list[[media]] = NormalizeData(seurat_list[[media]],verbose = T)
  copy_list[[media]] = CellCycleScoring(copy_list[[media]],  s.features=cc.genes.updated.2019$s.genes, g2m.features=cc.genes.updated.2019$g2m.genes)
  seurat_list[[media]]$Phase = copy_list[[media]]$Phase
  seurat_list[[media]]$S.Score = copy_list[[media]]$S.Score
  seurat_list[[media]]$G2M.Score = copy_list[[media]]$G2M.Score
}

## Normalizing the data ----

# Merge data again for joint analysis - might preserve some technical batch effects between the libraries
input_merged = merge(seurat_list$IGBN, list(seurat_list$IM, seurat_list$IMBN, seurat_list$ITGBN,
                                            seurat_list$ITM_ADMEM, seurat_list$ITMG),
                     add.cell.ids=media_names)
input_merged = SCTransform(input_merged, verbose = TRUE)


# Transforming individual libraries - for integration later
for(media in media_names){
  seurat_list[[media]] = SCTransform(seurat_list[[media]], verbose = TRUE)
  
}

### Perform dimensionality reduction by PCA and UMAP embedding

input_merged= calculate_PCA_UMAP_neighbors_clusters_merged(input_merged)


p = UMAPS_all_libs(input_merged)
# png("../../data/results/1.3.QC/media_merged_UMAP_nothingRegressed.png", 
#     width = 10, height = 12, units = "in", res = 400)
p # ITM ADMEM and ITMG seem more similar among each other
# dev.off()
p1 = DimPlot(input_merged, label = F,  group.by = "donor_id", shuffle = T)  + ggtitle("Donor")
p1

# removing underscores for plots
input_merged@meta.data$orig.ident = input_merged@meta.data$orig.ident %>%
  dplyr::recode("ITM_ADMEM" = "ITM ADMEM")

media_names = c("IGBN","IM","IMBN","ITGBN","ITM ADMEM","ITMG")

names(cols) = media_names


p1 =  SCpubr::do_DimPlot(sample = input_merged, 
                         plot.title = "Media",
                         dims = c(1, 2), 
                        group.by = "orig.ident",
                         colors.use = cols, 
                        shuffle = T)




png("../../data/results/1.3.QC/media_merged_UMAP_nothingRegressed_libraries_Fig5a.png", 
    width = 5.5, height = 4, units = "in", res = 400)
p1
dev.off()

p1 =  SCpubr::do_DimPlot(sample = input_merged, 
                         plot.title = "Clusters",
                         dims = c(1, 2),legend.ncol = 2)

p2 =  SCpubr::do_DimPlot(sample = input_merged, 
                         plot.title = "Cell cycle phase",
                         dims = c(1, 2),group.by = "Phase",
                         order = "G1", shuffle = F)

p3 =  SCpubr::do_DimPlot(sample = input_merged, 
                         plot.title = "Media",
                         reduction = "pca",
                         dims = c(1, 2), 
                         split.by = "orig.ident",
                         colors.use = cols, 
                         shuffle = T, legend = F)

png("../../data/results/1.3.QC/media_merged_UMAP_nothingRegressed_clusters_cellCycle_FigS5a.png", 
    width = 10, height = 9, units = "in", res = 400)
(p1 + p2) /  p3 +   plot_layout(heights = c(1, 2))
dev.off()

saveRDS(input_merged, "../../data/results/1.3.QC/media_merged_noRegression.rds")

rm(input_merged)
gc()
### Integrated - done here after removing merged to free up memory ----------

# Integrate data for joint analysis - might obfuscate some of the subtler biological differences between the experiments
# Most appropriate for per-cluster analysis

media_features = SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 3000)
seurat_list = PrepSCTIntegration(object.list = seurat_list, anchor.features = media_features, 
                                 verbose = TRUE)

plan("multiprocess", workers = 1) # To avoid serialize error

media_anchors = FindIntegrationAnchors(object.list = seurat_list, normalization.method = "SCT", 
                                       anchor.features = media_features, verbose = TRUE)
gc()
media_integrated = IntegrateData(anchorset = media_anchors, normalization.method = "SCT", 
                                 verbose = TRUE)


# PCA
media_integrated = calculate_PCA_UMAP_neighbors_clusters_integrated(media_integrated)

p = UMAPS_all_libs(media_integrated)
# png("../../data/results/1.3.QC/media_integrated_UMAP_nothingRegressed.png", 
#     width = 10, height = 12, units = "in", res = 400)
p
# dev.off()

media_integrated@meta.data$orig.ident = media_integrated@meta.data$orig.ident %>%
  dplyr::recode("ITM_ADMEM" = "ITM ADMEM")

p1 =  SCpubr::do_DimPlot(sample = media_integrated, 
                         plot.title = "Clusters",
                         dims = c(1, 2))

p2 =  SCpubr::do_DimPlot(sample = media_integrated, 
                         plot.title = "Cell cycle phase",
                         dims = c(1, 2),group.by = "Phase")


png("../../data/results/1.3.QC/media_integrated_UMAP_nothingRegressed_clusters_cellCycle.png", 
    width = 10, height = 4, units = "in", res = 400)
p1 + p2
dev.off()

p1 =  SCpubr::do_DimPlot(sample = media_integrated, 
                         plot.title = "Media",
                         dims = c(1, 2), 
                         group.by = "orig.ident",
                         colors.use = cols, shuffle = T)


png("../../data/results/1.3.QC/media_integrated_UMAP_nothingRegressed_libraries.png", 
    width = 5.5, height = 4, units = "in", res = 400)
p1
dev.off()


# checking donors
p1 = DimPlot(media_integrated, label = F,  group.by = "donor_id", shuffle = T)  + ggtitle("Donor")

p2 = DimPlot(media_integrated, label = F,  group.by = "donor_id", shuffle = T, reduction = "pca")  + ggtitle("Donor")

p3 = DimPlot(media_integrated, label = F,  group.by = "Phase", shuffle = T, reduction = "pca")  + ggtitle("Cell cycle")
# png("../../data/results/1.3.QC/media_integrated_UMAP_nothingRegressed_donors.png", 
#     width = 15, height = 5, units = "in", res = 400)
p1 + p2 + p3
# dev.off()

saveRDS(media_integrated, "../../data/results/1.3.QC/media_integrated_noRegression.rds")

# Info on packages used
writeLines(capture.output(sessionInfo()), "../../data/results/1.3.QC/sessionInfo.txt")

# R version 4.1.0 (2021-05-18)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 18.04.5 LTS
# 
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.7.1
# LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.7.1
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] SCpubr_0.1.0       future_1.25.0      patchwork_1.1.1    sp_1.4-7           SeuratObject_4.1.0 Seurat_4.1.1       forcats_0.5.1      stringr_1.4.0      dplyr_1.0.9       
# [10] purrr_0.3.4        readr_2.1.2        tidyr_1.2.0        tibble_3.1.7       ggplot2_3.3.6      tidyverse_1.3.1   
# 
# loaded via a namespace (and not attached):
#   [1] Rtsne_0.16            colorspace_2.0-3      deldir_1.0-6          ellipsis_0.3.2        ggridges_0.5.3        fs_1.5.2              spatstat.data_2.2-0   rstudioapi_0.13      
# [9] farver_2.1.0          leiden_0.4.2          listenv_0.8.0         ggrepel_0.9.1         RSpectra_0.16-1       fansi_1.0.3           lubridate_1.8.0       xml2_1.3.3           
# [17] codetools_0.2-18      splines_4.1.0         polyclip_1.10-0       jsonlite_1.8.0        broom_0.8.0           ica_1.0-2             cluster_2.1.3         dbplyr_2.1.1         
# [25] png_0.1-7             rgeos_0.5-9           uwot_0.1.11           spatstat.sparse_2.1-1 sctransform_0.3.3     shiny_1.7.1           compiler_4.1.0        httr_1.4.3           
# [33] backports_1.4.1       assertthat_0.2.1      Matrix_1.4-1          fastmap_1.1.0         lazyeval_0.2.2        cli_3.3.0             later_1.3.0           htmltools_0.5.2      
# [41] tools_4.1.0           igraph_1.3.1          gtable_0.3.0          glue_1.6.2            reshape2_1.4.4        RANN_2.6.1            Rcpp_1.0.8.3          scattermore_0.8      
# [49] cellranger_1.1.0      vctrs_0.4.1           nlme_3.1-157          progressr_0.10.0      lmtest_0.9-40         spatstat.random_2.2-0 globals_0.15.0        rvest_1.0.2          
# [57] mime_0.12             miniUI_0.1.1.1        lifecycle_1.0.1       irlba_2.3.5           goftest_1.2-3         MASS_7.3-57           zoo_1.8-10            scales_1.2.0         
# [65] spatstat.core_2.4-2   spatstat.utils_2.3-1  hms_1.1.1             promises_1.2.0.1      parallel_4.1.0        RColorBrewer_1.1-3    gridExtra_2.3         reticulate_1.25      
# [73] pbapply_1.5-0         rpart_4.1.16          stringi_1.7.6         rlang_1.0.2           pkgconfig_2.0.3       matrixStats_0.62.0    lattice_0.20-45       tensor_1.5           
# [81] ROCR_1.0-11           labeling_0.4.2        htmlwidgets_1.5.4     cowplot_1.1.1         tidyselect_1.1.2      parallelly_1.31.1     RcppAnnoy_0.0.19      plyr_1.8.7           
# [89] magrittr_2.0.3        R6_2.5.1              generics_0.1.2        DBI_1.1.2             mgcv_1.8-40           pillar_1.7.0          haven_2.5.0           withr_2.5.0          
# [97] fitdistrplus_1.1-8    abind_1.4-5           survival_3.3-1        future.apply_1.9.0    modelr_0.1.8          crayon_1.5.1          KernSmooth_2.23-20    utf8_1.2.2           
# [105] spatstat.geom_2.4-0   plotly_4.10.0         tzdb_0.3.0            grid_4.1.0            readxl_1.4.0          data.table_1.14.2     reprex_2.0.1          digest_0.6.29        
# [113] xtable_1.8-4          httpuv_1.6.5          munsell_0.5.0         viridisLite_0.4.0    
