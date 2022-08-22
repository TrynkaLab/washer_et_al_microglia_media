# Determine optimal clustering parameters for media integration seurat object
.libPaths("/software/teamtrynka/ma23/R4.1/libs")

library(clustree)
library(Seurat)
library(future)
options(future.globals.maxSize= 3097152000) # 3Gb for seurat
plan("multiprocess", workers = 14) # more workers
# check the current active plan
plan()

dir.create("../../data/results/1.2.determine_cluster_parameters", recursive = T)
###### functions --------
filter_seurat = function(seurat_object){
  
  message("Performing filter by number of genes and mitochondrial percentage.")
  seurat_object <- subset(seurat_object, subset = nFeature_RNA > 1000  & percent.mt < 10 & nCount_RNA < 40000 & vireo_doublets == "singlet")
  message("Now the object has ", dim(seurat_object)[1], " genes and ", dim(seurat_object)[2], " cells.")
  return(seurat_object)
}

calculate_PCA_UMAP_neighbors_clusters_integrated = function(seurat_object){
  seurat_object <- RunPCA(seurat_object, verbose = FALSE)
  
  seurat_object <- FindNeighbors(seurat_object, dims = 1:30, verbose = FALSE)
  seurat_object <- FindClusters(seurat_object, 
                                resolution = 0.3, # Selected in 1.1.Determine_clustering_params
                                verbose = FALSE)
  
  seurat_object <- RunUMAP(seurat_object, dims = 1:30, verbose = FALSE)
  
  
  
  return(seurat_object)
}

### Main ---------
# Load un-filtered data
media_names = c("IGBN","IM","IMBN","ITGBN","ITM_ADMEM","ITMG")

cellranger_data = list()
for(media in media_names){
  cellranger_data[[media]] =  Read10X(data.dir = paste0("../../data/",media,"/outs/filtered_feature_bc_matrix/"),
                                      strip.suffix = T)  ## remove trailing -1 from cell barcodes
  
}
str(cellranger_data)

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

copy_list = seurat_list
for(media in media_names){
  copy_list[[media]] = NormalizeData(seurat_list[[media]],verbose = T)
  copy_list[[media]] = CellCycleScoring(copy_list[[media]],  s.features=cc.genes.updated.2019$s.genes, g2m.features=cc.genes.updated.2019$g2m.genes)
  seurat_list[[media]]$Phase = copy_list[[media]]$Phase
  seurat_list[[media]]$S.Score = copy_list[[media]]$S.Score
  seurat_list[[media]]$G2M.Score = copy_list[[media]]$G2M.Score
}
# doublet
vireo_doublets = read.csv("../../data/results/1.1.Determine_doublets/vireo_doublets_all.csv")

for(media in media_names){
  seurat_list[[media]]@meta.data$original_cell = rownames(seurat_list[[media]]@meta.data)
  seurat_list[[media]]@meta.data = merge(seurat_list[[media]]@meta.data,vireo_doublets, by = "cell" )
  seurat_list[[media]]@meta.data$vireo_doublets = "singlet"
  seurat_list[[media]]@meta.data[seurat_list[[media]]@meta.data$donor_id == "doublet","vireo_doublets"] = "doublet"
  seurat_list[[media]]@meta.data[seurat_list[[media]]@meta.data$donor_id == "unassigned","vireo_doublets"] = "unassigned"
  rownames(seurat_list[[media]]@meta.data) = seurat_list[[media]]@meta.data$original_cell
}

#Filter

for(media in media_names){
  seurat_list[[media]] = filter_seurat(seurat_list[[media]])
}

# Normalise and regress
for(media in media_names){
  seurat_list[[media]] <- SCTransform(seurat_list[[media]], verbose = TRUE)
}


media_features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 3000)
seurat_list <- PrepSCTIntegration(object.list = seurat_list, anchor.features = media_features, 
                                  verbose = TRUE)

plan("multiprocess", workers = 1) # To avoid serialize error

media_anchors <- FindIntegrationAnchors(object.list = seurat_list, normalization.method = "SCT", 
                                        anchor.features = media_features, verbose = TRUE)
gc()
media_integrated <- IntegrateData(anchorset = media_anchors, normalization.method = "SCT", 
                                  verbose = TRUE)

# PCA
media_integrated = calculate_PCA_UMAP_neighbors_clusters_integrated(media_integrated)

## Clustree

media_integrated <- FindClusters(
  object = media_integrated ,
  reduction.type = "pca",
  resolution = c(0.1, 0.2, 0.3, 0.4, 0.6,0.8),
  dims.use = 1:30,
  save.SNN = TRUE
)


# Build cluster trees
clustree(media_integrated)

# png("../../data/results/1.2.determine_cluster_parameters/media_integrated_cellCycle_mitPercent_regressed_cluster_resolution_tree.png", 
#     width = 6, height = 10, units = "in", res = 400)
clustree(media_integrated)
# dev.off()

# png("../../data/results/1.2.determine_cluster_parameters/media_integrated_cellCycle_mitPercent_regressed_cluster_resolution_overlay_pca.png", 
#     width = 6, height = 4, units = "in", res = 400)
clustree_overlay(media_integrated, "pca1", "pca2", red_dim = "pca")

# dev.off()

# png("../../data/results/1.2.determine_cluster_parameters/media_integrated_cellCycle_mitPercent_regressed_cluster_resolution_overlay_umap.png", 
#     width = 6, height = 4, units = "in", res = 400)
clustree_overlay(media_integrated, "umap1", "umap2", red_dim = "umap")

# dev.off()


# Info on packages used
writeLines(capture.output(sessionInfo()), "../../data/results/1.2.determine_cluster_parameters/sessionInfo.txt")

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
# [10] purrr_0.3.4        readr_2.1.2        tidyr_1.2.0        tibble_3.1.7       ggplot2_3.3.6      tidyverse_1.3.1   clustree_0.4.3  
# 
# loaded via a namespace (and not attached):
#   [1] Rtsne_0.16            colorspace_2.0-3      deldir_1.0-6          ellipsis_0.3.2        ggridges_0.5.3        fs_1.5.2              spatstat.data_2.2-0   rstudioapi_0.13      
# [9] farver_2.1.0          leiden_0.4.2          listenv_0.8.0         ggrepel_0.9.1         fansi_1.0.3           lubridate_1.8.0       xml2_1.3.3            codetools_0.2-18     
# [17] splines_4.1.0         polyclip_1.10-0       jsonlite_1.8.0        broom_0.8.0           ica_1.0-2             cluster_2.1.3         dbplyr_2.1.1          png_0.1-7            
# [25] rgeos_0.5-9           uwot_0.1.11           spatstat.sparse_2.1-1 sctransform_0.3.3     shiny_1.7.1           compiler_4.1.0        httr_1.4.3            backports_1.4.1      
# [33] assertthat_0.2.1      Matrix_1.4-1          fastmap_1.1.0         lazyeval_0.2.2        cli_3.3.0             later_1.3.0           htmltools_0.5.2       tools_4.1.0          
# [41] igraph_1.3.1          gtable_0.3.0          glue_1.6.2            reshape2_1.4.4        RANN_2.6.1            Rcpp_1.0.8.3          scattermore_0.8       cellranger_1.1.0     
# [49] vctrs_0.4.1           nlme_3.1-157          progressr_0.10.0      lmtest_0.9-40         spatstat.random_2.2-0 globals_0.15.0        rvest_1.0.2           mime_0.12            
# [57] miniUI_0.1.1.1        lifecycle_1.0.1       irlba_2.3.5           goftest_1.2-3         MASS_7.3-57           zoo_1.8-10            scales_1.2.0          spatstat.core_2.4-2  
# [65] spatstat.utils_2.3-1  hms_1.1.1             promises_1.2.0.1      parallel_4.1.0        RColorBrewer_1.1-3    gridExtra_2.3         reticulate_1.25       pbapply_1.5-0        
# [73] rpart_4.1.16          stringi_1.7.6         rlang_1.0.2           pkgconfig_2.0.3       matrixStats_0.62.0    lattice_0.20-45       tensor_1.5            ROCR_1.0-11          
# [81] labeling_0.4.2        htmlwidgets_1.5.4     cowplot_1.1.1         tidyselect_1.1.2      parallelly_1.31.1     RcppAnnoy_0.0.19      plyr_1.8.7            magrittr_2.0.3       
# [89] R6_2.5.1              generics_0.1.2        DBI_1.1.2             mgcv_1.8-40           pillar_1.7.0          haven_2.5.0           withr_2.5.0           fitdistrplus_1.1-8   
# [97] abind_1.4-5           survival_3.3-1        future.apply_1.9.0    modelr_0.1.8          crayon_1.5.1          KernSmooth_2.23-20    utf8_1.2.2            spatstat.geom_2.4-0  
# [105] plotly_4.10.0         tzdb_0.3.0            grid_4.1.0            readxl_1.4.0          data.table_1.14.2     reprex_2.0.1          digest_0.6.29         xtable_1.8-4         
# [113] httpuv_1.6.5          munsell_0.5.0         viridisLite_0.4.0    
