# label transfer with SingleR
# https://bioconductor.org/packages/release/bioc/html/SingleR.html
# http://bioconductor.org/books/devel/SingleRBook/introduction.html#quick-start

# To run interactively on server's Rstudio
# .libPaths("/software/teamtrynka/ma23/R4.1/libs") # path to libraries


library(Seurat)
library(SingleR)
library(scuttle)
library(scater)
library(future)
library(tidyverse)
library(patchwork)
library(scran)
library(SCpubr)
options(stringsAsFactors = FALSE)

# change the current plan to access parallelization
plan("multiprocess", workers = 14)
options(future.globals.maxSize= 2097152000) # 2Gb

plan()

### main -------

media_names = c("IGBN","IM","IMBN","ITGBN","ITM_ADMEM","ITMG")
dir.create(file.path("../../data/results/3.singleR_labelTransfer"))

# Query data: my  data
media_merged= readRDS("../../data/results/1.3.QC/media_merged_noRegression.rds")
media_merged@meta.data$cell = rownames(media_merged@meta.data)
media_merged = as.SingleCellExperiment(media_merged)
media_merged = logNormCounts(media_merged)

# Reference data: foetal microglia #######################
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE133345
foetal_microglia_Bian = read.table("../../data/resources/Foetal_microglia_Bian_2020/GSE133345_Quality_controled_UMI_data_of_all_1231_embryonic_cells.txt")

foetal_microglia_Bian = CreateSeuratObject(foetal_microglia_Bian)

foetal_microglia_Bian_metadata = read.table("../../data/resources/Foetal_microglia_Bian_2020/GSE133345_Annotations_of_all_1231_embryonic_cells_updated_0620.txt")
foetal_microglia_Bian_metadata$cluster = dplyr::recode(foetal_microglia_Bian_metadata$cluster,
                                                       "Mac_1" = "Mac. precursor 1", "Mac_2" = "Mac. precursor 2",
                "Mac_3" = "Mac. precursor 3", "Mac_4" ="Microglial precursor", 
                "Mast_cell" = "Mast cell")
  
foetal_microglia_Bian = AddMetaData(foetal_microglia_Bian,metadata=foetal_microglia_Bian_metadata)
foetal_microglia_Bian@meta.data$cell = rownames(foetal_microglia_Bian@meta.data)

foetal_microglia_Bian = as.SingleCellExperiment(foetal_microglia_Bian)

# SingleR() expects reference datasets to be normalized and log-transformed.
foetal_microglia_Bian = logNormCounts(foetal_microglia_Bian)
pred.foetal = SingleR(test=media_merged, ref=foetal_microglia_Bian, 
                      labels=foetal_microglia_Bian$cluster)
table(pred.foetal$labels)

png("../../data/results/3.singleR_labelTransfer/foetal_microglia_Bian_score_heatmap.png", 
    width = 12, height = 3, units = "in", res = 400)
plotScoreHeatmap(pred.foetal)
dev.off()
png("../../data/results/3.singleR_labelTransfer/foetal_microglia_Bian_delta_distribution.png", 
    width = 9, height = 9, units = "in", res = 400)
plotDeltaDistribution(pred.foetal, ncol = 3)
dev.off()
summary(is.na(pred.foetal$pruned.labels))
all.markers = metadata(pred.foetal)$de.genes
media_merged$foetal_labels = pred.foetal$labels

# top 20 markers per prediction
all.markers = metadata(pred.foetal)$de.genes
empirical.markers = findMarkers(foetal_microglia_Bian,
                                 foetal_microglia_Bian$cluster, 
                                 direction="up")

collected = list()
for (lab in unique(pred.foetal$labels)) {
  lab.markers = unique(unlist(all.markers[[lab]]))
  m = match(lab.markers, rownames(empirical.markers[[lab]]))
  m = lab.markers[rank(m) <= 20]
  collected[[lab]] = plotHeatmap(media_merged, silent=TRUE, 
                                  order_columns_by="foetal_labels", main=lab, features=m)[[4]]
}

png("../../data/results/3.singleR_labelTransfer/foetal_microglia_Bian_heatmap_markers.png", 
    width = 18, height = 18, units = "in", res = 400)
do.call(gridExtra::grid.arrange, collected)
dev.off()

## Give always same colors per category
group.colors = c(`Mac. precursor 1`= "#73937E", `Mac. precursor 2` = "#218380",
                 `Mac. precursor 3` = "#40916c", `Microglial precursor` ="#84C318", 
                 Monocyte = "#8f84d8", 
                  CD7loP = "#EF798A" ,   Myeloblast = "#CEB992", ErP = "#ff9a17",
                  YSMP = "#514B23", `Mast cell` = "#F72C25", MkP = "#CC7178")
long = colData(media_merged) %>%
  as_tibble() %>%
  group_by(orig.ident) %>%
  dplyr::count(foetal_labels)  %>% mutate(Percent = n / sum(n)*100)


barplot1 = ggplot(data=long, aes(x=orig.ident, y=Percent, fill=foetal_labels)) +
  geom_bar(stat="identity", col="black") + theme_bw() + 
  xlab("Library") + scale_fill_manual(values=group.colors) + 
  theme(axis.text.x = element_text(face = "bold"), axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"), 
        legend.text = element_text(face = "bold"), legend.title = element_blank(),
        plot.title = element_text(face = "bold")) +
ggtitle("Label transfer: foetal macrophages") + 
  guides(fill=guide_legend(ncol=2))





# Reference data: Mancuso ex vivo (in mice) matured iPSC derived microglia human #######
ex_vivo_mancuso=readRDS("../../data/resources/Mancuso_Nat_Neurosci_2019/ex_vivo_mancuso_filtered_annotated_allSamples.rds")
ex_vivo_mancuso@meta.data$renamed_clusters = "Homeostatic microglia"
ex_vivo_mancuso@meta.data[ex_vivo_mancuso@meta.data$seurat_clusters == 1,"renamed_clusters"] = "Perivascular macrophages"
ex_vivo_mancuso@meta.data[ex_vivo_mancuso@meta.data$seurat_clusters == 4,"renamed_clusters"] = "Perivascular macrophages"
ex_vivo_mancuso@meta.data[ex_vivo_mancuso@meta.data$seurat_clusters == 3,"renamed_clusters"] = "Cytokine resp. microglia"
ex_vivo_mancuso@meta.data[ex_vivo_mancuso@meta.data$seurat_clusters == 5,"renamed_clusters"] = "Cycling cells"
ex_vivo_mancuso@meta.data[ex_vivo_mancuso@meta.data$seurat_clusters == 6,"renamed_clusters"] = "Amyloid resp. microglia"


ex_vivo_mancuso = as.SingleCellExperiment(ex_vivo_mancuso)

# SingleR() expects reference datasets to be normalized and log-transformed.
ex_vivo_mancuso = logNormCounts(ex_vivo_mancuso)
ped.mancuso = SingleR(test=media_merged, ref=ex_vivo_mancuso, 
                      labels=ex_vivo_mancuso$renamed_clusters)
table(ped.mancuso$labels)

png("../../data/results/3.singleR_labelTransfer/Mancuso_xenografts_score_heatmap.png", 
    width = 12, height = 3, units = "in", res = 400)
plotScoreHeatmap(ped.mancuso)
dev.off()
png("../../data/results/3.singleR_labelTransfer/Mancuso_xenografts_delta_distribution.png", 
    width = 9, height = 9, units = "in", res = 400)
plotDeltaDistribution(ped.mancuso, ncol = 3)
dev.off()
summary(is.na(ped.mancuso$pruned.labels))
all.markers = metadata(ped.mancuso)$de.genes
media_merged$xenograft_labels = ped.mancuso$labels


## Give always same colors per category

group.colors = c(`Homeostatic microglia` = "#73937E",  `Perivascular macrophages` = "#218380", `Cytokine resp. microglia` = "#8f84d8", 
                   `Cycling cells` = "#EF798A", `Amyloid resp. microglia` ="#0FA3B1")


long = colData(media_merged) %>%
  as_tibble() %>%
  group_by(orig.ident) %>%
  dplyr::count(xenograft_labels)  %>% mutate(Percent = n / sum(n)*100)


barplot2 = ggplot(data=long, aes(x=orig.ident, y=Percent, fill=xenograft_labels)) +
  geom_bar(stat="identity", col="black") + theme_bw() + 
  xlab("Library") + scale_fill_manual(values=group.colors) + 
  theme(axis.text.x = element_text(face = "bold"), axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"), 
        legend.text = element_text(face = "bold"), legend.title = element_blank(),
        plot.title = element_text(face = "bold")) +
  ggtitle("Label transfer: xenografted iPSC-microglia")

png("../../data/results/3.singleR_labelTransfer/label_transfer_barplots_per_library_Fig5c.png", 
    width = 8, height = 5, units = "in", res = 400)
barplot1 / barplot2
dev.off()

#both - combining inferences

#http://bioconductor.org/books/release/SingleRBook/using-multiple-references.html#combining-inferences-from-individual-references

pred.both = SingleR(test=media_merged, 
                    ref=list(foetal_microglia_Bian, ex_vivo_mancuso), 
                      labels=list(foetal_microglia_Bian$cluster,ex_vivo_mancuso$renamed_clusters))
                 
# tried using de.method=wilcox, the default and wilcox + changing de.n to 500 * (2/3) ^ log2(N unique labels) 
# the only parameter that changes results is de.n, which for wilcox method the default is 10. Seems way too low and pointless 
# This can't be used to distinguish which reference is more similar to our data though, as explained in the vignette
# but can chow which references share similar gene expression, particularly if they are subsets of the same cell type
table(pred.both$labels)
table(pred.both$reference)

png("../../data/results/3.singleR_labelTransfer/both_score_heatmap.png", 
    width = 12, height = 3, units = "in", res = 400)
plotScoreHeatmap(pred.both, 
                 annotation_col=as.data.frame(colData(media_merged)[,"donor_id",drop=FALSE]))
dev.off()
png("../../data/results/3.singleR_labelTransfer/both_delta_distribution.png", 
    width = 9, height = 9, units = "in", res = 400)
plotDeltaDistribution(ped.mancuso, ncol = 3)
dev.off()


# saving metadata
head(colData(media_merged))
metadata = colData(media_merged)
write.table(colData(media_merged), file = "../../data/results/3.singleR_labelTransfer/singleR_metadata.txt", 
            col.names = T,
            row.names = F, 
            quote = F, sep = "\t")

# Dimplots
media_merged_seurat= readRDS("../../data/results/1.3.QC/media_merged_noRegression.rds")
if(identical(media_merged_seurat@meta.data$cell, as.data.frame(metadata)$cell)){
  media_merged_seurat@meta.data$foetal_labels = metadata$foetal_labels
  media_merged_seurat@meta.data$xenograft_labels = metadata$xenograft_labels
  }

p1 = SCpubr::do_DimPlot(sample = media_merged_seurat, 
                        plot.title = "Label transfer: foetal macrophages",
                        dims = c(1, 2), 
                        group.by = "foetal_labels",
                        colors.use =   c(`Mac. precursor 1`= "#73937E", `Mac. precursor 2` = "#218380",
                                         `Mac. precursor 3` = "#40916c", `Microglial precursor` ="#84C318", 
                                         Monocyte = "#8f84d8", 
                                         CD7loP = "#EF798A" ,   Myeloblast = "#CEB992", ErP = "#ff9a17",
                                         YSMP = "#514B23", `Mast cell` = "#F72C25", MkP = "#CC7178"), 
                        order = c("CD7loP","Microglial precursor","Monocyte","Mac. precursor 2","Mac. precursor 3","Mac. precursor 1"), shuffle = F)


p2 = SCpubr::do_DimPlot(sample = media_merged_seurat, 
                        plot.title = "Label transfer: xenografted iPSC-microglia",
                        dims = c(1, 2), 
                        group.by = "xenograft_labels",
                        colors.use =   c(`Homeostatic microglia` = "#73937E",  `Perivascular macrophages` = "#218380", `Cytokine resp. microglia` = "#8f84d8", 
                                         `Cycling cells` = "#EF798A", `Amyloid resp. microglia` ="#0FA3B1"), 
                        order = c("Cycling cells","Amyloid resp. microglia","Perivascular macrophages"  ), shuffle = F)



p1 | p2 

png("../../data/results/3.singleR_labelTransfer/dimplots_merged_singleR_FigS5D.png", 
    width = 13, height = 5, units = "in", res = 400)
p1 | p2 
dev.off()
# with integrated data
media_integrated_seurat= readRDS("../../data/results/1.3.QC/media_integrated_noRegression.rds")

if(identical(media_integrated_seurat@meta.data$cell, as.data.frame(metadata)$cell)){
  media_integrated_seurat@meta.data$foetal_labels = metadata$foetal_labels
  media_integrated_seurat@meta.data$xenograft_labels = metadata$xenograft_labels
}

p1 = SCpubr::do_DimPlot(sample = media_integrated_seurat, 
                              plot.title = "Label transfer: foetal macrophages",
                              dims = c(1, 2), 
                              group.by = "foetal_labels",
                              colors.use =   c(`Mac. precursor 1`= "#73937E", `Mac. precursor 2` = "#218380",
                                               `Mac. precursor 3` = "#40916c", `Microglial precursor` ="#84C318", 
                                               Monocyte = "#8f84d8", 
                                               CD7loP = "#EF798A" ,   Myeloblast = "#CEB992", ErP = "#ff9a17",
                                               YSMP = "#514B23", `Mast cell` = "#F72C25", MkP = "#CC7178"), 
                        shuffle = F)



p2 = SCpubr::do_DimPlot(sample = media_integrated_seurat, 
                        plot.title = "Label transfer: xenografted iPSC-microglia",
                        dims = c(1, 2), 
                        group.by = "xenograft_labels",
                        colors.use =   c(`Homeostatic microglia` = "#73937E",  `Perivascular macrophages` = "#218380", `Cytokine resp. microglia` = "#8f84d8", 
                                         `Cycling cells` = "#EF798A", `Amyloid resp. microglia` ="#0FA3B1"), 
                        shuffle = F)




# png("../../data/results/3.singleR_labelTransfer/dimplots_integrated_singleR.png", 
#     width = 13, height = 4, units = "in", res = 400)
p1 | p2 
# dev.off()


# Info on packages used
writeLines(capture.output(sessionInfo()), "../../data/results/3.singleR_labelTransfer/sessionInfo.txt")

# R version 4.0.2 (2020-06-22)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Catalina 10.15.7
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8
# 
# attached base packages:
#   [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] UpSetR_1.4.0                SCpubr_0.1.0                scran_1.18.7                patchwork_1.1.1             forcats_0.5.1              
# [6] stringr_1.4.0               dplyr_1.0.9                 purrr_0.3.4                 readr_2.0.1                 tidyr_1.2.0                
# [11] tibble_3.1.4                tidyverse_1.3.1             future_1.22.1               scater_1.18.6               ggplot2_3.3.6              
# [16] scuttle_1.0.4               SingleCellExperiment_1.12.0 SingleR_1.4.1               SummarizedExperiment_1.20.0 Biobase_2.50.0             
# [21] GenomicRanges_1.42.0        GenomeInfoDb_1.26.7         IRanges_2.24.1              S4Vectors_0.28.1            BiocGenerics_0.36.1        
# [26] MatrixGenerics_1.2.1        matrixStats_0.60.1          sp_1.4-5                    SeuratObject_4.1.0          Seurat_4.1.1               
# [31] devtools_2.4.2              usethis_2.0.1              
# 
# loaded via a namespace (and not attached):
#   [1] utf8_1.2.2                reticulate_1.20           tidyselect_1.1.1          htmlwidgets_1.5.4         grid_4.0.2                BiocParallel_1.24.1      
# [7] Rtsne_0.15                munsell_0.5.0             codetools_0.2-16          ica_1.0-2                 statmod_1.4.36            miniUI_0.1.1.1           
# [13] withr_2.4.2               colorspace_2.0-2          progressr_0.10.0          rstudioapi_0.13           ROCR_1.0-11               ggsignif_0.6.3           
# [19] tensor_1.5                listenv_0.8.0             labeling_0.4.2            GenomeInfoDbData_1.2.4    polyclip_1.10-0           farver_2.1.0             
# [25] pheatmap_1.0.12           rprojroot_2.0.2           parallelly_1.28.1         vctrs_0.4.1               generics_0.1.0            colortools_0.1.5         
# [31] R6_2.5.1                  ggbeeswarm_0.6.0          rsvd_1.0.5                locfit_1.5-9.4            bitops_1.0-7              spatstat.utils_2.2-0     
# [37] cachem_1.0.6              DelayedArray_0.16.3       assertthat_0.2.1          promises_1.2.0.1          scales_1.2.0              rgeos_0.5-9              
# [43] beeswarm_0.4.0            gtable_0.3.0              beachmat_2.6.4            globals_0.14.0            processx_3.5.2            goftest_1.2-2            
# [49] rlang_1.0.2               splines_4.0.2             rstatix_0.7.0             lazyeval_0.2.2            spatstat.geom_2.2-2       broom_0.7.9              
# [55] reshape2_1.4.4            abind_1.4-5               modelr_0.1.8              backports_1.2.1           httpuv_1.6.3              tools_4.0.2              
# [61] ellipsis_0.3.2            spatstat.core_2.3-0       RColorBrewer_1.1-2        sessioninfo_1.1.1         ggridges_0.5.3            Rcpp_1.0.7               
# [67] plyr_1.8.6                sparseMatrixStats_1.2.1   zlibbioc_1.36.0           RCurl_1.98-1.4            ps_1.6.0                  prettyunits_1.1.1        
# [73] ggpubr_0.4.0              rpart_4.1-15              deldir_0.2-10             pbapply_1.4-3             viridis_0.6.2             cowplot_1.1.1            
# [79] zoo_1.8-9                 haven_2.4.3               ggrepel_0.9.1             cluster_2.1.0             fs_1.5.0                  magrittr_2.0.1           
# [85] data.table_1.14.0         scattermore_0.7           openxlsx_4.2.4            reprex_2.0.1              lmtest_0.9-38             RANN_2.6.1               
# [91] fitdistrplus_1.1-5        pkgload_1.2.2             hms_1.1.0                 mime_0.11                 xtable_1.8-4              rio_0.5.27               
# [97] readxl_1.3.1              gridExtra_2.3             testthat_3.0.4            compiler_4.0.2            KernSmooth_2.23-17        crayon_1.4.1             
# [103] htmltools_0.5.2           tzdb_0.1.2                mgcv_1.8-31               later_1.3.0               lubridate_1.7.10          DBI_1.1.1                
# [109] dbplyr_2.1.1              MASS_7.3-51.6             car_3.0-11                Matrix_1.4-1              cli_3.3.0                 igraph_1.2.6             
# [115] pkgconfig_2.0.3           foreign_0.8-80            plotly_4.9.4.1            spatstat.sparse_2.0-0     xml2_1.3.2                vipor_0.4.5              
# [121] dqrng_0.3.0               XVector_0.30.0            rvest_1.0.1               callr_3.7.0               digest_0.6.27             sctransform_0.3.3        
# [127] RcppAnnoy_0.0.19          spatstat.data_2.1-0       cellranger_1.1.0          leiden_0.3.9              edgeR_3.32.1              uwot_0.1.10              
# [133] DelayedMatrixStats_1.12.3 curl_4.3.2                shiny_1.6.0               lifecycle_1.0.1           nlme_3.1-148              jsonlite_1.7.2           
# [139] carData_3.0-4             BiocNeighbors_1.8.2       limma_3.46.0              desc_1.3.0                viridisLite_0.4.0         fansi_0.5.0              
# [145] pillar_1.6.2              lattice_0.20-41           fastmap_1.1.0             httr_1.4.2                pkgbuild_1.2.0            survival_3.1-12          
# [151] glue_1.6.2                remotes_2.4.0             zip_2.2.0                 png_0.1-7                 bluster_1.0.0             stringi_1.7.4            
# [157] BiocSingular_1.6.0        memoise_2.0.0             irlba_2.3.3               future.apply_1.8.1       
