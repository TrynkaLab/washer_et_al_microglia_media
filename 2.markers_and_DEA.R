# Plottingmarker gene lists, differential expression analysis
# this needs to be run locally because nebulosa won't install on linux on server
# .libPaths("/software/teamtrynka/ma23/R4.1/libs") # path to libraries

library(tidyverse)
library(Seurat)
library(patchwork)
library(SCpubr)
library(future)
library(UpSetR)
# change the current plan to access parallelization
plan("multiprocess", workers = 8)
plan()
options(future.globals.maxSize= 4097152000) # 2Gb
options(stringsAsFactors = FALSE)


dir.create(file.path("../../data/results/2.markers_and_DEA"))

media_merged= readRDS("../../data/results/1.3.QC/media_merged_noRegression.rds")

# removing underscores for plots
media_merged@meta.data$orig.ident = media_merged@meta.data$orig.ident %>%
  dplyr::recode("ITM_ADMEM" = "ITM ADMEM")

media_names = c("IGBN","IM","IMBN","ITGBN","ITM ADMEM","ITMG")
cols = c("#0B00F6","#12130F","#3C71C6","#20F8F6","#AA00FF","#F40006")


DefaultAssay(media_merged) = "SCT"

# Selected markers
micro = c( "CX3CR1",  
           "MERTK", 
           "OLFML3",
           "CXCL8", # IL8
           "TNF", # TNFa
           "TREM2", 
           "P2RY12", 
           "TMEM119", 
           "TGFB1",
           "GPR34",
           "IL1B",
           "CSF1R",
           "GAS6",
           "PROS1",
           "C1QA",
           "CD200R1")

peri= c("CD163",
        "COLEC12", 
        "LYVE1", 
        "F13A1")


p1 = SCpubr::do_FeaturePlot(sample = media_merged, 
                            features = micro) 
p2 = SCpubr::do_NebulosaPlot(sample = media_merged, 
                             features = micro) + 
  patchwork::plot_annotation(title = "Microglia markers") +
  patchwork::plot_layout(guides = 'collect') 



png("../../data/results/2.markers_and_DEA/nebulosa_plots_microglia_markers_FigS5b.png", 
    width = 11, height = 10, units = "in", res = 400)
p2
dev.off()



p1 = SCpubr::do_FeaturePlot(sample = media_merged, 
                            features = peri ) 
p2 = SCpubr::do_NebulosaPlot(sample = media_merged, 
                             features = peri) + 
  patchwork::plot_annotation(title = "Perivascular macrophage markers") +
  patchwork::plot_layout(guides = 'collect') 



png("../../data/results/2.markers_and_DEA/nebulosa_plots_pvM_markers_FigS5b.png", 
    width = 6, height = 5, units = "in", res = 400)
p2
dev.off()


# Differential expression analysis - against IGBN baseline

diff_expr_overall= function(seurat_object, media_names){
  DefaultAssay(seurat_object) = "RNA"
  #  DE needs the normalised counts from the RNA assay (all genes) without any correction that might introduce biases
  # Though it's unclear from Seurat forum and vignettes
  # https://github.com/satijalab/seurat/issues/3933
  # https://github.com/satijalab/seurat/issues/2180
  # https://satijalab.org/seurat/articles/integration_introduction.html
  # https://github.com/satijalab/seurat/issues/2067
  # Then you'll have to see which effects are caused by batch/technical issues by other methods (i.e. identifying and removing clusters of cycling cells)
  seurat_object = NormalizeData(object = seurat_object, assay = "RNA")
  combinations = matrix(c(rep(media_names[1],4),media_names), ncol = 5, byrow = T)
  Idents(seurat_object) = "orig.ident"
  
  diff_results = list()
  
  for(comb_col in 1:ncol(combinations)){
    message(paste0("Working on combination number ",comb_col))
    diff_results[[comb_col]] = FindMarkers(seurat_object, 
                                            ident.1 = combinations[2,comb_col], 
                                            ident.2 = combinations[1,comb_col],
                                            test.use = "wilcox",
                                            verbose = TRUE)
    
    # ident.2 is the baseline
    diff_results[[comb_col]]$gene = rownames(diff_results[[comb_col]])
    diff_results[[comb_col]]$contrast = paste0(combinations[2,comb_col], "/", combinations[1,comb_col]) 
  }
  diff_results = do.call("rbind",diff_results)
  return(diff_results)
}

diff_results_all = diff_expr_overall(media_merged, media_names)
diff_results_all_sign = diff_results_all[diff_results_all$p_val_adj<0.05,]

table(diff_results_all_sign$contrast)
mean(table(diff_results_all_sign$contrast))
genes = unique(diff_results_all_sign$gene)

write.csv(diff_results_all_sign, 
            file = "../../data/results/2.markers_and_DEA/DEA_table.csv",
            col.names = T, row.names = F, quote = F)

# prepare data for upsetr

upsetr_combined = purrr::reduce(list(data.frame(gene=diff_results_all_sign[diff_results_all_sign$contrast == "IM/IGBN", "gene"], "IM/IGBN"=1),
                        data.frame(gene=diff_results_all_sign[diff_results_all_sign$contrast == "IMBN/IGBN", "gene"], "IMBN/IGBN"=1),
                        data.frame(gene=diff_results_all_sign[diff_results_all_sign$contrast == "ITGBN/IGBN", "gene"], "ITGBN/IGBN"=1),
                                   data.frame(gene=diff_results_all_sign[diff_results_all_sign$contrast == "ITM ADMEM/IGBN", "gene"], "ITM ADMEM/IGBN" =1),
                                              data.frame(gene=diff_results_all_sign[diff_results_all_sign$contrast == "ITMG/IGBN", "gene"], "ITMG/IGBN"=1)),
                        full_join)
upsetr_combined[is.na(upsetr_combined)] = 0
names(upsetr_combined) = c("gene","IM/IGBN","IMBN/IGBN","ITGBN/IGBN","ITM ADMEM/IGBN", "ITMG/IGBN")

png("../../data/results/2.markers_and_DEA/upsetplot_diff_expr_FigS5c.png", 
    width = 6, height = 6, units = "in", res = 400)
upset(upsetr_combined, order.by = "freq", sets.bar.color = cols[c(5,6,2,3,4)])
dev.off()

# plot differentially expressed genes

DefaultAssay(media_merged) = "SCT"


names(cols) = media_names
p1 = Seurat::VlnPlot(media_merged, 
                     features = intersect(micro, genes), cols = cols,
                     split.by = "orig.ident", group.by = "orig.ident", 
                     pt.size = 0, combine = FALSE, stack = T, flip = T) + 
  ylab("SCT-normalised expression") + 
  theme(legend.position = "none", axis.title.x = element_blank()) +
  ggtitle("Microglia markers")

png("../../data/results/2.markers_and_DEA/violin_plots_microglia_markers_Fig5b.png", 
    width = 7, height = 10, units = "in", res = 400)
p1
dev.off()



p2 = Seurat::VlnPlot(media_merged, 
                     features = intersect(peri, genes), cols = cols,
                     split.by = "orig.ident", group.by = "orig.ident", 
                     pt.size = 0, combine = FALSE, stack = T, flip = T) + 
  ylab("SCT-normalised expression") + 
  theme(legend.position = "none", axis.title.x = element_blank()) +
  ggtitle("Perivascular macrophages markers")

png("../../data/results/2.markers_and_DEA/violin_plots_pvM_markers_Fig5b.png", 
    width = 6, height = 6, units = "in", res = 400)
p2
dev.off()

# Info on packages used
writeLines(capture.output(sessionInfo()), "../../data/results/2.markers_and_DEA/sessionInfo.txt")

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
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] UpSetR_1.4.0       SCpubr_0.1.0       future_1.22.1      patchwork_1.1.1    sp_1.4-5           SeuratObject_4.1.0 Seurat_4.1.1       forcats_0.5.1      stringr_1.4.0     
# [10] dplyr_1.0.9        purrr_0.3.4        readr_2.0.1        tidyr_1.2.0        tibble_3.1.4       ggplot2_3.3.6      tidyverse_1.3.1    devtools_2.4.2     usethis_2.0.1     
# 
# loaded via a namespace (and not attached):
#   [1] utf8_1.2.2                  ks_1.13.5                   reticulate_1.20             tidyselect_1.1.1            htmlwidgets_1.5.4           grid_4.0.2                 
# [7] Rtsne_0.15                  munsell_0.5.0               codetools_0.2-16            ica_1.0-2                   miniUI_0.1.1.1              withr_2.4.2                
# [13] colorspace_2.0-2            progressr_0.10.0            Biobase_2.50.0              rstudioapi_0.13             stats4_4.0.2                SingleCellExperiment_1.12.0
# [19] ROCR_1.0-11                 ggsignif_0.6.3              tensor_1.5                  listenv_0.8.0               MatrixGenerics_1.2.1        labeling_0.4.2             
# [25] GenomeInfoDbData_1.2.4      polyclip_1.10-0             farver_2.1.0                Nebulosa_1.0.2              rprojroot_2.0.2             parallelly_1.28.1          
# [31] vctrs_0.4.1                 generics_0.1.0              colortools_0.1.5            R6_2.5.1                    GenomeInfoDb_1.26.7         DelayedArray_0.16.3        
# [37] bitops_1.0-7                spatstat.utils_2.2-0        cachem_1.0.6                assertthat_0.2.1            promises_1.2.0.1            scales_1.2.0               
# [43] rgeos_0.5-9                 gtable_0.3.0                globals_0.14.0              processx_3.5.2              goftest_1.2-2               rlang_1.0.2                
# [49] splines_4.0.2               rstatix_0.7.0               lazyeval_0.2.2              spatstat.geom_2.2-2         broom_0.7.9                 reshape2_1.4.4             
# [55] abind_1.4-5                 modelr_0.1.8                backports_1.2.1             httpuv_1.6.3                tools_4.0.2                 ellipsis_0.3.2             
# [61] spatstat.core_2.3-0         RColorBrewer_1.1-2          BiocGenerics_0.36.1         sessioninfo_1.1.1           ggridges_0.5.3              Rcpp_1.0.7                 
# [67] plyr_1.8.6                  zlibbioc_1.36.0             RCurl_1.98-1.4              ps_1.6.0                    prettyunits_1.1.1           ggpubr_0.4.0               
# [73] rpart_4.1-15                deldir_0.2-10               pbapply_1.4-3               viridis_0.6.2               cowplot_1.1.1               S4Vectors_0.28.1           
# [79] zoo_1.8-9                   SummarizedExperiment_1.20.0 haven_2.4.3                 ggrepel_0.9.1               cluster_2.1.0               fs_1.5.0                   
# [85] magrittr_2.0.1              data.table_1.14.0           scattermore_0.7             openxlsx_4.2.4              lmtest_0.9-38               reprex_2.0.1               
# [91] RANN_2.6.1                  mvtnorm_1.1-3               fitdistrplus_1.1-5          matrixStats_0.60.1          pkgload_1.2.2               hms_1.1.0                  
# [97] mime_0.11                   xtable_1.8-4                rio_0.5.27                  mclust_5.4.7                readxl_1.3.1                IRanges_2.24.1             
# [103] gridExtra_2.3               testthat_3.0.4              compiler_4.0.2              KernSmooth_2.23-17          crayon_1.4.1                htmltools_0.5.2            
# [109] mgcv_1.8-31                 later_1.3.0                 tzdb_0.1.2                  lubridate_1.7.10            DBI_1.1.1                   dbplyr_2.1.1               
# [115] MASS_7.3-51.6               Matrix_1.4-1                car_3.0-11                  cli_3.3.0                   parallel_4.0.2              igraph_1.2.6               
# [121] GenomicRanges_1.42.0        pkgconfig_2.0.3             foreign_0.8-80              plotly_4.9.4.1              spatstat.sparse_2.0-0       xml2_1.3.2                 
# [127] XVector_0.30.0              rvest_1.0.1                 callr_3.7.0                 digest_0.6.27               pracma_2.3.3                sctransform_0.3.3          
# [133] RcppAnnoy_0.0.19            spatstat.data_2.1-0         cellranger_1.1.0            leiden_0.3.9                uwot_0.1.10                 curl_4.3.2                 
# [139] shiny_1.6.0                 lifecycle_1.0.1             nlme_3.1-148                jsonlite_1.7.2              carData_3.0-4               limma_3.46.0               
# [145] desc_1.3.0                  viridisLite_0.4.0           fansi_0.5.0                 pillar_1.6.2                lattice_0.20-41             fastmap_1.1.0              
# [151] httr_1.4.2                  pkgbuild_1.2.0              survival_3.1-12             glue_1.6.2                  remotes_2.4.0               zip_2.2.0                  
# [157] png_0.1-7                   stringi_1.7.4               memoise_2.0.0               irlba_2.3.3                 future.apply_1.8.1         
