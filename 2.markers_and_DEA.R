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

media_names = c("IGBN","IM","IMBN","ITGBN","ITM_RPMI","ITMG")

dir.create(file.path("../../data/results/2.markers_and_DEA"))

media_integrated =  readRDS(paste0("../../data/results/1.3.QC/media_integrated_noRegression.rds"))

DefaultAssay(media_integrated) = "SCT"

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


p1 = SCpubr::do_FeaturePlot(sample = media_integrated, 
                            features = micro) 
p2 = SCpubr::do_NebulosaPlot(sample = media_integrated, 
                             features = micro) + 
  patchwork::plot_annotation(title = "Microglia markers") +
  patchwork::plot_layout(guides = 'collect') 



png("../../data/results/2.markers_and_DEA/nebulosa_plots_microglia_markers_FigS6c.png", 
    width = 11, height = 10, units = "in", res = 400)
p2
dev.off()



p1 = SCpubr::do_FeaturePlot(sample = media_integrated, 
                            features = peri ) 
p2 = SCpubr::do_NebulosaPlot(sample = media_integrated, 
                             features = peri) + 
  patchwork::plot_annotation(title = "Perivascular macrophage markers") +
  patchwork::plot_layout(guides = 'collect') 



png("../../data/results/2.markers_and_DEA/nebulosa_plots_pvM_markers_FigS6c.png", 
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

diff_results_all = diff_expr_overall(media_integrated, media_names)
diff_results_all_sign = diff_results_all[diff_results_all$p_val_adj<0.05,]

table(diff_results_all_sign$contrast)
mean(table(diff_results_all_sign$contrast))
genes = unique(diff_results_all_sign$gene)

write.table(diff_results_all_sign, 
            file = "../../data/results/2.markers_and_DEA/DEA_table.txt",
            col.names = T, row.names = F, quote = F)

# prepare data for upsetr

upsetr_combined = purrr::reduce(list(data.frame(gene=diff_results_all_sign[diff_results_all_sign$contrast == "IM/IGBN", "gene"], "IM/IGBN"=1),
                        data.frame(gene=diff_results_all_sign[diff_results_all_sign$contrast == "IMBN/IGBN", "gene"], "IMBN/IGBN"=1),
                        data.frame(gene=diff_results_all_sign[diff_results_all_sign$contrast == "ITGBN/IGBN", "gene"], "ITGBN/IGBN"=1),
                                   data.frame(gene=diff_results_all_sign[diff_results_all_sign$contrast == "ITM_RPMI/IGBN", "gene"], "ITM_RPMI/IGBN" =1),
                                              data.frame(gene=diff_results_all_sign[diff_results_all_sign$contrast == "ITMG/IGBN", "gene"], "ITMG/IGBN"=1)),
                        full_join)
upsetr_combined[is.na(upsetr_combined)] = 0
names(upsetr_combined) = c("gene","IM/IGBN","IMBN/IGBN","ITGBN/IGBN","ITM_RPMI/IGBN", "ITMG/IGBN")

png("../../data/results/2.markers_and_DEA/upsetplot_diff_expr_FigS6d.png", 
    width = 6, height = 6, units = "in", res = 400)
upset(upsetr_combined, order.by = "freq", sets.bar.color = cols[c(4,3,2,6,5)])
dev.off()

# plot differentially expressed genes

DefaultAssay(media_integrated) = "SCT"


cols = c("#0B00F6","#12130F","#0D67F6","#20F8F6","#F46B08","#F40006")
names(cols) = media_names
p1 = Seurat::VlnPlot(media_integrated, 
                     features = intersect(micro, genes), cols = cols,
                     split.by = "orig.ident", group.by = "orig.ident", 
                     pt.size = 0, combine = FALSE, stack = T, flip = T) + 
  ylab("SCT-normalised expression") + 
  ggtitle("Microglia markers")

png("../../data/results/2.markers_and_DEA/violin_plots_microglia_markers_FigS6e.png", 
    width = 7, height = 10, units = "in", res = 400)
p1
dev.off()



p2 = Seurat::VlnPlot(media_integrated, 
                     features = intersect(peri, genes), cols = cols,
                     split.by = "orig.ident", group.by = "orig.ident", 
                     pt.size = 0, combine = FALSE, stack = T, flip = T) + 
  ylab("SCT-normalised expression") + 
  ggtitle("Perivascular macrophages markers")

png("../../data/results/2.markers_and_DEA/violin_plots_pvM_markers_FigS6e.png", 
    width = 6, height = 6, units = "in", res = 400)
p2
dev.off()

# Info on packages used
writeLines(capture.output(sessionInfo()), "../../data/results/2.markers_and_DEA/sessionInfo.txt")
