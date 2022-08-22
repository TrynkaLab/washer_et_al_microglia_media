# Check doublets and save the barcode info to remove them.
.libPaths("/software/teamtrynka/ma23/R4.1/libs") # path to libraries
library(Seurat)
library(tidyverse)
library(patchwork)

plan("multiprocess", workers = 14)
plan()
options(future.globals.maxSize= 2097152000) 
options(stringsAsFactors = FALSE)

dir.create(file.path("../../data/results/1.1.Determine_doublets"))

samples = c("IGBN","IM","IMBN","ITGBN","ITM_ADMEM","ITMG")

vireo = list()
donor_ids = list()
for (s in samples){
  donor_ids[[s]] = read.table(paste0("../../data/",s,"/Vireo_donor_deconvolution/vireoOutput/",
                                     s,"/donor_ids.tsv"), header = T)
  donor_ids[[s]]$category = s
  
  vireo[[s]] = read.table(paste0("../../data/",s,"/Vireo_donor_deconvolution/vireoOutput/",
                                 s,"/summary.tsv"), header = T)
  vireo[[s]]$category = s
}

donor_ids = do.call("rbind",donor_ids)
donor_ids$cell = gsub("\\-.*","",donor_ids$cell)
donor_ids$cell = paste(donor_ids$category,donor_ids$cell, sep = "_")

# save needed file to remove doublets in QC step
write.csv(donor_ids,"../../data/results/1.1.Determine_doublets/vireo_doublets_all.csv",row.names = F)

long = do.call("rbind",vireo)

colnames(long)[1] = "Donor"
long$Donor=factor(long$Donor,
                  levels = c("unassigned","doublet","aevs_1",     "aowh_2",     "coxy_33",   "eojr_2",
                             "hegp_3",     "hehd_1",  "iezw_2" ,   "letw_1",     "lizq_3",     "peoj_1",     "qolg_1",    
                             "sehp_2",     "tolg_4",     "ueah_4",      "uoxz_4",     "veku_2",     "woci_1",    
                             "xojn_3"   ,  "zaie_5"        ),ordered = T)

# Fill in missing categories

long = long %>% complete(Donor, nesting(category), fill = list(Freq = 0))
# Barplot with y axis in pseudo log scale (log after 0, linear at 0)
pseudolog_barplot_per_library = function(long_df){
  
  plot_list = list()
  for(library in unique(long$category)){
    plot_list[[library]]<- ggplot(long[long$category==library,], aes(x=Donor, y=Freq)) +
      geom_bar(stat="identity", color="black", position=position_dodge())+
      scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), breaks = c(0,1,10,50,100,250,500,1000,2500,5000,10000))+
      theme_bw() + 
      ylab("number of cells (log10 scale)") + 
      ggtitle(library) +
      theme(axis.text.x = element_text(face = "bold"), axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold"), 
            legend.text = element_text(face = "bold"), legend.title = element_blank(),
            plot.title = element_text(face = "bold"))
    
  }
  p = wrap_plots(plot_list,ncol = 1) + plot_layout(guides = "collect")
  return(p)
}


pseudolog_barplot_per_library(long)

cols = c("#0B00F6","#12130F","#0D67F6","#20F8F6","#F46B08","#F40006")

png("../../data/results/1.1.Determine_doublets/hipsci_ndonors_dots_fig_S5b.png", 
    width = 7, height = 4, units = "in", res = 400)
ggplot(long,aes(x=Donor,y=Freq,  color = category,fill=category)) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) + 
  scale_fill_manual(values = cols, 
                    breaks = samples,
                    labels = c("IGBN","IM","IMBN","ITGBN","ITM ADMEM", "ITMG"),
                    name = "Name") +
  scale_color_manual(values = cols, 
                     breaks = samples,
                     labels = c("IGBN","IM","IMBN","ITGBN","ITM ADMEM", "ITMG"),
                     name = "Name") +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), breaks = c(0,1,10,50,100,500,1000,2000, 3000, 5000))+
  theme_bw() + 
  theme(axis.text.x = element_text(colour = "black", angle = 90, vjust = 0.5, hjust =1),
        axis.title = element_text(face = "bold"), 
        legend.text = element_text(face = "bold"), legend.title = element_blank(),
        plot.title = element_text(face = "bold")) +
  ylab("Number of cells (log10 scale)") 

dev.off()

# Donor proportions

long = long %>% 
  group_by(category) %>%
  mutate(proportion = Freq / sum(Freq))

# write.table(long,"../../data/results/1.1.Determine_doublets/allmedia_donor_proportions.txt",
#             col.names = T, row.names = F, quote = F, sep = "\t")

# png("../../data/results/1.1.Determine_doublets/hipsci_donor_proportions_dots.png", 
#     width = 7, height = 4, units = "in", res = 400)
ggplot(long,aes(x=Donor,y=proportion, col=category, fill=category)) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) + 
  scale_fill_manual(values = cols, 
                    breaks = samples,
                    labels = c("IGBN","IM","IMBN","ITGBN","ITM ADMEM", "ITMG"),
                    name = "Name") +
  scale_color_manual(values = cols, 
                     breaks = samples,
                     labels = c("IGBN","IM","IMBN","ITGBN","ITM ADMEM", "ITMG"),
                     name = "Name") +
  theme_bw() + 
  theme(axis.text.x = element_text(colour = "black", angle = 90, vjust = 0.5, hjust =1)) +
  ylab("donor proportions") 

# dev.off()

writeLines(capture.output(sessionInfo()), "../../data/results/1.1.Determine_doublets/sessionInfo.txt")

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
