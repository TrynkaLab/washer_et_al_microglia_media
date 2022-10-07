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

png("../../data/results/1.1.Determine_doublets/hipsci_ndonors_dots_fig_S6b.png", 
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
