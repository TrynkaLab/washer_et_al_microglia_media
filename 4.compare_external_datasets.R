# compare our data to external datasets
# data not available on GEO was downloaded from stemformatics
library(Seurat)
library(future)
library(tidyverse)
library(SingleCellExperiment)
library(scater)
library(biomaRt)
library(limma)
library(edgeR)
library(ggrepel)
options(stringsAsFactors = FALSE)

dir.create(file.path("../../data/results/4.Compare_external_datasets"))

## functions ### ######
#to ensembl gene ids
geneName_to_ensembl = function(expr_df){
  
  
  ensembl38 = useMart("ensembl",dataset="hsapiens_gene_ensembl")
  
  
  gene_names =  getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                      filters = 'external_gene_name', values = rownames(expr_df), 
                      mart = ensembl38, useCache = FALSE)
  rownames(expr_df) = gene_names$ensembl_gene_id[match(rownames(expr_df),gene_names$external_gene_name)] 
  
  colSums(expr_df) # "library size" of log2norm counts
  expr_df = as.data.frame(expr_df)
  
  # remove NAs
  keep = !grepl("NA",rownames(expr_df))
  
  expr_df = expr_df[keep,]
  return(expr_df)
}
####main ######
# https://www.stemformatics.org/datasets/view?id=7268
Abud = read.table("../../../resources/Abud_2017_microglia_bulk/stemformatics_dataset_7268.raw.tsv")
#metadata from https://www.stemformatics.org/datasets/view?id=7268
colnames(Abud) = c("Abud_iMG_1","Abud_iMG_2","Abud_iMG_3","Abud_iMG_4","Abud_iMG_5","Abud_iMG_6",
                   "Abud_iHPC_1","Abud_iHPC_2","Abud_iHPC_3", # hematopoietic precursor
                   "Abud_iPSC_1","Abud_iPSC_2","Abud_iPSC_3","Abud_iPSC_4",
                   "Abud_monocyte_CD14_1","Abud_monocyte_CD14_2","Abud_monocyte_CD14_3","Abud_monocyte_CD14_4","Abud_monocyte_CD14_5",
                   "Abud_monocyte_CD14_CD16_1","Abud_monocyte_CD14_CD16_2","Abud_monocyte_CD14_CD16_3","Abud_monocyte_CD14_CD16_4",
                   "Abud_pFMGL_1","Abud_pFMGL_2","Abud_pFMGL_3", # foetal primary microglia
                   "Abud_pAMGL_1","Abud_pAMGL_2","Abud_pAMGL_3", # adult primary microglia
                   "Abud_DC_1","Abud_DC_2","Abud_DC_3", # dendritic cells
                   "Abud_iMG_neuron_1", "Abud_iMG_neuron_2","Abud_iMG_neuron_3","Abud_iMG_neuron_4","Abud_iMG_neuron_5", "Abud_iMG_neuron_6", # iPSC-derived microglia cocultured with neuron
                   "Abud_iMG_noCD200_noCX3CL1_1", "Abud_iMG_noCD200_noCX3CL1_2", "Abud_iMG_noCD200_noCX3CL1_3",  # microglial cell differentiated from induced pluripotent stem cell with CD200 and CX3CL1 withdrawl
                   "Abud_iMG_noTGFb_1","Abud_iMG_noTGFb_2","Abud_iMG_noTGFb_3") # microglial cell differentiated from induced pluripotent stem cell with TGFβ withdrawal

# ex vivo FACS-isolated microglia from fresh postmortem (Galatro_2017) or surgery-resected human brain (Gosselin_2017)
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE99074
Galatro = read.table("../../../resources/Galatro_2017/GSE99074_HumanMicrogliaBrainCounts.txt") # Grubman says these are raw counts but it's not true (there are decimal points)?

colnames(Galatro) = c( paste0("Galatro_exMGL_",1:39),
                       paste0("Galatro_exBrain_",1:16),
                       "S969" ,     "S947"  ,    "S861"   ,   "S810"  ,    "S805"  ,    "S755"    ,  "S291"  ,    "S262"    ,  "S243"    ,  "S147"   )

# Gosselin data is mouse, so not included

# Douvaras 2017
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97744 for metadata
# data from https://www.stemformatics.org/datasets/view?id=7193
Douvaras = read.table("../../../resources/Douvaras_2017/stemformatics_dataset_7193.raw.tsv")

colnames(Douvaras) = c(paste0("Douvaras_pMAC_",1:8), # primary macrophages
                       paste0("Douvaras_pAMGL_serum_",1:2), # primary microglia + serum
                       paste0("Douvaras_iMG_",1:10), # iPSC-derived microglia
                       paste0("Douvaras_pAMGL_",1:2), # primary microglia no serum
                       paste0("Douvaras_PMAC_hepatic_",1:2) # primary hepatic macrophages 
)


# Muffat 2016
# https://www.stemformatics.org/datasets/view?id=7062
Muffat = read.table("../../../resources/Muffat_2016/stemformatics_dataset_7062.raw.tsv")
colnames(Muffat) = c("Muffat_pFMGL_1","Muffat_pFMGL_2","Muffat_pFMGL_3",
                     paste0("Muffat_iMG_",1:5),
                     "Muffat_eMG",
                     paste0("Muffat_iMG_neuron_",1:3), # iPSC-microglia cultured in neural conditioned medium
                     paste0("Muffat_iNPC_",1:4) # iPSC-derived neuronal progenitor
)


# Mancuso data ex vivo
ex_vivo_mancuso=readRDS("../../data/resources/Mancuso_Nat_Neurosci_2019/ex_vivo_mancuso_filtered_annotated_allSamples.rds")
ex_vivo_mancuso@meta.data$renamed_clusters = "Mancuso_exMGL_1" # Homeostatic microglia
ex_vivo_mancuso@meta.data[ex_vivo_mancuso@meta.data$seurat_clusters == 1,"renamed_clusters"] = "Mancuso_exM_1" # perivascular macrophages
ex_vivo_mancuso@meta.data[ex_vivo_mancuso@meta.data$seurat_clusters == 4,"renamed_clusters"] = "Mancuso_exM_2" # perivascular macrophages
ex_vivo_mancuso@meta.data[ex_vivo_mancuso@meta.data$seurat_clusters == 3,"renamed_clusters"] = "Mancuso_exMGL_2" # Cytokine response microglia
ex_vivo_mancuso@meta.data[ex_vivo_mancuso@meta.data$seurat_clusters == 5,"renamed_clusters"] = "Mancuso_cycling"
ex_vivo_mancuso@meta.data[ex_vivo_mancuso@meta.data$seurat_clusters == 6,"renamed_clusters"] = "Mancuso_exMGL_3" # Amyloid response microglia

mg.sce = as.SingleCellExperiment(ex_vivo_mancuso)

# pseudobulk
ex_vivo_mancuso = aggregateAcrossCells(mg.sce, use_exprs_values = "counts",
                                       ids=colData(mg.sce)[,c( "renamed_clusters")]) 

ex_vivo_mancuso@colData
head(counts(ex_vivo_mancuso))



ex_vivo_mancuso = geneName_to_ensembl(counts(ex_vivo_mancuso))

# our data

media_merged= readRDS("../../data/results/1.3.QC/media_merged_noRegression.rds")
media_merged@meta.data$cell = rownames(media_merged@meta.data)
media_merged = as.SingleCellExperiment(media_merged)

media_merged = aggregateAcrossCells(media_merged, use_exprs_values = "counts",
                                    ids=colData(media_merged)[,c( "orig.ident")]) # sample

media_merged@colData
media_merged = geneName_to_ensembl(counts(media_merged))
colnames(media_merged) = c("IGBN_iMG","IM_iMG","IMBN_iMG","ITGBN_iMG","ITM_ADMEM_iMG","ITMG_iMG")

###################### merging and processing with limma
message("Merging")

combined_commongenes = merge(Abud, Douvaras, by = 0, all = TRUE)  # merge by row names (by=0 or by="row.names")
rownames(combined_commongenes) = combined_commongenes$Row.names  # give row names again
combined_commongenes = merge(combined_commongenes, Muffat, by = 0, all = TRUE)  # merge by row names (by=0 or by="row.names")
rownames(combined_commongenes) = combined_commongenes$Row.names  # give row names again
combined_commongenes = merge(combined_commongenes, Galatro, by = 0, all = TRUE)  # merge by row names (by=0 or by="row.names")
rownames(combined_commongenes) = combined_commongenes$Row.names  # give row names again

combined_commongenes = merge(combined_commongenes, media_merged, by = 0, all = TRUE)  # merge by row names (by=0 or by="row.names")
rownames(combined_commongenes) = combined_commongenes$Row.names  # give row names again
combined_commongenes = merge(combined_commongenes, ex_vivo_mancuso, by = 0, all = TRUE)  # merge by row names (by=0 or by="row.names")
rownames(combined_commongenes) = combined_commongenes$Row.names  # give row names again
combined_commongenes = combined_commongenes[,c(6:ncol(combined_commongenes))] # remove first column
combined_commongenes = na.omit(combined_commongenes)   # remove rows that contain NA values (reduces table to size of smallest table)

message(paste0("There are ",nrow(combined_commongenes)," genes remaining in the combined dataframe"))


combined_commongenes = DGEList(combined_commongenes, 
                               genes = rownames(combined_commongenes),
                               samples = colnames(combined_commongenes))
keep = edgeR::filterByExpr(combined_commongenes)

combined_commongenes = combined_commongenes[keep,]
nrow(combined_commongenes)

combined_commongenes = calcNormFactors(combined_commongenes)    # Calculate normalization factors. TMM by default


samples = c(colnames(Abud), colnames(Douvaras), 
            colnames(Muffat),  colnames(Galatro), 
            colnames(media_merged), colnames(ex_vivo_mancuso))   # create the design matrix
study = c(rep("Abud",ncol(Abud)),
          rep("Douvaras",ncol(Douvaras)),
          rep("Muffat",ncol(Muffat)),
          rep("Galatro",ncol(Galatro)),
          rep("Ours",ncol(media_merged)),
          rep("Mancuso",ncol(ex_vivo_mancuso)))

cell_type = c(rep("iMG",6), rep("iHPC",3),rep("iPSC",4),
              rep("exMonocyte",9), #rep("monocyte_CD14",5),rep("monocyte_CD14_CD16",4), # collapsing annotations to make plot clearer
              
              rep("pFMGL",3),rep("pAMGL",3),rep("exDC",3),rep("iMG_neuron",6),rep("iMG_noCD200_noCX3CL1",3),rep("iMG_noTGFb",3),
              rep("pMAC",8),  rep("pAMGL_serum",2),  rep("iMG",10),  rep("pAMGL",2), rep("PMAC_hepatic",2),
              rep("pFMGL",3),rep("iMG",5),"eMG",rep("iMG_neuron",3),  rep("iNPC",4),
              rep("exMGL",39),rep("exBrain",16),"S969" , "S947"  , "S861"  ,  "S810"  , "S805"  ,  "S755"  ,"S291"  ,  "S262"  ,  "S243"   , "S147" ,
              "iMG","iMG","iMG","iMG","iMG","iMG", # ours
              "cycling", rep("exMAC",2), rep("exMGL",3)
              
)

metadata = as.data.frame(cbind(samples,study,cell_type))
metadata$samples = as.factor(metadata$samples)
metadata$study = as.factor(metadata$study)
metadata$cell_type = as.factor(metadata$cell_type)


design = model.matrix(~study + cell_type )  

v = voom(combined_commongenes,design,plot = F) # voom normalize the read counts

v$E = removeBatchEffect(v$E,metadata$study)

plotPca = function(x, metadata. = metadata){ 
  pca1 = prcomp(t(x), retx = TRUE)
  
  percentVar = (pca1$sdev)^2 / sum(pca1$sdev^2)
  percentVar = round(100 * percentVar)
  pcs = as.data.frame(pca1$x)
  pcs = cbind(pcs,Samples = metadata$cell_type, Study=metadata$study, 
              study_samples = paste0(metadata$study,"_" ,metadata$cell_type),
              original_name = metadata$samples)
  
  
  ################# remove unwanted experiments & samples after calculating the pcs
  # 
  pcs = pcs[which(!(pcs$Samples  %in% c("S969" , "S947"  , "S861"  ,  "S810"  , "S805"  ,  "S755"  ,"S291"  ,  "S262"  ,  "S243"   , "S147" , 
                                        "cycling", "PMAC_hepatic", "iHPC","eMG","iMG_noCD200_noCX3CL1",
                                        "iMG_noTGFb"))),]
  pcs = pcs[which(!(pcs$original_name  %in% c("Douvaras_iMG_1","Douvaras_iMG_2"))),]
  pcs$labels = ""
  pcs$labels[which(!duplicated(pcs$study_samples))] = as.character(pcs$Samples[which(!duplicated(pcs$study_samples))])
  pcs$labels[which(pcs$original_name == "ITMG_iMG")] = "ours_ITMG_iMG"
  pcs$labels[which(pcs$original_name == "ITM_ADMEM_iMG")] = "ours_ITM_ADMEM_iMG"
  
  pcs$Study = relevel(pcs$Study,ref = "Ours")
  
  pcs$Samples = factor(pcs$Samples, ordered = T, levels = c("iPSC",
                                                            "exBrain","exDC","exMonocyte","iNPC",
                                                            "exMGL","exMAC",
                                                            "iMG","iMG_neuron",
                                                            "pMAC","pFMGL","pAMGL","pAMGL_serum"))
  colors = c("iPSC" = "#9A8F97" ,
             "exBrain" = "#BFCE96","exDC" = "#CAC4CE","exMonocyte" = "#FFDD4A","iNPC" = "#208AAE",
             "exMGL" = "#C4D88C","exMAC" = "#EBBAB9",
             "iMG" = "#000000","iMG_neuron" = "#4A7C59",
             "pMAC" = "#D58936","pFMGL" = "#A44200","pAMGL" = "#69140E","pAMGL_serum" = "#3C1518"
  )
  
  p = ggplot(pcs,aes(x = PC1, y = PC2, color = Samples, fill = Samples,
                     label = labels ,
                     shape = Study)) +
    geom_point(size = 5,
               alpha = 0.6,
               stroke = 1) +
    geom_text_repel(max.overlaps = Inf, box.padding = 0.5) +
    geom_point(data = pcs[pcs$labels != "",]) +
    scale_color_manual(values = colors) + 
    xlab(paste0( "PC1: " ,percentVar[ 1 ],"% variance")) +
    ylab(paste0( "PC2: ",percentVar[ 2 ],"% variance" ))
  
  
  p = p + theme(   panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   panel.background = element_blank(),
                   panel.border = element_blank(),
                   legend.key = element_blank(),
                   axis.title.y = element_text(  angle = 90, size = 24, vjust = 0.2),
                   axis.title.x = element_text(  size = 24, vjust = 0),
                   axis.text.x = element_text(  colour = "black", size = 24, vjust = 0.2, hjust = 1 ),
                   axis.text.y = element_text(  colour = "black",size = 24),
                   axis.ticks = element_line(colour = "black",size = 0.75),
                   axis.line = element_line(colour = "black",size = 0.30),
                   legend.text = element_text(  colour = "black",size = 24),
                   legend.title = element_text(  colour = "black",size = 24))
  
   #plot(p)
  
  return(p)
}

p1 = plotPca(v$E)
p1

png("../../data/results/4.Compare_external_datasets/PCA_external_datasets_FigS5a.png", 
    width = 13, height = 9, units = "in", res = 400)
plot(p1)
dev.off()

plotPca_highlights = function(x, metadata. = metadata){ 
  pca1 = prcomp(t(x), retx = TRUE)
  
  percentVar = (pca1$sdev)^2 / sum(pca1$sdev^2)
  percentVar = round(100 * percentVar)
  pcs = as.data.frame(pca1$x)
  pcs = cbind(pcs,Samples = metadata$cell_type, Study=metadata$study, 
              study_samples = paste0(metadata$study,"_" ,metadata$cell_type),
              original_name = metadata$samples)
  
  
  ################# remove unwanted experiments & samples after calculating the pcs
  # 
  pcs = pcs[which(!(pcs$Samples  %in% c("S969" , "S947"  , "S861"  ,  "S810"  , "S805"  ,  "S755"  ,"S291"  ,  "S262"  ,  "S243"   , "S147" , 
                                        "cycling", "PMAC_hepatic", "iHPC","eMG","iMG_noCD200_noCX3CL1",
                                        "iMG_noTGFb"))),]
  pcs = pcs[which(!(pcs$original_name  %in% c("Douvaras_iMG_1","Douvaras_iMG_2"))),]
  pcs$labels = ""
  pcs$labels[which(!duplicated(pcs$study_samples))] = as.character(pcs$Samples[which(!duplicated(pcs$study_samples))])
  pcs$labels[which(pcs$original_name == "ITMG_iMG")] = "ours_ITMG_iMG"
  pcs$labels[which(pcs$original_name == "ITM_ADMEM_iMG")] = "ours_ITM_ADMEM_iMG"
  pcs$labels[which(pcs$original_name == "IGBN_iMG")] = "ours_IGBN_iMG"
  pcs$labels[which(pcs$original_name == "IM_iMG")] = "ours_IM_iMG" 
  pcs$labels[which(pcs$original_name == "IMBN_iMG")] = "ours_IMBN_iMG"
  pcs$labels[which(pcs$original_name == "ITGBN_iMG")] = "ours_ITGBN_iMG"
  pcs$labels[which(pcs$labels == "exDC")] = ""
  pcs$labels[which(pcs$labels == "exMonocyte")] = ""
  pcs$labels[which(pcs$labels == "exMGL")] = ""
  pcs$labels[which(pcs$labels == "exMAC")] = ""
  pcs$labels[which(pcs$labels == "iMG_neuron")] = ""
  pcs$labels[which(pcs$labels == "exBrain")] = ""
  pcs$labels[which(pcs$labels == "iPSC")] = ""
  pcs$labels[which(pcs$labels == "pMAC")] = ""
  pcs$labels[which(pcs$labels == "iMG")] = ""
  pcs$labels[which(pcs$labels == "iNPC")] = ""
  

  pcs$Samples = factor(pcs$Samples, ordered = T, levels = c("iPSC",
                                                            "exBrain","exDC","exMonocyte","iNPC",
                                                            "exMGL","exMAC",
                                                            "iMG","iMG_neuron",
                                                            "pMAC","pFMGL","pAMGL","pAMGL_serum"))
  pcs$Study = relevel(pcs$Study,ref = "Ours")
  colors = c("iPSC" = "#9A8F97" ,
             "exBrain" = "#9A8F97","exDC" = "#9A8F97","exMonocyte" = "#9A8F97","iNPC" = "#9A8F97",
             "exMGL" = "#9A8F97","exMAC" = "#9A8F97",
             "iMG" = "#000000","iMG_neuron" = "#9A8F97",
             "pMAC" = "#D58936","pFMGL" = "#A44200","pAMGL" = "#69140E","pAMGL_serum" = "#3C1518"
  )
  
  p = ggplot(pcs,aes(x = PC1, y = PC2,
                     label = labels ,
                     shape = Study)) +
    geom_point(size = 5,
               alpha = 0.6,
               stroke = 1,color = "lightgrey") +
    geom_point(data = pcs[pcs$Study ==  "Ours",],color = "red", fill = "red",size = 5,) +
    geom_text_repel(max.overlaps = Inf, box.padding = 0.5) +
    geom_point(data = pcs[pcs$Sample %in%  c("pFMGL","pAMGL","pAMGL_serum"),],
               color = "black", fill = "black",size = 5,) +
    
    xlab(paste0( "PC1: " ,percentVar[ 1 ],"% variance")) +
    ylab(paste0( "PC2: ",percentVar[ 2 ],"% variance" ))
  
  
  p = p + theme(   panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   panel.background = element_blank(),
                   panel.border = element_blank(),
                   legend.key = element_blank(),
                   axis.title.y = element_text(  angle = 90, size = 24, vjust = 0.2),
                   axis.title.x = element_text(  size = 24, vjust = 0),
                   axis.text.x = element_text(  colour = "black", size = 24, vjust = 0.2, hjust = 1 ),
                   axis.text.y = element_text(  colour = "black",size = 24),
                   axis.ticks = element_line(colour = "black",size = 0.75),
                   axis.line = element_line(colour = "black",size = 0.30),
                   legend.text = element_text(  colour = "black",size = 24),
                   legend.title = element_text(  colour = "black",size = 24))
  
  #plot(p)
  
  return(p)
}

p2 = plotPca_highlights(v$E)
p2

png("../../data/results/4.Compare_external_datasets/PCA_external_datasets_highlights_FigS5b.png", 
    width = 13, height = 9, units = "in", res = 400)
plot(p2)
dev.off()
