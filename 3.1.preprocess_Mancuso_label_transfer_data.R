# QC and save seurat object for Mancuso iPSC-derived human microglia grown ex vivo in mouse

library(Seurat)
library(patchwork)
library(ggplot2)
library(future)
library(singleCellNet)
library(tidyverse)

options(stringsAsFactors = FALSE)

# change the current plan to access parallelization
plan("multiprocess", workers = 14)
options(future.globals.maxSize= 2097152000) # 2Gb

plan()
# Mancuso nat neuroscience 2019 ESC-derived microglia matured in mice (ex vivo) - raw gene expression counts per cell
# 10X 
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4107633
# RM37: microglia naive
# RM38: microglia + scrambled
# RM49: microglia + amyloid beta treatment

ex_vivo_mancuso = read.table("ESC_derived_microglia_in_mice/GSE137444_chimera_human_h9.tsv.gz")
ex_vivo_mancuso = CreateSeuratObject(counts = ex_vivo_mancuso, project = "ex_vivo_Mancuso")
table(ex_vivo_mancuso@meta.data$orig.ident) # 3 mice replicates
ex_vivo_mancuso[["percent.mt"]] <- PercentageFeatureSet(ex_vivo_mancuso, pattern = "^MT-")
dim(ex_vivo_mancuso@meta.data)

# how many cells per orig ident
table(ex_vivo_mancuso$orig.ident)

## Rename orig.ident

new.ids <- c("Naive", "Scrambled", "Amyloid_treatment")
names(new.ids) <- levels(ex_vivo_mancuso)
ex_vivo_mancuso <- RenameIdents(ex_vivo_mancuso, new.ids)
ex_vivo_mancuso$cell_type = Idents(ex_vivo_mancuso)

p1 = ex_vivo_mancuso@meta.data %>% 
  ggplot(aes(color=orig.ident, x=nCount_RNA, fill= cell_type)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  xlab("UMI counts") + 
  NoLegend()
p2 = ex_vivo_mancuso@meta.data %>% 
  ggplot(aes(color=cell_type, x=nFeature_RNA, fill= cell_type)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300) +
  ylab("Cell density") +
  xlab("Number of genes") + 
  NoLegend()



p3 = ex_vivo_mancuso@meta.data %>% 
  ggplot(aes(color=cell_type, x=percent.mt, fill=cell_type)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  geom_vline(xintercept = 10) +
  ylab("Cell density") +
  xlab("% of mitochondrial genes")


# Cells have been filtered for around 500 UMIs per cell min in each cell when creating the Seurat object, or does cellranger do that?
p1 + p2 + p3

filter_seurat = function(seurat_object){
  
  message("Performing filter by number of genes, mitochondrial percentage and UMI counts.")
  seurat_object <- subset(seurat_object, subset = nFeature_RNA > 300  & percent.mt < 10 & nCount_RNA > 1000)
  message("Now the object has ", dim(seurat_object)[1], " genes and ", dim(seurat_object)[2], " cells.")
  return(seurat_object)
}

## subset by gene expression
ex_vivo_mancuso = subset(ex_vivo_mancuso, subset = PLP1 < 1) # Remove oligodendrocytes
ex_vivo_mancuso = subset(ex_vivo_mancuso, subset = CD48 < 1) # Remove B cells
# Remove cells from cluster 7 - very far and few
cells_to_exclude =  c("RM37_AGCTTGATCACATAGC.1", "RM37_AGTCTTTTCCAAGTAC.1", "RM37_ATCACGACATTAACCG.1", "RM37_ATCTACTGTAGCTGCC.1",
                      "RM37_ATGTGTGCAAGAAGAG.1",
"RM37_CACTCCAGTGTCCTCT.1", "RM37_CGCGTTTAGAGTTGGC.1", "RM37_GACCTGGTCTGGAGCC.1", "RM37_GAGTCCGGTATAGGGC.1", "RM37_TTAGGCATCGTTACGA.1",
"RM38_CATGGCGAGACACTAA.1", "RM38_CTACATTCAATGGAGC.1" ,"RM38_CTACGTCTCTTAGAGC.1" ,"RM38_GACACGCTCCTCTAGC.1", "RM38_GGACATTCACCTTGTC.1",
"RM38_GTCTCGTGTCTAGAGG.1", "RM38_TTGGAACGTTTACTCT.1" ,"RM38_TTGTAGGCAAGAAGAG.1", "RM49_ACCTTTAGTCTTGATG.1", "RM49_ACGAGCCAGGAGCGTT.1",
"RM49_AGCTTGAGTCCATGAT.1", "RM49_ATCATGGAGAGTTGGC.1", "RM49_ATGTGTGGTCATCGGC.1", "RM49_CACACCTGTATAGTAG.1", "RM49_CACATAGGTCTACCTC.1",
"RM49_CCACCTAAGTACGTAA.1" ,"RM49_GACTACAGTGGTCTCG.1", "RM49_GCAGTTACAGTCACTA.1", "RM49_GGATTACAGTGCAAGC.1", "RM49_GGGATGAAGAACAACT.1",
"RM49_GTAGGCCCAGACGCTC.1", "RM49_TACGGGCGTGATGCCC.1", "RM49_TTTGTCAAGTGAACAT.1")
cells_to_include = rownames(ex_vivo_mancuso@meta.data[!rownames(ex_vivo_mancuso@meta.data) %in% cells_to_exclude,])
ex_vivo_mancuso = subset(ex_vivo_mancuso, cells =cells_to_include) 

ex_vivo_mancuso = filter_seurat(ex_vivo_mancuso)
cc.genes <- readLines("../../resources/regev_lab_cell_cycle_genes_Kowalczyk_2015.txt")
S.genes <- cc.genes[1:43]
G2M.genes <- cc.genes[44:97]

ex_vivo_mancuso <- CellCycleScoring(ex_vivo_mancuso, s.features=S.genes, g2m.features=G2M.genes, set.ident=F)
ex_vivo_mancuso <- SCTransform(ex_vivo_mancuso, verbose = TRUE)


calculate_PCA_UMAP_neighbors_clusters = function(seurat_object){
  seurat_object <- RunPCA(seurat_object, verbose = FALSE)
  
  seurat_object <- FindNeighbors(seurat_object, dims = 1:30, verbose = FALSE)
  seurat_object <- FindClusters(seurat_object,
                                verbose = FALSE,resolution = 0.3)
  
  seurat_object <- RunUMAP(seurat_object, dims = 1:30, verbose = FALSE)
  
  
  
  return(seurat_object)
}

ex_vivo_mancuso = calculate_PCA_UMAP_neighbors_clusters(ex_vivo_mancuso)

UMAPS_per_lib = function(seurat_object) {
  p1 = DimPlot(seurat_object, label = TRUE) + NoLegend() + ggtitle("Clusters")
  p12 = DimPlot(seurat_object, label = F,group.by = "cell_type")  + ggtitle("Sample")
  # Changing transparency
  p12[[1]]$layers[[1]]$aes_params$alpha = .4 
  
  p2 = FeaturePlot(seurat_object, label = F, features = "percent.mt", order = TRUE)  + ggtitle("Mit %") 
  
  p3 = FeaturePlot(seurat_object, label = F, features  = "nFeature_RNA", order = TRUE) + 
    ggtitle("Number of genes")
  p4 = FeaturePlot(seurat_object, label = F, features  = "nCount_RNA", order = TRUE) + 
    ggtitle("UMI counts")
  
  p5 = DimPlot(seurat_object, label = F,  group.by = "Phase")  + ggtitle("Cell cycle")
  p5[[1]]$layers[[1]]$aes_params$alpha = .4 
  

  p = (p1 | p12) / (p2 | p3) / (p4 | p5) 
  # / (p6 | p7) 
  
  return(p)
}
png("plots/GSE137444_chimera_human_h9_UMAPs.png",
    width = 7, height = 9, units = "in", res = 400)
p = UMAPS_per_lib(ex_vivo_mancuso)
p
dev.off()
## There's clustering by treatment within each of 2 separate big clusters
## remove cluster 7 because it's not clear what it is
# rownames(ex_vivo_mancuso@meta.data[ex_vivo_mancuso$seurat_clusters==7,])

## Investigate and remove other cell types as they do
# run this before removing cells over certain expression threshold
# other_markers = c("CCRL2", #neutrophils
#                   "MRC1", #CNS-associated macrophages (CAMs), they are perivascular macrophages but they identify by common marker (higher in pvM than microglia)
#                   "CD52", "CD48", # B cells
#                   "NKG7", "CD247", "CD7",  # natural killter T cells
#                   "CCR2", "LY6C11",  #monocytes 
#                   "CLU", #astrocytes 
#                   "MBP","PLP1", #oligodendrocytes 
#                   "TOP2A") #cycling cells 
# png("plots/GSE137444_chimera_human_h9_featureplot_other_markers.png",
#     width = 15, height = 15, units = "in", res = 400)
# FeaturePlot(ex_vivo_mancuso, label = F, features  = other_markers, order = TRUE) 
# dev.off()

top10_markers = function(seurat_object){
  seurat_object.markers <-
    FindAllMarkers(
      seurat_object,
      only.pos = FALSE,
      min.pct = 0.1,
      logfc.threshold = 0.25
    )
  top10 <-
    seurat_object.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  
  return(top10)
}


top10 = top10_markers(ex_vivo_mancuso)

png("plots/GSE137444_chimera_human_h9_heatmap_top10_markers.png",
    width = 7, height = 9, units = "in", res = 400)
DoHeatmap(subset(ex_vivo_mancuso, downsample = 100), features = top10$gene) + NoLegend()

dev.off()

## barplots of composition of cells (cell_type) per cluster
table(ex_vivo_mancuso$cell_type,ex_vivo_mancuso)
### txt files created from from supplementary data Table S3 - columns N, O, P

HM = readLines("SalaFrigerio_HomeostaticMicroglia.txt")
ARM = readLines("SalaFrigerio_ActivatedResponseMicroglia.txt")
CRM = readLines("Cytokine_response_microglia.txt")
macrophages = readLines("Prinz_CNS_associated_macrophages.txt")


microglia_markers = c( "P2RY12", # P2Y12, Microglia
                       "TMEM119", # TMEM119, Microglia
                       "FCRL1",# FCRLs, Microglia
                       "HEXB",  # Microglia
                       "SALL1", # SALL1, microglia
                       "SPI1") # Pu.1, microglia

microglia_macrophage_markers = c("ITGAM", # CD11B, Microglia &  macrophage marker
                                 "PTPRC", # CD45, higher in macrophage than microglia
                                 "MRC1", #CD206, higher in macrophage than microglia
                                 "CD68", # Microglia &  macrophage marker
                                 "CX3CR1", #Microglia &  macrophage marker
                                 "ADGRE1") # F4/80 

peripheral_macrophage_markers = c( "CD44", # CD44
                                   "SIGLEC1", # CD169
                                   "MYB") 


png( "plots/GSE137444_chimera_human_h9_heatmap_HM_per_cluster_heatmap_downsample100.png", 
     width = 14, height = 12, units = "in", res = 400)
DoHeatmap(subset(ex_vivo_mancuso, downsample = 100), features = HM) + NoLegend()
dev.off()

png( "plots/GSE137444_chimera_human_h9_heatmap_ARM_per_cluster_heatmap_downsample100.png", 
     width = 14, height = 7, units = "in", res = 400)
DoHeatmap(subset( ex_vivo_mancuso, downsample = 100), features = ARM) + NoLegend()
dev.off()

png( "plots/GSE137444_chimera_human_h9_heatmap_CRM_per_cluster_heatmap_downsample100.png", 
     width = 12, height = 3, units = "in", res = 400)
DoHeatmap(subset( ex_vivo_mancuso, downsample = 100), features = CRM) + NoLegend()
dev.off()

png( "plots/GSE137444_chimera_human_h9_heatmap_macrophages_per_cluster_heatmap_downsample100.png", 
     width = 12, height = 3, units = "in", res = 400)
DoHeatmap(subset( ex_vivo_mancuso, downsample = 100), features = macrophages) + NoLegend()
dev.off()

# by orig ident
png( "plots/GSE137444_chimera_human_h9_heatmap_HM_per_sample_heatmap_downsample100.png", 
     width = 14, height = 12, units = "in", res = 400)
DoHeatmap(subset(ex_vivo_mancuso, downsample = 100), features = HM, group.by = "cell_type") + NoLegend()
dev.off()

png( "plots/GSE137444_chimera_human_h9_heatmap_ARM_per_sample_heatmap_downsample100.png", 
     width = 14, height = 7, units = "in", res = 400)
DoHeatmap(subset( ex_vivo_mancuso, downsample = 100), features = ARM, group.by = "cell_type") + NoLegend()
dev.off()

png( "plots/GSE137444_chimera_human_h9_heatmap_CRM_per_sample_heatmap_downsample100.png", 
     width = 12, height = 3, units = "in", res = 400)
DoHeatmap(subset( ex_vivo_mancuso, downsample = 100), features = CRM, group.by = "cell_type") + NoLegend()
dev.off()

png( "plots/GSE137444_chimera_human_h9_heatmap_macrophages_per_sampler_heatmap_downsample100.png", 
     width = 12, height = 3, units = "in", res = 400)
DoHeatmap(subset( ex_vivo_mancuso, downsample = 100), features = macrophages, group.by = "cell_type") + NoLegend()
dev.off()


png( "plots/GSE137444_chimera_human_h9_heatmap_microglia_markers_literature_per_cluster_dotplot.png", 
     width = 6, height = 6, units = "in", res = 400)
DotPlot( ex_vivo_mancuso, features =microglia_markers) + RotatedAxis()
dev.off()

png( "plots/GSE137444_chimera_human_h9_heatmap_microglia_and_macrophage_markers_literature_per_cluster_dotplot.png", 
     width = 6, height = 6, units = "in", res = 400)
DotPlot( ex_vivo_mancuso, features =microglia_macrophage_markers) + RotatedAxis()
dev.off()

png( "plots/GSE137444_chimera_human_h9_heatmap_peripheral_macrophage_markers_literature_per_cluster_dotplot.png", 
     width = 6, height = 6, units = "in", res = 400)
DotPlot( ex_vivo_mancuso, features =peripheral_macrophage_markers) + RotatedAxis()
dev.off()


## module scoring

## Scoring cells according to microglia markers found in cluster 20
ex_vivo_mancuso <- AddModuleScore(
  object = ex_vivo_mancuso,
  features = list(peripheral_macrophage_markers),
  ctrl = 100,
  name = 'pvMacrophage_features'
)

ex_vivo_mancuso <- AddModuleScore(
  object = ex_vivo_mancuso,
  features = list(CRM),
  ctrl = 100,
  name = 'CRmicroglia_features'
)

p1 = FeaturePlot(ex_vivo_mancuso, label = F, features = "pvMacrophage_features1", order = TRUE)  + ggtitle("pvM scoring") 

p2 = FeaturePlot(ex_vivo_mancuso, label = F, features  = "CRmicroglia_features1", order = TRUE) + 
  ggtitle("CRM scoring")

png( "plots/GSE137444_chimera_human_h9_featureplot_pvM_CRM_scoring_UMAP.png", 
     width = 4, height = 6, units = "in", res = 400)
p1 + p2
dev.off()
# unclear

# Fibroblast markers. See Tabib et al., J Invest Dermatol, 2018
fibro = c( "FMO1",
           "SFRP2",  # FRP-2, FRP2, SARP1, SDF-5
           "CRABP1", "COL11A1", "FMO2", "PRG4", "C2ORF40",
           "COL3A1", "FBN1", "PDGFRA", "VIM", "DCN",
           "THY1", # also marker of hematopoietic stem cells, CD90
           "COL1A1","COL1A2",
           #"S100A4", #  FSP1, fibroblast specific protein 1, actuallt not good marker of fibroblast 
           #because it's detected in other cell types, including macrophages (see publication)
           "PDGFRA", "CD34",
           "CD73","CD146")  

png("plots/GSE137444_chimera_human_h9_heatmap_markers_fibroblast_per_cluster_dotplot.png", 
    width = 6, height = 6, units = "in", res = 400)
DotPlot( ex_vivo_mancuso, features =unique(fibro)) + RotatedAxis()
dev.off()

## Label transfer 

# this file is created in select_perivascular_markers/4.Inspect_30k_enriched_microglia.R
brain10x=readRDS("../../../../resources/1M_mice_brain_10X/seurat_30K_microglia_enriched_filtered_CCycleMouseRegression_humanOrthologs.rds")
head(brain10x@meta.data) 
brain10x

# Remove cell types that have less than 100 cells so I can do training on at least 50 cells
to_subset = names(which(table(brain10x@meta.data$renamed_clusters)>100))
brain10x = subset(brain10x, 
                  cells = brain10x@meta.data[brain10x@meta.data$renamed_clusters %in% to_subset,"Barcode"])


# to singleCellNet format
train_query_format = function(seurat_object){
  SCN_list = list()
  seuratfile_SCN= extractSeurat(seurat_object, exp_slot_name = "counts")
  SCN_list[["samp"]] =  seuratfile_SCN$sampTab
  SCN_list[["exp"]]= seuratfile_SCN$expDat
  return(SCN_list)
}
train = train_query_format(brain10x)

# Query data: my Hipsci data
query = train_query_format(ex_vivo_mancuso)
query$samp$cell=rownames(query$samp)
# subset the datasets based on common genes before proceeding to training
commonGenes=intersect(rownames(train$exp), rownames(query$exp))
length(commonGenes)
train$exp = train$exp[commonGenes, ]
query$exp = query$exp[commonGenes, ]

## Training
# Recommended minimum number of cells per category is 50

stList=singleCellNet::splitCommon(sampTab = train$samp, ncells = 50, dLevel = "renamed_clusters")
stTrain=stList[[1]]
expTrain = train$exp[,stTrain$Barcode]

stTest = stList[[2]]
expTest = train$exp[,stTest$Barcode]

sum(colnames(expTrain) %in% colnames(expTest)) # No shared cells in train and test set

# Train and asess the performance of the classifier
# If you increase nTopGenes and nTopGenePairs, you may get a even better classifier performance on query data!
class_info =
  scn_train(
    stTrain = stTrain,
    expTrain = expTrain,
    nTopGenes = 20,
    nRand = 70,
    nTrees = 1000,
    nTopGenePairs = 50,
    dLevel = "renamed_clusters",
    colName_samp = "Barcode"
  )
# There are  244  classification genes
# Finding top pairs
# nPairs =  190  for  Neuron_glutamanergic_c0 
# nPairs =  190  for  Neuron_glutamanergic_c7 
# nPairs =  190  for  Neuron_glutamanergic_c2 
# nPairs =  190  for  Neuron_glutamanergic_c5 
# nPairs =  190  for  Astrocytes 
# nPairs =  190  for  Endothelial 
# nPairs =  190  for  Smooth_muscle 
# nPairs =  190  for  Neuron_gabaergic_SST 
# nPairs =  190  for  Blood 
# nPairs =  190  for  Neuron_unknown_c9 
# nPairs =  190  for  Neuron_gabaergic_c4 
# nPairs =  190  for  Oligodendrocytes 
# nPairs =  190  for  Oligodendrocyte_precursors 
# nPairs =  190  for  Microglia_c1 
# nPairs =  190  for  Microglia_c14 
# There are 374 top gene pairs
# Finished pair transforming the data
# Number of missing genes  0 

message(" Here, we use the held-out data from the 10x Brain dataset to assess theTop-pair Random Forest (TP-RF) classifer.")
classRes_val_all = scn_predict(class_info[['cnProc']], expTest, nrand = 50)

tm_heldoutassessment = assess_comm(
  ct_scores = classRes_val_all,
  stTrain = stTrain,
  stQuery = stTest,
  dLevelSID = "Barcode",
  classTrain = "renamed_clusters",
  classQuery = "renamed_clusters"
)

message("Assessing the accuracy of classification of held-out test data with Precision-Recall curve")

head(tm_heldoutassessment)

png("plots/seurat_30K_microglia_enriched_filtered_CCycleMouseRegression_training_held-out_data_asessment_metrics.png", 
    width = 3, height = 4, units = "in", res = 400)
plot_metrics(tm_heldoutassessment)
dev.off()
png("plots/seurat_30K_microglia_enriched_filtered_CCycleMouseRegression_training_held-out_data_asessment_PRCurve.png", 
    width = 8, height = 5, units = "in", res = 400)
plot_PRs(tm_heldoutassessment)
dev.off()

# Method on query data
message("Predicting on query data (ex vivo Mancuso)")
cr_exvivo <- scn_predict(class_info[['cnProc']], query$exp, nrand = 50)

notmatch = cr_exvivo[ ,!(colnames(cr_exvivo) %in% query$samp$cell)]
head(notmatch)

# # Subset to only columns present
stQuery <-
  assign_cate(classRes = cr_exvivo[, (colnames(cr_exvivo) %in% query$samp$cell)],  
              sampTab =  query$samp,
              cThresh = 0.5)

sgrp<-as.vector(query$samp$seurat_clusters)
names(sgrp) = as.vector(query$samp$cell)

grpRand =rep("rand", 50)
names(grpRand) = paste("rand_", 1:50, sep='')
sgrp = append(sgrp, grpRand)  
# sgrp should be cell as name and category as character e.g.
# input_hipsci_lib1_AAACCTGAGAATAGGG input_hipsci_lib1_AAACCTGAGATGCGAC input_hipsci_lib1_AAACCTGAGCGATCCC 
# "c0"                    "c1"                    "c1" 
# heatmap classification result
## cr_exvivo should not contain rand (should be 50 rows smaller than sgrp)


png("plots/seurat_30K_microglia_enriched_filtered_CCycleMouseRegression_barplot_per_cluster.png", 
    width = 4, height = 3, units = "in", res = 400)
sc_hmClass(cr_exvivo, sgrp, max=5000, isBig=TRUE, cCol=F, font=8)
dev.off()

for_seurat = t(cr_exvivo[, (colnames(cr_exvivo) %in% query$samp$cell)])
for_seurat = cbind(for_seurat,stQuery)

ex_vivo_mancuso_annotated = AddMetaData(ex_vivo_mancuso,for_seurat)

DimPlot(ex_vivo_mancuso_annotated,group.by = "category", split.by = "category",ncol = 3)

## give categories the same color 
## Give always same colors per category
group.colors = c(Microglia_c1 = "#73937E", Perivascular_macrophages = "#218380", Neuron_glutamanergic_c0 = "#8f84d8", 
                 rand = "#d81159" )


long = ex_vivo_mancuso_annotated@meta.data %>%
  group_by(seurat_clusters) %>%
  count(category)  %>% mutate(Percent = n / sum(n)*100)

png("plots/seurat_30K_microglia_enriched_filtered_CCycleMouseRegression_barplot_per_cluster.png", 
    width = 5, height = 3, units = "in", res = 400)
ggplot(data=long, aes(x=seurat_clusters, y=Percent, fill=category)) +
  geom_bar(stat="identity") + theme_bw() + xlab("Cluster") + scale_fill_manual(values=group.colors)
dev.off()

long_lib_clusters = ex_vivo_mancuso_annotated@meta.data %>%
  group_by(orig.ident,seurat_clusters) %>%
  count(category)  %>% mutate(Percent = n / sum(n)*100)

p1 = ggplot(data=long_lib_clusters[long_lib_clusters$orig.ident == "RM37",], aes(x=seurat_clusters, y=n, fill=category)) +
  geom_bar(stat="identity") + theme_bw() + xlab("Cluster") + ggtitle("Mancuso naive microglia") + scale_fill_manual(values=group.colors)

p2 = ggplot(data=long_lib_clusters[long_lib_clusters$orig.ident == "RM38",], aes(x=seurat_clusters, y=n, fill=category)) +
  geom_bar(stat="identity") + theme_bw() + xlab("Cluster") + ggtitle("Mancuso scrambled microglia") + scale_fill_manual(values=group.colors)

p3 = ggplot(data=long_lib_clusters[long_lib_clusters$orig.ident == "RM49",], aes(x=seurat_clusters, y=n, fill=category)) +
  geom_bar(stat="identity") + theme_bw() + xlab("Cluster") + ggtitle("Mancuso Abeta microglia") + scale_fill_manual(values=group.colors)


png("plots/seurat_30K_microglia_enriched_filtered_CCycleMouseRegression_barplot_per_library_and_cluster.png", 
    width = 12, height = 4, units = "in", res = 400)
p1 + p2 + p3 + plot_layout(guides = "collect")
dev.off()


long = ex_vivo_mancuso_annotated@meta.data %>%
  group_by(cell_type) %>%
  count(category)  %>% mutate(Percent = n / sum(n)*100)

png("plots/seurat_30K_microglia_enriched_filtered_CCycleMouseRegression_barplot_per_library.png", 
    width = 6, height = 3, units = "in", res = 400)
ggplot(data=long, aes(x=cell_type, y=Percent, fill=category)) +
  geom_bar(stat="identity") + theme_bw() + xlab("Library")  + scale_fill_manual(values=group.colors)

dev.off()

## Table of percentages 
table = long[c(1,2,4)] %>% pivot_wider(names_from = cell_type,values_from = Percent, values_fill = 0)
write.table(table,"plots/table_from_labelTransfer_seurat_30K_microglia_enriched_filtered_CCycleMouseRegression_barplot_per_library.txt",
            sep = "\t",quote = F,row.names = F)


### just RM37 - microglia no treatment
ex_vivo_naive <- subset(ex_vivo_mancuso_annotated, subset = cell_type=="Naive")

p1 = DimPlot(ex_vivo_naive, label = TRUE) + NoLegend() + ggtitle("Clusters")
p12 = DimPlot(ex_vivo_naive, label = F,group.by = "category")  + ggtitle("Label transfer")
# Changing transparency
p12[[1]]$layers[[1]]$aes_params$alpha = .4 

p2 = FeaturePlot(ex_vivo_naive, label = F, features = "percent.mt", order = TRUE)  + ggtitle("Mit %") 

p5 = DimPlot(ex_vivo_naive, label = F,  group.by = "Phase")  + ggtitle("Cell cycle")
p5[[1]]$layers[[1]]$aes_params$alpha = .4 

png("plots/naive_UMAPs.png", 
    width = 10, height = 10, units = "in", res = 400)
p = (p1 | p12) / (p2 | p5) 
p
dev.off()

saveRDS(ex_vivo_mancuso_annotated,"ex_vivo_mancuso_filtered_annotated_allSamples.rds")
