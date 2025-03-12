#Import dependencies
library(ggpubr)
library(ggrepel)
library(Seurat)
library(scran)
library(harmony)
library(scuttle)
library(scater)
library(bluster)
library(tidyverse)
library(qs)
library(ggplot2)
library(readxl)
library(purrr)
library(kableExtra)
set.seed(1234)
#load obj
T_cells<-qread("/mnt/disks/disk/full_dataset_qs/T_cells_jan25_final.qs")
clusters_to_remove <- grep("^dbl", unique(Idents(T_cells)),value=TRUE)
clusters_to_remove<-unique(c(clusters_to_remove, "Plt", "lowquality", "Lowq", "lowq","Plt"))
Idents(T_cells)<-T_cells$sub_cluster
keep<-setdiff(unique(Idents(T_cells)), clusters_to_remove)
T_cells<- subset(T_cells, idents = keep)
meta<-read_excel("~/Immune_project/full_dataset_analysis/Meta_Immune_proj.xlsx")
meta<-distinct(meta, Sample_ID, .keep_all = TRUE)

T_cells@meta.data$new_barcodes<-rownames(T_cells@meta.data)
T_cells@meta.data<-T_cells@meta.data[-c(44:51)]
T_cells@meta.data<-left_join(T_cells@meta.data, meta, by="Sample_ID")
rownames(T_cells@meta.data)<-T_cells@meta.data$new_barcodes
unique(T_cells$Therapy)
#REMOVE LEN AND SCREENING, KEEP CART TALQ TEC HEALTHY
T_cells_filt<-subset(T_cells, subset=Timepoint%in%c("Post","Pre","Day14","Day30","Day100","Day56","6Months",   
                                                    "Healthy")&Therapy%in%c("TEC","TALQ","CART","Healthy"))

T_cells_filt<-RenameIdents(T_cells_filt, c('CD8+_GZMK_TEM'="CD8+_GZMK+_TEM"))
T_cells_filt$sub_cluster<-Idents(T_cells_filt)

PRE_PB_CD8<-subset(T_cells_filt, Timepoint=="Pre"&Therapy%in%c("TEC","TALQ","CART")&Tissue=="PB"& sub_cluster%in%c("CD8+_DR+_TEM", 
                                                                                                                  "CD8+_GZMB_TEM",
                                                                                                                   "CD8+_GZMK+_TEM",
                                                                                                                   "CD8+_ZEB2_TEM",
                                                                                                                   "CD8+_Tex",
                                                                                                                   "cycling_Tcells"))
PRE_PB_CD8$sub_cluster <- factor(PRE_PB_CD8$sub_cluster)
PRE_PB_CD8$sub_cluster <- droplevels(PRE_PB_CD8$sub_cluster)
df<-as.data.frame(table(PRE_PB_CD8$Sample_ID, PRE_PB_CD8$sub_cluster))
colnames(df)<-c("Sample_ID", "sub_cluster", "n")
df<-df%>%mutate(cluster_sample=paste0(Sample_ID,"_",sub_cluster))
df<-df%>%filter(n>10)
cs<-unique(df$cluster_sample)
# -------------------------------------------------------------------------
#pseudobulk
# -------------------------------------------------------------------------
# pseudobulk the counts based on donor-condition-celltype
pseudo_cd8 <- AggregateExpression(PRE_PB_CD8, assays = "RNA", return.seurat = TRUE, group.by = c("Sample_ID", "sub_cluster", "Disease"))
# each 'cell' is a donor-condition-celltype pseudobulk profile
tail(Cells(pseudo_cd8))

pseudo_cd8$celltype.dis <- paste(pseudo_cd8$sub_cluster, pseudo_cd8$Disease, sep = "_")
pseudo_cd8$cluster_sample <- paste( pseudo_cd8$Sample_ID,pseudo_cd8$sub_cluster, sep = "_")
unique(pseudo_cd8$cluster_sample)
pseudo_cd8$cluster_sample<-gsub("-","_", pseudo_cd8$cluster_sample)
#remove those observation with less than 10 cells
pseudo_cd8<-subset(pseudo_cd8, cluster_sample%in% cs)


output_dir<-"/mnt/disks/disk/full_dataset_qs/DEG/"
dir.create(output_dir)
clusters<-unique(pseudo_cd8$sub_cluster)
pseudo_cd8@assays$RNA$counts<-round(pseudo_cd8@assays$RNA$counts)

for (i in seq_along(clusters)){
  sub<-subset(pseudo_cd8, sub_cluster==clusters[[i]])
  Idents(sub) <- "celltype.dis"
  bulk.de <- FindMarkers(object = sub, 
                                          ident.1 = paste0(clusters[[i]],"_RRMM"), 
                                          ident.2 = paste0(clusters[[i]],"_HRSMM"),
                                          test.use = "DESeq2")
write.csv(bulk.de, paste0(output_dir,clusters[[i]],"PB_BL_RRMM_vs_HRSMM.csv"))
}


# -------------------------------------------------------------------------
#MAST
# -------------------------------------------------------------------------
library(SCOPfunctions)
clusters<-unique(PRE_PB_CD8$sub_cluster)
output_dir<-"/mnt/disks/disk/full_dataset_qs/DEG/MAST/"
for (i in seq_along(clusters)){
  sub<-subset(PRE_PB_CD8, sub_cluster==clusters[[i]])
  Idents(sub) <- "Disease"
  test<-DE_MAST_RE_seurat(object=sub, random_effect.vars = "Sample_ID", ident.1 = "RRMM", ident.2 = "HRSMM", base=2)
  write.csv(test, paste0(output_dir,clusters[[i]],"PB_BL_RRMM_vs_HRSMM.csv"))
}
