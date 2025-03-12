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

int<-qread("~/Immune_project/full_dataset_analysis/integrated_filtered.qs")
#plot markers
int[['ident']]<-Idents(int)
p1<-FeaturePlot(int, c("CD3D","NKG7","MS4A1","TNFRSF17","CD19","CD34","BEX2","HBB","AZU1","LYZ","CD1C","IL3RA"), order=T, label=T)
ggsave(plot=p1, "~/Immune_project/full_dataset_analysis/plots/features_integrated.pdf", width=20, height=14)
#look at batches
DimPlot(int, split.by = "Tissue_Cohort", ncol=4)
#assign major cluster
Idents(int)<-int$seurat_clusters
integrated_filt <- RenameIdents(integrated_filt, c('12'="NKcells",'4'="NKcells",'3' = "T-cells", '1' = "T-cells", '41' = "T-cells", 
                           '35' = "T-cells", '9' = "T-cells", '10' = "T-cells", 
                           '42' = "T-cells", '39' = "T-cells", '0' = "T-cells",'23' = "T-cells",'15' = "T-cells",'6' = "Mono", '16' = "dbl.PLT.Mono", '5' = "Mono", 
                           '34' = "Mono", '14' = "Mono", '7' = "Mono", 
                           '13' = "Mono", '24' = "Mono", '2' = "Mono", '40' = "Mono", '11'="cDC", '22'="pDC", '21'="EP", '18'="EP",'20'="HSC",
                           '19'="GMP",'30'="Other",'32'="Other",'31' = "Bcells", '8' = "Bcells", '29' = "Bcells", 
                           '33' = "Bcells", '28' = "Bcells", '37' = "Bcells", '25' = "Bcells",'17' = "Plasmacells", '38' = "Plasmacells", '36' = "Plasmacells",'26'="Other",'27'="Other"))
#visualize major cluster
int$major_clusters<-Idents(int)
int[['ident']]<-Idents(int)
DimPlot(int,label=T)

#subset celltypes
# Subset T-cells and save
Tcells <- subset(int, idents = "T-cells")
qsave(Tcells, "/mnt/disks/cellranger/full_dataset_qs/Tcells.qs")

# Subset Mono cells and save
Monos <- subset(int, idents = "Mono")
qsave(Monos, "/mnt/disks/cellranger/full_dataset_qs/Monos.qs")

# Subset DC cells and save
DCs <- subset(int, idents = c("cDC","pDC"))
qsave(DCs, "/mnt/disks/cellranger/full_dataset_qs/DCs.qs")

# Subset B-cells and save
Bcells <- subset(int, idents = "Bcells")
qsave(Bcells, "/mnt/disks/cellranger/full_dataset_qs/Bcells.qs")

# Subset NK cells and save
NK <- subset(int, idents = "NKcells")
qsave(NK, "/mnt/disks/cellranger/full_dataset_qs/ NK.qs")
