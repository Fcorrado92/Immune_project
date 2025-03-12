#Import dependencies
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

Bcells<-qread("/mnt/disks/disk/Bcells.qs")
#FindVariable Feat
Bcells <- FindVariableFeatures(object=Bcells, selection.method = "vst", nfeatures = 2000, verbose = FALSE) 
print("Scaling...")
genes_names <- VariableFeatures(Bcells)
Ig_genes <- genes_names[grep("^IGHV|^IGHJ|^IGKV|^IGKJ|^IGKD|^IGKC|^IGLV|^IGLJ|^IGLC|^IGLD|^IGHD", genes_names)]
genes_to_keep <- !(genes_names %in% Ig_genes)
genes_names <- genes_names[genes_to_keep]
VariableFeatures(Bcells)<-genes_names

Bcells <- ScaleData(object=Bcells, verbose = FALSE)
print("Running PCA...")
Bcells <- RunPCA(object=Bcells, npcs = 30, verbose = FALSE)
Bcells <- RunHarmony(Bcells, group.by.vars = "Pool",plot_convergence = TRUE)
print("Clustering corrected dataset...")
Bcells <- RunUMAP(Bcells, reduction = "harmony", dims=1:30)
Bcells <- FindNeighbors(Bcells, reduction = "harmony", dims=1:30)
Bcells <- FindClusters(Bcells, resolution = 1)
Bcells[['ident']]<-Idents(Bcells)

#FindMarkers
output_dir <- "/mnt/disks/disk/markers/"
output_file <- paste0(output_dir, "Bcells_markers.csv")

# Controlla se il file esiste
if (file.exists(output_file)) {
  file_content <- read.csv(output_file)
  print(file_content) } else {
    print("Performing Differential Expression for all clusters...")
    markers <- FindAllMarkers(Bcells)
    top_markers <- markers %>%
      group_by(cluster) %>%
      dplyr::filter(pct.1 > 0.3&avg_log2FC>0&p_val_adj<0.05) %>%
      arrange(desc(avg_log2FC)) %>%
      slice_head(n = 20)
    
    print(top_markers)
    write_csv(top_markers, output_file)
  }

Bcells[['ident']]<-Idents(Bcells)
DimPlot(Bcells, label=T)

# Rename identities for the specified clusters
Idents(Bcells)<-Bcells$seurat_clusters
Bcells <- RenameIdents(Bcells, 
                       '11' = "dbl.CD3.Bcells", 
                       '8' = "dbl.CD3.Bcells", 
                       '15' = "Plt", 
                       '6' = "dbl.Mono.Bcells", 
                       '17' = "dbl.Mono.Bcells", 
                       '18' = "dbl.Mono.Bcells", 
                       '12' = "lowquality", 
                       '16' = "lowquality", 
                       '9' = "lowquality", 
                       '20' = "lowquality")

# Get unique identities for the newly renamed clusters
renamed_clusters <- c("dbl.CD3.Bcells", "Plt", "dbl.Mono.Bcells", "lowquality")
# Find identities that are in Bcells or renamed clusters but not overlapping with other identities
keep <- setdiff(unique(Idents(Bcells)), renamed_clusters)

Bcells2<-subset(Bcells, idents=keep)
Bcells2[['ident']]<-Idents(Bcells2)
DimPlot(Bcells2)

Bcells2 <- FindVariableFeatures(object=Bcells2, selection.method = "vst", nfeatures = 2000, verbose = FALSE) 
print("Scaling...")
genes_names <- VariableFeatures(Bcells2)
Ig_genes <- genes_names[grep("^IGHV|^IGHJ|^IGKV|^IGKJ|^IGKD|^IGKC|^IGLV|^IGLJ|^IGLC|^IGLD|^IGHD", genes_names)]
genes_to_keep <- !(genes_names %in% Ig_genes)
genes_names <- genes_names[genes_to_keep]
VariableFeatures(Bcells2)<-genes_names

Bcells2 <- ScaleData(object=Bcells2, verbose = FALSE)
print("Running PCA...")
Bcells2 <- RunPCA(object=Bcells2, npcs = 30, verbose = FALSE)
Bcells2 <- RunHarmony(Bcells2, group.by.vars = "Pool",plot_convergence = TRUE)
print("Clustering corrected dataset...")
Bcells2 <- RunUMAP(Bcells2, reduction = "harmony", dims=1:30)
Bcells2 <- FindNeighbors(Bcells2, reduction = "harmony", dims=1:30)
Bcells2 <- FindClusters(Bcells2, resolution = 1)
Bcells2[['ident']]<-Idents(Bcells2)

#FindMarkers
output_dir <- "/mnt/disks/disk/markers/"
output_file <- paste0(output_dir, "Bcells2_markers.csv")

# Controlla se il file esiste
if (file.exists(output_file)) {
  file_content <- read.csv(output_file)
  print(file_content) } else {
    print("Performing Differential Expression for all clusters...")
    markers <- FindAllMarkers(Bcells2)
    top_markers <- markers %>%
      group_by(cluster) %>%
      dplyr::filter(pct.1 > 0.3&avg_log2FC>0&p_val_adj<0.05) %>%
      arrange(desc(avg_log2FC)) %>%
      slice_head(n = 20)
    
    print(top_markers)
    write_csv(top_markers, output_file)
  }

Bcells2[['ident']]<-Idents(Bcells2)
DimPlot(Bcells2, label = T)

FeaturePlot(Bcells2, c("CD19", "MME", "CD34", "IGLL5", "CD79A","MS4A1"), label = T, order = T)
FeaturePlot(Bcells2, c("percent.mt", "lognCount", "Doublet", "predicted_doublets", "scDblFinder.class"), label = T, order = T)
FeaturePlot(Bcells2, c("TCL1A", "IGHM", "CD1C", "MZB1", "IGHG2","IGHA","SOX5"), label = T, order = T)


Bcells2<-RenameIdents(Bcells2, c(
  '3'="Pre-B",
  '8'="lowq",
  '11'="Pro-B",
  '13'="dbl.Mono.Bcells",
  '14'="lowq",
  '16'="lowq",
  '18'="lowq",
  '19'="lowq",
  '20'="lowq",
  '21'="lowq"))
renamed_clusters <- c("Pre-B", "Pro-B", "dbl.Mono.Bcells", "lowq")
# Find identities that are in Bcells or renamed clusters but not overlapping with other identities
keep <- setdiff(unique(Idents(Bcells2)), renamed_clusters)

Bcells3<-subset(Bcells2, idents=keep)

Bcells3 <- FindVariableFeatures(object=Bcells3, selection.method = "vst", nfeatures = 2000, verbose = FALSE) 
print("Scaling...")
genes_names <- VariableFeatures(Bcells3)
Ig_genes <- genes_names[grep("^IGHV|^IGHJ|^IGKV|^IGKJ|^IGKD|^IGKC|^IGLV|^IGLJ|^IGLC|^IGLD|^IGHD", genes_names)]
genes_to_keep <- !(genes_names %in% Ig_genes)
genes_names <- genes_names[genes_to_keep]
VariableFeatures(Bcells3)<-genes_names

Bcells3 <- ScaleData(object=Bcells3, verbose = FALSE)
print("Running PCA...")
Bcells3 <- RunPCA(object=Bcells3, npcs = 30, verbose = FALSE)
Bcells3 <- RunHarmony(Bcells3, group.by.vars = "Pool",plot_convergence = TRUE)
print("Clustering corrected dataset...")
Bcells3 <- RunUMAP(Bcells3, reduction = "harmony", dims=1:30)
Bcells3 <- FindNeighbors(Bcells3, reduction = "harmony", dims=1:30)
Bcells3 <- FindClusters(Bcells3, resolution = 1)
Bcells3[['ident']]<-Idents(Bcells3)

#FindMarkers
output_dir <- "/mnt/disks/disk/markers/"
output_file <- paste0(output_dir, "Bcells3_markers.csv")

# Controlla se il file esiste
if (file.exists(output_file)) {
  file_content <- read.csv(output_file)
  print(file_content) } else {
    print("Performing Differential Expression for all clusters...")
    markers <- FindAllMarkers(Bcells3)
    top_markers <- markers %>%
      group_by(cluster) %>%
      dplyr::filter(pct.1 > 0.3&avg_log2FC>0&p_val_adj<0.05) %>%
      arrange(desc(avg_log2FC)) %>%
      slice_head(n = 20)
    
    print(top_markers)
    write_csv(top_markers, output_file)
  }
#TBC
FeaturePlot(Bcells3,c("CD19", "MS4A1", "CD38", "IGHM", "IGHD", "TCL1A", "IL4R", "SELL"), order = T, label=T)
#NBC
FeaturePlot(Bcells3,c("CD19", "MS4A1", "IGHM", "IGHD", "IL4R", "SELL", "CCR7", "FCER2", "ABCB1"), order = T, label=T)
#GCB
FeaturePlot(Bcells3,c("TNFSF9", "CXCR5", "CD83"), order = T, label=T)
#MZB
FeaturePlot(Bcells3,c("CD1C", "IGHM", "PLD4", "LINC01857"), order = T, label=T)
#ABC
FeaturePlot(Bcells3,c("TBX21", "ITGAX", "FCRL5", "ENC1", "TNFRSF1B", "SOX5", "MPP6"), order = T, label=T)

Bcells3[['ident']]<-Idents(Bcells3)
DimPlot(Bcells3, label = T)

Bcells$sub_cluster <- as.character(Idents(Bcells))
Bcells$sub_cluster[Cells(Bcells2)] <- paste(Idents(Bcells2))
Bcells$sub_cluster[Cells(Bcells3)] <- paste(Idents(Bcells3))

Bcells[['ident']]<-Bcells$sub_cluster
DimPlot(Bcells, label = T)
qsave(Bcells, "/mnt/disks/disk/ Bcells_annotated.qs")