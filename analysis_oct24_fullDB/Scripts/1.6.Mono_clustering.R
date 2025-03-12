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

#merge
Mono<-qread("/mnt/disks/disk/Monos.qs")
###############################
Mono <- FindVariableFeatures(object=Mono, selection.method = "vst", nfeatures = 2000, verbose = FALSE) 
print("Scaling...")
Mono <- ScaleData(object=Mono, verbose = FALSE)
#Identify the 15 most highly variable genes
ranked_variable_genes <- VariableFeatures(Mono)
top_genes <- ranked_variable_genes[1:15]
# Plot the average expression and variance of these genes
# With labels to indicate which genes are in the top 15
#p <- VariableFeaturePlot(Mono)
#LabelPoints(plot = p, points = top_genes, repel = TRUE)

print("Running PCA...")
Mono <- RunPCA(object=Mono, npcs = 30, verbose = FALSE)
Mono <- RunHarmony(Mono, group.by.vars = "Pool", plot=TRUE)
print("Clustering corrected dataset...")
Mono <- RunUMAP(Mono, reduction = "harmony", dims=1:30)
Mono <- FindNeighbors(Mono, reduction = "harmony", dims=1:30)
Mono <- FindClusters(Mono, resolution = 1)
Mono[['ident']]<-Idents(Mono)


#downsample before FindMarkers
output_file <- "/mnt/disks/disk/markers/Monos_markers.csv"
# Controlla se il file esiste
if (file.exists(output_file)) {
  file_content <- read.csv(output_file)
  print(file_content) } else {
    print("Performing Differential Expression for all clusters...")
    markers <- FindAllMarkers(Mono)
    top_markers <- markers %>%
      group_by(cluster) %>%
      dplyr::filter(pct.1 > 0.3) %>%
      arrange(desc(avg_log2FC)) %>%
      slice_head(n = 15)
    
    print(top_markers)
    write_csv(top_markers, output_file)
  }


Mono<-RenameIdents(Mono, c('0'="IL1B+Mono", '1'="CD14+Mono",'2'="IL1B+Mono", '3'="CD16+Mono",'4'="IFN+Mono",
                           '5'="CD14+Mono", '6'="CD14+Mono", '7'="CD14+Mono", '8'="DR+Mono", '9'="CD16+Mono",
                           '10'="dbl.Mono.CD3",'11'="dbl.Mono.CD3", '12'="dbl.Mono.CD3", '13'="C1Q+Macrophage",
                           '14'="dbl.Mono.CD3", '15'="CD16+Mono",'16'="CD16+Mono",'17'="CD14+Mono"))

Mono[['ident']]<-Idents(Mono)
Mono$sub_cluster<-Idents(Mono)
DimPlot(Mono, label = T)+NoLegend()
FeaturePlot(Mono, c("CD14","FCGR3A","MX1","IL1B","CXCL8","HLA-DRA"), label=T, order=T)
qsave(Mono, "/mnt/disks/disk/Mono_annotated.qs")
