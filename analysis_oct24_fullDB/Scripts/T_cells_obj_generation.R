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
library(kableExtra)
library(patchwork)

#load Tcell1
T_cells<-qread("/mnt/disks/cellranger/full_dataset_qs/T_cells.clustered.qs")
###############################
T_cells <- FindVariableFeatures(object=T_cells, selection.method = "vst", nfeatures = 2000, verbose = FALSE) 
#remove cilta and idecell from variable features
genes_names <- VariableFeatures(T_cells)
car <- genes_names[grep("ciltacabtagene|idecabtagene", genes_names)]
genes_to_keep <- !(genes_names %in% car)
genes_names <- genes_names[genes_to_keep]
VariableFeatures(T_cells)<-genes_names

print("Scaling...")
T_cells <- ScaleData(object=T_cells, verbose = FALSE)
# #Identify the 15 most highly variable genes
# ranked_variable_genes <- VariableFeatures(T_cells)
# top_genes <- ranked_variable_genes[1:15]
# # Plot the average expression and variance of these genes
# # With labels to indicate which genes are in the top 15
# #p <- VariableFeaturePlot(T_cells)
# #LabelPoints(plot = p, points = top_genes, repel = TRUE)

print("Running PCA...")
T_cells <- RunPCA(object=T_cells, npcs = 30, verbose = FALSE)
T_cells <- RunHarmony(T_cells, group.by.vars = "Pool", plot=TRUE)
print("Clustering corrected dataset...")
T_cells <- RunUMAP(T_cells, reduction = "harmony", dims=1:30)
T_cells <- FindNeighbors(T_cells, reduction = "harmony", dims=1:30)
T_cells <- FindClusters(T_cells, resolution = 1.2)
Idents(T_cells)<-T_cells$RNA_snn_res.1.2
T_cells[['ident']]<-Idents(T_cells)
p0<-DimPlot(T_cells, label=T)

keep<-setdiff(unique(Idents(T_cells)), c(21,15,19,25,27))
T_cells2<-subset(T_cells,idents=keep )

# -------------------------------------------------------------------------
T_cells2 <- FindVariableFeatures(object=T_cells2, selection.method = "vst", nfeatures = 2000, verbose = FALSE) 
#remove cilta and idecell from variable features
genes_names <- VariableFeatures(T_cells2)
car <- genes_names[grep("ciltacabtagene|idecabtagene", genes_names)]
genes_to_keep <- !(genes_names %in% car)
genes_names <- genes_names[genes_to_keep]
VariableFeatures(T_cells2)<-genes_names

print("Scaling...")
T_cells2 <- ScaleData(object=T_cells2, verbose = FALSE)
# #Identify the 15 most highly variable genes
# ranked_variable_genes <- VariableFeatures(T_cells2)
# top_genes <- ranked_variable_genes[1:15]
# # Plot the average expression and variance of these genes
# # With labels to indicate which genes are in the top 15
# #p <- VariableFeaturePlot(T_cells2)
# #LabelPoints(plot = p, points = top_genes, repel = TRUE)

print("Running PCA...")
T_cells2 <- RunPCA(object=T_cells2, npcs = 30, verbose = FALSE)
T_cells2 <- RunHarmony(T_cells2, group.by.vars = "Pool", plot=TRUE)
print("Clustering corrected dataset...")
T_cells2 <- RunUMAP(T_cells2, reduction = "harmony", dims=1:30)
T_cells2 <- FindNeighbors(T_cells2, reduction = "harmony", dims=1:30)
T_cells2 <- FindClusters(T_cells2, resolution = 1.2)
Idents(T_cells2)<-T_cells2$RNA_snn_res.1.2
T_cells2[['ident']]<-Idents(T_cells2)
p0<-DimPlot(T_cells2, label=T)




ggsave("~/Immune_project/analysis_oct24_fullDB/plots/Tcells/T_cells_UMAP1.pdf", plot = p0, width=7, height=4)

keep.cells <- T_cells2@meta.data |> tibble::rownames_to_column("cells") |> group_by(Sample_ID) |> 
  sample_frac(size=0.1) |> pull(cells)
T_cells3 <- T_cells2
T_cells2 <- subset(T_cells2, cells=keep.cells)

#downsample before FindMarkers
output_file <- "/mnt/disks/cellranger/full_dataset_qs/markers/T_cells_markers2.csv"
# Controlla se il file esiste
if (file.exists(output_file)) {
  file_content <- read.csv(output_file)
  print(file_content) } else {
    print("Performing Differential Expression for all clusters...")
    markers <- FindAllMarkers(T_cells2)
    top_markers <- markers %>%
      group_by(cluster) %>%
      dplyr::filter(pct.1 > 0.2) %>%
      arrange(desc(avg_log2FC)) %>%
      slice_head(n = 20)
    
    print(top_markers)
    write_csv(top_markers, output_file)
  }


#
FeaturePlot(T_cells, c("scDblFinder.class","predicted_doublets","Doublet","percent.mt","lognCount"), order=T, label=T)
Idents(T_cells3)<-T_cells3$seurat_clusters
T_cells3<-RenameIdents(T_cells3, c('0'="Gzmb+CD8_Tem", 
                                 '1'="CD4Naive", '2'="DR+CD8_Tem",
                                 '3'="CD4Tcm", '4'="Kir+CD8_Tem",
                                 '5'="Exhausted-like", 
                                 '6'="CD8Naive", '7'="Th17-like", 
                                 '8'="CD4Naive",
                                 '9'="Gzmk+CD8_Tem", 
                                 '10'="Tregs", 
                                 '11'="Th1-like",
                                 '12'="CD4Tcm",
                                 '13'="Zeb2+CD8_Tem",
                                 '14'="CD4Tcm",
                                 '15'="Mait",
                                 '16'="IFN+CD4Tcells",
                                 '17'="Th2-like",
                                 '18'="Cycling Tcells",
                                 '19'="dbl.granulo.CD3",
                                 '20'="CD4Naive", 
                                 '21'="Tgd",
                                 '22'="Cycling Tcells",
                                 '23'="Cycling Tcells",
                                 '24'="dbl.granulo.CD3",
                                 '25'="Cycling Tcells",
                                 '26'="Cycling Tcells",
                                 '27'="CD4Naive",
                                 '28'="dbl.granulo.CD3",
                                 '29'="CD8Naive", 
                                 '30'="dbl.B.Cd3",
                                 '31'="Cycling Tcells",
                                 '32'="dbl.B.Cd3",
                                 '33'="dbl.B.CD3",
                                 '34'="dbl.B.CD3"))

T_cells3[['ident']]<-Idents(T_cells3)
p0<-DimPlot(T_cells3, label=T)

qsave(T_cells3,"/mnt/disks/cellranger/full_dataset_qs/Tcells_annotated_final.qs")
FeaturePlot(T_cells2, c("ciltacabtagene", "idecabtagene"), order=T)

T_cells[['ident']]<-Idents(T_cells)
DimPlot(T_cells)
Idents(T_cells)<-T_cells$seurat_clusters
T_cells<-RenameIdents(T_cells, c('21'="dbl.Mono.CD3",'15'="dbl.Mono.CD3",'19'="dbl.Mono.CD3",'25'="dbl.Mono.CD3",'27'="dbl.Mono.CD3"))

T_cells$sub_cluster <- as.character(Idents(T_cells))

T_cells$sub_cluster[Cells(T_cells3)] <- paste(Idents(T_cells3))


