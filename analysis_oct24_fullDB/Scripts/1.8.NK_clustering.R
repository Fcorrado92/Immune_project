library(qs)
library(Seurat)
library(harmony)
NK_cells<-qread("/mnt/disks/cellranger/full_dataset_qs/NK.qs")
NK_cells <- FindVariableFeatures(object=NK_cells, selection.method = "vst", nfeatures = 2000, verbose = FALSE) 
print("Scaling...")
NK_cells <- ScaleData(object=NK_cells, verbose = FALSE)
#Identify the 15 most highly variable genes
# ranked_variable_genes <- VariableFeatures(NK_cells)
# top_genes <- ranked_variable_genes[1:15]
# # Plot the average expression and variance of these genes
# # With labels to indicate which genes are in the top 15
# p <- VariableFeaturePlot(NK_cells)
# LabelPoints(plot = p, points = top_genes, repel = TRUE)

print("Running PCA...")
NK_cells <- RunPCA(object=NK_cells, npcs = 30, verbose = FALSE)
ElbowPlot(NK_cells, ndims = 30)
NK_cells <- RunHarmony(NK_cells, group.by.vars = "Pool",plot_convergence = TRUE)
print("Clustering corrected dataset...")
NK_cells <- RunUMAP(NK_cells, reduction = "harmony", dims=1:20)
NK_cells <- FindNeighbors(NK_cells, reduction = "harmony", dims=1:20)
NK_cells <- FindClusters(NK_cells, resolution = 1)
NK_cells[['ident']]<-Idents(NK_cells)
DimPlot(NK_cells, label=T)+NoLegend()
FeaturePlot(NK_cells,"CD3D", label=T)
NK_cells<-qread("/mnt/disks/cellranger/full_dataset_qs/NK_annotated_final.qs")
T_cells_in_NK<-subset(NK_cells, idents=c(11,7,0,8))
qsave(T_cells_in_NK, "/mnt/disks/cellranger/full_dataset_qs/T_in_NK.qs")
# -------------------------------------------------------------------------

NK_cells2<-subset(NK_cells, idents=c(2,16,6,13,9,4,5,12,14,15,10,1,3))
NK_cells2 <- FindVariableFeatures(object=NK_cells2, selection.method = "vst", nfeatures = 2000, verbose = FALSE) 
print("Scaling...")
NK_cells2 <- ScaleData(object=NK_cells2, verbose = FALSE)
#Identify the 15 most highly variable genes
# ranked_variable_genes <- VariableFeatures(NK_cells2)
# top_genes <- ranked_variable_genes[1:15]
# # Plot the average expression and variance of these genes
# # With labels to indicate which genes are in the top 15
# p <- VariableFeaturePlot(NK_cells2)
# LabelPoints(plot = p, points = top_genes, repel = TRUE)

print("Running PCA...")
NK_cells2 <- RunPCA(object=NK_cells2, npcs = 30, verbose = FALSE)
ElbowPlot(NK_cells2, ndims = 30)
NK_cells2 <- RunHarmony(NK_cells2, group.by.vars = "Pool",plot_convergence = TRUE)
print("Clustering corrected dataset...")
NK_cells2 <- RunUMAP(NK_cells2, reduction = "harmony", dims=1:20)
NK_cells2 <- FindNeighbors(NK_cells2, reduction = "harmony", dims=1:20)
NK_cells2 <- FindClusters(NK_cells2, resolution = 1)
NK_cells2[['ident']]<-Idents(NK_cells2)
DimPlot(NK_cells2, label=T)+NoLegend()

# -------------------------------------------------------------------------


keep.cells <- NK_cells2@meta.data |> tibble::rownames_to_column("cells") |> group_by(Sample_ID) |> 
  sample_frac(size=0.1) |> pull(cells)
NK_cells3 <- NK_cells2
NK_cells2 <- subset(NK_cells2, cells=keep.cells)

markers<-FindAllMarkers(NK_cells2)

output_file <- paste0(output_dir, "/NK_cells_markers.csv")

# Controlla se il file esiste
if (file.exists(output_file)) {
  file_content <- read.csv(output_file)
  print(file_content) } else {
    print("Performing Differential Expression for all clusters...")
    markers <- FindAllMarkers(NK_cells)
    top_markers <- markers %>%
      group_by(cluster) %>%
      dplyr::filter(pct.1 > 0.3) %>%
      arrange(desc(avg_log2FC)) %>%
      slice_head(n = 15)
    
    print(top_markers)
    write_csv(top_markers, output_file)
  }
###
FeaturePlot(NK_cells3, c("percent.mt", "CD14", "CD3D","CD1C","CLEC9A","IL3RA", "Doublet","scDblFinder.class","predicted_doublets"), order=T)
Idents(NK_cells3)<-NK_cells3$seurat_clusters
NK_cells3<-RenameIdents(NK_cells3, c('0'="CD56dim_NK_cells", '1'="CD56br_NK_cells",
                                     '2'="CD56dim_NK_cells", '3'="CD56dim_NK_cells",
                                     '4'="aNK_cells", '5'="aNK_cells", '6'="CD56dim_NK_cells",
                       '7'="lowq", '8'="CD56dim_NK_cells", '9'="CD56dim_NK_cells", '10'="dbl.Mono.NK",
                       '11'="CD56dim_NK_cells",
                       '12'="CD56dim_NK_cells", '13'="CD56dim_NK_cells", '14'="dbl.Mono.NK", '15'="CD56dim_NK_cells", 
                       '16'="dbl.Mono.NK"))
NK_cells3[['ident']]<-Idents(NK_cells3)
DimPlot(NK_cells3, label=T)

qsave(NK_cells3, "/mnt/disks/disk/full_dataset_qs/NK_cells_annotated.qs")
