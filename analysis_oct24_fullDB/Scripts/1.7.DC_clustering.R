library(qs)
library(Seurat)
library(harmony)
DC<-qread("/mnt/disks/disk/DC.qs")
DC <- FindVariableFeatures(object=DC, selection.method = "vst", nfeatures = 2000, verbose = FALSE) 
print("Scaling...")
DC <- ScaleData(object=DC, verbose = FALSE)
#Identify the 15 most highly variable genes
# ranked_variable_genes <- VariableFeatures(DC)
# top_genes <- ranked_variable_genes[1:15]
# # Plot the average expression and variance of these genes
# # With labels to indicate which genes are in the top 15
# p <- VariableFeaturePlot(DC)
# LabelPoints(plot = p, points = top_genes, repel = TRUE)

print("Running PCA...")
DC <- RunPCA(object=DC, npcs = 30, verbose = FALSE)
DC <- RunHarmony(DC, group.by.vars = "Pool",plot_convergence = TRUE)
print("Clustering corrected dataset...")
DC <- RunUMAP(DC, reduction = "harmony", dims=1:30)
DC <- FindNeighbors(DC, reduction = "harmony", dims=1:30)
DC <- FindClusters(DC, resolution = 1)


output_dir <-"/mnt/disks/disk/markers"

output_file <- paste0(output_dir, "/DC_markers.csv")

# Controlla se il file esiste
if (file.exists(output_file)) {
  file_content <- read.csv(output_file)
  print(file_content) } else {
    print("Performing Differential Expression for all clusters...")
    markers <- FindAllMarkers(DC)
    top_markers <- markers %>%
      group_by(cluster) %>%
      dplyr::filter(pct.1 > 0.3) %>%
      arrange(desc(avg_log2FC)) %>%
      slice_head(n = 15)
    
    print(top_markers)
    write_csv(top_markers, output_file)
  }

###
DC[['ident']]<-Idents(DC)
DimPlot(DC, label=T)
FeaturePlot(DC, c("percent.mt", "CD14", "CD3D","CD1C","CLEC9A","IL3RA", "Doublet","scDblFinder.class","predicted_doublets"), order=T)
DC<-RenameIdents(DC, c('0'="pDC", '1'="cDC2", '2'="cDC2", '3'="IL1B+cDC2", '4'="Mono-DC", '5'="Mono-DC", '6'="Mono-DC",
                       '7'="IFN+cDC2", '8'="dbl.CD3.DC", '9'="dbl.CD3.DC", '10'="lowq",'11'="cDC2",
                       '12'="cDC1", '13'="lowq", '14'="dbl.CD3.DC", '15'="AS-DC", '16'="dbl.DC.Mono", '17'="Cycling cDC2", '18'="Cycling pDC", 
                       '19'="dbl.CD3.DC",'20'="Plt"))
DC[['ident']]<-Idents(DC)
DimPlot(DC, label=T)

qsave(DC, "/mnt/disks/disk/DC_annotated.qs")

