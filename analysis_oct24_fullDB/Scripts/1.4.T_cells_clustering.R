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
T_cells<-qread("/mnt/disks/cellranger/full_dataset_qs/Tcells.qs")
T_in_NK<-qread("/mnt/disks/cellranger/full_dataset_qs/T_in_NK.qs")
T_cells<-merge(T_cells, T_in_NK)
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
ElbowPlot(T_cells)
T_cells <- RunHarmony(T_cells, group.by.vars = "Pool", plot=TRUE)
print("Clustering corrected dataset...")
T_cells <- RunUMAP(T_cells, reduction = "harmony", dims=1:30)
T_cells <- FindNeighbors(T_cells, reduction = "harmony", dims=1:30)
T_cells <- FindClusters(T_cells, resolution = 2)
Idents(T_cells)<-T_cells$RNA_snn_res.1.2
T_cells[['ident']]<-Idents(T_cells)
p0<-DimPlot(T_cells, label=T)
FeaturePlot(T_cells, "CD3D", label=T)
keep<-setdiff(unique(Idents(T_cells)), c(24,32,28,14))
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
T_cells2 <- RunUMAP(T_cells2, reduction = "harmony", dims=1:20)
T_cells2 <- FindNeighbors(T_cells2, reduction = "harmony", dims=1:20)
T_cells2 <- FindClusters(T_cells2, resolution = 1.2)
Idents(T_cells2)<-T_cells2$RNA_snn_res.1.2
T_cells2[['ident']]<-Idents(T_cells2)
p0<-DimPlot(T_cells2, label=T)




ggsave("~/Immune_project/analysis_oct24_fullDB/plots/Tcells/T_cells_UMAP1.pdf", plot = p0, width=7, height=4)

keep.cells <- T_cells2@meta.data |> tibble::rownames_to_column("cells") |> group_by(Sample_ID) |> 
  sample_frac(size=0.1) |> pull(cells)
T_cells3 <- T_cells2
qsave(T_cells3,"/mnt/disks/disk/full_dataset_qs/T_cells_jan25.qs")

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
Idents(T_cells2)<-T_cells2$RNA_snn_res.1.2
T_cells2<-RenameIdents(T_cells2, c('0'="GZMB+_CD8+_TEM", '1'="Kir+CD8_Tem",'10'="CD4+_TEM",
                                   '11'="Exhausted-like", '12'="CD4+_TEM",
                                   '13'="CD4+Naive", '14'="Tregs",
                                   '15'="Mait", 
                                   '16'="dbl.NK.Tcells", '17'="lowq", 
                                   '18'="cycling-Tcells",
                                   '19'="cycling-Tcells", 
                                   '2'="DR+_GZMK+_CD8+_TEM", 
                                   '20'="CD4+Naive",
                                   '21'="lowq",
                                   '22'="lowq",
                                   '23'="lowq",
                                   '24'="lowq",
                                   '25'="lowq",
                                   '26'="lowq",
                                   '27'="lowq",
                                   '3'="CD4+Naive",
                                   '4'="CD4+_TEM", 
                                   '5'="CD4+_TCM", 
                                   '6'="CD8+Naive", 
                                   '7'="ZEB2+_CD8+_TEM",
                                   '8'="CD4+Naive",
                                   '9'="CD4+_TEM"))

T_cells2[['ident']]<-Idents(T_cells2)
p0<-DimPlot(T_cells2, label=T)

qsave(T_cells2,"/mnt/disks/disk/full_dataset_qs/T_cells_jan25_sub.qs")
FeaturePlot(T_cells2, c("ciltacabtagene", "idecabtagene"), order=T)

T_cells[['ident']]<-Idents(T_cells)
DimPlot(T_cells, label = T)
Idents(T_cells)<-T_cells$seurat_clusters
T_cells<-RenameIdents(T_cells, c('14'="dbl.Mono.CD3",'28'="dbl.Mono.CD3",'32'="dbl.Mono.CD3",
                                 '24'="dbl.Mono.CD3"))

T_cells$sub_cluster <- as.character(Idents(T_cells))

T_cells$sub_cluster[Cells(T_cells2)] <- paste(Idents(T_cells2))

T_cells[['ident']]<-T_cells$sub_cluster
qsave(T_cells,"/mnt/disks/disk/full_dataset_qs/T_cells_jan25_ann.qs")

T_cells<-qread("/mnt/disks/disk/full_dataset_qs/T_cells_jan25_ann.qs")
# try different resolution-tgd --------------------------------------------

DimPlot(T_cells2)
clusters_to_remove <- grep(
  "^dbl", 
  as.character(unique(Idents(T_cells2))),
  value = TRUE
)
clusters_to_remove<-unique(c(clusters_to_remove, "Plt", "lowquality", "Lowq", "lowq","Plt"))
keep<-setdiff(unique(Idents(T_cells2)), clusters_to_remove)
T_cells2<- subset(T_cells2, idents = keep)
T_cells2 <- FindClusters(T_cells2, resolution =c(1.2,1.5,1.8,2))
T_cells2$sub_cluster<-T_cells2$ident
qsave(T_cells2,"/mnt/disks/disk/full_dataset_qs/T_cells_jan25_sub.qs")
T_cells<-qread("/mnt/disks/disk/full_dataset_qs/T_cells_jan25_sub.qs")
# -------------------------------------------------------------------------
Idents(T_cells2)<-T_cells2$RNA_snn_res.2
T_cells2<-qread("/mnt/disks/disk/full_dataset_qs/T_cells_jan25_sub.qs")
markers<-FindAllMarkers(T_cells2)
write_csv(markers, "/mnt/disks/cellranger/full_dataset_qs/markers/T_cell_markers3.csv")
markers<-read_csv("/mnt/disks/cellranger/full_dataset_qs/markers/T_cells_markers3.csv")
top_markers <- markers %>%
  group_by(cluster) %>%
  dplyr::filter(pct.1 > 0.2) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 20)

T_cells2<-RenameIdents(T_cells2, c('0'="CD8+_GZMB+_TEM", 
                                   '1'="CD8+_KIR+_TEM",
                                   '10'="CD4+_TCM", '11'="CD4+_TCM",
                                   '12'="Treg", '13'="CD8+_TCM",
                                   '14'="CD8+_Tex", '15'="Th17",
                                   '16'="CD4+_TCM", '17'="CD8+_GZMB+_TEM",
                                   '18'="Mait",
                                   '19'="CD8+_KIR+_TEM", '2'="CD4+_Naive", 
                                   '20'="Th1", '21'="lowq",
                                   '22'="CD4+TfH", '23'="lowq", '24'="lowq", '25'="cycling_Tcells",
                                   '26'="CD8+_GZMK+_TEM", '27'="Th1",
                                   '28'="lowq", '29'="Tgd", '3'="CD8+_GZMB+_TEM",
                                   '30'="cycling_Tcells", '31'="cycling_Tcells", '32'="lowq", '33'="CD8+_Naive",
                                   '34'="lowq", '4'="CD8+_GZMK+_TEM", '5'="CD8+_Naive", '6'="CD4+_Naive", '7'="Th2",
                                   '8'="CD8+_DR+_TEM",'9'="CD4+_Naive"))

Idents(T_cells)<-T_cells$sub_cluster
T_cells$sub_cluster <- as.character(Idents(T_cells))
T_cells$sub_cluster[Cells(T_cells2)] <- paste(Idents(T_cells2))
T_cells[['ident']]<-T_cells$sub_cluster
Idents(T_cells)<-T_cells$sub_cluster
qsave(T_cells,"/mnt/disks/disk/full_dataset_qs/T_cells_jan25_final2.qs")
markers<-read_csv("/mnt/disks/cellranger/full_dataset_qs/markers/T_cells_markers2.csv")

HLA_DR_markers<-FindMarkers(draft, ident.1=8)
HLA_DR_markers$gene<-rownames(HLA_DR_markers)
HLA_DR_markers_sign<-HLA_DR_markers%>%filter(p_val_adj<0.05)



# -------------------------------------------------------------------------
#select baseline HD and pts CD8_TEM and compare expression of markers between exhausted and
#remaining TEM cells; only PB to avoid to have replicate samples PB and BM from some pts

# -------------------------------------------------------------------------


T_cells_filt_pre_pb<-  subset(T_cells_filt, Timepoint%in%c("Pre","Healthy"))

#select EM+EX
TEM<-c(unique(grep("GZMK|GZMB|KIR|DR|ex", T_cells_filt_pre_pb$sub_cluster, value=TRUE)))
EX_TEM<-subset(T_cells_filt_pre_pb, sub_cluster%in%TEM)
Idents(EX_TEM)<-EX_TEM$sub_cluster
# -------------------------------------------------------------------------
#rename TEMS
EX_TEM$sub_cluster<-Idents(EX_TEM)

#create metadata for filtering(sample+cluster)
EX_TEM@meta.data <- EX_TEM@meta.data %>%
  mutate(sample_cluster = paste(Sample_ID, sub_cluster, sep = "_"))

#look at how many cells x samples in the two clusters remove clusters with <10 cells
df<-as.data.frame(table(EX_TEM$Sample_ID, EX_TEM$sub_cluster))
colnames(df)<-c("Sample_ID","sub_cluster","n")
df<-df%>%mutate(sample_cluster=paste0(Sample_ID,"_",sub_cluster))
df<-df%>%group_by(sample_cluster)%>%filter(n>=10)
retain<-unique(df$sample_cluster)

EX_TEM_Filt<-subset(EX_TEM, subset=sample_cluster %in%retain)
# pseudobulk the counts based on donor-condition-celltype
pseudo<-AggregateExpression(EX_TEM_Filt, assays = "RNA",return.seurat = T, 
                            group.by = c("Sample_ID", "sub_cluster"))


# pseudobulk DEG with filtered dataset ------------------------------------
EX_TEM_Filt@meta.data%>%group_by(Disease, Therapy, Timepoint)%>%summarise((n=n_distinct(ID)))
EX_TEM_Filt@meta.data%>%group_by(Disease,Therapy, Timepoint)%>%summarise((n=n_distinct(Sample_ID)))

output_dir<-"~/Immune_project/full_dataset_analysis/plots/DEG/"
Idents(pseudo)<-pseudo$sub_cluster
pseudo@assays$RNA$counts<-round(pseudo@assays$RNA$counts)
# -------------------------------------------------------------------------

# Load DESeq2 library
library(DESeq2)

# Step 1: Extract raw counts from Seurat object
raw_counts <- GetAssayData(pseudo, assay = "RNA", slot = "counts")

# # Step 2: Round the counts to integers
# rounded_counts <- round(raw_counts)

# Step 3: Create DESeq2 dataset
coldata <- data.frame(row.names = colnames(raw_counts), Sample = colnames(raw_counts))
dds <- DESeqDataSetFromMatrix(countData = raw_counts, colData = coldata, design = ~ 1)

# Step 4: Normalize counts using DESeq2
dds <- estimateSizeFactors(dds)

sizeFactors(dds)

normalized_counts <- counts(dds, normalized = TRUE)
pseudo[["RNA"]] <- CreateAssayObject(data = as.matrix(normalized_counts))

# -------------------------------------------------------------------------
library(Seurat)
library(dplyr)
library(tidyr)

memory<-c("SELL","CD27","CD28","MYB")
exhaustion<- c("PDCD1","HAVCR2","LAG3","TIGIT","CTLA4","CD38","TOX2","TOX") 
cytotox1<-c("GZMK","GZMB","GZMH","GNLY","NKG7","PRF1", "IKZF2")
cytotox2<-c("ITGAM","FGFBP2","KLRF1", "TBX21","TYROBP")

# Step 1: Extract normalized counts
counts <- as.data.frame(GetAssayData(pseudo, slot = "data"))  # Use "data" for normalized counts
counts$gene <- rownames(counts)  # Add gene names as a column
# Step 2: Reshape the data to long format
counts_long <- counts %>%
  pivot_longer(
    cols = -c(gene),  # Keep gene and p_val_adj columns
    names_to = "cell",           # Column name for cell IDs
    values_to = "expression"     # Column name for expression values
  )

# Step 3: Add metadata (e.g., sub_cluster) to the counts data
metadata <- pseudo@meta.data
metadata$cell <- rownames(metadata)  # Ensure cell names are a column for merging

counts_long <- counts_long %>%
  left_join(metadata, by = "cell")  # Merge with metadata

Meta<-read_xlsx("~/Immune_project/full_dataset_analysis/Meta_Immune_proj_feb25.xlsx") 
meta<-Meta%>%distinct(Sample_ID, .keep_all=TRUE)
counts_long$Sample_ID<- gsub("-","_",counts_long$Sample_ID) 
counts_long<-left_join(counts_long, meta, by = "Sample_ID")

counts_long<-counts_long%>%mutate(Tissue_Cohort=paste0(Tissue,"_", Cohort))
counts_long$sub_cluster<-gsub("-"," ", counts_long$sub_cluster)
counts_long$sub_cluster<-factor(counts_long$sub_cluster, levels=c("CD8+ GZMK+ TEM","CD8+ DR+ TEM",  "CD8+ Tex"      , "CD8+ GZMB+ TEM" ,  "CD8+ KIR+ TEM" ))

# Esempio di definizione di una palette pastello scuro
dark_pastels <- c(
  "#6F9FD8", # blu pastello scuro
  "#F1B4B3", # rosa pastello scuro
  "#C3DE88", # verde pastello scuro
  "#F9CF87", # giallo pastello scuro
  "#BBA9D6", # viola pastello scuro
  "#EDA1A1"  # un altro rosa/rosso pastello scuro
)
median_iqr <- function(x) {
  qs <- quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  # Rinominiamo i valori come si aspetta stat_summary
  names(qs) <- c("ymin", "y", "ymax")
  return(qs)
}

counts_long_pb<-counts_long%>%filter(Tissue=="PB")
counts_long_bm<-counts_long%>%filter(Tissue=="BM")

# Define a function to create a plot for a given gene set
make_plot_pb <- function(gene_set, plot_title, text_offset) {
  dat <- counts_long_pb %>% filter(gene %in% gene_set)

  # Compute maximum expression per gene for positioning the q-value text
  max_expr <- dat %>% 
    group_by(gene) %>% 
    summarise(max_expr = max(expression))
  dat <- left_join(dat, max_expr, by = "gene")
  
  ggplot(dat, aes(x = sub_cluster, y = expression, color=sub_cluster)) +
    # (Esempio) punti jitter o violin
    geom_jitter(width = 0.2, size = 2) +
    scale_color_manual(values=dark_pastels)+
    # 1) Linea verticale dal Q1 al Q3 (senza “box”)
    stat_summary(
      fun.data = median_iqr,
      geom = "errorbar",      # disegna la “barra di errore” verticale
      width = 0,              # niente stanghette orizzontali in cima e in fondo
      color = "black",
      position = position_dodge(width = 0.9)  # se hai gruppi affiancati
    ) +
    
    # 2) Linea orizzontale alla mediana (senza scatola)
    stat_summary(
      fun.data = function(x) {
        m <- median(x, na.rm = TRUE)
        data.frame(y = m, ymin = m, ymax = m)  # y, ymin, ymax = mediana
      },
      geom = "crossbar",
      width = 0.5,            # lunghezza orizzontale della “linea di mediana”
      fatten = 0,             # toglie l’ispessimento tipico del “bar”
      color = "black",
      position = position_dodge(width = 0.9)
    ) +
    
    theme_classic()+
    facet_wrap(~gene, scales = "free", ncol = 6) +
    theme_classic() +
    theme(
      text=element_text(size=20),
      plot.title = element_text(hjust = 0.5),
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "bottom"
    )+labs(y="Normalized Expression", color="Cluster", title="Peripheral Blood")
  
  
}
make_plot_bm <- function(gene_set, plot_title, text_offset) {
  dat <- counts_long_bm %>% filter(gene %in% gene_set)
  
  # Compute maximum expression per gene for positioning the q-value text
  max_expr <- dat %>% 
    group_by(gene) %>% 
    summarise(max_expr = max(expression))
  dat <- left_join(dat, max_expr, by = "gene")
  
  ggplot(dat, aes(x = sub_cluster, y = expression, color=sub_cluster)) +
    # (Esempio) punti jitter o violin
    geom_jitter(width = 0.2, size = 2) +
    scale_color_manual(values=dark_pastels)+
    # 1) Linea verticale dal Q1 al Q3 (senza “box”)
    stat_summary(
      fun.data = median_iqr,
      geom = "errorbar",      # disegna la “barra di errore” verticale
      width = 0,              # niente stanghette orizzontali in cima e in fondo
      color = "black",
      position = position_dodge(width = 0.9)  # se hai gruppi affiancati
    ) +
    
    # 2) Linea orizzontale alla mediana (senza scatola)
    stat_summary(
      fun.data = function(x) {
        m <- median(x, na.rm = TRUE)
        data.frame(y = m, ymin = m, ymax = m)  # y, ymin, ymax = mediana
      },
      geom = "crossbar",
      width = 0.5,            # lunghezza orizzontale della “linea di mediana”
      fatten = 0,             # toglie l’ispessimento tipico del “bar”
      color = "black",
      position = position_dodge(width = 0.9)
    ) +
    
    theme_classic()+
    facet_wrap(~gene, scales = "free", ncol = 6) +
    theme_classic() +
    theme(
      text=element_text(size=20),
      plot.title = element_text(hjust = 0.5),
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "bottom"
    )+labs(y="Normalized Expression", color="Cluster", title="Bone Marrow")
  
  
}

# Generate plots for each gene set
p1_pb <- make_plot_pb(exhaustion, "Exhaustion Markers", 20)+geom_pwc(label.size = 4, method = "wilcox.test",
                                                               step.increase = 0.08,
                                                               size=0.1, tip.length = 0)
p1_pb<-ggadjust_pvalue(
  p1_pb, p.adjust.method = "BH",
  label = paste("q=","{p.adj.format}"),
  hide.ns = "p"
)
ggsave(filename=paste0(output_dir,"PB_TEM_exhaustion_subcluster_feb25.pdf"), plot=p1_pb, width=15, heigh=15)

# Generate plots for each gene set
p1_bm <- make_plot_bm(exhaustion, "Exhaustion Markers", 20)+geom_pwc(label.size = 4, method = "wilcox.test",
                                                                     step.increase = 0.08,
                                                                     size=0.1, tip.length = 0)
p1_bm<-ggadjust_pvalue(
  p1_bm, p.adjust.method = "BH",
  label = paste("q=","{p.adj.format}"),
  hide.ns = "p"
)
ggsave(filename=paste0(output_dir,"bm_TEM_exhaustion_subcluster_feb25.pdf"), plot=p1_bm, width=15, heigh=15)


# -------------------------------------------------------------------------
#memory markers
# -------------------------------------------------------------------------


# Generate plots for each gene set
p2_pb <- make_plot_pb(memory, "Memory|Stemness Markers", 20)+geom_pwc(label.size = 4, method = "wilcox.test",
                                                                     step.increase = 0.08,
                                                                     size=0.1, tip.length = 0)
p2_pb<-ggadjust_pvalue(
  p2_pb, p.adjust.method = "BH",
  label = paste("q=","{p.adj.format}"),
  hide.ns = "p"
)
ggsave(filename=paste0(output_dir,"PB_TEM_memory_subcluster_feb25.pdf"), plot=p2_pb, width=15, height=7)

# Generate plots for each gene set
p2_bm <- make_plot_bm(memory, "Memory|Stemness Markers", 20)+geom_pwc(label.size = 4, method = "wilcox.test",
                                                                     step.increase = 0.08,
                                                                     size=0.1, tip.length = 0)
p2_bm<-ggadjust_pvalue(
  p2_bm, p.adjust.method = "BH",
  label = paste("q=","{p.adj.format}"),
  hide.ns = "p"
)
ggsave(filename=paste0(output_dir,"bm_TEM_memory_subcluster_feb25.pdf"), plot=p2_bm, width=15, height=7)



# -------------------------------------------------------------------------
#cytotox
# -------------------------------------------------------------------------

cytotox<-c(cytotox1, cytotox2)
# Generate plots for each gene set
p3_pb <- make_plot_pb(cytotox,  20)+geom_pwc(label.size = 4, method = "wilcox.test",
                                                                      step.increase = 0.08,
                                                                      size=0.1, tip.length = 0)
p3_pb<-ggadjust_pvalue(
  p3_pb, p.adjust.method = "BH",
  label = paste("q=","{p.adj.format}"),
  hide.ns = "p"
)
ggsave(filename=paste0(output_dir,"PB_TEM_cytotox_subcluster_feb25.pdf"), plot=p3_pb, width=15, height=10)

# Generate plots for each gene set
p3_bm <- make_plot_bm(cytotox,  20)+geom_pwc(label.size = 4, method = "wilcox.test",
                                                                      step.increase = 0.08,
                                                                      size=0.1, tip.length = 0)
p3_bm<-ggadjust_pvalue(
  p3_bm, p.adjust.method = "BH",
  label = paste("q=","{p.adj.format}"),
  hide.ns = "p"
)
ggsave(filename=paste0(output_dir,"bm_TEM_cytotox_subcluster_feb25.pdf"), plot=p3_bm, width=15, height=10)



