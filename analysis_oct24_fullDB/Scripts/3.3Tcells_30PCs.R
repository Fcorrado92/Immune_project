#Look at T cells clustered with 30 PCs
T_cells<-qread("/mnt/disks/disk/full_dataset_qs/T_cells_jan25_final2.qs")
DimHeatmap(T_cells,reduction = "pca", dims = 1:30, cells = 500, balanced = TRUE)
print(T_cells[["pca"]], dims = 1:30, nfeatures = 10)
T_cells <- RunPCA(T_cells,dims=1:40)
T_cells <- RunUMAP(T_cells, reduction = "harmony", dims=1:30)
T_cells <- FindNeighbors(T_cells, reduction = "harmony", dims=1:30)
T_cells <- FindClusters(T_cells, resolution = c(1, 1.5, 2))
T_cells[['ident']]<-Idents(T_cells)
DimPlot(T_cells, label=T)
qsave(T_cells, "/mnt/disks/cellranger/full_dataset_qs/T_cells_jan25_final2_30pc.res1.qs")
# FindMarkers
output_file <- "/mnt/disks/cellranger/full_dataset_qs/markers/T_cells_markers_PCA30.csv"
markers30<-read_csv("/mnt/disks/cellranger/full_dataset_qs/markers/T_cells_markers_PCA30.csv")
# Controlla se il file esiste
if (file.exists(output_file)) {
  file_content <- read.csv(output_file)
  print(file_content) } else {
    print("Performing Differential Expression for all clusters...")
    markers <- FindAllMarkers(T_cells)
    top_markers <- markers %>%
      group_by(cluster) %>%
      dplyr::filter(pct.1 > 0.3) %>%
      arrange(desc(avg_log2FC)) %>%
      slice_head(n = 20)
    
    print(top_markers)
    write_csv(top_markers, output_file)
  }

T_cells<-qread( "/mnt/disks/cellranger/full_dataset_qs/T_cells_jan25_final2_30pc.res1.qs")
#remove low quality cells
keep<-setdiff(unique(Idents(T_cells)),c(9,16,12,13,18,22,27,24,26,29,30,31))
T_cells_filt<-subset(T_cells, idents=keep)
T_cells_filt<-FindVariableFeatures(object=T_cells_filt, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
T_cells_filt <- RunPCA(T_cells_filt,dims=1:50)
T_cells_filt <- RunUMAP(T_cells_filt, reduction = "harmony", dims=1:30)
T_cells_filt <- FindNeighbors(T_cells_filt, reduction = "harmony", dims=1:30)
T_cells_filt <- FindClusters(T_cells_filt, resolution = c(1))
T_cells_filt[['ident']]<-Idents(T_cells_filt)





meta<-read_excel("~/Immune_project/full_dataset_analysis/Meta_Immune_proj.xlsx")
meta<-distinct(meta, Sample_ID, .keep_all = TRUE)

T_cells_filt@meta.data$new_barcodes<-rownames(T_cells_filt@meta.data)
T_cells_filt@meta.data<-T_cells_filt@meta.data[-c(44:51)]
T_cells_filt@meta.data<-left_join(T_cells_filt@meta.data, meta, by="Sample_ID")
rownames(T_cells_filt@meta.data)<-T_cells_filt@meta.data$new_barcodes
#REMOVE LEN AND SCREENING, KEEP CART TALQ TEC HEALTHY
T_cells_filt<-subset(T_cells_filt, subset=Timepoint%in%c("Post","Pre","Day14","Day30","Day100","Day56","6Months",   
                                                         "Healthy")&Therapy%in%c("TEC","TALQ","CART","Healthy"))

qsave(T_cells_filt,"~/Immune_project/full_dataset_analysis/replicate_ASH/T_cells_filt_PC30.qs")
T_cells_filt<-qread("~/Immune_project/full_dataset_analysis/replicate_ASH/T_cells_filt_PC30.qs")
# pseudobulk DEG with filtered dataset ------------------------------------
#Look at expression of exhaustion markers in cluster sixteen
CD8_TEM<-subset(T_cells_filt, idents=c(17,2,23,20,21,9,3,19,5,4,16))
# 1) Extract current cluster identities into a vector
current_clusters <- Idents(CD8_TEM)

# 2) Replace cluster 16 with "EX", all others with "TEM"
new_clusters <- ifelse(current_clusters == 16, "EX", "TEM")

# 3) Assign these new labels back as the identities
Idents(CD8_TEM) <- new_clusters
CD8_TEM$sub_cluster<-Idents(CD8_TEM)
#create metadata for filtering(sample+cluster)
CD8_TEM@meta.data <- CD8_TEM@meta.data %>%
  mutate(sample_cluster = paste(Sample_ID, sub_cluster, sep = "_"))

#take only pre-treatment                    
EX_TEM_Pre<-subset(CD8_TEM, subset=Timepoint%in%c("Healthy","Pre"))
#look at how many cells x samples in the two clusters remove clusters with <10 cells

df<-as.data.frame(table(EX_TEM_Pre$Sample_ID, EX_TEM_Pre$sub_cluster))
colnames(df)<-c("Sample_ID","sub_cluster","n")
df<-df%>%mutate(sample_cluster=paste0(Sample_ID,"_",sub_cluster))
df<-df%>%group_by(sample_cluster)%>%filter(n>=10)
retain<-unique(df$sample_cluster)

EX_TEM_Pre_Filt<-subset(EX_TEM_Pre, subset=sample_cluster %in%retain)
# pseudobulk the counts based on donor-condition-celltype
pseudo<-AggregateExpression(EX_TEM_Pre_Filt, assays = "RNA",return.seurat = T, 
                            group.by = c("Sample_ID", "sub_cluster"))


# pseudobulk DEG with filtered dataset ------------------------------------
EX_TEM_Pre_Filt@meta.data%>%group_by(Disease, Therapy, Timepoint)%>%summarise((n=n_distinct(ID)))
EX_TEM_Pre_Filt@meta.data%>%group_by(Disease,Therapy, Timepoint)%>%summarise((n=n_distinct(Sample_ID)))

Idents(pseudo)<-pseudo$sub_cluster
pseudo@assays$RNA$counts<-round(pseudo@assays$RNA$counts)
DEG<-FindMarkers(pseudo, ident.1="EX", ident.2="TEM", test.use = "DESeq2")
DEG$gene<-rownames(DEG)
DEG<-DEG%>%mutate(significance=ifelse(p_val_adj<0.05, "q<0.05","q>0.05"))
DEG<-DEG%>%mutate(verse=ifelse(avg_log2FC<0, "downregulated","upregulated"))
sign<-DEG%>%filter(p_val_adj<0.05)

memory<-c("SELL","CD27","CD28","MYB")
exhaustion<- c("HAVCR2","CD38","TIGIT","CTLA4","PDCD1","TOX") 
cytotox1<-c("GZMK","GZMB","GZMH","GNLY","NKG7")
cytotox2<-c("ITGAM","FGFBP2","KLRG1", "KLRC2","TBX21")
show<-c(memory, exhaustion, cytotox1, cytotox2)



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



# Step 1: Extract normalized counts
counts <- as.data.frame(GetAssayData(pseudo, slot = "data"))  # Use "data" for normalized counts
counts$gene <- rownames(counts)  # Add gene names as a column
# Step 2: Reshape the data to long format
counts_long <- counts %>%
  pivot_longer(
    cols = -c(gene),  # Keep gene and p_val_adj columns
    names_to = "sample_cluster",           # Column name for cell IDs
    values_to = "expression"     # Column name for expression values
  )

# Step 3: Add metadata (e.g., sub_cluster) to the counts data
metadata <- pseudo@meta.data
metadata$sample_cluster <- rownames(metadata)  # Ensure cell names are a column for merging

counts_long <- counts_long %>%
  left_join(metadata, by = "sample_cluster")  # Merge with metadata

# Step 4: Filter for selected features (genes)
counts_filtered <- counts_long %>%
  filter(gene %in% exhaustion)  # Keep only the genes in the `show` vector

counts_filtered$sub_cluster<-factor(counts_filtered$sub_cluster, levels=c("TEM","EX"))
p1<-ggplot(counts_filtered, aes(x = gene, y = expression, fill = sub_cluster)) +# Overlay boxplot
  geom_point(
    aes(fill = sub_cluster, stroke = 1,alpha=0.8), color="black",  # Point colors by sub_cluster
    size = 4,  shape = 21,  # Shape 21 allows color fill with black outline
    position = position_jitterdodge(
      dodge.width = 0.9,  # Dodge for different sub_clusters
      jitter.width = 0.2  # Add jitter to points
    )
  )+
  geom_violin(alpha = 0.5) +
  geom_boxplot(
    size = 1, width = 0.5,alpha=0.5 ,outlier.shape = NA,
    position = position_dodge(width = 0.9)
  )  +  # Add jittered points
  labs(
    x = "",
    y = "",
    title = "Exhaustion Markers"
  ) +
  scale_fill_manual(values = c("gold", "orange")) +  # Fill colors for violin
  theme_classic() +
  facet_wrap(~gene, scales = "free", ncol = 6) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right"
  ) +
  guides(fill = FALSE,alpha=FALSE)


p1<-p1+ geom_pwc(label.size = 5, method = "wilcox.test",step.increase = 0.08, tip.length = 0.01, size=0.8)
p1<-ggadjust_pvalue(
  p1, p.adjust.method = "BH",
  label = paste("q=","{p.adj.format}"),
)
p1


# Suppose 'counts_filtered' has columns:
# 1) 'gene'         => Factor/character for gene name
# 2) 'expression'   => Numeric expression value
# 3) 'sub_cluster'  => Factor with two or more groups

# 1) Compute Wilcoxon test per gene
stat.test <- counts_filtered %>%
  group_by(gene) %>%
  wilcox_test(expression ~ sub_cluster) %>%  # e.g. two-group test
  adjust_pvalue(method = "BH")              # Add BH-adjusted p-values

# 2) Optionally rename 'p.adj' column to 'qvalue' for clarity
# stat.test <- stat.test %>%
#   rename(qvalue = p.adj)

# 3) Inspect the results
#    This gives you one row per gene (and sub-cluster comparison),
#    showing n1, n2, statistic, p-value, and q-value, etc.
stat.test


# -------------------------------------------------------------------------
#now only paired
# -------------------------------------------------------------------------
df2<-counts_filtered%>%group_by(Sample_ID)%>%summarise(n=n_distinct(sub_cluster))
df3<-df2%>%filter(n>1)
ids<-unique(df3$Sample_ID)
counts_filtered_paired<-counts_filtered%>%filter(Sample_ID%in%ids)
p1<-ggplot(counts_filtered_paired, aes(x = sub_cluster, y = expression, fill = sub_cluster)) +
  geom_line(aes(group = Sample_ID), color = "black", alpha = 0.6,
            position = position_dodge(width = 0.2))+# Add jittered points
  # Overlay boxplot
  geom_point(
    aes(fill = sub_cluster, stroke = 1,alpha=0.8), color="black",  # Point colors by sub_cluster
    size = 4,  shape = 21,  # Shape 21 allows color fill with black outline
    position = position_jitterdodge(
      dodge.width = 0.9,  # Dodge for different sub_clusters
      jitter.width = 0.2  # Add jitter to points
    )
  )+
  geom_violin(alpha = 0.5) +
  geom_boxplot(
    size = 1, width = 0.5,alpha=0.5 ,outlier.shape = NA,
    position = position_dodge(width = 0.9)
  )  + 
  labs(
    x = "",
    y = "",
    title = "Exhaustion Markers"
  ) +
  scale_fill_manual(values = c("gold", "orange")) +  # Fill colors for violin
  theme_classic() +
  facet_wrap(~gene, scales = "free", ncol = 6) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right"
  ) +
  guides(fill = FALSE,alpha=FALSE)

p1<-p1+ stat_compare_means(size = 3, method = "wilcox.test",label="p.format",paired=TRUE,step.increase = 0.08, tip.length = 0.01)

p1


# 1) Compute Wilcoxon test per gene
stat.test <- counts_filtered_paired %>%
  group_by(gene) %>%
  wilcox_test(expression ~ sub_cluster, paired=TRUE) %>%  # e.g. two-group test
  adjust_pvalue(method = "BH")              # Add BH-adjusted p-values

# 2) Optionally rename 'p.adj' column to 'qvalue' for clarity
stat.test <- stat.test %>%
  rename(qvalue = p.adj)

# 3) Inspect the results
#    This gives you one row per gene (and sub-cluster comparison),
#    showing n1, n2, statistic, p-value, and q-value, etc.
stat.test


# -------------------------------------------------------------------------
#Visualize difference in the percentage of exhausted cells
# -------------------------------------------------------------------------
Idents(T_cells_filt)<-T_cells_filt$RNA_snn_res.1
T_cells_filt$sub_cluster2<-Idents(T_cells_filt)
PBMC<-read_csv("~/Immune_project/full_dataset_analysis/compositional_analysis/SC_PBMC_cell_counts_over500.csv")
PBMC<-PBMC%>%group_by(Sample_ID)%>%summarise(total_PBMC=total)%>%distinct(Sample_ID, .keep_all = TRUE)
counts<-as.data.frame(table(T_cells_filt$Sample_ID, T_cells_filt$sub_cluster2))
colnames(counts)<-c("Sample_ID","sub_cluster","n")
counts<-left_join(counts, PBMC[c("Sample_ID","total_PBMC")], by="Sample_ID")

total<-counts%>% group_by(Sample_ID)%>% summarise(total=sum(n))
counts<-left_join(counts, total, by="Sample_ID")
counts<-counts%>% mutate(Freq= n/total)
counts<-counts%>% mutate(Freq_PBMC= n/total_PBMC)

meta<-read_excel("~/Immune_project/full_dataset_analysis/Meta_Immune_proj.xlsx")
meta<-meta%>%distinct(Sample_ID,.keep_all = TRUE)
counts<-left_join(counts, meta[,c("Sample_ID","BOR","Timepoint", "Tissue","Therapy","Disease","ID", "Cohort")], by="Sample_ID")
counts_filt<-counts%>%filter(Timepoint=="Pre"&Therapy%in%c("TEC")&Tissue%in%c("PB","BM"))
counts_filt <- counts_filt %>%
  mutate(response_disease = case_when(
    BOR %in% c("PR", "MR", "SD", "PD") & Disease == "RRMM" ~ "RRMM_NR",
    Disease == "RRMM" & !(BOR %in% c("PR", "MR", "SD", "PD")) ~ "RRMM_R",
    Disease == "HRSMM" & Therapy != "Len" ~ "HRSMM_R",
    Disease == "Healthy" ~ "Healthy_NA",
    TRUE ~ "SMM_NA"
  ))
counts_filt$response_disease<-factor(counts_filt$response_disease, levels=c("HRSMM_R","RRMM_R","RRMM_NR"))
ex<-counts_filt%>%filter(sub_cluster=="16")

p<-ggplot(ex, aes(x=response_disease, y=Freq_PBMC))+geom_boxplot(outlier.shape = NA)+geom_jitter(width=0.3)+
  facet_wrap(~Tissue ,scales="free")
p<-p+ geom_pwc(label.size = 5, method = "wilcox.test",step.increase = 0.08, tip.length = 0.01, size=0.8)
p<-ggadjust_pvalue(
  p, p.adjust.method = "BH",
  label = paste("p=","{p.format}"),
)
p



#Samples used for annotation of object in ASH
ASH<-qread("~/Immune_project/new_analysis_sep24/qs_files/T_cells_annotated_filtered.qs")
ASH<-RenameIdents(ASH,c('DR+Cd8_Tem'="DR+CD8_TEM", 'aCd4_Tem'="CD4_TEM", "Pro"="Cycling T cells",
                                'Gzmk+Cd8_Tem'="GZMK+CD8_TEM",'Kir+Cd8_Tem'="KIR+CD8_TEM",'Gzmb+Cd8_Tem'="GZMB+CD8_TEM",
                                'Cd8_Naive'="CD8_Naive",'Cd4_Naive'="CD4_Naive",'Cd4_Tem'="CD4_TEM",
                                'Exhausted_Cd8'="Exhausted-like CD8"))
ASH$sub_cluster<-Idents(ASH)
ASH_samples<-unique(ASH$Sample_ID)
ASH1@meta.data%>%group_by(Disease, Therapy, Timepoint)%>%summarise(n=n_distinct(Sample_ID))
# `summarise()` has grouped output by 'Disease', 'Therapy'. You can override using the `.groups` argument.
# # A tibble: 8 × 4
# # Groups:   Disease, Therapy [5]
# Disease Therapy Timepoint     n
# <chr>   <chr>   <chr>     <int>
#   1 Healthy HD      NA           10
# 2 RRMM    TALQ    Post          2
# 3 RRMM    TALQ    Pre           3
# 4 RRMM    TEC     Post          8
# 5 RRMM    TEC     Pre           8
# 6 SMM     Len     Pre           6
# 7 SMM     TEC     Post         11
# 8 SMM     TEC     Pre          19


# -------------------------------------------------------------------------
#select ASH_samples
# -------------------------------------------------------------------------
ASH1<-subset(T_cells_filt, Sample_ID%in% ASH_samples)

ASH1 <- FindVariableFeatures(object=ASH1, selection.method = "vst", nfeatures = 2000, verbose = FALSE) 

ASH1 <- ScaleData(object=ASH1, verbose = FALSE)
print("Running PCA...")
ASH1 <- RunPCA(object=ASH1, npcs = 30, verbose = FALSE)
ElbowPlot(ASH1, ndims = 30)
ASH1 <- RunHarmony(ASH1, group.by.vars = "Pool", plot=TRUE)
print("Clustering corrected dataset...")
ASH1 <- RunUMAP(ASH1, reduction = "harmony", dims=1:30)
ASH1 <- FindNeighbors(ASH1, reduction = "harmony", dims=1:30)
ASH1 <- FindClusters(ASH1, resolution = 1)
ASH1[['ident']]<-Idents(ASH1)
DimPlot(ASH1, label=T)

FeaturePlot(ASH1,c("GZMK","GZMB","CD8A"))

# pseudobulk DEG with filtered dataset ------------------------------------
#Look at expression of exhaustion markers in cluster sixteen
CD8_TEM<-subset(ASH1, idents=c(9,12,17,21,4,13,15,2,3,22,23,16))
# 1) Extract current cluster identities into a vector
current_clusters <- Idents(CD8_TEM)

# 2) Replace cluster 16 with "EX", all others with "TEM"
new_clusters <- ifelse(current_clusters == 17, "EX", "TEM")

# 3) Assign these new labels back as the identities
Idents(CD8_TEM) <- new_clusters
CD8_TEM$sub_cluster<-Idents(CD8_TEM)
#create metadata for filtering(sample+cluster)
CD8_TEM@meta.data <- CD8_TEM@meta.data %>%
  mutate(sample_cluster = paste(Sample_ID, sub_cluster, sep = "_"))

#take only pre-treatment                    
EX_TEM_Pre<-subset(CD8_TEM, subset=Timepoint%in%c("Healthy","Pre"))
#look at how many cells x samples in the two clusters remove clusters with <10 cells

df<-as.data.frame(table(EX_TEM_Pre$Sample_ID, EX_TEM_Pre$sub_cluster))
colnames(df)<-c("Sample_ID","sub_cluster","n")
df<-df%>%mutate(sample_cluster=paste0(Sample_ID,"_",sub_cluster))
df<-df%>%group_by(sample_cluster)%>%filter(n>=10)
retain<-unique(df$sample_cluster)

EX_TEM_Pre_Filt<-subset(EX_TEM_Pre, subset=sample_cluster %in%retain)
# pseudobulk the counts based on donor-condition-celltype
pseudo<-AggregateExpression(EX_TEM_Pre_Filt, assays = "RNA",return.seurat = T, 
                            group.by = c("Sample_ID", "sub_cluster"))


# pseudobulk DEG with filtered dataset ------------------------------------
EX_TEM_Pre_Filt@meta.data%>%group_by(Disease, Therapy, Timepoint)%>%summarise((n=n_distinct(ID)))
EX_TEM_Pre_Filt@meta.data%>%group_by(Disease,Therapy, Timepoint)%>%summarise((n=n_distinct(Sample_ID)))

Idents(pseudo)<-pseudo$sub_cluster
pseudo@assays$RNA$counts<-round(pseudo@assays$RNA$counts)
DEG<-FindMarkers(pseudo, ident.1="EX", ident.2="TEM", test.use = "DESeq2")
DEG$gene<-rownames(DEG)
DEG<-DEG%>%mutate(significance=ifelse(p_val_adj<0.05, "q<0.05","q>0.05"))
DEG<-DEG%>%mutate(verse=ifelse(avg_log2FC<0, "downregulated","upregulated"))
sign<-DEG%>%filter(p_val_adj<0.05)

memory<-c("SELL","CD27","CD28","MYB")
exhaustion<- c("HAVCR2","CD38","TIGIT","CTLA4","PDCD1","TOX") 
cytotox1<-c("GZMK","GZMB","GZMH","GNLY","NKG7")
cytotox2<-c("ITGAM","FGFBP2","KLRG1", "KLRC2","TBX21")
show<-c(memory, exhaustion, cytotox1, cytotox2)



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



# Step 1: Extract normalized counts
counts <- as.data.frame(GetAssayData(pseudo, slot = "data"))  # Use "data" for normalized counts
counts$gene <- rownames(counts)  # Add gene names as a column
# Step 2: Reshape the data to long format
counts_long <- counts %>%
  pivot_longer(
    cols = -c(gene),  # Keep gene and p_val_adj columns
    names_to = "sample_cluster",           # Column name for cell IDs
    values_to = "expression"     # Column name for expression values
  )

# Step 3: Add metadata (e.g., sub_cluster) to the counts data
metadata <- pseudo@meta.data
metadata$sample_cluster <- rownames(metadata)  # Ensure cell names are a column for merging

counts_long <- counts_long %>%
  left_join(metadata, by = "sample_cluster")  # Merge with metadata

# Step 4: Filter for selected features (genes)
counts_filtered <- counts_long %>%
  filter(gene %in% exhaustion)  # Keep only the genes in the `show` vector

counts_filtered$sub_cluster<-factor(counts_filtered$sub_cluster, levels=c("TEM","EX"))
p1<-ggplot(counts_filtered, aes(x = gene, y = expression, fill = sub_cluster)) +# Overlay boxplot
  geom_point(
    aes(fill = sub_cluster, stroke = 1,alpha=0.8), color="black",  # Point colors by sub_cluster
    size = 4,  shape = 21,  # Shape 21 allows color fill with black outline
    position = position_jitterdodge(
      dodge.width = 0.9,  # Dodge for different sub_clusters
      jitter.width = 0.2  # Add jitter to points
    )
  )+
  geom_violin(alpha = 0.5) +
  geom_boxplot(
    size = 1, width = 0.5,alpha=0.5 ,outlier.shape = NA,
    position = position_dodge(width = 0.9)
  )  +  # Add jittered points
  labs(
    x = "",
    y = "",
    title = "Exhaustion Markers"
  ) +
  scale_fill_manual(values = c("gold", "orange")) +  # Fill colors for violin
  theme_classic() +
  facet_wrap(~gene, scales = "free", ncol = 6) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right"
  ) +
  guides(fill = FALSE,alpha=FALSE)


p1<-p1+ geom_pwc(label.size = 5, method = "wilcox.test",step.increase = 0.08, tip.length = 0.01, size=0.8)
p1<-ggadjust_pvalue(
  p1, p.adjust.method = "BH",
  label = paste("q=","{p.adj.format}"),
)
p1


# Suppose 'counts_filtered' has columns:
# 1) 'gene'         => Factor/character for gene name
# 2) 'expression'   => Numeric expression value
# 3) 'sub_cluster'  => Factor with two or more groups

# 1) Compute Wilcoxon test per gene
stat.test <- counts_filtered %>%
  group_by(gene) %>%
  wilcox_test(expression ~ sub_cluster) %>%  # e.g. two-group test
  adjust_pvalue(method = "BH")              # Add BH-adjusted p-values

# 2) Optionally rename 'p.adj' column to 'qvalue' for clarity
# stat.test <- stat.test %>%
#   rename(qvalue = p.adj)

# 3) Inspect the results
#    This gives you one row per gene (and sub-cluster comparison),
#    showing n1, n2, statistic, p-value, and q-value, etc.
stat.test


# -------------------------------------------------------------------------
#now only paired
# -------------------------------------------------------------------------
df2<-counts_filtered%>%group_by(Sample_ID)%>%summarise(n=n_distinct(sub_cluster))
df3<-df2%>%filter(n>1)
ids<-unique(df3$Sample_ID)
counts_filtered_paired<-counts_filtered%>%filter(Sample_ID%in%ids)
p1<-ggplot(counts_filtered_paired, aes(x = sub_cluster, y = expression, fill = sub_cluster)) +
  geom_line(aes(group = Sample_ID), color = "black", alpha = 0.6,
            position = position_dodge(width = 0.2))+# Add jittered points
  # Overlay boxplot
  geom_point(
    aes(fill = sub_cluster, stroke = 1,alpha=0.8), color="black",  # Point colors by sub_cluster
    size = 4,  shape = 21,  # Shape 21 allows color fill with black outline
    position = position_jitterdodge(
      dodge.width = 0.9,  # Dodge for different sub_clusters
      jitter.width = 0.2  # Add jitter to points
    )
  )+
  geom_violin(alpha = 0.5) +
  geom_boxplot(
    size = 1, width = 0.5,alpha=0.5 ,outlier.shape = NA,
    position = position_dodge(width = 0.9)
  )  + 
  labs(
    x = "",
    y = "",
    title = "Exhaustion Markers"
  ) +
  scale_fill_manual(values = c("gold", "orange")) +  # Fill colors for violin
  theme_classic() +
  facet_wrap(~gene, scales = "free", ncol = 6) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right"
  ) +
  guides(fill = FALSE,alpha=FALSE)

p1<-p1+ stat_compare_means(size = 3, method = "wilcox.test",label="p.format",paired=TRUE,step.increase = 0.08, tip.length = 0.01)

p1


# 1) Compute Wilcoxon test per gene
stat.test <- counts_filtered_paired %>%
  group_by(gene) %>%
  wilcox_test(expression ~ sub_cluster, paired=TRUE) %>%  # e.g. two-group test
  adjust_pvalue(method = "BH")              # Add BH-adjusted p-values

# 2) Optionally rename 'p.adj' column to 'qvalue' for clarity
stat.test <- stat.test %>%
  rename(qvalue = p.adj)

# 3) Inspect the results
#    This gives you one row per gene (and sub-cluster comparison),
#    showing n1, n2, statistic, p-value, and q-value, etc.
stat.test


# -------------------------------------------------------------------------
#Visualize difference in the percentage of exhausted cells
# -------------------------------------------------------------------------
Idents(ASH1)<-ASH1$RNA_snn_res.1
ASH1$sub_cluster2<-Idents(ASH1)
PBMC<-read_csv("~/Immune_project/full_dataset_analysis/compositional_analysis/SC_PBMC_cell_counts_over500.csv")
PBMC<-PBMC%>%group_by(Sample_ID)%>%summarise(total_PBMC=total)%>%distinct(Sample_ID, .keep_all = TRUE)
counts<-as.data.frame(table(ASH1$Sample_ID, ASH1$sub_cluster2))
colnames(counts)<-c("Sample_ID","sub_cluster","n")
counts<-left_join(counts, PBMC[c("Sample_ID","total_PBMC")], by="Sample_ID")

total<-counts%>% group_by(Sample_ID)%>% summarise(total=sum(n))
counts<-left_join(counts, total, by="Sample_ID")
counts<-counts%>% mutate(Freq= n/total)
counts<-counts%>% mutate(Freq_PBMC= n/total_PBMC)

meta<-read_excel("~/Immune_project/full_dataset_analysis/Meta_Immune_proj.xlsx")
meta<-meta%>%distinct(Sample_ID,.keep_all = TRUE)
counts<-left_join(counts, meta[,c("Sample_ID","BOR","Timepoint", "Tissue","Therapy","Disease","ID", "Cohort")], by="Sample_ID")
counts_filt<-counts%>%filter(Timepoint=="Pre"&Therapy%in%c("TEC")&Tissue%in%c("PB","BM"))
counts_filt <- counts_filt %>%
  mutate(response_disease = case_when(
    BOR %in% c("PR", "MR", "SD", "PD") & Disease == "RRMM" ~ "RRMM_NR",
    Disease == "RRMM" & !(BOR %in% c("PR", "MR", "SD", "PD")) ~ "RRMM_R",
    Disease == "HRSMM" & Therapy != "Len" ~ "HRSMM_R",
    Disease == "Healthy" ~ "Healthy_NA",
    TRUE ~ "SMM_NA"
  ))
counts_filt$response_disease<-factor(counts_filt$response_disease, levels=c("HRSMM_R","RRMM_R","RRMM_NR"))
ex<-counts_filt%>%filter(sub_cluster=="17")

p<-ggplot(ex, aes(x=response_disease, y=Freq_PBMC))+geom_boxplot(outlier.shape = NA)+geom_jitter(width=0.3)+
  facet_wrap(~Tissue ,scales="free")
p<-p+ geom_pwc(label.size = 5, method = "wilcox.test",step.increase = 0.08, tip.length = 0.01, size=0.8)
p<-ggadjust_pvalue(
  p, p.adjust.method = "BH",
  label = paste("p=","{p.format}"),
)
p

qsave(ASH1,"~/Immune_project/full_dataset_analysis/replicate_ASH/ASH_subset.qs")





# -------------------------------------------------------------------------
#Look at PCs of my T cells
# -------------------------------------------------------------------------
T_cells<-qread("/mnt/disks/disk/full_dataset_qs/T_cells_jan25_final2.qs")
Idents(T_cells)<-T_cells$sub_cluster
clusters_to_remove <- grep("^dbl", unique(Idents(T_cells)),value=T)
clusters_to_remove<-unique(c(clusters_to_remove, "Plt", "lowquality", "Lowq", "lowq","Plt"))
Idents(T_cells)<-T_cells$sub_cluster
keep<-setdiff(unique(Idents(T_cells)), clusters_to_remove)
T_cells<- subset(T_cells, idents = keep)
meta<-read_excel("~/Immune_project/full_dataset_analysis/Meta_Immune_proj_feb25.xlsx")
meta<-distinct(meta, Sample_ID, .keep_all = TRUE)

T_cells@meta.data$new_barcodes<-rownames(T_cells@meta.data)
T_cells@meta.data<-T_cells@meta.data[-c(44:51)]
T_cells@meta.data<-left_join(T_cells@meta.data, meta, by="Sample_ID")
rownames(T_cells@meta.data)<-T_cells@meta.data$new_barcodes

# Carica le librerie necessarie
library(Seurat)
library(harmony)

# Definisci la classe base
base_class <- "DFCI_RRMM_PB_TEC"

# Crea il subset iniziale con la classe base
sobj_base <- subset(sobj, subset = class == base_class)

obj<-subset(T_cells, Sample_ID%in%ASH_samples)
# Funzione per eseguire il pipeline di analisi
run_pipeline <- function(obj) {
obj<-T_cells
  # Identificazione delle feature variabili (2000 features, metodo "vst")
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
  
  # Scaling dei dati
  obj <- ScaleData(obj)
  
  # Calcolo della PCA (30 componenti principali)
  obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
  
  # Integrazione con Harmony, raggruppando per "sample_id" usando i primi 30 PC
  obj <- RunHarmony(obj, group.by.vars = "Library", dims.use = 1:30)
  
  # Calcolo dell'UMAP partendo dalla riduzione Harmony
  obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30)
  
  # Calcolo dei vicini (neighbors) e clustering
  obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
  obj <- FindClusters(obj, resolution = 1.2)  # Il parametro resolution può essere adattato
  qsave(obj, "~/Immune_project/full_dataset_analysis/sobj_ASH_cohort.qs")
  
  return(obj)
}

# Esegui il pipeline sul subset base
sobj_base <- run_pipeline(sobj_base)

# Salva l'oggetto e le coordinate UMAP per il subset base
saveRDS(sobj_base, file = paste0("sobj_", base_class, ".rds"))
umap_base <- Embeddings(sobj_base, "umap")
write.csv(umap_base, file = paste0("umap_", base_class, ".csv"))

# Recupera tutte le classi presenti nell'oggetto e seleziona quelle diverse dalla classe base
all_classes <- unique(sobj@meta.data$class)
other_classes <- setdiff(all_classes, base_class)

# Loop: aggiungi le altre classi una alla volta alla classe base e analizza il clustering
for(cl in other_classes) {
  message("Analizzando aggiunta della classe: ", cl)
  
  # Crea il subset con la classe base + la nuova classe
  current_subset <- subset(sobj, subset = class %in% c(base_class, cl))
  
  # Esegui il pipeline di analisi sul subset corrente
  current_subset <- run_pipeline(current_subset)
  
  # Salva l'oggetto risultante
  filename_obj <- paste0("sobj_", base_class, "_", cl, ".rds")
  saveRDS(current_subset, file = filename_obj)
  
  # Salva le coordinate UMAP in formato CSV
  umap_data <- Embeddings(current_subset, "umap")
  filename_umap <- paste0("umap_", base_class, "_", cl, ".csv")
  write.csv(umap_data, file = filename_umap)
}



# -------------------------------------------------------------------------
#load ASH+ Calgary
obj<-qread("/mnt/disks/cellranger/replicate_ASH_qs/sobj_ASH_DFCI_RRMM_PB_CART.qs")
obj[['ident']]<-Idents(obj)
DimPlot(obj, label=T) +scale_color_manual(values = colors) 
FeaturePlot(obj,c("ciltacabtagene"))

# pseudobulk DEG with filtered dataset ------------------------------------
#Look at expression of exhaustion markers in cluster sixteen
CD8_TEM<-subset(obj, idents=c(15,5,3,12,17,0,4))
DimPlot(CD8_TEM)
# 1) Extract current cluster identities into a vector
current_clusters <- Idents(CD8_TEM)

# 2) Replace cluster 16 with "EX", all others with "TEM"
new_clusters <- ifelse(current_clusters == 15, "EX", "TEM")

# 3) Assign these new labels back as the identities
Idents(CD8_TEM) <- new_clusters
CD8_TEM$sub_cluster<-Idents(CD8_TEM)
#create metadata for filtering(sample+cluster)
CD8_TEM@meta.data <- CD8_TEM@meta.data %>%
  mutate(sample_cluster = paste(Sample_ID, sub_cluster, sep = "_"))

#take only pre-treatment                    
EX_TEM_Pre<-subset(CD8_TEM, subset=Timepoint%in%c("Healthy","Pre"))
#look at how many cells x samples in the two clusters remove clusters with <10 cells

df<-as.data.frame(table(EX_TEM_Pre$Sample_ID, EX_TEM_Pre$sub_cluster))
colnames(df)<-c("Sample_ID","sub_cluster","n")
df<-df%>%mutate(sample_cluster=paste0(Sample_ID,"_",sub_cluster))
df<-df%>%group_by(sample_cluster)%>%filter(n>=10)
retain<-unique(df$sample_cluster)

EX_TEM_Pre_Filt<-subset(EX_TEM_Pre, subset=sample_cluster %in%retain)
# pseudobulk the counts based on donor-condition-celltype
pseudo<-AggregateExpression(EX_TEM_Pre_Filt, assays = "RNA",return.seurat = T, 
                            group.by = c("Sample_ID", "sub_cluster"))


# pseudobulk DEG with filtered dataset ------------------------------------
EX_TEM_Pre_Filt@meta.data%>%group_by(Disease, Therapy, Timepoint,Tissue)%>%summarise((n=n_distinct(ID)))
EX_TEM_Pre_Filt@meta.data%>%group_by(Disease,Therapy, Timepoint)%>%summarise((n=n_distinct(Sample_ID)))

Idents(pseudo)<-pseudo$sub_cluster
pseudo@assays$RNA$counts<-round(pseudo@assays$RNA$counts)
# DEG<-FindMarkers(pseudo, ident.1="EX", ident.2="TEM", test.use = "DESeq2")
# DEG$gene<-rownames(DEG)
# DEG<-DEG%>%mutate(significance=ifelse(p_val_adj<0.05, "q<0.05","q>0.05"))
# DEG<-DEG%>%mutate(verse=ifelse(avg_log2FC<0, "downregulated","upregulated"))
# sign<-DEG%>%filter(p_val_adj<0.05)
# 
# memory<-c("SELL","CD27","CD28","MYB")
exhaustion<- c("HAVCR2","CD38","TIGIT","CTLA4","PDCD1","TOX") 
# cytotox1<-c("GZMK","GZMB","GZMH","GNLY","NKG7")
# cytotox2<-c("ITGAM","FGFBP2","KLRG1", "KLRC2","TBX21")
# show<-c(memory, exhaustion, cytotox1, cytotox2)



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



# Step 1: Extract normalized counts
counts <- as.data.frame(GetAssayData(pseudo, slot = "data"))  # Use "data" for normalized counts
counts$gene <- rownames(counts)  # Add gene names as a column
# Step 2: Reshape the data to long format
counts_long <- counts %>%
  pivot_longer(
    cols = -c(gene),  # Keep gene and p_val_adj columns
    names_to = "sample_cluster",           # Column name for cell IDs
    values_to = "expression"     # Column name for expression values
  )

# Step 3: Add metadata (e.g., sub_cluster) to the counts data
metadata <- pseudo@meta.data
metadata$sample_cluster <- rownames(metadata)  # Ensure cell names are a column for merging

counts_long <- counts_long %>%
  left_join(metadata, by = "sample_cluster")  # Merge with metadata

# Step 4: Filter for selected features (genes)
counts_filtered <- counts_long %>%
  filter(gene %in% exhaustion)  # Keep only the genes in the `show` vector

counts_filtered$sub_cluster<-factor(counts_filtered$sub_cluster, levels=c("TEM","EX"))
p1<-ggplot(counts_filtered, aes(x = gene, y = expression, fill = sub_cluster)) +# Overlay boxplot
  geom_point(
    aes(fill = sub_cluster, stroke = 1,alpha=0.8), color="black",  # Point colors by sub_cluster
    size = 4,  shape = 21,  # Shape 21 allows color fill with black outline
    position = position_jitterdodge(
      dodge.width = 0.9,  # Dodge for different sub_clusters
      jitter.width = 0.2  # Add jitter to points
    )
  )+
  geom_violin(alpha = 0.5) +
  geom_boxplot(
    size = 1, width = 0.5,alpha=0.5 ,outlier.shape = NA,
    position = position_dodge(width = 0.9)
  )  +  # Add jittered points
  labs(
    x = "",
    y = "",
    title = "Exhaustion Markers"
  ) +
  scale_fill_manual(values = c("gold", "orange")) +  # Fill colors for violin
  theme_classic() +
  facet_wrap(~gene, scales = "free", ncol = 6) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right"
  ) +
  guides(fill = FALSE,alpha=FALSE)


p1<-p1+ geom_pwc(label.size = 5, method = "wilcox.test",step.increase = 0.08, tip.length = 0.01, size=0.8)
p1<-ggadjust_pvalue(
  p1, p.adjust.method = "BH",
  label = paste("q=","{p.adj.format}"),
)
p1


# Suppose 'counts_filtered' has columns:
# 1) 'gene'         => Factor/character for gene name
# 2) 'expression'   => Numeric expression value
# 3) 'sub_cluster'  => Factor with two or more groups

# 1) Compute Wilcoxon test per gene
stat.test <- counts_filtered %>%
  group_by(gene) %>%
  wilcox_test(expression ~ sub_cluster) %>%  # e.g. two-group test
  adjust_pvalue(method = "BH")              # Add BH-adjusted p-values

# 2) Optionally rename 'p.adj' column to 'qvalue' for clarity
# stat.test <- stat.test %>%
#   rename(qvalue = p.adj)

# 3) Inspect the results
#    This gives you one row per gene (and sub-cluster comparison),
#    showing n1, n2, statistic, p-value, and q-value, etc.
stat.test


# -------------------------------------------------------------------------
#now only paired
# -------------------------------------------------------------------------
df2<-counts_filtered%>%group_by(Sample_ID)%>%summarise(n=n_distinct(sub_cluster))
df3<-df2%>%filter(n>1)
ids<-unique(df3$Sample_ID)
counts_filtered_paired<-counts_filtered%>%filter(Sample_ID%in%ids)
p1<-ggplot(counts_filtered_paired, aes(x = sub_cluster, y = expression, fill = sub_cluster)) +
  geom_line(aes(group = Sample_ID), color = "black", alpha = 0.6,
            position = position_dodge(width = 0.2))+# Add jittered points
  # Overlay boxplot
  geom_point(
    aes(fill = sub_cluster, stroke = 1,alpha=0.8), color="black",  # Point colors by sub_cluster
    size = 4,  shape = 21,  # Shape 21 allows color fill with black outline
    position = position_jitterdodge(
      dodge.width = 0.9,  # Dodge for different sub_clusters
      jitter.width = 0.2  # Add jitter to points
    )
  )+
  geom_violin(alpha = 0.5) +
  geom_boxplot(
    size = 1, width = 0.5,alpha=0.5 ,outlier.shape = NA,
    position = position_dodge(width = 0.9)
  )  + 
  labs(
    x = "",
    y = "",
    title = "Exhaustion Markers"
  ) +
  scale_fill_manual(values = c("gold", "orange")) +  # Fill colors for violin
  theme_classic() +
  facet_wrap(~gene, scales = "free", ncol = 6) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right"
  ) +
  guides(fill = FALSE,alpha=FALSE)

p1<-p1+ stat_compare_means(size = 3, method = "wilcox.test",label="p.format",paired=TRUE,step.increase = 0.08, tip.length = 0.01)

p1


# 1) Compute Wilcoxon test per gene
stat.test <- counts_filtered_paired %>%
  group_by(gene) %>%
  wilcox_test(expression ~ sub_cluster, paired=TRUE) %>%  # e.g. two-group test
  adjust_pvalue(method = "BH")              # Add BH-adjusted p-values

# # 2) Optionally rename 'p.adj' column to 'qvalue' for clarity
# stat.test <- stat.test %>%
#   rename(qvalue = p.adj)

# 3) Inspect the results
#    This gives you one row per gene (and sub-cluster comparison),
#    showing n1, n2, statistic, p-value, and q-value, etc.
stat.test


# -------------------------------------------------------------------------
#Visualize difference in the percentage of exhausted cells
# -------------------------------------------------------------------------
obj$sub_cluster2<-Idents(obj)
PBMC<-read_csv("~/Immune_project/full_dataset_analysis/compositional_analysis/SC_PBMC_cell_counts_over500.csv")
PBMC<-PBMC%>%group_by(Sample_ID)%>%summarise(total_PBMC=total)%>%distinct(Sample_ID, .keep_all = TRUE)
counts<-as.data.frame(table(obj$Sample_ID, obj$sub_cluster2))
colnames(counts)<-c("Sample_ID","sub_cluster","n")
counts<-left_join(counts, PBMC[c("Sample_ID","total_PBMC")], by="Sample_ID")

total<-counts%>% group_by(Sample_ID)%>% summarise(total=sum(n))
counts<-left_join(counts, total, by="Sample_ID")
counts<-counts%>% mutate(Freq= n/total)
counts<-counts%>% mutate(Freq_PBMC= n/total_PBMC)

meta<-read_excel("~/Immune_project/full_dataset_analysis/Meta_Immune_proj_feb25.xlsx")
meta$BOR[meta$ID=="617885"]<-"PR"

meta<-meta%>%distinct(Sample_ID,.keep_all = TRUE)
counts<-left_join(counts, meta[,c("Sample_ID","BOR","Timepoint", "Tissue","Therapy","Disease","ID", "Cohort")], by="Sample_ID")
counts_filt<-counts%>%filter(Timepoint=="Pre"&Therapy%in%c("TEC")&Tissue%in%c("PB","BM"))
counts_filt <- counts_filt %>%
  mutate(response_disease = case_when(
    BOR %in% c("PR", "MR", "SD", "PD") & Disease == "RRMM" ~ "RRMM_NR",
    Disease == "RRMM" & !(BOR %in% c("PR", "MR", "SD", "PD")) ~ "RRMM_R",
    Disease == "HRSMM" & Therapy != "Len" ~ "HRSMM_R",
    Disease == "Healthy" ~ "Healthy_NA",
    TRUE ~ "SMM_NA"
  ))
counts_filt$response_disease<-factor(counts_filt$response_disease, levels=c("HRSMM_R","RRMM_R","RRMM_NR"))
ex<-counts_filt%>%filter(sub_cluster=="15")

p<-ggplot(ex, aes(x=response_disease, y=Freq_PBMC))+geom_boxplot(outlier.shape = NA)+geom_jitter(width=0.3)+
  facet_wrap(~Tissue ,scales="free")
p<-p+ geom_pwc(label.size = 5, method = "wilcox.test",step.increase = 0.08, tip.length = 0.01, size=0.8)
p<-ggadjust_pvalue(
  p, p.adjust.method = "BH",
  label = paste("p=","{p.format}"),
)
p


















PBMC<-read_csv("~/Immune_project/full_dataset_analysis/compositional_analysis/SC_PBMC_cell_counts_over500.csv")
PBMC<-PBMC%>%group_by(Sample_ID)%>%summarise(total_PBMC=total)%>%distinct(Sample_ID, .keep_all = TRUE)
counts<-as.data.frame(table(obj$Sample_ID, obj$sub_cluster2))
colnames(counts)<-c("Sample_ID","sub_cluster","n")
counts<-left_join(counts, PBMC[c("Sample_ID","total_PBMC")], by="Sample_ID")

total<-counts%>% group_by(Sample_ID)%>% summarise(total=sum(n))
counts<-left_join(counts, total, by="Sample_ID")
counts<-counts%>% mutate(Freq= n/total)
counts<-counts%>% mutate(Freq_PBMC= n/total_PBMC)

counts<-counts%>%filter(!total<50)
meta<-read_excel("~/Immune_project/full_dataset_analysis/Meta_Immune_proj.xlsx")
meta<-meta%>%distinct(Sample_ID,.keep_all = TRUE)
counts<-left_join(counts, meta[,c("Sample_ID","BOR","Timepoint", "Tissue","Therapy","Disease","ID", "Cohort")], by="Sample_ID")
counts_filt<-counts%>%filter(Timepoint=="Pre"&Therapy%in%c("TEC")&Tissue=="PB")
counts_filt <- counts_filt %>%
  mutate(response_disease = case_when(
    BOR %in% c("PR", "MR", "SD", "PD") & Disease == "RRMM" ~ "RRMM_NR",
    Disease == "RRMM" & !(BOR %in% c("PR", "MR", "SD", "PD")) ~ "RRMM_R",
    Disease == "HRSMM" & Therapy != "Len" ~ "HRSMM_R",
    Disease == "Healthy" ~ "Healthy_NA",
    TRUE ~ "SMM_NA"
  ))
counts_filt$response_disease<-factor(counts_filt$response_disease, levels=c("HRSMM_R","RRMM_R","RRMM_NR"))
ex<-counts_filt%>%filter(sub_cluster=="17")
ex<-ex%>%filter(!Sample_ID=="XL1261ePB")

p<-ggplot(ex, aes(x=response_disease, y=Freq_PBMC))+geom_boxplot(outlier.shape = NA)+geom_jitter(width=0.3)+
  facet_wrap(~sub_cluster, scales="free")+geom_text(aes(label=Sample_ID))
p<-p+ geom_pwc(label.size = 5, method = "wilcox.test",step.increase = 0.08, tip.length = 0.01, size=0.8)
p<-ggadjust_pvalue(
  p, p.adjust.method = "BH",
  label = paste("p=","{p.format}"),
)
p




















#try with 30 PCs

