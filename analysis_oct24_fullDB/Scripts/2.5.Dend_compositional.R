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
DendCells<-qread("/mnt/disks/cellranger/full_dataset_qs/DC_annotated.qs")
clusters_to_remove <- grep("^dbl", unique(Idents(DendCells)),value=TRUE)
clusters_to_remove<-unique(c(clusters_to_remove, "Plt", "lowquality", "Lowq", "lowq","Plt"))
unique(Idents(DendCells))
keep<-setdiff(unique(Idents(DendCells)), clusters_to_remove)
DendCells<- subset(DendCells, idents = keep)
meta<-read_excel("~/Immune_project/full_dataset_analysis/Meta_Immune_proj.xlsx")
meta<-distinct(meta, Sample_ID, .keep_all = TRUE)

DendCells@meta.data$new_barcodes<-rownames(DendCells@meta.data)
DendCells@meta.data<-DendCells@meta.data[-c(44:51)]
DendCells@meta.data<-left_join(DendCells@meta.data, meta, by="Sample_ID")
rownames(DendCells@meta.data)<-DendCells@meta.data$new_barcodes
unique(DendCells$Therapy)
#REMOVE LEN AND SCREENING, KEEP CART TALQ TEC HEALTHY
DendCells_filt<-subset(DendCells, subset=Timepoint%in%c("Post","Pre","Day14","Day30","Day100","Day56","6Months",   
                                              "Healthy")&Therapy%in%c("TEC","TALQ","CART","Healthy"))

DendCells_filt$sub_cluster<-Idents(DendCells_filt)


# -------------------------------------------------------------------------
ncells<-DendCells_filt@meta.data%>%group_by(Disease,Sample_ID,Tissue, Timepoint ,Cohort, Therapy)%>%
  summarise(n=n())

N_cells <- ggplot(ncells, aes(x = n)) +
  geom_histogram(bins = 70, aes(y = ..density..), fill = "lightblue", color = "black", alpha = 0.7) + # Histogram with density scaling
  geom_density(color = "red", size = 1) + # Overlay density line
  theme_classic() +
  scale_x_continuous(breaks = c(0,10,50,500, 1000, 2000, 5000, 10000)) +
  labs(y = "Density", x = "Number of DendCells", title = "Distribution of number of DendCells per sample") +
  geom_vline(xintercept = 100, linetype = "dotted", size = 1)+
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1))

output_dir_DendCells<-"~/Immune_project/full_dataset_analysis/plots/DendCells/"
dir.create(output_dir_DendCells)
ggsave(paste0(output_dir_DendCells,"Ncells_baseline.pdf"), plot=N_cells, width=7, height=4)


# -------------------------------------------------------------------------
counts<-as.data.frame(table(DendCells_filt$Sample_ID, DendCells_filt$sub_cluster))
colnames(counts)<-c("Sample_ID","sub_cluster","n")
total<-counts%>% group_by(Sample_ID)%>% summarise(total=sum(n))
counts<-left_join(counts, total, by="Sample_ID")
counts<-counts%>% mutate(Freq= n/total)
meta<-meta%>%distinct(Sample_ID,.keep_all = TRUE)
counts<-left_join(counts, meta[,c("Sample_ID","BOR","Timepoint", "Tissue","Therapy","Disease","ID", "Cohort")], by="Sample_ID")
counts<-counts%>%filter(total>50)
write_csv(counts, "~/Immune_project/full_dataset_analysis/compositional_analysis/SC_DendCells_cell_counts_over50.csv")

#add levels
DendCells_counts<-read_csv("~/Immune_project/full_dataset_analysis/compositional_analysis/SC_DendCells_cell_counts_over50.csv")
DendCells_counts$Disease[DendCells_counts$Disease=="SMM"]<-"HRSMM"
DendCells_counts$Disease<-factor(DendCells_counts$Disease, levels=c("Healthy","HRSMM","RRMM"))
DendCells_counts$BOR<-factor(DendCells_counts$BOR, levels=c("Healthy","NA","CR","VGPR","PR","MR","SD","PD"))
#assign response_disease
DendCells_counts <- DendCells_counts %>%
  mutate(response_disease = case_when(
    BOR %in% c("PR", "MR","SD","PD") & Disease == "RRMM" ~ "RRMM_NR",            # RRMM Non-Responder
    Disease == "RRMM" & !BOR %in% c("PR", "MR","SD","PD") ~ "RRMM_R",            # RRMM Responder
    Disease == "HRSMM" & !Therapy== "Len" ~ "HRSMM_R",                     # SMM Responder
    Disease == "Healthy" ~ "Healthy_NA",                              # Healthy: NA response
    TRUE ~ "SMM_NA"                                                    # Default: SMM Non-Applicable
  ))
DendCells_counts%>%group_by(Disease, Therapy, Timepoint, Tissue, response_disease, BOR)%>%summarise(n=n_distinct(Sample_ID))
# -------------------------------------------------------------------------
#by filtering "Pre" I am selecting BL MM samples, counts here is DendCells_counts >1000
MM_BL_DendCells_counts<-DendCells_counts%>%filter(Timepoint%in%c("Pre")&Tissue=="PB")
#check number of samples x disease x timepoint
MM_BL_DendCells_counts%>%group_by(Disease, Timepoint)%>%summarise(n=n_distinct(Sample_ID))
# Disease Timepoint     n
# 1 HRSMM   Pre          30
# 2 RRMM    Pre          39



#do wilcox.test between RRMM and SMM for each celltype, store.pvalues, correct with BH
celltypes<-unique(MM_BL_DendCells_counts$sub_cluster)
results <- list()  
for (i in seq_along(celltypes)) {  
  sub <- MM_BL_DendCells_counts %>% filter(sub_cluster == celltypes[i])  
  test <- wilcox.test(Freq ~ Disease, data = sub)  
  results[[i]] <- test$p.value  
}
names(results) <- celltypes

adjusted_p_values <- as.data.frame(p.adjust(results, method = "BH"))
adjusted_p_values$sub_cluster <- rownames(adjusted_p_values)
colnames(adjusted_p_values)<-c("p_adj","sub_cluster")
results <- data.frame(
  sub_cluster = names(results),
  pval = unlist(results))
adjusted_p_values<-left_join(adjusted_p_values, results, by="sub_cluster")


# -------------------------------------------------------------------------

MM_BL_DendCells_counts$Disease<-factor(MM_BL_DendCells_counts$Disease, levels=c("HRSMM","RRMM"))

MM_BL_DendCells_counts%>%group_by(Disease, Therapy, Tissue, Timepoint)%>%summarise(n=n_distinct(Sample_ID))
#extract mean values for each subtype across disease groups
means<-MM_BL_DendCells_counts%>% group_by(Disease, sub_cluster)%>% summarise(mean_value=mean(Freq))
means_wide<-means%>% pivot_wider(names_from = Disease, values_from = mean_value)
means_wide<-means_wide%>%mutate(log2FC=log2(RRMM/HRSMM))
means_wide<-left_join(means_wide, adjusted_p_values, by="sub_cluster")
means_wide<-means_wide%>%mutate(enrichment= ifelse(log2FC>0,"RRMM","HRSMM"))
means_wide<-means_wide%>%mutate(significant= ifelse(p_adj<0.1,"q<0.1","q>0.1"))
#remove underscore for the graph but keep means_wide for further analysis
means_wide2<-means_wide
means_wide2$sub_cluster <- gsub("_", "", means_wide2$sub_cluster)


# -------------------------------------------------------------------------
logFCmeans <- ggplot(means_wide2, aes(x = log2FC, y = -log10(pval))) +
  geom_point(aes(shape = significant, 
                 alpha = ifelse(pval < 0.1, 1, 0.5)), size = 10) +  # Make points with pval<0.1 transparent
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "grey", linewidth=2) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey",linewidth=2) +
  scale_shape_manual(values = c("q<0.1" = 17, "q>0.1" = 19)) +  # Define shapes for significance
  scale_color_manual(values = c( "purple", "orange",  "darkblue", "darkgreen","red","gold"), 
                     name = "Lineage") +  # Set custom colors for major_cluster and rename the legend title to "Cell Type"
  geom_text_repel(aes(label = ifelse(pval < 0.05, sub_cluster, "")), 
                  size = 15, max.overlaps = Inf,force = 10,              # Increase force to repel labels
                  box.padding = 1,      # Add padding around boxes
                  point.padding = 1,
                  show.legend = FALSE) +  # Increase the text size
  labs(x = expression(Log[2]~Fold~Change), shape = "Significance",
       y = expression(-log[10](p~value)),
       fill = "Lineage" ) +  # Set legend titles
  theme_classic()+theme(
    text=element_text(size=40),
    legend.position = "bottom",
    # Position the legend at the top center of the plot
    legend.box = "horizontal",  # Arrange legends horizontally
    legend.spacing.x = unit(0.5, "cm"),  # Add space between legend items
    panel.background = element_blank(),  # Remove background from panel
    plot.background = element_blank(),  # Remove background from the plot
    legend.background = element_rect(fill = "white", color = "grey"),  # Add grey border around the legend box
    legend.key = element_blank()  # Remove the background and border from legend keys
  ) +
  guides(alpha = "none") +  # Remove the alpha legend
  annotate("text", x = -0.05, y = 0.1, label = "HRSMM n=30", hjust = 1, size = 14) +  # Left text
  annotate("text", x = 0.05, y = 0.1, label = "RRMM n=39", hjust = 0, size = 14)+
  xlim(c(-0.8,1.5))
pdf(paste0(output_dir_DendCells, "logFC_mean_diff_DendCells.pdf"), width=18, height = 14)
logFCmeans
dev.off()


# -------------------------------------------------------------------------

cols<-c("darkgreen","darkblue","darkred")
#do wilcox.test between RRMM and SMM for each celltype, store.pvalues, correct with BH
df<-means_wide%>%filter(p_adj<1)
subsets<-unique(df$sub_cluster)
PBMC_sign<-DendCells_counts%>%filter(Timepoint%in%c("Pre","Healthy") & Tissue=="PB" & sub_cluster%in%subsets)


PBMC_sign$sub_cluster<-gsub("_","",PBMC_sign$sub_cluster)
PBMC_sign$Disease<-factor(PBMC_sign$Disease, levels=c("Healthy","HRSMM","RRMM"))
p<-ggplot(PBMC_sign, aes(x=Disease, y=Freq, fill=Disease))+
  geom_violin(alpha=0.5)+geom_boxplot(size=3,alpha=0.6, outlier.shape = NA)+
  geom_jitter(size=10,width=0.1,aes(fill=Disease, stroke=2),color="black",shape=21, alpha=0.5)+
  labs(y="Percentage of DendCellscytes")+
  theme_classic()+guides(color=FALSE)+
  scale_fill_manual(values = cols)+
  facet_wrap(~sub_cluster, scales="free", ncol=6)+
  theme(plot.margin = unit(c(0.5, 0, 0.2, 0.1), "cm"),title=element_text(size=40),
        plot.title = element_text(hjust=0.5, face = "bold"), text=element_text(size=50,color = "black"),
        axis.text.x =element_blank(),legend.position = "bottom",
        axis.ticks.x = element_blank(),
        axis.text.y =element_text(size=40,color = "black"), axis.title.y =element_text(size=50,color = "black"))+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) 

p<-p+ geom_pwc(label.size = 12, method = "wilcox.test",step.increase = 0.08, tip.length = 0.01, size=1.2)
p<-ggadjust_pvalue(
  p, p.adjust.method = "BH",
  label = paste("q=","{p.adj.format}"),
  hide.ns = "p"
)

ggsave(paste0(output_dir_DendCells,"boxplots_sign_DendCells_cluster_disease_baseline.pdf"), 
       plot=p, width=30,height=20)

# -------------------------------------------------------------------------
#Look at DEG C1Q+ vs remaing clusters
# -------------------------------------------------------------------------
#pseudobulk
#tale DendCells in pts and HD at baseline
DendCells_PRE<-subset(DendCells_filt, subset=Timepoint%in%c("Healthy","Pre") & Tissue=="PB" &Therapy%in%c("TEC","TALQ","CART","Healthy"))

DendCells_PRE$sub_cluster <- factor(DendCells_PRE$sub_cluster)
DendCells_PRE$sub_cluster <- droplevels(DendCells_PRE$sub_cluster)
df<-as.data.frame(table(DendCells_PRE$Sample_ID, DendCells_PRE$sub_cluster))
colnames(df)<-c("Sample_ID", "sub_cluster", "n")
df<-df%>%mutate(cluster_sample=paste0(Sample_ID,"_",sub_cluster))
df<-df%>%filter(n>10)
df<-left_join(df, meta, by="Sample_ID")
cs<-unique(df$cluster_sample)

#take only clusters with more than 10 cells
DendCells_PRE@meta.data<-DendCells_PRE@meta.data%>%mutate(cluster_sample=paste0(Sample_ID,"_",sub_cluster))
DendCells_PRE_Filt<-subset(DendCells_PRE, subset=cluster_sample%in% cs)
# -------------------------------------------------------------------------
#pseudobulk
# -------------------------------------------------------------------------
# pseudobulk the counts based on donor-condition-celltype
pseudo_DendCells <- AggregateExpression(DendCells_PRE_Filt, assays = "RNA", return.seurat = TRUE, group.by = c("Sample_ID", "sub_cluster"))
# each 'cell' is a donor-condition-celltype pseudobulk profile
tail(Cells(pseudo_DendCells))

output_dir<-"/mnt/disks/disk/full_dataset_qs/DEG/DendCells/"
dir.create(output_dir)
pseudo_DendCells@assays$RNA$counts<-round(pseudo_DendCells@assays$RNA$counts)
Idents(pseudo_DendCells) <- "sub_cluster"
bulk.de <- FindMarkers(object = pseudo_DendCells, 
                       ident.1 = "C1Q+Macrophage", 
                       test.use = "DESeq2")
write.csv(bulk.de, paste0(output_dir,"PB_BL_RRMM_vs_HRSMM_C1Q+.csv"))

bulk.de$gene<-rownames(bulk.de)
bulk.de<-bulk.de%>%mutate(significance=ifelse(p_val_adj<0.05, "q<0.05","q>0.05"))
bulk.de<-bulk.de%>%mutate(verse=ifelse(avg_log2FC<0, "downregulated","upregulated"))
sign<-bulk.de%>%filter(p_val_adj<0.05)


# -------------------------------------------------------------------------
# However I am a bit concerned that I am comparing very unbalanced groups of cells in terms of disease and transcriptional differences between
# C1Q and remaining subset may instead reflect widespread differences between disease group
# restrict only to RRMM

# df<-left_join(df, meta, by="Sample_ID")
# df%>%group_by(Disease, sub_cluster)%>%summarise(n=n_distinct(Sample_ID))
# `summarise()` has grouped output by 'Disease'. You can override using the `.groups` argument.
# # A tibble: 18 × 3
# # Groups:   Disease [3]
# Disease sub_cluster        n
# <chr>   <fct>          <int>
#   1 HRSMM   IL1B+DendCells         30
# 2 HRSMM   CD14+DendCells         30
# 3 HRSMM   CD16+DendCells         30
# 4 HRSMM   IFN+DendCells          27
# 5 HRSMM   DR+DendCells           28
# 6 HRSMM   C1Q+Macrophage     6
# 7 Healthy IL1B+DendCells         10
# 8 Healthy CD14+DendCells         10
# 9 Healthy CD16+DendCells         10
# 10 Healthy IFN+DendCells           8
# 11 Healthy DR+DendCells           10
# 12 Healthy C1Q+Macrophage     2
# 13 RRMM    IL1B+DendCells         37
# 14 RRMM    CD14+DendCells         39
# 15 RRMM    CD16+DendCells         39
# 16 RRMM    IFN+DendCells          36
# 17 RRMM    DR+DendCells           37
# 18 RRMM    C1Q+Macrophage    21
# -------------------------------------------------------------------------

DendCells_PRE_Filt2<-subset(DendCells_PRE_Filt, subset=Disease=="RRMM")
# pseudobulk the counts based on donor-condition-celltype
pseudo_DendCells <- AggregateExpression(DendCells_PRE_Filt2, assays = "RNA", return.seurat = TRUE, group.by = c("Sample_ID", "sub_cluster"))
# each 'cell' is a donor-condition-celltype pseudobulk profile
tail(Cells(pseudo_DendCells))

pseudo_DendCells@assays$RNA$counts<-round(pseudo_DendCells@assays$RNA$counts)
Idents(pseudo_DendCells) <- "sub_cluster"
bulk.de <- FindMarkers(object = pseudo_DendCells, 
                       ident.1 = "C1Q+Macrophage", 
                       test.use = "DESeq2")
write.csv(bulk.de, paste0(output_dir,"PB_BL_RRMM_vs_HRSMM_C1Q+_ONLY_RR.csv"))

bulk.de$gene<-rownames(bulk.de)
bulk.de<-bulk.de%>%mutate(significance=ifelse(p_val_adj<0.05, "q<0.05","q>0.05"))
bulk.de<-bulk.de%>%mutate(verse=ifelse(avg_log2FC<0, "downregulated","upregulated"))
sign<-bulk.de%>%filter(p_val_adj<0.05)




# -------------------------------------------------------------------------
#FGSEA
# -------------------------------------------------------------------------
# Define mapping and org
mapping <- "org.Hs.eg.db"
org <- "human"
gene_sets <- msigdbr(species = org, category = "H") %>% dplyr::select(gs_name, entrez_gene, gene_symbol)
pathways <- split(gene_sets$entrez_gene, gene_sets$gs_name)

# Ensure pathways are named correctly and have unique genes
pathways <- lapply(pathways, unique)

res<-bulk.de
res <- as.data.frame(res)
# Map SYMBOL to ENTREZID
res$entrez <- mapIds(
  org.Hs.eg.db, 
  keys = rownames(res), 
  keytype = "SYMBOL",
  column = "ENTREZID",
  multiVals = "first"  # Use "list" if you prefer to keep all mappings
)

# Calculate ranks
ranks <- res$avg_log2FC * -log10(res$p_val_adj)
names(ranks) <- res$entrez

# Remove NAs and duplicates
ranks <- na.omit(ranks)
ranks <- ranks[!is.na(names(ranks))]
ranks <- ranks[!duplicated(names(ranks))]

# Ensure all ENTREZIDs are characters
ranks <- setNames(as.numeric(ranks), as.character(names(ranks)))

# Run fgsea
fgseaRes <- fgsea(
  pathways = pathways, 
  stats = ranks, 
  minSize = 10, 
  maxSize = 500  )


# -------------------------------------------------------------------------
# Volcano TAM vs Other
# -------------------------------------------------------------------------
DEG <- bulk.de
DEG <- DEG %>% 
  mutate(color = case_when(
    significance == "q<0.05" & verse == "upregulated" ~ "orange",
    significance == "q<0.05" & verse == "downregulated" ~ "gold",
    TRUE ~ "gray"  # Non-significant genes are gray
  ))

# Create the plot

p1 <- ggplot(DEG, aes(x = avg_log2FC, y = -log10(p_val), color = color)) +
  geom_point(size = 8, alpha = 0.5) +
  scale_color_identity() +  # Use colors from the 'color' column without modifying the scale
  labs(x = expression(Log[2]~Fold~Change), shape = "Significance",
       y = expression(-log[10](p~value)),
       title = "" ) +  # Set legend titles
  theme_classic() +
  geom_text_repel(
    data = DEG %>% filter(gene %in% tam),  # Only label genes in 'show'
    aes(label = gene),
    color = "black",
    size = 12,
    force=10,
    box.padding = 0.7,           # Increase padding around labels
    point.padding = 0.7,         # Increase padding around points
    min.segment.length = 0,      # Force connectors to display, even if short
    segment.color = "gray50",    # Connector line color
    segment.size = 0.6 ,
    fontface="bold"
  ) +
  theme(
    text = element_text(size = 50),
    plot.title = element_text(hjust = 0.5),
    axis.title.y = element_text(margin = margin(r = -5))
  )+annotate("text", x = -1, y = 0.2, label = "C1Q+Macrophages", hjust = 1, size = 15, fontface="bold") +  # Left text
  annotate("text", x = 0.5, y = 0.2, label = "DendCellscytes", hjust = 0, size =15 ,fontface="bold")  # Right text


p1

pdf(paste0(output_dir_DendCells, "/DEG_TAM.pdf"), width = 13, height = 12)
p1
dev.off()




# DendCells_PRE_Filt2<-subset(DendCells_PRE_Filt, subset=Disease=="RRMM")
# DendCells_PRE_Filt2@meta.data<-DendCells_PRE_Filt2@meta.data%>%mutate(TAM=ifelse(sub_cluster=="C1Q+Macrophage","C1Q+Macrophage","Other DendCellscytes"))
# # pseudobulk the counts based on donor-condition-celltype
# pseudo_DendCells <- AggregateExpression(DendCells_PRE_Filt2, assays = "RNA", return.seurat = TRUE, group.by = c("Sample_ID", "TAM"))
# # each 'cell' is a donor-condition-celltype pseudobulk profile
# tail(Cells(pseudo_DendCells))
# 
# tam<-c("C1QC","C1QA","C1QB","APOE","MRC1")
# # Load DESeq2 library
# library(DESeq2)
# 
# # Step 1: Extract raw counts from Seurat object
# raw_counts <- GetAssayData(pseudo_DendCells, assay = "RNA", slot = "counts")
# 
# # # Step 2: Round the counts to integers
# raw_counts <- round(raw_counts)
# 
# # Step 3: Create DESeq2 dataset
# coldata <- data.frame(row.names = colnames(raw_counts), Sample = colnames(raw_counts))
# dds <- DESeqDataSetFromMatrix(countData = raw_counts, colData = coldata, design = ~ 1)
# 
# # Step 4: Normalize counts using DESeq2
# dds <- estimateSizeFactors(dds)
# 
# sizeFactors(dds)
# 
# normalized_counts <- counts(dds, normalized = TRUE)
# pseudo_DendCells[["RNA"]] <- CreateAssayObject(data = as.matrix(normalized_counts))
# 
# 
# # -------------------------------------------------------------------------
# # Step 1: Extract normalized counts
# counts <- as.data.frame(GetAssayData(pseudo_DendCells, slot = "data"))  # Use "data" for normalized counts
# counts$gene <- rownames(counts)  # Add gene names as a column
# counts<-left_join(counts, sign[c("gene","p_val_adj")])
# # Step 2: Reshape the data to long format
# counts_long <- counts %>%
#   pivot_longer(
#     cols = -c(gene, p_val_adj),  # Keep gene and p_val_adj columns
#     names_to = "cell",           # Column name for cell IDs
#     values_to = "expression"     # Column name for expression values
#   )
# 
# # Step 3: Add metadata (e.g., sub_cluster) to the counts data
# metadata <- pseudo_DendCells@meta.data
# metadata$cell <- rownames(metadata)  # Ensure cell names are a column for merging
# 
# counts_long <- counts_long %>%
#   left_join(metadata, by = "cell")  # Merge with metadata
# 
# # Step 4: Filter for selected features (genes)
# counts_filtered <- counts_long %>%
#   filter(gene %in% tam)  # Keep only the genes in the `show` vector
# 
# 
# p1<-ggplot(counts_filtered, aes(x = gene, y = expression, fill = TAM)) +# Overlay boxplot
#   geom_point(
#     aes(fill = TAM, stroke = 1,alpha=0.8), color="black",  # Point colors by sub_cluster
#     size = 10,  shape = 21,  # Shape 21 allows color fill with black outline
#     position = position_jitterdodge(
#       dodge.width = 0.9,  # Dodge for different sub_clusters
#       jitter.width = 0.2  # Add jitter to points
#     )
#   )+
#   scale_fill_manual(values = c("gold","orange"))+
#   geom_violin(alpha = 0.5) +
#   geom_boxplot(
#     size = 1, alpha=0.5 ,outlier.shape = NA,
#     position = position_dodge(width = 0.9)
#   )  +  # Add jittered points
#   labs(
#     x = "",
#     y = "Normalized Expression",
#     title = ""
#   ) +
#   theme_classic() +
#   facet_wrap(~gene, scales = "free", ncol = 6) +
#   theme(plot.margin = unit(c(0.5, 0, 0.2, 0.1), "cm"),title=element_text(size=40),
#         plot.title = element_text(hjust=0.5, face = "bold"), text=element_text(size=50,color = "black"),
#         axis.text.x =element_blank(),legend.position = "bottom",
#         axis.ticks.x = element_blank(),
#         axis.text.y =element_text(size=40,color = "black"), axis.title.y =element_text(size=50,color = "black"))+
#   guides(alpha=FALSE)
# 
# # # Calculate the maximum expression for each gene (facet)
# # max_expr_per_gene <- counts_filtered %>%
# #   group_by(gene) %>%
# #   summarise(max_expr = max(expression))
# # 
# # # Merge the max_expr with your data
# # counts_filtered_with_max_expr <- counts_filtered %>%
# #   left_join(max_expr_per_gene, by = "gene")
# # 
# # # Plot with adjusted y position for geom_text
# # p1_with_q <- p1 + 
# #   # Add q-value text with formatting
# #   geom_text(
# #     data = counts_filtered_with_max_expr, color = "black", size = 6, 
# #     aes(
# #       x = gene,  # Add x aesthetic for gene
# #       y = max_expr + 50,  # Use the max expression for each gene
# #       label = ifelse(p_val_adj < 0.001, "q<0.001", paste("q =", formatC(p_val_adj, format = "f", digits = 3))),
# #       hjust = 0.5
# #     ), 
# #     inherit.aes = FALSE  # Prevent inheritance of aesthetics from the main plot
# #   )
# # 
# # # Display the plot
# # p1_with_q
# 
# ggsave(paste0(output_dir_DendCells,"TAM_markers.pdf"), plot=p1, width=18, height=15)








