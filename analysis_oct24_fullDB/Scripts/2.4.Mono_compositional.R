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
Mono<-qread("/mnt/disks/cellranger/full_dataset_qs/Mono_annotated.qs")
clusters_to_remove <- grep("^dbl", unique(Idents(Mono)),value=TRUE)
clusters_to_remove<-unique(c(clusters_to_remove, "Plt", "lowquality", "Lowq", "lowq","Plt"))
Idents(Mono)<-Mono$sub_cluster
keep<-setdiff(unique(Idents(Mono)), clusters_to_remove)
Mono<- subset(Mono, idents = keep)
meta<-read_excel("~/Immune_project/full_dataset_analysis/Meta_Immune_proj.xlsx")
meta<-distinct(meta, Sample_ID, .keep_all = TRUE)

Mono@meta.data$new_barcodes<-rownames(Mono@meta.data)
Mono@meta.data<-Mono@meta.data[-c(44:51)]
Mono@meta.data<-left_join(Mono@meta.data, meta, by="Sample_ID")
rownames(Mono@meta.data)<-Mono@meta.data$new_barcodes
unique(Mono$Therapy)
#REMOVE LEN AND SCREENING, KEEP CART TALQ TEC HEALTHY
Mono_filt<-subset(Mono, subset=Timepoint%in%c("Post","Pre","Day14","Day30","Day100","Day56","6Months",   
                                                    "Healthy")&Therapy%in%c("TEC","TALQ","CART","Healthy"))

Mono_filt$sub_cluster<-Idents(Mono_filt)


# -------------------------------------------------------------------------
ncells<-Mono_filt@meta.data%>%group_by(Disease,Sample_ID,Tissue, Timepoint ,Cohort, Therapy)%>%
  summarise(n=n())

N_cells <- ggplot(ncells, aes(x = n)) +
  geom_histogram(bins = 70, aes(y = ..density..), fill = "lightblue", color = "black", alpha = 0.7) + # Histogram with density scaling
  geom_density(color = "red", size = 1) + # Overlay density line
  theme_classic() +
  scale_x_continuous(breaks = c(500, 1000, 2000, 5000, 10000)) +
  labs(y = "Density", x = "Number of Mono", title = "Distribution of number of Mono per sample") +
  geom_vline(xintercept = 100, linetype = "dotted", size = 1)+
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1))

output_dir_mono<-"~/Immune_project/full_dataset_analysis/plots/Mono/"
dir.create(output_dir_mono)
ggsave(paste0(output_dir_mono,"Ncells_baseline.pdf"), plot=N_cells, width=7, height=4)


# -------------------------------------------------------------------------
counts<-as.data.frame(table(Mono_filt$Sample_ID, Mono_filt$sub_cluster))
colnames(counts)<-c("Sample_ID","sub_cluster","n")
total<-counts%>% group_by(Sample_ID)%>% summarise(total=sum(n))
counts<-left_join(counts, total, by="Sample_ID")
counts<-counts%>% mutate(Freq= n/total)
meta<-meta%>%distinct(Sample_ID,.keep_all = TRUE)
counts<-left_join(counts, meta[,c("Sample_ID","BOR","Timepoint", "Tissue","Therapy","Disease","ID", "Cohort")], by="Sample_ID")
counts<-counts%>%filter(total>50)
write_csv(counts, "~/Immune_project/full_dataset_analysis/compositional_analysis/SC_Mono_cell_counts_over50.csv")

#add levels
Mono_counts<-read_csv("~/Immune_project/full_dataset_analysis/compositional_analysis/SC_Mono_cell_counts_over50.csv")
Mono_counts$Disease[Mono_counts$Disease=="SMM"]<-"HRSMM"
Mono_counts$Disease<-factor(Mono_counts$Disease, levels=c("Healthy","HRSMM","RRMM"))
Mono_counts$BOR<-factor(Mono_counts$BOR, levels=c("Healthy","NA","CR","VGPR","PR","MR","SD","PD"))
#assign response_disease
Mono_counts <- Mono_counts %>%
  mutate(response_disease = case_when(
    BOR %in% c("PR", "MR","SD","PD") & Disease == "RRMM" ~ "RRMM_NR",            # RRMM Non-Responder
    Disease == "RRMM" & !BOR %in% c("PR", "MR","SD","PD") ~ "RRMM_R",            # RRMM Responder
    Disease == "HRSMM" & !Therapy== "Len" ~ "HRSMM_R",                     # SMM Responder
    Disease == "Healthy" ~ "Healthy_NA",                              # Healthy: NA response
    TRUE ~ "SMM_NA"                                                    # Default: SMM Non-Applicable
  ))
Mono_counts%>%group_by(Disease, Therapy, Timepoint, Tissue, response_disease, BOR)%>%summarise(n=n_distinct(Sample_ID))
# -------------------------------------------------------------------------
#by filtering "Pre" I am selecting BL MM samples, counts here is Mono_counts >1000
MM_BL_Mono_counts<-Mono_counts%>%filter(Timepoint%in%c("Pre")&Tissue=="PB")
#check number of samples x disease x timepoint
MM_BL_Mono_counts%>%group_by(Disease, Timepoint)%>%summarise(n=n_distinct(Sample_ID))
# Disease Timepoint     n
# 1 HRSMM   Pre          30
# 2 RRMM    Pre          39



#do wilcox.test between RRMM and SMM for each celltype, store.pvalues, correct with BH
celltypes<-unique(MM_BL_Mono_counts$sub_cluster)
results <- list()  
for (i in seq_along(celltypes)) {  
  sub <- MM_BL_Mono_counts %>% filter(sub_cluster == celltypes[i])  
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

MM_BL_Mono_counts$Disease<-factor(MM_BL_Mono_counts$Disease, levels=c("HRSMM","RRMM"))

MM_BL_Mono_counts%>%group_by(Disease, Therapy, Tissue, Timepoint)%>%summarise(n=n_distinct(Sample_ID))
#extract mean values for each subtype across disease groups
means<-MM_BL_Mono_counts%>% group_by(Disease, sub_cluster)%>% summarise(mean_value=mean(Freq))
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
pdf(paste0(output_dir_mono, "logFC_mean_diff_Mono.pdf"), width=18, height = 14)
logFCmeans
dev.off()


# -------------------------------------------------------------------------

cols<-c("darkgreen","darkblue","darkred")
#do wilcox.test between RRMM and SMM for each celltype, store.pvalues, correct with BH
df<-means_wide%>%filter(p_adj<0.15)
subsets<-unique(df$sub_cluster)
PBMC_sign<-Mono_counts%>%filter(Timepoint%in%c("Pre","Healthy") & Tissue=="PB" & sub_cluster%in%subsets)


PBMC_sign$sub_cluster<-gsub("_","",PBMC_sign$sub_cluster)
PBMC_sign$Disease<-factor(PBMC_sign$Disease, levels=c("Healthy","HRSMM","RRMM"))
p<-ggplot(PBMC_sign, aes(x=Disease, y=Freq, fill=Disease))+
  geom_violin(alpha=0.5)+geom_boxplot(size=3,alpha=0.6, outlier.shape = NA)+
  geom_jitter(size=10,width=0.1,aes(fill=Disease, stroke=2),color="black",shape=21, alpha=0.5)+
  labs(y="Percentage of Monocytes")+
  theme_classic()+guides(color=FALSE)+
  scale_fill_manual(values = cols)+
  facet_wrap(~sub_cluster, scales="free", ncol=4)+
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

ggsave(paste0(output_dir_mono,"boxplots_sign_Mono_cluster_disease_baseline.pdf"), 
       plot=p, width=15,height=13)

# -------------------------------------------------------------------------
# Volcano TAM vs Other
# -------------------------------------------------------------------------
#wilcoxon
#tale Mono in pts and HD at baseline
Mono_PRE<-subset(Mono_filt, subset=Timepoint%in%c("Pre") & Disease%in%c("HRSMM","RRMM")&Tissue=="PB" 
                 &Therapy%in%c("TEC","TALQ","CART"))
markers<-FindMarkers(Mono_PRE, ident.1 = "C1Q+Macrophage")
sign<-markers%>%filter(p_val_adj<0.05)
markers$gene<-rownames(markers)
markers<-markers%>%mutate(significance=ifelse(p_val_adj<0.05, "q<0.05","q>0.05"))
markers<-markers%>%mutate(verse=ifelse(avg_log2FC<=0, "downregulated","upregulated"))
sign<-markers%>%filter(p_val_adj<0.05)

DEG <- markers
DEG <- DEG %>% 
  mutate(color = case_when(
    significance == "q<0.05" & verse == "upregulated" ~ "orange",
    significance == "q<0.05" & verse == "downregulated" ~ "gold",
    TRUE ~ "gray"  # Non-significant genes are gray
  ))

# Create the plot
tam<-c("C1QC","C1QA","C1QB","MRC1","MSR1","HLA-DRA","HLA-DRB1",  "MAF","MAFB", "CEBPB")

# Define a maximum y-value, e.g., 300
max_y <- 300

# Replace Inf with max_y
DEG <- DEG %>%
  mutate(log_p_val = ifelse(is.infinite(-log10(p_val)), max_y, -log10(p_val)))

p1 <- ggplot(DEG, aes(x = avg_log2FC, y = log_p_val, color = color)) +
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
  )
pdf(paste0(output_dir_mono, "/DEG_TAM.pdf"), width = 13, height = 12)
p1
dev.off()



#rename tams and not tams
Mono_PRE@meta.data<-Mono_PRE@meta.data%>%mutate(TAM=ifelse(sub_cluster=="C1Q+Macrophage","C1Q+Macrophage","Other Monocytes"))
# pseudobulk the counts based on donor-condition-celltype
pseudo_mono <- AggregateExpression(Mono_PRE, assays = "RNA", return.seurat = TRUE, group.by = c("Sample_ID", "TAM"))
pseudo_mono@meta.data<-pseudo_mono@meta.data%>%mutate(cluster_sample=paste0(Sample_ID,"_",TAM))
df<-as.data.frame(table(Mono_PRE$Sample_ID, Mono_PRE$TAM))
colnames(df)<-c("Sample_ID", "sub_cluster", "n")
df<-df%>%mutate(cluster_sample=paste0(Sample_ID,"_",sub_cluster))
View(df)
df<-df%>%filter(n>10)
cs<-unique(df$cluster_sample)

#filter only samples with more than 10 cells
pseudo_mono_filt<-subset(pseudo_mono, subset=cluster_sample%in%cs)
# Step 1: Extract raw counts from Seurat object
raw_counts <- GetAssayData(pseudo_mono_filt, assay = "RNA", slot = "counts")
# # Step 2: Round the counts to integers
raw_counts <- round(raw_counts)

# Step 3: Create DESeq2 dataset
coldata <- data.frame(row.names = colnames(raw_counts), Sample = colnames(raw_counts))
dds <- DESeqDataSetFromMatrix(countData = raw_counts, colData = coldata, design = ~ 1)

# Step 4: Normalize counts using DESeq2
dds <- estimateSizeFactors(dds)

sizeFactors(dds)

normalized_counts <- counts(dds, normalized = TRUE)
pseudo_mono_filt[["RNA"]] <- CreateAssayObject(data = as.matrix(normalized_counts))

# # -------------------------------------------------------------------------
# Step 1: Extract normalized counts
counts <- as.data.frame(GetAssayData(pseudo_mono_filt, slot = "data"))  # Use "data" for normalized counts
counts$gene <- rownames(counts)  # Add gene names as a column
# Step 2: Reshape the data to long format
counts_long <- counts %>%
  pivot_longer(
    cols = -c(gene),  # Keep gene and p_val_adj columns
    names_to = "cell",           # Column name for cell IDs
    values_to = "expression"     # Column name for expression values
  )

# Step 3: Add metadata (e.g., sub_cluster) to the counts data
metadata <- pseudo_mono_filt@meta.data
metadata$cell <- rownames(metadata)  # Ensure cell names are a column for merging

counts_long <- counts_long %>%
  left_join(metadata, by = "cell")  # Merge with metadata

# Step 4: Filter for selected features (genes)
counts_filtered <- counts_long %>%
  filter(gene %in% tam)  # Keep only the genes in the `show` vector


p1<-ggplot(counts_filtered, aes(x = gene, y = expression, fill = TAM)) +# Overlay boxplot
  geom_point(
    aes(fill = TAM, stroke = 1,alpha=0.8), color="black",  # Point colors by sub_cluster
    size = 10,  shape = 21,  # Shape 21 allows color fill with black outline
    position = position_jitterdodge(
      dodge.width = 0.9,  # Dodge for different sub_clusters
      jitter.width = 0.2  # Add jitter to points
    )
  )+
  scale_fill_manual(values = c("gold","orange"))+
  geom_violin(alpha = 0.5) +
  geom_boxplot(
    size = 1, alpha=0.5 ,outlier.shape = NA,
    position = position_dodge(width = 0.9)
  )  +  # Add jittered points
  labs(
    x = "",
    y = "Normalized Expression",
    title = ""
  ) +
  theme_classic() +
  facet_wrap(~gene, scales = "free", ncol = 5) +
  theme(plot.margin = unit(c(0.5, 0, 0.2, 0.1), "cm"),title=element_text(size=40),
        plot.title = element_text(hjust=0.5, face = "bold"), text=element_text(size=50,color = "black"),
        axis.text.x =element_blank(),legend.position = "bottom",
        axis.ticks.x = element_blank(),
        axis.text.y =element_text(size=40,color = "black"), axis.title.y =element_text(size=50,color = "black"))+
  guides(alpha=FALSE)

p1<-p1+ geom_pwc(label.size = 8, method = "wilcox.test",step.increase = 0.08, tip.length = 0.01, size=1.2)
p1<-ggadjust_pvalue(
  p1, p.adjust.method = "BH",
  label = paste("q=","{p.adj.format}"),
  hide.ns = "p"
)

ggsave(paste0(output_dir_mono,"TAM_markers.pdf"), plot=p1, width=22, height=18)

#ABSOLUTE VALUES
Mono_PRE2<-subset(Mono_filt, subset=Timepoint%in%c("Pre","Healthy")& Tissue=="PB" 
                 &Therapy%in%c("TEC","TALQ","CART","Healthy")&sub_cluster=="C1Q+Macrophage")

df2<-as.data.frame(table(Mono_PRE2$Sample_ID))
colnames(df2)<-c("Sample_ID","n")
df2<-left_join(df2, meta, by="Sample_ID")

# 5) Convert PatientID to factor with sorted levels
df2$Disease <- factor(df2$Disease, levels = c("Healthy","HRSMM","RRMM"))
df2 <- df2 %>%
  arrange(Disease, n)
df2$Sample_ID <- factor(df2$Sample_ID, levels = unique(df2$Sample_ID))

# 6) Plot 1: bar plot of median expression
p <- ggplot(df2, aes(x = Sample_ID, y = n, fill = Disease)) +
  geom_col() +
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(
    y = paste0("Number of PB C1Q+ Macro"),
    x = "Samples",
    fill = "Diagnosis"
  ) +
  theme(legend.title = element_blank(), text = element_text(size = 50), legend.position = "bottom") +
  scale_fill_manual(
    values = c(
      "Healthy"          = "forestgreen",
      "HRSMM"           = "darkblue",
      "RRMM"          = "darkred"    )
  )

ggsave(paste0(output_dir_mono,"TAM_abs_n.pdf"), plot=p, width=18, height=15)
