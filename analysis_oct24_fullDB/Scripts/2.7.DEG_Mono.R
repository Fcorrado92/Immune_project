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
library(SCOPfunctions)
source("~/common_scripts/MAST_with_random_effect.R")
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
Mono_filt<-subset(Mono, subset=Timepoint%in%c("Pre")&Therapy%in%c("TEC","TALQ","CART"))
Mono_filt$sub_cluster<-Idents(Mono_filt)
output_dir<-"/mnt/disks/disk/full_dataset_qs/DEG/Pseudo_bulk/Mono/"
dir.create(output_dir)
clusters<-unique(Mono_filt$sub_cluster)

#pseudobulk and volcano PB
Mono_filt$sub_cluster <- factor(Mono_filt$sub_cluster)
Mono_filt$sub_cluster <- droplevels(Mono_filt$sub_cluster)
df<-as.data.frame(table(Mono_filt$Sample_ID, Mono_filt$sub_cluster))
colnames(df)<-c("Sample_ID", "sub_cluster", "n")
df<-df%>%mutate(cluster_sample=paste0(Sample_ID,"_",sub_cluster))
df<-df%>%filter(n>10)
df<-left_join(df,meta, by="Sample_ID")
pb<-df%>%filter(Tissue=="PB")
table(pb$sub_cluster, pb$Disease)
cs<-unique(df$cluster_sample) 
# -------------------------------------------------------------------------
#pseudobulk
# -------------------------------------------------------------------------
# pseudobulk the counts based on donor-condition-celltype
#For Mono only PB
Mono_filt<-subset(Mono_filt, Tissue=="PB")
pseudo_mono <- AggregateExpression(Mono_filt, assays = "RNA", return.seurat = TRUE, 
                                  group.by = c("Sample_ID","sub_cluster", "Disease"))
#remove those observation with less than 10 cells
pseudo_mono$cluster_sample <- paste(pseudo_mono$Sample_ID,pseudo_mono$sub_cluster, sep = "_")
pseudo_mono<-subset(pseudo_mono, cluster_sample%in% cs)
pseudo_mono$cluster_sample<-gsub("-","_", pseudo_mono$cluster_sample)

#sepecify method
met="Pseudo_bulk"
#rename sobj
sobj<-pseudo_mono
#specifiy Tissue
# Tis=unique(sobj$Tissue)
Cell_Type="Mono"
output_dir<-paste0("/mnt/disks/disk/full_dataset_qs/DEG/",met,"/", Cell_Type, "/")
clusters<-unique(sobj$sub_cluster)
sobj@assays$RNA$counts<-round(sobj@assays$RNA$counts)
for (j in seq_along(clusters)){
    sub<-subset(sobj, sub_cluster==clusters[[j]])
    Idents(sub) <- "Disease"
    bulk.de <- FindMarkers(object = sub, 
                           ident.1 ="RRMM", 
                           ident.2 = "HRSMM",
                           test.use = "DESeq2")
    write.csv(bulk.de, paste0(output_dir,clusters[[j]],"_BL_RRMM_vs_HRSMM.csv"))
  }
# -------------------------------------------------------------------------
#volcano PSEUDOBULK
# -------------------------------------------------------------------------
output_dir<-"/mnt/disks/disk/full_dataset_qs/DEG/Pseudo_bulk/Mono/"

# Initialize an empty list to store the results from each cluster
res_list <- vector("list", length(clusters))

# Loop over the clusters
for (i in seq_along(clusters)) {
    # Read the CSV file for the current cluster
    tmp <- read_csv(paste0(output_dir, clusters[[i]], "_BL_RRMM_vs_HRSMM.csv"))
    
    # Create a symbol column based on the row names (if row names exist)
    # Note: If tmp does not have row names, you might need to adjust this step.
    tmp$symbol <- tmp$...1
    
    # Add the current cluster as a new column
    tmp$cluster <- clusters[[i]]
    # Remove rows with NA in the p_val_adj column
    tmp <- tmp[!is.na(tmp$p_val_adj), ]
    
    # Order the data by p_val_adj
    tmp <- tmp[order(abs(tmp$p_val_adj)), ]
    
    
    # Convert to data.frame (if needed) and store in the list
    res_list[[paste0(clusters[[i]])]]  <- data.frame(tmp)
  }

# Combine the list of data.frames into one data.frame
res <- do.call(rbind, res_list)

res$cluster<-gsub("-","", res$cluster)

# Step A: Filter for only significant genes
sig_res <- res %>%
  filter(p_val_adj < 0.05)

# # Step B: Identify genes that are significant in >=2 distinct tissues for the same cluster
# overlap_genes <- sig_res %>%
#   group_by(cluster, symbol) %>%
#   summarise(
#     tissues_sig = n_distinct(Tissue),
#     .groups = "drop"
#   ) %>%
#   filter(tissues_sig >= 2)
# 
# Step C: Left join to set overlap=TRUE for those genes
# res <- res %>%
#   left_join(
#     overlap_genes %>% 
#       dplyr::select(cluster, symbol) %>% 
#       mutate(overlap = TRUE),
#     by = c("cluster", "symbol")
#   )
# 
# # Step D: Fill NA with FALSE
# res$overlap[is.na(res$overlap)] <- FALSE



  
  # Compute the number of significant genes (p_val_adj < 0.05) per cluster
  counts <- res%>%
    group_by(cluster) %>%
    summarise(n_signif = sum(p_val_adj < 0.05)) %>%
    ungroup()
  
  # Create a named vector of counts for use in the labeller
  counts_vec <- setNames(counts$n_signif, counts$cluster)
  
  # Define a custom labeller function that adds the count on top of each facet
  facet_labeller <- function(variable, value) {
    # value is a vector of cluster names; we use counts_vec to get the count for each
    paste0(value, "\n(n=", counts_vec[value], ")")
  }
  
  
  t<-res%>%group_by(cluster)%>%
    filter(p_val_adj<0.05)%>%
    arrange(desc(abs(avg_log2FC)))%>%slice_head(n=30)
  res <- res %>%
    left_join(
      t %>% 
        dplyr::select(cluster, symbol) %>% 
        mutate(show = TRUE),
      by = c("cluster", "symbol")
    )
  res$show[is.na(res$show)] <- FALSE
  res$gene <- ifelse(res$show == TRUE & !grepl("^ENSG", res$symbol),
                     res$symbol,
                     "")
  
  plot<-ggplot(res, aes(x=avg_log2FC, y=-log10(p_val), label=gene, col=p_val_adj<0.05)) + 
    
    geom_point(size = 8, alpha = 0.5) +
    
    geom_text_repel(  size = 12,
                            force=10,
                            box.padding = 0.7,           # Increase padding around labels
                            point.padding = 0.7,         # Increase padding around points
                            min.segment.length = 0,      # Force connectors to display, even if short
                            segment.color = "gray50",    # Connector line color
                            segment.size = 0.6, max.overlaps = Inf)+
    scale_color_manual(values=c("#CCCCCC20","#FF000080")) +
    
    facet_wrap(~cluster, ncol=7, scales = "free_x", labeller = facet_labeller) + theme_bw() + NoLegend()+
    labs(x = expression(Log[2]~Fold~Change), y = expression(-Log[10](p~value)), color="q<0.05")+
    theme(
      text = element_text(size = 50),
      plot.title = element_text(hjust = 0.5),
      axis.title.y = element_text(margin = margin(r = -5))
    )
  
  output_dir<-"~/Immune_project/full_dataset_analysis/plots/DEG/"
  pdf(paste0(output_dir,Cell_Type,"_volcano_pseudo_bulk.pdf"), width = 45, height = 25)
  print(plot)
  dev.off()





# -------------------------------------------------------------------------





