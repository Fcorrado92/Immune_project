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
T_cells<-qread("/mnt/disks/disk/full_dataset_qs/T_cells_jan25_final.qs")
clusters_to_remove <- grep("^dbl", unique(Idents(T_cells)),value=T)
clusters_to_remove<-unique(c(clusters_to_remove, "Plt", "lowquality", "Lowq", "lowq","Plt"))
Idents(T_cells)<-T_cells$sub_cluster
keep<-setdiff(unique(Idents(T_cells)), clusters_to_remove)
T_cells<- subset(T_cells, idents = keep)
meta<-read_excel("~/Immune_project/full_dataset_analysis/Meta_Immune_proj.xlsx")
meta<-distinct(meta, Sample_ID, .keep_all = TRUE)

T_cells@meta.data$new_barcodes<-rownames(T_cells@meta.data)
T_cells@meta.data<-T_cells@meta.data[-c(44:51)]
T_cells@meta.data<-left_join(T_cells@meta.data, meta, by="Sample_ID")
rownames(T_cells@meta.data)<-T_cells@meta.data$new_barcodes
unique(T_cells$Therapy)
#REMOVE LEN AND SCREENING, KEEP CART TALQ TEC HEALTHY
T_cells_filt<-subset(T_cells, subset=Timepoint%in%c("Post","Pre","Day14","Day30","Day100","Day56","6Months",   
                                                    "Healthy")&Therapy%in%c("TEC","TALQ","CART","Healthy"))

T_cells_filt<-RenameIdents(T_cells_filt, c('CD8+_GZMK_TEM'="CD8+_GZMK+_TEM"))
T_cells_filt$sub_cluster<-Idents(T_cells_filt)



# Find IDs with paired samples pre and post in Bispecifics ----------------
df<-T_cells_filt@meta.data%>%filter(Therapy%in%c("TEC", "TALQ","CART")&Tissue=="BM")%>%group_by(ID, Therapy, Disease)%>%summarise(n=n_distinct(Sample_ID))
paired_ids<-df%>%filter(n>1)
paired_ids<-unique(paired_ids$ID)
paired_counts<-Tcells_counts%>%filter(Tissue=="PB"&ID%in%paired_ids)
df2<-paired_counts%>%filter(Therapy%in%c("TEC", "TALQ")&Tissue=="PB")%>%group_by(ID)%>%summarise(n=n_distinct(Sample_ID))
#remove the ID with ony one sample with n>50
paired_counts<-paired_counts%>%filter(!ID=="844295")

#############HRSMM####################
#do wilcox.test between Pre and Post for each celltype, store.pvalues, correct with BH
paired_counts_SMM<-paired_counts%>%filter(Disease=="HRSMM")
paired_counts_SMM%>%group_by(Disease, Therapy, Tissue, Timepoint)%>%summarise(n=n_distinct(Sample_ID))

celltypes<-unique(paired_counts_SMM$sub_cluster)
results <- list()  
for (i in seq_along(celltypes)) {  
  sub <- paired_counts_SMM %>% filter(sub_cluster == celltypes[i])  
  test <- wilcox.test(Freq ~ Timepoint, data = sub, paired=TRUE)  
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


#extract mean values for each subtype across disease groups
input_dir<-"~/Immune_project/full_dataset_analysis/plots/"
means<-paired_counts_SMM%>% group_by(Disease,Timepoint, sub_cluster)%>% summarise(mean_value=mean(Freq))
means_wide<-means%>% pivot_wider(names_from = Timepoint, values_from = mean_value)
means_wide<-means_wide%>%mutate(log2FC=log2(Post/Pre))
means_wide<-left_join(means_wide, adjusted_p_values, by="sub_cluster")
means_wide<-means_wide%>%mutate(enrichment= ifelse(log2FC>0,"Post","Pre"))
means_wide<-means_wide%>%mutate(significant= ifelse(p_adj<0.1,"q<0.1","q>0.1"))
#remove underscore for the graph but keep means_wide for further analysis
means_wide2<-means_wide
means_wide2$sub_cluster <- gsub("_", "", means_wide2$sub_cluster)
means_wide2$sub_cluster[means_wide2$sub_cluster=="cyclingTcells"]<-"Cycling T-cells"
means_wide2$sub_cluster[means_wide2$sub_cluster=="CD8+GZMBTEM"]<-"CD8+GZMB+TEM"


# -------------------------------------------------------------------------
#assign lineage
CD4<-grep("CD4+|Th|Treg", unique(means_wide2$sub_cluster), value=TRUE)
Mait<-"Mait"
Tgd<-"Tgd"
CD8<-setdiff(unique(means_wide2$sub_cluster), c(Tgd,Mait,CD4))
means_wide2 <- means_wide2 %>%
  mutate(Lineage = case_when(
    sub_cluster %in% CD4 ~ "CD4",
    sub_cluster == Mait ~ "Mait",
    sub_cluster == Tgd ~ "Tgd",
    sub_cluster %in% CD8 ~ "CD8",
    TRUE ~ NA_character_  # Default case if none of the conditions are met.
  ))

logFCmeans <- ggplot(means_wide2, aes(x = log2FC, y = -log10(pval), color = Lineage)) +
  geom_point(aes(shape = significant, 
                 alpha = ifelse(pval < 0.05, 1, 0.5)), size = 10) + 
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
  annotate("text", x = -0.2, y = 0.1, label = "Pre n=11", hjust = 1, size = 14) +  # Left text
  annotate("text", x = 1, y = 0.1, label = "Post n=11", hjust = 0, size = 14)
pdf(paste0(input_dir, "logFC_mean_diff_Tcells_Timepoint_PB_SMM.pdf"), width=18, height = 14)
logFCmeans
dev.off()


# -------------------------------------------------------------------------

cols<-c("gray","purple")
#do wilcox.test between RRMM and SMM for each celltype, store.pvalues, correct with BH
df<-means_wide2%>%filter(p_adj<0.1)
subsets<-unique(df$sub_cluster)
df$sub_cluster[df$sub_cluster=="CD8+ZEB2TEM"]<-"CD8+ZEB2+TEM"

Tcells_counts$sub_cluster<-gsub("_","",Tcells_counts$sub_cluster)
PBMC_sign<-Tcells_counts%>%filter(Disease%in%c("HRSMM") & Therapy=="TEC"& Tissue=="PB" & sub_cluster%in%subsets&ID%in%paired_ids)
PBMC_sign<-PBMC_sign%>%filter(!ID=="844295")

PBMC_sign$sub_cluster[PBMC_sign$sub_cluster=="CD8+ZEB2TEM"]<-"CD8+ZEB2+TEM"

PBMC_sign$sub_cluster<-gsub("_","",PBMC_sign$sub_cluster)
PBMC_sign$sub_cluster<-factor(PBMC_sign$sub_cluster, levels=c("CD8+Naive", 
                                                              "CD8+DR+TEM", "CD8+ZEB2+TEM","Mait", "Th17" ,
                                                              "Th2","Treg"             
))

PBMC_sign$Timepoint<-factor(PBMC_sign$Timepoint, levels=c("Pre","Post"))

my_pvals <- df %>%
  mutate(
    group1 = "Pre",
    group2 = "Post",
    # rename p_adj column to p.adj for ggpubr convenience
    p.adj = round(p_adj,3)
  ) %>%
  select(sub_cluster, group1, group2, p.adj)
max_vals <- PBMC_sign %>%
  group_by(sub_cluster) %>%
  summarise(y.position = max(Freq, na.rm = TRUE) * 1.1)

my_pvals <- left_join(my_pvals, max_vals, by = "sub_cluster")

p <- ggplot(PBMC_sign, aes(x=Timepoint, y=Freq, fill=Timepoint))+
  geom_violin(alpha=0.5) +
  geom_boxplot(size=3, alpha=0.6, outlier.shape = NA) +
  geom_line(aes(group=ID),size=0.8) +
  geom_jitter(size=10, width=0.1,aes(fill=Timepoint, stroke=2),color="black",shape=21, alpha=0.5)+
  labs(y="Percentage of Tcells")+
  theme_classic()+guides(color=FALSE)+
  scale_fill_manual(values = c("gray","purple"))+
  facet_wrap(~sub_cluster, scales="free", ncol=4)+
  theme(
    plot.margin = unit(c(0.5, 0, 0.2, 0.1), "cm"), 
    title=element_text(size=40),
    plot.title = element_text(hjust=0.5, face = "bold"), 
    text=element_text(size=50,color = "black"),
    axis.text.x = element_blank(),
    legend.position = "bottom",
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size=40,color = "black"), 
    axis.title.y = element_text(size=50,color = "black")
  )+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))

p<-p+stat_pvalue_manual(
  data = my_pvals,
  label = "p.adj",
  xmin = "group1",
  xmax = "group2",
  size = 15,
  y.position = "y.position",
  facet.var = "sub_cluster",
  inherit.aes = FALSE
)
ggsave("~/Immune_project/full_dataset_analysis/plots/boxplots/longitudinal_T_cells_SMM.pdf", 
       plot=p, width=30, height=20)

############RRMM###############
#do wilcox.test between RRMM and SMM for each celltype, store.pvalues, correct with BH
paired_counts_RRMM<-paired_counts%>%filter(Disease=="RRMM"&Therapy=="TEC")
paired_counts_RRMM%>%group_by(Disease, Therapy, Tissue, Timepoint)%>%summarise(n=n_distinct(Sample_ID))

celltypes<-unique(paired_counts_RRMM$sub_cluster)
results <- list()  
for (i in seq_along(celltypes)) {  
  sub <- paired_counts_RRMM %>% filter(sub_cluster == celltypes[i])  
  test <- wilcox.test(Freq ~ Timepoint, data = sub, paired=TRUE)  
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


#extract mean values for each subtype across disease groups
input_dir<-"~/Immune_project/full_dataset_analysis/plots/"
means<-paired_counts_RRMM%>% group_by(Disease,Timepoint, sub_cluster)%>% summarise(mean_value=mean(Freq))
means_wide<-means%>% pivot_wider(names_from = Timepoint, values_from = mean_value)
means_wide<-means_wide%>%mutate(log2FC=log2(Post/Pre))
means_wide<-left_join(means_wide, adjusted_p_values, by="sub_cluster")
means_wide<-means_wide%>%mutate(enrichment= ifelse(log2FC>0,"Post","Pre"))
means_wide<-means_wide%>%mutate(significant= ifelse(p_adj<0.1,"q<0.1","q>0.1"))
#remove underscore for the graph but keep means_wide for further analysis
means_wide2<-means_wide
means_wide2$sub_cluster <- gsub("_", "", means_wide2$sub_cluster)
means_wide2$sub_cluster[means_wide2$sub_cluster=="cyclingTcells"]<-"Cycling T-cells"
means_wide2$sub_cluster[means_wide2$sub_cluster=="CD8+GZMBTEM"]<-"CD8+GZMB+TEM"


# -------------------------------------------------------------------------
#assign lineage
CD4<-grep("CD4+|Th|Treg", unique(means_wide2$sub_cluster), value=TRUE)
Mait<-"Mait"
Tgd<-"Tgd"
CD8<-setdiff(unique(means_wide2$sub_cluster), c(Tgd,Mait,CD4))
means_wide2 <- means_wide2 %>%
  mutate(Lineage = case_when(
    sub_cluster %in% CD4 ~ "CD4",
    sub_cluster == Mait ~ "Mait",
    sub_cluster == Tgd ~ "Tgd",
    sub_cluster %in% CD8 ~ "CD8",
    TRUE ~ NA_character_  # Default case if none of the conditions are met.
  ))

logFCmeans <- ggplot(means_wide2, aes(x = log2FC, y = -log10(pval), color = Lineage)) +
  geom_point(aes(shape = significant, 
                 alpha = ifelse(pval < 0.05, 1, 0.5)), size = 10) + 
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
  annotate("text", x = -0.2, y = 0.1, label = "Pre n=8", hjust = 1, size = 14) +  # Left text
  annotate("text", x = 0.2, y = 0.1, label = "Post n=8", hjust = 0, size = 14)
pdf(paste0(input_dir, "logFC_mean_diff_Tcells_Timepoint_PB_RRMM.pdf"), width=18, height = 14)
logFCmeans
dev.off()


# -------------------------------------------------------------------------
cols<-c("gray","purple")
#do wilcox.test between Post and PRe in subset significant for SMM for each celltype, store.pvalues, correct with BH
subsets<-c(subsets,"CD8+GZMB+TEM")
df<-means_wide2%>%filter(sub_cluster%in%subsets)
df$sub_cluster[df$sub_cluster=="CD8+ZEB2TEM"]<-"CD8+ZEB2+TEM"
subsets<-unique(df$sub_cluster)
Tcells_counts$sub_cluster[Tcells_counts$sub_cluster=="CD8+GZMBTEM"]<-"CD8+GZMB+TEM"
Tcells_counts$sub_cluster[Tcells_counts$sub_cluster=="CD8+ZEB2TEM"]<-"CD8+ZEB2+TEM"

Tcells_counts$sub_cluster<-gsub("_","",Tcells_counts$sub_cluster)
PBMC_sign<-Tcells_counts%>%filter(Disease%in%c("RRMM") & Therapy=="TEC"& Tissue=="PB" & sub_cluster%in%subsets&ID%in%paired_ids)


PBMC_sign$sub_cluster<-gsub("_","",PBMC_sign$sub_cluster)
PBMC_sign$sub_cluster<-factor(PBMC_sign$sub_cluster, levels=c("CD8+Naive", 
                                                              "CD8+DR+TEM", "CD8+ZEB2+TEM","CD8+GZMB+TEM","Mait", "Th17" ,
                                                              "Th2","Treg"             
))

PBMC_sign$Timepoint<-factor(PBMC_sign$Timepoint, levels=c("Pre","Post"))

my_pvals <- df %>%
  mutate(
    group1 = "Pre",
    group2 = "Post",
    # rename p_adj column to p.adj for ggpubr convenience
    p.adj = round(p_adj,3)
  ) %>%
  select(sub_cluster, group1, group2, p.adj)
max_vals <- PBMC_sign %>%
  group_by(sub_cluster) %>%
  summarise(y.position = max(Freq, na.rm = TRUE) * 1.1)

my_pvals <- left_join(my_pvals, max_vals, by = "sub_cluster")

p <- ggplot(PBMC_sign, aes(x=Timepoint, y=Freq, fill=Timepoint))+
  geom_violin(alpha=0.5) +
  geom_boxplot(size=3, alpha=0.6, outlier.shape = NA) +
  geom_line(aes(group=ID),size=0.8) +
  geom_jitter(size=10, width=0.1,aes(fill=Timepoint, stroke=2),color="black",shape=21, alpha=0.5)+
  labs(y="Percentage of Tcells")+
  theme_classic()+guides(color=FALSE)+
  scale_fill_manual(values = c("gray","purple"))+
  facet_wrap(~sub_cluster, scales="free", ncol=4)+
  theme(
    plot.margin = unit(c(0.5, 0, 0.2, 0.1), "cm"), 
    title=element_text(size=40),
    plot.title = element_text(hjust=0.5, face = "bold"), 
    text=element_text(size=50,color = "black"),
    axis.text.x = element_blank(),
    legend.position = "bottom",
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size=40,color = "black"), 
    axis.title.y = element_text(size=50,color = "black")
  )+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))

p<-p+stat_pvalue_manual(
  data = my_pvals,
  label = "p.adj",
  xmin = "group1",
  xmax = "group2",
  size = 15,
  y.position = "y.position",
  facet.var = "sub_cluster",
  inherit.aes = FALSE
)
ggsave("~/Immune_project/full_dataset_analysis/plots/boxplots/longitudinal_T_cells_RRMM.pdf", 
       plot=p, width=30, height=20)


# PBMC pre and post in Bispecifics ----------------
counts<-read_csv("~/Immune_project/full_dataset_analysis/compositional_analysis/SC_PBMC_cell_counts_over500.csv")
T<-grep("TEM|Tem|Naive|Mait|Tex|Treg|TCM|Tcells|Tgd|Th2|Th17|Th1|TCM|CD4+", unique(counts$sub_cluster), value = TRUE)
Bcells<-grep("*BC|MZB|GCB", unique(counts$sub_cluster), value = TRUE)
Dcells<-grep("*DC", unique(counts$sub_cluster), value = TRUE)
NKcells<-grep("*NK", unique(counts$sub_cluster), value = TRUE)
Mono<-grep("Mono|Macrophage", unique(counts$sub_cluster), value = TRUE)
PC<-grep("Plasmacells", unique(counts$sub_cluster), value = TRUE)
Other<-setdiff(unique(counts$sub_cluster), c(T,Bcells, Dcells, NKcells, Mono, PC))
counts <- counts %>%
  mutate(
    # Assign major clusters based on sub-cluster groups
    major_cluster = case_when(
      sub_cluster %in% T ~ "T cells",
      sub_cluster %in% Mono ~ "Monocytes",
      sub_cluster %in% Dcells ~ "Dendritic Cells",
      sub_cluster %in% Bcells ~ "B cells",
      sub_cluster %in% NKcells ~ "NK cells",
      TRUE ~ "Other"  # Catch-all for clusters not matching any of the above
    )
  )

df<-counts%>%filter(Therapy%in%c("TEC", "TALQ")&Tissue=="PB")%>%group_by(ID)%>%summarise(n=n_distinct(Sample_ID))
paired_ids<-df%>%filter(n>1)
paired_ids<-unique(paired_ids$ID)
paired_counts<-counts%>%filter(Tissue=="PB"&ID%in%paired_ids)
df2<-paired_counts%>%filter(Therapy%in%c("TEC", "TALQ")&Tissue=="PB")%>%group_by(ID)%>%summarise(n=n_distinct(Sample_ID))
#remove the ID with ony one sample with n>50
unique(paired_counts$total)

#############HRSMM####################
#do wilcox.test between Pre and Post for each celltype, store.pvalues, correct with BH
paired_counts_SMM<-paired_counts%>%filter(Disease=="HRSMM")
paired_counts_SMM%>%group_by(Disease, Therapy, Tissue, Timepoint)%>%summarise(n=n_distinct(Sample_ID))

celltypes<-unique(paired_counts_SMM$sub_cluster)
results <- list()  
for (i in seq_along(celltypes)) {  
  sub <- paired_counts_SMM %>% filter(sub_cluster == celltypes[i])  
  test <- wilcox.test(Freq ~ Timepoint, data = sub, paired=TRUE)  
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


#extract mean values for each subtype across disease groups
input_dir<-"~/Immune_project/full_dataset_analysis/plots/"
means<-paired_counts_SMM%>% group_by(Disease,Timepoint,major_cluster, sub_cluster)%>% summarise(mean_value=mean(Freq))
means_wide<-means%>% pivot_wider(names_from = Timepoint, values_from = mean_value)
means_wide<-means_wide%>%mutate(log2FC=log2(Post/Pre))
means_wide<-left_join(means_wide, adjusted_p_values, by="sub_cluster")
means_wide<-means_wide%>%mutate(enrichment= ifelse(log2FC>0,"Post","Pre"))
means_wide<-means_wide%>%mutate(significant= ifelse(p_adj<0.1,"q<0.1","q>0.1"))
#remove underscore for the graph but keep means_wide for further analysis
means_wide2<-means_wide
means_wide2$sub_cluster <- gsub("_", "", means_wide2$sub_cluster)
means_wide2$sub_cluster[means_wide2$sub_cluster=="cyclingTcells"]<-"Cycling T-cells"
means_wide2$sub_cluster[means_wide2$sub_cluster=="CD8+GZMBTEM"]<-"CD8+GZMB+TEM"


# -------------------------------------------------------------------------
means_wide2$sub_cluster <- gsub("_", "", means_wide2$sub_cluster)
min_value<-min(means_wide2$log2FC[is.finite(means_wide2$log2FC)], na.rm = TRUE)
means_wide2$log2FC[means_wide2$log2FC=="-Inf"]<-min_value
logFCmeans <- ggplot(means_wide2, aes(x = log2FC, y = -log10(pval), color = major_cluster)) +
  geom_point(aes(shape = significant, 
                 alpha = ifelse(pval < 0.05, 1, 0.5)), size = 10) +  # Make points with pval<0.1 transparent
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "grey", linewidth=2) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey",linewidth=2) +
  scale_shape_manual(values = c("q<0.1" = 17, "q>0.1" = 19)) +  # Define shapes for significance
  scale_color_manual(values = c( "purple", "orange",  "darkblue", "darkgreen","red","gold"), 
                     name = "Cell Type") +  # Set custom colors for major_cluster and rename the legend title to "Cell Type"
  geom_text_repel(aes(label = ifelse(pval < 0.05, sub_cluster, "")), 
                  size = 10, max.overlaps = Inf,force = 20,              # Increase force to repel labels
                  box.padding = 1,      # Add padding around boxes
                  point.padding = 1) +  # Increase the text size
  labs(x = expression(Log[2]~Fold~Change), shape = "Significance",
       y = expression(-log[10](p~value)),
       fill = "Cell Type" ) +  # Set legend titles
  
  # theme_classic(base_size = 18) +
  
  theme(
    text=element_text(size=40),
    legend.position = c(0.54, 0.8),  # Position the legend at the top center of the plot
    legend.box = "horizontal",
    legend.spacing.x = unit(0.5, "cm"),  # Add space between legend items
    panel.background = element_blank(),  # Remove background from panel
    plot.background = element_blank(),  # Remove background from the plot
    legend.background = element_rect(fill = "white", color = "grey"),  # Add grey border around the legend box
    legend.key = element_blank()  # Remove the background and border from legend keys
  ) +
  guides(alpha = "none") +  # Remove the alpha legend
  annotate("text", x = -1, y = 0.1, label = "Pre n=11", hjust = 1, size = 14) +  # Left text
  annotate("text", x = 1, y = 0.1, label = "Post n=11", hjust = 0, size = 14)
pdf(paste0(input_dir, "logFC_mean_diff_PBMC_Timepoint_PB_SMM.pdf"), width=18, height = 14)
logFCmeans
dev.off()


# -------------------------------------------------------------------------

cols<-c("gray","purple")
#do wilcox.test between RRMM and SMM for each celltype, store.pvalues, correct with BH
df<-means_wide2%>%filter(p_adj<0.2)
subsets<-unique(df$sub_cluster)

counts$sub_cluster<-gsub("_","",counts$sub_cluster)
PBMC_sign<-counts%>%filter(Disease%in%c("HRSMM") & Therapy=="TEC"& Tissue=="PB" & sub_cluster%in%subsets&ID%in%paired_ids)

# PBMC_sign$sub_cluster<-factor(PBMC_sign$sub_cluster, levels=c("TBC","NBC","IFN+NBC", "MZB","TCF7+MBC","S100A10+MBC","GCB","ABC",
#                                                              "pDC","Mono-DC","Mait","CD8+ DR+ TEM"))

PBMC_sign$Timepoint<-factor(PBMC_sign$Timepoint, levels=c("Pre","Post"))

my_pvals <- df %>%
  mutate(
    group1 = "Pre",
    group2 = "Post",
    # rename p_adj column to p.adj for ggpubr convenience
    p.adj = round(p_adj,3)
  ) %>%
  select(sub_cluster, group1, group2, p.adj)%>%
  ungroup()
max_vals <- PBMC_sign %>%
  group_by(sub_cluster) %>%
  summarise(y.position = max(Freq, na.rm = TRUE) * 1.1)%>%
  ungroup()

my_pvals <- left_join(my_pvals, max_vals, by = "sub_cluster")
my_pvals<-my_pvals%>%select(c("Disease", "sub_cluster",  "group1", "group2" ,"p.adj" ,"y.position"))

p <- ggplot(PBMC_sign, aes(x=Timepoint, y=Freq, fill=Timepoint))+
  geom_violin(alpha=0.5) +
  geom_boxplot(size=3, alpha=0.6, outlier.shape = NA) +
  geom_line(aes(group=ID),size=0.8) +
  geom_jitter(size=10, width=0.1,aes(fill=Timepoint, stroke=2),color="black",shape=21, alpha=0.5)+
  labs(y="Percentage of PBMCs")+
  theme_classic()+guides(color=FALSE)+
  scale_fill_manual(values = c("gray","purple"))+
  facet_wrap(~sub_cluster, scales="free", ncol=6)+
  theme(
    plot.margin = unit(c(0.5, 0, 0.2, 0.1), "cm"), 
    title=element_text(size=40),
    plot.title = element_text(hjust=0.5, face = "bold"), 
    text=element_text(size=50,color = "black"),
    axis.text.x = element_blank(),
    legend.position = "bottom",
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size=40,color = "black"), 
    axis.title.y = element_text(size=50,color = "black")
  )+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))

p<-p+stat_pvalue_manual(
  data = my_pvals,
  label = "p.adj",
  xmin = "group1",
  xmax = "group2",
  size = 15,
  y.position = "y.position",
  facet.var = "sub_cluster",
  inherit.aes = FALSE
)
ggsave("~/Immune_project/full_dataset_analysis/plots/boxplots/longitudinal_PBMC_SMM.pdf", 
       plot=p, width=35, height=25)

############RRMM###############
#do wilcox.test between RRMM and SMM for each celltype, store.pvalues, correct with BH
paired_counts_RRMM<-paired_counts%>%filter(Disease=="RRMM"&Therapy=="TEC")
paired_counts_RRMM%>%group_by(Disease, Therapy, Tissue, Timepoint)%>%summarise(n=n_distinct(Sample_ID))

celltypes<-unique(paired_counts_RRMM$sub_cluster)
results <- list()  
for (i in seq_along(celltypes)) {  
  sub <- paired_counts_RRMM %>% filter(sub_cluster == celltypes[i])  
  test <- wilcox.test(Freq ~ Timepoint, data = sub, paired=TRUE)  
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


#extract mean values for each subtype across disease groups
input_dir<-"~/Immune_project/full_dataset_analysis/plots/"
means<-paired_counts_RRMM%>% group_by(Disease,Timepoint,major_cluster,sub_cluster)%>% summarise(mean_value=mean(Freq))
means_wide<-means%>% pivot_wider(names_from = Timepoint, values_from = mean_value)
means_wide<-means_wide%>%mutate(log2FC=log2(Post/Pre))
means_wide<-left_join(means_wide, adjusted_p_values, by="sub_cluster")
means_wide<-means_wide%>%mutate(enrichment= ifelse(log2FC>0,"Post","Pre"))
means_wide<-means_wide%>%mutate(significant= ifelse(p_adj<0.1,"q<0.1","q>0.1"))
#remove underscore for the graph but keep means_wide for further analysis
means_wide2<-means_wide
means_wide2$sub_cluster <- gsub("_", "", means_wide2$sub_cluster)
means_wide2$sub_cluster[means_wide2$sub_cluster=="cyclingTcells"]<-"Cycling T-cells"
means_wide2$sub_cluster[means_wide2$sub_cluster=="CD8+GZMBTEM"]<-"CD8+GZMB+TEM"


# -------------------------------------------------------------------------
logFCmeans <- ggplot(means_wide2, aes(x = log2FC, y = -log10(pval), color = major_cluster)) +
  geom_point(aes(shape = significant, 
                 alpha = ifelse(pval < 0.05, 1, 0.5)), size = 10) + 
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
    legend.position = c(0.54, 0.8),
    # Position the legend at the top center of the plot
    legend.box = "horizontal",  # Arrange legends horizontally
    legend.spacing.x = unit(0.5, "cm"),  # Add space between legend items
    panel.background = element_blank(),  # Remove background from panel
    plot.background = element_blank(),  # Remove background from the plot
    legend.background = element_rect(fill = "white", color = "grey"),  # Add grey border around the legend box
    legend.key = element_blank()  # Remove the background and border from legend keys
  ) +
  guides(alpha = "none") +  # Remove the alpha legend
  annotate("text", x = -0.2, y = 0.1, label = "Pre n=8", hjust = 1, size = 14) +  # Left text
  annotate("text", x = 0.2, y = 0.1, label = "Post n=8", hjust = 0, size = 14)
pdf(paste0(input_dir, "logFC_mean_diff_PBMC_Timepoint_PB_RRMM.pdf"), width=18, height = 14)
logFCmeans
dev.off()


# -------------------------------------------------------------------------
cols<-c("gray","purple")
#do wilcox.test between Post and PRe in subset significant for SMM for each celltype, store.pvalues, correct with BH
df<-means_wide2%>%filter(sub_cluster%in%subsets)

counts$sub_cluster<-gsub("_","",counts$sub_cluster)
PBMC_sign<-counts%>%filter(Disease%in%c("RRMM") & Therapy=="TEC"& Tissue=="PB" & sub_cluster%in%subsets&ID%in%paired_ids)
unique(PBMC_sign$sub_cluster)

PBMC_sign$sub_cluster<-gsub("_","",PBMC_sign$sub_cluster)
# PBMC_sign$sub_cluster<-factor(PBMC_sign$sub_cluster, levels=c("TBC","NBC","IFN+NBC", "MZB","TCF7+MBC","S100A10+MBC","GCB","ABC",
#                                                               "pDC","Mono-DC","Mait","CD8+ DR+ TEM"))
# 
PBMC_sign$Timepoint<-factor(PBMC_sign$Timepoint, levels=c("Pre","Post"))

my_pvals <- df %>%
  mutate(
    group1 = "Pre",
    group2 = "Post",
    # rename p_adj column to p.adj for ggpubr convenience
    p.adj = round(p_adj,3)
  ) %>%
  select(sub_cluster, group1, group2, p.adj)%>%
  ungroup()
max_vals <- PBMC_sign %>%
  group_by(sub_cluster) %>%
  summarise(y.position = max(Freq, na.rm = TRUE) * 1.1)%>%
  ungroup()

my_pvals <- left_join(my_pvals, max_vals, by = "sub_cluster")
my_pvals<-my_pvals%>%select(c("Disease", "sub_cluster",  "group1", "group2" ,"p.adj" ,"y.position"))

p <- ggplot(PBMC_sign, aes(x=Timepoint, y=Freq, fill=Timepoint))+
  geom_violin(alpha=0.5) +
  geom_boxplot(size=3, alpha=0.6, outlier.shape = NA) +
  geom_line(aes(group=ID),size=0.8) +
  geom_jitter(size=10, width=0.1,aes(fill=Timepoint, stroke=2),color="black",shape=21, alpha=0.5)+
  labs(y="Percentage of PBMCs")+
  theme_classic()+guides(color=FALSE)+
  scale_fill_manual(values = c("gray","purple"))+
  facet_wrap(~sub_cluster, scales="free", ncol=6)+
  theme(
    plot.margin = unit(c(0.5, 0, 0.2, 0.1), "cm"), 
    title=element_text(size=40),
    plot.title = element_text(hjust=0.5, face = "bold"), 
    text=element_text(size=50,color = "black"),
    axis.text.x = element_blank(),
    legend.position = "bottom",
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size=40,color = "black"), 
    axis.title.y = element_text(size=50,color = "black")
  )+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))

p<-p+stat_pvalue_manual(
  data = my_pvals,
  label = "p.adj",
  xmin = "group1",
  xmax = "group2",
  size = 15,
  y.position = "y.position",
  facet.var = "sub_cluster",
  inherit.aes = FALSE
)
ggsave("~/Immune_project/full_dataset_analysis/plots/boxplots/longitudinal_PBMC_RRMM.pdf", 
       plot=p, width=35, height=25)

