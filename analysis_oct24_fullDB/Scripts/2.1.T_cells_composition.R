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
T_cells<-qread("/mnt/disks/disk/full_dataset_qs/T_cells_jan25_final2.qs")
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
unique(T_cells$Therapy)
#REMOVE LEN AND SCREENING, KEEP CART TALQ TEC HEALTHY
T_cells_filt<-subset(T_cells, subset=Timepoint%in%c("Post","Pre","Day14","Day30","Day100","Day56","6Months",   
                                                    "Healthy")&Therapy%in%c("TEC","TALQ","CART","Healthy"))

T_cells_filt2<-RenameIdents(T_cells_filt, c('CD8+_DR+_TEM'="CD8+_GZMK+_TEM",'CD8+_Tex'="CD8+_GZMK+_TEM",'CD8+_KIR+_TEM'="CD8+_GZMB+_TEM" ))
T_cells_filt$sub_cluster<-Idents(T_cells_filt)

unique(Idents(T_cells_filt2))
# -------------------------------------------------------------------------


#save metadata for diversity
remove <- unique(rownames(T_cells_filt@meta.data[is.na(T_cells_filt@meta.data$clonotypeID), ]))
keep<-setdiff(unique(rownames(T_cells_filt@meta.data)), remove)
T_cells_filt_clono<-subset(T_cells_filt, cells = keep)
nrow(T_cells_filt_clono@meta.data)
#visualize percentage of T_cells with associated clonotype
metadata<-T_cells_filt_clono@meta.data
write_csv(metadata, "/mnt/disks/disk/full_dataset_qs/diversity/metadata_Tcells_filt.csv")

# -------------------------------------------------------------------------
ncells<-T_cells_filt@meta.data%>%group_by(Disease,Sample_ID,Tissue, Timepoint ,Cohort, Therapy)%>%
  summarise(n=n())

N_cells <- ggplot(ncells, aes(x = n)) +
  geom_histogram(bins = 70, aes(y = ..density..), fill = "lightblue", color = "black", alpha = 0.7) + # Histogram with density scaling
  geom_density(color = "red", size = 1) + # Overlay density line
  theme_classic() +
  scale_x_continuous(breaks = c(500, 1000, 2000, 5000, 10000)) +
  labs(y = "Density", x = "Number of Tcells", title = "Distribution of number of Tcells per sample") +
  geom_vline(xintercept = 100, linetype = "dotted", size = 1)+
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1))

ggsave("~/Immune_project/full_dataset_analysis/plots/N_cells/Ncells_baseline_tcells.pdf", plot=N_cells, width=7, height=4)


# -------------------------------------------------------------------------
counts<-as.data.frame(table(T_cells_filt$Sample_ID, T_cells_filt$sub_cluster))
colnames(counts)<-c("Sample_ID","sub_cluster","n")
total<-counts%>% group_by(Sample_ID)%>% summarise(total=sum(n))
counts<-left_join(counts, total, by="Sample_ID")
counts<-counts%>% mutate(Freq= n/total)
meta<-meta%>%distinct(Sample_ID,.keep_all = TRUE)
counts<-left_join(counts, meta[,c("Sample_ID","BOR","Timepoint", "Tissue","Therapy","Disease","ID", "Cohort")], by="Sample_ID")
counts<-counts%>%filter(total>50)
write_csv(counts, "~/Immune_project/full_dataset_analysis/compositional_analysis/SC_T_cell_counts_over50.csv")

#add levels
Tcells_counts<-read_csv("~/Immune_project/full_dataset_analysis/compositional_analysis/SC_T_cell_counts_over50.csv")
Tcells_counts$Disease[Tcells_counts$Disease=="SMM"]<-"HRSMM"
Tcells_counts$Disease<-factor(Tcells_counts$Disease, levels=c("Healthy","HRSMM","RRMM"))
Tcells_counts$BOR<-factor(Tcells_counts$BOR, levels=c("Healthy","NA","CR","VGPR","PR","MR","SD","PD"))
#assign response_disease
Tcells_counts <- Tcells_counts %>%
  mutate(response_disease = case_when(
    BOR %in% c("PR", "MR","SD","PD") & Disease == "RRMM" ~ "RRMM_NR",            # RRMM Non-Responder
    Disease == "RRMM" & !BOR %in% c("PR", "MR","SD","PD") ~ "RRMM_R",            # RRMM Responder
    Disease == "HRSMM" & !Therapy== "Len" ~ "HRSMM_R",                     # SMM Responder
    Disease == "Healthy" ~ "Healthy_NA",                              # Healthy: NA response
    TRUE ~ "SMM_NA"                                                    # Default: SMM Non-Applicable
  ))
Tcells_counts%>%group_by(Disease, Therapy, Timepoint, Tissue, response_disease, BOR)%>%summarise(n=n_distinct(Sample_ID))
# -------------------------------------------------------------------------
#by filtering "Pre" I am selecting BL MM samples, counts here is Tcells_counts >1000
MM_BL_Tcells_counts<-Tcells_counts%>%filter(Timepoint%in%c("Pre")&Tissue=="PB")
#check number of samples x disease x timepoint
MM_BL_Tcells_counts%>%group_by(Disease, Timepoint)%>%summarise(n=n_distinct(Sample_ID))
# Disease Timepoint     n
# 1 HRSMM   Pre          29
# 2 RRMM    Pre          39



#do wilcox.test between RRMM and SMM for each celltype, store.pvalues, correct with BH
celltypes<-unique(MM_BL_Tcells_counts$sub_cluster)
results <- list()  
for (i in seq_along(celltypes)) {  
  sub <- MM_BL_Tcells_counts %>% filter(sub_cluster == celltypes[i])  
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

MM_BL_Tcells_counts$Disease<-factor(MM_BL_Tcells_counts$Disease, levels=c("HRSMM","RRMM"))

MM_BL_Tcells_counts%>%group_by(Disease, Therapy, Tissue, Timepoint)%>%summarise(n=n_distinct(Sample_ID))
#extract mean values for each subtype across disease groups
input_dir<-"~/Immune_project/full_dataset_analysis/plots/"
means<-MM_BL_Tcells_counts%>% group_by(Disease, sub_cluster)%>% summarise(mean_value=mean(Freq))
means_wide<-means%>% pivot_wider(names_from = Disease, values_from = mean_value)
means_wide<-means_wide%>%mutate(log2FC=log2(RRMM/HRSMM))
means_wide<-left_join(means_wide, adjusted_p_values, by="sub_cluster")
means_wide<-means_wide%>%mutate(enrichment= ifelse(log2FC>0,"RRMM","HRSMM"))
means_wide<-means_wide%>%mutate(significant= ifelse(p_adj<0.1,"q<0.1","q>0.1"))
#remove underscore for the graph but keep means_wide for further analysis
means_wide2<-means_wide
means_wide2$sub_cluster <- gsub("_", "", means_wide2$sub_cluster)
means_wide2$sub_cluster[means_wide2$sub_cluster=="cyclingTcells"]<-"Cycling T-cells"


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
  annotate("text", x = -0.2, y = 0.1, label = "HRSMM n=29", hjust = 1, size = 14) +  # Left text
  annotate("text", x = 1, y = 0.1, label = "RRMM n=39", hjust = 0, size = 14)
pdf(paste0(input_dir, "logFC_mean_diff_Tcells.pdf"), width=18, height = 14)
logFCmeans
dev.off()


# -------------------------------------------------------------------------

cols<-c("darkgreen","darkblue","darkred")
#do wilcox.test between RRMM and SMM for each celltype, store.pvalues, correct with BH
df<-means_wide%>%filter(p_adj<0.15)
subsets<-unique(df$sub_cluster)
PBMC_sign<-Tcells_counts%>%filter(Timepoint%in%c("Pre","Healthy") & Tissue=="PB" & sub_cluster%in%subsets)


PBMC_sign$sub_cluster<-gsub("_","",PBMC_sign$sub_cluster)
PBMC_sign$sub_cluster<-factor(PBMC_sign$sub_cluster, levels=c("CD8+Naive", 
                                                              "CD8+GZMK+TEM", "CD8+DR+TEM", "CD8+Tex" ,
                                                              "CD4+Naive","CD4+TCM" ,"Th1" ,"cyclingTcells"                 
                                                               ))

PBMC_sign$Disease<-factor(PBMC_sign$Disease, levels=c("Healthy","HRSMM","RRMM"))
p<-ggplot(PBMC_sign, aes(x=Disease, y=Freq, fill=Disease))+
  geom_violin(alpha=0.5)+geom_boxplot(size=3,alpha=0.6, outlier.shape = NA)+
  geom_jitter(size=10,width=0.1,aes(fill=Disease, stroke=2),color="black",shape=21, alpha=0.5)+
  labs(y="Percentage of Tcells")+
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

ggsave("~/Immune_project/full_dataset_analysis/plots/boxplots/boxplots_sign_Tcells_cluster_disease_baseline.pdf", 
       plot=p, width=30, height=20)


# BM---------------------------------------------------------
#by filtering "Pre" I am selecting BL MM samples, counts here is Tcells_counts >500
MM_BL_Tcells_counts_BM<-Tcells_counts%>%filter(Timepoint%in%c("Pre")&Tissue=="BM")
#check number of samples x disease x timepoint
MM_BL_Tcells_counts_BM%>%group_by(Disease, Timepoint)%>%summarise(n=n_distinct(ID))
# 1 HRSMM   Pre          14
# 2 RRMM    Pre           9



#do wilcox.test between RRMM and SMM for each celltype, store.pvalues, correct with BH
celltypes<-unique(MM_BL_Tcells_counts_BM$sub_cluster)
results <- list()  
for (i in seq_along(celltypes)) {  
  sub <- MM_BL_Tcells_counts_BM %>% filter(sub_cluster == celltypes[i])  
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

MM_BL_Tcells_counts_BM$Disease<-factor(MM_BL_Tcells_counts_BM$Disease, levels=c("HRSMM","RRMM"))

MM_BL_Tcells_counts_BM%>%group_by(Disease, Therapy, Tissue, Timepoint)%>%summarise(n=n_distinct(Sample_ID))
#extract mean values for each subtype across disease groups
input_dir<-"~/Immune_project/full_dataset_analysis/plots/"
means<-MM_BL_Tcells_counts_BM%>% group_by(Disease, sub_cluster)%>% summarise(mean_value=mean(Freq))
means_wide<-means%>% pivot_wider(names_from = Disease, values_from = mean_value)
means_wide<-means_wide%>%mutate(log2FC=log2(RRMM/HRSMM))
means_wide<-left_join(means_wide, adjusted_p_values, by="sub_cluster")
means_wide<-means_wide%>%mutate(enrichment= ifelse(log2FC>0,"RRMM","HRSMM"))
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
  annotate("text", x = -0.2, y = 0.1, label = "HRSMM n=14", hjust = 1, size = 14) +  # Left text
  annotate("text", x = 1, y = 0.1, label = "RRMM n=9", hjust = 0, size = 14)
pdf(paste0(input_dir, "logFC_mean_diff_Tcells_BM.pdf"), width=18, height = 14)
logFCmeans
dev.off()


# -------------------------------------------------------------------------

cols<-c("darkgreen","darkblue","darkred")
#do wilcox.test between RRMM and SMM for each celltype, store.pvalues, correct with BH
df<-means_wide2%>%filter(p_adj<0.1)
subsets<-unique(df$sub_cluster)
Tcells_counts$sub_cluster<-gsub("_","",Tcells_counts$sub_cluster)
Tcells_counts$sub_cluster[Tcells_counts$sub_cluster=="CD8+GZMBTEM"]<-"CD8+GZMB+TEM"
PBMC_sign<-Tcells_counts%>%filter(Timepoint%in%c("Pre","Healthy") & Tissue=="BM" & sub_cluster%in%subsets)
PBMC_sign$sub_cluster<-gsub("_","",PBMC_sign$sub_cluster)
PBMC_sign$sub_cluster<-factor(PBMC_sign$sub_cluster, levels=c("CD8+Naive", 
                                                              "CD8+GZMK+TEM", "CD8+DR+TEM","CD8+GZMB+TEM", "CD8+Tex" ,
                                                              "CD4+Naive","Th1"             
))

PBMC_sign$Disease<-factor(PBMC_sign$Disease, levels=c("Healthy","HRSMM","RRMM"))
p<-ggplot(PBMC_sign, aes(x=Disease, y=Freq, fill=Disease))+
  geom_violin(alpha=0.5)+geom_boxplot(size=3,alpha=0.6, outlier.shape = NA)+
  geom_jitter(size=10,width=0.1,aes(fill=Disease, stroke=2),color="black",shape=21, alpha=0.5)+
  labs(y="Percentage of Tcells")+
  theme_classic()+guides(color=FALSE)+
  scale_fill_manual(values = cols)+
  facet_wrap(~sub_cluster, scales="free", ncol=3)+
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

ggsave("~/Immune_project/full_dataset_analysis/plots/boxplots/BM_boxplots_sign_Tcells_cluster_disease_baseline.pdf", 
       plot=p, width=30, height=20)


# -------------------------------------------------------------------------
#Look at the proportions in expanded space
# -------------------------------------------------------------------------
T_cells_clono<-qread("/mnt/disks/disk/full_dataset_qs/diversity/T_cells_clono.qs")
T_cells_expanded<-subset(T_cells_clono, clonotype_size%in%c("Small", "Medium","Large"))
#Lokk at the number of cells x sample
ncells<-T_cells_expanded@meta.data%>%group_by(Disease,Sample_ID,Tissue, Timepoint ,Cohort, Therapy)%>%
  summarise(n=n())

N_cells <- ggplot(ncells, aes(x = n)) +
  geom_histogram(bins = 70, aes(y = ..density..), fill = "lightblue", color = "black", alpha = 0.7) + # Histogram with density scaling
  geom_density(color = "red", size = 1) + # Overlay density line
  theme_classic() +
  scale_x_continuous(breaks = c(500, 1000, 2000, 5000, 10000)) +
  labs(y = "Density", x = "Number of expanded Tcells", title = "Distribution of number of Tcells per sample") +
  geom_vline(xintercept = 100, linetype = "dotted", size = 1)+
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1))

ggsave("~/Immune_project/full_dataset_analysis/plots/N_cells/Ncells_baseline_tcells_expanded.pdf", plot=N_cells, width=7, height=4)

counts<-as.data.frame(table(T_cells_expanded$Sample_ID, T_cells_expanded$sub_cluster))
colnames(counts)<-c("Sample_ID","sub_cluster","n")
total<-counts%>% group_by(Sample_ID)%>% summarise(total=sum(n))
counts<-left_join(counts, total, by="Sample_ID")
counts<-counts%>% mutate(Freq= n/total)
meta<-meta%>%distinct(Sample_ID,.keep_all = TRUE)
counts<-left_join(counts, meta[,c("Sample_ID","BOR","Timepoint", "Tissue","Therapy","Disease","ID", "Cohort")], by="Sample_ID")
write_csv(counts, "~/Immune_project/full_dataset_analysis/compositional_analysis/SC_T_cell_expanded_counts.csv")

#add levels
Tcells_counts<-read_csv("~/Immune_project/full_dataset_analysis/compositional_analysis/SC_T_cell_expanded_counts.csv")
Tcells_counts$Disease[Tcells_counts$Disease=="SMM"]<-"HRSMM"
Tcells_counts$Disease<-factor(Tcells_counts$Disease, levels=c("Healthy","HRSMM","RRMM"))
Tcells_counts$BOR<-factor(Tcells_counts$BOR, levels=c("Healthy","NA","CR","VGPR","PR","MR","SD","PD"))
#assign response_disease
Tcells_counts <- Tcells_counts %>%
  mutate(response_disease = case_when(
    BOR %in% c("PR", "MR","SD","PD") & Disease == "RRMM" ~ "RRMM_NR",            # RRMM Non-Responder
    Disease == "RRMM" & !BOR %in% c("PR", "MR","SD","PD") ~ "RRMM_R",            # RRMM Responder
    Disease == "HRSMM" & !Therapy== "Len" ~ "HRSMM_R",                     # SMM Responder
    Disease == "Healthy" ~ "Healthy_NA",                              # Healthy: NA response
    TRUE ~ "SMM_NA"                                                    # Default: SMM Non-Applicable
  ))
Tcells_counts%>%group_by(Disease, Therapy, Timepoint, Tissue, response_disease, BOR)%>%summarise(n=n_distinct(Sample_ID))
# -------------------------------------------------------------------------
#by filtering "Pre" I am selecting BL MM samples, counts here is Tcells_counts >1000
MM_BL_Tcells_counts<-Tcells_counts%>%filter(Timepoint%in%c("Pre")&Tissue=="PB"&Therapy%in%c("TEC","TALQ","CART"))
#check number of samples x disease x timepoint
MM_BL_Tcells_counts<-MM_BL_Tcells_counts%>%filter(total>10)
MM_BL_Tcells_counts%>%group_by(Disease, Timepoint)%>%summarise(n=n_distinct(Sample_ID))
# Disease Timepoint     n
# 1 HRSMM   Pre          27
# 2 RRMM    Pre          35



#do wilcox.test between RRMM and SMM for each celltype, store.pvalues, correct with BH
celltypes<-unique(MM_BL_Tcells_counts$sub_cluster)
results <- list()  
for (i in seq_along(celltypes)) {  
  sub <- MM_BL_Tcells_counts %>% filter(sub_cluster == celltypes[i])  
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

MM_BL_Tcells_counts$Disease<-factor(MM_BL_Tcells_counts$Disease, levels=c("HRSMM","RRMM"))
MM_BL_Tcells_counts%>%group_by(Disease, Therapy, Tissue, Timepoint)%>%summarise(n=n_distinct(Sample_ID))
#extract mean values for each subtype across disease groups
input_dir<-"~/Immune_project/full_dataset_analysis/plots/"
means<-MM_BL_Tcells_counts%>% group_by(Disease, sub_cluster)%>% summarise(mean_value=mean(Freq))
means_wide<-means%>% pivot_wider(names_from = Disease, values_from = mean_value)
means_wide<-means_wide%>%mutate(log2FC=log2(RRMM/HRSMM))
means_wide<-left_join(means_wide, adjusted_p_values, by="sub_cluster")
means_wide<-means_wide%>%mutate(enrichment= ifelse(log2FC>0,"RRMM","HRSMM"))
means_wide<-means_wide%>%mutate(significant= ifelse(p_adj<0.1,"q<0.1","q>0.1"))
#remove underscore for the graph but keep means_wide for further analysis
means_wide2<-means_wide
means_wide2$sub_cluster <- gsub("_", "", means_wide2$sub_cluster)
means_wide2$sub_cluster[means_wide2$sub_cluster=="cyclingTcells"]<-"Cycling T-cells"
means_wide2$sub_cluster[means_wide2$sub_cluster=="CD8+GZMKTEM"]<-"CD8+GZMK+TEM"


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
                 alpha = ifelse(pval < 0.1, 1, 0.5)), size = 10) +  # Make points with pval<0.1 transparent
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "grey", linewidth=2) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey",linewidth=2) +
  scale_shape_manual(values = c("q<0.1" = 17, "q>0.1" = 19)) +  # Define shapes for significance
  scale_color_manual(values = c( "purple", "orange",  "darkblue", "darkgreen","red","gold"), 
                     name = "Lineage") +  # Set custom colors for major_cluster and rename the legend title to "Cell Type"
  geom_text_repel(aes(label = ifelse(pval < 0.1, sub_cluster, "")), 
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
  annotate("text", x = -0.2, y = 0.1, label = "HRSMM n=21", hjust = 1, size = 14) +  # Left text
  annotate("text", x = 1, y = 0.1, label = "RRMM n=34", hjust = 0, size = 14)+
  xlim(c(-2, 5))

pdf(paste0(input_dir, "logFC_mean_diff_Tcells_expanded.pdf"), width=18, height = 14)
logFCmeans
dev.off()


# -------------------------------------------------------------------------
#visualize in each sample contribution of various phenotype to expanded clones faceted for Disease
# -------------------------------------------------------------------------

cols<-c("darkgreen","darkblue","darkred")
#do wilcox.test between RRMM and SMM for each celltype, store.pvalues, correct with BH
df<-means_wide%>%filter(p_adj<0.1)
subsets<-unique(df$sub_cluster)
PBMC_sign<-Tcells_counts%>%filter(Timepoint%in%c("Pre","Healthy") & Tissue=="PB" & Therapy%in%c("Healthy","TEC","TALQ","CART")& sub_cluster%in%subsets)


PBMC_sign$sub_cluster<-gsub("_","",PBMC_sign$sub_cluster)
PBMC_sign$sub_cluster<-factor(PBMC_sign$sub_cluster, levels=c("CD8+Naive", 
                                                              "CD8+GZMK+TEM", "CD8+DR+TEM", "CD8+Tex" ,
                                                              "CD4+Naive","CD4+TCM" ,"Th1" ,"cyclingTcells"                 
))

PBMC_sign$Disease<-factor(PBMC_sign$Disease, levels=c("Healthy","HRSMM","RRMM"))
p<-ggplot(PBMC_sign, aes(x=Disease, y=Freq, fill=Disease))+
  geom_violin(alpha=0.5)+geom_boxplot(size=3,alpha=0.6, outlier.shape = NA)+
  geom_jitter(size=10,width=0.1,aes(fill=Disease, stroke=2),color="black",shape=21, alpha=0.5)+
  labs(y="Percentage of expanded Tcells(TCR>1%)")+
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

ggsave("~/Immune_project/full_dataset_analysis/plots/boxplots/boxplots_sign_Tcells_cluster_disease_baseline_expanded.pdf", 
       plot=p, width=30, height=20)




# -------------------------------------------------------------------------
#by filtering "Pre" I am selecting BL MM samples, counts here is Tcells_counts >1000
MM_BL_Tcells_counts<-Tcells_counts%>%filter(Timepoint%in%c("Pre")&Tissue=="BM"&Therapy%in%c("TEC","TALQ","CART"))
#check number of samples x disease x timepoint
MM_BL_Tcells_counts<-MM_BL_Tcells_counts%>%filter(total>10)
MM_BL_Tcells_counts%>%group_by(Disease, Timepoint)%>%summarise(n=n_distinct(Sample_ID))
# Disease Timepoint     n
# 1 HRSMM   Pre          7
# 2 RRMM    Pre          8



#do wilcox.test between RRMM and SMM for each celltype, store.pvalues, correct with BH
celltypes<-unique(MM_BL_Tcells_counts$sub_cluster)
results <- list()  
for (i in seq_along(celltypes)) {  
  sub <- MM_BL_Tcells_counts %>% filter(sub_cluster == celltypes[i])  
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

MM_BL_Tcells_counts$Disease<-factor(MM_BL_Tcells_counts$Disease, levels=c("HRSMM","RRMM"))
MM_BL_Tcells_counts%>%group_by(Disease, Therapy, Tissue, Timepoint)%>%summarise(n=n_distinct(Sample_ID))
#extract mean values for each subtype across disease groups
input_dir<-"~/Immune_project/full_dataset_analysis/plots/"
means<-MM_BL_Tcells_counts%>% group_by(Disease, sub_cluster)%>% summarise(mean_value=mean(Freq))
means_wide<-means%>% pivot_wider(names_from = Disease, values_from = mean_value)
means_wide<-means_wide%>%mutate(log2FC=log2(RRMM/HRSMM))
means_wide<-left_join(means_wide, adjusted_p_values, by="sub_cluster")
means_wide<-means_wide%>%mutate(enrichment= ifelse(log2FC>0,"RRMM","HRSMM"))
means_wide<-means_wide%>%mutate(significant= ifelse(p_adj<0.1,"q<0.1","q>0.1"))
#remove underscore for the graph but keep means_wide for further analysis
means_wide2<-means_wide
means_wide2$sub_cluster <- gsub("_", "", means_wide2$sub_cluster)
means_wide2$sub_cluster[means_wide2$sub_cluster=="cyclingTcells"]<-"Cycling T-cells"
means_wide2$sub_cluster[means_wide2$sub_cluster=="CD8+GZMKTEM"]<-"CD8+GZMK+TEM"


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
                 alpha = ifelse(pval < 0.1, 1, 0.5)), size = 10) +  # Make points with pval<0.1 transparent
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "grey", linewidth=2) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey",linewidth=2) +
  scale_shape_manual(values = c("q<0.1" = 17, "q>0.1" = 19)) +  # Define shapes for significance
  scale_color_manual(values = c( "purple", "orange",  "darkblue", "darkgreen","red","gold"), 
                     name = "Lineage") +  # Set custom colors for major_cluster and rename the legend title to "Cell Type"
  geom_text_repel(aes(label = ifelse(pval < 0.1, sub_cluster, "")), 
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
  annotate("text", x = -0.2, y = 0.1, label = "HRSMM n=7", hjust = 1, size = 14) +  # Left text
  annotate("text", x = 1, y = 0.1, label = "RRMM n=8", hjust = 0, size = 14)+
  xlim(c(-2, 5))

pdf(paste0(input_dir, "logFC_mean_diff_Tcells_expanded_BM.pdf"), width=18, height = 14)
logFCmeans
dev.off()


# -------------------------------------------------------------------------
#visualize in each sample contribution of various phenotype to expanded clones faceted for Disease
# -------------------------------------------------------------------------

cols<-c("darkgreen","darkblue","darkred")
#do wilcox.test between RRMM and SMM for each celltype, store.pvalues, correct with BH
df<-means_wide%>%filter(p_adj<0.1)
subsets<-unique(df$sub_cluster)
PBMC_sign<-Tcells_counts%>%filter(Timepoint%in%c("Pre","Healthy") & Tissue=="BM" & Therapy%in%c("Healthy","TEC","TALQ","CART")& sub_cluster%in%subsets)


PBMC_sign$sub_cluster<-gsub("_","",PBMC_sign$sub_cluster)
PBMC_sign$sub_cluster<-factor(PBMC_sign$sub_cluster, levels=c("CD8+Naive", 
                                                              "CD8+GZMK+TEM", "CD8+DR+TEM", "CD8+Tex" ,
                                                              "CD4+Naive","CD4+TCM" ,"Th1" ,"cyclingTcells"                 
))

PBMC_sign$Disease<-factor(PBMC_sign$Disease, levels=c("Healthy","HRSMM","RRMM"))
p<-ggplot(PBMC_sign, aes(x=Disease, y=Freq, fill=Disease))+
  geom_violin(alpha=0.5)+geom_boxplot(size=3,alpha=0.6, outlier.shape = NA)+
  geom_jitter(size=10,width=0.1,aes(fill=Disease, stroke=2),color="black",shape=21, alpha=0.5)+
  labs(y="Percentage of expanded Tcells(TCR>1%)")+
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

ggsave("~/Immune_project/full_dataset_analysis/plots/boxplots/boxplots_sign_Tcells_cluster_disease_baseline_expanded_BM.pdf", 
       plot=p, width=30, height=20)


# -------------------------------------------------------------------------
#Look at the correlation of PB and BM composition
# -------------------------------------------------------------------------
df<-Tcells_counts%>%mutate(ID_TP=paste0(ID,"_",Timepoint))%>%group_by(ID_TP, Disease, Therapy)%>%summarise(n=n_distinct(Tissue))%>%filter(n>1)
IDS<-unique(df$ID_TP)
Tcells_counts<-Tcells_counts%>%mutate(ID_TP=paste0(ID,"_",Timepoint))
Tcells_counts<-left_join(Tcells_counts, meta, by="Sample_ID")
paired_df<-Tcells_counts%>%filter(ID_TP%in%IDS & sub_cluster=="CD8+_Tex")%>%select("ID_TP","Freq","Tissue",)
paired_df<-paired_df%>%  pivot_wider(names_from = "Tissue", values_from = "Freq")

Tcells_counts2<-Tcells_counts%>%filter(Therapy=="TEC"&Timepoint=="Pre"&Tissue=="PM")
Tcells_counts2 <- Tcells_counts2 %>%
  mutate(response_disease = case_when(
    BOR %in% c("SD","MR", "PR","PD") & Disease == "RRMM" ~ "RRMM_NR",            # RRMM Non-Responder
    Disease == "RRMM" & !BOR %in% c( "SD","MR", "PR","PD") ~ "RRMM_R",            # RRMM Responder
    Disease == "HRSMM" & !Therapy== "Len" ~ "HRSMM_R",                     # SMM Responder
    Disease == "Healthy" ~ "Healthy_NA",                              # Healthy: NA response
    TRUE ~ "SMM_NA"                                                    # Default: SMM Non-Applicable
  ))
Tcells_counts2$response_disease<-factor(Tcells_counts2$response_disease, levels=c("HRSMM_R","RRMM_R","RRMM_NR"))
p<-ggplot(Tcells_counts2, aes(x=response_disease, y=Freq, fill=response_disease))+
  geom_violin(alpha=0.5)+geom_boxplot(alpha=0.5)+
  geom_jitter(size=1,width=0.1,aes(fill=response_disease, stroke=2),color="black",shape=21, alpha=0.5)+
  labs(y="Percentage of Tcells")+
  theme_classic()+guides(color=FALSE)+
  scale_fill_manual(values = cols)+
  facet_wrap(~sub_cluster, scales="free", ncol=4)
p<-p+ geom_pwc(label.size = 2, method = "wilcox.test",step.increase = 0.08, tip.length = 0.01, size=1.2)
p<-ggadjust_pvalue(
  p, p.adjust.method = "BH",
  label = paste("q=","{p.adj.format}"),
  hide.ns = "p"
)



# -------------------------------------------------------------------------
#Tregs
# -------------------------------------------------------------------------




CD4<-grep("CD4+|Th|Treg", unique(T_cells_filt$sub_cluster), value=TRUE)
Mait<-"Mait"
Tgd<-"Tgd"
CD8<-setdiff(unique(T_cells_filt$sub_cluster), c(Tgd,Mait,CD4))
CD4_Tcells<-subset(T_cells_filt,subset=sub_cluster%in%c("CD4+_TCM" ,  "Th2"      ,  "CD4+_Naive", "Th17"     ,  "Treg"    ,   "Th1"    ,    "CD4+_aTEM" , "CD4+TfH"))
CD4_baseline<-subset(CD4_Tcells, Timepoint=="Pre"& Disease%in%c("RRMM","HRSMM")&Therapy%in%c("TEC","TALQ"))
l<-CD4_baseline@meta.data%>%group_by(ID,Timepoint)%>%summarise(n=n_distinct(Tissue))
l<-l%>%filter(n>1)
l<-l%>%mutate(ID_Tissue=paste0(ID,"_","BM"))
ids<-unique(l$ID_Tissue)
keep<-setdiff(unique(CD4_baseline$ID_Tissue),ids)
CD4_baseline@meta.data<-CD4_baseline@meta.data%>%mutate(ID_Tissue=paste0(ID,"_",Tissue))
CD4_baseline<-subset(CD4_baseline, ID_Tissue%in%keep)
CD4_baseline$sub_cluster<-droplevels(CD4_baseline$sub_cluster)
counts<-as.data.frame(table(CD4_baseline$Sample_ID, CD4_baseline$sub_cluster))
colnames(counts)<-c("Sample_ID","sub_cluster","n")
total<-counts%>% group_by(Sample_ID)%>% summarise(total=sum(n))
counts<-left_join(counts, total, by="Sample_ID")
counts<-counts%>% mutate(Freq= n/total)
meta<-meta%>%distinct(Sample_ID,.keep_all = TRUE)
counts<-left_join(counts, meta[,c("Sample_ID","BOR","Timepoint", "Tissue","Therapy","Disease","ID", "Cohort")], by="Sample_ID")
Tcells_counts<-counts
Tcells_counts$Disease[Tcells_counts$Disease=="SMM"]<-"HRSMM"
Tcells_counts$Disease<-factor(Tcells_counts$Disease, levels=c("Healthy","HRSMM","RRMM"))
Tcells_counts$BOR<-factor(Tcells_counts$BOR, levels=c("Healthy","NA","CR","VGPR","PR","MR","SD","PD"))
#assign response_disease
Tcells_counts <- Tcells_counts %>%
  mutate(response_disease = case_when(
    BOR %in% c("PR", "MR","SD","PD") & Disease == "RRMM" ~ "RRMM_NR",            # RRMM Non-Responder
    Disease == "RRMM" & !BOR %in% c("PR", "MR","SD","PD") ~ "RRMM_R",            # RRMM Responder
    Disease == "HRSMM" & !Therapy== "Len" ~ "HRSMM_R",                     # SMM Responder
    Disease == "Healthy" ~ "Healthy_NA",                              # Healthy: NA response
    TRUE ~ "SMM_NA"                                                    # Default: SMM Non-Applicable
  ))
Tcells_counts%>%group_by(Disease, Therapy, Timepoint, Tissue, response_disease, BOR)%>%summarise(n=n_distinct(Sample_ID))
Tcells_counts$response_disease<-factor(Tcells_counts$response_disease, levels=c("HRSMM_R","RRMM_R","RRMM_NR"))
p<-ggplot(Tcells_counts, aes(x=response_disease, y=Freq, fill=response_disease))+
  geom_violin(alpha=0.5)+geom_boxplot(size=3,alpha=0.6, outlier.shape = NA)+
  geom_jitter(size=10,width=0.1,aes(fill=response_disease, stroke=2),color="black",shape=21, alpha=0.5)+
  labs(y="CD4+ T cells(%)")+
  theme_classic()+guides(color=FALSE)+
  scale_fill_manual(values = c("steelblue", "gold","orange"))+
  facet_wrap(~sub_cluster, scales="free", ncol=4)+
  theme(plot.margin = unit(c(0.5, 0, 0.2, 0.1), "cm"),title=element_text(size=40),
        plot.title = element_text(hjust=0.5, face = "bold"), text=element_text(size=50,color = "black"),
        axis.text.x =element_blank(),legend.position = "bottom",
        axis.ticks.x = element_blank(),
        axis.text.y =element_text(size=40,color = "black"), axis.title.y =element_text(size=50,color = "black"))+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))
p<-p+ geom_pwc(label.size = 12, method = "wilcox.test",step.increase = 0.08, tip.length = 0.01, size=1.2)


# -------------------------------------------------------------------------
CD8_Tcells<-subset(T_cells_filt,subset=sub_cluster%in%c("CD8+_TCM", "CD8+_GZMK+_TEM", "CD8+_Naive", "CD8+_KIR+_TEM",
                                                        "CD8+_Tex", "CD8+_GZMB_TEM", "CD8+_DR+_TEM",
                                                        "cycling_Tcells"))
CD8_baseline<-subset(CD8_Tcells, Timepoint=="Pre"& Disease%in%c("RRMM","HRSMM")&Therapy%in%c("TEC","TALQ"))
l<-CD8_baseline@meta.data%>%group_by(ID,Timepoint)%>%summarise(n=n_distinct(Tissue))
l<-l%>%filter(n>1)
l<-l%>%mutate(ID_Tissue=paste0(ID,"_","BM"))
ids<-unique(l$ID_Tissue)
CD8_baseline@meta.data<-CD8_baseline@meta.data%>%mutate(ID_Tissue=paste0(ID,"_",Tissue))
keep<-setdiff(unique(CD8_baseline$ID_Tissue),ids)
CD8_baseline<-subset(CD8_baseline, ID_Tissue%in%keep)
CD8_baseline$sub_cluster<-droplevels(CD8_baseline$sub_cluster)
counts<-as.data.frame(table(CD8_baseline$Sample_ID, CD8_baseline$sub_cluster))
colnames(counts)<-c("Sample_ID","sub_cluster","n")
total<-counts%>% group_by(Sample_ID)%>% summarise(total=sum(n))
counts<-left_join(counts, total, by="Sample_ID")
counts<-counts%>% mutate(Freq= n/total)
meta<-meta%>%distinct(Sample_ID,.keep_all = TRUE)
counts<-left_join(counts, meta[,c("Sample_ID","BOR","Timepoint", "Tissue","Therapy","Disease","ID", "Cohort")], by="Sample_ID")
counts<-counts%>%filter(n>10)
Tcells_counts<-counts
Tcells_counts$Disease[Tcells_counts$Disease=="SMM"]<-"HRSMM"
Tcells_counts$Disease<-factor(Tcells_counts$Disease, levels=c("Healthy","HRSMM","RRMM"))
Tcells_counts$BOR<-factor(Tcells_counts$BOR, levels=c("Healthy","NA","CR","VGPR","PR","MR","SD","PD"))
#assign response_disease
Tcells_counts <- Tcells_counts %>%
  mutate(response_disease = case_when(
    BOR %in% c("PR", "MR","SD","PD") & Disease == "RRMM" ~ "RRMM_NR",            # RRMM Non-Responder
    Disease == "RRMM" & !BOR %in% c("PR", "MR","SD","PD") ~ "RRMM_R",            # RRMM Responder
    Disease == "HRSMM" & !Therapy== "Len" ~ "HRSMM_R",                     # SMM Responder
    Disease == "Healthy" ~ "Healthy_NA",                              # Healthy: NA response
    TRUE ~ "SMM_NA"                                                    # Default: SMM Non-Applicable
  ))
Tcells_counts%>%group_by(Disease, Therapy, Timepoint, Tissue, response_disease, BOR)%>%summarise(n=n_distinct(Sample_ID))
Tcells_counts$response_disease<-factor(Tcells_counts$response_disease, levels=c("HRSMM_R","RRMM_R","RRMM_NR"))
p<-ggplot(Tcells_counts, aes(x=response_disease, y=Freq, fill=response_disease))+
  geom_violin(alpha=0.5)+geom_boxplot(size=3,alpha=0.6, outlier.shape = NA)+
  geom_jitter(size=10,width=0.1,aes(fill=response_disease, stroke=2),color="black",shape=21, alpha=0.5)+
  labs(y="CD8+ T cells(%)")+
  theme_classic()+guides(color=FALSE)+
  scale_fill_manual(values = c("steelblue", "gold","orange"))+
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
ggsave("~/Immune_project/full_dataset_analysis/plots/boxplots/CD8_baseline.pdf",
       plot=p, width=30, height=20)



# -------------------------------------------------------------------------
#rename gzmk and gzmb
# -------------------------------------------------------------------------
counts<-as.data.frame(table(T_cells_filt2$Sample_ID, T_cells_filt2$sub_cluster))
colnames(counts)<-c("Sample_ID","sub_cluster","n")
total<-counts%>% group_by(Sample_ID)%>% summarise(total=sum(n))
counts<-left_join(counts, total, by="Sample_ID")
counts<-counts%>% mutate(Freq= n/total)
meta<-meta%>%distinct(Sample_ID,.keep_all = TRUE)
counts<-left_join(counts, meta[,c("Sample_ID","BOR","Timepoint", "Tissue","Therapy","Disease","ID", "Cohort")], by="Sample_ID")
counts<-counts%>%filter(total>50)
Tcells_counts<-counts
#add levels
Tcells_counts$Disease[Tcells_counts$Disease=="SMM"]<-"HRSMM"
Tcells_counts$Disease<-factor(Tcells_counts$Disease, levels=c("Healthy","HRSMM","RRMM"))
Tcells_counts$BOR<-factor(Tcells_counts$BOR, levels=c("Healthy","NA","CR","VGPR","PR","MR","SD","PD"))
#assign response_disease
Tcells_counts <- Tcells_counts %>%
  mutate(response_disease = case_when(
    BOR %in% c("PR", "MR","SD","PD") & Disease == "RRMM" ~ "RRMM_NR",            # RRMM Non-Responder
    Disease == "RRMM" & !BOR %in% c("PR", "MR","SD","PD") ~ "RRMM_R",            # RRMM Responder
    Disease == "HRSMM" & !Therapy== "Len" ~ "HRSMM_R",                     # SMM Responder
    Disease == "Healthy" ~ "Healthy_NA",                              # Healthy: NA response
    TRUE ~ "SMM_NA"                                                    # Default: SMM Non-Applicable
  ))
Tcells_counts%>%group_by(Disease, Therapy, Timepoint, Tissue, response_disease, BOR)%>%summarise(n=n_distinct(Sample_ID))

Tcells_counts<-Tcells_counts%>%filter(Timepoint%in%c("Pre") & Therapy%in%c("TEC","TALQ"))
l<-Tcells_counts%>%group_by(ID,Timepoint)%>%summarise(n=n_distinct(Tissue))
l<-l%>%filter(n>1)
l<-l%>%mutate(ID_Tissue=paste0(ID,"_","BM"))
ids<-unique(l$ID_Tissue)

Tcells_counts<-Tcells_counts%>%mutate(ID_Tissue=paste0(ID,"_",Tissue))
Tcells_counts<-Tcells_counts%>%filter(!ID_Tissue%in%ids)


Tcells_counts$Disease<-factor(Tcells_counts$Disease, levels=c("Healthy","HRSMM","RRMM"))
p<-ggplot(Tcells_counts, aes(x=response_disease, y=Freq, fill=response_disease))+
  geom_violin(alpha=0.5)+geom_boxplot(size=3,alpha=0.6, outlier.shape = NA)+
  geom_jitter(size=10,width=0.1,aes(fill=response_disease, stroke=2),color="black",shape=21, alpha=0.5)+
  labs(y="Percentage of expanded Tcells(TCR>1%)")+
  theme_classic()+guides(color=FALSE)+
  scale_fill_manual(values = c("steelblue","gold","orange"))+
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

ggsave("~/Immune_project/full_dataset_analysis/plots/boxplots/response_v2.pdf", 
       plot=p, width=30, height=20)


