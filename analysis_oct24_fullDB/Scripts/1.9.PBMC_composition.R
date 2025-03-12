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
paired_pts<-c("383408",
              "514771",
              "565525",
              "617885",
              "634084",
              "674945",
              "763623",
              "771538",
              "778550",
              "844295",
              "929373",
              "935266",
              "977028",
              "1000524",
              "1032305",
              "1066795",
              "1073409",
              "1097876",
              "1100299")

#load all the objs
integrated<-qread("~/Immune_project/full_dataset_analysis/integrated_filtered.qs")
T_cells<-qread("/mnt/disks/disk/full_dataset_qs/T_cells_jan25_final2.qs")
meta_Tcells<-T_cells@meta.data
write_csv(meta_Tcells,"/mnt/disks/disk/full_dataset_qs/meta_T_cells_jan25_final2.csv")
Idents(T_cells)<-T_cells$sub_cluster
Mono<-qread("/mnt/disks/cellranger/full_dataset_qs/Mono_annotated.qs")
meta_Mono<-Mono@meta.data
write_csv(meta_Mono,"/mnt/disks/disk/full_dataset_qs/meta_Mono.csv")

NK_cells<-qread("/mnt/disks/disk/full_dataset_qs/NK_cells_annotated.qs")
Idents(NK_cells)<-NK_cells$ident
meta_NK<-NK_cells@meta.data
write_csv(meta_NK,"/mnt/disks/disk/full_dataset_qs/meta_NK.csv")

B_cells<-qread("/mnt/disks/cellranger/full_dataset_qs/Bcells_annotated.qs")
meta_B<-B_cells@meta.data
write_csv(meta_B,"/mnt/disks/disk/full_dataset_qs/meta_B.csv")

Idents(B_cells)<-B_cells$sub_cluster
DC<-qread("/mnt/disks/cellranger/full_dataset_qs/DC_annotated.qs")
meta_D<-DC@meta.data
write_csv(meta_D,"/mnt/disks/disk/full_dataset_qs/meta_DC.csv")

#assign clusters to integrated obj
Idents(integrated)<-integrated$sub_cluster
integrated$sub_cluster <- as.character(Idents(integrated))
integrated$sub_cluster[Cells(T_cells)] <- paste(Idents(T_cells))
integrated$sub_cluster[Cells(Mono)] <- paste(Idents(Mono))
integrated$sub_cluster[Cells(DC)] <- paste(Idents(DC))
integrated$sub_cluster[Cells(NK_cells)] <- paste(Idents(NK_cells))
integrated$sub_cluster[Cells(B_cells)] <- paste(Idents(B_cells))
integrated[['ident']]<-integrated$sub_cluster
Idents(integrated)<-integrated$sub_cluster
DimPlot(integrated, label=F, repel = T)
rm(T_cells)
rm(B_cells)
rm(NK_cells)
rm(Mono)
rm(DC)
gc()
integrated <- RenameIdents(integrated, c( '21'="EP", '18'="EP",'20'="HSC",
                           '19'="GMP",'30'="Pro-B",'32'="lowq",'17' = "Plasmacells", '38' = "Plasmacells", 
                           '36' = "Plasmacells",'26'="GMP",'27'="dbl.Mono.Precursors",'16'="Plt"))
integrated$sub_cluster<-Idents(integrated)
integrated[['ident']]<-integrated$sub_cluster
DimPlot(integrated, label=T, repel = T)+NoLegend()

Idents(integrated)<-integrated$sub_cluster


# -------------------------------------------------------------------------


#remove low quality and doublets cells
clusters_to_remove <- grep("^dbl", unique(Idents(integrated)),value=T)
clusters_to_remove<-unique(c(clusters_to_remove, "Plt", "lowquality", "Lowq", "lowq","Plt"))
keep<-setdiff(unique(Idents(integrated)), clusters_to_remove)
integrated<- subset(integrated, idents = keep)
#add metadata
meta<-read_excel("~/Immune_project/full_dataset_analysis/Meta_Immune_proj.xlsx")
meta<-distinct(meta, Sample_ID, .keep_all = TRUE)

View(integrated@meta.data)
integrated@meta.data$new_barcodes<-rownames(integrated@meta.data)
integrated@meta.data<-integrated@meta.data[-c(44:51)]
integrated@meta.data<-left_join(integrated@meta.data, meta, by="Sample_ID")
rownames(integrated@meta.data)<-integrated@meta.data$new_barcodes

# remove CAR_screening ----------------------------------------------------
table(integrated$Therapy)
keep<-setdiff(unique(integrated$Therapy), c("Len","TEC_SCREEN","CART_Screening"))
integrated<-subset(integrated, subset=Therapy%in%keep)
qsave(integrated, "/mnt/disks/disk/full_dataset_qs/integrated_only_Tx_filt_annotated_jan25.qs")
# -------------------------------------------------------------------------
Idents(integrated)<-integrated$sub_cluster
clusters_to_remove <- grep("^dbl", unique(Idents(integrated)),value=T)
clusters_to_remove<-unique(c(clusters_to_remove, "Plt", "lowquality", "Lowq", "lowq","Plt"))
keep<-setdiff(unique(Idents(integrated)), clusters_to_remove)
integrated<- subset(integrated, idents = keep)

# -------------------------------------------------------------------------
#metadata summary
metadata<-integrated@meta.data
t<-metadata%>%group_by(Disease, Therapy, Tissue, Timepoint, Cohort)%>%summarise(n_sample=n_distinct(Sample_ID))
t<-as.data.frame(t)

#UMAP
colors <- c(
  "#E41A1C",  # Red
  "#377EB8",  # Blue
  "#4DAF4A",  # Green
  "#984EA3",  # Purple
  "#FF7F00",  # Orange
  "#A65628",  # Brown
  "#F781BF",  # Light Pink
  "#A6CEE3",  # Light Blue
  "#1F78B4",  # Dark Blue
  "#B2DF8A",  # Light Green
  "#33A02C",  # Dark Green
  "#FB9A99",  # Pink
  "#FF69B4",  # Fuchsia
  "#FDBF6F",  # Yellow Orange
  "#8B4513",  # Saddle Brown
  "#FFD700",  # Gold
  "#7B68EE",  # Medium Slate Blue
  "#6A5ACD",  # Slate Blue
  "#8B0000",  # Dark Red
  "#00CED1",  # Dark Turquoise
  "#FF4500",  # Orange Red
  "#DA70D6",  # Orchid
  "#87CEEB",  # Sky Blue
  "#5F9EA0",  # Cadet Blue
  "#98FB98",  # Pale Green
  "#CD5C5C",  # Indian Red
  "#20B2AA",  # Light Sea Green
  "#4682B4",  # Steel Blue
  "#DB7093",  # Pale Violet Red
  "#66CDAA",  # Medium Aquamarine
  "#FF6347",  # Tomato
  "#40E0D0",  # Turquoise
  "#778899",  # Light Slate Gray
  "#32CD32",  # Lime Green
  "#FF1493",  # Deep Pink
  "#BDB76B",  # Dark Khaki
  "#6B8E23",  # Olive Drab
  "#9932CC",  # Dark Orchid
  "#DC143C",  # Crimson
  "#7FFF00",  # Chartreuse
  "#ADFF2F",  # Green Yellow
  "#8A2BE2",  # Blue Violet
  "#FF8C00",  # Dark Orange
  "#8FBC8F",  # Dark Sea Green
  "#FFE4C4",  # Bisque
  "#6495ED",  # Cornflower Blue
  "#BC8F8F",  # Rosy Brown
  "#FFDEAD",
  "gold",
  "darkorange",
  "blue"# Navajo White
)

integrated<-RenameIdents(integrated, 'CD8+_GZMK_TEM'="CD8+_GZMK+_TEM", 'CD8+_GZMB_TEM' ="CD8+_GZMB+_TEM",
                         'CD8+_ZEB2_TEM'="CD8+_ZEB2+_TEM"
)
integrated[['ident']]<-Idents(integrated)
#UMAP sub_cluster Fig1+ number of patients + number of samples+ number of cells
integrated$sub_cluster<-Idents(integrated)
integrated$sub_cluster<-gsub("_", " ", integrated$sub_cluster)
integrated[['ident']]<-integrated$sub_cluster

UMAP <- DimPlot(integrated) + 
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    text = element_text(size = 14),
    axis.ticks = element_blank(),
    axis.text = element_blank()
  ) + 
  labs(x = "UMAP 1", y = "UMAP 2") +  
  scale_color_manual(values = colors) + 
  NoLegend()

# Add labels with reduced overlap
UMAP <- LabelClusters(plot = UMAP, id = "ident", repel = TRUE, size = 4, nudge_y = -0.2, nudge_x=+0.2,
                      max.overlaps=getOption("ggrepel.max.overlaps", default = 15), force=10)

ggsave("~/Immune_project/full_dataset_analysis/plots/UMAP/UMAP_integrated.pdf",
       plot=UMAP, width=17, height=14, units="cm")

int_meta<-integrated_filt@meta.data
write_csv(int_meta,"/mnt/disks/disk/full_dataset_qs/integrated_meta.csv")

# look at number of cells x major clusters --------------------------------
T<-grep("TEM|Tem|Naive|Mait|Tex|Treg|TCM|Tcells|Tgd|Th2|Th17|Th1|TCM|CD4+", unique(integrated$sub_cluster), value = TRUE)
Bcells<-grep("*BC|MZB|GCB", unique(integrated$sub_cluster), value = TRUE)
Dcells<-grep("*DC", unique(integrated$sub_cluster), value = TRUE)
NKcells<-grep("*NK", unique(integrated$sub_cluster), value = TRUE)
Mono<-grep("Mono|Macrophage", unique(integrated$sub_cluster), value = TRUE)
PC<-grep("Plasmacells", unique(integrated$sub_cluster), value = TRUE)
Other<-setdiff(unique(integrated$sub_cluster), c(T,Bcells, Dcells, NKcells, Mono, PC))

integrated$major_cluster <- as.character(Idents(integrated))
integrated$major_cluster[c(integrated$sub_cluster%in%T)]<-"T cells"
integrated$major_cluster[c(integrated$sub_cluster%in%Bcells)]<-"B cells"
integrated$major_cluster[c(integrated$sub_cluster%in%Dcells)]<-"Dendritic cells"
integrated$major_cluster[c(integrated$sub_cluster%in%NKcells)]<-"NK cells"
integrated$major_cluster[c(integrated$sub_cluster%in%Mono)]<-"Monocytes"
integrated$major_cluster[c(integrated$sub_cluster%in%PC)]<-"Plasmacells"
integrated$major_cluster[c(integrated$sub_cluster%in%Other)]<-"Hematopoietic precursors"

integrated[['ident']]<-integrated$major_cluster
DimPlot(integrated)
table(integrated$major_cluster)

#look at Tissue differences
UMAP <- DimPlot(integrated, label = FALSE, repel = TRUE) + 
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    text = element_text(size = 14),
    axis.ticks = element_blank(),
    axis.text = element_blank()
  ) + 
  labs(x = "UMAP 1", y = "UMAP 2") +  
  scale_color_manual(values = colors) + 
  NoLegend()

# Add labels with reduced overlap
UMAP <- LabelClusters(plot = UMAP, id = "ident", repel = TRUE, size = 4, nudge_y = -0.2, nudge_x=+0.2,
                      max.overlaps=getOption("ggrepel.max.overlaps", default = 15), force=10)

ggsave("~/Immune_project/full_dataset_analysis/plots/UMAP/UMAP_MC.pdf",
       plot=UMAP, width=17, height=14, units="cm")



# Look at the distribution of PBMC number ---------------------------------
ncells<-int_meta%>%group_by(Library,Sample_ID)%>%
  summarise(n=n())

N_cells <- ggplot(ncells, aes(x = n)) +
  geom_histogram(bins = 70, aes(y = ..density..), fill = "lightblue", color = "black", alpha = 0.7) + # Histogram with density scaling
  geom_density(color = "red", size = 1) + # Overlay density line
  theme_classic() +
  scale_x_continuous(breaks = c(500, 1000, 2000, 5000, 10000)) +
  labs(y = "Density", x = "Number of PBMCs", title = "Distribution of number of PBMCs per sample") +
  geom_vline(xintercept = 500, linetype = "dotted", size = 1)+
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1))

ggsave("~/Immune_project/full_dataset_analysis/plots/N_cells/Ncells_baseline_pbmcs.pdf", plot=N_cells, width=7, height=4)


# calculate n x sub_cluster and total at sample level-------------------------------------------------------------------------
integrated<-subset(integrated, major_cluster%in%c("T cells","B cells", "Dendritic cells", "NK cells","Monocytes"))
integrated$major_cluster <- factor(integrated$major_cluster)
integrated$major_cluster <- droplevels(integrated$major_cluster)
counts<-as.data.frame(table(integrated$Sample_ID, integrated$major_cluster))
colnames(counts)<-c("Sample_ID","major_cluster","n")
total<-counts%>% group_by(Sample_ID)%>% summarise(total=sum(n))
counts<-left_join(counts, total, by="Sample_ID")
counts<-counts%>% mutate(Freq= n/total)
meta<-meta%>%distinct(Sample_ID,.keep_all = TRUE)
counts<-left_join(counts, meta[,c("Sample_ID","BOR","Timepoint", "Therapy","Tissue", "Disease","ID", "Cohort")], by="Sample_ID")
#remove less than 500cells
counts<-counts%>%filter(total>=500)
write_csv(counts, "~/Immune_project/full_dataset_analysis/compositional_analysis/MC_PBMC_cell_counts_over500.csv")

# -------------------------------------------------------------------------


cols<-c("darkgreen","darkblue","darkred")
#do wilcox.test between RRMM and SMM for each celltype, store.pvalues, correct with BH
counts$Disease<-factor(counts$Disease, levels=c("Healthy","HRSMM","RRMM"))
counts<-counts%>%filter(Tissue=="PB"&Timepoint%in%c("Healthy","Pre")&Therapy%in%c("TEC","TALQ","CART","Healthy"))
p<-ggplot(counts, aes(x=Disease, y=Freq, fill=Disease))+
  geom_violin(alpha=0.5)+geom_boxplot(size=3,alpha=0.6, outlier.shape = NA)+
  geom_jitter(size=10,width=0.1,aes(fill=Disease, stroke=2),color="black",shape=21, alpha=0.5)+
  labs(y="Percentage of PBMCs")+
  theme_classic()+guides(color=FALSE)+
  scale_fill_manual(values = cols)+
  facet_wrap(~major_cluster, scales="free", ncol=5)+
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

ggsave("~/Immune_project/full_dataset_analysis/plots/boxplots/boxplots_major_cluster_baseline_SMM_RRMM.pdf", 
       plot=p, width=25, height=15)

# Volcano Plot Cell types -------------------------------------------------
# unique(counts$Timepoint)
# [1] "Pre"        "Healthy"    "TEC_SCREEN"
# > unique(counts$Therapy)
# [1] "TEC"        "CART"       "Healthy"    "Len"        "TEC_SCREEN" "TALQ"      
# > unique(counts$Disease)
# [1] HRSMM   Healthy RRMM   
# unique(counts$BOR)
# [1] "CR"      "Healthy" "NA"      "PR"      "PD"      "VGPR" 
integrated<-read_csv("/mnt/disks/disk/full_dataset_qs/integrated_meta.csv")
integrated<-integrated%>%filter(!major_cluster%in%c("Hematopoietic precursors","Plasmacells"))
integrated$sub_cluster <- factor(integrated$sub_cluster)
integrated$sub_cluster <- droplevels(integrated$sub_cluster)
counts<-as.data.frame(table(integrated$Sample_ID, integrated$sub_cluster))
colnames(counts)<-c("Sample_ID","sub_cluster","n")
total<-counts%>% group_by(Sample_ID)%>% summarise(total=sum(n))
counts<-left_join(counts, total, by="Sample_ID")
counts<-counts%>% mutate(Freq= n/total)
meta<-meta%>%distinct(Sample_ID,.keep_all = TRUE)
counts<-left_join(counts, meta[,c("Sample_ID","BOR","Timepoint","Tissue", "Therapy","Disease","ID", "Cohort")], by="Sample_ID")
#remove less than 1000cells
counts<-counts%>%filter(total>=500)
write_csv(counts, "~/Immune_project/full_dataset_analysis/compositional_analysis/SC_PBMC_cell_counts_over500.csv")
counts<-read_csv("~/Immune_project/full_dataset_analysis/compositional_analysis/SC_PBMC_cell_counts_over500.csv")
#n of pts in the two groups
PBMC<-counts
#add levels
PBMC$Disease<-factor(PBMC$Disease, levels=c("Healthy","HRSMM","RRMM"))
PBMC$BOR<-factor(PBMC$BOR, levels=c("Healthy","NA","CR","VGPR","PR","PD"))
#assign response_disease
PBMC <- PBMC %>%
  mutate(response_disease = case_when(
    BOR %in% c("PR", "PD") & Disease == "RRMM" ~ "RRMM_NR",            # RRMM Non-Responder
    Disease == "RRMM" & !BOR %in% c("PR", "PD") ~ "RRMM_R",            # RRMM Responder
    Disease == "HRSMM" & !Therapy== "Len" ~ "SMM_R",                     # SMM Responder
    Disease == "Healthy" ~ "Healthy_NA",                              # Healthy: NA response
    TRUE ~ "SMM_NA"                                                    # Default: SMM Non-Applicable
  ))


table(PBMC$Therapy)
table(PBMC$Tissue)
table(PBMC$Timepoint)

#by filtering "Pre" I am selecting BL MM samples, counts here is PBMC >1000
MM_BL_PBMC<-PBMC%>%filter(Timepoint%in%c("Pre"), Therapy%in%c("CART","TEC","TALQ") & Tissue=="PB")
MM_BL_PBMC$sub_cluster<-factor(MM_BL_PBMC$sub_cluster)
MM_BL_PBMC$sub_cluster<-droplevels(MM_BL_PBMC$sub_cluster)
#check number of samples x disease x timepoint
MM_BL_PBMC%>%group_by(Disease, Timepoint)%>%summarise(n=n_distinct(Sample_ID))
#    Disease Timepoint     n
#   1 HRSMM   Pre          29
#   2 RRMM    Pre          38
#assign major clusters
# Create the major cluster column based on the sub-cluster groups
MM_BL_PBMC <- MM_BL_PBMC %>%
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


#do wilcox.test between RRMM and SMM for each celltype, store.pvalues, correct with BH
celltypes<-unique(MM_BL_PBMC$sub_cluster)
results <- list()  
for (i in seq_along(celltypes)) {  
  sub <- MM_BL_PBMC %>% filter(sub_cluster == celltypes[i])  
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
MM_BL_PBMC$Disease<-factor(MM_BL_PBMC$Disease, levels=c("HRSMM","RRMM"))
MM_BL_PBMC$sub_cluster<-factor(MM_BL_PBMC$sub_cluster)
MM_BL_PBMC$sub_cluster<-droplevels(MM_BL_PBMC$sub_cluster)


#extract mean values for each subtype across disease groups
input_dir<-"~/Immune_project/full_dataset_analysis/plots/"
means<-MM_BL_PBMC%>% group_by(Disease, sub_cluster, major_cluster)%>% summarise(mean_value=mean(Freq))
means_wide<-means%>% pivot_wider(names_from = Disease, values_from = mean_value)
means_wide<-means_wide%>%mutate(log2FC=log2(RRMM/HRSMM))
means_wide<-left_join(means_wide, adjusted_p_values, by="sub_cluster")
means_wide<-means_wide%>%mutate(enrichment= ifelse(log2FC>0,"RRMM","HRSMM"))
means_wide<-means_wide%>%mutate(significant= ifelse(p_adj<0.1,"q<0.1","q>0.1"))
#remove underscore for the graph but keep means_wide for further analysis
means_wide2<-means_wide
means_wide2$sub_cluster <- gsub("_", "", means_wide2$sub_cluster)
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
    legend.box = "horizontal",  # Arrange legends horizontally
    legend.spacing.x = unit(0.5, "cm"),  # Add space between legend items
    panel.background = element_blank(),  # Remove background from panel
    plot.background = element_blank(),  # Remove background from the plot
    legend.background = element_rect(fill = "white", color = "grey"),  # Add grey border around the legend box
    legend.key = element_blank()  # Remove the background and border from legend keys
  ) +
  guides(alpha = "none") +  # Remove the alpha legend
  annotate("text", x = -1, y = 0.1, label = "HRSMM n=29", hjust = 1, size = 14) +  # Left text
  annotate("text", x = 1, y = 0.1, label = "RRMM n=38", hjust = 0, size = 14)
# Display the plot
print(logFCmeans)

pdf(paste0(input_dir, "logFC_mean_diff_integrated.pdf"), width=20, height = 16)
logFCmeans
dev.off()


# -------------------------------------------------------------------------
#correlation with response
# -------------------------------------------------------------------------
counts<-read_csv("~/Immune_project/full_dataset_analysis/compositional_analysis/SC_T_cell_counts_over50.csv")
PBMC<-counts
#add levels
PBMC$Disease<-factor(PBMC$Disease, levels=c("Healthy","HRSMM","RRMM"))
PBMC$BOR<-factor(PBMC$BOR, levels=c("Healthy","NA","CR","VGPR","PR","PD"))
#assign response_disease
PBMC <- PBMC %>%
  mutate(response_disease = case_when(
    BOR %in% c("MR","SD","PR", "PD") & Disease == "RRMM" ~ "RRMM_NR",            # RRMM Non-Responder
    Disease == "RRMM" & !BOR %in% c("MR","SD","PR", "PD") ~ "RRMM_R",            # RRMM Responder
    Disease == "HRSMM" & !Therapy== "Len" ~ "HRSMM_R",                     # SMM Responder
    Disease == "Healthy" ~ "Healthy_NA",                              # Healthy: NA response
    TRUE ~ "SMM_NA"                                                    # Default: SMM Non-Applicable
  ))

table(PBMC$Therapy)
table(PBMC$Tissue)
table(PBMC$Timepoint)

#by filtering "Pre" I am selecting BL MM samples, counts here is PBMC >1000
MM_BL_PBMC<-PBMC%>%filter(Timepoint%in%c("Pre"), Therapy%in%c("TEC") & Tissue=="PB")
MM_BL_PBMC$sub_cluster<-factor(MM_BL_PBMC$sub_cluster)
MM_BL_PBMC$sub_cluster<-droplevels(MM_BL_PBMC$sub_cluster)
#check number of samples x disease x timepoint
MM_BL_PBMC%>%group_by(Disease, Timepoint)%>%summarise(n=n_distinct(Sample_ID))
MM_BL_PBMC$response_disease<-factor(MM_BL_PBMC$response_disease, c("HRSMM_R","RRMM_R", "RRMM_NR"))
p<-ggplot(MM_BL_PBMC, aes(x=response_disease, y=Freq, fill=Disease))+
  geom_boxplot()+facet_wrap(~sub_cluster,scales="free")+geom_point()
p<-p+ geom_pwc(method = "wilcox.test",step.increase = 0.08, tip.length = 0.01)
# p<-ggadjust_pvalue(
#   p, p.adjust.method = "BH",
#   label = paste("q=","{p.adj.format}"))

ggsave(filename = "~/Immune_project/full_dataset_analysis/plots/PBMC_response.pdf", plot=p, width=30, heigh=30)