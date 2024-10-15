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
integrated<-qread("~/Immune_project/new_analysis_sep24/integrated_annotated.qs")
T_cells<-qread("~/Immune_project/new_analysis_sep24/qs_files/T_cells_annotated.qs")
Mono<-qread("~/Immune_project/new_analysis_sep24/qs_files/Mono.qs")
NK_cells<-qread("~/Immune_project/new_analysis_sep24/qs_files/NKcells.qs")
B_cells<-qread("~/Immune_project/new_analysis_sep24/qs_files/Bcell.qs")
DC<-qread("~/Immune_project/new_analysis_sep24/qs_files/DCells.qs")

#assign clusters to integrated obj
Idents(integrated)<-integrated$seurat_clusters
integrated$sub_cluster <- as.character(Idents(integrated))
integrated$sub_cluster[Cells(T_cells)] <- paste(Idents(T_cells))
integrated$sub_cluster[Cells(Mono)] <- paste(Idents(Mono))
integrated$sub_cluster[Cells(DC)] <- paste(Idents(DC))
integrated$sub_cluster[Cells(NK_cells)] <- paste(Idents(NK_cells))
integrated$sub_cluster[Cells(B_cells)] <- paste(Idents(B_cells))
integrated[['ident']]<-integrated$sub_cluster
Idents(integrated)<-integrated$sub_cluster
DimPlot(integrated, label=T)


#rename idents
integrated<-RenameIdents(integrated, c('12'="Lowq", 
                                       '8'="Lowq", 
                                       '19'="Lowq",
                                       '21'="Plt", 
                                       '30'="HSCprec", 
                                       '34'="dbl.Mono.B",'Gzmk+Cd8_Tem'="DR+Cd8_Tem", 'Cd8_Tcm'="Gzmk+Cd8_Tem", 'lowq'="Lowq", 'PLT'="Plt", 
                                       'dbl.Mono.T'= "dbl.CD3.Mono",'CD163+CD14+Mono'="CD14+Mono",'Pro'="Cycling T-cells",
                                       'TIM3_TIGIT_Cd8_Tem'="Exhausted_Cd8",'IFN+NK'="CD56dim_NK"))

#assign to sub_clusters
integrated$sub_cluster<-Idents(integrated)
integrated[['ident']]<-Idents(integrated)
DimPlot(integrated)

#remove low quality and doublets cells
clusters_to_remove <- grep("^dbl", unique(Idents(integrated)),value=T)
clusters_to_remove<-unique(c(clusters_to_remove, "Plt", "HSCprec", "Lowq", "lowq","PLT", "Plasmacells"))
keep<-setdiff(unique(Idents(integrated)), clusters_to_remove)
integrated_filt <- subset(integrated, idents = keep)

colors <- c(
  "#E41A1C",  # Red
  "#377EB8",  # Blue
  "#4DAF4A",  # Green
  "#984EA3",  # Purple
  "#FF7F00",  # Orange
  "#FFFF33",  # Yellow
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
  "#FF8C00"   # Dark Orange
)

#UMAP sub_cluster Fig1+ number of patients + number of samples+ number of cells
UMAP <- DimPlot(integrated_filt, label = TRUE, repel = TRUE) + 
  theme(
    axis.ticks = element_blank(),
    axis.text = element_blank()
  ) + 
  labs(x = "Umap 1", y = "Umap 2") +  scale_color_manual(values = colors)+NoLegend()

#save
ggsave("~/Immune_project/new_analysis_sep24/plots/UMAP/UMAP_integrated.pdf", plot=UMAP, width=7, height=4)



#plot Ncells Supplementary Figure 1
ncells<-integrated_filt@meta.data%>%group_by(Sample_ID)%>%
  summarise(n=n())

N_cells<-ggplot(ncells, aes(x=n))+
  geom_histogram(bins=70)+
  theme_classic()+ scale_x_continuous(breaks = c(100,1000,2000,5000,10000))+
  labs(y="", x="", title="Distribution of number of PBMCs per sample")+
  geom_vline(xintercept = 1000, linetype="dotted", size=1)
ggsave("~/Immune_project/new_analysis_sep24/plots/N_cells/Ncells_integrated.pdf", plot=N_cells, width=7, height=4)


#check baseline differences in major cluster
integrated_filt$major_clusters[integrated_filt$major_clusters=="Unk"]<-"DC"
integrated_filt[['ident']]<-integrated_filt$major_clusters
#UMAP major_cluster Fig1+ number of patients + number of samples+ number of cells
integrated_filt[['ident']]<-integrated_filt$major_clusters
UMAP <- DimPlot(integrated_filt, label = TRUE, repel = TRUE) + 
  theme(
    axis.ticks = element_blank(),
    axis.text = element_blank()
  ) + 
  labs(x = "Umap 1", y = "Umap 2") +  scale_color_manual(values = colors)+NoLegend()

#save
ggsave("~/Immune_project/new_analysis_sep24/plots/UMAP/UMAP_MC_integrated.pdf", plot=UMAP, width=7, height=4)

#extract cell counts for MAJOR CLUSTERS 
meta <- read_excel("~/Immune_project/analysis_oct24_fullDB/Meta_Immune_proj.xlsx")
integrated_filt$major_clusters<-droplevels(integrated_filt$major_clusters)
counts<-as.data.frame(table(integrated_filt$Sample_ID, integrated_filt$major_clusters))
colnames(counts)<-c("Sample_ID","major_clusters","n")
total<-counts%>% group_by(Sample_ID)%>% summarise(total=sum(n))
counts<-left_join(counts, total, by="Sample_ID")
counts<-counts%>% mutate(Freq= n/total)
meta<-meta%>%distinct(Sample_ID,.keep_all = TRUE)
counts<-left_join(counts, meta[,c("Sample_ID","BOR","Timepoint", "Therapy","Disease","ID", "Cohort")], by="Sample_ID")
#remove less than 1000cells
counts<-counts%>%filter(total>=1000)
write_csv(counts, "~/Immune_project/new_analysis_sep24/compositional_analysis/MC_PBMC_cell_counts.csv")
#baseline differences between TEC SMM vs TEC RRMM
counts_pre<-counts%>%filter(Timepoint%in%c("Healthy","Pre") )
#19 vs 8
t<-counts_pre%>%group_by(Disease)%>%summarise(n=n_distinct(Sample_ID))
print(t)
#do wilcox.test between RRMM and SMM for each celltype, store.pvalues, correct with BH
counts_pre$Disease<-factor(counts_pre$Disease, levels=c("Healthy","SMM","RRMM"))
mycomparisons<-list(c("Healthy", "SMM"), c("SMM", "RRMM"),c("Healthy","RRMM"))
boxplots_major_cluster_baseline_SMM_RRMM<-ggplot(counts_pre, aes(x=Disease, fill=Disease, y=Freq*100))+geom_boxplot(alpha=0.5, outlier.shape = NA)+
  geom_jitter(width=0.1)+facet_wrap(~major_clusters, scales="free_y",ncol=5)+ labs(y="Percentage of PBMCs")+
  scale_fill_manual(values = c("darkgreen","darkblue", "darkred"))+theme_classic()+
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1))+
  guides(fill=FALSE)
boxplots_major_cluster_baseline_SMM_RRMM<-boxplots_major_cluster_baseline_SMM_RRMM+
  geom_pwc(label.size = 3, method = "wilcox.test",step.increase = 0.06, tip.length = 0.01)

boxplots_major_cluster_baseline_SMM_RRMM<-ggadjust_pvalue(
  boxplots_major_cluster_baseline_SMM_RRMM, p.adjust.method = "BH",
  label = paste("q=","{p.adj.format}"),
  hide.ns = "p",
  
)
ggsave("~/Immune_project/new_analysis_sep24/plots/boxplots/boxplots_major_cluster_baseline_SMM_RRMM.pdf", plot=boxplots_major_cluster_baseline_SMM_RRMM, width=7, height=4)

qsave(integrated_filt, "~/Immune_project/new_analysis_sep24/qs_files/integrated_filt.qs")
sessionInfo()
