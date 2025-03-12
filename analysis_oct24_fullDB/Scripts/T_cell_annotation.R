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

#load object
T_cells<-qread("/mnt/disks/cellranger/Accelerator_qs_files_oct24/Tcells2_merged_reduced_noCARfeat.qs")
#increase resolution
T_cells<-FindClusters(T_cells,resolution=c(1.5))
#use 1.5
#FindMarkers
output_file <- "~/Immune_project/analysis_oct24_fullDB/markers/Tcells_markers_noCAR_res1.5.csv"
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
      slice_head(n = 15)
    
    print(top_markers)
    write_csv(top_markers, output_file)
  }
qsave(T_cells, "/mnt/disks/cellranger/Accelerator_qs_files_oct24/Tcells_noCARfeat_1.5res.qs")


#load markers
Idents(T_cells2)<-T_cells2$seurat_clusters
#Rename T cells

T_cells2<-RenameIdents(T_cells2, c('0'="Gzmb+CD8_Tem", '1'="CD4Naive", '2'="CD4Tem",'3'="Gzmk+CD8_Tem", '4'="Exhausted-like", '5'="Kir+CD8_Tem", '6'="CD8Naive", '7'="CD4_Tem", '8'="CD4_Tem",'9'="Tregs",'10'="Cytotoxic_CD4", '11'="Mait", '12'="IFN+CD4", '13'="CD4Naive",'14'="lowq",'15'="dbl.Mono.CD3", '16'="Proliferating_Tcells",'17'="dbl.Mono.CD3",'18'="Cytotoxic_CD4",'19'="CD4Naive",'20'="lowq",'21'="Tgd",'22'="Proliferating_Tcells",'23'="Proliferating_Tcells",'24'="Gzmb+CD8_Tem",'25'="dbl.Mono.CD3",'26'="dbl.B.CD3",'27'="dbl.B.CD3", '28'="dbl.B.CD3", '29'="dbl.B.CD3",'30'="dbl.B.CD3"))
T_cells2[['ident']]<-Idents(T_cells2)
DimPlot(T_cells2, label=T)
# -------------------------------------------------------------------------
 eleven_vs_six<-FindMarkers(T_cells2, ident.1 = 11, ident.2 = 6, logfc.threshold = 1)

T_cells2[['ident']]<-Idents(T_cells2)
DimPlot(T_cells2)
qsave(T_cells2, "/mnt/disks/cellranger/Accelerator_qs_files_oct24/Tcells_annotated.qs")
#remove doublets
keep<-setdiff(unique(Idents(T_cells2)),grep("dbl*|lowq|31",unique(Idents(T_cells2)), value=TRUE))
T_cells_filt<-subset(T_cells2,idents=keep)

#look at distrubtion of CAR-Tcells[idecabtagene|ciltacabtagene>0]


#look at distribution of Tcell number
ncells<-T_cells_filt@meta.data%>%group_by(Sample_ID)%>%
  summarise(n=n())

N_cells<-ggplot(ncells, aes(x=n))+
  geom_histogram(bins=70)+
  theme_classic()+ scale_x_continuous(breaks = c(100,500,1000,2000,5000))+
  labs(y="", x="", title="Distribution of number of T-cells per sample")+
  geom_vline(xintercept = 100, linetype="dotted", size=1)


meta <- read_excel("~/Immune_project/analysis_oct24_fullDB/Meta_Immune_proj.xlsx")
# Remove duplicates based on the Sample_ID column
meta <- meta %>% distinct(Sample_ID, .keep_all = TRUE)

T_cells_filt$sub_cluster<-Idents(T_cells_filt)
counts<-as.data.frame(table(T_cells$Sample_ID, T_cells$seurat_clusters))
colnames(counts)<-c("Sample_ID","sub_cluster","n")
total<-counts%>% group_by(Sample_ID)%>% summarise(total=sum(n))
counts<-left_join(counts, total, by="Sample_ID")
counts<-counts%>% mutate(Freq= n/total)
counts<-left_join(counts, meta[,c("Sample_ID","BOR","Timepoint", "Therapy","Disease","ID","Tissue")], by="Sample_ID")
counts<-counts%>%filter(total>100)

counts$Disease<-factor(counts$Disease, levels=c("Healthy","SMM","RRMM"))
counts$Timepoint<-factor(counts$Timepoint, levels=c("Healthy","Pre","Post"))
counts$BOR<-factor(counts$BOR, levels=c("CR","VGPR","PR","PD"))
counts$response[counts$Disease=="SMM"&counts$Therapy%in%c("TEC","Len")]<-"R"
counts$response[counts$Disease=="RRMM"&counts$BOR%in%c("CR","VGPR")]<-"R"
counts$response[counts$Disease=="RRMM"&counts$BOR%in%c("PR","SD","MR","PD")]<-"NR"
counts$response[counts$Disease=="Healthy"]<-"Healthy"
counts$response<-factor(counts$response, levels=c("Healthy","R","NR"))
counts<-counts%>%mutate(response_disease=paste0(Disease,"_",response))
counts$response_disease[counts$response_disease=="Healthy_Healthy"]<-"Healthy"
counts$response_disease<-factor(counts$response_disease, levels=c("Healthy","SMM_R","RRMM_R","RRMM_NR"))




#baseline differences between TEC SMM vs TEC RRMM select only samples w/>1000 cells
Tcells_counts<-counts
#add levels
Tcells_counts$Disease<-factor(Tcells_counts$Disease, levels=c("Healthy","SMM","RRMM"))
Tcells_counts$Timepoint<-factor(Tcells_counts$Timepoint, levels=c("Healthy","Pre","Post"))
Tcells_counts$BOR<-factor(Tcells_counts$BOR, levels=c("CR","VGPR","PR","PD"))
Tcells_counts$response_disease<-factor(Tcells_counts$response_disease, levels=c("Healthy","SMM_R","RRMM_R","RRMM_NR"))

#by filtering "Pre" I am selecting BL MM samples, counts here is Tcells>100
MM_BL_TCELLS<-Tcells_counts%>%filter(Timepoint%in%c("Pre"))
#check number of samples x disease x timepoint
MM_BL_TCELLS%>%group_by(Disease, Timepoint)%>%summarise(n=n_distinct(Sample_ID))

#do wilcox.test between RRMM and SMM for each celltype, store.pvalues, correct with BH
celltypes<-unique(MM_BL_TCELLS$sub_cluster)
results <- list()  
for (i in seq_along(celltypes)) {  
  sub <- MM_BL_TCELLS %>% filter(sub_cluster == celltypes[i])  
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
MM_BL_TCELLS$Disease<-factor(MM_BL_TCELLS$Disease, levels=c("SMM","RRMM"))


#extract mean values for each subtype across disease groups
input_dir<-"~/Immune_project/new_analysis_sep24/plots/"
means<-MM_BL_TCELLS%>% group_by(Disease, sub_cluster)%>% summarise(mean_value=mean(Freq))
means_wide<-means%>% pivot_wider(names_from = Disease, values_from = mean_value)
means_wide<-means_wide%>%mutate(log2FC=log2(RRMM/SMM))
means_wide<-left_join(means_wide, adjusted_p_values, by="sub_cluster")
means_wide<-means_wide%>%mutate(enrichment= ifelse(log2FC>0,"RRMM","SMM"))
means_wide<-means_wide%>%mutate(significant= ifelse(p_adj<0.1,"q<0.1","q>0.1"))
#remove underscore for the graph but keep means_wide for further analysis
means_wide2<-means_wide
means_wide2$sub_cluster <- gsub("_", "", means_wide2$sub_cluster)
logFCmeans <- ggplot(means_wide2, aes(x=log2FC, y=-log10(pval))) +
  geom_point(aes(color=enrichment, shape=significant), size=3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "black") +
  geom_vline(xintercept=0, linetype="dotted", color="black")+
  scale_color_manual(values = c("darkblue", "darkred")) +
  scale_shape_manual(values = c("q<0.1" = 24, "q>0.1" = 21)) + # Adjust based on your data
  geom_text_repel(aes(label = ifelse(pval < 0.05, sub_cluster, "")), 
                  size = 3, max.overlaps = Inf) +
  labs(x="Log 2 Fold Change[RRMM/SMM]")+
  theme_classic()

pdf(paste0(input_dir, "logFC_mean_diff_Tcells.pdf"), width=7, height = 4)
logFCmeans
dev.off()


T_cells_clusters<-means_wide%>%filter(p_adj<0.1 &
                                        sub_cluster%in% grep("*Naive|Tgd|Exhausted*|*Tem|Cycling|Tregs|Tcm",means_wide$sub_cluster, value=TRUE))
T_cells_clusters<-unique(T_cells_clusters$sub_cluster)
print(T_cells_clusters)

BL_HD_MM_sign_Tcells<-Tcells_counts%>%filter(sub_cluster%in%T_cells_clusters & Timepoint%in%c("Pre", "Healthy"))

#Plot porportion of PBMC clusters showing diff between SMM and MM with absolute values and HD reference
Sign_T_cells_REL_BL<-ggplot(BL_HD_MM_sign_Tcells, aes(x=Disease, fill=Disease, y=Freq*100))+geom_boxplot(alpha=0.5, outlier.shape = NA)+
  geom_jitter(width=0.1)+facet_wrap(~sub_cluster, scales="free_y",ncol=7)+ labs(y="Percentage of T-cells")+
  scale_fill_manual(values = c("darkgreen","darkblue", "darkred"))+theme_classic()+guides(fill=FALSE)+
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1))
Sign_T_cells_REL_BL<-Sign_T_cells_REL_BL+ geom_pwc(label.size = 3, method = "wilcox.test",step.increase = 0.06, tip.length = 0.01)

Sign_T_cells_REL_BL<-ggadjust_pvalue(
  Sign_T_cells_REL_BL, p.adjust.method = "BH",
  label = paste("q=","{p.adj.format}"),
  hide.ns = "p"
)

BL_MM_sign_Tcells<-BL_HD_MM_sign_Tcells%>%filter(Disease%in%c("SMM","RRMM")&Therapy%in%c("TEC","Len"))
#Plot porportion of PBMC clusters showing diff between SMM and MM with absolute values and HD reference
Sign_T_cells_REL_BL_RX<-ggplot(BL_MM_sign_Tcells, aes(x=response_disease, fill=response_disease, y=Freq*100))+
  geom_boxplot(alpha=0.5, outlier.shape = NA)+
  geom_jitter(width=0.1)+facet_wrap(~sub_cluster, scales="free_y",ncol=7)+ labs(y="Percentage of Tcells")+
  scale_fill_manual(values = c("darkblue","gold","lightgrey"))+theme_classic()+guides(fill=FALSE)+
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1))
Sign_T_cells_REL_BL_RX<-Sign_T_cells_REL_BL_RX+ geom_pwc(label.size = 3, method = "wilcox.test",step.increase = 0.06, tip.length = 0.01)

Sign_T_cells_REL_BL_RX<-ggadjust_pvalue(
  Sign_T_cells_REL_BL_RX, p.adjust.method = "BH",
  label = paste("q=","{p.adj.format}"),
  hide.ns = "p"
)

