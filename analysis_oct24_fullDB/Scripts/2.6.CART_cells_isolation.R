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
clusters_to_remove <- grep("^dbl", unique(Idents(T_cells)),value=TRUE)
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

T_cells_filt$sub_cluster<-Idents(T_cells_filt)

#I want to extract how many cells have at least one copy of car-construct in each sample
# Add a metadata column for ciltacabtagene expression
T_cells_filt$express_cilta <- ifelse(
  GetAssayData(T_cells_filt, slot = "counts")["ciltacabtagene", ] > 0,
  1,
  0
)

# Add a metadata column for idecabtagene expression
T_cells_filt$express_ide<- ifelse(
  GetAssayData(T_cells_filt, slot = "counts")["idecabtagene", ] > 0,
  1,
  0
)

# Verify the new metadata columns
head(T_cells_filt@meta.data[, c("Sample_ID", "express_cilta", "express_ideca")])
df_cilta<-as.data.frame(table(T_cells_filt@meta.data$Sample_ID, T_cells_filt@meta.data$express_cilta))
df_cilta_wide<-df_cilta%>%pivot_wider(names_from = Var2, values_from = Freq)
colnames(df_cilta_wide)<-c("Sample_ID","CAR-","CAR+")
df_cilta_wide<-left_join(df_cilta_wide, meta, by="Sample_ID")
df_cilta_wide<-df_cilta_wide%>%filter(!Cohort=="Calgary_RRMM")
df_cilta_wide<-df_cilta_wide%>%mutate(`Post-Cilta Sample`=ifelse(Therapy=="CART"&!Timepoint%in%c("Pre","Healthy")&Product=="Cilta-cel","Yes","No"))
df_cilta_wide<-df_cilta_wide%>%mutate(`Cilta Positive Samples`=ifelse(`CAR+`>0,1,0))

df_cilta_wide$`Post-Cilta Sample` <- factor(df_cilta_wide$`Post-Cilta Sample`, levels = c("No","Yes"))
df_cilta_wide <- df_cilta_wide %>%
  arrange(`Post-Cilta Sample`, `CAR+` )
df_cilta_wide$Sample_ID <- factor(df_cilta_wide$Sample_ID, levels = unique(df_cilta_wide$Sample_ID))

p<-ggplot(df_cilta_wide, aes(x = `Post-Cilta Sample`, fill = factor(`Cilta Positive Samples`, labels = c("Negative", "Positive")))) +
  geom_bar(position = position_fill(), alpha=0.5) +
  labs(
    title = "Proportion of Samples\n with Cilta+cells",
    x = "Post-Cilta Samples",
    y = "Percentage",
    fill = "Cilta Positive Samples"
  ) +
  scale_y_continuous(labels = scales::percent) +
  theme_classic()+
  theme(text = element_text(size = 50), legend.position = "bottom") +
  scale_fill_manual(values = c("grey","purple"))
ggsave(paste0(output_dir_car, "false_pos_cilta.pdf"),plot = p ,width=15, height=13)

df_cilta_wide<-df_cilta_wide%>%filter(`CAR+`>0)
output_dir_car<-"~/Immune_project/full_dataset_analysis/plots/CAR/"
dir.create(output_dir_car)
# 6) Plot 1: bar plot of median expression
p <- ggplot(df_cilta_wide, aes(x = Sample_ID, y = `CAR+`, fill = `Post-Cilta Sample`)) +
  geom_col() +
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),text = element_text(size = 50), legend.position = "bottom") +
  labs(
    y = paste0("Number of CAR+"),
    x = "Samples",
    fill = "Post-Cilta Sample\n(previously untreated pts)"
  ) +
  scale_fill_manual(
    values = c(
      "No"          = "forestgreen",
      "Yes"           = "darkblue"  ))

ggsave(paste0(output_dir_car,"Cilta_abs.pdf"), plot=p, width=18, height=15)


# -IDE-------------------------------------------------------------------
df_ide<-as.data.frame(table(T_cells_filt@meta.data$Sample_ID, T_cells_filt@meta.data$express_ide))
df_ide_wide<-df_ide%>%pivot_wider(names_from = Var2, values_from = Freq)
colnames(df_ide_wide)<-c("Sample_ID","CAR-","CAR+")
df_ide_wide<-left_join(df_ide_wide, meta, by="Sample_ID")
df_ide_wide<-df_ide_wide%>%filter(!Cohort=="Calgary_RRMM")
df_ide_wide<-df_ide_wide%>%mutate(`Post-ide Sample`=ifelse(Therapy=="CART"&!Timepoint%in%c("Pre","Healthy")&Product=="Ide-celcel","Yes","No"))
df_ide_wide<-df_ide_wide%>%mutate(`ide Positive Samples`=ifelse(`CAR+`>0,1,0))

df_ide_wide$`Post-ide Sample` <- factor(df_ide_wide$`Post-ide Sample`, levels = c("No","Yes"))
df_ide_wide <- df_ide_wide %>%
  arrange(`Post-ide Sample`, `CAR+` )
df_ide_wide$Sample_ID <- factor(df_ide_wide$Sample_ID, levels = unique(df_ide_wide$Sample_ID))

p<-ggplot(df_ide_wide, aes(x = `Post-ide Sample`, fill = factor(`ide Positive Samples`, labels = c("Negative", "Positive")))) +
  geom_bar(position = position_fill(), alpha=0.5) +
  labs(
    title = "Proportion of Samples\n with Ide+cells",
    x = "Post-Ide Samples",
    y = "Percentage",
    fill = "Ide Positive Samples"
  ) +
  scale_y_continuous(labels = scales::percent) +
  theme_classic()+
  theme(text = element_text(size = 50), legend.position = "bottom") +
  scale_fill_manual(values = c("grey","purple"))
ggsave(paste0(output_dir_car, "false_pos_ide.pdf"), width=15, height=13)

df_ide_wide<-df_ide_wide%>%filter(`CAR+`>0)
output_dir_car<-"~/Immune_project/full_dataset_analysis/plots/CAR/"
dir.create(output_dir_car)
# 6) Plot 1: bar plot of median expression
p <- ggplot(df_ide_wide, aes(x = Sample_ID, y = `CAR+`, fill = `Post-ide Sample`)) +
  geom_col() +
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),text = element_text(size = 50), legend.position = "bottom") +
  labs(
    y = paste0("Number of CAR+"),
    x = "Samples",
    fill = "Post-Ide Sample\n(previously untreated pts)"
  ) +
  scale_fill_manual(
    values = c(
      "No"          = "forestgreen",
      "Yes"           = "darkblue"  ))

ggsave(paste0(output_dir_car,"ide_abs.pdf"), plot=p, width=18, height=15)


# -------------------------------------------------------------------------
#if Cilta>1 Positive Sample
# -------------------------------------------------------------------------
CAR<-subset(T_cells_filt, subset=Timepoint%in%c("Post","Day14","Day30","Day100","Day56","6Months")&Therapy=="CART")
#BUILD SWIMMER
RRMM<-read_xlsx("/mnt/disks/disk/full_dataset_qs/swimmer_RRMM_CART.xlsx")
RRMM<-RRMM%>%dplyr::select(c("MRN","Sample date","sample ID", "Timepoint"))
colnames(RRMM)<-c("ID","Date","Sample_ID", "Timepoint")
t<-RRMM%>%group_by(ID)%>%summarise(n=n_distinct(Timepoint))%>%filter(n==2)
ids<-unique(t$ID)
RRMM<-RRMM%>%filter(ID%in%ids)
df2<-RRMM%>%group_by(ID)%>%dplyr::select(c(ID,Date, Timepoint))%>%pivot_wider(names_from = Timepoint, values_from = Date )
df2<-df2%>%mutate(follow_up=difftime( Post, Pre, units = "days"))
df2$follow_up<-gsub(" days","",df2$follow_up)
df2$follow_up<-as.numeric(df2$follow_up)
df2<-df2%>%dplyr::select("ID","follow_up")%>%mutate(ID_followup=paste0(ID, follow_up))
df2$Tissue<-"PB"
df2$Timepoint<-df2$follow_up
df2$Disease<-"RRMM"
df2$ID<-as.character(df2$ID)
df2<-df2%>%mutate(Product=ifelse(ID=="756436","Ide-cel","Cilta-cel"))

#HRSMM
meta<-CAR@meta.data
meta<-meta%>%filter(Product=="Cilta-cel")
meta <- meta %>%filter(Disease=="HRSMM")%>%
  mutate(follow_up = case_when(
    Timepoint == "Day14"    ~ 14,
    Timepoint == "Day30"    ~ 30,
    Timepoint == "Day56"    ~ 60,
    Timepoint == "Day100"   ~ 100,
    Timepoint == "6Months"  ~ 180,  # Approssimazione di 6 mesi in giorni
    Timepoint == "Post"      ~ NA_real_,  # Se "Post" non ha un valore numerico specifico
    TRUE                     ~ NA_real_   # Per eventuali valori non previsti
  ))

df3<-meta%>%mutate(ID_follow_up=paste0(ID,follow_up))%>%dplyr::select(c(ID,follow_up, ID_follow_up, Tissue, Timepoint))%>%distinct(ID_follow_up,.keep_all = TRUE)
df3$Timepoint[df3$Timepoint=="6Months"]<-"Day180"
df3$Timepoint<-gsub("Day","",df3$Timepoint)
df3$Timepoint<-as.numeric(df3$Timepoint)
df3$Disease<-"HRSMM"
df3$Product<-"Cilta-cel"

meta<-T_cells_filt@meta.data
meta <- meta %>%filter(Cohort=="MXMERZ002A"&Timepoint=="Post")
rade_metadata<-read_xlsx("/mnt/disks/disk/full_dataset_qs/rade_metadata.xlsx")
rade_metadata<-rade_metadata%>%filter(follow_up>0)
meta$ID<-gsub("M","",meta$ID)
meta<-left_join(meta,rade_metadata, by=c("ID","Tissue"))
df4<-meta%>%mutate(ID_follow_up=paste0(ID,follow_up,Tissue))%>%dplyr::select(c(ID,follow_up, ID_follow_up, Tissue, Timepoint,Product))%>%distinct(ID_follow_up,.keep_all = TRUE)

df4$Timepoint<-df4$follow_up
df4$Timepoint<-as.numeric(df4$Timepoint)
df4$Disease<-"RRMM(Rade et al.)"
df4$Product[df4$Product=="Ide-celcel"]<-"Ide-cel"


df5<-rbind(df2, df3,df4)
# Assegnare un ordine agli eventi di follow-up per ogni paziente
library(dplyr)

# Assicurarsi che 'Disease' sia un fattore con ordine desiderato
df5 <- df5 %>%
  mutate(Disease = factor(Disease, levels = c("HRSMM", "RRMM", "RRMM(Rade et al.)")))

# Ordinare il dataframe
df5 <- df5 %>%
  arrange(Disease, ID, Timepoint) %>%
  group_by(ID) %>%
  mutate(EventOrder = row_number()) %>%
  ungroup()

# Calcolare Start e End per ogni paziente
patient_limits <- df5 %>%
  group_by(ID, Disease,Product, follow_up) %>%
  summarise(
    Start = min(Timepoint),
    End = max(Timepoint),
    .groups = 'drop'
  ) %>%
  arrange(Disease, desc(follow_up)) %>%
  mutate(ID = factor(ID, levels = unique(ID)))

# Unire i limiti temporali al dataframe originale
df5 <- df5 %>%
  left_join(patient_limits, by = c("ID", "Disease"))


library(ggplot2)

# Creare il Swimmer Plot
swimmer_car<-ggplot() +
  # Linee orizzontali per ogni paziente, colorate per Disease e linetype per Product
  geom_segment(
    data = patient_limits,
    aes(x = 0, xend = End, y = factor(ID), yend = factor(ID), 
        color = Disease, linetype = Product),
    size = 10,              # Ridurre lo spessore per aumentare la visibilità dei punti
    alpha = 0.8            # Aggiungere trasparenza per evitare sovrapposizioni pesanti
  ) +
  
  # Punti per ogni follow-up timepoint
  geom_point(
    data = df5,
    aes(x = Timepoint, y = factor(ID), fill = Tissue),
    size = 15, shape=21
  ) +
  
  # Aggiungere etichette e titolo
  labs(
    title = "",
    x = "Days from Infusion",
    y = "Patient ID",
    color = "Disease / Tissue",
    linetype = "Product"
  ) +
  
  # Temi e personalizzazioni
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 40),
    legend.text = element_text(size = 30),
    text = element_text(size=40)
  ) +
  
  # Ordinare l'asse y per Disease (HRSMM prima, poi RRMM)
  scale_y_discrete(limits = rev(levels(patient_limits$ID))) +
  
  # Definire i colori manuali per Disease e Tissue
  scale_color_manual(
    name = "Disease",
    values = c(
      "HRSMM" = "steelblue", 
      "RRMM" = "firebrick", 
      "RRMM(Rade et al.)" = "orange"
    ),
    breaks = c("HRSMM", "RRMM", "RRMM(Rade et al.)"),
    labels = c("HRSMM", "RRMM", "RRMM (Rade et al.)")
  ) +
  scale_fill_manual(
    name = "Tissue",
    values = c(
      "PB" = "skyblue", 
      "BM" = "purple"
    ),
    breaks = c( "PB", "BM"),
    labels = c( "PB", "BM")
  ) +
  
  # Definire le linetypes manuali per Product
  scale_linetype_manual(
    name = "Product",
    values = c("Ide-cel" = "dotted", "Cilta-cel" = "solid")
  ) +
  
  # Aggiungere una linea verticale all'origine (opzionale)
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  
  # Impostare i tick e le etichette sull'asse X a giorni specifici
  scale_x_continuous(
    breaks = c(30, 60, 100, 200,300),
    expand = expansion(mult = c(0, 0.1))
  ) +
  
  # Migliorare la visibilità delle legende
  guides(
    color = guide_legend(order = 1, title.position = "top", title.hjust = 0.5),
    linetype = guide_legend(order = 2, title.position = "top", title.hjust = 0.5)
  )
ggsave(filename="~/Immune_project/full_dataset_analysis/plots/CAR/swimmer.pdf", plot=swimmer_car, width=20, height=18)

View(patient_limits)

# -------------------------------------------------------------------------
#Proportion of CAR-T cells in T-cells at Day 30 HRSMM vs RRMM
# -------------------------------------------------------------------------

CAR<-subset(CAR, idents=keep)
#TAKE Post samples of SMM and make a plot of longitudinal CAR tracking
Cilta_SMM<-subset(CAR, Product=="Cilta-cel"&Disease=="HRSMM")
# Get unique IDs across both PB and BM and assign letters consistently
unique_ids <- unique(Cilta_SMM@meta.data$ID)  
id_mapping <- setNames(LETTERS[seq_along(unique_ids)], unique_ids)  # Assign A, B, C, D...

# Apply the same mapping to the full dataset
Cilta_SMM@meta.data <- Cilta_SMM@meta.data %>%
  mutate(ID2 = id_mapping[ID])
meta<-meta%>%
  mutate(ID2 = id_mapping[ID])

cilta_df<-as.data.frame(table(Cilta_SMM$Sample_ID, Cilta_SMM$express_cilta))
colnames(cilta_df)<-c("Sample_ID","CAR_status","Number")
cilta_wide<-cilta_df%>%pivot_wider(names_from = CAR_status, values_from = Number)
colnames(cilta_wide)<-c("Sample_ID","T_cells","CART")
cilta_wide<-cilta_wide%>%mutate(cilta_positive=ifelse(CART>1, CART,"0"))
cilta_wide<-left_join(cilta_wide, meta, by="Sample_ID")


cilta_wide_pb<-cilta_wide%>%filter(Tissue=="PB")
cilta_wide_pb$cilta_positive<-as.numeric(cilta_wide_pb$cilta_positive)
cilta_wide_pb<-cilta_wide_pb%>%mutate(total=T_cells+cilta_positive)%>%
  mutate(CART_prop=cilta_positive/total)%>%
  mutate(cilta_binary=ifelse(cilta_positive>0,"yes","no"))

#add day0 
# Create Day0 entries with CART_prop = 0 for each Sample_ID
day0_df <- cilta_wide_pb %>%
  distinct(ID2, .keep_all = TRUE) %>%  # Keep unique Sample_IDs
  mutate(Timepoint = "Day0", CART_prop = 0, cilta_positive = 0, total = T_cells)  # Set values for Day0

# Bind the new Day0 rows to the original dataframe
cilta_wide_pb <- bind_rows(cilta_wide_pb, day0_df) %>%
  arrange(ID2, Timepoint)  # Ensure correct ordering

# Function to compute binomial confidence intervals
compute_ci <- function(successes, total, conf.level = 0.95) {
  if (total == 0) return(c(0, 0))  # Avoid division by zero
  test <- binom.test(successes, total, conf.level = conf.level)
  return(c(test$conf.int[1], test$conf.int[2]))
}

# Add confidence intervals for CART proportion
cilta_wide_pb <- cilta_wide_pb %>%
  rowwise() %>%
  mutate(
    CI = list(compute_ci(cilta_positive, total)),
    CART_lower = CI[1],
    CART_upper = CI[2]
  ) %>%
  ungroup() %>%
  dplyr::select(-CI)  # Remove list column


p<-ggplot(cilta_wide_pb, aes(x = Timepoint, y = CART_prop, group = ID2, color = ID2)) +
  geom_line(size=3) +
  geom_point(size=10) +
  geom_errorbar(aes(ymin = CART_lower, ymax = CART_upper), width = 0.2, size=2) +
  theme_classic() +
  labs(title = "",
       y = "Tcells(%)",
       x = "Timepoint") +
  theme(text = element_text(size = 50), legend.position = "bottom")
ggsave(plot=p, filename="~/Immune_project/full_dataset_analysis/plots/CAR/long_SMM.pdf", width=15, height=13)


# -------------------------------------------------------------------------
#BM
# -------------------------------------------------------------------------


cilta_wide_bm<-cilta_wide%>%filter(Tissue=="BM")
cilta_wide_bm$cilta_positive<-as.numeric(cilta_wide_bm$cilta_positive)
cilta_wide_bm<-cilta_wide_bm%>%mutate(total=T_cells+cilta_positive)%>%
  mutate(CART_prop=cilta_positive/total)%>%
  mutate(cilta_binary=ifelse(cilta_positive>0,"yes","no"))

#add day0 
# Create Day0 entries with CART_prop = 0 for each Sample_ID
day0_df <- cilta_wide_bm %>%
  distinct(ID2, .keep_all = TRUE) %>%  # Keep unique Sample_IDs
  mutate(Timepoint = "Day0", CART_prop = 0, cilta_positive = 0, total = T_cells)  # Set values for Day0

# Bind the new Day0 rows to the original dataframe
cilta_wide_bm <- bind_rows(cilta_wide_bm, day0_df) %>%
  arrange(ID2, Timepoint)  # Ensure correct ordering

# Function to compute binomial confidence intervals
compute_ci <- function(successes, total, conf.level = 0.95) {
  if (total == 0) return(c(0, 0))  # Avoid division by zero
  test <- binom.test(successes, total, conf.level = conf.level)
  return(c(test$conf.int[1], test$conf.int[2]))
}

# Add confidence intervals for CART proportion
cilta_wide_bm <- cilta_wide_bm %>%
  rowwise() %>%
  mutate(
    CI = list(compute_ci(cilta_positive, total)),
    CART_lower = CI[1],
    CART_upper = CI[2]
  ) %>%
  ungroup() %>%
  dplyr::select(-CI)  # Remove list column

cilta_wide_bm$Timepoint<-factor(cilta_wide_bm$Timepoint, levels=c("Day0"   ,  "Day30" ,  "Day56", "Day100", "6Months"))    
p<-ggplot(cilta_wide_bm, aes(x = Timepoint, y = CART_prop, group = ID2, color = ID2)) +
  geom_line(size=3) +
  geom_point(size=10) +
  geom_errorbar(aes(ymin = CART_lower, ymax = CART_upper), width = 0.2, size=2) +
  theme_classic() +
  labs(title = "",
       y = "Tcells(%)",
       x = "Timepoint") +
  theme(text = element_text(size = 50), legend.position = "bottom")
ggsave(plot=p, filename="~/Immune_project/full_dataset_analysis/plots/CAR/long_SMM_BM.pdf", width=15, height=13)



# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
Cilta_SMM<-subset(Cilta_SMM, express_cilta=="1"&Disease=="HRSMM")
Cilta_SMM@meta.data<-Cilta_SMM@meta.data%>%dplyr::select(-c("ID2"))
Cilta_SMM@meta.data<-left_join(Cilta_SMM@meta.data, meta[c("ID","ID2")], by="ID")
rownames(Cilta_SMM@meta.data)<-Cilta_SMM@meta.data$new_barcodes
t<-as.data.frame(table(Cilta_SMM$Sample_ID))
t<-t%>%filter(Freq>1)
ids<-unique(t$Var1)
Cilta_SMM_filt<-subset(Cilta_SMM, Sample_ID%in%ids)
#look only at those with follow_up<70
# filt<-df5%>%filter(follow_up<120)
# ids<-unique(filt$ID)
Cilta_SMM_filt@meta.data<-Cilta_SMM_filt@meta.data%>%mutate(TP_Tissue=paste0(Tissue,"_",Timepoint))
Cilta_SMM_filt@meta.data$TP_Tissue<-factor(Cilta_SMM_filt@meta.data$TP_Tissue,levels=c("PB_Day14","PB_Day30","BM_Day30","BM_Day100"))


CD4<-grep("CD4+|Th|Treg", unique(Cilta_SMM_filt$sub_cluster), value=TRUE) 
CD8<-setdiff(unique(Cilta_SMM_filt$sub_cluster), CD4)
Cilta_SMM_filt@meta.data<-Cilta_SMM_filt@meta.data%>%mutate(lineage=ifelse(sub_cluster%in%CD4,"CD4+","CD8+"))

# Compute cell counts for each combination of TP_Tissue, ID, and lineage
cell_counts <- Cilta_SMM_filt@meta.data %>%
  group_by(TP_Tissue, ID) %>%
  summarise(count = n(), .groups = "drop")

p<-ggplot(Cilta_SMM_filt@meta.data, aes(x = TP_Tissue, fill = lineage)) +
  geom_bar(position = position_fill(), color = "black", alpha = 0.5) +
  geom_text(data = cell_counts, aes(x = TP_Tissue, y = 1, label = count), 
            inherit.aes = FALSE, vjust = 1, size = 5) +  # Fix aesthetics issue
  facet_wrap(~ID, scales = "free") +
  labs(y="Percentage of CART")+
  theme_classic()+scale_fill_manual(values=c("steelblue","darkorange"))
ggsave(paste0(output_dir_car, "CD4_CD8_cilta.pdf"),plot = p ,width=9, height=6)

ggplot(Cilta_SMM_filt@meta.data, aes(x=TP_Tissue, fill=sub_cluster))+
  geom_bar(position=position_fill(), color="black", alpha=0.5)+facet_wrap(~ID, scales="free")+theme_classic()+
  
  
  
ggplot(cilta_wide_pb, aes(x=Disease, y=CART_prop, color=Disease))+geom_point()+geom_boxplot()+stat_compare_means()
FeaturePlot(Cilta, "ciltacabtagene", split.by = "Disease")
Cilta_filt<-subset(Cilta, express_cilta=="1")
ggplot(Cilta_filt@meta.data, aes(x=Timepoint, fill=sub_cluster))+geom_bar(position=position_fill())+facet_wrap(~ID)
Cilta_filt<-subset(Cilta_filt, Timepoint%in%c("Post","Day30"))

Cilta_filt$sub_cluster<-droplevels(Cilta_filt$sub_cluster)
counts<-as.data.frame(table(Cilta_filt$Sample_ID, Cilta_filt$sub_cluster))
colnames(counts)<-c("Sample_ID","sub_cluster","n")
total<-counts%>% group_by(Sample_ID)%>% summarise(total=sum(n))
counts<-left_join(counts, total, by="Sample_ID")
counts<-counts%>% mutate(Freq= n/total)
meta<-meta%>%distinct(Sample_ID,.keep_all = TRUE)
counts<-left_join(counts, meta[,c("Sample_ID","BOR","Timepoint", "Tissue","Therapy","Disease","ID", "Cohort")], by="Sample_ID")

View(counts)










# -------------------------------------------------------------------------

CAR<-subset(CAR, express_cilta=="1")
df<-as.data.frame(table(CAR$Sample_ID))
df<-df%>%filter(Freq>10)
ids<-unique(df$Var1)
CAR_filt<-subset(CAR, Sample_ID%in%ids)
CAR_filt@meta.data%>%group_by(Disease, Therapy, Timepoint, Tissue)%>%summarise(n=n_distinct(Sample_ID))

CAR_filt<-subset(CAR_filt, Timepoint%in%c("Day30","Post"))
table(CAR_filt$ID,CAR_filt$Timepoint)
pseudo_cd8 <- AggregateExpression(CAR_filt, assays = "RNA", return.seurat = TRUE, 
                                  group.by = c("ID","Timepoint", "Disease"))

output_dir<-paste0("/mnt/disks/disk/full_dataset_qs/DEG/CART/")
dir.create(output_dir)
pseudo_cd8@assays$RNA$counts<-round(pseudo_cd8@assays$RNA$counts)
    Idents(pseudo_cd8) <- "Disease"
    bulk.de <- FindMarkers(object = pseudo_cd8, 
                           ident.1 ="RRMM", 
                           ident.2 = "HRSMM",
                           test.use = "DESeq2")




cilta_wide_bm<-cilta_wide%>%filter(Tissue=="BM")
cilta_wide_bm$cilta_positive<-as.numeric(cilta_wide_bm$cilta_positive)
cilta_wide_bm<-cilta_wide_bm%>%mutate(total=T_cells+cilta_positive)%>%
  mutate(CART_prop=cilta_positive/total)

#look only at those with follow_up<70
ggplot(cilta_wide_bm, aes(x=Timepoint, y=CART_prop, color=Disease))+geom_point()+geom_line(aes(group=ID))
FeaturePlot(Cilta, "ciltacabtagene", split.by = "Disease")
