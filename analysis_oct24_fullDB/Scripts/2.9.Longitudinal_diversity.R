
output_dir <- "/mnt/disks/disk/full_dataset_qs/diversity/"

#load meta
meta <- read_excel("~/Immune_project/full_dataset_analysis/Meta_Immune_proj.xlsx")
meta <- meta%>%
  mutate(response_disease = case_when(
    BOR %in% c("PR", "PD") & Disease == "RRMM" ~ "RRMM_NR",            # RRMM Non-Responder
    Disease == "RRMM" & !BOR %in% c("PR", "PD") ~ "RRMM_R",            # RRMM Responder
    Disease == "HRSMM" & !Therapy== "Len" ~ "HRSMM_R",                     # SMM Responder
    Disease == "Healthy" ~ "Healthy_NA",                              # Healthy: NA response
    TRUE ~ "SMM_NA"                                                    # Default: SMM Non-Applicable
  ))
meta<-meta%>%distinct(Sample_ID,.keep_all = TRUE)

out_mean<-read_csv(paste0(output_dir, "Diversity_Means_jan2025.csv"))
#ADD METADATA<- -2 sampples with less than 100 cells with clonotypes
out_mean<-left_join(out_mean, meta[,c("Sample_ID","response_disease","BOR","Timepoint","Tissue", "Therapy","Disease","ID")], by="Sample_ID")

df<-out_mean%>%filter(Therapy%in%c("TEC", "TALQ","CART")&Tissue=="PB")%>%group_by(ID, Therapy, Disease)%>%summarise(n=n_distinct(Sample_ID))
paired_ids<-df%>%filter(n>1)
paired_ids<-unique(paired_ids$ID)
paired_counts<-out_mean%>%filter(Tissue=="PB"&ID%in%paired_ids)
df2<-paired_counts%>%filter(Therapy%in%c("TEC", "TALQ")&Tissue=="PB")%>%group_by(ID)%>%summarise(n=n_distinct(Sample_ID))


# -------------------------------------------------------------------------
#HRSMM diversity longitudinal
# -------------------------------------------------------------------------
paired_counts_TEC<-paired_counts%>%filter(Therapy=="TEC")
paired_counts_TEC$Timepoint<-factor(paired_counts_TEC$Timepoint, levels=c("Pre","Post"))
chao_long<-ggplot(paired_counts_TEC, aes(x=Timepoint, y=1/Chao, fill=Timepoint))+
  geom_violin(alpha=0.5) +
  geom_boxplot(size=3, alpha=0.6, outlier.shape = NA) +
  geom_line(aes(group=ID),size=0.8) +
  geom_jitter(size=10, width=0.1,aes(fill=Timepoint, stroke=2),color="black",shape=21, alpha=0.5)+
  theme_classic()+guides(color=FALSE)+
  scale_fill_manual(values = c("gray","gold"))+
  facet_wrap(~Disease, scales="free")+
  theme(plot.margin = unit(c(0.5, 0, 0.2, 0.1), "cm"),title=element_text(size=40),
        plot.title = element_text(hjust=0.5, face = "bold"), text=element_text(size=50,color = "black"), axis.text.x =element_text(size=50,color = "black"),
        axis.text.y =element_text(size=40,color = "black"), axis.title.y =element_text(size=50,color = "black"), legend.position = "bottom")+
  labs(title="", y="1/Chao Index", x="")+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))+stat_compare_means(method = "wilcox.test", paired=TRUE, size=10, label="p.format")

output_dir<-"~/Immune_project/full_dataset_analysis/plots/diversity/"
ggsave(paste0(output_dir,"chao_long.pdf"),
       plot=chao_long, width=13, height=12 )


# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
T_cells_filt_clono <- qread("/mnt/disks/disk/full_dataset_qs/diversity/T_cells_clono.qs")

sub <- T_cells_filt_clono@meta.data
#TAKE ONLY BASELINE PB
df<-sub%>%filter(Therapy%in%c("TEC")&Tissue=="PB"&Disease=="HRSMM")%>%group_by(ID, Therapy, Disease)%>%summarise(n=n_distinct(Sample_ID))
paired_ids<-df%>%filter(n>1)
paired_ids<-unique(paired_ids$ID)

sub<-sub%>%filter(ID%in% paired_ids&Therapy=="TEC"&Disease=="HRSMM"&Tissue=="PB")
df2<-sub%>%group_by(ID)%>%summarise(n=n_distinct(Sample_ID))

set.seed(1234)


# 1) Genera tutte le combinazioni
combos <- expand.grid(
  Timepoint     = unique(sub$Timepoint),
  sub_cluster = unique(sub$sub_cluster),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

n      <- 100
n.iter <- 100

num_cores <- parallel::detectCores() - 10
cl        <- makeCluster(num_cores)
registerDoParallel(cl)

pt_results <- foreach(i = seq_len(nrow(combos)), 
                      .packages = c("dplyr")) %dopar% {
                        
                        current_dis <- combos$Timepoint[i]
                        current_sub <- combos$sub_cluster[i]
                        
                        # Filtra i dati
                        tmp <- sub[sub$Timepoint == current_dis & sub$sub_cluster == current_sub, ]
                        
                        # Se meno di 100 celle, restituisci data frame vuoto
                        if(nrow(tmp) < 100) {
                          return(
                            list(
                              out.iter = data.frame(),
                              pt_clonp = data.frame()
                            )
                          )
                        }
                        
                        # Prepara out.iter
                        out.iter <- data.frame(
                          Timepoint  = rep(current_dis, n.iter),
                          Cluster  = rep(current_sub, n.iter),
                          Freq_rare   = NA_real_,
                          Freq_small  = NA_real_,
                          Freq_medium = NA_real_,
                          Freq_large  = NA_real_
                        )
                        
                        pt_clonp_final <- data.frame()
                        
                        # Loop sugli n.iter
                        for(x in seq_len(n.iter)) {
                          smpl <- tmp[sample(nrow(tmp), size = n, replace = FALSE), ]
                          
                          tab_smpl <- prop.table(table(smpl$clonotype_size))
                          
                          # Crea un data frame per "pt_clonp" di questo campione
                          pt_clonp <- data.frame(
                            type  = names(tab_smpl),
                            Freq  = as.numeric(tab_smpl),
                            # Se vuoi, puoi mettere Sample_ID o Disease etc.
                            Timepoint  = current_dis,
                            Cluster  = current_sub,
                            iter     = x
                          )
                          
                          # Nel caso "Rare" non compaia, fallback a 0 (se vuoi 1, vedi discussione in altra risposta)
                          out.iter$Freq_rare[x]   <- if("Rare"   %in% names(tab_smpl)) tab_smpl[["Rare"]]   else 0
                          out.iter$Freq_small[x]  <- if("Small"  %in% names(tab_smpl)) tab_smpl[["Small"]]  else 0
                          out.iter$Freq_medium[x] <- if("Medium" %in% names(tab_smpl)) tab_smpl[["Medium"]] else 0
                          out.iter$Freq_large[x]  <- if("Large"  %in% names(tab_smpl)) tab_smpl[["Large"]]  else 0
                          
                          # Se vuoi salvare i pt_clonp di ogni iter:
                          pt_clonp_final <- rbind(pt_clonp_final, pt_clonp)
                        }
                        
                        return(list(
                          out.iter = out.iter,
                          pt_clonp = pt_clonp_final
                        ))
                      }

stopCluster(cl)


all_out.iter <- bind_rows(lapply(pt_results, `[[`, "out.iter"))
all_pt_clonp <- bind_rows(lapply(pt_results, `[[`, "pt_clonp"))

# ------------------------------------------------------------------
# 5) Combine results from parallel output
# ------------------------------------------------------------------

# Aggregate by Sample_ID: calculate mean & SD for each frequency type
out_iter_summary <- all_out.iter %>%
  group_by(Timepoint, Cluster) %>%
  summarise(
    Mean_Rare   = mean(Freq_rare, na.rm = TRUE),
    SD_Rare     = sd(Freq_rare, na.rm = TRUE),
    Mean_Small  = mean(Freq_small, na.rm = TRUE),
    SD_Small    = sd(Freq_small, na.rm = TRUE),
    Mean_Medium = mean(Freq_medium, na.rm = TRUE),
    SD_Medium   = sd(Freq_medium, na.rm = TRUE),
    Mean_Large  = mean(Freq_large, na.rm = TRUE),
    SD_Large    = sd(Freq_large, na.rm = TRUE)
  ) %>%
  pivot_longer(
    cols = -c(Timepoint, Cluster), 
    names_to = c("Metric", "Clone_Type"), 
    names_sep = "_", 
    values_to = "Value"
  )

# Check the transformed structure
print(out_iter_summary)


# Define custom colors for Clone Types
clone_colors <- c(
  "Large" = "red", 
  "Medium" = "darkorange", 
  "Small" = "gold", 
  "Rare" = "lightgrey"
)

out_iter_summary<-out_iter_summary%>%pivot_wider(names_from = Metric, values_from = Value)
out_iter_summary$Clone_Type<-factor(out_iter_summary$Clone_Type, levels=c("Large", "Medium", "Small","Rare"))

order_info <- out_iter_summary %>%
  filter(Clone_Type == "Large") %>%
  group_by(Timepoint, Cluster) %>%
  # Ordina in base a Mean in senso decrescente
  arrange(Mean) %>%
  # Crea un nuovo fattore per Sample_ID con livelli nell'ordine desiderato
  mutate(Cluster_ordered = factor(Cluster, levels = unique(Cluster))) %>%
  ungroup() %>%
  # Manteniamo solo le colonne che ci servono
  dplyr::select(Timepoint, Cluster, Cluster_ordered)

# 2) Attacchiamo quell'ordinamento al dataset completo
out_iter_summary <- out_iter_summary %>%
  left_join(order_info, by = c("Timepoint", "Cluster"))

cd8_clusters<-unique(grep("CD8+",out_iter_summary$Cluster, value=TRUE))
cd8_clusters<-cd8_clusters[-c(1,3)]
out_iter_summary_filt<-out_iter_summary%>%filter(Cluster%in%cd8_clusters)
out_iter_summary_filt$Cluster_ordered<-gsub("_","", out_iter_summary_filt$Cluster_ordered)
out_iter_summary_filt$Cluster_ordered<-factor(out_iter_summary_filt$Cluster_ordered, levels=c("CD8+GZMK+TEM",
                                                                                              "CD8+DR+TEM","CD8+Tex",    "CD8+GZMB+TEM" , "CD8+KIR+TEM"  ))
out_iter_summary_filt$Timepoint<-factor(out_iter_summary_filt$Timepoint, levels=c("Pre","Post"))
baseline_clone_size_sample_distribution<-ggplot(out_iter_summary_filt, aes(x = Timepoint, y = Mean, fill = Clone_Type)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8,alpha=0.5, color="black") +
  labs(x = NULL, y = "Percentage of Clonotypes", fill = "Clone Size") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = clone_colors)+
  facet_wrap(~ Cluster_ordered, ncol=5) +
  geom_text(aes(label = paste0(round((Mean*100),0),"%")), 
            position = position_stack(vjust = 0.5), 
            size = 12, 
            color = "black") +
  theme_classic() +
  theme(title=element_text(size=40),
        plot.title = element_text(hjust=0.5, face = "bold"), text=element_text(size=50,color = "black"), 
        axis.text.x =element_text(angle=45, vjust = 1, hjust = 1, color="black"),legend.position = "bottom",
        axis.text.y =element_text(size=40,color = "black"), axis.title.y =element_text(size=50,color = "black")) +
  labs(y="Clones(%)") 
output_dir<-"~/Immune_project/full_dataset_analysis/plots/diversity/"
ggsave(paste0(output_dir,"Long_clone_size_Cluster_distribution_PB_HRSMM.pdf"),
       plot=baseline_clone_size_sample_distribution, width=30, height=25 )
write_csv(out_iter_summary_filt,paste0(output_dir, "long_frequency_HRSMM_TEC.csv"))



# -------------------------------------------------------------------------
#Same for RRMM
# -------------------------------------------------------------------------

sub <- T_cells_filt_clono@meta.data
#TAKE ONLY BASELINE PB
df<-sub%>%filter(Therapy%in%c("TEC")&Tissue=="PB"&Disease=="RRMM")%>%group_by(ID, Therapy, Disease)%>%summarise(n=n_distinct(Sample_ID))
paired_ids<-df%>%filter(n>1)
paired_ids<-unique(paired_ids$ID)

sub<-sub%>%filter(ID%in% paired_ids&Therapy=="TEC"&Disease=="RRMM"&Tissue=="PB")
df2<-sub%>%group_by(ID)%>%summarise(n=n_distinct(Sample_ID))

set.seed(1234)


# 1) Genera tutte le combinazioni
combos <- expand.grid(
  Timepoint     = unique(sub$Timepoint),
  sub_cluster = unique(sub$sub_cluster),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

n      <- 100
n.iter <- 100

num_cores <- parallel::detectCores() - 10
cl        <- makeCluster(num_cores)
registerDoParallel(cl)

pt_results <- foreach(i = seq_len(nrow(combos)), 
                      .packages = c("dplyr")) %dopar% {
                        
                        current_dis <- combos$Timepoint[i]
                        current_sub <- combos$sub_cluster[i]
                        
                        # Filtra i dati
                        tmp <- sub[sub$Timepoint == current_dis & sub$sub_cluster == current_sub, ]
                        
                        # Se meno di 100 celle, restituisci data frame vuoto
                        if(nrow(tmp) < 100) {
                          return(
                            list(
                              out.iter = data.frame(),
                              pt_clonp = data.frame()
                            )
                          )
                        }
                        
                        # Prepara out.iter
                        out.iter <- data.frame(
                          Timepoint  = rep(current_dis, n.iter),
                          Cluster  = rep(current_sub, n.iter),
                          Freq_rare   = NA_real_,
                          Freq_small  = NA_real_,
                          Freq_medium = NA_real_,
                          Freq_large  = NA_real_
                        )
                        
                        pt_clonp_final <- data.frame()
                        
                        # Loop sugli n.iter
                        for(x in seq_len(n.iter)) {
                          smpl <- tmp[sample(nrow(tmp), size = n, replace = FALSE), ]
                          
                          tab_smpl <- prop.table(table(smpl$clonotype_size))
                          
                          # Crea un data frame per "pt_clonp" di questo campione
                          pt_clonp <- data.frame(
                            type  = names(tab_smpl),
                            Freq  = as.numeric(tab_smpl),
                            # Se vuoi, puoi mettere Sample_ID o Disease etc.
                            Timepoint  = current_dis,
                            Cluster  = current_sub,
                            iter     = x
                          )
                          
                          # Nel caso "Rare" non compaia, fallback a 0 (se vuoi 1, vedi discussione in altra risposta)
                          out.iter$Freq_rare[x]   <- if("Rare"   %in% names(tab_smpl)) tab_smpl[["Rare"]]   else 0
                          out.iter$Freq_small[x]  <- if("Small"  %in% names(tab_smpl)) tab_smpl[["Small"]]  else 0
                          out.iter$Freq_medium[x] <- if("Medium" %in% names(tab_smpl)) tab_smpl[["Medium"]] else 0
                          out.iter$Freq_large[x]  <- if("Large"  %in% names(tab_smpl)) tab_smpl[["Large"]]  else 0
                          
                          # Se vuoi salvare i pt_clonp di ogni iter:
                          pt_clonp_final <- rbind(pt_clonp_final, pt_clonp)
                        }
                        
                        return(list(
                          out.iter = out.iter,
                          pt_clonp = pt_clonp_final
                        ))
                      }

stopCluster(cl)


all_out.iter <- bind_rows(lapply(pt_results, `[[`, "out.iter"))
all_pt_clonp <- bind_rows(lapply(pt_results, `[[`, "pt_clonp"))

# ------------------------------------------------------------------
# 5) Combine results from parallel output
# ------------------------------------------------------------------

# Aggregate by Sample_ID: calculate mean & SD for each frequency type
out_iter_summary <- all_out.iter %>%
  group_by(Timepoint, Cluster) %>%
  summarise(
    Mean_Rare   = mean(Freq_rare, na.rm = TRUE),
    SD_Rare     = sd(Freq_rare, na.rm = TRUE),
    Mean_Small  = mean(Freq_small, na.rm = TRUE),
    SD_Small    = sd(Freq_small, na.rm = TRUE),
    Mean_Medium = mean(Freq_medium, na.rm = TRUE),
    SD_Medium   = sd(Freq_medium, na.rm = TRUE),
    Mean_Large  = mean(Freq_large, na.rm = TRUE),
    SD_Large    = sd(Freq_large, na.rm = TRUE)
  ) %>%
  pivot_longer(
    cols = -c(Timepoint, Cluster), 
    names_to = c("Metric", "Clone_Type"), 
    names_sep = "_", 
    values_to = "Value"
  )

# Check the transformed structure
print(out_iter_summary)


# Define custom colors for Clone Types
clone_colors <- c(
  "Large" = "red", 
  "Medium" = "darkorange", 
  "Small" = "gold", 
  "Rare" = "lightgrey"
)

out_iter_summary<-out_iter_summary%>%pivot_wider(names_from = Metric, values_from = Value)
out_iter_summary$Clone_Type<-factor(out_iter_summary$Clone_Type, levels=c("Large", "Medium", "Small","Rare"))

order_info <- out_iter_summary %>%
  filter(Clone_Type == "Large") %>%
  group_by(Timepoint, Cluster) %>%
  # Ordina in base a Mean in senso decrescente
  arrange(Mean) %>%
  # Crea un nuovo fattore per Sample_ID con livelli nell'ordine desiderato
  mutate(Cluster_ordered = factor(Cluster, levels = unique(Cluster))) %>%
  ungroup() %>%
  # Manteniamo solo le colonne che ci servono
  dplyr::select(Timepoint, Cluster, Cluster_ordered)

# 2) Attacchiamo quell'ordinamento al dataset completo
out_iter_summary <- out_iter_summary %>%
  left_join(order_info, by = c("Timepoint", "Cluster"))

cd8_clusters<-unique(grep("CD8+",out_iter_summary$Cluster, value=TRUE))
cd8_clusters<-cd8_clusters[-c(1,3)]
out_iter_summary_filt<-out_iter_summary%>%filter(Cluster%in%cd8_clusters)
out_iter_summary_filt$Cluster_ordered<-gsub("_","", out_iter_summary_filt$Cluster_ordered)
out_iter_summary_filt$Cluster_ordered<-factor(out_iter_summary_filt$Cluster_ordered, levels=c("CD8+GZMK+TEM",
                                                                                              "CD8+DR+TEM","CD8+Tex",    "CD8+GZMB+TEM" , "CD8+KIR+TEM"  ))
out_iter_summary_filt$Timepoint<-factor(out_iter_summary_filt$Timepoint, levels=c("Pre","Post"))
baseline_clone_size_sample_distribution<-ggplot(out_iter_summary_filt, aes(x = Timepoint, y = Mean, fill = Clone_Type)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8,alpha=0.5, color="black") +
  labs(x = NULL, y = "Percentage of Clonotypes", fill = "Clone Size") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = clone_colors)+
  facet_wrap(~ Cluster_ordered, ncol=5) +
  geom_text(aes(label = paste0(round((Mean*100),0),"%")), 
            position = position_stack(vjust = 0.5), 
            size = 12, 
            color = "black") +
  theme_classic() +
  theme(title=element_text(size=40),
        plot.title = element_text(hjust=0.5, face = "bold"), text=element_text(size=50,color = "black"), 
        axis.text.x =element_text(angle=45, vjust = 1, hjust = 1, color="black"),legend.position = "bottom",
        axis.text.y =element_text(size=40,color = "black"), axis.title.y =element_text(size=50,color = "black")) +
  labs(y="Clones(%)") 
output_dir<-"~/Immune_project/full_dataset_analysis/plots/diversity/"
ggsave(paste0(output_dir,"Long_clone_size_Cluster_distribution_PB_RRMM.pdf"),
       plot=baseline_clone_size_sample_distribution, width=30, height=25 )
write_csv(out_iter_summary_filt,paste0(output_dir, "long_frequency_RRMM_TEC.csv"))


# -------------------------------------------------------------------------
# Barplots showing the proportion of T cells in a given BL patient sample (P, n=14) or sample from a HD (n=11) that were determined 
# to belong to one of four clone size categories (Rare: ≤1%; Small: >1% and <5%; Medium: ≥5% and <10%; Large: ≥10%)
# through iterative (n=100) downsampling of 100 cells. The average proportion per clone size category was visualized and
# the standard deviation across iterations was depicted in solid-line error bars.
# -------------------------------------------------------------------------
# Load data
out_iter_combined<-read_csv("/mnt/disks/disk/full_dataset_qs/diversity/freq_clone_sizes_samples_100iter.csv")
# Aggregate by Sample_ID: calculate mean & SD for each frequency type
out_iter_summary <- out_iter_combined %>%
  group_by(Sample_ID) %>%
  summarise(
    Mean_Rare   = mean(Freq_rare, na.rm = TRUE),
    SD_Rare     = sd(Freq_rare, na.rm = TRUE),
    Mean_Small  = mean(Freq_small, na.rm = TRUE),
    SD_Small    = sd(Freq_small, na.rm = TRUE),
    Mean_Medium = mean(Freq_medium, na.rm = TRUE),
    SD_Medium   = sd(Freq_medium, na.rm = TRUE),
    Mean_Large  = mean(Freq_large, na.rm = TRUE),
    SD_Large    = sd(Freq_large, na.rm = TRUE)
  ) %>%
  pivot_longer(
    cols = -Sample_ID, 
    names_to = c("Metric", "Clone_Type"), 
    names_sep = "_", 
    values_to = "Value"
  )

# Check the transformed structure
print(out_iter_summary)


meta_filtered <- meta %>% dplyr::select(-Library)
meta_filtered<-meta_filtered%>%distinct(Sample_ID, .keep_all = TRUE)
out_iter_summary <- left_join(out_iter_summary, meta_filtered, by = "Sample_ID")

#TAKE ONLY BASELINE PB
df<-out_iter_summary%>%filter(Therapy%in%c("TEC")&Tissue=="PB"&Disease=="HRSMM")%>%group_by(ID, Therapy, Disease)%>%summarise(n=n_distinct(Timepoint))
paired_ids<-df%>%filter(n>1)
paired_ids<-unique(paired_ids$ID)

out_iter_summary_SMM_long<-out_iter_summary%>%filter(ID%in%paired_ids &Therapy%in%c("TEC")&Tissue=="PB"&Disease=="HRSMM")
# Define custom colors for Clone Types
clone_colors <- c(
  "Large" = "red", 
  "Medium" = "darkorange", 
  "Small" = "gold", 
  "Rare" = "lightgrey"
)

out_iter_summary_SMM_long<-out_iter_summary_SMM_long%>%pivot_wider(names_from = Metric, values_from = Value)
out_iter_summary_SMM_long$Clone_Type<-factor(out_iter_summary_SMM_long$Clone_Type, levels=c("Large", "Medium", "Small","Rare"))


out_iter_summary_SMM_long$Timepoint<-factor(out_iter_summary_SMM_long$Timepoint, levels=c("Pre","Post"))
baseline_clone_size_sample_distribution<-ggplot(out_iter_summary_SMM_long, aes(x = ID, y = Mean, fill = Clone_Type)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7,alpha=0.5,color="black") +
  labs(x = NULL, y = "Clones(%)", fill = "Clone Size") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = clone_colors)+
  facet_wrap(~ Timepoint, ncol=1) +
  theme_minimal() +
  theme(title=element_text(size=40),
        plot.title = element_text(hjust=0.5, face = "bold"), text=element_text(size=50,color = "black"), 
        legend.position = "bottom",
        axis.ticks.x =element_blank(),axis.text.x  =element_blank(),
        axis.text.y =element_text(size=40,color = "black"), axis.title.y =element_text(size=50,color = "black"))
  
ggsave(paste0(output_dir,"Long_clone_size_sample_distribution_PB_HRSMM.pdf"),
       plot=baseline_clone_size_sample_distribution, width=13, height=12 )


# -------------------------------------------------------------------------
#RRMM
# -------------------------------------------------------------------------
df<-out_iter_summary%>%filter(Therapy%in%c("TEC")&Tissue=="PB"&Disease=="RRMM")%>%group_by(ID, Therapy, Disease)%>%summarise(n=n_distinct(Timepoint))
paired_ids<-df%>%filter(n>1)
paired_ids<-unique(paired_ids$ID)

out_iter_summary_RRMM_long<-out_iter_summary%>%filter(ID%in%paired_ids &Therapy%in%c("TEC")&Tissue=="PB"&Disease=="RRMM")
# Define custom colors for Clone Types
clone_colors <- c(
  "Large" = "red", 
  "Medium" = "darkorange", 
  "Small" = "gold", 
  "Rare" = "lightgrey"
)

out_iter_summary_RRMM_long<-out_iter_summary_RRMM_long%>%pivot_wider(names_from = Metric, values_from = Value)
out_iter_summary_RRMM_long$Clone_Type<-factor(out_iter_summary_RRMM_long$Clone_Type, levels=c("Large", "Medium", "Small","Rare"))


out_iter_summary_RRMM_long$Timepoint<-factor(out_iter_summary_RRMM_long$Timepoint, levels=c("Pre","Post"))
baseline_clone_size_sample_distribution<-ggplot(out_iter_summary_RRMM_long, aes(x = ID, y = Mean, fill = Clone_Type)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7,alpha=0.5,color="black") +
  labs(x = NULL, y = "Clones(%)", fill = "Clone Size") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = clone_colors)+
  facet_wrap(~ Timepoint, ncol=1) +
  theme_minimal() +
  theme(title=element_text(size=40),
        plot.title = element_text(hjust=0.5, face = "bold"), text=element_text(size=50,color = "black"), 
        legend.position = "bottom",
        axis.ticks.x =element_blank(),axis.text.x  =element_blank(),
        axis.text.y =element_text(size=40,color = "black"), axis.title.y =element_text(size=50,color = "black"))

ggsave(paste0(output_dir,"Long_clone_size_sample_distribution_PB_RRMM.pdf"),
       plot=baseline_clone_size_sample_distribution, width=13, height=12 )



# -------------------------------------------------------------------------
#take all the clonotypes of expanded cells in HRSMM post samples for each sub_cluster and each sample
#then look in paired Pre-samples all T-cells to see if we find an intersection (using ID for pairing) and for the intersected one
#lets'see how many cells in that clone are representated by which cell types
# -------------------------------------------------------------------------
#First filter for HRSMM and Post and for type2 different from rare
sub<-subset(T_cells_filt_clono, Therapy=="TEC"&Disease=="RRMM"&Tissue=="PB"&Timepoint=="Post"&clonotype_size%in%c("Small","Medium","Large"))

df_post_fraction <- sub@meta.data %>%
  group_by(ID, clonotypeID) %>%
  summarise(
    n_total = n(),
    n_cd8TEM_Tex = sum(sub_cluster %in% c("CD8+_GZMB+_TEM","CD8+_GZMK+_TEM", "CD8+_KIR+_TEM","CD8+_DR+_TEM", "CD8+_Tex")),
    frac_cd8TEM_Tex = n_cd8TEM_Tex / n_total,
    .groups = "drop"
  ) %>%
  # Keep only clones with >50% in CD8+_TEM or Tex
  filter(frac_cd8TEM_Tex > 0.5)


# Raggruppa per ID e raccogli tutti i clonotipi presenti in ogni sub_cluster e campione
df <- df_post_fraction %>%
  group_by(ID) %>%
  summarise(clonotypes_post = list(unique(clonotypeID)), .groups = "drop")

#Take all pre (no restriction for clonotype size)
sub<-subset(T_cells_filt_clono, Therapy=="TEC"&Disease=="RRMM"&Tissue=="PB"&Timepoint=="Pre")
# Raggruppa per ID e raccogli tutti i clonotipi presenti in ogni sub_cluster e campione
df2 <- sub@meta.data %>%
  group_by(ID) %>%
  summarise(clonotypes_pre = list(unique(clonotypeID)), .groups = "drop")

df3<-left_join(df, df2, by="ID")


# Find the intersection of clonotypes for each ID
df3 <- df3 %>%
  rowwise() %>%
  mutate(
    clonotype_intersection = list(intersect(unlist(clonotypes_post), unlist(clonotypes_pre))),
    num_intersected = length(clonotype_intersection)
  ) %>%
  ungroup()

# Expand the dataframe so that each intersected clonotype gets its own row
df_expanded <- df3 %>%
  dplyr::select(ID, clonotype_intersection) %>%
  unnest(cols = clonotype_intersection)

# View results
df_expanded<-df_expanded%>%mutate(ID_clonotype=paste0(ID,clonotype_intersection))
keep<-unique(df_expanded$ID_clonotype)
T_cells_filt_clono@meta.data<-T_cells_filt_clono@meta.data%>%mutate(ID_clonotype=paste0(ID,clonotypeID))
sub2<-subset(T_cells_filt_clono, ID_clonotype%in%keep)
sub2@meta.data <- sub2@meta.data %>%
  mutate(
    sub_cluster2 = ifelse(
      !sub_cluster %in% c("CD8+_GZMB+_TEM", "CD8+_GZMK+_TEM", "CD8+_KIR+_TEM",
                          "CD8+_DR+_TEM",   "CD8+_Tex",       "Tgd", 
                          "cycling_Tcells"),
      "Other",
      as.character(sub_cluster)  # force to character
    )
  )
sub2@meta.data$sub_cluster2<-gsub("_", " ", sub2@meta.data$sub_cluster2)

library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(tidyr)

# 1. Ensure Timepoint is a factor with "Pre" before "Post"
df_for_plot <- sub2@meta.data %>%
  mutate(Timepoint = factor(Timepoint, levels = c("Pre","Post"))) %>%
  group_by(Timepoint, sub_cluster2) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(Timepoint) %>%
  mutate(prop = count / sum(count))

# 2. Donut Plot: Pre first, then Post
donut<-ggplot(df_for_plot, aes(x = 2, y = prop, fill = sub_cluster2)) +
  geom_bar(stat = "identity", color = "white", width = 1) +
  coord_polar(theta = "y", start = 0) +
  
  # Make 2 side-by-side donuts
  facet_wrap(~ Timepoint, ncol = 2) +
  
  # Create the "hole" by limiting x-range
  xlim(0.5, 2.5) +
  
  # Use a more saturated color palette
  scale_fill_brewer(palette = "Set1") +  
  # Try "Dark2" if you'd like an alternate look:
  # scale_fill_brewer(palette = "Dark2") +
  
  # Remove background/ticks, etc.
  theme_void() +
  theme(
    title = element_text(hjust=0.5),
    text=element_text(size=50),
    legend.position = "bottom",
  ) +
  
  # Labels
  labs(
    fill = "Sub-cluster",
    title = "Phenotype tracking of post-treatment expanded cytotoxic clones",
    subtitle = ""
  )
ggsave(plot=donut, filename=paste0(output_dir, "donut.pdf"), width=15, height=12)


meta_sub<-sub2@meta.data
ids<-unique(sub2$ID)
library(ggplot2)
library(dplyr)

for(i in seq_along(ids)){
  tmp <- meta_sub %>% filter(ID == ids[[i]])
  
  # Count cells per clonotype per timepoint
  counts <- tmp %>%
    group_by(clonotypeID, Timepoint) %>%
    summarise(cell_count = n(), .groups = "drop") %>%
    pivot_wider(names_from = Timepoint, values_from = cell_count, values_fill = 0)
  
  # Merge counts back into tmp for annotation
  tmp <- left_join(tmp, counts, by = "clonotypeID")
  
  # Determine number of columns dynamically
  n_col <- length(unique(tmp$clonotypeID))
  
  tmp$Timepoint<-factor(tmp$Timepoint,levels=c("Pre","Post"))
  p <- ggplot(tmp, aes(x = Timepoint, fill = sub_cluster)) +
    geom_bar(position = "fill") +
    facet_wrap(~clonotypeID, ncol = n_col) +
    theme_classic() +
    theme(legend.position = "bottom") +
    labs(y = "Proportion", x = "Timepoint") +
    geom_text(data = counts, aes(x = 1.5, y = 1.05, label = paste0("Pre: ", Pre, "\nPost: ", Post)), inherit.aes = FALSE, size = 3)
  
  ggsave(p, filename = paste0(output_dir, "clonotype_tracking_", ids[[i]], ".pdf"), width=12, height=10)
}




# -------------------------------------------------------------------------
#diversity x Response
# -------------------------------------------------------------------------

output_dir <- "/mnt/disks/disk/full_dataset_qs/diversity/"

#load meta
meta <- read_excel("~/Immune_project/full_dataset_analysis/Meta_Immune_proj.xlsx")
meta <- meta%>%
  mutate(response_disease = case_when(
    BOR %in% c("MR","SD","PR", "PD") & Disease == "RRMM" ~ "RRMM_NR",            # RRMM Non-Responder
    Disease == "RRMM" & !BOR %in% c("MR","SD","PR", "PD") ~ "RRMM_R",            # RRMM Responder
    Disease == "HRSMM" & !Therapy== "Len" ~ "HRSMM_R",                     # SMM Responder
    Disease == "Healthy" ~ "Healthy_NA",                              # Healthy: NA response
    TRUE ~ "SMM_NA"                                                    # Default: SMM Non-Applicable
  ))
meta<-meta%>%distinct(Sample_ID,.keep_all = TRUE)

out_mean<-read_csv(paste0(output_dir, "Diversity_Means_jan2025.csv"))
#ADD METADATA<- -2 sampples with less than 100 cells with clonotypes
out_mean<-left_join(out_mean, meta[,c("Sample_ID","response_disease","BOR","Timepoint","Tissue", "Therapy","Disease","ID")], by="Sample_ID")

df<-out_mean%>%filter(Therapy%in%c("TEC","TALQ"))%>%group_by(ID,Tissue, Therapy, Disease)%>%summarise(n=n_distinct(Timepoint))
paired_ids<-df%>%filter(n>1)
paired_ids<-paired_ids%>%mutate(Tissue_ID=paste0(Tissue, "_",ID))
paired_ids<-unique(paired_ids$Tissue_ID)
out_mean<-out_mean%>%mutate(Tissue_ID=paste0(Tissue, "_",ID))
paired_counts<-out_mean%>%filter(Tissue_ID%in%paired_ids)
df2<-paired_counts%>%filter(Therapy%in%c("TEC","TALQ"))%>%group_by(ID, Tissue,Disease)%>%summarise(n=n_distinct(Sample_ID))
paired_counts$Timepoint<-factor(paired_counts$Timepoint, levels=c("Pre","Post"))
chao_long<-ggplot(paired_counts, aes(x=Timepoint, y=1/Chao, fill=Timepoint))+
  geom_line(aes(group=ID),size=0.8) +
  geom_jitter(size=10, width=0.1,aes(fill=BOR,color=Therapy, stroke=2),shape=21, alpha=0.5)+
  theme_classic()+guides(color=FALSE)+
  facet_wrap(~response_disease, scales="free")+
  theme(plot.margin = unit(c(0.5, 0, 0.2, 0.1), "cm"),title=element_text(size=40),
        plot.title = element_text(hjust=0.5, face = "bold"), text=element_text(size=50,color = "black"), axis.text.x =element_text(size=50,color = "black"),
        axis.text.y =element_text(size=40,color = "black"), axis.title.y =element_text(size=50,color = "black"), legend.position = "bottom")+
  labs(title="", y="1/Chao Index", x="")+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))+stat_compare_means(method = "wilcox.test", paired=TRUE, size=10, label="p.format")

output_dir<-"~/Immune_project/full_dataset_analysis/plots/diversity/"
ggsave(paste0(output_dir,"response.pdf"),
       plot=chao_long, width=13, height=12 )

