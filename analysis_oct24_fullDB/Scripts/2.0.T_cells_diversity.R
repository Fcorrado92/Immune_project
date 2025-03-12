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
library(doParallel)
library(foreach)
output_dir <- "/mnt/disks/disk/full_dataset_qs/diversity/"

#load meta
meta <- read_excel("~/Immune_project/full_dataset_analysis/Meta_Immune_proj_feb25.xlsx")
meta <- meta%>%
  mutate(response_disease = case_when(
  BOR %in% c("PR", "PD") & Disease == "RRMM" ~ "RRMM_NR",            # RRMM Non-Responder
  Disease == "RRMM" & !BOR %in% c("PR", "PD") ~ "RRMM_R",            # RRMM Responder
  Disease == "HRSMM" & !Therapy== "Len" ~ "HRSMM_R",                     # SMM Responder
  Disease == "Healthy" ~ "Healthy_NA",                              # Healthy: NA response
  TRUE ~ "SMM_NA"                                                    # Default: SMM Non-Applicable
))
meta<-meta%>%distinct(Sample_ID,.keep_all = TRUE)
#load_object
T_cells_filt<-qread("/mnt/disks/disk/full_dataset_qs/T_cells_jan25_final.qs")

#add new Accelerator VDJ_TCR
##Do this before filtering souporcell singlets
vdj_folder<-"~/Immune_project/TCR/"
samples_VDJ<-  list.files(vdj_folder)
samples_VDJ <- grep("41|42|43|44|45|46|47|48", samples_VDJ, value=TRUE)
samples_VDJ<-samples_VDJ[c(1:8)]
setwd("~/Immune_project/TCR/")
#addTCRdata
getTCR <- function(samples) {
  #This function will return a dataframe with unique cell barcodes, all TRA/TRB chain clonotype information and MAIT/iNKT evidence in a single row.
  ##Load filtered_contig_annotations.csv files for all samples
  df_l <- lapply(samples, function(x) {
    read.csv(paste0(x, "/", "filtered_contig_annotations.csv"))
  })
  ## Load clonotypes.csv files for all samples
  cl_l <- lapply(samples, function(x) {
    read.csv(paste0(x, "/", "clonotypes.csv"))
  })
  ##Parse filtered_contig_annotations.csv files for each sample
  for (i in seq_along(df_l)){
    df <- df_l[[i]]
    df <- df[(df$chain %in% c("TRA","TRB")) & (df$productive %in% c("true","TRUE","True")),]
    ##For every chain per cell barcode, keep the row with the most UMIs (remove multiple rows per chain)
    ###Sorting by UMIs rather than reads.
    df <- df %>% group_by(barcode,chain) %>% arrange(desc(umis)) %>% slice_head(n=1)
    ##For every cell barcode, keep the two chains with the most UMIs (remove multiple chains per barcode)
    ###Using slice_head instead of top_n to handle ties in the number of UMIs per chain.
    df <- df %>% group_by(barcode) %>% arrange(desc(umis)) %>% slice_head(n=2)
    ##Add columns with TRA/TRB chain gene-level clonotype
    df <- df %>%
      mutate(TRAct = ifelse(chain == "TRA", paste(interaction(v_gene,  j_gene, c_gene)), NA)) %>%
      mutate(TRBct = ifelse(chain == "TRB", paste(interaction(v_gene,  j_gene, c_gene)), NA))
    ##Create dataframe with unique cell barcodes in rows and TRA/TRB chain gene-, amino acid- and nucleotide-level clonotype information
    out <- parseTCR(df)
    ##Add MAIT/iNKT evidence per cell
    cl <- cl_l[[i]]
    out$MAIT <- cl$mait_evidence[match(as.character(out$clonotypeID), as.character(cl$clonotype_id))]
    out$iNKT <- cl$inkt_evidence[match(as.character(out$clonotypeID), as.character(cl$clonotype_id))]
    out$clonotype_freq <- cl$proportion[match(as.character(out$clonotypeID), as.character(cl$clonotype_id))]
    sample_name <- samples[i]
    out$Library <- sub("Accelerator_(B[0-9]+).*", "\\1", sample_name)
    df_l[[i]] <- out
  }
  names(df_l) <- samples
  return(df_l)
}
parseTCR <- function(df){
  #This function will output a dataframe of unique cell barcodes and their TRA & TRB chain
  #clonotype information in a single row
  ##Create vector of unique cell barcodes in dataframe
  unique_cb <- unique(df$barcode)
  ##Create output dataframe where each cell barcode is matched to heavy and light chain clonotype information in a single row
  out <- data.frame(matrix(nrow=length(unique_cb),ncol=10))
  colnames(out) <- c("barcode","TRA_gene_ct","TRA_cdr3_aa","TRA_cdr3_nt","TRA_V_gene","TRB_gene_ct","TRB_cdr3_aa","TRB_cdr3_nt","TRB_V_gene","clonotypeID")
  out$barcode <- unique_cb
  ##For every unique cell barcode:
  for (i in seq_along(unique_cb)){
    barcode.i <- out$barcode[i]
    location.i <- which(df$barcode == barcode.i)
    ##Write in its clonotype ID
    ##Write in its clonotype ID
    if(length(location.i)>1){
      out[i,"clonotypeID"] <- df[location.i[1],"raw_clonotype_id"]
    } else {
      out[i,"clonotypeID"] <- df[location.i,"raw_clonotype_id"]
    }
    ##If there are two rows corresponding to it in the input data
    if (length(location.i) == 2){
      ##and the first row is not a TRA chain row
      if (is.na(df[location.i[1],c("TRAct")])) {
        ##then, write in the TRA chain clonotype info from the second row
        out[i,c("TRA_gene_ct","TRA_cdr3_aa","TRA_cdr3_nt","TRA_V_gene")]<-df[location.i[2], c("TRAct","cdr3","cdr3_nt","v_gene")]
        ##and write in the TRB chain clonotype info from the first row
        out[i,c("TRB_gene_ct","TRB_cdr3_aa","TRB_cdr3_nt","TRB_V_gene")]<-df[location.i[1], c("TRBct","cdr3","cdr3_nt","v_gene")]
      } else {
        ##otherwise, write in the TRA chain clonotype info from the first row
        out[i,c("TRA_gene_ct","TRA_cdr3_aa","TRA_cdr3_nt","TRA_V_gene")]<-df[location.i[1], c("TRAct","cdr3","cdr3_nt","v_gene")]
        ##and write in the TRB chain clonotype info from the second row
        out[i,c("TRB_gene_ct","TRB_cdr3_aa","TRB_cdr3_nt","TRB_V_gene")]<-df[location.i[2], c("TRBct","cdr3","cdr3_nt","v_gene")]
      }
      ##otherwise, if there is only one row for this cell barcode
    } else if (length(location.i) == 1) {
      chain.i <- df$chain[location.i]
      ##and it corresponds to a TRA chain
      if (chain.i == "TRA"){
        ##write in the TRA chain clonotype information from that row
        out[i,c("TRA_gene_ct","TRA_cdr3_aa","TRA_cdr3_nt","TRA_V_gene")]<-df[location.i, c("TRAct","cdr3","cdr3_nt","v_gene")]
        ##otherwise
      } else {
        ##write in the TRB chain clonotype information from that row
        out[i,c("TRB_gene_ct","TRB_cdr3_aa","TRB_cdr3_nt","TRB_V_gene")]<-df[location.i, c("TRBct","cdr3","cdr3_nt","v_gene")]
      }
    }
  }
  return(out)
}
vdj<-getTCR(samples_VDJ)
vdj_df <- bind_rows(vdj, .id = "SampleID")
vdj_df<-vdj_df%>%dplyr::select(c("Library","barcode","clonotypeID","TRB_cdr3_aa", "TRA_cdr3_aa"))

T_cells_filt@meta.data$Library<-T_cells_filt@meta.data$Pool
T_cells_filt$new_barcodes<-rownames(T_cells_filt@meta.data)
T_cells_filt@meta.data <- T_cells_filt@meta.data %>%
  left_join(vdj_df, by = c("barcode", "Library"), suffix = c("", ".new")) %>%
  mutate(
    # For each column we want to update, coalesce the old and the new
    clonotypeID = coalesce(clonotypeID, clonotypeID.new),
    TRB_cdr3_aa = coalesce(TRB_cdr3_aa, TRB_cdr3_aa.new),
    TRA_cdr3_aa = coalesce(TRA_cdr3_aa, TRA_cdr3_aa.new)
    
      ) %>%
  # Finally, remove the .new columns
  dplyr::select(-ends_with(".new"))

# Then assign back if this is a Seurat object
rownames(T_cells_filt@meta.data)<-T_cells_filt$new_barcodes
qsave(T_cells_filt,"/mnt/disks/disk/full_dataset_qs/T_cells_jan25_final_newTCR.qs")
#filter
T_cells_filt<-qread("/mnt/disks/disk/full_dataset_qs/T_cells_jan25_final_newTCR.qs")
clusters_to_remove <- grep("^dbl", unique(Idents(T_cells_filt)),value=T)
clusters_to_remove<-unique(c(clusters_to_remove, "Plt", "lowquality", "Lowq", "lowq","Plt"))
Idents(T_cells_filt)<-T_cells_filt$sub_cluster
keep<-setdiff(unique(Idents(T_cells_filt)), clusters_to_remove)
T_cells_filt<- subset(T_cells_filt, idents = keep)
DimPlot(T_cells_filt, label=T)
#remove rows with clonotypeID ==NA
remove <- unique(rownames(T_cells_filt@meta.data[is.na(T_cells_filt@meta.data$clonotypeID), ]))
keep<-setdiff(unique(rownames(T_cells_filt@meta.data)), remove)
T_cells_filt_clono<-subset(T_cells_filt, cells = keep)
nrow(T_cells_filt_clono@meta.data)

# # #visualize percentage of T_cells with associated clonotype
metadata<-T_cells_filt@meta.data
percentage_results <- list()
sample_ids<-unique(metadata$Sample_ID)
for (sample_id in sample_ids) {
  # Filtra i metadata per il campione corrente
  sample_metadata <- metadata %>%
    filter(Sample_ID == sample_id)

  # Calcola il valore percentuale per il campione corrente
  tcr_prop <- sum(!is.na(sample_metadata$clonotypeID)) / nrow(sample_metadata) * 100
  tcr_count <- sum(!is.na(sample_metadata$clonotypeID))
  n<-nrow(sample_metadata)
  # Conserva il risultato nella lista

  percentage_results[[length(percentage_results) + 1]] <- data.frame(
    Sample_ID = sample_id,
    MRN=unique(sample_metadata$ID),
    Pool_ID = unique(sample_metadata$Pool),
    Total_T_cells=n,
    tcr_prop = tcr_prop,
    tcr_count=tcr_count
  )
}
percentage_df <- do.call(rbind, percentage_results)
print(percentage_df)
#percentage of associated clonotypes
ggplot(percentage_df,aes(x=tcr_prop))+geom_histogram()


# -------------------------------------------------------------------------
# ------------------------------------------------------------------
# 1) Load libraries
# ------------------------------------------------------------------

meta<-read_excel("~/Immune_project/full_dataset_analysis/Meta_Immune_proj_feb25.xlsx")
meta<-distinct(meta, Sample_ID, .keep_all = TRUE)

# ------------------------------------------------------------------
# 4) Parallelize the Diversity Calculations
# ------------------------------------------------------------------
sub<-T_cells_filt_clono@meta.data
set.seed(1234)
n      <- 100       # sample size per iteration
n.iter <- 5000      # number of iterations (10K in your comment, but set to 1000 as in your code)
pt     <- unique(sub$Sample_ID)
sub$pt_clonotype <- sub$clonotypeID

# --- A) Set up parallel backend ---
num_cores <- parallel::detectCores() - 10
cl        <- makeCluster(num_cores)
registerDoParallel(cl)

# We will collect two pieces of information in parallel:
# 1) out.iter: the Richness, Shannon, Chao estimates across iterations
# 2) pt_clonp: the proportion of each clone (once per patient)

# We can return a list of length 2 for each patient:
#   - First data frame: out.iter
#   - Second data frame: pt_clonp
# Then we'll combine them at the end.

pt_results <- foreach(i = seq_along(pt), 
                      .packages = c("dplyr")) %dopar% {
                        
                        current_pt <- pt[i]
                        tmp <- sub[sub$Sample_ID %in% current_pt, ]
                        
                        # If less than 100 T cells with associated clonotype, skip
                        # Return empty data frames
                        if (nrow(tmp) < 100) {
                          return(list(out.iter = data.frame(),
                                      pt_clonp = data.frame()))
                        }
                        
                        # Prepare output dataframe for iterations
                        out.iter <- data.frame(
                          "Sample_ID" = rep(current_pt, n.iter), 
                          "Richness"  = NA_real_,
                          "Shannon"   = NA_real_,
                          "Chao"      = NA_real_
                        )
                        
                        # Summarize proportion of each clonotype once per patient
                        pt_clonp <- data.frame(prop.table(table(tmp$pt_clonotype)))
                        colnames(pt_clonp) <- c("pt_clonotype", "Freq")
                        pt_clonp$type <- ifelse(pt_clonp$Freq <= 0.01, "Rare",
                                                ifelse(pt_clonp$Freq > 0.01 & pt_clonp$Freq < 0.05, "Small",
                                                       ifelse(pt_clonp$Freq >= 0.05 & pt_clonp$Freq < 0.1, "Medium",
                                                              ifelse(pt_clonp$Freq >= 0.1, "Large", NA))))
                        pt_clonp$Sample_ID <- current_pt
                        
                        # --- B) Now the inner loop (x in 1:n.iter) ---
                        for (x in seq_len(n.iter)) {
                          # if(x %% 100 == 0) {
                          #   print(paste0("Patient: ", current_pt, " - Iteration ", x))
                          # }
                          
                          # Randomly sample n=100 T cells
                          smpl <- tmp[sample(x = nrow(tmp), size = n, replace = FALSE), ]
                          
                          # Calculate Shannon
                          clonp     <- data.frame(prop.table(table(smpl$pt_clonotype)))
                          clonp$logp <- log(clonp$Freq)
                          shannon   <- (-sum(clonp$Freq * clonp$logp, na.rm = TRUE)) / log(nrow(clonp))
                          
                          # Calculate Richness + Chao
                          clon  <- data.frame(table(smpl$pt_clonotype))
                          chao  <- nrow(clon) + (
                            (sum(clon$Freq == 1, na.rm = TRUE) * (sum(clon$Freq == 1, na.rm = TRUE) - 1)) /
                              (2 * (sum(clon$Freq == 2,  na.rm = TRUE) + 1))
                          )
                          
                          out.iter$Richness[x] <- nrow(clon)
                          out.iter$Shannon[x]  <- shannon
                          out.iter$Chao[x]     <- chao
                        }
                        
                        # Return the two data frames as a list
                        return(list(out.iter = out.iter, pt_clonp = pt_clonp))
                      }

# Shut down cluster
stopCluster(cl)

# ------------------------------------------------------------------
# 5) Combine results from parallel output
# ------------------------------------------------------------------

# pt_results is a list of length = length(pt). 
# Each element is a list: list(out.iter = ..., pt_clonp = ...).

# Extract out.iter from each patient
all_out_iters <- lapply(pt_results, `[[`, "out.iter")
# Extract pt_clonp from each patient
all_pt_clonp  <- lapply(pt_results, `[[`, "pt_clonp")

# Combine
out     <- do.call(rbind, all_out_iters)
sub.out <- do.call(rbind, all_pt_clonp)

# ------------------------------------------------------------------
# 6) Write outputs (same as your code)
# ------------------------------------------------------------------
write.csv(out,     paste0(output_dir, "Diversity_jan2025.csv"),  row.names = FALSE)
write.csv(sub.out, paste0(output_dir, "CloneSize_jan2025.csv"),  row.names = FALSE)

out_mean <- aggregate(out[, c("Richness", "Shannon", "Chao")],
                      list(out$Sample_ID),
                      mean, na.rm = TRUE)
colnames(out_mean)[1] <- "Sample_ID"


write.csv(out_mean, paste0(output_dir, "Diversity_Means_jan2025.csv"), row.names = FALSE)


# -------------------------------------------------------------------------
#AVERAGE DIVERSITY AND FREQUENCY AT THE SAMPLE LEVEL
out_mean<-read_csv(paste0(output_dir, "Diversity_Means_jan2025.csv"))

#ADD METADATA<- -2 sampples with less than 100 cells with clonotypes
out_mean<-left_join(out_mean, meta[,c("Sample_ID","response_disease","BOR","Timepoint","Tissue", "Therapy","Disease","ID")], by="Sample_ID")


#PLOT BASELINE DIFFERENCE IN DIVERSITY BETWEEN HD, RRMM, SMM -------------------------------------------------------------------------
#plot Chao&Richness at baseline
out_mean_pre<-out_mean%>%filter(Timepoint%in%c("Pre", "Healthy")&Therapy%in%c("TEC","Healthy","TALQ","CART")&Tissue=="PB")
out_mean_pre%>%group_by(Timepoint,Tissue, Therapy,Disease)%>%summarise(n=n_distinct(Sample_ID))
# Timepoint Tissue Therapy Disease     n
# 1 Healthy   PB     Healthy Healthy     9
# 2 Pre       PB     CART    HRSMM       6
# 3 Pre       PB     CART    RRMM       27
# 4 Pre       PB     TALQ    RRMM        3
# 5 Pre       PB     TEC     HRSMM      17
# 6 Pre       PB     TEC     RRMM        8
cols=c("darkgreen","darkblue", "darkred")
out_mean_pre$Disease[out_mean_pre$Disease=="SMM"]<-"HRSMM"
out_mean_pre$Disease<-factor(out_mean_pre$Disease, levels=c("Healthy","HRSMM","RRMM"))
disease_totals <- out_mean_pre %>%
  group_by(Disease) %>%
  summarise(Total = n(), .groups = "drop")
chao_baseline<-ggplot(out_mean_pre, aes(x=Disease, y=Chao, fill=Disease))+
  geom_violin(alpha=0.5)+geom_boxplot(size=3,alpha=0.6, outlier.shape = NA)+
  geom_jitter(size=15,width=0.3,aes(fill=Disease, stroke=2),color="black",shape=21, alpha=0.5)+
  labs(y="Chao Index Mean")+
  theme_classic()+guides(fill=FALSE,color=FALSE)+
  scale_fill_manual(values=cols)+
  theme(plot.margin = unit(c(0.5, 0, 0.2, 0.1), "cm"),title=element_text(size=40),
        plot.title = element_text(hjust=0.5, face = "bold"), text=element_text(size=50,color = "black"), axis.text.x =element_text(size=50,color = "black"),
        axis.text.y =element_text(size=40,color = "black"), axis.title.y =element_text(size=50,color = "black"))+
  labs(title="TCR Repertoire Diversity", y="Averaged Chao Index", x="")+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))+
geom_text(
  data = disease_totals, 
  aes(x = Disease, y = -300, label = paste0("n=", Total)), # Adjust y to -10%
  inherit.aes = TRUE, 
  size = 20, 
  hjust = 0.5
)
chao_baseline<-chao_baseline+ geom_pwc(label.size = 12, method = "wilcox.test",step.increase = 0.08, tip.length = 0.01, size=1.2)
chao_baseline<-ggadjust_pvalue(
  chao_baseline, p.adjust.method = "BH",
  label = paste("q=","{p.adj.format}"),
  hide.ns = "p"
)
output_dir<-"~/Immune_project/full_dataset_analysis/plots/diversity/"
ggsave(paste0(output_dir,"chao_baseline.pdf"),
       plot=chao_baseline, width=13, height=12 )


# BM ----------------------------------------------------------------------
#plot Chao&Richness at baseline
out_mean_pre<-out_mean%>%filter(Timepoint%in%c("Pre", "Healthy")&Therapy%in%c("TEC","Healthy","TALQ","CART")&Tissue=="BM")
out_mean_pre%>%group_by(Timepoint,Tissue, Therapy,Disease)%>%summarise(n=n_distinct(Sample_ID))
# 1 Healthy   BM     Healthy Healthy    11
# 2 Pre       BM     TEC     HRSMM       4
# 3 Pre       BM     TEC     RRMM        9
cols=c("darkgreen","darkblue", "darkred")
out_mean_pre$Disease[out_mean_pre$Disease=="SMM"]<-"HRSMM"
out_mean_pre$Disease<-factor(out_mean_pre$Disease, levels=c("Healthy","HRSMM","RRMM"))

# Calculate total counts for each Therapy group
disease_totals <- out_mean_pre %>%
  group_by(Disease) %>%
  summarise(Total = n(), .groups = "drop")

chao_baseline<-ggplot(out_mean_pre, aes(x=Disease, y=Chao, fill=Disease))+
  geom_violin(alpha=0.5)+geom_boxplot(size=3,alpha=0.6, outlier.shape = NA)+
  geom_jitter(size=15,width=0.3,aes(fill=Disease, stroke=2),color="black",shape=21, alpha=0.5)+
  labs(y="Chao Index Mean")+
  theme_classic()+guides(fill=FALSE,color=FALSE)+
  scale_fill_manual(values=cols)+
  theme(plot.margin = unit(c(0.5, 0, 0.2, 0.1), "cm"),title=element_text(size=40),
        plot.title = element_text(hjust=0.5, face = "bold"), text=element_text(size=50,color = "black"), axis.text.x =element_text(size=50,color = "black"),
        axis.text.y =element_text(size=40,color = "black"), axis.title.y =element_text(size=50,color = "black"))+
  labs(title="TCR Repertoire Diversity", y="Averaged Chao Index", x="")+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))+
  geom_text(
    data = disease_totals, 
    aes(x = Disease, y = -0.1, label = paste0("n=", Total)), # Adjust y to -10%
    inherit.aes = TRUE, 
    size = 20, 
    hjust = 0.5
  )

chao_baseline<-chao_baseline+ geom_pwc(label.size = 12, method = "wilcox.test",step.increase = 0.08, tip.length = 0.01, size=1.2)
chao_baseline<-ggadjust_pvalue(
  chao_baseline, p.adjust.method = "BH",
  label = paste("q=","{p.adj.format}"),
  hide.ns = "p"
)
output_dir<-"~/Immune_project/full_dataset_analysis/plots/diversity/"
ggsave(paste0(output_dir,"BM_chao_baseline.pdf"),
       plot=chao_baseline, width=13, height=12 )



# -------------------------------------------------------------------------
#assign stable clone size
# -------------------------------------------------------------------------
output_dir <- "/mnt/disks/disk/full_dataset_qs/diversity/"
sub.out<-read_csv(paste0(output_dir, "CloneSize_jan2025.csv"))

get_mode <- function(x) {
  uniqx <- unique(x)
  # tabulate(match(...)) finds how often each unique value appears
  uniqx[which.max(tabulate(match(x, uniqx)))]
}

# 2. Use the function in your dplyr pipeline
sub.out<-sub.out %>%
  group_by(pt_clonotype, Sample_ID) %>%
  summarise(mode_freq = get_mode(Freq))

#assign to T-cells_filt_clono
sub.out$clonotypeID<-sub.out$pt_clonotype
T_cells_filt_clono@meta.data<-left_join(T_cells_filt_clono@meta.data, sub.out, by=c("Sample_ID","clonotypeID"))
#now relabel based on the mode_freq
T_cells_filt_clono@meta.data <- T_cells_filt_clono@meta.data %>%
  mutate(clonotype_size = case_when(
    mode_freq <= 0.01                 ~ "Rare",
    mode_freq > 0.01  & mode_freq < 0.05  ~ "Small",
    mode_freq >= 0.05 & mode_freq < 0.1   ~ "Medium",
    mode_freq >= 0.1                    ~ "Large",
    TRUE                                ~ NA_character_
  ))         

rownames(T_cells_filt_clono@meta.data)<-T_cells_filt_clono@meta.data$new_barcodes
#remove samples with less than 100 cells
df<-as.data.frame(table(T_cells_filt_clono@meta.data$Sample_ID))
df<-df%>%filter(Freq>100)
keep<-unique(df$Var1)
T_cells_filt_clono<-subset(T_cells_filt_clono, Sample_ID%in%keep)
#save T_cells_filt_clono qithout samples with less than 100 cells
T_cells_filt_clono@meta.data<-T_cells_filt_clono@meta.data[-c(44:51)]
T_cells_filt_clono@meta.data<-left_join(T_cells_filt_clono@meta.data, meta, by="Sample_ID")
rownames(T_cells_filt_clono@meta.data)<-T_cells_filt_clono@meta.data$new_barcodes
View(T_cells_filt_clono@meta.data)
qsave(T_cells_filt_clono, paste0(output_dir, "T_cells_clono.qs"))
# -------------------------------------------------------------------------
#plot UMAP before treatment, Healthy, SMM, RRMM, TEC, CART, TALQ PB
# -------------------------------------------------------------------------
T_cells_filt_clono_pre<-subset(T_cells_filt_clono, Timepoint%in%c("Pre","Healthy")&Therapy%in%c("CART","Healthy","TEC","TALQ")&Tissue=="PB")
T_cells_filt_clono_pre@meta.data$clonotype_size<-factor(T_cells_filt_clono_pre@meta.data$clonotype_size, levels=c("Large", "Medium","Small","Rare"))
cols<-rev(c("lightgrey", "gold", "darkorange", "red"))

# Or in a tidy way:
group_counts <- T_cells_filt_clono_pre@meta.data %>%
  group_by(Disease) %>%
  tally(name = "num_cells")
group_counts

# 3. Determine the smallest group size (or pick your target size)
min_cells <- min(group_counts$num_cells)
min_cells

# 4. Downsample each disease group to the same size
#    Note: rownames(T_cells_filt_clono_pre@meta.data) == cell/barcode names
downsampled_barcodes <- T_cells_filt_clono_pre@meta.data %>%
  rownames_to_column(var = "cell_id") %>%  # Convert rownames to a column for cell IDs
  group_by(Disease) %>%
  sample_n(min_cells) %>%
  pull(cell_id)

# 5. Create a new Seurat object with only the sampled cells
T_cells_filt_clono_pre_ds <- subset(
  T_cells_filt_clono_pre,
  cells = downsampled_barcodes
)

T_cells_filt_clono_pre_ds@meta.data$Disease<-factor(T_cells_filt_clono_pre_ds@meta.data$Disease, levels=c("Healthy","HRSMM","RRMM"))
Clone_Size_UMAP <- DimPlot(
  T_cells_filt_clono_pre_ds, 
  group.by = "clonotype_size", 
  split.by = "Disease"
) +
  scale_color_manual(
    name = "Clone Size",   # Oppure = NULL per nascondere il titolo della legenda
    values = cols
  ) +
  labs(
    x = "UMAP 1", 
    y = "UMAP 2",
    title = NULL              # oppure title="" se vuoi stringa vuota
  ) +
  theme(
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    # legenda in basso, centrata
    legend.position = "bottom",
    legend.justification = "center",
    legend.box = "horizontal",  # assicura che la legenda sia orizzontale
    plot.title = element_blank() # rimuove il titolo generato di default da Seurat, se presente
  )

output_dir<-"~/Immune_project/full_dataset_analysis/plots/diversity/"
ggsave(paste0(output_dir, "UMAP_size_clones_pre_PB.pdf"), plot=Clone_Size_UMAP, width=7, height=4)


# -------------------------------------------------------------------------
# Barplots showing the proportion of T cells in a given BL patient sample (P, n=14) or sample from a HD (n=11) that were determined 
# to belong to one of four clone size categories (Rare: ≤1%; Small: >1% and <5%; Medium: ≥5% and <10%; Large: ≥10%)
# through iterative (n=100) downsampling of 100 cells. The average proportion per clone size category was visualized and
# the standard deviation across iterations was depicted in solid-line error bars.
# -------------------------------------------------------------------------
# Load data
T_cells_filt_clono <- qread("/mnt/disks/disk/full_dataset_qs/diversity/T_cells_clono.qs .qs")

sub <- T_cells_filt_clono@meta.data
set.seed(1234)

n      <- 100       # Sample size per iteration
n.iter <- 100       # Number of iterations
pt     <- unique(sub$Sample_ID)
sub$pt_clonotype <- sub$clonotypeID

# --- A) Set up parallel backend ---
num_cores <- parallel::detectCores() - 10
cl        <- makeCluster(num_cores)
registerDoParallel(cl)

# Run the foreach loop
pt_results <- foreach(i = seq_along(pt), 
                      .packages = c("dplyr")) %dopar% {
                        
                        current_pt <- pt[i]
                        tmp <- sub[sub$Sample_ID == current_pt, ]  # <-- Fix here (use dynamic Sample_ID)
                        
                        # If less than 100 T cells with associated clonotype, skip
                        if (nrow(tmp) < 100) {
                          return(list(out.iter = data.frame(),
                                      pt_clonp = data.frame()))
                        }
                        
                        # Prepare output dataframe for iterations
                        out.iter <- data.frame(
                          "Sample_ID"  = rep(current_pt, n.iter), 
                          "Freq_rare"  = NA_real_,
                          "Freq_small" = NA_real_,
                          "Freq_medium"= NA_real_,
                          "Freq_large" = NA_real_
                        )
                        
                        # --- B) Now the inner loop (x in 1:n.iter) ---
                        for (x in seq_len(n.iter)) {
                          # Randomly sample n=100 T cells
                          smpl <- tmp[sample(nrow(tmp), size = n, replace = FALSE), ]
                          
                          # Summarize proportion of each clonotype size
                          pt_clonp <- data.frame(prop.table(table(smpl$clonotype_size)))
                          colnames(pt_clonp) <- c("type", "Freq")
                          pt_clonp$Sample_ID <- current_pt  # Corrected variable name
                          
                          # Fill out.iter with calculated proportions
                          out.iter$Freq_rare[x]   <- ifelse("Rare"   %in% pt_clonp$type, pt_clonp$Freq[pt_clonp$type == "Rare"], 1)
                          out.iter$Freq_small[x]  <- ifelse("Small"  %in% pt_clonp$type, pt_clonp$Freq[pt_clonp$type == "Small"], 0)
                          out.iter$Freq_medium[x] <- ifelse("Medium" %in% pt_clonp$type, pt_clonp$Freq[pt_clonp$type == "Medium"], 0)
                          out.iter$Freq_large[x]  <- ifelse("Large"  %in% pt_clonp$type, pt_clonp$Freq[pt_clonp$type == "Large"], 0)
                        }
                        
                        # Return the two data frames as a list
                        return(list(out.iter = out.iter, pt_clonp = pt_clonp))
                      }

# Shut down cluster
stopCluster(cl)


# ------------------------------------------------------------------
# 5) Combine results from parallel output
# ------------------------------------------------------------------

out_iter_combined <- bind_rows(lapply(pt_results, `[[`, "out.iter"))
write_csv(out_iter_combined, "/mnt/disks/disk/full_dataset_qs/diversity/freq_clone_sizes_samples_100iter.csv")
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

out_iter_summary_PRE_PB<-out_iter_summary%>%filter(Timepoint%in%c("Pre","Healthy")&Therapy%in%c("CART","Healthy","TEC","TALQ")&Tissue=="PB")
# Define custom colors for Clone Types
clone_colors <- c(
  "Large" = "red", 
  "Medium" = "darkorange", 
  "Small" = "gold", 
  "Rare" = "lightgrey"
)

out_iter_summary_PRE_PB<-out_iter_summary_PRE_PB%>%pivot_wider(names_from = Metric, values_from = Value)
out_iter_summary_PRE_PB$Clone_Type<-factor(out_iter_summary_PRE_PB$Clone_Type, levels=c("Large", "Medium", "Small","Rare"))

order_info <- out_iter_summary_PRE_PB %>%
  filter(Clone_Type == "Large") %>%
  group_by(Disease) %>%
  # Ordina in base a Mean in senso decrescente
  arrange(Mean) %>%
  # Crea un nuovo fattore per Sample_ID con livelli nell'ordine desiderato
  mutate(Sample_ID_ordered = factor(Sample_ID, levels = unique(Sample_ID))) %>%
  ungroup() %>%
  # Manteniamo solo le colonne che ci servono
  dplyr::select(Disease, Sample_ID, Sample_ID_ordered)

# 2) Attacchiamo quell'ordinamento al dataset completo
out_iter_summary_PRE_PB <- out_iter_summary_PRE_PB %>%
  left_join(order_info, by = c("Disease", "Sample_ID"))


baseline_clone_size_sample_distribution<-ggplot(out_iter_summary_PRE_PB, aes(x = Sample_ID_ordered, y = Mean, fill = Clone_Type)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7,alpha=0.7) +
  labs(x = NULL, y = "Percentage of Clonotypes", fill = "Clone Size") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = clone_colors)+
  facet_wrap(~ Disease, scales = "free_x") +
  theme_minimal() +
  theme(title=element_text(size=40),
        plot.title = element_text(hjust=0.5, face = "bold"), text=element_text(size=50,color = "black"), 
        axis.text.x =element_blank(),legend.position = "bottom",
        axis.ticks.x =element_blank(),
        axis.text.y =element_text(size=40,color = "black"), axis.title.y =element_text(size=50,color = "black")) +
  labs(y="Clones(%)") 
ggsave(paste0(output_dir,"baseline_clone_size_sample_distribution_PB.pdf"),
       plot=baseline_clone_size_sample_distribution, width=13, height=12 )


# Create a lookup table of Rare clone type Mean values per sample
rare_means <- out_iter_summary_PRE_PB %>%
  filter(Clone_Type == "Rare") %>%
  dplyr::select(Sample_ID,Disease, Therapy, Timepoint, rare_Mean = Mean)%>%
  mutate(expanded_Mean=1-rare_Mean)
# Calculate total counts for each Therapy group
disease_totals <- rare_means %>%
  group_by(Disease) %>%
  summarise(Total = n(), .groups = "drop")

baseline_clone_size_disease_distribution<-ggplot(rare_means, aes(x = Disease, y = expanded_Mean, fill = Disease)) +
  geom_violin(alpha=0.5)+geom_boxplot(size=3,alpha=0.6, outlier.shape = NA)+
  geom_jitter(size=15,width=0.3,aes(fill=Disease, stroke=2),color="black",shape=21, alpha=0.5)+
  labs(y="Expanded Clones(%)")+
  scale_fill_manual(values=cols)+
  theme_classic()+guides(fill=FALSE,color=FALSE)+
    theme(plot.margin = unit(c(0.5, 0, 0.2, 0.1), "cm"),title=element_text(size=40),
        plot.title = element_text(hjust=0.5, face = "bold"), text=element_text(size=50,color = "black"),
        axis.text.x =element_text(size=50,color = "black"),
        legend.position="bottom",
        axis.text.y =element_text(size=40,color = "black"), axis.title.y =element_text(size=50,color = "black"))+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))+
  geom_text(
    data = disease_totals, 
    aes(x = Disease, y = -0.1, label = paste0("n=", Total)), # Adjust y to -10%
    inherit.aes = TRUE, 
    size = 20, 
    hjust = 0.5
  )

baseline_clone_size_disease_distribution<-baseline_clone_size_disease_distribution+ 
  geom_pwc(label.size = 12, method = "wilcox.test",step.increase = 0.08, tip.length = 0.01, size=1.2)
baseline_clone_size_disease_distribution<-ggadjust_pvalue(
  baseline_clone_size_disease_distribution, p.adjust.method = "BH",
  label = paste("q=","{p.adj.format}"),
  hide.ns = "p"
)


ggsave(paste0(output_dir,"boxplots_expanded_disease_distribution.pdf"),
       plot=baseline_clone_size_disease_distribution, width=13, height=12 )

# -------------------------------------------------------------------------
########BM

out_iter_summary_PRE_BM<-out_iter_summary%>%filter(Timepoint%in%c("Pre","Healthy")&Therapy%in%c("CART","Healthy","TEC","TALQ")&Tissue=="BM")
# Define custom colors for Clone Types
clone_colors <- c(
  "Large" = "red", 
  "Medium" = "darkorange", 
  "Small" = "gold", 
  "Rare" = "lightgrey"
)

out_iter_summary_PRE_BM<-out_iter_summary_PRE_BM%>%pivot_wider(names_from = Metric, values_from = Value)
out_iter_summary_PRE_BM$Clone_Type<-factor(out_iter_summary_PRE_BM$Clone_Type, levels=c("Large", "Medium", "Small","Rare"))

order_info <- out_iter_summary_PRE_BM %>%
  filter(Clone_Type == "Large") %>%
  group_by(Disease) %>%
  # Ordina in base a Mean in senso decrescente
  arrange(Mean) %>%
  # Crea un nuovo fattore per Sample_ID con livelli nell'ordine desiderato
  mutate(Sample_ID_ordered = factor(Sample_ID, levels = unique(Sample_ID))) %>%
  ungroup() %>%
  # Manteniamo solo le colonne che ci servono
  dplyr::select(Disease, Sample_ID, Sample_ID_ordered)

# 2) Attacchiamo quell'ordinamento al dataset completo
out_iter_summary_PRE_BM <- out_iter_summary_PRE_BM %>%
  left_join(order_info, by = c("Disease", "Sample_ID"))


baseline_clone_size_sample_distribution<-ggplot(out_iter_summary_PRE_BM, aes(x = Sample_ID_ordered, y = Mean, fill = Clone_Type)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7,alpha=0.7) +
  labs(x = NULL, y = "Percentage of Clonotypes", fill = "Clone Size") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = clone_colors)+
  facet_wrap(~ Disease, scales = "free_x") +
  theme_minimal() +
  theme(title=element_text(size=40),
        plot.title = element_text(hjust=0.5, face = "bold"), text=element_text(size=50,color = "black"), 
        axis.text.x =element_blank(),legend.position = "bottom",
        axis.ticks.x =element_blank(),
        axis.text.y =element_text(size=40,color = "black"), axis.title.y =element_text(size=50,color = "black")) +
  labs(y="Clones(%)") 
ggsave(paste0(output_dir,"baseline_clone_size_sample_distribution_BM.pdf"),
       plot=baseline_clone_size_sample_distribution, width=13, height=12 )


# -------------------------------------------------------------------------


# Create a lookup table of Rare clone type Mean values per sample
rare_means <- out_iter_summary_PRE_BM %>%
  filter(Clone_Type == "Rare") %>%
  dplyr::select(Sample_ID,Disease, Therapy, Timepoint, rare_Mean = Mean)%>%
  mutate(expanded_Mean=1-rare_Mean)

# Calculate total counts for each Therapy group
disease_totals <- rare_means %>%
  group_by(Disease) %>%
  summarise(Total = n(), .groups = "drop")

baseline_clone_size_disease_distribution<-ggplot(rare_means, aes(x = Disease, y = expanded_Mean, fill = Disease)) +
  geom_violin(alpha=0.5)+geom_boxplot(size=3,alpha=0.6, outlier.shape = NA)+
  geom_jitter(size=15,width=0.3,aes(fill=Disease, stroke=2),color="black",shape=21, alpha=0.5)+
  labs(y="Expanded Clones(%)")+
  scale_fill_manual(values=cols)+
  theme_classic()+guides(fill=FALSE,color=FALSE)+
  theme(plot.margin = unit(c(0.5, 0, 0.2, 0.1), "cm"),title=element_text(size=40),
        plot.title = element_text(hjust=0.5, face = "bold"), text=element_text(size=50,color = "black"),
        axis.text.x =element_text(size=50,color = "black"),
        legend.position="bottom",
        axis.text.y =element_text(size=40,color = "black"), axis.title.y =element_text(size=50,color = "black"))+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))+
  geom_text(
    data = disease_totals, 
    aes(x = Disease, y = -0.1, label = paste0("n=", Total)), # Adjust y to -10%
    inherit.aes = TRUE, 
    size = 20, 
    hjust = 0.5
  )
baseline_clone_size_disease_distribution<-baseline_clone_size_disease_distribution+ 
  geom_pwc(label.size = 12, method = "wilcox.test",step.increase = 0.08, tip.length = 0.01, size=1.2)
baseline_clone_size_disease_distribution<-ggadjust_pvalue(
  baseline_clone_size_disease_distribution, p.adjust.method = "BH",
  label = paste("q=","{p.adj.format}")
)


ggsave(paste0(output_dir,"BM_boxplots_expanded_disease_distribution.pdf"),
       plot=baseline_clone_size_disease_distribution, width=13, height=12 )


# -------------------------------------------------------------------------
# E) Barplots showing the proportion of clonotypes in a given T cell subtype across all patients (n=14) or HD (n=11) 
# that belonged to one of the four clone size categories. For each T cell subtype, 100 cells were randomly sampled 100 times
# from all patients or HD,
# and the proportion of expanded (1-Rare) clonotypes was compared between patients and HD using bootstrapping with 10,000 iterations.
# -------------------------------------------------------------------------
T_cells_filt_clono <- qread("/mnt/disks/disk/full_dataset_qs/diversity/T_cells_clono.qs")

sub <- T_cells_filt_clono@meta.data
#TAKE ONLY BASELINE PB
sub<-sub%>%filter(Timepoint%in%c("Pre","Healthy")&Therapy%in%c("CART","Healthy","TEC","TALQ")&Tissue=="PB")

set.seed(1234)


# 1) Genera tutte le combinazioni
combos <- expand.grid(
  Disease     = unique(sub$Disease),
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
                        
                        current_dis <- combos$Disease[i]
                        current_sub <- combos$sub_cluster[i]
                        
                        # Filtra i dati
                        tmp <- sub[sub$Disease == current_dis & sub$sub_cluster == current_sub, ]
                        
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
                          Disease  = rep(current_dis, n.iter),
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
                            Disease  = current_dis,
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
  group_by(Disease, Cluster) %>%
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
    cols = -c(Disease, Cluster), 
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
  group_by(Disease, Cluster) %>%
  # Ordina in base a Mean in senso decrescente
  arrange(Mean) %>%
  # Crea un nuovo fattore per Sample_ID con livelli nell'ordine desiderato
  mutate(Cluster_ordered = factor(Cluster, levels = unique(Cluster))) %>%
  ungroup() %>%
  # Manteniamo solo le colonne che ci servono
  dplyr::select(Disease, Cluster, Cluster_ordered)

# 2) Attacchiamo quell'ordinamento al dataset completo
out_iter_summary <- out_iter_summary %>%
  left_join(order_info, by = c("Disease", "Cluster"))

cd8_clusters<-unique(grep("CD8+",out_iter_summary$Cluster, value=TRUE))
cd8_clusters<-cd8_clusters[-c(1,3)]
out_iter_summary_filt<-out_iter_summary%>%filter(Cluster%in%cd8_clusters)
out_iter_summary_filt$Cluster_ordered<-gsub("_","", out_iter_summary_filt$Cluster_ordered)
out_iter_summary_filt$Cluster_ordered<-factor(out_iter_summary_filt$Cluster_ordered, levels=c("CD8+GZMK+TEM",
                                                                                              "CD8+DR+TEM","CD8+Tex",    "CD8+GZMB+TEM" , "CD8+KIR+TEM"  ))
baseline_clone_size_sample_distribution<-ggplot(out_iter_summary_filt, aes(x = Disease, y = Mean, fill = Clone_Type)) +
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
ggsave(paste0(output_dir,"baseline_clone_size_Cluster_distribution_PB_BL.pdf"),
       plot=baseline_clone_size_sample_distribution, width=30, height=25 )
write_csv(out_iter_summary_filt,paste0(output_dir, "baseline_frequency_exp_cluster_PB_DS.csv"))

# -------------------------------------------------------------------------
sub <- T_cells_filt_clono@meta.data
#TAKE ONLY BASELINE PB
sub<-sub%>%filter(Timepoint%in%c("Pre","Healthy")&Therapy%in%c("CART","Healthy","TEC","TALQ")&Tissue=="BM")

set.seed(1234)


# 1) Genera tutte le combinazioni
combos <- expand.grid(
  Disease     = unique(sub$Disease),
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
                        
                        current_dis <- combos$Disease[i]
                        current_sub <- combos$sub_cluster[i]
                        
                        # Filtra i dati
                        tmp <- sub[sub$Disease == current_dis & sub$sub_cluster == current_sub, ]
                        
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
                          Disease  = rep(current_dis, n.iter),
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
                            Disease  = current_dis,
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
  group_by(Disease, Cluster) %>%
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
    cols = -c(Disease, Cluster), 
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
  group_by(Disease, Cluster) %>%
  # Ordina in base a Mean in senso decrescente
  arrange(Mean) %>%
  # Crea un nuovo fattore per Sample_ID con livelli nell'ordine desiderato
  mutate(Cluster_ordered = factor(Cluster, levels = unique(Cluster))) %>%
  ungroup() %>%
  # Manteniamo solo le colonne che ci servono
  dplyr::select(Disease, Cluster, Cluster_ordered)

# 2) Attacchiamo quell'ordinamento al dataset completo
out_iter_summary <- out_iter_summary %>%
  left_join(order_info, by = c("Disease", "Cluster"))

cd8_clusters<-unique(grep("CD8+",out_iter_summary$Cluster, value=TRUE))
cd8_clusters<-cd8_clusters[-c(1,3)]
out_iter_summary_filt<-out_iter_summary%>%filter(Cluster%in%cd8_clusters)
out_iter_summary_filt$Cluster_ordered<-gsub("_", "", out_iter_summary_filt$Cluster_ordered)
out_iter_summary_filt$Cluster_ordered<-factor(out_iter_summary_filt$Cluster_ordered, levels=c("CD8+GZMK+TEM",
                                                                                              "CD8+DR+TEM","CD8+Tex",    "CD8+GZMB+TEM" ,   "CD8+KIR+TEM"  ))
baseline_clone_size_sample_distribution<-ggplot(out_iter_summary_filt, aes(x = Disease, y = Mean, fill = Clone_Type)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8,alpha=0.5, color="black") +
  labs(x = NULL, y = "Percentage of Clonotypes", fill = "Clone Size") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = clone_colors)+
  facet_wrap(~ Cluster_ordered, ncol=5) +
  geom_text(
    aes(label = ifelse(Mean > 0, paste0(round(Mean * 100, 0), "%"), NA)), 
    position = position_stack(vjust = 0.5), 
    size = 12, 
    color = "black", 
    na.rm = TRUE  # Ignora le etichette con NA
  ) +  theme_classic() +
  theme(title=element_text(size=40),
        plot.title = element_text(hjust=0.5, face = "bold"), text=element_text(size=50,color = "black"), 
        axis.text.x =element_text(angle=45, vjust = 1, hjust = 1, color="black"),legend.position = "bottom",
        axis.text.y =element_text(size=40,color = "black"), axis.title.y =element_text(size=50,color = "black")) +
  labs(y="Clones(%)") 
ggsave(paste0(output_dir,"baseline_clone_size_Cluster_distribution_BM_BL.pdf"),
       plot=baseline_clone_size_sample_distribution, width=30, height=25 )
write_csv(out_iter_summary_filt,paste0(output_dir, "baseline_frequency_exp_cluster_BM_DS.csv"))


# -------------------------------------------------------------------------
########BM

out_iter_summary_PRE_BM<-out_iter_summary%>%filter(Timepoint%in%c("Pre","Healthy")&Therapy%in%c("CART","Healthy","TEC","TALQ")&Tissue=="BM")
# Define custom colors for Clone Types
clone_colors <- c(
  "Large" = "red", 
  "Medium" = "darkorange", 
  "Small" = "gold", 
  "Rare" = "lightgrey"
)

out_iter_summary_PRE_BM<-out_iter_summary_PRE_BM%>%pivot_wider(names_from = Metric, values_from = Value)
out_iter_summary_PRE_BM$Clone_Type<-factor(out_iter_summary_PRE_BM$Clone_Type, levels=c("Large", "Medium", "Small","Rare"))

order_info <- out_iter_summary_PRE_BM %>%
  filter(Clone_Type == "Large") %>%
  group_by(Disease) %>%
  # Ordina in base a Mean in senso decrescente
  arrange(Mean) %>%
  # Crea un nuovo fattore per Sample_ID con livelli nell'ordine desiderato
  mutate(Sample_ID_ordered = factor(Sample_ID, levels = unique(Sample_ID))) %>%
  ungroup() %>%
  # Manteniamo solo le colonne che ci servono
  dplyr::select(Disease, Sample_ID, Sample_ID_ordered)

# 2) Attacchiamo quell'ordinamento al dataset completo
out_iter_summary_PRE_BM <- out_iter_summary_PRE_BM %>%
  left_join(order_info, by = c("Disease", "Sample_ID"))


baseline_clone_size_sample_distribution<-ggplot(out_iter_summary_PRE_BM, aes(x = Sample_ID_ordered, y = Mean, fill = Clone_Type)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7,alpha=0.7) +
  labs(x = NULL, y = "Percentage of Clonotypes", fill = "Clone Size") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = clone_colors)+
  facet_wrap(~ Disease, scales = "free_x") +
  theme_minimal() +
  theme(title=element_text(size=40),
        plot.title = element_text(hjust=0.5, face = "bold"), text=element_text(size=50,color = "black"), 
        axis.text.x =element_blank(),legend.position = "bottom",
        axis.ticks.x =element_blank(),
        axis.text.y =element_text(size=40,color = "black"), axis.title.y =element_text(size=50,color = "black")) +
  labs(y="Clones(%)") 
ggsave(paste0(output_dir,"baseline_clone_size_sample_distribution_BM.pdf"),
       plot=baseline_clone_size_sample_distribution, width=13, height=12 )




# -------------------------------------------------------------------------
#plot UMAP before treatment, Healthy, SMM, RRMM, TEC, CART, TALQ BM
# -------------------------------------------------------------------------
T_cells_filt_clono_pre<-subset(T_cells_filt_clono, Timepoint%in%c("Pre","Healthy")&Therapy%in%c("CART","Healthy","TEC","TALQ")&Tissue=="BM")
table(T_cells_filt_clono_pre$Disease)
T_cells_filt_clono_pre@meta.data$clonotype_size<-factor(T_cells_filt_clono_pre@meta.data$clonotype_size, levels=c("Large", "Medium","Small","Rare"))
cols<-rev(c("lightgrey", "gold", "darkorange", "red"))
#downsample
# Or in a tidy way:
group_counts <- T_cells_filt_clono_pre@meta.data %>%
  group_by(Disease) %>%
  tally(name = "num_cells")
group_counts

# 3. Determine the smallest group size (or pick your target size)
min_cells <- min(group_counts$num_cells)
min_cells

# 4. Downsample each disease group to the same size
#    Note: rownames(T_cells_filt_clono_pre@meta.data) == cell/barcode names
downsampled_barcodes <- T_cells_filt_clono_pre@meta.data %>%
  rownames_to_column(var = "cell_id") %>%  # Convert rownames to a column for cell IDs
  group_by(Disease) %>%
  sample_n(min_cells) %>%
  pull(cell_id)

# 5. Create a new Seurat object with only the sampled cells
T_cells_filt_clono_pre_ds <- subset(
  T_cells_filt_clono_pre,
  cells = downsampled_barcodes
)

Clone_Size_UMAP <- DimPlot(
  T_cells_filt_clono_pre_ds, 
  group.by = "clonotype_size", 
  split.by = "Disease"
) +
  scale_color_manual(
    name = "Clone Size",   # Oppure = NULL per nascondere il titolo della legenda
    values = cols
  ) +
  labs(
    x = "UMAP 1", 
    y = "UMAP 2",
    title = NULL              # oppure title="" se vuoi stringa vuota
  ) +
  theme(
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    # legenda in basso, centrata
    legend.position = "bottom",
    legend.justification = "center",
    legend.box = "horizontal",  # assicura che la legenda sia orizzontale
    plot.title = element_blank() # rimuove il titolo generato di default da Seurat, se presente
  )

output_dir<-"~/Immune_project/full_dataset_analysis/plots/diversity/"
ggsave(paste0(output_dir, "UMAP_size_clones_pre_BM.pdf"), plot=Clone_Size_UMAP, width=7, height=4)


# -------------------------------------------------------------------------
#UMAP T_cells
# -------------------------------------------------------------------------
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

T_cells_filt_clono_pre<-RenameIdents(T_cells_filt_clono_pre, 'CD8+_GZMK_TEM'="CD8+_GZMK+_TEM", 'CD8+_GZMB_TEM' ="CD8+_GZMB+_TEM",
                         'CD8+_ZEB2_TEM'="CD8+_ZEB2+_TEM"
)
T_cells_filt_clono_pre[['ident']]<-Idents(T_cells_filt_clono_pre)
#UMAP sub_cluster Fig1+ number of patients + number of samples+ number of cells
T_cells_filt_clono_pre$sub_cluster<-Idents(T_cells_filt_clono_pre)
T_cells_filt_clono_pre$sub_cluster<-gsub("_", " ", T_cells_filt_clono_pre$sub_cluster)
T_cells_filt_clono_pre[['ident']]<-T_cells_filt_clono_pre$sub_cluster
#UMAP sub_cluster Fig1+ number of patients + number of samples+ number of cells
T_cells_filt_clono_pre[['ident']]<-T_cells_filt_clono_pre$sub_cluster

UMAP <- DimPlot(T_cells_filt_clono_pre) + 
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    text = element_text(size = 14),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
  ) + 
  labs(x = "UMAP 1", y = "UMAP 2") +  
  scale_color_manual(values = colors) + 
  NoLegend()

# Add labels with reduced overlap
UMAP <- LabelClusters(plot = UMAP, id = "ident", repel = TRUE, size = 4, nudge_y = -0.2, nudge_x=+0.2,
                      max.overlaps=getOption("ggrepel.max.overlaps", default = 15), force=10)

ggsave("~/Immune_project/full_dataset_analysis/plots/UMAP/UMAP_T_cells.pdf",
       plot=UMAP, width=17, height=14, units="cm")


# -------------------------------------------------------------------------


#Theobserved difference will be represented by 
# To check for significant differences in expanded frequency
# Bootstrap SMM+RRMM cell clusters strtified by Disease(10000times)
# within bootstrapped sample calculate average frequency of expanded clones after downsampling(n=100 *100)
# calulate 95%ci for difference in average freq
T_cells_filt_clono <- qread("/mnt/disks/disk/full_dataset_qs/diversity/T_cells_clono.qs")
T_cells_filt_clono_pre<-subset(T_cells_filt_clono, subset=Timepoint=="Pre"&Therapy%in%c("CART","Healthy","TEC","TALQ")&Tissue=="PB"&sub_cluster%in% cd8_clusters)
unique(T_cells_filt_clono_pre$Therapy)
unique(T_cells_filt_clono_pre$Tissue)
unique(T_cells_filt_clono_pre$Disease)
unique(T_cells_filt_clono_pre$Timepoint)


# Estrai i metadati
metadata <- T_cells_filt_clono_pre@meta.data
metadata<-metadata%>%mutate(type2= ifelse(clonotype_size%in%c("Small","Medium","Large"),"Expanded","Rare"))
clusters <- cd8_clusters

results <- list()  # Initialize a list to store results

for (j in seq_along(clusters)) {
 
  # Suppose 'clusters' is a vector/list of cluster IDs and we're inside a loop:
  cluster <- clusters[j]
  
  # 1) Filter metadata for just one sub_cluster
  cluster_data <- metadata %>%
    filter(sub_cluster == cluster)
  
  # 2) Create a frequency table of Disease vs. clone size (type2)
  df <- cluster_data %>%
    dplyr::count(Disease, type2) %>%
    dplyr::rename(n = n)      # Rename for clarity
  
  # 3) Pivot data into wide format, with columns for Expanded and Rare
  df_wide <- df %>%
    pivot_wider(
      names_from  = type2, 
      values_from = n,
      values_fill = 0  # Fill missing combos with 0 instead of NA
    ) %>%
    mutate(
      Total         = Expanded + Rare,
      Freq_Expanded = Expanded / Total,
      Freq_Rare     = Rare / Total
    )
  
  # 4) Calculate the difference in Freq_Expanded between RRMM and HRSMM
  delta_expanded <- 
    df_wide[df_wide$Disease == "RRMM", "Freq_Expanded"] - 
    df_wide[df_wide$Disease == "HRSMM", "Freq_Expanded"]
  
  delta_expanded
  
  # Split cells based on disease
  group_smm <- cluster_data[cluster_data$Disease == "HRSMM", ]
  group_rrmm <- cluster_data[cluster_data$Disease == "RRMM", ]
  
  # Calculate group sizes
  n_smm <- nrow(group_smm)
  n_rrmm <- nrow(group_rrmm)
  n_min<-min(n_smm, n_rrmm)
  # Check if groups are non-empty to avoid sampling errors
  if (n_smm == 0 || n_rrmm == 0) {
    next  # Skip this cluster if any group is empty
  }
  
  # Calculate the observed frequency difference
  observed_diff <- delta_expanded$Freq_Expanded
  
  # Combine groups for sampling
  combined_group <- rbind(group_smm, group_rrmm)
  
  # Initialize a vector for null differences
  null_diffs <- numeric(10000)
  
  # Set seed for reproducibility
  set.seed(123)
  
  # Iterations to build the null distribution
  for (i in 1:10000) {
    # Randomly sample while maintaining original sizes
    smm_sample <- combined_group[sample(nrow(combined_group), n_min, replace = TRUE), ]
    rrmm_sample <- combined_group[sample(nrow(combined_group), n_min, replace = TRUE), ]
    
    # Calculate frequencies
    freq_expanded_sampled_rrmm <- prop.table(table(rrmm_sample$type2))["Expanded"]
    freq_expanded_sampled_smm <- prop.table(table(smm_sample$type2))["Expanded"]
    
    # Store the difference
    null_diffs[i] <- freq_expanded_sampled_smm - freq_expanded_sampled_rrmm
  }
  
  # Store the results for the current cluster
  results[[cluster]] <- list(observed_diff = observed_diff, null_diffs = null_diffs)
}

# Initialize an empty data frame to store all data for plotting
combined_df <- data.frame()

# Extract data for each cluster
for (cluster_name in names(results)) {
  # Get observed difference and null differences for the cluster
  observed_diff <- results[[cluster_name]]$observed_diff
  null_diffs <- results[[cluster_name]]$null_diffs
  
  # Create a data frame for the null distribution
  null_diffs_df <- data.frame(
    difference = null_diffs,
    cluster = cluster_name  )
  
  # Create a single-row data frame for the observed difference
  observed_df <- data.frame(
    difference = observed_diff,
    cluster = cluster_name  )
  
  # Left join the observed difference with the null distribution data
  cluster_df <- left_join(null_diffs_df, observed_df, by = "cluster", suffix = c("_null", "_observed"))
  
  # Stack the results for each cluster
  combined_df <- rbind(combined_df, cluster_df)
}

# Plot with ggplot2, using facet_wrap to separate clusters
combined_df$cluster<-gsub("_","",combined_df$cluster)
combined_df$cluster[combined_df$cluster == "CD8+ZEB2TEM"] <- "CD8+ZEB2+TEM"
combined_df$cluster[combined_df$cluster == "CD8+GZMKTEM"] <- "CD8+GZMK+TEM"
combined_df$cluster[combined_df$cluster == "CD8+GZMBTEM"] <- "CD8+GZMB+TEM"

p1<-ggplot(combined_df, aes(x = difference_null, fill=cluster)) +
  geom_histogram(alpha = 0.3) +
  geom_vline(data = combined_df, 
             aes(xintercept = difference_observed), color = "red", linetype = "dashed", size = 3) +
  facet_wrap(~ cluster, scales = "free") +
  labs(title = "",
       x = "Delta of Expanded Clones Frequency\nRRMM vs HRSMM",
       y = "Density") +
  theme_classic()+
theme(title=element_text(size=40),
      plot.title = element_text(hjust=0.5, face = "bold"), text=element_text(size=50,color = "black"), 
      axis.text.x =element_blank(),legend.position = "bottom",legend.title = element_blank(),
      axis.ticks.x =element_blank(),legend.text = element_text(size=30),
      axis.text.y =element_text(size=40,color = "black"), axis.title.y =element_text(size=50,color = "black"))

p_values <- combined_df %>%
  group_by(cluster) %>%
  summarize(
    observed_diff = unique(difference_observed),  # same for all rows in a cluster
    extreme_count = sum(difference_null >= unique(difference_observed)),
    p_value       = extreme_count / n()          # n() = total null draws for that cluster
  )

p_values
combined_df<-as.data.frame(combined_df)
write_csv(combined_df, paste0(output_dir,"null_distrib_PB_min_size.csv"))
ggsave("~/Immune_project/full_dataset_analysis/plots/diversity/BL_Disease_cluster_exp_NULL_DS_min.pdf",
       plot=p1, width=20, height=15)

# -------------------------------------------------------------------------
#Same for BM
# -------------------------------------------------------------------------
#Theobserved difference will be represented by 
# To check for significant differences in expanded frequency
# Bootstrap SMM+RRMM cell clusters strtified by Disease(10000times)
# within bootstrapped sample calculate average frequency of expanded clones after downsampling(n=100 *100)
# calulate 95%ci for difference in average freq
T_cells_filt_clono_pre<-subset(T_cells_filt_clono, subset=Timepoint=="Pre"&Therapy%in%c("CART","Healthy","TEC","TALQ")&
                                 Tissue=="BM"&sub_cluster%in% cd8_clusters)
unique(T_cells_filt_clono_pre$Therapy)
unique(T_cells_filt_clono_pre$Tissue)
unique(T_cells_filt_clono_pre$Disease)
unique(T_cells_filt_clono_pre$Timepoint)


# Estrai i metadati
metadata <- T_cells_filt_clono_pre@meta.data
metadata<-metadata%>%mutate(type2= ifelse(clonotype_size%in%c("Small","Medium","Large"),"Expanded","Rare"))
clusters <- cd8_clusters

results <- list()  # Initialize a list to store results

for (j in seq_along(clusters)) {
  
  # Suppose 'clusters' is a vector/list of cluster IDs and we're inside a loop:
  cluster <- clusters[j]
  
  # 1) Filter metadata for just one sub_cluster
  cluster_data <- metadata %>%
    filter(sub_cluster == cluster)
  
  # 2) Create a frequency table of Disease vs. clone size (type2)
  df <- cluster_data %>%
    count(Disease, type2) %>%
    rename(n = n)      # Rename for clarity
  
  # 3) Pivot data into wide format, with columns for Expanded and Rare
  df_wide <- df %>%
    pivot_wider(
      names_from  = type2, 
      values_from = n,
      values_fill = 0  # Fill missing combos with 0 instead of NA
    ) %>%
    mutate(
      Total         = Expanded + Rare,
      Freq_Expanded = Expanded / Total,
      Freq_Rare     = Rare / Total
    )
  
  # 4) Calculate the difference in Freq_Expanded between RRMM and HRSMM
  delta_expanded <- 
    df_wide[df_wide$Disease == "RRMM", "Freq_Expanded"] - 
    df_wide[df_wide$Disease == "HRSMM", "Freq_Expanded"]
  
  delta_expanded
  
  # Split cells based on disease
  group_smm <- cluster_data[cluster_data$Disease == "HRSMM", ]
  group_rrmm <- cluster_data[cluster_data$Disease == "RRMM", ]
  
  # Calculate group sizes
  n_smm <- nrow(group_smm)
  n_rrmm <- nrow(group_rrmm)
  
  # Check if groups are non-empty to avoid sampling errors
  if (n_smm == 0 || n_rrmm == 0) {
    next  # Skip this cluster if any group is empty
  }
  
  # Calculate the observed frequency difference
  observed_diff <- delta_expanded$Freq_Expanded
  
  # Combine groups for sampling
  combined_group <- rbind(group_smm, group_rrmm)
  
  # Initialize a vector for null differences
  null_diffs <- numeric(5000)
  
  # Set seed for reproducibility
  set.seed(123)
  
  # Iterations to build the null distribution
  for (i in 1:5000) {
    # Randomly sample while maintaining original sizes
    smm_sample <- combined_group[sample(nrow(combined_group), n_smm, replace = TRUE), ]
    rrmm_sample <- combined_group[sample(nrow(combined_group), n_rrmm, replace = TRUE), ]
    
    # Calculate frequencies
    freq_expanded_sampled_rrmm <- prop.table(table(rrmm_sample$type2))["Expanded"]
    freq_expanded_sampled_smm <- prop.table(table(smm_sample$type2))["Expanded"]
    
    # Store the difference
    null_diffs[i] <- freq_expanded_sampled_smm - freq_expanded_sampled_rrmm
  }
  
  # Store the results for the current cluster
  results[[cluster]] <- list(observed_diff = observed_diff, null_diffs = null_diffs)
}

# Initialize an empty data frame to store all data for plotting
combined_df <- data.frame()

# Extract data for each cluster
for (cluster_name in names(results)) {
  # Get observed difference and null differences for the cluster
  observed_diff <- results[[cluster_name]]$observed_diff
  null_diffs <- results[[cluster_name]]$null_diffs
  
  # Create a data frame for the null distribution
  null_diffs_df <- data.frame(
    difference = null_diffs,
    cluster = cluster_name  )
  
  # Create a single-row data frame for the observed difference
  observed_df <- data.frame(
    difference = observed_diff,
    cluster = cluster_name  )
  
  # Left join the observed difference with the null distribution data
  cluster_df <- left_join(null_diffs_df, observed_df, by = "cluster", suffix = c("_null", "_observed"))
  
  # Stack the results for each cluster
  combined_df <- rbind(combined_df, cluster_df)
}

# Plot with ggplot2, using facet_wrap to separate clusters
combined_df$cluster<-gsub("_","",combined_df$cluster)

p1<-ggplot(combined_df, aes(x = difference_null, fill=cluster)) +
  geom_histogram(alpha = 0.3) +
  geom_vline(data = combined_df, 
             aes(xintercept = difference_observed), color = "red", linetype = "dashed", size = 3) +
  facet_wrap(~ cluster, scales = "free") +
  labs(title = "",
       x = "Delta of Expanded Clones Frequency\nRRMM vs HRSMM",
       y = "Density") +
  theme_classic()+
  theme(title=element_text(size=40),
        plot.title = element_text(hjust=0.5, face = "bold"), text=element_text(size=50,color = "black"), 
        axis.text.x =element_blank(),legend.position = "bottom",legend.title = element_blank(),
        axis.ticks.x =element_blank(),legend.text = element_text(size=30),
        axis.text.y =element_text(size=40,color = "black"), axis.title.y =element_text(size=50,color = "black"))


p_values <- combined_df %>%
  group_by(cluster) %>%
  summarize(
    observed_diff = first(difference_observed),  # same for all rows in a cluster
    extreme_count = sum(difference_null >= first(difference_observed)),
    p_value       = extreme_count / n()          # n() = total null draws for that cluster
  )

p_values
combined_df<-as.data.frame(combined_df)
write_csv(combined_df, paste0(output_dir,"null_distrib_BM.csv"))
ggsave("~/Immune_project/full_dataset_analysis/plots/diversity/BL_BM_Disease_cluster_exp_NULL.pdf",
       plot=p1, width=20, height=15)




