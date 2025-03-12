#Import dependencies
library(ggpubr)
library(ggrepel)
library(dplyr)
library(stringr)
library(ggplot2)
library(patchwork)       # For wrap_plots
library(rstatix)         # For geom_pwc() and ggadjust_pvalue()

# 1) Read in your data (example)
output_dir <- "/mnt/disks/disk/full_dataset_qs/diversity/"

out_mean<-read_csv(paste0(output_dir, "Diversity_Means_jan2025.csv"))
meta <- read_excel("~/Immune_project/full_dataset_analysis/Meta_Immune_proj_feb25.xlsx")
meta<-meta%>%distinct(Sample_ID,.keep_all = TRUE)

#ADD METADATA<- -2 sampples with less than 100 cells with clonotypes
out_mean<-left_join(out_mean, meta[,c("Sample_ID","BOR","Timepoint","Tissue", "Therapy","Disease","ID")], by="Sample_ID")

out_mean<-out_mean%>%mutate(ter_Tis=paste0(Therapy,Tissue))
out_mean<-out_mean%>%filter(!ter_Tis=="CARTBM")
out_mean<-out_mean%>%filter(Timepoint%in%c("Pre","Post","Day30"))
out_mean$Timepoint[out_mean$Timepoint=="Day30"]<-"Post"
out_mean <- out_mean %>%
  mutate(ID_Tissue = paste0(ID, "_", Tissue))

df <- out_mean %>%
  filter(Therapy %in% c("TEC","TALQ","CART")) %>%
  group_by(ID_Tissue,ID, Therapy, Tissue, Disease) %>%
  summarise(n = n_distinct(Timepoint), .groups = "drop")

df3<- df %>%
  filter(n > 1) 
length(unique(df3$ID))

paired_ids <- df %>%
  filter(n > 1) %>%
  distinct(ID_Tissue) %>%
  pull()
length(unique(paired_ids))

paired_counts <- out_mean %>%
  filter(ID_Tissue %in% paired_ids)

df2 <- paired_counts %>%
  filter(Therapy %in% c("TEC", "TALQ", "CART")) %>%
  group_by(ID, Timepoint) %>%
  summarise(n = n_distinct(Tissue), .groups = "drop")


# ------------------------------------------------------------------------------
# 2) Define different "response classification" SCENARIOS.
#    For each scenario, you’ll specify how RRMM_R, RRMM_NR, SMM_R, etc., are defined.
# ------------------------------------------------------------------------------

classification_scenarios <- list(
  # Scenario 1: your original classification
  PR_NR = function(BOR, Disease, Therapy) {
    # Just copy your original logic:
    case_when(
      BOR %in% c("PR","MR","SD","PD") & Disease == "RRMM" ~ "RRMM_NR",   # RRMM Non-Responder
      Disease == "RRMM" & !BOR %in% c("PR","MR","SD","PD") ~ "RRMM_R",   # RRMM Responder
      Disease == "HRSMM" & Therapy != "Len"                ~ "HRSMM_R",  # SMM Responder
      Disease == "Healthy"                                 ~ "Healthy_NA",
      TRUE                                                 ~ "SMM_NA"
    )
  },
  
  # Scenario 2: for example, if now you want RRMM_NR to include "VGPR" as well
  VGPR_NR = function(BOR, Disease, Therapy) {
    case_when(
      BOR %in% c("VGPR","PR","MR","SD","PD") & Disease == "RRMM" ~ "RRMM_NR",  
      Disease == "RRMM" & !BOR %in% c("VGPR","PR","MR","SD","PD") ~ "RRMM_R",
      Disease == "HRSMM" & Therapy != "Len"                ~ "HRSMM_R",
      Disease == "Healthy"                                 ~ "Healthy_NA",
      TRUE                                                 ~ "SMM_NA"
    )
  },
  
  MR_NR = function(BOR, Disease, Therapy) {
    case_when(
      BOR %in% c("MR","SD","PD") & Disease == "RRMM" ~ "RRMM_NR",  
      Disease == "RRMM" & !BOR %in% c("MR","SD","PD") ~ "RRMM_R",
      Disease == "HRSMM" & Therapy != "Len"                ~ "HRSMM_R",
      Disease == "Healthy"                                 ~ "Healthy_NA",
      TRUE                                                 ~ "SMM_NA"
    )
  }
)

# ------------------------------------------------------------------------------
# 3) Define different therapy sets you want to iterate over
# ------------------------------------------------------------------------------
therapy_sets <- list(
  only_TEC          = c("TEC"),
  TEC_and_TALQ       = c("TEC","TALQ"),
  TEC_TALQ_and_CART  = c("TEC","TALQ","CART")
)

# ------------------------------------------------------------------------------
# 4) Loop over each combination of (scenario, therapy set), do your filtering,
#    classification, plotting, and saving
# ------------------------------------------------------------------------------
for (sc_name in names(classification_scenarios)) {
  
  # Extract the classification function for this scenario
  classify_fun <- classification_scenarios[[sc_name]]
  
  for (th_name in names(therapy_sets)) {
    selected_therapies <- therapy_sets[[th_name]]
    
    # Filter only the relevant therapy(ies)
    PBMC_sign <- out_mean %>%
      filter(Therapy %in% selected_therapies & ID_Tissue %in% paired_ids)
    
    # Factor ordering for Timepoint
    PBMC_sign$Timepoint <- factor(PBMC_sign$Timepoint, levels = c("Pre","Post"))
    
    # Apply scenario-specific classification
    PBMC_sign <- PBMC_sign %>%
      mutate(
        response_disease = classify_fun(BOR, Disease, Therapy)
      ) %>%
      mutate(
        response_disease = factor(response_disease, 
                                  levels=c("HRSMM_R","RRMM_R","RRMM_NR","Healthy_NA","SMM_NA"))
      )
    
      p <- ggplot(PBMC_sign, aes(x=Timepoint, y=1/Chao, fill=Timepoint)) +
        geom_boxplot(size=3, alpha=0.6, outlier.shape = NA) +
        geom_line(aes(group=ID)) +
        geom_point(
          size=10, 
          aes(fill=Timepoint, stroke=2),
          color="black",shape=21, alpha=0.5,
          position=position_dodge(width=0.5)
        ) +
        labs(
          y="1/Chao", 
        ) +
        theme_classic() + 
        guides(color=FALSE) +
        scale_fill_manual(values=c("gray","purple")) +
        facet_wrap(~ response_disease, scales="free") +
        theme(
          plot.margin = unit(c(0.5, 0, 0.2, 0.1), "cm"),
          plot.title = element_text(hjust=0.5, face="bold", size=20),
          text = element_text(size=30, color="black"),
          legend.position = "bottom",
          axis.ticks.x = element_blank()        ) +
        scale_y_continuous(expand = expansion(mult=c(0.1,0.1))) +
        
        # Pairwise comparisons: Wilcoxon
        geom_pwc(label.size = 4, method="wilcox.test",
                 step.increase=0.08, tip.length=0.01, size=0.8) 
      
      # Adjust p-values with BH, label them with "q = ..."
      p <- ggadjust_pvalue(
        p,
        p.adjust.method = "BH",
        label = paste("q=", "{p.adj.format}"),
        hide.ns = "p"
      )
      
    # 5) Save your combined plot. 
    #    Make the filename reflect both the scenario and the therapy set
      out_filename <- paste0(
        "~/Immune_project/full_dataset_analysis/plots/boxplots/longitudinal_Tcells_response/",
        "diversity_",
        sc_name,            # scenario
        "_", th_name,       # therapy set
        "_Tcells_boxplot.pdf"
      )
      
    ggsave(
      filename = out_filename,
      plot     = p,
      width    = 10,   # adjust as needed
      height   = 10
    )
    
    message("Saved plot: ", out_filename)
    
  } # end therapy loop

}






