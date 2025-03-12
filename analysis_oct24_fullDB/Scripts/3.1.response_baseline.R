library(tidyverse)
library(readxl)
library(ggpubr)      # For geom_pwc() and ggadjust_pvalue()
library(patchwork)   # For wrap_plots()

# 1) Load and process the data
Tcells_counts <- read_csv("~/Immune_project/full_dataset_analysis/compositional_analysis/SC_PBMC_cell_counts_over500.csv")
Tcells_counts <- Tcells_counts[-c(6:12)]
Meta <- read_xlsx("~/Immune_project/full_dataset_analysis/Meta_Immune_proj_feb25.xlsx")
meta <- Meta %>% distinct(Sample_ID, .keep_all = TRUE)
Tcells_counts <- left_join(Tcells_counts, meta, by = "Sample_ID")


# Create an ID_Tissue column for further grouping
Tcells_counts <- Tcells_counts %>%
  mutate(ID_Tissue = paste0(ID, "_", Tissue))

# Group counts per sample, therapy, tissue, and disease
df <- Tcells_counts %>%
  filter(Therapy %in% c("TEC", "TALQ", "CART")) %>%
  group_by(ID_Tissue, ID, Therapy, Tissue, Disease) %>%
  summarise(n = n_distinct(Timepoint), .groups = "drop")

# Filter to get only samples with more than one Timepoint
df3 <- df %>% filter(n > 1)
paired_ids <- df %>% filter(n > 1) %>% distinct(ID_Tissue) %>% pull()
paired_counts <- Tcells_counts %>% filter(ID_Tissue %in% paired_ids)

df2 <- paired_counts %>%
  filter(Therapy %in% c("TEC", "TALQ", "CART")) %>%
  group_by(ID,Therapy,Tissue, Timepoint,Disease) %>%
  summarise(n = n_distinct(Tissue), .groups = "drop")

# Get the sub-clusters of interest (those containing "CD8")
subsets <- unique(grep("CD8", paired_counts$sub_cluster, value = TRUE))

# ------------------------------------------------------------------------------
# 2) Define different "response classification" scenarios.
# ------------------------------------------------------------------------------
classification_scenarios <- list(
  PR_NR = function(BOR, Disease, Therapy) {
    case_when(
      BOR %in% c("PR", "MR", "SD", "PD") & Disease == "RRMM" ~ "RRMM_NR",
      Disease == "RRMM" & !BOR %in% c("PR", "MR", "SD", "PD") ~ "RRMM_R",
      Disease == "HRSMM" & Therapy != "Len" ~ "HRSMM_R",
      Disease == "Healthy" ~ "Healthy_NA",
      TRUE ~ "SMM_NA"
    )
  },
  VGPR_NR = function(BOR, Disease, Therapy) {
    case_when(
      BOR %in% c("VGPR", "PR", "MR", "SD", "PD") & Disease == "RRMM" ~ "RRMM_NR",
      Disease == "RRMM" & !BOR %in% c("VGPR", "PR", "MR", "SD", "PD") ~ "RRMM_R",
      Disease == "HRSMM" & Therapy != "Len" ~ "HRSMM_R",
      Disease == "Healthy" ~ "Healthy_NA",
      TRUE ~ "SMM_NA"
    )
  },
  MR_NR = function(BOR, Disease, Therapy) {
    case_when(
      BOR %in% c("MR", "SD", "PD") & Disease == "RRMM" ~ "RRMM_NR",
      Disease == "RRMM" & !BOR %in% c("MR", "SD", "PD") ~ "RRMM_R",
      Disease == "HRSMM" & Therapy != "Len" ~ "HRSMM_R",
      Disease == "Healthy" ~ "Healthy_NA",
      TRUE ~ "SMM_NA"
    )
  }
)

# ------------------------------------------------------------------------------
# 3) Define different therapy sets you want to iterate over.
# ------------------------------------------------------------------------------
therapy_sets <- list(
  only_TEC = c("TEC"),
  TEC_and_TALQ = c("TEC", "TALQ"),
  CART = c( "CART")
)

# ------------------------------------------------------------------------------
# 4) Loop over each combination of (scenario, therapy set), do your filtering,
#    classification, plotting, and saving.
# ------------------------------------------------------------------------------
for (sc_name in names(classification_scenarios)) {
  # Extract the classification function for this scenario
  classify_fun <- classification_scenarios[[sc_name]]
  for (th_name in names(therapy_sets)) {
    selected_therapies <- therapy_sets[[th_name]]
    # Filter only the relevant therapies and sub-clusters
    PBMC_sign <- paired_counts %>%
      filter(Therapy %in% selected_therapies & sub_cluster %in% subsets & ID_Tissue %in% paired_ids)
    
    # Set factor ordering for Timepoint
    PBMC_sign$Timepoint <- factor(PBMC_sign$Timepoint, levels = c("Pre", "Post","Day14", "Day30","Day56","Day100"))
    
    # Apply scenario-specific classification, forcing ID "PT17" to "RRMM_NR"
    PBMC_sign <- PBMC_sign %>%
      mutate(
        response_disease = classify_fun(BOR, Disease, Therapy)
      ) %>%
      mutate(
        response_disease = factor(response_disease,
                                  levels=c("HRSMM_R","RRMM_R","RRMM_NR","Healthy_NA","SMM_NA"))
      )
    
    # Create an additional grouping column if needed
    PBMC_sign <- PBMC_sign %>% mutate(ID_sub = paste0(sub_cluster, "_", ID))
    
    # Now create a plot for each sub-cluster in 'subsets'
    plot_list <- list()
    for (i in seq_along(subsets)) {
      sub_data <- PBMC_sign %>% filter(sub_cluster == subsets[[i]])
      p <- ggplot(sub_data, aes(x = Timepoint, y = Freq, fill = Timepoint)) +
        geom_boxplot(size = 3, alpha = 0.6, outlier.shape = NA) +
        geom_line(aes(group = ID)) +
        geom_point(
          size = 10,
          aes(fill = Timepoint, stroke = 2),
          color = "black", shape = 21, alpha = 0.5,
          position = position_dodge(width = 0.5)
        ) +
        labs(
          y = "Tcells(%)",
          title = paste0(subsets[[i]])
        ) +
        theme_classic() +
        guides(color = FALSE) +
        # Split the panels by Tissue (rows) and response_disease (columns)
        facet_grid(Tissue ~ response_disease, scales = "free") +
        theme(
          plot.margin = unit(c(0.5, 0, 0.2, 0.1), "cm"),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
          text = element_text(size = 30, color = "black"),
          legend.position = "bottom",
          axis.ticks.x = element_blank(),
          axis.text.x  = element_text(angle=45,vjust=1, hjust=1)
          
        ) +
        scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
        guides(fill = FALSE) +
        # Pairwise comparisons: Wilcoxon (Mann-Whitney test)
        geom_pwc(label.size = 10, method = "wilcox.test",
                 step.increase = 0.08, tip.length = 0.01, size = 1)
      
      p <- ggadjust_pvalue(
        p,
        p.adjust.method = "BH",
        label = paste("p=", "{p.format}")
      )
      
      plot_list[[i]] <- p
    }
    
    # Combine all plots into one
    combined_plot <- wrap_plots(plot_list, ncol = 4)
    
    # Save the combined plot with a filename reflecting scenario and therapy set
    out_filename <- paste0(
      "~/Immune_project/full_dataset_analysis/plots/boxplots/longitudinal_Tcells_response/",
      "longitudinal_",
      sc_name, "_", th_name, "_Tcells_boxplot.pdf"
    )
    ggsave(
      filename = out_filename,
      plot = combined_plot,
      width = 55,
      height = 20,limitsize = FALSE
    )
    message("Saved plot: ", out_filename)
  } # end therapy loop
} # end scenario loop
