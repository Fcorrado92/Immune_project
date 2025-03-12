library(qs)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)

# Load data
T_cells_filt_clono <- qread("/mnt/disks/disk/full_dataset_qs/diversity/T_cells_clono.qs")
T_cells_filt_clono_pre <- subset(
  T_cells_filt_clono,
  subset = Timepoint == "Pre" & Therapy %in% c("CART","Healthy","TEC","TALQ") &
    Tissue == "PB" & sub_cluster %in% cd8_clusters
)

# Check unique values
unique(T_cells_filt_clono_pre$Therapy)
unique(T_cells_filt_clono_pre$Tissue)
unique(T_cells_filt_clono_pre$Disease)
unique(T_cells_filt_clono_pre$Timepoint)

# Extract metadata and create new variable 'type2'
metadata <- T_cells_filt_clono_pre@meta.data %>%
  mutate(type2 = ifelse(clonotype_size %in% c("Small", "Medium", "Large"), "Expanded", "Rare"))

clusters <- cd8_clusters  # Your vector of clusters

# Define the nmin values you wish to try.
# (Make sure these values make sense given your data; e.g., if one cluster has only 200 cells,
# you might choose nmin values up to 200.)
nmin_values <- c(50, 100, 150, 200)

# Loop over each nmin value
for (n_min in nmin_values) {
  
  results <- list()  # Initialize a list to store results for this nmin
  
  for (j in seq_along(clusters)) {
    
    cluster <- clusters[j]
    
    # Filter metadata for this cluster
    cluster_data <- metadata %>%
      filter(sub_cluster == cluster)
    
    # Create frequency table (Disease vs. clone size)
    df <- cluster_data %>%
      dplyr::count(Disease, type2) %>%
      rename(n = n)
    
    # Pivot data to wide format
    df_wide <- df %>%
      pivot_wider(
        names_from  = type2, 
        values_from = n,
        values_fill = 0
      ) %>%
      mutate(
        Total         = Expanded + Rare,
        Freq_Expanded = Expanded / Total,
        Freq_Rare     = Rare / Total
      )
    
    # Calculate the difference in Freq_Expanded between RRMM and HRSMM
    delta_expanded <- df_wide[df_wide$Disease == "RRMM", "Freq_Expanded"] -
      df_wide[df_wide$Disease == "HRSMM", "Freq_Expanded"]
    
    # Split cells based on disease
    group_smm <- cluster_data %>% filter(Disease == "HRSMM")
    group_rrmm <- cluster_data %>% filter(Disease == "RRMM")
    
    n_smm <- nrow(group_smm)
    n_rrmm <- nrow(group_rrmm)
    
    # Skip cluster if any group is empty or has fewer cells than n_min
    if (n_smm < n_min | n_rrmm < n_min) {
      next  # Skip this cluster if the available cells are less than n_min
    }
    
    observed_diff <- delta_expanded$Freq_Expanded
    
    # Combine groups for sampling
    combined_group <- bind_rows(group_smm, group_rrmm)
    
    null_diffs <- numeric(10000)  # Number of bootstrap iterations
    set.seed(123)  # For reproducibility
    
    for (i in 1:10000) {
      # Sample n_min cells from each group (with replacement)
      smm_sample <- combined_group[sample(nrow(combined_group), n_min, replace = TRUE), ]
      rrmm_sample <- combined_group[sample(nrow(combined_group), n_min, replace = TRUE), ]
      
      # Calculate frequencies (handle the case where "Expanded" might be missing)
      freq_expanded_sampled_rrmm <- {
        tmp <- prop.table(table(rrmm_sample$type2))
        if ("Expanded" %in% names(tmp)) tmp["Expanded"] else 0
      }
      freq_expanded_sampled_smm <- {
        tmp <- prop.table(table(smm_sample$type2))
        if ("Expanded" %in% names(tmp)) tmp["Expanded"] else 0
      }
      
      null_diffs[i] <- freq_expanded_sampled_smm - freq_expanded_sampled_rrmm
    }
    
    # Store results for this cluster
    results[[cluster]] <- list(observed_diff = observed_diff, null_diffs = null_diffs)
  }  # end cluster loop
  
  # Combine results for all clusters into a single data frame
  combined_df <- data.frame()
  for (cluster_name in names(results)) {
    observed_diff <- results[[cluster_name]]$observed_diff
    null_diffs <- results[[cluster_name]]$null_diffs
    
    null_diffs_df <- data.frame(
      difference = null_diffs,
      cluster = cluster_name
    )
    
    observed_df <- data.frame(
      difference = observed_diff,
      cluster = cluster_name
    )
    
    cluster_df <- left_join(null_diffs_df, observed_df, by = "cluster", 
                            suffix = c("_null", "_observed"))
    combined_df <- bind_rows(combined_df, cluster_df)
  }
  
  # Clean up cluster names for plotting
  combined_df$cluster <- gsub("_", "", combined_df$cluster)
  combined_df$cluster[combined_df$cluster == "CD8+ZEB2TEM"] <- "CD8+ZEB2+TEM"
  combined_df$cluster[combined_df$cluster == "CD8+GZMKTEM"] <- "CD8+GZMK+TEM"
  combined_df$cluster[combined_df$cluster == "CD8+GZMBTEM"] <- "CD8+GZMB+TEM"
  
  # Generate the plot
  p1 <- ggplot(combined_df, aes(x = difference_null, fill = cluster)) +
    geom_histogram(alpha = 0.3, bins = 30) +
    geom_vline(data = combined_df, 
               aes(xintercept = difference_observed), color = "red", linetype = "dashed", size = 1.5) +
    facet_wrap(~ cluster, scales = "free") +
    labs(title = paste("Null Distribution (nmin =", n_min, ")"),
         x = "Delta of Expanded Clones Frequency\n(RRMM vs HRSMM)",
         y = "Density") +
    theme_classic() +
    theme(title = element_text(size = 40),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          text = element_text(size = 50, color = "black"),
          axis.text.x = element_blank(),
          legend.position = "bottom",
          legend.title = element_blank(),
          axis.ticks.x = element_blank(),
          legend.text = element_text(size = 30),
          axis.text.y = element_text(size = 40, color = "black"),
          axis.title.y = element_text(size = 50, color = "black"))
  
  # Save CSV and Plot for current n_min
  csv_filename <- paste0(output_dir, "null_distrib_PB_nmin", n_min, ".csv")
  plot_filename <- paste0("~/Immune_project/full_dataset_analysis/plots/diversity/BL_Disease_cluster_exp_NULL_nmin", n_min, ".pdf")
  
  write_csv(combined_df, csv_filename)
  ggsave(plot_filename, plot = p1, width = 20, height = 15)
  
  # Optionally, print p-values per cluster
  p_values <- combined_df %>%
    group_by(cluster) %>%
    summarize(
      observed_diff = unique(difference_observed),
      extreme_count = sum(difference_null >= unique(difference_observed)),
      p_value = extreme_count / n()  # n() is the number of bootstrap draws for that cluster
    )
  
  print(paste("nmin =", n_min))
  print(p_values)
}
