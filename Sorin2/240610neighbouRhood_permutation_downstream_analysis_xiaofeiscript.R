

#### Downstream analysis of neighbourhood permutation test ####

## Megan Cole
## 2021

## Scripts required to run previous to this:
# -Normalisation, concatenation and scaling of single cell data
# -Rphenograph clustering of normalised and scaled data
# -Preparation for neighbourhood permuation
# -NeighbouRhood permutation

# Clear working environment
rm(list = ls())

# Load required packages
library(dplyr)
library(ggplot2)
library(textshape)
library(tibble)

############################
##### GLOBAL VARIABLES #####
############################

base = getwd()

# Name of input and output directories:
input_path = file.path(base, "Data/enrichment_input//")
output_path = file.path(base, "Data/enrichment_results//")
input_path




###############################################################################################################

###########################################################
##Xiaofei test code :compare baseline and permutation together 
# data_baseline_perm_Whole = read.csv(paste(input_path, "permutation_results_300_meaned_280624.csv", sep=''))
# # Create a cellID column based on First Object Number and First Image Number
# data_baseline_perm_Whole$cellID = paste(data_baseline_perm_Whole$Patient_ID,
#                                         data_baseline_perm_Whole$source_cluster,
#                                         data_baseline_perm_Whole$target_cluster,
#                                         sep = "_")
# 
# permutation_results_300_dat_p=read.csv(paste(input_path, "permutation_results_300_pvals_280624.csv", sep=''))
# # gorup permutation_results_500_dat_p_mean  to get same observations
# permutation_results_300_dat_p_mean<-permutation_results_300_dat_p %>% group_by(source_cluster,target_cluster,Patient_ID) %>%
#   dplyr::summarise(mean_p = mean(p)) %>% mutate(sig = mean_p < 0.05)
# permutation_results_300_dat_p_mean$cellID = paste(permutation_results_300_dat_p_mean$Patient_ID,
#                                                   permutation_results_300_dat_p_mean$source_cluster,
#                                                   permutation_results_300_dat_p_mean$target_cluster,
#                                                   sep = "_")
# str(data_baseline_perm_Whole)
# # data_baseline_perm_Whole<- data_baseline_perm_Whole %>% mutate(across(c(ROI_ID,community,source_cluster,target_cluster,treatment),as.factor))
# # str(data_baseline_perm_Whole)
# data_baseline_perm_Whole<-data_baseline_perm_Whole %>% group_by(source_cluster,target_cluster,treatment,ROI_ID)
# Calculate log2 fold change between baseline and permutation ct values
# data_baseline_perm_Whole <- data_baseline_perm_Whole %>%
#   group_by(source_cluster, target_cluster, Patient_ID) %>%
#   dplyr::mutate(log2 = log2((mean_obs + 1) / (mean_perm + 1)) + 1)


# Calculate differences in CellID column
# diff_count <- length(unique(setdiff(data_baseline_perm_Whole$cellID, permutation_results_300_dat_p_mean$cellID))) +
#     length(unique(setdiff(permutation_results_300_dat_p_mean$cellID, data_baseline_perm_Whole$cellID)))
# 
# # Check diff_count and proceed accordingly
# if (diff_count == 0) {
#   # Proceed to next step
#   print("Proceeding to the next step...")
# } else {
#   # Print a message indicating diff_count is not 0
#   print(paste("diff_count is :", diff_count))
# }
# # Merge datasets based on CellID column and select specific columns, be careful  checking please. full_join allow variables have different names 
# merged_dat_p <- merge(data_baseline_perm_Whole, permutation_results_300_dat_p_mean, by = "cellID")
# merged_dat_p<- merged_dat_p %>% select(c(cellID,Patient_ID=Patient_ID.x,community=community,source_cluster=source_cluster.x,target_cluster=target_cluster.x,log2,mean_p,sig))
# merged_dat_p<- merged_dat_p %>% mutate(across(c(Patient_ID,source_cluster,target_cluster,community),as.factor))
# str(merged_dat_p)
# # merged_dat_p<- merged_dat_p %>% group_by(ROI_ID,treatment,source_cluster,target_cluster)
#   # Save data with calculations of difference between baseline and permutation ct values
#   write.csv(
#     merged_dat_p,
#     paste(
#       output_path,
#       "/neighbouRhood_dat_perm_mean_log2FC_300",
#       "Whole",
#       ".csv",
#       sep = ""
#     ),
#     row.names = F
  # )
##  Wrong need to check why target cluster are not selected 
#   for (i in 1:length(unique(merged_dat_p$source_cluster))) {
#     select_metacluster = merged_dat_p[which(merged_dat_p$source_cluster == unique(merged_dat_p$source_cluster)[i]),]
#     # select_metacluster <- select_metacluster %>% 
#     #   group_by(ROI_ID,treatment,source_cluster,target_cluster) 
#     #    summarise(mean(log2),mean(mean_p), na.rm = TRUE)
#   #  select_metacluster = merged_dat_p[merged_dat_p$source_cluster==unique(merged_dat_p$source_cluster)[1],]
#     print(unique(select_metacluster$source_cluster))
#     print(unique(select_metacluster$target_cluster))
#     print(unique(merged_dat_p$source_cluster)[i])
#     p = ggplot(select_metacluster,
#                aes(
#                  x = target_cluster,
#                  y = log2,
#                  # colour = factor(treatment, levels = c("MRTX+PD1", "MRTX+PD1+CTLA-4")),
#                  # fill = factor(treatment, levels = c("MRTX+PD1", "MRTX+PD1+CTLA-4"))
#                )) +
#       # geom_point(aes(shape = factor(ROI_ID)), size = 3, alpha = 0.6) +
#       geom_hline(yintercept = 0) +
#       geom_dotplot(
#         data = select_metacluster[which(select_metacluster$sig == FALSE),],
#         binaxis = 'y',
#         stackdir = 'center',
#         dotsize = 0.5,
#         fill = "white"
#       ) +
#       geom_dotplot(
#         data = select_metacluster[which(select_metacluster$sig == TRUE),],
#         binaxis = 'y',
#         stackdir = 'center',
#         dotsize = 0.5
#       ) +
#      
#       theme_minimal() +
#       labs(title = unique(select_metacluster$source_cluster)[i],
#            x = "",
#            y = "Log2FC enrichment") +
#       theme(
#         axis.text.x = element_text(
#           angle = 45,
#           hjust = 1,
#           size = 14
#         ),
#         axis.text.y = element_text(size = 14),
#         axis.title = element_text(size = 14),
#         plot.title = element_text(size = 16),
#         legend.title = element_blank(),
#         legend.text = element_text(size = 14)
#       ) +
#       theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
#       scale_x_discrete(
#         limits = c(
#           'Cancer', 'Cl_Mac', 'Th', 'B', 'Tc', 'Treg', 'Endothelial',
#           'Alt_Mac', 'Neutrophil', 'Cl_Mo', 'Unknown', 'NK', 'T_other',
#           'Non-Cl_Mo', 'Int_Mo', 'Mast', 'DC'
#         )
#       )
#     
#     # if (unique(merged_dat_p$source_cluster)[i] == "Macrophages type 1" |
#     #     unique(merged_dat_p$source_cluster)[i] == "Macrophages type 2") {
#     #   print("Setting y limits")
#     #   p = p + ylim(-4, 2.5)
#     # }
#     print(unique(merged_dat_p$source_cluster)[i])
#     filename = paste(
#       unique(merged_dat_p$source_cluster)[i],
#       "_neighbours_pointPlot_log2significance_perTreatment.png",
#       sep = ""
#     )
#      ggsave(
#       plot = p,
#       device = "png",
#       width = 30,
#       height = 20,
#       dpi = 300,
#       bg = 'White',
#       path = paste(
#         output_path,
#         "output_plots/point/new/",
#         "Whole",
#         "_v_wholeTissue/",
#         sep = ""
#       ),
#       filename = filename
#     )
#   }
# print("log2FC plotting complete")
# 

#######################
#### Run functions ####

# Cell types of interest
cellTypes = c(
  'Cancer', 'Cl_Mac', 'Th', 'B', 'Tc', 'Treg', 'Endothelial',
  'Alt_Mac', 'Neutrophil', 'Cl_Mo', 'Unknown', 'NK', 'T_other',
  'Non-Cl_Mo', 'Int_Mo', 'Mast', 'DC'
)

  
  
experiment_path = "Data/enrichment_input/"
path = experiment_path
# Call the function to prep the data & then run log2FC calculation and plotting
# data_prep(experiment_path, cellTypes)
######################################################################################

#### Comparing baseline and permutation neighbour data ##
data_prep = function(path, select) {
  # Load baseline and permutation data for each treatment group
  data_baseline_meaned = read.csv(
    paste(
      path,
      "permutation_results_300_meaned_280624.csv",
      sep = ""
    )
  )
  dat_perm_pvals = read.csv(
    paste(
      path,
      "permutation_results_300_pvals_280624.csv",
      sep = ""
    )
  )
  # dat_perm_ori = read.csv(
  #   paste(
  #     path,
  #     "permutation_300_original_values_280624.csv",
  #     sep = ""
  #   )
  # )
  
  # Create unique ID based on pairing neighbour combination, group and treatment for pvals and mean of the permutation runs
  dat_perm_pvals$ID = paste(
    dat_perm_pvals$Patient_ID,
    dat_perm_pvals$community,
    dat_perm_pvals$source_cluster,
    dat_perm_pvals$target_cluster,
    sep = "_"
  )
  data_baseline_meaned$ID = paste(
    data_baseline_meaned$Patient_ID,
    data_baseline_meaned$community,
    data_baseline_meaned$source_cluster,
    data_baseline_meaned$target_cluster,
    sep = "_"
  )
  
  #### Add p-value data to dat_perm for plotting ####
  # Create a unique ID based on first and second labels, image number and treatment
  data_meaned_pvals <- merge(data_baseline_meaned, dat_perm_pvals[,-c(1,2,3,4,5)], by = "ID")
  view(data_meaned_pvals)
  
  
  # Add significance values from dat_p to dat_perm_mean
# 
  dat_perm_mean = subset(dat_perm_mean, FirstLabel %in% select &
                           SecondLabel %in% select)
  dim(data_meaned_pvals)
  
  print("Data prep function complete")
  
  # Call log2FC function
  log2FC(path, data_meaned_pvals)
}

##################################
## Calculating log2 fold change and plotting volcanoes ##
data = data_meaned_pvals
data$log2 = log2(data$mean_obs / data$mean_perm)


# get rid of inf and -inf in log2 by removing 0 counts
celltypes_filter = c("Tc", "Th", 'B', 'Endothelial', 'Cl_Mac')
community_filter = c(26, 27)


data <- data[data$mean_obs != 0, ]
data <- data[data$mean_perm != 0, ]

# data <- data[data$source_cluster != 0, ]

data_tcells = subset(data, source_cluster %in% celltypes_filter &
                         target_cluster %in% celltypes_filter)
data_tcells = subset(data_tcells, community %in% community_filter)
data_tcells$community <- as.character(data_tcells$community)



# Modify 'sig' column to show "Significant" and "Non-significant"
data_tcells$sig[data_tcells$sig == "True"] <- "Significant"
data_tcells$sig[data_tcells$sig == "False"] <- "Not significant"

# Convert p-values to -log10 scale
data_tcells <- data_tcells %>%
  mutate(negLogP = -log10(p))

# view(data_tcells)

data_tcells$new_ID = paste(
  data_tcells$source_cluster,
  data_tcells$target_cluster,
  sep = " - "
)

# Define thresholds
x_threshold <- 0.3
y_threshold <- -log10(0.01)

# Calculate the counts above both thresholds for each facet and community
data_summary_pos <- data_tcells %>%
  group_by(new_ID, community) %>%
  summarise(
    count_above_pos = sum(log2 > x_threshold & negLogP > y_threshold)
  ) %>%
  ungroup() %>%
  mutate(label = paste("Pos. Hits:", count_above_pos))

# Calculate the counts above both thresholds for each facet and community
data_summary_neg <- data_tcells %>%
  group_by(new_ID, community) %>%
  summarise(
    count_above_neg = sum(-log2 > x_threshold & negLogP > y_threshold)
  ) %>%
  ungroup() %>%
  mutate(label = paste("Neg. Hits:", count_above_neg))


# Volcano plot with faceting by group
volcano_plot <- ggplot(data_tcells, aes(x = log2, y = negLogP)) +
  geom_point(aes(color = community), size = 0.1) +
  scale_color_manual(values = c('26' = "red", '27' = "blue")) +
  theme_minimal() +
  labs(
    title = "Volcano Plot",
    x = "Log2 Fold Change",
    y = "-Log10(P-value)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_blank()
  ) +
facet_wrap(~ new_ID, scales = "fixed") +
geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "black") + 
geom_vline(xintercept = 0.3, linetype = "dashed", color = "black") +
geom_vline(xintercept = -0.3, linetype = "dashed", color = "black")+
  geom_text(data = data_summary_pos[data_summary_pos$community == 26, ], aes(x = 1.5, y = 0.5, label = label),
            color = 'red', hjust = 1.1, vjust = 1.1, size = 2, inherit.aes = FALSE) +
  geom_text(data = data_summary_pos[data_summary_pos$community == 27, ], aes(x = 1.5, y = 0.2, label = label),
            color = 'blue', hjust = 1.1, vjust = 1.1, size = 2, inherit.aes = FALSE)+
  
  geom_text(data = data_summary_neg[data_summary_neg$community == 26, ], aes(x = -1.2, y = 0.5, label = label),
            color = 'red', hjust = 1.1, vjust = 1.1, size = 2, inherit.aes = FALSE) +
  geom_text(data = data_summary_neg[data_summary_neg$community == 27, ], aes(x = -1.2, y = 0.2, label = label),
            color = 'blue', hjust = 1.1, vjust = 1.1, size = 2, inherit.aes = FALSE)
# data_summary[data_summary$community == "27"]


# table(data_tcells$source_cluster)
# data_summary$community == "26"
ggsave(
  filename =   paste(
    output_path,
    "volcano_2627_cellfilter_hits",
    ".jpg",
    sep = ""
  ),
  plot = volcano_plot,
  width = 10,  # Width in inches
  height = 8,  # Height in inches
  dpi = 300    # Resolution in dots per inch
)
# write.csv(
#   data,
#   paste(
#     output_path,
#     "neighbouRhood_dat_perm_mean_log2FC",
#     ".csv",
#     sep = ""
#   ),
#   row.names = F
# )

# Volcano plot
ggplot(data, aes(x = log2, y = negLogP)) +
  geom_point(aes(color = sig)) +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey")) +
  theme_minimal() +
  labs(
    title = "Volcano Plot",
    x = "Log2 Fold Change",
    y = "-Log10(P-value)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_blank()
  )





##################################
## Calculating log2 fold change ##
path = experiment_path
log2FC = function(path, data) {
  # Calculate log2 fold change between baseline and permutation ct values
  data = data_meaned_pvals
  data$log2 = log2(data$mean_obs / data$mean_perm)
  
  # Save data with calculations of difference between baseline and permutation ct values
  # write.csv(
  #   data,
  #   paste(
  #     output_path,
  #     "neighbouRhood_dat_perm_mean_log2FC",
  #     ".csv",
  #     sep = ""
  #   ),
  #   row.names = F
  # )
  data = read.csv(
    paste(
      output_path,
      "neighbouRhood_dat_perm_mean_log2FC.csv",
      sep = ""
    )
  )

  # get rid of inf and -inf in log2 by removing 0 counts
  data <- data[data$mean_obs != 0, ]
  data <- data[data$mean_perm != 0, ]
  
  # Modify 'sig' column to show "Significant" and "Non-significant"
  data$sig[data$sig == "True"] <- "Significant"
  data$sig[data$sig == "False"] <- "Non-significant"
  
  
  # 
  # # Bin the values into ranges of 0.01
  # data$log2_binned <- floor(data$log2 / 0.1) * 0.1
  # 
  # value_counts <- table(data$log2_binned)
  # 
  # cat("\nValue Counts in mean_obs\n")
  # print(value_counts)
  # 

  for (i in 1:length(unique(data$source_cluster))) {
    select_metacluster = data[which(data$source_cluster == unique(data$source_cluster)[i]),]
    print(unique(data$source_cluster)[i])
    # 
    value_counts <- table(select_metacluster$sig)

    cat("\nValue Counts in mean_obs\n")
    print(value_counts)
    # 
    # print(select_metacluster[which(select_metacluster$sig == "True"),])
          
    p = ggplot(select_metacluster,
               aes(
                 x = target_cluster,
                 y = log2,
                 colour = factor(sig, levels = c("Significant", "Non-significant")),
                 # fill = factor(sig, levels = c("True", "False"))
               )) +
      geom_hline(yintercept = 0) +
      geom_violin(
        data = select_metacluster,
        binaxis = 'x',
        width = 1,
        stackdir = 'center',
        trim='False',
        bw = 0.5,
        alpha=0.8,
        # dotsize = 0.5,
        # colour = "blue",
        # alpha=0.1,
        scale="count"
      ) +
      # geom_violin(
      #   data = select_metacluster[which(select_metacluster$sig == "True"),],
      #   binaxis = 'y',
      #   stackdir = 'center',
      #   colour='red',
      #   alpha=0.1,
      #   scale="count"
      #   
      #   # dotsize = 0.5
      # ) +
    
      # geom_violin(
      #   data = select_metacluster[which(select_metacluster$sig == "False"),],
      #   binaxis = 'y',
      #   stackdir = 'center',
      #   # dotsize = 0.5,
      #   colour = "blue",
      #   alpha=0.1,
      #   scale="count"
      # ) +

      theme_minimal() +
    
      labs(title = unique(data$source_cluster)[i],
           x = "",
           y = "Log2FC enrichment", 
           color = "Significance") +
      theme(
        panel.spacing.x = unit(2, "cm"),  # increase horizontal space between violins
        
        axis.text.x = element_text(
          angle = 45,
          hjust = 1,
          size = 14
        ),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16),
        legend.title = element_blank(),
        legend.text = element_text(size = 14)
      ) +
      theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
      scale_x_discrete(
        limits = c(
          'Cancer', 'Cl_Mac', 'Th', 'B', 'Tc', 'Treg', 'Endothelial',
          'Alt_Mac', 'Neutrophil', 'Cl_Mo', 'Unknown', 'NK', 'T_other',
          'Non-Cl_Mo', 'Int_Mo', 'Mast', 'DC'
        )
      )
    # 
    # if (unique(data$source_cluster)[i] == "Macrophages type 1" |
    #     unique(data$source_cluster)[i] == "Macrophages type 2") {
    #   print("Setting y limits")
    #   p = p + ylim(-4, 2.5)
    # }
    
    filename = paste(
      unique(data$source_cluster)[i],
      "_neighbours_violinplot_log2significance.png",
      sep = ""
    )
    ggsave(
      plot = p,
      device = "png",
      width = 10,
      height = 5,
      dpi = 300,
      bg = 'White',
      path = paste(
        output_path,
        sep = ""
      ),
      filename = filename
    )
  }
  print("log2FC plotting complete")
}

#######################
#### Run functions ####

# Cell types of interest
cellTypes = c(
  'Cancer', 'Cl_Mac', 'Th', 'B', 'Tc', 'Treg', 'Endothelial',
  'Alt_Mac', 'Neutrophil', 'Cl_Mo', 'Unknown', 'NK', 'T_other',
  'Non-Cl_Mo', 'Int_Mo', 'Mast', 'DC'
)

# Call the function to prep the data & then run log2FC calculation and plotting
data_prep(experiment_path, cellTypes)
# Domain options: 'normal', 'interface', 'tumour', 'wholeTissue'


## END 
################################################################





