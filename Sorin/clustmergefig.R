# Cluster tree of neighbour clustering results 
neighb1$agglom30_average = agglom30$agglomerate_30[match(neighb1$cluster, agglom30$cluster)]

# Subset data and rename columns 
subset = neighb1 %>% 
  select(B_cells:cluster, agglom18_average, agglom30_average) %>% 
  dplyr::rename(Communities62 = cluster,
                Communities30 = agglom30_average,
                Communities18 = agglom18_average)

p = clustree(subset, prefix = "Communities")
ggsave(plot = p, device = "png", width=14, height=10, dpi=300, path = out_fig1,
       filename = paste(date, "_clustree_18_30_60_clusters.png", sep = ""))

