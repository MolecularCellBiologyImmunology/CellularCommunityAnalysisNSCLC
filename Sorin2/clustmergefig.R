library(ggplot2) 
library(dplyr)
library(Rtsne) 
library(clustree) 
library(ggdendro) 
library(tiff) 
library(reshape2) 

date = format(Sys.Date(), "%Y%m%d")

path = getwd()

#out_fig1 = paste0(path, "/Data/full_cd_clus_20_40_49.csv"))
#dir.create(out_dir)

clus = read.csv(paste0(path, "/Data/clustering/full_cd_clus_19_39_48.csv"))

# Cluster tree of neighbour clustering results 
#clus$temp = clus$agglomerateto_49[match(neighb1$cluster, clus$cluster)]

#neighb1$agglom30_average = agglom30$agglomerate_30[match(neighb1$cluster, agglom30$cluster)]

# Subset data and rename columns 
subset = clus %>% 
  select(cluster:agglomerateto_48) %>% 
  dplyr::rename(Communities275 = cluster,
                Communities48 = agglomerateto_48,
                Communities39 = agglomerateto_39,
                Communities19 = agglomerateto_19)

p = clustree(subset, prefix = "Communities", node_size=3, node_text_size = 2, node_text_angle = 60)
ggsave(plot = p, device = "png", width=20, height=10, dpi=300, path = paste0(path, "/Data/clustering_plots/"),
       filename = paste(date, "_clustree_19_39_48.png", sep = ""))

