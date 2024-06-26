library(ggplot2) 
library(dplyr)
library(Rtsne) 
library(clustree) 
library(ggdendro) 
library(tiff) 
library(reshape2) 

date = format(Sys.Date(), "%Y%m%d")

path = getwd()

#out_fig1 = paste0(path, "/Data/full_cd_clus_21_40_49.csv"))
#dir.create(out_dir)

clus = read.csv(paste0(path, "/data/agglo_clus/aggloclus_clus18.csv"))

# Cluster tree of neighbour clustering results 
#clus$temp = clus$agglomerateto_49[match(neighb1$cluster, clus$cluster)]

#neighb1$agglom30_average = agglom30$agglomerate_30[match(neighb1$cluster, agglom30$cluster)]

# Subset data and rename columns 
subset = clus %>% 
  select(cluster:agglomerateto_18) %>% 
  dplyr::rename(Communities31 = cluster,
                Communities18 = agglomerateto_18)

p = clustree(subset, prefix = "Communities", node_size=3, node_text_size = 2, node_text_angle = 60)
ggsave(plot = p, device = "png", width=20, height=10, dpi=300, path = paste0(path, "/Data/clustering_plots/"),
       filename = paste(date, "_clustree_21_40_49.png", sep = ""))

