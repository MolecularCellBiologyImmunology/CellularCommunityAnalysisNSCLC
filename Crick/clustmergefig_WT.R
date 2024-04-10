library(ggplot2) 
library(dplyr)
library(Rtsne) 
library(clustree) 
library(ggdendro) 
library(tiff) 
library(reshape2) 

date = format(Sys.Date(), "%Y%m%d")

path = getwd()

#dir.create(out_dir)

clus = read.csv(paste0(path, "/Data/clustering/full_cd_CRICK_WT_clus_29_40_59.csv"))

# Cluster tree of neighbour clustering results 
#clus$temp = clus$agglomerateto_49[match(neighb1$cluster, clus$cluster)]

#neighb1$agglom30_average = agglom30$agglomerate_30[match(neighb1$cluster, agglom30$cluster)]

# Subset data and rename columns 
subset = clus %>% 
  select(cluster:agglomerateto_59) %>% 
  dplyr::rename(Communities170 = cluster,
                Communities59 = agglomerateto_59,
                Communities40 = agglomerateto_40,
                Communities29 = agglomerateto_29)

p = clustree(subset, prefix = "Communities", node_size=3, node_text_size = 2, node_text_angle = 60)
ggsave(plot = p, device = "png", width=20, height=10, dpi=300, path = paste0(path, "/Data/plots/MOC2_WT/"),
       filename = paste(date, "_clustree_CRICK_WT_clusters_n3_w30.png", sep = ""))

