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

clus = read.csv(paste0(path, "/Data/clustering/full_cd_CRICK_KO_clus_27_33_72.csv"))

# Cluster tree of neighbour clustering results 
#clus$temp = clus$agglomerateto_49[match(neighb1$cluster, clus$cluster)]

#neighb1$agglom30_average = agglom30$agglomerate_30[match(neighb1$cluster, agglom30$cluster)]

# Subset data and rename columns 
subset = clus %>% 
  select(cluster:agglomerateto_72) %>% 
  dplyr::rename(Communities190 = cluster,
                Communities72 = agglomerateto_72,
                Communities33 = agglomerateto_33,
                Communities27 = agglomerateto_27)

p = clustree(subset, prefix = "Communities", node_size=3, node_text_size = 2, node_text_angle = 60)
ggsave(plot = p, device = "png", width=20, height=10, dpi=300, path = paste0(path, "/Data/plots/MOC2_CCR2KO/"),
       filename = paste(date, "_clustree_CRICK_KO_clusters_n3_w30.png", sep = ""))

