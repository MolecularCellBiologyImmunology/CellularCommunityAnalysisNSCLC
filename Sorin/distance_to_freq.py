# Import packages
import pandas as pd
import os 
from scipy.spatial.distance import cdist
from datetime import datetime 

# Get current directory
current_directory = os.getcwd()

# Find datafile
path = f"{current_directory}/distance_calcs/20240131_allROIs_25px_neighbours_distanceCalculation_allCells.csv"
dists = pd.read_csv(path)

# Convert datafile of distances to relative frequency of source and target cell per source_id and per patient
# Loopless!
count_df = dists.groupby(['ROI_name','source_ID','source_cluster', 'target_cluster']).size().reset_index(name='count')
total_counts_scid = count_df.groupby(['ROI_name','source_ID','source_cluster'])['count'].sum().reset_index(name='total_scid')


count_df = pd.merge(count_df, total_counts_scid, on=['ROI_name','source_ID' ,'source_cluster'])
totalscid = count_df['total_scid']
# print(count_df.head())
count_df['relative_frequency'] = (count_df['count'] / count_df['total_scid'])

# Convert dataframe to correct format for clustering
sc = dists['source_cluster'].unique()
freq = count_df.pivot(columns=['target_cluster'], index=['ROI_name','source_ID','source_cluster'], values=['relative_frequency'])
freq.reset_index(inplace=True)
freq.columns  = freq.columns.droplevel(0)
freq.columns = ['ROI_name','source_ID','source_cluster',*freq.columns[3:]]
freq.fillna(0, inplace=True)

print('Dataset converted, below is a preview:\n',freq.head())

# Save file
output_path = f"{current_directory}/freq_table/"
os.makedirs(output_path, exist_ok=True)
totalscid.to_csv("total_scid.csv", index = False)

# freq.to_csv(f"{output_path}20240205_freqtable_sorin.csv", index = False)

