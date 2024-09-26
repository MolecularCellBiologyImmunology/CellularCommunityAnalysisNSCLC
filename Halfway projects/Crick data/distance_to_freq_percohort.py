# Import packages
import pandas as pd
import os 
from scipy.spatial.distance import cdist
from datetime import datetime 

# Get current directory
cwd = os.getcwd()
input_path = f"{cwd}/Data/distances/"
output_path = f"{cwd}/Data/frequencies/"

groups = ['MOC1_MOCAF', 'MOC1_WT', 'MOC2_CCR2KO', 'MOC2_WT']
date = datetime.now().strftime("%Y%m%d")
cutoff = 25

for group in groups:
    # Find datafile
    dists = pd.read_csv(f"{input_path}dists_cellpairs_{group}_{cutoff}px_{date}.csv")

    # Convert datafile of distances to relative frequency of source and target cell per source_id and per patient
    count_df = dists.groupby(['Unique_ROI_ID','source_ID','source_cluster', 'target_cluster']).size().reset_index(name='count')
    total_counts_scid = count_df.groupby(['Unique_ROI_ID','source_ID','source_cluster'])['count'].sum().reset_index(name='total_scid')

    count_df = pd.merge(count_df, total_counts_scid, on=['Unique_ROI_ID','source_ID' ,'source_cluster'])
    totalscid = count_df['total_scid']
    count_df['relative_frequency'] = (count_df['count'] / count_df['total_scid'])

    # Convert dataframe to correct format for clustering
    sc = dists['source_cluster'].unique()
    freq = count_df.pivot(columns=['target_cluster'], index=['Unique_ROI_ID','source_ID','source_cluster'], values=['relative_frequency'])
    freq.reset_index(inplace=True)
    freq.columns  = freq.columns.droplevel(0)
    freq.columns = ['Unique_ROI_ID','source_ID','source_cluster',*freq.columns[3:]]
    freq.fillna(0, inplace=True)

    print('Dataset converted, below is a preview:\n',freq.head())

    # Save file
    freq.to_csv(f"{output_path}freqs_{group}_{cutoff}px_{date}.csv", index = False)


