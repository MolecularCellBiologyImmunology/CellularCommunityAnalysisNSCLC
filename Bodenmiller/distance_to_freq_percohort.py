# Import packages
import pandas as pd
import os 
from scipy.spatial.distance import cdist
from datetime import datetime 

# Get current directory
cwd = os.getcwd()
input_path = f"{cwd}/Data/distances/"
output_path = f"{cwd}/Data/frequencies/"

cohorts = ['86_A', '86_B', '86_C', '87_A', '87_B', '87_C', '88_A', '88_B', '88_C', '175_A', '175_B', '175_C', '176_A', '176_B', '176_C', '178_A', '178_B', '178_C']
date = datetime.now().strftime("%Y%m%d")
cutoff = 25

for cohort in cohorts:
    # Find datafile
    dists = pd.read_csv(f"{input_path}dists_cellpairs_{cohort}_{cutoff}px_{date}.csv")

    # Convert datafile of distances to relative frequency of source and target cell per source_id and per patient
    count_df = dists.groupby(['acID', 'Patient_ID','source_ID','source_cluster', 'target_cluster']).size().reset_index(name='count')
    total_counts_scid = count_df.groupby(['acID', 'Patient_ID','source_ID','source_cluster'])['count'].sum().reset_index(name='total_scid')

    count_df = pd.merge(count_df, total_counts_scid, on=['acID', 'Patient_ID','source_ID' ,'source_cluster'])
    totalscid = count_df['total_scid']
    count_df['relative_frequency'] = (count_df['count'] / count_df['total_scid'])

    # Convert dataframe to correct format for clustering
    sc = dists['source_cluster'].unique()
    freq = count_df.pivot(columns=['target_cluster'], index=['acID', 'Patient_ID','source_ID','source_cluster'], values=['relative_frequency'])
    freq.reset_index(inplace=True)
    freq.columns  = freq.columns.droplevel(0)
    freq.columns = ['acID', 'Patient_ID','source_ID','source_cluster',*freq.columns[4:]]
    freq.fillna(0, inplace=True)

    print('Dataset converted, below is a preview:\n',freq.head())

    # Save file
    freq.to_csv(f"{output_path}freqs_{cohort}_{cutoff}px_{date}.csv", index = False)


