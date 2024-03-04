
#### Script to run distance calculation using X- & Y-coordinates of each cell ####
## Optimised code for assigning_cellID function following distance calculation ##

## Heavily based on Megan Cole's script

## Jannes Roelink 29/2/2024

# Import packages
import pandas as pd
import os 
import scipy
from scipy.spatial.distance import cdist
from scipy.spatial import KDTree
from datetime import datetime 

#########################
##### SET FUNCTIONS #####
#########################

""" 
    #### COMPUTING DISTANCES BETWEEN TWO GROUPS OF CELLS - Megan Cole ####
    
    ** Input **
    cutoff     : distance cutoff between centre of cells (pixels)
    data       : input data frame 
    allData    : '0' --> not using whole tissue info, '1' --> using whole tissue info
    domain     : tissue domain of interest ('' if using whole tissue info)
    clustering : column for cell type information
    saveto     : filename for saving distance information for all ROIs at once
    *NOTE      :  Always use [] for lst1/lst2. If wanting to select all cell types keep [] empty* 
"""

def cell_distances(cutoff, data, date, group, output_path):
    
    #Create an empty data frame for later use
    merged = pd.DataFrame()
    #If allData = 0, take domain info. If allData = 1, skip taking domain info
        
    ROIs = set(data.Unique_ROI_ID)
    ROIs = list(ROIs)
    # Run distance calculation on each acID separately

    for ROI in ROIs:
        print(ROI)
        cells_A = data[data['Unique_ROI_ID'] == ROI]
        cells_B = data[data['Unique_ROI_ID'] == ROI]
        # Patient_ID = cells_A['Patient_ID'].unique()
        print("cells_A and cells_B created")
            
        # Call distance calculator function 
        distances  = distance_matrix(cutoff, cells_A[['Location_Center_X', 'Location_Center_Y']], cells_B[['Location_Center_X', 'Location_Center_Y']])
        print('distances function complete')
        if distances.empty == True:
            print(f"No neighbours within {cutoff} pixels were identified for {ROI}")
            continue
        print("Tested for presence of neighbours in distances dictionary")
        
        # Call assinging cellID function 
        neighbours = assigning_cellID(cells_A, cells_B, distances, ROI, cutoff)
        print('assigning cell ID function complete')
        merged = pd.concat([merged, neighbours], ignore_index=True)
        # merged = merged.append(neighbours)

    # Reset the index 
    merged = merged.reset_index(drop = True)

    # Save the concatenated dataset
    merged.to_csv(f"{output_path}dists_cellpairs_{group}_{cutoff}px_{date}.csv", index = False)
    print("Function complete")


#### Distance calculation function ####
def distance_matrix(cutoff, points1, points2):
    
    tree1 = scipy.spatial.cKDTree(points1, leafsize=16)
    if points2 is None:
        points2 = points1
    tree2 = scipy.spatial.cKDTree(points2,leafsize=16)
    
    distances = tree1.sparse_distance_matrix(tree2, cutoff, output_type='dict')
        
    # CONVERT DISTANCES DIRECTORY INTO NEIGHBOURS DATAFRAME
    # Only carry on with next steps if distances dictionary has neighbour values 
    if bool(distances):
        print(f"distances found for ROI")
        keys = pd.DataFrame.from_dict(distances.keys())
        values = pd.DataFrame.from_dict(distances.values())
        # Give name to values dataframe 
        values.columns = ['distance']
        # Concatenate keys and values dataframes 
        neighbours = pd.concat([keys, values], axis = 1)
        # Sort data frame based on ascending order of first column 
        neighbours.sort_values([0,1], inplace = True) #Check this is sorting and not new values
        # Reset the index 
        neighbours = neighbours.reset_index(drop = True)
        # Rename column names 0 --> 'source' and 1 --> 'target'
        neighbours.rename(columns={neighbours.columns[0]: 'source', neighbours.columns[1]: 'target'}, inplace=True)
        # Subset for distances > 0 
        neighbours = neighbours[((neighbours['distance'] > 0))]
    else:
        pd.DataFrame(distances)
    
    return neighbours


#### Assigning information to cells identified as neighbours ####
def assigning_cellID(cells_A, cells_B, distances, ROI, cutoff):
    
    # Link distances.source values to cellIDs in subset1
    cellsA_id = pd.DataFrame(cells_A.CellID.unique(), columns = ['source_cellID'])
    cellsA_id = cellsA_id[cellsA_id.index.isin(distances.source)]

    # Link distances.target values to cellIDs in subset2
    cellsB_id = pd.DataFrame(cells_B.CellID.unique(), columns = ['CellID'])
    cellsB_id = cellsB_id[cellsB_id.index.isin(distances.target)]

    # Add cellID info for source and target cells
    distances = distances.set_index(['source'])
    distj = distances.join(cellsA_id)
    distj = distj.set_index(['target'])
    distj = distj.join(cellsB_id)

    # Add marker info for source cells 
    distj = distj.set_index(['source_cellID'])
    cells_A = cells_A.set_index(['CellID'])
    distj = distj.join(cells_A)

    distj = distj.reset_index()

    # Rename new columns linking them to source cells
    distj = distj.rename(columns = {'source_cellID':'source_ID', 'cluster_num':'source_cluster', 'Location_Center_X':'source_X', 'Location_Center_Y':'source_Y'})

    # Add marker info for target cells 
    distj = distj.set_index(['CellID'])
    cells_B = cells_B.set_index(['CellID'])
    cells_B = cells_B.rename(columns = {'source_cluster': 'target_cluster'})
    cells_B.drop(columns=['Unique_ROI_ID'], inplace=True)
    distj = distj.join(cells_B)

    # Order dataframe by source_ID
    distj.sort_values('source_ID')
    distj = distj.reset_index()

    # Rename and re-order columns 
    distj = distj.rename(columns = {'CellID':'target_ID', 'cluster_num':'target_cluster', 'Location_Center_X':'target_X', 'Location_Center_Y':'target_Y'})
    distj = distj[['source_ID', 'target_ID', 'distance', 'source_X', 'source_Y', 'target_X', 'target_Y', 'source_cluster',
            'target_cluster','Unique_ROI_ID']]


    # distj.to_csv(f"{output_path}20240228_{acID}_{cutoff}px_{save_ind}", index = False)
    print(f'{ROI}, Assigning CellID complete')
    return distj

######################################################################################################

# Set the working directory 
# path = "/Users/vanmalf/Documents/Amsterdam/Collaborations/Iris Miedema/Results/outputs_Febe/neighbours/"
cwd = os.getcwd()
input_path = f"{cwd}/Data/celldata/"

output_path = f"{cwd}/Data/distances/"

groups = ['MOC1_MOCAF', 'MOC1_WT', 'MOC2_CCR2KO', 'MOC2_WT']

today = datetime.now()
date = today.strftime("%Y%m%d")

for group in groups:
    celldata = pd.read_csv(f"{input_path}{group}_celldata_{date}.csv")
    # Only take cell columns of interest
    celldata = celldata[['Location_Center_X', 'Location_Center_Y', 'cluster_num', 'Unique_ROI_ID', 'CellID']]
    # Set maximum pixel/micron distance for cells to be considered neighbours
    # Note Distance is calculated cell center to cell center, so take radius into account
    max_distance = 25
    neighbours = cell_distances(25, celldata, date, group, output_path = output_path)

