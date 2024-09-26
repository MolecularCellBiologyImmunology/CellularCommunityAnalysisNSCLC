
#### Script to run distance calculation using X- & Y-coordinates of each cell ####
## Optimised code for assigning_cellID function following distance calculatiom ##

## Megan Cole 
## November 2021

## Altered by Jannes Roelink
## March 2024

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
    saveto     : filename for saving distance information for all ROIs at once
"""

def cell_distances(cutoff, data, output_path):
    
    #Create an empty data frame for later use
    merged = pd.DataFrame()
    
    # get unique patients
    pats = set(data.Patient_ID)
    pats = list(pats)

    # Run distance calculation on each patient separately
    for pat in pats:
        print(f'Currently working on patient: {pat}')

        # take all cells of the current patient twice: source and target cell for distance calc 
        cells_A = data[data['Patient_ID'].str.match(pat)]

        cells_B = data[data['Patient_ID'].str.match(pat)]

        print("cells_A and cells_B created")
            
        # Call distance calculator function 
        distances  = distance_matrix(cutoff, cells_A[['Location_Center_X', 'Location_Center_Y']], cells_B[['Location_Center_X', 'Location_Center_Y']], pat)
        print('distances function complete')
        
        if distances.empty == True:
            print(f"No neighbours within {cutoff} pixels were identified for {pat}")
            continue
        print("Tested for presence of neighbours in distances dictionary")
        
        # Call assinging cellID function 
        neighbours = assigning_cellID(cells_A, cells_B, distances, pat, cutoff)
        print('assigning cell ID function complete')
        merged = pd.concat([merged, neighbours], ignore_index=True)

    # Reset the index 
    merged = merged.reset_index(drop = True)

    # Save the concatenated dataset
    merged.to_csv(f"{output_path}cell_distances_{cutoff}px_{date}_allpatients.csv", index = False)
    print("Function complete")


#### Distance calculation function ####
def distance_matrix(cutoff, points1, points2, pat):
    la, lb = len(points1), len(points2)
    tree1 = scipy.spatial.cKDTree(points1, leafsize=16)
    if points2 is None:
        points2 = points1
    tree2 = scipy.spatial.cKDTree(points2,leafsize=16)
    
    distances = tree1.sparse_distance_matrix(tree2, cutoff, output_type='dict')
        
    # CONVERT DISTANCES DIRECTORY INTO NEIGHBOURS DATAFRAME
    # Only carry on with next steps if distances dictionary has neighbour values 
    if bool(distances) == True:
        print("distances found for {pat}")
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

        if neighbours.iloc[:,0].nunique() != la:
            print('acells missing')
        # Rename column names 0 --> 'source' and 1 --> 'target'
        neighbours.rename(columns={neighbours.columns[0]: 'source', neighbours.columns[1]: 'target'}, inplace=True)
        # Subset for distances > 0 
        neighbours = neighbours[((neighbours['distance'] > 0))]
    else:
        pd.DataFrame(distances)
    
    return neighbours


#### Assigning information to cells identified as neighbours ####
def assigning_cellID(cells_A, cells_B, distances, pat, cutoff):
    
    # Link distances.source values to cellIDs in subset1
    cellsA_id = pd.DataFrame(cells_A.cellID.unique(), columns = ['source_cellID'])
    cellsA_id = cellsA_id[cellsA_id.index.isin(distances.source)]

    # Link distances.target values to cellIDs in subset2
    cellsB_id = pd.DataFrame(cells_B.cellID.unique(), columns = ['cellID'])
    cellsB_id = cellsB_id[cellsB_id.index.isin(distances.target)]

    # Add cellID info for source and target cells
    distances = distances.set_index(['source'])
    distj = distances.join(cellsA_id)
    distj = distj.set_index(['target'])
    distj = distj.join(cellsB_id)

    # Add marker info for source cells 
    distj = distj.set_index(['source_cellID'])
    cells_A = cells_A.set_index(['cellID'])
    distj = distj.join(cells_A)
    distj = distj.reset_index()

    # Rename new columns linking them to source cells
    distj = distj.rename(columns = {'source_cellID':'source_ID', "celltype":'source_cluster', 'Location_Center_X':'source_X', 'Location_Center_Y':'source_Y'})

    # Add marker info for target cells 
    distj = distj.set_index(['cellID', 'Patient_ID'])
    cells_B = cells_B.set_index(['cellID', 'Patient_ID'])
    distj = distj.join(cells_B)

    # Order dataframe by source_ID
    distj.sort_values('source_ID')
    distj = distj.reset_index()

    # Rename and re-order columns 
    distj = distj.rename(columns = {'cellID':'target_ID', "celltype":'target_cluster', 'Location_Center_X':'target_X', 'Location_Center_Y':'target_Y'})
    distj = distj[['source_ID', 'target_ID', 'distance', 'source_X', 'source_Y', 'target_X', 'target_Y', 'source_cluster',
            'target_cluster','Patient_ID']]

    print('Assigning cellID complete')
    return distj

######################################################################################################

# Set the working directory 
date = datetime.now().strftime("%Y%m%d")

path = os.getcwd()

output_path = f"{path}/Data/distances/"
input_path = f"{path}/Data/celldata/"

celldata = pd.read_csv(f"{input_path}celldata_{date}.csv")

# Only take cell columns of interest
celldata = celldata[['Patient_ID', 'cellID', 'celltype','Location_Center_X', 'Location_Center_Y']]

# run distance calculating function
# specify cutoff in pixels (1 pixel : 1 micron)
neighbours = cell_distances(25, celldata, output_path = output_path)

