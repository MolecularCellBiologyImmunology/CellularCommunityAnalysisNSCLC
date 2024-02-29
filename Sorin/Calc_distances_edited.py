
#### Script to run distance calculation using X- & Y-coordinates of each cell ####
## Optimised code for assinging_cellID function following distance calculatiom ##

## Megan Cole
## November 2021

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
    lst1       : 'source' cell type(s) of interest (cluster/metacluster) 
    lst2       : 'target' neighbour cell type(s) (cluster/metacluster)
    saveto     : filename for saving distance information for all ROIs at once
    *NOTE      :  Always use [] for lst1/lst2. If wanting to select all cell types keep [] empty* 
"""

def cell_distances(cutoff, data, allData, domain, clustering, lst1, lst2, save_ind, save_all, output_path):
    
    #Create an empty data frame for later use
    merged = pd.DataFrame()
    #If allData = 0, take domain info. If allData = 1, skip taking domain info
    if allData == "0":
        # Remove any rows that have an NA in the domain column
        data = data[data['domain'].notna()]
        # Take data just from the selected domain
        data = data[data['domain'].str.match(domain)]
        
    ROIs = set(data.ROI_name)
    ROIs = list(ROIs)
    # Run distance calculation on each ROI separately
    for ROI in ROIs:
        print(ROI)
        # Subset data for cell types of interest 
        if not lst1:
            cells_A = data[data['ROI_name'].str.match(ROI)]
        else: 
            cells_A = data[data['ROI_name'].str.match(ROI) & data[clustering].isin(lst1)]
        if not lst2:
            cells_B = data[data['ROI_name'].str.match(ROI)]
        else: 
            cells_B = data[data['ROI_name'].str.match(ROI) & data[clustering].isin(lst2)]
        print("cells_A and cells_B created")
            
        # Call distance calculator function 
        distances  = distance_matrix(cutoff, cells_A[['Location_Center_X', 'Location_Center_Y']], cells_B[['Location_Center_X', 'Location_Center_Y']])
        print('distances function complete')
        
        if distances.empty == True:
            print(f"No neighbours within {cutoff} pixels were identified for {ROI}")
            continue
        print("Tested for presence of neighbours in distances dictionary")
        
        # Call assinging cellID function 
        neighbours = assigning_cellID(cells_A, cells_B, distances, ROI, cutoff, clustering, save_ind)
        print('assigning cell ID function complete')
        # print(neighbours, 'uhigrh')
        merged = pd.concat([merged, neighbours], ignore_index=True)
        # merged = merged.append(neighbours)

    # Reset the index 
    merged = merged.reset_index(drop = True)

    # Save the concatenated dataset
    merged = merged.drop(columns=['treatment'])
    merged.to_csv(f"{output_path}20240217_allROIs_{cutoff}px_{save_all}", index = False)
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
    if bool(distances) == True:
        print("distances found for {ROI}")
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
def assigning_cellID(cells_A, cells_B, distances, ROI, cutoff, clustering, save_ind):
    
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
    # changed this line cause think its a typo
    #     distj = distj.rename(columns = {'index':'source_ID', f"{clustering}":'source_cluster', 'Location_Center_X':'source_X', 'Location_Center_Y':'source_Y'})

    distj = distj.rename(columns = {'source_cellID':'source_ID', f"{clustering}":'source_cluster', 'Location_Center_X':'source_X', 'Location_Center_Y':'source_Y'})

    # Add marker info for target cells 
    distj = distj.set_index(['cellID', 'treatment', 'ROI_name'])
    cells_B = cells_B.set_index(['cellID', 'treatment', 'ROI_name'])
    distj = distj.join(cells_B)
    # Order dataframe by source_ID
    distj.sort_values('source_ID')
    distj = distj.reset_index()
    # Rename and re-order columns 
    distj = distj.rename(columns = {'cellID':'target_ID', f"{clustering}":'target_cluster', 'Location_Center_X':'target_X', 'Location_Center_Y':'target_Y'})
    distj = distj[['source_ID', 'target_ID', 'distance', 'source_X', 'source_Y', 'target_X', 'target_Y', 'source_cluster',
            'target_cluster','ROI_name', 'treatment']]

    # Return/save neighbours 
    distj_notreat = distj.drop(columns=['treatment'])

    distj_notreat.to_csv(f"{output_path}20240217_{ROI}_{cutoff}px_{save_ind}", index = False)
    print('Assinging cellID complete')
    return distj

######################################################################################################

# Set the working directory 
# path = "/Users/vanmalf/Documents/Amsterdam/Collaborations/Iris Miedema/Results/outputs_Febe/neighbours/"
current_directory = os.getcwd()

output_path = f"{current_directory}/distance_calcs/"
os.makedirs(output_path, exist_ok=True)

# Load the data 
input_path = f"{current_directory}/celldata/"

celldata = pd.read_csv(f"{input_path}20240208_celldata.csv")

# Only take cell columns of interest
# list(celldata)
celldata = celldata[['ROI_name', 'cellID', 'study', 'class','Location_Center_X', 'Location_Center_Y']]
celldata = celldata.rename(columns={ "study":"treatment"})

# Please note: Do you want to set a threshold? This will filter out all other distances!
# Only use threshold if this table will be reused to build a neighbours table
# If the distance as measurement is used for generating plots, then set a very large threshold.
# Call cell distances function to get dataset with distances calculated 
neighbours = cell_distances(35, celldata, allData = "1", domain = "", clustering = "class", lst1 = [], lst2 = [],
save_ind = "neighbours_distance_calculation.csv", save_all = "neighbours_distanceCalculation_allCells.csv", output_path = output_path)

