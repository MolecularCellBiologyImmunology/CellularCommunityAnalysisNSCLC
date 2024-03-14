
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

def cell_distances(cutoff, data, output_path):
    
    #Create an empty data frame for later use
    merged = pd.DataFrame()
        
    pats = set(data.Patient_ID)
    pats = list(pats)
    # Run distance calculation on each ROI separately
    for pat in pats:
        print(pat)
        # Subset data for cell types of interest 

        cells_A = data[data['Patient_ID'].str.match(pat)]

        cells_B = data[data['Patient_ID'].str.match(pat)]
        if  cells_A.isna().sum().sum():

            print(cells_A.isna().sum())
            print('bye, cellsa, ya')
        print("cells_A and cells_B created")

        # Call distance calculator function 
        distances  = distance_matrix(cutoff, cells_A[['Location_Center_X', 'Location_Center_Y']], cells_B[['Location_Center_X', 'Location_Center_Y']], pat)
        print('distances function complete')
        # print(distances.head())

        if distances.empty == True:
            print(f"No neighbours within {cutoff} pixels were identified for {pat}")
            continue
        print("Tested for presence of neighbours in distances dictionary")
        
        # Call assigning cellID function 
        neighbours = assigning_cellID(cells_A, cells_B, distances, pat, cutoff)
        print('assigning cell ID function complete')
        # print(neighbours, 'uhigrh')
        merged = pd.concat([merged, neighbours], ignore_index=True)
        # merged = merged.append(neighbours)

    # Reset the index 
    merged = merged.reset_index(drop = True)

    # Save the concatenated dataset
    merged.to_csv(f"{output_path}dists_cellpairs_{cutoff}px_{date}_subset.csv", index = False)
    print("Function complete")


#### Distance calculation function ####
def distance_matrix(cutoff, points1, points2, pat):
    # print(points1.head(), 'points1')
    tree1 = scipy.spatial.cKDTree(points1, leafsize=16)
    if points2 is None:
        points2 = points1
    tree2 = scipy.spatial.cKDTree(points2,leafsize=16)
    
    distances = tree1.sparse_distance_matrix(tree2, cutoff, output_type='dict')
    # CONVERT DISTANCES DIRECTORY INTO NEIGHBOURS DATAFRAME
    # Only carry on with next steps if distances dictionary has neighbour values 
    if bool(distances) == True:
        print(f"distances found for {pat}")
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
        # print(neighbours.head())
        # Rename column names 0 --> 'source' and 1 --> 'target'
        neighbours.rename(columns={neighbours.columns[0]: 'source_num', neighbours.columns[1]: 'target_num'}, inplace=True)
        # Subset for distances > 0 
        neighbours = neighbours[((neighbours['distance'] > 0))]
    else:
        pd.DataFrame(distances)
    
    return neighbours


#### Assigning information to cells identified as neighbours ####
def assigning_cellID(cells_A, cells_B, distances, pat, cutoff):
    
    # Link distances.source values to cellIDs in subset1
    cellsA_id = pd.DataFrame(cells_A.cellID.unique(), columns = ['source_ID'])
    print('dists\n', distances)
    cellsA_id = cellsA_id[cellsA_id.index.isin(distances.source_num)]
    # print(distj[distj.isnull().any(axis=1)])
    # print('cela \n',cellsA_id.index.max())

    # if  cellsA_id.isna().sum().sum():

    #     print(cellsA_id.isna().sum())
    #     print('bye, cella, ya')
    # Link distances.target values to cellIDs in subset2
    cellsB_id = pd.DataFrame(cells_B.cellID.unique(), columns = ['target_ID'])
    cellsB_id = cellsB_id[cellsB_id.index.isin(distances.target_num)]
    # Add cellID info for source and target cells
    distances = distances.set_index(['source_num'])
    distj = distances.join(cellsA_id)
    # print(distj[distj.isnull().any(axis=1)])

    # print(distances.head())

    distj = distj.set_index(['target_num'])
    distj = distj.join(cellsB_id)

    # Add marker info for source cells 
    distj = distj.set_index(['source_ID'])
    cells_A = cells_A.set_index(['cellID'])
    if  distj.isna().sum().sum():

        print(distj.isna().sum())
        print('hey')


    # print(cells_A.head())
    # print(distj.head())

    # if  cells_A.isna().sum().sum():

    #     print(cells_A.isna().sum())
    #     print('bye, cella')
    distj = distj.join(cells_A)

    # print(distj[distj.isnull().any(axis=1)])
    # print('after', distj.head())
    # if  distj.isna().sum().sum():

    #     print(distj.isna().sum())
    #     print('bye')
    distj = distj.reset_index()

    # Rename new columns linking them to source cells
    # changed this line cause think its a typo
    #     distj = distj.rename(columns = {'index':'source_ID', f"{clustering}":'source_cluster', 'Location_Center_X':'source_X', 'Location_Center_Y':'source_Y'})

    distj = distj.rename(columns = {'cellID':'source_ID', "target_ID":"cellID", "class":'source_cluster', 'Location_Center_X':'source_X', 'Location_Center_Y':'source_Y'})


    # Add marker info for target cells 
    distj = distj.set_index(['cellID', 'Patient_ID'])
    cells_B = cells_B.set_index(['cellID', 'Patient_ID'])
    print('check1')
    # print(cells_B.head())

    distj = distj.join(cells_B)
    # Order dataframe by source_ID
    # if  distj.isna().sum().sum():

    #     print(distj.isna().sum())
    print('check2')

    distj.sort_values('source_ID')
    distj = distj.reset_index()
    # print(distj.head())

    # Rename and re-order columns 
    distj = distj.rename(columns = {'cellID':'target_ID', "class":'target_cluster', 'Location_Center_X':'target_X', 'Location_Center_Y':'target_Y'})
    distj = distj[['source_ID', 'target_ID', 'distance', 'source_X', 'source_Y', 'target_X', 'target_Y', 'source_cluster',
            'target_cluster','Patient_ID']]

    # Return/save neighbours 

    # distj_notreat.to_csv(f"{output_path}20240217_{ROI}_{cutoff}px_{save_ind}", index = False)
    print('Assigning cellID complete')
    return distj

######################################################################################################

# Set the working directory 
# path = "/Users/vanmalf/Documents/Amsterdam/Collaborations/Iris Miedema/Results/outputs_Febe/neighbours/"
cwd = os.getcwd()

date = datetime.now().strftime("%Y%m%d")
cutoff = 25
output_path = f"{cwd}/Data/distances/"

# Load the data 
input_path = f"{cwd}/Data/celldata/"

celldata = pd.read_csv(f"{input_path}celldata_{date}.csv")
celldata = celldata.loc[celldata['Patient_ID'] == 'LUAD_D001']
print('data nans \n',celldata.isna().sum())
print(celldata['class'].unique())
# Only take cell columns of interest
# list(celldata)
celldata = celldata[['Patient_ID', 'cellID', 'class','Location_Center_X', 'Location_Center_Y']]

# Please note: Do you want to set a threshold? This will filter out all other distances!
# Only use threshold if this table will be reused to build a neighbours table
# If the distance as measurement is used for generating plots, then set a very large threshold.
# Call cell distances function to get dataset with distances calculated 
neighbours = cell_distances(cutoff, celldata, output_path = output_path)

