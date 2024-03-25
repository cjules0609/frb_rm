import h5py
import numpy as np

def read(path):
    data = []
    with h5py.File(path, "r") as f:
        # Print all root level object names (aka keys) 
        # these can be group or dataset names 
        # print("Keys: %s" % f.keys())
        
        meshblock_size = f.attrs['MeshBlockSize']
        num_meshblocks = f.attrs['NumMeshBlocks']
        meshblock_loc = f['LogicalLocations'][()]

        # get first object name/key; may or may NOT be a group
        # a_group_key = list(f.keys())[0] # 'B'

        # get the object type for a_group_key: usually group or dataset
        # print(type(f[a_group_key])) 

        # If a_group_key is a group name, 
        # this gets the object names in the group and returns as a list
        B = f['B'][()]

        # B in x-direction
        B_x = np.zeros([
            int(np.sqrt(len(B[0]))*len(B[0][0][0])),        # x-coord
            int(np.sqrt(len(B[0]))*len(B[0][0][0][0]))      # y-coord
            ])
        # B in y-direction
        B_y = np.zeros([
            int(np.sqrt(len(B[1]))*len(B[0][0][0])),        # x-coord
            int(np.sqrt(len(B[1]))*len(B[0][0][0][0]))      # y-coord
            ])
        # B in z-direction
        B_z = np.zeros([
            int(np.sqrt(len(B[2]))*len(B[0][0][0])),        # x-coord
            int(np.sqrt(len(B[2]))*len(B[0][0][0][0]))      # y-coord
            ])

        for i,meshblock in enumerate(B[0]):
            B_x[int((meshblock_loc[i][1]*meshblock_size[0])):int((meshblock_loc[i][1]+1)*meshblock_size[0]), 
                int((meshblock_loc[i][0]*meshblock_size[0])):int((meshblock_loc[i][0]+1)*meshblock_size[0]),
                ] = meshblock
            
        for i,meshblock in enumerate(B[1]):
            B_y[int((meshblock_loc[i][1]*meshblock_size[0])):int((meshblock_loc[i][1]+1)*meshblock_size[0]), 
                int((meshblock_loc[i][0]*meshblock_size[0])):int((meshblock_loc[i][0]+1)*meshblock_size[0]),
                ] = meshblock
            
        for i,meshblock in enumerate(B[2]):
            B_z[int((meshblock_loc[i][1]*meshblock_size[0])):int((meshblock_loc[i][1]+1)*meshblock_size[0]), 
                int((meshblock_loc[i][0]*meshblock_size[0])):int((meshblock_loc[i][0]+1)*meshblock_size[0]),
                ] = meshblock
    return B_z



    # If a_group_key is a dataset name, 
    # this gets the dataset values and returns as a list
    # data = list(f[a_group_key])
    # preferred methods to get dataset values:
    # ds_obj = f[a_group_key]      # returns as a h5py dataset object
    # ds_arr = f[a_group_key][()]  # returns as a numpy array
