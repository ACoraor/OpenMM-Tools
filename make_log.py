import os
import mdtraj as md
# from subtool import *
import pandas as pd
import numpy as np
# from clean_thermo import clean_all


def main():
    """Create a logfile with basic data from state_data. Remove data clashes
    And properly filter array.
    """
    
    #Read state_data data
    state_data = read_state('state_data.txt',recursive=False)
    #Use recursive if state_datas are buried in subdirectories    

    print("Reading hdf5 files.")
    #Read full_output
    if os.path.isfile('full_output.h5'):
        full_out = md.load_hdf5('full_output.h5')
    else:
        full_out = None

    #Read output.h5
    if os.path.isfile('output.h5'):
        out = md.load_hdf5('output.h5')
    else:
        out = full_out
        full_out = None #Swap full_out to out

    #Create new logfile
    print("Writing logfile.")
    write_logfile(state_data,full_out,out)
    print("Success. Exiting...")
    
def read_state(fn='state_data.txt',recursive=False):
    """Read all state data. If there exists state_data in subdirectories,
    load those first. Delete invalid entries.
    
    Parameters:
        fn: *str*
            Path to state data file.
    Returns:
        data: *np.array*, shape: (n_timesteps)
    """
    statefiles = []
    if recursive:
        dirs = [ele for ele in os.listdir('.') if os.path.isdir(ele)]
        #Check subdirs for state_data. Load them.
        for d in dirs:
            joined_name = os.path.join(d,fn)
            if os.path.isfile(joined_name):
                statefiles.append(pd.read_csv(joined_name))
    if not recursive:
        statefiles.append(pd.read_csv(fn))

    #Join all statefiles
    if len(statefiles) > 1:
        data = pd.concat(statefiles)
    else:
        data = statefiles[0]
    #print("pre-drop length:",len(data))
    #Remove non-numerical data
    #Get indices
    data.index = np.arange(len(data.index))
    #print('index:',data.index)
    drop_inds = np.argwhere(pd.to_numeric(data['#"Step"'],errors='coerce').isnull()).ravel()
    data = data.drop(index=drop_inds)
    #print("Final length:",len(data2))

    #Drop duplicates of thermo indices
    #data = data.drop_duplicates(subset='#"Step"')
    return data
    


def write_logfile(state_data,full_out,out):
    """Create a logfile with an explicit row for every valid timestep.
    Fill missing data with NULL. 
    
    Parameters:
        fn: *str*
            Path to state data file.
    Returns:
        data: *np.array*, shape: (n_timesteps)
    """
    pass
    #If number of state_data frames and output frames don't match, crash
    if full_out is None:
        out_len = len(out.xyz)
    else:
        out_len = len(full_out.xyz) + len(out.xyz)
    if out_len != len(state_data):
        raise ValueError("Output frames and state data do not match in length:%s vs %s" % (
            out_len,len(state_data)))

    renamer = {'#"Step"':'timestep', "Potential Energy (kJ/mole)":'pe',
                "Temperature (K)":"temp"}
    data = state_data.rename(columns=renamer)
    data.to_csv('log.dframe',index=False,na_rep='NULL')
 

if __name__ == "__main__":
    main()
