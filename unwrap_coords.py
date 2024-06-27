#Module conv_h5
#Created by Aria Coraor
#Written 3/28/24

import numpy as np
import mdtraj as md
import argparse


def main():
    """Unwrap the trajectory."""
    a = md.load_hdf5(args.file)
    xyz = a.xyz

    chains = [ch for ch in a.topology.chains]
    atoms_1 = np.array([at.index for at in chains[0].atoms])
    atoms_2 = np.array([at.index for at in chains[1].atoms])

    #Calculate COMS for each chain
    sense_com = np.average(xyz[:,atoms_1],axis=1) #Shape: n_ts, 3
    antisense_com = np.average(xyz[:,atoms_2],axis=1)


    #find the frame with the median distance between COMS
    sa_dist = antisense_com - sense_com

    dist_lens = np.linalg.norm(sa_dist,axis=1)
    #Get index of median
    approx_median = np.argsort(dist_lens)[len(dist_lens)//2]
    
    #Create a rolling corrector, nullifying any hard delta on any frame.
    correctors = np.zeros((2,sense_com.shape[0],3)) #First row is sense/anti,
        #Then timestep, then xyz displacement
    sense_diff = sense_com[:-1] - sense_com[1:]
    antisense_diff = antisense_com[:-1] - antisense_com[:-1]
    
    bounds = (20.0,5.0,5.0)

    #i = approx_median
    #Move forward from the median and check for deltas
    for i in range(approx_median,len(sense_diff)):
        #Check if moved off into each direction
        for j in range(len(bounds)):
            if sense_diff[i][j] > 0.9*bounds[j]:
                correctors[0,i:,j] -= bounds[j]
            elif sense_diff[i][j] < -0.9*bounds[j]:
                correctors[0,i:,j] += bounds[j]
        #Now check antisense strand
        for j in range(len(bounds)):
            if antisense_diff[i][j] > 0.9*bounds[j]:
                correctors[1,i:,j] -= bounds[j]
            elif antisense_diff[i][j] < -0.9*bounds[j]:
                correctors[1,i:,j] += bounds[j]

    #Move backward from the median and check for deltas
    for i in reversed([x for x in range(approx_median)]):
        #Check if moved off into each direction
        for j in range(len(bounds)):
            if sense_diff[i][j] > 0.9*bounds[j]:
                correctors[0,:i+1,j] -= bounds[j]
            elif sense_diff[i][j] < -0.9*bounds[j]:
                correctors[0,:i+1,j] += bounds[j]
        #Now check antisense strand
        for j in range(len(bounds)):
            if antisense_diff[i][j] > 0.9*bounds[j]:
                correctors[1,:i+1,j] -= bounds[j]
            elif antisense_diff[i][j] < -0.9*bounds[j]:
                correctors[1,:i+1,j] += bounds[j]
    
    print("Correctors constructed.")
    #Construct new trajectory with corrections in place
    uw_xyz = np.copy(xyz)
    uw_xyz[:,atoms_1] += correctors[0]
    uw_xyz[:,atoms_2] += correctors[1]
    
    print("Corrected xyz created. Validating...")
    sense_com = np.average(uw_xyz[:,atoms_1],axis=1) #Shape: n_ts, 3
    antisense_com = np.average(uw_xyz[:,atoms_2],axis=1)
    sa_dist = antisense_com - sense_com
    dist_lens = np.linalg.norm(sa_dist,axis=1)
    print("Max distances:",np.max(dist_lens,axis=0))
    print("Min distances:",np.min(dist_lens,axis=0))
    
    traj = md.Trajectory(xyz = uw_xyz, topology = a.topology, time = a.time,
            unitcell_lengths = bounds)
    print("Saving trajectory to uw_output.h5")
    traj.save_hdf5('uw_output.h5')

    
    print("Saved converted traj.")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f','--file', type=str, default='output.h5',
        help='Path to hdf5 output')
    args = parser.parse_args()

    main()

