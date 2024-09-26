#Module conv_h5
#Created by Aria Coraor
#Written 3/28/24

import numpy as np
import mdtraj as md
import argparse
#import zfit
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
#import sklearn.cluster


def main():
    """Unwrap the trajectory."""

    #Load the file
    a = md.load_hdf5(args.file)
    xyz = a.xyz

    chains = [ch for ch in a.topology.chains]
    atoms_1 = np.array([at.index for at in chains[0].atoms])
    atoms_2 = np.array([at.index for at in chains[1].atoms])

    #Calculate COMS for each chain
    sense_com = np.average(xyz[:,atoms_1],axis=1) #Shape: n_ts, 3
    antisense_com = np.average(xyz[:,atoms_2],axis=1)

    print("Max sense_com:",np.max(sense_com,axis=0))
    print("Min sense_com:",np.min(sense_com,axis=0))
    print("Max antisense_com:",np.max(antisense_com,axis=0))
    print("Min antisense_com:",np.min(antisense_com,axis=0))
    
    #Wrap into usable box:
    #sense_x_corrections = np.argwhere(sense_com[:,0] > 27.0)
    
    

    #quit()
    #find the frame with the median distance between COMS
    sa_dist = antisense_com - sense_com

    dist_lens = np.linalg.norm(sa_dist,axis=1)
    max_diff = np.argmax(dist_lens)
    break_start = max_diff - 10
    break_end = max_diff + 10
    print("distances, %s:%s:" % (break_start,break_end),dist_lens[break_start:break_end])
    print("sense_coms, %s:%s:" % (break_start,break_end),sense_com[break_start:break_end])
    print("antisense_coms, %s:%s:" % (break_start,break_end),antisense_com[break_start:break_end])
    #Get index of median
    
    #approx_median = np.argsort(dist_lens)[len(dist_lens)//2]
    approx_median = np.argmin(dist_lens)
    
    #Create a rolling corrector, nullifying any hard delta on any frame.
    correctors = np.zeros((2,sense_com.shape[0],3)) #First row is sense/anti,
        #Then timestep, then xyz displacement
    sense_diff = sense_com[1:] - sense_com[:-1]
    antisense_diff = antisense_com[1:] - antisense_com[:-1]
    #print("sense_diff, 590:610:",sense_diff[590:610])
    #print("antisense_diff, 590:610:",antisense_diff[590:610])
    print("sense_diff, %s:%s:" % (break_start,break_end),sense_diff[break_start:break_end])
    print("antisense_diff, %s:%s:" % (break_start,break_end),antisense_diff[break_start:break_end])

    #Plot pooled xdiff, ydiff, zdiff kdes
    plot_kdes(sense_diff,antisense_diff)
    print("Plotted kdes.") 
    
    bounds = (20.0,5.0,5.0)

    #Set parameters for unwrapping
    partial_multiple = 0.7
    partial_threshold = 0.45
    full_multiple = 1.0657
    full_threshold = 1.0

    #Scan through and find all these deltas. Cluster, and subtract.
    diffs = np.concatenate((sense_diff,antisense_diff),axis=0)
    km_x = sklearn.cluster.KMeans(n_cluster = 3).fit(diffs[:,0])
    km_y = sklearn.cluster.KMeans(n_cluster = 3).fit(diffs[:,1])
    km_z = sklearn.cluster.KMeans(n_cluster = 3).fit(diffs[:,2])
    

    
    #i = approx_median
    #Move forward from the median and check for deltas
    for i in range(approx_median,len(sense_diff)):
        #Check if moved off into each direction
        for j in range(len(bounds)):
            #Check funky oscillation with frequency of 1 timestep
            #if np.abs(sense_diff[i][j] > 0.25*bounds[j]):
                #See if adding next delta and this delta solves problems
                #oscillation_ratio = np.abs((sense_diff[i][j])/sense_diff[i][j])

            if sense_diff[i][j] > full_threshold*bounds[j]:
                correctors[0,i+1:,j] -= full_multiple*bounds[j]
            elif sense_diff[i][j] > partial_threshold*bounds[j]:
                correctors[0,i+1:,j] -= partial_multiple*bounds[j]
            elif sense_diff[i][j] < -full_threshold*bounds[j]:
                correctors[0,i+1:,j] += full_multiple*bounds[j]
            elif sense_diff[i][j] < -partial_threshold*bounds[j]:
                correctors[0,i+1:,j] += partial_multiple*bounds[j]


        #Now check antisense strand
        for j in range(len(bounds)):
            if antisense_diff[i][j] > full_threshold*bounds[j]:
                correctors[1,i+1:,j] -= full_multiple*bounds[j]
            elif antisense_diff[i][j] > partial_threshold*bounds[j]:
                correctors[1,i+1:,j] -= partial_multiple*bounds[j]
            elif antisense_diff[i][j] < -full_threshold*bounds[j]:
                correctors[1,i+1:,j] += full_multiple*bounds[j]
            elif antisense_diff[i][j] < -partial_threshold*bounds[j]:
                correctors[1,i+1:,j] += partial_multiple*bounds[j]

    #Move backward from the median and check for deltas
    for i in reversed([x for x in range(approx_median)]):
        #Check if moved off into each direction
        for j in range(len(bounds)):
            if sense_diff[i][j] > full_threshold*bounds[j]:
                correctors[0,i+1:,j] -= full_multiple*bounds[j]
            elif sense_diff[i][j] > partial_threshold*bounds[j]:
                correctors[0,i+1:,j] -= partial_multiple*bounds[j]
            elif sense_diff[i][j] < -full_threshold*bounds[j]:
                correctors[0,i+1:,j] += full_multiple*bounds[j]
            elif sense_diff[i][j] < -partial_threshold*bounds[j]:
                correctors[0,i+1:,j] += partial_multiple*bounds[j]
        #Now check antisense strand
        for j in range(len(bounds)):
            if antisense_diff[i][j] > full_threshold*bounds[j]:
                correctors[1,i+1:,j] -= full_multiple*bounds[j]
            elif antisense_diff[i][j] > partial_threshold*bounds[j]:
                correctors[1,i+1:,j] -= partial_multiple*bounds[j]
            elif antisense_diff[i][j] < -full_threshold*bounds[j]:
                correctors[1,i+1:,j] += full_multiple*bounds[j]
            elif antisense_diff[i][j] < -partial_threshold*bounds[j]:
                correctors[1,i+1:,j] += partial_multiple*bounds[j]
        '''
        for j in range(len(bounds)):
            if sense_diff[i][j] > 0.45*bounds[j]:
                correctors[0,:i+1,j] -= 0.5*bounds[j] + 1
            elif sense_diff[i][j] < -0.45*bounds[j]:
                correctors[0,:i+1,j] += 0.5*bounds[j] + 1
        #Now check antisense strand
        for j in range(len(bounds)):
            if antisense_diff[i][j] > 0.45*bounds[j]:
                correctors[1,:i+1,j] -= 0.5*bounds[j] + 1
            elif antisense_diff[i][j] < -0.45*bounds[j]:
                correctors[1,:i+1,j] += 0.5*bounds[j] + 1
        '''

    print("Correctors constructed.")
    #Construct new trajectory with corrections in place
    uw_xyz = np.copy(xyz)
    uw_xyz = np.swapaxes(uw_xyz,0,1)
    uw_xyz[atoms_1] += correctors[0]
    uw_xyz[atoms_2] += correctors[1]
    uw_xyz = np.swapaxes(uw_xyz,0,1)
    np.save('correctors.npy',correctors)    


    print("Corrected xyz created. Validating...")
    sense_com = np.average(uw_xyz[:,atoms_1],axis=1) #Shape: n_ts, 3
    antisense_com = np.average(uw_xyz[:,atoms_2],axis=1)
    sa_dist = antisense_com - sense_com
    dist_lens = np.linalg.norm(sa_dist,axis=1)
    print("Corrected Max distances:",np.max(dist_lens,axis=0))
    print("Corrected Min distances:",np.min(dist_lens,axis=0))
    print("Corrected distances, %s:%s:" % (break_start,break_end),dist_lens[break_start:break_end])
    
    sense_diff = sense_com[1:] - sense_com[:-1]
    antisense_diff = antisense_com[1:] - antisense_com[:-1]
    
    print("Corrected sense_diff, %s:%s:" % (break_start,break_end),sense_diff[break_start:break_end])
    print("Corrected antisense_diff, %s:%s:" % (break_start,break_end),antisense_diff[break_start:break_end])
   
    #ucl = a.unitcell_lengths
    ucl = np.zeros((len(a.xyz),3))
    ucl[:,0] = bounds[0]
    ucl[:,1] = bounds[1]
    ucl[:,2] = bounds[2] 
    traj = md.Trajectory(xyz = uw_xyz, topology = a.topology, time = a.time)
    print("Saving trajectory to uw_output.h5")
    traj.save_hdf5('uw_output.h5')

    
    print("Saved converted traj.")


def plot_kdes(sense_diff,antisense_diff):
    """Create kdes of x, y, and z diffs using ISJ.

    Parameters:
        sense_diff: *np.array*, shape: (n_ts,3)
            Frame-by-frame displacement of sense strand
        antisense_diff: *np.array*, shape: (n_ts,3)
            Frame-by-frame displacement of antisense strand

    """
    pooled_diffs = np.concatenate((sense_diff,antisense_diff),axis=0)

    if not os.path.isdir('outputs'):
        os.mkdir('outputs')
    #Plot X
    plt.clf()
    fig,ax = plt.subplots(layout='constrained')    

    kde_mins = np.min(pooled_diffs,axis = 0)
    kde_mins = kde_mins - np.abs(kde_mins) # lower than the min

    kde_maxes = np.max(pooled_diffs,axis = 0)
    kde_maxes = kde_maxes + np.abs(kde_maxes) # lower than the min

    obs_x = zfit.Space("X_diffs", limits=(kde_mins[0], kde_maxes[0]))
    obs_y = zfit.Space("Y_diffs", limits=(kde_mins[1], kde_maxes[1]))
    obs_z = zfit.Space("Z_diffs", limits=(kde_mins[2], kde_maxes[2]))
    kde_x = zfit.pdf.KDE1DimISJ(pooled_diffs[:,0], obs=obs_x)
    x = np.linspace(kde_mins[0],kde_maxes[0],1000)
    un_color = '#009E73'
    kde_vals = kde_x.pdf(x)
    ax.plot(x, kde_vals, linewidth=0.5, color=un_color)
    ax.fill_between(x, kde_vals, color=un_color, alpha=0.25)
    ax.set_xlabel("X differences (nm)")
    ax.set_ylabel("Marginal probability")
    ax.set_yscale('log')
    fig.savefig('outputs/x_diff_kde.png',dpi=600)
    fig.savefig('outputs/x_diff_kde.svg')
    
    #Plot Y
    plt.clf()
    fig,ax = plt.subplots(layout='constrained')    
    
    y = np.linspace(kde_mins[1],kde_maxes[1],1000)
    un_color = '#E69F00'
    kde_y = zfit.pdf.KDE1DimISJ(pooled_diffs[:,1], obs=obs_y)
    kde_vals = kde_y.pdf(y)
    ax.plot(y, kde_vals, linewidth=0.5, color=un_color)
    ax.fill_between(y, kde_vals, color=un_color, alpha=0.25)
    ax.set_xlabel("Y differences (nm)")
    ax.set_ylabel("Marginal probability")
    ax.set_yscale('log')
    fig.savefig('outputs/y_diff_kde.png',dpi=600)
    fig.savefig('outputs/y_diff_kde.svg')
    
    #Plot Z
    plt.clf()
    fig,ax = plt.subplots(layout='constrained')    
    z = np.linspace(kde_mins[2],kde_maxes[2],1000)
    kde_z = zfit.pdf.KDE1DimISJ(pooled_diffs[:,2], obs=obs_z)
    un_color = '#D55E00'
    kde_vals = kde_z.pdf(z)
    ax.plot(z, kde_vals, linewidth=0.5, color=un_color)
    ax.fill_between(z, kde_vals, color=un_color, alpha=0.25)
    ax.set_xlabel("Z differences (nm)")
    ax.set_ylabel("Marginal probability")
    ax.set_yscale('log')
    fig.savefig('outputs/z_diff_kde.png',dpi=600)
    fig.savefig('outputs/z_diff_kde.svg')
    
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f','--file', type=str, default='output.h5',
        help='Path to hdf5 output')
    args = parser.parse_args()

    main()

