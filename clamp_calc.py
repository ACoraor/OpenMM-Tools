# Module clamp_calc
# Created by Aria Coraor, 4/29/24

import numpy as np
import os
import sys

# from merlin import files
# from nucl_load import nucl_load,np_conv_lammps
# from subtool import inb4
import pandas as pd
import mdtraj as md

# from recalc import calc_re
# import rvect_calc
import argparse


def main():
    """Calculate the inter-phosphate distance in the XY plane."""

    # Load hdf5s
    if os.path.isfile('output.h5'):
        out = md.load_hdf5("output.h5")
    else:
        out = None
    
    if os.path.isfile('full_output.h5'):
        full_out = md.load_hdf5("full_output.h5")
        if out is not None:
            xyz = np.concatenate((full_out.xyz, out.xyz), axis=0)
            top = out.topology
        else:
            xyz = full_out.xyz
            top = full_out.topology
    else:
        xyz = out.xyz
        top = out.topology

    chains = [ch for ch in top.chains]
    ch0 = chains[0]
    res = [r for r in ch0.residues]
    end_1A_phos = [at.index for at in res[1].atoms if "P" == at.name][0]
    end_2A_phos = [at.index for at in res[-1].atoms if "P" == at.name][0]
    bound_ends = np.array(((end_1A_phos,),(end_2A_phos,)))
    print("Bound_ends:",bound_ends)
    #Get atom inds using the procedure from run_openmm
    atom_inds = bound_ends
    # Calculate angles
    angs = calculate_clamp_dist(xyz, atom_inds)
    # Write to pandas dframe
    save_to_dframe(angs,args.name)

    print("Saved distances in %s. Exiting..." % os.getcwd())

def save_to_dframe(angs,name):
    """Save the fiber angles to log.dframe

    *Parameters*
            angs:  *np.array of floats*, shape: (n_ts)
                    Angle, in radians, of the fiber bending.
    """
    df = pd.read_csv("log.dframe")
    df[name] = angs
    df.to_csv("log.dframe", index=False, na_rep="NULL")


def calculate_clamp_dist(xyz, atom_inds):
    """Calculate distance between atom inds only in xy plane.

    *Parameters*
            xyz:  *np.array of floats*, shape: (n_ts,n_atoms,3)
                    Atomic positions.
            atom_inds: *array-like of int*, shape:3
                    Indices to calculate angles between
    """
    try:
        pos1 = np.average(xyz[:,atom_inds[0]],axis=1)
        pos2 = np.average(xyz[:,atom_inds[1]],axis=1)
        vec1 = pos2 - pos1
    except:
        new_xyz = xyz[:,atom_inds]
        vec1 = new_xyz[:,0] - new_xyz[:,1]
            
    #Account for periodic conditions in y, z
    #vec1[:,1] = np.where(vec1[:,1] < 2.5, vec1[:,1], 5.0 - vec1[:,1])    
    #vec1[:,2] = np.where(vec1[:,2] < 2.5, vec1[:,2], 5.0 - vec1[:,2])
    print("vec1.shape:",vec1.shape)
    
    #Calculate distances excluding z axis
    dists = np.linalg.norm(vec1[:,:-1],axis=1)
    return dists
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-d",
        "--dir",
        type=str,
        default=".",
        help="Path to simulation Processing directory.",
    )
    parser.add_argument(
        "-n",
        "--name",
        type=str,
        default="clamp_dist",
        help=("Name of distance cv. Options include 'clamp_dist'" ),
    )
    # parser.add_argument('-f','--file',type=str,help="Path to lammpstrj to analyze.")
    args = parser.parse_args()

    main()

