# Module distance
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
    """Calculate the angle of the fiber bending and write to pandas dframe."""

    # Load hdf5s
    full_out = md.load_hdf5("full_output.h5")
    out = md.load_hdf5("output.h5")
    # Calculate angles
    xyz = np.concatenate((full_out.xyz, out.xyz), axis=0)
    if args.name == "cy5_N+1":
        atom_inds = (2383, 918)  # Fiber angle at Cy5 for soft_041
    elif args.name == "cy5_interc":
        atom_inds = (2383,888)
    elif args.name == "cy5_5-stack":
        atom_inds = (2383,2348)
    angs = calculate_dist(xyz, atom_inds)
    # Write to pandas dframe
    save_to_dframe(angs,args.name)

    print("Saved angles in %s. Exiting..." % args.dir)

def save_to_dframe(angs,name):
    """Save the fiber angles to log.dframe

    *Parameters*
            angs:  *np.array of floats*, shape: (n_ts)
                    Angle, in radians, of the fiber bending.
    """
    df = pd.read_csv("log.dframe")
    df[name] = angs
    df.to_csv("log.dframe", index=False, na_rep="NULL")


def calculate_dist(xyz, atom_inds):
    """Calculate angle between atomic coordinates.

    *Parameters*
            xyz:  *np.array of floats*, shape: (n_ts,n_atoms,3)
                    Atomic positions.
            atom_inds: *array-like of int*, shape:3
                    Indices to calculate angles between
    """
    new_xyz = xyz[:,atom_inds]
    vec1 = new_xyz[:,0] - new_xyz[:,1]
    dists = np.linalg.norm(vec1,axis=1)
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
        default="cy5_N+1",
        help="Name of distance cv. Options include 'cy5_N+1', 'cy5_interc',
            cy5_5-stack",
    )
    # parser.add_argument('-f','--file',type=str,help="Path to lammpstrj to analyze.")
    args = parser.parse_args()

    main()

