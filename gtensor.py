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
    if os.path.isfile('output.h5'):
        out = md.load_hdf5("output.h5")
        top = out.topology
    else:
        out = None
    
    if os.path.isfile('full_output.h5'):
        full_out = md.load_hdf5("full_output.h5")
        if out is not None:
            xyz = np.concatenate((full_out.xyz, out.xyz), axis=0)
        else:
            xyz = full_out.xyz
            top = full_out.topology
    else:
        xyz = out.xyz
        top = out.topology
    
    #Get indices of all atoms in our DNA
    dnas_resnames = set(["DA","DC","DT","DG","C5N","C3N"])
    atom_inds = []
    #HACK: Only use chain1
    chains = [ch for ch in top.chains]
    
    for at in chains[0].atoms:
        if at.residue.name in dnas_resnames:
            atom_inds.append(at.index)
    atom_inds = np.array(atom_inds)

    # Calculate angles
    print("Calculating gtensor components.")
    gtensor_comps = calculate_gtensor(xyz, atom_inds)
    gtensor_comps = np.transpose(gtensor_comps)
    print("gtensor_comps shape:",gtensor_comps.shape)
    print("Calculation done. Saving...")
    names = ['gyration_lz2','gyration_ly2','gyration_lx2',
            'gyration_b_reduced','gyration_c_reduced','gyration_k2_reduced',
            'gyration_rg2']
    # Write to pandas dframe
    for data,name in zip(gtensor_comps,names):
        print("data shape for %s:" % name,data.shape)
        #print(data)
        save_to_dframe(data,name)

    #save_to_dframe(gtensor_comps[0],'gtensor_lx2')
    #save_to_dframe(gtensor_comps[1],'gtensor_ly2')
    #save_to_dframe(gtensor_comps[2],'gtensor_lz2')

    print("Saved gtensor comps in %s. Exiting..." % os.getcwd())

def save_to_dframe(angs,name):
    """Save the fiber angles to log.dframe

    *Parameters*
            angs:  *np.array of floats*, shape: (n_ts)
                    Angle, in radians, of the fiber bending.
    """
    df = pd.read_csv("log.dframe")
    df[name] = angs
    df.to_csv("log.dframe", index=False, na_rep="NULL")


def calculate_gtensor(xyz, atom_inds):
    """Calculate the principal moments of the gyration tensor.

    *Parameters*
            xyz:  *np.array of floats*, shape: (n_ts,n_atoms,3)
                    Atomic positions.
            atom_inds: *array-like of int*
                    Indices to include in gtensor calculation
    *Returns*
            gtensor_comps: *np.array*, shape: 3
                    Each principal moment of the gtensor, sorted by size
                    descending
    """
    #STUB
    gtensor_comps = np.zeros(3)
    xyz = xyz[:,atom_inds]
    COG = np.average(xyz, axis=1)
    # print("COG.shape:",COG.shape)
    centered_pos = np.zeros(xyz.shape)
    n_atoms = xyz.shape[1]

    boxdim = (20.0,5.0,5.0)

    # print("Centering.")
    for i in range(xyz.shape[0]):
        centered_pos[i, :, :] = xyz[i] - COG[i]
    # print("Centered_pos:",centered_pos)

    # print("Computing tensor.")
    gyration_tensor = np.zeros((centered_pos.shape[0], 3, 3))
    for ts_i in range(centered_pos.shape[0]):
        for nuc_j in range(centered_pos.shape[1]):
            current_pos = centered_pos[ts_i, nuc_j]
            gyration_tensor[ts_i] += np.matmul(
                np.reshape(current_pos, (3, 1)), np.reshape(current_pos, (1, 3))
            )
        gyration_tensor[ts_i] = gyration_tensor[ts_i] / n_atoms

    print("Calculating eigendecomposition.")
    # Diagonalize gyration tensor
    #gyration_tensor = np.sum(gyration_tensor,axis=0)
    print("Gyration tensor shape:",gyration_tensor.shape)
    gyr_eigenvalues, gyr_eigenvectors = np.linalg.eig(gyration_tensor)
    # print("gyr_eigenvalues.shape:",gyr_eigenvalues.shape)
    # print("gyr_eigenvectors.shape:",gyr_eigenvectors.shape)
    # Sort eigenvalues
    principal_moments = np.sort(gyr_eigenvalues, axis=1)
    print("Principal moments shape:",principal_moments.shape)

    rg2 = np.sum(principal_moments, axis=1)
    # Calculate b, c, k2
    # print("Last frame eigenvalues:",gyr_eigenvalues[-1])

    b = principal_moments[:, 2] - 0.5 * (
        principal_moments[:, 0] + principal_moments[:, 1]
    )
    c = principal_moments[:, 1] - principal_moments[:, 0]
    # Calculate reduced parameters by Rg2
    b_reduced = b / rg2
    c_reduced = c / rg2

    k2 = (
        1.5
        * (
            np.sum(principal_moments**2, axis=1)
            / (np.sum(principal_moments, axis=1) ** 2)
        )
        - 0.5
    )

    # Return results
    gtensor_comps = np.concatenate(
        (
            principal_moments,
            np.reshape(b_reduced, (-1, 1)),
            np.reshape(c_reduced, (-1, 1)),
            np.reshape(k2, (-1, 1)),
            np.reshape(rg2, (-1, 1)),
        ),
        axis=1,
    )

    return gtensor_comps
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-d",
        "--dir",
        type=str,
        default=".",
        help="Path to simulation Processing directory.",
    )
    # parser.add_argument('-f','--file',type=str,help="Path to lammpstrj to analyze.")
    args = parser.parse_args()

    main()

