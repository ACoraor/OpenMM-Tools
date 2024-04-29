#Module bending_angle
#Created by Aria Coraor, 4/29/24

import numpy as np
import os
import sys
#from merlin import files
#from nucl_load import nucl_load,np_conv_lammps
#from subtool import inb4
import pandas as pd
import mdtraj as md
#from recalc import calc_re
#import rvect_calc
import argparse

def main():
	"""Calculate the angle of the fiber bending and write to pandas dframe."""

	#Load hdf5s
        full_out = md.load_hdf5('full_output.h5')
        out = md.load_hdf5('output.h5')
        #Calculate angles
        xyz = np.concatenate((full_out.xyz,out.xyz),axis=0)
        atom_inds = (55,896,1580) #Fiber angle at Cy5 for soft_041
        angs = calculate_angles(xyz,atom_inds)
        #Write to pandas dframe
        save_to_dframe(angs)


	print("Saved tetras in %s. Exiting..." % d)

def save_to_dframe(angs):
	"""Save the fiber angles to log.dframe
		
	*Parameters*
		angs:  *np.array of floats*, shape: (n_ts)
			Angle, in radians, of the fiber bending.
	"""
        df = pd.read_csv('log.dframe')
        df['fiber_angle'] = angs
        df.to_csv('log.dframe',index=False,na_rep='NULL')

def calculate_angles(xyz,atom_inds):
	"""Calculate angle between atomic coordinates.
		
	*Parameters*
		xyz:  *np.array of floats*, shape: (n_ts,n_atoms,3)
			Atomic positions.
                atom_inds: *array-like of int*, shape:3
                        Indices to calculate angles between
	"""
        df = pd.read_csv('log.dframe')
        df['fiber_angle'] = angs
        df.to_csv('log.dframe',index=False,na_rep='NULL')

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('-d','--dir', type=str, default=".",
		help='Path to simulation Processing directory.')
	#parser.add_argument('-f','--file',type=str,help="Path to lammpstrj to analyze.")
	args = parser.parse_args()
	
	main()
