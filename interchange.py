#Module interchange
#Created by Aria Coraor
#Written 6/7/23

import copy
from openff.interchange import Interchange
from openff.toolkit import Molecule,ForceField
from openff.units import unit
import numpy as np
#import mdtraj as md
import argparse
import warnings
import os


def main():
	"""Interchange GROMACS-style AMBER-dyes FF into OpenMM"""

        fn = ("/project2/depablo/coraor/soft_002/amber_dyes/" 
                + "amber99sb_dyes.ff/ffnonbonded.itp")

        amber = ForceField(fn)

	#Read the raw text into a list. Split into columns, and strip each element.
	xyz, top = read_raw_pdb(args.file)

	

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-f','--file', type=str, default='input.pdb',
		help='PDB file to snip.')
	parser.add_argument('-o','--output', type=str, default='output_pdb.txt',
		help='Path to output pdb file.')
	args = parser.parse_args()

	main()

