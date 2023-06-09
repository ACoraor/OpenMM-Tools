#Module construct_dna
#Created by Aria Coraor
#Written 6/6/23

import copy
import numpy as np
import mdtraj as md
import argparse
import warnings
import os


def main():
        """Create a construct PDB by attaching Cy3 and Cy5 to the specified
        NA residues."""
    
    #Load input pdb
    pdb = md.load_pdb(args.file)
    xyz = pdb.xyz
    top = pdb.top




def foo(bar):
    """Init.

    Parameters:
        bar: *np.array*
            description
    Returns:
        biz: *np.arra*
            description
    """
    pass
	


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f','--file', type=str, default='input.pdb',
        help='Path to PDB of ideal structure, with correct endings.')
    parser.add_argument('-o','--output', type=str, default='labelled.pdb',
        help='Path to output labelled PDB.')
    parser.add_argument('-d','--donor', type=int, 
        default=15, help='Index of Cy3 in donor (Ashort) strand. Default=15')
    parser.add_argument('-a','--acceptor', type=int, 
        default=23, help='Index of Cy5 in Acceptor (B) strand. Default=23 (B10)')
    args = parser.parse_args()

    main()

