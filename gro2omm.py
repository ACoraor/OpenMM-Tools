#Module construct_dna
#Created by Aria Coraor
#Written 6/6/23

import copy
import numpy as np
import mdtraj as md
import argparse
import warnings
import os
import xml.etree.ElementTree as ET

MASS = {"C":12.0107, "N": 14.0067, "O": 15.999, "P": 30.973762}

def main():
    """Modify the existing Cy3 and Cy5 pdbs to have a phosphate linker,
    for direct inclusion in openMM structure. Phosphate group is coming
    from the 5' side.

    Truncate the lysine linkage and add an additional alkyl linkage.
    """
    
    #Load base DNA xml
    dna_xml = ET.parse

    print("PDBs loaded. Modifying...")


def foo(bar):
    """Init.

    Parameters:
        bar: *np.array*
            description
    Returns:
        biz: *np.array*
            description
    """
    pass
	


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-x','--xml', type=str, 
        default='test.xml', help='Path to test xml file for format.')
    parser.add_argument('-a','--acceptor', type=str, 
        default='C5N_L1R.pdb', help='Path to basic Cy5 pdb.')
    parser.add_argument('-a','--acceptor', type=str, 
        default='C5N_L1R.pdb', help='Path to basic Cy5 pdb.')
    args = parser.parse_args()

    main()

