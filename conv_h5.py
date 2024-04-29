#Module conv_h5
#Created by Aria Coraor
#Written 3/28/24

import numpy as np
import mdtraj as md
import argparse


def main():
    """Convert hdf5 trajectory into pdb trajectory."""
    a = md.load_hdf5(args.file)
    dot_index = args.file.find('.')
    fn = args.file[:dot_index]+".pdb"
    a.save_pdb(fn)
    print("Saved converted %s." % fn)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f','--file', type=str, default='output.h5',
        help='Path to hdf5 output')
    args = parser.parse_args()

    main()

