#Module make_traj
#Created by Aria Coraor
#Written 6/6/23

import copy
import numpy as np
import mdtraj as md
import argparse
import warnings
import os


def main():
        """Concatenate all xmls into one trajectory."""

        #Generate filenames
        fs = [ele for ele in os.listdir(args.dir) if ('state' in ele
            and '.xml' in ele)]

        n_states = len(fs)
        fns = [os.path.join(args.dir,"state_%s.xml" % n) for n in 
            range(n_states)]
        
        #Load initial pdbs
        init_pdb = md.load_pdb(args.init)
        top = init_pdb.top

        #Actually load data
        print("Loading data.")
        xmls = []
        for i, fn in enumerate(fns):
            if i % 50 == 0:
                print("On traj",i)
                xmls.append(md.load_xml(fn,top=top)) 
        #xmls = [md.load_xml(fn,top=top) for fn in fns]
        traj = md.join(xmls)
        #traj = md.load_xml(fns[0],top=top)
        #for i,fn in enumerate(fns[1:]):
        #    if i % 10 == 0:
        #        print("On traj",i)
        #    traj = md.join((traj,md.load_xml(fn,top=top)))
        
        traj.save_lammpstrj(args.output)
        print("Data saved. Exiting.")

	


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-d','--dir', type=str, default='states',
        help='Path to savestate directory')
    parser.add_argument('-i','--init', type=str, default='init.pdb',
        help='Path to PDB of initialized structure.')
    parser.add_argument('-o','--output', type=str, default='cat_state.lammpsrj',
        help='Name for final lammpstrj')
    args = parser.parse_args()

    main()

