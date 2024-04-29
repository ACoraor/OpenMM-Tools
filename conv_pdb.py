#Module conv_pdb
#Created by Aria Coraor
#Written 9/7/23

import copy
import numpy as np
#import mdtraj as md
import argparse
import warnings
import os
from subtool import psub


def main():
        """Concatenate all xmls into one trajectory."""
        #Generate filenames
        cmd = "sed -n "

        n_steps = 2000
        starts = [4 + 74718*i for i in range(n_steps)]
        ends = [4209 + 74718*i for i in range(n_steps)]
        cmd += "'"
        for (s,e) in zip(starts,ends):
            cmd += "%s,%sp;" % (s,e)
        outfile = "data.xyz"
        cmd = cmd[:-1]
        cmd += "' %s > %s" % (args.file,outfile)
        print(cmd)
        psub(cmd)
        print("Done.")
        quit()
        

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
        
        traj.save_xyz(args.output)
        no_solvent = traj.remove_solvent()
        no_solvent.save_xyz(args.ns_output)
        print("Data saved. Exiting.")

	


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f','--file', type=str, default='output.pdb',
        help='Path to broken pdb')
    args = parser.parse_args()

    main()

