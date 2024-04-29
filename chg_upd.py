#Module chg_upd
#Created by Aria Coraor
#Written 6/6/23

import copy
import numpy as np
import mdtraj as md
import pandas as pd
import argparse
import warnings
import datetime
import os
from io import StringIO
import xml.etree.ElementTree as ET
from xml.dom import minidom


MASS = {"C":12.0107, "N": 14.0067, "O": 15.999, "P": 30.973762}

def main():
    """Update gaff xml with Multiwfn charges.
    """

    #Load existing xml
    xml = load_ff_xml(args.gaff)

    #Load charges
    cy3_labels, cy3_chg = load_mwfn_charges(args.charges_don)
    cy5_labels, cy5_chg = load_mwfn_charges(args.charges_acc)
    
    #Write ff xml
    write_ff_xml(xml,cy3_labels,cy3_chg,cy5_labels,cy5_chg,args.output)

    print("DNA dye pdbs written! Exiting...")


def load_ff_xml(gaff_path):
    """Load the default gaff forcefield as plaintext

    Parameters:
        gaff_path: *str*
            Path to previously produced gaff output.
    Returns:
        curr_gaff: *str*
            Fully loaded current element tree, split by line
    """
    with open(gaff_path,'r') as f:
        curr_gaff = f.readlines()
    return curr_gaff

def load_mwfn_charges(charge_path):
    """Load the Multiwfn charges.

    Parameters:
        charge_path: *str*
            Path to previously produced Multiwfn output.
    Returns:
        labels: *np.array*, shape: (n_atoms)
            Elements for all atoms, in sequence
        charges: *np.array*, shape: (n_atoms,4)
            Charges for all atoms, in sequence
    """
    charges = np.loadtxt(charge_path,dtype=str)
    labels = charges[:,0]
    charges = charges[:,1:].astype(float)
    return labels,charges


def write_ff_xml(curr_gaff,cy3_labels,cy3_chg,cy5_labels,cy5_chg,output):
    """Create OpenMM-formatted forcefield file based on AMBER-dyes
    parameters with updated charges

    Parameters:
        curr_gaff: *str*
            Fully loaded current element tree, split by line
        cy3_labels: *np.array*, shape: (n_atoms)
            elements for all atoms, in Cy3 sequence
        cy3_chg: *np.array*, shape: (n_atoms,4)
            Charges for all atoms, in Cy3 sequence
        cy5_labels: *np.array*, shape: (n_atoms)
            elements for all atoms, in Cy5 sequence
        cy5_chg: *np.array*, shape: (n_atoms,4)
            Charges for all atoms, in Cy5 sequence
    """
    #Combine hydrogens with bonded oxygen
    #For Cy3: Hydrogen 76 attached to Oxygen 63
    #For Cy3: Hydrogen 77 attached to oxygen 73
    try:
        cy3_chg[63] = cy3_chg[63] + cy3_chg[76]
        cy3_chg[73] = cy3_chg[73] + cy3_chg[77]
        cy3_chg = cy3_chg[:-2]
        cy3_labels = cy3_labels[:-2]
    except:
        print("Not combining hydrogen charges.")

    #For Cy5: Hydrogen 80 attached to Oxygen 79
    #For Cy5: Hydrogen 81 attached to Oxygen 78
    try:
        cy5_chg[79] = cy5_chg[79] + cy5_chg[80]
        cy5_chg[78] = cy5_chg[78] + cy5_chg[81]
        cy5_chg = cy5_chg[:-2]
        cy5_labels = cy5_labels[:-2]
    except:
        print("Not combining hydrogen charges.")
    
    #Force-parse the plaintext, replacing charges line by line
    atom_i = 0
    molecule_i = 0
    
    new_gaff = []

    for j,line in enumerate(curr_gaff):
        if "charge" in line:
            charge_index = line.find('charge')
            charge_end = line.find('"/>')
            #Get new charge
            if molecule_i == 0:
                new_charge = cy3_chg[atom_i,3]
            else:
                new_charge = cy5_chg[atom_i,3]
            new_line = line[:charge_index] + 'charge="%s"/>\n' % (new_charge)
            atom_i = atom_i + 1
        else:
            new_line = line

        #If finished cy3, confirm this
        if molecule_i == 0 and atom_i >= len(cy3_chg):
            atom_i = 0
            molecule_i = 1         
            
        new_gaff.append(new_line)


    with open(output,'w') as f:
        for line in new_gaff:
            f.write(line)


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
    parser.add_argument('-d','--charges_don', type=str, 
        default='orca.cy3.chg', help='Path to Multiwfn RESP output for Cy3')
    parser.add_argument('-a','--charges_acc', type=str, 
        default='orca.cy5.chg', help='Path to Multiwfn RESP output for Cy5')
    parser.add_argument('-g','--gaff',type=str,default = "gaff-isolated.xml",
	help="Path to GAFF parameter file.")
    parser.add_argument('-o','--output',type=str,default = "gaff-isolated_v2.xml",
	help="Path for output gaff xml.")
    args = parser.parse_args()

    main()

