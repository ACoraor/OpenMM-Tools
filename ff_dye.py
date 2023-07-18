#Module construct_dna
#Created by Aria Coraor
#Written 6/6/23

import copy
import numpy as np
import mdtraj as md
import pandas as pd
import argparse
import warnings
import os
from io import StringIO
import xml.etree.ElementTree as ET
from xml.dom import minidom


MASS = {"C":12.0107, "N": 14.0067, "O": 15.999, "P": 30.973762}

def main():
    """Create OpenMM Forcefield files based on the ff files from AMBER
    dyes.
    """
    #Load basic pdb
    cy3_amber_top, cy5_amber_top = load_amber_pdbs(args.donor,args.acceptor)

    #Load new pdbs
    cy3_top, cy5_top = load_label_pdbs()

    #Load AMBER-dyes forcefield
    amber_dyes = load_amber_ffs(args.ff,cy3_amber_top,cy5_amber_top,
        cy3_top,cy5_top)
    
    #Write ff xml
    write_ff_xml(amber_dyes,cy3_top,cy5_top,cy3_amber_top,cy5_amber_top)

    print("DNA dye pdbs written! Exiting...")


def load_amber_ffs(dye_path,cy3_amber_top,cy5_amber_top,cy3_top,cy5_top):
    """Load AMBER-dyes gromacs-style forcefield file.

    Parameters:
        dye_path: *str*
            Path to AMBER dye parameter root directory, typically
            amber99sb_dyes.ff.
        cy3_amber_top: *md.Topology*
            Loaded Cy3_L1N Topology.
        cy5_amber_top: *md.Topology*
            Loaded Cy5_L1N Topology.
        cy3_top: *md.Topology*
            Loaded cy3 label topology.
        cy5_top: *md.Topology*
            Loaded cy5 label topology.
    Returns:
        amber_charges: *list of pd.DataFrame*, shape: 2
            DataFrame mapping atom names to charges for amber donor
            acceptor
        amber_dyes: *list of str*, shape: ????
            Parsed AMBER-dyes parameters. Details undescribed.
    """
    #Notes from AMBER-dyes:
    #Charge parameters and residue atom information can be found in the 
    #amber99sb_dyes.ff/aminoacids.rtp.
    #Bonded and nonbonded parameters are located in 
    #amber99sb_dyes.ff/ffbonded.itp and amber99sb_dyes.ff/ffnonbonded.itp, 
    #respectively.
    #Atom types are located in amber99sb_dyes.ff/atomtypes.atp

    #aminoacids.arn: atom renaming."
    #aminoacids.hdb: unparseable residue details."
    #aminoacids.vsd: Virtual site details
    #aminoacids.c.tdb and n.tdb: Empty
    #aminoacids.r2b: Residue to buildingblock table?

    #Create paths to parameter files
    #charge_res_path = os.path.join(dye_path,'aminoacids.rtp')
    cy3_rtp_path = os.path.join(dye_path,'cy3.rtp')
    cy5_rtp_path = os.path.join(dye_path,'cy5.rtp')

    #Get names of cy3, cy5 atoms
    cy3_names = list(set([at.name for at in cy3_top.atoms]))
    cy5_names = list(set([at.name for at in cy5_top.atoms]))
    cy3_names.sort()
    cy5_names.sort()

   
    #Parse amino acids for pandas
    with open(cy3_rtp_path,'r') as f:
        cy3_rtp = f.readlines()
    with open(cy5_rtp_path,'r') as f:
        cy5_rtp = f.readlines()
   
    #Manually parse the rtps
    cy3_charges = cy3_rtp[5:74]
    cy3_bonds = cy3_rtp[75:147]
    cy3_impropers = cy3_rtp[148:167]

    cy5_charges = cy5_rtp[5:78]
    cy5_bonds = cy5_rtp[79:155]
    cy5_impropers = cy5_rtp[156:177]
     
    #Split and process:
    cy3_charges = [ele.strip().split() for ele in cy3_charges]
    cy3_bonds = [ele.strip().split() for ele in cy3_bonds]
    cy3_impropers = [ele.strip().split() for ele in cy3_impropers]
    cy5_charges = [ele.strip().split() for ele in cy5_charges]
    cy5_bonds = [ele.strip().split() for ele in cy5_bonds]
    cy5_impropers = [ele.strip().split() for ele in cy5_impropers]
    
    cy3_atom1 = [row[0] for row in cy3_charges]
    cy3_atom2 = [row[1] for row in cy3_charges]
    cy3_resp = [float(row[2]) for row in cy3_charges]

    cy3_charge_df = pd.DataFrame()
    cy3_charge_df['atom1'] = cy3_atom1
    cy3_charge_df['atom2'] = cy3_atom2
    cy3_charge_df['QqmmmRESP'] = cy3_resp
    cy3_charge_df.to_csv('cy3_charge.csv',index=False)

    cy3_bond_a1 = [row[0] for row in cy3_bonds]
    cy3_bond_a2 = [row[1] for row in cy3_bonds]
    cy3_bond_df = pd.DataFrame()
    cy3_bond_df['a1'] = cy3_bond_a1
    cy3_bond_df['a2'] = cy3_bond_a2
    cy3_bond_df.to_csv('cy3_bond.csv',index=False)
    

    cy3_improper_a1 = [row[0] for row in cy3_impropers]
    cy3_improper_a2 = [row[1] for row in cy3_impropers]
    cy3_improper_a3 = [row[2] for row in cy3_impropers]
    cy3_improper_a4 = [row[3] for row in cy3_impropers]
    cy3_improper_df = pd.DataFrame()
    cy3_improper_df['a1'] = cy3_improper_a1
    cy3_improper_df['a2'] = cy3_improper_a2
    cy3_improper_df['a3'] = cy3_improper_a3
    cy3_improper_df['a4'] = cy3_improper_a4
    cy3_improper_df.to_csv('cy3_improper.csv',index=False)
    #Do the same for cy5
    cy5_atom1 = [row[0] for row in cy5_charges]
    cy5_atom2 = [row[1] for row in cy5_charges]
    cy5_resp = [float(row[2]) for row in cy5_charges]

    cy5_charge_df = pd.DataFrame()
    cy5_charge_df['atom1'] = cy5_atom1
    cy5_charge_df['atom2'] = cy5_atom2
    cy5_charge_df['QqmmmRESP'] = cy5_resp
    cy5_charge_df.to_csv('cy5_charge.csv',index=False)

    cy5_bond_a1 = [row[0] for row in cy5_bonds]
    cy5_bond_a2 = [row[1] for row in cy5_bonds]
    cy5_bond_df = pd.DataFrame()
    cy5_bond_df['a1'] = cy5_bond_a1
    cy5_bond_df['a2'] = cy5_bond_a2
    cy5_bond_df.to_csv('cy5_bond.csv',index=False)
    

    cy5_improper_a1 = [row[0] for row in cy5_impropers]
    cy5_improper_a2 = [row[1] for row in cy5_impropers]
    cy5_improper_a3 = [row[2] for row in cy5_impropers]
    cy5_improper_a4 = [row[3] for row in cy5_impropers]
    cy5_improper_df = pd.DataFrame()
    cy5_improper_df['a1'] = cy5_improper_a1
    cy5_improper_df['a2'] = cy5_improper_a2
    cy5_improper_df['a3'] = cy5_improper_a3
    cy5_improper_df['a4'] = cy5_improper_a4
    cy5_improper_df.to_csv('cy5_improper.csv',index=False)

    print("STUB!")
    return dye_path


def load_amber_pdbs(donor_path,acceptor_path):
    """Load the pdbs at each path.

    Parameters:
        donor_path: *str*
            Path to Cy3_L1N.pdb
        acceptor_path: *str*
            Path to Cy5_L1N.pdb
    Returns:
        cy3_amber_top: *md.Topology*
            Loaded Cy3_L1N Topology.
        cy5_amber_top: *md.Topology*
            Loaded Cy5_L1N Topology.
    """
    cy3 = md.load_pdb(args.donor)
    cy5 = md.load_pdb(args.acceptor)

    return cy3.top, cy5.top


def load_label_pdbs(donor_path="cy3_dna.pdb",acceptor_path="cy5_dna.pdb"):
    """Load the processed dna label pdbs at each path.

    Parameters:
        donor_path: *str*
            Path to label-processed cy3 dye. Default: cy3_dna.pdb
        acceptor_path: *str*
            Path to label-processed cy5 dye. Default: cy5_dna.pdb
    Returns:
        cy3_top: *md.Topology*
            Loaded cy3 label topology.
        cy5_top: *md.Topology*
            Loaded cy5 label topology.
    """
    cy3 = md.load_pdb(args.donor)
    cy5 = md.load_pdb(args.acceptor)

    return cy3.top, cy5.top

def write_ff_xml(amber_dyes,cy3_top,cy5_top,cy3_amber_top,
        cy5_amber_top,out='dyes.xml'):
    """Create OpenMM-formatted forcefield file based on AMBER-dyes
    parameters for the cy3 and cy5 residues.

    Parameters:
        amber_dyes: *list of str*, shape: ????
            Parsed AMBER-dyes parameters. Details undescribed.
        cy3_top: *md.Topology*
            Loaded cy3 label topology.
        cy5_top: *md.Topology*
            Loaded cy5 label topology.
        cy3_amber_top: *md.Topology*
            Loaded Cy3_L1N Topology.
        cy5_amber_top: *md.Topology*
            Loaded Cy5_L1N Topology.
        out: *str*
            Path to outputfile. Default: 'dyes.xml'
    """
    xml_shell = make_xml_shell()

    #Populate xml_shell
    print("STUB!")



def make_xml_shell():
    """Create empty xml file in OpenMM Forcefield format.

    Returns:
        xml_shell: *xml.etree.ElementTree*
            Empty ff file with appropriate tags.
    """
    print("STUB!")
    root = minidom.Document()
    ff = root.createElement('ForceField')
    root.appendChild(ff)

    atom_types = root.createElement('AtomTypes')
    residues = root.createElement('Residues')
    ff.appendChild(atom_types)
    ff.appendChild(residues)

    
    xml_str = root.toprettyxml(indent='\t')
    with open("test.xml",'w') as f:
        f.write(xml_str)
    return root    

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
    parser.add_argument('-f','--ff', type=str, 
        default='ff', help='Path to AMBER-dyes forcefield file.')
    parser.add_argument('-d','--donor', type=str, 
        default='C3N_L1R.pdb', help='Path to basic Cy3 pdb.')
    parser.add_argument('-a','--acceptor', type=str, 
        default='C5N_L1R.pdb', help='Path to basic Cy5 pdb.')
    args = parser.parse_args()

    main()

