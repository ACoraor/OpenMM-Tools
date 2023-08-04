#Module construct_dna
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
    """Create OpenMM Forcefield files based on the ff files from AMBER
    dyes.
    """
    #Load basic pdb
    cy3_amber_top, cy5_amber_top = load_amber_pdbs(args.donor,args.acceptor)

    #Load new pdbs
    cy3_top, cy5_top = load_label_pdbs()

    #Load AMBER-dyes forcefield
    amber_charges, amber_bonds, amber_impropers = load_amber_ffs(args.ff
        ,cy3_amber_top,cy5_amber_top,cy3_top,cy5_top)
    
    #Write ff xml
    write_ff_xml(amber_charges,amber_bonds,amber_impropers,cy3_top,cy5_top,
        cy3_amber_top,cy5_amber_top,args.gaff,args.output)

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
            DataFrame mapping atom names to GAFF classes & charges for 
            amber donor, acceptor
        amber_bonds: *list of pd.DataFrame*, shape: 2
            DataFrame with bonds between adjacent atoms.
        amber_impropers: *list of pd.DataFrame*, shape: 2
            DataFrame with improper dihedrals.
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

    return ([cy3_charge_df,cy5_charge_df], [cy3_bond_df,cy5_bond_df],
        [cy3_improper_df, cy5_improper_df])


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
    cy3 = md.load_pdb(donor_path)
    cy5 = md.load_pdb(acceptor_path)

    return cy3.top, cy5.top

def write_ff_xml(amber_charges,amber_bonds,amber_impropers,cy3_top,cy5_top,
        cy3_amber_top,cy5_amber_top,gaff,output):
    """Create OpenMM-formatted forcefield file based on AMBER-dyes
    parameters for the cy3 and cy5 residues.

    Parameters:
        amber_charges: *list of pd.DataFrame*, shape: 2
            DataFrame mapping atom names to GAFF classes & charges for 
            amber donor, acceptor
        amber_bonds: *list of pd.DataFrame*, shape: 2
            DataFrame with bonds between adjacent atoms.
        amber_impropers: *list of pd.DataFrame*, shape: 2
            DataFrame with improper dihedrals.
        cy3_top: *md.Topology*
            Loaded cy3 label topology.
        cy5_top: *md.Topology*
            Loaded cy5 label topology.
        cy3_amber_top: *md.Topology*
            Loaded Cy3_L1N Topology.
        cy5_amber_top: *md.Topology*
            Loaded Cy5_L1N Topology.
        gaff: *str*
            Path to GAFF xml.
        output: *str*
            Path to output xml.
    """
    #Read gaff xml
    xml = ET.parse(gaff)
    xml_root = xml.getroot()    
    
    #Get names of cy3, cy5 atoms
    cy3_names = set([at.name for at in cy3_top.atoms])
    cy5_names = set([at.name for at in cy5_top.atoms])
    #cy3_names.sort()
    #cy5_names.sort()
    
    
    #Create map from atom name to class, charge.
    cy3_charges = amber_charges[0]
    cy5_charges = amber_charges[1]
    cy3_bonds = amber_bonds[0]
    cy5_bonds = amber_bonds[1]
    cy3_charge_dict = {}
    for i in range(len(cy3_charges['atom1'])):
        at1 = cy3_charges['atom1'][i]
        at2 = cy3_charges['atom2'][i]
        qres = cy3_charges['QqmmmRESP'][i]
        cy3_charge_dict[at1] = [at2,qres]
        
    cy5_charge_dict = {}
    for i in range(len(cy5_charges['atom1'])):
        at1 = cy5_charges['atom1'][i]
        at2 = cy5_charges['atom2'][i]
        qres = cy5_charges['QqmmmRESP'][i]
        cy5_charge_dict[at1] = [at2,qres]

    #Create new xml
    new_root = ET.Element('ForceField')
    new_xml = ET.ElementTree(new_root)
    info = ET.SubElement(new_root,'Info')
    date_gen = ET.SubElement(info,'DateGenerated')
    date_gen.set('date',str(datetime.date.today()))
    atomtypes = ET.SubElement(new_root,'AtomTypes')
    for name in cy3_names:
        if name not in cy3_charge_dict:
            if name[0] == "C" and int(name[1:]) > 85:
                #Aliphatic carbon
                aliph_carbon = ET.SubElement(atomtypes,'Type')
                aliph_carbon.set('element','C')
                aliph_carbon.set('name',name)
                aliph_carbon.set('class','c3')
                aliph_carbon.set('mass','12.01')
                #aliph_carbon.set('charge','-0.073')
            if name[0] == "H" and int(name[1:]) > 85:
                #Assume hcg, q=0.0229
                aliph_hydrogen = ET.SubElement(atomtypes,'Type')
                aliph_hydrogen.set('element','H')
                aliph_hydrogen.set('name',name)
                aliph_hydrogen.set('class','hc')
                aliph_hydrogen.set('mass','1.008')
                #aliph_hydrogen.set('charge','0.0229')
    res_root = ET.SubElement(new_root,'Residues')

    #Create residue template
    cy3_res = ET.SubElement(res_root,"Residue")
    cy3_res.set('name','C3N')
    cy5_res = ET.SubElement(res_root,'Residue')
    cy5_res.set('name','C5N')

    #Loop through Cy3, cy5, creating residue
    for at in cy3_top.atoms:
        new_at = ET.SubElement(cy3_res,"Atom")
        new_at.set('name',at.name)
        if at.name in cy3_charge_dict:
            charge_data = cy3_charge_dict[at.name]
        elif "C" in at.name:
            charge_data = ['c3g','-0.073']
        elif "H" in at.name:
            charge_data = ['hcg','0.0229']
        elif "O" in at.name:
            charge_data = ['og','-0.6458']
        elif "P" in at.name:
            charge_data = ['DNA-Pg','1.1659']
        new_at.set('type',charge_data[0][:-1])
        new_at.set('charge', str(np.float32(charge_data[1])))
    #Add bonds
    atom1 = cy3_bonds['a1']
    atom2 = cy3_bonds['a2']
    for i in range(len(atom1)):
        if atom1[i] in cy3_names and atom2[i] in cy3_names:
            new_bond = ET.SubElement(cy3_res,"Bond")
            new_bond.set('atomName1',atom1[i])
            new_bond.set('atomName2',atom2[i])

    #Add external bonds
    ext_p = ET.SubElement(cy3_res,'ExternalBond')
    ext_p.set('atomName',"P")
    ext_o = ET.SubElement(cy3_res,'ExternalBond')
    ext_o.set('atomName',"O3'")

    for at in cy5_top.atoms:
        new_at = ET.SubElement(cy5_res,"Atom")
        new_at.set('name',at.name)
        if at.name in cy3_charge_dict:
            charge_data = cy5_charge_dict[at.name]
        elif "C" in at.name:
            charge_data = ['c3g','-0.073']
        elif "H" in at.name:
            charge_data = ['hcg','0.0229']
        elif "O" in at.name:
            charge_data = ['og','-0.6458']
        elif "P" in at.name:
            charge_data = ['DNA-Pg','1.1659']
        new_at.set('type',charge_data[0][:-1])
        new_at.set('charge', str(np.float32(charge_data[1])))
    
    #Add bonds
    atom1 = cy5_bonds['a1']
    atom2 = cy5_bonds['a2']
    for i in range(len(atom1)):
        if atom1[i] in cy5_names and atom2[i] in cy5_names:
            new_bond = ET.SubElement(cy5_res,"Bond")
            new_bond.set('atomName1',atom1[i])
            new_bond.set('atomName2',atom2[i])

    #Add external bonds
    ext_p = ET.SubElement(cy5_res,'ExternalBond')
    ext_p.set('atomName',"P")
    ext_o = ET.SubElement(cy5_res,'ExternalBond')
    ext_o.set('atomName',"O3'")
    #Populate xml_shell
 
    basic_str = ET.tostring(new_root,encoding='unicode')
    #with open(output,'w') as f:   
    #    f.write(basic_str)
    #basic_str = str(basic_str)#.replace('\n','')
    reparsed = minidom.parseString(basic_str)
    xml_str = reparsed.toprettyxml(newl='\n')
    #Replace silly strings
    xml_str = xml_str.replace('\n\t\t\n','\n\t\t')
    xml_str = xml_str.replace('\n\t\t    \n\t\t','\n\t\t')
    xml_str = xml_str.replace('\n\t\t    \n\t\t','\n\t\t')
    xml_str = xml_str.replace('\n\t\t  \n\t','\n\t')
    with open(output,'w') as f:
        f.write(xml_str)


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
        default='amber_subset', help='Path to AMBER-dyes forcefield file.')
    parser.add_argument('-d','--donor', type=str, 
        default='C3N_L1R.pdb', help='Path to basic Cy3 pdb.')
    parser.add_argument('-a','--acceptor', type=str, 
        default='C5N_L1R.pdb', help='Path to basic Cy5 pdb.')
    parser.add_argument('-g','--gaff',type=str,default = "gaff.xml",
	help="Path to GAFF parameter file.")
    parser.add_argument('-o','--output',type=str,default = "gaff-aec.xml",
	help="Path for output gaff xml.")
    args = parser.parse_args()

    main()

