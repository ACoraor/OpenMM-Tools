#Module construct_dna
#Created by Aria Coraor
#Written 6/6/23

import copy
import numpy as np
import mdtraj as md
import argparse
import warnings
from vect_quat_util import *
import os


def main():
    """Create a construct PDB by attaching Cy3 and Cy5 to the specified
    NA residues."""
    #Assume cy3_dna.pdb and cy5_dna.pdb are in this directory!    

    #Load input pdb
    pdb = md.load_pdb(args.file)
    xyz = pdb.xyz[0]
    top = pdb.top

    #Load cy3 and cy5 pdbs
    d_pdb = md.load_pdb("cy3_dna.pdb")
    d_xyz = d_pdb.xyz[0]
    d_top = d_pdb.top
    
    a_pdb = md.load_pdb("cy5_dna.pdb")
    a_xyz = a_pdb.xyz[0]
    a_top = a_pdb.top

    #Add dyes in the correct locations
    print("STUB")
    xyz, top = add_dyes(xyz, top, d_xyz, d_top, a_xyz, a_top, args.donor,
        args.acceptor)

    #Save new pdb
    pdb = md.formats.PDBTrajectoryFile(args.output,mode='w')
    pdb.write(10.*xyz,top)

    print("Labelled pdb written. Exiting...")

def add_dyes(xyz, top, donor_xyz, donor_top, acceptor_xyz, acceptor_top,
         donor_ind, acceptor_ind):
    """Remove the appropriate atoms of the DNA residue to replace. Reorient
    the donor and acceptor to the appropriate residue and add to the pdb.

    Parameters:
        xyz: *np.array*, shape: (n_atoms,3)
            Positions of particles in the input DNA pdb.
        top: *mdtraj.topology*
            Topology of the input DNA pdb.
        donor_xyz: *np.array*, shape: (n_atoms,3)
            Positions of particles in the Cy3 dye with linkers.
        donor_top: *mdtraj.topology*
            Topology of the Cy3 dye with linkers.
        acceptor_xyz: *np.array*, shape: (n_atoms,3)
            Positions of particles in the Cy5 dye with linkers.
        acceptor_top: *mdtraj.topology*
            Topology of the Cy5 dye with linkers.
        donor_ind: *int*
            Index of the donor with respect to its strand, 5' direction.
        acceptor_ind: *int*
            Index of the acceptor with respect to its strand, 5' direction.
    Returns:
        labelled_xyz: *np.array*, shape: (n_atoms,3)
            Positions of particles in the labelled DNA pdb.
        labelled_top: *mdtraj.topology*
            Topology of the labelled DNA pdb.
    """
    #First: create a truncated version with all the appropriate residue 
    #Atoms removed 
    trunc_xyz, trunc_top = remove_dna_atoms(xyz, top, args.donor, 
        args.acceptor)

    #Next: get the inlet and outlet phosphorus and oxygen positions
    d_phos_in, d_ox_in, d_ox_out = get_phos_ox(xyz, top, 
        donor_ind, is_donor=True) 
    a_phos_in, a_ox_in, a_ox_out = get_phos_ox(xyz, top, 
        acceptor_ind, is_donor=False)

    #Rotate and translate the dye molecules to line up the phosphates
    affine_donor_xyz, affine_donor_top = affine_translate(d_phos_in,
        d_ox_in, d_ox_out, donor_xyz, donor_top)
    affine_acceptor_xyz, affine_acceptor_top = affine_translate(a_phos_in,
        a_ox_in, a_ox_out, acceptor_xyz,acceptor_top)

    #Combine xyzs and tops:
    labelled_xyz = np.concatenate((trunc_xyz, affine_donor_xyz, affine_acceptor_xyz))
    labelled_top = trunc_top + affine_donor_top + affine_acceptor_top

    return labelled_xyz, labelled_top

def remove_dna_atoms(xyz, top, donor_ind, acceptor_ind):
    """Remove relevant atoms from each donor/acceptor residue in the
    pdb structure.

    Parameters:
        xyz: *np.array*, shape: (n_atoms,3)
            Positions of particles in the input DNA pdb.
        top: *mdtraj.topology*
            Topology of the input DNA pdb.
        donor_ind: *int*
            Index of the donor with respect to its strand, 5' direction.
        acceptor_ind: *int*
            Index of the acceptor with respect to its strand, 5' direction.
    Returns:
        trunc_xyz: *np.array*, shape: (n_atoms,3)
            Positions of particles in the truncated DNA pdb.
        trunc_top: *mdtraj.topology*
            Topology of the truncated DNA pdb.
    """
    rm_inds = []

    all_atoms = [at for at in top.atoms]

    chains = [c for c in top.chains]
    ch1_res = [r for r in chains[0].residues]
    donor_residue = ch1_res[donor_ind]
    ch2_res = [r for r in chains[1].residues]
    acceptor_residue = ch2_res[acceptor_ind]

    for at in donor_residue.atoms:
        rm_inds.append(at.index)
    for at in acceptor_residue.atoms:
        rm_inds.append(at.index)

    keep_inds = [i for i in range(len(all_atoms)) if i not in set(rm_inds)]
    trunc_top = top.subset(keep_inds)
    trunc_xyz = xyz[keep_inds]
   
    pdb = md.formats.PDBTrajectoryFile('rmdna.pdb',mode='w')
    pdb.write(10.*trunc_xyz,trunc_top)

    return trunc_xyz, trunc_top

def get_phos_ox(xyz, top,index, is_donor):
    """Get the positions of the inlet and outlet phosphates and O?3? to 
    determine orientation of dyes.

    Parameters:
        xyz: *np.array*, shape: (n_atoms,3)
            Positions of particles in the input DNA pdb.
        top: *mdtraj.topology*
            Topology of the input DNA pdb.
        index: *int*
            Index of the dye with respect to its strand, 5' direction.
        is_donor: *bool*
            Whether or not this is on the donor or the acceptor strand.
    Returns:
        phos_in: *np.array*, shape: (3)
            Position of the inlet phosphorus from the 3' direction.
        ox_in: *np.array*, shape: (3)
            Position of the inlet O3' from the 3' direction.
        ox_out: *np.array*, shape: (3)
            Position of the outlet O5' from the 5' direction.
    """
    #We seek to exactly overlay the incoming 3' phosphate, and its OP1

    #We rotate the molecule until this is matched, and then rotate until 
    #The outlet O5' position is matched.
    chains = [ch for ch in top.chains]
    if is_donor:
        chain = chains[0]
    else:
        chain = chains[1]

    residues = [r for r in chain.residues]
    incoming_res = residues[index-1]
    dye_res = residues[index]
    incoming_atoms = [at for at in incoming_res.atoms]
    dye_res_atoms = [at for at in dye_res.atoms]
    for at in incoming_atoms:
        if at.name == "P":
            phos_in = xyz[at.index]
    for at in dye_res_atoms:
        if at.name == "O3'":
            ox_in = xyz[at.index]
        elif at.name == "O5'":
            ox_out = xyz[at.index]

    return phos_in, ox_in, ox_out

def affine_translate(phos_in, ox_in, ox_out,dye_xyz,dye_top):
    """Rotate the dye until it aligns with the appropriate phosphate in/out

    Parameters:
        phos_in: *np.array*, shape: (3)
            Position of the inlet phosphorus from the 3' direction.
        ox_in: *np.array*, shape: (3)
            Position of the inlet OP1 from the 3' direction.
        ox_out: *np.array*, shape: (3)
            Position of the outlet OP2 from the 5' direction.
        dye_xyz: *np.array*, shape: (n_atoms,3)
            Positions of particles in the dye with linkers.
        dye_top: *mdtraj.topology*
            Topology of the dye pdb.
    Returns:
        dye_xyz: *np.array*, shape: (n_atoms,3)
            Positions of particles in the rotated dye with linkers.
    """
    fvu0 = np.eye(3)
    quat0 = tu2rotquat(1e-5,np.array([1,0,0]))
    quat = np.copy(quat0)

    #Get relevant dye atoms
    dye_ats = [at for at in dye_top.atoms]
    for at in dye_ats:
        if at.name == "P":
            dye_P = dye_xyz[at.index]
        elif at.name == "O5'":
            dye_o5 = dye_xyz[at.index]
        elif at.name == "O3'":
            dye_o3 = dye_xyz[at.index]

    #Calculate angle to align P-O vectors
    in_po_vector = ox_in - phos_in
    in_po_vector = in_po_vector/np.linalg.norm(in_po_vector)
    current_po_vector = (dye_o5 - dye_P)/np.linalg.norm(dye_o5-dye_P)

    theta = np.arccos(np.dot(current_po_vector,in_po_vector))
    rot_vector = np.cross(in_po_vector,current_po_vector)
    rot_vector /= np.linalg.norm(rot_vector)
    
    q = tu2rotquat(theta,rot_vector)
    quat = quat_multiply(q,quat)
    fvu = quat_fvu_rot(fvu0,quat)
    
    #Rotate by quat
    dye_xyz = np.array([quat_vec_rot(row, quat) for row in dye_xyz])
    dye_ats = [at for at in dye_top.atoms]
    for at in dye_ats:
        if at.name == "P":
            dye_P = dye_xyz[at.index]
        elif at.name == "O5'":
            dye_o5 = dye_xyz[at.index]
        elif at.name == "O3'":
            dye_o3 = dye_xyz[at.index]

    #Now: Calculate angle to rotate along PO vector to align O5
    p_O5_vect = dye_o5 - dye_P
    dna_p_O5_vect = ox_out - phos_in

    #Project dna vector onto the plan of rotation
    p_O5_proj_vect = dna_p_O5_vect - np.dot(dna_p_O5_vect,
        in_po_vector)*in_po_vecotr
    
    theta = np.arccos(np.dot(p_O5_proj_vect,p_O5_vect))
    q = tu2rotquat(theta,in_po_vector)
    
    #Rotate dye into place   
    dye_xyz = np.array([quat_vec_rot(row, quat) for row in dye_xyz])
    dye_ats = [at for at in dye_top.atoms]
    for at in dye_ats:
        if at.name == "P":
            dye_P = dye_xyz[at.index]
        elif at.name == "O5'":
            dye_o5 = dye_xyz[at.index]
        elif at.name == "O3'":
            dye_o3 = dye_xyz[at.index]
    
    #Finally: translate dye into place
    disp = phos_in - dye_P
    
    dye_xyz = dye_xyz + disp
    return dye_xyz

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
        help='Path to PDB of ideal DNA structure, with correct endings.')
    parser.add_argument('-o','--output', type=str, default='labelled.pdb',
        help='Path to output labelled PDB.')
    parser.add_argument('-d','--donor', type=int, 
        default=15, help='Index of Cy3 in donor (Ashort) strand. Default=15')
    parser.add_argument('-a','--acceptor', type=int, 
        default=23, help='Index of Cy5 in Acceptor (B) strand. Default=23 (B10)')
    args = parser.parse_args()

    main()

