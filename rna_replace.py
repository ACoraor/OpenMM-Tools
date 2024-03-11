#Module rna_replace
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
    """Create a chimeric PDB by attaching construct to RNA"""
    #Assume cy3_dna.pdb and cy5_dna.pdb are in this directory!    

    #Load input RNA pdb
    pdb = md.load_pdb(args.file)
    xyz = pdb.xyz[0]
    top = pdb.top

    #Load construct pdbs
    construct_pdb = md.load_pdb(args.label)
    construct_xyz = construct_pdb.xyz[0]
    construct_top = construct_pdb.top

    

    #Add construct in the correct locations
    xyz, top = add_dyes(xyz, top, construct_xyz,construct_top,args.start,args.end)

    #Save new pdb
    pdb = md.formats.PDBTrajectoryFile(args.output,mode='w')
    pdb.write(10.*xyz,top)

    print("Labelled pdb written. Exiting...")

def add_dyes(xyz, top, construct_xyz, construct_top,start,end):
    """Remove the appropriate atoms of the RNA residues to replace. Reorient
    the construct to the appropriate positions and orientations

    Parameters:
        xyz: *np.array*, shape: (n_atoms,3)
            Positions of particles in the input DNA pdb.
        top: *mdtraj.topology*
            Topology of the input DNA pdb.
        construct_xyz: *np.array*, shape: (n_atoms,3)
            Positions of particles in the Cy3 dye with linkers.
        construct_top: *mdtraj.topology*
            Topology of the Cy3 dye with linkers.
        start: *int*
            Index of residue start for label insertion
        end: *int*
            Index of residue end for label insertion
    Returns:
        xyz: *np.array*, shape: (n_atoms,3)
            Positions of particles in the chimeric pdb.
        top: *mdtraj.topology*
            Topology of the Chimeric pdb.
    """
    #First: create a truncated version with all the appropriate residue 
    #Atoms removed 
    trunc_xyz, trunc_top = remove_rna_atoms(xyz, top, start, 
        end)

    #Get inlet and outlet phosphorus positions
    phos_in, c2_pos,o5_pos, phos_out, start_com = get_phos_rna(xyz, top, 
        start,end) 
    
    #Rotate and translate the dye molecules to line up the phosphates
    affine_construct_xyz, affine_construct_top = affine_translate_rna(phos_in,
        c2_pos,o5_pos, phos_out, construct_xyz, construct_top,start_com)

    #Combine xyzs and tops:
    #labelled_xyz = np.concatenate((trunc_xyz, affine_donor_xyz, affine_acceptor_xyz))
    
    #Manually construct new topology obeying proper sequencing
    labelled_top = md.Topology()
    labelled_xyz = []#Concatenate-able list of np.array

    #loop over all residues
    res_index = 0
    old_chains = [ch for ch in trunc_top.chains]
    for chain_i, chain in enumerate(old_chains):
        labelled_top.add_chain()
        curr_ch = labelled_top.chain(-1)
        for res_i, res in enumerate(chain.residues):
            #If you hit the dye residue, first add the dye
            if chain_i == 0 and res_i == start:
                donor_res = [r for r in affine_construct_top.residues]
                for construct_res in donor_res:
                    labelled_top.add_residue(name=construct_res.name,chain=curr_ch)
                    curr_res = labelled_top.residue(-1)
                    for d_at in construct_res.atoms:
                        labelled_top.add_atom(name=d_at.name,element=d_at.element,
                            residue = curr_res
                labelled_xyz += [affine_construct_xyz]
            
            #Now, add the normal residue and atoms
            labelled_top.add_residue(name=res.name,chain=curr_ch)
            curr_res = labelled_top.residue(-1)
            for at in res.atoms:
                labelled_top.add_atom(name=at.name,element=at.element,
                    residue=curr_res)
                labelled_xyz += [np.reshape(trunc_xyz[at.index],(1,-1))]
    labelled_xyz = np.concatenate(labelled_xyz)
    return labelled_xyz, labelled_top
    print("Unreachable code! legacy version of same code.")
    #Donor is on chain 0, acceptor is on chain 1
    labelled_top = copy.deepcopy(trunc_top)
    labelled_chains = [ch for ch in labelled_top.chains]
    
    #Loop through residues to add donor
    for don_r in affine_donor_top.residues:
        #Add residue
        labelled_top.add_residue("C3N",labelled_chains[0])
        current_r = labelled_top.residue(-1)
        #Loop through atoms
        for don_atom in don_r.atoms:
            labelled_top.add_atom(don_atom.name,element=don_atom.element,
                residue=current_r)

    #Loop through residues to add acceptor
    for acc_r in affine_acceptor_top.residues:
        #Add residue
        labelled_top.add_residue(acc_r.name,labelled_chains[1])
        next_r = labelled_top.residue(-1)
        #Loop through atoms
        for acc_atom in acc_r.atoms:
            labelled_top.add_atom(acc_atom.name,element=acc_atom.element,
                residue=next_r)
    
    #Old-style joining, makes chain ends:
    #labelled_top = trunc_top.join(affine_donor_top)
    #labelled_top = labelled_top.join(affine_acceptor_top)

    return labelled_xyz, labelled_top

def remove_rna_atoms(xyz, top, start, end):
    """Remove relevant atoms from each donor/acceptor residue in the
    pdb structure.

    Parameters:
        xyz: *np.array*, shape: (n_atoms,3)
            Positions of particles in the input DNA pdb.
        top: *mdtraj.topology*
            Topology of the input DNA pdb.
        start: *int*
            Index of residue start for label insertion
        end: *int*
            Index of residue end for label insertion
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
    #donor_residue = ch1_res[donor_ind]
    #ch2_res = [r for r in chains[1].residues]
    #acceptor_residue = ch2_res[acceptor_ind]
    for rm_ind in range(start,end+1):
        rm_res = ch1_res[rm_ind]
        for at in rm_res.atoms:
            rm_inds.append(at.index)

    keep_inds = [i for i in range(len(all_atoms)) if i not in set(rm_inds)]
    trunc_top = top.subset(keep_inds)
    trunc_xyz = xyz[keep_inds]
   
    pdb = md.formats.PDBTrajectoryFile('rmdna.pdb',mode='w')
    pdb.write(10.*trunc_xyz,trunc_top)

    return trunc_xyz, trunc_top

def get_phos_rna(xyz, top,start,end):
    """Get the positions of the inlet and outlet phosphates and O?3? to 
    determine orientation of dyes.

    Parameters:
        xyz: *np.array*, shape: (n_atoms,3)
            Positions of particles in the input DNA pdb.
        top: *mdtraj.topology*
            Topology of the input DNA pdb.
        start: *int*
            Index of residue start for label insertion
        end: *int*
            Index of residue end for label insertion
    Returns:
        phos_st: *np.array*, shape: (3)
            Position of the inlet phosphorus from the 5' direction.
        c2_pos: *np.array*, shape: (3)
            Position of the c2' carbon at start ind.
        o5_pos: *np.array*, shape: (3)
            Position of the O5' oxyen direction.
        phos_out: *np.array*, shape: (3)
            Position of the outlet P at end.
        rna_com: *np.array*, shape: (3)
            Average position of the labelled DNA residue.
    """
    #We seek to exactly overlay the incoming 3' phosphate, and its OP1

    #We rotate the molecule until this is matched, and then rotate until 
    #The outlet O5' position is matched.
    chains = [ch for ch in top.chains]
    chain = chains[0]

    residues = [r for r in chain.residues]
    #outgoing_res = residues[index+1]
    dye_res = residues[start]
    #outgoing_atoms = [at for at in outgoing_res.atoms]
    dye_res_atoms = [at for at in dye_res.atoms]
    #for at in outgoing_atoms:
    #    if at.name == "P":
    #        phos_out = xyz[at.index]
    for at in dye_res_atoms:
        if at.name == "P":
            phos_in = xyz[at.index]
        if at.name == "O5'":
            o5_pos = xyz[at.index]
        if at.name == "C2'":
            c2_pos = xyz[at.index]

    #Get phos out
    for at in [at for at in residues[end].atoms]
        if at.name == "P":
            phos_out = xyz[at.index]

    rna_com = np.average(xyz[[at.index for at in dye_res_atoms]],axis=0)

    return phos_in, c2_pos, o5_pos, phos_out,rna_com

def affine_translate_rna(phos_in, c2_pos, o5_pos, phos_out,construct_xyz,construct_top,rna_com):
    """Rotate the dye until it aligns with the appropriate phosphate in/out

    Parameters:
        phos_in: *np.array*, shape: (3)
            Position of the inlet phosphorus from the 3' direction.
        c2_pos: *np.array*, shape: (3)
            Position of the c2' carbon at start ind.
        o5_pos: *np.array*, shape: (3)
            Position of the O5' oxyen direction.
        phos_out: *np.array*, shape: (3)
            Position of the outlet OP2 from the 5' direction.
        construct_xyz: *np.array*, shape: (n_atoms,3)
            Positions of particles in the dye with linkers.
        construct_top: *mdtraj.topology*
            Topology of the dye pdb.
        rna_com: *np.array*, shape: (3)
            Average position of the labelled DNA residue.
    Returns:
        affine_xyz: *np.array*, shape: (n_atoms,3)
            Positions of particles in the rotated dye with linkers.
        affine_top: *mdtraj.topology*
            Topology of the dye pdb.
    """
    fvu0 = np.eye(3)
    quat0 = tu2rotquat(1e-5,np.array([1,0,0]))
    quat = np.copy(quat0) #Running accumulator of quaternion

    #Translate last opposing-strand phosphate to 0
    conjugate_res = [res for res in construct_top.residues]
    conjugate_inds = [at.index for at in conjugate_res[-1].atoms]
    phos_ind = [at.index for at in conjugate_res[-1].atoms if at.name == "P"][0]

    construct_xyz = construct_xyz - construct_xyz[phos_ind]

    #Get relevant dye atoms
    construct_ats = [at for at in conjugate_res.atoms]
    for at in dye_ats:
        if at.name == "P":
            construct_P = construct_xyz[at.index]
            construct_phos_ind = at.index
        elif at.name == "C2'":
            construct_c2 = construct_xyz[at.index]
            construct_c2_ind = at.index
        elif at.name == "O5'":
            construct_o5 = construct_xyz[at.index]
            construct_o5_ind = at.index

    #Idea: Bring P --> C2, P --> O5, and C2 --> O5 vectors into alignment
    p_c2_vect = construct_c2 - construct_P
    dna_p_c2_vect = c2_pos - phos_in
    axis_vect = np.cross(p_c2_vect,dna_p_c2_vect)
    axis_vect /= np.linalg.norm(axis_vect)
    first_angle = vector_angle(p_c2_vect,dna_p_c2_vect)
    
    q = tu2rotquat(first_angle,axis_vect)
    dye_xyz = np.array([quat_vec_rot(row, q) for row in construct_xyz])
    #Now: P-->C2 axes are matched

    #Bring P --> O5 vectors into alignment
    construct_o5 = dye_xyz[construct_o5_ind]
    construct_P = dye_xyz[construct_phos_ind]
    construct_c2 = dye_xyz[construct_c2_ind]
    p_o5_vect = construct_o5 - construct_P
    dna_p_o5_vect = o5_pos - phos_in
    axis_vect = dna_p_c2_vect#Rotate along dna_p_c2_vect
    #Verify that dots are zero:
    print("dna-axis dot (should be 0):",np.dot(dna_p_o5_vect,axis_vect))    
    print("construct-axis dot (should be 0):",np.dot(p_o5_vect,axis_vect))    
    second_angle = vector_angle(p_o5_vect,dna_p_o5_vect)
    q = tu2rotquat(second_angle,axis_vect)

    dye_xyz = np.array([quat_vec_rot(row, q) for row in dye_xyz])
    construct_o5 = dye_xyz[construct_o5_ind]
    construct_P = dye_xyz[construct_phos_ind]
    construct_c2 = dye_xyz[construct_c2_ind]

    c2_o5_vect = construct_o5 - construct_P
    dna_c2_o5_vect = o5_pos - c2_pos

 
    print("Final dot between c2_O5 and DNA c2_o5 vectors (should be ~1):",
            np.dot(c2_o5_vect,dna_c2_o5_vect))
    fn = "affine_construct.pdb"
    pdb = md.formats.PDBTrajectoryFile(fn,mode='w')
    pdb.write(10.*dye_xyz,construct_top)

    return dye_xyz, construct_top

    print("Unreachable code!")
    #Dye is now rotated so phos_in - ox_out axis is correct. now rotate
    #the COM to the outside.
    dye_com = np.average(dye_xyz,axis=0)
    p_dna_vect = (dna_com - phos_in)/np.linalg.norm(dna_com-phos_in)
    p_dna_perp = p_dna_vect - np.dot(p_dna_vect,dna_p_O3_vect)*dna_p_O3_vect
    p_dna_perp /= np.linalg.norm(p_dna_perp)
    
    p_dye_vect = (dye_com - dye_P)/np.linalg.norm(dye_com - dye_P)
    p_dye_perp = p_dye_vect - np.dot(p_dye_vect,dna_p_O3_vect)*dna_p_O3_vect
    p_dye_perp /= np.linalg.norm(p_dye_perp)
    
    dotp = np.dot(p_dye_perp,p_dna_perp)
    if np.abs(dotp) > 1:
        dotp = dotp/np.abs(dotp)
    theta = np.arccos(dotp)
    #Rotate p_dye_perp by theta+pi/2
    q = tu2rotquat(theta + np.pi,dna_p_O3_vect)

    dye_xyz = np.array([quat_vec_rot(row, q) for row in dye_xyz])
    dye_ats = [at for at in dye_top.atoms]
    for at in dye_ats:
        if at.name == "P":
            dye_P = dye_xyz[at.index]
        elif at.name == "O5'":
            dye_o5 = dye_xyz[at.index]
        elif at.name == "O3'":
            dye_o3 = dye_xyz[at.index]
    
    p_dye_vect = (dye_com - dye_P)/np.linalg.norm(dye_com - dye_P)
    p_dye_perp = p_dye_vect - np.dot(p_dye_vect,dna_p_O3_vect)*dna_p_O3_vect
    p_dye_perp /= np.linalg.norm(p_dye_perp)
    dotp = np.dot(p_dye_perp,p_dna_perp)
    if np.abs(dotp) > 1:
        dotp = dotp/np.abs(dotp)
    print("Final dot between p_dye_perp and p_dna_perp vectors (should be -1):",
            np.dot(p_dye_perp,p_dna_perp))
    #If dotp > 0.5, rotate by pi
    if dotp > 0.9:
        print("Strange rotation error? Flipping dye COM.")
        q = tu2rotquat(np.pi,dna_p_O3_vect)
        dye_xyz = np.array([quat_vec_rot(row, q) for row in dye_xyz])
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

    #Write just dye to PDB
    n_atoms = len([at for at in dye_top.atoms])
    is_donor = n_atoms < 73
    if is_donor:
        fn = "affine_cy3.pdb"
    else:
        fn = "affine_cy5.pdb"
    pdb = md.formats.PDBTrajectoryFile(fn,mode='w')
    pdb.write(10.*dye_xyz,dye_top)
    return dye_xyz, dye_top

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
    parser.add_argument('-f','--file', type=str, default='EGFP_key.pdb',
        help='Path to PDB of ideal RNA structure.')
    parser.add_argument('-l','--label', type=str, default='82bp_labelled.pdb',
        help='Path to PDB of labelled DNA.')
    parser.add_argument('-o','--output', type=str, default='labelled_dna.pdb',
        help='Path to output chimeric PDB.')
    parser.add_argument('-s','--start', type=int, default=100,
        help='Index of first removed RNA base, from 0')
    parser.add_argument('-e','--end', type=int, default=111,
        help='Index of last removed RNA base, from 0')
    args = parser.parse_args()

    main()

