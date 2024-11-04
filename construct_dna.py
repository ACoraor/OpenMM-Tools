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
    xyz, top = add_dyes(xyz, top, d_xyz, d_top, a_xyz, a_top, args.donor,
        args.acceptor,args.swap_strands)

    #Save new pdb
    pdb = md.formats.PDBTrajectoryFile(args.output,mode='w')
    pdb.write(10.*xyz,top)

    print("Labelled pdb written. Exiting...")

def add_dyes(xyz, top, donor_xyz, donor_top, acceptor_xyz, acceptor_top,
         donor_ind, acceptor_ind,swap_strands):
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
        swap_strands: *bool*
            If true, 'donor' is Cy5 and 'acceptor' is Cy3.
    Returns:
        labelled_xyz: *np.array*, shape: (n_atoms,3)
            Positions of particles in the labelled DNA pdb.
        labelled_top: *mdtraj.topology*
            Topology of the labelled DNA pdb.
    """
    #Swap out Cy3 and Cy5, if called for.
    if swap_strands:
        tmp_xyz = np.copy(donor_xyz)
        tmp_top = donor_top
        donor_xyz = acceptor_xyz
        donor_top = acceptor_top
        acceptor_xyz = tmp_xyz
        acceptor_top = tmp_top

    #First: create a truncated version with all the appropriate residue 
    #Atoms removed 
    trunc_xyz, trunc_top = remove_dna_atoms(xyz, top, args.donor, 
        args.acceptor)

    #Next: get the inlet and outlet phosphorus and oxygen positions
    d_phos_in, d_ox_in, d_ox_out, d_dna_com = get_phos_ox(xyz, top, 
        donor_ind, is_donor=True) 
    a_phos_in, a_ox_in, a_ox_out, a_dna_com = get_phos_ox(xyz, top, 
        acceptor_ind, is_donor=False)


    #Rotate and translate the dye molecules to line up the phosphates
    affine_donor_xyz, affine_donor_top = affine_translate(d_phos_in,
        d_ox_in, d_ox_out, donor_xyz, donor_top,d_dna_com)
    affine_acceptor_xyz, affine_acceptor_top = affine_translate(a_phos_in,
        a_ox_in, a_ox_out, acceptor_xyz,acceptor_top,a_dna_com)

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
            if chain_i == 0 and res_i == donor_ind:
                donor_res = [r for r in affine_donor_top.residues][0]
                if not swap_strand:
                    labelled_top.add_residue(name='C3N',chain=curr_ch)
                else:
                    labelled_top.add_residue(name='C5N',chain=curr_ch)
                curr_res = labelled_top.residue(-1)
                for d_at in donor_res.atoms:
                    labelled_top.add_atom(name=d_at.name,
                        element=d_at.element, residue=curr_res)
                labelled_xyz += [affine_donor_xyz]
            elif chain_i == 1 and res_i == acceptor_ind:
                acc_res = [r for r in affine_acceptor_top.residues][0]
                if not swap_strand:
                    labelled_top.add_residue(name='C5N',chain=curr_ch)
                else:
                    labelled_top.add_residue(name='C3N',chain=curr_ch)
                curr_res = labelled_top.residue(-1)
                for acc_at in acc_res.atoms:
                    labelled_top.add_atom(name=acc_at.name,
                        element=acc_at.element, residue=curr_res)
                labelled_xyz += [affine_acceptor_xyz]
            
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

def remove_dna_atoms(xyz, top, donor_ind, acceptor_ind):
    """Remove relevant atoms from each donor/acceptor residue in the
    pdb structure.

    Parameters:
        xyz: *np.array*, shape: (n_atoms,3)
            Positions of particles in the input DNA pdb.
        top: *mdtraj.topology*
            Topology of the input DNA pdb.
        donor_ind: *int*
            Index of the donor with respect to its strand, 5' direction. If 
			negative, do not remove.
        acceptor_ind: *int*
            Index of the acceptor with respect to its strand, 5' direction. If
			negative, do not remove.
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
    if donor_ind >= 0:
        donor_residue = ch1_res[donor_ind]
        for at in donor_residue.atoms:
            rm_inds.append(at.index)
    ch2_res = [r for r in chains[1].residues]
    if acceptor_ind >= 0:
        acceptor_residue = ch2_res[acceptor_ind]
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
        dna_com: *np.array*, shape: (3)
            Average position of the labelled DNA residue.
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
    #outgoing_res = residues[index+1]
    dye_res = residues[index]
    #outgoing_atoms = [at for at in outgoing_res.atoms]
    dye_res_atoms = [at for at in dye_res.atoms]
    #for at in outgoing_atoms:
    #    if at.name == "P":
    #        phos_out = xyz[at.index]
    for at in dye_res_atoms:
        if at.name == "P":
            phos_in = xyz[at.index]
        if at.name == "OP1":
            ox_in = xyz[at.index]
        elif at.name == "O3'":
            ox_out = xyz[at.index]

    dna_com = np.average(xyz[[at.index for at in dye_res_atoms]],axis=0)

    return phos_in, ox_in, ox_out, dna_com

def affine_translate(phos_in, ox_in, ox_out,dye_xyz,dye_top,dna_com):
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
        dna_com: *np.array*, shape: (3)
            Average position of the labelled DNA residue.
    Returns:
        dye_xyz: *np.array*, shape: (n_atoms,3)
            Positions of particles in the rotated dye with linkers.
    """
    fvu0 = np.eye(3)
    quat0 = tu2rotquat(1e-5,np.array([1,0,0]))
    quat = np.copy(quat0)

    #Translate dye to 0
    dye_xyz = dye_xyz - np.average(dye_xyz,axis=0)

    #Get relevant dye atoms
    dye_ats = [at for at in dye_top.atoms]
    for at in dye_ats:
        if at.name == "P":
            dye_P = dye_xyz[at.index]
        elif at.name == "O5'":
            dye_o5 = dye_xyz[at.index]
        elif at.name == "O3'":
            dye_o3 = dye_xyz[at.index]

    #Now: Calculate angle to rotate along PO vector to align O5
    p_O3_vect = dye_o3 - dye_P
    dna_p_O3_vect = ox_out - phos_in

    #Project dna vector onto the plane of rotation
    p_O3_vect /= np.linalg.norm(p_O3_vect)
    dna_p_O3_vect /= np.linalg.norm(dna_p_O3_vect)
    
    rot_axis = np.cross(p_O3_vect,dna_p_O3_vect)
    rot_axis /= np.linalg.norm(rot_axis)

    dotp = np.dot(p_O3_vect,dna_p_O3_vect)
    if np.abs(dotp) > 1:
        dotp = dotp/np.abs(dotp)
    theta = np.arccos(dotp)
    q = tu2rotquat(theta,rot_axis)

    #Rotate dye into place   
    dye_xyz = np.array([quat_vec_rot(row, q) for row in dye_xyz])
    dye_ats = [at for at in dye_top.atoms]
    for at in dye_ats:
        if at.name == "P":
            dye_P = dye_xyz[at.index]
        elif at.name == "O5'":
            dye_o5 = dye_xyz[at.index]
        elif at.name == "O3'":
            dye_o3 = dye_xyz[at.index]
    
    p_O3_vect = dye_o3 - dye_P
    p_O3_vect /= np.linalg.norm(p_O3_vect)
    print("Final dot between p_O3 and DNA p_o3 vectors (should be 1):",
            np.dot(p_O3_vect,dna_p_O3_vect))

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
    parser.add_argument('-f','--file', type=str, default='input.pdb',
        help='Path to PDB of ideal DNA structure, with correct endings.')
    parser.add_argument('-o','--output', type=str, default='labelled_dna.pdb',
        help='Path to output labelled PDB.')
    parser.add_argument('-d','--donor', type=int, 
        default=-1, help='Index of Cy3 in donor (Ashort) strand. Pass -1 for no donor. Default=15')
    parser.add_argument('-a','--acceptor', type=int, 
        default=-1, help='Index of Cy5 in Acceptor (B) strand. Pass -1 for no acceptor. Default=23 (B10)')
    parser.add_argument('-s','--swap_strand', action='store_true', 
            help="If true, 'donor' is Cy5 and 'acceptor' is Cy3.")
    args = parser.parse_args()

    main()

