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
    phos_in, c2_pos,o5_pos, phos_out, start_com,c1_pos,n4_pos, o3_pos, o3_out, c1_out = get_phos_rna(xyz, top, 
        start,end) 
    
    #Rotate and translate the dye molecules to line up the phosphates
    affine_construct_xyz, affine_construct_top = affine_translate_rna(phos_in,
        c2_pos,o5_pos, phos_out, construct_xyz, construct_top,start_com,
            c1_pos,n4_pos,o3_pos,o3_out,c1_out)

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
                            residue = curr_res)
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
    for rm_ind in range(start-1,end+1):
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
        c1_pos: *np.array*, shape: (3)
            STUB
        n4_pos: *np.array*, shape: (3)
            STUB
        o3_pos: *np.array*, shape: (3)
            STUB
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
        if at.name == "C1'":
            c1_pos = xyz[at.index]
        if at.name == "N4":
            n4_pos = xyz[at.index]
        if at.name == "O5'":
            o5_pos = xyz[at.index]
        if at.name == "C2'":
            c2_pos = xyz[at.index]
        if at.name == "O3'":
            o3_pos = xyz[at.index]

    #Get next phosphate?
    #for at in [at for at in residues[start+1].atoms]:
    #    if at.name == "P":
    #        phos_in = xyz[at.index]#Account for weird Openmm phosphate issue?
    #Get phos out
    for at in [at for at in residues[end].atoms]:
        if at.name == "P":
            phos_out = xyz[at.index]
        if at.name == "O3'":
            o3_out = xyz[at.index]
        if at.name == "C1'":
            c1_out = xyz[at.index]

    rna_com = np.average(xyz[[at.index for at in dye_res_atoms]],axis=0)

    return phos_in, c2_pos, o5_pos, phos_out,rna_com,c1_pos,n4_pos,o3_pos,o3_out, c1_out

def affine_translate_rna(phos_in, c2_pos, o5_pos, phos_out,construct_xyz,
        construct_top,rna_com, c1_pos, n4_pos,o3_pos,o3_out,c1_out):
    """Rotate the dye until it aligns with the appropriate phosphate in/out

    Parameters:
        phos_in: *np.array*, shape: (3)
            Position of the inlet phosphorus from the 3' direction.
        c2_pos: *np.array*, shape: (3)
            Position of the c2' carbon at start ind.
        o5_pos: *np.array*, shape: (3)
            Position of the O5' oxyen direction.
        phos_out: *np.array*, shape: (3)
            Position of the end-1 residue's P.
        construct_xyz: *np.array*, shape: (n_atoms,3)
            Positions of particles in the dye with linkers.
        construct_top: *mdtraj.topology*
            Topology of the dye pdb.
        rna_com: *np.array*, shape: (3)
            Average position of the labelled DNA residue.
        c1_pos: *np.array*, shape: (3)
            STUB
        n4_pos: *np.array*, shape: (3)
            STUB
        o3_pos: *np.array*, shape: (3)
            STUB
        o3_out: *np.array*, shape: (3)
            Phosphate position of end - 1 rna base.
        c1_out: *np.array*, shape: (3)
            Phosphate position of end - 1 rna base.
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
    start_ind = len(conjugate_res)//2 - 17
    conjugate_inds = [at.index for at in conjugate_res[start_ind].atoms]
    phos_ind = [at.index for at in conjugate_res[start_ind].atoms if at.name == "P"][0]

    construct_xyz = construct_xyz - construct_xyz[phos_ind]

    #Get relevant dye atoms
    construct_ats = [at for at in conjugate_res[start_ind].atoms]
    for at in construct_ats:
        if at.name == "P":
            construct_P = construct_xyz[at.index]
            construct_phos_ind = at.index
        elif at.name == "C2'":
            construct_c2 = construct_xyz[at.index]
            construct_c2_ind = at.index
        elif at.name == "C1'":
            construct_c1 = construct_xyz[at.index]
            construct_c1_ind = at.index
        elif at.name == "N4":
            construct_n4 = construct_xyz[at.index]
            construct_n4_ind = at.index
        elif at.name == "O5'":
            construct_o5 = construct_xyz[at.index]
            construct_o5_ind = at.index
        elif at.name == "O3'":
            construct_o3 = construct_xyz[at.index]
            construct_o3_ind = at.index

    #Idea: bring C1' --> N4 in line
    '''
    c1_n4_vect =construct_n4 - construct_c1
    dna_c1_n4_vect = n4_pos - c1_pos
    axis_vect = np.cross(c1_n4_vect,dna_c1_n4_vect)
    axis_vect /= np.linalg.norm(axis_vect)
    first_angle = vector_angle(c1_n4_vect,dna_c1_n4_vect)
    
    q = tu2rotquat(first_angle,axis_vect)
    dye_xyz = np.array([quat_vec_rot(row, q) for row in construct_xyz])
    '''
    
    #Idea: Bring P --> C1, P --> O5, and C1 --> O5 vectors into alignment
    p_c1_vect = construct_c1 - construct_P
    dna_p_c1_vect = c1_pos - phos_in
    axis_vect = np.cross(p_c1_vect,dna_p_c1_vect)
    axis_vect /= np.linalg.norm(axis_vect)
    first_angle = vector_angle(p_c1_vect,dna_p_c1_vect)
    
    q = tu2rotquat(first_angle,axis_vect)
    dye_xyz = np.array([quat_vec_rot(row, q) for row in construct_xyz])
    #Now: P-->C1 axes are matched, i.e. 5' side --> sugar aligned
    '''
    #Bring C1' --> P vectors into alignment
    construct_o5 = dye_xyz[construct_o5_ind]
    construct_P = dye_xyz[construct_phos_ind]
    construct_c2 = dye_xyz[construct_c2_ind]
    construct_c1 = dye_xyz[construct_c1_ind]
    construct_n4 = dye_xyz[construct_n4_ind]

    c1_p_vect = construct_P - construct_c1
    dna_c1_p_vect = phos_in - c1_pos
    axis_vect = dna_c1_n4_vect # Rotate along previous axis
    #Remove components along this vector to determine angle
    c1_p_vect_comp = c1_p_vect - np.dot(c1_p_vect,axis_vect)
    dna_c1_p_vect_comp = dna_c1_p_vect - dna_c1_p_vect - np.dot(dna_c1_p_vect,axis_vect)
    
    #Now calculate remaining angle in-plane
    second_angle = vector_angle(c1_p_vect_comp,dna_c1_p_vect_comp)
    q = tu2rotquat(np.pi+second_angle + 37.5/180*np.pi,axis_vect)
    dye_xyz = np.array([quat_vec_rot(row, q) for row in dye_xyz])
    '''
    
    #Now bring P--> O3 axis in alignment 
    construct_o3 = dye_xyz[construct_o3_ind]
    construct_P = dye_xyz[construct_phos_ind]
    construct_c2 = dye_xyz[construct_c2_ind]
    construct_c1 = dye_xyz[construct_c1_ind]
    construct_n4 = dye_xyz[construct_n4_ind]
    axis_vect = dna_p_c1_vect # Rotate along previous axis
    axis_vect /= np.linalg.norm(axis_vect)
    p_o3_vect = construct_o3 - construct_P
    dna_p_o3_vect = o3_pos - phos_in
    
    #Remove components along this vector to determine angle
    p_o3_vect_comp = p_o3_vect - np.dot(p_o3_vect,axis_vect)*axis_vect
    dna_p_o3_vect_comp = dna_p_o3_vect - np.dot(dna_p_o3_vect,axis_vect)*axis_vect
    
    #Verify that dots are zero:
    print("dna-axis dot (should be 0):",np.dot(dna_p_o3_vect_comp,axis_vect))    
    print("construct-axis dot (should be 0):",np.dot(p_o3_vect_comp,axis_vect))    
    second_angle = vector_angle(p_o3_vect_comp,dna_p_o3_vect_comp)
    #q = tu2rotquat(np.pi-second_angle,axis_vect) #Hacked for v1
    q = tu2rotquat(second_angle,-axis_vect)
    dye_xyz = np.array([quat_vec_rot(row, q) for row in dye_xyz])
    
    #Final checking and saving
    construct_o3 = dye_xyz[construct_o3_ind]
    construct_P = dye_xyz[construct_phos_ind]
    construct_c2 = dye_xyz[construct_c2_ind]
    construct_c1 = dye_xyz[construct_c1_ind]
    construct_n4 = dye_xyz[construct_n4_ind]

    c2_o3_vect = construct_o3 - construct_c2
    dna_c2_o3_vect = o3_pos - c2_pos

    c2_o3_vect /= np.linalg.norm(c2_o3_vect)
    dna_c2_o3_vect /= np.linalg.norm(dna_c2_o3_vect)
    
    print("Final dot between c2_O3 and DNA c2_o3 vectors (should be ~1):",
            np.dot(c2_o3_vect,dna_c2_o3_vect))
    fn = "affine_construct.pdb"
    pdb = md.formats.PDBTrajectoryFile(fn,mode='w')
    pdb.write(10.*dye_xyz,construct_top)
    
    dye_xyz = dye_xyz + phos_in
    #Now, bend the RNA to tie in to the end
    end_ind = len(conjugate_res)//2 - 1
    conjugate_inds = [at.index for at in conjugate_res[end_ind].atoms]
    end_phos_ind = [at.index for at in conjugate_res[end_ind].atoms if at.name == "P"][0]

    #Reset to start = 0 for alignment
    #dye_xyz = dye_xyz - dye_xyz[phos_ind]

    for _ in range(5):   
     
        sheared_dye_xyz = np.copy(dye_xyz)
        
        #Calculate vector pointing from current pos to correct pos
        p_correction_vector = phos_out - dye_xyz[end_phos_ind]

        axis_vector = dye_xyz[end_phos_ind] - dye_xyz[phos_ind]
        axis_vector /= np.linalg.norm(axis_vector)

        #Positions along axis
        #axis_positions = np.array([np.dot(dye_xyz[ind],axis_vector) for ind in bridge_conj_inds])
        n_tilts = 16

        #Force allow the manipulated base to match!
        new_phos_res = conjugate_res[len(conjugate_res)//2-1-n_tilts]
        for at in new_phos_res.atoms:
            if at.name == "P":
                phos_ind = at.index
        
        #Get parabola parameters for p_correction_vector
        start_p_axis_pos = np.dot(sheared_dye_xyz[phos_ind],axis_vector)
        #axis_positions = axis_positions - start_p_axis_pos
        end_p_axis_pos = np.dot(dye_xyz[end_phos_ind],axis_vector) - start_p_axis_pos
        parabola_a = p_correction_vector/(end_p_axis_pos)

        deflection_axis = np.cross(axis_vector,p_correction_vector)
        deflection_axis /= np.linalg.norm(deflection_axis)

        #slopes = 2*np.linalg.norm(parabola_a) *axis_positions
        #deflections = np.arctan(slopes)
        res_rots = []
        for i in range(n_tilts):
            res_A = conjugate_res[len(conjugate_res)//2 - 1 - i]
            res_B = conjugate_res[len(conjugate_res)//2 + i]
            atoms = [at for at in res_A.atoms] + [at for at in res_B.atoms]
            #atoms = [at for at in conjugate_res[len(conjugate_res)//2-1-i]]+ 
            com = np.average([sheared_dye_xyz[at.index] for at in atoms],axis=0)
            com_axis_pos = np.dot(com,axis_vector) - start_p_axis_pos
            slope = np.linalg.norm(parabola_a)
            theta = np.arctan(slope)    
            print("Deflection angle:",theta*180/np.pi)
            q = tu2rotquat(theta,deflection_axis)
            for index in [at.index for at in atoms]:
                loc_vect = sheared_dye_xyz[index] - com
                new_vect = quat_vec_rot(loc_vect,q)
                com_shear = parabola_a * com_axis_pos
                sheared_dye_xyz[index] = new_vect + com + com_shear
                
        dye_xyz = np.copy(sheared_dye_xyz)    

    #for i,ind in enumerate(bridge_conj_inds):
    #    sheared_dye_xyz[ind] += parabola_a*axis_positions[i]
    #sheared_dye_xyz = sheared_dye_xyz + dye_xyz[phos_ind]
    fn = "sheared_construct.pdb"
    pdb = md.formats.PDBTrajectoryFile(fn,mode='w')
    pdb.write(10.*sheared_dye_xyz,construct_top)
    
    #construct_xyz = dye_xyz - dye_xyz[phos_ind]

    #Calculate phosphate displacement to get end to correct position
    
    #Calculate vector along start --> end phos to calculate addition

    return sheared_dye_xyz, construct_top
    '''
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
    '''


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
    parser.add_argument('-o','--output', type=str, default='chimera.pdb',
        help='Path to output chimeric PDB.')
    parser.add_argument('-s','--start', type=int, default=100,
        help='Index of first removed RNA base, from 0')
    parser.add_argument('-e','--end', type=int, default=111,
        help='Index of last removed RNA base, from 0')
    args = parser.parse_args()

    main()

