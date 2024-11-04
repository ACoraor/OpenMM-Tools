# Module distance
# Created by Aria Coraor, 4/29/24

import numpy as np
import os
import sys

# from merlin import files
# from nucl_load import nucl_load,np_conv_lammps
# from subtool import inb4
import pandas as pd
import mdtraj as md

# from recalc import calc_re
# import rvect_calc
import argparse


def main():
    """Calculate the average Darboux vector distribution for C3N and C5N 
    residues relative to the opposing and N+-1 basepairs."""

    # Load hdf5s
    if os.path.isfile('output.h5'):
        out = md.load_hdf5("output.h5")
    else:
        out = None
    
    if os.path.isfile('full_output.h5'):
        full_out = md.load_hdf5("full_output.h5")
        if out is not None:
            xyz = np.concatenate((full_out.xyz, out.xyz), axis=0)
            top = full_out.topology
        else:
            xyz = full_out.xyz
            top = full_out.topology

    else:
        xyz = out.xyz
        top = out.topology

    #Calculate the orientation of each dye residue in the Darboux frame
    t_vects, planar_t_vects, names = calc_dye_darboux(xyz,top)
    
    #Combine t_vects and planar t_vects before saving
    t_vects = t_vects + planar_t_vects

    # Write to pandas dframe
    save_to_dframe(t_vects,names)

    print("Saved distances in %s. Exiting..." % os.getcwd())

def save_to_dframe(t_vects,names):
    """Save the Darboux orientation vectors to log.dframe

    *Parameters*
            angs:  *np.array of floats*, shape: (n_ts)
                    Angle, in radians, of the fiber bending.
    """
    df = pd.read_csv("log.dframe")

    #Loop through names and add each column of t_vects in
    for i, name in enumerate(names):
        df[name+"_t1"] = t_vects[i][:,0]
        df[name+"_t2"] = t_vects[i][:,1]
        df[name+"_t3"] = t_vects[i][:,2]
    df.to_csv("log.dframe", index=False, na_rep="NULL",float_format='%.8f')


def calc_dye_darboux(xyz, top):
    """Calculate the darboux orientation vector for each dye residue
    in the topology for all frames.

    *Parameters*
            xyz:  *np.array of floats*, shape: (n_ts,n_atoms,3)
                    Atomic positions.
            top: *mdtraj.Topology*
                    Topology of the system
    *Returns*
            t_vects: *list of np.array*, shape: [n_dyes,(n_ts,3)]
                    Components of the dye's orientation in the Darboux frame.
            planar_t_vects: *list of np.array*, shape: [3*n_dyes,(n_ts,3)]
                    Components of the planar vectors in the Darboux frame.
            names: *list of str*
                    Names of the dyes.
    """
    
    #Determine the residues of dyes
    chains = [ch for ch in top.chains]
    
    cy3_res = []
    cy5_res = []
    residues = [res for res in top.residues]
    for r in residues:
        if r.name == "C3N":
            cy3_res.append(r)
        elif r.name == "C5N":
            cy5_res.append(r)

    all_dyes = cy3_res + cy5_res
    dye_names = ["cy3_%s" % i for i in range(len(cy3_res))]+ [
            "cy5_%s" % i for i in range(len(cy5_res))]


    res_5p = [] #List of N+1 residues in the 5' direction
    res_3p = [] #List of N-1 residues in the 3' direction
    res_5p_opp = [] #List of N+1 residues complementary to the 5p
    res_3p_opp = [] #List of N-1 residues complementary to the 3p

    #Loop through dyes, identifying corresponding N+1 and N-1 residues
    for i, dye in enumerate(all_dyes):
        res_5p.append(residues[dye.index-1])
        res_3p.append(residues[dye.index+1])

        #Calculate complementary basepairs
        n_bp = len([r for r in chains[0].residues])
        opp_bp = lambda x: (2*n_bp - x) - 1 #Calculate index of complementary bp
        res_5p_opp.append(residues[opp_bp(res_5p[-1].index)])
        res_3p_opp.append(residues[opp_bp(res_3p[-1].index)])
        
    #Calculate dye orientation vectors in real space
    dye_vects = [calc_dye_vects(xyz,r) for r in all_dyes]

    #t1 vector is the average between the cross-strand sugar axis in the
    #N+1 and N-1 planes, minus its projection onto t3
    #t3 vector is the average normal to the N+1 and N-1 planes, pointing
    #3' to 5'
    t1_vects,t3_vects = calc_t1_t3_vects(xyz,res_5p,res_3p,res_5p_opp,res_3p_opp)

    print("T1 and t3 vects calculated.")
    print("Calculating T2 vect...")

    #t2 = t3 cross t1
    if type(t1_vects) != list:
        t1_vects = [t1_vects]
        t3_vects = [t3_vects]
    t2_vects = [np.cross(t3_vects[i],t1_vects[i]) for i in range(len(t1_vects))]

    renorm_2d = lambda x: x / (np.linalg.norm(x,axis=1)[:,np.newaxis])
    
    t_vects = []
    #Calculate components in the Darboux frame
    for i, dv in enumerate(dye_vects):
        tv = np.zeros((len(dv),3), dtype='float')
        tv[:,0] = np.sum(t1_vects[i]*dv, axis=-1)
        tv[:,1] = np.sum(t2_vects[i]*dv, axis=-1)
        tv[:,2] = np.sum(t3_vects[i]*dv, axis=-1)
        #Renormalize
        tv = renorm_2d(tv)
        t_vects.append(np.copy(tv)) 

    #Calculate the DNA ribbon planar vectors
    p3_planar, p5_planar, p35_planar = calc_planar_vects(xyz,res_5p,res_3p,
            res_5p_opp,res_3p_opp)
    #Add these names:
    new_names = []
    for i, name in enumerate(dye_names):
        new_names += [name + ele for ele in ["_p3_planar","_p5_planar","_p35_planar"]]
    dye_names = dye_names + new_names
    
    #Change basis from xyz to t1,t2,t3 and renormalize
    planar_vects = [p3_planar,p5_planar,p35_planar]
    planar_t_vects = []
    
    #Planar_vects should have shape: [3,n_dyes,(n_ts,3)]
    #We want to recast to [n_dyes,(n_ts,3)]
    flat_planar_vects = []
    n_dyes = len(p3_planar)
    for i in range(n_dyes):
        flat_planar_vects += [p3_planar[i],p5_planar[i],p35_planar[i]]


    #Calculate components in the Darboux frame
    for i, pv in enumerate(flat_planar_vects):
        tv = np.zeros((len(pv),3), dtype='float')
        #Use the t vectors from the corresponding dye with index i//3
        tv[:,0] = np.sum(t1_vects[i//3]*pv, axis=-1)
        tv[:,1] = np.sum(t2_vects[i//3]*pv, axis=-1)
        tv[:,2] = np.sum(t3_vects[i//3]*pv, axis=-1)
        #Renormalize
        tv = renorm_2d(tv)
        planar_t_vects.append(np.copy(tv)) 
        


    return t_vects, planar_t_vects, dye_names
    
def calc_dye_vects(xyz, res):
    """Calculate the real-space orientation vector of the dye res, 5' to 3'

    *Parameters*
            xyz:  *np.array of floats*, shape: (n_ts,n_atoms,3)
                    Atomic positions.
            res: *mdtraj.residue*
                    Residue of the dye
    *Returns*
            dye_vects: *list of np.array*, shape: [(n_ts,3)]
                    Components of the dye's orientation in the xyz frame.
    """
    
    #Determine the residues of dyes
    indole_5p_inds = []
    indole_3p_inds = []

    #Determine name of all atoms in dye
    if res.name == "C5N":
        indole_5p_names = set(['N2', 'C15', 'C14','C13','C12','C17','C16',
                'C18','C19','C20','C30'])
        indole_3p_names = set(['N1','C1','C2','C3','C4','C5','C6','C7','C29',
                'C28','C8'])
    elif res.name == "C3N":
        indole_5p_names = set(['N2', 'C15', 'C14','C13','C12','C17','C16',
                'C18','C19','C20','C30'])
        indole_3p_names = set(['N1','C1','C2','C3','C4','C5','C6','C7','C29',
                'C28','C8'])
    else:
        raise ValueError("Incorrect name of dye residue: %s" % res.name)
    
    indole_5p_inds = np.array([at.index for at in res.atoms if at.name in indole_5p_names])
    indole_3p_inds = np.array([at.index for at in res.atoms if at.name in indole_3p_names])

    #Calculate COMs
    indole_5p_com = np.average(xyz[:,indole_5p_inds],axis=1)
    indole_3p_com = np.average(xyz[:,indole_3p_inds],axis=1)

    dye_vects = indole_3p_com - indole_5p_com
    #Normalize
    #print("Dye_vects shape:",dye_vects.shape)
    dye_vects = np.transpose(np.transpose(dye_vects) / 
            np.linalg.norm(dye_vects,axis=1) )

    return dye_vects

def calc_t1_t3_vects(xyz, res_5p,res_3p,res_5p_opp,res_3p_opp):
    """Calculate the real-space orientation vector t3, the principal tangent.
    t3 vector is the average normal to the N+1 and N-1 planes.

    *Parameters*
            xyz:  *np.array of floats*, shape: (n_ts,n_atoms,3)
                    Atomic positions.
            res_5p: *mdtraj.residue or list of residues*
                    Residue of DNA in the same-strand 5'-ward of the dye.
            res_3p: *mdtraj.residue or list of residues*
                    Residue of DNA in the same-strand 3'-ward of the dye.
            res_5p_opp: *mdtraj.residue or list of residues*
                    Residue of DNA in the opposite-strand complementary
                    to res_5p.
            res_3p_opp: *mdtraj.residue or list of residues*
                    Residue of DNA in the opposite-strand complementary
                    to res_3p.
    *Returns*
            t1_vects: *list of np.array*, shape: [n_dyes,(n_ts,3)]
                    Normalized components of the t1 orientation in the xyz frame.
            t3_vects: *list of np.array*, shape: [n_dyes,(n_ts,3)]
                    Normalized components of the t3 orientation in the xyz frame.
    """
    #If called on a strand with multiple dyes, call recursively
    if type(res_5p) == list:
        t1_vects = []
        t3_vects = []
        for i in range(len(res_5p)):
            t1v, t3v = calc_t1_t3_vects(xyz,res_5p[i],res_3p[i],res_5p_opp[i],
                    res_3p_opp[i])
            t1_vects.append(np.copy(t1v))
            t3_vects.append(np.copy(t3v))
        return t1_vects, t3_vects


    #Determine name of all atoms in dye
    vect_dict = {} #Key: Resname. Value: "startpoint atom name set", "
    #                  endpoint atom name set" 
    
    #Set atom names for every possible residue name
    vect_dict["C5N"] = [set(['N2', 'C15', 'C14','C13','C12','C17','C16',  
                'C18','C19','C20','C30']),set(['N1','C1','C2','C3','C4','C5','C6','C7','C29',
                'C28','C8'])]
    vect_dict["C3N"] = [set(['N2', 'C15', 'C14','C13','C12','C17','C16',  
                'C18','C19','C20','C30']),set(['N1','C1','C2','C3','C4','C5','C6','C7','C29',
                'C28','C8'])]
    
    vect_dict["DG"] = [set(["C4'","O4'","C1'","C2'","C3'"]),set(['N1','C2','N2','N3','C4','C5',
                    'C6','O6','N7','C8','N9'])]
    vect_dict["DC"] = [set(["C4'","O4'","C1'","C2'","C3'"]),set(['N1','C2',
                    'O2','N3','C4','N4','C5','C6'])]
    vect_dict["DA"] = [set(["C4'","O4'","C1'","C2'","C3'"]),set(['N1','C2','N3','C4','C5',
                    'C6','N6','N7','C8','N9'])]
    vect_dict["DT"] = [set(["C4'","O4'","C1'","C2'","C3'"]),set(['N1','C2',
                    'O2','N3','C4','O4','C5','C6','C7'])]

    #Replace residues with atom indices
    res_labels = ["res_5p","res_3p","res_5p_opp","res_3p_opp"]
    all_res = [res_5p,res_3p,res_5p_opp,res_3p_opp]
    all_start_inds = [] #Shape: 4, n_atoms
    all_end_inds = [] #Shape: 4,n_atoms
    for i, res in enumerate(all_res):
        all_start_inds.append(np.array([at.index for at in res.atoms if at.name in vect_dict[res.name][0]]))
        all_end_inds.append(np.array([at.index for at in res.atoms if at.name in vect_dict[res.name][1]]))


    #Calculate COMs
    all_start_coms = np.array([np.average(xyz[:,start_inds],axis=1) for start_inds in all_start_inds])
    all_end_coms = np.array([np.average(xyz[:,end_inds],axis=1) for end_inds in all_end_inds])

    #Calculate explicit orientation vectors
    dye_vects = all_end_coms - all_start_coms
    #print("Dye_vects shape, should be (4,n_ts,3):",dye_vects.shape)


    #Split into individual variables
    vect_5p = dye_vects[0]
    vect_3p = dye_vects[1]
    vect_5p_opp = dye_vects[2]
    vect_3p_opp = dye_vects[3]

    renorm_2d = lambda x: x / (np.linalg.norm(x,axis=1)[:,np.newaxis])
    renorm_3d = lambda x: x / (np.linalg.norm(x,axis=2)[:,:,np.newaxis])

    #print("vect_5p shape:",vect_5p.shape)

    #Normalize each vector
    vect_5p = renorm_2d(vect_5p)
    vect_3p = renorm_2d(vect_3p)
    vect_5p_opp = renorm_2d(vect_5p_opp)
    vect_3p_opp = renorm_2d(vect_3p_opp)

    #These vectors should be roughly negative of each other. We can
    #Take an appropriate average by subtracting them and renormalizing!

    #print("Finished renormalizing. Subtracting...")
    t1_5p = vect_5p - vect_5p_opp
    t1_5p = renorm_2d(t1_5p)
    t1_3p = vect_3p - vect_3p_opp
    t1_3p = renorm_2d(t1_3p)

    t1_vects = np.average((t1_5p,t1_3p),axis=0)
    t1_vects = renorm_2d(t1_vects)
    
    #calculate t3 vectors as the average same- and opposite strand vectors
    #Pointing between sugars.
    t3_approx = all_start_coms[0] - all_start_coms[1]
    t3_opp_approx = all_start_coms[2] - all_start_coms[3]
    t3_vects = renorm_2d(t3_approx + t3_opp_approx)

    #print("Performing cross product...")
    #We do the same thing with a crossproduct to calculate t3.
    #t3_5p = np.cross(vect_5p,vect_5p_opp)
    #t3_3p = np.cross(vect_3p,vect_3p_opp)
    #These should point in the same direction. Renormalize, average, then
    #renormalize.
    #t3_5p = renorm_2d(t3_5p)
    #t3_3p = renorm_2d(t3_3p)
    
    #print("Removing projection...")
    #t3_vects = renorm_2d(t3_5p + t3_3p)
    #Remove projection of t3 onto t1
    t3_vects = t3_vects - t1_vects * np.sum(t3_vects * t1_vects,axis=-1)[:,np.newaxis]
    t3_vects = renorm_2d(t3_vects)
    #print("returning.")
    return t1_vects, t3_vects


def calc_planar_vects(xyz, res_5p,res_3p,res_5p_opp,res_3p_opp):
    """Calculate the real-space vectors p3, p5, and p35. These point from the
    4-residue COM to the 3' phosphate, 5' phosphate, and 3' --> 5' phosphates
    respectively.

    *Parameters*
            xyz:  *np.array of floats*, shape: (n_ts,n_atoms,3)
                    Atomic positions.
            res_5p: *mdtraj.residue or list of residues*
                    Residue of DNA in the same-strand 5'-ward of the dye.
            res_3p: *mdtraj.residue or list of residues*
                    Residue of DNA in the same-strand 3'-ward of the dye.
            res_5p_opp: *mdtraj.residue or list of residues*
                    Residue of DNA in the opposite-strand complementary
                    to res_5p.
            res_3p_opp: *mdtraj.residue or list of residues*
                    Residue of DNA in the opposite-strand complementary
                    to res_3p.
    *Returns*
            p3_vects: *list of np.array*, shape: [n_dyes,(n_ts,3)]
                    Normalized components of the t1 orientation in the xyz frame.
            p5_vects: *list of np.array*, shape: [n_dyes,(n_ts,3)]
                    Normalized components of the t3 orientation in the xyz frame.
            p35_vects: *list of np.array*, shape: [n_dyes,(n_ts,3)]
                    Normalized components of the t3 orientation in the xyz frame.
    """
    #If called on a strand with multiple dyes, call recursively
    if type(res_5p) == list:
        p3_vects = []
        p5_vects = []
        p35_vects = []
        for i in range(len(res_5p)):
            p3v, p5v, p35v = calc_planar_vects(xyz,res_5p[i],res_3p[i],res_5p_opp[i],
                    res_3p_opp[i])
            p3_vects.append(np.copy(p3v))
            p5_vects.append(np.copy(p5v))
            p35_vects.append(np.copy(p35v))
        return p3_vects, p5_vects, p35_vects

    #Replace residues with atom indices
    res_labels = ["res_5p","res_3p","res_5p_opp","res_3p_opp"]
    all_res = [res_5p,res_3p,res_5p_opp,res_3p_opp]
    all_inds = [] #Shape: 4, n_atoms
    all_phos_inds = []

    for i, res in enumerate(all_res):
        all_inds.append(np.array([at.index for at in res.atoms]))
        all_phos_inds.append(np.array([at.index for at in res.atoms if at.name == "P"]))

    #Calculate COMs, first within all 4 residues, then across all 4
    res_coms = np.array([np.average(xyz[:,inds],axis=1) for inds in all_inds])
    phos_pos = np.array([xyz[:,ind] for ind in all_phos_inds])
    com = np.average(res_coms,axis=0)
    #print("Planar COM:",com)
    
    #print("phos_pos[1] shape:",phos_pos[1].shape)
    #Calculate vectors and flatten appropriately
    p3_vects = phos_pos[1][:,0] - com
    p5_vects = phos_pos[0][:,0] - com
    p35_vects = phos_pos[0][:,0] - phos_pos[1][:,0]


    renorm_2d = lambda x: x / (np.linalg.norm(x,axis=1)[:,np.newaxis])
    renorm_3d = lambda x: x / (np.linalg.norm(x,axis=2)[:,:,np.newaxis])

    #print("p3_vects shape:",p3_vects.shape)

    #Normalize each vector
    p3_vects = renorm_2d(p3_vects)
    p5_vects = renorm_2d(p5_vects)
    p35_vects = renorm_2d(p35_vects)

    return p3_vects, p5_vects, p35_vects

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-d",
        "--dir",
        type=str,
        default=".",
        help="Path to simulation Processing directory.",
    )
    # parser.add_argument('-f','--file',type=str,help="Path to lammpstrj to analyze.")
    args = parser.parse_args()

    main()

