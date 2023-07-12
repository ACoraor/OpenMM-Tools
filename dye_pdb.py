#Module construct_dna
#Created by Aria Coraor
#Written 6/6/23

import copy
import numpy as np
import mdtraj as md
import argparse
import warnings
import os

MASS = {"C":12.0107, "N": 14.0067, "O": 15.999, "P": 30.973762}

def main():
    """Modify the existing Cy3 and Cy5 pdbs to have a phosphate linker,
    for direct inclusion in openMM structure. Phosphate group is coming
    from the 5' side.

    Truncate the lysine linkage and add an additional alkyl linkage.
    """
    
    #Load input pdb
    cy3 = md.load_pdb(args.donor)
    cy3_xyz = cy3.xyz[0]
    cy3_top = cy3.top

    cy5 = md.load_pdb(args.acceptor)
    cy5_xyz = cy5.xyz[0]
    cy5_top = cy5.top

    print("PDBs loaded. Modifying...")

    #Remove lysine tail
    cy3_nt_xyz, cy3_nt_top = remove_tail(cy3_xyz,cy3_top)
    cy5_nt_xyz, cy5_nt_top = remove_tail(cy5_xyz,cy5_top)

    #Add phosphate and new alkyl linker
    cy3_built_xyz, cy3_built_top = add_links(cy3_nt_xyz,cy3_nt_top)
    cy5_built_xyz, cy5_built_top = add_links(cy5_nt_xyz,cy5_nt_top)

    #Save new pdbs
    pdb = md.formats.PDBTrajectoryFile('cy3_dna.pdb',mode='w')
    pdb.write(10.*cy3_built_xyz,cy3_built_top)
    
    pdb = md.formats.PDBTrajectoryFile('cy5_dna.pdb',mode='w')
    pdb.write(10.*cy5_built_xyz,cy5_built_top)

    print("DNA dye pdbs written! Exiting...")


def remove_tail(xyz,top):
    """Remove lysine tail from pdb.

    Parameters:
        xyz: *np.array*, shape: (n_atoms,3)
            Positions of particles in the raw AMBER dye.
        top: *mdtraj.topology*
            Topology of the amber dye.
    Returns:
        nt_xyz: *np.array*, shape: (n_atoms,3)
            Positions of particles in the new AMBER dye without tail atoms.
        nt_top: *mdtraj.topology*
            A new Topology of the amber dye without tail atoms.
    """
    #Remove lysine
    res = [r for r in top.residues]

    rm_inds = []
    #Remove all lysine atoms
    for at in res[1].atoms:
        rm_inds.append(at.index)

    rm_names = ["C25","C26","C99","O91", "H27","H28","H29","H30"]
    for at in res[0].atoms:
        if at.name in rm_names:
            rm_inds.append(at.index)

    all_ats = [at for at in top.atoms]
    keep_inds = [n for n in range(len(all_ats)) if n not in set(rm_inds)]
    keep_inds = np.array(keep_inds,dtype=int)

    nt_xyz = xyz[keep_inds]
    nt_top = top.subset(keep_inds)

    return nt_xyz,nt_top


def add_links(xyz,top):
    """Add initial phosphate and terminal alkyl chain.

    Parameters:
        xyz: *np.array*, shape: (n_atoms,3)
            Positions of particles in the AMBER dye.
        top: *mdtraj.topology*
            Topology of the amber dye.
    Returns:
        built_xyz: *np.array*, shape: (n_atoms,3)
            Positions of particles in the new AMBER dye with linkers.
        built_top: *mdtraj.topology*
            A new Topology of the amber dye with linkers.
    """
    res1 = [r for r in top.residues][0]
    
    #Construct phosphate
    sym = "P"
    top.add_atom(sym,md.element.Element.getBySymbol(sym),residue=res1)

    sym = "O"
    top.add_atom("O5'",md.element.Element.getBySymbol(sym),residue=res1)
    sym = "O"
    top.add_atom("OP1",md.element.Element.getBySymbol(sym),residue=res1)
    sym = "O"
    top.add_atom("OP2",md.element.Element.getBySymbol(sym),residue=res1)

    #Construct phosphate positions
    #Treat the vector from C to O as the same as the prior leg
    c0_i = top.select('name C24')[0]
    c1_i = top.select('name C23')[0]
    c2_i = top.select('name C22')[0]
    #print("Indices:",c0_i,c1_i,c2_i)

    #Works for 
    C0 = xyz[c0_i]
    C1 = xyz[c1_i]
    C2 = xyz[c2_i]
    
    d1 = C1 - C2
    d2 = C0 - C1
    o5 = C0 + d1 #Oxygen 1 position
    P = o5 + d2 #Phosphorus core position

    d3 = d2 - d1 # Bonus vector to place OP3
    h1_i = top.select('name H25')[0]
    h2_i = top.select('name H26')[0]
    H1 = xyz[h1_i]
    H2 = xyz[h2_i]
    dh3 = H1 - C0
    dh4 = H2 - C0
    op1 = P + 1.5*dh3
    op2 = P + 1.5*dh4
    


    new_xyz = np.stack((P,o5,op1,op2))
    xyz = np.concatenate((xyz,new_xyz))


    #Add alkyl linker
    #Remove hydrogens from other side
    c2p_i = top.select('name C21')[0]
    h0p_i = top.select('name H18')[0]
    h1p_i = top.select('name H19')[0]
    h2p_i = top.select('name H20')[0]

    #Move H18, H19, prep to delete H20
    c2_h1 = top.select('name H21')[0]
    c2_h2 = top.select('name H22')[0]
    d1p = xyz[c2_h1] - C2
    d2p = xyz[c2_h2] - C2

    xyz[h0p_i] = xyz[c2p_i] + d1p
    xyz[h1p_i] = xyz[c2p_i] + d2p




    #Add carbon linkers
    sym = "C"
    top.add_atom("CG",md.element.Element.getBySymbol(sym),residue=res1)
    sym = "H"
    top.add_atom("HG",md.element.Element.getBySymbol(sym),residue=res1)
    sym = "H"
    top.add_atom("HI",md.element.Element.getBySymbol(sym),residue=res1)
    sym = "C"
    top.add_atom("CI",md.element.Element.getBySymbol(sym),residue=res1)
    sym = "H"
    top.add_atom("HN",md.element.Element.getBySymbol(sym),residue=res1)
    sym = "H"
    top.add_atom("HN",md.element.Element.getBySymbol(sym),residue=res1)
    #All carbons added
    C1p = xyz[c2p_i]+d1
    C0p = C1p + d2

    #Add hydrogens
    h1p1 = C1p - d1p
    h2p1 = C1p - d2p
    h1p2 = C0p + d1p
    h2p2 = C0p + d2p

    new_xyz = np.stack((C1p,h1p1,h2p1,C0p,h1p2,h2p2))
    xyz = np.concatenate((xyz,new_xyz))


    xyz = np.concatenate((xyz[:h2p_i],xyz[h2p_i+1:]))
    ats = [at for at in top.atoms]
    inds = [i for i in range(len(ats)) if i != h2p_i]
    top = top.subset(inds)

    return xyz, top


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
    parser.add_argument('-d','--donor', type=str, 
        default='C3N_L1R.pdb', help='Path to basic Cy3 pdb.')
    parser.add_argument('-a','--acceptor', type=str, 
        default='C5N_L1R.pdb', help='Path to basic Cy5 pdb.')
    args = parser.parse_args()

    main()

