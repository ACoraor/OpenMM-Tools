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
    top.add_atom("OP1",md.element.Element.getBySymbol(sym),residue=res1)
    sym = "O"
    top.add_atom("OP2",md.element.Element.getBySymbol(sym),residue=res1)

    #Construct phosphate positions

    #Add alkyl linker
    #Remove hydrogen


    print("STUB!")

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

