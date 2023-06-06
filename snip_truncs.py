#Module snip_truncs
#Created by Aria Coraor
#Written 4/11/22

import copy
import numpy as np
import mdtraj as md
import argparse
import warnings
import os

MASS = {"C":12.0107, "N": 14.0067, "O": 15.999, "P": 30.973762}
ELE_MAP = {"ARG":"Ar","THR":"Th","LYS":'VS','GLN':'VS','ALA':'Al','SER':'Se',
		"PRO":"Pr","VAL":"V","ASP":"As","ILE":"I"}
TYPE_MAP = {'A':0,
'ALA'   :1,
'ARG'   :2,
'ASN'   :3,
'ASP'   :4,
'C'	 :5,
'CYS'   :6,
'G'	 :7,
'GLN'   :8,
'GLU'   :9,
'GLY'   :10,
'HIS'   :11,
'ILE'   :12,
'LEU'   :13,
'LYS'   :14,
'MET'   :15,
'P'	 :16,
'PHE'   :17,
'PRO'   :18,
'S'	 :19,
'SER'   :20,
'T'	 :21,
'THR'   :22,
'TYR'   :23,
'VAL'   :24}


def main():
	"""Change 3' and 5' truncation pattern to match parmBsc1."""

	#Read the raw text into a list. Split into columns, and strip each element.
	xyz, top = read_raw_pdb(args.file)

	#Identify removable DNA inds
	chs = [ch for ch in top.chains]
	#Remove phosphate on 5' ends
	ats = [at for at in top.atoms]
	ch1_at = [at for at in chs[0].atoms]
	ch2_at = [at for at in chs[1].atoms]
	rm_ats = []
	for i in range(3):
		rm_ats.append(ch1_at[i].index)
		rm_ats.append(ch2_at[i].index)
	#rm_ats = [ele for ele in range(3)] + [ele for ele in range(len(ats)-3,len(ats))]
	keep_inds = [i for i in range(len(ats)) if i not in set(rm_ats)]
	
	new_top = top.subset(keep_inds)
	new_xyz = xyz[np.array(keep_inds,dtype=int)]
	
	

	tf = md.formats.PDBTrajectoryFile(args.output,mode='w')
	tf.write(10.*new_xyz,new_top)
	print("Snipped.")		
	

def read_raw_pdb(fn):
	""" Read the pdb file, parse out the all-atom representation.
	Clean it, turning relevant fields into floats/ints/etc.

	*Parameters*
			fn: *str*
					Path to pdb file.
	*Returns*
			xyz: *3d np.array*, <N_timesteps>, n_atoms, xyz
					xyz position of each particle.
			top: *mdtraj.Topology*
					Nucleosome topology.
	"""
	a = md.load(fn)
	top = a.topology
	return a.xyz[0], top

def make_MTpdb_rep(xyx,top,at_xyz,at_top,outfile):
	""" Create a MTpdb representation (3/13/23) of the atomistic pdb.

	*Parameters*
		xyz: *3d np.array*, <N_timesteps>, n_atoms, xyz
				xyz position of each particle.
		top: *mdtraj.Topology*
				Nucleosome topology.
		at_xyz: *3d np.array*, <N_timesteps>, n_atoms, xyz
				xyz position of each atom.
		at_top: *mdtraj.Topology*
				Atomistic nucleosome topology.
		outfile: *str*
			Path to write to file.
	"""

	cg_aas = [ch for ch in top.chains][2:]

	atomistic_chains = [ch for ch in at_top.chains]
	atomistic_aas = atomistic_chains[2:]
	
	#unpack all residues
	atomistic_residues = [r for ch in atomistic_aas for r in ch.residues]
	cg_residues = [r for ch in cg_aas for r in ch.residues]

	if not len(atomistic_residues) == len(cg_residues):
		raise ValueError("Different numbers of atomistic vs. cg residues: %s vs. %s" % (
			len(atomistic_residues),len(cg_residues)))

	#Parse and write to file
	with open(outfile,'w') as f:
		#Looping over all residues
		for i, (aa_r,cg_r) in enumerate(zip(atomistic_residues,cg_residues)):
			atomtype = str(aa_r.name)
			if str(aa_r.name) != str(cg_r.name):
				raise ValueError("Name mismatch between aa and cg: " +
					"%s vs. %s" % (aa_r.name,cg_r.name))
			f.write("# %s %s \n" % (atomtype,str(i) ))
			#Create atom lists
			aa_ats = [at for at in aa_r.atoms]
			cg_ats = [at for at in cg_r.atoms]
			#Atomnumber line
			atom_number = TYPE_MAP[cg_r.name]
			f.write("%s\n" % atom_number)
			n = len(aa_ats)
			#n line
			f.write("%s\n" % n)
			
			atom_ids = [at.name for at in aa_ats]
			inds = [at.index for at in aa_ats]
			atom_pos = at_xyz[inds]
			#print("atom_pos.shape:",atom_pos.shape)
			for atomid, pos in zip(atom_ids,atom_pos):
				f.write("%s %s %s %s\n" % (atomid, pos[0], pos[1],
						pos[2]) )
	return	


	

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-f','--file', type=str, default='input.pdb',
		help='PDB file to snip.')
	parser.add_argument('-o','--output', type=str, default='output_pdb.txt',
		help='Path to output pdb file.')
	args = parser.parse_args()

	main()

