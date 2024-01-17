#Module run_openmm
#Created by Aria Coraor 6/2/23


import numpy as np
import argparse
#Following openMM docs
from openmm import *
import openmm
import sys
import openmm.app as app
import openmm.unit as unit
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mdtraj as md
import pysages
from openmm.vec3 import Vec3

def main():
    """Run an OpenMM simulation on the FRET constructs.
    """
    #Set defaults
    #pdb = "input.pdb"
    #amber = "DNA.bsc1.xml"
    #water = "spce_standard.xml"
    print("Starting system.")

    #Load system
    pdb, ff = load_system(args.pdb,args.amber,args.water,args.dyes,args.gaff)

    #Create GAFF templates for OpenMM
    #pdb,ff = create_dye_templates(pdb,ff)    

    #Set up simulation environment
    simulation = generate_simulation(pdb,ff)
    print("Simulation generated. Running...")
    #Run the system
    results = run_simulation(simulation,args.timesteps)

    print("Simulation complete.")
    print("results:",repr(results))

def load_system(pdb,amber,water,dyes,gaff):
    """Load and parse pdb, and amber/water parameter files.

    Parameters:
	pdb: *str*
	    Path to input pdb.
	amber: *str*
	    Path to amber parameter xml.
	water: *str*
                Path to SPCE xml.
	dyes: *str*
                Path to gaff-aec xml.
	gaff: *str*
	    Path to GAFF xml.
    Returns:
	pdb: *openmm.app.PDBFile*
	    Loaded DNA construct pdb
	ff: *openmm.app.Forcefield*
	    Loaded forcefield file.
    """
    pdb = app.PDBFile(pdb)
    

    ff = app.ForceField(amber,water,dyes,gaff)
    return pdb,ff


def generate_simulation(pdb, ff):
    """Generate an openmm simulation.

    Parameters:
	pdb: *openmm.app.PDBFile*
	    Loaded DNA construct pdb
	ff: *openmm.app.Forcefield*
	    Loaded forcefield file.
    Returns:
	simulation: *openmm.app.simulation*
	    Preset simulation with langevin integrator.
	pdb: *openmm.app.PDBFile*
	    Loaded DNA construct pdb
	ff: *openmm.app.Forcefield*
	    Loaded forcefield file.
    """
    #Set upp system constants
    T = 298.15 * unit.kelvin
    dt = 2.0 * unit.femtoseconds
    cutoff_distance = 1.0 * unit.nanometer
    topology = pdb.topology
    
    pos = pdb.positions

    #pos *= 0.1 # Convert from angstroms to nanometers
    
    #Create standard bonds
    topology.loadBondDefinitions(args.amber)
    topology.createStandardBonds()
    
    #Manually add bonds for C3N, C5N residues
    rs = [res for res in topology.residues()]

    dye_names = ["C3N","C5N"]
    dye_res = [r for r in rs if r.name in dye_names]
    print("dye_res:",dye_res)

    #Manually check: if any bonds exist between non-dye residues and dye res,
    #Delete them.
    all_dye_atoms = []
    for dye_r in dye_res:
        all_dye_atoms = all_dye_atoms + [at for at in dye_r.atoms()]
    all_dye_atoms = set(all_dye_atoms)
    
    bond_rms = []
    for bond in topology.bonds():
        at1_in_dye = bond[0] in all_dye_atoms
        at2_in_dye = bond[1] in all_dye_atoms
        if at1_in_dye != at2_in_dye:
            bond_rms.append(bond)

    #Rebuild topology without bond_rms
    copy_top = app.topology.Topology()
    for chain in topology.chains():
        copy_top.addChain(id=chain.id)
        curr_ch = [ch for ch in copy_top.chains()][-1]
        for res in chain.residues():
            copy_top.addResidue(name=res.name,chain=curr_ch)
            curr_res = [r for r in copy_top.residues()][-1]
            for atom in res.atoms():
                copy_top.addAtom(name=atom.name,element=atom.element,
                    residue = curr_res)
    bond_rms = set(bond_rms)
    copy_ats = [at for at in copy_top.atoms()]
    for bond in topology.bonds():
        if bond not in bond_rms:
            at1 = copy_ats[bond[0].index]
            at2 = copy_ats[bond[1].index]
            copy_top.addBond(at1,at2)

    print("Old topology:",topology)
    print("New topology:",copy_top)
    print("bond_rms:",list(bond_rms))
    topology = copy_top

    rs = [res for res in topology.residues()]
    dye_res = [r for r in rs if r.name in dye_names]

    bond_dist = 0.18 #nanometers
    for dye_r in dye_res:
        dye_ats = [at for at in dye_r.atoms()]

        for i in range(len(dye_ats)):
            for j in range(i+1,len(dye_ats)):
                index_i = dye_ats[i].index
                index_j = dye_ats[j].index
                #Manually cast Vec3 to numpy array
                dist = (np.array([float(pos[index_i].x),float(pos[index_i].y),
                    float(pos[index_i].z)]) - np.array([float(pos[index_j].x),
                    float(pos[index_j].y),float(pos[index_j].z)]))
                if np.linalg.norm(dist) < bond_dist:
                    topology.addBond(dye_ats[i],dye_ats[j])
        #Manually add bonds to the prior and later residue atoms
        dye_phos = [at for at in dye_ats if at.name=="P"][0]
        dye_O3 = [at for at in dye_ats if at.name=="O3'"][0]
        #Bond dye O3' to phosphate from earlier residue
        prev_res = rs[dye_r.index-1]
        prev_res_O3 = [at for at in prev_res.atoms() if at.name=="O3'"][0]
        topology.addBond(dye_phos,prev_res_O3)
        #Bond dye phosphate to O3' from next residue
        next_res = rs[dye_r.index+1]
        next_res_phos = [at for at in next_res.atoms() if at.name=="P"][0]
        topology.addBond(dye_O3,next_res_phos)
   

    
    #res = [r for r in topology.residues()]
    #res[0].name = "DG5" 
    #Add spce solvent molecules
    modeller = app.modeller.Modeller(topology,pos)
    
    bounds = Vec3(5.0,5.0,32.0)*unit.nanometer #Change
    #dna_ff = app.ForceField(args.amber) 
    modeller.loadHydrogenDefinitions(args.amber)
    #print("Adding hydrogens...")
    modeller.addHydrogens(ff)
    print("Adding water...")
    modeller.addSolvent(ff,boxSize=bounds,model='spce')
    print("Creating system.")
    
    #Write true initial state
    with open("init.pdb",'w') as f:
       app.PDBFile.writeFile(modeller.topology,modeller.positions,f) 

    #Create system
    system = ff.createSystem(
        modeller.topology, constraints = app.HBonds, 
        nonbondedMethod = app.PME, nonbondedCutoff = cutoff_distance)
   
    #NPT at 1bar
    #system.addForce(MonteCarloAnisotropicBarostat((1,1,1)*unit.bar,T)) 

    #Create Simulation

    integrator = openmm.LangevinIntegrator(T, 1 / unit.picosecond, dt)

    simulation = app.Simulation(modeller.topology, system, integrator)
    simulation.context.setPositions(modeller.positions)

    #Set PBC ?
    #a = Vec3(10.0,0.0,0.0)
    #b = Vec3(0.0,10.0,0.0)
    #c = Vec3(0.0,0.0,32.0)
    #simulation.context.setPeroidicBoxVectors

    return simulation 

def run_simulation(simulation,n_ts=10000):
    """Run simulation given current setup.

    Parameters:
	simulation: *openmm.app.simulation*
	    Preset simulation with langevin integrator.
	n_ts: *int*
	    Number of timesteps to run. Default=10,000
    Returns:
	simulation: *openmm.app.simulation*
	    Preset simulation with langevin integrator.
    """
    #Minimize initial condition
    print("Minimizing initial energy")
    #simulation.minimizeEnergy(maxIterations=1000)
    simulation.minimizeEnergy()
    print("Energy minimized.")

    with open("minstate.xml",'w') as f:
        simulation.saveState(f)

    traj_int = 10000

    #Add reporters
    simulation.reporters.append(app.StateDataReporter(sys.stdout,5000,
        step=True,potentialEnergy=True,temperature=True))
    
    relevant_inds = np.arange(4205) # Manually counted)
    simulation.reporters.append(md.reporters.HDF5Reporter('output.h5',5000,
        atomSubset=relevant_inds))
    #simulation.reporters.append(app.PDBReporter('output.pdb',5000))



    #Run simulation
    print("Running!")

    
    #Old style running:
    for i in range(n_ts//traj_int):
        simulation.step(traj_int)
        pdb_fn = 'structure_%s.pdb' % i
        #app.PDBFile.writeFile(simulation.topology,simulation.getpositions,pdb_fn) 
        #print("simulation.topology:",simulation.topology())
        #print("simulation positions:",simulation.getPositions())
        
        fn = "restart_%s.xml" % (i % 2)
        with open(fn,'w') as f:
            simulation.saveState(f)
    
    print("Finished!")
    return simulation

def foo():
    """Do foo.

    Parameters:
	bar: *np.array*
	    Biz.
    """
    pass


if __name__ == '__main__':
    #gen_probs()
    #gen_dih()
    parser = argparse.ArgumentParser()
    parser.add_argument('-p','--pdb',type=str,default = "input.pdb",
	help="Path to solvent-free DNA construct PDB.")
    parser.add_argument('-a','--amber',type=str,default = "DNA.bsc1.xml",
	help="Path to Amber parameter file.")
    parser.add_argument('-w','--water',type=str,default = "spce_standard.xml",
	help="Path to Water/SPCE parameter file.")
    parser.add_argument('-d','--dyes',type=str,default = "gaff-aec.xml",
	help="Path to dyes parameter file.")
    parser.add_argument('-g','--gaff',type=str,default = "gaff.xml",
	help="Path to GAFF parameter file.")
    parser.add_argument('-t','--timesteps',type=int,default = 10000000,
	help="Number of timesteps for simulation.")
	
    args = parser.parse_args()
    main()
