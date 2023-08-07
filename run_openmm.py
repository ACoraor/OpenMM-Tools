#Module run_openmm
#Created by Aria Coraor 6/2/23


import numpy as np
import argparse
#Following openMM docs
from openmm import *
import sys
import openmm.app as app
import openmm.unit as unit
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pysages
from openmm.vec3 import Vec3

def main():
    """Run an OpenMM simulation on the FRET constructs.
    """
    #Set defaults
    #pdb = "input.pdb"
    #amber = "DNA.bsc1.xml"
    #water = "spce_standard.xml"


    #Load system
    pdb, ff = load_system(args.pdb,args.amber,args.water,args.gaff)

    #Create GAFF templates for OpenMM
    #pdb,ff = create_dye_templates(pdb,ff)    

    #Set up simulation environment
    simulation = generate_simulation(pdb,ff)

    #Run the system
    results = run_simulation(simulation,args.timesteps)

    print("Simulation complete.")
    print("results:",repr(results))

def load_system(pdb,amber,water,gaff):
    """Load and parse pdb, and amber/water parameter files.

    Parameters:
	pdb: *str*
	    Path to input pdb.
	amber: *str*
	    Path to amber parameter xml.
	water: *str*
                Path to SPCE xml.
	gaff: *str*
	    Path to GAFF xml.
    Returns:
	pdb: *openmm.app.PDBFile*
	    Loaded DNA construct pdb
	ff: *openmm.app.Forcefield*
	    Loaded forcefield file.
    """
    pdb = app.PDBFile(pdb)
    

    ff = app.ForceField(amber,water,gaff)
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
    rs = [res for res in topology.residues]
    dye_names = ["C3N","C5N"]
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
                print("dist:",dist)
                if np.linalg.norm(dist) < bond_dist:
                    topology.addBond(dye_ats[i],dye_ats[j])
   
    

    #bonds = [b for b in topology.bonds()]
    #print("All bonds:",bonds)
    
    #res = [r for r in topology.residues()]
    #res[0].name = "DG5" 
    #Add spce solvent molecules
    modeller = app.modeller.Modeller(topology,pos)
    
    bounds = Vec3(10.0,10.0,32.0)*unit.nanometer #Change
    #dna_ff = app.ForceField(args.amber) 
    modeller.loadHydrogenDefinitions(args.amber)
    print("Adding hydrogens...")
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

    simulation = app.Simulation(topology, system, integrator)
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
    #simulation.reporters.append(app.PDBReporter('output.pdb',50))
	    

    #Run simulation
    print("Running!")
	
    for i in range(n_ts//traj_int):
        simulation.step(traj_int)
        fn = "state_%s.xml" % i
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
    parser.add_argument('-g','--gaff',type=str,default = "gaff-aec.xml",
	help="Path to GAFF parameter file.")
    parser.add_argument('-t','--timesteps',type=int,default = 10000000,
	help="Number of timesteps for simulation.")
	
    args = parser.parse_args()
    main()
