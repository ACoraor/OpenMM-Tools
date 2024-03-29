#Module run_openmm
#Created by Aria Coraor 6/2/23


import numpy as np
import argparse
#Following openMM docs
from openmm import *
import os
import openmm
import sys
import openmm.app as app
import openmm.unit as unit
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mdtraj as md
from shutil import copy2
from openmm.vec3 import Vec3
import pickle

#Import pysages things
import pysages
from pysages.grids import Grid
from pysages.colvars import DihedralAngle
from pysages.methods import SpectralABF,ABF


def main():
    """Run an OpenMM simulation on the FRET constructs.
    """
    #Set defaults
    #pdb = "input.pdb"
    #amber = "DNA.bsc1.xml"
    #water = "spce_standard.xml"
    print("Starting system.")

    
    #Load system
    if args.restart:
        pdbfn = "init.pdb"
    else:
        pdbfn = args.pdb
    pdb, ff = load_system(pdbfn,args.amber,args.water,args.dyes,args.gaff)

    #Create GAFF templates for OpenMM
    #pdb,ff = create_dye_templates(pdb,ff)    

    #Set up simulation environment
    if not args.restart:
        simulation = generate_simulation(pdb,ff)
    else:
        T = 298.15 * unit.kelvin
        dt = 2.0 * unit.femtoseconds
        cutoff_distance = 1.0 * unit.nanometer
        system = ff.createSystem(
            pdb.topology, constraints = app.HBonds, 
            nonbondedMethod = app.PME, nonbondedCutoff = cutoff_distance)
        integrator = openmm.LangevinIntegrator(T, 1 / unit.picosecond, dt)
        simulation = app.Simulation(pdb.topology, system, integrator)
    print("Simulation generated.")

    #Add reporters
    simulation = add_reporters(simulation)

    #Run the system
    print("Running simulation.")
    results = run_simulation(simulation,args.timesteps,args.restart,args.enhanced)

    print("Simulation complete.")

    print("Plotting data.")
    plot_data()

def autogen_simulation():
    """Automatically generate the entire simulation for pysages.
    
    Returns:
	simulation: *openmm.app.simulation*
	    Preset simulation with langevin integrator.
    """
    #pdb,ff = load_system('input.pdb','DNA.bsc1.xml','spce_standard.xml',
    #    'gaff-isolated.xml','gaff.xml')
    pdb,ff = load_system('input.pdb','DNA.bsc1.xml','spce_standard.xml',
        'gaff-aec.xml','gaff.xml')
    try:
        simulation = generate_simulation(pdb,ff,is_restart=True)
    except:
        print("Failed to restart.")
        simulation = generate_simulation(pdb,ff)
    simulation.minimizeEnergy()
    add_reporters(simulation)
    return simulation

def add_reporters(simulation,is_restart=False):
    """Add structure and thermodynamic mdtraj reporters to simulation.

    Parameters:
	simulation: *openmm.app.simulation*
	    Preset simulation with langevin integrator.
	is_restart: *bool*
	    Is this simulation a restart
    Returns:
	simulation: *openmm.app.simulation*
	    Preset simulation with langevin integrator.
    """
    relevant_inds = np.arange(3171) # Manually counted)

    #Set up reporters
    if os.path.isfile('output.h5'):        
        if is_restart:
            print("Concatenating outputs...")
            if not os.path.isfile('full_output.h5'):
                copy2('output.h5','full_output.h5')
            else:
                out_data = md.load_hdf5('output.h5')
                full_data = md.load_hdf5('full_output.h5')
                total_xyz = np.concatenate((full_data._xyz,out_data._xyz),axis=0)
                new_traj = md.Trajectory(total_xyz,full_data._topology)
                new_traj.save_hdf5('full_output.h5',force_overwrite=True)
        
        #All existing save-worthy data is now saved at full_output.h5
        os.remove('output.h5')
    odf = md.reporters.HDF5Reporter('output.h5',5000,atomSubset=relevant_inds)
    simulation.reporters.append(odf)

    #Use try/finally to handle datafile closure
    #try:
    state_data_file = open('state_data.txt','a')
    simulation.reporters.append(app.StateDataReporter(state_data_file,5000,
            step=True,potentialEnergy=True,temperature=True))
    
    return simulation

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


def generate_simulation(pdb, ff,is_restart=False):
    """Generate an openmm simulation.

    Parameters:
	pdb: *openmm.app.PDBFile*
	    Loaded DNA construct pdb
	ff: *openmm.app.Forcefield*
	    Loaded forcefield file.
	is_restart: *bool*
	    Is this simulation a restart
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
    amber = "DNA.bsc1.xml"
    topology.loadBondDefinitions(amber)#args.amber)
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
    
    bounds = Vec3(32.0,5.0,5.0)*unit.nanometer #Change
    #dna_ff = app.ForceField(args.amber) 
    modeller.loadHydrogenDefinitions(amber)
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
    system.addForce(MonteCarloAnisotropicBarostat(Vec3(1,1,1)*unit.bar,T)) 

    #Create Simulation

    integrator = openmm.LangevinIntegrator(T, 1 / unit.picosecond, dt)
    
    #STUB
    is_restart=False
    if is_restart:
        try:
            simulation = app.Simulation(modeller.topology, system, integrator,
                state='restart_1.xml')
        except:
            simulation = app.Simulation(modeller.topology, system, integrator,
                state='restart_0.xml')
    else:
        simulation = app.Simulation(modeller.topology, system, integrator)
        simulation.context.setPositions(modeller.positions)

    #Set PBC ?
    #a = Vec3(10.0,0.0,0.0)
    #b = Vec3(0.0,10.0,0.0)
    #c = Vec3(0.0,0.0,32.0)
    #simulation.context.setPeroidicBoxVectors

    return simulation 

def run_simulation(simulation,n_ts=10000,is_restart=False,is_enhanced=False):
    """Run simulation given current setup.

    Parameters:
	simulation: *openmm.app.simulation*
	    Preset simulation with langevin integrator.
	n_ts: *int*
	    Number of timesteps to run. Default=10,000
	is_restart: *bool*
	    Is this simulation a restart
	is_enhanced: *bool*
	    Is this an enhanced sampling simulation
    Returns:
	simulation: *openmm.app.simulation*
	    Preset simulation with langevin integrator.
    """

    if is_restart:
        #Attempt to load restart 1, if failed, try 0
        print("Restarting simulation...")
        try:
            simulation.loadCheckpoint("restart_1.xml")
        except:
            simulation.loadCheckpoint("restart_0.xml")
    else:
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
    
    #Set up pysages cvs
    cvs = [DihedralAngle([26,27,28,22]),DihedralAngle([27,28,22,11])]
    grid = pysages.Grid(lower=(-np.pi, -np.pi), upper=(np.pi, np.pi), shape=(64, 64), periodic=True)

    method = SpectralABF(cvs, grid) 
    
    
    #Run simulation
    try:
        print("Running!")
        if is_enhanced:
            # Create histogram logger
            histogram_logger = HistogramLogger(period=100)

            run_result = pysages.run(method,autogen_simulation,n_ts,
                    callback=histogram_logger)
            print("Simulation completed! Dumping...")
           
            with open("run_res.pickle",'wb') as f:
                pickle.dump(run_result,f) 
            with open("run_hist.pickle",'wb') as g:
                pickle.dump(histogram_logger,g)
            print("Pysages data dumped to 'run_res.pickle'.")
            print("Histograms dumped to 'run_hist.pickle'.")
            print("Analyzing...")
            return simulation
        else:
            #Old style running:
            for i in range(n_ts//traj_int):
                simulation.step(traj_int)
                pdb_fn = 'structure_%s.pdb' % i
                #app.PDBFile.writeFile(simulation.topology,simulation.getpositions,pdb_fn) 
                #print("simulation.topology:",simulation.topology())
                #print("simulation positions:",simulation.getPositions())
                
                fn = "restart_%s.xml" % (i % 2)
                with open(fn,'wb') as f:
                    simulation.saveCheckpoint(f)
    finally:
        for reporter in simulation.reporters:
            try:
                reporter.file.close()
            except:
                print("Couldn't close reporter %s" % str(reporter))
    
    print("Finished!")
    return simulation
def plot_data(pickle_fn='run_res.pickle'):
    """Plot the output of a crashed pysages run.

    Parameters:
	pickle_fn: *str*
	    Path to pickled pysages result.
    """
    with open(pickle_fn,'rb') as f:
        res = pickle.load(f)
    print("Training on data...")
    result = pysages.analyze(res)
    print("Network trained.")
    A = result['free_energy']
    A = A.max() - A#A.max() - A
    grid = pysages.Grid(lower=(-np.pi, -np.pi), upper=(np.pi, np.pi), shape=(64, 64), periodic=True)
    A = A.reshape(grid.shape)
    if type(A) == np.array:
        np.savetxt("free_eng.txt",A)
    else:
        print("Type of A:",type(A))
    fig, ax = plt.subplots(dpi=120)
    A_kT = A/(2.479) # Convert from kJ/mol to kT
    A = A_kT
    im = ax.imshow(
        A, interpolation="bicubic", origin="lower", 
            extent=[-np.pi, np.pi, -np.pi, np.pi])
    ax.contour(
        A, levels=list(range(0,int(A.max()),10)), linewidths=0.75, colors="k", 
            extent=[-np.pi, np.pi, -np.pi, np.pi])
    ax.set_xlabel(r"$\phi_2$ (rad)")
    ax.set_ylabel(r"$\phi_1$ (rad)")

    cbar = plt.colorbar(im)
    cbar.ax.set_ylabel(r"$A~[kT]$", rotation=270, labelpad=20)


    #plt.contour(A,linewidths=0.75,colors='black',extent=[-np.pi,np.pi,-np.pi,np.pi])
    #plt.set_xlabel("$\\phi_2$ (deg)")
    #plt.set_xlabel("$\\phi_3$ (deg)")
    plt.savefig("free_eng.pdf")
    plt.savefig("free_eng.png",dpi=600)
    save_energy_forces(result)
    plot_histogram(result)
    plot_forces(result)

def save_energy_forces(result):
    """Energy/force saving function from Ludwig/Pablo

    Parameters:
	result: *pysages.analyze.result*
	    Result of analyzed pysages simulation.
    """
    energy = numpy.asarray(result["free_energy"])
    forces = numpy.asarray(result["mean_force"])
    grid = numpy.asarray(result["mesh"])
    numpy.savetxt("FES.csv", numpy.hstack([grid, energy.reshape(-1, 1)]))
    numpy.savetxt("Forces.csv", numpy.hstack([grid, forces.reshape(-1, grid.shape[1])]))

def plot_histogram(result):
    surface = numpy.asarray(result["histogram"]) / numpy.nanmax(numpy.asarray(result["histogram"]))
    fig, ax = plt.subplots()
    im = ax.imshow(
        surface, interpolation="bicubic", origin="lower", extent=[-np.pi, np.pi, -np.pi, np.pi], aspect=1
    )
    ax.contour(surface, levels=15, linewidths=0.75, colors="k", extent=[-np.pi, np.pi, -np.pi, np.pi])
    plt.colorbar(im)
    fig.savefig("histogram.pdf")

def plot_forces(result):
    fig, ax = plt.subplots()

    ax.set_xlabel("CV")
    ax.set_ylabel("Forces $[\\epsilon]$")

    forces = numpy.asarray(result["mean_force"])
    x = numpy.asarray(result["mesh"])
    plt.quiver(x, forces, width=(0.0002 * (x[x.shape[0] - 1, 0] - x[0, 0])), headwidth=3)

    fig.savefig("forces.pdf")

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
    parser.add_argument('-r','--restart',action='store_true',
	help="Restart existing simulation.")
    parser.add_argument('-e','--enhanced',action='store_true',
	help="Enhanced sampling simulation.")
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
