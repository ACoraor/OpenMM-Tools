#Module run_openmm
#Created by Aria Coraor 6/2/23


import numpy as np
import argparse
#Following openMM docs
from openmm import *
import os
import openmm
#import openmmtools.forces as ommforces
import sys
import openmm.app as app
import openmm.unit as unit
import openmmtools.forces as ommforces
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mdtraj as md
from shutil import copy2
from openmm.vec3 import Vec3
import pickle
from jax.lib import xla_bridge

#Import pysages things
import pysages
from pysages.grids import Grid
from pysages.colvars import DihedralAngle,Distance,Angle,Acylindricity,Asphericity, Component
from pysages.methods import SpectralABF,ABF,HarmonicBias
from pysages.methods.utils import HistogramLogger

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
    simulation = generate_simulation(pdb,ff,args.restart)
    '''
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
    '''
    print("Simulation generated.")

    #Add reporters
    if not args.enhanced:
        simulation = add_reporters(simulation,args.restart)

    #Run the system
    print("Running simulation.")
    results = run_simulation(simulation,args.timesteps,args.restart,args.enhanced)

    print("Simulation complete.")

    if args.enhanced:
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
    #pdb,ff = load_system('init.pdb','DNA.bsc1.xml','spce_standard.xml',
    #    'gaff-aec.xml','gaff.xml')
    pdb,ff = load_system('init.pdb','DNA.bsc1.xml','opc_standard.xml',
        'gaff-aec_v3.xml','gaff.xml')

    '''#Old style:
    try:
        T = 298.15 * unit.kelvin
        dt = 2.0 * unit.femtoseconds
        cutoff_distance = 1.0 * unit.nanometer
        system = ff.createSystem(
            pdb.topology, constraints = app.HBonds, 
            nonbondedMethod = app.PME, nonbondedCutoff = cutoff_distance)
        integrator = openmm.LangevinIntegrator(T, 1 / unit.picosecond, dt)
        simulation = app.Simulation(pdb.topology, system, integrator)
    except:
        print("Failed to restart.")
        simulation = generate_simulation(pdb,ff)
    '''
    simulation = generate_simulation(pdb,ff)
    try:
        simulation.loadCheckpoint("restart_m.chk")
    except:
        print("Failed standard restart.")
        try:
            simulation.loadCheckpoint("restart_1.xml")
        except:
            simulation.loadCheckpoint("restart_0.xml")
    #simulation.minimizeEnergy()
    add_reporters(simulation,is_restart=True)
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

    n_atoms = 1600
    relevant_inds = np.arange(n_atoms) # Manually counted)
    #print("simulation.context.__dict__:",simulation.context.__dict__.keys())
    try:
        top = simulation.topology
        chains = [ch for ch in top.chains()]
        dna_inds = []
        for ch in chains[:2]:
            for at in ch.atoms():
                dna_inds.append(at.index)
        dna_inds = np.array(dna_inds,dtype=int)
        relevant_inds = dna_inds
    except Exception as e:
        print("Failed dna inds setting:",str(e))
        

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
    simulation.reporters.append(app.StateDataReporter(sys.stdout,5000,
        step=True,potentialEnergy=True,temperature=True))
    chk_file = open('restart_m.chk','wb')
    simulation.reporters.append(app.CheckpointReporter(chk_file,10000))

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
    molar = 1.0* unit.molar 
    pos = pdb.positions

    #pos *= 0.1 # Convert from angstroms to nanometers
    if not is_restart: 
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
            

        #Set position at center
        com = pos[0]
        for i in range(len(pos)-1):
            com += pos[i+1]
        com = com / len(pos)
        #Loop through pos to maintain correct typing
        for i in range(len(pos)):
            pos[i] = pos[i] - com
        
        #res = [r for r in topology.residues()]
        #res[0].name = "DG5" 
        #Add spce solvent molecules
        modeller = app.modeller.Modeller(topology,pos)
        
        bounds = Vec3(6.0,6.0,11.0)*unit.nanometer #Change
        #dna_ff = app.ForceField(args.amber) 
        modeller.loadHydrogenDefinitions(amber)
        #print("Adding hydrogens...")
        modeller.addHydrogens(ff)
        print("Adding water...")
        modeller.addSolvent(ff,boxSize=bounds,model='tip4pew',ionicStrength=0.15*molar)
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

    #Add custom restraint forces
    k_phos = 3.075 * unit.kilojoules_per_mole/(unit.nanometer**2)

    #bound_ends = [(1489,30),(788,725)] #Pair of tuples to be bound
    #print("Should be:",bound_ends)
    for end in bound_ends:
        harm_rest_force = ommforces.HarmonicRestraintBondForce(
            spring_constant = k_phos, restrained_atom_index1 = end[0],
            restrained_atom_index2 = end[1])
        system.addForce(harm_rest_force)
    #Counter-rotation forces
    k_rot = k_phos*400#400
    restraint1 = CustomExternalForce('k_rot*periodicdistance(x,y,0.0,0.0,0.0,0.0)^2')
    restraint1.addGlobalParameter('k_rot',k_rot)
    restraint1.addParticle(bound_ends[0][1])
    restraint1.addParticle(bound_ends[1][1])
    system.addForce(restraint1)

    #Set up pysages cvs
    
    #Create Simulation

    #Angle,distance
    #For soft_073 - soft_076, determine end phosphates:
    #bound_ends = [(1109,31),(537,600)] #Pair of tuples to be bound

    ksprings = [3.075,3.075] #1/10th current distribution constant
    #in kJ/mol(nm^2)
    centers = [1.96,1.96]

    #Angle,distance
    ang = Angle((845,1164,1484))# Manually identified 11/4/24

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

    
    

    #Create harmonic restraints to prevent rotation and dehybridization
    sys = simulation.context.getSystem()
    forces = sys.getForces()
    print("Forces:")
    for f in forces:
        print(f) 

    #dehybridization forces
    #k_phos = 3.075 * unit.kilojoules_per_mole/(unit.nanometer**2)
    #Determine bound end indices at runtime from phosphates
    bound_ends = [[0,0],[0,0]]
    chains = [ch for ch in simulation.topology.chains()]
    ch0 = chains[0]
    res = [r for r in ch0.residues()]
    end_1A_phos = [at.index for at in res[1].atoms() if "P" == at.name][0]
    end_2A_phos = [at.index for at in res[-1].atoms() if "P" == at.name][0]
    ch1 = chains[1]
    res = [r for r in ch1.residues()]
    end_2B_phos = [at.index for at in res[1].atoms() if "P" == at.name][0]
    end_1B_phos = [at.index for at in res[-1].atoms() if "P" == at.name][0]
    bound_ends = [(end_1B_phos,end_1A_phos),(end_2B_phos,end_2A_phos)]
    print("Bound_ends:",bound_ends)
    #Only add these forces if not a restart, otherwise they are already there
    '''
    if not is_restart:
        k_phos = 3.075 * unit.kilojoules_per_mole/(unit.nanometer**2)

        #bound_ends = [(1489,30),(788,725)] #Pair of tuples to be bound
        #print("Should be:",bound_ends)
        for end in bound_ends:
            harm_rest_force = ommforces.HarmonicRestraintBondForce(
                spring_constant = k_phos, restrained_atom_index1 = end[0],
                restrained_atom_index2 = end[1])
            sys.addForce(harm_rest_force)
        #Counter-rotation forces
        k_rot = k_phos*400#400
        restraint1 = CustomExternalForce('k_rot*periodicdistance(x,y,0.0,0.0,0.0,0.0)^2')
        #restraint1 = CustomExternalForce('k_rot*((x-x0)^2+(y-y0)^2)')
        #restraint1.addPerParticleParameter("x0")
        #restraint1.addPerParticleParameter("y0")
        restraint1.addGlobalParameter('k_rot',k_rot)
        #restraint1.addParticle(bound_ends[0][1],(0.0,0.0)*unit.nanometer)
        #restraint1.addParticle(bound_ends[1][1],(0.0,0.0)*unit.nanometer)
        restraint1.addParticle(bound_ends[0][1])
        restraint1.addParticle(bound_ends[1][1])
        sys.addForce(restraint1)

        #Set up pysages cvs
        
        #Rebuild the simulation using everything here
        old_integ = simulation.integrator
        integrator = openmm.LangevinIntegrator(old_integ.getTemperature(), 
            old_integ.getFriction(), old_integ.getStepSize())
        biased_sim = app.Simulation(simulation.topology,sys,integrator)
        pos = simulation.context.getState(getPositions=True).getPositions()
        biased_sim.context.setPositions(pos) 

        for reporter in simulation.reporters:
            biased_sim.reporters.append(reporter)

        #Replace original simulation with biased simulation
        tmp_sim = simulation
        simulation = biased_sim
    '''
    #Angle,distance
    #For soft_073 - soft_076, determine end phosphates:
    #bound_ends = [(1109,31),(537,600)] #Pair of tuples to be bound

    ksprings = [3.075,3.075] #1/10th current distribution constant
    #in kJ/mol(nm^2)
    centers = [1.96,1.96]

    #Angle,distance
    ang = Angle((845,1164,1484))# Manually identified 11/4/24
    
    #dna_inds = np.arange(3216,dtype=int)
    #acyl = Acylindricity(dna_inds,"xy")
    #asphr = Asphericity(dna_inds)
    cvs = [ang]
    #cvs = [acyl,ang]
    #grid = pysages.Grid(lower=(0.0, np.pi/2), upper=(60.0,np.pi), shape=(32, 32), periodic=False)
    #cvs = [asphr]
    grid = pysages.Grid(lower=(np.pi/2,),upper=(np.pi,),shape=(64,),periodic=False) 

    #method = SpectralABF(cvs, grid) 
    method = ABF(cvs, grid) 
    
    
    #Run simulation
    try:
        print("Running!")
        if is_enhanced:
            # Create histogram logger
            histogram_logger = HistogramLogger(period=100)
            try:
                print("Attempting to restart from run_res.pickle")
                result = pysages.load('run_res.pickle')
                run_result = pysages.run(result,autogen_simulation,n_ts)#,
                    #callback=histogram_logger)
            except Exception as e:
                print("Failed to restart from run_res.pickle.")        
                print("error:",e)
                run_result = pysages.run(method,autogen_simulation,n_ts,
                    callback=histogram_logger)
            print("Simulation completed! Dumping...")
           
            pysages.save(run_result,"run_res.pickle")
            #with open("run_res.pickle",'wb') as f:
            #    pickle.dump(run_result,f) 
            #with open("run_hist.pickle",'wb') as g:
            #    pickle.dump(histogram_logger,g)
            print("Pysages data dumped to 'run_res.pickle'.")
            #print("Histograms dumped to 'run_hist.pickle'.")
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
    result = pysages.analyze(res,topology=(4,))
    #result = pysages.analyze(res,topology=(14,))
    print("Network trained.")
    A = result['free_energy']
    #A = A.max() - A#A.max() - A
    A = A - A.min()
    A_kT = A/(2.479) # Convert from kJ/mol to kT
    A = A_kT
    ####NOTE!!!!
    #Convert a to degrees, r to angstroms
    grid_bounds = [np.pi/2,np.pi] # low_1, high_1, low_2, high_2
    grid = pysages.Grid(lower=(grid_bounds[0],), upper=(grid_bounds[1],), shape=(64), periodic=False)
    #grid = pysages.Grid(lower=(-np.pi, -np.pi), upper=(np.pi, np.pi), shape=(64, 64), periodic=True)
    A = A.reshape(grid.shape)
    if type(A) == np.array:
        np.savetxt("free_eng.txt",A)
    else:
        print("Type of A:",type(A))
    
    #Begin plotting
    if len(grid.shape) == 2:
        fig, ax = plt.subplots(dpi=600)
        im = ax.imshow(
            A, interpolation="bicubic", origin="lower", 
                extent=[grid_bounds[0],grid_bounds[1]])
        ax.contour(
            A, levels=list(range(0,int(A.max()),10)), linewidths=0.75, colors="k", 
                extent=[grid_bounds[0],grid_bounds[1]])
        ax.set_xlabel(r"$a$ (deg)")
        ax.set_ylabel(r"$\theta_{bending}$ (deg)")

        #####
        cbar = plt.colorbar(im)
        cbar.ax.set_ylabel(r"$A~[kT]$", rotation=270, labelpad=20)


        #plt.contour(A,linewidths=0.75,colors='black',extent=[-np.pi,np.pi,-np.pi,np.pi])
        #plt.set_xlabel("$\\phi_2$ (deg)")
        #plt.set_xlabel("$\\phi_3$ (deg)")
        plt.gca().set_aspect('auto')
        plt.savefig("free_eng.pdf")
        plt.savefig("free_eng.png",dpi=600)
    else:
        print("Not a 2D FES. Saving...")
    save_energy_forces(result)
    try:
        plot_histogram(result)
    except Exception as e:
            print("Histogram plotting failed under conditions:")
            print(e)
    try:
        plot_forces(result)
    except Exception as e:
            print("Forces plotting failed under conditions:")
            print(e)

def save_energy_forces(result):
    """Energy/force saving function from Ludwig/Pablo

    Parameters:
	result: *pysages.analyze.result*
	    Result of analyzed pysages simulation.
    """
    energy = numpy.asarray(result["free_energy"])
    forces = numpy.asarray(result["mean_force"])
    grid = numpy.asarray(result["mesh"])
    print("grid.shape:",grid.shape)
    print("energy.shape:",energy.shape)
    #reshape grid to be column-like
    if len(grid.shape) == 1:
        grid = np.reshape(grid,(-1,1))
    numpy.savetxt("FES.csv", numpy.hstack([grid, energy.reshape(-1, 1)]))
    numpy.savetxt("Forces.csv", numpy.hstack([grid, forces.reshape(-1, grid.shape[1])]))

def plot_histogram(result):
    np.savetxt('histogram.txt',np.array(result["histogram"]))
    np.savetxt('histogram_mesh.txt',np.array(result["mesh"]))
    surface = numpy.asarray(result["histogram"]) / numpy.nanmax(numpy.asarray(result["histogram"]))
    fig, ax = plt.subplots()
    grid_bounds = [np.pi/2,np.pi] # low_1, high_1, low_2, high_2
    im = ax.imshow(
        surface, interpolation="bicubic", origin="lower", 
        extent=[grid_bounds[0],grid_bounds[1]], aspect=1
    )
    ax.contour(surface, levels=15, linewidths=0.75, colors="k", 
        extent=[grid_bounds[0],grid_bounds[1]])
    plt.colorbar(im)
    ax.set_xlabel(r"$a$ (deg)")
    ax.set_ylabel(r"$\theta_{bending}$ (deg)")
    plt.gca().set_aspect('auto')
    fig.savefig("histogram.pdf")
    fig.savefig("histogram.png",dpi=600)

def plot_forces(result):
    fig, ax = plt.subplots()
    ax.set_xlabel(r"$a$ (deg)")
    ax.set_ylabel(r"$\theta_{bending}$ (deg)")
    #ax.set_xlabel(r"$a$ (rad)")
    #ax.set_ylabel(r"$r$ (nm)")

    #ax.set_xlabel("CV")
    #ax.set_ylabel("Forces $[\\epsilon]$")

    forces = numpy.asarray(result["mean_force"])
    x = numpy.asarray(result["mesh"])
    forces = np.reshape(forces,x.shape)
    plt.quiver(x[:,0],x[:,1]*(180.0/np.pi),forces[:,0],forces[:,1])#, width=(0.0002 * (x[x.shape[0] - 1, 0] - x[0, 0])), headwidth=3)

    plt.gca().set_aspect('auto')
    fig.savefig("forces.pdf")
    fig.savefig("forces.png",dpi=600)

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
    parser.add_argument('-w','--water',type=str,default = "opc_standard.xml",
	help="Path to Water/SPCE parameter file.")
    parser.add_argument('-d','--dyes',type=str,default = "gaff-aec_v3.xml",
	help="Path to dyes parameter file.")
    parser.add_argument('-g','--gaff',type=str,default = "gaff.xml",
	help="Path to GAFF parameter file.")
    parser.add_argument('-t','--timesteps',type=int,default = 10000000,
	help="Number of timesteps for simulation.")
	
    args = parser.parse_args()
    print(xla_bridge.get_backend().platform)
    main()
