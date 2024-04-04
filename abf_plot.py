#Module abf_plot
#Created by Aria Coraor
#Written 4/4/24

import numpy as np
import mdtraj as md
import argparse
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm


def main():
    """Plot the raw ABF integration FES properly."""
    fes = np.loadtxt(args.file)

    make_plot(fes)
    print("Saved converted %s." % fn)

def make_plot(fes):
    '''Create the plot.

    Parameters:
        fes: *np.array*, shape: (len(x)*len(y),3)
            First two columns are x and y, last column is free energy
                in kJ/mol
    '''

    plt.clf()
    energies = fes[:,2]
    energies -= np.min(energies)
    energies /= 2.479 # Convert kJ/mol to kT
    cmap = cm.Viridis_r
    norm = matplotlib.colors.Normalize(vmin=0,vmax=20)
    sm = cm.ScalarMappable(norm=norm,cmap=cmap)

    plt.contourf(fes[:,0]*180/np.pi,fes[:,1]*10.0,energies,levels = np.arange(20),
        norm=norm,cmap=cmap)
    cbar = plt.colorbar(mappable=sm,label = "Free energy (kT)")
    plt.xlabel("Angle (deg)")
    plt.ylabel("Intercalation distance ($\\AA$)")
    plt.tight_layout()
    plt.savefig('abf_fes.png',dpi=600)
    plt.savefig('abf_fes.pdf')



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f','--file', type=str, default='output.h5',
        help='Path to ABF integrated output')
    args = parser.parse_args()

    main()

