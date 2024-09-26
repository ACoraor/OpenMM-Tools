#Module abf_plot
#Created by Aria Coraor
#Written 4/4/24

import numpy as np
import os
import argparse
import matplotlib
matplotlib.use("Agg")
from matplotlib import rc
import matplotlib.pyplot as plt
import matplotlib.cm as cm


def main():
    """Plot the raw ABF integration FES properly."""
    fes = np.loadtxt(args.file)
    
    #Check for 1D FES
    if fes.shape[1] == 3:
        make_plot(fes,args.max)
    elif fes.shape[1] == 2:
        make_1d_plot(fes,args.max)
        make_1d_hist(fes)
    else:
        raise Exception("Number of columns in %s can't be processed." % args.file)
    fn = 'abf_fes.png'
    print("Saved converted %s." % fn)

def make_plot(fes,max_eng=20):
    '''Create the plot.

    Parameters:
        fes: *np.array*, shape: (len(x)*len(y),3)
            First two columns are x and y, last column is free energy
                in kJ/mol
        max_eng: *int or float*
            Maximum energy to plot up to. Default = 20.
    '''
    
    plt.clf()
    #rc('font',**{'family':'sans-serif','sans-serif':['Arial']})
    #rc('text',usetex=True)
    #plt.rcParams['font.family']='sans-serif'
    #plt.rcParams['font.sans-serif']=['Arial']
    fig,ax = plt.subplots(layout='constrained')
    energies = fes[:,2]
    energies -= np.min(energies)
    energies /= 2.479 # Convert kJ/mol to kT
    cmap = cm.get_cmap('viridis',int(max_eng))

    

    norm = matplotlib.colors.Normalize(vmin=0,vmax=max_eng)
    sm = cm.ScalarMappable(norm=norm,cmap=cmap)

    #Reshape into meshgrid for plotting
    XX = np.reshape(fes[:,0]*180.0/np.pi,(200,200))
    #YY = np.reshape(fes[:,1]*10.0,(200,200))
    YY = np.reshape(fes[:,1]*180.0/np.pi,(200,200))
    ZZ = np.reshape(energies,(200,200))

    #For soft_041: flip xx and yy
    #tmp = np.copy(XX)
    #XX = YY
    #YY = tmp

    ax.contourf(XX,YY,ZZ,levels = np.arange(max_eng),
        norm=norm,cmap=cmap)
    ax.contour(XX,YY,ZZ,levels = np.arange(max_eng),
        colors = "black",linewidths=0.5)
    cbar = plt.colorbar(mappable=sm,label = "Free energy (kT)",ax=ax,
        orientation='vertical')
    ax.set_xlabel("Backbone Angle (deg)")
    #ax.set_xlabel("Intercalation Angle (deg)") #soft_041
    #ax.set_ylabel("Intercalation distance ($\\AA$)")
    #ax.set_xlabel("Intercalation distance ($\\AA$)")
    ax.set_ylabel(r"$\theta_{bending}$ (deg)") #soft_041
    #ax.set_xlabel(r"$\theta_{bending}$ (deg)")
    plt.savefig('abf_fes.png',dpi=600)
    plt.savefig('abf_fes.pdf')
    plt.savefig('abf_fes.eps')


def make_1d_plot(fes,max_eng=20):
    '''Create a 1d fes plot.

    Parameters:
        fes: *np.array*, shape: (len(x)*len(y),3)
            First two columns are x and y, last column is free energy
                in kJ/mol
        max_eng: *int or float*
            Maximum energy to plot up to. Default = 20.
    '''
    
    plt.clf()
    fig,ax = plt.subplots(layout='constrained')
    energies = fes[:,1]
    energies = energies / 2.479 # Convert kJ/mol to kT
    energies = energies - np.min(energies)
    XX = fes[:,0]*(180.0/np.pi)
    ZZ = energies

    plt.plot(XX,ZZ)

    ax.set_xlabel(r"Bending angle $\theta$ (deg)")
    ax.set_ylabel(r"Free energy (kT)") #soft_041
    ax.set_ylim(bottom=0.0,top=max_eng)
    plt.savefig('abf_fes.png',dpi=600)
    plt.savefig('abf_fes.pdf')
    plt.savefig('abf_fes.eps')

    #Now create probability plot
    plt.clf()
    fig,ax = plt.subplots(layout='constrained')
    #energies = fes[:,1]
    #energies -= np.min(energies)
    #energies /= 2.479 # Convert kJ/mol to kT
    XX = fes[:,0]*(180.0/np.pi)
    #Calculate probability distribution, reweighted to sum to 1
    prob = np.exp(-energies)
    prob = prob / np.sum(prob)
    #print("prob:",prob)

    ZZ = prob

    plt.plot(XX,ZZ)
    ax.set_xlabel(r"Bending angle $\theta$ (deg)")
    ax.set_ylabel(r"Marginal probability") #soft_041
    max_prob = np.max(prob)
    ax.set_ylim(bottom=0.0,top=max_prob*1.1)
    plt.savefig('abf_prob.png',dpi=600)
    plt.savefig('abf_prob.pdf')
    plt.savefig('abf_prob.eps')
    
def make_1d_hist(fes):
    '''Create a 1d histogram plot to check convergence

    Parameters:
        fes: *np.array*, shape: (len(x)*len(y),3)
            First two columns are x and y, last column is free energy
                in kJ/mol
    '''
    #Load histograms: Most recent, less recent.
    
    curr_hist = np.loadtxt('histogram.txt')

    dirs = [fn for fn in os.path.listdir('.') if "output_" in fn]
    dirs.sort(key = lambda x: int(x[7:]))
    #attempt to load hists
    hists = []
    for d in dirs:
        fn = os.path.join(d,'histogram.txt')
        hists.append(np.loadtxt(fn))
    
    burn_in_hist = hists[len(hists)//3]
    second_last_hist = hists[2*len(hists)//3]
    
    #Check if histogram from 1/3rd --> 2/3rds is the same as 2/3rds --> on.
    curr_hist = curr_hist - second_last_hist
    second_last_hist = second_last_hist - burn_in_hist
    


    plt.clf()
    fig,ax = plt.subplots(layout='constrained')
    XX = fes[:,0]*(180.0/np.pi)
    x_grid = np.linspace(XX[0],XX[-1],len(curr_hist))

    plt.plot(x_grid,second_last_hist,color='#E69F00',label = "Middle third")
    plt.plot(x_grid,curr_hist,color='#56B4E9',label = "Last third")

    ax.set_xlabel(r"Bending angle $\theta$ (deg)")
    ax.set_ylabel(r"Histogram probability") #soft_041
    #ax.set_ylim(bottom=0.0,top=max_eng)
    plt.savefig('conv_hist.png',dpi=600)
    plt.savefig('conv_hist.pdf')
    plt.savefig('conv_hist.eps')



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f','--file', type=str, default='abf_fes_integrated.txt',
        help='Path to ABF integrated output')
    parser.add_argument('-m','--max', type=float, default=40.0,
        help='Maximum free energy to plot up to. Default = 40.0.')
    args = parser.parse_args()

    main()

