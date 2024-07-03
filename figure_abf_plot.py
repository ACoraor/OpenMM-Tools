#Module figure_abf_plot
#Created by Aria Coraor
#Written  7/3/24

import numpy as np
import mdtraj as md
import argparse
import matplotlib
matplotlib.use("Agg")
from matplotlib import rc
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scipy.integrate as integrate


def main():
    """Plot a faux 1-D ABF run"""
    
    #Generate the samples and the corresponding forces
    data, forces = make_samples(args.num)

    #Add random noise to forces
    forces = add_noise(forces)
    
    #Create ABF plot
    plot_abf(data,forces)


def make_samples(num):
    '''Create samples from a trimodal distribution, and forces.

    Parameters:
        num: *int*
            Number of samples to produce
    Returns:
        data: *np.array*, shape: (num,)
            Raw samples from distribution.
        forces: *np.array of float*, shape: (num,)
            Negative local Gradient of pdf.
    '''
    
    #Set up hyperparameters
    means = np.array([2.0,5.0,7.5])
    sigs = np.array([1.0,3.5,0.5])
    mags = np.array([0.25,0.6,0.15])
    mag_cumsum = np.cumsum(mags)
    
    #Generate rands to choose distributions
    rand_ns = np.random.uniform(size=num)
    #Identify totals for each dist
    totals = np.zeros(3,dtype=int)
    totals[0] = np.array([len(np.argwhere(rand_ns < mag_cumsum[0]).ravel())
    totals[1] = np.array([len(np.argwhere(rand_ns < mag_cumsum[1]).ravel())
    totals[0] = np.array([len(np.argwhere(rand_ns < mag_cumsum[2]).ravel())
    totals[1] = totals[1] - totals[0]
    totals[2] = totals[2] - totals[1] - totals[0]

    #Generate totals-many random samples

    data = []
    forces = []
    for i in range(len(means)):
        new_x = np.random.normal(loc=means[i],scale=sigs[i],size=totals[i])
        new_forces = g_pdf_deriv(new_x,means[i],sigs[i])
        data.append(new_x)
        forces.append(new_forces)
    data = np.concatenate(data)
    forces = np.concatenate(forces)  

    print("Data:",data)
    print("Forces:",forces)
    
    return data, forces

def add_noise(forces,mag = 0.05,sig = 0.05):
    '''Add random noise to the forces.

    Parameters:
        forces: *np.array of float*, shape: (num,)
            Negative local Gradient of pdf.
    Returns:
        forces: *np.array of float*, shape: (num,)
            Negative local Gradient of pdf, plus noise .
    '''
    
    return forces + np.random.normal(loc=mag,scale=sig,size=forces.shape)

def g_pdf(x,mu,sig):
    """Calculate gaussian pdf value at x, with parameters mu and sigma."""
    a = 1/(np.sqrt(2*np.pi*sig**2))
    b = mu
    c = (2*sig**2)
    return 1/(np.sqrt(2*np.pi*sig**2)) * np.exp((-(x-mu)**2)/(2*sig**2))

def g_pdf_deriv(x,mu,sig):
    """Calculate the derivative of pdf value at x, with parameters mu and sigma."""
    a = 1/(np.sqrt(2*np.pi*sig**2))
    b = mu
    c = (2*sig**2)
    return -(2*a*(x-b))/c * np.exp(-((x-b)**2)/c)

def plot_abf(data, forces):
    '''Create a plot without the ABF biasing added, and a second plot
        with the ABF biasing added.

    Parameters:
        data: *np.array*, shape: (num,)
            Raw samples from distribution.
        forces: *np.array of float*, shape: (num,)
            Negative local Gradient of pdf.
    '''
    
    plt.clf()
    fig,ax = plt.subplots(layout='constrained')
    
    means = np.array([2.0,5.0,7.5])
    sigs = np.array([1.0,3.5,0.5])
    mags = np.array([0.25,0.6,0.15])

    x = np.linspace(0,10,1000)
    baseline_pdf = np.zeros(x.shape)
    for i in range(len(means)):
        baseline_pdf = baseline_pdf + mags[i]*g_pdf(x,means[i],sigs[i])

    baseline_pdf = baseline_pdf / np.sum(baseline_pdf)

    #Integrate forces
    n_bins = 32
    edges = np.linspace(0,10,n_bins+1)
    medians = (edges[:-1] + edges[1:])/2
    #x_sorted = np.sort(data)
    bin_pops = np.zeros(n_bins)
    bias_forces = np.zeros(n_bins)
    
    #Loop through all forces, adding to histogram bins
    for i, force in enumerate(forces):
        bin_ind = np.searchsorted(medians,force)
        bin_pops[bin_ind] = bin_pops[bin_ind] + 1
        bias_forces[bin_ind] = bias_forces[bin_ind] + force
    bias_forces = bias_forces / bin_pops

    bias = np.cumsum(bias_forces * (edges[1]-edges[0]))
    bias = bias - np.min(bias)

    biased_pdf = np.exp(-(-np.log(baseline_pdf) + bias))
    biased_pdf = biased_pdf / np.sum(biased_pdf)
    
    ax.hist(data,bins = edges,density=True,color='red',edgecolor='black',
        linewidth='1',label = 'histogram')

    ax2 = ax.twinx()
    ax2.plot(x,baseline_pdf,color='blue',label="Unbiased pdf")
    ax2.plot(x,biased_pdf,color='blue',label="Biased pdf")
    ax2.set_ylabel("Marginal probability")

    ax.set_xlabel("Collective variable")
    ax.set_ylabel("Probability")

    fig.tight_layout()
 
    fig.savefig('abf_diagram.svg')
    fig.savefig('abf_diagram.png',dpi=600)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-n','--num', type=int, default=250,
        help='Number of samples to include in abf generation. Default = 250.')
    args = parser.parse_args()

    main()

