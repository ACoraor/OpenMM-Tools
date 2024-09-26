#Module figure_abf_plot
#Created by Aria Coraor
#Written  7/3/24

import os
import numpy as np
#import mdtraj as md
import argparse
import matplotlib
matplotlib.use("Agg")
from matplotlib import rc
import matplotlib.pyplot as plt
import matplotlib.cm as cm
#import scipy.integrate as integrate


def main():
    """Plot a faux 1-D ABF run"""
    
    #Generate the samples and the corresponding forces
    data, forces = make_samples(args.num)
    print("Data generated.")
    
    #Add random noise to forces
    forces = add_noise(forces,args.mag)
    print("Noise added.")    

    #Create ABF plot
    plot_abf(data,forces,args.grid,args.mag,args.switch)

    print("Figures generated. Exiting...")

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
    sigs = np.array([1.0,1.5,0.5])
    mags = np.array([0.25,0.6,0.15])
    #means = np.array([5.0])
    #sigs = np.array([2.0])
    #mags = np.array([1.0])
    mag_cumsum = np.cumsum(mags)
    
    #Generate rands to choose distributions
    np.random.seed(12345)
    rand_ns = np.random.uniform(size=num)
    #Identify totals for each dist
    totals = np.zeros(len(mags),dtype=int)
    
    for i in range(len(totals)):
        totals[i] = len(np.argwhere(rand_ns < mag_cumsum[i]).ravel())
        totals[i] = totals[i] - np.sum(totals[:i])
    #totals[1] = len(np.argwhere(rand_ns < mag_cumsum[1]).ravel())
    #totals[2] = len(np.argwhere(rand_ns < mag_cumsum[2]).ravel())
    #print("totals:",totals)
    #totals[1] = totals[1] - totals[0]
    #totals[2] = totals[2] - totals[1] - totals[0]

    #Generate totals-many random samples
    print("totals:",totals)
    data = []
    forces = []
    for i in range(len(means)):
        new_x = np.random.normal(loc=means[i],scale=sigs[i],size=totals[i])
        new_forces = g_pdf_deriv(new_x,means[i],sigs[i])
        data.append(np.copy(new_x))
        forces.append(np.copy(new_forces))
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
        mag: *float*
            Magnitude of gaussian noise to add to forces
        sig: *float*
            Stdev of gaussian noise to add to forces
    Returns:
        forces: *np.array of float*, shape: (num,)
            Negative local Gradient of pdf, plus noise .
    '''
    
    return forces + np.random.normal(loc=0.0,scale=mag,size=forces.shape)

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
    #return -(2*a*(x-b))/c * np.exp(-((x-b)**2)/c)
    prob = g_pdf(x,mu,sig)
    #Wolframalpha: d(E^U(x),x) = E^U(x) * dU/dx = F * E^U(x)
    return (-(2*a*(x-b))/c * np.exp(-((x-b)**2)/c))/prob

def plot_abf(data, forces,grid=32,mag=.05,switch=2500):
    '''Create a plot without the ABF biasing added, and a second plot
        with the ABF biasing added.

    Parameters:
        data: *np.array*, shape: (num,)
            Raw samples from distribution.
        forces: *np.array of float*, shape: (num,)
            Negative local Gradient of pdf.
        mag: *float*
            Magnitude of gaussian noise to add to forces
    '''
    
    plt.clf()
    font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : 16}

    matplotlib.rc('font', **font)
    fig,ax = plt.subplots(layout='constrained')
    
    #means = np.array([2.0,5.0,7.5])
    #sigs = np.array([1.0,2.5,0.5])
    #mags = np.array([0.25,0.6,0.15])
    means = np.array([2.0,5.0,7.5])
    sigs = np.array([1.0,1.5,0.5])
    mags = np.array([0.25,0.6,0.15])
    #means = np.array([5.0])
    #sigs = np.array([2.0])
    #mags = np.array([1.0])

    x = np.linspace(0,10,1000)
    baseline_pdf = np.zeros(x.shape)
    for i in range(len(means)):
        baseline_pdf = baseline_pdf + mags[i]*g_pdf(x,means[i],sigs[i])

    baseline_pdf = baseline_pdf / np.sum(baseline_pdf)

    #Integrate forces
    n_bins = grid
    edges = np.linspace(0,10,n_bins+1)
    medians = (edges[:-1] + edges[1:])/2
    #x_sorted = np.sort(data)
    bin_pops = np.zeros(n_bins)
    bias_forces = np.zeros(n_bins)
    
    #Loop through all forces, adding to histogram bins
    for i, force in enumerate(forces):
        bin_ind = np.searchsorted(edges,data[i]) - 1
        if bin_ind < 0 or bin_ind >= len(bin_pops):
            continue
        bin_pops[bin_ind] = bin_pops[bin_ind] + 1
        bias_forces[bin_ind] = bias_forces[bin_ind] + force
    print("Bin pops:",bin_pops)
    bias_forces = bias_forces / np.where(bin_pops > 0,bin_pops,1.0)

    print("medians shape:",medians.shape)
    print("Bias forces shape:",bias_forces.shape)
    saver_biases = np.zeros((len(bias_forces),2))
    saver_biases[:,0] = medians
    saver_biases[:,1] = bias_forces
    np.savetxt("figure_abf_forces.txt",saver_biases)


    #Fudge factor
    #bias_forces = 2*bias_forces

    print("biases integrated:",bias_forces)

    bias = np.cumsum(-bias_forces * (edges[1]-edges[0]))
    bias = bias - np.max(bias)

    #Interpolate to create same size
    interp_bias = np.interp(x,medians,bias)

    
    baseline_eng = -np.log(baseline_pdf)
    #print("interp bias:",interp_bias)
    baseline_eng = baseline_eng - baseline_eng[len(baseline_eng)//2]
    interp_bias = interp_bias - interp_bias[len(interp_bias)//2]
    if os.path.isfile('figure_abf_fes_integrated.txt'):
        print("Loading pre-saved data.")
        int_data = np.loadtxt('figure_abf_fes_integrated.txt')
        interp_bias = int_data[:,1] - int_data[:,1][int_data.shape[1]//2]
    else:
        print("No saved real abf integration found.")
    total_eng = baseline_eng - interp_bias
    biased_pdf = np.exp(-total_eng)
    biased_pdf = biased_pdf / np.sum(biased_pdf)*len(x)/args.grid 
    
    #print("Biased pdf:",biased_pdf)
    #print("Base pdf:",baseline_pdf)
    
    baseline_pdf = baseline_pdf * len(x) * args.num / args.grid

    #Skew the samples to imply a uniform distribution occurring
    #Above 2500 samples, for half of datapoints
    unif_samples = max(len(data) - switch,0)
    subsampled_data = np.random.choice(data,size=min(switch,len(data)))
    
    unif_data = np.concatenate((subsampled_data,np.linspace(0.0,10.0,unif_samples)))


    #ln1 = ax.hist(data,bins = edges,color='red',edgecolor='black',
    #    linewidth='1',label = 'Samples')
    
    ln1 = ax.hist(unif_data,bins = edges,color='red',edgecolor='black',
        linewidth='1',label = 'Samples')

    #Plot just the bias distribution
    abf_bias_pdf = np.exp(-interp_bias)
    abf_bias_pdf = abf_bias_pdf / np.sum(abf_bias_pdf)*len(x)/args.grid

    ax2 = ax.twinx()
    ln2 = ax.plot(x,baseline_pdf,color='blue',label="True probability density")
    ln5 = ax2.plot(int_data[:,0],interp_bias,color='#009E73',label="ABF bias")
    #ln3 = ax2.plot(x,biased_pdf,color='black',label="Biased pdf")
    #ln4 = ax2.plot(x,abf_bias_pdf,color='orange',label="Direct ABF pdf")
    #ln2 = ax2.plot(x,baseline_eng,color='#56B4E9',label="True free energy")
    #if os.path.isfile('figure_abf_fes_integrated.txt'):
    #    ln5 = ax2.plot(int_data[:,0],interp_bias,color='#009E73',label="ABF bias")
    #else:
    #    ln4 = ax2.plot(x,-interp_bias,color='#E69F00',label="Direct ABF eng")
    #print("total_eng:",total_eng)
    #ln6 = ax2.plot(x,total_eng,color='#CC79A7',label="Total energy")
    #sampling_dist = np.exp(-total_eng)/np.sum(np.exp(-total_eng))
    sampling_dist = np.exp(-total_eng)
    #print("Raw sampling probs:",sampling_dist)
    sampling_dist = sampling_dist / np.sum(sampling_dist) * args.num *len(interp_bias) / args.grid
    #print("Reweighted sampling probs:",sampling_dist)
    ln3 = ax.plot(x,sampling_dist,color='black',label="Sampling distribution")
    
    #Add green quiver arrows to the top
    
    X = medians
    #Y = np.ones(medians.shape) * -2.5
    #Set Y exactly 2 units below ABF bias curve
    #Y = interp_bias - 2.0
    Y = np.interp(medians,x,interp_bias) - 1.0


    U = -bias_forces
    V = np.random.normal(0.0,mag,size=medians.shape)/np.sqrt(bin_pops)
    print("X.shape:",X.shape)
    print("Y.shape:",Y.shape)
    print("U.shape:",U.shape)
    print("V.shape:",V.shape)


    ln6 = ax2.quiver(X,Y,U,V,color='#009E73'),
        #label="ABF Forces")
    
    print("RMS biasing force:",np.sqrt(np.average(bias_forces**2)))

    ax2.set_ylabel("Free energy ($k_{B}T$)")

    ax.set_xlabel("Collective variable")
    ax.set_ylabel("Relative Probability")

    ax.set_ylim(0,np.max(baseline_pdf)*4)
    ax2.set_ylim(8*np.min(interp_bias),1.1*np.max(interp_bias))

    #ax.legend([ln2,ln5,ln3],['True Probability','ABF bias',
    #            'Sampling distribution'])
    ax.legend(loc='center left',edgecolor='grey')#CC79A7.legend(loc='best',edgecolor='grey')
    ax2.legend(loc='upper center',edgecolor='grey')#CC79A7.legend(loc='best',edgecolor='grey')

    fig.tight_layout()
 
    fig.savefig('abf_diagram.svg')
    #fig.savefig('abf_diagram.png',dpi=600)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-n','--num', type=int, default=250,
        help='Number of samples to include in abf generation. Default = 250.')
    parser.add_argument('-g','--grid', type=int, default=32,
        help='Size of integration grid. Default = 32.')
    parser.add_argument('-m','--mag', type=float, default=1.0,
        help='Magnitude of Gaussian force noise Default=1.0.')
    parser.add_argument('-s','--switch', type=int, default=2500,
        help='Fake number of timesteps at which ABF converges. Default=2500.')
    args = parser.parse_args()

    main()

