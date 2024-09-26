#Module plot_converge
#Created by Aria Coraor
#Written 4/30/24

import numpy as np
import argparse
import matplotlib
matplotlib.use("Agg")
from matplotlib import rc
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pandas as pd
import os
import zfit

def main():
    """Plot each cv versus time in a bastardization for poster."""
    #Load pandas dataframe
    df = pd.read_csv('log.dframe')
    cvs = df.columns[1:]

    #Cvs of choice: pe, rg2

    #Plot all pairs of cvs
    #for i, cv1 in enumerate(cvs[:-1]):
    #    for j, cv2 in enumerate(cvs[i+1:]):
    cv1 = 'pe'
    cv2 = "gyration_rg2"

    #cv1 = 'cy5_5-stack'
    #cv1 = "gyration_c_reduced"
    #cv2 = "angle_d15"
    

    double_plot(df['timestep'],df[cv1],df[cv2],cv1,cv2)
        
    print("Created all correlation plots. Exiting...")


def double_plot(timestep,data1,data2,name1,name2):
    '''Create a line plot for both datasets, and a comparison correlation plot
    normalized to PC space.

    Parameters:
        timestep: *arraylike*, shape: (n_ts)
                Timestep for each datapoint.
        data1: *arraylike*, shape: (n_ts)
            First column of data to be plotted.
        data2: *arraylike*, shape: (n_ts)
            Second column of data to be plotted.
        name1: *str*
            Name of the first column
        name2: *str*
            Name of the first column
    '''
    #Co-plot both lines
    plt.clf()
    font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : 16}

    matplotlib.rc('font', **font)
    ts_sec = 2e-15 # Timestep in seconds
    plot_times = timestep * ts_sec / (1e-9) #Plot in nanoseconds
    #data1 = data1[10000:]#[:10001]
    #data2 = data2[10000:]#[:10001]
    data1 = data1[:10001]
    data2 = data2[:10001]
    #data2 = data2*180.0/np.pi
    data2 = np.sqrt(data2)
    plot_times= plot_times[:10001]    
    
    fig,ax = plt.subplots(layout='constrained')
    fig.set_size_inches(10.3,3.2)
    #ax.set_aspect(3.2/10.3)
    secax = ax.twinx()
    #secax.set_aspect(3.2/10.3)
    ax.set_xlabel("Simulation time (ns)")
    ax.set_ylabel("Potential Energy (kJ/mol)") #soft_041
    #ax.set_ylabel("Cy5 5'-ward stacking distance (nm)") #soft_041
    #ax.set_ylabel("Gyration tensor acylindricity") #soft_041

    #n_running = 10
    n_running = 100
    running_filter = np.ones(n_running)/float(n_running)
    running_data1 = np.convolve(data1,running_filter,mode='valid')
    data1 = running_data1
    running_times = np.convolve(plot_times,running_filter,mode='valid')
    plot_times = running_times
    ax.plot(plot_times,data1,color='#E69F00',label="PE")
    ax.spines['left'].set_color('#E69F00')
    ax.yaxis.label.set_color('#E69F00')

    running_data2 = np.convolve(data2,running_filter,mode='valid')
    data2 = running_data2
    secax.plot(plot_times,data2,color='#56B4E9',label="$R_G$")
    secax.spines['right'].set_color('#56B4E9')
    secax.spines['left'].set_color('#E69F00')
    secax.tick_params(axis='y',colors='#56B4E9')
    secax.yaxis.label.set_color('#56B4E9')
    ax.tick_params(axis='y',colors='#E69F00')
    
    secax.set_ylabel("Radius of Gyration (nm)")
    #secax.set_ylabel("Fiber bending angle (deg)")
    
    #Set ylims to separate values
    IQD_pe = np.max(data1) - np.min(data1)
    IQD_rg = np.max(data2) - np.min(data2)

    pe_range = (np.min(data1) - 0.1*IQD_pe, np.max(data1) + 1.5*IQD_pe)
    rg_range = (np.min(data2) - 1.5*IQD_rg, np.max(data2) + 0.1*IQD_rg)
    #print("pe_range:",pe_range)
    #pe_range = (0.0,0.1)
    ax.set_ylim(pe_range[0],pe_range[1])
    #print("Ax ylim:",ax.get_ylim())
    secax.set_ylim(rg_range[0],rg_range[1])
    #plt.gca().set_aspect(3.2/10.3)

    if not os.path.isdir('outputs'):
        os.mkdir('outputs')
    #name = '%s_vs_%s' % (name1,name2)
    name = "convergence_plot_run"
    #name = "fangle_gyration_b"
    #plt.legend(loc='best',edgecolor='grey')
    fn = os.path.join('outputs',name+'.png')
    plt.savefig(fn,dpi=600)
    
    #fn = os.path.join('outputs',name+'.pdf')
    #plt.savefig(fn)
    fn = os.path.join('outputs',name+'.svg')
    plt.savefig(fn)
    quit()
    
    #Calculate cross-correlation
    white_data1 = (data1 - np.average(data1))/np.std(data1)
    white_data2 = (data2 - np.average(data2))/np.std(data2)
    cross_corr = white_data1 * white_data2
    plt.clf()
    fig,ax = plt.subplots(layout='constrained')

    ax.plot(plot_times,cross_corr,color='blue')
    ax.set_xlabel("Simulation time (ns)")
    ax.set_ylabel('Cross Correlation') #soft_041
    plt.title("%s vs %s" % (name1,name2))
    if not os.path.isdir('outputs'):
        os.mkdir('outputs')
    name = 'cross_corr_%s_%s' % (name1,name2)
    fn = os.path.join('outputs',name+'.png')
    plt.savefig(fn,dpi=600)
    #fn = os.path.join('outputs',name+'.pdf')
    #plt.savefig(fn)
    #fn = os.path.join('outputs',name+'.eps')
    #plt.savefig(fn)
    print("%s vs %s <cross_correlation>:" % (name1,name2),np.average(cross_corr))
    plt.clf()
    fig,ax = plt.subplots(layout='constrained')
    ax.scatter(white_data1,white_data2,marker='.',s=0.1)
    ax.set_xlabel(name1)
    ax.set_ylabel(name2)
    plt.tight_layout()
    name = 'comparison_%s_%s' % (name1,name2)
    fn = os.path.join('outputs',name+'.png')
    plt.savefig(fn,dpi=600)
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f','--file', type=str, default='output.h5',
        help='Path to ABF integrated output')
    parser.add_argument('-m','--max', type=float, default=40.0,
        help='Maximum free energy to plot up to. Default = 40.0.')
    args = parser.parse_args()

    main()

