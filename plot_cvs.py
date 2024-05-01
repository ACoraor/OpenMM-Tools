#Module plot_cvs
#Created by Aria Coraor
#Written 4/30/24

import numpy as np
import mdtraj as md
import argparse
import matplotlib
matplotlib.use("Agg")
from matplotlib import rc
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pandas as pd

def main():
    """Plot each cv versus time"""
    #Load pandas dataframe
    df = pd.read_csv('log.dframe')
    cvs = df.columns[3:]

    #Plot all single cvs
    for cv in cvs:
        single_plot(df['timestep'],df[cv],cv)

    #Plot all pairs of cvs
    for i, cv1 in enumerate(cvs[:-1]):
        for j, cv2 in enumerate(cvs[i+1:]):
            double_plot(df['timestep'],df[cv1],df[cv2],cv1,cv2)
    
    print("Created all correlation plots. Exiting...")


def single_plot(timestep,data,name):
    '''Create a single plot for the cv.

    Parameters:
        timestep: *arraylike*, shape: (n_ts)
                Timestep for each datapoint.
        data: *arraylike*, shape: (n_ts)
            The column of data to be plotted.
        name: *str*
            Name of the column
    '''
    
    plt.clf()
    ts_sec = 2e-15 # Timestep in seconds
    plot_times = timestep * ts_sec / (1e-9) #Plot in nanoseconds
    

    fig,ax = plt.subplots(layout='constrained')
    ax.plot(timestep,data,color='blue')
    
    
    ax.set_xlabel("Simulation time (ns)")
    ax.set_ylabel(name) #soft_041
    try:
        plt.tight_layout()
    except:
        print("Error: cannot do tight_layout.")
    if not os.path.isdir('outputs'):
        os.mkdir('outputs')
    fn = os.path.join('outputs',name+'.png')
    plt.savefig(fn,dpi=600)
    fn = os.path.join('outputs',name+'.pdf')
    plt.savefig(fn)
    fn = os.path.join('outputs',name+'.eps')
    plt.savefig(fn)

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
    ts_sec = 2e-15 # Timestep in seconds
    plot_times = timestep * ts_sec / (1e-9) #Plot in nanoseconds
    
    fig,ax = plt.subplots(layout='constrained')

    ax.plot(timestep,data1,color='blue',label=name1)
    secax = ax.secondary_yaxis('right')
    secax.plot(timestep,data2,color='red',label=name2)
    
    ax.set_xlabel("Simulation time (ns)")
    ax.set_ylabel(name1) #soft_041
    secax.set_ylabel(name2)
    if not os.path.isdir('outputs'):
        os.mkdir('outputs')
    name = '%s_vs_%s' % (name1,name2)
    plt.legend(loc='best',edgecolor='grey')
    fn = os.path.join('outputs',name+'.png')
    plt.savefig(fn,dpi=600)
    fn = os.path.join('outputs',name+'.pdf')
    plt.savefig(fn)
    fn = os.path.join('outputs',name+'.eps')
    plt.savefig(fn)

    #Calculate cross-correlation
    white_data1 = (data1 - np.average(data1))/np.std(data1)
    white_data2 = (data2 - np.average(data2))/np.std(data2)
    cross_corr = white_data1 * white_data2
    plt.clf()
    fig,ax = plt.subplots(layout='constrained')

    ax.plot(timestep,cross_corr,color='blue')
    ax.set_xlabel("Simulation time (ns)")
    ax.set_ylabel('Cross Correlation') #soft_041
    ax.title("%s vs %s" % (name1,name2))
    if not os.path.isdir('outputs'):
        os.mkdir('outputs')
    name = 'cross_corr_%s_%s' % (name1,name2)
    fn = os.path.join('outputs',name+'.png')
    plt.savefig(fn,dpi=600)
    fn = os.path.join('outputs',name+'.pdf')
    plt.savefig(fn)
    fn = os.path.join('outputs',name+'.eps')
    plt.savefig(fn)
    print("%s vs %s <cross_correlation>:" % (name1,name2),np.average(cross_corr))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f','--file', type=str, default='output.h5',
        help='Path to ABF integrated output')
    parser.add_argument('-m','--max', type=float, default=40.0,
        help='Maximum free energy to plot up to. Default = 40.0.')
    args = parser.parse_args()

    main()

