#Module plot_cvs
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
    """Compare the kdes of 3 plots on the same graph."""
    #Load pandas dataframes
    #dirs = ["soft_047/0","soft_047/1","soft_047/2","soft_048/0","soft_048/1",
    #    "soft_048/2","soft_049/0","soft_049/1","soft_049/2"]
    dirs = ["soft_047/0","soft_048/0","soft_049/0"]
    #dirs = ["soft_038/0","soft_048/0","soft_049/0"]
    cvs = "fiber_angle"

    fns = [os.path.join(ele,"log.dframe") for ele in dirs]
    dfs = [pd.read_csv(fn) for fn in fns]
    #dfs_real = [pd.concat((dfs[0],dfs[1],dfs[2])),pd.concat((dfs[3],dfs[4],
    #        dfs[5])),pd.concat((dfs[6],dfs[7],dfs[8]))]

    compare_plots(dfs,cvs)

    print("Created all correlation plots. Exiting...")

def compare_plots(dfs,name):
    '''Create a joint kde plot for all data.

    Parameters:
        dfs: *list of pandas.dframe*
            Dataframes containing the data to plot.
        name: *str*
            Name of the column
    '''
    
    #Plot zfit KDE
    plt.clf()
    font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : 16}

    matplotlib.rc('font', **font)
    fig,ax = plt.subplots(layout='constrained')
    bounds = [40,180]
    kde_min = bounds[0]
    kde_max = bounds[1]
    obs = zfit.Space([name], limits=(kde_min, kde_max))
    data = [df[name].to_numpy()*180.0/np.pi for df in dfs]

    #Create kdes
    kde_d20c0 = zfit.pdf.KDE1DimISJ(data[0], obs=obs)
    kde_d0c0 = zfit.pdf.KDE1DimISJ(data[1], obs=obs)
    kde_d20a5c0 = zfit.pdf.KDE1DimISJ(data[2], obs=obs)
    x = np.linspace(kde_min,kde_max,1000)
    colors = ["#E69F00","#56B4E9","#009E73"]
    kde_d20c0_vals = kde_d20c0.pdf(x)
    kde_d0c0_vals = kde_d0c0.pdf(x)
    kde_d20a5c0_vals = kde_d20a5c0.pdf(x)
    
    kde_vals = [kde_d20c0_vals,kde_d0c0_vals,kde_d20a5c0_vals]
    averages = [np.sum(kde_vals[i]*x)/np.sum(kde_vals[i]) for i in range(len(dfs))]
    #averages = [np.average(ele) for ele in kde_vals]
    x_inds = [np.argmax(x > kde_ave) for kde_ave in averages]
    print("x_inds:",x_inds)
    y_ends = [kde_vals[i][x_inds[i]] for i in range(len(dfs))]
    

    #Actually add to plot
    ax.plot(x, kde_d20c0_vals, linewidth=0.5, color=colors[0],label='D20C0')
    ax.fill_between(x, kde_d20c0_vals, color=colors[0], alpha=0.25)
    ax.vlines(averages[0],ymin=0.0,ymax=float(y_ends[0]),linestyle='--',linewidth=0.5,
            color = colors[0])

    ax.plot(x, kde_d0c0_vals, linewidth=0.5, color=colors[1],label='D0C0')
    ax.fill_between(x, kde_d0c0_vals, color=colors[1], alpha=0.25)
    ax.vlines(averages[1],ymin=0.0,ymax=float(y_ends[1]),linestyle='--',linewidth=0.5,
            color = colors[1])

    ax.plot(x, kde_d20a5c0_vals, linewidth=0.5, color=colors[2],label='D20A5C0')
    ax.fill_between(x, kde_d20a5c0_vals, color=colors[2], alpha=0.25)
    ax.vlines(averages[2],ymin=0.0,ymax=float(y_ends[2]),linestyle='--',linewidth=0.5,
            color = colors[2])
    ax.legend(loc='best',edgecolor='grey')
    ax.set_xlabel("Fiber angle (deg)")
    ax.set_ylabel("Marginal probability")
    ax.set_xlim(90,180)
    fn = os.path.join('outputs',name+"_kde.png")
    fig.savefig(fn,dpi=600)

    fn = os.path.join('outputs',name+"_kde.svg")
    fig.savefig(fn)


    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f','--file', type=str, default='output.h5',
        help='Path to ABF integrated output')
    parser.add_argument('-m','--max', type=float, default=40.0,
        help='Maximum free energy to plot up to. Default = 40.0.')
    args = parser.parse_args()

    main()

