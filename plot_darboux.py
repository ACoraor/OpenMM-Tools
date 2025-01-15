#Module plot_darboux
#Created by Aria Coraor
#Written 9/13/24

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
import KDEpy
import matplotlib.transforms as transforms
import cartopy.crs as ccrs
from scipy.interpolate import RegularGridInterpolator as RGI

def main():
    """Plot each cv versus time"""
    #Load pandas dataframe
    df = pd.read_csv('log.dframe')

    if not os.path.isdir('outputs'):
        os.mkdir('outputs')

    #Find columns and plot
    cols = df.columns
    darb_cols = [col for col in cols if col[-3:] in set(["_t1","_t2","_t3"])]
    darboux_plot(df,darb_cols)

        
    print("Created globe plot. Exiting...")


def darboux_plot(df,names):
    '''Create a 2d-mapped plot of the probability density of the Darboux
    vector using Cartopy.

    Parameters:
        df: *pd.Dataframe*
                Full log data.
        names: *list of str*
            Names of all Darboux columns
    '''

    #Extract the names and number of each dye
    dye_names = list(set([name[:-3] for name in names]))
    dye_names.sort()

    #Iterate over each separate dye name and create plot
    for i, dye in enumerate(dye_names):

        #Calculate phi, theta from t1 / t2 / t3 components
        t1 = df["%s_t1" % dye].to_numpy()
        t2 = df["%s_t2" % dye].to_numpy()
        t3 = df["%s_t3" % dye].to_numpy()
        phis = np.arccos(t3) # Azimuthal angle, radians, 0-Pi
        thetas = np.arctan2(t2,t1) #Longitudinal angle, 0-2Pi
        
        #Transform to plate carree, degrees
        thetas *= 180.0/np.pi
        phis = 90 - (phis*180.0/np.pi)

        print("Theta range:",[np.min(thetas),np.max(thetas)])
        print("Phis range:",[np.min(phis),np.max(phis)])


        #Create numpy histogram for pcolormesh
        xedges = np.linspace(-180.0,180.0,91)
        yedges = np.linspace(-90.0,90.0,91)
        xcent = (xedges[:-1] + xedges[1:])/2
        ycent = (yedges[:-1] + yedges[1:])/2
        #If 2d KDEs are working, use that:
        try:    
            H, xed, yed = calc_darboux_kde(t1,t2,t3,name=dye)
        except Exception as e:
            print("Failed in kde:",e)
            
            print("Continuing with standard histogram.")
            H, xed, yed = np.histogram2d(thetas,phis,bins=[xedges,yedges])

        #Finally, plot:
        plt.clf()
        fig = plt.figure()
        gs = fig.add_gridspec(2,2,width_ratios = (4,1),height_ratios = (1,4))
        ax1 = fig.add_subplot(gs[0,0])
        ax2 = fig.add_subplot(gs[0,1])
        ax3 = fig.add_subplot(gs[1,0],projection=ccrs.EckertIV())
        ax4 = fig.add_subplot(gs[1,1])
        '''
        fig, axs = plt.subplots(
            nrows=2,
            ncols=2,
            sharex="col",
            sharey="row",
            layout="constrained",
            width_ratios=(4, 1),
            height_ratios=(1, 4),
            subplot_kw = {'projection': 'rectilinear'}
        )
        ax1 = axs[0, 0]
        ax2 = axs[0, 1]
        ax3 = axs[1, 0]
        ax4 = axs[1, 1]
        '''
        ax2.set_visible(False)
        # Define matplotlib transforms
        base = ax4.transData
        scale = transforms.Affine2D().scale(sx=-1, sy=1)
        rot = scale.rotate_deg(270)
        
        #ax3 = 
        #ax3.projection=ccrs.EckertIV()
        ax3.pcolormesh(xcent,ycent,H,cmap='viridis',shading='auto',
            transform=ccrs.PlateCarree())
        ax3.gridlines(linestyle=':')
        #fig.savefig('%s_darboux.png'% dye,dpi=600)
        #fig.savefig('%s_darboux.svg'% dye)
        #fig.savefig('%s_darboux.pdf'% dye)
        #Pasted from joined_pca:

def calc_darboux_kde(t1,t2,t3,gridpoints=512,name='darboux_scat'):
    '''Create a proper KDE of the darboux vectors in PlateCarRee coordinates,

    Parameters:
        t1: *np.array*, shape: (n_ts)
                Components of the orientation in the t1 direction (Stiffest).
        t2: *np.array*, shape: (n_ts)
                Components of the orientation in the t2 direction. (Floppy)
        t3: *np.array*, shape: (n_ts)
                Components of the orientation in the t3 direction. (Principal)
        gridpoints: *int*
            Number of gridpoints for the evaluation of the KDE.
        name: *str*
            Name for saving 2d plot.

    Returns:
        H: *2d np.array*, shape: (gridpoints,gridpoints)
                KDE density in theta-phi space.
        xcent: *np.array*, shape: (gridpoints)
                Center of 'histogram' bins for kde density evaluations, 
                in longitudinal degrees.
        ycent: *np.array*, shape: (gridpoints)
                Center of 'histogram' bins for kde density evaluations, 
                in latitudinal degrees.
    '''

    #If plotting planar_35, switch to planar_53
    if name[6:9] == "p35":
        name = name[:6] + "p53" + name[9:]
        t1 = -t1
        t2 = -t2
        t3 = -t3
    
    #Perform the KDE in Darboux space, then map the pdf to PlateCaree space,
        #Then project onto EckertIV space
    #phi = zfit.Space("Phi", lower = -90, upper = 90)
    #theta = zfit.Space("Theta", lower = -180, upper = 180)
    t1_obs = zfit.Space("t1", lower = -1, upper = 1)
    t2_obs = zfit.Space("Theta", lower = -1, upper = 1)
    t3_obs = zfit.Space("Theta", lower = -1, upper = 1)

    print("Fitting kdes...")
    #phi_data = zfit.Data(phi,obs=phi
    kde_t1 = zfit.pdf.KDE1DimISJ(t1, obs=t1_obs)
    kde_t2 = zfit.pdf.KDE1DimISJ(t2, obs=t2_obs)
    kde_t3 = zfit.pdf.KDE1DimISJ(t3, obs=t3_obs)
    #kde_phi = zfit.pdf.KDE1DimISJ(phi, obs=phi)
    #kde_theta = zfit.pdf.KDE1DimISJ(theta, obs=theta)
    print("KDEs fit.")
    #H = kde_eval()
            

    print("Performing fft kdes.")
    #Calculate geometric mean bandwidth
    geom_mean_bandwidth = (float(kde_t1._bandwidth.numpy()) *
            float(kde_t2._bandwidth.numpy())*
            float(kde_t3._bandwidth.numpy()))**(1/3)
    

    #Convert joined_pca code into 3d KDE
    fft_kde = KDEpy.FFTKDE(kernel='gaussian',
            bw=geom_mean_bandwidth)
    data_all = np.column_stack((t1.ravel(),t2.ravel(),t3.ravel()))
    fft_kde.fit(data_all)
    
    #Calculate density for 1v2 and 1v3
    fft_points = gridpoints
    gridpoints, density = fft_kde.evaluate(fft_points)
    density = density / np.max(density) # normalize

    #Convert to RGI formatting
    grid_x = np.sort(np.unique(gridpoints[:,0]))
    grid_y = np.sort(np.unique(gridpoints[:,1]))
    grid_z = np.sort(np.unique(gridpoints[:,2]))
    density_mesh = np.reshape(density,(fft_points,fft_points,fft_points))
    print("grid_x.shape:",grid_x.shape)

    #Density mesh: normalized pdf in t1, t2, t3 space!
    
    #Calculate scatter color for all points in data
    interpolator = RGI((grid_x,grid_y,grid_z),density_mesh)
    colors = interpolator(data_all)
    #colors_1v3 = interpolator_1v3(data_1v3)

    #Now: from density_mesh, calculate phi and theta!
    #phi_theta_mesh = np.zeros(shape=(fft_points,fft_points))
    #Using formulas for the pdf transformation from <https://tinyurl.com/bdjkpjsk>
    phis = np.arccos(data_all[:,2]) # Azimuthal angle, radians, 0-Pi
    thetas = np.arctan2(data_all[:,1],data_all[:,0]) #Longitudinal angle, 0-2Pi
    
    #Transform to plate carree, degrees
    thetas *= 180.0/np.pi
    phis = 90 - (phis*180.0/np.pi)

    #Create scatterplot on ax3
    plt.clf()
    fig = plt.figure()
    gs = fig.add_gridspec(7,2,width_ratios = (4,1),
        height_ratios = (6,4,4,4,4,4,4))
    ax1 = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[0,1])
    ax3 = fig.add_subplot(gs[1:,0],projection=ccrs.EckertIV())
    ax4 = fig.add_subplot(gs[2:-1,1])

    ax2.set_visible(False)
    base = ax4.transData
    scale = transforms.Affine2D().scale(sx=-1, sy=1)
    rot = scale.rotate_deg(270)
    cmap = cm.viridis
    
    #Sort so that brightest colors are scattered last
    color_shape = colors.shape
    color_sorter = np.argsort(colors.ravel()).ravel()
    colors_sorted = colors.ravel()[color_sorter]
    thetas_sorted = thetas.ravel()[color_sorter]
    phis_sorted = phis.ravel()[color_sorter]

    ax3.scatter(thetas_sorted,phis_sorted, marker=".", s=1, c=cmap(colors_sorted),
        transform=ccrs.PlateCarree())
    ax3.gridlines(linestyle=':')


    #ax3.set_xlabel("PC 1")
    #ax3.set_ylabel("PC 2")
    
    #bounds_1v2 = (-4.0,7.0)
    #bounds_1v3 = (-3.0, 4.0)
    ax3.set_extent([-180,180,-90,90])
    #ax3.set_ylim(-90,90)
    #ax3.set_xlim(-180,180)

    #If pandas dataframe already created, scatter the modes:
    if os.path.isfile("aves.csv"):
        aves = pd.read_csv('aves.csv')
        #Get p3, p5, and p35 axes of this name
        

        #Add these axes to the scatter

        #Calculate plane in the great circle using the formulas at
        # <https://tinyurl.com/ycxsmw95>

    
    #Plot other kdes
    phi_sp = zfit.Space("Phi", lower = -90, upper = 90)
    theta_sp = zfit.Space("Theta", lower = -180, upper = 180)
    kde_phi = zfit.pdf.KDE1DimISJ(phis, obs=phi_sp) #True azimuthal angle
    kde_theta = zfit.pdf.KDE1DimISJ(thetas, obs=theta_sp)
    #xlim = ax3.get_xlim()
    #ylim = ax3.get_ylim()
    x = np.linspace(-180.0,180.0,1000)
    y_eval = np.linspace(-90.0,90.0,1000)
    kde_horiz = kde_theta.pdf(x)
    kde_vert = kde_phi.pdf(y_eval)
    y_plot = np.linspace(0.0,180.0,1000)

    lw = 0.5
    color = '#0072B2'
    ax1.plot(x, kde_horiz, linewidth=lw, color=color)
    ax1.fill_between(x, kde_horiz, color=color, alpha=0.25)
    ax1.set_xlim((-180.0,180.0))
    theta_labels = ["-180","-120","-60","0","60","120","180"]
    ax1.set_xticks([float(ele) for ele in theta_labels])
    ax1.set_xticklabels(theta_labels)
    #ax4.set_xlim((0.0,180.0))
    ax4.plot(y_plot, kde_vert, linewidth=lw, color=color, 
        transform = rot+base)
    xlim = ax4.get_xlim()
    ax4.fill_between(y_plot, kde_vert, color=color, alpha=0.25, 
            transform  = rot+base)
    ax4.set_xlim(xlim)
    ax4.set_ylim(0.0,180.0)
    phi_labels = ["180","150","120","90","60","30","0"]
    print("xlim:",xlim)
    print("ylim:",ax4.get_ylim())
    ax4.set_yticks([float(ele) for ele in phi_labels])
    ax4.set_yticklabels([ele for ele in reversed(phi_labels)])
    #ax4.set_xlim(0,np.max(kde_vert)*1.1)

    #Save modes to dataframe
    if os.path.isfile("aves.csv"):
        df = pd.read_csv('aves.csv')
        cols = df.columns
    else:
        df = pd.DataFrame()
        cols = ['name','mode_theta','ave_theta','mode_phi','ave_phi']
        #for col in cols:
        #    df[col] = None
        #df.columns = cols
    mode = np.argmax(colors_sorted).flatten()
    mode_theta = thetas_sorted[mode]
    mode_phi = phis_sorted[mode]
    ave_theta = np.average(thetas_sorted)
    ave_phi = np.average(phis_sorted)
    #Reformat for ease of DataFrame entry
    data_list = [name,mode_theta,ave_theta,mode_phi,ave_phi]
    data_dict = {col:ele for col,ele in zip(cols,data_list)}
    if len(df.columns) == 0:
        df = pd.DataFrame(data_dict)
    else:
        df.loc[-1] = [name,mode_theta,ave_theta,mode_phi,ave_phi]
    #df.columns = cols
    df.to_csv('aves.csv',index=False,na_rep='NULL')
    
    
    base = "outputs/" + name + "_scatt"
    plt.savefig(base + ".png", dpi=600)
    plt.savefig(base + ".pdf")
    plt.savefig(base + ".svg")
    plt.savefig(base + ".eps")
    print("Saved %s.png" % base)

    return #Haven't calculated these...

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    args = parser.parse_args()

    main()

