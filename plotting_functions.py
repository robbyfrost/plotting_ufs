# --------------------------------
# Name: plotting_functions.py
# Author: Robert M. Frost
# NOAA Global Systems Laboratory
# Created: 16 June 2023
# Purpose: A collection of functions for
# plotting output from the UFS short range
# weather application.
# --------------------------------
import matplotlib.pyplot as plt
import pygrib
from metpy.plots import ctables
from matplotlib import rc
import cartopy.crs as ccrs
import cartopy.feature as cpf

def read_grib(hr, dgrib, nat_prs, mesg_num, array_only=False):
    """
    Purpose: Function to read in grib output from the UFS SRW app.
    :param int hr: forecast hour to be read
    :param str dgrib: directory where forecast output are located
    :param str nat_prs: reading in natlev or prslev
    :param int mesg_num: message number for the variable to be read in
    :param array_only bool: returns just variable array if true, returns variable array, lat lon, and valid_time if false
    """
    # set filenames
    if hr < 10:
        dgrib_h = f"{dgrib}rrfs.t00z.{nat_prs}.f00{hr}.rrfs_conuscompact_3km.grib2"
    else:
        dgrib_h = f"{dgrib}rrfs.t00z.{nat_prs}.f0{hr}.rrfs_conuscompact_3km.grib2"

    # open hrrr and rap output
    print("Reading in grib output")
    grbs = pygrib.open(dgrib_h)
    # extract variable of interest
    grb = grbs[mesg_num]
    # extract latitude and longitude arrays
    lat, lon = grb.latlons()
    # extract datatime
    valid_time = grb.validDate
    
    if array_only:
        return grb
    else:
        return grb, lat, lon, valid_time
# --------------------------------
def plot_precip_acc(hr, clevs, lon, lat, tot_accu1, tot_accu2, 
                    nonc_accu1, nonc_accu2, figdir, valid_date,
                    plt_area=[-104, -94, 30, 39]):
    """
    Adapted from plot_precip_accu.py
    Purpose: Plots comparison of total precipitation accumilation between two model runs and cumulus schemes
    :param int hr: forecast hour being plot
    :param list or array clevs: Determines the number and positions of the contour regions
    :param array lon: array containing longitude values for plotting
    :param array lat: array containing latitude values for plotting
    :param array tot_accu1: first array of total precipitation accumulation
    :param array tot_accu2: second array of total precipitation accumulation
    :param array nonc_accu1: first array of non-convective precipitation accumulation
    :param array nonc_accu2: second array of non-convective precipitation accumulation
    :param str figdir: directory for figure to be saved
    :param datetime valid_date: valid date and time that corresponds to forecast hour
    :param list plt_area: mapping bounds in longitude and latitude [W,E,S,N] (centered over southern plains by default)
    """
    # set font size
    rc("font",weight="normal",size=15)
    # precipitation colortable
    colors = ctables.registry.get_colortable('precipitation')

    # create plot
    fig, ax = plt.subplots(ncols=2, nrows=2, subplot_kw={'projection': ccrs.PlateCarree()}, 
                        figsize=(18,16), constrained_layout=True)

    # mapping
    for i, iax in enumerate(ax[:,0]):
        iax.coastlines()
        iax.add_feature(cpf.BORDERS)
        iax.add_feature(cpf.STATES)
        iax.set_extent(plt_area)
    for i, iax in enumerate(ax[:,1]):
        iax.coastlines()
        iax.add_feature(cpf.BORDERS)
        iax.add_feature(cpf.STATES)
        iax.set_extent(plt_area)
        
    # plot
    print(f"Plotting F0{hr}")
    c0 = ax[0,0].contourf(lon, lat, tot_accu1.values, 
                        clevs, transform=ccrs.PlateCarree(), 
                        cmap=colors, extend='max')
    c1 = ax[0,1].contourf(lon, lat, tot_accu2.values, 
                        clevs, transform=ccrs.PlateCarree(), 
                        cmap=colors, extend='max')
    c2 = ax[1,0].contourf(lon, lat, tot_accu1.values - nonc_accu1.values, 
                        clevs, transform=ccrs.PlateCarree(), 
                        cmap=colors, extend='max')
    c3 = ax[1,1].contourf(lon, lat, tot_accu2.values - nonc_accu2.values, 
                        clevs, transform=ccrs.PlateCarree(), 
                        cmap=colors, extend='max')

    # pretty up
    ax[0,0].set_title(f"HRRR F0{hr},  Valid {valid_date} UTC")
    ax[0,1].set_title(f"RAP F0{hr},  Valid {valid_date} UTC")
    ax[1,0].set_title(f"HRRR F0{hr}, Valid {valid_date} UTC")
    ax[1,1].set_title(f"RAP F0{hr}, Valid {valid_date} UTC")

    # Add colorbar
    cbar = fig.colorbar(c0, ax=ax[0,:], orientation='vertical', 
                        extend=True, pad=0.02, fraction=0.013, 
                        aspect=30)
    cbar.set_label('Total Precipitation Accumilation [mm]')
    cbar2 = fig.colorbar(c3, ax=ax[1,:], orientation='vertical', 
                        extend=True, pad=0.02, fraction=0.013, 
                        aspect=30)
    cbar2.set_label('Un-Resolved Precipitation Accumilation [mm]')

    # save fig
    dout = f"{figdir}precip_accum_f{hr}.png"
    plt.savefig(dout)
    plt.close()
    print(f"Figure saved to {dout}")

    return