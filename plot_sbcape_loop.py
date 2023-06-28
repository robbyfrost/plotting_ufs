# --------------------------------
# Name: plot_sbcape_loop.py
# Author: Robert M. Frost
# NOAA Global Systems Laboratory
# Created: 23 June 2023
# Purpose: Loop to plot surface based
# cape comparisons at different times 
# during forecast runs
# --------------------------------
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cpf
from datetime import datetime
from matplotlib import rc
import seaborn
import numpy as np
# --------------------------------
for hr in range(7, 37):
    # read in rap and hrrr data
    if hr < 10:
        dnc_h = f"/scratch2/BMC/fv3lam/Robby.Frost/expt_dirs/2019052000_3km_hrrrphys/2019052000/postprd/rrfs.t00z.natlev.f00{hr}.rrfs_conuscompact_3km.grib2"
        dnc_r = f"/scratch2/BMC/fv3lam/Robby.Frost/expt_dirs/2019052000_3km_rapphys/2019052000/postprd/rrfs.t00z.natlev.f00{hr}.rrfs_conuscompact_3km.grib2"
    else:
        dnc_h = f"/scratch2/BMC/fv3lam/Robby.Frost/expt_dirs/2019052000_3km_hrrrphys/2019052000/postprd/rrfs.t00z.natlev.f0{hr}.rrfs_conuscompact_3km.grib2"
        dnc_r = f"/scratch2/BMC/fv3lam/Robby.Frost/expt_dirs/2019052000_3km_rapphys/2019052000/postprd/rrfs.t00z.natlev.f0{hr}.rrfs_conuscompact_3km.grib2"
    hrrr = xr.open_dataset(dnc_h,
                        engine="cfgrib",
                        filter_by_keys={'stepType': 'instant', 'typeOfLevel': 'surface'})
    rap = xr.open_dataset(dnc_r,
                        engine="cfgrib",
                        filter_by_keys={'stepType': 'instant', 'typeOfLevel': 'surface'})
    
    # set lats and lons
    lat = hrrr.latitude
    lon = hrrr.longitude

    # take numpy.datatime64 and put into datetime format
    valid_time_str = str(hrrr.time.values + hrrr.time.step.values)
    valid_time_str = valid_time_str[:19]
    valid_time = datetime.strptime(valid_time_str, '%Y-%m-%dT%H:%M:%S')

    # plotting setup
    rc('font',weight='normal',size=15)#,family='serif',serif='Times New Roman')
    # rc('text',usetex='True')
    rc('figure',facecolor='white')

    # --------------------------------
    # side by side sbcape plot

    # Define your custom colorbar bounds
    cbar_min = 100
    cbar_max = 5000
    clevs = np.arange(cbar_min, cbar_max, 100)
    # color palette
    colors = "CMRmap_r"

    # create plot
    fig, ax = plt.subplots(ncols=2, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(18,10))
    # mapping
    plt_area = [-104, -94, 30, 39] # W, E, S, N
    for i, iax in enumerate(ax):
        iax.coastlines()
        iax.add_feature(cpf.BORDERS)
        iax.add_feature(cpf.STATES)
        iax.set_extent(plt_area)
    # plot
    c0 = ax[0].contourf(lon, lat, hrrr.cape, clevs, 
                        transform=ccrs.PlateCarree(), 
                        cmap=colors,
                        extend="both")
    c1 = ax[1].contourf(lon, lat, rap.cape, clevs,
                        transform=ccrs.PlateCarree(), 
                        cmap=colors,
                        extend="both")
    # pretty up
    ax[0].set_title(f"HRRR F0{hr},  Valid {valid_time} UTC")
    ax[1].set_title(f"RAP F0{hr},  Valid {valid_time} UTC")
    # Add colorbar
    cbar = fig.colorbar(c1, ax=ax, orientation='horizontal', extend=True, pad=0.05, aspect=50)
    cbar.set_label('SBCAPE [J kg$^{-1}$]')
    # save and close figure
    plt.savefig(f"/scratch2/BMC/fv3lam/Robby.Frost/figures/20190520/sbcape_sidebyside_f{hr}.png")
    plt.close()
     # --------------------------------
    # sbcape difference plot

    # Define your custom colorbar bounds
    cbar_min = -3500
    cbar_max = 3500
    clevs = np.arange(cbar_min, cbar_max, 100)

    # Create a normalized colormap
    cdiff = seaborn.color_palette("seismic", as_cmap=True)

    # create plot
    fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(10,10))

    # mapping
    plt_area = [-104, -94, 30, 39] # W, E, S, N
    ax.coastlines()
    ax.add_feature(cpf.BORDERS)
    ax.add_feature(cpf.STATES)
    ax.set_extent(plt_area)
    # plot
    c0 = ax.contourf(lon, lat, hrrr.cape - rap.cape, clevs, 
                    transform=ccrs.PlateCarree(), 
                    cmap=cdiff, 
                    vmin=cbar_min, vax=cbar_max,
                    extend="both")

    # pretty up
    ax.set_title(f"HRRR - RAP F0{hr},  Valid {valid_time} UTC")

    # Add colorbar
    cbar = fig.colorbar(c0, ax=ax, orientation='horizontal', extend=True, pad=0.05, aspect=50)
    cbar.set_label('SBCAPE [J kg$^{-1}$]')

    plt.savefig(f"/scratch2/BMC/fv3lam/Robby.Frost/figures/20190520/sbcape_diff_f{hr}.png")
    plt.close()