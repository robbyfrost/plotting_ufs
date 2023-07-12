# --------------------------------
# Name: plot_cin_loop.py
# Author: Robert M. Frost
# NOAA Global Systems Laboratory
# Created: 03 July 2023
# Purpose: Plot CIN comparisons of 
# SRW output
# --------------------------------
import matplotlib.pyplot as plt
from matplotlib import rc
import cartopy.crs as ccrs
import cartopy.feature as cpf
import numpy as np
from metpy.plots import ctables
from UFSutils import read_grib
import geopandas as gpd
import seaborn
# --------------------------------
# Important parameters

# directory where hrrr grib data are located
dgrib_h = "/scratch2/BMC/fv3lam/Robby.Frost/expt_dirs/2023041900_3km_hrrrphys/2023041900/postprd/"
# directory where rap grib data are located
dgrib_r = "/scratch2/BMC/fv3lam/Robby.Frost/expt_dirs/2023041900_3km_rapphys/2023041900/postprd/"
# natlev or prslev
nat_prs = "natlev"
# message number of SBCIN
mn_sbcin = 1409
# directory for figure to be output
figdir = "/scratch2/BMC/fv3lam/Robby.Frost/figures/20230419/sbcin/"
# --------------------------------
# Read in SBCIN

# loop over time
for hr in range(37):
    print(f"Hour {hr}")
    hrrr, sbcin_h, lat, lon, valid_date = read_grib(hr, dgrib_h, nat_prs, mn_sbcin)
    sbcin_r = read_grib(hr, dgrib_r, nat_prs, mn_sbcin, ret_type=1)
    # --------------------------------
    # plotting setup
    rc('font',weight='normal',size=12.5)
    # rc('text',usetex='True')
    rc('figure',facecolor='white')
    # --------------------------------
    # Side by side SBCIN plot
    print(f"Creating 1 x 2 SBCIN Plot")

    # Define your custom colorbar bounds
    cbar_min = -1000
    cbar_max = -9
    # levels for sbcape to be plot
    clevs = np.arange(cbar_min, cbar_max, 10)
    # color palette
    colors = "CMRmap"

    # create plot
    fig, ax = plt.subplots(ncols=2, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(16,6), constrained_layout=True)

    # plot HRRR
    c0 = ax[0].contourf(lon, lat, sbcin_h.values, clevs, 
                        transform=ccrs.PlateCarree(), 
                        cmap=colors, extend="both")
    # plot RAP
    c1 = ax[1].contourf(lon, lat, sbcin_r.values, clevs, 
                        transform=ccrs.PlateCarree(), 
                        cmap=colors, extend="both")

    # mapping
    plt_area = [-101, -94, 33.5, 37.5] # W, E, S, N
    for i, iax in enumerate(ax):
        iax.coastlines()
        iax.add_feature(cpf.BORDERS)
        iax.add_feature(cpf.STATES)
        iax.set_extent(plt_area)
        # Load the json file with county coordinates
        geoData = gpd.read_file('https://raw.githubusercontent.com/holtzy/The-Python-Graph-Gallery/master/static/data/US-counties.geojson')
        geoData.plot(ax=iax, color="none", lw=0.3, aspect=1)

    # set title
    ax[0].set_title(f"HRRR F0{hr},  Valid {valid_date} UTC")
    ax[1].set_title(f"RAP F0{hr},  Valid {valid_date} UTC")

    # Add colorbar
    cbar = fig.colorbar(c1, ax=ax, orientation='horizontal', extend=True, pad=0.03, aspect=50)
    cbar.set_label('Surface Based CIN [J kg$^{-1}$]')
    cbar.set_ticks(np.arange(cbar_min,0.1,100))

    # save and close figure
    figdir_full = f"{figdir}sbcin_sidebyside_f{hr}.png"
    print(f"Saving figure to {figdir_full}")
    plt.savefig(figdir_full)
    plt.close()
    print("Finished plotting 1 x 2 SBCIN!")
    # --------------------------------
    # SBCIN Difference Plot
    print("Creating SBCIN Difference Plot")

    # Define your custom colorbar bounds
    cbar_min = -500
    cbar_max = 500.1
    clevs = np.linspace(cbar_min, cbar_max, 51)
    # color palette
    colors = seaborn.color_palette("seismic", as_cmap=True)

    # create plot
    fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()}, 
                        figsize=(10,6.5), constrained_layout=True)

    # plot HRRR - RAP
    c0 = ax.contourf(lon, lat, sbcin_h.values - sbcin_r.values,
                    clevs, transform=ccrs.PlateCarree(),
                    cmap=colors, extend="both")

    # mapping
    plt_area = [-101, -94, 33.5, 37.5] # W, E, S, N
    ax.coastlines()
    ax.add_feature(cpf.BORDERS)
    ax.add_feature(cpf.STATES)
    ax.set_extent(plt_area)
    # Load the json file with county coordinates
    geoData = gpd.read_file('https://raw.githubusercontent.com/holtzy/The-Python-Graph-Gallery/master/static/data/US-counties.geojson')
    geoData.plot(ax=ax, color="none", lw=0.3, aspect=1)

    # set title
    ax.set_title(f"HRRR - RAP F0{hr},  Valid {valid_date} UTC")

    # Add colorbar
    cbar = fig.colorbar(c0, ax=ax, orientation='horizontal', extend=True, pad=0.03, aspect=50)
    cbar.set_label('Surface Based CIN [J kg$^{-1}$]')
    cbar.set_ticks(np.arange(cbar_min, cbar_max, 250))

    # save and close figure
    figdir_full = f"{figdir}sbcin_diff_f{hr}.png"
    print(f"Saving figure to {figdir_full}")
    plt.savefig(figdir_full)
    plt.close()
    print(f"Finished with hour {hr}! \n")