# --------------------------------
# Name: plot_sbcape.ipynb
# Author: Robert M. Frost
# NOAA Global Systems Laboratory
# Created: 03 July 2023
# Purpose: Plot surface based cape 
# utilizing read_grib from 
# plotting_functions.py
# --------------------------------
import matplotlib.pyplot as plt
from matplotlib import rc
import cartopy.crs as ccrs
import cartopy.feature as cpf
import numpy as np
from metpy.plots import ctables
from plotting_functions import read_grib
import geopandas as gpd
import seaborn
# --------------------------------
# important parameters

# first and last hour to be plotted
hri = 0
hrf = 37 # last hour + 1
# directory where hrrr grib data are located
dgrib_h = "/scratch2/BMC/fv3lam/Robby.Frost/expt_dirs/2023041900_3km_hrrrphys/2023041900/postprd/"
# directory where rap grib data are located
dgrib_r = "/scratch2/BMC/fv3lam/Robby.Frost/expt_dirs/2023041900_3km_rapphys/2023041900/postprd/"
# natlev or prslev
nat_prs = "natlev"
# message number for SBCAPE
mn_refc = 1408
# directory for figure to be output
figdir = "/scratch2/BMC/fv3lam/Robby.Frost/figures/20230419/sbcape/"
# --------------------------------
# plotting setup
rc('font',weight='normal',size=12.5)
# rc('text',usetex='True')
rc('figure',facecolor='white')
# --------------------------------
for hr in range(hri,hrf):
    print(f"Hour {hr}")
    # read in data
    sbcape_h, lat, lon, valid_date = read_grib(hr, dgrib_h, nat_prs, mn_refc)
    sbcape_r = read_grib(hr, dgrib_r, nat_prs, mn_refc, array_only=True)
    # --------------------------------
    # Side by side sbcape plot
    print(f"Creating 1 x 2 SBCAPE Plot")

    # Define your custom colorbar bounds
    cbar_min = 0
    cbar_max = 4501
    clevs = np.arange(cbar_min, cbar_max, 100)
    # color palette
    colors = "CMRmap_r"

    # create plot
    fig, ax = plt.subplots(ncols=2, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(16,6), constrained_layout=True)

    # plot HRRR
    c0 = ax[0].contourf(lon, lat, sbcape_h.values, clevs, 
                        transform=ccrs.PlateCarree(), 
                        cmap=colors)
    # plot RAP
    c1 = ax[1].contourf(lon, lat, sbcape_r.values, clevs, 
                        transform=ccrs.PlateCarree(), 
                        cmap=colors)

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
    cbar.set_label('SBCAPE [J kg$^{-1}$]')

    # save and close figure
    figdir_full = f"{figdir}sbcape_sidebyside_f{hr}.png"
    print(f"Saving figure to {figdir_full}")
    plt.savefig(figdir_full)
    plt.close()
    print("Finished plotting 1 x 2 SBCAPE!")
    # --------------------------------
    # SBCAPE Difference Plot
    print("Creating SBCAPE Difference Plot")

    # Define your custom colorbar bounds
    cbar_min = -1000.1
    cbar_max = 1000.1
    clevs = np.linspace(cbar_min, cbar_max, 51)
    # color palette
    colors = seaborn.color_palette("seismic", as_cmap=True)

    # create plot
    fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()}, 
                        figsize=(10,6.5), constrained_layout=True)

    # plot HRRR - RAP
    c0 = ax.contourf(lon, lat, sbcape_h.values - sbcape_r.values,
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
    cbar.set_label('SBCAPE [J kg$^{-1}$]')
    cbar.set_ticks(np.arange(cbar_min, cbar_max, 250))

    # save and close figure
    figdir_full = f"{figdir}sbcape_diff_f{hr}.png"
    print(f"Saving figure to {figdir_full}")
    plt.savefig(figdir_full)
    plt.close()
    print(f"Finished with hour {hr}! \n")
