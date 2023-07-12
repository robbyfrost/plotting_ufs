# --------------------------------
# Name: plot_refc_loop.py
# Author: Robert M. Frost
# NOAA Global Systems Laboratory
# Created: 02 July 2023
# Purpose: Loop to plot composite
# reflectivity over Oklahoma and 
# Texas using HRRR and RAP output
# and looping over the forecast
# --------------------------------
import matplotlib.pyplot as plt
from matplotlib import rc
import cartopy.crs as ccrs
import cartopy.feature as cpf
import numpy as np
from metpy.plots import ctables
from UFSutils import read_grib
import geopandas as gpd
# --------------------------------
# important parameters

# directory where hrrr grib data are located
dgrib_h = "/scratch2/BMC/fv3lam/Robby.Frost/expt_dirs/2023041900_3km_hrrrphys/2023041900/postprd/"
# directory where rap grib data are located
dgrib_r = "/scratch2/BMC/fv3lam/Robby.Frost/expt_dirs/2023041900_3km_rapphys/2023041900/postprd/"
# natlev or prslev
nat_prs = "prslev"
# message number for composite reflectivity
mn_refc = 40
# message number for updraft helicity
mn_heli = 936
# directory for figure to be output
figdir = "/scratch2/BMC/fv3lam/Robby.Frost/figures/20230419/refc/"
# flag to plot 4 hour refc
four_hr = False
# flag to plot single hour refc
single_hr = True
# --------------------------------
# plotting setup
rc('font',weight='normal',size=15)
# rc('text',usetex='True')
rc('figure',facecolor='white')
# --------------------------------
# Side by side 4 hour reflectivity plot
if four_hr:
    # read in data

    # create lists to store data
    refch_all = []
    refcr_all = []
    # start and end hours
    hri = 1 # first hour to be plot
    hrf = 5 # enter last hour + 1
    # loop over forecast hours of interest
    for hr in range(hri,hrf):
        hrrr, refc_h, lat, lon, valid_date = read_grib(hr, dgrib_h, nat_prs, mn_refc)
        refc_r = read_grib(hr, dgrib_r, nat_prs, mn_refc, ret_type=1)
        # append hourly refc array to list
        refch_all.append(refc_h)
        refcr_all.append(refc_r)
    # --------------------------------
    # Plot

    # Define your custom colorbar bounds
    cbar_min = 0
    cbar_max = 75.1
    # set reflectivity levels to be plotted
    clevs = np.arange(cbar_min, cbar_max, 5)
    # define color table using metpy colortables
    colors = ctables.registry.get_colortable('NWSReflectivity')

    # create plot
    print(f"Creating 4 x 2 Reflectivity plot")
    fig, ax = plt.subplots(ncols=2, nrows=4, subplot_kw={'projection': ccrs.PlateCarree()}, 
                        figsize=(9,16.7), constrained_layout=True)
    # mapping
    geoData = gpd.read_file('https://raw.githubusercontent.com/holtzy/The-Python-Graph-Gallery/master/static/data/US-counties.geojson')
    plt_area = [-101, -94, 33.5, 37.5] # W, E, S, N
    for i, iax in enumerate(ax[:,0]):
        iax.coastlines()
        iax.add_feature(cpf.BORDERS)
        iax.add_feature(cpf.STATES)
        iax.set_extent(plt_area)
        geoData.plot(ax=iax, color="none", lw=0.3, aspect=1)
    for i, iax in enumerate(ax[:,1]):
        iax.coastlines()
        iax.add_feature(cpf.BORDERS)
        iax.add_feature(cpf.STATES)
        iax.set_extent(plt_area)
        # Load the json file with county coordinates
        geoData.plot(ax=iax, color="none", lw=0.3, aspect=1)

    # plot HRRR
    c0 = ax[0,0].contourf(lon, lat, refch_all[0].values, clevs, 
                        transform=ccrs.PlateCarree(), cmap=colors)
    c2 = ax[1,0].contourf(lon, lat, refch_all[1].values, clevs, 
                        transform=ccrs.PlateCarree(), cmap=colors)
    c4 = ax[2,0].contourf(lon, lat, refch_all[2].values, clevs, 
                        transform=ccrs.PlateCarree(), cmap=colors)
    c6 = ax[3,0].contourf(lon, lat, refch_all[3].values, clevs, 
                        transform=ccrs.PlateCarree(), cmap=colors)
    # plot RAP
    c1 = ax[0,1].contourf(lon, lat, refcr_all[0].values, clevs,
                        transform=ccrs.PlateCarree(), cmap=colors)
    c3 = ax[1,1].contourf(lon, lat, refcr_all[1].values, clevs,
                        transform=ccrs.PlateCarree(), cmap=colors)
    c5 = ax[2,1].contourf(lon, lat, refcr_all[2].values, clevs,
                        transform=ccrs.PlateCarree(), cmap=colors)
    c7 = ax[3,1].contourf(lon, lat, refcr_all[3].values, clevs,
                        transform=ccrs.PlateCarree(), cmap=colors)

    # add axes titles
    for i, iax in enumerate(ax[:,0]):
        iax.set_title(f"HRRR Valid at {refch_all[i].validDate} UTC")
    for i, iax in enumerate(ax[:,1]):
        iax.set_title(f"RAP Valid at {refcr_all[i].validDate} UTC")

    # Add colorbar
    cbar = fig.colorbar(c0, ax=ax, orientation='horizontal', extend=True, pad=0.01, fraction=0.017, aspect=30)
    cbar.set_label('Simulated Composite Reflectivty [dBZ]')

    # tight_layout
    plt.tight_layout

    # save figure
    figdir_full = f"{figdir}refc_4x2_f{hri}-f{hrf-1}.png"
    print(f"Saving figure to {figdir_full}")
    plt.savefig(figdir_full)
    plt.close()
    print("Finished plotting! \n")
# --------------------------------
# Side by side 1 hour reflectivity

if single_hr:
    # loop over forecast hours of interest
    for hr in range(37):
        print(f"Creating Reflectivity Plot for Hour {hr}")
        # read in reflectivity
        grbs_h, refc_h, lat, lon, valid_date = read_grib(hr, dgrib_h, nat_prs, mn_refc)
        grbs_r, refc_r, lat, lon, valid_date = read_grib(hr, dgrib_r, nat_prs, mn_refc)
        # read in updraft helicity
        heli_h = grbs_h[mn_heli]
        heli_r = grbs_r[mn_heli]

        # Define your custom colorbar bounds
        cbar_min = 0
        cbar_max = 75.1
        # set reflectivity levels to be plotted
        clevs = np.arange(cbar_min, cbar_max, 5)
        # define color table using metpy colortables
        colors = ctables.registry.get_colortable('NWSReflectivity')

        # create plot
        fig, ax = plt.subplots(ncols=2, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(16,10.2), constrained_layout=True)

        # plot HRRR
        c0 = ax[0].contourf(lon, lat, refc_h.values, clevs, 
                            transform=ccrs.PlateCarree(), 
                            cmap=colors)
        ax[0].contourf(lon, lat, heli_h.values, [25, 500],
                      transform=ccrs.PlateCarree(), colors="black", alpha=0.5)
        ax[0].contour(lon, lat, heli_h.values, [25, 500],
                      transform=ccrs.PlateCarree(), colors="black")
        # plot RAP
        c1 = ax[1].contourf(lon, lat, refc_r.values, clevs,
                            transform=ccrs.PlateCarree(), 
                            cmap=colors)
        ax[1].contourf(lon, lat, heli_r.values, [25, 500],
                      transform=ccrs.PlateCarree(), colors="black", alpha=0.5)
        ax[1].contour(lon, lat, heli_r.values, [25, 500],
                      transform=ccrs.PlateCarree(), colors="black")

        # mapping
        plt_area = [-101, -94, 30, 37.5] # W, E, S, N
        for i, iax in enumerate(ax):
            iax.coastlines()
            iax.add_feature(cpf.BORDERS)
            iax.add_feature(cpf.STATES)
            iax.set_extent(plt_area)
            # Load the json file with county coordinates
            geoData = gpd.read_file('https://raw.githubusercontent.com/holtzy/The-Python-Graph-Gallery/master/static/data/US-counties.geojson')
            geoData.plot(ax=iax, color="none", lw=0.3, aspect=1)

        # add axes titles
        plt.suptitle("Composite Reflectivity [dBZ], UH > 25 [m$^{2}$ s$^{-2}$]")
        ax[0].set_title(f"HRRR F0{hr},  Valid {valid_date} UTC")
        ax[1].set_title(f"RAP F0{hr},  Valid {valid_date} UTC")

        # Add colorbar
        cbar = fig.colorbar(c0, ax=ax, orientation='horizontal', 
                            extend=True, pad=0.03, aspect=50)

        # save figure
        figdir_full = f"{figdir}refc_sidebyside_f{hr}.png"
        print(f"Saving figure to {figdir_full}")
        plt.savefig(figdir_full)
        plt.close()
        print("Finished plotting! \n")