# --------------------------------
# Name: plot_sbcape_loop.py
# Author: Robert M. Frost
# NOAA Global Systems Laboratory
# Created: 26 June 2023
# Purpose: Loop to plot 2m dew point
# and wind barb comparisons at 
# different times during forecast runs
# --------------------------------
from UFSutils import read_grib
import matplotlib.pyplot as plt
from matplotlib import rc
import cartopy.crs as ccrs
import cartopy.feature as cpf
import seaborn
import numpy as np
import geopandas as gpd
# --------------------------------
# settings

# date being plot
date = "2023041912" #YYYMMDDHH
# hour forecast was initialized (UTC)
init = 12
# directory where hrrr grib data are located
dgrib_h = f"/scratch2/BMC/fv3lam/Robby.Frost/expt_dirs/{date}_3km_hrrrphys/{date}/postprd/"
# directory where rap grib data are located
dgrib_r = f"/scratch2/BMC/fv3lam/Robby.Frost/expt_dirs/{date}_3km_rapphys/{date}/postprd/"
# natlev or prslev
nat_prs = "natlev"
# message number for dew point
mn_td2m = 1358
# message number for u at 10m
mn_u10 = 1364
# message number for v at 10m
mn_v10 = 1365
# directory for figure to be output
figdir = f"/scratch2/BMC/fv3lam/Robby.Frost/figures/{date}/td2m/"
# --------------------------------
# plotting setup
rc('font',weight='normal',size=12.5)
# rc('text',usetex='True')
rc('figure',facecolor='white')
# --------------------------------
# NWS dew point colorbar

import matplotlib.colors as colors
a = np.array([0,10,20,30,40,45,50,55,60,65,70,75,80])
# Normalize the bin between 0 and 1 (uneven bins are important here)
norm = [(float(i)-min(a))/(max(a)-min(a)) for i in a]
# Color tuple for every bin
C = np.array([[59,34,4],
            [84,48,5],
            [140,82,10],
            [191,129,45],
            [204,168,84],
            [223,194,125],
            [230,217,181],
            [211,235,231],
            [169,219,211],
            [114,184,173],
            [49,140,133],
            [1,102,95],
            [0,60,48],
            [0,41,33]])
# Create a tuple for every color indicating the normalized position on the colormap and the assigned color.
COLORS = []
for i, n in enumerate(norm):
    COLORS.append((n, np.array(C[i])/255.))
# Create the colormap
cmap = colors.LinearSegmentedColormap.from_list("dewpoint", COLORS)
# --------------------------------
# loop over time

for hr in range(0,37):
    print(f"Hour {hr}")

    # read in dew point
    hrrr, td2m_h, lat, lon, valid_date = read_grib(init, hr, dgrib_h, nat_prs, mn_td2m, ret_type=0)
    rap, td2m_r, lat, lon, valid_date = read_grib(init, hr, dgrib_r, nat_prs, mn_td2m, ret_type=0)
    # convert to fahrenheit (superior unit of temperature)
    td2m_h = (td2m_h.values - 273.15) * (9/5) + 32
    td2m_r = (td2m_r.values - 273.15) * (9/5) + 32

    # read in 10m wind
    u10_h = hrrr[mn_u10].values
    v10_h = hrrr[mn_v10].values
    u10_r = rap[mn_u10].values
    v10_r = rap[mn_v10].values

    # convert 10m wind to knots
    u10_h = u10_h * 1.944
    v10_h = v10_h * 1.944
    u10_r = u10_r * 1.944
    v10_r = v10_r * 1.944
    # --------------------------------
    # Plot dew point comparison
    print(f"Creating 1 x 2 Td2m Plot")

    # Define your custom colorbar bounds
    cbar_min = 0
    cbar_max = 80.1
    # levels for sbcape to be plot
    clevs = np.arange(cbar_min, cbar_max, 2)

    # create plot
    fig, ax = plt.subplots(ncols=2, subplot_kw={'projection': ccrs.PlateCarree()}, 
                        figsize=(16,10), constrained_layout=True)

    # plot HRRR
    c0 = ax[0].contourf(lon, lat, td2m_h, clevs, 
                        transform=ccrs.PlateCarree(), 
                        cmap=cmap, extend="both")
    # plot RAP
    c1 = ax[1].contourf(lon, lat, td2m_r, clevs, 
                        transform=ccrs.PlateCarree(), 
                        cmap=cmap, extend="both")

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

    # set title
    ax[0].set_title(f"No-GF F0{hr},  Valid {valid_date} UTC")
    ax[1].set_title(f"GF F0{hr},  Valid {valid_date} UTC")

    # Add colorbar
    cbar = fig.colorbar(c1, ax=ax, orientation='horizontal', extend=True, pad=0.03, aspect=50)
    cbar.set_label('2m Dew Point Temperature [$^{\circ}$F]')
    cbar.set_ticks(np.arange(cbar_min, cbar_max, 10))

    # Wind barbs
    spacing=25 #barbspacing (smaller if zoomed in)
    ax[0].barbs(lon[::spacing,::spacing], lat[::spacing,::spacing],
                u10_h[::spacing,::spacing], v10_h[::spacing,::spacing], 
                length=6)
    ax[1].barbs(lon[::spacing,::spacing], lat[::spacing,::spacing],
                u10_r[::spacing,::spacing], v10_r[::spacing,::spacing], 
                length=6)
    
    # save and close figure
    figdir_full = f"{figdir}td2m_sidebyside_f{hr}.png"
    print(f"Saving figure to {figdir_full}")
    plt.savefig(figdir_full)
    plt.close() 
    print("Finished plotting 1 x 2 Td2m!")
    # --------------------------------
    # Plot dew point comparison
    print("Creating Td2m Difference Plot!")

    # Define your custom colorbar bounds
    cbar_min = -30
    cbar_max = 30.1
    # contour levels
    clevs = np.linspace(cbar_min, cbar_max, 50)
    # color palette
    colors = seaborn.color_palette("seismic", as_cmap=True)

    # create plot
    fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()}, 
                        figsize=(10,6.5), constrained_layout=True)

    # plot HRRR - RAP
    c0 = ax.contourf(lon, lat, td2m_h - td2m_r,
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
    cbar.set_label('HRRR - RAP 2m Dew Point Temperature [$^{\circ}$F]')
    cbar.set_ticks(np.arange(cbar_min, cbar_max, 5))

    # save and close figure
    figdir_full = f"{figdir}td2m_diff_f{hr}.png"
    print(f"Saving figure to {figdir_full}")
    plt.savefig(figdir_full)
    plt.close()
    print(f"Finished with hour {hr}! \n")