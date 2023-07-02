# --------------------------------
# Name: plot_sbcape_loop.py
# Author: Robert M. Frost
# NOAA Global Systems Laboratory
# Created: 26 June 2023
# Purpose: Loop to plot 2m dew point
# and wind barb comparisons at 
# different times during forecast runs
# --------------------------------
from plotting_functions import read_grib
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cpf
import seaborn
import numpy as np
# --------------------------------
# settings

# directory where hrrr grib data are located
dgrib_h = "/scratch2/BMC/fv3lam/Robby.Frost/expt_dirs/2019052000_3km_hrrrphys/2019052000/postprd/"
# directory where rap grib data are located
dgrib_r = "/scratch2/BMC/fv3lam/Robby.Frost/expt_dirs/2019052000_3km_rapphys/2019052000/postprd/"
# natlev or prslev
nat_prs = "natlev"
# message number for dew point
mn_td2m = 1358
# message number for u at 10m
mn_u10 = 1364
# message number for v at 10m
mn_v10 = 1365
# directory for figure to be output
figdir = "/scratch2/BMC/fv3lam/Robby.Frost/figures/20190520/td2m/"

for hr in range(0,37):
    # read in dew point
    td2m_h, lat, lon, valid_date = read_grib(hr, dgrib_h, nat_prs, mn_td2m)
    td2m_r = read_grib(hr, dgrib_r, nat_prs, mn_td2m, array_only=True)
    # convert to fahrenheit (superior unit of temperature)
    td2m_h = (td2m_h.values - 273.15) * (9/5) + 32
    td2m_r = (td2m_r.values - 273.15) * (9/5) + 32

    # read in 10m wind
    u10_h = read_grib(hr, dgrib_h, nat_prs, mn_u10, array_only=True).values
    v10_h = read_grib(hr, dgrib_h, nat_prs, mn_v10, array_only=True).values
    u10_r = read_grib(hr, dgrib_r, nat_prs, mn_u10, array_only=True).values
    v10_r = read_grib(hr, dgrib_r, nat_prs, mn_v10, array_only=True).values
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
    # Plot dew point comparison

    # Define your custom colorbar bounds
    cbar_min = 0
    cbar_max = 80
    # contour levels
    clevs = np.arange(cbar_min, cbar_max, 2)
    # set NWS dew point colormap
    colors = cmap

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
    c0 = ax[0].contourf(lon, lat, td2m_h, clevs, 
                        transform=ccrs.PlateCarree(), 
                        cmap=colors,
                        extend="both")
    c1 = ax[1].contourf(lon, lat, td2m_r, clevs,
                        transform=ccrs.PlateCarree(), 
                        cmap=colors,
                        extend="both")

    # pretty up
    ax[0].set_title(f"HRRR F0{hr},  Valid {valid_date} UTC")
    ax[1].set_title(f"RAP F0{hr},  Valid {valid_date} UTC")

    # Add colorbar
    cbar = fig.colorbar(c1, ax=ax, orientation='horizontal', extend=True, pad=0.05, aspect=50)
    cbar.set_label('2m Dew Point Temperature [$^{\circ}$F]')

    # Wind barbs
    spacing=25 #barbspacing (smaller if zoomed in)
    ax[0].barbs(lon[::spacing,::spacing], lat[::spacing,::spacing],
                u10_h[::spacing,::spacing], v10_h[::spacing,::spacing], 
                length=6)
    ax[1].barbs(lon[::spacing,::spacing], lat[::spacing,::spacing],
                u10_r[::spacing,::spacing], v10_r[::spacing,::spacing], 
                length=6)

    plt.savefig(f"{figdir}td2m_sidebyside_f{hr}.png")
    plt.close() 
    # --------------------------------
    # Plot dew point comparison

    # Define your custom colorbar bounds
    cbar_min = -30
    cbar_max = 30
    # contour levels
    clevs = np.arange(cbar_min, cbar_max, 2)
    # color palette
    colors = seaborn.color_palette("seismic", as_cmap=True)

    # create plot
    fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(10,10))

    # mapping
    plt_area = [-104, -94, 30, 39] # W, E, S, N
    ax.coastlines()
    ax.add_feature(cpf.BORDERS)
    ax.add_feature(cpf.STATES)
    ax.set_extent(plt_area)

    # plot
    c0 = ax.contourf(lon, lat, td2m_h - td2m_r, clevs, 
                        transform=ccrs.PlateCarree(), 
                        cmap=colors,
                        extend="both")

    # pretty up
    ax.set_title(f"F0{hr},  Valid {valid_date} UTC")

    # Add colorbar
    cbar = fig.colorbar(c0, ax=ax, orientation='horizontal', extend=True, pad=0.05, aspect=50)
    cbar.set_label('HRRR - RAP 2m Dew Point Temperature [$^{\circ}$F]')

    # save figure
    plt.savefig(f"{figdir}td2m_diff_f{hr}.png")
    plt.close()