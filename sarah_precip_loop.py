# Name: sarah_precip_loop.py
# 
# Author: Robert M. Frost
# 
# NOAA Global Systems Laboratory
# 
# Created: 28 July 2023
# 
# Purpose: Plot precipitation accumulation for each hour of model run
# --------------------------------
import matplotlib.pyplot as plt
from matplotlib import rc
import cartopy.crs as ccrs
import cartopy.feature as cpf
import numpy as np
from UFSutils import read_grib
import geopandas as gpd
import seaborn
import matplotlib.colors as mcolors
import seaborn
# --------------------------------
# Important parameters

# hour at which the forecast was initialized (UTC)
init = 12
# directory where hrrr grib data are located
dgrib_h = "/scratch2/BMC/fv3lam/Robby.Frost/expt_dirs/2023041912_3km_hrrrphys/2023041912/postprd/"
# directory where rap grib data are located
dgrib_r = "/scratch2/BMC/fv3lam/Robby.Frost/expt_dirs/2023041912_3km_rapphys/2023041912/postprd/"
# natlev or prslev
nat_prs = "natlev"
# message number for total precip accu
mn_precip = 1374
# message number for non-convective precip accu
mn_nonc = 1376
# directory for figure to be output
figdir = "/scratch2/BMC/fv3lam/Robby.Frost/figures/sarah_test/"
# --------------------------------
# plotting setup
rc('font',weight='normal',size=12.5)
# rc('text',usetex='True')
rc('figure',facecolor='white')
# custom colorbar using NWS precip accumulation
nws_precip_colors = [
    "#fdfdfd",
    "#a9f5f4",
    "#33aff2",
    "#0300f4",
    "#02fd02",
    "#01c501",
    "#008e00",
    "#fdf802",
    "#e5bc00",
    "#fd9500",
    "#fd0000",
    "#d40000",
    "#bc0000",
    "#f800fd",
    "#9854c6",
    "#fdfdfd" 
]
# --------------------------------
# Read in precip

for hr in range(0,2):
    hrrr, precip_h, lat, lon, valid_date = read_grib(init, hr, dgrib_h, nat_prs, mn_precip, ret_type=0)
    rap, precip_r, lat, lon, valid_date = read_grib(init, hr, dgrib_r, nat_prs, mn_precip, ret_type=0)
    nonc_h = hrrr[mn_nonc]
    nonc_r = rap[mn_nonc]
    # --------------------------------
    # Plot precip accumulation in Oklahoma and Texas

    colors = mcolors.ListedColormap(nws_precip_colors)
    color_diff = seaborn.color_palette("seismic", as_cmap=True)
    # custom colorbar values
    clevs = [-0.1, 0.0, 2.5, 5, 10, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100]
   # Define your custom colorbar bounds
    cbar_min = -50
    cbar_max = 50.1
    clevs_diff = np.linspace(cbar_min, cbar_max, 50)
    # normalize colors around clevs
    norm = mcolors.BoundaryNorm(clevs, 16)

    # create plot
    fig, ax = plt.subplots(nrows=3, subplot_kw={'projection': ccrs.PlateCarree()}, 
                        figsize=(10,15), constrained_layout=True)

    # plot HRRR
    c0 = ax[0].contourf(lon, lat, precip_h.values, clevs, 
                        transform=ccrs.PlateCarree(), 
                        cmap=colors, norm=norm, extend="max")
    # plot RAP
    c1 = ax[1].contourf(lon, lat, precip_r.values, clevs, 
                        transform=ccrs.PlateCarree(), 
                        cmap=colors, norm=norm, extend="max")
    # plot difference
    c2 = ax[2].contourf(lon, lat, precip_r.values - precip_h.values,
                        clevs_diff, transform = ccrs.PlateCarree(), cmap=color_diff)

    # mapping
    for i, iax in enumerate(ax):
        iax.coastlines()
        iax.add_feature(cpf.BORDERS)
        iax.add_feature(cpf.STATES)

    # set title
    ax[0].set_title(f"HRRR F0{hr},  Valid 2023-04-20 12:00:00 UTC")
    ax[1].set_title(f"RAP F0{hr},  Valid 2023-04-20 12:00:00 UTC")

    # Add colorbar
    cbar = fig.colorbar(c1, ax=ax[1], orientation='horizontal', extend=True, pad=0.03, aspect=50)
    cbar.set_label('Total Forecast Precipitation Accumulation [mm]')
    cbar.set_ticks(clevs)

    cbar_diff = fig.colorbar(c2,ax=ax[2], orientation='horizontal', extend=True, pad=0.03, aspect=50)
    cbar_diff.set_label('RAP - HRRR Total Precipitation Accumulation [mm]')
    cbar_diff.set_ticks(np.arange(cbar_min, cbar_max, 10))

    # save and close figure
    plt.savefig(f"{figdir}sarah_test_f{hr}.png")
    plt.close()