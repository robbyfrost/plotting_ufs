# --------------------------------
# Name: plot_refc_sidebyside.py
# Author: Robert M. Frost
# NOAA Global Systems Laboratory
# Created: 28 June 2023
# Purpose: Loop to plot composite
# reflectivity over Oklahoma and 
# Texas using HRRR and RAP output
# and looping over the forecast
# --------------------------------
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cpf
import numpy as np
from metpy.plots import ctables
import pygrib
# --------------------------------
# define directory location of grib2 output
dhrrr = "/scratch2/BMC/fv3lam/Robby.Frost/expt_dirs/2019052000_3km_hrrrphys/2019052000/postprd/"
drap = "/scratch2/BMC/fv3lam/Robby.Frost/expt_dirs/2019052000_3km_rapphys/2019052000/postprd/"
# define prslev or natlev
nat_prs = "prslev"

# loop over forecast hours of interest
for hr in range(37):
    # set full filepaths
    if hr < 10:
        dgrib_h = f"{dhrrr}rrfs.t00z.{nat_prs}.f00{hr}.rrfs_conuscompact_3km.grib2"
        dgrib_r = f"{drap}rrfs.t00z.{nat_prs}.f00{hr}.rrfs_conuscompact_3km.grib2"
    else:
        dgrib_h = f"{dhrrr}rrfs.t00z.{nat_prs}.f0{hr}.rrfs_conuscompact_3km.grib2"
        dgrib_r = f"{drap}rrfs.t00z.{nat_prs}.f0{hr}.rrfs_conuscompact_3km.grib2"

    # read in hrrr and rap output
    print(f"Reading in F0{hr} HRRR and RAP output.")
    grbs_h = pygrib.open(dgrib_h)
    grbs_r = pygrib.open(dgrib_r)
    # extract composite reflectivity arrays
    refc_h = grbs_h[40]
    refc_r = grbs_r[40]
    # set arrays containing latitude and longitude values
    lat, lon = refc_h.latlons()
    # extract date in datetime format
    valid_time = refc_h.validDate

    # --------------------------------
    # side by side reflectivity plot

    # Define your custom colorbar bounds
    cbar_min = 0
    cbar_max = 75.1
    # set reflectivity levels to be plotted
    clevs = np.arange(cbar_min, cbar_max, 5)
    # define color table using metpy colortables
    colors = ctables.registry.get_colortable('NWSReflectivity')

    # create plot
    print(f"Creating Reflectivity plot for F0{hr}")
    fig, ax = plt.subplots(ncols=2, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(12,6.7), constrained_layout=True)
    # mapping
    plt_area = [-104, -94, 30, 39] # W, E, S, N
    for i, iax in enumerate(ax):
        iax.coastlines()
        iax.add_feature(cpf.BORDERS)
        iax.add_feature(cpf.STATES)
        iax.set_extent(plt_area)

    # plot HRRR
    c0 = ax[0].contourf(lon, lat, refc_h.values, clevs, 
                        transform=ccrs.PlateCarree(), 
                        cmap=colors)
    # plot RAP
    c1 = ax[1].contourf(lon, lat, refc_r.values, clevs,
                        transform=ccrs.PlateCarree(), 
                        cmap=colors)

    # add axes titles
    ax[0].set_title(f"HRRR F0{hr},  Valid {valid_time} UTC")
    ax[1].set_title(f"RAP F0{hr},  Valid {valid_time} UTC")

    # Add colorbar
    cbar = fig.colorbar(c0, ax=ax, orientation='horizontal', extend=True, pad=0.05, aspect=50)
    cbar.set_label('Simulated Composite Reflectivty [dBZ]')

    dout = f"/scratch2/BMC/fv3lam/Robby.Frost/figures/20190520/refc/refc_sidebyside_f{hr}.png"
    print(f"Saving figure here: {dout}")
    plt.savefig(dout)
    plt.close()
    print(f"Finished plotting F0{hr}!")