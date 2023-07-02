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
import cartopy.crs as ccrs
import cartopy.feature as cpf
import numpy as np
from metpy.plots import ctables
from plotting_functions import read_grib
# --------------------------------
# important parameters

# directory where hrrr grib data are located
dgrib_h = "/scratch2/BMC/fv3lam/Robby.Frost/expt_dirs/2019052000_3km_hrrrphys/2019052000/postprd/"
# directory where rap grib data are located
dgrib_r = "/scratch2/BMC/fv3lam/Robby.Frost/expt_dirs/2019052000_3km_rapphys/2019052000/postprd/"
# natlev or prslev
nat_prs = "prslev"
# message number for composite reflectivity
mn_refc = 40
# directory for figure to be output
figdir = "/scratch2/BMC/fv3lam/Robby.Frost/figures/20190520/refc/"
# --------------------------------
# read in data

# create lists to store data
refch_all = []
refcr_all = []
# start and end hours
hri = 33 # first hour to be plot
hrf = 37 # enter last hour + 1
# loop over forecast hours of interest
for hr in range(hri,hrf):
    refc_h, lat, lon, valid_date = read_grib(hr, dgrib_h, nat_prs, mn_refc)
    refc_r = read_grib(hr, dgrib_r, nat_prs, mn_refc, array_only=True)
    # append hourly refc array to list
    refch_all.append(refc_h)
    refcr_all.append(refc_r)
# --------------------------------
# Side by side reflectivity plot

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
plt_area = [-104, -94, 30, 39] # W, E, S, N
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
plt.savefig(f"{figdir_full}")
plt.close()
print(f"Finished plotting!")