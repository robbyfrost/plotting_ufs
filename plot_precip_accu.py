# --------------------------------
# Name: plot_refc_sidebyside.py
# Author: Robert M. Frost
# NOAA Global Systems Laboratory
# Created: 29 June 2023
# Purpose: Plot forecast precipitation
# accumilation totals and un-resolved
# precipitation to see the impacts of
# the GF convective scheme
# --------------------------------
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cpf
from metpy.plots import ctables
import pygrib
import matplotlib.colors as colors
from matplotlib import rc
# --------------------------------
# hour to be plot
hr = 36
# define directory location of grib output
dhrrr = "/scratch2/BMC/fv3lam/Robby.Frost/expt_dirs/2019052000_3km_hrrrphys/2019052000/postprd/"
drap = "/scratch2/BMC/fv3lam/Robby.Frost/expt_dirs/2019052000_3km_rapphys/2019052000/postprd/"
# natlev of prslev
nat_prs = "natlev"
# set filepaths
if hr < 10:
    dgrib_h = f"{dhrrr}rrfs.t00z.{nat_prs}.f00{hr}.rrfs_conuscompact_3km.grib2"
    dgrib_r = f"{drap}rrfs.t00z.{nat_prs}.f00{hr}.rrfs_conuscompact_3km.grib2"
else:
    dgrib_h = f"{dhrrr}rrfs.t00z.{nat_prs}.f0{hr}.rrfs_conuscompact_3km.grib2"
    dgrib_r = f"{drap}rrfs.t00z.{nat_prs}.f0{hr}.rrfs_conuscompact_3km.grib2"

# open hrrr and rap output
print("Reading in grib output")
grbs_h = pygrib.open(dgrib_h)
grbs_r = pygrib.open(dgrib_r)
# extract total precip arrays (f000 - f036)
tot_accu_h = grbs_h[1374]
tot_accu_r = grbs_r[1374]
# extract non-convective precip arrays (f000 - f036)
nonc_accu_h = grbs_h[1376]
nonc_accu_r = grbs_r[1376]

# extract latitude and longitude arrays
lat, lon = tot_accu_h.latlons()
# extract datatime TODO: Figure out why this doesn't work LOL
valid_time = tot_accu_h.validDate
# --------------------------------
# precip accumilation 2x2 plot

# set font size
rc("font",weight="normal",size=15)

# Define your custom colorbar bounds
# cbar_min = 0
# cbar_max = 400
# levels to be plot
clevs = [0, 5, 10, 15, 25, 40, 60, 80, 100, 125, 150, 200, 250, 300, 400]
# precipitation colortable
colors = ctables.registry.get_colortable('precipitation')

# create plot
fig, ax = plt.subplots(ncols=2, nrows=2, subplot_kw={'projection': ccrs.PlateCarree()}, 
                       figsize=(18,16), constrained_layout=True)

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
    
# plot
print(f"Plotting F0{hr}")
c0 = ax[0,0].contourf(lon, lat, tot_accu_h.values, 
                    clevs, transform=ccrs.PlateCarree(), 
                    cmap=colors, extend='max')
c1 = ax[0,1].contourf(lon, lat, tot_accu_r.values, 
                    clevs, transform=ccrs.PlateCarree(), 
                    cmap=colors, extend='max')
c2 = ax[1,0].contourf(lon, lat, tot_accu_h.values - nonc_accu_h.values, 
                    clevs, transform=ccrs.PlateCarree(), 
                    cmap=colors, extend='max')
c3 = ax[1,1].contourf(lon, lat, tot_accu_r.values - nonc_accu_r.values, 
                    clevs, transform=ccrs.PlateCarree(), 
                    cmap=colors, extend='max')

# pretty up
ax[0,0].set_title(f"HRRR F0{hr},  Valid 2019-05-21 12:00:00 UTC")
ax[0,1].set_title(f"RAP F0{hr},  Valid 2019-05-21 12:00:00 UTC")
ax[1,0].set_title(f"HRRR F0{hr}, Valid 2019-05-21 12:00:00 UTC")
ax[1,1].set_title(f"RAP F0{hr}, Valid 2019-05-21 12:00:00 UTC")

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
dout = f"/scratch2/BMC/fv3lam/Robby.Frost/figures/20190520/precip/precip_accum_f{hr}.png"
plt.savefig(dout)
plt.close()
print(f"Figure saved to {dout}")