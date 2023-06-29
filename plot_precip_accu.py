# This program is used for developing new plotting scripts. Once scripts work properly and can be iterated over, this notebook is imported to a script, where it will be run for most cases. Therefore, code here is not normally used and up to date.

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cpf
import numpy as np
from metpy.plots import ctables
import pygrib
import matplotlib.colors as colors

hr = 36
dhrrr = "/scratch2/BMC/fv3lam/Robby.Frost/expt_dirs/2019052000_3km_hrrrphys/2019052000/postprd/"
drap = "/scratch2/BMC/fv3lam/Robby.Frost/expt_dirs/2019052000_3km_rapphys/2019052000/postprd/"
nat_prs = "natlev"
# set filepaths
if hr < 10:
    dgrib_h = f"{dhrrr}rrfs.t00z.{nat_prs}.f00{hr}.rrfs_conuscompact_3km.grib2"
    dgrib_r = f"{drap}rrfs.t00z.{nat_prs}.f00{hr}.rrfs_conuscompact_3km.grib2"
else:
    dgrib_h = f"{dhrrr}rrfs.t00z.{nat_prs}.f0{hr}.rrfs_conuscompact_3km.grib2"
    dgrib_r = f"{drap}rrfs.t00z.{nat_prs}.f0{hr}.rrfs_conuscompact_3km.grib2"

grbs_h = pygrib.open(dgrib_h)
grbs_r = pygrib.open(dgrib_r)

tot_accu_h = grbs_h[1374]
tot_accu_r = grbs_r[1374]
nonc_accu_h = grbs_h[1376]
nonc_accu_r = grbs_r[1376]

lat, lon = tot_accu_h.latlons()

valid_time = tot_accu_h.validDate

# %%
# Define your custom colorbar bounds
cbar_min = 10
cbar_max = 400

clevs = np.linspace(cbar_min, cbar_max, 50)

clevs_diff = np.linspace(10, 20, 50)

colors = ctables.registry.get_colortable('NWSReflectivity')

# create plot
fig, ax = plt.subplots(ncols=2, nrows=2, subplot_kw={'projection': ccrs.PlateCarree()}, 
                       figsize=(18,16))#, constrained_layout=True)
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
c0 = ax[0,0].contourf(lon, lat, tot_accu_h.values, 
                    clevs, transform=ccrs.PlateCarree(), 
                    cmap=colors, extend='max')
c1 = ax[0,1].contourf(lon, lat, tot_accu_r.values, 
                    clevs, transform=ccrs.PlateCarree(), 
                    cmap=colors, extend='max')
c2 = ax[1,0].contourf(lon, lat, tot_accu_h.values - nonc_accu_h.values, 
                    clevs_diff, transform=ccrs.PlateCarree(), 
                    cmap=colors, extend='max')
c3 = ax[1,1].contourf(lon, lat, tot_accu_r.values - nonc_accu_r.values, 
                    clevs_diff, transform=ccrs.PlateCarree(), 
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
plt.savefig(f"/scratch2/BMC/fv3lam/Robby.Frost/figures/20190520/precip_accum_f{hr}.png")
plt.close()

# %%



