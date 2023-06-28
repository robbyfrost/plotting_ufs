# --------------------------------
# Name: plot_sbcape_loop.py
# Author: Robert M. Frost
# NOAA Global Systems Laboratory
# Created: 26 June 2023
# Purpose: Loop to plot 2m dew point
# and wind barb comparisons at 
# different times during forecast runs
# --------------------------------
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cpf
from datetime import datetime
from matplotlib import rc
import seaborn
import numpy as np
# --------------------------------
for hr in range(0,37):
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
    # read in dew point data
    hrrr = xr.open_dataset(dgrib_h, engine="cfgrib", 
                        filter_by_keys={'stepType': 'instant', 'typeOfLevel': 'heightAboveGround'})
    rap = xr.open_dataset(dgrib_r, engine="cfgrib",
                        filter_by_keys={'stepType': 'instant', 'typeOfLevel': 'heightAboveGround'})
    # convert dew point to fahrenheit
    hrrr_td2m = ((hrrr.d2m - 273.15) * 9 / 5) + 32
    rap_td2m = ((rap.d2m - 273.15) * 9 / 5) + 32
    # read in wind barb data
    h = xr.open_dataset(dgrib_h, engine="cfgrib", 
                       filter_by_keys={'typeOfLevel': 'hybrid'})
    r = xr.open_dataset(dgrib_r, engine="cfgrib", 
                        filter_by_keys={'typeOfLevel': 'hybrid'})
    # arrays of u and v
    hu = h.u
    hv = h.v
    ru = r.u
    rv = r.v

    # latitude and longitude arrays
    lat = hrrr.latitude
    lon = hrrr.longitude

    # convert numpy.datetime to datetime value
    valid_time_str = str(hrrr.time.values + hrrr.time.step.values)
    valid_time_str = valid_time_str[:19]
    valid_time = datetime.strptime(valid_time_str, '%Y-%m-%dT%H:%M:%S')

    # plotting setup
    rc('font',weight='normal',size=15)#,family='serif',serif='Times New Roman')
    # rc('text',usetex='True')
    rc('figure',facecolor='white')
    # cmap = seaborn.color_palette("ColorBrewer", as_cmap=True)
    colors = seaborn.color_palette("YlGnBu", as_cmap=True)

    # --------------------------------
    # side by side td2m plot

    # dew point colorbar
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

    # Define your custom colorbar bounds
    cbar_min = 0
    cbar_max = 80

    clevs = np.arange(cbar_min, cbar_max, 2)

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
    c0 = ax[0].contourf(lon, lat, hrrr_td2m, clevs, 
                        transform=ccrs.PlateCarree(), 
                        cmap=colors,
                        extend="both")
    c1 = ax[1].contourf(lon, lat, rap_td2m, clevs,
                        transform=ccrs.PlateCarree(), 
                        cmap=colors,
                        extend="both")

    # pretty up
    ax[0].set_title(f"HRRR F0{hr},  Valid {valid_time} UTC")
    ax[1].set_title(f"RAP F0{hr},  Valid {valid_time} UTC")

    # Add colorbar
    cbar = fig.colorbar(c1, ax=ax, orientation='horizontal', extend=True, pad=0.05, aspect=50)
    cbar.set_label('2m Dew Point Temperature [$^{\circ}$F]')

    # Wind barbs
    spacing=25 #barbspacing (smaller if zoomed in)
    ax[0].barbs(lon[::spacing,::spacing], lat[::spacing,::spacing],
                hu[0,::spacing,::spacing], hv[0,::spacing,::spacing], 
                length=6)
    ax[1].barbs(lon[::spacing,::spacing], lat[::spacing,::spacing],
                ru[0,::spacing,::spacing], rv[0,::spacing,::spacing], 
                length=6)

    plt.savefig(f"/scratch2/BMC/fv3lam/Robby.Frost/figures/20190520/td2m_sidebyside_f{hr}.png")
    plt.close()

    # --------------------------------
    # comparison td2m plot

    # Define your custom colorbar bounds
    cbar_min = -30
    cbar_max = 30

    clevs = np.arange(cbar_min, cbar_max, 2)

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
    c0 = ax.contourf(lon, lat, hrrr_td2m - rap_td2m, clevs, 
                        transform=ccrs.PlateCarree(), 
                        cmap=colors,
                        extend="both")

    # pretty up
    ax.set_title(f"F0{hr},  Valid {valid_time} UTC")

    # Add colorbar
    cbar = fig.colorbar(c0, ax=ax, orientation='horizontal', extend=True, pad=0.05, aspect=50)
    cbar.set_label('HRRR - RAP 2m Dew Point Temperature [$^{\circ}$F]')

    plt.savefig(f"/scratch2/BMC/fv3lam/Robby.Frost/figures/20190520/td2m_diff_f{hr}.png")
    plt.close()