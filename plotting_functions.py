# --------------------------------
# Name: plotting_functions.py
# Author: Robert M. Frost
# NOAA Global Systems Laboratory
# Created: 16 June 2023
# Purpose: A collection of functions for
# plotting output from the UFS short range
# weather application.
# --------------------------------
import xarray as xr
import matplotlib.pyplot as plt
import scipy
import cartopy.crs as ccrs
import cartopy.feature as cpf
from datetime import datetime
import imageio
# --------------------------------
# begin defining functinons


def plot_conus():
    """
    Input desired forecast hour and simulation grib output, 
    outputs contour plot of desired variable across the CONUS
    :param str hr: desired forecast hour to be plotted
    :param str filepath: director in which grib output lies
    """
    # forecast hour to be plotted
    hr = int(input("Forecast hour to be plot: "))
    # directory where grib output is stored
    filepath = input("Directory where grib output are stored: ")
    # forecast initialization time
    init_str = input("Forecast start time (UTC): ")
    # natlev or prslev
    level_type = input("natlev or prslev: ")
    # grid 
    grid = input("Grid used: ")
    # full file path
    if hr < 10:
        file_name = f"rrfs.t{init_str}z.{level_type}.f00{hr}.{grid}.grib2"
    else:
        file_name = f"rrfs.t{init_str}z.{level_type}.f0{hr}.{grid}.grib2"
    # stepType
    if level_type == "prslev":
        typeOfLevel = input("What type of level are you using? Here are your options: \n meanSea \n hybrid \n atmosphereSingleLayer \n surface \n isothermal \n planetaryBoundaryLayer \n isobaricInhPa \n isobaricLayer \n heightAboveGround \n heightAboveGround \n depthBelowLandLayer \n depthBelowLand \n lowestLevelWetBulb0 \n nominalTop \n unknown \n lowCloudLayer \n middleCloudLayer \n highCloudLayer \n cloudBase \n cloudCeiling \n gridScaleCloudBottom \n cloudTop \n gridScaleCloudTop \n tropopause \n maxWind \n heightAboveSea \n isothermZero \n highestTroposphericFreezing \n pressureFromGroundLayer \n adiabaticCondensation \n \n Answer: ")
        var_name = input("What variable would you like plotted? See the following link for the list of variables available (https://ufs-srweather-app.readthedocs.io/en/develop/tables/SRW_PRSLEV_table.html) \n Answer: ")
    if level_type == "natlev":
        typeOfLevel = input("What type of level are you using? Here are your options: \n meanSea \n hybrid \n atmosphereSingleLayer \n surface \n isothermal \n planetaryBoundaryLayer \n isobaricInhPa \n isobaricLayer \n heightAboveGround \n heightAboveGroundLayer \n atmosphere \n depthBelowLandLayer \n depthBelowLand \n lowestLevelWetBulb0 \n nominalTop \n hybridLayer \n lowCloudLayer \n middleCloudLayer \n highCloudLayer \n cloudBase \n gridScaleCloudBottom \n cloudTop \n gridScaleCloudTop \n tropopause \n maxWind \n heightAboveSea \n isothermZero \n highestTroposphericFreezing \n pressureFromGroundLayer \n \n Answer: ")
        var_name = input("Variable being plotted? See the following link for the list of variables available (https://ufs-srweather-app.readthedocs.io/en/develop/tables/SRW_NATLEV_table.html) \n \n Answer: ")
    # read in grib output
    ds = xr.open_dataset(f"{filepath}/{file_name}", engine="cfgrib", filter_by_keys={'typeOfLevel': f'{typeOfLevel}'})
    # isolate needed variables/dimensions
    lat = ds.latitude
    lon = ds.longitude
    var = ds[var_name]
    # create datetime strings
    date_str = input("Date being plotted (YYYYMMDD): ")
    datetime_str = date_str + str(hr)
    date = datetime.strptime(datetime_str, "%Y%m%d%H")

    # create plot
    fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(10,10))
    # mapping
    ax.coastlines()
    ax.add_feature(cpf.BORDERS)
    ax.add_feature(cpf.STATES)
    # plot
    ax.contourf(lon, lat, var[0,:,:], transform=ccrs.PlateCarree())
    # variable name for plot title
    var_title = input("Full name of the variable being plotted: ")
    # variable units
    var_units = input("latex code for the units of the variable being plotted (including $ symbols): ")
    # forecast hour title
    fcst_hr_title = input("Forecast hour (0HH): ")
    plt.title(f"{var_title} [{var_units}], F0{fcst_hr_title} Valid at: {date} UTC")
    # save plot
    output_path = input("Path for figure to be output: ")
    fig_path = f"{output_path}/{var_name}.F{fcst_hr_title}.png"
    plt.savefig(fig_path)
    plt.close(fig)
# --------------------------------
