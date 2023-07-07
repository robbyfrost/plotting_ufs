# --------------------------------
# Name: plotting_functions.py
# Author: Robert M. Frost
# NOAA Global Systems Laboratory
# Created: 16 June 2023
# Purpose: Function to read in
# grib model output and extract
# desired variable
# --------------------------------
import matplotlib.pyplot as plt
import pygrib
from metpy.plots import ctables
from matplotlib import rc
import cartopy.crs as ccrs
import cartopy.feature as cpf

def read_grib(hr, dgrib, nat_prs, mesg_num, array_only=False):
    """
    Purpose: Function to read in grib output from the UFS SRW app.
    :param int hr: forecast hour to be read
    :param str dgrib: directory where forecast output are located
    :param str nat_prs: reading in natlev or prslev
    :param int mesg_num: message number for the variable to be read in
    :param array_only bool: returns just variable array if true, returns variable array, lat lon, and valid_time if false
    """
    # set filenames
    if hr < 10:
        dgrib = f"{dgrib}rrfs.t00z.{nat_prs}.f00{hr}.rrfs_conuscompact_3km.grib2"
    else:
        dgrib = f"{dgrib}rrfs.t00z.{nat_prs}.f0{hr}.rrfs_conuscompact_3km.grib2"

    # open hrrr and rap output
    print(f"Reading in {dgrib}")
    grbs = pygrib.open(dgrib)
    # extract variable of interest
    grb = grbs[mesg_num]
    # extract latitude and longitude arrays
    lat, lon = grb.latlons()
    # extract datatime
    valid_time = grb.validDate
    
    print(f"Finished reading in {grb.name}")
    if array_only:
        return grb
    else:
        return grb, lat, lon, valid_time