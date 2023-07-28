# --------------------------------
# Name: run_skew_to_nc.py
# Author: Robert M. Frost
# NOAA Global Systems Laboratory
# Created: 12 July 2023
# Purpose: Script to run skew_to_nc()
# from UFSutils.py
# --------------------------------
from UFSutils import skew_to_nc
# --------------------------------
# directory where hrrr grib data are located
dgrib = "/scratch2/BMC/fv3lam/Robby.Frost/expt_dirs/2023041912_3km_rapphys/2023041912/postprd/"
# directory for dataset to be output
dout = "/scratch2/BMC/fv3lam/Robby.Frost/skewt_data/2023041912/rap/"

# loop over forecast hours
for hr in range(25):
    skew_to_nc(hr, dgrib, dout)