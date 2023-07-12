# --------------------------------
# Name: plotting_functions.py
# Author: Robert M. Frost
# NOAA Global Systems Laboratory
# Created: 16 June 2023
# Purpose: Function to read in
# grib model output and extract
# desired variable
# --------------------------------
import numpy as np
import xarray as xr
import pygrib

def read_grib(hr, dgrib, nat_prs, mesg_num, ret_type=0):
    """
    Purpose: Function to read in grib output from the UFS SRW app.
    :param int hr: forecast hour to be read
    :param str dgrib: directory where forecast output are located
    :param str nat_prs: reading in natlev or prslev
    :param int mesg_num: message number for the variable to be read in
    :param ret_type int: If 0, returns grbs, grb, lat, lon, and valid_time,
    if 1, returns just grbs, if 2, returns just valid_time
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
    if (ret_type == 0):
        return grbs, grb, lat, lon, valid_time
    if (ret_type == 1):
        return grb
    if (ret_type == 2):
        return valid_time
# --------------------------------
def skew_to_nc(hr, dgrib, dout):
    """
    Adapted from skew_to_new.py
    Inputs grb output for specified hour and outputs dataset 
    with 3d temperature, specific humidity, pressure, and u and v 
    wind components in base units for plotting skew-t's
    :param int hr: forecast hour of interest
    :param str dgrib: directory where grib output are found
    :param str dout: director for dataset to be output
    """
    print(f"Starting hour {hr}")
    # message number for T at isobaric levels in hPa
    mn_t = np.arange(14,794.1,20).astype(int)
    # message number for specific humidity
    mn_q = mn_t + 1
    # message number for u wind component
    mn_u = mn_t + 2
    # message number for v wind component
    mn_v = mn_t + 3
    # message number for pressure
    mn_p = np.arange(2,782.1,20)

    # read in output files
    grbs, T, lat, lon, valid_date = read_grib(hr, dgrib, "natlev", int(mn_t[0]))

    # create arrays for temperature
    T = np.empty((mn_t.size, lat[:,0].size, lon[0,:].size))
    # create arrays for specific humidity
    q = np.empty((mn_t.size, lat[:,0].size, lon[0,:].size))
    # create arrays for pressure
    p = np.empty((mn_t.size, lat[:,0].size, lon[0,:].size))
    # create arrays for u wind
    u = np.empty((mn_t.size, lat[:,0].size, lon[0,:].size))
    # create arrays for v wind
    v = np.empty((mn_t.size, lat[:,0].size, lon[0,:].size))

    # loop over hybrid levels
    for i in range(mn_t.size):
        # read in temperature
        T_2d = grbs[int(mn_t[i])]
        # store in array
        T[i] = T_2d.values

        # read in specific humidity
        q_2d = grbs[int(mn_q[i])]
        # store in array
        q[i] = q_2d.values

        # read in pressure
        p_2d = grbs[int(mn_p[i])]
        # store in array
        p[i] = p_2d.values

        # read in u component
        u_2d = grbs[int(mn_u[i])]
        # store in array
        u[i] = u_2d.values

        # read in v component
        v_2d = grbs[int(mn_v[i])]
        # store in array
        v[i] = v_2d.values

        print(f"Finished with hybrid level {i}")

    # output to xarray dataset
    ds = xr.Dataset()
    # add arrays to dataset
    ds['T'] = (['hybrid_level', 'latitude', 'longitude'], T)
    ds['q'] = (['hybrid_level', 'latitude', 'longitude'], q)
    ds['p'] = (['hybrid_level', 'latitude', 'longitude'], p)
    ds['u'] = (['hybrid_level', 'latitude', 'longitude'], u)
    ds['v'] = (['hybrid_level', 'latitude', 'longitude'], v)
    # add coordinates
    ds['hybrid_level'] = np.arange(mn_t.size)
    ds['north_south'] = (["latitude", "longitude"], lat)
    ds['west_east'] = (["latitude", "longitude"], lon)
    # add attributes
    ds['T'].attrs['units'] = 'K'
    ds['q'].attrs['units'] = 'kg/kg'
    ds['p'].attrs['units'] = 'Pa'
    ds['u'].attrs['units'] = 'm/s'
    ds['v'].attrs['units'] = 'm/s'
    ds['hybrid_level'].attrs['units'] = 'level'
    ds['latitude'].attrs['units'] = 'degrees_north'
    ds['longitude'].attrs['units'] = 'degrees_east'

    # save dataset
    fsave = f"{dout}skew_f{hr}.nc"
    print(f"Saving file: {fsave}")
    ds.to_netcdf(fsave)
    print(f"Finished with hour {hr} \n!")
    return