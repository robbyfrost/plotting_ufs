import xarray as xr
import numpy as np
from plotting_functions import read_grib

def skew_to_nc(hr, dgrib, dout):
    """
    Inputs grb output for specified hour and outputs dataset 
    with 3d temperature, dew point, pressure, and u and v 
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


# directory where hrrr grib data are located
dgrib = "/scratch2/BMC/fv3lam/Robby.Frost/expt_dirs/2023041900_3km_hrrrphys/2023041900/postprd/"
# directory for dataset to be output
dout = "/scratch2/BMC/fv3lam/Robby.Frost/skewt_data/20230419/hrrr/"

# loop over forecast hours
for hr in range(18,37):
    skew_to_nc(hr, dgrib, dout)