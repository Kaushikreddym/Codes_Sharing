import numpy as np
import xarray as xray
from glob import glob
from datetime import datetime
from netCDF4 import Dataset, MFDataset
from wrf import getvar, interplevel, to_np, latlon_coords, get_cartopy,     cartopy_xlim, cartopy_ylim, ALL_TIMES, geo_bounds
import cmor
import pandas as pd
import datetime as dt
from dateutil.relativedelta import relativedelta
import subprocess as sp
import sys
import numpy as np

def make_bounds(data):
    data_b = np.zeros(data.shape[0] + 1)
    if data.shape[0] == 1:
        data_b[0] = 0.0
        data_b[1] = 2 * data[0]
    else:

        # print(data[1], data[0])

        deldata = (float(data[1]) - float(data[0])) / 2

        # print(deldata)

        data_b[0] = data[0] - deldata
        deldata = (data[data.shape[0] - 1] - data[data.shape[0] - 2])             / 2
        data_b[data.shape[0]] = data[data.shape[0] - 1] + deldata
        for i in range(data.shape[0] - 1):
            data_b[i + 1] = (data[i] + data[i + 1]) / 2
    return data_b


def create_ilon(lon):  # creates longitude axis
    lonb = make_bounds(lon)
    ilon = cmor.axis(table_entry='longitude', units='degrees_east',
                     coord_vals=lon, cell_bounds=lonb)
    return ilon


def create_ilat(lat):  # creates latitude axis
    latb = make_bounds(lat)
    ilat = cmor.axis(table_entry='latitude', units='degrees_north',
                     coord_vals=lat, cell_bounds=latb)
    return ilat


def create_itime(time):  # creates time axis
    timeb = make_bounds(time)
    itime = cmor.axis('time', units='days since 2015')
    return itime


def create_ilevel(level):  # creates pressure axis
    levelb = make_bounds(level)
    ilevel = cmor.axis('plev27', coord_vals=level, units='Pa',
                       cell_bounds=levelb)
    return ilevel


def cmorize_2d(
    lat,
    lon,
    time,
    data,
    units,
    scale,
    cmorname,
    table,
    user_input='common_user_input.json'
    ):

    # 'cmorize_2d' takes arrays and corresponding dimensions to convert data into CF compliant format
    #

    cmor.setup(inpath='../cmor/TestTables',
               netcdf_file_action=cmor.CMOR_REPLACE_4)
    cmor.dataset_json('../cmor/'+user_input)
    cmor.load_table(table)  # 'table corresponds to the json table containing variable information'

    ilat = create_ilat(lat)
    ilon = create_ilon(lon)
    itime = create_itime(time)

    axes = [itime, ilat, ilon]
    ivar = cmor.variable(cmorname, units, axes,positive="up")

    timeb = make_bounds(time)

    cmor.write(ivar, data * scale, time_vals=time, time_bnds=timeb)
    filename = cmor.close(ivar, file_name=True)
    print ('stored in :', filename)
    cmor.close()


def cmorize_3d(
    lat,
    lon,
    time,
    level,
    data,
    units,
    scale,
    cmorname,
    table,
    user_input='common_user_input.json'
    ):

    # similar to the 'cmorize_2d' function but for 3d variables with added pressure dimension

    cmor.setup(inpath='../cmor/TestTables',
               netcdf_file_action=cmor.CMOR_REPLACE_4)
    cmor.dataset_json('../cmor/'+user_input)
    cmor.load_table(table)

    ilat = create_ilat(lat)
    ilon = create_ilon(lon)
    ilevel = create_ilevel(level)
    itime = create_itime(time)

    axes = [itime, ilevel, ilat, ilon]

    timeb = make_bounds(time)

    ivar = cmor.variable(cmorname, units, axes)
    cmor.write(ivar, data * scale, time_vals=time, time_bnds=timeb)
    filename = cmor.close(ivar, file_name=True)
    print ('stored in :', filename)
    cmor.close()


def read_cmor():
    flist = sorted(glob(REL_PATH + DATE + '/' + fname_prefix + '*'))[1:]
    ncfile = [Dataset(f) for f in flist]
    dset = getvar(ncfile, aerlist.loc[i]['varname'], timeidx=ALL_TIMES,
                  method='cat')


