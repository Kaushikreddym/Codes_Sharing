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

from cmor_func import *

#REL_PATH='/mnt/stime/scratch/PartclAeroResLab/kaushik.reddy.m/3.9.1/sam/run_cpu80/fc/'
fold = 'only_WRF10'
REL_PATH='/mnt/nas/DATA/kaushik/raw/'+fold+'/fc/'
fname_prefix = 'wrfout_d01_'
plevs = [
    100000.,
    97500.,
    95000.,
    92500.,
    90000.,
    87500.,
    85000.,
    82500.,
    80000.,
    77500.,
    75000.,
    70000.,
    65000.,
    60000.,
    55000.,
    50000.,
    45000.,
    40000.,
    35000.,
    30000.,
    25000.,
    22500.,
    20000.,
    17500.,
    15000.,
    12500.,
    10000.,
    ]


# In[ ]:


for mon_no in ['01','02']:#[ '01','02','03','04','05','06','07','08','09','10','11','12']:
    datelist = [sorted(glob(REL_PATH + '/2015' + mon_no + '*'
                ))[i].split('/')[-1] for i in
                range(len(sorted(glob(REL_PATH + '/2015' + mon_no + '*'
                ))))][:6]  # get the folders inside the fc folder

    flag_met = 0
    flag_rain = 0
    flag_aod = 0
    flag_soa = 0
    flag_aerosol = 0
    flag_cloud = 0
    flag_rad = 0
    flag_cdnc = 0
    flag_cf = 0
    flag_cloudvi = 1
    flag_tau= 0
    flag_pres = 0
    flag_gas = 0
    table = 'CMIP6_1hr.json'
    
    def get_attrs(dset):
        lon = dset.XLONG[0, :].data
        lat = dset.XLAT[:, 0].data
        time = np.array(dset.Time.data - np.datetime64('2015-01-01T00'
                        ), dtype=np.float) / (1e9 * 3600 * 24)
        level = np.array(plevs)
        return (lat, lon, level, time)


    if flag_rain == 1:
        dset_PR = []
        dset_PRCSH = []
        dset_PRRA = []
        dset_PRC = []
        for DATE in datelist:
            print(DATE)

            # flist = sorted(glob(REL_PATH+DATE+'/'+fname_prefix+'*'))

            flist = sorted(glob(REL_PATH + DATE + '/' + fname_prefix
                           + '*'))
            ncfile = [Dataset(f) for f in flist]
            dset_C = getvar(ncfile, 'RAINC', timeidx=ALL_TIMES,
                            method='cat')
            dset_SH = getvar(ncfile, 'RAINNC', timeidx=ALL_TIMES,
                             method='cat')
            dset_SC = getvar(ncfile, 'RAINSH', timeidx=ALL_TIMES,
                             method='cat')
            dset = dset_C + dset_SH.data + dset_SC.data

            dset_pr = dset[1:] - dset[:-1].data
            dset_sh = dset_SH[1:] - dset_SH[:-1].data
            dset_sc = dset_SC[1:] - dset_SC[:-1].data
            dset_c = dset_C[1:] - dset_C[:-1].data

            dset_PR.append(dset_pr)
            dset_PRCSH.append(dset_sc)
            dset_PRRA.append(dset_sh)
            dset_PRC.append(dset_c)

        dset_PR = xray.concat(dset_PR, 'Time')
        dset_PRCSH = xray.concat(dset_PRCSH, 'Time')
        dset_PRC = xray.concat(dset_PRC, 'Time')
        dset_PRRA = xray.concat(dset_PRRA, 'Time')
        # level = np.array(plevs)

        (lat, lon, level, time) = get_attrs(dset_PR)
        cmorize_2d(
            lat,
            lon,
            time,
            dset_PR.data * 3600 ** -1,
            'kg m-2 s-1',
            1,
            'pr',
            table,
            )
        cmorize_2d(
            lat,
            lon,
            time,
            dset_PRCSH.data * 3600 ** -1,
            'kg m-2 s-1',
            1,
            'prsh',
            table,
            )
        cmorize_2d(
            lat,
            lon,
            time,
            dset_PRC.data * 3600 ** -1,
            'kg m-2 s-1',
            1,
            'prrc',
            table,
            )
        cmorize_2d(
            lat,
            lon,
            time,
            dset_PRRA.data * 3600 ** -1,
            'kg m-2 s-1',
            1,
            'prcsh',
            table,
            )

    if flag_aod == 1:
        dset_DA = []
        dset_z = []
        for DATE in datelist:
            print(DATE)

            # flist = sorted(glob(REL_PATH+DATE+'/'+fname_prefix+'*'))

            flist = sorted(glob(REL_PATH + DATE + '/' + fname_prefix
                           + '*'))[1:]
            ncfile = [Dataset(f) for f in flist]
            dset = getvar(ncfile, 'EXTCOF55', timeidx=ALL_TIMES,
                          method='cat')
            dset_zstag = getvar(ncfile, 'zstag', timeidx=ALL_TIMES,
                                method='cat')
            thickness = dset_zstag[:, 1:].data - dset_zstag[:, :-1].data
            dset_AOD = dset * thickness / 1e3  # .sum(['bottom_top'])
            dset_DA.append(dset_AOD)
        dset_DA = xray.concat(dset_DA, 'Time')

        (lat, lon, level, time) = get_attrs(dset_DA)
        cmorize_2d(
            lat,
            lon,
            time,
            dset_DA.sum('bottom_top').data,
            '1',
            1,
            'od550aer',
            table,
            )
    if flag_met == 1:

        #"""

        aerlist = pd.read_csv('../cmor/MET6hourly_3d_varlist')
        for i in range(aerlist.shape[0]):
            dset_DA = []
            dset_z = []
            for DATE in datelist:

                # print(DATE)
                # flist = sorted(glob(REL_PATH+DATE+'/'+fname_prefix+'*'))

                flist = sorted(glob(REL_PATH + DATE + '/'
                               + fname_prefix + '*'))[1:]
                ncfile = [Dataset(f) for f in flist]
                dset = getvar(ncfile, aerlist.loc[i]['varname'],
                              timeidx=ALL_TIMES, method='cat')
                pres = getvar(ncfile, 'pressure', timeidx=ALL_TIMES,
                              method='cat')

                # tk = getvar(ncfile,"tk", timeidx=ALL_TIMES, method="cat")

                dset_interp = interplevel(dset, pres, np.array(plevs)
                        / 100)
                dset_DA.append(dset_interp)
            dset_DA = xray.concat(dset_DA, 'Time')
            if aerlist.loc[i]['cmorname'] == 'zg':
                units = 'm'
            else:
                units = dset.units

            (lat, lon, level, time) = get_attrs(dset_DA)
            cmorize_3d(
                lat,
                lon,
                time,
                level,
                dset_DA.data,
                units,
                aerlist.loc[i]['scale'],
                aerlist.loc[i]['cmorname'],
                table,
                )

        #"""

        aerlist = pd.read_csv('../cmor/MET6hourly_2d_varlist')
        for i in range(aerlist.shape[0]):
            dset_DA = []
            for DATE in datelist:

                print(DATE)
                # flist = sorted(glob(REL_PATH+DATE+'/'+fname_prefix+'*'))

                flist = sorted(glob(REL_PATH + DATE + '/'
                               + fname_prefix + '*'))[1:]
                ncfile = [Dataset(f) for f in flist]
                if aerlist.loc[i]['cmorname'] == 'uas':
                    dset = getvar(ncfile, aerlist.loc[i]['varname'],
                                  timeidx=ALL_TIMES, method='cat')[0]
                elif aerlist.loc[i]['cmorname'] == 'vas':
                    dset = getvar(ncfile, aerlist.loc[i]['varname'],
                                  timeidx=ALL_TIMES, method='cat')[1]
                elif aerlist.loc[i]['cmorname'] == 'sfcWind':
                    dset = getvar(ncfile, aerlist.loc[i]['varname'],
                                  timeidx=ALL_TIMES, method='cat')[0]
                else:
                    dset = getvar(ncfile, aerlist.loc[i]['varname'],
                                  timeidx=ALL_TIMES, method='cat')

                # pres = getvar(ncfile,"pressure", timeidx=ALL_TIMES, method="cat")
                # tk = getvar(ncfile,"tk", timeidx=ALL_TIMES, method="cat")
                # dset_interp = interplevel(dset,pres,np.array(plevs)/100)

                dset_DA.append(dset)
            dset_DA = xray.concat(dset_DA, 'Time')
            (lat, lon, level, time) = get_attrs(dset_DA)
            cmorize_2d(
                lat,
                lon,
                time,
                dset_DA.data,
                aerlist.loc[i]['units'],
                aerlist.loc[i]['scale'],
                aerlist.loc[i]['cmorname'],
                table,
                )

    if flag_aerosol == 1:
        #'''
        aerlist = pd.read_csv('../cmor/AER6hourly_3d_varlist')
        for i in range(aerlist.shape[0]):
            dset_DA = []
            dset_z = []
            for DATE in datelist:
                print(DATE)

                # flist = sorted(glob(REL_PATH+DATE+'/'+fname_prefix+'*'))

                flist = sorted(glob(REL_PATH + DATE + '/'
                               + fname_prefix + '*'))[1:]
                ncfile = [Dataset(f) for f in flist]
                dset = getvar(ncfile, aerlist.loc[i]['varname'] + 'i',
                              timeidx=ALL_TIMES, method='cat') \
                    + getvar(ncfile, aerlist.loc[i]['varname'] + 'j',
                             timeidx=ALL_TIMES, method='cat').data
                pres = getvar(ncfile, 'pressure', timeidx=ALL_TIMES,
                              method='cat')

                # tk = getvar(ncfile,"tk", timeidx=ALL_TIMES, method="cat")

                dset_interp = interplevel(dset, pres, np.array(plevs)
                        / 100)
                dset_DA.append(dset_interp)
            dset_DA = xray.concat(dset_DA, 'Time')
            (lat, lon, level, time) = get_attrs(dset_DA)
            cmorize_3d(
                lat,
                lon,
                time,
                level,
                dset_DA.data,
                'kg kg-1',
                aerlist.loc[i]['scale'],
                aerlist.loc[i]['cmorname'],
                table,
                )

        aerlist = pd.read_csv('../cmor/AER6hourly_2d_varlist')
        for i in range(aerlist.shape[0]):
            print(aerlist.loc[i]['cmorname'])
            dset_DA = []
            for DATE in datelist:

                # print(DATE)
                # flist = sorted(glob(REL_PATH+DATE+'/'+fname_prefix+'*'))

                flist = sorted(glob(REL_PATH + DATE + '/'
                               + fname_prefix + '*'))[1:]
                ncfile = [Dataset(f) for f in flist]
                dset = getvar(ncfile, aerlist.loc[i]['varname'] + 'i',
                              timeidx=ALL_TIMES, method='cat') \
                    + getvar(ncfile, aerlist.loc[i]['varname'] + 'j',
                             timeidx=ALL_TIMES, method='cat').data
                pres = getvar(ncfile, 'pressure', timeidx=ALL_TIMES,
                              method='cat')
                tk = getvar(ncfile, 'tk', timeidx=ALL_TIMES,
                            method='cat')
                rho = pres.data * 1e2 * 28.979 * 1e-3 / (8.314
                        * tk.data)
                dset_zstag = getvar(ncfile, 'zstag', timeidx=ALL_TIMES,
                                    method='cat')
                thickness = dset_zstag[:, 1:].data - dset_zstag[:, :
                        -1].data
                dset_sum = (dset * thickness * rho).sum(['bottom_top'])
                dset_DA.append(dset_sum)
            dset_DA = xray.concat(dset_DA, 'Time')
            (lat, lon, level, time) = get_attrs(dset_DA)
            cmorize_2d(
                lat,
                lon,
                time,
                dset_DA.data,
                'kg m-2',
                aerlist.loc[i]['scale'],
                aerlist.loc[i]['cmorname'],
                table,
                )

        aerlist = pd.read_csv('../cmor/AER6hourly_surf_varlist')

        for i in range(aerlist.shape[0]):
            print(aerlist.loc[i]['cmorname'])
            dset_DA = []
            for DATE in datelist:

                # print(DATE)
                # flist = sorted(glob(REL_PATH+DATE+'/'+fname_prefix+'*'))

                flist = sorted(glob(REL_PATH + DATE + '/'
                               + fname_prefix + '*'))[1:]
                ncfile = [Dataset(f) for f in flist]
                dset = getvar(ncfile, aerlist.loc[i]['varname'] + 'i',
                              timeidx=ALL_TIMES, method='cat') \
                    + getvar(ncfile, aerlist.loc[i]['varname'] + 'j',
                             timeidx=ALL_TIMES, method='cat').data
                pres = getvar(ncfile, 'pressure', timeidx=ALL_TIMES,
                              method='cat')
                tk = getvar(ncfile, 'tk', timeidx=ALL_TIMES,
                            method='cat')
                rho = pres.data * 1e2 * 28.979 * 1e-3 / (8.314
                        * tk.data)
                dset_zstag = getvar(ncfile, 'zstag', timeidx=ALL_TIMES,
                                    method='cat')
                thickness = dset_zstag[:, 1:].data - dset_zstag[:, :
                        -1].data
                dset = (dset * rho).sel(bottom_top=0)
                dset_DA.append(dset)
            dset_DA = xray.concat(dset_DA, 'Time')
            (lat, lon, level, time) = get_attrs(dset_DA)
            cmorize_2d(
                lat,
                lon,
                time,
                dset_DA.data,
                'kg m-3',
                aerlist.loc[i]['scale'],
                aerlist.loc[i]['cmorname'],
                table,
                )
        #'''
        dset_DA = []
        dset_surf = []
        dset_load = []
        for DATE in datelist:
            print(DATE)

            # flist = sorted(glob(REL_PATH+DATE+'/'+fname_prefix+'*'))

            flist = sorted(glob(REL_PATH + DATE + '/' + fname_prefix
                           + '*'))[1:]
            ncfile = [Dataset(f) for f in flist]
            dset = getvar(ncfile, 'PM2_5_DRY', timeidx=ALL_TIMES,
                          method='cat')
            pres = getvar(ncfile, 'pressure', timeidx=ALL_TIMES,
                          method='cat')
            dset_zstag = getvar(ncfile, 'zstag', timeidx=ALL_TIMES,
                          method='cat')
            thickness = dset_zstag[:, 1:].data - dset_zstag[:, :
                        -1].data
            
            dset_interp = interplevel(dset, pres, np.array(plevs) / 100)
            
            dset_DA.append(dset_interp)
            dset_surf.append(dset.sel(bottom_top=0))
            dset_load.append((dset*thickness).sum('bottom_top'))
        
        dset_DA = xray.concat(dset_DA, 'Time')
        (lat, lon, level, time) = get_attrs(dset_DA)

        cmorize_3d(
            lat,
            lon,
            time,
            level,
            dset_DA.data,
            'kg m-3',
            1e-9,
            'mcpm2p5',
            table,
            )

        dset_surf = xray.concat(dset_surf, 'Time')
        cmorize_2d(
            lat,
            lon,
            time,
            dset_surf.data,
            'kg m-3',
            1e-9,
            'sconcpm2p5',
            table,
            )
        dset_load = xray.concat(dset_load, 'Time')
        cmorize_2d(
            lat,
            lon,
            time,
            dset_load.data,
            'kg m-2',
            1e-9,
            'loadpm2p5',
            table,
            )
        #'''
        dset_DA = []
        dset_surf = []
        dset_load = []
        for DATE in datelist:
            print(DATE)

            # flist = sorted(glob(REL_PATH+DATE+'/'+fname_prefix+'*'))

            flist = sorted(glob(REL_PATH + DATE + '/' + fname_prefix
                           + '*'))[1:]
            ncfile = [Dataset(f) for f in flist]
            dset = getvar(ncfile, 'PM2_5_WATER', timeidx=ALL_TIMES,
                          method='cat')
            pres = getvar(ncfile, 'pressure', timeidx=ALL_TIMES,
                          method='cat')
            dset_zstag = getvar(ncfile, 'zstag', timeidx=ALL_TIMES,
                          method='cat')
            thickness = dset_zstag[:, 1:].data - dset_zstag[:, :
                        -1].data
            
            dset_interp = interplevel(dset, pres, np.array(plevs) / 100)
            
            dset_DA.append(dset_interp)
            dset_surf.append(dset.sel(bottom_top=0))
            dset_load.append((dset*thickness).sum('bottom_top'))
        
        dset_DA = xray.concat(dset_DA, 'Time')
        (lat, lon, level, time) = get_attrs(dset_DA)

        cmorize_3d(
            lat,
            lon,
            time,
            level,
            dset_DA.data,
            'kg m-3',
            1e-9,
            'mcaerh2o',
            table,
            )

        dset_surf = xray.concat(dset_surf, 'Time')
        cmorize_2d(
            lat,
            lon,
            time,
            dset_surf.data,
            'kg m-3',
            1e-9,
            'sconcaerh2o',
            table,
            )
        dset_load = xray.concat(dset_load, 'Time')
        cmorize_2d(
            lat,
            lon,
            time,
            dset_load.data,
            'kg m-2',
            1e-9,
            'loadaerh2o',
            table,
            )
    #'''

    if flag_soa == 1:

        dset_DA = []
        dset_z = []
        dset_surf = []
        for DATE in datelist:
            print(DATE)

            # flist = sorted(glob(REL_PATH+DATE+'/'+fname_prefix+'*'))

            flist = sorted(glob(REL_PATH + DATE + '/' + fname_prefix
                           + '*'))[1:]
            ncfile = [Dataset(f) for f in flist]
            soa = []
            #for aero in ['asoa1', 'asoa2', 'asoa3', 'asoa4','bsoa1', 'bsoa2', 'bsoa3', 'bsoa4']:
            for aero in ['orgaro1','orgaro2','orgalk1','orgole1','orgba1','orgba2','orgba3','orgba4']:
                dset = getvar(ncfile, aero + 'i', timeidx=ALL_TIMES,
                              method='cat') + getvar(ncfile, aero + 'j'
                        , timeidx=ALL_TIMES, method='cat').data
                soa.append(dset)
            soa = xray.concat(soa, 'soa').sum('soa')
            pres = getvar(ncfile, 'pressure', timeidx=ALL_TIMES,
                          method='cat')
            tk = getvar(ncfile, 'tk', timeidx=ALL_TIMES, method='cat')
            rho = pres.data * 1e2 * 28.979 * 1e-3 / (8.314 * tk.data)
            dset_zstag = getvar(ncfile, 'zstag', timeidx=ALL_TIMES,
                                method='cat')
            thickness = dset_zstag[:, 1:].data - dset_zstag[:, :-1].data
            dset_sum = (soa * thickness * rho).sum(['bottom_top'])

            # tk = getvar(ncfile,"tk", timeidx=ALL_TIMES, method="cat")

            dset_interp = interplevel(soa, pres, np.array(plevs) / 100)

            dset_DA.append(dset_sum)
            dset_z.append(dset_interp)
            dset_surf.append((soa * rho).sel(bottom_top=0))

        
        dset_z = xray.concat(dset_z, 'Time')
        (lat, lon, level, time) = get_attrs(dset_z)
        cmorize_3d(
            lat,
            lon,
            time,
            level,
            dset_z.data,
            'kg kg-1',
            1e-9,
            'mmrsoa',
            table,
            )

        dset_DA = xray.concat(dset_DA, 'Time')
        cmorize_2d(
            lat,
            lon,
            time,
            dset_DA.data,
            'kg m-2',
            1e-9,
            'loadsoa',
            table,
            )

        dset_surf = xray.concat(dset_surf, 'Time')
        cmorize_2d(
            lat,
            lon,
            time,
            dset_surf.data,
            'kg m-3',
            1e-9,
            'sconcsoa',
            table,
            )
    if flag_cloud == 1:

        aerlist = pd.read_csv('../cmor/CLOUD6hourly_3d_varlist')

        for i in range(aerlist.shape[0]):
            dset_DA = []
            dset_z = []
            for DATE in datelist:
                print(DATE)
                flist = sorted(glob(REL_PATH + DATE + '/'
                               + fname_prefix + '*'))

                # flist = sorted([glob(REL_PATH+DATE+'/'+fname_prefix+'*'+i+':00:00*') for i in ['00','06','12','18']])[1:]

                ncfile = [Dataset(f) for f in flist[1:]]
                dset = getvar(ncfile, aerlist.loc[i]['varname'],
                              timeidx=ALL_TIMES, method='cat')
                pres = getvar(ncfile, 'pressure', timeidx=ALL_TIMES,
                              method='cat')

                # tk = getvar(ncfile,"tk", timeidx=ALL_TIMES, method="cat")

                dset_interp = interplevel(dset, pres, np.array(plevs)
                        / 100)
                dset_DA.append(dset_interp)
            dset_DA = xray.concat(dset_DA, 'Time')

            (lat, lon, level, time) = get_attrs(dset_DA)

            data_cloud = dset_DA.data
            scale = 1
            cmorname = aerlist.loc[i]['cmorname']
            units = 'kg kg-1'
            table = 'CMIP6_cloud_1hr.json'
            cmorize_3d(
                lat,
                lon,
                time,
                level,
                data_cloud,
                units,
                scale,
                cmorname,
                table,
                )

        #aerlist = pd.read_csv('../cmor/CLOUD6hourly_2d_varlist')

        #for i in range(aerlist.shape[0]):
        #print(aerlist.loc[i]['cmorname'])
    if flag_cloudvi == 1:
        dset_DA = []
        dset_z = []
        for DATE in datelist:
            print(DATE)
            flist = sorted(glob(REL_PATH + DATE + '/'
                           + fname_prefix + '*'))

            # flist = sorted([glob(REL_PATH+DATE+'/'+fname_prefix+'*'+i+':00:00*') for i in ['00','06','12','18']])[1:]

            ncfile = [Dataset(f) for f in flist[1:]]
            dsetvar=[]
            for varname in ['QVAPOR','QCLOUD','QRAIN','QICE','QSNOW','QGRAUP']:
                dset = getvar(ncfile,varname,timeidx=ALL_TIMES, method='cat')
                dsetvar.append(dset)
            dset = xray.concat(dsetvar,pd.Index(['QVAPOR','QCLOUD','QRAIN','QICE','QSNOW','QGRAUP'],name='var'))
            
            dset_prw = dset.sel(var='QVAPOR').drop('var')
            dset_clwvi = dset.sel(var=['QCLOUD','QRAIN']).sum('var')
            dset_clivi = dset.sel(var=['QICE','QSNOW','QGRAUP']).sum('var')
            
            dset = xray.concat([dset_prw,dset_clwvi,dset_clivi],pd.Index(['prw','clwvi','clivi'],name='var'))
            
            #METHOD 1
            #pres = getvar(ncfile, 'pressure', timeidx=ALL_TIMES,
            #              method='cat')
            #tk = getvar(ncfile, 'tk', timeidx=ALL_TIMES,
            #            method='cat')
            #rho = pres.data * 1e2 * 28.979 * 1e-3 / (8.314
            #        * tk.data)
            #dset_zstag = getvar(ncfile, 'zstag', timeidx=ALL_TIMES,
            #                    method='cat')
            #thickness = dset_zstag[:, 1:].data - dset_zstag[:, :
            #        -1].data
            #dset_sum = (dset * thickness * rho).sum(['bottom_top'])
            
            #METHOD2

            mass = getvar(ncfile,'MU',timeidx=ALL_TIMES,method='cat')+getvar(ncfile,'MUB',timeidx=ALL_TIMES,method='cat')
            dnw = -1*getvar(ncfile,'DNW',timeidx=ALL_TIMES,method='cat')
            dset_sum = ((dset*dnw).sum('bottom_top'))*(mass/9.81)

            dset_DA.append(dset_sum)
        dset_DA = xray.concat(dset_DA, 'Time')
        (lat, lon, level, time) = get_attrs(dset_DA)
        scale = 1
        for cmorname in ['prw','clwvi','clivi']:
            units = "kg m-2"
            table = 'CMIP6_cloud_1hr.json'
            cmorize_2d(
                lat,
                lon,
                time,
                dset_DA.sel(var=cmorname).data,
                units,
                scale,
                cmorname,
                table,
                )
    if flag_rad == 1:
        aerlist = pd.read_csv('../cmor/RAD6hourly_2d_varlist')
        for i in range(aerlist.shape[0]):
            dset_DA = []
            for DATE in datelist:
                print(DATE)
                flist = sorted(glob(REL_PATH + DATE + '/'
                               + fname_prefix + '*'))[1:]
                ncfile = [Dataset(f) for f in flist]
                dset = getvar(ncfile, aerlist.loc[i]['varname'],
                              timeidx=ALL_TIMES, method='cat')
                dset_DA.append(dset)
            dset_DA = xray.concat(dset_DA, 'Time')
            (lat, lon, level, time) = get_attrs(dset_DA)
            scale = 1
            cmorname = aerlist.loc[i]['cmorname']
            units = aerlist.loc[i]['units']
            table = 'CMIP6_1hr_radiation.json'
            cmorize_2d(
                lat,
                lon,
                time,
                dset_DA.data,
                aerlist.loc[i]['units'],
                scale,
                aerlist.loc[i]['cmorname'],
                table,
                )
    if flag_cdnc == 1: 
        dset_DA = []
        dset_z = []
        for DATE in datelist:
            print(DATE)
            flist = sorted(glob(REL_PATH + DATE + '/'
                           + fname_prefix + '*'))

            # flist = sorted([glob(REL_PATH+DATE+'/'+fname_prefix+'*'+i+':00:00*') for i in ['00','06','12','18']])[1:]

            ncfile = [Dataset(f) for f in flist[1:]]
            dset = getvar(ncfile, 'QNDROP',
                          timeidx=ALL_TIMES, method='cat')
            pres = getvar(ncfile, 'pressure', timeidx=ALL_TIMES,
                          method='cat')
            tk = getvar(ncfile, 'tk', timeidx=ALL_TIMES,
                        method='cat')
            rho = pres.data * 1e2 * 28.979 * 1e-3 / (8.314
                    * tk.data)
            

            dset_interp = interplevel(dset*rho.data, pres, np.array(plevs)
                    / 100)
            dset_DA.append(dset_interp)
        dset_DA = xray.concat(dset_DA, 'Time')

        (lat, lon, level, time) = get_attrs(dset_DA)

        data_cloud = dset_DA.data
        scale = 1
        cmorname = 'cdnc'
        units = 'm-3'
        table = 'CMIP6_cloud_1hr.json'
        cmorize_3d(
            lat,
            lon,
            time,
            level,
            data_cloud,
            units,
            scale,
            cmorname,
            table,
            ) 
    if flag_tau == 1:
        dset_DA = []
        dset_z = []
        for DATE in datelist:
            print(DATE)
            flist = sorted(glob(REL_PATH + DATE + '/'
                           + fname_prefix + '*'))

            # flist = sorted([glob(REL_PATH+DATE+'/'+fname_prefix+'*'+i+':00:00*') for i in ['00','06','12','18']])[1:]
            #from scipy.special import Gamma
            ncfile = [Dataset(f) for f in flist[1:]]
            dset_Nd = getvar(ncfile, 'QNDROP',
                          timeidx=ALL_TIMES, method='cat')
            dset_qc = getvar(ncfile, 'QCLOUD',
                          timeidx=ALL_TIMES, method='cat')
            pres = getvar(ncfile, 'pressure', timeidx=ALL_TIMES,
                          method='cat')
            tk = getvar(ncfile, 'tk', timeidx=ALL_TIMES,
                        method='cat')
            rho = pres.data * 1e2 * 28.979 * 1e-3 / (8.314
                    * tk.data)
            

            dset_interp = interplevel(dset*rho.data, pres, np.array(plevs)
                    / 100)
            dset_DA.append(dset_interp)
        dset_DA = xray.concat(dset_DA, 'Time')

        (lat, lon, level, time) = get_attrs(dset_DA)

        data_cloud = dset_DA.data
        scale = 1
        cmorname = 'cdnc'
        units = 'm-3'
        table = 'CMIP6_cloud_1hr.json'
        cmorize_3d(
            lat,
            lon,
            time,
            level,
            data_cloud,
            units,
            scale,
            cmorname,
            table,
            )
    if flag_cf == 1: 
        dset_DA = []
        dset_z = []
        for DATE in datelist:
            print(DATE)
            flist = sorted(glob(REL_PATH + DATE + '/'
                           + fname_prefix + '*'))

            # flist = sorted([glob(REL_PATH+DATE+'/'+fname_prefix+'*'+i+':00:00*') for i in ['00','06','12','18']])[1:]

            ncfile = [Dataset(f) for f in flist[1:]]
            dset = getvar(ncfile, 'CLDFRA',
                          timeidx=ALL_TIMES, method='cat')
            pres = getvar(ncfile, 'pressure', timeidx=ALL_TIMES,
                          method='cat')

            dset_interp = interplevel(dset, pres, np.array(plevs)
                    / 100)
            dset_DA.append(dset_interp)
        dset_DA = xray.concat(dset_DA, 'Time')

        (lat, lon, level, time) = get_attrs(dset_DA)

        scale = 1
        cmorname = 'cldfra'
        units = '1'
        table = 'CMIP6_cloud_1hr.json'
        cmorize_3d(
            lat,
            lon,
            time,
            level,
            dset_DA.data,
            units,
            scale,
            cmorname,
            table,
            )
        dset_frac = []
        for DATE in datelist:
            print(DATE)
            # flist = sorted(glob(REL_PATH+DATE+'/'+fname_prefix+'*'))
            flist = sorted(glob(REL_PATH + DATE + '/' + fname_prefix
                           + '*'))[1:]
            ncfile = [Dataset(f) for f in flist]
            dset = getvar(ncfile, 'cloudfrac', timeidx=ALL_TIMES,
                                method='cat',vert_type='pres',low_thresh=97000,mid_thresh=68000,high_thresh=40000)
            dset_frac.append(dset)
        dset_frac = xray.concat(dset_frac, 'Time')

        scale = 1
        units = '%'
        table = 'CMIP6_cloud_1hr.json'

        (lat, lon, level, time) = get_attrs(dset_frac)

        cmorname = 'clh'
        cldlev = 'high'
        cmorize_2d(
            lat,
            lon,
            time,
            dset_frac.sel(low_mid_high = cldlev).data,
            units,
            scale,
            cmorname,
            table,
            )
        cmorname = 'clm'
        cldlev = 'mid'
        cmorize_2d(
            lat,
            lon,
            time,
            dset_frac.sel(low_mid_high = cldlev).data,
            units,
            scale,
            cmorname,
            table,
            )
        cmorname = 'cll'
        cldlev = 'low'
        cmorize_2d(
            lat,
            lon,
            time,
            dset_frac.sel(low_mid_high = cldlev).data,
            units,
            scale,
            cmorname,
            table,
            )
    if flag_pres == 1: 
        dset_DA = []
        dset_z = []
        for DATE in datelist:
            print(DATE)
            flist = sorted(glob(REL_PATH + DATE + '/'
                           + fname_prefix + '*'))

            # flist = sorted([glob(REL_PATH+DATE+'/'+fname_prefix+'*'+i+':00:00*') for i in ['00','06','12','18']])[1:]

            ncfile = [Dataset(f) for f in flist[1:]]
            pres = getvar(ncfile, 'pressure', timeidx=ALL_TIMES,
                          method='cat')
            dset_interp = interplevel(pres, pres, np.array(plevs)
                    / 100)
            dset_DA.append(dset_interp)
        dset_DA = xray.concat(dset_DA, 'Time')

        (lat, lon, level, time) = get_attrs(dset_DA)

        data_cloud = dset_DA.data
        scale = 1
        cmorname = 'pfull'
        units = 'Pa'
        table = 'CMIP6_1hr.json'
        cmorize_3d(
            lat,
            lon,
            time,
            level,
            data_cloud,
            units,
            scale,
            cmorname,
            table,
            )
    if flag_gas == 1:
        varname_plev=pd.read_csv('../cmor/varname_plev.csv')
        varname_sconc=pd.read_csv('../cmor/varname_sconc.csv')
        for i in range(varname_sconc.shape[0]):
            dset_load = []
            dset_sconc = []
            dset_plev = []
            
            if varname_plev.loc[i]['varname']!='nan':

                for DATE in datelist: # put datelist[:n] to process first n files
                    print(DATE) 
                    flist = sorted(glob(REL_PATH+DATE+'/'+fname_prefix+'*'))
                    ncfile = [Dataset(f) for f in flist[1:]] # 1 neglects the 1st file of each 12 hour segment
                    dset = getvar(ncfile,varname_plev.loc[i]['varname'], timeidx=ALL_TIMES, method="cat")
                    lat = getvar(ncfile,"latitude", timeidx=ALL_TIMES, method="cat")
                    lon = getvar(ncfile,"longitude", timeidx=ALL_TIMES, method="cat")

                    pres = getvar(ncfile,"pressure", timeidx=ALL_TIMES, method="cat")
                    tk = getvar(ncfile,"tk", timeidx=ALL_TIMES, method="cat")
                    mod_height=getvar(ncfile,"zstag", timeidx=ALL_TIMES, method="cat")

                    mod_thick = mod_height[:,1:].data-mod_height[:,:-1].data #calculating thickness of each level
                    molec_vol = 8.314*tk.data/(pres.data*100)  # m3 mol-1    (molar volume) 

                    dset_col = dset*6.022e23*mod_thick/molec_vol.data 
                    dset_load.append(dset_col.sum('bottom_top'))       

                    dset_interp = interplevel(dset,pres,np.array(plevs)/100)          
                    dset_plev.append(dset_interp) 

                    dset_surf = dset*varname_sconc.loc[i]['mweight']*1e-3/molec_vol
                    dset_surf = dset_surf.sel(bottom_top=0)
                    dset_sconc.append(dset_surf)

                dset_load=xray.concat(dset_load,'Time')    
                dset_plev=xray.concat(dset_plev,'Time')
                dset_sconc=xray.concat(dset_sconc,'Time')
                (lat, lon, level, time) = get_attrs(dset_plev)
                scale = 1e-6
                table='CMIP6_1hr_gaseous.json'
                units_plev='mol mol-1'
                cmorname_plev = varname_plev.loc[i]['cmorname']
                cmorize_3d(lat,lon,time,level,dset_plev.data,units_plev,scale,cmorname_plev,table)

                units_load='molec m-2'
                cmorname_load = 'load' + varname_plev.loc[i]['cmorname']
                cmorize_2d(lat,lon,time,dset_load.data,units_load,scale,cmorname_load,table)

                units_sconc ='kg m-3'
                cmorname_sconc = varname_sconc.loc[i]['cmorname']
                cmorize_2d(lat,lon,time,dset_sconc.data,units_sconc,scale,cmorname_sconc,table)
            else:
                print(varname_plev.loc[i]['varname']+' is not there in RACM mechanism')

    import time as tmod
    seconds = tmod.time()
    local_time = tmod.ctime(seconds)
    print('Local time:', local_time)


# In[ ]:




