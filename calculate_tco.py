
# coding: utf-8
#!/usr/bin/python

'''
This routine calculates total column from TOMCAT output files
author: Sandip Dhomse (sandip.dhomse@gmail.com)

Here we use O3 and air mass in a given grid box to get total column 
As this simulation does not have tropospheric chemistry- 30 DU is added everywhere

'''

import netCDF4 as nc
import numpy as np
import xarray as xr
import pandas as pd
import math
from datetime import datetime
import os

# Convert numpy array to Xarray
def make_xarray(odata,lats,lons,ntimes,vname,sdate,unit1):
    ## get the middle of the month date but does not work all the time
    time1 = pd.date_range(sdate, periods=ntimes, freq="1M")+ pd.DateOffset(days=-16)
    ds = xr.Dataset({vname: (['time', 'lat', 'lon'], odata.astype('float32'), {'units':unit1})},
                 coords={'lon': (['lon'], lons, {'units':'degrees_east'}),
                         'lat': (['lat'], lats, {'units':'degrees_north'}),
                         'time': (['time'], time1)
                })

    ds.lon.attrs['long_name'] = 'longitude'
    ds.lat.attrs['long_name'] = 'latitude'


    ds.attrs['Conventions'] = 'CF-1.7'
    ds.attrs['title'] = 'TOMCAT MPC741'
    ds.attrs['nc.institution'] = 'University of Leeds'
    ds.attrs['source'] = 'TOMCAT with ERA5'
    ds.attrs['history'] = str(datetime.utcnow()) + ' Python'
    ds.attrs['references'] = 'Dhomse et al., 2021, Dhomse et al., 2022'
    ds.attrs['comment'] = ''
    return ds

## Get surface area of each grid box
def get_area(nlons,nlats,lats):
    re=6.37e6
    dlon= math.radians(360./nlons)
    dlat = math.radians(180./nlats)
    sarea = np.zeros((nlats,nlons))
    for l in range(0,nlats):
        sarea[l,:] =re*re*math.cos(math.radians(lats[l]))*dlon*dlat
    return sarea
    


def get_tco(afile):
    COLFAC = 2.132E20
    grav = 9.80616#[ms-2] 
    vn1 = nc.MFDataset(afile)
   
    lons = vn1.variables["lon"][:]
    lats = vn1.variables["lat"][:]
    niv = vn1.variables["niv"][:]
    
    tracer = vn1.variables["O3_mm"]
    mass = vn1.variables["sm_mm"][:]
    time = vn1.variables["time"][:]
    
    nlons  =len(lons)
    nlats  =len(lats)
    ntimes = len(time)
    nlev = len(niv)
    sarea =get_area( nlons, nlats, lats)

    mass1 = mass*grav*COLFAC
    mass2 = 0.*mass
    
    for t in range(ntimes):
        for l in range(0,nlev):
            mass2[t,l,:,:] = (mass1[t,l,:,:]*tracer[t,l,:,:])/sarea[:,:]
    tdata =np.sum(mass2,axis=1)/(2.69*1.e16)
    ## Add 30 DU as tropospheric column
    tdata = tdata +30

    ds = xr.open_dataset(afile)
    sdate = str(ds.time.values[0])[0:10]
    vname = "tco"
    unit= 'DU'
    ## create a xarray
    ds= make_xarray(tdata,lats,lons, ntimes, vname, sdate,unit)
    return ds
   
## note that input file is in the same directory
## If input files are somewhere use glob module

fname = "TOMCAT_RUN741_testfile.nc"## input file
fname = "/nfs/a176/lecmc/MPC741/MPC741_t042_201901mm.nc"
## check if files exist and size is greater than 0

if os.path.exists(fname) and os.path.getsize(fname) > 0:
    xdata =get_tco(fname)
    #write in netcdf format
    xdata.to_netcdf("tco.nc")


# if you want to have quick look
xdata1 = xdata["tco"]
xdata1[1,:,:].plot(cmap='Reds')
