
# coding: utf-8
#!/usr/bin/python

'''
A Code to interpolate TOMCAT tracers on MLS pressure levels
author: Sandip Dhomse (sandip.dhomse@gmail.com)
'''

import netCDF4 as nc
import numpy as np
import xarray as xr
import pandas as pd
from datetime import datetime 
import os

# Convert numpy array to Xarray
def make_xarray(odata,levs,lats,lons,ntimes,vname,sdate,unit1):
    
    time1 = pd.date_range(sdate, periods=ntimes, freq="1M")+ pd.DateOffset(days=-16)
    ds = xr.Dataset({vname: (['time', 'plev', 'lat', 'lon'], odata.astype('float32'), {'units':unit1})},
                 coords={'lon': (['lon'], lons, {'units':'degrees_east'}),
                         'lat': (['lat'], lats, {'units':'degrees_north'}),
                         'plev': (['plev'], levs, {'units':'hPa'},{'positive':'down'}),
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


## sometimes model grid point pressure are larger than 1000hPa or lesser than 0.1hPa

def revise_ax(nlev,ax,ay):
    nax = np.zeros(nlev+2)
    nay = np.zeros(nlev+2)
    nax[0]= 0.001## lets fix very small pressure for top of the atmosphere
    nax[1:nlev+1] = ax[:]
    nax[nlev+1] = 1200.## some high pressure at the bottom
    
    nay[0]= 0.8*ay[0]
    nay[1:nlev+1] = ay[:]
    nay[nlev+1] = 0.001*ay[nlev-1]
    return nax,nay



def read_and_interpolate(afile,mlsp,vname,unit):
    vn1 = nc.MFDataset(afile)
    #print(vn1)
    lons = vn1.variables["lon"][:]
    lats = vn1.variables["lat"][:]
    tracer = vn1.variables[vname+"_mm"][:]
    pres = vn1.variables["p_mm"][:]
    time = vn1.variables["time"][:]
    #clo =0
    nlon  =len(lons)
    nlat  =len(lats)
    ntimes = np.size(time)
    nlev = np.size(tracer[0,:,0,0])
    
    adata = np.zeros((ntimes,np.size(mlsp),nlat,nlon))
    
    for itime in range(ntimes):
        for ilat in range(0,nlat):
            for ilon in range(0,nlon):
                ay = tracer[itime,:,ilat,ilon]
                ax = pres[itime,:,ilat,ilon]/100.## hPa from Pa
                
                nax,nay = revise_ax(nlev,ax,ay)
                newy =np.interp(np.log(mlsp),np.log(nax),nay)
                adata[itime,:,ilat,ilon] = newy*1.e9## convert ppb 
    
    ds = xr.open_dataset(afile)
    sdate = str(ds.time.values[0])[0:10]
    ## create a xarray
    ds = make_xarray(adata,mlsp, lats, lons, ntimes, vname, sdate, unit)
    return ds

## create mls level array
ll = np.arange(0,25)## 25 pressure levels 1000 hPa to 0.1hPa
mlsp = np.round(1000*10**(-ll/6.),6)

## test file has only three tracers 
## normal file would have many more..
## give tracer names below

vnames = ["O3","HCL","HNO3"]
units = ["ppb","ppb","ppb"]

for vname,unit in zip(vnames,units):
    fname = "TOMCAT_RUN741_testfile.nc"## input file
    ## check if files exist and size is greater than 0

    if os.path.exists(fname) and os.path.getsize(fname) > 0:
        xdata =read_and_interpolate(fname,mlsp, vname,unit)
        ## write in netcdf format
        xdata.to_netcdf(vname+"_test_plev.nc")

# if you want to have quick look
#xdata1 = xdata[vname]
#xdata1[0,:,:,1].plot(yscale='log',ylim=(300,0.1))
