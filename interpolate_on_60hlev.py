
# coding: utf-8
#!/usr/bin/python


'''
A Code to interpolate TOMCAT tracers on altitude levels
author: Sandip Dhomse (sandip.dhomse@gmail.com)
'''

import netCDF4 as nc
import numpy as np
import xarray as xr
import pandas as pd
from datetime import datetime 
import os

# Convert numpy array to Xarray



# Convert numpy array to Xarray
def make_xarray(odata,levs,lats,lons,ntimes,vname,sdate,unit1):
    
    time1 = pd.date_range(sdate, periods=ntimes, freq="1M")+ pd.DateOffset(days=-16)
    ds = xr.Dataset({vname: (['time', 'height', 'lat', 'lon'], odata.astype('float32'), {'units':unit1})},
                 coords={'lon': (['lon'], lons, {'units':'degrees_east'}),
                         'lat': (['lat'], lats, {'units':'degrees_north'}),
                         'height': (['height'], levs, {'units':'km'},{'positive':'up'}),
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



## sometimes model grid point height are greater than 1km or less than 60km

def revise_ax(nlev,ax,ay):
    nax = np.zeros(nlev+2)
    nay = np.zeros(nlev+2)
    nax[0]=70## lets top height of the model atmosphere
    nax[1:nlev+1] = ax[:]
    nax[nlev+1] = -1## set bottom height 
    
    nay[0]= 0.8*ay[0]
    nay[1:nlev+1] = ay[:]
    nay[nlev+1] = 0.001*ay[nlev-1]
    return nax[::-1],nay[::-1]



def read_and_interpolate(afile,newht,vname,unit):
    vn1 = nc.MFDataset(afile)
    #print(vn1)
    lons = vn1.variables["lon"][:]
    lats = vn1.variables["lat"][:]
    tracer = vn1.variables[vname+"_mm"][:]
    gg = vn1.variables["g_mm"][:]
    time = vn1.variables["time"][:]
    #clo =0
    nlon  =len(lons)
    nlat  =len(lats)
    ntimes = np.size(time)
    #ntimes =4
    nlev = np.size(tracer[0,:,0,0])
    
    adata = np.zeros((ntimes,np.size(newht),nlat,nlon))
    ht =gg/980## convert geopotential to height
    for itime in range(ntimes):
        for ilat in range(0,nlat):
            for ilon in range(0,nlon):
                ay = tracer[itime,:,ilat,ilon]
                ax = ht[itime,:,ilat,ilon]
                
                nax,nay = revise_ax(nlev,ax,ay)
                newy =np.interp(newht,nax,nay)
                adata[itime,:,ilat,ilon] = newy*1.e9## convert ppb 
               
    
    ds = xr.open_dataset(afile)
    sdate = str(ds.time.values[0])[0:10]
    ## create a xarray
    ds = make_xarray(adata,newht, lats, lons, ntimes, vname, sdate,unit)
    return ds[vname]


###new height levels
newht = np.arange(1,61)
## test file has only one tracer ## normal file would have many more..
## give their names below



vnames = ["O3","HCL","HNO3"]
units = ["ppb","ppb","ppb"]

for vname,unit in zip(vnames,units):
    fname = "TOMCAT_RUN741_testfile.nc"
    
    if os.path.exists(fname) and os.path.getsize(fname) > 0:
        xdata =read_and_interpolate(fname,newht, vname,unit)
        ## write in netcdf format
        xdata.to_netcdf(vname+"_test_ht.nc")

# if you want to have quick look
#xdata1= xdata[vname]
#xdata[0,:,:,0].plot(ylim=(10,60),cmap = 'Reds')
