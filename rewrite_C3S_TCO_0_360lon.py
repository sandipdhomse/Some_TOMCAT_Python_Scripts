'''
C3S data file longitudes range from -180 to 180
this code combines C3S total column data files
and reorganises them as 0 to 360
Author: Sandip Dhomse (sandip.dhomse@gmail.com)
'''
# coding: utf-8
#!/usr/bin/python
import netCDF4 as nc
import numpy as np
from scipy import interpolate
import datetime
import os


# Function to create netcdf files
def writenc( ncfile,mlons,mlats,mtimes,tracer,tracername, expname, runtitle ):
    dataset = nc.Dataset(ncfile, 'w', format='NETCDF3_CLASSIC')
    nlon = len(mlons)
    nlat = len(mlats)
    ntime = len(mtimes)

    lon=dataset.createDimension('lon', nlon)
    lat=dataset.createDimension('lat', nlat)
    time=dataset.createDimension('time', None)

    times=dataset.createVariable('time',np.float64,('time',))
    lons=dataset.createVariable('lon',np.float32,('lon',))
    lats=dataset.createVariable('lat',np.float32,('lat',))

    temp1=dataset.createVariable(tracername,np.float32,('time','lat','lon',))

    time_num  = datetime.datetime.now().strftime("%I:%M%p on %B %d, %Y")
    
    ##  Add global attributes

    dataset.branch_Time ='November 2015'
    dataset.contact = 'Sandip Dhomse and Martyn Chipperfield (s.s.dhomse@leeds.ac.uk)'
    dataset.Convensions ='CF-1.6'
    dataset.creation_date = str(time_num)
    dataset.frequency = 'mon'
    dataset.frequency = 'mon'
    dataset.institution = "University of Leeds, UK"
    dataset.product ="total ozone"
    import uuid
    fid= uuid.uuid4()
    dataset.tracking_id =str(fid)

    lats.units = 'degrees_north'
    lons.units = 'degrees_east'
    times.units = "months since 1979-01-15 00:00:00"
    lons[:] = mlons[:]
    lats[:] = mlats[:]

    times[:] = mtimes[:]
    temp1[:,:,:] = (tracer[:,:,:])
    dataset.close()
    return

nyrs=39
nmons=nyrs*12
bdata=np.zeros((nmons,361,720))
## Give path for input files
inpath  ="/nfs/../.../../"
tt=0
for year in range(1979,1979+nyrs):
    for month in range(1,13):
        mm=str(month)
        yy=str(year)
        if month<10:
            mm='0'+mm
        print(yy+mm)
        vn1 = nc.Dataset(inpath+"C3S_OZONE-L4-TC-ASSIM_MSR-"+yy+mm+"-fv0020.nc")

        mlons = vn1.variables["longitude"][:]
        mlats = vn1.variables["latitude"][:]
        o3 = vn1.variables["total_ozone_column"][:]
        times = vn1.variables["time"][:]

        o3=np.squeeze(o3)
        nlon  =len(mlons)
        nlat  =len(mlats)
        bdata[tt,:,0:361]=o3[:,359:721]
        bdata[tt,:,361:721]=o3[:,0:359]
        tt=tt+1


bdata[bdata<80]=np.nan
mlevs =[1000.]
tname="toz"

ofile = 'TCO_C3S_0_360lon.nc'
tt=np.arange(0,nmons,1)
print(ofile)
times=np.arange(0,nmons,1.)
newlons=np.arange(0,360,.5)
writenc( ofile,newlons,mlats,times,bdata,tname, "C3S","C3S" )
