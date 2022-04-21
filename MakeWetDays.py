import numpy as np
import xarray as xr

pcip_path = '/home/swenson/era5_data_direct_1.0degree_Version3/tp/'
flist = []
for year in range(1980,2011):
    flist.append(pcip_path+'ERA5_'+str(year)+'_tp.nc')
pcip_dat = xr.open_mfdataset(flist,concat_dim='time')
print("data loaded")

pcip = pcip_dat.resample(time='3H').sum()
print("time resampled")

WetDays = np.where(pcip['tp']>0.001/8,pcip['tp'],np.nan)
print("wet days made")

pcip['WetDays_lt_125e-6'] = (['time','latitude','longitude'], WetDays)

pcip.to_netcdf(pcip_path+'ERA5_1980-2010_tp_3H.nc')

print("DONE")
