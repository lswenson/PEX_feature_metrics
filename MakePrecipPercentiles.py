import numpy as np
import xarray as xr
from tqdm import tqdm

pcip_path = '/home/swenson/era5_data_direct_1.0degree_Version3/tp/'
flist = []
for year in tqdm(range(1980,2011)):
    flist.append(pcip_path+'ERA5_'+str(year)+'_tp_3H.nc')
pcip = xr.open_mfdataset(flist,concat_dim='time')

WD = pcip['WetDays_lt_125e-6']
WD_per = np.zeros((100,WD.shape[1],WD.shape[2]))
for per in tqdm(range(1,100)):
    WD_per[per,:,:] = np.nanpercentile(WD,q=per,axis=0)

tp = pcip['tp']
tp_per = np.zeros((100,tp.shape[1],tp.shape[2]))
for per in tqdm(range(1,100)):
    tp_per[per,:,:] = np.nanpercentile(tp,q=per,axis=0)

pcip.assign_coords(coords={'percentile':np.arange(0,100)})

pcip['WD_percentile'] = (['percentile','latitude','longitude'], WD_per)
pcip['tp_percentile'] = (['percentile','latitude','longitude'], tp_per)

pcip.to_netcdf(pcip_path+'ERA5_1980-2010_tp_3H_percentiles.nc')

print("DONE")