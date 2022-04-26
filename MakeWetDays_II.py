import numpy as np
import xarray as xr
from tqdm import tqdm

pcip_path = '/home/swenson/era5_data_direct_1.0degree_Version3/tp/'
flist = []
for year in tqdm(range(1980,2003)):
    pcip = xr.open_dataset(pcip_path+'ERA5_'+str(year)+'_tp.nc')

    pcip1 = pcip.resample(time='3H').sum()

    WetDays = np.where(pcip1['tp']>0.001/8,pcip1['tp'],np.nan)

    pcip1['WetDays_lt_125e-6'] = (['time','latitude','longitude'], WetDays)

    pcip1.to_netcdf(pcip_path+'ERA5_'+str(year)+'_tp_3H.nc')

print("DONE")
