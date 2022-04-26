import MetricCalculatorFunctions as mcf
import numpy as np
import xarray as xr
import metpy.calc as mpcalc
import metpy.constants as mpconstants
import metpy
from metpy.units import units

hours = 3
dt = hours*60*60

CAPE = xr.open_mfdataset('/home/swenson/era5_data_direct_1.0degree_Version3/cape/*.nc',concat_dim='time',combine='by_coords')
CAPE = CAPE.resample(time=str(hours)+'H').nearest()

ddt = mcf.dCAPEdt(CAPE,dt)

# mask positive ddt values (CAPE Generation isn't interesting)
NEGddt = np.where(ddt.magnitude<=0,ddt,np.nan)

# run nanpercenitle on the masked ddt array
PERddt = np.zeros((101,CAPE.shape[-2],CAPE.shape[-1]))
for p in range(101):
    PERddt[p,:,:] = np.nanpercentile(-1*NEGddt,q=p,axis=0)

# make a time,lat,lon array of the negative dCAPEdt percentile values.
THEmetric = np.zeros(NEGddt.shape)
for p in range(1,101):
    THEmetric = np.where(-1*NEGddt.magnitude>np.broadcast_to(PERddt[p,:,:],NEGddt.shape),p,THEmetric)

coords = dict(

    time = (['time'],CAPE['time']),
    latitude = (['latitude'],CAPE['latitude']),
    longitude = (['longitude'],CAPE['longitude'])

)

data_vars=dict(

    CAPE_consumption_percentile = (['time','latitude','longitude'], THEmetric)

    )

plot_var_out = xr.Dataset(

    data_vars=data_vars,
    coords=coords,

    )

plot_var_out.to_netcdf('/home/swenson/projects/PEX_feature_metrics/CAPE_consumption_percentile.nc',mode='w')

print('DONE')
