import MetricCalculatorFunctions as mcf
import numpy as np
import xarray as xr
import metpy.calc as mpcalc
import metpy.constants as mpconstants
import metpy
from metpy.units import units
from tqdm import tqdm

hours = 3

level = 500

RelVort = xr.open_mfdataset('/home/swenson/era5_data_direct_1.0degree_Version3/Vorticity_'+str(level)+'hPa/*.nc',concat_dim='time',combine='by_coords')
RelVort = RelVort.resample(time=str(hours)+'H').nearest()
U = xr.open_mfdataset('/home/swenson/era5_data_direct_1.0degree_Version3/U_'+str(level)+'hPa/*.nc',concat_dim='time',combine='by_coords')
U = U.resample(time=str(hours)+'H').nearest()
V = xr.open_mfdataset('/home/swenson/era5_data_direct_1.0degree_Version3/V_'+str(level)+'hPa/*.nc',concat_dim='time',combine='by_coords')
V = V.resample(time=str(hours)+'H').nearest()

RelVortAdv = mcf.RelativeVorticityAdvection(RelVort['vo'],U['u'],V['v'])

RelVortAdv_mask = np.where(RelVortAdv>=0,RelVortAdv,np.nan)

# run nanpercenitle on the masked RelVortAdv array
PER = np.zeros((101,RelVortAdv.shape[-2],RelVortAdv.shape[-1]))
for p in tqdm(range(101)):
    PER[p,:,:] = np.nanpercentile(RelVortAdv_mask,q=p,axis=0)

# make a time,lat,lon array of the positive RelVortAdv percentile values.
THEmetric = np.zeros(RelVortAdv.shape)
for p in tqdm(range(1,101)):
    THEmetric = np.where(RelVortAdv.values>np.broadcast_to(PER[p,:,:],RelVortAdv.shape),p,THEmetric)

coords = dict(

    time = (['time'],U['time']),
    latitude = (['latitude'],U['latitude']),
    longitude = (['longitude'],U['longitude'])

)

data_vars=dict(

    PositiveRelativeVorticityAdvection_Percentile = (['time','latitude','longitude'], THEmetric)

    )

plot_var_out = xr.Dataset(

    data_vars=data_vars,
    coords=coords,

    )

plot_var_out.to_netcdf('/home/swenson/projects/PEX_feature_metrics/PositiveRelativeVorticityAdvection_Percentile.nc',mode='w')

print('DONE')