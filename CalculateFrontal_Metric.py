import MetricCalculatorFunctions as mcf
import numpy as np
import xarray as xr
import metpy.calc as mpcalc
import metpy.constants as mpconstants
import metpy
from metpy.units import units
from tqdm import tqdm

y1 = 1979
y2 = 2010

hours=3

level = 850

T = xr.open_mfdataset('/home/swenson/era5_data_direct_1.0degree_Version3/Temperature_'+str(level)+'hPa/*.nc',concat_dim='time',combine='by_coords')
T = T.resample(time=str(hours)+"H").nearest()
SPH = xr.open_mfdataset('/home/swenson/era5_data_direct_1.0degree_Version3/specific_humidity_'+str(level)+'hPa/*.nc',concat_dim='time',combine='by_coords')
SPH = SPH.resample(time=str(hours)+"H").nearest()

gradEPT = mcf.EquivalentPotentialTemperatureGradient(SPH['q'],T['t'],level * units('hPa'))

gradEPT_per = np.zeros((101,T['t'].shape[-2],T['t'].shape[-1]))
for p in tqdm(range(101)):
    gradEPT_per[p,:,:] = np.nanpercentile(gradEPT,q=p,axis=0)

# save gradEPT_per

# mask every gradEPT value below the 50th percentile
gradEPT_mask = np.where(gradEPT.magnitude<np.broadcast_to(gradEPT_per[50,:,:],gradEPT.shape),np.nan,gradEPT)

# run nanpercentile loop on the masked array to create the scores that translate to the ranks of the actual metric
gradEPT_mask_per = np.zeros((101,T['t'].shape[-2],T['t'].shape[-1]))
for p in tqdm(range(101)):
    gradEPT_mask_per[p,:,:] = np.nanpercentile(gradEPT_mask,q=p,axis=0)

# create a time,lat,lon array with the metric scores 
THEmetric = np.zeros(gradEPT.shape)
for p in tqdm(range(1,101)):
    THEmetric = np.where(gradEPT_mask.magnitude>np.broadcast_to(gradEPT_mask_per[p,:,:],gradEPT.shape),gradEPT_mask_per[p,:,:],THEmetric)

# create a time,lat,lon array with the raw percentile scores.
RAWper = np.zeros(gradEPT.shape)
for p in tqdm(range(1,101)):
    RAWper = np.where(gradEPT.magnitude>np.broadcast_to(gradEPT_per[p,:,:],gradEPT.shape),gradEPT_per[p,:,:],RAWper)

coords = dict(

    time = (['time'],T['time']),
    latitude = (['latitude'],T['latitude']),
    longitude = (['longitude'],T['longitude'])

)

data_vars=dict(

    gradEPT_percentile_AboveMedian = (['time','latitude','longitude'], THEmetric),
    gradEPT_percentile = (['time','latitude','longitude'], RAWper)

    )

plot_var_out = xr.Dataset(

    data_vars=data_vars,
    coords=coords,

    )

plot_var_out.to_netcdf('/home/swenson/projects/PEX_feature_metrics/gradEPT_Frontal_percentile.nc',mode='w')

print('DONE')