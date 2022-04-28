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

OMEGA = 7.2921159e-5  #rotation rate of the earth in radians / second
f = 2 * OMEGA * np.sin(np.radians(U.latitude.values))

_,PlanetaryVort,_ = np.meshgrid(np.arange(len(U.time)),f,U.longitude.values,indexing = 'ij')# * units('s-1')

TotVortAdv = mcf.RelativeVorticityAdvection(RelVort['vo']+PlanetaryVort,U['u'],V['v'])

TotVortAdv_mask = np.where(TotVortAdv>=0,TotVortAdv,np.nan)

# run nanpercenitle on the masked RelVortAdv array
PER = np.zeros((101,TotVortAdv.shape[-2],TotVortAdv.shape[-1]))
for p in tqdm(range(101)):
    PER[p,:,:] = np.nanpercentile(TotVortAdv_mask,q=p,axis=0)

# make a time,lat,lon array of the positive RelVortAdv percentile values.
THEmetric = np.zeros(TotVortAdv.shape)
for p in tqdm(range(1,101)):
    THEmetric = np.where(TotVortAdv.values>np.broadcast_to(PER[p,:,:],TotVortAdv.shape),p,THEmetric)

# THEmetric_within_range = np.zeros(THEmetric.shape)
# for tt in range(THEmetric_within_range.shape[0]):
#     if tt == 0:
#         t1=0
#         t2=1
#     elif tt == THEmetric_within_range.shape[0]-1:
#         t1=-2
#         t2=-1
#     else:
#         t1=tt-1
#         t2=tt+2
#     for jj in range(THEmetric_within_range.shape[1]):
#         if jj == 0:
#             j1=0
#             j2=1
#         elif jj == THEmetric_within_range.shape[1]-1:
#             j1=-2
#             j2=-1
#         else:
#             j1=jj-1
#             j2=jj+2
#         for ii in range(THEmetric_within_range.shape[2]):
#             if ii == 0:
#                 i1=0
#                 i2=1
#             elif ii == THEmetric_within_range.shape[2]-1:
#                 i1=-2
#                 i2=-1
#             else:
#                 i1=ii-1
#                 i2=ii+2
            
#             if 
#             within_range = THEmetric
        

coords = dict(

    time = (['time'],U['time']),
    latitude = (['latitude'],U['latitude']),
    longitude = (['longitude'],U['longitude'])

)

data_vars=dict(

    PositiveTotalVorticityAdvection_Percentile = (['time','latitude','longitude'], THEmetric)

    )

plot_var_out = xr.Dataset(

    data_vars=data_vars,
    coords=coords,

    )

plot_var_out.to_netcdf('/home/swenson/projects/PEX_feature_metrics/PositiveTotalVorticityAdvection_Percentile.nc',mode='w')

print('DONE')