import numpy as np
import xarray as xr

season = 'Summer'
region = 4
path2somIN = '/home/swenson/ResearchUpdates/20220112_Composites4matlab/'+season+'/Region'+str(region)+'/CenteredByGrid_HovMoller_Data_Region'+str(region)+'_GausianAverage3x3_lbc.nc'
path2somOUT = '/home/swenson/ResearchUpdates/20220421_MatlabSOM_2Dlattice/PlusMinus_6_TimeSteps/P_lt_5/'

# somINname = 'Summer_SOMpats_Region4_GausianAverage3x3_NodeShape_2x2_Kis4.nc'
# Summer_SOMpats_Region5_GausianAverage3x3_Kis4.nc
# somINname = 'Winter_SOMpats_Region4_GausianAverage3x3_NodeShape_3x2_noTC_Kis6.nc'
# Winter_SOMpats_Region4_GausianAverage3x3_Kis6.nc

# Make the hovmoller plots

# Make the horizontal composites