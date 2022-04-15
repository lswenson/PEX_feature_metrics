import numpy as np
import xarray as xr
import metpy.calc as mpcalc
import metpy.constants as mpconstants
import metpy
from metpy.units import units

# relative_humidity_from_specific_humidity(pressure, temperature, specific_humidity)
# dewpoint_from_relative_humidity(temperature, relative_humidity)
# equivalent_potential_temperature(pressure, temperature, dewpoint)

# mpcalc.gradient(th,axes=[0],deltas=[1*60*60])[0] *units('s-1')

# mpcalc.first_derivative(u*th, axis=-1, delta=[wrf_xr.attrs['DX']])

# mpcalc.advection(th, u=u, dx=wrf_xr.attrs['DX']*units('m'), x_dim=- 1)
# mpcalc.advection(th, v=v, dy=wrf_xr.attrs['DY']*units('m'), y_dim=- 2)
# mpcalc.advection(th, w=w, dz=np.diff(pres,axis=1)*units('Pa'), vertical_dim=- 3)

def EquivalentPotentialTemperatureGradient(sph,T,p):
    
    dewpoint = mpcalc.dewpoint_from_specific_humidity(p, T, sph)
    
    EQtheta = mpcalc.equivalent_potential_temperature(p, T, dewpoint)

    DX,DY = mpcalc.lat_lon_grid_deltas(T['longitude'], T['latitude'], x_dim=- 1, y_dim=- 2, geod=None)
    ddx = np.broadcast_to(DX,(EQtheta.shape[0],DX.shape[0],DX.shape[1]))
    ddy = np.broadcast_to(DY,(EQtheta.shape[0],DY.shape[0],DY.shape[1]))

    grad1,grad2 = mpcalc.gradient(EQtheta,axes=[-2,-1],deltas=[ddy,ddx])
    
    return np.power(np.power(grad1,2)+np.power(grad2,2),1/2)

def ThermalFrontalParameter(T,lat,lon):
    
    dTdy,dTdx = mpcalc.gradient(T, axes=[1,2])
    
    magGradT = np.power(np.power(dTdy,2)+np.power(dTdx,2),1/2)

    DX,DY = mpcalc.lat_lon_grid_deltas(T['longitude'], T['latitude'], x_dim=- 1, y_dim=- 2, geod=None)
    ddx = np.broadcast_to(DX,(T.shape[0],DX.shape[0],DX.shape[1]))
    ddy = np.broadcast_to(DY,(T.shape[0],DY.shape[0],DY.shape[1]))
    
    dMGTdy,dMGTdx = mpcalc.gradient(magGradT, axes=[1,2], deltas=[ddy,ddx])
    
    Yside = np.multiply(dMGTdy,np.divide(dTdy,magGradT))
    Xside = np.multiply(dMGTdx,np.divide(dTdx,magGradT))
    
    return -1 * (Yside + Xside)

def dCAPEdt(CAPE,dt):

    ddt = mpcalc.gradient(CAPE,axes=[0],deltas=[dt])

    return ddt

def RelativeVorticityAdvection(vorticity,u,v):

    DX,DY = mpcalc.lat_lon_grid_deltas(u['longitude'], u['latitude'], x_dim=- 1, y_dim=- 2, geod=None)
    ddx = np.broadcast_to(DX,(u.shape[0],DX.shape[0],DX.shape[1]))
    ddy = np.broadcast_to(DY,(v.shape[0],DY.shape[0],DY.shape[1]))

    VA = mpcalc.advection(vorticity, u=u, v=v, x_dim=-1, y_dim=-2, dx=ddx, dy=ddy)

    return VA

def FindCutOffLow(GHT200,GHT300,U,TFP):
    # MAYBE WE DONT NEED TO USE THIS

    # step zero
    IDarray = np.zeros(GHT200.shape)
    # step one
    for tt in np.arange(0,GHT200.shape[0]):
        for jj in np.arange(1,GHT200.shape[-2]):
            for ii in np.arange(1,GHT200.shape[-1]):
                IsLower = GHT200[tt,jj-1:jj+2,ii-1:ii+2]
                IsLower = np.where(IsLower.flatten()<IsLower[1,1])
                if len(IsLower) >= 6:
                    IDarray[tt,jj,ii] = 1

    # step two
    ID_s1 = np.where(IDarray==1)
    for lm in range(len(ID_s1[0])):
        tt = ID_s1[0][lm]
        jj = ID_s1[1][lm]
        ii = ID_s1[2][lm]
        thk = GHT200[tt,jj,ii] - GHT300[tt,jj,ii]
        thk_east = GHT200[tt,jj,ii+1] - GHT300[tt,jj,ii+1]
        if thk > thk_east:
            IDarray[tt,jj,ii] = 0
    # step three
        elif TFP[tt,jj,ii] > TFP[tt,jj,ii+1]:
            IDarray[tt,jj,ii] = 0

    return IDarray

def MakeUpslopeFlow(U,V,TOPO):

    DX,DY = mpcalc.lat_lon_grid_deltas(U['longitude'], U['latitude'], x_dim=- 1, y_dim=- 2, geod=None)
    UdotTopo = mpcalc.advection(TOPO, u=U, v=V, dx=DX, dy=DY, x_dim=- 1, y_dim=-2)
    
    return UdotTopo

def great_circle_array(clon,clat,XLON,XLAT):
    """
    clon and clat are the target point (in degrees) to be compared to every point in XLON and XLAT
    XLON and XLAT should be the result of meshgrid meaning they are both 2D (in degrees)
    Out put is a 2D array of Great Circle Distances (GCD) from each point to the target point in meters
    """
    rearth = 6371000 #meters
    clon,clat = map(radians,[clon,clat])
    XLON,XLAT = map(np.radians,[XLON,XLAT])
    
    return rearth * ( np.arccos(sin(clat) * np.sin(XLAT) + np.cos(clat) * np.cos(XLAT) * np.cos(clon - XLON)) )

def GCD_m2deg(GCD):
    rearth = 6371000 #meters
    return (GCD / (2*np.pi*rearth)) * 360

def GCD_deg2m(GCD):
    rearth = 6371000 #meters
    return (GCD / 360) * (2*np.pi*rearth)