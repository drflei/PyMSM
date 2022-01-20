import sys, os

currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
#sys.path.append(parentdir)
sys.path.insert(0,parentdir)

import pymsm as pm
import plottools as pt

import numpy as np

def global_map(alti=36000.,date='2015-01-01T12:00:00', kp=2):
    '''

    '''
    res = 10  # grid size
    nlon = int(360/res) + 1
    nlat = int(180/res) + 1
    ndata = nlon*nlat
    print(nlat, nlon)
#
    kps = np.empty(ndata,dtype=int)
    kps.fill(kp)

    times = np.empty(ndata,dtype='object')   
    times.fill(date)
#      
    coords = []
    xi = []
    yi = []
    for lat, lon in [(lat, lon) for lat in np.linspace(-87,87,nlat) for lon in np.linspace(0,360,nlon)]:
        if lat < -88: lat = -88
        if lat > 88: lat = 88
        # cannot deal with the extreme latitudes!  
        
        coords.append([alti,lat,lon])  # [alti, lati, longi]
        yi.append(lat)
        xi.append(lon)
                        
#  
    msm = pm.PyMSM(times=times,positions=coords,kps=kps)
    lm, bm, mlats, rcv, es, tf = msm.getTransmissionFunctions()

    #plot the map
    pt.plotmap_contour(np.array(xi),np.array(yi),rcv)

if __name__ == '__main__':
    global_map() 