import sys, os

currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
#sys.path.append(parentdir)
sys.path.insert(0,parentdir)

import pymsm as pm
import plottools as pt

import numpy as np

def global_map(alti=500.,date='2015-01-01T12:00:00', kp=2):
    '''

    '''
    res = 15  # grid size
    nlon = int(360/res) + 1
    nlat = int(180/res) + 1
    ndata = nlon*nlat

#
    kps = np.empty(ndata,dtype=int)
    kps.fill(kp)

    times = np.empty(ndata,dtype='object')   
    times.fill(date)
#      
    coords = []
    xi = []
    yi = []
    for j in range(90,-95,-res): 
        for i in range(0,365,res):
            coords.append([alti,i,j])  # [alti, lati, longi]
            yi.append(j)
            xi.append(i)
    xi = np.array(xi).reshape(nlat,nlon).transpose()
    yi = np.array(yi).reshape(nlat,nlon).transpose()
#  
    msm = pm.PyMSM(times,np.array(coords),kps)
    rc = msm.getRc()
    zi = np.array(rc).reshape(nlat,nlon).transpose()

    #plot the map
    pt.plotmap_contour(xi,yi,zi)

if __name__ == '__main__':
    global_map() 