import sys, os

currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
#sys.path.append(parentdir)
sys.path.insert(0,parentdir)

import src.pymsm as ms
import numpy as np
import datetime

def test1():
    '''

    '''
    nlon = 73
    nlat = 37
    ndata = nlon*nlat
    res = 5  # grid size
#
    kps = np.empty(ndata,dtype=int)
    kps.fill(2)
    # rc = [0.1, 0.2, 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, \
    #   6., 6.5, 7., 7.5, 8., 9., 10., 11., 12., 13, 14., 15., 16., 17., \
    #   20., 25., 20., 30., 40., 50., 60.] 
#
    times = np.empty(ndata,dtype='object')   
    #dat = datetime.datetime(2012,7,1,12,0,0)
    dat = datetime.datetime.fromisoformat('2015-01-01T12:00:00')
    
    times.fill(dat)
#      
    coord = []
    alti = 500.
    for i in range(90,-91,-res): 
        for j in range(0,361,res):
            coord.append([alti,i,j])  # [alti, lati, longi]

#  
    coord = np.array(coord)
#    print (coord.shape)
#    print(coord[1000,1])
    pm = ms.PyMSM(times,np.array(coord),kps)
    rc = pm.getRc()
    zi = np.array(rc).reshape(nlat,nlon)
    #zi = zi.transpose()
    print(zi.size, zi.shape)
    
    print (zi[18, :] )

    ms.plotmap_c(zi)

#    lm, rc, rclm2 = mapMgr.getMap('2000','9','12')
#    plotscatter(lm,rc)
#    plotmap_c(rc)

if __name__ == '__main__':
    test1() 