import spacepy.time as spt
import spacepy.coordinates as spc
import spacepy.irbempy as ib
import spacepy.omni as om
import datetime

import numpy as np

import pymsm as pm
import plottools as pt


def testpymsm():
    '''

    '''    
    # map size
    nlon = 73
    nlat = 37
    ndata = nlon*nlat
    # make up the grid for the plotting tests
    xi = []
    yi = []
    for j in range(90,-95,-5):   # Y grid
        for i in range(0, 365, 5):    # X grid
            yi.append(j)
            xi.append(i)
    xi = np.array(xi).reshape(nlat,nlon)
    yi = np.array(yi).reshape(nlat,nlon)
    xi = xi.transpose()
    yi = yi.transpose()

    mapMgr = pm.MapDB()
    
    # test plotting: plot a pre-calculated map
    lm, rc, rclm2 = mapMgr.getMap('2025','1','12')
    pt.plotscatter(lm,rc, xtit='Lm', ytit='Rc (GV)')
    pt.plotmap(xi,yi,rc)
    pt.plotmap_contour(xi,yi,rc)
    pt.plotmap_basemap(xi,yi,rc)
    pt.plotmapfile("./MAPS/2000/AVKP9T12.AVG")
        
    print("Plotting test Completed!")
    
    #
    N = 4*5 
    times = np.empty(N,dtype='object')
    kps = np.empty(N,dtype=int)
    mjd = 55000.5
    for i in range(N):
        mjd += 0.00000001*i
        times[i] = mjd  
#    times.fill('2002-02-02T12:00:00.0')
    kps.fill(1)
#       
    coords =[]
    alti = 1000.
    for i in range(-60,60,30): 
        for j in range(0,360,72):
            coords.append([alti, i, j ]) 
    t = spt.Ticktock(times,'MJD')
    y = spc.Coords(coords,'GDZ','sph')
    y.ticks = t
    pmsm = pm.PyMSM(times=t,positions=y, kps=kps)
        
#     t = spt.Ticktock(['2002-02-02T12:00:00','2012-02-02T12:00:00','2019-02-02T12:00:00'])
#     y = spc.Coords([[1000,45,0],[1000,-45,90],[500,45,-100]],'GDZ','sph')
#     y.ticks = t
#     kp = [0,1,2]
#     rc = [0., 0.2, 1, 10., 30.]
#     pm = PyMSM(t,y, kp,rc)

    lm, bm, mlats, rcv, es, tf = pmsm.getTransmissionFunctions()
    
    pt.plotscatter(times,lm, xtit='MJD', ytit='Lm (Re)', title = ' ')
    pt.plotscatter(times,bm, xtit='MJD', ytit='Bm ', title = ' ')
    pt.plotscatter(times,mlats,xtit='MJD', ytit='MLat (rad)', title = ' ')
    pt.plotscatter(times,rcv, xtit='MJD', ytit='Rc (Gv)', title = ' ')
    pt.plotscatter(times,es, xtit='MJD', ytit='ES ()', title = ' ')
    
    colname = [0.1, 0.2, 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, \
        6., 6.5, 7., 7.5, 8., 9., 10., 11., 12., 13, 14., 15., 16., 17., \
        20., 25., 20., 30., 40., 50., 60.]
    cn = []
    for i in range(len(colname)):
        if i%4 == 0: cn.append(colname[i]) 
        
    rowname = times
    rn = []
    for i in range(len(rowname)):
        if i%2 == 0: rn.append(rowname[i])  
    
    pt.plot3D(tf, xgrid = colname,xtit='Rigidity (GV)', ytit=' Time (MJD)',ztit='Transmission', title='Transmission Functions' )
    pt.plotbar3d(tf, colname = cn, rowname= rn, xtit='Rigidity (GV)', ytit=' Time (MJD)',ztit='Transmission', title='Transmission Functions')
    
    print ('Test completed')
    

if __name__ == '__main__':
    testpymsm() 
