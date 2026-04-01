"""Test script for generating global rigidity cutoff maps."""
import os
import sys

currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

import numpy as np
import pymsm as pm
import plottools as pt


def global_map(alti=36000., date='2015-01-01T12:00:00', kp=2):
    """

    Generate a global vertical cutoff rigidity map.
    
    Args:
        alti: Altitude in km (default: 36000).
        date: Date-time in ISO format (default: '2015-01-01T12:00:00').
        kp: Kp index (default: 2).
    """
    res = 10  # Grid size
    nlon = int(360 / res) + 1
    nlat = int(180 / res) + 1
    ndata = nlon * nlat
    print(nlat, nlon)
    
    kps = np.empty(ndata, dtype=int)
    kps.fill(kp)

    times = np.empty(ndata, dtype='object')
    times.fill(date)
    
    coords = []
    xi = []
    yi = []
    for lat, lon in [(lat, lon) for lat in np.linspace(-87, 87, nlat) for lon in np.linspace(0, 360, nlon)]:
        if lat < -88:
            lat = -88
        if lat > 88:
            lat = 88
        # Cannot deal with extreme latitudes
        coords.append([alti, lat, lon])  # [alti, lati, longi]
        yi.append(lat)
        xi.append(lon)
    
    msm = pm.PyMSM(times=times, positions=coords, kps=kps)
    lm, bm, mlats, rcv, es, tf = msm.getTransmissionFunctions()

    # Plot the map
    pt.plotmap_contour(np.array(xi), np.array(yi), rcv)


if __name__ == '__main__':
    global_map() 