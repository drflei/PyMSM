# IRBEM version of PyMSM

## Installation:

Easy one step installation using `pip`:

`pip install --upgrade https://github.com/drflei/PyMSM/tarball/irbem`


Note: This version of PyMSM uses the IRBEM package. See [https://github.com/PRBEM/IRBEM/tree/main/python](https://github.com/PRBEM/IRBEM/tree/main/python)

What's used in the installation is actually a forked version of IRBEM which has the TS89 model extended to Kp> 6. 
[https://github.com/drflei/IRBEM/tree/ts89c-ext/python](https://github.com/drflei/IRBEM/tree/ts89c-ext/python)

## Usage:
PyMSM should be imported and instanciated from another python code, e.g.:

```
...    

import pymsm as pm

#   
#     t = ['2002-02-02T12:00:00','2012-02-02T12:00:00','2019-02-02T12:00:00']
#       in ISO string format
#     y = [[1000,45,0],[1000,-45,90],[500,45,-100]] 
#       in GDZ coordinates
#     kp = [0,1,2]
#     rc = [0., 0.2, 1, 10., 30.]
#       in GV units
#     pmsm = pm.PyMSM(t,y, kp,rc)
...

	pmsm = pm.PyMSM(times=t,positions=y, kps=kps)

#   
	lm, bm, mlats, rcv, es, tf = pmsm.getTransmissionFunctions()

#        Return all relevant results for the specified (times, locations) series:
#        Lm: the McIlwain's L-parameter, in np.array
#        Bm: the magnetic field intensity at the location, in np.array
#        Mlat: the magnetic latitude, in np.array
#        ES: the Earth's shadowing factor, in np.array
#        TF: the transmission function, in 2D np.array [len(times) x len(rc)]. The default rc is of the size 34.  

```

## Tests:

There is a build-in test folder, and the main test is the testpymsm.py which can executed directly: 

` $> python testpymsm.py ` 

testpymsm.py uses the plotting utilities from plotools.py to display the various results. 

The other script testglobalmap.py produces the vertical cut-off rigidity map at 36000km altitude: 

` $> python testglobalmap.py ` 
