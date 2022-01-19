# Spacepy version of PyMSM

## Usage:

This version of PyMSM requires the spacepy package and it can be installled via pip:

`$> pip install spacepy` 

PyMSM should be imported and instanciated from another python code, e.g.:

```
...    

import PyMSM as pm

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

#     Returns are list correspond to the specified times and positions series:
#        Lm: the McIlwain's L-parameter, in a list
#        Bm: the magnetic field intensity at the mirror point, in a list
#        Mlat: the magnetic latitude, in a list
#        ES: the Earth's shadowing factor, in a list
#        TF: the transmission function, in 2D array [len(times) x len(rc)]. The default rc is of the 
#            size 34. 

```

## Test:

There is a build-in test folder, and the main test is the testpymsm.py which can executed directly: 

` $> python testpymsm.py ` 

testpymsm.py uses the plotting utilities from plotools.py to display the various results. 
 