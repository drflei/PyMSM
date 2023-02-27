# IRBEM version of PyMSM

## Installation:

Easy one step installation using `pip`:

~~`pip install --upgrade https://github.com/drflei/PyMSM/tarball/irbem`~~

The above is not working properly. Has to be done in mutiple steps:
```
$ git clone --recurse-submodules https://github.com/drflei/PyMSM.git
$ cd PyMSM
$ pip install -r requirements.txt -e ./
```

==Note 1:== This version of PyMSM uses the IRBEM package and it must have been installed separately. See [https://github.com/PRBEM/IRBEM/tree/main/python](https://github.com/PRBEM/IRBEM/tree/main/python). However, a forked version of IRBEM which has the TS89 model extended to Kp> 6 is best to be used together with PyMSM [https://github.com/drflei/IRBEM/tree/ts89c-ext/python](https://github.com/drflei/IRBEM/tree/ts89c-ext/python) 
and it can be install following these steps:

```
$  git clone -b ts89c-ext --single-branch https://github.com/drflei/IRBEM.git
$  cd IRBEM/
$  make OS=linux64 ENV=gfortran64 all
$  make install
$  cd python/
$  pip install -r requirements.txt
```
==Note 2:==
The above installation is not working in a Ubuntu20.04 docker container! The reason is still unknow. One work-around solution is to install it in a virtual env instead:

```
cd /opt/ # or any tmp folder
git clone --recurse-submodules https://github.com/drflei/PyMSM.git
cd /opt
git clone -b ts89c-ext --single-branch https://github.com/drflei/IRBEM.git
cd IRBEM/
make OS=linux64 ENV=gfortran64 all
make install
#
apt update
apt install python3-virtualenv
cd /opt
virtualenv -p python3 msmtest
# msmtest is the virtual enve
source msmtest/bin/activate
cd PyMSM
pip install -r requirements.txt -e ./
cd ../IRBEM/python/
pip install -r requirements.txt
cd /opt
# Installation completed and run the tests
# you may needed to set your DISPLAY pointing to you x11 server before run below
python3 /opt/PyMSM/test/testglobalmap.py
wget http://distfiles.macports.org/py-matplotlib-basemap/basemap-1.2.2.tar.gz
python3.8 -m pip install basemap-1.2.2.tar.gz
python3 /opt/PyMSM/test/testpymsm.py
```

or

```
cd /opt/ # or any tmp folder
git clone --recurse-submodules https://github.com/drflei/PyMSM.git
cd PyMSM
source installondocker.sh
```

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
