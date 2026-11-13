# IRBEM version of PyMSM
## Introduction
Please see the technical report for a description of the model and all the technical details
[Technical Report](https://github.com/drflei/PyMSM/blob/irbem/MSM%20Technical%20Report.pdf)
## Installation:

### Pre-requisite ###
This version of PyMSM uses the IRBEM package and it must have been installed separately already. See [https://github.com/PRBEM/IRBEM/tree/main/python](https://github.com/PRBEM/IRBEM/tree/main/python).

### Quick Installation (Recommended)

Install directly from GitHub using pip:

```bash
pip install git+https://github.com/drflei/PyMSM
```

### Development Installation

For development or to get the latest code:

```bash
git clone --recurse-submodules https://github.com/drflei/PyMSM.git
cd PyMSM
./scripts/update_maps.sh
pip install -r requirements.txt -e ./
```

## Usage:
PyMSM should be imported and instanciated from another python code, e.g.:

PyMSM uses IRBEM's Kp-based Dst estimate for the T89c Boberg model; it does not
accept measured Dst as an input.
Map values are linearly interpolated in time between five-year epochs from 1955
through 2025; dates outside that range use the nearest endpoint map.

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

## Tests

Run these commands from the repository root, with PyMSM installed, the `MAPS`
submodule initialized (`git submodule update --init -- MAPS`), and an IRBEM
library that supports `kext=15` installed. The automated tests also require
`pytest`, which is not included in `requirements.txt`.

### Automated tests

```bash
python -m pytest -q test/test_boberg_missing_dst.py \
  test/test_epoch_interpolation.py test/test_kp_cutoff.py
```

For a headless environment, prefix the command with `MPLBACKEND=Agg`.

- `test/test_boberg_missing_dst.py`: checks low and high Kp, mixed-Kp
  inputs, and the Boberg model's Kp-based estimate when Dst is missing.
- `test/test_epoch_interpolation.py`: checks map-epoch interpolation at
  endpoints, between epochs, across leap years, and outside the map range.
- `test/test_kp_cutoff.py`: calculates and plots vertical cutoff rigidity for
  Kp 0–9 at 50° N, 0° E, 400 km on 2020-01-01 00:00 UT; it also checks
  individual versus batch calculations and PNG output.

### Plotting scripts

```bash
python test/testpymsm.py
python test/testglobalmap.py
python test/test_kp_cutoff.py [output.png]
```

- `testpymsm.py` plots sample maps and an orbit trajectory's cutoff, magnetic
  parameters, shadowing factor, and transmission functions. Its Basemap plot
  requires the separately installed `basemap` package.
- `testglobalmap.py` plots a global vertical cutoff map at 36,000 km for
  2015-01-01 12:00 UT and Kp 2.
- `test_kp_cutoff.py` prints the Kp sweep and saves its plot to the specified
  path, or to `kp_cutoff_2020.png` in the current directory by default.
