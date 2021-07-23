# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#     MODULE:           example1.py
#
#     Version:          0.B
#     Date:             16/08/2019
#     Author:           Pete Truscott
#     Organisation:     Kallisto Consultancy Ltd, UK
#     Customer:         ESA/ESTEC, NOORDWIJK
#     Contract:         4000117974/16/NL/LF (VALIRENE Project)
#
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#     LICENCE AND COPYRIGHT/DISTRIBUTION CONDITIONS
#
#     (To be confirmed with Hugh Evans (ESA) and Daniel Heynderickx (DHC))
#
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#     DESCRIPTION
#
#     Very first-order example of a Python version of UNILIB Example 1 program.
#     Evaluation of the magnetic field vector at a given geographic location.
#     The magnetic field is described by the combination of the DGRF/IGRF 1995
#     geomagnetic model and the Olson & Pfitzer dynamic model. The geographic
#     location is provided in geodetic coordinates.
#
#     This example is intended to be pure Python, without any Cython.
#
#     NOTE(1) This example does not open a file with unit kunit, so the
#     output will be directed by default to a file named "fort.1".  Obviously,
#     it should be relatively easy to create a further interface to open and
#     close formatted files with specific unit numbers and names/paths within
#     Python through similar wrappers shown in PyUNILIB.
#
#     NOTE(2) This example must be run from the top-level PyUNILIB directory or
#     a subdirectory.  Obviously more permanent implementations can make this
#     more general.
#
#     NOTE(3) Ensure you have the location of the libUNILIB.so library in the
#     LD_LIBRARY_PATH environment variable.
#
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ==============================================================================
import sys
import os

example = "example2"
currdir = os.getcwd() + "/"
i       = currdir.rfind("/UNILIB/")
if ( i == -1 ) :
  print("The example cannot be run from outside of the UNLIB directory")
  print("or a subdirectory")
else :
  libdir = currdir[:i+8] + "lib"
  sys.path.insert(0, libdir)
  print("Library location is ", libdir)

from PyUNILIB import PyUT990, PyUM510, PyUT540, PyUM520, PyUM522, PyUM524, PyUM536, PyUM530, PyUL220, PyUL225
from PyUNILIB import PyZDAT, PyZGEO, PyZVEC

def main():
#
#
# Initialize the library; the output will be to "fort.1", but you can open a
# suitable unit with number kunit.
#
  kunit  = 1
  kinit  = 1
  ifail  = PyUT990(kunit, kinit)
#  print ("Result from UT990 call = ", ifail)
  if ( ifail < 0 ) :
    print("Error in call to PyUT990.  ifail = ",ifail)
    return
#
#
# Set the geomagnetic field model
# (kint=0, DGRF/IGRF)
#
  kint   = 0
  year   = 2008.0
  lbint, ifail  = PyUM510(kint, year, kunit)
#  print ("Result from UM510 call = ", lbint, ifail)
  if ( ifail < 0 ) :
    print("Error in call to PyUM510.  ifail = ",ifail)
    return
#
#
# Set date and time
#  
  mdate         = PyZDAT()
  mdate.iyear   = 2008
  mdate.imonth  =    1
  mdate.iday    =    1
  mdate.ihour   =    0
  mdate.imin    =    0
  mdate.secs    =    0.0
#
#
# Compute the modified julian day
# (based on January 1, 1950)
#
  mdate1 = PyUT540(mdate)
#  print ("Result from UT540 call = ", mdate1.amjd)
  if ( ifail < 0 ) :
    print("Error in call to PyUT540.  ifail = ",ifail)
    return
#
#
# Set the external magnetic field model
# (kext=6, Olson et Pfizer dynamic)
#
  kext  = 6
  param = [0.0, -30.0, 0.0, 25.0, 300.0, 0.0, 0.0, 0.0, 0.0, 0.0]
  amjd  = mdate1.amjd
  lbext, ifail  = PyUM520(kext, amjd, param, kunit)
#  print ("Result from UM520 call = ", lbext)
  if ( ifail < 0 ) :
    print("Error in call to PyUM520.  ifail = ",ifail)
    return
#
#
# Set the geographic location
#  
  mgeod         = PyZGEO()
  mgeod.radius  = 300.0 + 6371.2
  mgeod.colat   = 120.0
  mgeod.elong   =  30.0
#
#
# Convert from geodetic to geocentric
#
  mpos = PyUM536 (mgeod)
#
#
# Evaluate the magnetic field vector
#
  print(" -------------------------")
  print(" RESULT: magnetic field_shell")
  alpha = [90.,]
  nfbm = 1
  for i in range(24):
        mdate.ihour = i
        mdate1 = PyUT540 (mdate)
        PyUM522 (mdate1.amjd, kunit)
        PyUM524 (kunit)
       
        dbm, dlm, dkm, dsm, fbeq, ds, ifail = PyUL220(mpos, alpha, nfbm)
        
        if ( ifail < 0 ) :
            print("Error in call to PyUL220.  ifail = ",ifail)
            return
#
        frd, fla, ifail = PyUL225(dbm, dlm, 1)
        if ( ifail < 0 ) :
            print("Error in call to PyUL225.  ifail = ",ifail)
            return
#
#     Write the result but to stdout, not kunit in this case.
#
        print("                          amjd = %22.15E " %mdate1.amjd)
        print("              magetic shell Lm = %22.15E " %dlm[0])
        print("              invarant radius  = %22.15E Re" %frd[0])
        print("              invarant latitude= %22.15E deg" %fla[0])
        print(" -------------------------")
#  print("              along colatitude = %22.15E Gauss" %mb.theta)
#  print("               along longitude = %22.15E Gauss" %mb.phi)
#  print(" -------------------------")
  
  return
#
# ==============================================================================
#
if __name__ == "__main__":
  main()
#
# ==============================================================================
#
