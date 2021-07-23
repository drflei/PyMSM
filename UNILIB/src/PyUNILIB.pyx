# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#     MODULE:           PyUNILIB.py
#
#     Version:          0.A
#     Date:             28/04/2018
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
#     Example/demonstration Python wrapper for some of the UNILIB subroutines, 
#     specifically those needed for UNILIB example1.
#
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ==============================================================================
#from libcpp cimport bool
from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free
from astropy.units import ds, na
#
#
# Declare the interfaces to the C class/Fortran.
#
cdef extern from "../include/PyUNILIB_Interface.h" :
  cdef struct zgeo_ :
    double radius
    double colat
    double elong
    
  cdef struct zvec_ :
    double dnrm
    double rho
    double theta
    double phi

  cdef struct zdat_ :
    double secs
    double amjd
    int    iyear, imonth, iday
    int    ihour, imin, idummy

  cdef void ut990_ ( int *kunit, int *kinit, int *ivar )
  cdef void um510_ ( int *kint, double *year, char *lbint, int *kunit, int *ifail, int lblen )
  cdef void ut540_ ( zdat_ *mdate )
  cdef void um520_ ( int *kext, double *amjd, double param[20], char *lbext, int *kunit, int* ifail, int lblen)
  cdef void um522_ ( double *amjd, int *kunit)
  cdef void um524_ ( int *kunit)
  cdef void um536_ ( zgeo_ *mgde, zgeo_ *mgeo)
  cdef void um530_ ( zgeo_ *mpos, zvec_ *mb, int *ifail)
  cdef void ul220_ ( zgeo_ *mpos, double *alpha, int *nfbm, double *fbm, double *flm, double*, double *fkm, double *fsm, double*fbeq, int *ifail)
  cdef void ul225_ ( double *fbm, double *flm, int *nfbm, double *frd, double *fla, int *ifail);
#
# ==============================================================================
#
# The following defines Python equivalents for the Fortran structs used in
# UNILIB.  Probably a better approach is to use the Python recordtype module
# which is a closer match to the C/Fortran struct.  However, I wanted to avoid
# too much dependence on additional Python modules for this demonstration -
# the objective was to make it as close as possible to pure-Python.
#
# Define struct PyZGEO
#
class PyZGEO (object) :
  def __init__ (self) :
   self.radius = 0.0
   self.colat  = 0.0
   self.elong  = 0.0
#
#
# Define struct PyZVEC
#
class PyZVEC (object) :
  def __init__ (self) :
   self.dnrm  = 0.0
   self.rho   = 0.0
   self.theta = 0.0
   self.phi   = 0.0
#
#
# Define struct PyZDAT
#
class PyZDAT (object) :
  def __init__ (self) :
    self.amjd   =    0.0
    self.iyear  = 1950
    self.imonth =    1
    self.iday   =    0
    self.ihour  =    0
    self.imin   =    0
    self.secs   =    0.0
    self.idummy =    0
#
# ==============================================================================
#
# PyUT990
# -------
#
# Define the wrapper to the Fortran subroutine UT990.
# Initialize the UNILIB libray.
#
cpdef int PyUT990 ( int kunit, int kinit ) :
  cdef int ivar = 0
  ut990_ (&kunit, &kinit, &ivar)
  return ivar
#
# ==============================================================================
#
# PyUM510
# -------
#
# Define the wrapper to the Fortran subroutine UM510.
# Select a geomagnetic field model.
#
cpdef PyUM510 ( int kint, double year, int kunit ) :
  cdef int lblen   = 32
  cdef char *lbint = <char *>PyMem_Malloc(lblen)
  cdef int ifail   = 0
  um510_ ( &kint, &year, lbint, &kunit, &ifail, lblen)
  if lbint == NULL :
    py_lbint = ""
  else :
    py_lbint = lbint.decode("UTF-8","replace")[:lblen].rstrip()
  PyMem_Free(lbint)
    
  return py_lbint, ifail
#
# ==============================================================================
#
# PyUT540
# -------
#
# Define the wrapper to the Fortran subroutine UT540.
# Compute modified Julian Day from date.
#
cpdef PyUT540 ( mdate ) :
  cdef zdat_ s
  s.secs   = mdate.secs
  s.amjd   = mdate.amjd
  s.iyear  = mdate.iyear
  s.imonth = mdate.imonth
  s.iday   = mdate.iday
  s.ihour  = mdate.ihour
  s.imin   = mdate.imin
  s.idummy = mdate.idummy
  ut540_ (&s)
  mdate.secs   = s.secs
  mdate.amjd   = s.amjd
  mdate.iyear  = s.iyear
  mdate.imonth = s.imonth
  mdate.iday   = s.iday
  mdate.ihour  = s.ihour
  mdate.imin   = s.imin
  mdate.idummy = s.idummy
  
  return mdate
#
# ==============================================================================
#
# PyUM520
# -------
#
# Define the wrapper to the Fortran subroutine UM520.
# Select an external magnetic field model.
#
cpdef PyUM520 ( int kext, double amjd, list param, int kunit ) :
  cdef int lblen   = 32
  cdef char *lbext = <char *>PyMem_Malloc(lblen)
  cdef int ifail   = 0
  cdef double param1[20]
  for i in range (20) :
    if ( i < len(param) ) :
      param1[i] = param[i]
    else:
      param1[i] = 0.0
  um520_ ( &kext, &amjd, param1, lbext, &kunit, &ifail, lblen)
  if lbext == NULL :
    py_lbext = ""
  else :
    py_lbext = lbext.decode("UTF-8","replace")[:lblen].rstrip()
  PyMem_Free(lbext)
  
  return py_lbext, ifail

#
# ==============================================================================
#
# PyUM522
# -------
#
# Define the wrapper to the Fortran subroutine UM522.
# Compute the position of the Sun.
#
cpdef PyUM522 ( double amjd, int kunit ) :
  um522_ (&amjd, &kunit)
    
  return

#
# ==============================================================================
#
# PyUM524
# -------
#
# Define the wrapper to the Fortran subroutine UM524.
# Transformation matrix from GEO to SM
#
cpdef PyUM524 ( int kunit ) :
  um524_ (&kunit)
    
  return

# ==============================================================================
#
# PyUM536
# -------
#
# Define the wrapper to the Fortran subroutine UM536.
# Geodetic to geocentric transformation.
#
cpdef PyUM536 ( mgde ) :
  cdef zgeo_ s
  s.radius = mgde.radius
  s.colat  = mgde.colat
  s.elong  = mgde.elong
  cdef zgeo_ q
  um536_ (&s, &q)
  mgeo        = PyZGEO()
  mgeo.radius = q.radius
  mgeo.colat  = q.colat
  mgeo.elong  = q.elong
  
  return mgeo
#
# ==============================================================================
#
# PyUM530
# -------
#
# Define the wrapper to the Fortran subroutine UM530.
# Evaluate the magnetic field vector.
#
cpdef PyUM530 ( mgeo ) :
  cdef zgeo_ s
  s.radius = mgeo.radius
  s.colat  = mgeo.colat
  s.elong  = mgeo.elong
  cdef zvec_ v
  cdef int ifail = 0
  um530_ (&s, &v, &ifail)
  mb       = PyZVEC()
  mb.dnrm  = v.dnrm
  mb.rho   = v.rho
  mb.theta = v.theta
  mb.phi   = v.phi
  
  return mb, ifail
#
# ==============================================================================
#
# PyUL220
# -------
#
# Define the wrapper to the Fortran subroutine UL220.
# Get information on a magnetic field line segment.
#
cpdef PyUL220 ( mgeo, list alpha, int nfbm ) :
  cdef zgeo_ s
  s.radius = mgeo.radius
  s.colat  = mgeo.colat
  s.elong  = mgeo.elong
  
  cdef int na  = nfbm
  cdef double al[20]
  
  for i in range (na) :
    if ( i < len(alpha) ) :
      al[i] = alpha[i]
    else:
      al[i] = 0.0

  cdef int ifail = 0
  cdef double fbm[20]
  cdef double flm[20]
  cdef double fkm[20]
  cdef double fsm[20]
  cdef double fs[20]
  cdef double fbeq
  
  ul220_ (&s, al, &na, fbm, flm, fkm, fsm, &fbeq, fs, &ifail)
  
  cdef list dbm = []
  cdef list dlm = []
  cdef list dkm = []
  cdef list dsm = []
  cdef list ds = []
  
  for i in range(na) :
      dbm.append(fbm[i])
      dlm.append(flm[i])
      dkm.append(fkm[i])
      dsm.append(fsm[i])
      ds.append(fs[i])

  return dbm, dlm, dkm, dsm, fbeq, ds, ifail
#
# ==============================================================================
#
# PyUL225
# -------
#
# Define the wrapper to the Fortran subroutine UL225.
# evaluate the invariant radius and latitude
#
cpdef PyUL225 ( list fbm, list flm, int nfbm ) :

  cdef int na  = nfbm
  cdef double bm[20]
  cdef double lm[20]
  
  for i in range (na) :
    if ( i < len(flm) ) :
      bm[i] = fbm[i]
      lm[i] = flm[i]
    else:
      bm[i] = 0.0
      lm[i] = 0.0
  
  cdef int ifail = 0
  cdef double frd[20]
  cdef double fla[20]
  
  ul225_ (bm, lm, &na, frd, fla, &ifail)
  
  cdef list drd = []
  cdef list dla = []
  
  for i in range(na) :
      drd.append(frd[i])
      dla.append(fla[i])

  return drd, dla, ifail
