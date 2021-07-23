      PROGRAM example1
C
      INCLUDE 'structure.h'
C
      INTEGER*4    kunit, kinit, ifail, kint, kext
      CHARACTER*32 lbint, lbext
      REAL*8       year, param(10)
      TYPE(ZDAT) mdate
      TYPE(ZGEO) mgeod, mpos
      TYPE(ZVEC) mb
C
C     initialize variables
C
      DATA kunit, kinit/ 6, 1/
      DATA kint, kext/ 0, 6/
      DATA year, param/ 1995.d0, 0.d0, -30.d0, 0.d0, 25.d0, 300.d0, 
     :                           0.d0, 0.d0, 0.d0, 0.d0, 0.d0/
C
C     initialize the library
C
C
      CALL UT990 (kunit, kinit, ifail)
      IF( ifail .LT. 0 )STOP
C
C     set the geomagnetic field model
C     (kint=0, DGRF/IGRF)
C
      CALL UM510 (kint, year, lbint, kunit, ifail)
      IF( ifail .LT. 0 )STOP
C
C     set date and time
C
      mdate%iyear   = 1995
      mdate%imonth  =    1
      mdate%iday    =    1
      mdate%ihour   =    0
      mdate%imin    =    0
      mdate%secs    =    0.0d0
C
C     compute the modified julian day
C     (based on January 1, 1950)
C
      CALL UT540 (mdate)          
C
C     set the external magnetic field model
C     (kext=6, Olson et Pfizer dynamic)
C
      CALL UM520 (kext, mdate%amjd, param,
     :           lbext, kunit, ifail)
      IF( ifail .LT. 0 )STOP
C
C     set the geographic location
C
      mgeod%radius  =  300.0d0 + 6371.2d0
      mgeod%colat   =  120.0d0
      mgeod%elong   =   30.0d0
C
C     convert from geodetic to geocentric
C
      CALL UM536 (mgeod, mpos)
C
C     evaluate the magnetic field vector
C
      CALL UM530 (mpos, mb, ifail)
      IF( ifail .LT. 0 )STOP
C
C     write the result
C
      WRITE(kunit,1000) mb%dnrm, mb%rho, mb%theta, mb%phi
C
 1000 FORMAT(
     :  /25('-'),/'RESULT: magnetic field',/25x,
     :               'norm =', es24.15, ' gauss',/13x,
     :   'radial component =', es24.15, ' gauss',/13x,
     :   'along colatitude =', es24.15, ' gauss',/14x,
     :    'along longitude =', es24.15, ' gauss',/25('-') )
C
      END 
