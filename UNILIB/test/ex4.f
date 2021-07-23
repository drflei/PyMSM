      PROGRAM example4
C
      INCLUDE 'structure.h'
C
      INTEGER*4    kunit, kinit, ifail, kint, kext
      CHARACTER*32 lbint, lbext
      REAL*8       year, param(10), amjd
      INTEGER*4    knfl, ktyplus
      REAL*8       fbm0, flm0, falt
      TYPE(ZGEO) mpn, mps
C
C     initialize variables
C
      DATA kunit, kinit/ 6, 1/
      DATA kint, kext/ 2, 0/
      DATA year, param, amjd/ 1970.d0, 10*0.d0, 0.d0/
C
C     initialize the library
C
      CALL UT990 (kunit, kinit, ifail)
      IF( ifail .LT. 0 )STOP
C
C     set the geomagnetic field model
C     (kint=2, GFSC 12/66)
C
      CALL UM510 (kint, year, lbint, kunit, ifail)
      IF( ifail .LT. 0 )STOP
C
C     set no external magnetic field model
C
      CALL UM520 (kext, amjd, param,
     :           lbext, kunit, ifail)
      IF( ifail .LT. 0 )STOP
C
C     trace the drift shell
C
      fbm0          =    0.19d0
      flm0          =    2.0d0
      falt          = -999.9d0
      knfl          =   60
      ktyplus       =    3
C
      CALL UD310 (fbm0, flm0, falt, knfl, 
     :           ktyplus, ifail)
      IF( ifail .LT. 0 )STOP
C
      CALL UD315 (mpn, mps, ifail)
      IF( ifail .LT. 0 )STOP
C
C     write the result
C         
      WRITE (kunit,1000) mpn%radius, mpn%colat, mpn%elong,
     :                   mps%radius, mps%colat, mps%elong
C
 1000 FORMAT(
     :  /25('-'),/'RESULT: coordinates (radius, colatitude, longitude)',
     :   ' of the points with the lowest altitude',
     :  /'North : ', 3f20.10,
     :  /'South : ', 3f20.10,
     :  /25('-') )
C
      END

