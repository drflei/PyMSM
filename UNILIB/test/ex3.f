      PROGRAM example3
C
      INCLUDE 'structure.h'
C
      INTEGER*4    kunit, kinit, ifail, kint, kext
      CHARACTER*32 lbint, lbext
      REAL*8       year, param(10), amjd
      INTEGER*4    knfl, ktyplus
      REAL*8       fbm0, flm0, falt, phi, star
C
C     initialize variables
C
      DATA kunit, kinit/ 6, 1/
      DATA kint, kext/ 0, 0/
      DATA year, param, amjd/ 1995.d0, 10*0.0d0, 0.0d0/
C
C     initialize the library
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
C     set no external magnetic field model
C
C
      CALL UM520 (kext, amjd, param,
     :           lbext, kunit, ifail)
      IF( ifail .LT. 0 )STOP
C
C     trace the drift shell
C
      fbm0          =    0.19d0
      flm0          =    2.0d0
      falt          =    0.0d0
      knfl          =   60
      ktyplus       =    3
C
      CALL UD310 (fbm0, flm0, falt, knfl, 
     :           ktyplus, ifail)
      IF( ifail .LT. 0 )STOP
C
C     evaluate the third invariant
C
      CALL UD330 (phi, star, ifail) 
      IF( ifail .LT. 0 )STOP
C
C     write the result
C
      WRITE(kunit,1000) phi,star
C
 1000 FORMAT(
     :  /25('-'),/'RESULT:',
     :  /'Third invariant =', es24.14, ' ; Lstar =', es24.14,
     :  /25('-') )
C
      END

