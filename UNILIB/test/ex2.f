      PROGRAM example2
C
      INCLUDE 'structure.h'
C
      INTEGER*4    kunit, kinit, ifail, kint, kext, nfbm, i, noprint
      CHARACTER*32 lbint, lbext
      REAL*8       year, param(10), alpha
      REAL*8       fbm, flm, fkm, fsm, fbeq, fs
      TYPE(ZDAT) mdate
      TYPE(ZGEO) mgeod, mpos
C
C     initialisation
C
      DATA kunit, kinit, kint, kext, nfbm, noprint/ 6, 1, 0, 5, 1, -1/
      DATA year, param, alpha/ 1995.0d0, 10*0.0d0, 90.0d0/
C     Note: param is not used when kext=5
C
C     date and time
      mdate%iyear   = 1995
      mdate%imonth  =    1
      mdate%iday    =    1
      mdate%ihour   =    0
      mdate%imin    =    0
      mdate%secs    =    0.0d0
C
C     geodetic location
      mgeod%radius  =20000.0d0 + 6371.2d0
      mgeod%colat   =  120.0d0
      mgeod%elong   =   30.0d0        
C
      CALL UT990 (kunit, kinit, ifail)
      IF( ifail .LT. 0 )STOP
C
      CALL UM510 (kint, year, lbint, kunit, ifail)
      IF( ifail .LT. 0 )STOP
C
      CALL UT540 (mdate)          
C
      CALL UM520 (kext, mdate%amjd, param,
     :           lbext, kunit, ifail)
      IF( ifail .LT. 0 )STOP
C
C     convert the location from geodetic to geocentric coordinates
C
      CALL UM536 (mgeod, mpos)
C
C     loop on the hours
C
      WRITE(kunit,1000)
C
      DO i=0, 23
C
C       update the Sun position and the SM coordinates
C
        mdate%ihour =    i
        CALL UT540 (mdate)          
C
        CALL UM522 (mdate%amjd, noprint)
        CALL UM524 (noprint)
C
C       evaluate the (B, L) coordinates
C
        CALL UL220 (mpos, alpha, nfbm, fbm, flm, fkm, fsm, fbeq, 
     :              fs, ifail)
        IF( ifail .LT. 0 )STOP
C
C       write the result
C
        WRITE(kunit,1001) i, fbm, flm
C
      ENDDO
C
      WRITE(kunit,1002)
C
 1000 FORMAT(
     :  /25('-'),/'RESULT: Magnetic field intensity at mirror points ',
     :   ' and McIlwain shell parameter',
     :  /4x,'hour',15x,'intensity',7x,'shell parameter')
 1001 FORMAT( 4x,i4,2es24.14 )
 1002 FORMAT( /25('-') )
C
      END

