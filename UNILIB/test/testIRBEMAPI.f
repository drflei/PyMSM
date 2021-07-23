      IMPLICIT NONE

      INCLUDE 'structure.h'

      INCLUDE 'IRBEMAPI.inc'

      INTEGER*4 NALP
      PARAMETER (NALP=25)
      REAL*8 SECS, UT, XIN(3), XOUT(3), MLT, MAGINPUT(25), BXGEO(3), BL
      REAL*8 XGDZ(3), XGSM(3), XGEI(3), XGEO(3), ALTI, LATI, LONGI, AMJD
      INTEGER*4 KUNIT /6/, KINIT /1/, IFAIL, VERSION, NTIME
      INTEGER*4 SYSAXESIN, SYSAXESOUT, IYEAR, IMONTH, IDAY, IHOUR
      INTEGER*4 IMINUTE, IDOY, KEXT, OPTIONS(5)
      CHARACTER*80 RELEASE
      CHARACTER*32 LBINT

      TYPE(zgeo) msph
      TYPE(zxyz) mxyz

      REAL*8 UTM(NTIME_MAX), XIN1(NTIME_MAX), XIN2(NTIME_MAX)
      REAL*8 XIN3(NTIME_MAX), MAGINPUTM(25,NTIME_MAX), LM(NTIME_MAX)
      REAL*8 LSTAR(NTIME_MAX), BLOCAL(NTIME_MAX), BMIN(NTIME_MAX)
      REAL*8 XJ(NTIME_MAX), MLTM(NTIME_MAX), PLM(NTIME_MAX,NALP)
      REAL*8 PLSTAR(NTIME_MAX,NALP), PBLOCAL(NTIME_MAX,NALP)
      REAL*8 PXJ(NTIME_MAX,NALP), PHI(NTIME_MAX), ALPHA(NALP)
      INTEGER*4 YEARM(NTIME_MAX), DOYM(NTIME_MAX), NIPA

C Functions
      REAL*8 YDOYUT_TO_JULDAY
      INTEGER*4 GET_DOY
      EXTERNAL YDOYUT_TO_JULDAY, GET_DOY

      CALL UT990(KUNIT, KINIT, IFAIL)
      WRITE(6, *)

C UNILIB release
      CALL IRBEM_FORTRAN_RELEASE1(RELEASE)
      WRITE(6, '(2A,/)') 'UNILIB release: ', TRIM(ADJUSTL(RELEASE))

C UNILIB version
      CALL IRBEM_FORTRAN_VERSION1(VERSION)
      WRITE(6, '(A,I6,/)') 'UNILIB version: ', version

C Size of pre-allocated arrays
      CALL GET_IRBEM_NTIME_MAX1(NTIME)
      WRITE(6, '(A,I8,/)') 'Maximum array size: ', NTIME

C Set a date and time (2010-05-17 06:00:25)
      IYEAR = 2010
      IMONTH = 5
      IDAY = 17
      IHOUR = 6
      IMINUTE = 0
      SECS = 25.0D0
      UT = IHOUR * 3600.0D0 + IMINUTE * 60.0D0 + SECS

C Day of year
      IDOY = GET_DOY(IYEAR, IMONTH, IDAY)
      WRITE(6, '(A,I4,/)') 'Day of year:', IDOY

C Coordinate transformation
C GDZ -> GEI
      SYSAXESIN = 0
      SYSAXESOUT = 5
      XGDZ(1) = 1000.0D0
      XGDZ(2) = 50.88D0
      XGDZ(3) = 4.7D0
      WRITE(6, '(A,F6.1,2(1H,,F7.2),A)')
     &  'GDZ coordinates: (', XGDZ, ')'
      CALL COORD_TRANS1(SYSAXESIN, SYSAXESOUT, IYEAR, IDOY, UT, XGDZ,
     &                  XGEI)
      WRITE(6, '(A,1P,2(E13.6,1H,),E13.6,0P,A)')
     &  'GEI coordinates (Cartesian): (', XGEI, ')'
      mxyz%X = XGEI(1)
      mxyz%Y = XGEI(2)
      mxyz%Z = XGEI(3)
      CALL UT546(mxyz, msph)
      WRITE(6, '(A,F6.1,2(1H,,F7.2),A)')
     &  'GEI coordinates (spherical): (', msph%RADIUS-6371.2D0,
     &  90.0D0-msph%COLAT, MOD(msph%ELONG+360.0D0, 360.0D0), ')'
      WRITE(6, '(A,F6.1,2(1H,,F7.2),A)')
     &  'GEI coordinates (spherical, SPENVIS): (', 7365.33D0-6371.2D0,
     &  90.0D0-39.28D0, MOD(-30.36D0+360.0D0, 360.0D0), ')'
C GEI -> GDZ
      SYSAXESIN = 5
      SYSAXESOUT = 0
      CALL COORD_TRANS1(SYSAXESIN, SYSAXESOUT, IYEAR, IDOY, UT, XGEI,
     &                  XOUT)
      WRITE(6, '(A,F6.1,2(1H,,F7.2),A)')
     &  'GDZ coordinates: (', XOUT, ')'
C GEI -> GSM
      CALL UM510(0, 2010.0D0, LBINT, 0, IFAIL)
      AMJD = YDOYUT_TO_JULDAY(IYEAR, IDOY, UT) - 2433282.5D0
      CALL UM522(AMJD, 0)
      SYSAXESIN = 5
      SYSAXESOUT = 2
      CALL COORD_TRANS1(SYSAXESIN, SYSAXESOUT, IYEAR, IDOY, UT, XGEI,
     &                  XGSM)
      mxyz%X = XGSM(1)
      mxyz%Y = XGSM(2)
      mxyz%Z = XGSM(3)
      CALL UT546(mxyz, msph)
      WRITE(6, '(A,F7.1,2(1H,,F7.2),A)')
     &  'GSM coordinates (spherical): (', msph%RADIUS-6371.2D0,
     &  90.0D0-msph%COLAT, MOD(msph%ELONG+360.0D0, 360.0D0), ')'
      WRITE(6, '(A,F7.1,2(1H,,F7.2),A)')
     &  'GSM coordinates (spherical, SPENVIS): (', 7365.33D0-6371.2D0,
     &  90.0D0-41.91D0, MOD(-61.82D0+360.0D0, 360.0D0), ')'
C GEI -> GEO
      SYSAXESIN = 5
      SYSAXESOUT = 1
      CALL COORD_TRANS1(SYSAXESIN, SYSAXESOUT, IYEAR, IDOY, UT, XGEI,
     &                  XGEO)
      mxyz%X = XGEO(1)
      mxyz%Y = XGEO(2)
      mxyz%Z = XGEO(3)
      CALL UT546(mxyz, msph)
      WRITE(6, '(A,F7.1,2(1H,,F7.2),A)')
     &  'GEO coordinates (spherical): (', msph%RADIUS-6371.2D0,
     &  90.0D0-msph%COLAT, MOD(msph%ELONG+360.0D0, 360.0D0), ')'
      WRITE(6, '(A,F7.1,2(1H,,F7.2),A)')
     &  'GEO coordinates (spherical, SPENVIS): (', 7365.33D0-6371.2D0,
     &  90.0D0-39.28D0, MOD(4.70D0+360.0D0, 360.0D0), ')'

C Magnetic local time
      CALL GET_MLT1(IYEAR, IDOY, UT, XGEO, MLT)
      WRITE(6, '(A,F5.2)') 'Magnetic local time: ', MLT

C Magnetic field vector
      SYSAXESIN = 1
      OPTIONS(2) = 0
      OPTIONS(5) = 0
      KEXT = 5
      MAGINPUT(1) = 0.0D0
      CALL GET_FIELD1(KEXT, OPTIONS, SYSAXESIN, IYEAR, IDOY, UT,
     &                XGEO(1), XGEO(2), XGEO(3), MAGINPUT, BXGEO, BL)
      WRITE(6, '(A,F8.1,3(1H,,F8.1))')
     &  'B strength and vector (GEO Cart., OPQ, nT): ', BL, BXGEO

C GSM -> GEO, using GET_COORDINATES
      SYSAXESIN = 2
      CALL GET_COORDINATES(SYSAXESIN, XGSM(1), XGSM(2), XGSM(3), ALTI,
     &                     LATI, LONGI, XOUT)
      mxyz%X = XOUT(1)
      mxyz%Y = XOUT(2)
      mxyz%Z = XOUT(3)
      CALL UT546(mxyz, msph)
c      WRITE(6, '(A,F7.1,2(1H,,F7.2),A)')
c     &  'GEO coordinates (spherical): (', msph%RADIUS-6371.2D0,
c     &  90.0D0-msph%COLAT, MOD(msph%ELONG+360.0D0, 360.0D0), ')'
c      WRITE(6, '(A,F6.1,2(1H,,F7.2),A)')
c     &  'GDZ coordinates: (', ALTI, LATI, LONGI, ')'

C L* for 90deg pitch angle
      YEARM(1) = IYEAR
      DOYM(1) = IDOY
      UTM(1) = UT
      XIN1(1) = XGEO(1)
      XIN2(1) = XGEO(2)
      XIN3(1) = XGEO(3)
      MAGINPUTM(1:25,1) = MAGINPUT
      CALL MAKE_LSTAR1(1, KEXT, OPTIONS, 1, YEARM, DOYM, UTM,
     &                 XIN1, XIN2, XIN3, MAGINPUTM, LM, LSTAR, BLOCAL,
     &                 BMIN, XJ, MLTM)
      WRITE(6, '(/,A)') 'MAKE_LSTAR1'
      WRITE(6, '(A,F6.3,1H,,F7.3)') 'L, L*: ', LM(1), LSTAR(1)
      WRITE(6, '(A,F8.1,1H,,F9.1)') 'B, Beq: ', BLOCAL(1), BMIN(1)
      WRITE(6, '(A,F6.3)') 'I: ', XJ(1)
      WRITE(6, '(A,F5.2)') 'Magnetic local time: ', MLTM(1)

C Phi <> L*
      CALL LSTAR_PHI1(1, 1, OPTIONS, IYEAR, IDOY, LSTAR, PHI)
      WRITE(6, '(/,A,F7.3,A,1P,E11.3,0P)') 'L* -> Phi: ', LSTAR(1),
     &  ' -> ', PHI(1)
      CALL LSTAR_PHI1(1, 2, OPTIONS, IYEAR, IDOY, LSTAR, PHI)
      WRITE(6, '(/,A,1P,E11.3,0P,A,F7.3)') 'Phi -> L*: ', PHI(1),
     &  ' -> ', LSTAR(1)

C L* for multiple pitch angles
      YEARM(1) = IYEAR
      DOYM(1) = IDOY
      UTM(1) = UT
      XIN1(1) = XGEO(1)
      XIN2(1) = XGEO(2)
      XIN3(1) = XGEO(3)
      MAGINPUTM(1:25,1) = MAGINPUT
      NIPA = 4
      ALPHA(1) = 150.0D0
      ALPHA(2) = 75.0D0
      ALPHA(3) = 90.0D0
      ALPHA(4) = 85.0D0
      CALL MAKE_LSTAR_SHELL_SPLITTING1(1, NIPA, KEXT, OPTIONS, 1, YEARM,
     &                                 DOYM, UTM, XIN1, XIN2, XIN3,
     &                                 ALPHA, MAGINPUTM, PLM, PLSTAR,
     &                                 PBLOCAL, BMIN, PXJ, MLTM)
      WRITE(6, '(/,A)') 'MAKE_LSTAR_SHELL_SPLITTING1'
      WRITE(6, '(A,10(F6.3:,1H,))') 'L:  ', PLM(1,1:NIPA)
      WRITE(6, '(A,10(F6.3:,1H,))') 'L*: ', PLSTAR(1,1:NIPA)
      WRITE(6, '(A,10(F9.1:,1H,))') 'B: ', PBLOCAL(1,1:NIPA)
      WRITE(6, '(A,F8.1)') 'Beq: ', BMIN(1)
      WRITE(6, '(A,10(F6.3:,1H,))') 'I: ', PXJ(1,1:NIPA)
      WRITE(6, '(A,F5.2)') 'Magnetic local time: ', MLTM(1)


      END
