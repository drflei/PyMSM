


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Get UNILIB release string
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE IRBEM_FORTRAN_RELEASE1(RELEASE)

      IMPLICIT NONE

      CHARACTER*80 RELEASE
      CALL myversion(RELEASE)

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Get UNILIB version string
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE IRBEM_FORTRAN_VERSION1(VERSION)

      IMPLICIT NONE

      INTEGER*4 VERSION
      VERSION = 2230

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Get maximum array size
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE GET_IRBEM_NTIME_MAX1(NTIME)

      IMPLICIT NONE

      INCLUDE 'IRBEMAPI.inc'

      INTEGER*4 NTIME

      NTIME = NTIME_MAX

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Day of year
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER*4 FUNCTION GET_DOY(IYEAR, IMONTH, IDAY)

      IMPLICIT NONE

      INCLUDE 'structure.h'

      INTEGER*4 IYEAR, IMONTH, IDAY

      TYPE(ZDAT) MDATE, MDATE1

      MDATE%IYEAR = IYEAR
      MDATE%IMONTH = IMONTH
      MDATE%IDAY = IDAY
      MDATE%IHOUR = 0
      MDATE%IMIN = 0
      MDATE%SECS = 0.0D0
      CALL UT540(MDATE)
      MDATE1%IYEAR = IYEAR
      MDATE1%IMONTH = 1
      MDATE1%IDAY = 1
      MDATE1%IHOUR = 0
      MDATE1%IMIN = 0
      MDATE1%SECS = 0.0D0
      CALL UT540(MDATE1)

      GET_DOY = MDATE%AMJD - MDATE1%AMJD + 1

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Convert year, day of year, UT (in seconds) to Julian day
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      REAL*8 FUNCTION YDOYUT_TO_JULDAY(YEAR, DOY, UT)

      IMPLICIT NONE

      INCLUDE 'structure.h'

      REAL*8 UT
      INTEGER*4 YEAR, DOY

      TYPE(zdat) AMJD

      AMJD%IYEAR = YEAR
      AMJD%IMONTH = 1
      AMJD%IDAY = 1
      AMJD%IHOUR = 0
      AMJD%IMIN = 0
      AMJD%SECS = 0.0D0
      CALL UT540(AMJD)
      YDOYUT_TO_JULDAY = AMJD%AMJD + DOY - 1 + UT / 86400.0D0 +
     &                   2433282.5D0

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Convert Julian day to year, day of year, UT (in seconds)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE JULDAY_TO_YDOYUT(JULDAY, YEAR, DOY, UT)

      IMPLICIT NONE

      INCLUDE 'structure.h'

      REAL*8 JULDAY, UT
      INTEGER*4 YEAR, DOY

C Functions
      INTEGER*4 GET_DOY
      EXTERNAL GET_DOY

      TYPE(zdat) AMJD

      AMJD%AMJD = JULDAY - 2433282.5D0
      CALL UT545(AMJD)

      YEAR = AMJD%IYEAR
      DOY = GET_DOY(YEAR, AMJD%IMONTH, AMJD%IDAY)
      UT = AMJD%IHOUR * 3600.0D0 + AMJD%IMIN * 60.0D0 + AMJD%SECS

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Coordinate transformations (single location)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE COORD_TRANS1(SYSAXESIN, SYSAXESOUT, YEAR, DOY, UT, XIN,
     &                        XOUT)

      IMPLICIT NONE

      INCLUDE 'structure.h'

      INCLUDE 'IRBEMAPI.inc'

      REAL*8 UT, XIN(3), XOUT(3), AMJD, TRANS(3,3)
      INTEGER*4 SYSAXESIN, SYSAXESOUT, YEAR, DOY
      INTEGER*4 KFROM, KTO, IFAIL, KUNIT /0/
C Correspondence table between IRBEM and UNILIB coordinate systems
      INTEGER*4 IRBEMTOUNILIB(6) /1,5,6,4,2,3/

      TYPE(zgeo) msph, mgde, mgeo, mout
      TYPE(zxyz) mxyz

C Functions
      REAL*8 YDOYUT_TO_JULDAY
      EXTERNAL YDOYUT_TO_JULDAY

C Exit for IRBEM coordinate systems not supported in UNILIB
      IF ((SYSAXESIN .LT. 0) .OR. (SYSAXESIN .GT. 8) .OR.
     &    (SYSAXESOUT .LT. 0) .OR. (SYSAXESOUT .GT. 8)) THEN
        XOUT = BADDATA
        RETURN
      END IF

      IF (SYSAXESOUT .EQ. SYSAXESIN) THEN
        XOUT = XIN
        RETURN
      END IF

      AMJD = YDOYUT_TO_JULDAY(YEAR, DOY, UT) - 2433282.5D0
      CALL UM522(AMJD, KUNIT)

C Input coordinates
C     Geodetic
      IF (SYSAXESIN .EQ. 0) THEN
        mgde%RADIUS = XIN(1) + 6371.2D0
        mgde%COLAT = 90.0D0 - XIN(2)
        mgde%ELONG = XIN(3)
        CALL UM536(mgde, mgeo)
C       GEO Cartesian
        IF (SYSAXESOUT .EQ. 1) THEN
          CALL UT541(mgeo, mxyz)
          XOUT(1) = mxyz%X
          XOUT(2) = mxyz%Y
          XOUT(3) = mxyz%Z
C       GEO spherical
        ELSE IF ((SYSAXESOUT .EQ. 7) .OR. (SYSAXESOUT .EQ. 8)) THEN
          XOUT(1) = mgeo%RADIUS / 6371.2D0
          XOUT(2) = 90.0D0 - mgeo%COLAT
          XOUT(3) = mgeo%ELONG
        ELSE
          KFROM = 1
          KTO = IRBEMTOUNILIB(SYSAXESOUT)
          CALL UT550(KFROM, KTO, TRANS, IFAIL)
          CALL UT555(mgeo, mout, TRANS)
          CALL UT541(mout, mxyz)
          XOUT(1) = mxyz%X
          XOUT(2) = mxyz%Y
          XOUT(3) = mxyz%Z
        END IF
      ELSE
C       GEO spherical
        IF ((SYSAXESIN .EQ. 7) .OR. (SYSAXESIN .EQ. 8)) THEN
          msph%RADIUS = XIN(1) * 6371.2D0
          msph%COLAT = 90.0D0 - XIN(2)
          msph%ELONG = XIN(3)
        ELSE
          mxyz%X = XIN(1)
          mxyz%Y = XIN(2)
          mxyz%Z = XIN(3)
          CALL UT546(mxyz, msph)
        END IF
        KFROM = IRBEMTOUNILIB(SYSAXESIN)
        IF ((SYSAXESOUT .EQ. 0) .OR.
     &      (SYSAXESOUT .EQ. 7) .OR. (SYSAXESOUT .EQ. 8)) THEN
          CALL UT550(KFROM, 1, TRANS, IFAIL)
          CALL UT555(msph, mgeo, TRANS)
          IF (SYSAXESOUT .EQ. 0) THEN
            CALL UM535(mgeo, mgde)
            XOUT(1) = mgde%RADIUS - 6371.2D0
            XOUT(2) = 90.0D0 - mgde%COLAT
            XOUT(3) = mgde%ELONG
          ELSE IF ((SYSAXESOUT .EQ. 7) .OR. (SYSAXESOUT .EQ. 8)) THEN
            XOUT(1) = mgeo%RADIUS / 6371.2D0
            XOUT(2) = 90.0D0 - mgeo%COLAT
            XOUT(3) = mgeo%ELONG
          END IF
        ELSE
          KTO = IRBEMTOUNILIB(SYSAXESOUT)
          CALL UT550(KFROM, KTO, TRANS, IFAIL)
          CALL UT555(msph, mout, TRANS)
          CALL UT541(mout, mxyz)
          XOUT(1) = mxyz%X
          XOUT(2) = mxyz%Y
          XOUT(3) = mxyz%Z
        END IF
      END IF

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Coordinate transformations (multiple locations)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE COORD_TRANS_VEC1(NTIME, SYSAXESIN, SYSAXESOUT, YEAR,
     &                            DOY, UT, XINV, XOUTV)

      IMPLICIT NONE

      INCLUDE 'IRBEMAPI.inc'

      REAL*8 UT(NTIME_MAX), XINV(3,NTIME_MAX), XOUTV(3,NTIME_MAX)
      INTEGER*4 NTIME, SYSAXESIN, SYSAXESOUT, YEAR(NTIME_MAX)
      INTEGER*4 DOY(NTIME_MAX), I

      IF (SYSAXESOUT .EQ. SYSAXESIN) THEN
        XOUTV(1:3,1:NTIME) = XINV(1:3,1:NTIME)
      ELSE
        DO I=1,NTIME
          CALL COORD_TRANS1(SYSAXESIN, SYSAXESOUT, YEAR(I), DOY(I),
     &                      UT(I), XINV(1:3,I), XOUTV(1:3,I))
        END DO
      END IF

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Get GEO spherical and Cartesian coordinates for a single location
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE GET_COORDINATES(SYSAXES, XIN1, XIN2, XIN3, ALTI, LATI,
     &                           LONGI, XGEO)

      IMPLICIT NONE

      INCLUDE 'IRBEMAPI.inc'

      REAL*8 XIN1, XIN2, XIN3, ALTI, LATI, LONGI, XGEO(3), XIN(3)
      REAL*8 XGEOD(3), UT /0.0D0/
      INTEGER*4 SYSAXES, YEAR

      XIN(1) = XIN1
      XIN(2) = XIN2
      XIN(3) = XIN3
C GEO Cartesian
      CALL COORD_TRANS1(SYSAXES, 1, 2010, 137, 21625.0D0, XIN, XGEO)
C Geodetic
      CALL COORD_TRANS1(1, 0, 0, 0, 0.0D0, XGEO, XGEOD)
      ALTI = XGEOD(1)
      LATI = XGEOD(2)
      LONGI = XGEOD(3)

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Calculate magnetic local time
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE GET_MLT1(YEAR, DOY, UT, XGEO, MLT)

      IMPLICIT NONE

      INCLUDE 'structure.h'

      INCLUDE 'IRBEMAPI.inc'

      REAL*8 UT, XGEO(3), MLT, MLAT, AMJD
      INTEGER*4 YEAR, DOY, IFAIL, KUNIT /0/

      TYPE(zgeo) SGEO

      TYPE(zxyz) CGEO

C Functions
      REAL*8 YDOYUT_TO_JULDAY
      EXTERNAL YDOYUT_TO_JULDAY

      CGEO%X = XGEO(1)
      CGEO%Y = XGEO(2)
      CGEO%Z = XGEO(3)
      CALL UT546(CGEO, SGEO)
      AMJD = YDOYUT_TO_JULDAY(YEAR, DOY, UT) - 2433282.5D0
      CALL UM522(AMJD, KUNIT)
      CALL UM538(SGEO, AMJD, MLT, MLAT, IFAIL)
      IF ((IFAIL .LT. 0) .OR. ISNAN(MLT))
     &  MLT = BADDATA

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Initialise fields
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE INITIALISE_FIELDS(OPTIONS, YEAR, DOY, UT, KEXT,
     &                             MAGINPUT)

      IMPLICIT NONE

      REAL*8 UT, MAGINPUT(25), PARAM(20), AMJD
      INTEGER*4 OPTIONS(5), YEAR, DOY, KEXT, IFAIL, KINT, KEXTU
      INTEGER*4 KUNIT /0/, OYEAR, I
      CHARACTER*32 LBINT, LBEXT

      SAVE OYEAR

C Functions
      REAL*8 YDOYUT_TO_JULDAY
      EXTERNAL YDOYUT_TO_JULDAY

C Internal field model selection
      IF (OPTIONS(5) .EQ. 0) THEN
        KINT = 0
      ELSE IF (OPTIONS(5) .EQ. 2) THEN
        KINT = 1
      ELSE IF (OPTIONS(5) .EQ. 3) THEN
        KINT = 2
      ELSE IF (OPTIONS(5) .EQ. 5) THEN
        KINT = 3
      ELSE
        RETURN
      END IF
C Set field epoch and initialise internal field
      IF ((OPTIONS(2) .EQ. 0) .AND. (YEAR .NE. OYEAR)) THEN
        OYEAR = YEAR
        CALL UM510(KINT, YEAR+0.5D0, LBINT, KUNIT, IFAIL)
      END IF
      IF (OPTIONS(2) .NE. 0)
     &  CALL UM510(KINT, YEAR+(DOY+UT/86400.0D0)/365.25D0, LBINT,
     &             KUNIT, IFAIL)
      IF (IFAIL .LT. 0)
     &  RETURN

C External field model selection
      IF ((KEXT .GE. 0) .AND. (KEXT .LE. 9)) THEN
        KEXTU = KEXT
      ELSE IF ((KEXT .EQ. 11) .OR. (KEXT .EQ. 12)) THEN
        KEXTU = KEXT - 1
      ELSE
        RETURN
      END IF
C Magnetic field parameters
      PARAM(1) = MAGINPUT(1) / 10.0D0
      PARAM(2) = MAGINPUT(2)
      PARAM(3) = MAGINPUT(5)
      PARAM(4) = MAGINPUT(3)
      PARAM(5) = MAGINPUT(4)
      PARAM(7) = MAGINPUT(6)
      PARAM(8) = MAGINPUT(7)
      PARAM(10) = MAGINPUT(17)
      PARAM(11) = MAGINPUT(8)
      PARAM(12) = MAGINPUT(9)
      DO I=13,18
        PARAM(I) = MAGINPUT(I-2)
      END DO
C Set field epoch and initialise external field
      AMJD = YDOYUT_TO_JULDAY(YEAR, DOY, UT) - 2433282.5D0
      CALL UM520(KEXTU, AMJD, PARAM, LBEXT, KUNIT, IFAIL)

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Calculate magnetic field vector components and strength
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE GET_FIELD1(KEXT, OPTIONS, SYSAXES, YEAR, DOY, UT, XIN1,
     &                      XIN2, XIN3, MAGINPUT, BXGEO, BL)

      IMPLICIT NONE

      INCLUDE 'structure.h'

      INCLUDE 'IRBEMAPI.inc'

      REAL*8 UT, XIN1, XIN2, XIN3, BXGEO(3), BL, XIN(3), XOUT(3)
      REAL*8 MAGINPUT(25)
      INTEGER*4 KEXT, OPTIONS(5), SYSAXES, YEAR, DOY, IFAIL, KUNIT /0/

      TYPE(zgeo) MPOS

      TYPE(zvec) MBSPH

      TYPE(zxyz) MXYZ, MBCAR

C Functions
      REAL*8 YDOYUT_TO_JULDAY
      EXTERNAL YDOYUT_TO_JULDAY

      BXGEO = BADDATA
      BL = BADDATA

C Not all IRBEM magnetic field models and options are supported
C by UNILIB
      IF ((OPTIONS(5) .LT. 0) .OR. (OPTIONS(5) .EQ. 1) .OR.
     &    (OPTIONS(5) .EQ. 4) .OR. (OPTIONS(5) .GT. 5) .OR.
     &    (KEXT .LT. 0) .OR. (KEXT .GT. 12))
     &  RETURN

      CALL INITIALISE_FIELDS(OPTIONS, YEAR, DOY, UT, KEXT, MAGINPUT)
      XIN(1) = XIN1
      XIN(2) = XIN2
      XIN(3) = XIN3
      CALL COORD_TRANS1(SYSAXES, 1, YEAR, DOY, UT, XIN, XOUT)
      MXYZ%X = XOUT(1)
      MXYZ%Y = XOUT(2)
      MXYZ%Z = XOUT(3)
      CALL UT546(MXYZ, MPOS)
      CALL UM530(MPOS, MBSPH, IFAIL)
      IF (IFAIL .LT. 0) THEN
        RETURN
      ELSE
        BL = MBSPH%DNRM * 1.0D5
        CALL UT542(MPOS, MBSPH, MBCAR)
        BXGEO(1)= MBCAR%X * 1.0D5
        BXGEO(2)= MBCAR%Y * 1.0D5
        BXGEO(3)= MBCAR%Z * 1.0D5
      END IF

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Calculate magnetic field vector components and strength
C for multiple locations
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE GET_FIELD_MULTI(NTIME, KEXT, OPTIONS, SYSAXES, YEAR,
     &                           DOY, UT, XIN1, XIN2, XIN3, MAGINPUT,
     &                           BXGEO, BL)

      IMPLICIT NONE

      INCLUDE 'IRBEMAPI.inc'

      REAL*8 UT(NTIME_MAX), XIN1(NTIME_MAX), XIN2(NTIME_MAX)
      REAL*8 XIN3(NTIME_MAX), MAGINPUT(25,NTIME_MAX), BXGEO(3,NTIME_MAX)
      REAL*8 BL(NTIME_MAX), XIN(3), XOUT(3)
      INTEGER*4 NTIME, KEXT, OPTIONS(5), SYSAXES, YEAR(NTIME_MAX)
      INTEGER*4 DOY(NTIME_MAX), I

      DO I=1,NTIME
        CALL GET_FIELD1(KEXT, OPTIONS, SYSAXES, YEAR(I), DOY(I), UT(I),
     &                  XIN1(I), XIN2(I), XIN3(I), MAGINPUT(1,I),
     &                  BXGEO(1,I), BL(I))
      END DO

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Calculate L* for pitch angle 90deg
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE MAKE_LSTAR1(NTIME, KEXT, OPTIONS, SYSAXES, YEAR, DOY,
     &                       UT, XIN1, XIN2, XIN3, MAGINPUT, LM, LSTAR,
     &                       BLOCAL, BMIN, XJ, MLT)

      IMPLICIT NONE

      INCLUDE 'structure.h'

      INCLUDE 'IRBEMAPI.inc'

      REAL*8 UT(NTIME_MAX), XIN1(NTIME_MAX), XIN2(NTIME_MAX)
      REAL*8 XIN3(NTIME_MAX), MAGINPUT(25,NTIME_MAX), LM(NTIME_MAX)
      REAL*8 LSTAR(NTIME_MAX), BLOCAL(NTIME_MAX), BMIN(NTIME_MAX)
      REAL*8 XJ(NTIME_MAX), MLT(NTIME_MAX), AMJD, PHI
      REAL*8 XIN(3), XOUT(3), FKM, FSM, FS, MLAT, FALT /0.0D0/, ALTMIN
      INTEGER*4 NTIME, KEXT, OPTIONS(5), SYSAXES, YEAR(NTIME_MAX)
      INTEGER*4 DOY(NTIME_MAX), IFAIL, KTYPLUS /3/, I

      TYPE(zgeo) SGEO

      TYPE(zxyz) CGEO

      TYPE(zlbl) MLAB

C Functions
      REAL*8 YDOYUT_TO_JULDAY
      EXTERNAL YDOYUT_TO_JULDAY

      DO I=1,NTIME
        CALL INITIALISE_FIELDS(OPTIONS, YEAR(I), DOY(I), UT(I), KEXT,
     &                         MAGINPUT(1,I))
C Convert to GEO spherical
        XIN(1) = XIN1(I)
        XIN(2) = XIN2(I)
        XIN(3) = XIN3(I)
        CALL COORD_TRANS1(SYSAXES, 1, YEAR(I), DOY(I), UT(I), XIN, XOUT)
        CGEO%X = XOUT(1)
        CGEO%Y = XOUT(2)
        CGEO%Z = XOUT(3)
        CALL UT546(CGEO, SGEO)
        CALL UL220(SGEO, 90.0D0, 1, BLOCAL(I), LM(I), FKM, FS, BMIN(I),
     &             FS, IFAIL)
        IF ((BLOCAL(I) .LE. 0.0D0) .OR. ISNAN(BLOCAL(I)))
     &    BLOCAL(I) = BADDATA
        IF ((BLOCAL(I) .GT. 0.0D0) .AND. (FKM .GT. 0.0D0))
     &    XJ(I) = FKM / SQRT(BLOCAL(I)) / 6371.2D0
        IF ((LM(I) .LE. 0.0D0) .OR. ISNAN(LM(I)))
     &    LM(I) = BADDATA
        IF ((BLOCAL(I) .GT. 0.0D0) .AND. (LM(I) .GT. 0.0D0)) THEN
          MLAB%LINV = .FALSE.
          MLAB%LBMP = .TRUE.
          MLAB%LKAUF = .FALSE.
          MLAB%LLMI = .TRUE.
          MLAB%LALP0 = .FALSE.
          MLAB%LPHI = .FALSE.
          MLAB%LTIM = .FALSE.
          MLAB%FBMP = BLOCAL(I)
          MLAB%FLMI = LM(I)
          IFAIL = 0
          CALL UD317(MLAB, FALT, KTYPLUS, ALTMIN, IFAIL)
          IF (IFAIL .GE. 0) THEN
            IFAIL = 0
            CALL UD330(PHI, LSTAR(I), IFAIL)
            IF ((LSTAR(I) .LE. 0.0D0) .OR. ISNAN(LSTAR(I)))
     &        LSTAR(I) = BADDATA
          ELSE
            LSTAR(I) = BADDATA
          END IF
        END IF
C Scale to nT for output
        IF (BLOCAL(I) .GT. 0.0D0)
     &    BLOCAL(I) = BLOCAL(I) * 1.0D5
C Scale to nT for output
        IF ((BMIN(I) .LE. 0.0D0) .OR. ISNAN(BMIN(I))) THEN
          BMIN(I) = BADDATA
        ELSE
          BMIN(I) = BMIN(I) * 1.0D5
        END IF
        IFAIL = 0
        AMJD = YDOYUT_TO_JULDAY(YEAR(I), DOY(I), UT(I)) - 2433282.5D0
        CALL UM538(SGEO, AMJD, MLT(I), MLAT, IFAIL)
        IF ((IFAIL .LT. 0) .OR. ISNAN(MLT(I)))
     &    MLT(I) = BADDATA
      END DO

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Calculate L* for a range of pitch angles
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE MAKE_LSTAR_SHELL_SPLITTING1(NTIME, NIPA, KEXT, OPTIONS,
     &                                       SYSAXES, YEAR, DOY, UT,
     &                                       XIN1, XIN2, XIN3, ALPHA,
     &                                       MAGINPUT, LM, LSTAR,
     &                                       BLOCAL, BMIN, XJ, MLT)

      IMPLICIT NONE

      INCLUDE 'structure.h'

      INCLUDE 'IRBEMAPI.inc'

      INTEGER*4 NALP
      PARAMETER (NALP=25)
      REAL*8 UT(NTIME_MAX), XIN1(NTIME_MAX), XIN2(NTIME_MAX)
      REAL*8 XIN3(NTIME_MAX), ALPHA(NALP), MAGINPUT(25,NTIME_MAX)
      REAL*8 LM(NTIME_MAX,NALP), LSTAR(NTIME_MAX,NALP)
      REAL*8 BLOCAL(NTIME_MAX,NALP), BMIN(NTIME_MAX), XJ(NTIME_MAX,NALP)
      REAL*8 MLT(NTIME_MAX), ALPHAS(NIPA), FBM(NIPA), FLM(NIPA)
      REAL*8 FKM(NIPA), FSM(NIPA), FS(NIPA), DALPHA(NIPA), FIM(NIPA)
      REAL*8 FLSTAR(NIPA), XIN(3), XOUT(3), MLAT, FALT /0.0D0/, ALTMIN
      REAL*8 AMJD, PHI
      INTEGER*4 NTIME, NIPA, KEXT, OPTIONS(5), SYSAXES, YEAR(NTIME_MAX)
      INTEGER*4 DOY(NTIME_MAX), IFAIL, NSORT(NIPA), NU, KTYPLUS /3/
      INTEGER*4 I, J

      TYPE(zgeo) SGEO

      TYPE(zxyz) CGEO

      TYPE(zlbl) MLAB

C Functions
      REAL*8 YDOYUT_TO_JULDAY
      EXTERNAL YDOYUT_TO_JULDAY

C Sort pitch angles into ascending order
      DO I=1,NIPA
        IF (ALPHA(I) .GT. 90.0D0) THEN
          DALPHA(I) = 180.0D0 - ALPHA(I)
        ELSE
          DALPHA(I) = ALPHA(I)
        END IF
      END DO
      IF (NIPA .GT. 1) THEN
        CALL ALPHA_SORT(NIPA, DALPHA, ALPHAS, NSORT, NU)
      ELSE
        ALPHAS(1) = DALPHA(1)
        NU = 1
        NSORT(1) = 1
      END IF

      DO I=1,NTIME
        CALL INITIALISE_FIELDS(OPTIONS, YEAR(I), DOY(I), UT(I), KEXT,
     &                         MAGINPUT(1,I))
C Convert to GEO spherical
        XIN(1) = XIN1(I)
        XIN(2) = XIN2(I)
        XIN(3) = XIN3(I)
        CALL COORD_TRANS1(SYSAXES, 1, YEAR(I), DOY(I), UT(I), XIN, XOUT)
        CGEO%X = XOUT(1)
        CGEO%Y = XOUT(2)
        CGEO%Z = XOUT(3)
        CALL UT546(CGEO, SGEO)
        CALL UL220(SGEO, ALPHAS, NU, FBM, FLM, FKM, FSM, BMIN(I), FS,
     &             IFAIL)
        DO J=1,NU
          IF ((FBM(J) .LE. 0.0D0) .OR. ISNAN(FBM(J)))
     &      FBM(J) = BADDATA
          IF ((FBM(J) .GT. 0.0D0) .AND. (FKM(J) .GT. 0.0D0)) THEN
            FIM(J) = FKM(J) / SQRT(FBM(J)) / 6371.2D0
          ELSE
            FIM(J) = BADDATA
          END IF
          IF ((FLM(J) .LE. 0.0D0) .OR. ISNAN(FLM(J)))
     &      FLM(J) = BADDATA
          IF ((FBM(J) .GT. 0.0D0) .AND. (FLM(J) .GT. 0.0D0)) THEN
            MLAB%LINV = .FALSE.
            MLAB%LBMP = .TRUE.
            MLAB%LKAUF = .FALSE.
            MLAB%LLMI = .TRUE.
            MLAB%LALP0 = .FALSE.
            MLAB%LPHI = .FALSE.
            MLAB%LTIM = .FALSE.
            MLAB%FBMP = FBM(J)
            MLAB%FLMI = FLM(J)
            IFAIL = 0
            CALL UD317(MLAB, FALT, KTYPLUS, ALTMIN, IFAIL)
            IF (IFAIL .GE. 0) THEN
              IFAIL = 0
              CALL UD330(PHI, FLSTAR(J), IFAIL)
              IF ((FLSTAR(J) .LE. 0.0D0) .OR. ISNAN(FLSTAR(J)))
     &          FLSTAR(J) = BADDATA
            ELSE
              FLSTAR(J) = BADDATA
            END IF
          ELSE
            FLSTAR(J) = BADDATA
          END IF
C Scale to nT for output
          IF (FBM(J) .GT. 0.0D0)
     &      FBM(J) = FBM(J) * 1.0D5
        END DO
C Scale to nT for output
        IF ((BMIN(I) .LE. 0.0D0) .OR. ISNAN(BMIN(I))) THEN
          BMIN(I) = BADDATA
        ELSE
          BMIN(I) = BMIN(I) * 1.0D5
        END IF
        IFAIL = 0
        AMJD = YDOYUT_TO_JULDAY(YEAR(I), DOY(I), UT(I)) - 2433282.5D0
        CALL UM538(SGEO, AMJD, MLT(I), MLAT, IFAIL)
        IF ((IFAIL .LT. 0) .OR. ISNAN(MLT(I)))
     &    MLT(I) = BADDATA
        DO J=1,NIPA
          LM(I,J) = FLM(NSORT(J))
          LSTAR(I,J) = FLSTAR(NSORT(J))
          BLOCAL(I,J) = FBM(NSORT(J))
          XJ(I,J) = FIM(NSORT(J))
        END DO
      END DO

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Calculate L* for pitch angle 90deg
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE LSTAR_PHI1(NTIME, WHICHINV, OPTIONS, YEAR, DOY, LSTAR,
     &                      PHI)

      IMPLICIT NONE

      INCLUDE 'structure.h'

      INCLUDE 'IRBEMAPI.inc'

      REAL*8 LSTAR(NTIME_MAX), PHI(NTIME_MAX), MAGINPUT(25), GM
      INTEGER*4 NTIME, WHICHINV, OPTIONS(5), YEAR(NTIME_MAX), TWOPI
      INTEGER*4 DOY(NTIME_MAX), I

      TYPE(zimf) mint
      TYPE(zemf) mext 
      TYPE(zsun) msun

      COMMON /UC140/ mint, mext, msun

      REAL*8 pi, deg, re, gmagmo, eclipt, geoid(3), uma(30)

      COMMON /UC160/ pi, deg, re, gmagmo, eclipt, geoid, uma

      REAL*8 prop, stepx, stpmin, umsq, upsq, uk2, uk3, epskm, epsrel
      REAL*8 stplst, xclat
      INTEGER*4 kmflg, nxstp

      COMMON /UC190/ prop, stepx, stpmin, umsq, upsq, uk2, uk3, 
     &               epskm, epsrel, stplst, xclat, kmflg, nxstp

      DO I=1,NTIME
        CALL INITIALISE_FIELDS(OPTIONS, YEAR(I), DOY(I), 0.0D0, 0,
     &                         MAGINPUT)
C Set the magnetic field moment
        IF (MOD(kmflg, 10) .EQ. 0) THEN
          GM = gmagmo
        ELSE
          GM = mint%gmmo
        END IF

C L* -> Phi
        IF (WHICHINV .EQ. 1) THEN
          IF (LSTAR(I) .EQ. BADDATA) THEN
            PHI(I) = BADDATA
          ELSE
            PHI(I) = 2.0D0 * PI * GM * RE**2 / LSTAR(I)
          END IF
C Phi -> L*
        ELSE
          IF (PHI(I) .EQ. BADDATA) THEN
            LSTAR(I) = BADDATA
          ELSE
            LSTAR(I) = 2.0D0 * PI * GM * RE**2 / PHI(I)
          END IF
        END IF
      END DO

      RETURN
      END

C
C-----------------------------------------------------------------------
C
C Sorting routine used by card players.
C Pick out the second card and put it in order with respect to the
C first; then pick out the third card and insert it into the sequence
C among the first two; and so on until the last card has been picked
C out and inserted.
      SUBROUTINE ALPHA_SORT(N, X, XS, NSORT, NU)

      IMPLICIT NONE

      INTEGER*4 N, NSORT(N), NU, NS(N), I, J
      REAL*8 X(N), XS(N), XT

      XS = X
      NS(1) = 1
      DO J=2,N
        XT = XS(J)
        NS(J) = J
        DO I=J-1,1,-1
          IF (XS(I) .LE. XT) GO TO 10
          XS(I+1) = XS(I)
          NS(I+1) = NS(I)
        END DO
        I = 0
 10     XS(I+1) = XT
        NS(I+1) = J
      END DO

      DO I=1,N
        DO J=1,N
          IF (NS(J) .EQ. I) THEN
            NSORT(I) = J
            EXIT
          END IF
        END DO
      END DO

C Pitch angles too close together can cause problems with line tracing.
C Return unique pitch angle values only (more than 0.5° apart).
      NU = 1
      NSORT(1) = 1
      DO I=2,N
        IF (ABS(XS(I)-XS(I-1)) .GT. 0.5D0) THEN
          NU = NU + 1
          XS(NU) = XS(I)
        END IF
        NSORT(NS(I)) = NU
      END DO

      RETURN
      END
