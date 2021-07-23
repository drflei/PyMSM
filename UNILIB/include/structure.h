C     INCLUDE FILE
C
C      IMPLICIT NONE
C
C     VARIABLE DEFAULT TYPE
C
C        A-H, P-Y      REAL*8
C        I, J, K, N    INTEGER*4
C        L             CHARACTER*(*), LOGICAL*1
C        M             TYPE
C
C        Z             TYPE definition
C
C
C
C
C     GLOBAL PARAMETERS
C
      INTEGER*4  nx120, nx130, nx140, nx170
C
      PARAMETER( nx120 =    500 )
      PARAMETER( nx130 =  20000 )
      PARAMETER( nx140 =     16 )
      PARAMETER( nx170 =    500 )      
C
C
C
C     TYPE DEFINITION
C
C     (X,Y,Z) Coordinate [by default, km]
C
      TYPE ZXYZ
        SEQUENCE
        REAL*8       X, Y, Z
      END TYPE ZXYZ
C
C     Geocentric coordinate [km, deg]
C
      TYPE ZGEO
        SEQUENCE
        REAL*8       RADIUS
        REAL*8       COLAT
        REAL*8       ELONG
      END TYPE ZGEO
C
C     Vector magnitude and spherical component [by default Gauss]
C
      TYPE ZVEC
        SEQUENCE
        REAL*8       DNRM
        REAL*8       RHO
        REAL*8       THETA
        REAL*8       PHI
      END TYPE ZVEC
C
C     Drift shell labels
C
      TYPE ZLBL
        SEQUENCE
        REAL*8       FINV
        REAL*8       FBMP
        REAL*8       FKAUF
        REAL*8       FLMI
        REAL*8       FALP0
        REAL*8       FPHI
        REAL*8       FTIM
        CHARACTER*1  LABEL
        LOGICAL*1    LINV, LBMP, LKAUF, LLMI, LALP0, LPHI, LTIM
      END TYPE ZLBL
C
C     Point on a field line
C
      TYPE ZPNT
        SEQUENCE
        TYPE(ZGEO) :: COORD
        TYPE(ZVEC) :: B
        REAL*8       RCURV
      END TYPE ZPNT
C
C     Elementary segment along a field line
C
      TYPE ZSEG
        SEQUENCE
        TYPE(ZPNT) :: BEG
        REAL*8       ARCL
        REAL*8       CSALP
        REAL*8       DTBND
        REAL*8       RKSTP(3)
      END TYPE ZSEG
C
C     Indices of a field line segment
C
      TYPE ZIND
        SEQUENCE
        INTEGER*4    JBEG
        INTEGER*4    JEND
        INTEGER*4    JMIRPN
        INTEGER*4    JMIRPS
      END TYPE ZIND
C
C     Description of a field line segment which is part of a drift shell
C
      TYPE ZFLN
        SEQUENCE
        TYPE(ZPNT) :: EQUAT
        TYPE(ZGEO) :: FOOTPN
        TYPE(ZGEO) :: FOOTPS
        TYPE(ZVEC) :: DRIFT
        REAL*8       DTDFT
        INTEGER*4    KWEST, KEAST
        TYPE(ZIND) :: IND
      END TYPE ZFLN
C
C     Geomagnetic field model
C
      TYPE ZIMF
        SEQUENCE
        CHARACTER*32 LABEL
        INTEGER*4    KINNER
        INTEGER*4    NORDER
        REAL*8       COEF(NX140,NX140)
        REAL*8       COLAT, ELONG
        REAL*8       GMMO
        REAL*8       TZERO, EPOCH, SAAROT
      END TYPE ZIMF
C
C     Sun position
C
      TYPE ZSUN
        SEQUENCE
        REAL*8       UTDEG
        REAL*8       GHA
        TYPE(ZXYZ) :: DIR
      END TYPE ZSUN
C
C     External magnetic field
C
      TYPE ZEMF
        SEQUENCE
        REAL*8       VDST
        REAL*8       WDENS, WVEL
        REAL*8       VKP
        REAL*8       VAL
        REAL*8       PDYN, BXIMF, BYIMF, BZIMF, STDOFF
        REAL*8       G1, G2, W1, W2, W3, W4, W5, W6
        REAL*8       TRANS(3,3)
        REAL*8       TILT
        CHARACTER*20 LABEL
        INTEGER*4    KOUTER
        INTEGER*4    IKP
        CHARACTER*4  LBLTNS
      END TYPE ZEMF
C
C     Date and Time
C
      TYPE ZDAT
        SEQUENCE
        REAL*8       SECS
        REAL*8       AMJD
        INTEGER*4    IYEAR, IMONTH, IDAY
        INTEGER*4    IHOUR, IMIN, IDUMMY
      END TYPE ZDAT
C
C     Atmospheric model
C
      TYPE ZATM
        SEQUENCE
        REAL*8       UT
        REAL*8       RZSS
        REAL*8       F107A, F107
        REAL*8       APIND(7)
        REAL*8       FKPX
        INTEGER*4    KATM
        INTEGER*4    KION
        INTEGER*4    KYEAR, KDAY
      END TYPE ZATM
C
C
C
C     COMMON LIST
C
C     COMMON /UC110/ mlbl, k1st, klmp
C
C            INTEGER*4     k1st, klmp
C            TYPE(zlbl) :: mlbl
C
C     COMMON /UC120/ nbrfl, kurfl, mfl
C
C            INTEGER*4     nbrfl, kurfl
C            TYPE(zfln) :: mfl(nx120) 
C
C     COMMON /UC130/ nbrsg, kursg, mseg
C
C            INTEGER*4     nbrsg, kursg
C            TYPE(zseg) :: mseg(nx130)
C
C     COMMON /UC140/ mint, mext, msun
C
C            TYPE(zimf) :: mint
C            TYPE(zsun) :: msun
C            TYPE(zemf) :: mext          
C
C     COMMON /UC150/matm, ntspec, nnspec, kspec, kflag
C
C            TYPE(zatm) :: matm
C            INTEGER*4     ntspec, nnspec, kspec(30)
C            INTEGER*4     kflag(50)
C                                      
C     COMMON /UC160/ pi, deg, re, gmagmo, eclipt, geoid, uma
C
C            REAL*8        pi, deg, re
C            REAL*8        gmagmo
C            REAL*8        eclipt, geoid(3)
C            REAL*8        uma(30)
C
C     COMMON /UC170/ nsg, kgp, mlab, mlin, mele
C
C            INTEGER*4     nsg, kgp
C            TYPE(zlbl) :: mlab
C            TYPE(zfln) :: mlin
C            TYPE(zseg) :: mele(nx170)
C
C     COMMON /UC190/ prop, stepx, stpmin, umsq, upsq, uk2, uk3, 
C    :               epskm, epsrel, stplst, xclat, kmflg, nxstp       
C
C            REAL*8        prop, stepx, stpmin
C            REAL*8        umsq, upsq, uk2, uk3
C            REAL*8        epskm, epsrel, stplst, xclat
C            INTEGER*4     kmflg, nxstp       
C
C     COMMON /UC192/ xrmin, xbmin, xtmin, xbmax, epslon, epsfl,
C    :               fvet, pvet, epsomeg, dltlat
C
C            REAL*8        xrmin, xbmin, xtmin, xbmax, epslon, epsfl
C            REAL*8        fvet, pvet, epsomeg, dltla                   
C
