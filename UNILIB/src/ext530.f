# 1 "ext530.f"
      subroutine ext530
      end
      SUBROUTINE MEAD(xx, yy, zz, tilt, kp, abx, aby, abz)

C CALCULATES INCREMENTAL HIGH ALTITUDE GEOMAGNETIC FIELD IN GAUSS
C xx, yy, zz   SM coordinates in Re
C tilt      in degrees
C abx, aby, abz in gauss

      IMPLICIT NONE

C INPUT:
      REAL*8  xx, yy, zz, tilt
      INTEGER*4 KP
C OUTPUT:
      REAL*8  abx, aby, abz
C LOCAL:
      REAL*8  a(7,4), b(3,4), c(7,4), x, y, z, t
      INTEGER*4 i, j, k, l, m, n 

C COEFFICIENTS OF MODEL 1 (MEAD-FAIRFIELD 1975)

      DATA ((A(I,J),I=1,7),J=1,4)/17.93D0,-5.79D0,2.98D0,-2.57D0,
     1-0.30D0,-1.47D0,1.05D0,21.79D0,-7.03D0,3.02D0,-2.99D0,-0.62D0,
     2-1.22D0,0.95D0,33.16D0,-6.39D0,4.30D0,-3.25D0,-0.44D0,-1.27D0,
     30.45D0,39.48D0,-2.91D0,5.17D0,-3.86D0,-1.04D0,-1.29D0,-1.14D0/,
     4((B(K,L),K=1,3),L=1,4)/-10.11D0,-1.98D0,0.09D0,-11.84D0,-2.57D0,
     5-0.28D0,-16.54D0,-3.08D0,0.22D0,-19.10D0,-3.50D0,0.23D0/,
     6((C(M,N),M=1,7),N=1,4)/-9.41D0,15.07D0,13.16D0,8.36D0,7.95D0,
     74.55D0,0.51D0,-11.96D0,17.87D0,15.88D0,9.77D0,9.43D0,5.57D0,
     81.53D0,-19.88D0,20.23D0,22.72D0,13.23D0,11.46D0,6.33D0,0.67D0,
     9-22.90D0,22.70D0,26.50D0,15.54D0,11.00D0,7.36D0,1.85D0/

C ROTATE AXES ABOUT Z BY 4 DEGREES
      X = 0.99756405D0 * XX - 0.06975647D0 * YY
      Y = 0.06975647D0 * XX + 0.99756405D0 * YY
C SCALE COORDINATES AND TILT ANGLE
      x = x / 10.0D0
      y = y / 10.0D0
      z = zz / 10.0D0
      T = TILT / 10.0D0
C CALCULATE EXTERNAL FIELD COMPONENTS
      ABX = A(1,KP) * Z + A(2,KP) * X * Z + T * (A(3,KP) + A(4,KP) * X
     &      + A(5,KP) * X * X + A(6,KP) * Y * Y + A(7,KP) * Z * Z)
      ABY = B(1,KP) * Y * Z + T * (B(2,KP) * Y + B(3,KP) * X * Y)
      ABZ = C(1,KP) + C(2,KP) * X + C(3,KP) * X * X + C(4,KP) * Y * Y
     &      + C(5,KP) * Z * Z + T * (C(6,KP) * Z + C(7,KP) * X * Z)

      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE TSY87S(IOPT, PS, X, Y, Z, BX, BY, BZ)
C-----------------------------------------------------------------------
C  Computes GSM components of the magnetic field of extraterrestrial 
C  current systems up to geocentric radial distances about 30 Re. 
C  Corresponds to magnetospheric field model (N.A.Tsyganenko, 1986). 
C  Based on IMP-A,C,D,E,F,G,H,I,J (years 1966-1980) and HEOS-1,2 
C  (1969-1974) satellite merged data set.
C
C  INPUT:   IOPT  specifies model version, according to ground 
C           disturbance level, as specified below
C        PS geodipole tilt angle in radians
C        X,Y,Z    GSM coordinates of the point in Earth radii
C
C  OUTPUT:  BX,BY,BZ    GSM components of the external
C              magnetic field in nanotesla
C
C  IOPT    1     2      3       4      5        6         7      8
C   KP    0,0+  1-,1   1+,2-   2,2+  3-,3,3+   4-,4,4+   >=5-   >=5+
C-----------------------------------------------------------------------
      
      IMPLICIT NONE

C Input
      REAL*8 x, y, z, ps
      INTEGER*4 IOPT
C Output
      REAL*8 bx, by, bz
C Local
      REAL*8 GA(24,8), PA(24), HPI, FC, RT
      REAL*8 C1, RRC2, DSTR, XN, RH, X1, DY, B0, B1, XN21, SPS, CPS, RPS
      REAL*8 ZS, ZP, ZM, FY, XNX, XNX2, XC1, XC12, B20, B2P, B2M, B, BP
      REAL*8 BM, XA1, XAP1, XAM1, XNA, XNAP, XNAM, XLN1, XLNP1, XLNM1
      REAL*8 ALN, S0, S0P, S0M, S1, S1P, S1M, G1, G1P, G1M, EX, Y2, Z2
      REAL*8 YZ, XSM, ZSM, RR, RR2, ZN, BRSM, BZSM, BY1, BXSM
      INTEGER*4 I, IP
      SAVE IP,PA,C1,RRC2,DSTR,XN,RH,X1,DY,B0,B1,XN21
C
C    &12.72,-.00867,-.001953,-.3437,-.002903,-.000999,18.41,-270.3,
C
C errata 1.29.90 D.P. Stern
C
      DATA GA /1.126d0,26.66d0,-.077d0,-.06102d0,-.06197d0,-2.048d0,
     :.00327d0,.008473d0,12.72d0,-.00867d0,-.01953d0,-.3437d0,
     :-.002903d0,-.000999d0,18.41d0,-270.3d0,-25.94d0,5.21d0,-6.2d0,
     :2.29d0,11.96d0,8.315d0,44.22d0,11.15d0,1.403d0,29.24d0,-.0693d0,
     :-.0864d0,-.07202d0,-2.068d0,.00286d0,.007438d0,16.37d0,-.02705d0,
     :-.0281d0,-.604d0,-.002256d0,.000152d0,20.2d0,-140.1d0,-29.65d0,
     :5.62d0,-5.52d0,2.02d0,14.66d0,8.06d0,27.76d0,10.94d0,1.589d0,
     :31.07d0,-.06527d0,-.07447d0,-.07632d0,-2.413d0,.002719d0,
     :.01098d0,16.2d0,-.02355d0,-.03475d0,-.4377d0,-.002169d0,
     :-.001383d0,18.7d0,-292.6d0,-35.25d0,5.29d0,-5.18d0,2.21d0,
     :14.03d0,7.66d0,17.56d0,10.9d0,1.699d0,36.28d0,-.07514d0,-.1448d0,
     :-.08049d0,-2.209d0,.000919d0,.01084d0,17.38d0,-.03516d0,
     :-.03886d0,-1.169d0,.004239d0,.000881d0,21.79d0,-162.d0,-41.87d0,
     :5.15d0,-3.62d0,2.35d0,17.26d0,7.61d0,17.99d0,10.74d0,2.141d0,
     :41.51d0,-.1518d0,-.1857d0,-.1015d0,-2.929d0,.004584d0,.01589d0,
     :18.29d0,-.02514d0,-.05927d0,-1.336d0,.00185d0,.001066d0,21.31d0,
     :-358.8d0,-47.91d0,5.13d0,-3.74d0,2.07d0,17.23d0,6.33d0,32.51d0,
     :9.73d0,2.252d0,39.35d0,-.04525d0,-.2062d0,-.1491d0,-3.059d0,
     :-.000183d0,.02614d0,15.48d0,-.02144d0,-.06608d0,-1.855d0,
     :.006199d0,-.00013d0,23.91d0,-161.d0,-51.48d0,4.61d0,-3.32d0,
     :1.68d0,15.22d0,6.68d0,.6765d0,8.007d0,2.773d0,40.95d0,.00667d0,
     :-.133d0,-.1304d0,-5.187d0,.004623d0,.03651d0,20.d0,-.03765d0,
     :-.09066d0,.5838d0,-.01462d0,-.007189d0,24.87d0,-186.1d0,-74.81d0,
     :4.57d0,-4.03d0,1.7d0,12.15d0,6.87d0,-1.746d0,8.9d0,2.919d0,
     :34.96d0,2*0.d0,-.1609d0,-5.077d0,2*0.d0,22.1d0,-.05915d0,
     :-.1051d0,.6321d0,2*0.d0,28.11d0,-330.1d0,-86.82d0,4.d0,-3.d0,
     :1.73d0,12.56d0,5.11d0,4.d0,7.866d0/

      DATA IP, FC, RT /100, 0.3183099031D0, 30.0D0/

      HPI = 2.0D0 * DATAN(1.0D0)

      IF (IOPT .NE. IP) THEN
        IP = IOPT
        DO I=1,24
          PA(I) = GA(I,IP)         
        END DO
        C1 = PA(20)**2
        RRC2 = PA(18)**2
        DSTR = PA(17) / RRC2 * 4.0D0
        XN = PA(19)
        RH = PA(22)
        X1 = PA(23)
        DY = PA(21)
        B0 = PA(15)
        B1 = PA(16)
        XN21 = (XN-X1)**2
      END IF
      SPS = DSIN(PS)
      CPS = DCOS(PS)
      RPS = RH * SPS
C
C   COMPUTATION BEGINS HERE IF PARAMETER IOPT REMAINED UNCHANGED AFTER THE 
C   PRECEDING CALL OF THIS SUBROUTINE.
C
      ZS = Z - RPS
      ZP = Z - RT
      ZM = Z + RT
      FY = FC / (1.0D0+(Y/DY)**2)
      XNX = XN - X
      XNX2 = XNX**2
      XC1 = X - X1
      XC12 = XC1**2
      B20 = ZS**2 + C1
      B2P = ZP**2 + C1
      B2M = ZM**2 + C1
      B = DSQRT(B20)
      BP = DSQRT(B2P)
      BM = DSQRT(B2M)
      XA1 = XC12 + B20
      XAP1 = XC12 + B2P
      XAM1 = XC12 + B2M
      XNA = XNX2 + B20
      XNAP = XNX2 + B2P
      XNAM = XNX2 + B2M
      XLN1 = DLOG(XN21/XNA)
      XLNP1 = DLOG(XN21/XNAP)
      XLNM1 = DLOG(XN21/XNAM)
      ALN = 0.25D0 * (XLNP1+XLNM1-2.0D0*XLN1)
      S0 = (DATAN(XNX/B)+HPI) / B
      S0P = (DATAN(XNX/BP)+HPI) / BP
      S0M = (DATAN(XNX/BM)+HPI) / BM
      S1 = (XLN1*0.5D0+XC1*S0) / XA1
      S1P = (XLNP1*0.5D0+XC1*S0P) / XAP1
      S1M = (XLNM1*0.5D0+XC1*S0M) / XAM1
      G1 = (B20*S0-0.5D0*XC1*XLN1) / XA1
      G1P = (B2P*S0P-0.5D0*XC1*XLNP1) / XAP1
      G1M = (B2M*S0M-0.5D0*XC1*XLNM1) / XAM1
      BX = FY * (B0*(ZS*S0-0.5D0*(ZP*S0P+ZM*S0M))
     &     +B1*(ZS*S1-0.5D0*(ZP*S1P+ZM*S1M)))
      BY = 0.0D0
      BZ = FY * (B0*ALN+B1*(G1-0.5D0*(G1P+G1M)))
C
C    CALCULATION OF THE MAGNETOTAIL CURRENT CONTRIBUTION IS FINISHED
C
      EX = DEXP(X/PA(24))
      Y2 = Y**2
      Z2 = Z**2
      YZ = Y * Z
      BX = BX + EX * (CPS*PA(1)*Z+SPS*(PA(2)+PA(3)*YZ+PA(4)*Z2))
      BY = EX * (CPS*PA(5)*YZ+SPS*Y*(PA(6)+PA(7)*Y2+PA(8)*Z2))
      BZ = BZ + EX * (CPS*(PA(9)+PA(10)*Y2+PA(11)*Z2)+SPS*Z*(PA(12)
     &     +PA(13)*Y2+PA(14)*Z2))
C
C   DCF AND FAC CONTRIBUTION HAS BEEN ADDED TO BX, BY, AND BZ
C
      XSM = X * CPS - Z * SPS
      ZSM = X * SPS + Z * CPS
      Z2 = ZSM**2
      RR2 = XSM**2 + Y**2
      RR = DSQRT(RR2)
      ZN = ((RR2+Z2)/RRC2+4.0D0)**2.5D0
      BRSM = DSTR * 3.0D0 * ZSM / ZN
      BZSM = DSTR * (2.0D0*Z2-RR2+8.0D0*RRC2) / ZN
      BY1 = BRSM * Y
      BXSM = BRSM * XSM
      BX = BX + BXSM * CPS + BZSM * SPS
      BZ = BZ - BXSM * SPS + BZSM * CPS
      BY = BY + BY1
C
C   RING CURRENT FIELD HAS BEEN TAKEN INTO ACCOUNT
C
      RETURN
      END
C
      SUBROUTINE TSY87L(IOPT, PS, X, Y, Z, BX, BY, BZ)
C-----------------------------------------------------------------------
C  'LONG' version of the magnetospheric magnetic field model.
C  Computes GSM components of the magnetic field of extraterrestrial
C  current systems up to geocentric radial distances about 70 Re.
C  Corresponds to magnetospheric field model (N.A.Tsyganenko, 1986)
C  based on IMP-A,C,D,E,F,G,H,I,J (years 1966-1974) and HEOS-1,2 
C  (1969-1974) satellite merged data set.
C
C  INPUT:   IOPT  specifies model version, according to ground 
C           disturbance level, as specified below
C        PS geodipole tilt angle in radians
C        X,Y,Z    GSM coordinates of the point in Earth radii
C
C  OUTPUT:  BX,BY,BZ    GSM components of the external
C              magnetic field in nanotesla
C
C  IOPT     1      2         3        4         5        6
C   KP    0,0+  1-,1,1+   2-,2,2+  3-,3,3+   4-,4,4+   >=5-  
C-----------------------------------------------------------------------

      IMPLICIT NONE

C Input
      REAL*8 ps, x, y, z
      INTEGER*4 IOPT
C Output
      REAL*8 bx, by, bz
C Local
      REAL*8 GA(32,6), PA(32), HPI, FC, RT, X1, X2
      REAL*8 C1, RRC2, DSTR, XN, RH, DY, B0, B1, B2, XN21, XN2, XNR
      REAL*8 XN22, ADLN, SPS, CPS, RPS, ZS, ZP, ZM, FY, XNX, XNX2
      REAL*8 XC1, XC2, XC22, XR2, XC12, B20, B2P, B2M, B, BP, BM, XA1
      REAL*8 XAP1, XAM1, XA2, XAP2, XAM2, XNA, XNAP, XNAM, F, FP, FM
      REAL*8 XLN1, XLNP1, XLNM1, XLN2, XLNP2, XLNM2, ALN, S0, S0P, S0M
      REAL*8 S1, S1P, S1M, S2, S2P, S2M, G1, G1P, G1M, G2, G2P, G2M
      REAL*8 EX1, EX2, Y2, Z2, YZ, XSM, ZSM, RR, RR2, ZN, BRSM, BZSM
      REAL*8 BY1, BXSM
      INTEGER*4 I, IP
      SAVE IP,PA,C1,RRC2,DSTR,XN,RH,X1,DY,B0,B1,XN21,XN2,XNR,XN22,ADLN
C
C    &-.1058,-3.221,-.00114,-.02166,-30.43,.04049,.05464,.008884,42.,
C
C errata 1.29.90 D.P. Stern
C
      DATA GA /-.09673d0,-10.63d0,1.21d0,34.57d0,-.04502d0,-.06553d0,
     :-.02952d0,.3852d0,-.03665d0,-2.084d0,.001795d0,.00638d0,-23.49d0,
     :.06082d0,.01642d0,-.02137d0,32.21d0,-.04373d0,-.02311d0,-.2832d0,
     :-.002303d0,-.000631d0,-6.397d0,-967.d0,-8650.d0,-20.55d0,5.18d0,
     :-2.796d0,2.715d0,13.58d0,8.038d0,29.21d0,-.485d0,-12.84d0,
     :1.856d0,40.06d0,-.0294d0,-.09071d0,-.02993d0,.5465d0,-.04928d0,
     :-2.453d0,.001587d0,.007402d0,-29.41d0,.08101d0,.02322d0,-.1091d0,
     :40.75d0,-.07995d0,-.03859d0,-.2755d0,-.002759d0,-.000408d0,
     :-6.189d0,-957.8d0,-7246.d0,-25.51d0,5.207d0,-4.184d0,2.641d0,
     :16.56d0,7.795d0,29.36d0,-1.132d0,-18.05d0,2.625d0,48.55d0,
     :-.004868d0,-.1087d0,-.03824d0,.8514d0,-.0522d0,-2.881d0,
     :-.000295d0,.009055d0,-29.48d0,.06394d0,.03864d0,-.2288d0,41.77d0,
     :-.05849d0,-.06443d0,-.4683d0,.001222d0,-.000519d0,-3.696d0,
     :-991.1d0,-6955.d0,-31.43d0,4.878d0,-3.151d0,3.277d0,19.19d0,
     :7.248d0,28.99d0,-1.003d0,-16.98d0,3.14d0,52.81d0,-.08625d0,
     :-.1478d0,-.03501d0,.55d0,-.07778d0,-2.97d0,.002086d0,.01275d0,
     :-26.79d0,.06328d0,.03622d0,.08345d0,39.72d0,-.06009d0,-.07825d0,
     :-.9698d0,.000178d0,-.000573d0,-.9328d0,-872.5d0,-5851.d0,
     :-39.68d0,4.902d0,-3.848d0,2.79d0,20.91d0,6.193d0,26.81d0,
     :-1.539d0,-14.29d0,3.479d0,53.36d0,-.004201d0,-.2043d0,-.03932d0,
     :.6409d0,-.1058d0,-3.221d0,-.00114d0,.02166d0,-30.43d0,.04049d0,
     :.05464d0,.008884d0,42.d0,-.01035d0,-.1053d0,-1.63d0,.003802d0,
     :-.001029d0,4.204d0,-665.6d0,-1011.d0,-43.49d0,4.514d0,-2.948d0,
     :2.99d0,21.59d0,6.005d0,22.d0,-2.581d0,-7.726d0,5.045d0,53.31d0,
     :.02262d0,-.1972d0,-.01981d0,.428d0,-.1055d0,-5.075d0,.002762d0,
     :.03277d0,-27.35d0,.04986d0,.06119d0,-.1211d0,47.48d0,-.0502d0,
     :-.1477d0,.838d0,-.01008d0,-.0057d0,9.231d0,-674.3d0,-900.d0,
     :-74.43d0,4.658d0,-3.245d0,3.39d0,21.8d0,5.62d0,25.17d0/

      DATA IP, FC, RT, X1, X2 /100, 0.3183099031D0, 30.0D0, 4.0D0,5.0D0/

      HPI = 2.0D0 * DATAN(1.0D0)

      IF (IOPT .NE. IP) THEN
        IP = IOPT
        DO I=1,32
          PA(I) = GA(I,IP)
        END DO
        C1 = PA(29)**2
        RRC2 = PA(27)**2
        DSTR = PA(26) / RRC2 * 4.0D0
        XN = PA(28)
        RH = PA(31)
        DY = PA(30)
        B0 = PA(23)
        B1 = PA(24)
        B2 = PA(25)
        XN21 = (XN-X1)**2
        XN2 = XN - X2
        XNR = 1.0 / XN2
        XN22 = XN2**2
        ADLN = DLOG(XN22/XN21)
      END IF
      SPS = DSIN(PS)
      CPS = DCOS(PS)
      RPS = RH * SPS
C
C   COMPUTATION BEGINS HERE IF PARAMETER IOPT REMAINED UNCHANGED AFTER THE 
C   PRECEDING CALL OF THIS SUBROUTINE.
C
      ZS = Z - RPS
      ZP = Z - RT
      ZM = Z + RT
      FY = FC / (1.0D0+(Y/DY)**2.0D0)
      XNX = XN - X
      XNX2 = XNX**2.0D0
      XC1 = X - X1
      XC2 = X - X2
      XC22 = XC2**2.0D0
      XR2 = XC2 * XNR
      XC12 = XC1**2.0D0
      B20 = ZS**2 + C1
      B2P = ZP**2 + C1
      B2M = ZM**2 + C1
      B = DSQRT(B20)
      BP = DSQRT(B2P)
      BM = DSQRT(B2M)
      XA1 = XC12 + B20
      XAP1 = XC12 + B2P
      XAM1 = XC12 + B2M
      XA2 = 1.0D0 / (XC22+B20)
      XAP2 = 1.0D0 / (XC22+B2P)
      XAM2 = 1.0D0 / (XC22+B2M)
      XNA = XNX2 + B20
      XNAP = XNX2 + B2P
      XNAM = XNX2 + B2M
      F = B20 - XC22
      FP = B2P - XC22
      FM = B2M - XC22
      XLN1 = DLOG(XN21/XNA)
      XLNP1 = DLOG(XN21/XNAP)
      XLNM1 = DLOG(XN21/XNAM)
      XLN2 = XLN1 + ADLN
      XLNP2 = XLNP1 + ADLN
      XLNM2 = XLNM1 + ADLN
      ALN = 0.25D0 * (XLNP1+XLNM1-2.0D0*XLN1)
      S0 = (DATAN(XNX/B)+HPI) / B
      S0P = (DATAN(XNX/BP)+HPI) / BP
      S0M = (DATAN(XNX/BM)+HPI) / BM
      S1 = (XLN1*0.5D0+XC1*S0) / XA1
      S1P = (XLNP1*0.5D0+XC1*S0P) / XAP1
      S1M = (XLNM1*0.5D0+XC1*S0M) / XAM1
      S2 = (XC2*XA2*XLN2-XNR-F*XA2*S0) * XA2
      S2P = (XC2*XAP2*XLNP2-XNR-FP*XAP2*S0P) * XAP2
      S2M = (XC2*XAM2*XLNM2-XNR-FM*XAM2*S0M) * XAM2
      G1 = (B20*S0-0.5D0*XC1*XLN1) / XA1
      G1P = (B2P*S0P-0.5D0*XC1*XLNP1) / XAP1
      G1M = (B2M*S0M-0.5D0*XC1*XLNM1) / XAM1
      G2 = ((0.5D0*F*XLN2+2.0D0*S0*B20*XC2)*XA2+XR2) * XA2
      G2P = ((0.5D0*FP*XLNP2+2.0D0*S0P*B2P*XC2)*XAP2+XR2) * XAP2
      G2M = ((0.5D0*FM*XLNM2+2.0D0*S0M*B2M*XC2)*XAM2+XR2) * XAM2
      BX = FY * (B0*(ZS*S0-0.5D0*(ZP*S0P+ZM*S0M))+B1*
     &     (ZS*S1-0.5D0*(ZP*S1P+ZM*S1M))+B2*
     &     (ZS*S2-0.5D0*(ZP*S2P+ZM*S2M)))
      BY = 0.0D0
      BZ = FY * (B0*ALN+B1*(G1-0.5D0*(G1P+G1M))+B2*(G2-0.5D0*(G2P+G2M)))
C
C    CALCULATION OF THE MAGNETOTAIL CURRENT CONTRIBUTION IS FINISHED.
C
      EX1 = DEXP(X/PA(32))
      EX2 = EX1**2
      Y2 = Y**2
      Z2 = Z**2
      YZ = Y * Z
      BX = BX + (EX1*PA(1)+EX2*PA(3)) * Z * CPS +
     &     (EX1*PA(2)+EX2*(PA(4)+PA(5)*Y2+PA(6)*Z2)) * SPS
      BY = (EX1*PA(7)+EX2*PA(9)) *YZ * CPS +
     &     (EX1*PA(8)+EX2*(PA(10)+PA(11)*Y2+PA(12)*Z2)) * Y * SPS
      BZ = BZ + (EX1*(PA(13)+PA(14)*Y2+PA(15)*Z2)+EX2*(PA(17)+PA(18)*Y2+
     &     PA(19)*Z2))*CPS+(EX1*PA(16)+EX2*(PA(20)+PA(21)*Y2+PA(22)*Z2))
     &     * Z * SPS
C
C   DCF AND FAC CONTRIBUTION HAS BEEN ADDED TO BX, BY, AND BZ.
C
      XSM = X * CPS - Z * SPS
      ZSM = X * SPS + Z * CPS
      Z2 = ZSM**2
      RR2 = XSM**2 + Y**2
      RR = DSQRT(RR2)
      ZN = ((RR2+Z2)/RRC2+4.0D0)**2.5D0
      BRSM = DSTR * 3.0D0 * ZSM / ZN
      BZSM = DSTR * (2.0D0*Z2-RR2+8.0D0*RRC2) / ZN
      BY1 = BRSM * Y
      BXSM = BRSM * XSM
      BX = BX + BXSM * CPS + BZSM * SPS
      BZ = BZ - BXSM * SPS + BZSM * CPS
      BY = BY + BY1
C
C   RING CURRENT FIELD HAS BEEN TAKEN INTO ACCOUNT.
C
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE T89C(IOPT,PARMOD,PS,X,Y,Z,BX,BY,BZ)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C
C
C   COMPUTES GSM COMPONENTS OF THE MAGNETIC FIELD PRODUCED BY EXTRA-
C   TERRESTRIAL CURRENT SYSTEMS IN THE GEOMAGNETOSPHERE. THE MODEL IS
C   VALID UP TO GEOCENTRIC DISTANCES OF 70 RE AND IS BASED ON THE MER-
C   GED IMP-A,C,D,E,F,G,H,I,J (1966-1974), HEOS-1 AND -2 (1969-1974),
C   AND ISEE-1 AND -2  SPACECRAFT DATA SET.
C
C   THIS IS A MODIFIED VERSION (T89c), WHICH REPLACED THE ORIGINAL ONE
C     IN 1992 AND DIFFERS FROM IT IN THE FOLLOWING:
C
C   (1)  ISEE-1,2 DATA WERE ADDED TO THE ORIGINAL IMP-HEOS DATASET
C   (2)  TWO TERMS WERE ADDED TO THE ORIGINAL TAIL FIELD MODES, ALLOWING
C          A MODULATION OF THE CURRENT BY THE GEODIPOLE TILT ANGLE
C
C
C  REFERENCE FOR THE ORIGINAL MODEL: N.A. TSYGANENKO, A MAGNETOSPHERIC MAGNETIC
C       FIELD MODEL WITH A WARPED TAIL CURRENT SHEET: PLANET.SPACE SCI., V.37,
C         PP.5-20, 1989.
C
C----INPUT PARAMETERS: IOPT - SPECIFIES THE GROUND DISTURBANCE LEVEL:
C
C   IOPT= 1       2        3        4        5        6      7
C                  CORRESPOND TO:
C    KP= 0,0+  1-,1,1+  2-,2,2+  3-,3,3+  4-,4,4+  5-,5,5+  > =6-
C
C    PS - GEODIPOLE TILT ANGLE IN RADIANS
C    X, Y, Z  - GSM COORDINATES OF THE POINT IN EARTH RADII
C
C----OUTPUT PARAMETERS: BX,BY,BZ - GSM COMPONENTS OF THE MODEL MAGNETIC
C                        FIELD IN NANOTESLAS
c
c   THE PARAMETER PARMOD(10) IS A DUMMY ARRAY.  IT IS NOT USED IN THIS
C        SUBROUTINE AND IS PROVIDED JUST FOR MAKING IT COMPATIBLE WITH THE
C           NEW VERSION (4/16/96) OF THE GEOPACK SOFTWARE.
C
C   THIS RELEASE OF T89C IS DATED  FEB 12, 1996;
C--------------------------------------------------------------------------
C
C
C              AUTHOR:     NIKOLAI A. TSYGANENKO
C                          HSTX CORP./NASA GSFC
C
       DIMENSION XI(4),F(3),DER(3,30),PARAM(30,7),A(30),PARMOD(10)
c-mk-c       DOUBLE PRECISION F,DER
       real*8 F,DER
        DATA PARAM/-116.53,-10719.,42.375,59.753,-11363.,1.7844,30.268,
     * -0.35372E-01,-0.66832E-01,0.16456E-01,-1.3024,0.16529E-02,
     * 0.20293E-02,20.289,-0.25203E-01,224.91,-9234.8,22.788,7.8813,
     * 1.8362,-0.27228,8.8184,2.8714,14.468,32.177,0.01,0.0,
     * 7.0459,4.0,20.0,-55.553,-13198.,60.647,61.072,-16064.,
     * 2.2534,34.407,-0.38887E-01,-0.94571E-01,0.27154E-01,-1.3901,
     * 0.13460E-02,0.13238E-02,23.005,-0.30565E-01,55.047,-3875.7,
     * 20.178,7.9693,1.4575,0.89471,9.4039,3.5215,14.474,36.555,
     * 0.01,0.0,7.0787,4.0,20.0,-101.34,-13480.,111.35,12.386,-24699.,
     * 2.6459,38.948,-0.34080E-01,-0.12404,0.29702E-01,-1.4052,
     * 0.12103E-02,0.16381E-02,24.49,-0.37705E-01,-298.32,4400.9,18.692,
     * 7.9064,1.3047,2.4541,9.7012,7.1624,14.288,33.822,0.01,0.0,6.7442,
     * 4.0,20.0,-181.69,-12320.,173.79,-96.664,-39051.,3.2633,44.968,
     * -0.46377E-01,-0.16686,0.048298,-1.5473,0.10277E-02,0.31632E-02,
     * 27.341,-0.50655E-01,-514.10,12482.,16.257,8.5834,1.0194,3.6148,
     * 8.6042,5.5057,13.778,32.373,0.01,0.0,7.3195,4.0,20.0,-436.54,
     * -9001.0,323.66,-410.08,-50340.,3.9932,58.524,-0.38519E-01,
     * -0.26822,0.74528E-01,-1.4268,-0.10985E-02,0.96613E-02,27.557,
     * -0.56522E-01,-867.03,20652.,14.101,8.3501,0.72996,3.8149,9.2908,
     *  6.4674,13.729,28.353,0.01,0.0,7.4237,4.0,20.0,-707.77,-4471.9,
     * 432.81,-435.51,-60400.,4.6229,68.178,-0.88245E-01,-0.21002,
     * 0.11846,-2.6711,0.22305E-02,0.10910E-01,27.547,-0.54080E-01,
     * -424.23,1100.2,13.954,7.5337,0.89714,3.7813,8.2945,5.174,14.213,
     * 25.237,0.01,0.0,7.0037,4.0,20.0,-1190.4,2749.9,742.56,-1110.3,
     * -77193.,7.6727,102.05,-0.96015E-01,-0.74507,0.11214,-1.3614,
     * 0.15157E-02,0.22283E-01,23.164,-0.74146E-01,-2219.1,48253.,
     * 12.714,7.6777,0.57138,2.9633,9.3909,9.7263,11.123,21.558,0.01,
     * 0.0,4.4518,4.0,20.0/

        DATA IOP/10/
C
         id=0
         IF (IOP.NE.IOPT) THEN
C
            ID=1
            IOP=IOPT
            DO 1 I=1,30
   1        A(I)=PARAM(I,IOPT)
C
         ENDIF
C
        XI(1)=X
        XI(2)=Y
        XI(3)=Z
        XI(4)=PS
         CALL T89(ID,A,XI,F,DER)
          IF (ID.EQ.1) ID=2
        BX=F(1)
        BY=F(2)
        BZ=F(3)
        RETURN
       END
C-------------------------------------------------------------------
C
          SUBROUTINE  T89 (ID, A, XI, F, DER)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C
C        ***  N.A. Tsyganenko ***  8-10.12.1991  ***
C
C      Calculates dependent model variables and their deriva-
C  tives for given independent variables and model parame-
C  ters.  Specifies model functions with free parameters which
C  must be determined by means of least squares fits (RMS
C  minimization procedure).
C
C      Description of parameters:
C
C  ID  - number of the data point in a set (initial assignments are performed
c        only for ID=1, saving thus CPU time)
C  A   - input vector containing model parameters;
C  XI  - input vector containing independent variables;
C  F   - output double precision vector containing
C        calculated values of dependent variables;
C  DER   - output double precision vector containing
C        calculated values for derivatives of dependent
C        variables with respect to model parameters;
C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
C      T89 represents external magnetospheric magnetic field
C  in Cartesian SOLAR MAGNETOSPHERIC coordinates (Tsyganenko N.A.,
C  Planet. Space Sci., 1989, v.37, p.5-20; the "T89 model" with the warped
c  tail current sheet) + A MODIFICATION ADDED IN APRIL 1992 (SEE BELOW)
C
C      Model formulas for the magnetic field components contain in total
c  30 free parameters (17 linear and 13 nonlinear parameters).
C      First 2 independent linear parameters A(1)-A(2) correspond to contribu-
c  tion from the tail current system, then follow A(3) and A(4) which are the
c  amplitudes of symmetric and antisymmetric terms in the contribution from
c  the closure currents; A(5) is the ring current amplitude. Then follow the
c coefficients A(6)-A(15) which define Chapman-Ferraro+Birkeland current field.
c    The coefficients c16-c19  (see Formula 20 in the original paper),
c   due to DivB=0 condition, are expressed through A(6)-A(15) and hence are not
c    independent ones.
c  A(16) AND A(17) CORRESPOND TO THE TERMS WHICH YIELD THE TILT ANGLE DEPEN-
C    DENCE OF THE TAIL CURRENT INTENSITY (ADDED ON APRIL 9, 1992)
C
C      Nonlinear parameters:
C
C    A(18) : DX - Characteristic scale of the Chapman-Ferraro field along the
c        X-axis
C    A(19) : ADR (aRC) - Characteristic radius of the ring current
c    A(20) : D0 - Basic half-thickness of the tail current sheet
C    A(21) : DD (GamRC)- defines rate of thickening of the ring current, as
c             we go from night- to dayside
C    A(22) : Rc - an analog of "hinging distance" entering formula (11)
C    A(23) : G - amplitude of tail current warping in the Y-direction
C    A(24) : aT - Characteristic radius of the tail current
c    A(25) : Dy - characteristic scale distance in the Y direction entering
c                 in W(x,y) in (13)
c    A(26) : Delta - defines the rate of thickening of the tail current sheet
c                 in the Y-direction (in T89 it was fixed at 0.01)
c    A(27) : Q - this parameter was fixed at 0 in the final version of T89;
c              initially it was introduced for making Dy to depend on X
c    A(28) : Sx (Xo) - enters in W(x,y) ; see (13)
c    A(29) : Gam (GamT) - enters in DT in (13) and defines rate of tail sheet
c              thickening on going from night to dayside; in T89 fixed at 4.0
c    A(30) : Dyc - the Dy parameter for closure current system; in T89 fixed
c               at 20.0
c  - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
C
         REAL*8  A(*), XI(*)
C
         DIMENSION  F(3), DER(3,30)
C
         INTEGER*4  ID, I, L
         save
         DATA A02,XLW2,YN,RPI,RT/25.D0,170.D0,30.D0,0.31830989D0,30.D0/
         DATA XD,XLD2/0.D0,40.D0/
C
C   The last four quantities define variation of tail sheet thickness along X
C
      DATA SXC,XLWC2/4.D0,50.D0/
C
C   The two quantities belong to the function WC which confines tail closure
c    current in X- and Y- direction
C
      DATA DXL/20.D0/
C
C
         IF (ID.NE.1)  GOTO  3
           DO  2  I = 1, 30
             DO  1  L = 1, 3
  1            DER(L,I) = 0.0D0
  2        CONTINUE
C
       DYC=A(30)
       DYC2=DYC**2
       DX=A(18)
       HA02=0.5D0*A02
       RDX2M=-1.D0/DX**2
       RDX2=-RDX2M
       RDYC2=1.D0/DYC2
       HLWC2M=-0.5D0*XLWC2
       DRDYC2=-2.D0*RDYC2
       DRDYC3=2.D0*RDYC2*DSQRT(RDYC2)
       HXLW2M=-0.5D0*XLW2
       ADR=A(19)
       D0=A(20)
       DD=A(21)
       RC=A(22)
       G=A(23)
       AT=A(24)
       DT=D0
       DEL=A(26)
       P=A(25)
       Q=A(27)
       SX=A(28)
       GAM=A(29)
       HXLD2M=-0.5D0*XLD2
       ADSL=0.D0
       XGHS=0.D0
       H=0.D0
       HS=0.D0
       GAMH=0.D0
       W1=-0.5D0/DX
       DBLDEL=2.D0*DEL
       W2=W1*2.D0
       W4=-1.D0/3.D0
       W3=W4/DX
       W5=-0.5D0
       W6=-3.D0
       AK1=A(1)
       AK2=A(2)
       AK3=A(3)
       AK4=A(4)
       AK5=A(5)
       AK6=A(6)
       AK7=A(7)
       AK8=A(8)
       AK9=A(9)
       AK10=A(10)
       AK11=A(11)
       AK12=A(12)
       AK13=A(13)
       AK14=A(14)
       AK15=A(15)
        AK16=A(16)
        AK17=A(17)
       SXA=0.D0
       SYA=0.D0
       SZA=0.D0
       AK610=AK6*W1+AK10*W5
       AK711=AK7*W2-AK11
       AK812=AK8*W2+AK12*W6
       AK913=AK9*W3+AK13*W4
       RDXL=1.D0/DXL
       HRDXL=0.5D0*RDXL
       A6H=AK6*0.5D0
       A9T=AK9/3.D0
       YNP=RPI/YN*0.5D0
       YND=2.D0*YN
C
  3      CONTINUE
C
           X  = XI(1)
           Y  = XI(2)
           Z  = XI(3)
           TILT=XI(4)
               TLT2=TILT**2
           SPS = DSIN(TILT)
           CPS = DSQRT (1.0D0 - SPS ** 2)
C
       X2=X*X
       Y2=Y*Y
       Z2=Z*Z
       TPS=SPS/CPS
       HTP=TPS*0.5D0
       GSP=G*SPS
       XSM=X*CPS-Z*SPS
       ZSM=X*SPS+Z*CPS
C
C   CALCULATE THE FUNCTION ZS DEFINING THE SHAPE OF THE TAIL CURRENT SHEET
C    AND ITS SPATIAL DERIVATIVES:
C
       XRC=XSM+RC
       XRC16=XRC**2+16.D0
       SXRC=DSQRT(XRC16)
       Y4=Y2*Y2
       Y410=Y4+1.D4
       SY4=SPS/Y410
       GSY4=G*SY4
       ZS1=HTP*(XRC-SXRC)
       DZSX=-ZS1/SXRC
       ZS=ZS1-GSY4*Y4
       D2ZSGY=-SY4/Y410*4.D4*Y2*Y
       DZSY=G*D2ZSGY
C
C   CALCULATE THE COMPONENTS OF THE RING CURRENT CONTRIBUTION:
C
       XSM2=XSM**2
       DSQT=DSQRT(XSM2+A02)
       FA0=0.5D0*(1.D0+XSM/DSQT)
       DDR=D0+DD*FA0
       DFA0=HA02/DSQT**3
       ZR=ZSM-ZS
       TR=DSQRT(ZR**2+DDR**2)
       RTR=1.D0/TR
       RO2=XSM2+Y2
       ADRT=ADR+TR
       ADRT2=ADRT**2
       FK=1.D0/(ADRT2+RO2)
       DSFC=DSQRT(FK)
       FC=FK**2*DSFC
       FACXY=3.0D0*ADRT*FC*RTR
       XZR=XSM*ZR
       YZR=Y*ZR
       DBXDP=FACXY*XZR
       DER(2,5)=FACXY*YZR
       XZYZ=XSM*DZSX+Y*DZSY
       FAQ=ZR*XZYZ-DDR*DD*DFA0*XSM
       DBZDP=FC*(2.D0*ADRT2-RO2)+FACXY*FAQ
       DER(1,5)=DBXDP*CPS+DBZDP*SPS
       DER(3,5)=DBZDP*CPS-DBXDP*SPS
C
C  CALCULATE THE TAIL CURRENT SHEET CONTRIBUTION:
C
       DELY2=DEL*Y2
       D=DT+DELY2
       IF (DABS(GAM).LT.1.D-6) GOTO 8
       XXD=XSM-XD
       RQD=1.D0/(XXD**2+XLD2)
       RQDS=DSQRT(RQD)
       H=0.5D0*(1.D0+XXD*RQDS)
       HS=-HXLD2M*RQD*RQDS
       GAMH=GAM*H
       D=D+GAMH
       XGHS=XSM*GAM*HS
       ADSL=-D*XGHS
   8   D2=D**2
       T=DSQRT(ZR**2+D2)
       XSMX=XSM-SX
       RDSQ2=1.D0/(XSMX**2+XLW2)
       RDSQ=DSQRT(RDSQ2)
       V=0.5D0*(1.D0-XSMX*RDSQ)
       DVX=HXLW2M*RDSQ*RDSQ2
       OM=DSQRT(DSQRT(XSM2+16.D0)-XSM)
       OMS=-OM/(OM*OM+XSM)*0.5D0
       RDY=1.D0/(P+Q*OM)
       OMSV=OMS*V
       RDY2=RDY**2
       FY=1.D0/(1.D0+Y2*RDY2)
       W=V*FY
       YFY1=2.D0*FY*Y2*RDY2
       FYPR=YFY1*RDY
       FYDY=FYPR*FY
       DWX=DVX*FY+FYDY*Q*OMSV
       YDWY=-V*YFY1*FY
       DDY=DBLDEL*Y
       ATT=AT+T
       S1=DSQRT(ATT**2+RO2)
       F5=1.D0/S1
       F7=1.D0/(S1+ATT)
       F1=F5*F7
       F3=F5**3
       F9=ATT*F3
       FS=ZR*XZYZ-D*Y*DDY+ADSL
       XDWX=XSM*DWX+YDWY
       RTT=1.D0/T
       WT=W*RTT
       BRRZ1=WT*F1
       BRRZ2=WT*F3
       DBXC1=BRRZ1*XZR
       DBXC2=BRRZ2*XZR
       DER(2,1)=BRRZ1*YZR
       DER(2,2)=BRRZ2*YZR
          DER(2,16)=DER(2,1)*TLT2
          DER(2,17)=DER(2,2)*TLT2
       WTFS=WT*FS
       DBZC1=W*F5+XDWX*F7+WTFS*F1
       DBZC2=W*F9+XDWX*F1+WTFS*F3
       DER(1,1)=DBXC1*CPS+DBZC1*SPS
       DER(1,2)=DBXC2*CPS+DBZC2*SPS
       DER(3,1)=DBZC1*CPS-DBXC1*SPS
       DER(3,2)=DBZC2*CPS-DBXC2*SPS
          DER(1,16)=DER(1,1)*TLT2
          DER(1,17)=DER(1,2)*TLT2
          DER(3,16)=DER(3,1)*TLT2
          DER(3,17)=DER(3,2)*TLT2
C
C  CALCULATE CONTRIBUTION FROM THE CLOSURE CURRENTS
C
       ZPL=Z+RT
       ZMN=Z-RT
       ROGSM2=X2+Y2
       SPL=DSQRT(ZPL**2+ROGSM2)
       SMN=DSQRT(ZMN**2+ROGSM2)
       XSXC=X-SXC
       RQC2=1.D0/(XSXC**2+XLWC2)
       RQC=DSQRT(RQC2)
       FYC=1.D0/(1.D0+Y2*RDYC2)
       WC=0.5D0*(1.D0-XSXC*RQC)*FYC
       DWCX=HLWC2M*RQC2*RQC*FYC
       DWCY=DRDYC2*WC*FYC*Y
       SZRP=1.D0/(SPL+ZPL)
       SZRM=1.D0/(SMN-ZMN)
       XYWC=X*DWCX+Y*DWCY
       WCSP=WC/SPL
       WCSM=WC/SMN
       FXYP=WCSP*SZRP
       FXYM=WCSM*SZRM
       FXPL=X*FXYP
       FXMN=-X*FXYM
       FYPL=Y*FXYP
       FYMN=-Y*FXYM
       FZPL=WCSP+XYWC*SZRP
       FZMN=WCSM+XYWC*SZRM
       DER(1,3)=FXPL+FXMN
       DER(1,4)=(FXPL-FXMN)*SPS
       DER(2,3)=FYPL+FYMN
       DER(2,4)=(FYPL-FYMN)*SPS
       DER(3,3)=FZPL+FZMN
       DER(3,4)=(FZPL-FZMN)*SPS
C
C   NOW CALCULATE CONTRIBUTION FROM CHAPMAN-FERRARO SOURCES + ALL OTHER
C
           EX=DEXP(X/DX)
           EC=EX*CPS
           ES=EX*SPS
           ECZ=EC*Z
           ESZ=ES*Z
           ESZY2=ESZ*Y2
           ESZZ2=ESZ*Z2
           ECZ2=ECZ*Z
           ESY=ES*Y
C
           DER(1,6)=ECZ
           DER(1,7)=ES
           DER(1,8)=ESY*Y
           DER(1,9)=ESZ*Z
           DER(2,10)=ECZ*Y
           DER(2,11)=ESY
           DER(2,12)=ESY*Y2
           DER(2,13)=ESY*Z2
           DER(3,14)=EC
           DER(3,15)=EC*Y2
           DER(3,6)=ECZ2*W1
           DER(3,10)=ECZ2*W5
           DER(3,7)=ESZ*W2
           DER(3,11)=-ESZ
           DER(3,8)=ESZY2*W2
           DER(3,12)=ESZY2*W6
           DER(3,9)=ESZZ2*W3
           DER(3,13)=ESZZ2*W4
C
C  FINALLY, CALCULATE NET EXTERNAL MAGNETIC FIELD COMPONENTS,
C    BUT FIRST OF ALL THOSE FOR C.-F. FIELD:
C
      SX1=AK6*DER(1,6)+AK7*DER(1,7)+AK8*DER(1,8)+AK9*DER(1,9)
      SY1=AK10*DER(2,10)+AK11*DER(2,11)+AK12*DER(2,12)+AK13*DER(2,13)
      SZ1=AK14*DER(3,14)+AK15*DER(3,15)+AK610*ECZ2+AK711*ESZ+AK812
     * *ESZY2+AK913*ESZZ2
       BXCL=AK3*DER(1,3)+AK4*DER(1,4)
       BYCL=AK3*DER(2,3)+AK4*DER(2,4)
       BZCL=AK3*DER(3,3)+AK4*DER(3,4)
       BXT=AK1*DER(1,1)+AK2*DER(1,2)+BXCL +AK16*DER(1,16)+AK17*DER(1,17)
       BYT=AK1*DER(2,1)+AK2*DER(2,2)+BYCL +AK16*DER(2,16)+AK17*DER(2,17)
       BZT=AK1*DER(3,1)+AK2*DER(3,2)+BZCL +AK16*DER(3,16)+AK17*DER(3,17)
       F(1)=BXT+AK5*DER(1,5)+SX1+SXA
       F(2)=BYT+AK5*DER(2,5)+SY1+SYA
       F(3)=BZT+AK5*DER(3,5)+SZ1+SZA
C
       RETURN
       END
c
C-----------------------------------------------------------------------
      SUBROUTINE BXYZMU(x, y , z, TILT, abx, aby, abz)
C
C  VERSION 11/01/76
C
C     PURPOSE
C       TO CALCULATE THE CONTRIBUTION TO THE EARTHS MAGNETIC FIELD BY SOURCES 
C       EXTERNAL TO THE EARTH. NO INTERNAL FIELD IS INCLUDED IN THIS ROUTINE.
C
C     METHOD
C        THE ROUTINE INCLUDES THE FIELD CONTRIBUTIONS FROM THE MAGNETOPAUSE 
C        CURRENTS, AND CURRENTS DISTRIBUTED THROUGHOUT THE MAGNETOSPHERE 
C        (THE TAIL AND RING CURRENTS). IT IS VALID FOR ALL TILTS OF THEEARTH'S  
C        DIPOLE AXIS AND IS VALID DURING QUIET MAGNETIC CONDITIONS.
C        A GENERALIZED ORTHONORMAL LEAST-SQUARES PROGRAM WAS USED TO FIT THE 
C        COEFFICIENTS OF A POWER SERIES (INCLUDING EXPONENTIAL TERMS) THROUGH 
C        FOURTH ORDER IN SPACE AND THIRD ORDER IN TILT. THIS EXPANSION HAS
C        BEEN OPTIMIZED FOR THE NEAR EARTH REGION AND IS VALID TO 15 EARTH
C        RADII. OUTSIDE OF THIS REGION THE FIELD DIVERGES RAPIDLY AND A
C        TEMPLATE SETS THE FIELD TO ZERO. IN ORDER TO IMPROVE COMPUTATIONAL 
C        SPEED THE FIELD IS SET TO ZERO BELOW 2 EARTH RADII. (IN THIS REGION 
C        THE EARTH'S INTERNAL FIELD DOMINATES AND THE VARIATIONS EXPRESSED 
C        BY THIS EXPANSION ARE NOT SUFFICIENTLY ACCURATE TO PREDICT VARIATIONS 
C        ON THE EARTH'S SURFACE).
C
C        THE POWER SERIES REPRESENTING THE MAGNETIC FIELD IS
C        BX=SUM OVER I,J,K OF ( A(I,J,K)*X**(I-1)*Y**(2*J-2)*Z**(K-1)
C               + B(I,J,K)*X**(I-1)*Y**(2*J-2)*Z**(K-1)*EXP(-.06*R**2))
C            I GOES FROM 1 TO 5,J FROM 1 TO 3, K FROM 1 TO 5
C            THE SUM OF I + 2*J + K IS LESS THAN OR EQUAL TO 9
C        BY=SUM OVER I,J,K OF ( C(I,J,K)*X**(I-1)*Y**(2*J-1)*Z**(K-1)
C               + D(I,J,K)*X**(I-1)*Y**(2*J-1)*Z**(K-1)*EXP(-.06*R**2))
C            I GOES FROM 1 TO 5, J FROM 1 TO 3, K FROM 1 TO 5
C            THE SUM OF I + 2*J+1 + K IS LESS THAN OR EQUAL TO 9
C        BZ=SUM OVER I,J,K OF ( E(I,J,K)*X**(I-1)*Y**(2*J-2)*Z**(K-1)
C               + F(I,J,K)*X**(I-1)*Y**(2*J-2)*Z**(K-1)*EXP(-.06*R**2))
C            I GOES FROM 1 TO 5, J FROM 1 TO 3, K FROM 1 TO 5
C            THE SUM OF I + 2*J + K IS LESS THAN OR EQUAL TO 9
C        THE COEFFICIENTS A-F ARE DEPENDENT ONLY ON POSITION AND ARE 
C        RECALCULATED EACH TIME THE TILT OF THE DIPOLE IS CHANGED.
C        THE COEEFICIENTS A-F ARE DETERMINED FROM THE TILT DEPENDENT
C        CONSTANTS AA-FF BY THE FOLLOWING EXPRESSIONS
C        A(I,J,K)=AA(I,J,K,1)*TILT**(K-1-(K-1)/2*2)
C                 +AA(I,J,K,2)*TILT**(K+1-(K-1)/2*2)
C        B(I,J,K)=BB......
C        C(I,J,K)=CC......
C        D(I,J,K)=DD......
C        E(I,J,K)=EE(I,J,K,1)*TILT**(K-(K)/2*2)
C                 +EE(I,J,K,2)*TILT**(K+2-(K)/2*2)
C        F(I,J,K)=FF......
C
C     INPUT -- CALLING SEQUENCE
C        X,Y,Z  REALS GIVING THE POSITION WHERE THE MAGNETIC FIELD 
C               IS TO BE EVALUATED. XX(1), XX(2), XX(3) ARE RESPECTIVELY 
C               THE X, Y, AND Z SOLAR MAGNETIC COORDINATES IN EARTH RADII. 
C               Z IS ALONG THE EARTH'S NORTH DIPOLE AXIS, X IS PERPENDICULAR 
C               TO Z AND IN THE PLANE CONTAINING THE Z AXIS AND THE SUN-EARTH
C               LINE, Y IS PERPENDICULAR TO X AND Z FORMING A RIGHT HANDED
C               COORDINATE SYSTEM. X IS POSITIVE IN THE SOLAR DIRECTION.
C        TILT   IS THE TILT OF THE DIPOLE AXIS IN DEGREES. IT IS THE COMPLEMENT 
C               OF THE ANGLE BETWEEN THE NORTH DIPOLE AXIS AND THE SOLAR 
C               DIRECTION (PSI). TILT=90-PSI.
C
C     OUTPUT -- CALLING SEQUENCE
C        ABX,ABY,ABZ  THE X, Y, AND Z COMPONENTS OF THE MAGNETOSPHERIC 
C               MAGNETIC FIELD IN GAMMA. 
C
C     CONSTANTS
C        AA,BB,CC,DD,EE,FF ARE REAL ARRAYS CONTAINING THE TILT DEPENDENT
C               COEFFICIENTS. AA(I,J,K,L) ARE STORED SUCH THAT L VARIES MOST
C               RAPIDLY, FOLLOWED IN ORDER BY K, J AND I. I VARIES THE SLOWEST.
C               THE ARRAY IS CLOSE PACKED AND ALL COEFFICIENTS THAT ARE ZERO 
C               BECAUSE OF SYMMETRY OR BECAUSE THE CROSS TERM POWER IS TOO 
C               LARGE ARE DELETED.
C
C     VARIABLES
C        A,B,C,D,E,F  THE TILT INDEPENDENT COEFFICIENTS. THEIR USE IS 
C               DESCRIBED UNDER METHOD.
C        ITA    A REAL ARRAY WHICH CONTAINS THE SYMMETRY OF THE TILT
C               DEPENDENCE FOR EACH OF THE A AND B COEFFICIENTS
C               ITA(1) HAS THE SYMMETRY INFORMATION FOR A(1,1,1,1)
C                      AND A(1,1,1,2)
C               ITA(2) HAS THE SYMMETRY INFORMATION FOR A(1,1,2,1)
C                      AND A(1,1,2,2) ETC.
C               IF ITA = 1 TILT SYMMETRY IS EVEN WITH RESPECT TO Z SYM.
C               IF ITA = 2 TILT SYMMETRY IS ODD WITH RESPECT TO Z SYM.
C        ITB    SYMMETRY POINTER FOR C AND D ARRAYS
C        ITC    SYMMETRY POINTER FOR E AND F ARRAYS
C        X      X COMPONENT OF POSITION
C        Y      Y COMPONENT OF POSITION
C        Z      Z COMPONENT OF POSITION
C        Y2     Y**2
C        Z2     Z**2
C        R2     X**2 + Y**2 + Z**2
C        R      SQRT(R2)
C        I      DO LOOP VARIABLE. IN THE FIELD EXPANSION LOOP IT REPRESENTS 
C               THE POWER TO WHICH X IS CURRENTLY RAISED, I.E. X**(I-1)
C        J      DO LOOP VARIABLE. ALSO Y**(2*J-2)
C        K      DO LOOP VARIABLE. ALSO Z**(K-1)
C        XB     X**(I-1)
C        YEXB   X**(I-1)*Y**(2*J-2)
C        ZEYEXB X**(I-1)*Y**(2*J-2)*Z**(K-1)
C        IJK    I + 2*J + K
C        II     POINTS TO THE ARRAY LOCATION WHERE THE CURRENT POWER SERIES 
C               COEFFICIENT FOR BX IS LOCATED.
C        JJ     BY COEFFICIENT LOCATION POINTER
C        KK     BZ COEFFICIENT LOCATION POINTER
C        BX,BY,BZ ARE USED TO CONSTRUCT THE MAGNETIC FIELD WITHIN THE POWER 
C               SERIES LOOP.
C        EXPR   EXP(-.06*R2)
C        TILTL  HOLDS THE LAST VALUE OF THE TILT FOR WHICH THE TILT
C               INDEPENDENT COEFFICIENTS A-F WERE CALCULATED.
C        TT     A REAL ARRAY HOLDING THE POWERS OF THE TILT.
C               TT(1)=TILT**0, TT(2)=TILT**1, ETC.
C        CON    =0 FOR R LESS THAN 2
C               =1 FOR R GREATER THAN 2.5
C               GOES FROM 0 TO 1 IN THE REGION 2 TO 2.5.
C
C     FOR MORE INFORMATION CALL OR WRITE K.A. PFITZER OR W.P. OLSON AT 
C     MCDONNEL DOUGLAS ASTRONAUTICS CO. 5301 BOLSA AVE, HUNTINGTON
C     CALIF., PHONE (714) 896-3231.
C

      IMPLICIT NONE

C Input
      REAL*8 X, Y, Z, TILT
C Output
      REAL*8 ABX, ABY, ABZ
C Local
      REAL*8 AA(64), BB(64), CC(44), DD(44), EE(64), FF(64), A(32)
      REAL*8 B(32), C(22), D(22), E(32), F(32), TT(4), TILTL
      REAL*8 Y2, Z2, R2, BX, BY, BZ, CON, EXPR, XB, YEXB, ZEYEXB
      INTEGER*4 ITA(32), ITB(22), ITC(32), I, J, K, II, JJ, KK, IJK
      SAVE TILTL,TT,A,B,E,F,C,D

      DATA ITA /2,1,2,1,2,2,1,2,1,2,1,2,1,2,1,2,2,1,2,2,2,1,
     *          2,1,2,1,2,1,2,2,2,1/
      DATA ITB /2,1,2,1,2,2,1,2,2,2,1,2,1,2,1,2,1,2,2,2,1,2/
      DATA ITC /1,2,1,2,1,1,2,1,2,1,2,1,2,1,2,1,1,2,1,1,1,2,
     *          1,2,1,2,1,2,1,1,1,2/

      DATA AA/-2.26836d-02,-1.01863d-04,3.42986d+00,
     *-3.12195d-04, 9.50629d-03,-2.91512d-06,-1.57317d-03, 8.62856d-08,
     *-4.26478d-05, 1.62924d-08,-1.27549d-04, 1.90732d-06,-1.65983d-02,
     * 8.46680d-09,-5.55850d-05, 1.37404d-08, 9.91815d-05, 1.59296d-08,
     * 4.52864d-07,-7.17669d-09, 4.98627d-05, 3.33662d-10,-5.97620d-02,
     * 1.60669d-05,-2.29457d-01,-1.43777d-04, 1.09403d-03,-9.15606d-07,
     * 1.60658d-03,-4.01198d-07,-3.15064d-06, 2.03125d-09, 4.92887d-04,
     *-1.80676d-07,-1.12022d-03, 5.98568d-07,-5.90009d-06, 5.16504d-09,
     *-1.48737d-06, 4.83477d-10,-7.44379d-04, 3.82472d-06, 7.41737d-04,
     *-1.31468d-05,-1.24729d-04, 1.92930d-08,-1.91764d-04,-5.30371d-08,
     * 1.38186d-05,-2.81594d-08, 7.46386d-06, 2.64404d-08, 2.45049d-04,
     *-1.81802d-07,-1.00278d-03, 1.98742d-06,-1.16425d-05, 1.17556d-08,
     *-2.46079d-06,-3.45831d-10, 1.02440d-05,-1.90716d-08,-4.00855d-05,
     * 1.25818d-07/

      DATA BB/ 9.47753d-02, 1.45981d-04,-1.82933d+00,
     * 5.54882d-04, 5.03665d-03,-2.07698d-06, 1.10959d-01,-3.45837d-05,
     *-4.40075d-05, 5.06464d-07,-1.20112d-03, 3.64911d-06, 1.49849d-01,
     *-7.44929d-05, 2.46382d-04, 9.65870d-07,-9.54881d-04, 2.43647d-07,
     * 3.06520d-04, 3.07836d-07, 6.48301d-03, 1.26251d-06,-7.09548d-03,
     *-1.55596d-05, 3.06465d+00,-7.84893d-05, 4.95145d-03, 3.71921d-06,
     *-1.52002d-01, 6.81988d-06,-8.55686d-05,-9.01230d-08,-3.71458d-04,
     * 1.30476d-07,-1.82971d-01, 1.51390d-05,-1.45912d-04,-2.22778d-07,
     * 6.49278d-05,-3.72758d-08,-1.59932d-03, 8.04921d-06, 5.38012d-01,
     *-1.43182d-04, 1.50000d-04, 5.88020d-07,-1.59000d-02, 1.60744d-06,
     * 3.17837d-04, 1.78959d-07,-8.93794d-03, 6.37549d-06, 1.27887d-03,
     *-2.45878d-07,-1.93210d-01, 6.91233d-06,-2.80637d-04,-2.57073d-07,
     * 5.78343d-05, 4.52128d-10, 1.89621d-04,-4.84911d-08,-1.50058d-02,
     * 6.21772d-06/

      DATA CC/-1.88177d-02,-1.92493d-06,-2.89064d-01,
     *-8.49439d-05,-4.76380d-04,-4.52998d-08, 1.61086d-03, 3.18728d-07,
     * 1.29159d-06, 5.52259d-10, 3.95543d-05, 5.61209d-08, 1.38287d-03,
     * 5.74237d-07, 1.86489d-06, 7.10175d-10, 1.45243d-07,-2.97591d-10,
     *-2.43029d-03,-6.70000d-07,-2.30624d-02,-6.22193d-06,-2.40815d-05,
     * 2.01689d-08, 1.76721d-04, 3.78689d-08, 9.88496d-06, 7.33820d-09,
     * 7.32126d-05, 8.43986d-08, 8.82449d-06,-6.11708d-08, 1.78881d-04,
     * 8.62589d-07, 3.43724d-06, 2.53783d-09,-2.04239d-07, 8.16641d-10,
     * 1.68075d-05, 7.62815d-09, 2.26026d-04, 3.66341d-08, 3.44637d-07,
     * 2.25531d-10/

      DATA DD/ 2.50143d-03, 1.01200d-06, 3.23821d+00,                   
     * 1.08589d-05,-3.39199d-05,-5.27052d-07,-9.46161d-02,-1.95413d-09,
     *-4.23614d-06, 1.43153d-08,-2.62948d-04, 1.05138d-07,-2.15784d-01,
     *-2.20717d-07,-2.65687d-05, 1.26370d-08, 5.88917d-07,-1.13658d-08,
     * 1.64385d-03, 1.44263d-06,-1.66045d-01,-1.46096d-05, 1.22811d-04,
     * 3.43922d-08, 9.66760d-05,-6.32150d-07,-4.97400d-05,-2.78578d-08,
     * 1.77366d-02, 2.05401d-07,-1.91756d-03,-9.49392d-07,-1.99488d-01,
     *-2.07170d-06,-5.40443d-05, 1.59289d-08, 7.30914d-05, 3.38786d-08,
     *-1.59537d-04,-1.65504d-07, 1.90940d-02, 2.03238d-06, 1.01148d-04,
     * 5.20815d-08/

      DATA EE/-2.77924d+01,-1.01457d-03, 9.21436d-02,                   
     *-8.52177d-06, 5.19106d-01, 8.28881d-05,-5.59651d-04, 1.16736d-07,
     *-2.11206d-03,-5.35469d-07, 4.41990d-01,-1.33679d-05,-7.18642d-04,
     * 6.17358d-08,-3.51990d-03,-5.29070d-07, 1.88443d-06,-6.60696d-10,
     *-1.34708d-03, 1.02160d-07, 1.58219d-06, 2.05040d-10, 1.18039d+00,
     * 1.58903d-04, 1.86944d-02,-4.46477d-06, 5.49869d-02, 4.94690d-06,
     *-1.18335d-04, 6.95684d-09,-2.73839d-04,-9.17883d-08, 2.79126d-02,
     *-1.02567d-05,-1.25427d-04, 3.07143d-08,-5.31826d-04,-2.98476d-08,
     *-4.89899d-05, 4.91480d-08, 3.85563d-01, 4.16966d-05, 6.74744d-04,
     *-2.08736d-07,-3.42654d-03,-3.13957d-06,-6.31361d-06,-2.92981d-09,
     *-2.63883d-03,-1.32235d-07,-6.19406d-06, 3.54334d-09, 6.65986d-03,
     *-5.81949d-06,-1.88809d-04, 3.62055d-08,-4.64380d-04,-2.21159d-07,
     *-1.77496d-04, 4.95560d-08,-3.18867d-04,-3.17697d-07,-1.05815d-05,
     * 2.22220d-09/

      DATA FF/-5.07092d+00, 4.71960d-03,-3.79851d-03,                   
     *-3.67309d-06,-6.02439d-01, 1.08490d-04, 5.09287d-04, 5.62210d-07,
     * 7.05718d-02, 5.13160d-06,-2.85571d+00,-4.31728d-05, 1.03185d-03,
     * 1.05332d-07, 1.04106d-02, 1.60749d-05, 4.18031d-05, 3.32759d-08,
     * 1.20113d-01, 1.40486d-05,-3.37993d-05, 5.48340d-09, 9.10815d-02,
     *-4.00608d-04, 3.75393d-03,-4.69939d-07,-2.48561d-02, 1.31836d-04,
     *-2.67755d-04,-7.60285d-08, 3.04443d-03,-3.28956d-06, 5.82367d-01,
     * 5.39496d-06,-6.15261d-04, 4.05316d-08, 1.13546d-02,-4.26493d-06,
     *-2.72007d-02, 5.72523d-08,-2.98576d+00, 3.07325d-05, 1.51645d-03,
     * 1.25098d-06, 4.07213d-02, 1.05964d-05, 1.04232d-04, 1.77381d-08,
     * 1.92781d-01, 2.15734d-05,-1.65741d-05,-1.88683d-09, 2.44803d-01,
     * 1.51316d-05,-3.01157d-04, 8.47006d-08, 1.86971d-02,-6.94074d-06,
     * 9.13198d-03,-2.38052d-07, 1.28552d-01, 6.92595d-06,-8.36792d-05,
     *-6.10021d-08/

      DATA TILTL /99.0D0/
C
C  SET UP SOME OF THE INITIAL POSITION VARIABLES.
      Y2 = Y**2
      Z2 = Z**2
      R2 = X**2 + Y2 + Z2
C
C  SET MAGNETIC FIELD VARIABLES TO ZERO.
      BX = 0.0D0
      BY = 0.0D0
      BZ = 0.0D0
C
C  CHECK TO SEE IF POSITION IS WITHIN REGION OF VALIDITY.
      CON = 1.0D0
C  IF DISTANCE TOO LARGE TAKE ERROR EXIT.
      IF (R2 .GT. 225.0D0) GO TO 50
C  IF DISTANCE TOO SMALL SET FIELD TO ZERO AND EXIT.
      IF (R2 .LT. 4.0D0) GO TO 40
      IF (R2 .LT. 6.25D0) CON = CON * (R2-4.0D0) / 2.25D0
C
C  IF TILT HAS NOT CHANGED, GO DIRECTLY TO FIELD CALCULATION.
      IF (TILTL .EQ. TILT) GO TO 6
C  SET UP POWERS OF TILT.
      TILTL = TILT
      TT(1) = 1.0D0
      TT(2) = TILTL
      TT(3) = TILTL**2
      TT(4) = TILT * TT(3)
C
C  SET UP THE X AND Z COMPONENT TILT INDEPENDENT COEFFICIENTS.
      DO 1 I=1,32
        J = (I-1) * 2 + 1
        K = ITA(I)
        A(I) = AA(J) * TT(K) + AA(J+1) * TT(K+2)
        B(I) = BB(J) * TT(K) + BB(J+1) * TT(K+2)
        K = ITC(I)
        E(I) = EE(J) * TT(K) + EE(J+1) * TT(K+2)
        F(I) = FF(J) * TT(K) + FF(J+1) * TT(K+2)
   1  CONTINUE
C
C  SET UP THE Y COMPONENT TILT INDEPENDENT COEFFICIENTS.
      DO 2 I=1,22
        J = (I-1) * 2 + 1
        K = ITB(I)
        C(I) = CC(J) * TT(K) + CC(J+1) * TT(K+2)
        D(I) = DD(J) * TT(K) + DD(J+1) * TT(K+2)
   2  CONTINUE
   6  EXPR = DEXP(-0.06D0*R2)
C
C  INITIALIZE THE POINTERS.
      II=1
      JJ=1
      KK=1
      XB=1.0D0
C
C  BEGIN SUM OVER X.
      DO 30 I=1,5
        YEXB = XB
C
C  BEGIN SUM OVER Y.
        DO 20 J=1,3
          IF (I+2*J .GT. 8) GO TO 25
            ZEYEXB = YEXB
            IJK = I + 2 * J + 1
            K = 1
C
C  Z LOOP STARTS HERE.
  10  BZ = BZ + (E(KK)+F(KK)*EXPR) * ZEYEXB
      KK = KK + 1
      BX = BX + (A(II)+B(II)*EXPR) * ZEYEXB
      II = II + 1
      IF (IJK .GT. 8) GO TO 15
      BY = BY + (C(JJ)+D(JJ)*EXPR) * ZEYEXB * Y
      JJ = JJ + 1
      ZEYEXB = ZEYEXB * Z
      IJK = IJK + 1
      K = K + 1
      IF ((IJK .LE. 9) .AND. (K .LE. 5)) GO TO 10
  15  YEXB = YEXB * Y2
  20  CONTINUE
  25  XB = XB * X
  30  CONTINUE
C
C  SET UP THE OUTPUT ARRAY, MULTIPLY BY CON. CON IS NORMALLY ONE
C  BUT INSIDE OF R=2.5 IT GOES TO ZERO. INSIDE R=2 IT IS ZERO.
  40  ABX = BX * CON
      ABY = BY * CON
      ABZ = BZ * CON
      RETURN

C  ERROR EXIT IF OUTSIDE OF R = 15.
  50  CONTINUE
c     WRITE(6,*)' OUTSIDE THE VALID REGION--POWER SERIES DIVERGES 
c    &BFIELD IS SET TO ZERO'
  60  FORMAT(4H X= ,E10.3,4H Y= ,E10.3,4H Z= ,E10.3,76H 'IS OUTSIDE THE'
     &'VALID REGION--POWER SERIES DIVERGES BFIELD IS SET TO ZERO')
      GO TO 40

      END
C-----------------------------------------------------------------------
      SUBROUTINE BDYN(den, vel, dst, X, Y, Z, ABX, ABY, ABZ) 

      IMPLICIT NONE

C Input
      REAL*8 den, vel, dst, x, y, z
C Output
      REAL*8 abx, aby, abz
C Local
      REAL*8 XX(3), B(3), SOFFD, SRING, STAIL, STDOFF, RINGST
C
      XX(1) = X
      XX(2) = Y
      XX(3) = Z
C
      SOFFD = STDOFF(VEL, DEN)
      SRING = RINGST(SOFFD, DST)
      STAIL = 1.0D0
    
C     GET THE EXTERNAL CONTRIBUTION 
      CALL BDYNAM(XX, B, SOFFD, SRING, STAIL)
      ABX = B(1)
      ABY = B(2)
      ABZ = B(3)

      RETURN
      END 
C-----------------------------------------------------------------------
      SUBROUTINE BDYNAM(XX, BB, SOFFD, SRING, STAIL)
C  VERSION 5/13/88
C  DEVELOPED MCDONNELL DOUGLAS
C  FOR INFORMATION CALL KARL PFITZER (714) 896-3231
C     
C     PURPOSE 
C        CALCULATE THE TOTAL EXTERNAL MAGNETIC FIELD DURING DISTURBED TIMES.
C
C     METHOD
C        CALLS THE EXTERNAL QUIET TIME SUBROUTINES AND COMBINES THEM
C        ACCORDING TO THE DYNAMIC SCALING ALGORITHMS.
C     
C     INPUT -- ARGUMENT LIST
C        XX     A REAL ARRAY GIVING THE POSITION WHERE THE MAGNETIC FIELD 
C               IS TO BE DETERMINED. XX(1), XX(2), XX(3) ARE RESPECTIVELY  
C               THE X, Y, AND Z SOLAR MAGNETIC COORDINATES IN EARTH RADII.
C               Z IS ALONG THE EARTH'S NORTH DIPOLE AXIS. X IS PERPENDICULAR 
C               TO Z AND IN THE PLANE CONTAINING THE Z AXIS AND THE SUN-EARTH  
C               LINE (X IS POSITIVE IN THE SOLAR DIRECTION). Y IS PERPENDICULAR 
C               TO X AND Z AND X Y Z FORM A RIGHT HANDED COORDINATE SYSTEM. 
C     
C        SOFFD  THE STANDOFF DISTANCE OF THE MAGNETOPAUSE. THE QUIET STANDOFF 
C               DISTANCE IS 10.5 EARTH RADII. ACCEPTABLE VALUES RANGE BETWEEN 
C               6 AND 11. THIS VALUE IS USED TO CALCULATE THE STRENGTH OF
C               THE MAGNETOPAUSE CURRENTS AND TO SCALE THE SIZE OF THE 
C               MAGNETOPAUSE. THIS VALUE ALSO SCALES THE SIZE OF THE TAIL
C               CURRENT SYSTEM. THE RING SYSTEM IS NOT SCALED, SINCE ITS 
C               SOURCE IS PRIMARILY AT RADIAL DISTANCES.
C
C        SRING  RELATIVE STRENGTH OF THE RING CURRENT. A VALUE OF ONE
C               UTILIZES THE NOMINAL QUIET RING VALUES BUILT INTO THE BASIC 
C               MODEL. THIS BASIC MODEL HAS A MAXIMUM RING DEPRESSION OF
C               40 NT AT L=4 RE. IF SRING IS SET TO 2 THE RING DEPRESSION
C               WOULD BE 80 NT.
C
C        STAIL  A TAIL CURRENT STRENGTH MULTIPLIER. WHEN STAIL IS EQUALTO 
C               1.0 THEN THE TAIL SCALES WITH THE STRENGTH OF THE MAGNETOPAUSE 
C               CURRENTS. TO WEAKEN THE TAIL FROM THIS VALUE USE VALUES 
C               LESS THAN 1.0, TO STRENGTHEN USE VALUES GREATER THAN 1.0.
C
C     OUTPUT -- ARGUMENT LIST 
C        BB     A REAL ARRAY CONTAINING THE X, Y, AND Z COMPONENTS OF THE 
C               EARTH'S TOTAL MAGNETIC FIELD IN SOLAR MAGNETIC COORDINATES.  
C               BB(1), BB(2) AND BB(3) ARE THE BX, BY, AND BZ COMPONENTS.  
C               THE UNITS ARE GAUSS. 
C

      IMPLICIT NONE

C Input
      REAL*8 XX(3), SOFFD, SRING, STAIL
C Output
      REAL*8 BB(3) 
C Local
      REAL*8 XXX(3), BM(3), BR(3), BT(3), xxxx(3)
      REAL*8 SCL, STRMAG, STRRIN, STRTAI 
      INTEGER*4 I

C   CALCULATE MAGNETOPAUSE SCALE FACTOR.
      SCL = 10.5D0 / SOFFD
C   CALCULATE STRENGTH OF MAGNETOPAUSE CURRENTS.
      STRMAG = SCL**3.0D0 
C   SET STRENGTH OF RING AND TAIL.
      STRRIN = Sring
      STRTAI = STAIL * STRMAG 
C   CALCULATE SCALED DISTANCES
      DO I=1,3 
        XXX(I) = XX(I) * SCL
        XXXX(I) = XX(I) * 1.0D0
      END DO
C   CALL THE QUIET TIME SUBROUTINE.
      CALL BFMAGP(XXX, BM) 
      CALL BFRING(XXXX, BR)
      CALL BFTAIL(XXX, BT) 
C   COMBINE THE COMPONENTS OF THE MAGNETIC FIELD ACCORDING TO THEIR
C   RELATIVE STRENGTHS.
      DO I=1,3 
        BB(I) = STRMAG * BM(I) + STRRIN * BR(I) + STRTAI * BT(I)
      END DO

      RETURN
      END 
C-----------------------------------------------------------------------
      SUBROUTINE BFMAGP(XX, BB)
C  VERSION 5/13/88
C  DEVELOPED MCDONNELL DOUGLAS
C  FOR INFORMATION CALL KARL PFITZER (714) 896-3231
C
C     INPUT -- ARGUMENT LIST
C        XX     A REAL ARRAY GIVING THE POSITION WHERE THE MAGNETIC FIELD 
C               IS TO BE DETERMINED. XX(1), XX(2), XX(3) ARE RESPECTIVELY 
C               THE X, Y, AND Z SOLAR MAGNETIC COORDINATES IN EARTH RADII.
C               Z IS ALONG THE EARTH'S NORTH DIPOLE AXIS. X IS PERPENDICULAR 
C               TO Z AND IN THE PLANE CONTAINING THE Z AXIS AND THE SUN-EARTH
C               LINE (X IS POSITIVE IN THE SOLAR DIRECTION). Y IS PERPENDICULAR 
C               TO X AND Z AND X Y Z FORM A RIGHT HANDED COORDINATE SYSTEM. 
C     
C     OUTPUT -- ARGUMENT LIST 
C        BB     A REAL ARRAY CONTAINING THE X, Y, AND Z COMPONENTS OF THE 
C               EARTH'S TOTAL MAGNETIC FIELD IN SOLAR MAGNETIC COORDINATES.  
C               BB(1), BB(2) AND BB(3) ARE THE BX, BY, AND BZ COMPONENTS.
C               THE UNITS ARE NANOTESLA.
C
C THIS SUBROUTINE CONTAINS NEW, REFITTED COEFFICIENTS FOR COMPUTING ALL THE 
C B-MAGNETOPAUSE COMPONENTS.  
C THE FORM OF THE EXPANSION IS GIVEN IN THE NEXT FEW STATEMENTS.
C     
C BMPX=(1 + X + X**2 + X**3 + X**4 + 1/(10-(X-2)**2))*(1 + Y**2 + Y**4)*
C      (Z + Z**3 + Z**5). 
C BMPY=(1 + X + X**2 + X**3 + X**4)*(Y + Y**3 + Y**5)*(Z + Z**3 + Z**5).
C BMPZ=(1 + X + X**2 + X**3 + X**4 + 1/(15-X) + 1/((30-X)**2))* 
C      (1 + Y**2 + Y**4)*(1 + Z**2 + Z**4). 
C   
C COEFFICIENTS COMPUTED FROM COMBINED OLSON DATA-SETS BOUNDARY.DAT AND
C BOUND.DAT, FOR Y>0 AND Z>0 ONLY, I.E. A TOTAL OF 1009 DATA POINTS. 
C 2 EXTRA EXTRAPOLATED POINTS ADDED FOR Z-COEFF, TO IMPROVE FIT, 
C NAMELY X=10,Y=Z=0,BZ=29 AND X=11,Y=Z=0,BZ=30.25.
C 
C***PRELIMINARY ROUTINE 
C***VALID TO APPROXIMATELY -60 RE 
C***MAGNETIC FIELD FROM MAGNETOPAUSE CURRENTS ONLY
C

      IMPLICIT NONE

C Input
      REAL*8 XX(3)
C Output
      REAL*8 BB(3)
C Local
      REAL*8 A(54), B(45), C(63), C1(30), C2(33), X(5), Y(5), Z(5)
      REAL*8 XN, YN, ZN, FX1, XF1, XF2
      INTEGER*4 I

      EQUIVALENCE (C(1),C1(1)), (C(31),C2(1))

      DATA A/ 
     *  0.113275039d+01, 0.354408138d-01,-0.152252289d-02,
     * -0.683306571d-04,-0.642841428d-06,-0.121504674d-01,
     * -0.839622808d-03,-0.167520029d-04,-0.385962942d-07,
     *  0.107674747d-08, 0.558984066d-04, 0.551508083d-05,
     *  0.206288036d-06, 0.335316730d-08, 0.198413126d-10,
     * -0.545824692d-02,-0.264107861d-03, 0.143533146d-05,
     *  0.195177861d-06, 0.207546358d-08, 0.211199178d-03,
     *  0.220245929d-04, 0.860991804d-06, 0.145349395d-07,
     *  0.886173426d-10,-0.949615014d-06,-0.110830563d-06,
     * -0.477998707d-08,-0.873645670d-10,-0.569051859d-12,
     *  0.271760982d-04, 0.266707661d-05, 0.994617153d-07,
     *  0.167023062d-08, 0.104617062d-10,-0.989193381d-06,
     * -0.113236254d-06,-0.482686247d-08,-0.880319914d-10,
     * -0.575385009d-12, 0.487020380d-08, 0.586310778d-09,
     *  0.260182431d-10, 0.488435735d-12, 0.326678627d-14,
     *  0.193470073d+01, 0.402453184d+00,-0.193471275d-02,
     *  0.682263300d-01,-0.576195028d-02, 0.237557251d-04,
     * -0.529665092d-03, 0.255710365d-04,-0.120115033d-06/
C
      DATA B/ 
     * -0.519952811d-01,-0.230140495d-02, 0.146173188d-03,
     *  0.809832090d-05, 0.888401672d-07,-0.370911323d-03,
     * -0.101231737d-03,-0.742647399d-05,-0.196170248d-06,
     * -0.165503899d-08, 0.150949325d-05, 0.308240260d-06,
     *  0.195390104d-07, 0.472441419d-09, 0.375989214d-11,
     * -0.422217818d-04,-0.621468353d-04,-0.620102765d-05,
     * -0.189322407d-06,-0.172039538d-08, 0.445292017d-05,
     *  0.118324999d-05, 0.855768008d-07, 0.223059815d-08,
     *  0.183677951d-10,-0.550030643d-08,-0.150351465d-08,
     * -0.107031245d-09,-0.268793755d-11,-0.205845354d-13,
     *  0.813519478d-06, 0.279971147d-06, 0.227601529d-07,
     *  0.643000209d-09, 0.561745876d-11,-0.983297266d-08,
     * -0.265465072d-08,-0.194798427d-09,-0.513382522d-11,
     * -0.420117906d-13,-0.469398392d-11,-0.543405219d-12,
     * -0.121854998d-13,-0.483310746d-16,-0.429469692d-17/
C
      DATA C1/
     *  0.406363373d+02, 0.291153884d+01, 0.991215929d-01,
     *  0.161603605d-02, 0.994476977d-05,-0.566497850d+01,
     * -0.346289247d+00,-0.102486340d-01,-0.153071058d-03,
     * -0.892381365d-06, 0.182735808d-01, 0.106282183d-02,
     *  0.311990625d-04, 0.464014079d-06, 0.269492229d-08,
     * -0.102119482d+01,-0.649643913d-01,-0.205774955d-02,
     * -0.323610875d-04,-0.195236396d-06, 0.531459488d-01,
     *  0.324825896d-02, 0.991819543d-04, 0.152400162d-05,
     *  0.907312536d-08,-0.132267553d-03,-0.871756401d-05,
     * -0.262251859d-06,-0.395617938d-08,-0.232419934d-10/

      DATA C2/
     *  0.144323579d-02, 0.799393092d-04, 0.322526876d-05,
     *  0.596131713d-07, 0.395406097d-09,-0.839159111d-05,
     *  0.564246250d-06,-0.212045990d-07,-0.866837990d-09,
     * -0.746255575d-11,-0.685688633d-06,-0.523054773d-07,
     * -0.130326583d-08,-0.157964718d-10,-0.759061461d-13,
     *  0.836994324d+02,-0.609500999d+02, 0.100208335d+00,
     * -0.688268995d+01, 0.397136599d+00,-0.250137411d-02,
     * -0.594024621d-01, 0.457714684d-02,-0.449951913d-04,
     * -0.273244004d+05, 0.875882129d+04,-0.227706509d+02,
     *  0.129124341d+04,-0.715722046d+02, 0.266965359d+00,
     *  0.240404391d+01,-0.269608498d+00, 0.332747493d-02/
C
      XN = XX(1)
      YN = XX(2)
      ZN = XX(3)
      DO 1 I=1,5  
        X(I) = XN 
        Y(I) = YN 
        Z(I) = ZN 
        XN = XN * XX(1) 
        YN = YN * XX(2) 
        ZN = ZN * XX(3) 
    1 CONTINUE
      FX1 = 1.0D0 / (10.0D0+(X(1)-2.0D0)**2.0D0) 
      XF1 = 1.0D0 / (15.0D0-X(1)) 
      XF2 = 1.0D0 / ((30.0D0-X(1))**2.0D0)
C
      BB(1)=Z(1)*(A(1)+A(2)*X(1)+A(3)*X(2)+A(4)*X(3)+A(5)*X(4))+
     * Z(3)*(A(6)+A(7)*X(1)+A(8)*X(2)+A(9)*X(3)+A(10)*X(4))+
     * Z(5)*(A(11)+A(12)*X(1)+A(13)*X(2)+A(14)*X(3)+A(15)*X(4))+
     * Y(2)*Z(1)*(A(16)+A(17)*X(1)+A(18)*X(2)+A(19)*X(3)+A(20)*X(4))+ 
     * Y(2)*Z(3)*(A(21)+A(22)*X(1)+A(23)*X(2)+A(24)*X(3)+A(25)*X(4))+ 
     * Y(2)*Z(5)*(A(26)+A(27)*X(1)+A(28)*X(2)+A(29)*X(3)+A(30)*X(4))+ 
     * Y(4)*Z(1)*(A(31)+A(32)*X(1)+A(33)*X(2)+A(34)*X(3)+A(35)*X(4))+ 
     * Y(4)*Z(3)*(A(36)+A(37)*X(1)+A(38)*X(2)+A(39)*X(3)+A(40)*X(4))+ 
     * Y(4)*Z(5)*(A(41)+A(42)*X(1)+A(43)*X(2)+A(44)*X(3)+A(45)*X(4))+ 
     * FX1*(A(46)*Z(1)+A(47)*Z(3)+A(48)*Z(5))+FX1*Y(2)*(A(49)*Z(1)+ 
     * A(50)*Z(3)+A(51)*Z(5))+FX1*Y(4)*(A(52)*Z(1)+A(53)*Z(3)+A(54)*
     * Z(5))
C
      BB(2)=Z(1)*Y(1)*(B(1)+B(2)*X(1)+B(3)*X(2)+B(4)*X(3)+B(5)*X(4))+ 
     * Z(3)*Y(1)*(B(6)+B(7)*X(1)+B(8)*X(2)+B(9)*X(3)+B(10)*X(4))+ 
     * Z(5)*Y(1)*(B(11)+B(12)*X(1)+B(13)*X(2)+B(14)*X(3)+B(15)*X(4))+ 
     * Y(3)*Z(1)*(B(16)+B(17)*X(1)+B(18)*X(2)+B(19)*X(3)+B(20)*X(4))+ 
     * Y(3)*Z(3)*(B(21)+B(22)*X(1)+B(23)*X(2)+B(24)*X(3)+B(25)*X(4))+ 
     * Y(3)*Z(5)*(B(26)+B(27)*X(1)+B(28)*X(2)+B(29)*X(3)+B(30)*X(4))+ 
     * Y(5)*Z(1)*(B(31)+B(32)*X(1)+B(33)*X(2)+B(34)*X(3)+B(35)*X(4))+ 
     * Y(5)*Z(3)*(B(36)+B(37)*X(1)+B(38)*X(2)+B(39)*X(3)+B(40)*X(4))+ 
     * Y(5)*Z(5)*(B(41)+B(42)*X(1)+B(43)*X(2)+B(44)*X(3)+B(45)*X(4))
C
       BB(3)=C(1)+C(2)*X(1)+C(3)*X(2)+C(4)*X(3)+C(5)*X(4)+
     * Z(2)*(C(6)+C(7)*X(1)+C(8)*X(2)+C(9)*X(3)+C(10)*X(4))+
     * Z(4)*(C(11)+C(12)*X(1)+C(13)*X(2)+C(14)*X(3)+C(15)*X(4))+
     * Y(2)*(C(16)+C(17)*X(1)+C(18)*X(2)+C(19)*X(3)+
     * C(20)*X(4))+Y(2)*Z(2)*(C(21)+C(22)*X(1)+C(23)*X(2)+
     * C(24)*X(3)+C(25)*X(4))+Y(2)*Z(4)*(C(26)+C(27)*X(1)+
     * C(28)*X(2)+C(29)*X(3)+C(30)*X(4))+Y(4)*(C(31)+ 
     * C(32)*X(1)+C(33)*X(2)+C(34)*X(3)+C(35)*X(4)) 
      BB(3)=BB(3)+
     * Y(4)*Z(2)*(C(36)+C(37)*X(1)+C(38)*X(2)+C(39)*X(3)+C(40)*X(4))+ 
     * Y(4)*Z(4)*(C(41)+C(42)*X(1)+C(43)*X(2)+C(44)*X(3)+ 
     * C(45)*X(4))+XF1*(C(46)+C(47)*Z(2)+C(48)*Z(4)+
     * C(49)*Y(2)+C(50)*Y(2)*Z(2)+C(51)*Y(2)*Z(4)+C(52)*Y(4)+ 
     * C(53)*Y(4)*Z(2)+C(54)*Y(4)*Z(4))+XF2*(C(55)+C(56)*Z(2)+
     * C(57)*Z(4)+C(58)*Y(2)+C(59)*Y(2)*Z(2)+C(60)*Y(2)*Z(4)+C(61)*Y(4)+
     * C(62)*Y(4)*Z(2)+C(63)*Y(4)*Z(4)) 

       RETURN 
       END
C-----------------------------------------------------------------------
      SUBROUTINE BFTAIL(XX, BB)
C  VERSION 5/13/88
C  DEVELOPED MCDONNELL DOUGLAS
C  FOR INFORMATION CALL KARL PFITZER (714) 896-3231
C
C     INPUT -- ARGUMENT LIST
C        XX     A REAL ARRAY GIVING THE POSITION WHERE THE MAGNETIC FIELD 
C               IS TO BE DETERMINED. XX(1), XX(2), XX(3) ARE RESPECTIVELY 
C               THE X, Y, AND Z SOLAR MAGNETIC COORDINATES IN EARTH RADII.
C               Z IS ALONG THE EARTHS NORTH DIPOLE AXIS. X IS PERPENDICULAR
C               TO Z AND IN THE PLANE CONTAINING THE Z AXIS AND THE SUN-EARTH  
C               LINE (X IS POSITIVE IN THE SOLAR DIRECTION). Y IS PERPENDICULAR 
C               TO X AND Z AND X Y Z FORM A RIGHT HANDED COORDINATE SYSTEM. 
C     
C     OUTPUT -- ARGUMENT LIST 
C        BB     A REAL ARRAY CONTAINING THE X, Y, AND Z COMPONENTS OF THE 
C               EARTH'S TOTAL MAGNETIC FIELD IN SOLAR MAGNETIC COORDINATES.  
C               BB(1), BB(2) AND BB(3) ARE THE BX, BY, AND BZ COMPONENTS.
C               THE UNITS ARE NANOTESLA.
C
C  THIS IS THE EXPANSION FOR THE TAIL CURRENT SYSTEM. THE EXPANSION IS 
C  VALID FROM THE SUBSOLAR POIN TO -60 RE. THE EXPANSION IS BASED ON A 
C  FIT TO VALUES CALCULATED USING THE WIRE LOOP TAIL SYSTEM.

      IMPLICIT NONE

C Input
      REAL*8 XX(3)
C Output
      REAL*8 BB(3)
C Local
      REAL*8 A(65), B(17), C(39), X(5), Y(5), Z(5)
      REAL*8 XN, YN, ZN, R, R2, R22, EXPC, TANZR, EXPR
      INTEGER*4 I

      DATA A/ 
     *-.118386794d-12, .260137167d+01, .408016277d-12,-.306063863d+00,
     * .852659791d-13, .848404600d-14,-.568097241d-02,-.601368497d-14,
     *-.336276159d-13,-.676779936d-15,-.110762251d-02,-.150912058d-15,
     *-.477506548d-14,-.805245718d-02,-.130105300d-14, .442299435d-16,
     *-.432185140d-04,-.520612496d-01,-.918209408d-04,-.686114562d-03,
     * .275041492d-04, .235864029d-15,-.628394374d-04,-.236539414d-16,
     * .379298441d-18,-.109452698d-14,-.163675727d-16,-.766199004d-04,
     *-.110519916d-15,-.111417355d-17, .311215382d-17,-.605957952d-06,
     *-.609414361d+01, .207037106d-12, .130315144d+00,-.250115110d-13,
     * .325228977d+00, .169606672d-01,-.131084126d-14, .232305257d-03,
     * .254138418d-01,-.585580678d-03, .344211139d-16, .268904941d-05,
     * .561115936d-01,-.855121118d-15, .577135898d-03,-.389637036d-04,
     * .531094438d-18, .517250317d-14, .163439821d-17, .280008382d-15,
     * .311491125d-17, .165293989d-02,-.149174308d-16, .406457779d-05,
     *-.415855886d-06, .127866736d-03,-.106070848d-04, .105524883d-17,
     * .293942950d-05,-.417367450d-06, .134032750d-04,-.139506296d-18,
     * 0.0d0/ 

      DATA B/ 
     *-.323149328d-01, .430535014d-02, .115661689d-03,-.486002660d-04,
     *-.102777234d-04,-.489864422d-05,-.356884232d-04,-.334316125d-07,
     * .122456608d+00, .202317315d-01,-.487990709d-03, .338684854d-04,
     *-.511755985d-04, .119096933d-04, .609353153d-03,-.243627124d-05,
     * 0.0d0/ 

      DATA C/ 
     * .318422091d+00, .154017442d+00, .337581827d-01, .436882397d-01,
     *-.153732787d-03, .362817457d-02, .179382198d-03,-.394772816d-05,
     *-.193942567d-01,-.263603775d-04,-.314364082d-04,-.103110548d-02,
     * .386165884d-06,-.301272556d-06,-.102838611d-03,-.725608973d-04,
     *-.893564810d-05,-.200670765d-05,-.805631807d-05,-.217861072d+02,
     *-.219688864d+01, .178558432d+00, .144137907d-01,-.293171667d-04,
     * .178727330d-01, .846703874d-02, .292860242d-04,-.583591628d+00,
     * .177991433d-02, .253212943d-02,-.629907297d-01, .669977751d-04,
     * .141706101d-03,-.334067698d-03, .122648694d-03,-.259383966d-07,
     * .252027517d-04,-.212223753d-02,
     * 0.0d0/ 

      XN = XX(1)
      YN = XX(2)
      ZN = XX(3)
      R2 = XN * XN + YN * YN + ZN * ZN
      R = DSQRT(R2)
      DO 1 I=1,5  
        X(I) = XN 
        Y(I) = YN 
        Z(I) = ZN 
        XN = XN * XX(1) 
        YN = YN * XX(2) 
        ZN = ZN * XX(3) 
    1 CONTINUE
      R22 = DSQRT((22.0D0-X(1))**2+Y(2)+Z(2)) 
      EXPC = DEXP(X(1)/15.0D0)
      TANZR = DTANH(Z(1)) * (1.0D0-DTANH((8.0D0-R)/5.0D0)) 
      EXPR = DEXP(-(R22-29.0D0)**2/60.0D0)
      BB(1)=(A( 2)*Z( 1)+A( 4)*X( 1)*Z( 1)
     *+A( 7)*Y( 2)*Z( 1)
     *+A(11)*X( 1)*Y( 2)*Z( 1)
     *+A(17)*X( 2)*Y( 2)*Z( 1)+A(18)*Z( 3)+A(19)*Y( 2)*Z( 3)+A(20)*X( 1)
     **Z( 3)+A(21)*X( 2)*Z( 3)+A(23)*X( 3)*Z( 1)
     *+A(28)*Y( 4)*Z( 1)
     *+A(32)*X( 4)*Z( 1))*EXPC
     *+(0.0D0 +A(33)+A(35)*X( 1)+A(37)*Z( 2)
     *+A(38)*Y( 2)+A(40)*Y( 2)*Z( 2)+A(41)*X( 1)*Z( 2)
     *+A(42)*X( 1)*Y( 2)+A(44)*X( 1)*Y( 2)*Z( 2)
     *+A(45)*X( 2)+A(47)*X( 2)*Z( 2)+A(48)*X( 2)*Y( 2)
     *+A(54)*X( 3)+A(56)*X( 3)
     **Z( 2)+A(57)*X( 3)*Y( 2)+A(58)*Z( 4)+A(59)*Y( 4)
     *+A(61)*X( 1)*Z( 4)+A(62)*X( 1)*Y( 4)+A(63)*X( 4))*TANZR 
      BB(2)=(B( 1)*Y( 1)*Z( 1)+B( 2)*X( 1)*Y( 1)*Z( 1)+B( 3)*Y( 1)*Z( 3)
     *+B( 4)*Y( 3)*Z( 1)+B( 5)*X( 1)*Y( 1)*Z( 3)+B( 6)*X( 1)*Y( 3)*Z( 1)
     *+B( 7)*X( 2)*Y( 1)*Z( 1)+B( 8)*X( 3)*Y( 1)*Z( 1))*EXPC
     *+(0.0D0 +B( 9)*Y( 1)
     **Z( 1)+B(10)*X( 1)*Y( 1)*Z( 1)+B(11)*Y( 1)*Z( 3)+B(12)*Y( 3)*Z( 1)
     *+B(13)*X( 1)*Y( 1)*Z( 3)+B(14)*X( 1)*Y( 3)*Z( 1)+B(15)*X( 2)*Y( 1)
     **Z( 1)+B(16)*X( 3)*Y( 1)*Z( 1))*EXPR
      BB(3)=(C( 1)+C( 2)*X( 1)+C( 3)*Z( 2)+C( 4)*Y( 2)+C( 5)*Y( 2)*Z( 2)
     *+C( 6)*X( 1)*Z( 2)+C( 7)*X( 1)*Y( 2)+C( 8)*X( 1)*Y( 2)*Z( 2)+C( 9)
     **X( 2)+C(10)*X( 2)*Z( 2)+C(11)*X( 2)*Y( 2)+C(12)*X( 3)+C(13)*X( 3)
     **Z( 2)+C(14)*X( 3)*Y( 2)+C(15)*Z( 4)+C(16)*Y( 4)+C(17)*X( 1)*Z( 4)
     *+C(18)*X( 1)*Y( 4)+C(19)*X( 4))*EXPC
     *+(0.0D0 +C(20)+C(21)*X( 1)+C(22)*Z( 2)
     *+C(23)*Y( 2)+C(24)*Y( 2)*Z( 2)+C(25)*X( 1)*Z( 2)+C(26)*X( 1)*Y( 2)
     *+C(27)*X( 1)*Y( 2)*Z( 2)+C(28)*X( 2)+C(29)*X( 2)*Z( 2)+C(30)*X( 2)
     **Y( 2)+C(31)*X( 3)+C(32)*X( 3)*Z( 2)+C(33)*X( 3)*Y( 2)+C(34)*Z( 4)
     *+C(35)*Y( 4)+C(36)*X( 1)*Z( 4)+C(37)*X( 1)*Y( 4)+C(38)*X( 4))*EXPR

      RETURN
      END 
C-----------------------------------------------------------------------
      SUBROUTINE BFRING(XX, BB)
C  VERSION 5/13/88
C  DEVELOPED MCDONNELL DOUGLAS
C  FOR INFORMATION CALL KARL PFITZER (714) 896-3231
C
C     INPUT -- ARGUMENT LIST
C        XX     A REAL ARRAY GIVING THE POSITION WHERE THE MAGNETIC FIELD 
C               IS TO BE DETERMINED. XX(1), XX(2), XX(3) ARE RESPECTIVELY 
C               THE X, Y, AND Z SOLAR MAGNETIC COORDINATES IN EARTH RADII.
C               Z IS ALONG THE EARTHS NORTH DIPOLE AXIS. X IS PERPENDICULAR
C               TO Z AND IN THE PLANE CONTAINING THE Z AXIS AND THE SUN-EARTH 
C               LINE (X IS POSITIVE IN THE SOLAR DIRECTION). Y IS PERPENDICULAR 
C               TO X AND Z AND X,Y,Z FORM A RIGHT HANDED COORDINATE SYSTEM. 
C     
C     OUTPUT -- ARGUMENT LIST 
C        BB     A REAL ARRAY CONTAINING THE X, Y, AND Z COMPONENTS OF THE 
C               EARTH'S TOTAL MAGNETIC FIELD IN SOLAR MAGNETIC COORDINATES.  
C               BB(1), BB(2), AND BB(3) ARE THE BX, BY, AND BZ COMPONENTS.
C               THE UNITS ARE NANOTESLA.
C
C  THIS SUBROUTINE CALCULATES THE FIELD FROM THE RING CURRENT SYSTEM.
C  THE EXPANSION IS A FIT TO VALUES CALCULATED FROM THE WIRE RING
C  CURRENT MODEL. THE EXPANSION IS VALID FROM THE SUBSOLAR POINT.
C  TO -60 RE.

      IMPLICIT NONE

C Input
      REAL*8 XX(3)
C Output
      REAL*8 BB(3)
C Local
      REAL*8 A(29), B(17), C(39), X(5), Y(5), Z(5)
      REAL*8 XN, YN, ZN, R, R2, EXPC, EXPR
      INTEGER*4 I

      DATA A/ 
     * .937029737d+00,-.734269078d+00,-.125896726d-01,-.843388063d-02,
     * .756104711d-04, .294507011d-02,-.719118601d-03,-.177154663d-01,
     * .104113319d-03,-.339745485d-04, .324439655d-03, .492786378d-04,
     *-.100821105d-04, .109966887d-04, .119616338d+00, .403556177d+01,
     *-.363651494d-01,-.337286459d-01,-.908902973d-04,-.980450316d-01,
     *-.220988518d+00,-.244671475d+00,-.977974501d-03, .311933785d-01,
     *-.249204900d+00, .825058070d-03, .464195892d-02, .223651513d-01,
     * 0.0d0/ 

      DATA B/ 
     *-.908641389d+00,-.249680217d-01, .443512048d-02,-.124215709d-03,
     * .211679921d-03,-.368134800d-04, .547288643d-03, .164845371d-04,
     * .407818714d+01,-.129156231d+00,-.940633654d-01,-.220684438d+00,
     * .878070158d-04, .174193445d-01,-.223040987d+00, .151981648d-01,
     * 0.0d0/ 

      DATA C/ 
     *-.381390073d+02,-.362173083d+01,-.410551306d+00, .532760526d+00,
     *-.151227645d-02, .182345800d-01, .358417761d-01,-.103889316d-03,
     * .395514004d+00, .100299786d-02, .138275245d-03, .288046807d-01,
     *-.127951613d-05,-.177797800d-04, .239511803d-02,-.284121147d-03,
     * .939796129d-04,-.101830861d-04, .504629929d-03, .105982946d+02,
     * .265464860d+01,-.157855689d+01,-.548140707d+01,-.181759612d-01,
     * .653535097d-01, .405331254d+00,-.726064092d-02,-.554702622d+01,
     *-.652021402d-02, .802389538d-01, .167926792d+00,-.384118806d-02,
     * .872021714d-02, .474604567d-01, .772720393d-01, .144274860d-02,
     *-.179837707d-01, .871619151d-01,
     * 0.0d0/ 
      XN = XX(1)
      YN = XX(2)
      ZN = XX(3)
      R2 = XN * XN + YN * YN + ZN * ZN
      R = DSQRT(R2)
      DO I=1,5  
        X(I) = XN 
        Y(I) = YN 
        Z(I) = ZN 
        XN = XN * XX(1) 
        YN = YN * XX(2) 
        ZN = ZN * XX(3) 
      END DO
      EXPC = DEXP(-R/5.20D0)
      IF (R2 .GT. 900.0D0) R2 = 900.0D0 
      EXPR = DEXP(-0.060D0*R2) 
      BB(1)=(A( 1)*Z( 1)+A( 2)*X( 1)*Z( 1)+A( 3)*Z( 3)+A( 4)*Y( 2)*Z( 1)
     *+A( 5)*Y( 2)*Z( 3)+A( 6)*X( 1)*Z( 3)+A( 7)*X( 1)*Y( 2)*Z( 1)+A( 8)
     **X( 2)*Z( 1)+A( 9)*X( 2)*Z( 3)+A(10)*X( 2)*Y( 2)*Z( 1)+A(11)*X( 3)
     **Z( 1)+A(12)*Z( 5)+A(13)*Y( 4)*Z( 1)+A(14)*X( 4)*Z( 1))*EXPC
     *+(0.0 +A(15)
     **Z( 1)+A(16)*X( 1)*Z( 1)+A(17)*Z( 3)+A(18)*Y( 2)*Z( 1)+A(19)*Y( 2)
     **Z( 3)+A(20)*X( 1)*Z( 3)+A(21)*X( 1)*Y( 2)*Z( 1)+A(22)*X( 2)*Z( 1)
     *+A(23)*X( 2)*Z( 3)+A(24)*X( 2)*Y( 2)*Z( 1)+A(25)*X( 3)*Z( 1)+A(26)
     **Z( 5)+A(27)*Y( 4)*Z( 1)+A(28)*X( 4)*Z( 1))*EXPR
      BB(2)=(B( 1)*Y( 1)*Z( 1)+B( 2)*X( 1)*Y( 1)*Z( 1)+B( 3)*Y( 1)*Z( 3)
     *+B( 4)*Y( 3)*Z( 1)+B( 5)*X( 1)*Y( 1)*Z( 3)+B( 6)*X( 1)*Y( 3)*Z( 1)
     *+B( 7)*X( 2)*Y( 1)*Z( 1)+B( 8)*X( 3)*Y( 1)*Z( 1))*EXPC
     *+(0.0 +B( 9)*Y( 1)
     **Z( 1)+B(10)*X( 1)*Y( 1)*Z( 1)+B(11)*Y( 1)*Z( 3)+B(12)*Y( 3)*Z( 1)
     *+B(13)*X( 1)*Y( 1)*Z( 3)+B(14)*X( 1)*Y( 3)*Z( 1)+B(15)*X( 2)*Y( 1)
     **Z( 1)+B(16)*X( 3)*Y( 1)*Z( 1))*EXPR
      BB(3)=(C( 1)+C( 2)*X( 1)+C( 3)*Z( 2)+C( 4)*Y( 2)+C( 5)*Y( 2)*Z( 2)
     *+C( 6)*X( 1)*Z( 2)+C( 7)*X( 1)*Y( 2)+C( 8)*X( 1)*Y( 2)*Z( 2)+C( 9)
     **X( 2)+C(10)*X( 2)*Z( 2)+C(11)*X( 2)*Y( 2)+C(12)*X( 3)+C(13)*X( 3)
     **Z( 2)+C(14)*X( 3)*Y( 2)+C(15)*Z( 4)+C(16)*Y( 4)+C(17)*X( 1)*Z( 4)
     *+C(18)*X( 1)*Y( 4)+C(19)*X( 4))*EXPC
     *+(0.0 +C(20)+C(21)*X( 1)+C(22)*Z( 2)
     *+C(23)*Y( 2)+C(24)*Y( 2)*Z( 2)+C(25)*X( 1)*Z( 2)+C(26)*X( 1)*Y( 2)
     *+C(27)*X( 1)*Y( 2)*Z( 2)+C(28)*X( 2)+C(29)*X( 2)*Z( 2)+C(30)*X( 2)
     **Y( 2)+C(31)*X( 3)+C(32)*X( 3)*Z( 2)+C(33)*X( 3)*Y( 2)+C(34)*Z( 4)
     *+C(35)*Y( 4)+C(36)*X( 1)*Z( 4)+C(37)*X( 1)*Y( 4)+C(38)*X( 4))*EXPR

      RETURN
      END 
C-----------------------------------------------------------------------
      REAL*8 FUNCTION RINGST(SOFFD, DST)
C
C  THIS FUNCTION CALCULATES THE STRENGTH OF THE RING CURRENT FROM THE 
C  STANDOFF DISTANCE AND THE DST.
C
C  THIS FUNCTION CAN BE USED TO CALCULATE A VALUE FOR SRING, ONE OF THE
C  REQUIRED PARAMETERS FOR CALCUATING THE DYNAMIC MAGNETIC FIELD.
C
C  IT CALCULATES THE CONTRIBUTION OF THE MAGNETOPAUSE CURRENTS TO GROUND 
C  BASED SIGNATURE AND SUBTRACTS THAT COMPONENT FROM THE OBSERVED VALUE 
C  OF DST. IT ATTRIBUTES THE REMAINDER TO THE RING CURRENT.
C
C  INPUT PARAMETERS
C
C        SOFFD  THE STANDOFF DISTANCE OF THE MAGNETOPAUSE. THE QUIET STANDOFF 
C               DISTANCE IS 10.5 EARTH RADII. ACCEPTABLE VALUES RANGE BETWEEN 
C               6 AND 11. THIS VALUE IS USED TO CALCULATE THE STRENGTH OF
C               THE MAGNETOPAUSE CURRENTS AND TO SCALE THE SIZE OF THE 
C               MAGNETOPAUSE. THIS VALUE ALSO SCALES THE SIZE OF THE TAIL 
C               CURRENT SYSTEM. THE RING SYSTEM IS NOT SCALED, SINCE ITS 
C               SOURCE IS PRIMARILY AT RADIAL DISTANCES.
C
C        DST    DST IS THE STANDARD PUBLISHED DST VALUE IN NANOTESLA, THE 
C               STORMTIME DISTURBED EQUATORIAL FIELD.
C
C  CON    SCALES THE EFFECT OF THE DST (ITS VALUE OF .03 IS STILL
C         SOMEWHAT UNCERTAIN).
C
      IMPLICIT NONE

C Input
      REAL*8 SOFFD, DST
C Local
      REAL*8 SCL, SCM, DSTMOD, CON / 0.03D0 /

      SCL = 10.5D0 / SOFFD
      SCM = SCL**3.0D0
      DSTMOD = (SCM-1.0D0) * 15.0D0 - DST
      RINGST = 1.0D0 + DSTMOD * CON

      RETURN
      END
C-----------------------------------------------------------------------
      REAL*8 FUNCTION STDOFF(VEL, DEN)
C
C  THIS FUNCTION CALCULATES THE STANOFF DISTANCE FROM THE SOLAR WIND VELOCITY 
C  AND DENSITY. IT CAN BE USED TO EVALUATE THE PARAMETER SOFFD. SOFFD IS 
C  REQUIRED FOR ALL SCALING OPERATIONS.
C
C  INPUT PARAMETERS
C
C     VEL    THE SOLAR WIND VELOCITY IN KM/SEC NEAR THE SUBSOLAR POINT.
C            TYPICAL VALUES ARE 300 TO 500.
C
C     DEN    THE NUMBER DENSITY OF THE SOLAR WIND IN NUMBER PER CC.
C            TYPICAL VALUES ARE 5 TO 50.
C
C  OUTPUT
C
C     STDOFF THE DISTANCE TO THE SUBSOLAR POINT IN RE.
C
      IMPLICIT NONE

C Input
      REAL*8 VEL, DEN

      STDOFF = 98.0D0 / ((DEN*VEL**2)**(1.0D0/6.0D0))

      RETURN
      END
C-----------------------------------------------------------------------
C----------------------------------------------------------------------
c
      SUBROUTINE T96_01 (IOPT,PARMOD,PS,X,Y,Z,BX,BY,BZ)
C
c     RELEASE DATE OF THIS VERSION:   JUNE 22, 1996.
C----------------------------------------------------------------------
C
C  WITH TWO CORRECTIONS, SUGGESTED BY T.SOTIRELIS' COMMENTS (APR.7, 1997)
C
C  (1) A "STRAY "  CLOSING PARENTHESIS WAS REMOVED IN THE S/R   R2_BIRK
C  (2) A 0/0 PROBLEM ON THE Z-AXIS WAS SIDESTEPPED (LINES 44-46 OF THE
c       DOUBLE PRECISION FUNCTION XKSI)
c--------------------------------------------------------------------
C DATA-BASED MODEL CALIBRATED BY (1) SOLAR WIND PRESSURE PDYN (NANOPASCALS),
C           (2) DST (NANOTESLA),  (3) BYIMF, AND (4) BZIMF (NANOTESLA).
c THESE INPUT PARAMETERS SHOULD BE PLACED IN THE FIRST 4 ELEMENTS
c OF THE ARRAY PARMOD(10).
C
C   THE REST OF THE INPUT VARIABLES ARE: THE GEODIPOLE TILT ANGLE PS (RADIANS),
C AND   X,Y,Z -  GSM POSITION (RE)
C
c   IOPT  IS JUST A DUMMY INPUT PARAMETER, NECESSARY TO MAKE THIS SUBROUTINE
C COMPATIBLE WITH THE NEW RELEASE (APRIL 1996) OF THE TRACING SOFTWARE
C PACKAGE (GEOPACK). IOPT VALUE DOES NOT AFFECT THE OUTPUT FIELD.
c
C
c OUTPUT:  GSM COMPONENTS OF THE EXTERNAL MAGNETIC FIELD (BX,BY,BZ, nanotesla)
C            COMPUTED AS A SUM OF CONTRIBUTIONS FROM PRINCIPAL FIELD SOURCES
C
c  (C) Copr. 1995, 1996, Nikolai A. Tsyganenko, Hughes STX, Code 695, NASA GSFC
c      Greenbelt, MD 20771, USA
c
C                            REFERENCES:
C
C               (1) N.A. TSYGANENKO AND D.P. STERN, A NEW-GENERATION GLOBAL
C           MAGNETOSPHERE FIELD MODEL  , BASED ON SPACECRAFT MAGNETOMETER DATA,
C           ISTP NEWSLETTER, V.6, NO.1, P.21, FEB.1996.
C
c              (2) N.A.TSYGANENKO,  MODELING THE EARTH'S MAGNETOSPHERIC
C           MAGNETIC FIELD CONFINED WITHIN A REALISTIC MAGNETOPAUSE,
C           J.GEOPHYS.RES., V.100, P. 5599, 1995.
C
C              (3) N.A. TSYGANENKO AND M.PEREDO, ANALYTICAL MODELS OF THE
C           MAGNETIC FIELD OF DISK-SHAPED CURRENT SHEETS, J.GEOPHYS.RES.,
C           V.99, P. 199, 1994.
C
c----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 PDYN,DST,BYIMF,BZIMF,PS,X,Y,Z,BX,BY,BZ,QX,QY,QZ,PARMOD(10),
     *   A(9)
c
      DATA PDYN0,EPS10 /2.d0,3630.7d0/
C
      DATA A/1.162d0,22.344d0,18.50d0,2.602d0,6.903d0,5.287d0,0.5790d0,
     :       0.4462d0,0.7850d0/
C
      DATA  AM0,S0,X00,DSIG/70.d0,1.08d0,5.48d0,0.005d0/
      DATA  DELIMFX,DELIMFY /20.d0,10.d0/
C
       idummy=iopt
       PDYN=PARMOD(1)
       DST=PARMOD(2)
       BYIMF=PARMOD(3)
       BZIMF=PARMOD(4)
C
       SPS=SIN(PS)
       PPS=PS
C
       DEPR=0.8*DST-13.*SQRT(PDYN)  
c!  DEPR is an estimate of total near-Earth
c                                         depression, based on DST and Pdyn
c                                             (usually, DEPR &lt 0 )
C
C  CALCULATE THE IMF-RELATED QUANTITIES:
C
       Bt=SQRT(BYIMF**2+BZIMF**2)

       IF (BYIMF.EQ.0.d0.AND.BZIMF.EQ.0.d0) THEN
          THETA=0.
          GOTO 1
       ENDIF
C
       THETA=ATAN2(BYIMF,BZIMF)
       IF (THETA.LE.0.D0) THETA=THETA+atan(1.0d0)*8
  1    CT=COS(THETA)
       ST=SIN(THETA)
       EPS=718.5d0*SQRT(Pdyn)*Bt*SIN(THETA/2)
C
       FACTEPS=EPS/EPS10-1
       FACTPD=SQRT(PDYN/PDYN0)-1
C
       RCAMPL=-A(1)*DEPR     
c!   RCAMPL is the amplitude of the ring current
c                  (positive and equal to abs.value of RC depression at origin)
C
       TAMPL2=A(2)+A(3)*FACTPD+A(4)*FACTEPS
       TAMPL3=A(5)+A(6)*FACTPD
       B1AMPL=A(7)+A(8)*FACTEPS
       B2AMPL=20.*B1AMPL  
c! IT IS EQUIVALENT TO ASSUMING THAT THE TOTAL CURRENT
C                           IN THE REGION 2 SYSTEM IS 40% OF THAT IN REGION 1
       RECONN=A(9)
C
       XAPPA=(PDYN/PDYN0)**0.14d0
       XAPPA3=XAPPA**3
       YS=Y*CT-Z*ST
       ZS=Z*CT+Y*ST
C
       FACTIMF=EXP(X/DELIMFX-(YS/DELIMFY)**2)
C
C  CALCULATE THE "IMF" COMPONENTS OUTSIDE THE LAYER  (HENCE BEGIN WITH "O")
C
       OIMFX=0
       OIMFY=RECONN*BYIMF*FACTIMF
       OIMFZ=RECONN*BZIMF*FACTIMF
C
       RIMFAMPL=RECONN*Bt
C
       PPS=PS
       XX=X*XAPPA
       YY=Y*XAPPA
       ZZ=Z*XAPPA
C
C  SCALE AND CALCULATE THE MAGNETOPAUSE PARAMETERS FOR THE INTERPOLATION ACROSS
C   THE BOUNDARY LAYER (THE COORDINATES XX,YY,ZZ  ARE ALREADY SCALED)
C
       X0=X00/XAPPA
       AM=AM0/XAPPA
       RHO2=Y**2+Z**2
       ASQ=AM**2
       XMXM=AM+X-X0
       IF (XMXM.LT.0) XMXM=0 
c! THE BOUNDARY IS A CYLINDER TAILWARD OF X=X0-AM
       AXX0=XMXM**2
       ARO=ASQ+RHO2
       SIGMA=SQRT((ARO+AXX0+SQRT((ARO+AXX0)**2-4*ASQ*AXX0))/(2*ASQ))
C
C   NOW, THERE ARE THREE POSSIBLE CASES:
C    (1) INSIDE THE MAGNETOSPHERE
C    (2) IN THE BOUNDARY LAYER
C    (3) OUTSIDE THE MAGNETOSPHERE AND B.LAYER
C       FIRST OF ALL, CONSIDER THE CASES (1) AND (2):
C
      IF (SIGMA.LT.S0+DSIG) THEN  
c!  CALCULATE THE T95_06 FIELD (WITH THE
C                                POTENTIAL "PENETRATED" INTERCONNECTION FIELD):

       CALL DIPSHLD(PPS,XX,YY,ZZ,CFX,CFY,CFZ)
       CALL TAILRC96(SPS,XX,YY,ZZ,BXRC,BYRC,BZRC,BXT2,BYT2,BZT2,
     *   BXT3,BYT3,BZT3)
       CALL BIRK1TOT_02(PPS,XX,YY,ZZ,R1X,R1Y,R1Z)
       CALL BIRK2TOT_02(PPS,XX,YY,ZZ,R2X,R2Y,R2Z)
       CALL INTERCON(XX,YS*XAPPA,ZS*XAPPA,RIMFX,RIMFYS,RIMFZS)
       RIMFY=RIMFYS*CT+RIMFZS*ST
       RIMFZ=RIMFZS*CT-RIMFYS*ST
C
      FX=CFX*XAPPA3+RCAMPL*BXRC +TAMPL2*BXT2+TAMPL3*BXT3
     *  +B1AMPL*R1X +B2AMPL*R2X +RIMFAMPL*RIMFX
      FY=CFY*XAPPA3+RCAMPL*BYRC +TAMPL2*BYT2+TAMPL3*BYT3
     *  +B1AMPL*R1Y +B2AMPL*R2Y +RIMFAMPL*RIMFY
      FZ=CFZ*XAPPA3+RCAMPL*BZRC +TAMPL2*BZT2+TAMPL3*BZT3
     *  +B1AMPL*R1Z +B2AMPL*R2Z +RIMFAMPL*RIMFZ
C
C  NOW, LET US CHECK WHETHER WE HAVE THE CASE (1). IF YES - WE ARE DONE:
C
       IF (SIGMA.LT.S0-DSIG) THEN
         BX=FX
         BY=FY
         BZ=FZ
                        ELSE  
c!  THIS IS THE MOST COMPLEX CASE: WE ARE INSIDE
C                                         THE INTERPOLATION REGION
       FINT=0.5d0*(1-(SIGMA-S0)/DSIG)
       FEXT=0.5d0*(1+(SIGMA-S0)/DSIG)
C
       CALL DIPOLE_96(PS,X,Y,Z,QX,QY,QZ)
       BX=(FX+QX)*FINT+OIMFX*FEXT -QX
       BY=(FY+QY)*FINT+OIMFY*FEXT -QY
       BZ=(FZ+QZ)*FINT+OIMFZ*FEXT -QZ
c
        ENDIF  
c!   THE CASES (1) AND (2) ARE EXHAUSTED; THE ONLY REMAINING
C                      POSSIBILITY IS NOW THE CASE (3):
         ELSE
                CALL DIPOLE_96(PS,X,Y,Z,QX,QY,QZ)
                BX=OIMFX-QX
                BY=OIMFY-QY
                BZ=OIMFZ-QZ
         ENDIF
C
       RETURN
       END
C=====================================================================

      SUBROUTINE DIPSHLD(PS,X,Y,Z,BX,BY,BZ)
C
C   CALCULATES GSM COMPONENTS OF THE EXTERNAL MAGNETIC FIELD DUE TO
C    SHIELDING OF THE EARTH'S DIPOLE ONLY
C
       IMPLICIT REAL*8 (A-H,O-Z)
       DIMENSION A1(12),A2(12)
      DATA A1 /.24777d0,-27.003d0,-.46815d0,7.0637d0,-1.5918d0,
     :   -.90317d-01,57.522d0,13.757d0,2.0100d0,10.458d0,4.5798d0,
     :   2.1695d0/
      DATA A2/-.65385d0,-18.061d0,-.40457d0,-5.0995d0,1.2846d0,
     :   .78231d-01,39.592d0,13.291d0,1.9970d0,10.062d0,4.5140d0,
     :   2.1558d0/
C
          CPS=DCOS(PS)
          SPS=DSIN(PS)
          CALL CYLHARM(A1,X,Y,Z,HX,HY,HZ)
          CALL CYLHAR1(A2,X,Y,Z,FX,FY,FZ)
C
          BX=HX*CPS+FX*SPS
          BY=HY*CPS+FY*SPS
          BZ=HZ*CPS+FZ*SPS
        RETURN
       END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C  THIS CODE YIELDS THE SHIELDING FIELD FOR THE PERPENDICULAR DIPOLE
C
         SUBROUTINE  CYLHARM( A, X,Y,Z, BX,BY,BZ)
C
C
C   ***  N.A. Tsyganenko ***  Sept. 14-18, 1993; revised March 16, 1994 ***
C
C   An approximation for the Chapman-Ferraro field by a sum of 6 cylin-
c   drical harmonics (see pp. 97-113 in the brown GSFC notebook #1)
c
C      Description of parameters:
C
C  A   - input vector containing model parameters;
C  X,Y,Z   -  input GSM coordinates
C  BX,BY,BZ - output GSM components of the shielding field
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  The 6 linear parameters A(1)-A(6) are amplitudes of the cylindrical harmonic
c       terms.
c  The 6 nonlinear parameters A(7)-A(12) are the corresponding scale lengths
C       for each term (see GSFC brown notebook).
c
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
       IMPLICIT  REAL * 8  (A - H, O - Z)
C
       DIMENSION  A(12)
C
           RHO=DSQRT(Y**2+Z**2)
            IF (RHO.LT.1.D-8) THEN
               SINFI=1.D0
               COSFI=0.D0
               RHO=1.D-8
               GOTO 1
            ENDIF
C
           SINFI=Z/RHO
           COSFI=Y/RHO
  1        SINFI2=SINFI**2
           SI2CO2=SINFI2-COSFI**2
C
             BX=0.D0
             BY=0.D0
             BZ=0.D0
C
         DO 11 I=1,3
             DZETA=RHO/A(I+6)
             XJ0=BES(DZETA,0)
             XJ1=BES(DZETA,1)
             XEXP=DEXP(X/A(I+6))
             BX=BX-A(I)*XJ1*XEXP*SINFI
             BY=BY+A(I)*(2.D0*XJ1/DZETA-XJ0)*XEXP*SINFI*COSFI
             BZ=BZ+A(I)*(XJ1/DZETA*SI2CO2-XJ0*SINFI2)*XEXP
   11        CONTINUE
c
         DO 12 I=4,6
             DZETA=RHO/A(I+6)
             XKSI=X/A(I+6)
             XJ0=BES(DZETA,0)
             XJ1=BES(DZETA,1)
             XEXP=DEXP(XKSI)
             BRHO=(XKSI*XJ0-(DZETA**2+XKSI-1.D0)*XJ1/DZETA)*XEXP*SINFI
             BPHI=(XJ0+XJ1/DZETA*(XKSI-1.D0))*XEXP*COSFI
             BX=BX+A(I)*(DZETA*XJ0+XKSI*XJ1)*XEXP*SINFI
             BY=BY+A(I)*(BRHO*COSFI-BPHI*SINFI)
             BZ=BZ+A(I)*(BRHO*SINFI+BPHI*COSFI)
   12        CONTINUE
C
c
         RETURN
       END
C
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C  THIS CODE YIELDS THE SHIELDING FIELD FOR THE PARALLEL DIPOLE
C
         SUBROUTINE  CYLHAR1(A, X,Y,Z, BX,BY,BZ)
C
C
C   ***  N.A. Tsyganenko ***  Sept. 14-18, 1993; revised March 16, 1994 ***
C
C   An approximation of the Chapman-Ferraro field by a sum of 6 cylin-
c   drical harmonics (see pages 97-113 in the brown GSFC notebook #1)
c
C      Description of parameters:
C
C  A   - input vector containing model parameters;
C  X,Y,Z - input GSM coordinates,
C  BX,BY,BZ - output GSM components of the shielding field
C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
C      The 6 linear parameters A(1)-A(6) are amplitudes of the cylindrical
c  harmonic terms.
c      The 6 nonlinear parameters A(7)-A(12) are the corresponding scale
c  lengths for each term (see GSFC brown notebook).
c
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
       IMPLICIT  REAL * 8  (A - H, O - Z)
C
       DIMENSION  A(12)
C
           RHO=DSQRT(Y**2+Z**2)
            IF (RHO.LT.1.D-10) THEN
               SINFI=1.D0
               COSFI=0.D0
               GOTO 1
            ENDIF
C
           SINFI=Z/RHO
           COSFI=Y/RHO
C
   1      BX=0.D0
          BY=0.D0
          BZ=0.D0
C
             DO 11 I=1,3
             DZETA=RHO/A(I+6)
             XKSI=X/A(I+6)
             XJ0=BES(DZETA,0)
             XJ1=BES(DZETA,1)
             XEXP=DEXP(XKSI)
             BRHO=XJ1*XEXP
             BX=BX-A(I)*XJ0*XEXP
             BY=BY+A(I)*BRHO*COSFI
             BZ=BZ+A(I)*BRHO*SINFI
   11        CONTINUE
c
         DO 12 I=4,6
             DZETA=RHO/A(I+6)
             XKSI=X/A(I+6)
             XJ0=BES(DZETA,0)
             XJ1=BES(DZETA,1)
             XEXP=DEXP(XKSI)
             BRHO=(DZETA*XJ0+XKSI*XJ1)*XEXP
             BX=BX+A(I)*(DZETA*XJ1-XJ0*(XKSI+1.D0))*XEXP
             BY=BY+A(I)*BRHO*COSFI
             BZ=BZ+A(I)*BRHO*SINFI
   12        CONTINUE
C
         RETURN
       END

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
      REAL*8 FUNCTION BES(X,K)
      IMPLICIT REAL*8 (A-H,O-Z)
C
      IF (K.EQ.0) THEN
        BES=BES0(X)
        RETURN
      ENDIF
C
      IF (K.EQ.1) THEN
        BES=BES1(X)
        RETURN
      ENDIF
C
      IF (X.EQ.0.D0) THEN
        BES=0.D0
        RETURN
      ENDIF
C
      G=2.D0/X
      IF (X.LE.K) GOTO 10
C
      N=1
      XJN=BES1(X)
      XJNM1=BES0(X)
C
  1   XJNP1=G*N*XJN-XJNM1
      N=N+1
      IF (N.LT.K) GOTO 2
      BES=XJNP1
      RETURN
C
 2    XJNM1=XJN
      XJN=XJNP1
      GOTO 1
C
 10   N=24
      XJN=1.D0
      XJNP1=0.D0
      SUM=0.D0
C
  3   IF (MOD(N,2).EQ.0) SUM=SUM+XJN
      XJNM1=G*N*XJN-XJNP1
      N=N-1
C
      XJNP1=XJN
      XJN=XJNM1
      IF (N.EQ.K) BES=XJN
C
      IF (DABS(XJN).GT.1.D5) THEN
        XJNP1=XJNP1*1.D-5
        XJN=XJN*1.D-5
        SUM=SUM*1.D-5
        IF (N.LE.K) BES=BES*1.D-5
      ENDIF
C
      IF (N.EQ.0) GOTO 4
      GOTO 3
C
  4   SUM=XJN+2.D0*SUM
      BES=BES/SUM
      RETURN
      END
c-------------------------------------------------------------------
c
       REAL*8 FUNCTION BES0(X)
C
        IMPLICIT REAL*8 (A-H,O-Z)
C
        IF (DABS(X).LT.3.D0) THEN
          X32=(X/3.D0)**2
          BES0=1.D0-X32*(2.2499997D0-X32*(1.2656208D0-X32*
     *    (0.3163866D0-X32*(0.0444479D0-X32*(0.0039444D0
     *     -X32*0.00021D0)))))
        ELSE
        XD3=3.D0/X
        F0=0.79788456D0-XD3*(0.00000077D0+XD3*(0.00552740D0+XD3*
     *   (0.00009512D0-XD3*(0.00137237D0-XD3*(0.00072805D0
     *   -XD3*0.00014476D0)))))
        T0=X-0.78539816D0-XD3*(0.04166397D0+XD3*(0.00003954D0-XD3*
     *   (0.00262573D0-XD3*(0.00054125D0+XD3*(0.00029333D0
     *   -XD3*0.00013558D0)))))
        BES0=F0/DSQRT(X)*DCOS(T0)
       ENDIF
       RETURN
       END
c
c--------------------------------------------------------------------------
c
       REAL*8 FUNCTION BES1(X)
C
        IMPLICIT REAL*8 (A-H,O-Z)
C
       IF (DABS(X).LT.3.D0) THEN
        X32=(X/3.D0)**2
        BES1XM1=0.5D0-X32*(0.56249985D0-X32*(0.21093573D0-X32*
     *  (0.03954289D0-X32*(0.00443319D0-X32*(0.00031761D0
     *  -X32*0.00001109D0)))))
         BES1=BES1XM1*X
       ELSE
        XD3=3.D0/X
        F1=0.79788456D0+XD3*(0.00000156D0+XD3*(0.01659667D0+XD3*
     *   (0.00017105D0-XD3*(0.00249511D0-XD3*(0.00113653D0
     *   -XD3*0.00020033D0)))))
        T1=X-2.35619449D0+XD3*(0.12499612D0+XD3*(0.0000565D0-XD3*
     *   (0.00637879D0-XD3*(0.00074348D0+XD3*(0.00079824D0
     *   -XD3*0.00029166D0)))))
        BES1=F1/DSQRT(X)*DCOS(T1)
       ENDIF
       RETURN
       END
C------------------------------------------------------------
C
         SUBROUTINE INTERCON(X,Y,Z,BX,BY,BZ)
C
C      Calculates the potential interconnection field inside the magnetosphere,
c  corresponding to  DELTA_X = 20Re and DELTA_Y = 10Re (NB#3, p.90, 6/6/1996).
C  The position (X,Y,Z) and field components BX,BY,BZ are given in the rotated
c   coordinate system, in which the Z-axis is always directed along the BzIMF
c   (i.e. rotated by the IMF clock angle Theta)
C   It is also assumed that the IMF Bt=1, so that the components should be
c     (i) multiplied by the actual Bt, and
c     (ii) transformed to standard GSM coords by rotating back around X axis
c              by the angle -Theta.
c
C      Description of parameters:
C
C     X,Y,Z -   GSM POSITION
C      BX,BY,BZ - INTERCONNECTION FIELD COMPONENTS INSIDE THE MAGNETOSPHERE
C        OF A STANDARD SIZE (TO TAKE INTO ACCOUNT EFFECTS OF PRESSURE CHANGES,
C         APPLY THE SCALING TRANSFORMATION)
C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
C     The 9 linear parameters are amplitudes of the "cartesian" harmonics
c     The 6 nonlinear parameters are the scales Pi and Ri entering
c    the arguments of exponents, sines, and cosines in the 9 "Cartesian"
c       harmonics (3+3)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
       IMPLICIT  REAL * 8  (A - H, O - Z)
C
        DIMENSION A(15),RP(3),RR(3),P(3),R(3)
C
      SAVE M,RP,RR
C
      DATA A/-8.411078731d0,5932254.951d0,-9073284.93d0,-11.68794634d0,
     * 6027598.824d0,-9218378.368d0,-6.508798398d0,-11824.42793d0,
     : 18015.66212d0,
     * 7.99754043d0,13.9669886d0,90.24475036d0,16.75728834d0,
     : 1015.645781d0,
     * 1553.493216d0/
C
        DATA M/0/
C
        IF (M.NE.0) GOTO 111
        M=1
C
         P(1)=A(10)
         P(2)=A(11)
         P(3)=A(12)
         R(1)=A(13)
         R(2)=A(14)
         R(3)=A(15)
C
C
           DO 11 I=1,3
            RP(I)=1.D0/P(I)
  11        RR(I)=1.D0/R(I)
C
  111   CONTINUE
C
            L=0
C
               BX=0.d0
               BY=0.d0
               BZ=0.d0
C
c        "PERPENDICULAR" KIND OF SYMMETRY ONLY
C
               DO 2 I=1,3
                  CYPI=DCOS(Y*RP(I))
                  SYPI=DSIN(Y*RP(I))
C
                DO 2 K=1,3
                   SZRK=DSIN(Z*RR(K))
                   CZRK=DCOS(Z*RR(K))
                     SQPR=DSQRT(RP(I)**2+RR(K)**2)
                      EPR=DEXP(X*SQPR)
C
                     HX=-SQPR*EPR*CYPI*SZRK
                     HY=RP(I)*EPR*SYPI*SZRK
                     HZ=-RR(K)*EPR*CYPI*CZRK
             L=L+1
c
          BX=BX+A(L)*HX
          BY=BY+A(L)*HY
          BZ=BZ+A(L)*HZ
  2   CONTINUE
C
      RETURN
      END

C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE TAILRC96(SPS,X,Y,Z,BXRC,BYRC,BZRC,BXT2,BYT2,BZT2,
     *   BXT3,BYT3,BZT3)
c
c  COMPUTES THE COMPONENTS OF THE FIELD OF THE MODEL RING CURRENT AND THREE
c                   TAIL MODES WITH UNIT AMPLITUDES
C      (FOR THE RING CURRENT, IT MEANS THE DISTURBANCE OF Bz=-1nT AT ORIGIN,
C   AND FOR THE TAIL MODES IT MEANS MAXIMAL BX JUST ABOVE THE SHEET EQUAL 1 nT.
C
         IMPLICIT REAL*8 (A-H,O-Z)
         DIMENSION ARC(48),ATAIL2(48),ATAIL3(48)
         COMMON /WARP/ CPSS,SPSS,DPSRR,RPS,WARP,D,XS,ZS,DXSX,DXSY,DXSZ,
     *   DZSX,DZSY,DZSZ,DZETAS,DDZETADX,DDZETADY,DDZETADZ,ZSWW
C
         DATA ARC/-3.087699646d0,3.516259114d0,18.81380577d0,
     :  -13.95772338d0,-5.497076303d0,0.1712890838d0,2.392629189d0,
     :  -2.728020808d0,-14.79349936d0,11.08738083d0,4.388174084d0,
     :  0.2492163197d-01,0.7030375685d0,-.7966023165d0,
     :  -3.835041334d0,2.642228681d0,-0.2405352424d0,-0.7297705678d0,
     :  -0.3680255045d0,0.1333685557d0,2.795140897d0,-1.078379954d0,
     :  0.8014028630d0,0.1245825565d0,0.6149982835d0,-0.2207267314d0,
     :  -4.424578723d0,1.730471572d0,-1.716313926d0,-0.2306302941d0,
     :  -0.2450342688d0,0.8617173961d-01,1.54697858d0,
     :  -0.6569391113d0,-0.6537525353d0,0.2079417515d0,12.75434981d0,
     :  11.37659788d0,636.4346279d0,1.752483754d0,3.604231143d0,
     :  12.83078674d0,7.412066636d0,9.434625736d0,676.7557193d0,
     :  1.701162737d0,3.580307144d0,14.64298662d0/
C
        DATA ATAIL2/.8747515218d0,-.9116821411d0,2.209365387d0,
     : -2.159059518d0,-7.059828867d0,5.924671028d0,-1.916935691d0,
     : 1.996707344d0,-3.877101873d0,3.947666061d0,11.38715899d0,
     : -8.343210833d0,1.194109867d0,-1.244316975d0,3.73895491d0,
     : -4.406522465d0,-20.66884863d0,3.020952989d0,.2189908481d0,
     : -.09942543549d0,-.927225562d0,.1555224669d0,.6994137909d0,
     : -.08111721003d0,-.7565493881d0,.4686588792d0,4.266058082d0,
     : -.3717470262d0,-3.920787807d0,.02298569870d0,.7039506341d0,
     : -.5498352719d0,-6.675140817d0,.8279283559d0,-2.234773608d0,
     : -1.622656137d0,5.187666221d0,6.802472048d0,39.13543412d0,
     : 2.784722096d0,6.979576616d0,25.71716760d0,4.495005873d0,
     : 8.068408272d0,93.47887103d0,4.158030104d0,9.313492566d0,
     :57.18240483d0/
C
        DATA ATAIL3/-19091.95061d0,-3011.613928d0,20582.16203d0,
     : 4242.918430d0,-2377.091102d0,-1504.820043d0,19884.04650d0,
     : 2725.150544d0,-21389.04845d0,-3990.475093d0,2401.610097d0,
     : 1548.171792d0,-946.5493963d0,490.1528941d0,986.9156625d0,
     : -489.3265930d0,-67.99278499d0,8.711175710d0,-45.15734260d0,
     : -10.76106500d0,210.7927312d0,11.41764141d0,-178.0262808d0,
     : .7558830028d0,339.3806753d0,9.904695974d0,69.50583193d0,
     : -118.0271581d0,22.85935896d0,45.91014857d0,-425.6607164d0,
     : 15.47250738d0,118.2988915d0,65.58594397d0,-201.4478068d0,
     : -14.57062940d0,19.69877970d0,20.30095680d0,86.45407420d0,
     : 22.50403727d0,23.41617329d0,48.48140573d0,24.61031329d0,
     : 123.5395974d0,223.5367692d0,39.50824342d0,65.83385762d0,
     : 266.2948657d0/
C
       DATA RH,DR,G,D0,DELTADY/9.d0,4.d0,10.d0,2.d0,10.d0/
C
C   TO ECONOMIZE THE CODE, WE FIRST CALCULATE COMMON VARIABLES, WHICH ARE
C      THE SAME FOR ALL MODES, AND PUT THEM IN THE COMMON-BLOCK /WARP/
C
       DR2=DR*DR
       C11=DSQRT((1.D0+RH)**2+DR2)
       C12=DSQRT((1.D0-RH)**2+DR2)
       C1=C11-C12
       SPSC1=SPS/C1
       RPS=0.5d0*(C11+C12)*SPS  
c!  THIS IS THE SHIFT OF OF THE SHEET WITH RESPECT
C                            TO GSM EQ.PLANE FOR THE 3RD (ASYMPTOTIC) TAIL MODE
C
        R=DSQRT(X*X+Y*Y+Z*Z)
        SQ1=DSQRT((R+RH)**2+DR2)
        SQ2=DSQRT((R-RH)**2+DR2)
        C=SQ1-SQ2
        CS=(R+RH)/SQ1-(R-RH)/SQ2
        SPSS=SPSC1/R*C
        CPSS=DSQRT(1.D0-SPSS**2)
        DPSRR=SPS/(R*R)*(CS*R-C)/DSQRT((R*C1)**2-(C*SPS)**2)
C
        WFAC=Y/(Y**4+1.D4)   
c!   WARPING
        W=WFAC*Y**3
        WS=4.D4*Y*WFAC**2
        WARP=G*SPS*W
        XS=X*CPSS-Z*SPSS
        ZSWW=Z*CPSS+X*SPSS  
c! "WW" MEANS "WITHOUT Y-Z WARPING" (IN X-Z ONLY)
        ZS=ZSWW +WARP

        DXSX=CPSS-X*ZSWW*DPSRR
        DXSY=-Y*ZSWW*DPSRR
        DXSZ=-SPSS-Z*ZSWW*DPSRR
        DZSX=SPSS+X*XS*DPSRR
        DZSY=XS*Y*DPSRR  +G*SPS*WS  
c!  THE LAST TERM IS FOR THE Y-Z WARP
        DZSZ=CPSS+XS*Z*DPSRR        
c!      (TAIL MODES ONLY)

        D=D0+DELTADY*(Y/20.D0)**2   
c!  SHEET HALF-THICKNESS FOR THE TAIL MODES
        DDDY=DELTADY*Y*0.005D0      
c!  (THICKENS TO FLANKS, BUT NO VARIATION
C                                         ALONG X, IN CONTRAST TO RING CURRENT)
C
        DZETAS=DSQRT(ZS**2+D**2)  
c!  THIS IS THE SAME SIMPLE WAY TO SPREAD
C                                        OUT THE SHEET, AS THAT USED IN T89
        DDZETADX=ZS*DZSX/DZETAS
        DDZETADY=(ZS*DZSY+D*DDDY)/DZETAS
        DDZETADZ=ZS*DZSZ/DZETAS
C
        CALL SHLCAR3X3_96(ARC,X,Y,Z,SPS,WX,WY,WZ)
        CALL RINGCURR96(X,Y,Z,HX,HY,HZ)
        BXRC=WX+HX
        BYRC=WY+HY
        BZRC=WZ+HZ
C
        CALL SHLCAR3X3_96(ATAIL2,X,Y,Z,SPS,WX,WY,WZ)
        CALL TAILDISK_96(X,Y,Z,HX,HY,HZ)
        BXT2=WX+HX
        BYT2=WY+HY
        BZT2=WZ+HZ
C
        CALL SHLCAR3X3_96(ATAIL3,X,Y,Z,SPS,WX,WY,WZ)
        CALL TAIL87(X,Z,HX,HZ)
        BXT3=WX+HX
        BYT3=WY
        BZT3=WZ+HZ
C
        RETURN
        END
C
c********************************************************************
C
        SUBROUTINE RINGCURR96(X,Y,Z,BX,BY,BZ)
c
c       THIS SUBROUTINE COMPUTES THE COMPONENTS OF THE RING CURRENT FIELD,
C        SIMILAR TO THAT DESCRIBED BY TSYGANENKO AND PEREDO (1994).  THE
C          DIFFERENCE IS THAT NOW WE USE SPACEWARPING, AS DESCRIBED IN THE
C           PAPER ON MODELING BIRKELAND CURRENTS (TSYGANENKO AND STERN, 1996),
C            INSTEAD OF SHEARING IT IN THE SPIRIT OF THE T89 TAIL MODEL.
C
C          IN  ADDITION, INSTEAD OF 7 TERMS FOR THE RING CURRENT MODEL, WE USE
C             NOW ONLY 2 TERMS;  THIS SIMPLIFICATION ALSO GIVES RISE TO AN
C                EASTWARD RING CURRENT LOCATED EARTHWARD FROM THE MAIN ONE,
C                  IN LINE WITH WHAT IS ACTUALLY OBSERVED
C
C             FOR DETAILS, SEE NB #3, PAGES 70-73
C
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION F(2),BETA(2)
        COMMON /WARP/ CPSS,SPSS,DPSRR, XNEXT(3),XS,ZSWARPED,DXSX,DXSY,
     *   DXSZ,DZSX,DZSYWARPED,DZSZ,OTHER(4),ZS  
c!  ZS HERE IS WITHOUT Y-Z WARP
C

      DATA D0,DELTADX,XD,XLDX /2.d0,0.d0,0.,4.d0/  
c!  ACHTUNG !!  THE RC IS NOW
C                                            COMPLETELY SYMMETRIC (DELTADX=0)

C
        DATA F,BETA /569.895366D0,-1603.386993D0,2.722188D0,3.766875D0/
C
C  THE ORIGINAL VALUES OF F(I) WERE MULTIPLIED BY BETA(I) (TO REDUCE THE
C     NUMBER OF MULTIPLICATIONS BELOW)  AND BY THE FACTOR -0.43, NORMALIZING
C      THE DISTURBANCE AT ORIGIN  TO  B=-1nT
C
           DZSY=XS*Y*DPSRR  
c! NO WARPING IN THE Y-Z PLANE (ALONG X ONLY), AND
C                         THIS IS WHY WE DO NOT USE  DZSY FROM THE COMMON-BLOCK
           XXD=X-XD
           FDX=0.5D0*(1.D0+XXD/DSQRT(XXD**2+XLDX**2))
           DDDX=DELTADX*0.5D0*XLDX**2/DSQRT(XXD**2+XLDX**2)**3
           D=D0+DELTADX*FDX

           DZETAS=DSQRT(ZS**2+D**2)  
c!  THIS IS THE SAME SIMPLE WAY TO SPREAD
C                                        OUT THE SHEET, AS THAT USED IN T89
           RHOS=DSQRT(XS**2+Y**2)
           DDZETADX=(ZS*DZSX+D*DDDX)/DZETAS
           DDZETADY=ZS*DZSY/DZETAS
           DDZETADZ=ZS*DZSZ/DZETAS
         IF (RHOS.LT.1.D-5) THEN
            DRHOSDX=0.D0
            DRHOSDY=DSIGN(1.D0,Y)
            DRHOSDZ=0.D0
          ELSE
           DRHOSDX=XS*DXSX/RHOS
           DRHOSDY=(XS*DXSY+Y)/RHOS
           DRHOSDZ=XS*DXSZ/RHOS
         ENDIF
C
           BX=0.D0
           BY=0.D0
           BZ=0.D0
C
           DO 1 I=1,2
C
           BI=BETA(I)
C
           S1=DSQRT((DZETAS+BI)**2+(RHOS+BI)**2)
           S2=DSQRT((DZETAS+BI)**2+(RHOS-BI)**2)
           DS1DDZ=(DZETAS+BI)/S1
           DS2DDZ=(DZETAS+BI)/S2
           DS1DRHOS=(RHOS+BI)/S1
           DS2DRHOS=(RHOS-BI)/S2
C
           DS1DX=DS1DDZ*DDZETADX+DS1DRHOS*DRHOSDX
           DS1DY=DS1DDZ*DDZETADY+DS1DRHOS*DRHOSDY
           DS1DZ=DS1DDZ*DDZETADZ+DS1DRHOS*DRHOSDZ
C
           DS2DX=DS2DDZ*DDZETADX+DS2DRHOS*DRHOSDX
           DS2DY=DS2DDZ*DDZETADY+DS2DRHOS*DRHOSDY
           DS2DZ=DS2DDZ*DDZETADZ+DS2DRHOS*DRHOSDZ
C
           S1TS2=S1*S2
           S1PS2=S1+S2
           S1PS2SQ=S1PS2**2
           FAC1=DSQRT(S1PS2SQ-(2.D0*BI)**2)
           AS=FAC1/(S1TS2*S1PS2SQ)
           TERM1=1.D0/(S1TS2*S1PS2*FAC1)
           FAC2=AS/S1PS2SQ
           DASDS1=TERM1-FAC2/S1*(S2*S2+S1*(3.D0*S1+4.D0*S2))
           DASDS2=TERM1-FAC2/S2*(S1*S1+S2*(3.D0*S2+4.D0*S1))
C
           DASDX=DASDS1*DS1DX+DASDS2*DS2DX
           DASDY=DASDS1*DS1DY+DASDS2*DS2DY
           DASDZ=DASDS1*DS1DZ+DASDS2*DS2DZ
C
      BX=BX+F(I)*((2.D0*AS+Y*DASDY)*SPSS-XS*DASDZ
     *   +AS*DPSRR*(Y**2*CPSS+Z*ZS))
      BY=BY-F(I)*Y*(AS*DPSRR*XS+DASDZ*CPSS+DASDX*SPSS)
  1   BZ=BZ+F(I)*((2.D0*AS+Y*DASDY)*CPSS+XS*DASDX
     *   -AS*DPSRR*(X*ZS+Y**2*SPSS))
C
       RETURN
       END
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
         SUBROUTINE TAILDISK_96(X,Y,Z,BX,BY,BZ)
C
c
c       THIS SUBROUTINE COMPUTES THE COMPONENTS OF THE TAIL CURRENT FIELD,
C        SIMILAR TO THAT DESCRIBED BY TSYGANENKO AND PEREDO (1994).  THE
C          DIFFERENCE IS THAT NOW WE USE SPACEWARPING, AS DESCRIBED IN OUR
C           PAPER ON MODELING BIRKELAND CURRENTS (TSYGANENKO AND STERN, 1996)
C            INSTEAD OF SHEARING IT IN THE SPIRIT OF T89 TAIL MODEL.
C
C          IN  ADDITION, INSTEAD OF 8 TERMS FOR THE TAIL CURRENT MODEL, WE USE
C           NOW ONLY 4 TERMS
C
C             FOR DETAILS, SEE NB #3, PAGES 74-
C
         IMPLICIT REAL*8 (A-H,O-Z)
         DIMENSION F(4),BETA(4)
         COMMON /WARP/ CPSS,SPSS,DPSRR,XNEXT(3),XS,ZS,DXSX,DXSY,DXSZ,
     *    OTHER(3),DZETAS,DDZETADX,DDZETADY,DDZETADZ,ZSWW
C
         DATA XSHIFT /4.5d0/
C
         DATA F,BETA
     * / -745796.7338D0,1176470.141D0,-444610.529D0,-57508.01028D0,
     *   7.9250000D0,8.0850000D0,8.4712500D0,27.89500D0/
c
c  here original F(I) are multiplied by BETA(I), to economize
c    calculations
C
           RHOS=DSQRT((XS-XSHIFT)**2+Y**2)
         IF (RHOS.LT.1.D-5) THEN
            DRHOSDX=0.D0
            DRHOSDY=DSIGN(1.D0,Y)
            DRHOSDZ=0.D0
         ELSE
           DRHOSDX=(XS-XSHIFT)*DXSX/RHOS
           DRHOSDY=((XS-XSHIFT)*DXSY+Y)/RHOS
           DRHOSDZ=(XS-XSHIFT)*DXSZ/RHOS
         ENDIF
C
           BX=0.D0
           BY=0.D0
           BZ=0.D0
C
           DO 1 I=1,4
C
           BI=BETA(I)
C
           S1=DSQRT((DZETAS+BI)**2+(RHOS+BI)**2)
           S2=DSQRT((DZETAS+BI)**2+(RHOS-BI)**2)
           DS1DDZ=(DZETAS+BI)/S1
           DS2DDZ=(DZETAS+BI)/S2
           DS1DRHOS=(RHOS+BI)/S1
           DS2DRHOS=(RHOS-BI)/S2
C
           DS1DX=DS1DDZ*DDZETADX+DS1DRHOS*DRHOSDX
           DS1DY=DS1DDZ*DDZETADY+DS1DRHOS*DRHOSDY
           DS1DZ=DS1DDZ*DDZETADZ+DS1DRHOS*DRHOSDZ
C
           DS2DX=DS2DDZ*DDZETADX+DS2DRHOS*DRHOSDX
           DS2DY=DS2DDZ*DDZETADY+DS2DRHOS*DRHOSDY
           DS2DZ=DS2DDZ*DDZETADZ+DS2DRHOS*DRHOSDZ
C
           S1TS2=S1*S2
           S1PS2=S1+S2
           S1PS2SQ=S1PS2**2
           FAC1=DSQRT(S1PS2SQ-(2.D0*BI)**2)
           AS=FAC1/(S1TS2*S1PS2SQ)
           TERM1=1.D0/(S1TS2*S1PS2*FAC1)
           FAC2=AS/S1PS2SQ
           DASDS1=TERM1-FAC2/S1*(S2*S2+S1*(3.D0*S1+4.D0*S2))
           DASDS2=TERM1-FAC2/S2*(S1*S1+S2*(3.D0*S2+4.D0*S1))
C
           DASDX=DASDS1*DS1DX+DASDS2*DS2DX
           DASDY=DASDS1*DS1DY+DASDS2*DS2DY
           DASDZ=DASDS1*DS1DZ+DASDS2*DS2DZ
C
      BX=BX+F(I)*((2.D0*AS+Y*DASDY)*SPSS-(XS-XSHIFT)*DASDZ
     *   +AS*DPSRR*(Y**2*CPSS+Z*ZSWW))
C
      BY=BY-F(I)*Y*(AS*DPSRR*XS+DASDZ*CPSS+DASDX*SPSS)
  1   BZ=BZ+F(I)*((2.D0*AS+Y*DASDY)*CPSS+(XS-XSHIFT)*DASDX
     *   -AS*DPSRR*(X*ZSWW+Y**2*SPSS))

       RETURN
       END

C-------------------------------------------------------------------------
C
      SUBROUTINE TAIL87(X,Z,BX,BZ)

      IMPLICIT REAL*8 (A-H,O-Z)

      COMMON/WARP/ FIRST(3), RPS,WARP,D, OTHER(13)
C
C      'LONG' VERSION OF THE 1987 TAIL MAGNETIC FIELD MODEL
C              (N.A.TSYGANENKO, PLANET. SPACE SCI., V.35, P.1347, 1987)
C
C      D   IS THE Y-DEPENDENT SHEET HALF-THICKNESS (INCREASING TOWARDS FLANKS)
C      RPS  IS THE TILT-DEPENDENT SHIFT OF THE SHEET IN THE Z-DIRECTION,
C           CORRESPONDING TO THE ASYMPTOTIC HINGING DISTANCE, DEFINED IN THE
C           MAIN SUBROUTINE (TAILRC96) FROM THE PARAMETERS RH AND DR OF THE
C           T96-TYPE MODULE, AND
C      WARP  IS THE BENDING OF THE SHEET FLANKS IN THE Z-DIRECTION, DIRECTED
C           OPPOSITE TO RPS, AND INCREASING WITH DIPOLE TILT AND |Y|
C

        DATA DD/3.d0/
C
      DATA HPI,RT,XN,X1,X2,B0,B1,B2,XN21,XNR,ADLN
     * /1.5707963d0,40.d0,-10.d0,
     * -1.261d0,-0.663d0,0.391734d0,5.89715d0,24.6833d0,76.37d0,
     : -0.1071d0,0.13238005d0/
C                !!!   THESE ARE NEW VALUES OF  X1, X2, B0, B1, B2,
C                       CORRESPONDING TO TSCALE=1, INSTEAD OF TSCALE=0.6
C
C  THE ABOVE QUANTITIES WERE DEFINED AS FOLLOWS:------------------------
C       HPI=PI/2
C       RT=40.      !  Z-POSITION OF UPPER AND LOWER ADDITIONAL SHEETS
C       XN=-10.     !  INNER EDGE POSITION
C
C       TSCALE=1  !  SCALING FACTOR, DEFINING THE RATE OF INCREASE OF THE
C                       CURRENT DENSITY TAILWARDS
C
c  ATTENTION !  NOW I HAVE CHANGED TSCALE TO:  TSCALE=1.0, INSTEAD OF 0.6
c                  OF THE PREVIOUS VERSION
c
C       B0=0.391734
C       B1=5.89715 *TSCALE
C       B2=24.6833 *TSCALE**2
C
C    HERE ORIGINAL VALUES OF THE MODE AMPLITUDES (P.77, NB#3) WERE NORMALIZED
C      SO THAT ASYMPTOTIC  BX=1  AT X=-200RE
C
C      X1=(4.589  -5.85) *TSCALE -(TSCALE-1.)*XN ! NONLINEAR PARAMETERS OF THE
C                                                         CURRENT FUNCTION
C      X2=(5.187  -5.85) *TSCALE -(TSCALE-1.)*XN
c
c
C      XN21=(XN-X1)**2
C      XNR=1./(XN-X2)
C      ADLN=-DLOG(XNR**2*XN21)
C
C---------------------------------------------------------------
C
      ZS=Z -RPS +WARP
      ZP=Z-RT
      ZM=Z+RT
C
      XNX=XN-X
      XNX2=XNX**2
      XC1=X-X1
      XC2=X-X2
      XC22=XC2**2
      XR2=XC2*XNR
      XC12=XC1**2
      D2=DD**2    
c!  SQUARE OF THE TOTAL HALFTHICKNESS (DD=3Re for this mode)
      B20=ZS**2+D2
      B2P=ZP**2+D2
      B2M=ZM**2+D2
      B=DSQRT(B20)
      BP=DSQRT(B2P)
      BM=DSQRT(B2M)
      XA1=XC12+B20
      XAP1=XC12+B2P
      XAM1=XC12+B2M
      XA2=1.d0/(XC22+B20)
      XAP2=1.d0/(XC22+B2P)
      XAM2=1.d0/(XC22+B2M)
      XNA=XNX2+B20
      XNAP=XNX2+B2P
      XNAM=XNX2+B2M
      F=B20-XC22
      FP=B2P-XC22
      FM=B2M-XC22
      XLN1=DLOG(XN21/XNA)
      XLNP1=DLOG(XN21/XNAP)
      XLNM1=DLOG(XN21/XNAM)
      XLN2=XLN1+ADLN
      XLNP2=XLNP1+ADLN
      XLNM2=XLNM1+ADLN
      ALN=0.25d0*(XLNP1+XLNM1-2.d0*XLN1)
      S0=(DATAN(XNX/B)+HPI)/B
      S0P=(DATAN(XNX/BP)+HPI)/BP
      S0M=(DATAN(XNX/BM)+HPI)/BM
      S1=(XLN1*.5d0+XC1*S0)/XA1
      S1P=(XLNP1*.5d0+XC1*S0P)/XAP1
      S1M=(XLNM1*.5d0+XC1*S0M)/XAM1
      S2=(XC2*XA2*XLN2-XNR-F*XA2*S0)*XA2
      S2P=(XC2*XAP2*XLNP2-XNR-FP*XAP2*S0P)*XAP2
      S2M=(XC2*XAM2*XLNM2-XNR-FM*XAM2*S0M)*XAM2
      G1=(B20*S0-0.5d0*XC1*XLN1)/XA1
      G1P=(B2P*S0P-0.5d0*XC1*XLNP1)/XAP1
      G1M=(B2M*S0M-0.5d0*XC1*XLNM1)/XAM1
      G2=((0.5d0*F*XLN2+2.d0*S0*B20*XC2)*XA2+XR2)*XA2
      G2P=((0.5d0*FP*XLNP2+2.d0*S0P*B2P*XC2)*XAP2+XR2)*XAP2
      G2M=((0.5d0*FM*XLNM2+2.d0*S0M*B2M*XC2)*XAM2+XR2)*XAM2
      BX=B0*(ZS*S0-0.5d0*(ZP*S0P+ZM*S0M))
     : +B1*(ZS*S1-0.5d0*(ZP*S1P+ZM*S1M))
     : +B2*(ZS*S2-0.5d0*(ZP*S2P+ZM*S2M))
      BZ=B0*ALN+B1*(G1-0.5d0*(G1P+G1M))+B2*(G2-0.5d0*(G2P+G2M))
C
C    CALCULATION OF THE MAGNETOTAIL CURRENT CONTRIBUTION IS FINISHED
C
      RETURN
      END

C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C THIS CODE RETURNS THE SHIELDING FIELD REPRESENTED BY  2x3x3=18 "CARTESIAN"
C    HARMONICS
C
         SUBROUTINE  SHLCAR3X3_96(A,X,Y,Z,SPS,HX,HY,HZ)
C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  The 36 coefficients enter in pairs in the amplitudes of the "cartesian"
c    harmonics (A(1)-A(36).
c  The 12 nonlinear parameters (A(37)-A(48) are the scales Pi,Ri,Qi,and Si
C   entering the arguments of exponents, sines, and cosines in each of the
C   18 "Cartesian" harmonics
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
       IMPLICIT  REAL * 8  (A - H, O - Z)
C
         DIMENSION A(48)
C
          CPS=DSQRT(1.D0-SPS**2)
          S3PS=4.D0*CPS**2-1.D0   
c!  THIS IS SIN(3*PS)/SIN(PS)
C
           HX=0.D0
           HY=0.D0
           HZ=0.D0
           L=0
C
           DO 1 M=1,2     
c!    M=1 IS FOR THE 1ST SUM ("PERP." SYMMETRY)
C                           AND M=2 IS FOR THE SECOND SUM ("PARALL." SYMMETRY)
             DO 2 I=1,3
                  P=A(36+I)
                  Q=A(42+I)
                  CYPI=DCOS(Y/P)
                  CYQI=DCOS(Y/Q)
                  SYPI=DSIN(Y/P)
                  SYQI=DSIN(Y/Q)
C
              DO 3 K=1,3
                   R=A(39+K)
                   S=A(45+K)
                   SZRK=DSIN(Z/R)
                   CZSK=DCOS(Z/S)
                   CZRK=DCOS(Z/R)
                   SZSK=DSIN(Z/S)
                     SQPR=DSQRT(1.D0/P**2+1.D0/R**2)
                     SQQS=DSQRT(1.D0/Q**2+1.D0/S**2)
                        EPR=DEXP(X*SQPR)
                        EQS=DEXP(X*SQQS)
C
                   DO 4 N=1,2  
c! N=1 IS FOR THE FIRST PART OF EACH COEFFICIENT
C                                  AND N=2 IS FOR THE SECOND ONE
C
                    L=L+1
                     IF (M.EQ.1) THEN
                       IF (N.EQ.1) THEN
                         DX=-SQPR*EPR*CYPI*SZRK
                         DY=EPR/P*SYPI*SZRK
                         DZ=-EPR/R*CYPI*CZRK
                         HX=HX+A(L)*DX
                         HY=HY+A(L)*DY
                         HZ=HZ+A(L)*DZ
                                   ELSE
                         DX=DX*CPS
                         DY=DY*CPS
                         DZ=DZ*CPS
                         HX=HX+A(L)*DX
                         HY=HY+A(L)*DY
                         HZ=HZ+A(L)*DZ
                                   ENDIF
                     ELSE
                       IF (N.EQ.1) THEN
                         DX=-SPS*SQQS*EQS*CYQI*CZSK
                         DY=SPS*EQS/Q*SYQI*CZSK
                         DZ=SPS*EQS/S*CYQI*SZSK
                         HX=HX+A(L)*DX
                         HY=HY+A(L)*DY
                         HZ=HZ+A(L)*DZ
                                   ELSE
                         DX=DX*S3PS
                         DY=DY*S3PS
                         DZ=DZ*S3PS
                         HX=HX+A(L)*DX
                         HY=HY+A(L)*DY
                         HZ=HZ+A(L)*DZ
                       ENDIF
                 ENDIF
c
  4   CONTINUE
  3   CONTINUE
  2   CONTINUE
  1   CONTINUE
C
         RETURN
       END

C
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C
       SUBROUTINE BIRK1TOT_02(PS,X,Y,Z,BX,BY,BZ)
C
C  THIS IS THE SECOND VERSION OF THE ANALYTICAL MODEL OF THE REGION 1 FIELD
C   BASED ON A SEPARATE REPRESENTATION OF THE POTENTIAL FIELD IN THE INNER AND
C   OUTER SPACE, MAPPED BY MEANS OF A SPHERO-DIPOLAR COORDINATE SYSTEM (NB #3,
C   P.91).   THE DIFFERENCE FROM THE FIRST ONE IS THAT INSTEAD OF OCTAGONAL
C   CURRENT LOOPS, CIRCULAR ONES ARE USED IN THIS VERSION FOR APPROXIMATING THE
C   FIELD IN THE OUTER REGION, WHICH IS FASTER.
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION D1(3,26),D2(3,79),XI(4),C1(26),C2(79)

         COMMON /COORD11/ XX1(12),YY1(12)
         COMMON /RHDR/ RH,DR
         COMMON /LOOPDIP1/ TILT,XCENTRE(2),RADIUS(2), DIPX,DIPY
C
         COMMON /COORD21/ XX2(14),YY2(14),ZZ2(14)
         COMMON /DX1/ DX,SCALEIN,SCALEOUT
C
      DATA C1/-0.911582d-03,-0.376654d-02,-0.727423d-02,-0.270084d-02,
     * -0.123899d-02,-0.154387d-02,-0.340040d-02,-0.191858d-01,
     * -0.518979d-01,0.635061d-01,0.440680d0,-0.396570d0,0.561238d-02,
     *  0.160938d-02,-0.451229d-02,-0.251810d-02,-0.151599d-02,
     * -0.133665d-02,-0.962089d-03,-0.272085d-01,-0.524319d-01,
     *  0.717024d-01,0.523439d0,-0.405015d0,-89.5587d0,23.2806d0/

C
      DATA C2/6.04133d0,.305415d0,.606066d-02,.128379d-03,-.179406d-04,
     * 1.41714d0,-27.2586d0,-4.28833d0,-1.30675d0,35.5607d0,8.95792d0,
     : .961617d-03,
     * -.801477d-03,-.782795d-03,-1.65242d0,-16.5242d0,-5.33798d0,
     : .424878d-03,
     * .331787d-03,-.704305d-03,.844342d-03,.953682d-04,.886271d-03,
     * 25.1120d0,20.9299d0,5.14569d0,-44.1670d0,-51.0672d0,-1.87725d0,
     : 20.2998d0,
     * 48.7505d0,-2.97415d0,3.35184d0,-54.2921d0,-.838712d0,-10.5123d0,
     : 70.7594d0,
     * -4.94104d0,.106166d-03,.465791d-03,-.193719d-03,10.8439d0,
     : -29.7968d0,
     *  8.08068d0,.463507d-03,-.224475d-04,.177035d-03,-.317581d-03,
     * -.264487d-03,.102075d-03,7.71390d0,10.1915d0,-4.99797d0,
     : -23.1114d0,
     * -29.2043d0,12.2928d0,10.9542d0,33.6671d0,-9.3851d0,.174615d-03,
     : -.789777d-06,
     * .686047d-03,.460104d-04,-.345216d-02,.221871d-02,.110078d-01,
     * -.661373d-02,.249201d-02,.343978d-01,-.193145d-05,.493963d-05,
     * -.535748d-04,.191833d-04,-.100496d-03,-.210103d-03,-.232195d-02,
     * .315335d-02,-.134320d-01,-.263222d-01/
c
      DATA TILT,XCENTRE,RADIUS,DIPX,DIPY /1.00891d0,2.28397d0,
     : -5.60831d0,
     * 1.86106d0,7.83281d0,1.12541d0,0.945719d0/

      DATA DX,SCALEIN,SCALEOUT /-0.16D0,0.08D0,0.4D0/
      DATA XX1/-11.D0,2*-7.D0,2*-3.D0,3*1.D0,2*5.D0,2*9.D0/
      DATA YY1/2.D0,0.D0,4.D0,2.D0,6.D0,0.D0,4.D0,8.D0,2.D0,6.D0,0.D0,
     *  4.D0/
      DATA XX2/-10.D0,-7.D0,2*-4.D0,0.D0,2*4.D0,7.D0,10.D0,5*0.D0/
      DATA YY2/3.D0,6.D0,3.D0,9.D0,6.D0,3.D0,9.D0,6.D0,3.D0,5*0.D0/
      DATA ZZ2/2*20.D0,4.D0,20.D0,2*4.D0,3*20.D0,2.D0,3.D0,4.5D0,
     *  7.D0,10.D0/
C
      DATA RH,DR /9.D0,4.D0/   
c!  RH IS THE "HINGING DISTANCE" AND DR IS THE
C                                TRANSITION SCALE LENGTH, DEFINING THE
C                                CURVATURE  OF THE WARPING (SEE P.89, NB #2)
C
      DATA XLTDAY,XLTNGHT /78.D0,70.D0/  
c!  THESE ARE LATITUDES OF THE R-1 OVAL
C                                             AT NOON AND AT MIDNIGHT
      DATA DTET0 /0.034906d0/   
c!   THIS IS THE LATITUDINAL HALF-THICKNESS OF THE
C                                  R-1 OVAL (THE INTERPOLATION REGION BETWEEN
C                                    THE HIGH-LAT. AND THE PLASMA SHEET)
C
        TNOONN=(90.D0-XLTDAY)*0.01745329D0
        TNOONS=3.141592654D0-TNOONN     
c! HERE WE ASSUME THAT THE POSITIONS OF
C                                          THE NORTHERN AND SOUTHERN R-1 OVALS
C                                          ARE SYMMETRIC IN THE SM-COORDINATES
        DTETDN=(XLTDAY-XLTNGHT)*0.01745329D0
        DR2=DR**2
C
      SPS=DSIN(PS)
      R2=X**2+Y**2+Z**2
      R=DSQRT(R2)
      R3=R*R2
C
      RMRH=R-RH
      RPRH=R+RH
      SQM=DSQRT(RMRH**2+DR2)
      SQP=DSQRT(RPRH**2+DR2)
      C=SQP-SQM
      Q=DSQRT((RH+1.D0)**2+DR2)-DSQRT((RH-1.D0)**2+DR2)
      SPSAS=SPS/R*C/Q
      CPSAS=DSQRT(1.D0-SPSAS**2)
       XAS = X*CPSAS-Z*SPSAS
       ZAS = X*SPSAS+Z*CPSAS
        IF (XAS.NE.0.D0.OR.Y.NE.0.D0) THEN
          PAS = DATAN2(Y,XAS)
                                      ELSE
          PAS=0.D0
        ENDIF
C
      TAS=DATAN2(DSQRT(XAS**2+Y**2),ZAS)
      STAS=DSIN(TAS)
      F=STAS/(STAS**6*(1.D0-R3)+R3)**0.1666666667D0
C
      TET0=DASIN(F)
      IF (TAS.GT.1.5707963D0) TET0=3.141592654D0-TET0
      DTET=DTETDN*DSIN(PAS*0.5D0)**2
      TETR1N=TNOONN+DTET
      TETR1S=TNOONS-DTET
C
C NOW LET'S DEFINE WHICH OF THE FOUR REGIONS (HIGH-LAT., NORTHERN PSBL,
C   PLASMA SHEET, SOUTHERN PSBL) DOES THE POINT (X,Y,Z) BELONG TO:
C
       IF (TET0.LT.TETR1N-DTET0.OR.TET0.GT.TETR1S+DTET0)  LOC=1 
c! HIGH-LAT.
       IF (TET0.GT.TETR1N+DTET0.AND.TET0.LT.TETR1S-DTET0) LOC=2 
c! PL.SHEET
       IF (TET0.GE.TETR1N-DTET0.AND.TET0.LE.TETR1N+DTET0) LOC=3 
c! NORTH PSBL
       IF (TET0.GE.TETR1S-DTET0.AND.TET0.LE.TETR1S+DTET0) LOC=4 
c! SOUTH PSBL
C
       IF (LOC.EQ.1) THEN   
c! IN THE HIGH-LAT. REGION USE THE SUBROUTINE DIPOCT
C
C      print *, '  LOC=1 (HIGH-LAT)'    
c!  (test printout; disabled now)
         XI(1)=X
         XI(2)=Y
         XI(3)=Z
         XI(4)=PS
         CALL  DIPLOOP1(XI,D1)
          BX=0.D0
          BY=0.D0
          BZ=0.D0
            DO 1 I=1,26
              BX=BX+C1(I)*D1(1,I)
              BY=BY+C1(I)*D1(2,I)
  1           BZ=BZ+C1(I)*D1(3,I)
       ENDIF                                           
c!  END OF THE CASE 1
C
       IF (LOC.EQ.2) THEN
C           print *, '  LOC=2 (PLASMA SHEET)'  
c!  (test printout; disabled now)
C
         XI(1)=X
         XI(2)=Y
         XI(3)=Z
         XI(4)=PS
         CALL  CONDIP1(XI,D2)
          BX=0.D0
          BY=0.D0
          BZ=0.D0
            DO 2 I=1,79
              BX=BX+C2(I)*D2(1,I)
              BY=BY+C2(I)*D2(2,I)
  2           BZ=BZ+C2(I)*D2(3,I)
       ENDIF                                           
c!   END OF THE CASE 2
C
       IF (LOC.EQ.3) THEN
C       print *, '  LOC=3 (north PSBL)'  
c!  (test printout; disabled now)
C
         T01=TETR1N-DTET0
         T02=TETR1N+DTET0
          SQR=DSQRT(R)
          ST01AS=SQR/(R3+1.D0/DSIN(T01)**6-1.D0)**0.1666666667
          ST02AS=SQR/(R3+1.D0/DSIN(T02)**6-1.D0)**0.1666666667
          CT01AS=DSQRT(1.D0-ST01AS**2)
          CT02AS=DSQRT(1.D0-ST02AS**2)
         XAS1=R*ST01AS*DCOS(PAS)
         Y1=  R*ST01AS*DSIN(PAS)
         ZAS1=R*CT01AS
         X1=XAS1*CPSAS+ZAS1*SPSAS
         Z1=-XAS1*SPSAS+ZAS1*CPSAS 
c! X1,Y1,Z1 ARE COORDS OF THE NORTHERN
c                                                      BOUNDARY POINT
         XI(1)=X1
         XI(2)=Y1
         XI(3)=Z1
         XI(4)=PS
         CALL  DIPLOOP1(XI,D1)
          BX1=0.D0
          BY1=0.D0
          BZ1=0.D0
            DO 11 I=1,26
              BX1=BX1+C1(I)*D1(1,I) 
c!   BX1,BY1,BZ1  ARE FIELD COMPONENTS
              BY1=BY1+C1(I)*D1(2,I)  
c!  IN THE NORTHERN BOUNDARY POINT
 11           BZ1=BZ1+C1(I)*D1(3,I)  
c!
C
         XAS2=R*ST02AS*DCOS(PAS)
         Y2=  R*ST02AS*DSIN(PAS)
         ZAS2=R*CT02AS
         X2=XAS2*CPSAS+ZAS2*SPSAS
         Z2=-XAS2*SPSAS+ZAS2*CPSAS 
c! X2,Y2,Z2 ARE COORDS OF THE SOUTHERN
C                                        BOUNDARY POINT
         XI(1)=X2
         XI(2)=Y2
         XI(3)=Z2
         XI(4)=PS
         CALL  CONDIP1(XI,D2)
          BX2=0.D0
          BY2=0.D0
          BZ2=0.D0
            DO 12 I=1,79
              BX2=BX2+C2(I)*D2(1,I)
c!  BX2,BY2,BZ2  ARE FIELD COMPONENTS
              BY2=BY2+C2(I)*D2(2,I) 
c!  IN THE SOUTHERN BOUNDARY POINT
  12          BZ2=BZ2+C2(I)*D2(3,I)
C
C  NOW INTERPOLATE:
C
         SS=DSQRT((X2-X1)**2+(Y2-Y1)**2+(Z2-Z1)**2)
         DS=DSQRT((X-X1)**2+(Y-Y1)**2+(Z-Z1)**2)
         FRAC=DS/SS
         BX=BX1*(1.D0-FRAC)+BX2*FRAC
         BY=BY1*(1.D0-FRAC)+BY2*FRAC
         BZ=BZ1*(1.D0-FRAC)+BZ2*FRAC
C
        ENDIF                                              
c! END OF THE CASE 3
C
        IF (LOC.EQ.4) THEN
C       print *, '  LOC=4 (south PSBL)'  
c!  (test printout; disabled now)
C
         T01=TETR1S-DTET0
         T02=TETR1S+DTET0
          SQR=DSQRT(R)
          ST01AS=SQR/(R3+1.D0/DSIN(T01)**6-1.D0)**0.1666666667
          ST02AS=SQR/(R3+1.D0/DSIN(T02)**6-1.D0)**0.1666666667
          CT01AS=-DSQRT(1.D0-ST01AS**2)
          CT02AS=-DSQRT(1.D0-ST02AS**2)
         XAS1=R*ST01AS*DCOS(PAS)
         Y1=  R*ST01AS*DSIN(PAS)
         ZAS1=R*CT01AS
         X1=XAS1*CPSAS+ZAS1*SPSAS
         Z1=-XAS1*SPSAS+ZAS1*CPSAS 
c! X1,Y1,Z1 ARE COORDS OF THE NORTHERN
C                                               BOUNDARY POINT
         XI(1)=X1
         XI(2)=Y1
         XI(3)=Z1
         XI(4)=PS
         CALL  CONDIP1(XI,D2)
          BX1=0.D0
          BY1=0.D0
          BZ1=0.D0
            DO 21 I=1,79
              BX1=BX1+C2(I)*D2(1,I) 
c!  BX1,BY1,BZ1  ARE FIELD COMPONENTS
              BY1=BY1+C2(I)*D2(2,I)  
c!  IN THE NORTHERN BOUNDARY POINT
 21           BZ1=BZ1+C2(I)*D2(3,I)  
c!
C
         XAS2=R*ST02AS*DCOS(PAS)
         Y2=  R*ST02AS*DSIN(PAS)
         ZAS2=R*CT02AS
         X2=XAS2*CPSAS+ZAS2*SPSAS
         Z2=-XAS2*SPSAS+ZAS2*CPSAS 
c! X2,Y2,Z2 ARE COORDS OF THE SOUTHERN
C                                          BOUNDARY POINT
         XI(1)=X2
         XI(2)=Y2
         XI(3)=Z2
         XI(4)=PS
         CALL  DIPLOOP1(XI,D1)
          BX2=0.D0
          BY2=0.D0
          BZ2=0.D0
            DO 22 I=1,26
              BX2=BX2+C1(I)*D1(1,I) 
c!  BX2,BY2,BZ2  ARE FIELD COMPONENTS
              BY2=BY2+C1(I)*D1(2,I) 
c!     IN THE SOUTHERN BOUNDARY POINT
  22          BZ2=BZ2+C1(I)*D1(3,I)
C
C  NOW INTERPOLATE:
C
         SS=DSQRT((X2-X1)**2+(Y2-Y1)**2+(Z2-Z1)**2)
         DS=DSQRT((X-X1)**2+(Y-Y1)**2+(Z-Z1)**2)
         FRAC=DS/SS
         BX=BX1*(1.D0-FRAC)+BX2*FRAC
         BY=BY1*(1.D0-FRAC)+BY2*FRAC
         BZ=BZ1*(1.D0-FRAC)+BZ2*FRAC
C
        ENDIF                                        
c! END OF THE CASE 4
C
C   NOW, LET US ADD THE SHIELDING FIELD
C
        CALL  BIRK1SHLD(PS,X,Y,Z,BSX,BSY,BSZ)
        BX=BX+BSX
        BY=BY+BSY
        BZ=BZ+BSZ
         RETURN
         END
C
C------------------------------------------------------------------------------
C
C
         SUBROUTINE  DIPLOOP1(XI,D)
C
C
C      Calculates dependent model variables and their deriva-
C  tives for given independent variables and model parame-
C  ters.  Specifies model functions with free parameters which
C  must be determined by means of least squares fits (RMS
C  minimization procedure).
C
C      Description of parameters:
C
C  XI  - input vector containing independent variables;
C  D   - output double precision vector containing
C        calculated values for derivatives of dependent
C        variables with respect to LINEAR model parameters;
C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
c  The  26 coefficients are moments (Z- and X-components) of 12 dipoles placed
C    inside the  R1-shell,  PLUS amplitudes of two octagonal double loops.
C     The dipoles with nonzero  Yi appear in pairs with equal moments.
c                  (see the notebook #2, pp.102-103, for details)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
         IMPLICIT  REAL * 8  (A - H, O - Z)
C
         COMMON /COORD11/ XX(12),YY(12)
         COMMON /LOOPDIP1/ TILT,XCENTRE(2),RADIUS(2),  DIPX,DIPY
         COMMON /RHDR/RH,DR
         DIMENSION XI(4),D(3,26)
C
           X = XI(1)
         Y = XI(2)
         Z = XI(3)
           PS= XI(4)
           SPS=DSIN(PS)
C
         DO 1 I=1,12
           R2=(XX(I)*DIPX)**2+(YY(I)*DIPY)**2
           R=DSQRT(R2)
             RMRH=R-RH
             RPRH=R+RH
             DR2=DR**2
             SQM=DSQRT(RMRH**2+DR2)
             SQP=DSQRT(RPRH**2+DR2)
             C=SQP-SQM
             Q=DSQRT((RH+1.D0)**2+DR2)-DSQRT((RH-1.D0)**2+DR2)
             SPSAS=SPS/R*C/Q
             CPSAS=DSQRT(1.D0-SPSAS**2)
         XD= (XX(I)*DIPX)*CPSAS
         YD= (YY(I)*DIPY)
         ZD=-(XX(I)*DIPX)*SPSAS
      CALL DIPXYZ(X-XD,Y-YD,Z-ZD,BX1X,BY1X,BZ1X,BX1Y,BY1Y,BZ1Y,
     *  BX1Z,BY1Z,BZ1Z)
        IF (DABS(YD).GT.1.D-10) THEN
      CALL DIPXYZ(X-XD,Y+YD,Z-ZD,BX2X,BY2X,BZ2X,BX2Y,BY2Y,BZ2Y,
     *  BX2Z,BY2Z,BZ2Z)
                                   ELSE
        BX2X=0.D0
        BY2X=0.D0
        BZ2X=0.D0
C
        BX2Z=0.D0
        BY2Z=0.D0
        BZ2Z=0.D0
                                   ENDIF
C
            D(1,I)=BX1Z+BX2Z
            D(2,I)=BY1Z+BY2Z
            D(3,I)=BZ1Z+BZ2Z
            D(1,I+12)=(BX1X+BX2X)*SPS
            D(2,I+12)=(BY1X+BY2X)*SPS
            D(3,I+12)=(BZ1X+BZ2X)*SPS
  1   CONTINUE
c
           R2=(XCENTRE(1)+RADIUS(1))**2
           R=DSQRT(R2)
             RMRH=R-RH
             RPRH=R+RH
             DR2=DR**2
             SQM=DSQRT(RMRH**2+DR2)
             SQP=DSQRT(RPRH**2+DR2)
             C=SQP-SQM
             Q=DSQRT((RH+1.D0)**2+DR2)-DSQRT((RH-1.D0)**2+DR2)
             SPSAS=SPS/R*C/Q
             CPSAS=DSQRT(1.D0-SPSAS**2)
         XOCT1= X*CPSAS-Z*SPSAS
         YOCT1= Y
         ZOCT1= X*SPSAS+Z*CPSAS
C
      CALL CROSSLP(XOCT1,YOCT1,ZOCT1,BXOCT1,BYOCT1,BZOCT1,XCENTRE(1),
     *        RADIUS(1),TILT)
            D(1,25)=BXOCT1*CPSAS+BZOCT1*SPSAS
            D(2,25)=BYOCT1
            D(3,25)=-BXOCT1*SPSAS+BZOCT1*CPSAS
C
           R2=(RADIUS(2)-XCENTRE(2))**2
           R=DSQRT(R2)
             RMRH=R-RH
             RPRH=R+RH
             DR2=DR**2
             SQM=DSQRT(RMRH**2+DR2)
             SQP=DSQRT(RPRH**2+DR2)
             C=SQP-SQM
             Q=DSQRT((RH+1.D0)**2+DR2)-DSQRT((RH-1.D0)**2+DR2)
             SPSAS=SPS/R*C/Q
             CPSAS=DSQRT(1.D0-SPSAS**2)
         XOCT2= X*CPSAS-Z*SPSAS -XCENTRE(2)
         YOCT2= Y
         ZOCT2= X*SPSAS+Z*CPSAS
            CALL CIRCLE(XOCT2,YOCT2,ZOCT2,RADIUS(2),BX,BY,BZ)
            D(1,26) =  BX*CPSAS+BZ*SPSAS
            D(2,26) =  BY
            D(3,26) = -BX*SPSAS+BZ*CPSAS
C
            RETURN
            END
c-------------------------------------------------------------------------
C
        SUBROUTINE CIRCLE(X,Y,Z,RL,BX,BY,BZ)
C
C  RETURNS COMPONENTS OF THE FIELD FROM A CIRCULAR CURRENT LOOP OF RADIUS RL
C  USES THE SECOND (MORE ACCURATE) APPROXIMATION GIVEN IN ABRAMOWITZ AND STEGUN

        IMPLICIT REAL*8 (A-H,O-Z)
        REAL*8 K
cccccccc        DATA PI/3.141592654D0/
        pi=atan(1.0d0)*4
C
        RHO2=X*X+Y*Y
        RHO=DSQRT(RHO2)
        R22=Z*Z+(RHO+RL)**2
        R2=DSQRT(R22)
        R12=R22-4.D0*RHO*RL
        R32=0.5D0*(R12+R22)
        XK2=1.D0-R12/R22
        XK2S=1.D0-XK2
        DL=DLOG(1.D0/XK2S)
        K=1.38629436112d0+XK2S*(0.09666344259D0+XK2S*(0.03590092383d0+
     *     XK2S*(0.03742563713d0+XK2S*0.01451196212d0))) +DL*
     *     (0.5D0+XK2S*(0.12498593597D0+XK2S*(0.06880248576D0+
     *     XK2S*(0.03328355346D0+XK2S*0.00441787012D0))))
        E=1.D0+XK2S*(0.44325141463D0+XK2S*(0.0626060122D0+XK2S*
     *      (0.04757383546D0+XK2S*0.01736506451D0))) +DL*
     *     XK2S*(0.2499836831D0+XK2S*(0.09200180037D0+XK2S*
     *       (0.04069697526D0+XK2S*0.00526449639D0)))

        IF (RHO.GT.1.D-6) THEN
           BRHO=Z/(RHO2*R2)*(R32/R12*E-K) 
c!  THIS IS NOT EXACTLY THE B-RHO COM-
                           ELSE           
c!   PONENT - NOTE THE ADDITIONAL
           BRHO=PI*RL/R2*(RL-RHO)/R12*Z/(R32-RHO2)  
c!      DIVISION BY RHO
        ENDIF

        BX=BRHO*X
        BY=BRHO*Y
        BZ=(K-E*(R32-2.D0*RL*RL)/R12)/R2
        RETURN
        END
C-------------------------------------------------------------
C
        SUBROUTINE CROSSLP(X,Y,Z,BX,BY,BZ,XC,RL,AL)
C
c   RETURNS FIELD COMPONENTS OF A PAIR OF LOOPS WITH A COMMON CENTER AND
C    DIAMETER,  COINCIDING WITH THE X AXIS. THE LOOPS ARE INCLINED TO THE
C    EQUATORIAL PLANE BY THE ANGLE AL (RADIANS) AND SHIFTED IN THE POSITIVE
C     X-DIRECTION BY THE DISTANCE  XC.
c
        IMPLICIT REAL*8 (A-H,O-Z)
C
            CAL=DCOS(AL)
            SAL=DSIN(AL)
C
        Y1=Y*CAL-Z*SAL
        Z1=Y*SAL+Z*CAL
        Y2=Y*CAL+Z*SAL
        Z2=-Y*SAL+Z*CAL
        CALL CIRCLE(X-XC,Y1,Z1,RL,BX1,BY1,BZ1)
        CALL CIRCLE(X-XC,Y2,Z2,RL,BX2,BY2,BZ2)
        BX=BX1+BX2
        BY= (BY1+BY2)*CAL+(BZ1-BZ2)*SAL
        BZ=-(BY1-BY2)*SAL+(BZ1+BZ2)*CAL
C
        RETURN
        END

C*******************************************************************

       SUBROUTINE DIPXYZ(X,Y,Z,BXX,BYX,BZX,BXY,BYY,BZY,BXZ,BYZ,BZZ)
C
C       RETURNS THE FIELD COMPONENTS PRODUCED BY THREE DIPOLES, EACH
C        HAVING M=Me AND ORIENTED PARALLEL TO X,Y, and Z AXIS, RESP.
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      X2=X**2
      Y2=Y**2
      Z2=Z**2
      R2=X2+Y2+Z2

      XMR5=30574.D0/(R2*R2*DSQRT(R2))
      XMR53=3.D0*XMR5
      BXX=XMR5*(3.D0*X2-R2)
      BYX=XMR53*X*Y
      BZX=XMR53*X*Z
C
      BXY=BYX
      BYY=XMR5*(3.D0*Y2-R2)
      BZY=XMR53*Y*Z
C
      BXZ=BZX
      BYZ=BZY
      BZZ=XMR5*(3.D0*Z2-R2)
C
      RETURN
      END
C
C------------------------------------------------------------------------------
         SUBROUTINE  CONDIP1(XI,D)
C
C      Calculates dependent model variables and their derivatives for given
C  independent variables and model parameters.  Specifies model functions with
C  free parameters which must be determined by means of least squares fits
C  (RMS minimization procedure).
C
C      Description of parameters:
C
C  XI  - input vector containing independent variables;
C  D   - output double precision vector containing
C        calculated values for derivatives of dependent
C        variables with respect to LINEAR model parameters;
C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
c  The  79 coefficients are (1) 5 amplitudes of the conical harmonics, plus
c                           (2) (9x3+5x2)x2=74 components of the dipole moments
c              (see the notebook #2, pp.113-..., for details)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
         IMPLICIT  REAL * 8  (A - H, O - Z)
C
      COMMON/DX1/ DX,SCALEIN,SCALEOUT
      COMMON/COORD21/ XX(14),YY(14),ZZ(14)
c
         DIMENSION XI(4),D(3,79),CF(5),SF(5)
C
           X = XI(1)
         Y = XI(2)
         Z = XI(3)
           PS= XI(4)
           SPS=DSIN(PS)
           CPS=DCOS(PS)
C
      XSM=X*CPS-Z*SPS  - DX
      ZSM=Z*CPS+X*SPS
      RO2=XSM**2+Y**2
      RO=SQRT(RO2)
C
      CF(1)=XSM/RO
      SF(1)=Y/RO
C
      CF(2)=CF(1)**2-SF(1)**2
      SF(2)=2.d0*SF(1)*CF(1)
      CF(3)=CF(2)*CF(1)-SF(2)*SF(1)
      SF(3)=SF(2)*CF(1)+CF(2)*SF(1)
      CF(4)=CF(3)*CF(1)-SF(3)*SF(1)
      SF(4)=SF(3)*CF(1)+CF(3)*SF(1)
      CF(5)=CF(4)*CF(1)-SF(4)*SF(1)
      SF(5)=SF(4)*CF(1)+CF(4)*SF(1)
C
      R2=RO2+ZSM**2
      R=DSQRT(R2)
      C=ZSM/R
      S=RO/R
      CH=DSQRT(0.5D0*(1.D0+C))
      SH=DSQRT(0.5D0*(1.D0-C))
      TNH=SH/CH
      CNH=1.D0/TNH
C
      DO 1 M=1,5
       BT=M*CF(M)/(R*S)*(TNH**M+CNH**M)
       BF=-0.5D0*M*SF(M)/R*(TNH**(M-1)/CH**2-CNH**(M-1)/SH**2)
       BXSM=BT*C*CF(1)-BF*SF(1)
       BY=BT*C*SF(1)+BF*CF(1)
       BZSM=-BT*S
C
       D(1,M)=BXSM*CPS+BZSM*SPS
       D(2,M)=BY
  1    D(3,M)=-BXSM*SPS+BZSM*CPS
C
      XSM = X*CPS-Z*SPS
      ZSM = Z*CPS+X*SPS
C
        DO 2 I=1,9
C
        IF (I.EQ.3.OR.I.EQ.5.OR.I.EQ.6) THEN
                XD =  XX(I)*SCALEIN
                YD =  YY(I)*SCALEIN
                                         ELSE
                XD =  XX(I)*SCALEOUT
                YD =  YY(I)*SCALEOUT
        ENDIF
C
         ZD =  ZZ(I)
C
      CALL DIPXYZ(XSM-XD,Y-YD,ZSM-ZD,BX1X,BY1X,BZ1X,BX1Y,BY1Y,BZ1Y,
     *  BX1Z,BY1Z,BZ1Z)
      CALL DIPXYZ(XSM-XD,Y+YD,ZSM-ZD,BX2X,BY2X,BZ2X,BX2Y,BY2Y,BZ2Y,
     *  BX2Z,BY2Z,BZ2Z)
      CALL DIPXYZ(XSM-XD,Y-YD,ZSM+ZD,BX3X,BY3X,BZ3X,BX3Y,BY3Y,BZ3Y,
     *  BX3Z,BY3Z,BZ3Z)
      CALL DIPXYZ(XSM-XD,Y+YD,ZSM+ZD,BX4X,BY4X,BZ4X,BX4Y,BY4Y,BZ4Y,
     *  BX4Z,BY4Z,BZ4Z)
C
      IX=I*3+3
      IY=IX+1
      IZ=IY+1
C
      D(1,IX)=(BX1X+BX2X-BX3X-BX4X)*CPS+(BZ1X+BZ2X-BZ3X-BZ4X)*SPS
      D(2,IX)= BY1X+BY2X-BY3X-BY4X
      D(3,IX)=(BZ1X+BZ2X-BZ3X-BZ4X)*CPS-(BX1X+BX2X-BX3X-BX4X)*SPS
C
      D(1,IY)=(BX1Y-BX2Y-BX3Y+BX4Y)*CPS+(BZ1Y-BZ2Y-BZ3Y+BZ4Y)*SPS
      D(2,IY)= BY1Y-BY2Y-BY3Y+BY4Y
      D(3,IY)=(BZ1Y-BZ2Y-BZ3Y+BZ4Y)*CPS-(BX1Y-BX2Y-BX3Y+BX4Y)*SPS
C
      D(1,IZ)=(BX1Z+BX2Z+BX3Z+BX4Z)*CPS+(BZ1Z+BZ2Z+BZ3Z+BZ4Z)*SPS
      D(2,IZ)= BY1Z+BY2Z+BY3Z+BY4Z
      D(3,IZ)=(BZ1Z+BZ2Z+BZ3Z+BZ4Z)*CPS-(BX1Z+BX2Z+BX3Z+BX4Z)*SPS
C
      IX=IX+27
      IY=IY+27
      IZ=IZ+27
C
      D(1,IX)=SPS*((BX1X+BX2X+BX3X+BX4X)*CPS+(BZ1X+BZ2X+BZ3X+BZ4X)*SPS)
      D(2,IX)=SPS*(BY1X+BY2X+BY3X+BY4X)
      D(3,IX)=SPS*((BZ1X+BZ2X+BZ3X+BZ4X)*CPS-(BX1X+BX2X+BX3X+BX4X)*SPS)
C
      D(1,IY)=SPS*((BX1Y-BX2Y+BX3Y-BX4Y)*CPS+(BZ1Y-BZ2Y+BZ3Y-BZ4Y)*SPS)
      D(2,IY)=SPS*(BY1Y-BY2Y+BY3Y-BY4Y)
      D(3,IY)=SPS*((BZ1Y-BZ2Y+BZ3Y-BZ4Y)*CPS-(BX1Y-BX2Y+BX3Y-BX4Y)*SPS)
C
      D(1,IZ)=SPS*((BX1Z+BX2Z-BX3Z-BX4Z)*CPS+(BZ1Z+BZ2Z-BZ3Z-BZ4Z)*SPS)
      D(2,IZ)=SPS*(BY1Z+BY2Z-BY3Z-BY4Z)
      D(3,IZ)=SPS*((BZ1Z+BZ2Z-BZ3Z-BZ4Z)*CPS-(BX1Z+BX2Z-BX3Z-BX4Z)*SPS)
  2   CONTINUE
C
      DO 3 I=1,5
      ZD=ZZ(I+9)
      CALL DIPXYZ(XSM,Y,ZSM-ZD,BX1X,BY1X,BZ1X,BX1Y,BY1Y,BZ1Y,BX1Z,BY1Z,
     *  BZ1Z)
      CALL DIPXYZ(XSM,Y,ZSM+ZD,BX2X,BY2X,BZ2X,BX2Y,BY2Y,BZ2Y,BX2Z,BY2Z,
     *  BZ2Z)
       IX=58+I*2
       IZ=IX+1
      D(1,IX)=(BX1X-BX2X)*CPS+(BZ1X-BZ2X)*SPS
      D(2,IX)=BY1X-BY2X
      D(3,IX)=(BZ1X-BZ2X)*CPS-(BX1X-BX2X)*SPS
C
      D(1,IZ)=(BX1Z+BX2Z)*CPS+(BZ1Z+BZ2Z)*SPS
      D(2,IZ)=BY1Z+BY2Z
      D(3,IZ)=(BZ1Z+BZ2Z)*CPS-(BX1Z+BX2Z)*SPS
C
      IX=IX+10
      IZ=IZ+10
      D(1,IX)=SPS*((BX1X+BX2X)*CPS+(BZ1X+BZ2X)*SPS)
      D(2,IX)=SPS*(BY1X+BY2X)
      D(3,IX)=SPS*((BZ1X+BZ2X)*CPS-(BX1X+BX2X)*SPS)
C
      D(1,IZ)=SPS*((BX1Z-BX2Z)*CPS+(BZ1Z-BZ2Z)*SPS)
      D(2,IZ)=SPS*(BY1Z-BY2Z)
  3   D(3,IZ)=SPS*((BZ1Z-BZ2Z)*CPS-(BX1Z-BX2Z)*SPS)
C
            RETURN
            END
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
         SUBROUTINE  BIRK1SHLD(PS,X,Y,Z,BX,BY,BZ)
C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
C  The 64 linear parameters are amplitudes of the "box" harmonics.
c The 16 nonlinear parameters are the scales Pi, and Qk entering the arguments
C  of sines/cosines and exponents in each of  32 cartesian harmonics
c  N.A. Tsyganenko, Spring 1994, adjusted for the Birkeland field Aug.22, 1995
c    Revised  June 12, 1996.
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
      IMPLICIT  REAL * 8  (A - H, O - Z)
C
      DIMENSION A(80)
      DIMENSION P1(4),R1(4),Q1(4),S1(4),RP(4),RR(4),RQ(4),RS(4)
C
      EQUIVALENCE (P1(1),A(65)),(R1(1),A(69)),(Q1(1),A(73)),
     * (S1(1),A(77))
C
      DATA A/1.174198045d0,-1.463820502d0,4.840161537d0,-3.674506864d0,
     * 82.18368896d0,-94.94071588d0,-4122.331796d0,4670.278676d0,
     : -21.54975037d0,26.72661293d0,-72.81365728d0,44.09887902d0,
     : 40.08073706d0,-51.23563510d0,1955.348537d0,-1940.971550d0,
     : 794.0496433d0,-982.2441344d0,1889.837171d0,-558.9779727d0,
     : -1260.543238d0,1260.063802d0,-293.5942373d0,344.7250789d0,
     * -773.7002492d0,957.0094135d0,-1824.143669d0,520.7994379d0,
     : 1192.484774d0,-1192.184565d0,89.15537624d0,-98.52042999d0,
     : -0.8168777675d-01,0.4255969908d-01,0.3155237661d0,
     : -0.3841755213d0,2.494553332d0,-0.6571440817d-01,-2.765661310d0,
     : 0.4331001908d0,0.1099181537d0,-0.6154126980d-01,-0.3258649260d0,
     : 0.6698439193d0,-5.542735524d0,0.1604203535d0,5.854456934d0,
     : -0.8323632049d0,3.732608869d0,-3.130002153d0,107.0972607d0,
     : -32.28483411d0,-115.2389298d0,54.45064360d0,-0.5826853320d0,
     * -3.582482231d0,-4.046544561d0,3.311978102d0,-104.0839563d0,
     : 30.26401293d0,97.29109008d0,-50.62370872d0,-296.3734955d0,
     : 127.7872523d0,5.303648988d0,10.40368955d0,69.65230348d0,
     : 466.5099509d0,1.645049286d0,3.825838190d0,11.66675599d0,
     : 558.9781177d0,1.826531343d0,2.066018073d0,25.40971369d0,
     * 990.2795225d0,2.319489258d0,4.555148484d0,9.691185703d0,
     : 591.8280358d0/
C
         BX=0.D0
         BY=0.D0
         BZ=0.D0
         CPS=DCOS(PS)
         SPS=DSIN(PS)
         S3PS=4.D0*CPS**2-1.D0
C
         DO 11 I=1,4
          RP(I)=1.D0/P1(I)
          RR(I)=1.D0/R1(I)
          RQ(I)=1.D0/Q1(I)
 11       RS(I)=1.D0/S1(I)
C
          L=0
C
           DO 1 M=1,2     
c!    M=1 IS FOR THE 1ST SUM ("PERP." SYMMETRY)
C                           AND M=2 IS FOR THE SECOND SUM ("PARALL." SYMMETRY)
             DO 2 I=1,4
                  CYPI=DCOS(Y*RP(I))
                  CYQI=DCOS(Y*RQ(I))
                  SYPI=DSIN(Y*RP(I))
                  SYQI=DSIN(Y*RQ(I))
C
                DO 3 K=1,4
                   SZRK=DSIN(Z*RR(K))
                   CZSK=DCOS(Z*RS(K))
                   CZRK=DCOS(Z*RR(K))
                   SZSK=DSIN(Z*RS(K))
                     SQPR=DSQRT(RP(I)**2+RR(K)**2)
                     SQQS=DSQRT(RQ(I)**2+RS(K)**2)
                        EPR=DEXP(X*SQPR)
                        EQS=DEXP(X*SQQS)
C
                    DO 4 N=1,2  
c! N=1 IS FOR THE FIRST PART OF EACH COEFFICIENT
C                                  AND N=2 IS FOR THE SECOND ONE
                     IF (M.EQ.1) THEN
                       IF (N.EQ.1) THEN
                         HX=-SQPR*EPR*CYPI*SZRK
                         HY=RP(I)*EPR*SYPI*SZRK
                         HZ=-RR(K)*EPR*CYPI*CZRK
                                   ELSE
                         HX=HX*CPS
                         HY=HY*CPS
                         HZ=HZ*CPS
                                   ENDIF
                     ELSE
                       IF (N.EQ.1) THEN
                         HX=-SPS*SQQS*EQS*CYQI*CZSK
                         HY=SPS*RQ(I)*EQS*SYQI*CZSK
                         HZ=SPS*RS(K)*EQS*CYQI*SZSK
                                   ELSE
                         HX=HX*S3PS
                         HY=HY*S3PS
                         HZ=HZ*S3PS
                       ENDIF
                 ENDIF
       L=L+1
c
       BX=BX+A(L)*HX
       BY=BY+A(L)*HY
  4    BZ=BZ+A(L)*HZ
  3   CONTINUE
  2   CONTINUE
  1   CONTINUE
C
         RETURN
       END
C
C##########################################################################
C
         SUBROUTINE BIRK2TOT_02(PS,X,Y,Z,BX,BY,BZ)
C
          IMPLICIT REAL*8 (A-H,O-Z)
C
          CALL BIRK2SHL(X,Y,Z,PS,WX,WY,WZ)
          CALL R2_BIRK(X,Y,Z,PS,HX,HY,HZ)
         BX=WX+HX
         BY=WY+HY
         BZ=WZ+HZ
         RETURN
         END
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C THIS CODE IS FOR THE FIELD FROM  2x2x2=8 "CARTESIAN" HARMONICS
C
         SUBROUTINE  BIRK2SHL(X,Y,Z,PS,HX,HY,HZ)
C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C    The model parameters are provided to this module via common-block /A/.
C  The 16 linear parameters enter in pairs in the amplitudes of the
c       "cartesian" harmonics.
c    The 8 nonlinear parameters are the scales Pi,Ri,Qi,and Si entering the
c  arguments of exponents, sines, and cosines in each of the 8 "Cartesian"
c   harmonics
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
       IMPLICIT  REAL * 8  (A - H, O - Z)
C
         DIMENSION P(2),R(2),Q(2),S(2)
         DIMENSION A(24)
C
         EQUIVALENCE(P(1),A(17)),(R(1),A(19)),(Q(1),A(21)),(S(1),A(23))
         DATA A/-111.6371348d0,124.5402702d0,110.3735178d0,
     : -122.0095905d0,111.9448247d0,-129.1957743d0,-110.7586562d0,
     : 126.5649012d0,-0.7865034384d0,-0.2483462721d0,0.8026023894d0,
     : 0.2531397188d0,10.72890902d0,0.8483902118d0,-10.96884315d0,
     : -0.8583297219d0,13.85650567d0,14.90554500d0,10.21914434d0,
     * 10.09021632d0,6.340382460d0,14.40432686d0,12.71023437d0,
     : 12.83966657d0/
C
            CPS=DCOS(PS)
            SPS=DSIN(PS)
            S3PS=4.D0*CPS**2-1.D0   
c!  THIS IS SIN(3*PS)/SIN(PS)
C
           HX=0.D0
           HY=0.D0
           HZ=0.D0
           L=0
C
           DO 1 M=1,2     
c!    M=1 IS FOR THE 1ST SUM ("PERP." SYMMETRY)
C                           AND M=2 IS FOR THE SECOND SUM ("PARALL." SYMMETRY)
             DO 2 I=1,2
                  CYPI=DCOS(Y/P(I))
                  CYQI=DCOS(Y/Q(I))
                  SYPI=DSIN(Y/P(I))
                  SYQI=DSIN(Y/Q(I))
C
               DO 3 K=1,2
                   SZRK=DSIN(Z/R(K))
                   CZSK=DCOS(Z/S(K))
                   CZRK=DCOS(Z/R(K))
                   SZSK=DSIN(Z/S(K))
                     SQPR=DSQRT(1.D0/P(I)**2+1.D0/R(K)**2)
                     SQQS=DSQRT(1.D0/Q(I)**2+1.D0/S(K)**2)
                        EPR=DEXP(X*SQPR)
                        EQS=DEXP(X*SQQS)
C
                   DO 4 N=1,2  
c! N=1 IS FOR THE FIRST PART OF EACH COEFFICIENT
C                                  AND N=2 IS FOR THE SECOND ONE
C
                    L=L+1
                     IF (M.EQ.1) THEN
                       IF (N.EQ.1) THEN
                         DX=-SQPR*EPR*CYPI*SZRK
                         DY=EPR/P(I)*SYPI*SZRK
                         DZ=-EPR/R(K)*CYPI*CZRK
                         HX=HX+A(L)*DX
                         HY=HY+A(L)*DY
                         HZ=HZ+A(L)*DZ
                                   ELSE
                         DX=DX*CPS
                         DY=DY*CPS
                         DZ=DZ*CPS
                         HX=HX+A(L)*DX
                         HY=HY+A(L)*DY
                         HZ=HZ+A(L)*DZ
                                   ENDIF
                     ELSE
                       IF (N.EQ.1) THEN
                         DX=-SPS*SQQS*EQS*CYQI*CZSK
                         DY=SPS*EQS/Q(I)*SYQI*CZSK
                         DZ=SPS*EQS/S(K)*CYQI*SZSK
                         HX=HX+A(L)*DX
                         HY=HY+A(L)*DY
                         HZ=HZ+A(L)*DZ
                                   ELSE
                         DX=DX*S3PS
                         DY=DY*S3PS
                         DZ=DZ*S3PS
                         HX=HX+A(L)*DX
                         HY=HY+A(L)*DY
                         HZ=HZ+A(L)*DZ
                       ENDIF
                 ENDIF
c
  4   CONTINUE
  3   CONTINUE
  2   CONTINUE
  1   CONTINUE
C
         RETURN
       END

c********************************************************************
C
       SUBROUTINE R2_BIRK(X,Y,Z,PS,BX,BY,BZ)
C
C  RETURNS THE MODEL FIELD FOR THE REGION 2 BIRKELAND CURRENT/PARTIAL RC
C    (WITHOUT SHIELDING FIELD)
C
       IMPLICIT REAL*8 (A-H,O-Z)
       SAVE PSI,CPS,SPS
       DATA DELARG/0.030D0/,DELARG1/0.015D0/,PSI/10.D0/
C
       IF (DABS(PSI-PS).GT.1.D-10) THEN
         PSI=PS
         CPS=DCOS(PS)
         SPS=DSIN(PS)
       ENDIF
C
       XSM=X*CPS-Z*SPS
       ZSM=Z*CPS+X*SPS
C
       XKS=XKSI(XSM,Y,ZSM)
      IF (XKS.LT.-(DELARG+DELARG1)) THEN
        CALL R2OUTER(XSM,Y,ZSM,BXSM,BY,BZSM)
         BXSM=-BXSM*0.02d0      
c!  ALL COMPONENTS ARE MULTIPLIED BY THE
       BY=-BY*0.02d0          
c!  FACTOR -0.02, IN ORDER TO NORMALIZE THE
         BZSM=-BZSM*0.02d0      
c!  FIELD (SO THAT Bz=-1 nT at X=-5.3 RE, Y=Z=0)
      ENDIF
      IF (XKS.GE.-(DELARG+DELARG1).AND.XKS.LT.-DELARG+DELARG1) THEN
        CALL R2OUTER(XSM,Y,ZSM,BXSM1,BY1,BZSM1)
        CALL R2SHEET(XSM,Y,ZSM,BXSM2,BY2,BZSM2)
        F2=-0.02d0*TKSI(XKS,-DELARG,DELARG1)
        F1=-0.02d0-F2
        BXSM=BXSM1*F1+BXSM2*F2
        BY=BY1*F1+BY2*F2
        BZSM=BZSM1*F1+BZSM2*F2
      ENDIF

      IF (XKS.GE.-DELARG+DELARG1.AND.XKS.LT.DELARG-DELARG1) THEN
       CALL R2SHEET(XSM,Y,ZSM,BXSM,BY,BZSM)
         BXSM=-BXSM*0.02d0
         BY=-BY*0.02d0
         BZSM=-BZSM*0.02d0
      ENDIF
      IF (XKS.GE.DELARG-DELARG1.AND.XKS.LT.DELARG+DELARG1) THEN
        CALL R2INNER(XSM,Y,ZSM,BXSM1,BY1,BZSM1)
        CALL R2SHEET(XSM,Y,ZSM,BXSM2,BY2,BZSM2)
        F1=-0.02d0*TKSI(XKS,DELARG,DELARG1)
        F2=-0.02d0-F1
        BXSM=BXSM1*F1+BXSM2*F2
        BY=BY1*F1+BY2*F2
        BZSM=BZSM1*F1+BZSM2*F2
      ENDIF
      IF (XKS.GE.DELARG+DELARG1) THEN
         CALL R2INNER(XSM,Y,ZSM,BXSM,BY,BZSM)
         BXSM=-BXSM*0.02d0
         BY=-BY*0.02d0
         BZSM=-BZSM*0.02d0
      ENDIF
C
        BX=BXSM*CPS+BZSM*SPS
        BZ=BZSM*CPS-BXSM*SPS
C
        RETURN
        END
C
C****************************************************************

c
      SUBROUTINE R2INNER (X,Y,Z,BX,BY,BZ)
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION CBX(5),CBY(5),CBZ(5)
C
      DATA PL1,PL2,PL3,PL4,PL5,PL6,PL7,PL8/154.185d0,-2.12446d0,
     : .601735d-01,-.153954d-02,.355077d-04,29.9996d0,262.886d0,
     :  99.9132d0/
      DATA PN1,PN2,PN3,PN4,PN5,PN6,PN7,PN8/-8.1902d0,6.5239d0,5.504d0,
     : 7.7815d0,.8573d0,3.0986d0,.0774d0,-.038d0/
C
      CALL BCONIC(X,Y,Z,CBX,CBY,CBZ,5)
C
C   NOW INTRODUCE  ONE  4-LOOP SYSTEM:
C
       CALL LOOPS4(X,Y,Z,DBX8,DBY8,DBZ8,PN1,PN2,PN3,PN4,PN5,PN6)
C
       CALL DIPDISTR(X-PN7,Y,Z,DBX6,DBY6,DBZ6,0)
       CALL DIPDISTR(X-PN8,Y,Z,DBX7,DBY7,DBZ7,1)

C                           NOW COMPUTE THE FIELD COMPONENTS:

      BX=PL1*CBX(1)+PL2*CBX(2)+PL3*CBX(3)+PL4*CBX(4)+PL5*CBX(5)
     * +PL6*DBX6+PL7*DBX7+PL8*DBX8
      BY=PL1*CBY(1)+PL2*CBY(2)+PL3*CBY(3)+PL4*CBY(4)+PL5*CBY(5)
     * +PL6*DBY6+PL7*DBY7+PL8*DBY8
      BZ=PL1*CBZ(1)+PL2*CBZ(2)+PL3*CBZ(3)+PL4*CBZ(4)+PL5*CBZ(5)
     * +PL6*DBZ6+PL7*DBZ7+PL8*DBZ8
C
      RETURN
      END
C-----------------------------------------------------------------------

      SUBROUTINE BCONIC(X,Y,Z,CBX,CBY,CBZ,NMAX)
C
c   "CONICAL" HARMONICS
c
       IMPLICIT REAL*8 (A-H,O-Z)
C
       DIMENSION CBX(NMAX),CBY(NMAX),CBZ(NMAX)

       RO2=X**2+Y**2
       RO=SQRT(RO2)
C
       CF=X/RO
       SF=Y/RO
       CFM1=1.D0
       SFM1=0.D0
C
      R2=RO2+Z**2
      R=DSQRT(R2)
      C=Z/R
      S=RO/R
      CH=DSQRT(0.5D0*(1.D0+C))
      SH=DSQRT(0.5D0*(1.D0-C))
      TNHM1=1.D0
      CNHM1=1.D0
      TNH=SH/CH
      CNH=1.D0/TNH
C
      DO 1 M=1,NMAX
        CFM=CFM1*CF-SFM1*SF
        SFM=CFM1*SF+SFM1*CF
        CFM1=CFM
        SFM1=SFM
        TNHM=TNHM1*TNH
        CNHM=CNHM1*CNH
       BT=M*CFM/(R*S)*(TNHM+CNHM)
       BF=-0.5D0*M*SFM/R*(TNHM1/CH**2-CNHM1/SH**2)
         TNHM1=TNHM
         CNHM1=CNHM
       CBX(M)=BT*C*CF-BF*SF
       CBY(M)=BT*C*SF+BF*CF
  1    CBZ(M)=-BT*S
C
       RETURN
       END

C-------------------------------------------------------------------
C
       SUBROUTINE DIPDISTR(X,Y,Z,BX,BY,BZ,MODE)
C
C   RETURNS FIELD COMPONENTS FROM A LINEAR DISTRIBUTION OF DIPOLAR SOURCES
C     ON THE Z-AXIS.  THE PARAMETER MODE DEFINES HOW THE DIPOLE STRENGTH
C     VARIES ALONG THE Z-AXIS:  MODE=0 IS FOR A STEP-FUNCTION (Mx=const &gt 0
c         FOR Z &gt 0, AND Mx=-const &lt 0 FOR Z &lt 0)
C      WHILE MODE=1 IS FOR A LINEAR VARIATION OF THE DIPOLE MOMENT DENSITY
C       SEE NB#3, PAGE 53 FOR DETAILS.
C
C
C INPUT: X,Y,Z OF A POINT OF SPACE, AND MODE
C
        IMPLICIT REAL*8 (A-H,O-Z)
        X2=X*X
        RHO2=X2+Y*Y
        R2=RHO2+Z*Z
        R3=R2*DSQRT(R2)

        IF (MODE.EQ.0) THEN
         BX=Z/RHO2**2*(R2*(Y*Y-X2)-RHO2*X2)/R3
         BY=-X*Y*Z/RHO2**2*(2.D0*R2+RHO2)/R3
         BZ=X/R3
        ELSE
         BX=Z/RHO2**2*(Y*Y-X2)
         BY=-2.D0*X*Y*Z/RHO2**2
         BZ=X/RHO2
        ENDIF
         RETURN
         END

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE R2OUTER (X,Y,Z,BX,BY,BZ)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DATA PL1,PL2,PL3,PL4,PL5/-34.105d0,-2.00019d0,628.639d0,
     : 73.4847d0,12.5162d0/
      DATA PN1,PN2,PN3,PN4,PN5,PN6,PN7,PN8,PN9,PN10,PN11,PN12,PN13,PN14,
     *  PN15,PN16,PN17 /.55d0,.694d0,.0031d0,1.55d0,2.8d0,.1375d0,
     : -.7d0,.2d0,.9625d0,-2.994d0,2.925d0,-1.775d0,4.3d0,-.275d0,
     : 2.7d0,.4312d0,1.55d0/
c
C    THREE PAIRS OF CROSSED LOOPS:
C
      CALL CROSSLP(X,Y,Z,DBX1,DBY1,DBZ1,PN1,PN2,PN3)
      CALL CROSSLP(X,Y,Z,DBX2,DBY2,DBZ2,PN4,PN5,PN6)
      CALL CROSSLP(X,Y,Z,DBX3,DBY3,DBZ3,PN7,PN8,PN9)
C
C    NOW AN EQUATORIAL LOOP ON THE NIGHTSIDE
C
      CALL CIRCLE(X-PN10,Y,Z,PN11,DBX4,DBY4,DBZ4)
c
c   NOW A 4-LOOP SYSTEM ON THE NIGHTSIDE
c

      CALL LOOPS4(X,Y,Z,DBX5,DBY5,DBZ5,PN12,PN13,PN14,PN15,PN16,PN17)

C---------------------------------------------------------------------

C                           NOW COMPUTE THE FIELD COMPONENTS:

      BX=PL1*DBX1+PL2*DBX2+PL3*DBX3+PL4*DBX4+PL5*DBX5
      BY=PL1*DBY1+PL2*DBY2+PL3*DBY3+PL4*DBY4+PL5*DBY5
      BZ=PL1*DBZ1+PL2*DBZ2+PL3*DBZ3+PL4*DBZ4+PL5*DBZ5

       RETURN
       END
C
C--------------------------------------------------------------------
C
       SUBROUTINE LOOPS4(X,Y,Z,BX,BY,BZ,XC,YC,ZC,R,THETA,PHI)
C
C   RETURNS FIELD COMPONENTS FROM A SYSTEM OF 4 CURRENT LOOPS, POSITIONED
C     SYMMETRICALLY WITH RESPECT TO NOON-MIDNIGHT MERIDIAN AND EQUATORIAL
C      PLANES.
C  INPUT: X,Y,Z OF A POINT OF SPACE
C        XC,YC,ZC (YC &gt 0 AND ZC &gt 0) - POSITION OF THE CENTER OF THE
C                                         1ST-QUADRANT LOOP
C        R - LOOP RADIUS (THE SAME FOR ALL FOUR)
C        THETA, PHI  -  SPECIFY THE ORIENTATION OF THE NORMAL OF THE 1ST LOOP
c      -----------------------------------------------------------

        IMPLICIT REAL*8 (A-H,O-Z)
C
          CT=DCOS(THETA)
          ST=DSIN(THETA)
          CP=DCOS(PHI)
          SP=DSIN(PHI)
C------------------------------------1ST QUADRANT:
        XS=(X-XC)*CP+(Y-YC)*SP
        YSS=(Y-YC)*CP-(X-XC)*SP
        ZS=Z-ZC
        XSS=XS*CT-ZS*ST
        ZSS=ZS*CT+XS*ST

        CALL CIRCLE(XSS,YSS,ZSS,R,BXSS,BYS,BZSS)
          BXS=BXSS*CT+BZSS*ST
          BZ1=BZSS*CT-BXSS*ST
          BX1=BXS*CP-BYS*SP
          BY1=BXS*SP+BYS*CP
C-------------------------------------2nd QUADRANT:
        XS=(X-XC)*CP-(Y+YC)*SP
        YSS=(Y+YC)*CP+(X-XC)*SP
        ZS=Z-ZC
        XSS=XS*CT-ZS*ST
        ZSS=ZS*CT+XS*ST

        CALL CIRCLE(XSS,YSS,ZSS,R,BXSS,BYS,BZSS)
          BXS=BXSS*CT+BZSS*ST
          BZ2=BZSS*CT-BXSS*ST
          BX2=BXS*CP+BYS*SP
          BY2=-BXS*SP+BYS*CP
C-------------------------------------3RD QUADRANT:
        XS=-(X-XC)*CP+(Y+YC)*SP
        YSS=-(Y+YC)*CP-(X-XC)*SP
        ZS=Z+ZC
        XSS=XS*CT-ZS*ST
        ZSS=ZS*CT+XS*ST

        CALL CIRCLE(XSS,YSS,ZSS,R,BXSS,BYS,BZSS)
          BXS=BXSS*CT+BZSS*ST
          BZ3=BZSS*CT-BXSS*ST
          BX3=-BXS*CP-BYS*SP
          BY3=BXS*SP-BYS*CP
C-------------------------------------4TH QUADRANT:
        XS=-(X-XC)*CP-(Y-YC)*SP
        YSS=-(Y-YC)*CP+(X-XC)*SP
        ZS=Z+ZC
        XSS=XS*CT-ZS*ST
        ZSS=ZS*CT+XS*ST

        CALL CIRCLE(XSS,YSS,ZSS,R,BXSS,BYS,BZSS)
          BXS=BXSS*CT+BZSS*ST
          BZ4=BZSS*CT-BXSS*ST
          BX4=-BXS*CP+BYS*SP
          BY4=-BXS*SP-BYS*CP

        BX=BX1+BX2+BX3+BX4
        BY=BY1+BY2+BY3+BY4
        BZ=BZ1+BZ2+BZ3+BZ4

         RETURN
         END
C
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C
      SUBROUTINE R2SHEET(X,Y,Z,BX,BY,BZ)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DATA PNONX1,PNONX2,PNONX3,PNONX4,PNONX5,PNONX6,PNONX7,PNONX8,
     *     PNONY1,PNONY2,PNONY3,PNONY4,PNONY5,PNONY6,PNONY7,PNONY8,
     *     PNONZ1,PNONZ2,PNONZ3,PNONZ4,PNONZ5,PNONZ6,PNONZ7,PNONZ8
     */-19.0969D0,-9.28828D0,-0.129687D0,5.58594D0,22.5055D0,
     *  0.483750D-01,0.396953D-01,0.579023D-01,-13.6750D0,-6.70625D0,
     *  2.31875D0,11.4062D0,20.4562D0,0.478750D-01,0.363750D-01,
     * 0.567500D-01,-16.7125D0,-16.4625D0,-0.1625D0,5.1D0,23.7125D0,
     * 0.355625D-01,0.318750D-01,0.538750D-01/
C
C
      DATA A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15,A16,A17,
     *  A18,A19,A20,A21,A22,A23,A24,A25,A26,A27,A28,A29,A30,A31,A32,A33,
     *  A34,A35,A36,A37,A38,A39,A40,A41,A42,A43,A44,A45,A46,A47,A48,A49,
     *  A50,A51,A52,A53,A54,A55,A56,A57,A58,A59,A60,A61,A62,A63,A64,A65,
     *  A66,A67,A68,A69,A70,A71,A72,A73,A74,A75,A76,A77,A78,A79,A80
     * /8.07190D0,-7.39582D0,-7.62341D0,0.684671D0,-13.5672D0,11.6681D0,
     * 13.1154,-0.890217D0,7.78726D0,-5.38346D0,-8.08738D0,0.609385D0,
     * -2.70410D0, 3.53741D0,3.15549D0,-1.11069D0,-8.47555D0,0.278122D0,
     *  2.73514D0,4.55625D0,13.1134D0,1.15848D0,-3.52648D0,-8.24698D0,
     * -6.85710D0,-2.81369D0, 2.03795D0, 4.64383D0,2.49309D0,-1.22041D0,
     * -1.67432D0,-0.422526D0,-5.39796D0,7.10326D0,5.53730D0,-13.1918D0,
     *  4.67853D0,-7.60329D0,-2.53066D0, 7.76338D0, 5.60165D0,5.34816D0,
     * -4.56441D0,7.05976D0,-2.62723D0,-0.529078D0,1.42019D0,-2.93919D0,
     *  55.6338D0,-1.55181D0,39.8311D0,-80.6561D0,-46.9655D0,32.8925D0,
     * -6.32296D0,19.7841D0,124.731D0,10.4347D0,-30.7581D0,102.680D0,
     * -47.4037D0,-3.31278D0,9.37141D0,-50.0268D0,-533.319D0,110.426D0,
     *  1000.20D0,-1051.40D0, 1619.48D0,589.855D0,-1462.73D0,1087.10D0,
     *  -1994.73D0,-1654.12D0,1263.33D0,-260.210D0,1424.84D0,1255.71D0,
     *  -956.733D0, 219.946D0/
C
C
      DATA B1,B2,B3,B4,B5,B6,B7,B8,B9,B10,B11,B12,B13,B14,B15,B16,B17,
     *  B18,B19,B20,B21,B22,B23,B24,B25,B26,B27,B28,B29,B30,B31,B32,B33,
     *  B34,B35,B36,B37,B38,B39,B40,B41,B42,B43,B44,B45,B46,B47,B48,B49,
     *  B50,B51,B52,B53,B54,B55,B56,B57,B58,B59,B60,B61,B62,B63,B64,B65,
     *  B66,B67,B68,B69,B70,B71,B72,B73,B74,B75,B76,B77,B78,B79,B80
     */-9.08427D0,10.6777D0,10.3288D0,-0.969987D0,6.45257D0,-8.42508D0,
     * -7.97464D0,1.41996D0,-1.92490D0,3.93575D0,2.83283D0,-1.48621D0,
     *0.244033D0,-0.757941D0,-0.386557D0,0.344566D0,9.56674D0,-2.5365D0,
     * -3.32916D0,-5.86712D0,-6.19625D0,1.83879D0,2.52772D0,4.34417D0,
     * 1.87268D0,-2.13213D0,-1.69134D0,-.176379D0,-.261359D0,.566419D0,
     * 0.3138D0,-0.134699D0,-3.83086D0,-8.4154D0,4.77005D0,-9.31479D0,
     * 37.5715D0,19.3992D0,-17.9582D0,36.4604D0,-14.9993D0,-3.1442D0,
     * 6.17409D0,-15.5519D0,2.28621D0,-0.891549D-2,-.462912D0,2.47314D0,
     * 41.7555D0,208.614D0,-45.7861D0,-77.8687D0,239.357D0,-67.9226D0,
     * 66.8743D0,238.534D0,-112.136D0,16.2069D0,-40.4706D0,-134.328D0,
     * 21.56D0,-0.201725D0,2.21D0,32.5855D0,-108.217D0,-1005.98D0,
     * 585.753D0,323.668D0,-817.056D0,235.750D0,-560.965D0,-576.892D0,
     * 684.193D0,85.0275D0,168.394D0,477.776D0,-289.253D0,-123.216D0,
     * 75.6501D0,-178.605D0/
C
      DATA C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,C13,C14,C15,C16,C17,
     *  C18,C19,C20,C21,C22,C23,C24,C25,C26,C27,C28,C29,C30,C31,C32,C33,
     *  C34,C35,C36,C37,C38,C39,C40,C41,C42,C43,C44,C45,C46,C47,C48,C49,
     *  C50,C51,C52,C53,C54,C55,C56,C57,C58,C59,C60,C61,C62,C63,C64,C65,
     *  C66,C67,C68,C69,C70,C71,C72,C73,C74,C75,C76,C77,C78,C79,C80
     * / 1167.61D0,-917.782D0,-1253.2D0,-274.128D0,-1538.75D0,1257.62D0,
     * 1745.07D0,113.479D0,393.326D0,-426.858D0,-641.1D0,190.833D0,
     * -29.9435D0,-1.04881D0,117.125D0,-25.7663D0,-1168.16D0,910.247D0,
     * 1239.31D0,289.515D0,1540.56D0,-1248.29D0,-1727.61D0,-131.785D0,
     * -394.577D0,426.163D0,637.422D0,-187.965D0,30.0348D0,0.221898D0,
     * -116.68D0,26.0291D0,12.6804D0,4.84091D0,1.18166D0,-2.75946D0,
     * -17.9822D0,-6.80357D0,-1.47134D0,3.02266D0,4.79648D0,0.665255D0,
     * -0.256229D0,-0.857282D-1,-0.588997D0,0.634812D-1,0.164303D0,
     * -0.15285D0,22.2524D0,-22.4376D0,-3.85595D0,6.07625D0,-105.959D0,
     * -41.6698D0,0.378615D0,1.55958D0,44.3981D0,18.8521D0,3.19466D0,
     *  5.89142D0,-8.63227D0,-2.36418D0,-1.027D0,-2.31515D0,1035.38D0,
     *  2040.66D0,-131.881D0,-744.533D0,-3274.93D0,-4845.61D0,482.438D0,
     * 1567.43D0,1354.02D0,2040.47D0,-151.653D0,-845.012D0,-111.723D0,
     * -265.343D0,-26.1171D0,216.632D0/
C
c------------------------------------------------------------------
C
       XKS=XKSI(X,Y,Z)    
c!  variation across the current sheet
       T1X=XKS/DSQRT(XKS**2+PNONX6**2)
       T2X=PNONX7**3/DSQRT(XKS**2+PNONX7**2)**3
       T3X=XKS/DSQRT(XKS**2+PNONX8**2)**5 *3.493856D0*PNONX8**4
C
       T1Y=XKS/DSQRT(XKS**2+PNONY6**2)
       T2Y=PNONY7**3/DSQRT(XKS**2+PNONY7**2)**3
       T3Y=XKS/DSQRT(XKS**2+PNONY8**2)**5 *3.493856D0*PNONY8**4
C
       T1Z=XKS/DSQRT(XKS**2+PNONZ6**2)
       T2Z=PNONZ7**3/DSQRT(XKS**2+PNONZ7**2)**3
       T3Z=XKS/DSQRT(XKS**2+PNONZ8**2)**5 *3.493856D0*PNONZ8**4
C
      RHO2=X*X+Y*Y
      R=DSQRT(RHO2+Z*Z)
      RHO=DSQRT(RHO2)
C
      C1P=X/RHO
      S1P=Y/RHO
      S2P=2.D0*S1P*C1P
      C2P=C1P*C1P-S1P*S1P
      S3P=S2P*C1P+C2P*S1P
      C3P=C2P*C1P-S2P*S1P
      S4P=S3P*C1P+C3P*S1P
      CT=Z/R
      ST=RHO/R
C
      S1=FEXP(CT,PNONX1)
      S2=FEXP(CT,PNONX2)
      S3=FEXP(CT,PNONX3)
      S4=FEXP(CT,PNONX4)
      S5=FEXP(CT,PNONX5)
C
C                   NOW COMPUTE THE GSM FIELD COMPONENTS:
C
C
      BX=S1*((A1+A2*T1X+A3*T2X+A4*T3X)
     *        +C1P*(A5+A6*T1X+A7*T2X+A8*T3X)
     *        +C2P*(A9+A10*T1X+A11*T2X+A12*T3X)
     *        +C3P*(A13+A14*T1X+A15*T2X+A16*T3X))
     *    +S2*((A17+A18*T1X+A19*T2X+A20*T3X)
     *        +C1P*(A21+A22*T1X+A23*T2X+A24*T3X)
     *        +C2P*(A25+A26*T1X+A27*T2X+A28*T3X)
     *        +C3P*(A29+A30*T1X+A31*T2X+A32*T3X))
     *    +S3*((A33+A34*T1X+A35*T2X+A36*T3X)
     *        +C1P*(A37+A38*T1X+A39*T2X+A40*T3X)
     *        +C2P*(A41+A42*T1X+A43*T2X+A44*T3X)
     *        +C3P*(A45+A46*T1X+A47*T2X+A48*T3X))
     *    +S4*((A49+A50*T1X+A51*T2X+A52*T3X)
     *        +C1P*(A53+A54*T1X+A55*T2X+A56*T3X)
     *        +C2P*(A57+A58*T1X+A59*T2X+A60*T3X)
     *        +C3P*(A61+A62*T1X+A63*T2X+A64*T3X))
     *    +S5*((A65+A66*T1X+A67*T2X+A68*T3X)
     *        +C1P*(A69+A70*T1X+A71*T2X+A72*T3X)
     *        +C2P*(A73+A74*T1X+A75*T2X+A76*T3X)
     *        +C3P*(A77+A78*T1X+A79*T2X+A80*T3X))
C
C
      S1=FEXP(CT,PNONY1)
      S2=FEXP(CT,PNONY2)
      S3=FEXP(CT,PNONY3)
      S4=FEXP(CT,PNONY4)
      S5=FEXP(CT,PNONY5)
C
      BY=S1*(S1P*(B1+B2*T1Y+B3*T2Y+B4*T3Y)
     *      +S2P*(B5+B6*T1Y+B7*T2Y+B8*T3Y)
     *      +S3P*(B9+B10*T1Y+B11*T2Y+B12*T3Y)
     *      +S4P*(B13+B14*T1Y+B15*T2Y+B16*T3Y))
     *  +S2*(S1P*(B17+B18*T1Y+B19*T2Y+B20*T3Y)
     *      +S2P*(B21+B22*T1Y+B23*T2Y+B24*T3Y)
     *      +S3P*(B25+B26*T1Y+B27*T2Y+B28*T3Y)
     *      +S4P*(B29+B30*T1Y+B31*T2Y+B32*T3Y))
     *  +S3*(S1P*(B33+B34*T1Y+B35*T2Y+B36*T3Y)
     *      +S2P*(B37+B38*T1Y+B39*T2Y+B40*T3Y)
     *      +S3P*(B41+B42*T1Y+B43*T2Y+B44*T3Y)
     *      +S4P*(B45+B46*T1Y+B47*T2Y+B48*T3Y))
     *  +S4*(S1P*(B49+B50*T1Y+B51*T2Y+B52*T3Y)
     *      +S2P*(B53+B54*T1Y+B55*T2Y+B56*T3Y)
     *      +S3P*(B57+B58*T1Y+B59*T2Y+B60*T3Y)
     *      +S4P*(B61+B62*T1Y+B63*T2Y+B64*T3Y))
     *  +S5*(S1P*(B65+B66*T1Y+B67*T2Y+B68*T3Y)
     *      +S2P*(B69+B70*T1Y+B71*T2Y+B72*T3Y)
     *      +S3P*(B73+B74*T1Y+B75*T2Y+B76*T3Y)
     *      +S4P*(B77+B78*T1Y+B79*T2Y+B80*T3Y))
C
      S1=FEXP1(CT,PNONZ1)
      S2=FEXP1(CT,PNONZ2)
      S3=FEXP1(CT,PNONZ3)
      S4=FEXP1(CT,PNONZ4)
      S5=FEXP1(CT,PNONZ5)
C
      BZ=S1*((C1+C2*T1Z+C3*T2Z+C4*T3Z)
     *      +C1P*(C5+C6*T1Z+C7*T2Z+C8*T3Z)
     *      +C2P*(C9+C10*T1Z+C11*T2Z+C12*T3Z)
     *      +C3P*(C13+C14*T1Z+C15*T2Z+C16*T3Z))
     *   +S2*((C17+C18*T1Z+C19*T2Z+C20*T3Z)
     *      +C1P*(C21+C22*T1Z+C23*T2Z+C24*T3Z)
     *      +C2P*(C25+C26*T1Z+C27*T2Z+C28*T3Z)
     *      +C3P*(C29+C30*T1Z+C31*T2Z+C32*T3Z))
     *   +S3*((C33+C34*T1Z+C35*T2Z+C36*T3Z)
     *      +C1P*(C37+C38*T1Z+C39*T2Z+C40*T3Z)
     *      +C2P*(C41+C42*T1Z+C43*T2Z+C44*T3Z)
     *      +C3P*(C45+C46*T1Z+C47*T2Z+C48*T3Z))
     *   +S4*((C49+C50*T1Z+C51*T2Z+C52*T3Z)
     *      +C1P*(C53+C54*T1Z+C55*T2Z+C56*T3Z)
     *      +C2P*(C57+C58*T1Z+C59*T2Z+C60*T3Z)
     *      +C3P*(C61+C62*T1Z+C63*T2Z+C64*T3Z))
     *   +S5*((C65+C66*T1Z+C67*T2Z+C68*T3Z)
     *      +C1P*(C69+C70*T1Z+C71*T2Z+C72*T3Z)
     *      +C2P*(C73+C74*T1Z+C75*T2Z+C76*T3Z)
     *      +C3P*(C77+C78*T1Z+C79*T2Z+C80*T3Z))
C
       RETURN
       END
C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      REAL*8 FUNCTION XKSI(X,Y,Z)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
C   A11 - C72, R0, and DR below  ARE STRETCH PARAMETERS (P.26-27, NB# 3),
C
      DATA A11A12,A21A22,A41A42,A51A52,A61A62,B11B12,B21B22,C61C62,
     *  C71C72,R0,DR /0.305662d0,-0.383593d0,0.2677733d0,-0.097656d0,
     :  -0.636034d0,-0.359862d0,0.424706d0,-0.126366d0,0.292578d0,
     :  1.21563d0,7.50937d0/

      DATA TNOON,DTETA/0.3665191d0,0.09599309d0/ 
c! Correspond to noon and midnight
C                                         latitudes 69 and 63.5 degs, resp.
       DR2=DR*DR
C
       X2=X*X
       Y2=Y*Y
       Z2=Z*Z
       XY=X*Y
       XYZ=XY*Z
       R2=X2+Y2+Z2
       R=DSQRT(R2)
       R3=R2*R
       R4=R2*R2
       XR=X/R
       YR=Y/R
       ZR=Z/R
C
       IF (R.LT.R0) THEN
         PR=0.D0
       ELSE
         PR=DSQRT((R-R0)**2+DR2)-DR
       ENDIF
C
      F=X+PR*(A11A12+A21A22*XR+A41A42*XR*XR+A51A52*YR*YR+
     *        A61A62*ZR*ZR)
      G=Y+PR*(B11B12*YR+B21B22*XR*YR)
      H=Z+PR*(C61C62*ZR+C71C72*XR*ZR)
      G2=G*G
C
      FGH=F**2+G2+H**2
      FGH32=DSQRT(FGH)**3
      FCHSG2=F**2+G2

      IF (FCHSG2.LT.1.D-5) THEN
         XKSI=-1.D0               
c!  THIS IS JUST FOR ELIMINATING PROBLEMS
         RETURN                    
c!  ON THE Z-AXIS
      ENDIF

      SQFCHSG2=DSQRT(FCHSG2)
      ALPHA=FCHSG2/FGH32
      THETA=TNOON+0.5D0*DTETA*(1.D0-F/SQFCHSG2)
      PHI=DSIN(THETA)**2
C
      XKSI=ALPHA-PHI
C
      RETURN
      END
C
C--------------------------------------------------------------------
C
        FUNCTION FEXP(S,A)
         IMPLICIT REAL*8 (A-H,O-Z)
          DATA E/2.718281828459D0/
          IF (A.LT.0.D0) FEXP=DSQRT(-2.D0*A*E)*S*DEXP(A*S*S)
          IF (A.GE.0.D0) FEXP=S*DEXP(A*(S*S-1.D0))
         RETURN
         END
C
C-----------------------------------------------------------------------
        FUNCTION FEXP1(S,A)
         IMPLICIT REAL*8 (A-H,O-Z)
         IF (A.LE.0.D0) FEXP1=DEXP(A*S*S)
         IF (A.GT.0.D0) FEXP1=DEXP(A*(S*S-1.D0))
         RETURN
         END
C
C************************************************************************
C
         REAL*8 FUNCTION TKSI(XKSI,XKS0,DXKSI)
         IMPLICIT REAL*8 (A-H,O-Z)
         SAVE M,TDZ3
         DATA M/0/
C
         IF (M.EQ.0) THEN
         TDZ3=2.d0*DXKSI**3
         M=1
         ENDIF
C
         IF (XKSI-XKS0.LT.-DXKSI) TKSII=0.d0
         IF (XKSI-XKS0.GE.DXKSI)  TKSII=1.d0
C
         IF (XKSI.GE.XKS0-DXKSI.AND.XKSI.LT.XKS0) THEN
           BR3=(XKSI-XKS0+DXKSI)**3
           TKSII=1.5d0*BR3/(TDZ3+BR3)
         ENDIF
C
         IF (XKSI.GE.XKS0.AND.XKSI.LT.XKS0+DXKSI) THEN
           BR3=(XKSI-XKS0-DXKSI)**3
           TKSII=1.d0+1.5d0*BR3/(TDZ3-BR3)
         ENDIF
           TKSI=TKSII
         END
C
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
       SUBROUTINE DIPOLE_96(PS,X,Y,Z,BX,BY,BZ)
C
C  CALCULATES GSM COMPONENTS OF GEODIPOLE FIELD WITH THE DIPOLE MOMENT
C  CORRESPONDING TO THE EPOCH OF 1980.
C------------INPUT PARAMETERS:
C   PS - GEODIPOLE TILT ANGLE IN RADIANS, X,Y,Z - GSM COORDINATES IN RE
C------------OUTPUT PARAMETERS:
C   BX,BY,BZ - FIELD COMPONENTS IN GSM SYSTEM, IN NANOTESLA.
C
C
C                   AUTHOR: NIKOLAI A. TSYGANENKO
C                           INSTITUTE OF PHYSICS
C                           ST.-PETERSBURG STATE UNIVERSITY
C                           STARY PETERGOF 198904
C                           ST.-PETERSBURG
C                           RUSSIA
C
      IMPLICIT NONE
C
      REAL*8 PS,X,Y,Z,BX,BY,BZ,PSI,SPS,CPS,P,U,V,T,Q
      INTEGER M
      SAVE M,PSI,SPS,CPS

      DATA M,PSI/0,5.d0/
      IF(M.EQ.1.AND.ABS(PS-PSI).LT.1.d-5) GOTO 1
      SPS=SIN(PS)
      CPS=COS(PS)
      PSI=PS
      M=1
  1   P=X**2
      U=Z**2
      V=3.d0*Z*X
      T=Y**2
      Q=30574.d0/SQRT(P+T+U)**5
      BX=Q*((T+U-2.d0*P)*SPS-V*CPS)
      BY=-3.d0*Y*Q*(X*SPS+Z*CPS)
      BZ=Q*((P+T-2.d0*U)*CPS-V*SPS)
      RETURN
      END
c
c
c
c
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE BAZ_T(X,Y,Z,N,B)
C     BAZIS FOR TILTE-MODEL(S) MAGNETIC FIELD
C     X,Y,Z = X,Y,Z/10.D0 !  IN NORMALIZATION UNITS
      INTEGER I,J,K,L,N
      REAL*8 X,Y,Z,B(3,1),SYM(3,10),ASY(3,10),JWORK(3,10)
      K=0
C 1-3:******** 1.Usym(1,3,5  ) FOR N=17 ******
C 1-4:******** 1.Usym(1,3,5,7) FOR N=29 ******
      L=3
      IF(N.EQ.29) L=4
      CALL PTNCL(X,Y,Z,7,SYM,ASY)   ! CALL FIRST TIME ONLY
      DO I=1,L
         DO J=1,3
            B(J,I+K)=SYM(J,2*I-1)
         ENDDO
      ENDDO
      K=K+L !  K=3/4
C 4- 6:******* 2.JETsym(1:3) FOR N=17 *********************
C 5-10:******* 2.JETsym(1:6) FOR N=29 *********************
      L=3
      IF(N.EQ.29) L=6
      CALL JETSYM(X,Y,Z,L,JWORK)
      DO I=1,L
         DO J=1,3
            B(J,I+K)=JWORK(J,I)
         ENDDO
      ENDDO
      K=K+L !  K=6/10
C   7,8:***** 3.Uasy(2,4)   FOR N=17 *************
C 11-13:***** 3.Uasy(2,4,6) FOR N=29 *************
      L=2
      IF(N.EQ.29) L=3
      DO I=1,L
         DO J=1,3
            B(J,I+K)=ASY(J,2*I)
         ENDDO
      ENDDO
      K=K+L !  K=8/13
C  9-12:***** 4.JETSasy(1:4) FOR N=17 *******************
C 14-22:***** 4.JETSasy(1:9) FOR N=29 *******************
      L=4
      IF(N.EQ.29) L=9
      CALL JETASY(X,Y,Z,L,JWORK)
      DO I=1,L
         DO J=1,3
            B(J,I+K)=JWORK(J,I)
         ENDDO
      ENDDO
      K=K+L !  K=12/22
C 13,14:**** 5.Usym(2,4)*sin(PSI)   FOR N=17 ************
C 23-25:**** 5.Usym(2,4,6)*sin(PSI) FOR N=29 ************
      L=2
      IF(N.EQ.29) L=3
      DO I=1,L
         DO J=1,3
            B(J,I+K)=SYM(J,2*I)
         ENDDO
      ENDDO
      K=K+L !  K=14/25
C 15-17:*** 6.Uasy(1,3,5)*sin(PSI)   FOR N=17 **********
C 26-29:*** 6.Uasy(1,3,5,7)*sin(PSI) FOR N=29 **********
      L=3
      IF(N.EQ.29) L=4
      DO I=1,L
         DO J=1,3
            B(J,I+K)=ASY(J,2*I-1)
         ENDDO
      ENDDO
      K=K+L !  K=17/29
      IF(K.EQ.N) RETURN
      WRITE(*,*) ' ERROR GENERATED IN EXT530, OSTAPENKO-MALTSEV 1997'
      WRITE(*,*) ' ERROR IN MODULE BAZ_T K#MF:',K,'#',N
      STOP
      END
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE JETSYM(X,Y,Z,L,B)
C     SYMMETRIC JET
      REAL*8 X,Y,Z,B(3,*),Z2,P2
      INTEGER I,L,J
      DO J=1,L
         DO I=1,3
            B(I,J)=0.D0
         ENDDO
      ENDDO
      P2=X*X+Y*Y
      Z2=Z*Z
      B(3,1)=P2
      B(3,2)=P2*P2
      B(2,3)=-2.D0*Z*Z2    !  #3: [-2pz3,z4]
      B(3,3)=Z2*Z2
      IF(L.EQ.6) THEN
         B(3,4)=P2**3
         B(2,5)=-Z*Z2*(P2-2.D0*Z2/5.D0) !       revised 28 april 1997
         B(3,5)=Z2*Z2*(P2-2.D0*Z2/15.D0)  !     revised 28 april 1997
         B(2,6)=-3.D0*Z**5 !  #6: [-3pz5,z6]
         B(3,6)=Z2**3
      ENDIF
      DO J=1,L
         B(1,J)=X*B(2,J)
         B(2,J)=Y*B(2,J)
      ENDDO
      RETURN
      END
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE JETASY(X,Y,Z,L,B)
C     ASYMMETRIC JET ~R3 - REVISED 06/06/95
      REAL*8 X,Y,Z,B(3,*),P2,Z2,WORK
      INTEGER L,I,J
      DO J=1,L
         DO I=1,3
            B(I,J)=0.D0
         ENDDO
      ENDDO
      P2=X*X+Y*Y
      Z2=Z*Z
      B(1,1)=Z
      B(3,2)=P2*X
      B(1,3)=Z2*Z
      B(1,4)=-Y*Y*Z
      B(2,4)= X*Y*Z
      B(3,4)=-X*Z2/2.D0
      IF(L.EQ.4) RETURN
C     ASYMMETRIC JET ~R5 - REVISED 16/06/95, NORMALIZATION 14/05/97
      B(3,5)=P2*P2*X       !  NORMALIZATION = 1
      WORK=3.D0*P2*Z       !  NORMALIZATION = 3
      B(1,6)=-Y*Y*WORK
      B(2,6)= X*Y*WORK
      B(3,6)=-X*Z*WORK/2.D0
      WORK=10.D0*P2*X*Z    !  NORMALIZATION = 10
      B(1,7)= X*WORK/5.D0
      B(2,7)= Y*WORK/5.D0
      B(3,7)=-Z*WORK/2.D0
      WORK=5.D0*Z2*Z       !  NORMALIZATION = 5
      B(1,8)=-Y*Y*WORK
      B(2,8)= X*Y*WORK
      B(3,8)=-X*Z*WORK/4.D0
      WORK=10.D0*X*Z*Z2
      B(1,9)= X*WORK/3.D0  !  NORMALIZATION = 10
      B(2,9)= Y*WORK/3.D0
      B(3,9)=-Z*WORK/4.D0
      RETURN
      END
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE PTNCL(X,Y,Z,N,SYM,ASY)
C     DERIVES grad [ P(i,j)*r^i]
      REAL*8 X,Y,Z,SYM(3,7),ASY(3,7),P(0:8,0:8),WORK,R2,C,R(-1:10),RN
      INTEGER N,I,J
C
      IF(N.GT.7.OR.N.LT.1) THEN
         WRITE(*,*) ' ERROR GENERATED IN EXT530, OSTAPENKO-MALTSEV 1997'
         WRITE(*,'(26H PTNCL: WRONG PARAMETER N=I2)') N
         STOP
      ENDIF
      DO J=1,N
         DO I=1,3
            SYM(I,J)=0.D0
            ASY(I,J)=0.D0
         ENDDO
      ENDDO
      SYM(3,1)=1.D0
      ASY(1,1)=1.D0
      R2=X*X+Y*Y+Z*Z
      IF(R2.EQ.0.D0) RETURN
      R(0)=1.D0
      R(1)=SQRT(R2)
      R(2)=R2
      R(-1)=1.D0/R(1)
      DO I=3,N+3
         R(I)=R(I-1)*R(1)
      ENDDO
      C=Z/R(1)
      CALL LEGNDR(C,N,P)
      DO I=2,N
C ZERO GARMONICS
         WORK=R(I-2)*(I*P(I,0)-C*P(I,1))
         SYM(1,I)=X*WORK                  ! d/dx
         SYM(2,I)=Y*WORK                  ! d/dy
         SYM(3,I)=Z*WORK+R(I-1)*P(I,1)    ! d/dz
C FIRST GARMONICS
         WORK=X*R(I-3)*((I-1)*P(I,1)-C*P(I,2))
         ASY(1,I)=X*WORK+R(I-1)*P(I,1)    ! d/dx
         ASY(2,I)=Y*WORK                  ! d/dy
         ASY(3,I)=Z*WORK+X*R(I-2)*P(I,2)  ! d/dz
      ENDDO
C SCHMIDT NORMALIZATION FOR M=1,N=2,3,...
      WORK=1.D0
      DO I=2,N
         WORK=(I-1)*WORK/(I+1)
         RN=SQRT(WORK)
         DO J=1,3
            ASY(J,I)=ASY(J,I)*RN
         ENDDO
      ENDDO
      RETURN
      END
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE LEGNDR(X,N,P)
C     P(i,0) = LEGENDRE'S POLYNOMS Pi(x)
C     AND ITS DERIVES: P(i,j)=(d/dx)^j (Pi(x))
      REAL*8 X,P(0:8,0:8)
      INTEGER N,I,J
      IF(N.GT.7.OR.N.LT.1) THEN
         WRITE(*,*) ' ERROR GENERATED IN EXT530, OSTAPENKO-MALTSEV 1997'
         WRITE(*,'(37H LEGENDRE POLYNOM: WRONG PARAMETER N=I2)') N
         STOP
      ENDIF
C
      P(0,0)=1.D0
      P(0,1)=0.D0
      P(1,0)=X
      P(1,1)=1.D0
      DO I=1,N-1
         P(I,I+1)=0.D0
         DO J=I+1,1,-1
            P(I+1,J)=X*P(I,J)+(I+J)*P(I,J-1)
         ENDDO
         P(I+1,0)=((2*I+1)*X*P(I,0)-I*P(I-1,0))/(I+1)
      ENDDO
      RETURN
      END
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE BOM97(RE,B)
C     EXTERNAL MODEL MAGNETIC FIELD (WITH TILTE)
      INTEGER MF,I,J,NA
      PARAMETER (MF=29)
      REAL*8 A(MF),RE(3),R(3),X,Y,Z,B(3),BAZIS(3,MF)
      COMMON/COEFOM97/A,NA
C
      DO J=1,3
         R(J)=RE(J)/10.D0 !  TO NORMALIZATION UNITS
         B(J)=0.D0
      ENDDO
      X=R(1)
      Y=R(2)
      Z=R(3)
      CALL BAZ_T(X,Y,Z,NA,BAZIS)
      DO I=1,NA
         DO J=1,3
            B(J)=B(J)+BAZIS(J,I)*A(I)
         ENDDO
      ENDDO
      RETURN
      END
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE SET_A(Dst,Pdyn,Kp,IMFz,SN)
      INTEGER MF,I,J,K,NC,L,NA
      PARAMETER (MF=29,NC=4)
      REAL*8 AA(85),C(4),SN,A(MF),Dst,Pdyn,Kp,IMFz
      COMMON/COEFOM97/A,NA
C Input parameters:
C Dst (nT), Pdyn (nPa),Kp (numeric),IMFz (nT), SN=sin(tilt of dipole)
C Kp (as key )=  0,  0+,  1-,   1,  1+,...
C Kp (numeric)=0.0, 0.3, 0.7, 1.0, 1.3,...
      DATA AA/
     ,-4.3392E+01, 2.0364E+01,-7.4532E-01,-3.2941E+00,-1.6294E+00,
     , 4.0308E+01,-1.3734E+01, 6.4463E+00, 3.4054E+00, 2.8944E+00,
     , 7.0475E+00,-3.2920E+00, 2.2886E+00, 1.0694E-01, 8.0280E-01,
     , 1.3347E+02,-6.1188E+01, 2.5064E+01, 2.5815E+00, 1.1370E+01,
     ,-3.6968E+01, 2.4557E+01,-1.2768E+01, 4.1262E+00,-4.2696E+00,
     ,-1.3010E+02, 4.6964E+01,-2.5640E+01,-1.4840E+01,-9.6417E+00,
     , 2.3632E+01, 4.9102E+00, 7.9419E+00, 4.3164E+00, 2.6925E+00,
     ,-1.7139E+00,-7.4376E-01, 2.8973E-01,-9.6157E-01,-5.0774E-01,
     , 1.2742E+00,-1.5849E+01, 7.1720E+00,-1.6217E+00,-1.0925E+01,
     ,-2.1437E+01,-1.0040E+01,-7.1910E+00,-9.9341E+00,-5.0171E+00,
     , 4.4863E+00, 3.1827E+00,-9.6490E+00,-6.9168E-01, 3.4412E+00,
     , 2.3315E+01,-6.8824E+00, 1.2323E+01, 3.5851E+00,-6.1224E+00,
     , 2.3179E+01,-2.0127E+00, 9.0600E+00,-1.8949E-01, 1.5000E-01,
     ,-1.9672E+00,-1.3847E+00, 2.9655E+00,-3.4092E+00, 1.4149E+00,
     , 1.2970E+01,-4.3018E+00, 6.4866E+00,-1.2664E-01,-2.1103E+00,
     , 5.7215E+00,-3.1930E+00, 2.2242E+00, 3.5120E-01, 9.1624E-02,
     , 4.9085E+00,-1.0754E+00, 3.5705E+00,-9.9892E-01,-6.2140E-01/
C     Input parameters: Dst (nT), Pdyn (nPa),Kp (numeric),IMFz (nT)
C     Convert from values to normalized parameters:
      C(1)=(Dst+16.9367D0)/25.2834D0
      C(2)=(Pdyn-2.278138D0)/1.882804D0
      C(3)=(Kp-2.30896D0)/1.35401D0
      C(4)=(IMFz-.0180D0)/3.7051D0
      NA=17
C     WRITE(*,'(24H NUMBER OF COEFFICIENTS=I3)') NA
      DO I=1,NA
         K=(I-1)*(NC+1)+1           !!!! REVISED 24/12/97
         A(I)=AA(K)                 !!!! REVISED 24/12/97
         DO J=1,NC
            A(I)=A(I)+AA(K+J)*C(J)  !!!! REVISED 24/12/97
         ENDDO
      ENDDO
      if( na.eq.17 )then
         L=13
      elseif( na.eq.29 )then
         L=23
      else
         WRITE(*,*) ' ERROR GENERATED IN EXT530, OSTAPENKO-MALTSEV 1997'
         WRITE(*,*) 'WRONG PARAMETER NA'
         STOP
      ENDif
      DO I=L,NA
         A(I)=A(I)*SN
      ENDDO
      RETURN
      END
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c
c
      SUBROUTINE T01_01 (IOPT,PARMOD,PS,X,Y,Z,BX,BY,BZ)
c
c     RELEASE DATE OF THIS VERSION:   AUGUST 8, 2001.
C
c--------------------------------------------------------------------
C   A DATA-BASED MODEL OF THE EXTERNAL (I.E., WITHOUT EARTH'S CONTRIBUTION) PART OF THE
C   MAGNETOSPHERIC MAGNETIC FIELD, CALIBRATED BY
C    (1) SOLAR WIND PRESSURE PDYN (NANOPASCALS),
C    (2) DST (NANOTESLA),
C    (3) BYIMF,
C    (4) BZIMF (NANOTESLA)
C    (5) G1-INDEX
C    (6) G2-INDEX  (SEE TSYGANENKO [2001] FOR AN EXACT DEFINITION OF THESE TWO INDICES)

c   THESE INPUT PARAMETERS SHOULD BE PLACED IN THE FIRST 6 ELEMENTS
c   OF THE ARRAY PARMOD(10).
C
C   THE REST OF THE INPUT VARIABLES ARE: THE GEODIPOLE TILT ANGLE PS (RADIANS),
C     AND   X,Y,Z -  GSM POSITION (RE)
C
c   IOPT  IS JUST A DUMMY INPUT PARAMETER, NECESSARY TO MAKE THIS SUBROUTINE
C   COMPATIBLE WITH THE TRACING SOFTWARE PACKAGE (GEOPACK). IN THIS MODEL
C   IT DOES NOT AFFECT THE OUTPUT FIELD.
c
C*******************************************************************************************
c** ATTENTION:  THE MODEL IS BASED ON DATA TAKEN SUNWARD FROM X=-15Re, AND HENCE BECOMES   *
C**              INVALID AT LARGER TAILWARD DISTANCES !!!                                  *
C*******************************************************************************************
C
c   OUTPUT:  GSM COMPONENTS OF THE EXTERNAL MAGNETIC FIELD (BX,BY,BZ, nanotesla)
C            COMPUTED AS A SUM OF CONTRIBUTIONS FROM PRINCIPAL FIELD SOURCES
C
c  (C) Copr. 2001, Nikolai A. Tsyganenko, USRA, Code 690.2, NASA GSFC
c      Greenbelt, MD 20771, USA
c
C                            REFERENCE:
C
C    N. A. Tsyganenko, A new data-based model of the near magnetosphere magnetic field:
c       1. Mathematical structure.
c       2. Parameterization and fitting to observations.
c
c             (submitted to JGR, July 2001; available online in the PDF format
c              from anonymous ftp-area www-istp.gsfc.nasa.gov,  /pub/kolya/T01)
c
c----------------------------------------------------------------------
c
      IMPLICIT  REAL * 8  (A - H, O - Z)
C
      REAL*8 PARMOD(10),PS,X,Y,Z,BX,BY,BZ
      REAL*8 A(43),PDYN,DST_AST,BYIMF,BZIMF,G1,G2,PSS,XX,YY,ZZ,
     *  BXCF,BYCF,BZCF,BXT1,BYT1,BZT1,BXT2,BYT2,BZT2,
     *  BXSRC,BYSRC,BZSRC,BXPRC,BYPRC,BZPRC, BXR11,BYR11,BZR11,
     *  BXR12,BYR12,BZR12,BXR21,BYR21,BZR21,BXR22,BYR22,BZR22,HXIMF,
     *  HYIMF,HZIMF,BBX,BBY,BBZ
C
      DATA A /1.00000,2.48341,.58315,.31917,-.08796,-1.17266,3.57478,
     * -.06143,-.01113,.70924,-.01675,-.46056,-.87754,-.03025,.18933,
     * .28089,.16636,-.02932,.02592,-.23537,-.07659,.09117,-.02492,
     * .06816,.55417,.68918,-.04604,2.33521,3.90147,1.28978,.03139,
     * .98751,.21824,41.60182,1.12761,.01376,1.02751,.02969,.15790,
     * 8.94335,28.31280,1.24364,.38013/
C
      IF (X.LT.-20.) THEN
      PRINT *,
     * '  ATTENTION:  THE MODEL IS VALID SUNWARD FROM X=-15 Re ONLY,'
      PRINT *,'              WHILE YOU ARE TRYING TO USE IT AT X=', X
      ENDIF
C
      PDYN=PARMOD(1)
      DST_AST=PARMOD(2)*0.8-13.*SQRT(PDYN)
      BYIMF=PARMOD(3)
      BZIMF=PARMOD(4)
C
      G1=PARMOD(5)
      G2=PARMOD(6)
      PSS=PS
      XX=X
      YY=Y
      ZZ=Z
C
      CALL EXTALL (0,0,0,0,A,43,PDYN,DST_AST,BYIMF,BZIMF,G1,G2,
     *  PSS,XX,YY,ZZ,BXCF,BYCF,BZCF,BXT1,BYT1,BZT1,BXT2,BYT2,BZT2,
     *  BXSRC,BYSRC,BZSRC,BXPRC,BYPRC,BZPRC, BXR11,BYR11,BZR11,
     *  BXR12,BYR12,BZR12,BXR21,BYR21,BZR21,BXR22,BYR22,BZR22,HXIMF,
     *  HYIMF,HZIMF,BBX,BBY,BBZ)
C
      BX=BBX
      BY=BBY
      BZ=BBZ
C
      RETURN
      END
C
C================================================================
      SUBROUTINE EXTALL (IOPGEN,IOPT,IOPB,IOPR,A,NTOT,
     *  PDYN,DST,BYIMF,BZIMF,VBIMF1,VBIMF2,PS,X,Y,Z,
     *  BXCF,BYCF,BZCF,BXT1,BYT1,BZT1,BXT2,BYT2,BZT2,
     *  BXSRC,BYSRC,BZSRC,BXPRC,BYPRC,BZPRC, BXR11,BYR11,BZR11,
     *  BXR12,BYR12,BZR12,BXR21,BYR21,BZR21,BXR22,BYR22,BZR22,HXIMF,
     *  HYIMF,HZIMF,BX,BY,BZ)
C
C   IOPGEN - GENERAL OPTION FLAG:  IOPGEN=0 - CALCULATE TOTAL FIELD
C                                  IOPGEN=1 - DIPOLE SHIELDING ONLY
C                                  IOPGEN=2 - TAIL FIELD ONLY
C                                  IOPGEN=3 - BIRKELAND FIELD ONLY
C                                  IOPGEN=4 - RING CURRENT FIELD ONLY
C                                  IOPGEN=5 - INTERCONNECTION FIELD ONLY
C
C   IOPT -  TAIL FIELD FLAG:       IOPT=0  -  BOTH MODES
C                                  IOPT=1  -  MODE 1 ONLY
C                                  IOPT=2  -  MODE 2 ONLY
C
C   IOPB -  BIRKELAND FIELD FLAG:  IOPB=0  -  ALL 4 TERMS
C                                  IOPB=1  -  REGION 1, MODES 1 AND 2
C                                  IOPB=2  -  REGION 2, MODES 1 AND 2
C
C   IOPR -  RING CURRENT FLAG:     IOPR=0  -  BOTH SRC AND PRC
C                                  IOPR=1  -  SRC ONLY
C                                  IOPR=2  -  PRC ONLY
C
      IMPLICIT  REAL * 8  (A - H, O - Z)
C
      DIMENSION A(NTOT)
C
      COMMON/TAIL/ DXSHIFT1,DXSHIFT2,D,DELTADY  ! THE COMMON BLOCKS FORWARD NONLINEAR PARAMETERS
      COMMON/BIRKPAR/ XKAPPA1,XKAPPA2
      COMMON/RCPAR/ SC_SY,SC_AS,PHI
      COMMON/TSY01_G/ G
      COMMON/TSY01_RH0/ RH0
C
      DATA A0_A,A0_S0,A0_X0 /34.586D0,1.1960D0,3.4397D0/   !   SHUE ET AL. PARAMETERS
      DATA DSIG /0.003D0/, RH0,RH2 /8.0D0,-5.2D0/
c
      XAPPA=(PDYN/2.)**A(39)   !  NOW THIS IS A VARIABLE PARAMETER
      RH0=A(40)
      G=A(41)

      XAPPA3=XAPPA**3

      XX=X*XAPPA
      YY=Y*XAPPA
      ZZ=Z*XAPPA
C
      SPS=DSIN(PS)
c
      X0=A0_X0/XAPPA
      AM=A0_A/XAPPA
      S0=A0_S0
c
      BPERP=DSQRT(BYIMF**2+BZIMF**2)
C
C   CALCULATE THE IMF CLOCK ANGLE:
C
        IF (BYIMF.EQ.0.D0.AND.BZIMF.EQ.0.D0) THEN
            THETA=0.D0
         ELSE
            THETA=DATAN2(BYIMF,BZIMF)
            IF (THETA.LE.0.D0) THETA=THETA+6.283185307D0
        ENDIF
c
       CT=COS(THETA)
       ST=SIN(THETA)
       YS=Y*CT-Z*ST
       ZS=Z*CT+Y*ST

       STHETAH=SIN(THETA/2.)**2
C
C  CALCULATE "IMF" COMPONENTS OUTSIDE THE MAGNETOPAUSE LAYER (HENCE BEGIN WITH "O")
C  THEY ARE NEEDED ONLY IF THE POINT (X,Y,Z) IS WITHIN THE TRANSITION MAGNETOPAUSE LAYER
C  OR OUTSIDE THE MAGNETOSPHERE:
C
      FACTIMF=A(24)+A(25)*STHETAH
c
      OIMFX=0.D0
      OIMFY=BYIMF*FACTIMF
      OIMFZ=BZIMF*FACTIMF
c
      R=SQRT(X**2+Y**2+Z**2)
      XSS=X
      ZSS=Z

  1   XSOLD=XSS      !   BEGIN ITERATIVE SEARCH OF UNWARPED COORDS (TO FIND SIGMA)
      ZSOLD=ZSS

      RH=RH0+RH2*(ZSS/R)**2
      SINPSAS=SPS/(1.D0+(R/RH)**3)**0.33333333D0
      COSPSAS=DSQRT(1.D0-SINPSAS**2)
      ZSS=X*SINPSAS+Z*COSPSAS
      XSS=X*COSPSAS-Z*SINPSAS
      DD=DABS(XSS-XSOLD)+DABS(ZSS-ZSOLD)
      IF (DD.GT.1.D-6) GOTO 1
C                                END OF ITERATIVE SEARCH
      RHO2=Y**2+ZSS**2
      ASQ=AM**2
      XMXM=AM+XSS-X0
      IF (XMXM.LT.0.) XMXM=0. ! THE BOUNDARY IS A CYLINDER TAILWARD OF X=X0-AM
      AXX0=XMXM**2
      ARO=ASQ+RHO2
      SIGMA=DSQRT((ARO+AXX0+SQRT((ARO+AXX0)**2-4.*ASQ*AXX0))/(2.*ASQ))
C
C   NOW, THERE ARE THREE POSSIBLE CASES:
C    (1) INSIDE THE MAGNETOSPHERE   (SIGMA
C    (2) IN THE BOUNDARY LAYER
C    (3) OUTSIDE THE MAGNETOSPHERE AND B.LAYER
C       FIRST OF ALL, CONSIDER THE CASES (1) AND (2):
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF (SIGMA.LT.S0+DSIG) THEN  !  CASES (1) OR (2); CALCULATE THE MODEL FIELD
C                              (WITH THE POTENTIAL "PENETRATED" INTERCONNECTION FIELD):
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      IF (IOPGEN.LE.1) THEN
         CALL SHLCAR3X3(XX,YY,ZZ,PS,CFX,CFY,CFZ)         !  DIPOLE SHIELDING FIELD
         BXCF=CFX*XAPPA3
         BYCF=CFY*XAPPA3
         BZCF=CFZ*XAPPA3
      ELSE
         BXCF=0.D0
         BYCF=0.D0
         BZCF=0.D0
      ENDIF

      IF (IOPGEN.EQ.0.OR.IOPGEN.EQ.2) THEN
         DXSHIFT1=A(26)+A(27)*VBIMF2
         DXSHIFT2=0.D0
         D=A(28)
         DELTADY=A(29)
         CALL DEFORMED (IOPT,PS,XX,YY,ZZ,                !  TAIL FIELD (THREE MODES)
     *    BXT1,BYT1,BZT1,BXT2,BYT2,BZT2)
      ELSE
         BXT1=0.D0
         BYT1=0.D0
         BZT1=0.D0
         BXT2=0.D0
         BYT2=0.D0
         BZT2=0.D0
      ENDIF

      IF (IOPGEN.EQ.0.OR.IOPGEN.EQ.3) THEN
         XKAPPA1=A(35)+A(36)*VBIMF2
         XKAPPA2=A(37)+A(38)*VBIMF2
         CALL BIRK_TOT_01 (IOPB,PS,XX,YY,ZZ,BXR11,BYR11,BZR11,BXR12,
     *   BYR12,BZR12,BXR21,BYR21,BZR21,BXR22,BYR22,BZR22)    !   BIRKELAND FIELD (TWO MODES FOR R1 AND TWO MODES FOR R2)
      ELSE
         BXR11=0.D0
         BYR11=0.D0
         BZR11=0.D0
         BXR12=0.D0
         BYR12=0.D0
         BZR12=0.D0
         BXR21=0.D0
         BYR21=0.D0
         BZR21=0.D0
         BXR22=0.D0
         BYR22=0.D0
         BZR22=0.D0
      ENDIF

      IF (IOPGEN.EQ.0.OR.IOPGEN.EQ.4) THEN
         PHI=1.5707963D0*DTANH(DABS(DST)/A(34))
         ZNAM=DABS(DST)
         IF (ZNAM.LT.20.D0) ZNAM=20.D0
         SC_SY=A(30)*(20.D0/ZNAM)**A(31) *XAPPA    !
         SC_AS=A(32)*(20.D0/ZNAM)**A(33) *XAPPA
         CALL FULL_RC(IOPR,PS,XX,YY,ZZ,BXSRC,BYSRC,BZSRC,BXPRC,BYPRC,
     *                                        BZPRC)  !  SHIELDED RING CURRENT (SRC AND PRC)
      ELSE
         BXSRC=0.D0
         BYSRC=0.D0
         BZSRC=0.D0
         BXPRC=0.D0
         BYPRC=0.D0
         BZPRC=0.D0
      ENDIF
C
      IF (IOPGEN.EQ.0.OR.IOPGEN.EQ.5) THEN
         HXIMF=0.D0
         HYIMF=BYIMF
         HZIMF=BZIMF   ! THESE ARE COMPONENTS OF THE PENETRATED FIELD PER UNIT OF THE PENETRATION COEFFICIENT.
C                        IN OTHER WORDS, THESE ARE DERIVATIVES OF THE PENETRATION FIELD COMPONENTS WITH RESPECT
C                        TO THE PENETRATION COEFFICIENT.   WE ASSUME THAT ONLY THE TRANSVERSE COMPONENT OF THE
C                        FIELD PENETRATES INSIDE.
       ELSE
         HXIMF=0.D0
         HYIMF=0.D0
         HZIMF=0.D0
       ENDIF
C
C-----------------------------------------------------------
C
C    NOW, ADD UP ALL THE COMPONENTS:

      DLP1=(PDYN/2.D0)**A(42)
      DLP2=(PDYN/2.D0)**A(43)

      TAMP1=A(2)+A(3)*DLP1+A(4)*VBIMF1+A(5)*DST
      TAMP2=A(6)+A(7)*DLP2+A(8)*VBIMF1+A(9)*DST
      A_SRC=A(10)+A(11)*DST+A(12)*DSQRT(PDYN)
      A_PRC=A(13)+A(14)*DST+A(15)*DSQRT(PDYN)
      A_R11=A(16)+A(17)*VBIMF2
      A_R12=A(18)+A(19)*VBIMF2
      A_R21=A(20)+A(21)*VBIMF2
      A_R22=A(22)+A(23)*VBIMF2

      BBX=A(1)*BXCF+TAMP1*BXT1+TAMP2*BXT2+A_SRC*BXSRC+A_PRC*BXPRC
     * +A_R11*BXR11+A_R12*BXR12+A_R21*BXR21+A_R22*BXR22
     *   +A(24)*HXIMF+A(25)*HXIMF*STHETAH

      BBY=A(1)*BYCF+TAMP1*BYT1+TAMP2*BYT2+A_SRC*BYSRC+A_PRC*BYPRC
     * +A_R11*BYR11+A_R12*BYR12+A_R21*BYR21+A_R22*BYR22
     *   +A(24)*HYIMF+A(25)*HYIMF*STHETAH

      BBZ=A(1)*BZCF+TAMP1*BZT1+TAMP2*BZT2+A_SRC*BZSRC+A_PRC*BZPRC
     * +A_R11*BZR11+A_R12*BZR12+A_R21*BZR21+A_R22*BZR22
     *   +A(24)*HZIMF+A(25)*HZIMF*STHETAH

C
C   AND WE HAVE THE TOTAL EXTERNAL FIELD.
C
C  NOW, LET US CHECK WHETHER WE HAVE THE CASE (1). IF YES - WE ARE DONE:
C
      IF (SIGMA.LT.S0-DSIG) THEN    !  (X,Y,Z) IS INSIDE THE MAGNETOSPHERE
C-------------------------------------------------------------------------
       BX=BBX
       BY=BBY
       BZ=BBZ
C-------------------------------------------------------------------------
                     ELSE           !  THIS IS THE MOST COMPLEX CASE: WE ARE INSIDE
C                                             THE INTERPOLATION REGION
       FINT=0.5*(1.-(SIGMA-S0)/DSIG)
       FEXT=0.5*(1.+(SIGMA-S0)/DSIG)
C
       CALL DIPOLE_01 (PS,X,Y,Z,QX,QY,QZ)
       BX=(BBX+QX)*FINT+OIMFX*FEXT -QX
       BY=(BBY+QY)*FINT+OIMFY*FEXT -QY
       BZ=(BBZ+QZ)*FINT+OIMFZ*FEXT -QZ
c
        ENDIF  !   THE CASES (1) AND (2) ARE EXHAUSTED; THE ONLY REMAINING
C                      POSSIBILITY IS NOW THE CASE (3):
C--------------------------------------------------------------------------
        ELSE
                CALL DIPOLE_01 (PS,X,Y,Z,QX,QY,QZ)
                BX=OIMFX-QX
                BY=OIMFY-QY
                BZ=OIMFZ-QZ
        ENDIF
C
      END
c
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
         SUBROUTINE  SHLCAR3X3(X,Y,Z,PS,BX,BY,BZ)
C
C   THIS S/R RETURNS THE SHIELDING FIELD FOR THE EARTH'S DIPOLE,
C   REPRESENTED BY  2x3x3=18 "CARTESIAN" HARMONICS, tilted with respect
C   to the z=0 plane  (NB#4, p.74)
C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  The 36 coefficients enter in pairs in the amplitudes of the "cartesian"
c    harmonics (A(1)-A(36).
c  The 14 nonlinear parameters (A(37)-A(50) are the scales Pi,Ri,Qi,and Si
C   entering the arguments of exponents, sines, and cosines in each of the
C   18 "Cartesian" harmonics  PLUS TWO TILT ANGLES FOR THE CARTESIAN HARMONICS
C       (ONE FOR THE PSI=0 MODE AND ANOTHER FOR THE PSI=90 MODE)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
      IMPLICIT  REAL * 8  (A - H, O - Z)
C
      DIMENSION A(50)
      DATA A/-901.2327248,895.8011176,817.6208321,-845.5880889,
     *-83.73539535,86.58542841,336.8781402,-329.3619944,-311.2947120,
     *308.6011161,31.94469304,-31.30824526,125.8739681,-372.3384278,
     *-235.4720434,286.7594095,21.86305585,-27.42344605,-150.4874688,
     *2.669338538,1.395023949,-.5540427503,-56.85224007,3.681827033,
     *-43.48705106,5.103131905,1.073551279,-.6673083508,12.21404266,
     *4.177465543,5.799964188,-.3977802319,-1.044652977,.5703560010,
     *3.536082962,-3.222069852,9.620648151,6.082014949,27.75216226,
     *12.44199571,5.122226936,6.982039615,20.12149582,6.150973118,
     *4.663639687,15.73319647,2.303504968,5.840511214,.8385953499E-01,
     *.3477844929/
C
       P1=A(37)
       P2=A(38)
       P3=A(39)
       R1=A(40)
       R2=A(41)
       R3=A(42)
       Q1=A(43)
       Q2=A(44)
       Q3=A(45)
       S1=A(46)
       S2=A(47)
       S3=A(48)

       T1 =A(49)
       T2 =A(50)
C
       CPS=DCOS(PS)
       SPS=DSIN(PS)
       S2PS=2.D0*CPS      !   MODIFIED HERE (SIN(2*PS) INSTEAD OF SIN(3*PS))
C
       ST1=DSIN(PS*T1)
       CT1=DCOS(PS*T1)
       ST2=DSIN(PS*T2)
       CT2=DCOS(PS*T2)

       X1=X*CT1-Z*ST1
       Z1=X*ST1+Z*CT1
       X2=X*CT2-Z*ST2
       Z2=X*ST2+Z*CT2
C
C
c  MAKE THE TERMS IN THE 1ST SUM ("PERPENDICULAR" SYMMETRY):
C
C       I=1:
C
       SQPR= DSQRT(1.D0/P1**2+1.D0/R1**2)
       CYP = DCOS(Y/P1)
       SYP = DSIN(Y/P1)
       CZR = DCOS(Z1/R1)
       SZR = DSIN(Z1/R1)
       EXPR= DEXP(SQPR*X1)
       FX1 =-SQPR*EXPR*CYP*SZR
       HY1 = EXPR/P1*SYP*SZR
       FZ1 =-EXPR*CYP/R1*CZR
       HX1 = FX1*CT1+FZ1*ST1
       HZ1 =-FX1*ST1+FZ1*CT1

       SQPR= DSQRT(1.D0/P1**2+1.D0/R2**2)
       CYP = DCOS(Y/P1)
       SYP = DSIN(Y/P1)
       CZR = DCOS(Z1/R2)
       SZR = DSIN(Z1/R2)
       EXPR= DEXP(SQPR*X1)
       FX2 =-SQPR*EXPR*CYP*SZR
       HY2 = EXPR/P1*SYP*SZR
       FZ2 =-EXPR*CYP/R2*CZR
       HX2 = FX2*CT1+FZ2*ST1
       HZ2 =-FX2*ST1+FZ2*CT1

       SQPR= DSQRT(1.D0/P1**2+1.D0/R3**2)
       CYP = DCOS(Y/P1)
       SYP = DSIN(Y/P1)
       CZR = DCOS(Z1/R3)
       SZR = DSIN(Z1/R3)
       EXPR= DEXP(SQPR*X1)
       FX3 =-EXPR*CYP*(SQPR*Z1*CZR+SZR/R3*(X1+1.D0/SQPR))
       HY3 = EXPR/P1*SYP*(Z1*CZR+X1/R3*SZR/SQPR)
       FZ3 =-EXPR*CYP*(CZR*(1.D0+X1/R3**2/SQPR)-Z1/R3*SZR)
       HX3 = FX3*CT1+FZ3*ST1
       HZ3 =-FX3*ST1+FZ3*CT1
C
C       I=2:
C
       SQPR= DSQRT(1.D0/P2**2+1.D0/R1**2)
       CYP = DCOS(Y/P2)
       SYP = DSIN(Y/P2)
       CZR = DCOS(Z1/R1)
       SZR = DSIN(Z1/R1)
       EXPR= DEXP(SQPR*X1)
       FX4 =-SQPR*EXPR*CYP*SZR
       HY4 = EXPR/P2*SYP*SZR
       FZ4 =-EXPR*CYP/R1*CZR
       HX4 = FX4*CT1+FZ4*ST1
       HZ4 =-FX4*ST1+FZ4*CT1

       SQPR= DSQRT(1.D0/P2**2+1.D0/R2**2)
       CYP = DCOS(Y/P2)
       SYP = DSIN(Y/P2)
       CZR = DCOS(Z1/R2)
       SZR = DSIN(Z1/R2)
       EXPR= DEXP(SQPR*X1)
       FX5 =-SQPR*EXPR*CYP*SZR
       HY5 = EXPR/P2*SYP*SZR
       FZ5 =-EXPR*CYP/R2*CZR
       HX5 = FX5*CT1+FZ5*ST1
       HZ5 =-FX5*ST1+FZ5*CT1

       SQPR= DSQRT(1.D0/P2**2+1.D0/R3**2)
       CYP = DCOS(Y/P2)
       SYP = DSIN(Y/P2)
       CZR = DCOS(Z1/R3)
       SZR = DSIN(Z1/R3)
       EXPR= DEXP(SQPR*X1)
       FX6 =-EXPR*CYP*(SQPR*Z1*CZR+SZR/R3*(X1+1.D0/SQPR))
       HY6 = EXPR/P2*SYP*(Z1*CZR+X1/R3*SZR/SQPR)
       FZ6 =-EXPR*CYP*(CZR*(1.D0+X1/R3**2/SQPR)-Z1/R3*SZR)
       HX6 = FX6*CT1+FZ6*ST1
       HZ6 =-FX6*ST1+FZ6*CT1
C
C      I=3:
C
       SQPR= DSQRT(1.D0/P3**2+1.D0/R1**2)
       CYP = DCOS(Y/P3)
       SYP = DSIN(Y/P3)
       CZR = DCOS(Z1/R1)
       SZR = DSIN(Z1/R1)
       EXPR= DEXP(SQPR*X1)
       FX7 =-SQPR*EXPR*CYP*SZR
       HY7 = EXPR/P3*SYP*SZR
       FZ7 =-EXPR*CYP/R1*CZR
       HX7 = FX7*CT1+FZ7*ST1
       HZ7 =-FX7*ST1+FZ7*CT1

       SQPR= DSQRT(1.D0/P3**2+1.D0/R2**2)
       CYP = DCOS(Y/P3)
       SYP = DSIN(Y/P3)
       CZR = DCOS(Z1/R2)
       SZR = DSIN(Z1/R2)
       EXPR= DEXP(SQPR*X1)
       FX8 =-SQPR*EXPR*CYP*SZR
       HY8 = EXPR/P3*SYP*SZR
       FZ8 =-EXPR*CYP/R2*CZR
       HX8 = FX8*CT1+FZ8*ST1
       HZ8 =-FX8*ST1+FZ8*CT1

       SQPR= DSQRT(1.D0/P3**2+1.D0/R3**2)
       CYP = DCOS(Y/P3)
       SYP = DSIN(Y/P3)
       CZR = DCOS(Z1/R3)
       SZR = DSIN(Z1/R3)
       EXPR= DEXP(SQPR*X1)
       FX9 =-EXPR*CYP*(SQPR*Z1*CZR+SZR/R3*(X1+1.D0/SQPR))
       HY9 = EXPR/P3*SYP*(Z1*CZR+X1/R3*SZR/SQPR)
       FZ9 =-EXPR*CYP*(CZR*(1.D0+X1/R3**2/SQPR)-Z1/R3*SZR)
       HX9 = FX9*CT1+FZ9*ST1
       HZ9 =-FX9*ST1+FZ9*CT1
C
       A1=A(1)+A(2)*CPS
       A2=A(3)+A(4)*CPS
       A3=A(5)+A(6)*CPS
       A4=A(7)+A(8)*CPS
       A5=A(9)+A(10)*CPS
       A6=A(11)+A(12)*CPS
       A7=A(13)+A(14)*CPS
       A8=A(15)+A(16)*CPS
       A9=A(17)+A(18)*CPS
       BX=A1*HX1+A2*HX2+A3*HX3+A4*HX4+A5*HX5+A6*HX6+A7*HX7+A8*HX8+A9*HX9
       BY=A1*HY1+A2*HY2+A3*HY3+A4*HY4+A5*HY5+A6*HY6+A7*HY7+A8*HY8+A9*HY9
       BZ=A1*HZ1+A2*HZ2+A3*HZ3+A4*HZ4+A5*HZ5+A6*HZ6+A7*HZ7+A8*HZ8+A9*HZ9


c  MAKE THE TERMS IN THE 2ND SUM ("PARALLEL" SYMMETRY):
C
C       I=1
C
       SQQS= DSQRT(1.D0/Q1**2+1.D0/S1**2)
       CYQ = DCOS(Y/Q1)
       SYQ = DSIN(Y/Q1)
       CZS = DCOS(Z2/S1)
       SZS = DSIN(Z2/S1)
       EXQS= DEXP(SQQS*X2)
       FX1 =-SQQS*EXQS*CYQ*CZS *SPS
       HY1 = EXQS/Q1*SYQ*CZS   *SPS
       FZ1 = EXQS*CYQ/S1*SZS   *SPS
       HX1 = FX1*CT2+FZ1*ST2
       HZ1 =-FX1*ST2+FZ1*CT2

       SQQS= DSQRT(1.D0/Q1**2+1.D0/S2**2)
       CYQ = DCOS(Y/Q1)
       SYQ = DSIN(Y/Q1)
       CZS = DCOS(Z2/S2)
       SZS = DSIN(Z2/S2)
       EXQS= DEXP(SQQS*X2)
       FX2 =-SQQS*EXQS*CYQ*CZS *SPS
       HY2 = EXQS/Q1*SYQ*CZS   *SPS
       FZ2 = EXQS*CYQ/S2*SZS   *SPS
       HX2 = FX2*CT2+FZ2*ST2
       HZ2 =-FX2*ST2+FZ2*CT2

       SQQS= DSQRT(1.D0/Q1**2+1.D0/S3**2)
       CYQ = DCOS(Y/Q1)
       SYQ = DSIN(Y/Q1)
       CZS = DCOS(Z2/S3)
       SZS = DSIN(Z2/S3)
       EXQS= DEXP(SQQS*X2)
       FX3 =-SQQS*EXQS*CYQ*CZS *SPS
       HY3 = EXQS/Q1*SYQ*CZS   *SPS
       FZ3 = EXQS*CYQ/S3*SZS   *SPS
       HX3 = FX3*CT2+FZ3*ST2
       HZ3 =-FX3*ST2+FZ3*CT2
C
C       I=2:
C
       SQQS= DSQRT(1.D0/Q2**2+1.D0/S1**2)
       CYQ = DCOS(Y/Q2)
       SYQ = DSIN(Y/Q2)
       CZS = DCOS(Z2/S1)
       SZS = DSIN(Z2/S1)
       EXQS= DEXP(SQQS*X2)
       FX4 =-SQQS*EXQS*CYQ*CZS *SPS
       HY4 = EXQS/Q2*SYQ*CZS   *SPS
       FZ4 = EXQS*CYQ/S1*SZS   *SPS
       HX4 = FX4*CT2+FZ4*ST2
       HZ4 =-FX4*ST2+FZ4*CT2

       SQQS= DSQRT(1.D0/Q2**2+1.D0/S2**2)
       CYQ = DCOS(Y/Q2)
       SYQ = DSIN(Y/Q2)
       CZS = DCOS(Z2/S2)
       SZS = DSIN(Z2/S2)
       EXQS= DEXP(SQQS*X2)
       FX5 =-SQQS*EXQS*CYQ*CZS *SPS
       HY5 = EXQS/Q2*SYQ*CZS   *SPS
       FZ5 = EXQS*CYQ/S2*SZS   *SPS
       HX5 = FX5*CT2+FZ5*ST2
       HZ5 =-FX5*ST2+FZ5*CT2

       SQQS= DSQRT(1.D0/Q2**2+1.D0/S3**2)
       CYQ = DCOS(Y/Q2)
       SYQ = DSIN(Y/Q2)
       CZS = DCOS(Z2/S3)
       SZS = DSIN(Z2/S3)
       EXQS= DEXP(SQQS*X2)
       FX6 =-SQQS*EXQS*CYQ*CZS *SPS
       HY6 = EXQS/Q2*SYQ*CZS   *SPS
       FZ6 = EXQS*CYQ/S3*SZS   *SPS
       HX6 = FX6*CT2+FZ6*ST2
       HZ6 =-FX6*ST2+FZ6*CT2
C
C       I=3:
C
       SQQS= DSQRT(1.D0/Q3**2+1.D0/S1**2)
       CYQ = DCOS(Y/Q3)
       SYQ = DSIN(Y/Q3)
       CZS = DCOS(Z2/S1)
       SZS = DSIN(Z2/S1)
       EXQS= DEXP(SQQS*X2)
       FX7 =-SQQS*EXQS*CYQ*CZS *SPS
       HY7 = EXQS/Q3*SYQ*CZS   *SPS
       FZ7 = EXQS*CYQ/S1*SZS   *SPS
       HX7 = FX7*CT2+FZ7*ST2
       HZ7 =-FX7*ST2+FZ7*CT2

       SQQS= DSQRT(1.D0/Q3**2+1.D0/S2**2)
       CYQ = DCOS(Y/Q3)
       SYQ = DSIN(Y/Q3)
       CZS = DCOS(Z2/S2)
       SZS = DSIN(Z2/S2)
       EXQS= DEXP(SQQS*X2)
       FX8 =-SQQS*EXQS*CYQ*CZS *SPS
       HY8 = EXQS/Q3*SYQ*CZS   *SPS
       FZ8 = EXQS*CYQ/S2*SZS   *SPS
       HX8 = FX8*CT2+FZ8*ST2
       HZ8 =-FX8*ST2+FZ8*CT2

       SQQS= DSQRT(1.D0/Q3**2+1.D0/S3**2)
       CYQ = DCOS(Y/Q3)
       SYQ = DSIN(Y/Q3)
       CZS = DCOS(Z2/S3)
       SZS = DSIN(Z2/S3)
       EXQS= DEXP(SQQS*X2)
       FX9 =-SQQS*EXQS*CYQ*CZS *SPS
       HY9 = EXQS/Q3*SYQ*CZS   *SPS
       FZ9 = EXQS*CYQ/S3*SZS   *SPS
       HX9 = FX9*CT2+FZ9*ST2
       HZ9 =-FX9*ST2+FZ9*CT2

       A1=A(19)+A(20)*S2PS
       A2=A(21)+A(22)*S2PS
       A3=A(23)+A(24)*S2PS
       A4=A(25)+A(26)*S2PS
       A5=A(27)+A(28)*S2PS
       A6=A(29)+A(30)*S2PS
       A7=A(31)+A(32)*S2PS
       A8=A(33)+A(34)*S2PS
       A9=A(35)+A(36)*S2PS

       BX=BX+A1*HX1+A2*HX2+A3*HX3+A4*HX4+A5*HX5+A6*HX6+A7*HX7+A8*HX8
     *   +A9*HX9
       BY=BY+A1*HY1+A2*HY2+A3*HY3+A4*HY4+A5*HY5+A6*HY6+A7*HY7+A8*HY8
     *   +A9*HY9
       BZ=BZ+A1*HZ1+A2*HZ2+A3*HZ3+A4*HZ4+A5*HZ5+A6*HZ6+A7*HZ7+A8*HZ8
     *   +A9*HZ9
C
       RETURN
       END
c
c############################################################################
c
C
      SUBROUTINE DEFORMED (IOPT,PS,X,Y,Z,BX1,BY1,BZ1,BX2,BY2,BZ2)
C
C   IOPT - TAIL FIELD MODE FLAG:   IOPT=0 - THE TWO TAIL MODES ARE ADDED UP
C                                  IOPT=1 - MODE 1 ONLY
C                                  IOPT=2 - MODE 2 ONLY
C
C   CALCULATES GSM COMPONENTS OF TWO UNIT-AMPLITUDE TAIL FIELD MODES,
C    TAKING INTO ACCOUNT BOTH EFFECTS OF DIPOLE TILT:
C    WARPING IN Y-Z (DONE BY THE S/R WARPED) AND BENDING IN X-Z (DONE BY THIS SUBROUTINE)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/TSY01_RH0/ RH0
      DATA RH2,IEPS /-5.2D0,3/
C
C  RH0,RH1,RH2, AND IEPS CONTROL THE TILT-RELATED DEFORMATION OF THE TAIL FIELD
C
      SPS=DSIN(PS)
      CPS=DSQRT(1.D0-SPS**2)
      R2=X**2+Y**2+Z**2
      R=SQRT(R2)
      ZR=Z/R
      RH=RH0+RH2*ZR**2
      DRHDR=-ZR/R*2.D0*RH2*ZR
      DRHDZ= 2.D0*RH2*ZR/R
C
      RRH=R/RH
      F=1.D0/(1.D0+RRH**IEPS)**(1.D0/IEPS)
      DFDR=-RRH**(IEPS-1)*F**(IEPS+1)/RH
      DFDRH=-RRH*DFDR
c
      SPSAS=SPS*F
      CPSAS=DSQRT(1.D0-SPSAS**2)
C
      XAS=X*CPSAS-Z*SPSAS
      ZAS=X*SPSAS+Z*CPSAS
C
      FACPS=SPS/CPSAS*(DFDR+DFDRH*DRHDR)/R
      PSASX=FACPS*X
      PSASY=FACPS*Y
      PSASZ=FACPS*Z+SPS/CPSAS*DFDRH*DRHDZ
C
      DXASDX=CPSAS-ZAS*PSASX
      DXASDY=-ZAS*PSASY
      DXASDZ=-SPSAS-ZAS*PSASZ
      DZASDX=SPSAS+XAS*PSASX
      DZASDY=XAS*PSASY
      DZASDZ=CPSAS+XAS*PSASZ
      FAC1=DXASDZ*DZASDY-DXASDY*DZASDZ
      FAC2=DXASDX*DZASDZ-DXASDZ*DZASDX
      FAC3=DZASDX*DXASDY-DXASDX*DZASDY
C
C     DEFORM:
C
      CALL WARPED(IOPT,PS,XAS,Y,ZAS,BXAS1,BYAS1,BZAS1,BXAS2,BYAS2,BZAS2)
C
      BX1=BXAS1*DZASDZ-BZAS1*DXASDZ +BYAS1*FAC1
      BY1=BYAS1*FAC2
      BZ1=BZAS1*DXASDX-BXAS1*DZASDX +BYAS1*FAC3

      BX2=BXAS2*DZASDZ-BZAS2*DXASDZ +BYAS2*FAC1
      BY2=BYAS2*FAC2
      BZ2=BZAS2*DXASDX-BXAS2*DZASDX +BYAS2*FAC3

      RETURN
      END
C
C------------------------------------------------------------------
C
      SUBROUTINE WARPED (IOPT,PS,X,Y,Z,BX1,BY1,BZ1,BX2,BY2,BZ2)
C
C   CALCULATES GSM COMPONENTS OF THE WARPED FIELD FOR TWO TAIL UNIT MODES.
C   THE WARPING DEFORMATION IS IMPOSED ON THE UNWARPED FIELD, COMPUTED
C   BY THE S/R "UNWARPED".  THE WARPING PARAMETER G WAS OBTAINED BY LEAST
C   SQUARES FITTING TO THE ENTIRE DATASET.
C
C   IOPT - TAIL FIELD MODE FLAG:   IOPT=0 - THE TWO TAIL MODES ARE ADDED UP
C                                  IOPT=1 - MODE 1 ONLY
C                                  IOPT=2 - MODE 2 ONLY
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON/TSY01_G/ G
      DGDX=0.D0
      XL=20.D0
      DXLDX=0.D0

      SPS=DSIN(PS)
      RHO2=Y**2+Z**2
      RHO=DSQRT(RHO2)

      IF (Y.EQ.0.D0.AND.Z.EQ.0.D0) THEN
       PHI=0.D0
       CPHI=1.D0
       SPHI=0.D0
      ELSE
       PHI=DATAN2(Z,Y)
       CPHI=Y/RHO
       SPHI=Z/RHO
      ENDIF

      RR4L4=RHO/(RHO2**2+XL**4)

      F=PHI+G*RHO2*RR4L4*CPHI*SPS
      DFDPHI=1.D0-G*RHO2*RR4L4*SPHI*SPS
      DFDRHO=G*RR4L4**2*(3.D0*XL**4-RHO2**2)*CPHI*SPS
      DFDX=RR4L4*CPHI*SPS*(DGDX*RHO2-G*RHO*RR4L4*4.D0*XL**3*DXLDX)

      CF=DCOS(F)
      SF=DSIN(F)
      YAS=RHO*CF
      ZAS=RHO*SF

      CALL UNWARPED (IOPT,X,YAS,ZAS,BX_AS1,BY_AS1,BZ_AS1,
     *  BX_AS2,BY_AS2,BZ_AS2)

      BRHO_AS =  BY_AS1*CF+BZ_AS1*SF      !   DEFORM THE 1ST MODE
      BPHI_AS = -BY_AS1*SF+BZ_AS1*CF

      BRHO_S = BRHO_AS*DFDPHI
      BPHI_S = BPHI_AS-RHO*(BX_AS1*DFDX+BRHO_AS*DFDRHO)
      BX1    = BX_AS1*DFDPHI

      BY1    = BRHO_S*CPHI-BPHI_S*SPHI
      BZ1    = BRHO_S*SPHI+BPHI_S*CPHI    !   DONE

      BRHO_AS =  BY_AS2*CF+BZ_AS2*SF      !   DEFORM THE 2ND MODE
      BPHI_AS = -BY_AS2*SF+BZ_AS2*CF

      BRHO_S = BRHO_AS*DFDPHI
      BPHI_S = BPHI_AS-RHO*(BX_AS2*DFDX+BRHO_AS*DFDRHO)
      BX2    = BX_AS2*DFDPHI

      BY2    = BRHO_S*CPHI-BPHI_S*SPHI
      BZ2    = BRHO_S*SPHI+BPHI_S*CPHI    !   DONE

      RETURN
      END
C
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
       SUBROUTINE UNWARPED (IOPT,X,Y,Z,BX1,BY1,BZ1,BX2,BY2,BZ2)

C   IOPT - TAIL FIELD MODE FLAG:   IOPT=0 - THE TWO TAIL MODES ARE ADDED UP
C                                  IOPT=1 - MODE 1 ONLY
C                                  IOPT=2 - MODE 2 ONLY
C
C    CALCULATES GSM COMPONENTS OF THE SHIELDED FIELD OF TWO TAIL MODES WITH UNIT
C    AMPLITUDES,  WITHOUT ANY WARPING OR BENDING.  NONLINEAR PARAMETERS OF THE MODES
C    ARE FORWARDED HERE VIA A COMMON BLOCK /TAIL/.
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION A1(60),A2(60)  !   TAIL SHIELDING FIELD PARAMETERS FOR THE MODES #1 & #2

      COMMON/TAIL/ DXSHIFT1,DXSHIFT2,D0,DELTADY
C
      DATA DELTADX1,ALPHA1,XSHIFT1 /1.D0,1.1D0,6.D0/
      DATA DELTADX2,ALPHA2,XSHIFT2 /0.D0,.25D0,4.D0/

      DATA A1/-25.45869857,57.35899080,317.5501869,-2.626756717,
     *-93.38053698,-199.6467926,-858.8129729,34.09192395,845.4214929,
     *-29.07463068,47.10678547,-128.9797943,-781.7512093,6.165038619,
     *167.8905046,492.0680410,1654.724031,-46.77337920,-1635.922669,
     *40.86186772,-.1349775602,-.9661991179E-01,-.1662302354,
     *.002810467517,.2487355077,.1025565237,-14.41750229,-.8185333989,
     *11.07693629,.7569503173,-9.655264745,112.2446542,777.5948964,
     *-5.745008536,-83.03921993,-490.2278695,-1155.004209,39.08023320,
     *1172.780574,-39.44349797,-14.07211198,-40.41201127,-313.2277343,
     *2.203920979,8.232835341,197.7065115,391.2733948,-18.57424451,
     *-437.2779053,23.04976898,11.75673963,13.60497313,4.691927060,
     *18.20923547,27.59044809,6.677425469,1.398283308,2.839005878,
     *31.24817706,24.53577264/

      DATA A2/-287187.1962,4970.499233,410490.1952,-1347.839052,
     *-386370.3240,3317.983750,-143462.3895,5706.513767,171176.2904,
     *250.8882750,-506570.8891,5733.592632,397975.5842,9771.762168,
     *-941834.2436,7990.975260,54313.10318,447.5388060,528046.3449,
     *12751.04453,-21920.98301,-21.05075617,31971.07875,3012.641612,
     *-301822.9103,-3601.107387,1797.577552,-6.315855803,142578.8406,
     *13161.93640,804184.8410,-14168.99698,-851926.6360,-1890.885671,
     *972475.6869,-8571.862853,26432.49197,-2554.752298,-482308.3431,
     *-4391.473324,105155.9160,-1134.622050,-74353.53091,-5382.670711,
     *695055.0788,-916.3365144,-12111.06667,67.20923358,-367200.9285,
     *-21414.14421,14.75567902,20.75638190,59.78601609,16.86431444,
     *32.58482365,23.69472951,17.24977936,13.64902647,68.40989058,
     *11.67828167/

      DATA XM1,XM2/2*-12.D0/

      IF (IOPT.EQ.2) GOTO 1

      XSC1=(X-XSHIFT1-DXSHIFT1)*ALPHA1-XM1*(ALPHA1-1.D0)
      YSC1=Y*ALPHA1
      ZSC1=Z*ALPHA1
      D0SC1=D0*ALPHA1   ! HERE WE USE A SINGLE VALUE D0 OF THE THICKNESS FOR BOTH MODES

      CALL TAILDISK(D0SC1,DELTADX1,DELTADY,XSC1,YSC1,ZSC1,FX1,FY1,FZ1)
      CALL SHLCAR5X5(A1,X,Y,Z,DXSHIFT1,HX1,HY1,HZ1)

      BX1=FX1+HX1
      BY1=FY1+HY1
      BZ1=FZ1+HZ1

      IF (IOPT.EQ.1) THEN
        BX2=0.D0
        BY2=0.D0
        BZ2=0.D0
        RETURN
      ENDIF

 1    XSC2=(X-XSHIFT2-DXSHIFT2)*ALPHA2-XM2*(ALPHA2-1.D0)
      YSC2=Y*ALPHA2
      ZSC2=Z*ALPHA2
      D0SC2=D0*ALPHA2   ! HERE WE USE A SINGLE VALUE D0 OF THE THICKNESS FOR BOTH MODES

      CALL TAILDISK(D0SC2,DELTADX2,DELTADY,XSC2,YSC2,ZSC2,FX2,FY2,FZ2)
      CALL SHLCAR5X5(A2,X,Y,Z,DXSHIFT2,HX2,HY2,HZ2)

      BX2=FX2+HX2
      BY2=FY2+HY2
      BZ2=FZ2+HZ2

      IF (IOPT.EQ.2) THEN
        BX1=0.D0
        BY1=0.D0
        BZ1=0.D0
        RETURN
      ENDIF

      RETURN
      END
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      SUBROUTINE TAILDISK(D0,DELTADX,DELTADY,X,Y,Z,BX,BY,BZ)
c
c       THIS SUBROUTINE COMPUTES THE COMPONENTS OF THE TAIL CURRENT FIELD,
C       SIMILAR TO THAT DESCRIBED BY TSYGANENKO AND PEREDO (1994).  THE
C       DIFFERENCE IS THAT NOW WE USE SPACEWARPING, AS DESCRIBED IN OUR
C       PAPER ON MODELING BIRKELAND CURRENTS (TSYGANENKO AND STERN, 1996)
C       INSTEAD OF SHEARING IT IN THE SPIRIT OF THE T89 TAIL MODEL.
C
      IMPLICIT REAL*8 (A-H,O-Z)
c
      DIMENSION F(5),B(5),C(5)
C
      DATA F /-71.09346626D0,-1014.308601D0,-1272.939359D0,
     *        -3224.935936D0,-44546.86232D0/
      DATA B /10.90101242D0,12.68393898D0,13.51791954D0,14.86775017D0,
     *          15.12306404D0/
      DATA C /.7954069972D0,.6716601849D0,1.174866319D0,2.565249920D0,
     *          10.01986790D0/
C
      RHO=DSQRT(X**2+Y**2)
      DRHODX=X/RHO
      DRHODY=Y/RHO

      DEX=DEXP(X/7.D0)
      D=D0+DELTADY*(Y/20.D0)**2  +DELTADX*DEX !   THE LAST TERM (INTRODUCED 10/11/2000) MAKES THE SHEET
      DDDY=DELTADY*Y*0.005D0                  !   THICKEN SUNWARD, TO AVOID PROBLEMS IN THE SUBSOLAR REGION
      DDDX=DELTADX/7.D0*DEX
C
      DZETA=DSQRT(Z**2+D**2)  !  THIS IS THE SAME SIMPLE WAY TO SPREAD
C                                        OUT THE SHEET, AS THAT USED IN T89
      DDZETADX=D*DDDX/DZETA
      DDZETADY=D*DDDY/DZETA
      DDZETADZ=Z/DZETA

C
      DBX=0.D0
      DBY=0.D0
      DBZ=0.D0
C
      DO 1 I=1,5
C
      BI=B(I)
      CI=C(I)
C
      S1=DSQRT((RHO+BI)**2+(DZETA+CI)**2)
      S2=DSQRT((RHO-BI)**2+(DZETA+CI)**2)

      DS1DRHO=(RHO+BI)/S1
      DS2DRHO=(RHO-BI)/S2
      DS1DDZ=(DZETA+CI)/S1
      DS2DDZ=(DZETA+CI)/S2
C
      DS1DX=DS1DRHO*DRHODX  +DS1DDZ*DDZETADX
      DS1DY=DS1DRHO*DRHODY  +   DS1DDZ*DDZETADY
      DS1DZ=                      DS1DDZ*DDZETADZ
C
      DS2DX=DS2DRHO*DRHODX  +DS2DDZ*DDZETADX
      DS2DY=DS2DRHO*DRHODY  +   DS2DDZ*DDZETADY
      DS2DZ=                    DS2DDZ*DDZETADZ
C
      S1TS2=S1*S2
      S1PS2=S1+S2
      S1PS2SQ=S1PS2**2

      FAC1=DSQRT(S1PS2SQ-(2.D0*BI)**2)
      AS=FAC1/(S1TS2*S1PS2SQ)
      DASDS1=(1.D0/(FAC1*S2)-AS/S1PS2*(S2*S2+S1*(3.D0*S1+4.D0*S2)))
     *          /(S1*S1PS2)
      DASDS2=(1.D0/(FAC1*S1)-AS/S1PS2*(S1*S1+S2*(3.D0*S2+4.D0*S1)))
     *          /(S2*S1PS2)
C
      DASDX=DASDS1*DS1DX+DASDS2*DS2DX
      DASDY=DASDS1*DS1DY+DASDS2*DS2DY
      DASDZ=DASDS1*DS1DZ+DASDS2*DS2DZ
C
      DBX=DBX-F(I)*X*DASDZ
      DBY=DBY-F(I)*Y*DASDZ
  1   DBZ=DBZ+F(I)*(2.D0*AS+X*DASDX+Y*DASDY)

      BX=DBX
      BY=DBY
      BZ=DBZ

      RETURN
      END
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C THIS CODE RETURNS THE SHIELDING FIELD REPRESENTED BY  5x5=25 "CARTESIAN"
C    HARMONICS
C
         SUBROUTINE  SHLCAR5X5(A,X,Y,Z,DSHIFT,HX,HY,HZ)
C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  The NLIN coefficients are the amplitudes of the "cartesian"
c    harmonics (A(1)-A(NLIN).
c  The NNP nonlinear parameters (A(NLIN+1)-A(NTOT) are the scales Pi and Ri
C   entering the arguments of exponents, sines, and cosines in each of the
C   NLIN "Cartesian" harmonics
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
         IMPLICIT  REAL * 8  (A - H, O - Z)
C
         DIMENSION A(60)
C
         DHX=0.D0
         DHY=0.D0
         DHZ=0.D0

         L=0
C
         DO 2 I=1,5
         RP=1.D0/A(50+I)
         CYPI=DCOS(Y*RP)
         SYPI=DSIN(Y*RP)
C
         DO 2 K=1,5
         RR=1.D0/A(55+K)
         SZRK=DSIN(Z*RR)
         CZRK=DCOS(Z*RR)
         SQPR=DSQRT(RP**2+RR**2)
         EPR=DEXP(X*SQPR)
C
         DBX=-SQPR*EPR*CYPI*SZRK
         DBY= RP*EPR*SYPI*SZRK
         DBZ=-RR*EPR*CYPI*CZRK

         L=L+2
         COEF=A(L-1)+A(L)*DSHIFT

         DHX=DHX+COEF*DBX
         DHY=DHY+COEF*DBY
         DHZ=DHZ+COEF*DBZ
c
  2      CONTINUE

         HX=DHX
         HY=DHY
         HZ=DHZ
C
      RETURN
      END
c
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
      SUBROUTINE BIRK_TOT_01 (IOPB,PS,X,Y,Z,BX11,BY11,BZ11,BX12,BY12,
     *                          BZ12,BX21,BY21,BZ21,BX22,BY22,BZ22)
C
C      IOPB -  BIRKELAND FIELD MODE FLAG:
C         IOPB=0 - ALL COMPONENTS
C         IOPB=1 - REGION 1, MODES 1 & 2
C         IOPB=2 - REGION 2, MODES 1 & 2
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION SH11(86),SH12(86),SH21(86),SH22(86)
      COMMON/BIRKPAR/ XKAPPA1,XKAPPA2   !  INPUT PARAMETERS, SPECIFIED FROM S/R EXTALL
      COMMON/DPHI_B_RHO0/ DPHI,B,RHO_0,XKAPPA ! PARAMETERS, CONTROLLING THE DAY-NIGHT ASYMMETRY OF F.A.C.

      DATA SH11/46488.84663,-15541.95244,-23210.09824,-32625.03856,
     *-109894.4551,-71415.32808,58168.94612,55564.87578,-22890.60626,
     *-6056.763968,5091.368100,239.7001538,-13899.49253,4648.016991,
     *6971.310672,9699.351891,32633.34599,21028.48811,-17395.96190,
     *-16461.11037,7447.621471,2528.844345,-1934.094784,-588.3108359,
     *-32588.88216,10894.11453,16238.25044,22925.60557,77251.11274,
     *50375.97787,-40763.78048,-39088.60660,15546.53559,3559.617561,
     *-3187.730438,309.1487975,88.22153914,-243.0721938,-63.63543051,
     *191.1109142,69.94451996,-187.9539415,-49.89923833,104.0902848,
     *-120.2459738,253.5572433,89.25456949,-205.6516252,-44.93654156,
     *124.7026309,32.53005523,-98.85321751,-36.51904756,98.88241690,
     *24.88493459,-55.04058524,61.14493565,-128.4224895,-45.35023460,
     *105.0548704,-43.66748755,119.3284161,31.38442798,-92.87946767,
     *-33.52716686,89.98992001,25.87341323,-48.86305045,59.69362881,
     *-126.5353789,-44.39474251,101.5196856,59.41537992,41.18892281,
     *80.86101200,3.066809418,7.893523804,30.56212082,10.36861082,
     *8.222335945,19.97575641,2.050148531,4.992657093,2.300564232,
     *.2256245602,-.05841594319/

      DATA SH12/210260.4816,-1443587.401,-1468919.281,281939.2993,
     *-1131124.839,729331.7943,2573541.307,304616.7457,468887.5847,
     *181554.7517,-1300722.650,-257012.8601,645888.8041,-2048126.412,
     *-2529093.041,571093.7972,-2115508.353,1122035.951,4489168.802,
     *75234.22743,823905.6909,147926.6121,-2276322.876,-155528.5992,
     *-858076.2979,3474422.388,3986279.931,-834613.9747,3250625.781,
     *-1818680.377,-7040468.986,-414359.6073,-1295117.666,-346320.6487,
     *3565527.409,430091.9496,-.1565573462,7.377619826,.4115646037,
     *-6.146078880,3.808028815,-.5232034932,1.454841807,-12.32274869,
     *-4.466974237,-2.941184626,-.6172620658,12.64613490,1.494922012,
     *-21.35489898,-1.652256960,16.81799898,-1.404079922,-24.09369677,
     *-10.99900839,45.94237820,2.248579894,31.91234041,7.575026816,
     *-45.80833339,-1.507664976,14.60016998,1.348516288,-11.05980247,
     *-5.402866968,31.69094514,12.28261196,-37.55354174,4.155626879,
     *-33.70159657,-8.437907434,36.22672602,145.0262164,70.73187036,
     *85.51110098,21.47490989,24.34554406,31.34405345,4.655207476,
     *5.747889264,7.802304187,1.844169801,4.867254550,2.941393119,
     *.1379899178,.06607020029/

      DATA SH21/162294.6224,503885.1125,-27057.67122,-531450.1339,
     *84747.05678,-237142.1712,84133.61490,259530.0402,69196.05160,
     *-189093.5264,-19278.55134,195724.5034,-263082.6367,-818899.6923,
     *43061.10073,863506.6932,-139707.9428,389984.8850,-135167.5555,
     *-426286.9206,-109504.0387,295258.3531,30415.07087,-305502.9405,
     *100785.3400,315010.9567,-15999.50673,-332052.2548,54964.34639,
     *-152808.3750,51024.67566,166720.0603,40389.67945,-106257.7272,
     *-11126.14442,109876.2047,2.978695024,558.6019011,2.685592939,
     *-338.0004730,-81.99724090,-444.1102659,89.44617716,212.0849592,
     *-32.58562625,-982.7336105,-35.10860935,567.8931751,-1.917212423,
     *-260.2023543,-1.023821735,157.5533477,23.00200055,232.0603673,
     *-36.79100036,-111.9110936,18.05429984,447.0481000,15.10187415,
     *-258.7297813,-1.032340149,-298.6402478,-1.676201415,180.5856487,
     *64.52313024,209.0160857,-53.85574010,-98.52164290,14.35891214,
     *536.7666279,20.09318806,-309.7349530,58.54144539,67.45226850,
     *97.92374406,4.752449760,10.46824379,32.91856110,12.05124381,
     *9.962933904,15.91258637,1.804233877,6.578149088,2.515223491,
     *.1930034238,-.02261109942/

      DATA SH22/-131287.8986,-631927.6885,-318797.4173,616785.8782,
     *-50027.36189,863099.9833,47680.20240,-1053367.944,-501120.3811,
     *-174400.9476,222328.6873,333551.7374,-389338.7841,-1995527.467,
     *-982971.3024,1960434.268,297239.7137,2676525.168,-147113.4775,
     *-3358059.979,-2106979.191,-462827.1322,1017607.960,1039018.475,
     *520266.9296,2627427.473,1301981.763,-2577171.706,-238071.9956,
     *-3539781.111,94628.16420,4411304.724,2598205.733,637504.9351,
     *-1234794.298,-1372562.403,-2.646186796,-31.10055575,2.295799273,
     *19.20203279,30.01931202,-302.1028550,-14.78310655,162.1561899,
     *.4943938056,176.8089129,-.2444921680,-100.6148929,9.172262228,
     *137.4303440,-8.451613443,-84.20684224,-167.3354083,1321.830393,
     *76.89928813,-705.7586223,18.28186732,-770.1665162,-9.084224422,
     *436.3368157,-6.374255638,-107.2730177,6.080451222,65.53843753,
     *143.2872994,-1028.009017,-64.22739330,547.8536586,-20.58928632,
     *597.3893669,10.17964133,-337.7800252,159.3532209,76.34445954,
     *84.74398828,12.76722651,27.63870691,32.69873634,5.145153451,
     *6.310949163,6.996159733,1.971629939,4.436299219,2.904964304,
     *.1486276863,.06859991529/
C
      XKAPPA=XKAPPA1        !  FORWARDED IN BIRK_1N2
      X_SC=XKAPPA1-1.1D0    !  FORWARDED IN BIRK_SHL

      IF (IOPB.EQ.0.OR.IOPB.EQ.1) THEN

      CALL BIRK_1N2_01 (1,1,PS,X,Y,Z,FX11,FY11,FZ11)           !  REGION 1, MODE 1
      CALL BIRK_SHL (SH11,PS,X_SC,X,Y,Z,HX11,HY11,HZ11)
      BX11=FX11+HX11
      BY11=FY11+HY11
      BZ11=FZ11+HZ11

      CALL BIRK_1N2_01 (1,2,PS,X,Y,Z,FX12,FY12,FZ12)           !  REGION 1, MODE 2
      CALL BIRK_SHL (SH12,PS,X_SC,X,Y,Z,HX12,HY12,HZ12)
      BX12=FX12+HX12
      BY12=FY12+HY12
      BZ12=FZ12+HZ12

      ENDIF

      XKAPPA=XKAPPA2        !  FORWARDED IN BIRK_1N2
      X_SC=XKAPPA2-1.0D0    !  FORWARDED IN BIRK_SHL

      IF (IOPB.EQ.0.OR.IOPB.EQ.2) THEN

      CALL BIRK_1N2_01 (2,1,PS,X,Y,Z,FX21,FY21,FZ21)           !  REGION 2, MODE 1
      CALL BIRK_SHL (SH21,PS,X_SC,X,Y,Z,HX21,HY21,HZ21)
      BX21=FX21+HX21
      BY21=FY21+HY21
      BZ21=FZ21+HZ21

      CALL BIRK_1N2_01 (2,2,PS,X,Y,Z,FX22,FY22,FZ22)           !  REGION 2, MODE 2
      CALL BIRK_SHL (SH22,PS,X_SC,X,Y,Z,HX22,HY22,HZ22)
      BX22=FX22+HX22
      BY22=FY22+HY22
      BZ22=FZ22+HZ22

      ENDIF

      RETURN
      END
C
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c
      SUBROUTINE BIRK_1N2_01 (NUMB,MODE,PS,X,Y,Z,BX,BY,BZ)        !   NB# 6, P.60
C
C  CALCULATES COMPONENTS  OF REGION 1/2 FIELD IN SPHERICAL COORDS.  DERIVED FROM THE S/R DIPDEF2C (WHICH
C    DOES THE SAME JOB, BUT INPUT/OUTPUT THERE WAS IN SPHERICAL COORDS, WHILE HERE WE USE CARTESIAN ONES)
C
C   INPUT:  NUMB=1 (2) FOR REGION 1 (2) CURRENTS
C           MODE=1 YIELDS SIMPLE SINUSOIDAL MLT VARIATION, WITH MAXIMUM CURRENT AT DAWN/DUSK MERIDIAN
C     WHILE MODE=2 YIELDS THE SECOND HARMONIC.
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A11(31),A12(31),A21(31),A22(31)
      COMMON/MODENUM/ M
      COMMON/DTHETA/ DTHETA

      COMMON/DPHI_B_RHO0/ DPHI,B,RHO_0,XKAPPA ! THESE PARAMETERS CONTROL DAY-NIGHT ASYMMETRY OF F.A.C., AS FOLLOWS:

C  (1) DPHI:   HALF-DIFFERENCE (IN RADIANS) BETWEEN DAY AND NIGHT LATITUDE OF FAC OVAL AT IONOSPHERIC ALTITUDE;
C              TYPICAL VALUE: 0.06
C  (2) B:      AN ASYMMETRY FACTOR AT HIGH-ALTITUDES;  FOR B=0, THE ONLY ASYMMETRY IS THAT FROM DPHI
C              TYPICAL VALUES: 0.35-0.70
C  (3) RHO_0:  A FIXED PARAMETER, DEFINING THE DISTANCE RHO, AT WHICH THE LATITUDE SHIFT GRADUALLY SATURATES AND
C              STOPS INCREASING
C              ITS VALUE WAS ASSUMED FIXED, EQUAL TO 7.0.
C  (4) XKAPPA: AN OVERALL SCALING FACTOR, WHICH CAN BE USED FOR CHANGING THE SIZE OF THE F.A.C. OVAL
C
      DATA BETA,RH,EPS/0.9D0,10.D0,3.D0/ ! parameters of the tilt-dependent deformation of the untilted F.A.C. field

      DATA A11/.1618068350,-.1797957553,2.999642482,-.9322708978,
     *-.6811059760,.2099057262,-8.358815746,-14.86033550,.3838362986,
     *-16.30945494,4.537022847,2.685836007,27.97833029,6.330871059,
     *1.876532361,18.95619213,.9651528100,.4217195118,-.08957770020,
     *-1.823555887,.7457045438,-.5785916524,-1.010200918,.01112389357,
     *.09572927448,-.3599292276,8.713700514,.9763932955,3.834602998,
     *2.492118385,.7113544659/
      DATA A12/.7058026940,-.2845938535,5.715471266,-2.472820880,
     *-.7738802408,.3478293930,-11.37653694,-38.64768867,.6932927651,
     *-212.4017288,4.944204937,3.071270411,33.05882281,7.387533799,
     *2.366769108,79.22572682,.6154290178,.5592050551,-.1796585105,
     *-1.654932210,.7309108776,-.4926292779,-1.130266095,-.009613974555,
     *.1484586169,-.2215347198,7.883592948,.02768251655,2.950280953,
     *1.212634762,.5567714182/
      DATA A21/.1278764024,-.2320034273,1.805623266,-32.37241440,
     *-.9931490648,.3175085630,-2.492465814,-16.21600096,.2695393416,
     *-6.752691265,3.971794901,14.54477563,41.10158386,7.912889730,
     *1.258297372,9.583547721,1.014141963,.5104134759,-.1790430468,
     *-1.756358428,.7561986717,-.6775248254,-.04014016420,.01446794851,
     *.1200521731,-.2203584559,4.508963850,.8221623576,1.779933730,
     *1.102649543,.8867880020/
      DATA A22/.4036015198,-.3302974212,2.827730930,-45.44405830,
     *-1.611103927,.4927112073,-.003258457559,-49.59014949,.3796217108,
     *-233.7884098,4.312666980,18.05051709,28.95320323,11.09948019,
     *.7471649558,67.10246193,.5667096597,.6468519751,-.1560665317,
     *-1.460805289,.7719653528,-.6658988668,.2515179349E-05,
     *.02426021891,.1195003324,-.2625739255,4.377172556,.2421190547,
     *2.503482679,1.071587299,.7247997430/

      B=0.5
      RHO_0=7.0

      M=MODE
      IF (NUMB.EQ.1) THEN
          DPHI=0.055D0
          DTHETA=0.06D0
      ENDIF

      IF (NUMB.EQ.2) THEN
          DPHI=0.030D0
          DTHETA=0.09D0
      ENDIF

      Xsc=X*XKAPPA
      Ysc=Y*XKAPPA
      Zsc=Z*XKAPPA
      RHO=DSQRT(Xsc**2+Zsc**2)

      Rsc=DSQRT(Xsc**2+Ysc**2+Zsc**2)                                 !  SCALED
      RHO2=RHO_0**2

      IF (Xsc.EQ.0.D0.AND.Zsc.EQ.0.D0) THEN
         PHI=0.D0
      ELSE
         PHI=DATAN2(-Zsc,Xsc)  !  FROM CARTESIAN TO CYLINDRICAL (RHO,PHI,Y)
      ENDIF

      SPHIC=DSIN(PHI)
      CPHIC=DCOS(PHI)  !  "C" means "CYLINDRICAL", TO DISTINGUISH FROM SPHERICAL PHI

      BRACK=DPHI+B*RHO2/(RHO2+1.D0)*(RHO**2-1.D0)/(RHO2+RHO**2)
      R1RH=(Rsc-1.D0)/RH
      IF (R1RH.LT.0.D0) R1RH=0.D0  !  AVOID NEGATIVE VALUES OF R1RH, WHICH MAY OCCUR IF THE
C                                     POINT (X,Y,Z) LIES CLOSE TO EARTH'S SURFACE AND THE S.W.
C                                     PRESSURE IS ABNORMALLY LOW

      PSIAS=BETA*PS/(1.D0+R1RH**EPS)**(1.D0/EPS)

      PHIS=PHI-BRACK*DSIN(PHI) -PSIAS
      DPHISPHI=1.D0-BRACK*DCOS(PHI)
      DPHISRHO=-2.D0*B*RHO2*RHO/(RHO2+RHO**2)**2 *DSIN(PHI)
     *   +BETA*PS*R1RH**(EPS-1.D0)*RHO/(RH*Rsc*
     *   (1.D0+R1RH**EPS)**(1.D0/EPS+1.D0))
      DPHISDY= BETA*PS*R1RH**(EPS-1.D0)*Ysc/(RH*Rsc*
     *   (1.D0+R1RH**EPS)**(1.D0/EPS+1.D0))

      SPHICS=DSIN(PHIS)
      CPHICS=DCOS(PHIS)

      XS= RHO*CPHICS
      ZS=-RHO*SPHICS

      IF (NUMB.EQ.1) THEN
        IF (MODE.EQ.1) CALL TWOCONES (A11,XS,Ysc,ZS,BXS,BYAS,BZS)
        IF (MODE.EQ.2) CALL TWOCONES (A12,XS,Ysc,ZS,BXS,BYAS,BZS)
      ELSE
        IF (MODE.EQ.1) CALL TWOCONES (A21,XS,Ysc,ZS,BXS,BYAS,BZS)
        IF (MODE.EQ.2) CALL TWOCONES (A22,XS,Ysc,ZS,BXS,BYAS,BZS)
      ENDIF

      BRHOAS=BXS*CPHICS-BZS*SPHICS
      BPHIAS=-BXS*SPHICS-BZS*CPHICS

      BRHO_S=BRHOAS*DPHISPHI                             *XKAPPA        ! SCALING
      BPHI_S=(BPHIAS-RHO*(BYAS*DPHISDY+BRHOAS*DPHISRHO)) *XKAPPA
      BY_S=BYAS*DPHISPHI                                 *XKAPPA

      BX=BRHO_S*CPHIC-BPHI_S*SPHIC
      BY=BY_S
      BZ=-BRHO_S*SPHIC-BPHI_S*CPHIC

      RETURN
      END
c
C=========================================================================
c
      SUBROUTINE TWOCONES (A,X,Y,Z,BX,BY,BZ)
C
C  ADDS FIELDS FROM TWO CONES (NORTHERN AND SOUTHERN), WITH A PROPER SYMMETRY
C  OF THE CURRENT AND FIELD, CORRESPONDING TO THE REGION 1 BIRKELAND CURRENTS. (NB #6, P.58).
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(31)

      CALL ONE_CONE (A,X,Y,Z,BXN,BYN,BZN)
      CALL ONE_CONE (A,X,-Y,-Z,BXS,BYS,BZS)
      BX=BXN-BXS
      BY=BYN+BYS
      BZ=BZN+BZS

      RETURN
      END
c
C-------------------------------------------------------------------------
C
      SUBROUTINE ONE_CONE(A,X,Y,Z,BX,BY,BZ)
c
c  RETURNS FIELD COMPONENTS FOR A DEFORMED CONICAL CURRENT SYSTEM, FITTED TO A BIOSAVART FIELD
c  HERE ONLY THE NORTHERN CONE IS TAKEN INTO ACCOUNT.
c
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(31)

      COMMON/DTHETA/ DTHETA
      COMMON/MODENUM/ M

      DATA DR,DT/1.D-6,1.D-6/  !   JUST FOR NUMERICAL DIFFERENTIATION

      THETA0=A(31)

      RHO2=X**2+Y**2
      RHO=DSQRT(RHO2)
      R=DSQRT(RHO2+Z**2)
      THETA=DATAN2(RHO,Z)
      PHI=DATAN2(Y,X)
C
C   MAKE THE DEFORMATION OF COORDINATES:
C
       RS=R_S(A,R,THETA)
       THETAS=THETA_S(A,R,THETA)
       PHIS=PHI
C
C   CALCULATE FIELD COMPONENTS AT THE NEW POSITION (ASTERISKED):
C
       CALL FIALCOS (RS,THETAS,PHIS,BTAST,BFAST,M,THETA0,DTHETA)    !   MODE #M
C
C   NOW TRANSFORM B{R,T,F}_AST BY THE DEFORMATION TENSOR:
C
C      FIRST OF ALL, FIND THE DERIVATIVES:
C
       DRSDR=(R_S(A,R+DR,THETA)-R_S(A,R-DR,THETA))/(2.D0*DR)
       DRSDT=(R_S(A,R,THETA+DT)-R_S(A,R,THETA-DT))/(2.D0*DT)
       DTSDR=(THETA_S(A,R+DR,THETA)-THETA_S(A,R-DR,THETA))/(2.D0*DR)
       DTSDT=(THETA_S(A,R,THETA+DT)-THETA_S(A,R,THETA-DT))/(2.D0*DT)

       STSST=DSIN(THETAS)/DSIN(THETA)
       RSR=RS/R

       BR     =-RSR/R*STSST*BTAST*DRSDT                 !   NB#6, P.43    BRAST DOES NOT ENTER HERE
       BTHETA = RSR*STSST*BTAST*DRSDR                  !          (IT IS IDENTICALLY ZERO IN OUR CASE)
       BPHI   = RSR*BFAST*(DRSDR*DTSDT-DRSDT*DTSDR)

       S=RHO/R
       C=Z/R
       SF=Y/RHO
       CF=X/RHO

       BE=BR*S+BTHETA*C

       BX=A(1)*(BE*CF-BPHI*SF)
       BY=A(1)*(BE*SF+BPHI*CF)
       BZ=A(1)*(BR*C-BTHETA*S)

       RETURN
       END
C
C=====================================================================================
      DOUBLE PRECISION FUNCTION R_S(A,R,THETA)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(31)
C
      R_S=R+A(2)/R+A(3)*R/DSQRT(R**2+A(11)**2)+A(4)*R/(R**2+A(12)**2)
     *+(A(5)+A(6)/R+A(7)*R/DSQRT(R**2+A(13)**2)+A(8)*R/(R**2+A(14)**2))*
     * DCOS(THETA)
     *+(A(9)*R/DSQRT(R**2+A(15)**2)+A(10)*R/(R**2+A(16)**2)**2)
     * *DCOS(2.D0*THETA)
C
      RETURN
      END
C
C-----------------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION THETA_S(A,R,THETA)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(31)
c
      THETA_S=THETA+(A(17)+A(18)/R+A(19)/R**2
     *                +A(20)*R/DSQRT(R**2+A(27)**2))*DSIN(THETA)
     * +(A(21)+A(22)*R/DSQRT(R**2+A(28)**2)
     *                +A(23)*R/(R**2+A(29)**2))*DSIN(2.D0*THETA)
     * +(A(24)+A(25)/R+A(26)*R/(R**2+A(30)**2))*DSIN(3.D0*THETA)
C
      RETURN
      END
C
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
      SUBROUTINE FIALCOS(R,THETA,PHI,BTHETA,BPHI,N,THETA0,DT)
C
C  CONICAL MODEL OF BIRKELAND CURRENT FIELD; BASED ON THE OLD S/R FIALCO (OF 1990-91)
C  NB OF 1985-86-88, NOTE OF MARCH 5, BUT HERE BOTH INPUT AND OUTPUT ARE IN SPHERICAL CDS.

C  BTN, AND BPN ARE THE ARRAYS OF BTHETA AND BPHI (BTN(i), BPN(i) CORRESPOND TO i-th MODE).
C   ONLY FIRST  N  MODE AMPLITUDES ARE COMPUTED (N<=10).
C    THETA0 IS THE ANGULAR HALF-WIDTH OF THE CONE, DT IS THE ANGULAR H.-W. OF THE CURRENT LAYER

C   NOTE:  BR=0  (BECAUSE ONLY RADIAL CURRENTS ARE PRESENT IN THIS MODEL)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION  BTN(10),BPN(10),CCOS(10),SSIN(10)

      SINTE=DSIN(THETA)
      RO=R*SINTE
      COSTE=DCOS(THETA)
      SINFI=DSIN(PHI)
      COSFI=DCOS(PHI)
      TG=SINTE/(1.D0+COSTE)   !        TAN(THETA/2)
      CTG=SINTE/(1.D0-COSTE)  !        COT(THETA/2)
C
C
      TETANP=THETA0+DT
      TETANM=THETA0-DT
      IF(THETA.LT.TETANM) GOTO 1
      TGP=DTAN(TETANP*0.5D0)
      TGM=DTAN(TETANM*0.5D0)
      TGM2=TGM*TGM
      TGP2=TGP*TGP
  1   CONTINUE

      COSM1=1.D0
      SINM1=0.D0
      TM=1.D0
      TGM2M=1.D0
      TGP2M=1.D0

      DO 2 M=1,N
      TM=TM*TG
      CCOS(M)=COSM1*COSFI-SINM1*SINFI
      SSIN(M)=SINM1*COSFI+COSM1*SINFI
      COSM1=CCOS(M)
      SINM1=SSIN(M)
      IF(THETA.LT.TETANM) THEN
      T=TM
      DTT=0.5D0*M*TM*(TG+CTG)
      DTT0=0.D0
      ELSE IF(THETA.LT.TETANP) THEN
      TGM2M=TGM2M*TGM2
      FC=1.D0/(TGP-TGM)
      FC1=1.D0/(2*M+1)
      TGM2M1=TGM2M*TGM
      TG21=1.D0+TG*TG
      T=FC*(TM*(TGP-TG)+FC1*(TM*TG-TGM2M1/TM))
      DTT=0.5D0*M*FC*TG21*(TM/TG*(TGP-TG)-FC1*(TM-TGM2M1/(TM*TG)))
      DTT0=0.5D0*FC*((TGP+TGM)*(TM*TG-FC1*(TM*TG-TGM2M1/TM))+
     * TM*(1.D0-TGP*TGM)-(1.D0+TGM2)*TGM2M/TM)
      ELSE
      TGP2M=TGP2M*TGP2
      TGM2M=TGM2M*TGM2
      FC=1.D0/(TGP-TGM)
      FC1=1.D0/(2*M+1)
      T=FC*FC1*(TGP2M*TGP-TGM2M*TGM)/TM
      DTT=-T*M*0.5D0*(TG+CTG)
      ENDIF

      BTN(M)=M*T*CCOS(M)/RO
  2   BPN(M)=-DTT*SSIN(M)/R

      BTHETA=BTN(N) *800.
      BPHI  =BPN(N) *800.

      RETURN
      END
C
C-------------------------------------------------------------------------
C
C
         SUBROUTINE BIRK_SHL (A,PS,X_SC,X,Y,Z,BX,BY,BZ)
C
       IMPLICIT  REAL * 8  (A - H, O - Z)
         DIMENSION A(86)
C
         CPS=DCOS(PS)
         SPS=DSIN(PS)

         S3PS=2.D0*CPS
C
         PST1=PS*A(85)
         PST2=PS*A(86)

         ST1=DSIN(PST1)
         CT1=DCOS(PST1)
         ST2=DSIN(PST2)
         CT2=DCOS(PST2)

         X1=X*CT1-Z*ST1
         Z1=X*ST1+Z*CT1
         X2=X*CT2-Z*ST2
         Z2=X*ST2+Z*CT2
C
         L=0
         GX=0.D0
         GY=0.D0
         GZ=0.D0
C
         DO 1 M=1,2     !    M=1 IS FOR THE 1ST SUM ("PERP." SYMMETRY)
C                          AND M=2 IS FOR THE SECOND SUM ("PARALL." SYMMETRY)
             DO 2 I=1,3
                  P=A(72+I)
                  Q=A(78+I)
                  CYPI=DCOS(Y/P)
                  CYQI=DCOS(Y/Q)
                  SYPI=DSIN(Y/P)
                  SYQI=DSIN(Y/Q)
C
                DO 3 K=1,3
                   R=A(75+K)
                   S=A(81+K)
                   SZRK=DSIN(Z1/R)
                   CZSK=DCOS(Z2/S)
                   CZRK=DCOS(Z1/R)
                   SZSK=DSIN(Z2/S)
                     SQPR=DSQRT(1.D0/P**2+1.D0/R**2)
                     SQQS=DSQRT(1.D0/Q**2+1.D0/S**2)
                        EPR=DEXP(X1*SQPR)
                        EQS=DEXP(X2*SQQS)
C
                  DO 4 N=1,2  ! N=1 IS FOR THE FIRST PART OF EACH COEFFICIENT
C                                AND N=2 IS FOR THE SECOND ONE

                    DO 5 NN=1,2 !   NN = 1,2 FURTHER SPLITS THE COEFFICIENTS INTO 2 PARTS,
C                                         TO TAKE INTO ACCOUNT THE SCALE FACTOR DEPENDENCE

                    IF (M.EQ.1) THEN
                         FX=-SQPR*EPR*CYPI*SZRK
                         FY=EPR*SYPI*SZRK/P
                         FZ=-EPR*CYPI*CZRK/R
                       IF (N.EQ.1) THEN
                         IF (NN.EQ.1) THEN
                          HX=FX
                          HY=FY
                          HZ=FZ
                         ELSE
                          HX=FX*X_SC
                          HY=FY*X_SC
                          HZ=FZ*X_SC
                         ENDIF
                       ELSE
                         IF (NN.EQ.1) THEN
                          HX=FX*CPS
                          HY=FY*CPS
                          HZ=FZ*CPS
                         ELSE
                          HX=FX*CPS*X_SC
                          HY=FY*CPS*X_SC
                          HZ=FZ*CPS*X_SC
                         ENDIF
                       ENDIF

                     ELSE                            !   M.EQ.2
                         FX=-SPS*SQQS*EQS*CYQI*CZSK
                         FY=SPS/Q*EQS*SYQI*CZSK
                         FZ=SPS/S*EQS*CYQI*SZSK
                       IF (N.EQ.1) THEN
                        IF (NN.EQ.1) THEN
                          HX=FX
                          HY=FY
                          HZ=FZ
                        ELSE
                          HX=FX*X_SC
                          HY=FY*X_SC
                          HZ=FZ*X_SC
                        ENDIF
                       ELSE
                        IF (NN.EQ.1) THEN
                         HX=FX*S3PS
                         HY=FY*S3PS
                         HZ=FZ*S3PS
                        ELSE
                         HX=FX*S3PS*X_SC
                         HY=FY*S3PS*X_SC
                         HZ=FZ*S3PS*X_SC
                        ENDIF
                       ENDIF
                  ENDIF
       L=L+1

       IF (M.EQ.1) THEN
       HXR=HX*CT1+HZ*ST1
       HZR=-HX*ST1+HZ*CT1
       ELSE
       HXR=HX*CT2+HZ*ST2
       HZR=-HX*ST2+HZ*CT2
       ENDIF

       GX=GX+HXR*A(L)
       GY=GY+HY *A(L)
  5    GZ=GZ+HZR*A(L)

  4   CONTINUE
  3   CONTINUE
  2   CONTINUE
  1   CONTINUE

      BX=GX
      BY=GY
      BZ=GZ

      RETURN
      END

C
C************************************************************************************
C
      SUBROUTINE FULL_RC (IOPR,PS,X,Y,Z,BXSRC,BYSRC,BZSRC,BXPRC,BYPRC,
     *  BZPRC)
C
C   CALCULATES GSM FIELD COMPONENTS OF THE SYMMETRIC (SRC) AND PARTIAL (PRC) COMPONENTS OF THE RING CURRENT
C   SRC  PROVIDES A DEPRESSION OF -28 nT AT EARTH
C   PRC  CORRESPONDS TO THE PRESSURE DIFFERENCE OF 2 nPa BETWEEN MIDNIGHT AND NOON RING CURRENT
C             PARTICLE PRESSURE AND YIELDS A DEPRESSION OF -17 nT AT X=-6Re
C
C   SC_SY AND SC_PR ARE SCALING FACTORS FOR THE SYMMETRIC AND PARTIAL COMPONENTS:
C          VALUES LARGER THAN 1 RESULT IN SPATIALLY LARGER CURRENTS
C
C   PHI IS THE ROTATION ANGLE IN RADIANS OF THE PARTIAL RING CURRENT (MEASURED FROM MIDNIGHT TOWARD DUSK)
C
C     IOPR -  A RING CURRENT CALCULATION FLAG (FOR LEAST-SQUARES FITTING ONLY):
C             IOPR=0 - BOTH SRC AND PRC FIELDS ARE CALCULATED
C             IOPR=1 - SRC ONLY
C             IOPR=2 - PRC ONLY
C
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION C_SY(86),C_PR(86)
        COMMON /RCPAR/ SC_SY,SC_PR,PHI
C
        DATA C_SY/1675.694858,1780.006388,-961.6082149,-1668.914259,
     *-27.40437029,-107.4169670,27.76189943,92.89740503,-43.92949274,
     *-403.6444072,6.167161865,298.2779761,-1680.779044,-1780.933039,
     *964.1861088,1670.988659,27.48864650,107.7809519,-27.84600972,
     *-93.20691865,44.28496784,404.4537249,-6.281958730,-298.6050952,
     *-7.971914848,2.017383761,-1.492230168,-1.957411655,-.08525523181,
     *-.3811813235,.08446716725,.3215044399,-.7141912767,-.9086294596,
     *.2966677742,-.04736679933,-11.38731325,.1719795189,1.356233066,
     *.8613438429,-.09143823092,-.2593979098,.04244838338,.06318383319,
     *-.5861372726,-.03368780733,-.07104470269,-.06909052953,
     *-60.18659631,-32.87563877,11.76450433,5.891673644,2.562360333,
     *6.215377232,-1.273945165,-1.864704763,-5.394837143,-8.799382627,
     *3.743066561,-.7649164511,57.09210569,32.61236511,-11.28688017,
     *-5.849523392,-2.470635922,-5.961417272,1.230031099,1.793192595,
     *5.383736074,8.369895153,-3.611544412,.7898988697,7.970609948,
     *7.981216562,35.16822497,12.45651654,1.689755359,3.678712366,
     *23.66117284,6.987136092,6.886678677,20.91245928,1.650064156,
     *3.474068566,.3474715765,.6564043111/

        DATA C_PR/-64820.58481,-63965.62048,66267.93413,135049.7504,
     *-36.56316878,124.6614669,56.75637955,-87.56841077,5848.631425,
     *4981.097722,-6233.712207,-10986.40188,68716.52057,65682.69473,
     *-69673.32198,-138829.3568,43.45817708,-117.9565488,-62.14836263,
     *79.83651604,-6211.451069,-5151.633113,6544.481271,11353.03491,
     *23.72352603,-256.4846331,25.77629189,145.2377187,-4.472639098,
     *-3.554312754,2.936973114,2.682302576,2.728979958,26.43396781,
     *-9.312348296,-29.65427726,-247.5855336,-206.9111326,74.25277664,
     *106.4069993,15.45391072,16.35943569,-5.965177750,-6.079451700,
     *115.6748385,-35.27377307,-32.28763497,-32.53122151,93.74409310,
     *84.25677504,-29.23010465,-43.79485175,-6.434679514,-6.620247951,
     *2.443524317,2.266538956,-43.82903825,6.904117876,12.24289401,
     *17.62014361,152.3078796,124.5505289,-44.58690290,-63.02382410,
     *-8.999368955,-9.693774119,3.510930306,3.770949738,-77.96705716,
     *22.07730961,20.46491655,18.67728847,9.451290614,9.313661792,
     *644.7620970,418.2515954,7.183754387,35.62128817,19.43180682,
     *39.57218411,15.69384715,7.123215241,2.300635346,21.90881131,
     *-.01775839370,.3996346710/

        CALL SRC_PRC (IOPR,SC_SY,SC_PR,PHI,PS,X,Y,Z,HXSRC,HYSRC,HZSRC,
     *      HXPRC,HYPRC,HZPRC)

        X_SC=SC_SY-1.D0
        IF (IOPR.EQ.0.OR.IOPR.EQ.1) THEN
          CALL RC_SHIELD (C_SY,PS,X_SC,X,Y,Z,FSX,FSY,FSZ)
        ELSE
           FSX=0.D0
           FSY=0.D0
           FSZ=0.D0
        ENDIF

        X_SC=SC_PR-1.D0
        IF (IOPR.EQ.0.OR.IOPR.EQ.2) THEN
          CALL RC_SHIELD (C_PR,PS,X_SC,X,Y,Z,FPX,FPY,FPZ)
        ELSE
           FPX=0.D0
           FPY=0.D0
           FPZ=0.D0
        ENDIF

        BXSRC=HXSRC+FSX
        BYSRC=HYSRC+FSY
        BZSRC=HZSRC+FSZ

        BXPRC=HXPRC+FPX
        BYPRC=HYPRC+FPY
        BZPRC=HZPRC+FPZ

        RETURN
        END
C---------------------------------------------------------------------------------------
C
       SUBROUTINE SRC_PRC (IOPR,SC_SY,SC_PR,PHI,PS,X,Y,Z,BXSRC,BYSRC,
     *    BZSRC,BXPRC,BYPRC,BZPRC)
C
C   RETURNS FIELD COMPONENTS FROM A MODEL RING CURRENT, INCLUDING ITS SYMMETRIC PART
C     AND A PARTIAL RING CURRENT, CLOSED VIA BIRKELAND CURRENTS. BASED ON RESULTS, DESCRIBED
C     IN A PAPER "MODELING THE INNER MAGNETOSPHERE: ASYMMETRIC RING CURRENT AND REGION 2
C     BIRKELAND CURRENTS REVISITED" (JGR, DEC.2000).
C
C     IOPR -  A RING CURRENT CALCULATION FLAG (FOR LEAST-SQUARES FITTING ONLY):
C             IOPR=0 - BOTH SRC AND PRC FIELDS ARE CALCULATED
C             IOPR=1 - SRC ONLY
C             IOPR=2 - PRC ONLY
C
C     SC_SY &  SC_PR ARE SCALE FACTORS FOR THE ABOVE COMPONENTS;  TAKING SC<1 OR SC>1 MAKES THE CURRENTS
C                      SHRINK OR EXPAND, RESPECTIVELY.
C
C   PHI IS THE ROTATION ANGLE (RADIANS) OF THE PARTIAL RING CURRENT (MEASURED FROM MIDNIGHT TOWARD DUSK)
C
        IMPLICIT REAL*8 (A-H,O-Z)
c
c   1.  TRANSFORM TO TILTED COORDINATES (i.e., SM coordinates):
C
        CPS=DCOS(PS)
        SPS=DSIN(PS)

        XT=X*CPS-Z*SPS
        ZT=Z*CPS+X*SPS
C
C   2.  SCALE THE COORDINATES FOR THE SYMMETRIC AND PARTIAL RC COMPONENTS:
C
        XTS=XT/SC_SY    !  SYMMETRIC
        YTS=Y /SC_SY
        ZTS=ZT/SC_SY

        XTA=XT/SC_PR    !  PARTIAL
        YTA=Y /SC_PR
        ZTA=ZT/SC_PR
C
C   3.  CALCULATE COMPONENTS OF THE TOTAL FIELD IN THE TILTED (SOLAR-MAGNETIC) COORDINATE SYSTEM:
C
C==========   ONLY FOR LEAST SQUARES FITTING:
        BXS=0.D0
        BYS=0.D0
        BZS=0.D0
        BXA_S=0.D0
        BYA_S=0.D0
        BZA_S=0.D0
        BXA_QR=0.D0
        BYA_QR=0.D0
        BZA_Q=0.D0
C============================================
C
C    3a. SYMMETRIC FIELD:
C
        IF (IOPR.LE.1) CALL RC_SYMM(XTS,YTS,ZTS,BXS,BYS,BZS)
        IF (IOPR.EQ.0.OR.IOPR.EQ.2)
     *                 CALL PRC_SYMM(XTA,YTA,ZTA,BXA_S,BYA_S,BZA_S)

C    3b. ROTATE THE SCALED SM COORDINATES BY PHI AROUND ZSM AXIS AND CALCULATE QUADRUPOLE PRC FIELD
C         IN THOSE COORDS:

        CP=DCOS(PHI)
        SP=DSIN(PHI)
        XR=XTA*CP-YTA*SP
        YR=XTA*SP+YTA*CP

        IF (IOPR.EQ.0.OR.IOPR.EQ.2)
     *                 CALL PRC_QUAD(XR,YR,ZTA,BXA_QR,BYA_QR,BZA_Q)

C    3c. TRANSFORM THE QUADRUPOLE FIELD COMPONENTS BACK TO THE SM COORDS:
C
        BXA_Q= BXA_QR*CP+BYA_QR*SP
        BYA_Q=-BXA_QR*SP+BYA_QR*CP

C    3d. FIND THE TOTAL FIELD OF PRC (SYMM.+QUADR.) IN THE SM COORDS:
C
        BXP=BXA_S+BXA_Q
        BYP=BYA_S+BYA_Q
        BZP=BZA_S+BZA_Q
C
C   4.  TRANSFORM THE FIELDS OF BOTH PARTS OF THE RING CURRENT BACK TO THE GSM SYSTEM:
C
        BXSRC=BXS*CPS+BZS*SPS   !    SYMMETRIC RC
        BYSRC=BYS
        BZSRC=BZS*CPS-BXS*SPS
C
        BXPRC=BXP*CPS+BZP*SPS   !    PARTIAL RC
        BYPRC=BYP
        BZPRC=BZP*CPS-BXP*SPS
C
        RETURN
        END
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
      SUBROUTINE RC_SYMM (X,Y,Z,BX,BY,BZ)
      IMPLICIT  REAL * 8  (A - H, O - Z)
      DATA DS,DC/1.D-2,0.99994999875D0/, D/1.D-4/,DRD/5.D3/  ! DS=SIN(THETA) AT THE BOUNDARY OF THE LINEARITY
                                                                        REGION; DC=SQRT(1-DS**2);  DRD=1/(2*D)
      RHO2=X**2+Y**2
      R2=RHO2+Z**2
      R=DSQRT(R2)
      RP=R+D
      RM=R-D
      SINT=DSQRT(RHO2)/R
      COST=Z/R

      IF (SINT.LT.DS) THEN  !  TOO CLOSE TO THE Z-AXIS; USING A LINEAR APPROXIMATION A_PHI~SINT,
C                                    TO AVOID THE SINGULARITY PROBLEM
        A=AP(R,DS,DC)/DS
        DARDR=(RP*AP(RP,DS,DC)-RM*AP(RM,DS,DC))*DRD
        FXY=Z*(2.D0*A-DARDR)/(R*R2)
        BX=FXY*X
        BY=FXY*Y
        BZ=(2.D0*A*COST**2+DARDR*SINT**2)/R

       ELSE

        THETA=DATAN2(SINT,COST)
        TP=THETA+D
        TM=THETA-D
        SINTP=DSIN(TP)
        SINTM=DSIN(TM)
        COSTP=DCOS(TP)
        COSTM=DCOS(TM)
        BR=(SINTP*AP(R,SINTP,COSTP)-SINTM*AP(R,SINTM,COSTM))
     *       /(R*SINT)*DRD
        BT=(RM*AP(RM,SINT,COST)-RP*AP(RP,SINT,COST))/R*DRD
        FXY=(BR+BT*COST/SINT)/R
        BX=FXY*X
        BY=FXY*Y
        BZ=BR*COST-BT*SINT

      ENDIF

      RETURN
      END
c
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
      DOUBLE PRECISION FUNCTION AP(R,SINT,COST)
C
C      Calculates azimuthal component of the vector potential of the symmetric
c  part of the model ring current.
C
      IMPLICIT  REAL * 8  (A - H, O - Z)
      LOGICAL PROX   !  INDICATES WHETHER WE ARE TOO CLOSE TO THE AXIS OF SYMMETRY, WHERE THE INVERSION
C                                                             OF DIPOLAR COORDINATES BECOMES INACCURATE
      DATA A1,A2,RRC1,DD1,RRC2,DD2,P1,R1,DR1,DLA1,P2,R2,DR2,DLA2,P3,
     *R3,DR3/-563.3722359,425.0891691,4.150588549,2.266150226,
     * 3.334503403,3.079071195,.02602428295,8.937790598,3.327934895,
     *.4487061833,.09125832351,6.243029867,1.750145910,.4181957162,
     *.06106691992,2.079908581,.6828548533/

      PROX=.FALSE.
      SINT1=SINT
      COST1=COST
      IF (SINT1.LT.1.D-2) THEN  !  TOO CLOSE TO Z-AXIS;  USE LINEAR INTERPOLATION BETWEEN SINT=0 & SINT=0.01
        SINT1=1.D-2
        COST1=.99994999875
        PROX=.TRUE.
      ENDIF

         ALPHA=SINT1**2/R         !  R,THETA -> ALPHA,GAMMA
         GAMMA=COST1/R**2

         ARG1=-((R-R1)/DR1)**2-(COST1/DLA1)**2
         ARG2=-((R-R2)/DR2)**2-(COST1/DLA2)**2
         ARG3=-((R-R3)/DR3)**2

         IF (ARG1.LT.-500.D0) THEN        !   TO PREVENT "FLOATING UNDERFLOW" CRASHES
           DEXP1=0.D0
         ELSE
           DEXP1=DEXP(ARG1)
         ENDIF

         IF (ARG2.LT.-500.D0) THEN
           DEXP2=0.D0
         ELSE
           DEXP2=DEXP(ARG2)
         ENDIF

         IF (ARG3.LT.-500.D0) THEN
           DEXP3=0.D0
         ELSE
           DEXP3=DEXP(ARG3)
         ENDIF


         ALPHA_S=ALPHA*(1.D0+P1*DEXP1+P2*DEXP2+P3*DEXP3)     !  ALPHA -> ALPHA_S  (DEFORMED)

         GAMMA_S=GAMMA
         GAMMAS2=GAMMA_S**2


         ALSQH=ALPHA_S**2/2.D0            !  ALPHA_S,GAMMA_S -> RS,SINTS,COSTS
         F=64.D0/27.D0*GAMMAS2+ALSQH**2
         Q=(DSQRT(F)+ALSQH)**(1.D0/3.D0)
         C=Q-4.D0*GAMMAS2**(1.D0/3.D0)/(3.D0*Q)
         IF (C.LT.0.D0) C=0.D0
         G=DSQRT(C**2+4.D0*GAMMAS2**(1.D0/3.D0))
         RS=4.D0/((DSQRT(2.D0*G-C)+DSQRT(C))*(G+C))
         COSTS=GAMMA_S*RS**2
         SINTS=DSQRT(1.D0-COSTS**2)
         RHOS=RS*SINTS
         RHOS2=RHOS**2
         ZS=RS*COSTS
C
c  1st loop:

         P=(RRC1+RHOS)**2+ZS**2+DD1**2
         XK2=4.D0*RRC1*RHOS/P
         XK=SQRT(XK2)
         XKRHO12=XK*SQRT(RHOS)     !   SEE NB#4, P.3
C
      XK2S=1.D0-XK2
      DL=DLOG(1.D0/XK2S)
      ELK=1.38629436112d0+XK2S*(0.09666344259D0+XK2S*(0.03590092383+
     *     XK2S*(0.03742563713+XK2S*0.01451196212))) +DL*
     *     (0.5D0+XK2S*(0.12498593597D0+XK2S*(0.06880248576D0+
     *     XK2S*(0.03328355346D0+XK2S*0.00441787012D0))))
      ELE=1.D0+XK2S*(0.44325141463D0+XK2S*(0.0626060122D0+XK2S*
     *      (0.04757383546D0+XK2S*0.01736506451D0))) +DL*
     *     XK2S*(0.2499836831D0+XK2S*(0.09200180037D0+XK2S*
     *       (0.04069697526D0+XK2S*0.00526449639D0)))
C
      APHI1=((1.D0-XK2*0.5D0)*ELK-ELE)/XKRHO12
c
c  2nd loop:

         P=(RRC2+RHOS)**2+ZS**2+DD2**2
         XK2=4.D0*RRC2*RHOS/P
         XK=SQRT(XK2)
         XKRHO12=XK*SQRT(RHOS)     !   SEE NB#4, P.3
C
      XK2S=1.D0-XK2
      DL=DLOG(1.D0/XK2S)
      ELK=1.38629436112d0+XK2S*(0.09666344259D0+XK2S*(0.03590092383+
     *     XK2S*(0.03742563713+XK2S*0.01451196212))) +DL*
     *     (0.5D0+XK2S*(0.12498593597D0+XK2S*(0.06880248576D0+
     *     XK2S*(0.03328355346D0+XK2S*0.00441787012D0))))
      ELE=1.D0+XK2S*(0.44325141463D0+XK2S*(0.0626060122D0+XK2S*
     *      (0.04757383546D0+XK2S*0.01736506451D0))) +DL*
     *     XK2S*(0.2499836831D0+XK2S*(0.09200180037D0+XK2S*
     *       (0.04069697526D0+XK2S*0.00526449639D0)))
C
       APHI2=((1.D0-XK2*0.5D0)*ELK-ELE)/XKRHO12

       AP=A1*APHI1+A2*APHI2
       IF (PROX) AP=AP*SINT/SINT1   !   LINEAR INTERPOLATION, IF TOO CLOSE TO THE Z-AXIS
C
       RETURN
       END
c
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C
      SUBROUTINE PRC_SYMM (X,Y,Z,BX,BY,BZ)
      IMPLICIT  REAL * 8  (A - H, O - Z)
      DATA DS,DC/1.D-2,0.99994999875D0/, D/1.D-4/,DRD/5.D3/  ! DS=SIN(THETA) AT THE BOUNDARY OF THE LINEARITY
                                                                        REGION; DC=SQRT(1-DS**2);  DRD=1/(2*D)
      RHO2=X**2+Y**2
      R2=RHO2+Z**2
      R=DSQRT(R2)
      RP=R+D
      RM=R-D
      SINT=DSQRT(RHO2)/R
      COST=Z/R

      IF (SINT.LT.DS) THEN  !  TOO CLOSE TO THE Z-AXIS; USING A LINEAR APPROXIMATION A_PHI~SINT,
C                                    TO AVOID THE SINGULARITY PROBLEM
        A=APPRC(R,DS,DC)/DS
        DARDR=(RP*APPRC(RP,DS,DC)-RM*APPRC(RM,DS,DC))*DRD
        FXY=Z*(2.D0*A-DARDR)/(R*R2)
        BX=FXY*X
        BY=FXY*Y
        BZ=(2.D0*A*COST**2+DARDR*SINT**2)/R

       ELSE

        THETA=DATAN2(SINT,COST)
        TP=THETA+D
        TM=THETA-D
        SINTP=DSIN(TP)
        SINTM=DSIN(TM)
        COSTP=DCOS(TP)
        COSTM=DCOS(TM)
        BR=(SINTP*APPRC(R,SINTP,COSTP)-SINTM*APPRC(R,SINTM,COSTM))
     *       /(R*SINT)*DRD
        BT=(RM*APPRC(RM,SINT,COST)-RP*APPRC(RP,SINT,COST))/R*DRD
        FXY=(BR+BT*COST/SINT)/R
        BX=FXY*X
        BY=FXY*Y
        BZ=BR*COST-BT*SINT

      ENDIF

      RETURN
      END
C
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C
C
         SUBROUTINE PRC_QUAD (X,Y,Z,BX,BY,BZ)
C
C  CALCULATES COMPONENTS OF THE FIELD FROM THE "QUADRUPOLE" COMPONENT OF THE PRC
C
         IMPLICIT  REAL * 8  (A - H, O - Z)

         DATA D,DD/1.D-4,2.D-4/, DS/1.D-2/,DC/0.99994999875D0/

         RHO2=X**2+Y**2
         R=DSQRT(RHO2+Z**2)
         RHO=DSQRT(RHO2)
         SINT=RHO/R
         COST=Z/R
         RP=R+D
         RM=R-D

         IF (SINT.GT.DS) THEN
           CPHI=X/RHO
           SPHI=Y/RHO
           BR=BR_PRC_Q(R,SINT,COST)
           BT=BT_PRC_Q(R,SINT,COST)
           DBRR=(BR_PRC_Q(RP,SINT,COST)-BR_PRC_Q(RM,SINT,COST))/DD
           THETA=DATAN2(SINT,COST)
           TP=THETA+D
           TM=THETA-D
           SINTP=DSIN(TP)
           COSTP=DCOS(TP)
           SINTM=DSIN(TM)
           COSTM=DCOS(TM)
           DBTT=(BT_PRC_Q(R,SINTP,COSTP)-BT_PRC_Q(R,SINTM,COSTM))/DD
           BX=SINT*(BR+(BR+R*DBRR+DBTT)*SPHI**2)+COST*BT
           BY=-SINT*SPHI*CPHI*(BR+R*DBRR+DBTT)
           BZ=(BR*COST-BT*SINT)*CPHI
         ELSE
           ST=DS
           CT=DC
           IF (Z.LT.0.D0) CT=-DC
           THETA=DATAN2(ST,CT)
           TP=THETA+D
           TM=THETA-D
           SINTP=DSIN(TP)
           COSTP=DCOS(TP)
           SINTM=DSIN(TM)
           COSTM=DCOS(TM)
           BR=BR_PRC_Q(R,ST,CT)
           BT=BT_PRC_Q(R,ST,CT)
           DBRR=(BR_PRC_Q(RP,ST,CT)-BR_PRC_Q(RM,ST,CT))/DD
           DBTT=(BT_PRC_Q(R,SINTP,COSTP)-BT_PRC_Q(R,SINTM,COSTM))/DD
           FCXY=R*DBRR+DBTT
           BX=(BR*(X**2+2.D0*Y**2)+FCXY*Y**2)/(R*ST)**2+BT*COST
           BY=-(BR+FCXY)*X*Y/(R*ST)**2
           BZ=(BR*COST/ST-BT)*X/R
         ENDIF

         RETURN
         END
c
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
      DOUBLE PRECISION FUNCTION BR_PRC_Q (R,SINT,COST)
C
Calculates the radial component of the "quadrupole" part of the model partial ring current.
C
      IMPLICIT  REAL * 8  (A - H, O - Z)

      DATA A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15,A16,A17,   ! ALL LINEAR PARAMETERS HERE
     * A18,XK1,AL1,DAL1,B1,BE1,XK2,AL2,DAL2,B2,BE2,XK3,XK4,AL3,DAL3,B3,  ! WERE MULTIPLIED BY 0.1,
     * BE3,AL4,DAL4,DG1,AL5,DAL5,DG2,C1,C2,C3,AL6,DAL6,DRM/-21.2666329,  ! SO THAT THEY CORRESPOND TO P_0=1 nPa,
     *32.24527521,-6.062894078,7.515660734,233.7341288,-227.1195714,     ! RATHER THAN THE ORIGINAL VALUE OF 10 nPa
     *8.483233889,16.80642754,-24.63534184,9.067120578,-1.052686913,     ! ASSUMED IN THE BIOT-SAVART INTEGRAL
     *-12.08384538,18.61969572,-12.71686069,47017.35679,-50646.71204,
     *7746.058231,1.531069371,2.318824273,.1417519429,.6388013110E-02,
     *5.303934488,4.213397467,.7955534018,.1401142771,.2306094179E-01,
     *3.462235072,2.568743010,3.477425908,1.922155110,.1485233485,
     *.2319676273E-01,7.830223587,8.492933868,.1295221828,.01753008801,
     *.01125504083,.1811846095,.04841237481,.01981805097,6.557801891,
     *6.348576071,5.744436687,.2265212965,.1301957209,.5654023158/

        SINT2=SINT**2
        COST2=COST**2
        SC=SINT*COST
        ALPHA=SINT2/R
        GAMMA=COST/R**2

        CALL FFS(ALPHA,AL1,DAL1,F,FA,FS)
        D1=SC*F**XK1/((R/B1)**BE1+1.D0)
        D2=D1*COST2

        CALL FFS(ALPHA,AL2,DAL2,F,FA,FS)
        D3=SC*FS**XK2/((R/B2)**BE2+1.D0)
        D4=D3*COST2

        CALL FFS(ALPHA,AL3,DAL3,F,FA,FS)
        D5=SC*(ALPHA**XK3)*(FS**XK4)/((R/B3)**BE3+1.D0)
        D6=D5*COST2

        ARGA=((ALPHA-AL4)/DAL4)**2+1.D0
        ARGG=1.D0+(GAMMA/DG1)**2

        D7=SC/ARGA/ARGG
        D8=D7/ARGA
        D9=D8/ARGA
        D10=D9/ARGA

        ARGA=((ALPHA-AL5)/DAL5)**2+1.D0
        ARGG=1.D0+(GAMMA/DG2)**2

        D11=SC/ARGA/ARGG
        D12=D11/ARGA
        D13=D12/ARGA
        D14=D13/ARGA


        D15=SC/(R**4+C1**4)
        D16=SC/(R**4+C2**4)*COST2
        D17=SC/(R**4+C3**4)*COST2**2

        CALL FFS(ALPHA,AL6,DAL6,F,FA,FS)
        D18=SC*FS/(1.D0+((R-1.2D0)/DRM)**2)

        BR_PRC_Q=A1*D1+A2*D2+A3*D3+A4*D4+A5*D5+A6*D6+A7*D7+A8*D8+A9*D9+
     *  A10*D10+A11*D11+A12*D12+A13*D13+A14*D14+A15*D15+A16*D16+A17*D17+
     *   A18*D18
C
        RETURN
        END
c
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
        DOUBLE PRECISION FUNCTION BT_PRC_Q (R,SINT,COST)
C
Calculates the Theta component of the "quadrupole" part of the model partial ring current.
C
        IMPLICIT  REAL * 8  (A - H, O - Z)

      DATA A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15,A16,A17,  ! ALL LINEAR PARAMETERS HERE
     *XK1,AL1,DAL1,B1,BE1,XK2,AL2,DAL2,BE2,XK3,XK4,AL3,DAL3,B3,BE3,AL4, ! WERE MULTIPLIED BY 0.1,
     *DAL4,DG1,AL5,DAL5,DG2,C1,C2,C3/12.74640393,-7.516393516,          ! SO THAT THEY CORRESPOND TO P_0=1 nPa,
     *-5.476233865,3.212704645,-59.10926169,46.62198189,-.01644280062,  ! RATHER THAN THE ORIGINAL VALUE OF 10 nPa
     *.1234229112,-.08579198697,.01321366966,.8970494003,9.136186247,   ! ASSUMED IN THE BIOT-SAVART INTEGRAL
     *-38.19301215,21.73775846,-410.0783424,-69.90832690,-848.8543440,
     *1.243288286,.2071721360,.05030555417,7.471332374,3.180533613,
     *1.376743507,.1568504222,.02092910682,1.985148197,.3157139940,
     *1.056309517,.1701395257,.1019870070,6.293740981,5.671824276,
     *.1280772299,.02189060799,.01040696080,.1648265607,.04701592613,
     *.01526400086,12.88384229,3.361775101,23.44173897/

        SINT2=SINT**2
        COST2=COST**2
        SC=SINT*COST
        ALPHA=SINT2/R
        GAMMA=COST/R**2

        CALL FFS(ALPHA,AL1,DAL1,F,FA,FS)
        D1=F**XK1/((R/B1)**BE1+1.D0)
        D2=D1*COST2

        CALL FFS(ALPHA,AL2,DAL2,F,FA,FS)
        D3=FA**XK2/R**BE2
        D4=D3*COST2

        CALL FFS(ALPHA,AL3,DAL3,F,FA,FS)
        D5=FS**XK3*ALPHA**XK4/((R/B3)**BE3+1.D0)
        D6=D5*COST2

        CALL FFS(GAMMA,0.D0,DG1,F,FA,FS)
        FCC=(1.D0+((ALPHA-AL4)/DAL4)**2)
        D7 =1.D0/FCC*FS
        D8 =D7/FCC
        D9 =D8/FCC
        D10=D9/FCC

        ARG=1.D0+((ALPHA-AL5)/DAL5)**2
        D11=1.D0/ARG/(1.D0+(GAMMA/DG2)**2)
        D12=D11/ARG
        D13=D12/ARG
        D14=D13/ARG

        D15=1.D0/(R**4+C1**2)
        D16=COST2/(R**4+C2**2)
        D17=COST2**2/(R**4+C3**2)
C
        BT_PRC_Q=A1*D1+A2*D2+A3*D3+A4*D4+A5*D5+A6*D6+A7*D7+A8*D8+A9*D9+
     *   A10*D10+A11*D11+A12*D12+A13*D13+A14*D14+A15*D15+A16*D16+A17*D17
C
       RETURN
       END
c
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

       SUBROUTINE FFS(A,A0,DA,F,FA,FS)
       IMPLICIT  REAL * 8  (A - H, O - Z)
       SQ1=DSQRT((A+A0)**2+DA**2)
       SQ2=DSQRT((A-A0)**2+DA**2)
       FA=2.D0/(SQ1+SQ2)
       F=FA*A
       FS=0.5D0*(SQ1+SQ2)/(SQ1*SQ2)*(1.D0-F*F)
       RETURN
       END
C
C||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
C
C
         SUBROUTINE RC_SHIELD (A,PS,X_SC,X,Y,Z,BX,BY,BZ)
C
C   COMPUTES THE COMPONENTS OF THE SHIELDING FIELD FOR THE RING CURRENT
C       (EITHER PARTIAL OR AXISYMMETRICAL)
C   INPUT:   A - AN ARRAY CONTAINING THE HARMONIC COEFFICIENTS AND NONLINEAR PARAMETERS
C            PS - GEODIPOLE TILT ANGLE IN RADIANS
C            X_SC - SCALING FACTOR ( X_SC> 1 AND X_SC< 1 CORRESPOND TO LARGER/SMALLER
C                  RING CURRENT, RESP.)
C            X,Y,Z - POSITION IN RE (GSM COORDS)
C   OUTPUT:  BX,BY,BZ - SHIELDING FIELD COMPONENTS (GSM)
C
       IMPLICIT  REAL * 8  (A - H, O - Z)
         DIMENSION A(86)
C
         FAC_SC=(X_SC+1.D0)**3
C
         CPS=DCOS(PS)
         SPS=DSIN(PS)

         S3PS=2.D0*CPS
C
         PST1=PS*A(85)
         PST2=PS*A(86)

         ST1=DSIN(PST1)
         CT1=DCOS(PST1)
         ST2=DSIN(PST2)
         CT2=DCOS(PST2)

         X1=X*CT1-Z*ST1
         Z1=X*ST1+Z*CT1
         X2=X*CT2-Z*ST2
         Z2=X*ST2+Z*CT2
C
         L=0
         GX=0.D0
         GY=0.D0
         GZ=0.D0
C
         DO 1 M=1,2     !    M=1 IS FOR THE 1ST SUM ("PERP." SYMMETRY)
C                           AND M=2 IS FOR THE SECOND SUM ("PARALL." SYMMETRY)
             DO 2 I=1,3
                  P=A(72+I)
                  Q=A(78+I)
                  CYPI=DCOS(Y/P)
                  CYQI=DCOS(Y/Q)
                  SYPI=DSIN(Y/P)
                  SYQI=DSIN(Y/Q)
C
                DO 3 K=1,3
                   R=A(75+K)
                   S=A(81+K)
                   SZRK=DSIN(Z1/R)
                   CZSK=DCOS(Z2/S)
                   CZRK=DCOS(Z1/R)
                   SZSK=DSIN(Z2/S)
                     SQPR=DSQRT(1.D0/P**2+1.D0/R**2)
                     SQQS=DSQRT(1.D0/Q**2+1.D0/S**2)
                        EPR=DEXP(X1*SQPR)
                        EQS=DEXP(X2*SQQS)
C
                  DO 4 N=1,2  ! N=1 IS FOR THE FIRST PART OF EACH COEFFICIENT
C                                AND N=2 IS FOR THE SECOND ONE

                    DO 5 NN=1,2 !   NN = 1,2 FURTHER SPLITS THE COEFFICIENTS INTO 2 PARTS,
C                                         TO TAKE INTO ACCOUNT THE SCALE FACTOR DEPENDENCE

                    IF (M.EQ.1) THEN
                         FX=-SQPR*EPR*CYPI*SZRK  *FAC_SC
                         FY=EPR*SYPI*SZRK/P   *FAC_SC
                         FZ=-EPR*CYPI*CZRK/R  *FAC_SC
                       IF (N.EQ.1) THEN
                         IF (NN.EQ.1) THEN
                          HX=FX
                          HY=FY
                          HZ=FZ
                         ELSE
                          HX=FX*X_SC
                          HY=FY*X_SC
                          HZ=FZ*X_SC
                         ENDIF
                       ELSE
                         IF (NN.EQ.1) THEN
                          HX=FX*CPS
                          HY=FY*CPS
                          HZ=FZ*CPS
                         ELSE
                          HX=FX*CPS*X_SC
                          HY=FY*CPS*X_SC
                          HZ=FZ*CPS*X_SC
                         ENDIF
                       ENDIF

                     ELSE                            !   M.EQ.2
                         FX=-SPS*SQQS*EQS*CYQI*CZSK  *FAC_SC
                         FY=SPS/Q*EQS*SYQI*CZSK   *FAC_SC
                         FZ=SPS/S*EQS*CYQI*SZSK   *FAC_SC
                       IF (N.EQ.1) THEN
                        IF (NN.EQ.1) THEN
                          HX=FX
                          HY=FY
                          HZ=FZ
                        ELSE
                          HX=FX*X_SC
                          HY=FY*X_SC
                          HZ=FZ*X_SC
                        ENDIF
                       ELSE
                        IF (NN.EQ.1) THEN
                         HX=FX*S3PS
                         HY=FY*S3PS
                         HZ=FZ*S3PS
                        ELSE
                         HX=FX*S3PS*X_SC
                         HY=FY*S3PS*X_SC
                         HZ=FZ*S3PS*X_SC
                        ENDIF
                       ENDIF
                  ENDIF
       L=L+1

       IF (M.EQ.1) THEN
       HXR=HX*CT1+HZ*ST1
       HZR=-HX*ST1+HZ*CT1
       ELSE
       HXR=HX*CT2+HZ*ST2
       HZR=-HX*ST2+HZ*CT2
       ENDIF

       GX=GX+HXR*A(L)
       GY=GY+HY *A(L)
  5    GZ=GZ+HZR*A(L)

  4   CONTINUE
  3   CONTINUE
  2   CONTINUE
  1   CONTINUE

      BX=GX
      BY=GY
      BZ=GZ

      RETURN
      END
C
c===========================================================================
c
      SUBROUTINE DIPOLE_01 (PS,X,Y,Z,BX,BY,BZ)
C
C     THIS IS A DOUBLE PRECISION ROUTINE, OTHERWISE IDENTICAL TO THE S/R DIP OF GEOPACK
C
C     CALCULATES GSM COMPONENTS OF A GEODIPOLE FIELD WITH THE DIPOLE MOMENT
C     CORRESPONDING TO THE EPOCH OF 2000.
C
C------INPUT PARAMETERS:
C       PS - GEODIPOLE TILT ANGLE IN RADIANS,
C       X,Y,Z - GSM COORDINATES IN RE (1 RE = 6371.2 km)
C
C----OUTPUT PARAMETERS:
C     BX,BY,BZ - FIELD COMPONENTS IN GSM SYSTEM, IN NANOTESLA.
C
C    LAST MODIFICATION: JAN. 5, 2001. THE VALUE OF THE DIPOLE MOMENT WAS UPDATED TO 2000.
C      AND A "SAVE" STATEMENT HAS BEEN ADDED, TO AVOID POTENTIAL PROBLEMS WITH SOME
C      FORTRAN COMPILERS
C
C    WRITTEN BY: N. A. TSYGANENKO
C
      IMPLICIT REAL*8 (A-H,O-Z)
      SAVE M,PSI
      DATA M,PSI/0,5.D0/
      IF(M.EQ.1.AND.DABS(PS-PSI).LT.1.D-5) GOTO 1   !   THIS IS TO AVOID MULTIPLE CALCULATIONS
      SPS=DSIN(PS)                                  !   OF SIN(PS) AND COS(PS), IF THE ANGLE PS
      CPS=DCOS(PS)                                  !   REMAINS UNCHANGED
      PSI=PS
      M=1
  1   P=X**2
      U=Z**2
      V=3.D0*Z*X
      T=Y**2
      Q=30115.D0/DSQRT(P+T+U)**5
      BX=Q*((T+U-2.D0*P)*SPS-V*CPS)
      BY=-3.D0*Y*Q*(X*SPS+Z*CPS)
      BZ=Q*((P+T-2.D0*U)*CPS-V*SPS)
      RETURN
      END



c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


c====================================================================================
c
c
      SUBROUTINE T04_s (IOPT,PARMOD,PS,X,Y,Z,BX,BY,BZ)
c
c     ASSEMBLED:  MARCH 25, 2004; UPDATED:  AUGUST 2 & 31, DECEMBER 27, 2004.
C
c--------------------------------------------------------------------
C   A DATA-BASED MODEL OF THE EXTERNAL (I.E., WITHOUT EARTH'S CONTRIBUTION) PART OF THE
C   MAGNETOSPHERIC MAGNETIC FIELD, CALIBRATED BY
C    (1) SOLAR WIND PRESSURE PDYN (NANOPASCALS),
C    (2) DST (NANOTESLA),
C    (3) BYIMF,
C    (4) BZIMF (NANOTESLA)
C    (5-10)   INDICES W1 - W6, CALCULATED AS TIME INTEGRALS FROM THE BEGINNING OF A STORM
c               SEE THE REFERENCE (3) BELOW, FOR A DETAILED DEFINITION OF THOSE VARIABLES
C
c   THE ABOVE 10 INPUT PARAMETERS SHOULD BE PLACED IN THE ELEMENTS
c   OF THE ARRAY PARMOD(10).
C
C   THE REST OF THE INPUT VARIABLES ARE: THE GEODIPOLE TILT ANGLE PS (RADIANS),
C        X,Y,Z -  GSM POSITION (RE)
C
c   IOPT IS A DUMMY INPUT PARAMETER, INCLUDED TO MAKE THIS SUBROUTINE
C   COMPATIBLE WITH THE TRACING SOFTWARE PACKAGE (GEOPACK). IN THIS MODEL,
C   THE PARAMETER IOPT DOES NOT AFFECT THE OUTPUT FIELD.
c
C*******************************************************************************************
c** ATTENTION:  THE MODEL IS BASED ON DATA TAKEN SUNWARD FROM X=-15Re, AND HENCE BECOMES   *
C**              INVALID AT LARGER TAILWARD DISTANCES !!!                                  *
C*******************************************************************************************
C
c   OUTPUT:  GSM COMPONENTS OF THE EXTERNAL MAGNETIC FIELD (BX,BY,BZ, nanotesla)
C            COMPUTED AS A SUM OF CONTRIBUTIONS FROM PRINCIPAL FIELD SOURCES
C
c  (C) Copr. 2004, Nikolai A. Tsyganenko, USRA/Code 695.1, NASA GSFC
c      Greenbelt, MD 20771, USA
c
C                            REFERENCES:
C
C  (1)   N. A. Tsyganenko, A new data-based model of the near magnetosphere magnetic field:
c       1. Mathematical structure.
c       2. Parameterization and fitting to observations.  JGR v. 107(A8), 1176/1179, doi:10.1029/2001JA000219/220, 2002.
c
c  (2)  N. A. Tsyganenko, H. J. Singer, J. C. Kasper, Storm-time distortion of the
c           inner magnetosphere: How severe can it get ?  JGR v. 108(A5), 1209, doi:10.1029/2002JA009808, 2003.

c  (3)  N. A. Tsyganenko and M. I. Sitnov, Modeling the dynamics of the inner
c           magnetosphere during strong geomagnetic storms,  JGR v. 110, 2005, in press.
c----------------------------------------------------------------------
c
      IMPLICIT  REAL * 8  (A - H, O - Z)
C
      REAL*8 PARMOD(10),PS,X,Y,Z,BX,BY,BZ
      REAL*8 A(69),PDYN,DST_AST,BXIMF,BYIMF,BZIMF,W1,W2,W3,W4,W5,W6,
     *  PSS,XX,YY,ZZ,BXCF,BYCF,BZCF,BXT1,BYT1,BZT1,BXT2,BYT2,BZT2,
     *  BXSRC,BYSRC,BZSRC,BXPRC,BYPRC,BZPRC, BXR11,BYR11,BZR11,
     *  BXR12,BYR12,BZR12,BXR21,BYR21,BZR21,BXR22,BYR22,BZR22,HXIMF,
     *  HYIMF,HZIMF,BBX,BBY,BBZ
C
      DATA A/1.00000,5.19884,0.923524,8.68111,0.00000,-6.44922,11.3109,
     * -3.84555,0.00000,0.558081,0.937044,0.00000,0.772433,0.687241,
     * 0.00000,0.320369,1.22531,-0.432246E-01,-0.382436,0.457468,
     * 0.741917,0.227194,0.154269,5.75196,22.3113,10.3526,64.3312,
     * 1.01977,-0.200859E-01,0.971643,0.295525E-01,1.01032,0.215561,
     * 1.50059,0.730898E-01,1.93625,1.74545,1.29533,0.714744,0.391687,
     * 3.31283,75.0127,6.36283,4.43561,0.387801,0.699661,0.305352E-01,
     * 0.581002,1.14671,0.876060,0.386060,0.801831,0.874315,0.463634,
     * 0.175077,0.673053,0.388341,2.32074,1.32373,0.419800,1.24968,
     * 1.28903,.409286,1.57622,.690036,1.28836,2.4054,.528557,.564247/

      DATA IOPGEN,IOPTT,IOPB,IOPR/0,0,0,0/
C
      PDYN=PARMOD(1)
      DST_AST=PARMOD(2)*0.8-13.*SQRT(PDYN)
      BYIMF=PARMOD(3)
      BZIMF=PARMOD(4)
C
      W1=PARMOD (5)
      W2=PARMOD (6)
      W3=PARMOD (7)
      W4=PARMOD (8)
      W5=PARMOD (9)
      W6=PARMOD(10)

      PSS=PS
      XX=X
      YY=Y
      ZZ=Z
C
      CALL TSY04_EXTERN (IOPGEN,IOPTT,IOPB,IOPR,A,79,PDYN,DST_AST,
     + BXIMF,BYIMF,
     + BZIMF,W1,W2,W3,W4,W5,W6,PSS,XX,YY,ZZ,BXCF,BYCF,BZCF,BXT1,BYT1,
     + BZT1,BXT2,BYT2,BZT2,BXSRC,BYSRC,BZSRC,BXPRC,BYPRC,BZPRC, BXR11,
     + BYR11,BZR11,BXR12,BYR12,BZR12,BXR21,BYR21,BZR21,BXR22,BYR22,
     + BZR22,HXIMF,HYIMF,HZIMF,BBX,BBY,BBZ)
C
      BX=BBX
      BY=BBY
      BZ=BBZ
C
      RETURN
      END

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
      SUBROUTINE TSY04_EXTERN (IOPGEN,IOPT,IOPB,IOPR,A,NTOT,
     *  PDYN,DST,BXIMF,BYIMF,BZIMF,W1,W2,W3,W4,W5,W6,PS,X,Y,Z,
     *  BXCF,BYCF,BZCF,BXT1,BYT1,BZT1,BXT2,BYT2,BZT2,
     *  BXSRC,BYSRC,BZSRC,BXPRC,BYPRC,BZPRC, BXR11,BYR11,BZR11,
     *  BXR12,BYR12,BZR12,BXR21,BYR21,BZR21,BXR22,BYR22,BZR22,HXIMF,
     *  HYIMF,HZIMF,BX,BY,BZ)
C
C   IOPGEN - GENERAL OPTION FLAG:  IOPGEN=0 - CALCULATE TOTAL FIELD
C                                  IOPGEN=1 - DIPOLE SHIELDING ONLY
C                                  IOPGEN=2 - TAIL FIELD ONLY
C                                  IOPGEN=3 - BIRKELAND FIELD ONLY
C                                  IOPGEN=4 - RING CURRENT FIELD ONLY
C                                  IOPGEN=5 - INTERCONNECTION FIELD ONLY
C
C   IOPT -  TAIL FIELD FLAG:       IOPT=0  -  BOTH MODES
C                                  IOPT=1  -  MODE 1 ONLY
C                                  IOPT=2  -  MODE 2 ONLY
C
C   IOPB -  BIRKELAND FIELD FLAG:  IOPB=0  -  ALL 4 TERMS
C                                  IOPB=1  -  REGION 1, MODES 1 AND 2
C                                  IOPB=2  -  REGION 2, MODES 1 AND 2
C
C   IOPR -  RING CURRENT FLAG:     IOPR=0  -  BOTH SRC AND PRC
C                                  IOPR=1  -  SRC ONLY
C                                  IOPR=2  -  PRC ONLY
C
      IMPLICIT  REAL * 8  (A - H, O - Z)
C
      DIMENSION A(NTOT)
C
      COMMON/TAIL/ DXSHIFT1,DXSHIFT2,D,DELTADY  ! THE COMMON BLOCKS FORWARD NONLINEAR PARAMETERS
      COMMON/BIRKPAR/ XKAPPA1,XKAPPA2
      COMMON/RCPAR/ SC_SY,SC_AS,PHI
      COMMON/TSY04_G/ G
      COMMON/TSY04_RH0/ RH0
C
C
      DATA A0_A,A0_S0,A0_X0 /34.586D0,1.1960D0,3.4397D0/   !   SHUE ET AL. PARAMETERS
      DATA DSIG /0.005D0/, RH0,RH2 /7.5D0,-5.2D0/
c
      XAPPA=(PDYN/2.)**A(23)   !  OVERALL SCALING PARAMETER
      RH0=7.5                  !  TAIL HINGING DISTANCE
c
      G=  35.0                 !  TAIL WARPING PARAMETER

      XAPPA3=XAPPA**3

      XX=X*XAPPA
      YY=Y*XAPPA
      ZZ=Z*XAPPA
C
      SPS=DSIN(PS)
c
      X0=A0_X0/XAPPA
      AM=A0_A/XAPPA
      S0=A0_S0
c
C  CALCULATE "IMF" COMPONENTS OUTSIDE THE MAGNETOPAUSE LAYER (HENCE BEGIN WITH "O")
C  THEY ARE NEEDED ONLY IF THE POINT (X,Y,Z) IS WITHIN THE TRANSITION MAGNETOPAUSE LAYER
C  OR OUTSIDE THE MAGNETOSPHERE:
C
      FACTIMF=A(20)
c
      OIMFX=0.D0
      OIMFY=BYIMF*FACTIMF
      OIMFZ=BZIMF*FACTIMF
c
      R=DSQRT(X**2+Y**2+Z**2)
      XSS=X
      ZSS=Z

  1   XSOLD=XSS      !   BEGIN ITERATIVE SEARCH OF UNWARPED COORDS (TO FIND SIGMA)
      ZSOLD=ZSS

      RH=RH0+RH2*(ZSS/R)**2
      SINPSAS=SPS/(1.D0+(R/RH)**3)**0.33333333D0
      COSPSAS=DSQRT(1.D0-SINPSAS**2)
      ZSS=X*SINPSAS+Z*COSPSAS
      XSS=X*COSPSAS-Z*SINPSAS
      DD=DABS(XSS-XSOLD)+DABS(ZSS-ZSOLD)
      IF (DD.GT.1.D-6) GOTO 1
C                                END OF ITERATIVE SEARCH
      RHO2=Y**2+ZSS**2
      ASQ=AM**2
      XMXM=AM+XSS-X0
      IF (XMXM.LT.0.) XMXM=0. ! THE BOUNDARY IS A CYLINDER TAILWARD OF X=X0-AM
      AXX0=XMXM**2
      ARO=ASQ+RHO2
      SIGMA=DSQRT((ARO+AXX0+SQRT((ARO+AXX0)**2-4.*ASQ*AXX0))/(2.*ASQ))
C
C   NOW, THERE ARE THREE POSSIBLE CASES:
C    (1) INSIDE THE MAGNETOSPHERE   (SIGMA
C    (2) IN THE BOUNDARY LAYER
C    (3) OUTSIDE THE MAGNETOSPHERE AND B.LAYER
C       FIRST OF ALL, CONSIDER THE CASES (1) AND (2):
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF (SIGMA.LT.S0+DSIG) THEN  !  CASES (1) OR (2); CALCULATE THE MODEL FIELD
C                                   (WITH THE POTENTIAL "PENETRATED" INTERCONNECTION FIELD):
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      IF (IOPGEN.LE.1) THEN
         CALL SHLCAR3X3(XX,YY,ZZ,PS,CFX,CFY,CFZ)         !  DIPOLE SHIELDING FIELD
         BXCF=CFX*XAPPA3
         BYCF=CFY*XAPPA3
         BZCF=CFZ*XAPPA3
      ELSE
         BXCF=0.D0
         BYCF=0.D0
         BZCF=0.D0
      ENDIF                                              !  DONE

      IF (IOPGEN.EQ.0.OR.IOPGEN.EQ.2) THEN
          DSTT=-20.
          IF (DST.LT.DSTT) DSTT=DST
          ZNAM=DABS(DSTT)**0.37
         DXSHIFT1=A(24)-A(25)/ZNAM
         DXSHIFT2=A(26)-A(27)/ZNAM
         D=A(36)*DEXP(-W1/A(37))  +A(69)
         DELTADY=4.7

         CALL DEFORMED (IOPT,PS,XX,YY,ZZ,                !  TAIL FIELD (THREE MODES)
     *    BXT1,BYT1,BZT1,BXT2,BYT2,BZT2)
      ELSE
         BXT1=0.D0
         BYT1=0.D0
         BZT1=0.D0
         BXT2=0.D0
         BYT2=0.D0
         BZT2=0.D0
      ENDIF

      IF (IOPGEN.EQ.0.OR.IOPGEN.EQ.3) THEN

          ZNAM=DABS(DST)
          IF (DST.GE.-20.D0) ZNAM=20.D0
          XKAPPA1=A(32)*(ZNAM/20.D0)**A(33)
          XKAPPA2=A(34)*(ZNAM/20.D0)**A(35)

         CALL BIRK_TOT_04 (IOPB,PS,XX,YY,ZZ,BXR11,BYR11,BZR11,BXR12,
     *   BYR12,BZR12,BXR21,BYR21,BZR21,BXR22,BYR22,BZR22)    !   BIRKELAND FIELD (TWO MODES FOR R1 AND TWO MODES FOR R2)
      ELSE
         BXR11=0.D0
         BYR11=0.D0
         BZR11=0.D0
         BXR21=0.D0
         BYR21=0.D0
         BZR21=0.D0
      ENDIF

      IF (IOPGEN.EQ.0.OR.IOPGEN.EQ.4) THEN
         PHI=A(38)

          ZNAM=DABS(DST)
          IF (DST.GE.-20.D0) ZNAM=20.D0
          SC_SY=A(28)*(20.D0/ZNAM)**A(29) *XAPPA    !
          SC_AS=A(30)*(20.D0/ZNAM)**A(31) *XAPPA    !  MULTIPLICATION  BY XAPPA IS MADE IN ORDER TO MAKE THE SRC AND PRC
                                                    !     SCALING COMPLETELY INDEPENDENT OF THE GENERAL SCALING DUE TO THE
C                                                         MAGNETOPAUSE COMPRESSION/EXPANSION                             !
C
         CALL FULL_RC(IOPR,PS,XX,YY,ZZ,BXSRC,BYSRC,BZSRC,BXPRC,BYPRC,
     *                                        BZPRC)  !  SHIELDED RING CURRENT (SRC AND PRC)
      ELSE
         BXSRC=0.D0
         BYSRC=0.D0
         BZSRC=0.D0
         BXPRC=0.D0
         BYPRC=0.D0
         BZPRC=0.D0
      ENDIF
C
      IF (IOPGEN.EQ.0.OR.IOPGEN.EQ.5) THEN
         HXIMF=0.D0
         HYIMF=BYIMF
         HZIMF=BZIMF   ! THESE ARE COMPONENTS OF THE PENETRATED FIELD PER UNIT OF THE PENETRATION COEFFICIENT.
C                             IN OTHER WORDS, THESE ARE DERIVATIVES OF THE PENETRATION FIELD COMPONENTS WITH RESPECT
C                             TO THE PENETRATION COEFFICIENT.   WE ASSUME THAT ONLY TRANSVERSE COMPONENT OF THE
C                             FIELD PENETRATES INSIDE.
       ELSE
         HXIMF=0.D0
         HYIMF=0.D0
         HZIMF=0.D0
       ENDIF
C
C-----------------------------------------------------------
C
C    NOW, ADD UP ALL THE COMPONENTS:
c
      DLP1=(PDYN/2.D0)**A(21)
      DLP2=(PDYN/2.D0)**A(22)

      TAMP1=A(2)+A(3)*DLP1+A(4)*A(39)*W1/DSQRT(W1**2+A(39)**2)+A(5)*DST
      TAMP2=A(6)+A(7)*DLP2+A(8)*A(40)*W2/DSQRT(W2**2+A(40)**2)+A(9)*DST
      A_SRC=A(10)+A(11)*A(41)*W3/DSQRT(W3**2+A(41)**2)
     *  +A(12)*DST
      A_PRC=A(13)+A(14)*A(42)*W4/DSQRT(W4**2+A(42)**2)
     *  +A(15)*DST
      A_R11=A(16)+A(17)*A(43)*W5/DSQRT(W5**2+A(43)**2)
      A_R21=A(18)+A(19)*A(44)*W6/DSQRT(W6**2+A(44)**2)

      BBX=A(1)*BXCF+TAMP1*BXT1+TAMP2*BXT2+A_SRC*BXSRC+A_PRC*BXPRC
     * +A_R11*BXR11+A_R21*BXR21+A(20)*HXIMF

      BBY=A(1)*BYCF+TAMP1*BYT1+TAMP2*BYT2+A_SRC*BYSRC+A_PRC*BYPRC
     * +A_R11*BYR11+A_R21*BYR21+A(20)*HYIMF

      BBZ=A(1)*BZCF+TAMP1*BZT1+TAMP2*BZT2+A_SRC*BZSRC+A_PRC*BZPRC
     * +A_R11*BZR11+A_R21*BZR21+A(20)*HZIMF
C
C   AND WE HAVE THE TOTAL EXTERNAL FIELD.
C
C
C  NOW, LET US CHECK WHETHER WE HAVE THE CASE (1). IF YES - ALL DONE:
C
      IF (SIGMA.LT.S0-DSIG) THEN    !  (X,Y,Z) IS INSIDE THE MAGNETOSPHERE

       BX=BBX
       BY=BBY
       BZ=BBZ
                     ELSE           !  THIS IS THE MOST COMPLEX CASE: WE ARE INSIDE
C                                             THE INTERPOLATION REGION
       FINT=0.5*(1.-(SIGMA-S0)/DSIG)
       FEXT=0.5*(1.+(SIGMA-S0)/DSIG)
C
       CALL DIPOLE_04 (PS,X,Y,Z,QX,QY,QZ)
       BX=(BBX+QX)*FINT+OIMFX*FEXT -QX
       BY=(BBY+QY)*FINT+OIMFY*FEXT -QY
       BZ=(BBZ+QZ)*FINT+OIMFZ*FEXT -QZ
c
        ENDIF  !   THE CASES (1) AND (2) ARE EXHAUSTED; THE ONLY REMAINING
C                      POSSIBILITY IS NOW THE CASE (3):
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ELSE
                CALL DIPOLE_04 (PS,X,Y,Z,QX,QY,QZ)
                BX=OIMFX-QX
                BY=OIMFY-QY
                BZ=OIMFZ-QZ
        ENDIF
C
      END
c
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
      SUBROUTINE BIRK_TOT_04 (IOPB,PS,X,Y,Z,BX11,BY11,BZ11,BX12,BY12,
     *                          BZ12,BX21,BY21,BZ21,BX22,BY22,BZ22)
C
C      IOPB -  BIRKELAND FIELD MODE FLAG:
C         IOPB=0 - ALL COMPONENTS
C         IOPB=1 - REGION 1, MODES 1 & 2
C         IOPB=2 - REGION 2, MODES 1 & 2
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION SH11(86),SH12(86),SH21(86),SH22(86)
      COMMON/BIRKPAR/ XKAPPA1,XKAPPA2   !  INPUT PARAMETERS, SPECIFIED FROM A MAIN PROGRAM
      COMMON/DPHI_B_RHO0/ DPHI,B,RHO_0,XKAPPA ! PARAMETERS, CONTROL DAY-NIGHT ASYMMETRY OF F.A.C.

      DATA SH11/46488.84663,-15541.95244,-23210.09824,-32625.03856,
     *-109894.4551,-71415.32808,58168.94612,55564.87578,-22890.60626,
     *-6056.763968,5091.368100,239.7001538,-13899.49253,4648.016991,
     *6971.310672,9699.351891,32633.34599,21028.48811,-17395.96190,
     *-16461.11037,7447.621471,2528.844345,-1934.094784,-588.3108359,
     *-32588.88216,10894.11453,16238.25044,22925.60557,77251.11274,
     *50375.97787,-40763.78048,-39088.60660,15546.53559,3559.617561,
     *-3187.730438,309.1487975,88.22153914,-243.0721938,-63.63543051,
     *191.1109142,69.94451996,-187.9539415,-49.89923833,104.0902848,
     *-120.2459738,253.5572433,89.25456949,-205.6516252,-44.93654156,
     *124.7026309,32.53005523,-98.85321751,-36.51904756,98.88241690,
     *24.88493459,-55.04058524,61.14493565,-128.4224895,-45.35023460,
     *105.0548704,-43.66748755,119.3284161,31.38442798,-92.87946767,
     *-33.52716686,89.98992001,25.87341323,-48.86305045,59.69362881,
     *-126.5353789,-44.39474251,101.5196856,59.41537992,41.18892281,
     *80.86101200,3.066809418,7.893523804,30.56212082,10.36861082,
     *8.222335945,19.97575641,2.050148531,4.992657093,2.300564232,
     *.2256245602,-.05841594319/

      DATA SH12/210260.4816,-1443587.401,-1468919.281,281939.2993,
     *-1131124.839,729331.7943,2573541.307,304616.7457,468887.5847,
     *181554.7517,-1300722.650,-257012.8601,645888.8041,-2048126.412,
     *-2529093.041,571093.7972,-2115508.353,1122035.951,4489168.802,
     *75234.22743,823905.6909,147926.6121,-2276322.876,-155528.5992,
     *-858076.2979,3474422.388,3986279.931,-834613.9747,3250625.781,
     *-1818680.377,-7040468.986,-414359.6073,-1295117.666,-346320.6487,
     *3565527.409,430091.9496,-.1565573462,7.377619826,.4115646037,
     *-6.146078880,3.808028815,-.5232034932,1.454841807,-12.32274869,
     *-4.466974237,-2.941184626,-.6172620658,12.64613490,1.494922012,
     *-21.35489898,-1.652256960,16.81799898,-1.404079922,-24.09369677,
     *-10.99900839,45.94237820,2.248579894,31.91234041,7.575026816,
     *-45.80833339,-1.507664976,14.60016998,1.348516288,-11.05980247,
     *-5.402866968,31.69094514,12.28261196,-37.55354174,4.155626879,
     *-33.70159657,-8.437907434,36.22672602,145.0262164,70.73187036,
     *85.51110098,21.47490989,24.34554406,31.34405345,4.655207476,
     *5.747889264,7.802304187,1.844169801,4.867254550,2.941393119,
     *.1379899178,.06607020029/

      DATA SH21/162294.6224,503885.1125,-27057.67122,-531450.1339,
     *84747.05678,-237142.1712,84133.61490,259530.0402,69196.05160,
     *-189093.5264,-19278.55134,195724.5034,-263082.6367,-818899.6923,
     *43061.10073,863506.6932,-139707.9428,389984.8850,-135167.5555,
     *-426286.9206,-109504.0387,295258.3531,30415.07087,-305502.9405,
     *100785.3400,315010.9567,-15999.50673,-332052.2548,54964.34639,
     *-152808.3750,51024.67566,166720.0603,40389.67945,-106257.7272,
     *-11126.14442,109876.2047,2.978695024,558.6019011,2.685592939,
     *-338.0004730,-81.99724090,-444.1102659,89.44617716,212.0849592,
     *-32.58562625,-982.7336105,-35.10860935,567.8931751,-1.917212423,
     *-260.2023543,-1.023821735,157.5533477,23.00200055,232.0603673,
     *-36.79100036,-111.9110936,18.05429984,447.0481000,15.10187415,
     *-258.7297813,-1.032340149,-298.6402478,-1.676201415,180.5856487,
     *64.52313024,209.0160857,-53.85574010,-98.52164290,14.35891214,
     *536.7666279,20.09318806,-309.7349530,58.54144539,67.45226850,
     *97.92374406,4.752449760,10.46824379,32.91856110,12.05124381,
     *9.962933904,15.91258637,1.804233877,6.578149088,2.515223491,
     *.1930034238,-.02261109942/

      DATA SH22/-131287.8986,-631927.6885,-318797.4173,616785.8782,
     *-50027.36189,863099.9833,47680.20240,-1053367.944,-501120.3811,
     *-174400.9476,222328.6873,333551.7374,-389338.7841,-1995527.467,
     *-982971.3024,1960434.268,297239.7137,2676525.168,-147113.4775,
     *-3358059.979,-2106979.191,-462827.1322,1017607.960,1039018.475,
     *520266.9296,2627427.473,1301981.763,-2577171.706,-238071.9956,
     *-3539781.111,94628.16420,4411304.724,2598205.733,637504.9351,
     *-1234794.298,-1372562.403,-2.646186796,-31.10055575,2.295799273,
     *19.20203279,30.01931202,-302.1028550,-14.78310655,162.1561899,
     *.4943938056,176.8089129,-.2444921680,-100.6148929,9.172262228,
     *137.4303440,-8.451613443,-84.20684224,-167.3354083,1321.830393,
     *76.89928813,-705.7586223,18.28186732,-770.1665162,-9.084224422,
     *436.3368157,-6.374255638,-107.2730177,6.080451222,65.53843753,
     *143.2872994,-1028.009017,-64.22739330,547.8536586,-20.58928632,
     *597.3893669,10.17964133,-337.7800252,159.3532209,76.34445954,
     *84.74398828,12.76722651,27.63870691,32.69873634,5.145153451,
     *6.310949163,6.996159733,1.971629939,4.436299219,2.904964304,
     *.1486276863,.06859991529/

      XKAPPA=XKAPPA1        !  FORWARDED IN BIRK_1N2
      X_SC=XKAPPA1-1.1D0    !  FORWARDED IN BIRK_SHL

      IF (IOPB.EQ.0.OR.IOPB.EQ.1) THEN

      CALL BIRK_1N2_04 (1,1,PS,X,Y,Z,FX11,FY11,FZ11)           !  REGION 1, MODE 1
      CALL BIRK_SHL (SH11,PS,X_SC,X,Y,Z,HX11,HY11,HZ11)
      BX11=FX11+HX11
      BY11=FY11+HY11
      BZ11=FZ11+HZ11

      CALL BIRK_1N2_04 (1,2,PS,X,Y,Z,FX12,FY12,FZ12)           !  REGION 1, MODE 2
      CALL BIRK_SHL (SH12,PS,X_SC,X,Y,Z,HX12,HY12,HZ12)
      BX12=FX12+HX12
      BY12=FY12+HY12
      BZ12=FZ12+HZ12

      ENDIF

      XKAPPA=XKAPPA2        !  FORWARDED IN BIRK_1N2
      X_SC=XKAPPA2-1.0D0    !  FORWARDED IN BIRK_SHL

      IF (IOPB.EQ.0.OR.IOPB.EQ.2) THEN

      CALL BIRK_1N2_04 (2,1,PS,X,Y,Z,FX21,FY21,FZ21)           !  REGION 2, MODE 1
      CALL BIRK_SHL (SH21,PS,X_SC,X,Y,Z,HX21,HY21,HZ21)
      BX21=FX21+HX21
      BY21=FY21+HY21
      BZ21=FZ21+HZ21

      CALL BIRK_1N2_04 (2,2,PS,X,Y,Z,FX22,FY22,FZ22)           !  REGION 2, MODE 2
      CALL BIRK_SHL (SH22,PS,X_SC,X,Y,Z,HX22,HY22,HZ22)
      BX22=FX22+HX22
      BY22=FY22+HY22
      BZ22=FZ22+HZ22

      ENDIF

      RETURN
      END
C
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c
      SUBROUTINE BIRK_1N2_04 (NUMB,MODE,PS,X,Y,Z,BX,BY,BZ)
C
C  CALCULATES COMPONENTS  OF REGION 1/2 FIELD IN SPHERICAL COORDS.  DERIVED FROM THE S/R DIPDEF2C (WHICH
C    DOES THE SAME JOB, BUT INPUT/OUTPUT THERE WAS IN SPHERICAL COORDS, WHILE HERE WE USE CARTESIAN ONES)
C
C   INPUT:  NUMB=1 (2) FOR REGION 1 (2) CURRENTS
C           MODE=1 YIELDS SIMPLE SINUSOIDAL MLT VARIATION, WITH MAXIMUM CURRENT AT DAWN/DUSK MERIDIAN
C     WHILE MODE=2 YIELDS THE SECOND HARMONIC.
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A11(31),A12(31),A21(31),A22(31)
      COMMON/MODENUM/ M
      COMMON/DTHETA/ DTHETA

      COMMON/DPHI_B_RHO0/ DPHI,B,RHO_0,XKAPPA ! THESE PARAMETERS CONTROL DAY-NIGHT ASYMMETRY OF F.A.C., AS FOLLOWS:

C  (1) DPHI:   HALF-DIFFERENCE (IN RADIANS) BETWEEN DAY AND NIGHT LATITUDE OF FAC OVAL AT IONOSPHERIC ALTITUDE;
C              TYPICAL VALUE: 0.06
C  (2) B:      AN ASYMMETRY FACTOR AT HIGH-ALTITUDES;  FOR B=0, THE ONLY ASYMMETRY IS THAT FROM DPHI
C              TYPICAL VALUES: 0.35-0.70
C  (3) RHO_0:  A FIXED PARAMETER, DEFINING THE DISTANCE RHO, AT WHICH THE LATITUDE SHIFT GRADUALLY SATURATES AND
C              STOPS INCREASING
C              ITS VALUE WAS ASSUMED FIXED, EQUAL TO 7.0.
C  (4) XKAPPA: AN OVERALL SCALING FACTOR, WHICH CAN BE USED FOR CHANGING THE SIZE OF THE F.A.C. OVAL
C
      DATA BETA,RH,EPS/0.9D0,10.D0,3.D0/ ! parameters of the tilt-dependent deformation of the untilted F.A.C. field

      DATA A11/.1618068350,-.1797957553,2.999642482,-.9322708978,
     *-.6811059760,.2099057262,-8.358815746,-14.86033550,.3838362986,
     *-16.30945494,4.537022847,2.685836007,27.97833029,6.330871059,
     *1.876532361,18.95619213,.9651528100,.4217195118,-.08957770020,
     *-1.823555887,.7457045438,-.5785916524,-1.010200918,.01112389357,
     *.09572927448,-.3599292276,8.713700514,.9763932955,3.834602998,
     *2.492118385,.7113544659/
      DATA A12/.7058026940,-.2845938535,5.715471266,-2.472820880,
     *-.7738802408,.3478293930,-11.37653694,-38.64768867,.6932927651,
     *-212.4017288,4.944204937,3.071270411,33.05882281,7.387533799,
     *2.366769108,79.22572682,.6154290178,.5592050551,-.1796585105,
     *-1.654932210,.7309108776,-.4926292779,-1.130266095,-.009613974555,
     *.1484586169,-.2215347198,7.883592948,.02768251655,2.950280953,
     *1.212634762,.5567714182/
      DATA A21/.1278764024,-.2320034273,1.805623266,-32.37241440,
     *-.9931490648,.3175085630,-2.492465814,-16.21600096,.2695393416,
     *-6.752691265,3.971794901,14.54477563,41.10158386,7.912889730,
     *1.258297372,9.583547721,1.014141963,.5104134759,-.1790430468,
     *-1.756358428,.7561986717,-.6775248254,-.04014016420,.01446794851,
     *.1200521731,-.2203584559,4.508963850,.8221623576,1.779933730,
     *1.102649543,.8867880020/
      DATA A22/.4036015198,-.3302974212,2.827730930,-45.44405830,
     *-1.611103927,.4927112073,-.003258457559,-49.59014949,.3796217108,
     *-233.7884098,4.312666980,18.05051709,28.95320323,11.09948019,
     *.7471649558,67.10246193,.5667096597,.6468519751,-.1560665317,
     *-1.460805289,.7719653528,-.6658988668,.2515179349E-05,
     *.02426021891,.1195003324,-.2625739255,4.377172556,.2421190547,
     *2.503482679,1.071587299,.7247997430/

      B=0.5
      RHO_0=7.0

      M=MODE
      IF (NUMB.EQ.1) THEN
          DPHI=0.055D0
          DTHETA=0.06D0
      ENDIF

      IF (NUMB.EQ.2) THEN
          DPHI=0.030D0
          DTHETA=0.09D0
      ENDIF

      Xsc=X*XKAPPA
      Ysc=Y*XKAPPA
      Zsc=Z*XKAPPA
      RHO=DSQRT(Xsc**2+Zsc**2)

      Rsc=DSQRT(Xsc**2+Ysc**2+Zsc**2)                                 !  SCALED
      RHO2=RHO_0**2

      IF (Xsc.EQ.0.D0.AND.Zsc.EQ.0.D0) THEN
         PHI=0.D0
      ELSE
         PHI=DATAN2(-Zsc,Xsc)  !  FROM CARTESIAN TO CYLINDRICAL (RHO,PHI,Y)
      ENDIF

      SPHIC=DSIN(PHI)
      CPHIC=DCOS(PHI)  !  "C" means "CYLINDRICAL", TO DISTINGUISH FROM SPHERICAL PHI

      BRACK=DPHI+B*RHO2/(RHO2+1.D0)*(RHO**2-1.D0)/(RHO2+RHO**2)
      R1RH=(Rsc-1.D0)/RH
      PSIAS=BETA*PS/(1.D0+R1RH**EPS)**(1.D0/EPS)

      PHIS=PHI-BRACK*DSIN(PHI) -PSIAS
      DPHISPHI=1.D0-BRACK*DCOS(PHI)
      DPHISRHO=-2.D0*B*RHO2*RHO/(RHO2+RHO**2)**2 *DSIN(PHI)
     *   +BETA*PS*R1RH**(EPS-1.D0)*RHO/(RH*Rsc*
     *   (1.D0+R1RH**EPS)**(1.D0/EPS+1.D0))
      DPHISDY= BETA*PS*R1RH**(EPS-1.D0)*Ysc/(RH*Rsc*
     *   (1.D0+R1RH**EPS)**(1.D0/EPS+1.D0))

      SPHICS=DSIN(PHIS)
      CPHICS=DCOS(PHIS)

      XS= RHO*CPHICS
      ZS=-RHO*SPHICS

      IF (NUMB.EQ.1) THEN
        IF (MODE.EQ.1) CALL TWOCONES (A11,XS,Ysc,ZS,BXS,BYAS,BZS)
        IF (MODE.EQ.2) CALL TWOCONES (A12,XS,Ysc,ZS,BXS,BYAS,BZS)
      ELSE
        IF (MODE.EQ.1) CALL TWOCONES (A21,XS,Ysc,ZS,BXS,BYAS,BZS)
        IF (MODE.EQ.2) CALL TWOCONES (A22,XS,Ysc,ZS,BXS,BYAS,BZS)
      ENDIF

      BRHOAS=BXS*CPHICS-BZS*SPHICS
      BPHIAS=-BXS*SPHICS-BZS*CPHICS

      BRHO_S=BRHOAS*DPHISPHI                             *XKAPPA        ! SCALING
      BPHI_S=(BPHIAS-RHO*(BYAS*DPHISDY+BRHOAS*DPHISRHO)) *XKAPPA
      BY_S=BYAS*DPHISPHI                                 *XKAPPA

      BX=BRHO_S*CPHIC-BPHI_S*SPHIC
      BY=BY_S
      BZ=-BRHO_S*SPHIC-BPHI_S*CPHIC

      RETURN
      END
c
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C
      DOUBLE PRECISION FUNCTION APPRC(R,SINT,COST)
C
C      Calculates azimuthal component of the vector potential of the symmetric
c  part of the model PARTIAL ring current.
C
      IMPLICIT  REAL * 8  (A - H, O - Z)
      LOGICAL PROX
      DATA A1,A2,RRC1,DD1,RRC2,DD2,P1,ALPHA1,DAL1,BETA1,DG1,P2,ALPHA2,
     * DAL2,BETA2,DG2,BETA3,P3,ALPHA3,DAL3,BETA4,DG3,BETA5,Q0,Q1,ALPHA4,
     * DAL4,DG4,Q2,ALPHA5,DAL5,DG5,BETA6,BETA7
     * /-80.11202281,12.58246758,6.560486035,1.930711037,3.827208119,
     *.7789990504,.3058309043,.1817139853,.1257532909,3.422509402,
     *.04742939676,-4.800458958,-.02845643596,.2188114228,2.545944574,
     *.00813272793,.35868244,103.1601001,-.00764731187,.1046487459,
     *2.958863546,.01172314188,.4382872938,.01134908150,14.51339943,
     *.2647095287,.07091230197,.01512963586,6.861329631,.1677400816,
     *.04433648846,.05553741389,.7665599464,.7277854652/

      PROX=.FALSE.
      SINT1=SINT
      COST1=COST
      IF (SINT1.LT.1.D-2) THEN  !  TOO CLOSE TO Z-AXIS;  USE LINEAR INTERPOLATION BETWEEN SINT=0 & SINT=0.01
        SINT1=1.D-2
        COST1=.99994999875
        PROX=.TRUE.
      ENDIF

         ALPHA=SINT1**2/R         !  R,THETA -> ALPHA,GAMMA
         GAMMA=COST1/R**2

         ARG1=-(GAMMA/DG1)**2
         ARG2=-((ALPHA-ALPHA4)/DAL4)**2-(GAMMA/DG4)**2

         IF (ARG1.LT.-500.D0) THEN        !   TO PREVENT "FLOATING UNDERFLOW" CRASHES
           DEXP1=0.D0
         ELSE
           DEXP1=DEXP(ARG1)
         ENDIF

         IF (ARG2.LT.-500.D0) THEN        !   TO PREVENT "FLOATING UNDERFLOW" CRASHES
           DEXP2=0.D0
         ELSE
           DEXP2=DEXP(ARG2)
         ENDIF

         ALPHA_S=ALPHA*(1.D0+P1/(1.D0+((ALPHA-ALPHA1)/DAL1)**2)**BETA1
     * *DEXP1+P2*(ALPHA-ALPHA2)/(1.D0+((ALPHA-ALPHA2)/DAL2)**2)**BETA2
     */(1.D0+(GAMMA/DG2)**2)**BETA3
     *+P3*(ALPHA-ALPHA3)**2/(1.D0+((ALPHA-ALPHA3)/DAL3)**2)**BETA4
     */(1.D0+(GAMMA/DG3)**2)**BETA5)     !  ALPHA -> ALPHA_S  (DEFORMED)

         GAMMA_S=GAMMA*(1.D0+Q0+Q1*(ALPHA-ALPHA4)*DEXP2              !  GAMMA -> GAMMA_  (DEFORMED)
     * +Q2*(ALPHA-ALPHA5)/(1.D0+((ALPHA-ALPHA5)/DAL5)**2)**BETA6
     * /(1.D0+(GAMMA/DG5)**2)**BETA7)

         GAMMAS2=GAMMA_S**2

         ALSQH=ALPHA_S**2/2.D0                            !  ALPHA_S,GAMMA_S -> RS,SINTS,COSTS
         F=64.D0/27.D0*GAMMAS2+ALSQH**2
         Q=(DSQRT(F)+ALSQH)**(1.D0/3.D0)
         C=Q-4.D0*GAMMAS2**(1.D0/3.D0)/(3.D0*Q)
         IF (C.LT.0.D0) C=0.D0
         G=DSQRT(C**2+4.D0*GAMMAS2**(1.D0/3.D0))
         RS=4.D0/((DSQRT(2.D0*G-C)+DSQRT(C))*(G+C))
         COSTS=GAMMA_S*RS**2
         SINTS=DSQRT(1.D0-COSTS**2)
         RHOS=RS*SINTS
         RHOS2=RHOS**2
         ZS=RS*COSTS
C
c  1st loop:

         P=(RRC1+RHOS)**2+ZS**2+DD1**2
         XK2=4.D0*RRC1*RHOS/P
         XK=SQRT(XK2)
         XKRHO12=XK*SQRT(RHOS)
C
      XK2S=1.D0-XK2
      DL=DLOG(1.D0/XK2S)
      ELK=1.38629436112d0+XK2S*(0.09666344259D0+XK2S*(0.03590092383+
     *     XK2S*(0.03742563713+XK2S*0.01451196212))) +DL*
     *     (0.5D0+XK2S*(0.12498593597D0+XK2S*(0.06880248576D0+
     *     XK2S*(0.03328355346D0+XK2S*0.00441787012D0))))
      ELE=1.D0+XK2S*(0.44325141463D0+XK2S*(0.0626060122D0+XK2S*
     *      (0.04757383546D0+XK2S*0.01736506451D0))) +DL*
     *     XK2S*(0.2499836831D0+XK2S*(0.09200180037D0+XK2S*
     *       (0.04069697526D0+XK2S*0.00526449639D0)))
C
      APHI1=((1.D0-XK2*0.5D0)*ELK-ELE)/XKRHO12
c
c  2nd loop:

         P=(RRC2+RHOS)**2+ZS**2+DD2**2
         XK2=4.D0*RRC2*RHOS/P
         XK=SQRT(XK2)
         XKRHO12=XK*SQRT(RHOS)
C
      XK2S=1.D0-XK2
      DL=DLOG(1.D0/XK2S)
      ELK=1.38629436112d0+XK2S*(0.09666344259D0+XK2S*(0.03590092383+
     *     XK2S*(0.03742563713+XK2S*0.01451196212))) +DL*
     *     (0.5D0+XK2S*(0.12498593597D0+XK2S*(0.06880248576D0+
     *     XK2S*(0.03328355346D0+XK2S*0.00441787012D0))))
      ELE=1.D0+XK2S*(0.44325141463D0+XK2S*(0.0626060122D0+XK2S*
     *      (0.04757383546D0+XK2S*0.01736506451D0))) +DL*
     *     XK2S*(0.2499836831D0+XK2S*(0.09200180037D0+XK2S*
     *       (0.04069697526D0+XK2S*0.00526449639D0)))
C
      APHI2=((1.D0-XK2*0.5D0)*ELK-ELE)/XKRHO12

      APPRC=A1*APHI1+A2*APHI2
      IF (PROX) APPRC=APPRC*SINT/SINT1   !   LINEAR INTERPOLATION, IF TOO CLOSE TO THE Z-AXIS
C
      RETURN
      END
C
c===========================================================================
c
       SUBROUTINE DIPOLE_04 (PS,X,Y,Z,BX,BY,BZ)
C
C      A DOUBLE PRECISION ROUTINE
C
C  CALCULATES GSM COMPONENTS OF A GEODIPOLE FIELD WITH THE DIPOLE MOMENT
C  CORRESPONDING TO THE EPOCH OF 2000.
C
C----INPUT PARAMETERS:
C     PS - GEODIPOLE TILT ANGLE IN RADIANS,
C     X,Y,Z - GSM COORDINATES IN RE (1 RE = 6371.2 km)
C
C----OUTPUT PARAMETERS:
C     BX,BY,BZ - FIELD COMPONENTS IN GSM SYSTEM, IN NANOTESLA.
C
      IMPLICIT REAL*8 (A-H,O-Z)
      SPS=DSIN(PS)
      CPS=DCOS(PS)
      P=X**2
      U=Z**2
      V=3.D0*Z*X
      T=Y**2
      Q=30115.D0/DSQRT(P+T+U)**5
      BX=Q*((T+U-2.D0*P)*SPS-V*CPS)
      BY=-3.D0*Y*Q*(X*SPS+Z*CPS)
      BZ=Q*((P+T-2.D0*U)*CPS-V*SPS)
      RETURN
      END



C
C(((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((



***************************************************************************
*                                					   
*									   
*                     DYNAMIC PARABOLOID MODEL OF			   
*                      THE EARTH'S MAGNETOSPHERE			   
*                           (March 2000)				   
*	              (Revision on 19 June 2002)			   
*									   
* The dynamic paraboloid model  allows to calculate the 		   
* magnetic field inside the Earth's magnetosphere. In the base of the      
* model lies CPMOD (closed paraboloid model) software. Dynamic paraboloid  
* model enable calculation the time variations of each large scale source  
* of magnetospheric magnetic field by empirical data. The possibilities    
* are provided to "switch on" and "switch off" separate sources of the     
* magnetic field, and to change the parameters defining their intensity.   
*									   
* Contacts: Igor Alexeev						   
*           Scobeltsyn Inst. of Nuclear Physics 			   
*           Moscow State University,            			   
*           Moscow, 119899, Russia					   
*           alexeev@dec1.sinp.msu.ru					   
***************************************************************************

      subroutine a2000(ut,iy,mo,id,ro,v,bimf,dst,al,x,bm,bb)
*--------------------------------------------------------------------------
*  Calculation of magnetic field vector Bm in point x(3) 		   
*              inside the Earth's magnetosphere 			   
*  INPUT: 								   
*           time moment UT (Universal Time) by				   
*           year=iy;							   
*           month=mo;							   
*           day in month=id.						   
*           ro, V, bz - solar wind density and velocity, IMF Bz.	   
*           dst - value of Dst index;					   
*           AL  - value 0f al index. 					   
*  OUTPUT: magnetic field components at the point x0(3) 		   
*                      in GSM coordinates, nT:				   
*         bm(i) - total magnetic field (i=1,3); 			   
*         bb(1,i) - geomagnetic dipole magnetic field;			   
*         bb(2,i) - ring current magnetic field;			   
*         bb(3,i) - geomagnetic tail currents magnetic field;		   
*         bb(4,i) - magnetic field of CF currents shielding dipole;	   
*         bb(5,i) - magnetic field of CF currents shielding ring current;  
*         bb(6,i) - magnetic field of Region 1 FAC;     		   
*         bb(7,i) - IMF penetrated into the magnetosphere.      	   
*									   
* WARNING: Because of the paraboloid coordinates singularity, avoid	   
*          the magnetic field calculations at the Ox axis.		   
* 									   
*          Written by I. Alexeev					   
*--------------------------------------------------------------------------

      IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION x(3), bm(3), par(10), bb(7,3)
      call submod(ut,iy,mo,id,ro,v,bimf,dst,al,par)
      r1=par(6)
       IF(x(1)+0.5/R1*(x(2)**2+x(3)**2).GT.R1) then
       do i=1,3
       bm(i)=0.
       do j=1,7
       bb(j,i)=0.
       end do
       end do
       return
       end if
      call a_field(x,par,bm,bb)
c      do 4 i=1,3
c4     bm(i)=bm(i)-bb(1,i)
      return
      END
C
      subroutine submod(ut,iy,mo,id,ro,v,bimf,dst,al,par)
*---------------------------------------------------------------------*
*  Calculation of the paraboloid model input parameters 	      *
*  by empirical data						      *
*  INPUT (empirical data):					      *
*           time moment UT (Universal Time)			      *
*           year=iy;						      *
*           month=mo;						      *
*           day in month=id.					      *
*           ro, V, bimf- solar wind density and velocity, IMF.	      *
*           dst - value of Dst index;				      *
*           AL  - value of al index. 				      *
*  OUTPUT (model input parameters): 				      *
*   	      par(1) - geomagnetic dipole tilt angle, degrees;        *
*   	      par(2) - dipole magnetic field at equator, nT;	      *
*   	      par(3) - magnetic flux through the tail lobes, Wb;      *
*   	      par(4) - maximum ring current intensity, nT;	      *
*   	      par(5) - the total current of Region 1 FAC, MA;	      *
*   	      par(6) - magnetopause stand-off distance, Re;	      *
*   	      par(7) - distance to the inner edge of geotail	      *
*   		       current sheet;				      *
*   	      par(8-10) -IMF penetrated components in GSM coord., nT. *
*								      *
*    Written by V. Kalegaev					      *
*---------------------------------------------------------------------*
      IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION x(3), bm(3), bimf(3), par(10), bb(7,3)
      iday=IDD(iy,mo,id)                ! the day number in year
*
***>  calculation of the tilt angle, tpsi [degrees],
***>  and the dipole magnetic field at the Earth's equator, BD [nT] 
*
      call trans(ut,iday,iy,tpsi,bd)

*
***>  calculation of the magnetopause stand-off distance, R1 [Re], 
*     by Shue et al. 1997
*
      Psw=(1.67e-6)*ro*V*V              ! Solar wind pressure [nT]          
      bz=bimf(3)
      if (bz.ge.0.) then
      r1=(11.4+0.013*bz)*(psw**(-1./6.6))
      else
      r1=(11.4+0.14*bz)*(psw**(-1./6.6))
      end if
*
***>  calculation of the distance to the geotail inner edge: R2 [Re]
*
      if(dst.ge.-10.0)then
      R2=0.7*r1
      goto 3
      end if
      fie=74.9-8.6*log10(-dst)
      fi=fie*3.1416/180.0
      R2=1.0/cos(fi)/cos(fi)
 3     al0=sqrt(1.+2*R2/R1)
*
***>  calculation of the geotail lobe magnetic flux: FLUX [Wb] 
*
      BT=-AL/7                                                                  
      FLUX=1.5708*R1*R1*AL0*BT*6.37816*6.37816*1.E3                            
      flux=flux+395.98278e6                                                      
*
***>  calculation of the ring current mag. f. in the Earth's
*     center: BR [nT] 
*
      BR=dst-10.
      if(dst.ge.0.0)br=-10.0

*
***>  calculation of the total FAC intensity: AJ0 [MA] 
*
      AJ0=0.327744
      if(bz.le.-1.61133) AJ0=-1.017*bz/5.
      AJ0=2.*AJ0*sqrt(400./v)*((5/ro)**0.125)

       par(1)=tpsi
       par(2)=BD
       par(3)=FLUX
       par(4)=BR
       par(5)=AJ0
       par(6)=R1
       par(7)=R2
      do 2 i=1,3
2     par(i+7)=bimf(i)
      return
      END
C
      SUBROUTINE TRANS (UT,IDAY,iyear,tpsi,BD)
***************************************************************
*   Calculation of the geomagnetic dipole tilt angle and      .
*   matrices of transition from the geographic coordinates    .
*   (Cartesian) to the GSM-coordinates (G2GSM(3,3) in the     .
*   COMMON BLOCK /TRAN/).                                     .
*                _            _                               .
*                Xgsm = G2GSM*Xgeogr                          .
*                                                             .
*    ALPHA1 is the angle between the Earth's axis and the     .
*           dipole moment;                                    .
*    ALPHA2 is the angle between the Earth's axis and         .
*           the normal to the ecliptic plane;                 .
*    PHIM   is the angle between the midnight meridian plane  .
*           and the meridional plane;                         .
*    PHISE  is the angle between the Earth-Sun line and       .
*           the projection of the Earth's axis on the         .
*           ecliptic plane;                                   .
*    PSI    is the tilt angle of the geomagnetic dipole;      .
*    B      is the Sun's deflection;                          .
*    UT     is universal time;                                .
*    IDAY   is the day number;                                .
*    B1     is western longitude of the noon meridian;        .
*    B2     is the angle between the noon geogr. meridian     .
*           and the Y=0 plane in the GSM-coordinates.         .
*       Written by I. Alexeev                                 .
*							      .
*   Acknowledgment:                                           .  
*       Program SUN written by Gilbert D. Mead is used	      .
*       to determine some parameters dependent on the Earth   .
*       and Sun mutual position (SLONG, sind, and cosd from   .
*       COMMON BLOCK /ddd/).                                  .
***************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)

      COMMON/TRAN/g2gsm(3,3)
      COMMON/ddd/ sind,cosd
      dimension gauss(3)         
      Ihour=int(ut)
      dmin=(ut-ihour)*60.
      min=int(dmin)
      isec=int((dmin-min)*60.)
      call ALX00_SUN(IYeaR,IDAY,IHOUR,MIN,ISEC,GST,SLONG,
     +               SRASN,SDEC)
      P= 3.1415926/180.
      pi2=3.1415926/2
          call dipgarm(iyear, gauss)
          G1=GAUSS(2)                                                           
          H1=GAUSS(3)                                                           
          PD=G1**2+H1**2                                                        
          G10=GAUSS(1)                                                           
          BD=-SQRT(G10*G10+PD)                                                    
          ALPHA1= ATAN(-1.*SQRT(PD)/G10)                                         
          PHINP=  ATAN(H1/G1)                                                   
        ALPHA2= 23.4419*P
        PHIM= P*(UT*15+phinp/p)
        PHISE= pi2-slong+9.924E-5
         
        B1=P*15.*(UT-12)
        SB= SIN(ALPHA2)*COS(PHISE)
        CB= SQRT(1-SB*SB)
        SB1=SIN(B1)
        CB1=COS(B1)
          SPSI= -SB*COS(ALPHA1) + CB*SIN(ALPHA1)*COS(PHIM)
          CPSI=  SQRT(1-SPSI*SPSI)
          psi=asin(spsi)/p
        CB2=(COS(ALPHA1)+SB*SPSI)/CPSI/CB
        SB2=SQRT(abs(1-CB2*CB2))
        IF(PHIM.LE.0..OR.PHIM.GE.3.1415926) SB2=-SB2
          tpsi=psi
      g2gsm(1,1)=cb1*cb
      g2gsm(1,2)=-sb1*cb
      g2gsm(1,3)=sb
      g2gsm(2,1)=sb1*cb2-cb1*sb*sb2
      g2gsm(2,2)=cb1*cb2+sb1*sb*sb2
      g2gsm(2,3)=cb*sb2
      g2gsm(3,1)=-sb1*sb2-cb1*sb*cb2
      g2gsm(3,2)=-cb1*sb2+sb1*sb*cb2
      g2gsm(3,3)=cb*cb2
         RETURN
         END
  
      FUNCTION IDD(iy,mo,id)
********************************************************************
*  Calculation of the day number in a year			   *
*  INPUT PARAMETERS: year (IY), month (MO), day in the month (ID)  *
*  Written by V. Kalegaev					   *
********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)

      II=0
      IF (MO.EQ.1) GOTO 4
 5    DO 10 M=1,MO-1
      GOTO (1,2,1,3,1,3,1,1,3,1,3) M
 1    II=II+31
      GOTO 10
 2    II=II+28
      GOTO 10
 3    II=II+30
 10   CONTINUE
 4    II=II+ID
      if (mod(iy,100).eq.0.and.mod(iy,400).ne.0) goto 6
      IF (MOD(IY,4).EQ.0.AND.II.GT.59.AND.MO.GT.2) II=II+1
 6    IDD=II
      RETURN
      END

C     
      subroutine dipgarm(iyear,gauss)
*------------------------------------------------------------------------*
*   Calculation of the first three Gaussian coefficients for given year  *
*   year>=1900  							 *
*   Last IGRF values are for 1995.					 *
*   Written by V. Kalegaev						 *
*------------------------------------------------------------------------*
      IMPLICIT REAL*8 (A-H,O-Z)

      dimension gauss(3), g(60), sg95(3)
      data g/
     *-31543., -2298.,  5922.,
     *-31464., -2298.,  5909.,
     *-31354., -2297.,  5898.,
     *-31212., -2306.,  5875.,
     *-31060., -2317.,  5845.,
     *-30926., -2318.,  5817.,
     *-30805., -2316.,  5808.,
     *-30715., -2306.,  5812.,
     *-30654., -2292.,  5821.,
     *-30594., -2285.,  5810.,
     *-30554., -2250.,  5815.,
     *-30500., -2215.,  5820.,
     *-30421., -2169.,  5791.,
     *-30334., -2119.,  5776.,
     *-30220., -2068.,  5737.,
     *-30100., -2013.,  5675.,
     *-29992., -1956.,  5604.,
     *-29873., -1905.,  5500.,
     *-29775., -1848.,  5406.,
     *-29682., -1789.,  5318./
      data sg95/17.6, 13.0, -18.3/
      i1=(iyear-1900)
      i2=i1/5
      i3=i1-5*i2
      n=(i2)*3+1
      if(i2.gt.18) goto 1
      do i=1,3
      a=g(n-1+i)
      b=g(n+2+i)
      gauss(i)=a+(b-a)/5.*i3
      end do
      return
 1    continue
      do i=1,3
      a=g(57+i)
      gauss(i)=a+sg95(i)*(iyear-1995)
      end do
      return
      end

      subroutine A_field(x0,par,bm,bdd)
******************************************************************************
*  Calculation of the magnetic field in the magnetosphere                    *
*  inthe point x0(3) by model input parameters par(10).	                     *
*      Ver. 3 on March 2000.                                                 *
*  (The model's parameters par(1-10) are moved to  COMMON/T2/ )              *
*									     *
*  INPUT PARAMETERS: x0(3) is a point where the magnetic field is being      *
*                          calculated, in GSM coordinates, Re;		     *
*                    par(1) - geomagnetic dipole tilt angle, degrees;	     *
*                    par(2) - dipole magnetic field at equator, nT;	     *
*                    par(3) - magnetic flux through the tail lobes, Wb;      *
*                    par(4) - maximum ring current intensity, nT;	     *
*                    par(5) - the total current of Region 1 FAC, MA;	     *
*                    par(6) - magnetopause stand-off distance, Re;	     *
*                    par(7) - distance to the inner edge of geotail 	     *
*                             current sheet;				     *
*                    par(8-10) -IMF penetrated components in GSM coord., nT. *
*                                                  			     *
*  OUTPUT - magnetic field components at the point x0(3) in GSM  coord., nT. *
*         bm(i) - total magnetic field (i=1,3); 			     *
*         bdd(1,i) - geomagnetic dipole magnetic field; 		     *
*         bdd(2,i) - ring current magnetic field;			     *
*         bdd(3,i) - geomagnetic tail currents magnetic field;		     *
*         bdd(4,i) - magnetic field of CF currents shielding dipole;	     *
*         bdd(5,i) - magnetic field of CF currents shielding ring current;   *
*         bdd(6,i) - magnetic field of Region 1 FAC;     		     *
*         bdd(7,i) - IMF penetrated into the magnetosphere.     	     *
*									     *
* WARNING: Because of the paraboloid coordinates singularity, avoid	     *
*          the magnetic field calculations at the Ox axis.		     *
*									     *
* Written by  V.Kalegaev	         				     *
******************************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)

      COMMON/TK/B1(3),B2(3),B1A(3),B1R(3),B2A(3),B2R(3),
     *BA(3),DB(3),AB,B(3)
      COMMON/BEGFC/B1CF(3),B1CFD(3),B1CFR(3),B1D(3),B1RC(3)
      COMMON/T2/PI,R1,BETA0,AL0,C0,E5,AL1,BT,CPSI,SPSI,PSI,Z0,B0,bd,bd0
      COMMON/IMFd/Bimf(3),b3(3)
      COMMON/fac12/ami1,ami2,tm1,tm2
      COMMON/bfac12/bfac1(3),bfac2(3)
      DIMENSION X0(3), FF(3), bm(3), par(10), bdd(7,3)

      psi=par(1)
      bd=par(2)
      flux=par(3) ! magnetic flux through the tail lobes
      BR=par(4)   ! maximum ring current intensity
      ami1=par(5) ! maximum of Region 1 FAC intensity
      R1=par(6)   ! magnetopause stand-off distance
      R2=par(7)   ! distance to the inner edge of geotail current sheet
      do 2 i=1,3
2     bimf(i)=par(i+7) ! IMF components
        call mas1d (flux,BR,R1,R2)
        call field(X0,FF)

         do 1 i=1,3
         bm(i)=b(i)         ! total magnetic field
         bdd(1,i)=b1d(i)    ! geomagnetic dipole contribution
         bdd(2,i)=b1rc(i)   ! ring current contribution
         bdd(3,i)=b2(i)     ! geomagnetic tail contribution
         bdd(4,i)=b1cfd(i)  ! contribution of CF currents shielding dip.
         bdd(5,i)=b1cfr(i)  ! contribution of CF currents shielding RC
         bdd(6,i)=bfac1(i)  !  contribution of Region 1 FAC. 
 1       bdd(7,i)=b3(i)     !  contribution of IMF. 
      return
      END


**********************************************************************
*                             C P M O D 			     *
*								     *
*                     CLOSED PARABOLOID MODEL OF		     *
*             MAGNETIC FIELD IN THE EARTH'S MAGNETOSPHERE	     *
*                     (ver. 3, February 2000)			     *
*								     *
* The CPMOD (closed paraboloid model) software was written based     *
* on the I.I.Alexeev paraboloid model A78. The paraboloid	     *
* model allows to calculate the magnetic field inside the Earth's    *
* magnetosphere. The magnetospheric magnetic field is described by   *
* sum of magnetic field of the following sources:                    *
* (1) the geomagnetic dipole, 	                                     *
* (2) the ring current;                                              *   
* (3) the current system of magnetotail including the dawn-dusk      *
* currents across the magnetotail current sheet and the closure	     *
* currents on the magnetopause;	                                     *
* magnetopause;                                                      *   
* (4) the Chapmen-Ferraro currents on the  magnetopause (the dipole  *
* magnetopause (the dipole screening field);                         *
* (5) the currents on the magnetopause screening the ring current;   *
* (6) Region 1 FAC;                                                  *
* (7) the IMF penetrated from magnetosheath.  			     *
*								     *
*  Written by I. Alexeev.					     *
**********************************************************************

      SUBROUTINE FIELD(UF,FF)                                           
***************************************************************************
* Program calculating the magnetic field in the magnetosphere.		  *
*									  *
* INPUT PARAMETERS:  x0(3) is a point where the magnetic field is being   *
*                          calculated, in GSM coordinates, Re;		  *
* OUTPUT PARAMETERS: ff(3) is a normalized vector of the magnetic	  *
*                          field at the point x0(3).			  *
* NOTE: The total magnetic field and the magnetic fields of the 	  *
*       magnetospheric current systems are stored in COMMON BLOCKs	  *
*       /TK/ and /BEGFC/ .						  *
* WARNING: Because of the paraboloid coordinates singularity, avoid	  *
*          the magnetic field calculations at the Ox axis.		  *
*									  *
* Written by I. Alexeev. 						  *
***************************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)

      COMMON/TK/B1(3),B2(3),B1A(3),B1R(3),B2A(3),B2R(3),                
     *BA(3),DB(3),AB,B(3)
      COMMON/SM/SSCP,SSP,simf,smd,ssd,ssr,smr,sbt,ss1,ss2
      COMMON/COR1/AL,BE,SQ,PQ,QA                                        
      COMMON/COR2/CFI,SFI                                               
      COMMON/COR3/R,CT,ST                                               
      COMMON/GN/V2(3)                                                   
      COMMON/T2/PI,R1,BETA0,AL0,C0,E5,AL1,BT,CPSI,SPSI,PSI,Z0,b0,bd,bd0
      COMMON/IMFd/Bimf(3),b3(3)
      COMMON/T21/BD1,R0,RKM,BK1,BKA,BKB,bkc
      COMMON/BEGF/UFCF(3)
      COMMON/BEGFC/B1CF(3),B1CFD(3),B1CFR(3),B1D(3),B1RC(3)
      COMMON/fac12/ami1,ami2,tm1,tm2
      COMMON/bfac12/bfac1(3),bfac2(3)
      DIMENSION UF(3),FF(3),V1(3),V3(3),                                
     *UZ(3,3),B1IJ(3,3),B2AA(3,3),ZU(3,3)                               
     *,EZ(3,3),EA(3,3),V4(3),V5(3)                                      
     *,V6(3),V7(3),V8(3),V9(3)                                          
     *,A2X(3,3),X2A(3,3)                                                
      X=UF(1)                                                           
      YX=UF(2)                                                          
      Z=UF(3)                                                           
      T=SQRT(YX*YX+Z*Z)                                                 
      R=SQRT(X*X+T*T)                                                   
      CT=X/R                                                            
      ST=T/R                                                            
      Z=Z+Z0                                                            
      K=1                                                               
2     K=K+1                                                             
      R=SQRT(X*X+YX*YX+Z*Z)                                             
      RX=R/R1                                                           
      X1=X/R1-0.5                                                       
      RY=RX**2-X1-0.25                                                  
      RY=SQRT(ABS(RY))                                                  
      BE=RY+X1                                                          
      BE=SQRT(ABS(BE))                                                  
      AL=RY-X1                                                          
      AL=SQRT(ABS(AL))                                                  
      PQ=AL*BE                                                          
c      T=R1*PQ                                                           
      T=SQRT(YX*YX+Z*Z)                                                 
      QA=AL*AL+BE*BE                                                    
      SQ=SQRT(QA)                                                       
      CFI=Z/T 
      SFI=YX/T                                                          
      
                                                                
      IF(K.EQ.3)GOTO3                                                   
      CALL DERY4D(B2A,B2AA)                                             
c      CALL COM(B2A,EA)                                                  
      CALL PRIS(A2X,X2A)                                                
      Z=Z-Z0                                                            
      IF(K.EQ.2) GOTO2                                                  
3     IF(AL1-AL)4,4,5                                                   
4     CALL FLYD(V1,B1IJ)                                                
      CALL PRIS(ZU,UZ)                                                  
      GOTO 6                                                            
5     CALL BEG(V1,B1IJ)                                                
12    EZ(1,1)=0.                                                        
      EZ(1,2)=-ST*V1(1)/R-CT*V1(2)/R                                    
      EZ(1,3)=0.                                                        
      EZ(2,1)=0.                                                        
      EZ(2,2)=SFI/R*(CT*V1(1)-ST*V1(2))                                 
      EZ(2,3)=(CFI*V1(1)+(CT*CFI*V1(2)-SFI*V1(3))/ST)/R                 
      EZ(3,1)=0.                                                        
      EZ(3,2)=CFI/R*(CT*V1(1)-ST*V1(2))                                 
      EZ(3,3)=-(SFI*V1(1)+(CT*SFI*V1(2)+CFI*V1(3))/ST)/R                
      ZU(1,1)=CT                                                        
      ZU(1,2)=-ST                                                       
      ZU(1,3)=0.                                                        
      ZU(2,1)=ST*SFI                                                    
      ZU(2,2)=CT*SFI                                                    
      ZU(2,3)=CFI                                                       
      ZU(3,1)=ST*CFI                                                    
      ZU(3,2)=CT*CFI                                                    
      ZU(3,3)=-SFI                                                      
      DO 7 I=1,3                                                        
      DO 7 J=1,3                                                        
      UZ(I,J)=ZU(J,I)                                                   
7     CONTINUE                                                          
6     CALL PERE(V1,B1,ZU) 
        call bdipc(uf,b0,b1d)
        if (al.ge.al1) then
        do i=1,3
        b1cf(i)=b1(i)-b1d(i)
        end do
        else      
         do i=1,3
        b1cf(i)=b1(i)
        end do
        end if
      call bring(ff)
        CALL PERE(ff,B1rc,ZU)
        do i=1,3
        b1d(i)=ssd*b1d(i)*bd0
        b1cfd(i)=smr*b1cf(i)*bd0
        b1cfr(i)=b1cf(i)-b1cfd(i)
        end do

         if (bk1.eq.0.) then
        do i=1,3
         b1rc(i)=0.
         b1cfr(i)=0.
        end do
         end if
      
        do i=1,3
        b1cf(i)=b1cfr(i)+b1cfd(i)
        b1(i)=b1cf(i)+b1d(i)
        end do

      CALL PERE(B2A,B2,A2X)                                             
      CONTINUE                                                          
      CALL PERE(B1,V2,X2A)                                              
      DO 8 I=1,3                                                        
      BA(I)=B2A(I)+V2(I)                                                
8     CONTINUE                                                          
      CALL PERE(BA,B,A2X)                                               

      b3(1)=0.0116*Bimf(1)*simf
      b3(2)=0.0978*Bimf(2)*simf
      b3(3)=0.0978*Bimf(3)*simf

      call FAC(uf,Bfac1,bfac2)    ! BFAC2=0 in this version

      do 11 j=1,3
c11    b(j)=b(j)+B3(j)+bfac1(j)+bfac2(j)
11    b(j)=b1(j)+b2(j)+B3(j)+bfac1(j)+bfac2(j)+b1rc(j)
       AB=SQRT(B(1)*B(1)+B(2)*B(2)+B(3)*B(3))                            

      DO 9 I=1,3                                                        
      FF(I)=B(I)/AB                                                     
9     CONTINUE                                                          
      RETURN                                                            
      END                                                               
C
      SUBROUTINE mas1d (flux,Br,Rs,R2)
*****************************************************************************
* Program calculating the internal coefficients for the Bessel series,	    *
* describing the magnetic field of the magnetospheric current systems in    *
* the paraboloid model. 						    *
*									    *
* INPUT PARAMETERS: 							    *
*  FLUX is the magnetic flux through the tail lobes, Wb;		    *
*  BR	is the maximum intensity of ring current, nT;			    *
*  R1	is the distance to the subsolar point of the magnetosphere, Re;     *
*  R2	is the distance to the earthward edge of geotail current sheet, Re. *
* WARNING: You should call MAS1D after each changes in the model input      *
*          parameters.							    *
* Written by I. Alexeev. 						    *
*****************************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)

      COMMON/SM/SSCP,SSP,simf,smd,ssd,ssr,smr,sbt,ss1,ss2
      COMMON/AA/BM,ZN,HN,ON
     *,CP,V7
      COMMON/T1/A1(12)
      COMMON/A/Y(3),F(5),V(3),U(3),YF(3)
      COMMON/T2/PI,R1,BETA0,AL0,C0,E5,AL1,BT,CPSI,SPSI,PSI,Z0,b0,bd,bd0
      COMMON/T21/BD1,R0,RKM,BK1,BKA,BKB,bkc
      COMMON/fac12/ami1,ami2,tm1,tm2
      dimension d1(12)
      DATA D1/0.64972264,0.21646207,
     +0.043429128,-0.000846358,-0.004917225,-0.002224403,0.94028094,
     +0.4649891,0.12928167,-0.014765534,-0.016942754,-0.022559739/
      re=6.37816
      r1=rs
      ZN = 1.
      BM = -2.*bd
      AL0=SQRT(1.+2*R2/r1)
      bt= FLUX/(1.5708*R1*R1*AL0*re*re*1.E3)          
      HN=0.06875*SQRT(HN)
      tpsi=psi
      PSI=PI/180.*PSI
      CPSI=COS(PSI)
      SPSI=SIN(PSI)
      psi=tpsi
      Z0=R1*2.*SPSI*CPSI/(3.+SPSI**2)


C******************************************************************
C      definition of the TETAM, STM, and CTMf for Region 1 FAC 
C      3,912353=31.2/7.97474; 7974.74 MWb=2*31200*10**(-9)*PI*RE**2 
C******************************************************************
          STM2=3.912353*FLUX/abs(BD)*1.e-6 
        CTM=SQRT(1.-STM2)
        STM=SQRT(STM2)
          TM1=ASIN(STM)*180./pi
          ami1=ami1/2.*(1.+CTM)
        ami2=0.75*ami1
        tm2=tm1+3.                                                   

      r0=4.
      bk1=br
      If (r2.GE.6.) then
      rkm=6.
      else 
      rkm=r2
      end if
      al1=sqrt(1.+(2.*rkm+1.)/r1)

         call masring
         b0=bd+bd1*bka
         bd0=bd/b0
      CALL MAS2D
      P=B0/R1/R1
      DO 6 I=1,6
      P=P/R1
      P1=P1/R1
      A1(I)=CPSI*(D1(I)*P)
      A1(I+6)=SPSI*(D1(I+6)*P)
    6 CONTINUE
      return
      end

      SUBROUTINE BDIPC(x,BM,B)                                           
****************************************************************************
* Program calculating the dipole magnetic field in Cartesian coordinates.  *
*									   *
* INPUT PARAMETERS:  x(3) is a coordinates of the point where the magnetic *
*                    field is being calculated, in GSM coordinates, Re;	   *
*                    BM   is dipole moment, nT*Re^3.			   *
*                    SPSI,CPSI are Sin and Cos of dipole tilt angle. 	   *
* OUTPUT PARAMETERS: B(3) is the magnetic field in GSM coordinates, nT.    *
* Written by V. Kalegaev. 						   * 
****************************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)

      dimension x(3),b(3)
      COMMON/T2/PI,R1,BETA0,AL0,C0,E5,AL1,BT,CPSI,SPSI,PSI,Z0,b0,bd,bd0
      x1=x(1)
      x2=x(2)
      x3=x(3)
      r2=(x1*x1+x2*x2+x3*x3)
      r=sqrt(r2)      
      r5=r2*r2*r
      p=x3*cpsi-x1*spsi
      br=bm/r5
      b(1)=-Br*(-r*r*spsi-3.*x1*p)
      b(2)= Br*(3.*x2*p)
      b(3)=-Br*(r*r*cpsi-3.*x3*p)
      RETURN                                                            
      END                                                               

      SUBROUTINE BDIP (P,bdp)
***********************************************************************
* Program calculating the magnetic field of the geomagnetic dipole    *
* in spherical coordinates (OX is polar axes) in the current point,   *
* defined by /cor2/, /cor3/ common blocks.			      *
*								      *
* INPUT PARAMETERS: 						      *
*     BDp  is the dipole magnetic moment, nT*Re^3;		      *
*     SPSI, CPSI - are Sin and Cos of dipole tilt angle. 	      *
* OUTPUT PARAMETERS: 						      *
*     P(3) is the dipole magnetic field, nT.			      *
* Written by I. Alexeev 					      *
***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)

      COMMON/COR2/CFI,SFI
      COMMON/COR3/R,CT,ST
      COMMON/T2/PI,R1,BETA0,AL0,C0,E5,AL1,BT,CPSI,SPSI,PSI,Z0,B0,bd,bd0
      DIMENSION P(3)
      RR=R*R
      T=Bdp/R/RR
      TC=T*CPSI
      TS=T*SPSI
      CPF=TC*CFI
       P(1)= 2*(CPF*ST-TS*CT)
       P(2)= -CPF*CT-TS*ST
       P(3)= TC*SFI
      RETURN
      END

      SUBROUTINE MASRING
***************************************************************************
* Program calculating the magnetic moment and internal coefficients 	  *
* describing the magnetic field of the ring current in  		  *
* the paraboloid model. 						  *
*									  *
* INPUT PARAMETERS: 							  *
*     BK1  is the maximum intensity of ring current, nT;		  *
*     R0   is the distance to the ring current maximum, Re;		  *
*     RKM  is the distance to the ring current edge, Re.		  *
* OUTPUT PARAMETERS: 							  *
*     BD1     is the ring current magnetic moment, nT*Re^3;		  *
*     BD1*BKA is the dipole, screening ring current, magnetic moment, 	  *
*             nT*Re^3;							  *
* WARNING: You should call MASRING after each changes in the model input  *
*          parameters.							  *
* Written by I. Alexeev 					      	  *
***************************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)

      COMMON/T21/BD1,R0,RKM,BK1,BKA,BKB,bkc
      RK =RKM
      RK3 = RK*RK*RK
      T = RK/R0/2
      T2 = T*T
      T5 = T2*T2*T
      V = 1+T2
      TS = SQRT(V)
      A = T5/V/V/TS
      BKA = A
      BKB = A/RK3/T2
      BKC = 4*R0*R0
      BD1 = BK1*RK3*T2/2/(T5-A)
      RETURN
      END

                              SUBROUTINE BRING (P)
C*********************************************************************
C Calculation of the ring current field                              *
C New version 29.08.2001					     *
C                                                                    *
C    Written by Igor I. Alexeev 				     *
C*********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)

      COMMON/T2/PI,R1,R2,BETA0,AL0,C0,AL1,BT,CPSI,SPSI,PSI,Z0,B0,BD
      COMMON/T21/BD1,R0,RKM,BK1,BKA,BKB,BKC
      COMMON/COR2/CFI,SFI
      COMMON/COR3/R,CT,ST
      DIMENSION P(3),PD(3),UFR(3)
      T=BD1*BKA/R/R/R
      P(1)=2.*T*(CPSI*CFI*ST-SPSI*CT)
      P(2)=-T*(CPSI*CT*CFI+SPSI*ST)
      P(3)=SFI*T*CPSI
      IF (R.GT.RKM) then
      RETURN
      else
      CALL BDIP (PD,BD1)
       RR=R*R
       RKT2= RR+BKC
       RKT = SQRT(RKT2)
       T2=RR/RKT2
       T3=T2*R/RKT
       TB=BKB*RR*R
       TC=T3*BKC/RKT2

        F1 = T3-TB-BKA
        F2 = F1+3*(TB-TC)

       P(1)=P(1)+PD(1)*F1
       P(2)=P(2)+PD(2)*F2
       P(3)=P(3)+PD(3)*F2
       endif
                             RETURN
                              END


      SUBROUTINE BRING1 (P,f1,f2)
***********************************************************************
* Program calculating the magnetic field of    Bring-BDR*Bdip         *
* in spherical coordinates (BDR=Bd1*BKA=B0-BD)			      *
*								      *
* INPUT PARAMETERS: 						      *
*     R0  - distance to ring current maximum;			      *
*     BK1 - maximum ring current magnetic field; 		      *
*     BD1,RKM,BKA,BKB,bkc - are calculated in MASRING subroutine.     *
* OUTPUT PARAMETERS: 						      *
*     P(3) is the  "Bring-BDR*Bdip"  magnetic field, nT.	      *
* Written by I. Alexeev 					      *    
***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)

      COMMON/T21/BD1,R0,RKM,BK1,BKA,BKB,bkc
      COMMON/COR3/R,CT,ST
      DIMENSION P(3),PD(3)

***Calculation of the Magnetic Fields of Geodipole and Ring Current 

       CALL BDIP (PD,BD1)
       RR=R*R
       RKT2= RR+BKC
       RKT = SQRT(RKT2)
       T2=RR/RKT2
       T3=T2*R/RKT
       TB=BKB*RR*R
       TC=T3*BKC/RKT2
        F1 = T3-TB-BKA
        F2 = F1+3*(TB-TC)
       P(1)= PD(1)*F1
       P(2)= PD(2)*F2
       P(3)= PD(3)*F2
      RETURN
      END

      SUBROUTINE BEG(UF,VV)
****************************************************************
* Calculation of the summary dipole, ring current and Chapman- *
* Ferraro currents at the magnetopause magnetic field in       *
* the inner magnetosphere (al<al1)			       *
* Written by I. Alexeev 				       *
****************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)

      COMMON/COR2/CFI,SFI	
      COMMON/SM/SSCP,SSP,simf,smd,ssd,ssr,smr,sbt,ss1,ss2
      COMMON/COR3/R,CT,ST
      COMMON/T2/PI,R1,BETA0,AL0,C0,E5,AL1,BT,CPSI,SPSI,PSI,Z0,b0,bd,bd0
      COMMON/T1/A1(12)
      COMMON/AA/BM,ZN,HN,ON
     *,CP,V7
      COMMON/begf/UFCF(3)
      COMMON/T21/BD1,R0,RKM,BK1,BKA,BKB,bkc
      DIMENSION UF(3),VV(3,3),ufr(3)
      DIMENSION EL(7),EL1(7),EL2(7)
      EL(1)=1.
      EL1(1)=1.
      EL(2)=CT
      EL1(2)=3.*CT
      EL2(1)=0.
      EL2(2)=3.
      EL2(3)=15.*CT
      DO 19 I=2,6
      P=1.-1./I
      EL(I+1)=(1.+P)*CT*EL(I)-P*EL(I-1)
19    EL1(I+1)=(3.-P)*CT*EL1(I)-(2.-P)*EL1(I-1)
      DO 20 I=3,6
      P=1.+3./(I-1)
      EL2(I+1)=(1.+P)*CT*EL2(I)-P*EL2(I-1)
   20 CONTINUE

***Calculation of the Magnetic field of Magnetopause Screening Currents 
      T=A1(6)
      U3=6.*T*EL1(7)
      V=EL1(6)
      U2=T*V
      U1=6.*U2
      U4=36.*T*V
      U5=6.*U3
      U6=T*EL2(7)
      U7=T*EL2(6)
      T=A1(12)
      V1=6.*T*EL(7)
      V2=T*V
      V3=6.*V1
      V4=6.*V2
      V5=6.*T*EL1(7)
      I=5

21    CONTINUE
      T=A1(I)
      U3=U3*R+I*T*EL1(I+1)
      U1=U1*R+I*T*EL1(I)
      U2=U2*R+T*EL1(I)
      U4=U4*R+I*T*I*EL1(I)
      U5=U5*R+I*I*T*EL1(I+1)
      U6=U6*R+T*EL2(I+1)
      U7=U7*R+T*EL2(I)
      T=A1(I+6)
      V1=V1*R+I*T*EL(I+1)
      V2=V2*R+T*EL1(I)
      V3=V3*R+I*T*EL(I+1)*I
      V4=V4*R+I*T*EL1(I)
      V5=V5*R+I*T*EL1(I+1)
      I=I-1
      IF(I.GT.0)GOTO 21
      UFCF(1)=-ST*CFI*U1+V1
      UFCF(2)=CFI*(CT*(U1+U2)-U3)-V2*ST
      UFCF(3)=SFI*U2
      do 24 i=1,3
 24   uf(i)=ufcf(i)
      RETURN
      END

      SUBROUTINE FLYD(P,BB)
******************************************************************
* Calculation of the summary dipole, ring current and Chapman-	 *
* Ferraro currents at the magnetopause magnetic field in	 *
* the geomagnetic tail (al>al1) 				 *
* Written by I. Alexeev 					 *     
******************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)

      COMMON/S2/ CF0(5),CF1(5),CF2(5),CF3(5),CF4(5)
      REAL*8 L,L0
      COMMON/T3/ L(6,5),L0(5)
      COMMON/COR1/AL,BE,SQ,PQ,QA
      COMMON/COR2/CFI,SFI
      COMMON/SM/SSCP,SSP,simf,smd,ssd,ssr,smr,sbt,ss1,ss2
      COMMON/T2/PI,R1,BETA0,AL0,C0,E5,AL1,BT,CPSI,SPSI,PSI,Z0,b0,bd,bd0
      DIMENSION P(3),BB(3,3)
      A2=AL*AL
      B2=BE*BE
      Y1=0.
      Y2=0.
      Y3=0.
      Y4=0.
      Y5=0.
      Y6=0.
      Y7=0.
      Y8=0.
      Y9=0.
      DO 2 N=1,5
      Z=L(1,N)
      CALL BESS(1,Z*BE,U,DU)
      CALL BESK(1,Z*AL,V,DV)
      X=CF2(N)
      Y1=Y1+X*U*DV
      Y2=Y2+X*DU*V
      Y3=Y3+CF1(N)*U*V
      X=CF3(N)
      Y4=Y4+X*U*V
      Y5=Y5+X*DU*DV
      Z=L0(N)
      X=Z*BE
      Z=Z*AL
      DU=CF0(N)
      DV=CF4(N)
      U=BESK0(Z)
      Z=BESK1(Z)
      V=BESJ0(X)
      X=BESJ1(X)
      Y6=Y6+DU*Z*V
      Y7=Y7+DU*U*X
      Y8=Y8+DV*V*U
      Y9=Y9+DV*Z*X
    2 CONTINUE
      P(1)=(-Y1*CFI+Y6)/SQ
      P(2)=(-Y2*CFI+Y7)/SQ
      P(3)=Y3*SFI/PQ
      RETURN
      END
C
      SUBROUTINE BFAC (P)
C*************************************************************
C Calculation of the magnetic field of Region 1 FAC (P(3)).
C R,CTE,STE,CFIE,SFIE are the geomagnetic spherical coordinates 
C (polar axis is directed on the Morth magnetic pole)
C STM,CTM are sin(tetam), cos(tetam),
C BFACP0=0,000098*AJ0/STM , 
C AJ0 is total field aligned current in hemisphere
C BFACP1=BFACP0*(1-CTM).
C Written by I. Alexeev
C*************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)

      COMMON/COR3/R,CT,ST
      COMMON/COR4/CTE,STE,CFIE,SFIE
      COMMON/TFAC/STM,CTM,BFAC0,BFAC1,TETAM,AJ0
      COMMON/T2/PI,R1,R2,BETA0,AL0,C0,E5,AL1,BT,CPSI,SPSI,PSI,Z0,B0,BD
      DIMENSION P(3)

        U= TETAM*PI/180.
        STM=SIN(U)
        CTM=COS(U)
        BFAC0=0.000098*AJ0/STM
        BFAC1=BFAC0*(1-CTM)


      IF (R.GT.1.AND.R.LT.R1) GOTO 3
       P(1)= 0.
       P(2)= 0.
       P(3)= 0.
                              RETURN
3      CONTINUE

      U=BFAC0/R
      U1=BFAC1/R
      IF (STE.GT.STM)        GOTO 2
      IF (CTE.GT.0)          GOTO 1
C***************************************************
C               South polar cap	 
C***************************************************
       U3=U/(1-CTE)
       P(1)= 0.
       P(2)= U3*CFIE
       P(3)= U3*SFIE

                              RETURN
C***************************************************
C            North polar cap 
C***************************************************
1      CONTINUE
       U3=U/(1+CTE)
       P(1)= 0.
       P(2)= U3*CFIE
       P(3)= -U3*SFIE

                              RETURN
C***************************************************
C              Inner magnetosphere
C***************************************************
2      CONTINUE
       U3=U1/STE/STE
       P(1)= 0.
       P(2)= U3*CFIE
       P(3)= U3*SFIE*CTE
                              RETURN
                              END
      		    

      SUBROUTINE DERY4D(P,BB)
***************************************************************
* Calculation of the magnetic field of geomagnetic tail
* current system
* Written by I. Alexeev 					 *     
***************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)

      REAL*8 L,L0
      COMMON/COR1/AL,BE,SQ,PQ,QA
      COMMON/COR2/CFI,SFI
      COMMON/SM/SSCP,SSP,simf,smd,ssd,ssr,smr,sbt,ss1,ss2
      COMMON/T2/PI,R1,BETA0,AL0,C0,E5,AL1,BT,CPSI,SPSI,PSI,Z0,b0,bd,bd0
      COMMON/T3/ L(6,5),L0(5)
      COMMON/S1/CB(6,5),CB2(6,5),
     +CD(6,5),CB3(6,5),CD2(6,5),CD3(6,5)
      DIMENSION P(3),BB(3,3),U(6,5),DU(6,5)

      if (sbt.eq.1.)then 
      IF(AL-AL0)2,2,3
2     CONTINUE
      IF(AL-174.673)40,40,41
41    PRINT 42,AL
42    FORMAT(19H GRAND EXP-DERY,AL=,E12.5)
      AL=174.670
40    CONTINUE
      E4=EXP(AL)
      DO 21 K=1,6
      M=2*K-1
      DO 21 N=1,5
      X =AL *L(K,N)
      Y=X -16.118095651
      E4=EXP(Y)
      CALL BESM(M,X,Z,DZ)
      U(K,N)=Z*E4
      DU(K,N)=DZ*E4
   21 CONTINUE
      W=0.
      GO TO 4
3     DO 22 K=1,6
      M=2*K-1
      DO 22 N=1,5
      X=AL*L(K,N)
      CALL BESK(M,X,Z,DZ)
      U(K,N)=Z*E5
      DU(K,N)=DZ*E5
22    CONTINUE
      W=+SIGN(1.D0,CFI)*C0/AL**2
      IF (CFI.EQ.0.) W=0.
4     R=+1.
      V2=0.
      V1=0.
      V3=0.
      V4=0.
      V5=0.
      V6=0.
      V7=0.
      V8=0.
      V9=0.
      BE1=BE/BETA0
      CVI=CFI
      SVI=-SFI
      SVS=2.*SFI*CFI
      CVS=2.*CFI**2-1.
      DO 23 K=1,6
      M=2*K-1
      Y1=0.
      Y2=0.
      Y3=0.
      Y4=0.
      Y5=0.
      CMFI=CVI*CVS-SVI*SVS
      SMFI=SVI*CVS+CVI*SVS
      CVI=CMFI
      SVI=SMFI
      DO 24 N=1,5
      X=L(K,N)*BE1
      CALL BESS(M,X,Z,DZ)
      IF(AL-AL0)5,5,6
    5 X1=CB(K,N)
      X2=CB2(K,N)
      X3=CB3(K,N)
      GO TO 7
6     X1=CD(K,N)
      X2=CD2(K,N)
      X3=CD3(K,N)
7     Y1=Y1+X2*Z*DU(K,N)
      Y2=Y2+X2*DZ*U(K,N)
      Y3=Y3+X1*Z*U(K,N)
      Y4=Y4+X3*Z*U(K,N)
      Y5=Y5+X3*DZ*DU(K,N)
   24 CONTINUE
      X1=CMFI/M*R
      X2=CMFI*M*R
      X3=SMFI*R
      V1=V1+Y1*X1
      V2=V2+Y2*X1
      V3=V3+Y3*X3
      V4=V4+Y4*X1
      V5=V5+Y3*X2
      V6=V6+Y2*X1
      V7=V7+Y1*X3
      V8=V8+Y2*X3
      R=-R
   23 CONTINUE
      IF(AL-AL0)60,60,70
60    V1=V1*E5
      V2=V2*E5
      V3=V3*E5
      V4=V4*E5
      V5=V5*E5
      V6=V6*E5
      V7=V7*E5
      V8=V8*E5
70    CONTINUE
      P(1)=(-V1-W*AL)/SQ
      P(2)=-V2/SQ
      P(3)=V3/PQ
      
      else
      do 1 i=1,3
      p(i)=0.
      do 1 j=1,3
1     bb(j,i)=0.
      end if 
      RETURN
      END
C
      SUBROUTINE PERE(UF,VV,DET)
      IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION UF(3),VV(3),DET(3,3)
      DO 1 I=1,3
      P=0.
      DO 2 J=1,3
      P=P+DET(I,J)*UF(J)
    2 CONTINUE
      VV(I)=P
    1 CONTINUE
      RETURN
      END
C
      SUBROUTINE PRIS(UF,VV)
      IMPLICIT REAL*8 (A-H,O-Z)

      COMMON/COR1/AL,BE,SQ,PQ,QA
      COMMON/COR2/CFI,SFI
      DIMENSION UF(3,3),VV(3,3)
      UF(1,1)=-AL/SQ
      UF(1,2)=BE/SQ
      UF(1,3)=0.
      Z=SFI/SQ
      UF(2,1)=BE*Z
      UF(2,2)=AL*Z
      UF(2,3)=CFI
      Z=CFI/SQ
      UF(3,1)=BE*Z
      UF(3,2)=AL*Z
      UF(3,3)=-SFI
      DO 1 I=1,3
      DO 1 J=1,3
      VV(I,J)=UF(J,I)
    1 CONTINUE
      RETURN
      END

      BLOCK DATA
      IMPLICIT REAL*8 (A-H,O-Z)

      REAL*8 L,L0
      COMMON/SM/SSCP,SSP,simf,smd,ssd,ssr,smr,sbt,ss1,ss2
      COMMON/T3/L(6,5),L0(5)
      COMMON/T2/PI,R1,BETA0,AL0,C0,E5,AL1,BT,CPSI,SPSI,PSI,Z0,b0,bd,bd0
      COMMON/S2/ CF0(5),CF1(5),CF2(5),CF3(5),CF4(5)
      COMMON/S5/CI0(9),CI1(9)
      COMMON/AA/BM,ZN,HN,ON
     *,CP,V7
      DATA
     *L0/3.83170597,7.01558667,10.17346814,13.3236919,16.47063005/,
     *L/1.84118390,4.2011889412,6.4156163752,8.5778364889,10.711433969,
     +12.826491226,5.3314427000,8.0152365984,10.519860874,12.932386237,
     *15.286737667,17.600266557,8.5363163000,11.345924311,13.987188630,
     *16.529365884,19.004593538,21.430854238,11.706005000,14.585848286,
     *17.312842488,19.941853366,22.501398726,25.008518704,14.863588700,
     *17.788747866,20.575514521,23.268052926,25.891277276,28.460857279
     */
      DATA
     *HN/100./,R1/10./,ON/10.4372/,B0/-31200./,PI/3.1415926/,BT/40./
      DATA SSCP,SSP,simf,smd,ssd,ssr,smr,sbt,ss1,ss2
     */1.,0.,1.,1.,1.,1.,1.,1.,1.,0./
      DATA CI0/
     *+0.003923767,-0.016476329,+0.026355372,
     *-0.020577063,+0.009162808,-0.001575649,
     *+0.002253187,+0.013285917,+0.398942280/
      DATA CI1/
     *-0.004200587,+0.017876535,-0.028953121,
     *0.022929673,-0.010315550,+0.001638014,
     *-0.003620183,-0.039880242,+0.398942280/
      END

      SUBROUTINE MAS2D
* Written by I. Alexeev 					 *     
      IMPLICIT REAL*8 (A-H,O-Z)

      REAL*8 L,L0
      COMMON/T3/L(6,5),L0(5)
      COMMON/T2/PI,R1,BETA0,AL0,C0,E5,AL1,BT,CPSI,SPSI,PSI,Z0,b0,bd,bd0
      COMMON/S1/ CB(6,5),CB2(6,5),
     *CD(6,5),CB3(6,5),CD2(6,5),
     +CD3(6,5)
      COMMON/S2/ CF0(5),CF1(5),CF2(5),CF3(5),CF4(5)
      COMMON/S5/CI0(9),CI1(9)
      COMMON/SM/SSCP,SSP,simf,smd,ssd,ssr,smr,sbt,ss1,ss2
      DIMENSION CF00(5),CF01(5),
     *cb0(6,5),cd0(6,5)
        DATA                      CB0/
     *7.493663380000,3.119532240E-1,1.706022280E-2,1.03359770E-3,
     *6.633343340E-5,4.421213830E-6,5.777707910E-3,2.103707660E-4,
     *7.538502750E-6,3.011267810E-7,1.313975350E-8,6.13862642E-10,
     *1.867244830E-5,5.965184960E-7,1.712466940E-8,5.43331737E-10,
     *1.89762765E-11,7.17571687E-13,8.545308210E-8,2.571561030E-9,
     *6.45145558E-11,1.75563377E-12,5.24527385E-14,1.70163718E-15,
     *4.43194331E-10,1.26958726E-11,2.88578820E-13,6.99784044E-15,
     *1.85385885E-16,0.
     */,

     *                             CD0/
     *5.40613258E-5,1.17552607E-3,2.01041928E-2,0.31415060400,
     *4.67288038000,67.3370338000,2.50590102E-3,0.202158093,
     *7.58793123,216.779913,5328.09218,118727.936,
     *1.72159326E-1,21.3628656,1157.19802,45465.0611,
     *1482947.1,42673478.6,14.734532, 2354.2227,
     *161804.316,7916120.89,316102914.,10974469500.,
     *1367.45791,254342.53,20562943.2,1179904930.,
     *54901114700.,2205064660000./
      DATA CF00/-94630.534652,-10798057.06,-652953011.
     +92,-30158380850.,-1198163932400.
     */
      DATA
     *CF01/-1318.693367,-190097.98824,-96033
     +82.4047,-369794441.86,-12502247266.
     */
      P=((10./R1)**3)
      P=P*B0/(-31200.)
      DO 1 N=1,5
      CF0(N)=CF00(N)*SPSI*P
      CF1(N)=CF01(N)*CPSI*P
      CF2(N)=CF1(N)*L(1,N)
      CF3(N)=CF2(N)*L(1,N)
      CF4(N)=CF0(N)*L0(N)
    1 CONTINUE
      AL01=SQRT(2.4)
    		      
cc	print *, 'mas',
cc     *cb(1,1),cb0(1,1),ska,sk,al0,al01,bt
      	      
      DO 2 K=1,6
      M=2*K-1
      DO 2 N=1,5
      ZZ=L(k,N)*AL01
      ZA=L(k,N)*AL0
      CALL BESK(M,ZZ,SK ,DSK )
      CALL BESM(M,ZZ,SI ,DSI )
      CALL BESK(M,ZA,SKA,DSKA)
      CALL BESM(M,ZA,SIA,DSIA)
      SIII=ZA-ZZ
      SIA=SIA/SI*EXP(SIII)
      CB(K,N)=BT/40.*CB0(K,N)*AL0*SKA/AL01/SK
      CD(K,N)=BT/40.*CD0(K,N)*AL0*SIA/AL01
      CB2(K,N)=L(K,N)*CB(K,N)
      CB3(K,N)=L(K,N)*CB2(K,N)
      CD2(K,N)=L(K,N)*CD(K,N)
      CD3(K,N)=L(K,N)*CD2(K,N)
    2 CONTINUE
      AL1=SQRT(2.4)
      C0=   61.96773354*BT/40.*al0/al01
      E5=10.**7
      BETA0=1.0+0.2*(PSI/30.)**2
*--------------------------------------------------------------
*********> place to insert block "test.for"<*******************
*--------------------------------------------------------------
      RETURN
      END
C
      SUBROUTINE BESK(M,X,V,DV)
* Written by I. Alexeev 					 *     
      IMPLICIT REAL*8 (A-H,O-Z)

      U=BESK0(X)
      V=BESK1(X)
      L=M-1
      IF(L)3,4,3
3     DO 2 K=1,L
      S=2.*K*V/X+U
      U=V
2     V=S
4     DV=-(M*V/X+U)
      RETURN
      END
C
      FUNCTION RINT(DZ,FX,FZ)
      DX=FX*DZ/FZ
      RINT=DX
      return
      END
C
      SUBROUTINE BESM(M,X,V,DV)
* Written by I. Alexeev 					 *     
      IMPLICIT REAL*8 (A-H,O-Z)

      COMMON/S5/CI0(9),CI1(9)
      IF(X.GT.3.75)GO TO 5
      XA=-X
      IF(XA-174.673)50,50,51
 51   PRINT 52,X
52    FORMAT(18H GRAND EXP-BESM,X=,E12.5)
      X=-174.670
50    CONTINUE
      E=EXP(-X)
      V=BSI(M,X)*E
      DV=0.5*E*(BSI(M-1,X)+BSI(M+1,X))
      RETURN
5     CONTINUE
      IF(X)90,91,91
90    PRINT 92,X
92    FORMAT(19H NEGATIVE X-BESM,X=,E12.5)
   91 X=ABS(X)
      U=UG(CI0(1),X)/SQRT(X)
      V=UG(CI1(1),X)/SQRT(X)
      L=M-1
      IF(L)3,4,3
3     DO 2 K=1,L
      S=U-2.*K*V/X
      U=V
2     V=S
4     DV=U-M*V/X
      RETURN
      END
C
      SUBROUTINE BESS(M,X,V,DV)
* Written by I. Alexeev 					 *     
      IMPLICIT REAL*8 (A-H,O-Z)

      IF(X.GT.3.75)GO TO 5
      V=BSJ(M,X)
      DV=0.5*(BSJ(M-1,X)-BSJ(M+1,X))
      RETURN
5     U=BESJ0(X)
      V=BESJ1(X)
      L=M-1
      IF(L)3,4,3
3     DO 2 K=1,L
      S=2*K*V/X-U
      U=V
2      V=S
4     DV=U-M*V/X
      RETURN
      END
C
      FUNCTION BSJ(N,X)
* Written by I. Alexeev 					 *     
      IMPLICIT REAL*8 (A-H,O-Z)

      SUM=1.
      P=1.
      Z=-X*X/4.
      DO 2 K=1,7
      P=P*Z/K/(K+N)
2     SUM=SUM+P
      IF(N.LE.0)GO TO 3
      DO 1 K=1,N
1     SUM=SUM/K
  3   CONTINUE
7     FORMAT(16H EXP NEGATIVE,N=,I3,2HX=,E12.5)
      IF(X)4,5,4
5     CONTINUE
      BSJ=0.
      PRINT 7,N,X
      GO TO 6
4     CONTINUE
      BSJ=SUM*(X/2.)**N
6     RETURN
      END
C
      FUNCTION BSI(N,X)
* Written by I. Alexeev 					 *     
      IMPLICIT REAL*8 (A-H,O-Z)

      SUM=1.
      P=1.
      Z=X*X/4.
      DO 2 K=1,7
      P=P*Z/K/(K+N)
2     SUM=SUM+P
      IF(N.LE.0)GO TO 3
      DO 1 K=1,N
1     SUM=SUM/K
  3   CONTINUE
7     FORMAT(16H EXP NEGATIVE,N=,I3,2HX=,E12.5)
      IF(X)4,5,4
5     CONTINUE
      PRINT 7,N,X
      BSI=0.
      GO TO 6
4     CONTINUE
      BSI=SUM*(X/2.)**N
6     RETURN
      END
C
      FUNCTION UG(V,X)
* Written by I. Alexeev 					 *     
      IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION V(9)
      UG=V(1)
      DO 7 I=2,9
      UG=UG*(3.75/X)+V(I)
7     CONTINUE
      RETURN
      END
C
C
      FUNCTION BESJY(X)
C*************************************************************************
C Calculation of the Bessel functions J0(x), J1(x), Y0(x) or Y1(x)       *
C   INPUT:                                                               * 
C        x - argument of the Bessel function                             * 
C                                                                        * 
C   OUTPUT: 					         		 *
C	BESJ0(X)= J0(x)                                           	 *
C	BESJ1(X)= J1(x)                                                  *  
C	BESY0(X)= Y0(x)                                           	 *
C	BESY1(X)= Y1(x)                                                  *  
C   Special subroutine from MSU Computer Center library  		 *
C*************************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)

      LOGICAL L
C
      ENTRY BESJ0(X)
C
      L=.TRUE.
      V=ABS(X)
      IF(V .GE. 8.0) GO TO 4
    8 F=0.0625*X**2-2.0
      A =           - 0.00000 00000 000008
      B = F * A     + 0.00000 00000 000413
      A = F * B - A - 0.00000 00000 019438
      B = F * A - B + 0.00000 00000 784870
      A = F * B - A - 0.00000 00026 792535
      B = F * A - B + 0.00000 00760 816359
      A = F * B - A - 0.00000 17619 469078
      B = F * A - B + 0.00003 24603 288210
      A = F * B - A - 0.00046 06261 662063
      B = F * A - B + 0.00481 91800 694676
      A = F * B - A - 0.03489 37694 114089
      B = F * A - B + 0.15806 71023 320973
      A = F * B - A - 0.37009 49938 726498
      B = F * A - B + 0.26517 86132 033368
      A = F * B - A - 0.00872 34423 528522
      A = F * A - B + 0.31545 59429 497802
      BESJY=0.5*(A-B)
      IF(L) RETURN
C
      A =           + 0.00000 00000 000016
      B = F * A     - 0.00000 00000 000875
      A = F * B - A + 0.00000 00000 040263
      B = F * A - B - 0.00000 00001 583755
      A = F * B - A + 0.00000 00052 487948
      B = F * A - B - 0.00000 01440 723327
      A = F * B - A + 0.00000 32065 325377
      B = F * A - B - 0.00005 63207 914106
      A = F * B - A + 0.00075 31135 932578
      B = F * A - B - 0.00728 79624 795521
      A = F * B - A + 0.04719 66895 957634
      B = F * A - B - 0.17730 20127 811436
      A = F * B - A + 0.26156 73462 550466
      B = F * A - B + 0.17903 43140 771827
      A = F * B - A - 0.27447 43055 297453
      A = F * A - B - 0.06629 22264 065699
      BESJY=0.636619772367581*LOG(X)*BESJY+0.5*(A-B)
      RETURN
C
    4 F=256.0/X**2-2.0
      B =           + 0.00000 00000 000007
      A = F * B     - 0.00000 00000 000051
      B = F * A - B + 0.00000 00000 000433
      A = F * B - A - 0.00000 00000 004305
      B = F * A - B + 0.00000 00000 051683
      A = F * B - A - 0.00000 00000 786409
      B = F * A - B + 0.00000 00016 306465
      A = F * B - A - 0.00000 00517 059454
      B = F * A - B + 0.00000 30751 847875
      A = F * B - A - 0.00053 65220 468132
      A = F * A - B + 1.99892 06986 950373
      P=A-B
      B =           - 0.00000 00000 000006
      A = F * B     + 0.00000 00000 000043
      B = F * A - B - 0.00000 00000 000334
      A = F * B - A + 0.00000 00000 003006
      B = F * A - B - 0.00000 00000 032067
      A = F * B - A + 0.00000 00000 422012
      B = F * A - B - 0.00000 00007 271916
      A = F * B - A + 0.00000 00179 724572
      B = F * A - B - 0.00000 07414 498411
      A = F * B - A + 0.00006 83851 994261
      A = F * A - B - 0.03111 17092 106740
      Q=8.0*(A-B)/V
      F=V-0.785398163397448
      A=COS(F)
      B=SIN(F)
      F=0.398942280401432/SQRT(V)
      IF(L) GO TO 6
      BESJY=F*(Q*A+P*B)
      RETURN
    6 BESJY=F*(P*A-Q*B)
      RETURN
C
      ENTRY BESJ1(X)
C
      L=.TRUE.
      V=ABS(X)
      IF(V .GE. 8.0) GO TO 5
    3 F=0.0625*X**2-2.0
      B =           + 0.00000 00000 000114
      A = F * B     - 0.00000 00000 005777
      B = F * A - B + 0.00000 00000 252812
      A = F * B - A - 0.00000 00009 424213
      B = F * A - B + 0.00000 00294 970701
      A = F * B - A - 0.00000 07617 587805
      B = F * A - B + 0.00001 58870 192399
      A = F * B - A - 0.00026 04443 893486
      B = F * A - B + 0.00324 02701 826839
      A = F * B - A - 0.02917 55248 061542
      B = F * A - B + 0.17770 91172 397283
      A = F * B - A - 0.66144 39341 345433
      B = F * A - B + 1.28799 40988 576776
      A = F * B - A - 1.19180 11605 412169
      A = F * A - B + 1.29671 75412 105298
      BESJY=0.0625*(A-B)*X
      IF(L) RETURN
C
      B =           - 0.00000 00000 000244
      A = F * B     + 0.00000 00000 012114
      B = F * A - B - 0.00000 00000 517212
      A = F * B - A + 0.00000 00018 754703
      B = F * A - B - 0.00000 00568 844004
      A = F * B - A + 0.00000 14166 243645
      B = F * A - B - 0.00002 83046 401495
      A = F * B - A + 0.00044 04786 298671
      B = F * A - B - 0.00513 16411 610611
      A = F * B - A + 0.04231 91803 533369
      B = F * A - B - 0.22662 49915 567549
      A = F * B - A + 0.67561 57807 721877
      B = F * A - B - 0.76729 63628 866459
      A = F * B - A - 0.12869 73843 813500
      A = F * A - B + 0.04060 82117 718685
      BESJY=0.636619772367581*LOG(X)*BESJY-0.636619772367581/X
     1     +0.0625*(A-B)*X
      RETURN
C
    5 F=256.0/X**2-2.0
      B =           - 0.00000 00000 000007
      A = F * B     + 0.00000 00000 000055
      B = F * A - B - 0.00000 00000 000468
      A = F * B - A + 0.00000 00000 004699
      B = F * A - B - 0.00000 00000 057049
      A = F * B - A + 0.00000 00000 881690
      B = F * A - B - 0.00000 00018 718907
      A = F * B - A + 0.00000 00617 763396
      B = F * A - B - 0.00000 39872 843005
      A = F * B - A + 0.00089 89898 330859
      A = F * A - B + 2.00180 60817 200274
      P=A-B
      B =           + 0.00000 00000 000007
      A = F * B     - 0.00000 00000 000046
      B = F * A - B + 0.00000 00000 000360
      A = F * B - A - 0.00000 00000 003264
      B = F * A - B + 0.00000 00000 035152
      A = F * B - A - 0.00000 00000 468636
      B = F * A - B + 0.00000 00008 229193
      A = F * B - A - 0.00000 00209 597814
      B = F * A - B + 0.00000 09138 615258
      A = F * B - A - 0.00009 62772 354916
      A = F * A - B + 0.09355 55741 390707
      Q=8.0*(A-B)/V
      F=V-2.356194490192345
      A=COS(F)
      B=SIN(F)
      F=0.398942280401432/SQRT(V)
      IF(L) GO TO 7
      BESJY=F*(Q*A+P*B)
      RETURN
    7 BESJY=F*(P*A-Q*B)
      IF(X .LT. 0.0) BESJY=-BESJY
      RETURN
C
      ENTRY BESY0(X)
C
      IF(X .LE. 0.0) GO TO 9
      L=.FALSE.
      V=X
      IF(V .GE. 8.0) GO TO 4
      GO TO 8
C
      ENTRY BESY1(X)
C
      IF(X .LE. 0.0) GO TO 9
      L=.FALSE.
      V=X
      IF(V .GE. 8.0) GO TO 5
      GO TO 3
C
    9 BESJY=0.
      PRINT 100,X
      RETURN
100   FORMAT(1X,32HBESJY...NON-POSITIVE ARGUMENT X=,E12.5)
C
      END


      FUNCTION BESIK(X)
C*************************************************************************
C Calculation of the Bessel functions K0(x), K1(x), I0(x) or I1(x)       *
C   INPUT:                                                               * 
C        x - argument of the Bessel function                             * 
C                                                                        * 
C   OUTPUT: 					         		 *
C	BESI0(X)= I0(x)                                           	 *
C	BESI1(X)= I1(x)                                                  *  
C	BESK0(X)= K0(x)                                           	 *
C	BESK1(X)= K1(x)                                                  *  
C	EBESI0(X)= exp(-|x|)*I0(x)                               	 *
C	EBESI1(X)= exp(-|x|)*I1(x)                                       *  
C	EBESK0(X)= exp(x)*K0(x)                                 	 *
C	EBESK1(X)= exp(x)*K1(x)                                          *  
C   Special subroutine from MSU Computer Center library         	 *
C*************************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)

      LOGICAL L,E
C
      ENTRY EBESI0(X)
C
      E=.TRUE.
      GO TO 1
      ENTRY BESI0(X)
      E=.FALSE.
    1 L=.TRUE.
      V=ABS(X)
      IF(V .GE. 8.0) GO TO 4
    8 F=0.0625*X**2-2.0
      A =               0.00000 00000 00002
      B = F * A     +   0.00000 00000 00120
      A = F * B - A +   0.00000 00000 06097
      B = F * A - B +   0.00000 00002 68828
      A = F * B - A +   0.00000 00101 69727
      B = F * A - B +   0.00000 03260 91051
      A = F * B - A +   0.00000 87383 15497
      B = F * A - B +   0.00019 24693 59688
      A = F * B - A +   0.00341 63317 66012
      B = F * A - B +   0.04771 87487 98174
      A = F * B - A +   0.50949 33654 39983
      B = F * A - B +   4.01167 37601 79349
      A = F * B - A +  22.27481 92424 62231
      B = F * A - B +  82.48903 27440 24100
      A = F * B - A + 190.49432 01727 42844
      A = F * A - B + 255.46687 96243 62167
      BESIK=0.5*(A-B)
      IF(L .AND. E) BESIK=EXP(-V)*BESIK
      IF(L) RETURN
      A =           +   0.00000 00000 00003
      B = F * A     +   0.00000 00000 00159
      A = F * B - A +   0.00000 00000 07658
      B = F * A - B +   0.00000 00003 18588
      A = F * B - A +   0.00000 00112 81211
      B = F * A - B +   0.00000 03351 95256
      A = F * B - A +   0.00000 82160 25940
      B = F * A - B +   0.00016 27083 79043
      A = F * B - A +   0.00253 63081 88086
      B = F * A - B +   0.03008 07224 20512
      A = F * B - A +   0.25908 44324 34900
      B = F * A - B +   1.51153 56760 29228
      A = F * B - A +   5.28363 28668 73920
      B = F * A - B +   8.00536 88687 00334
      A = F * B - A -   4.56343 35864 48395
      A = F * A - B -  21.05766 01774 02440
      BESIK=-LOG(0.125*X)*BESIK+0.5*(A-B)
      IF(E) BESIK=EXP(X)*BESIK
      RETURN
    4 F=32.0/V-2.0
      B =           - 0.00000 00000 00001
      A = F * B     - 0.00000 00000 00001
      B = F * A - B + 0.00000 00000 00004
      A = F * B - A + 0.00000 00000 00010
      B = F * A - B - 0.00000 00000 00024
      A = F * B - A - 0.00000 00000 00104
      B = F * A - B + 0.00000 00000 00039
      A = F * B - A + 0.00000 00000 00966
      B = F * A - B + 0.00000 00000 01800
      A = F * B - A - 0.00000 00000 04497
      B = F * A - B - 0.00000 00000 33127
      A = F * B - A - 0.00000 00000 78957
      B = F * A - B + 0.00000 00000 29802
      A = F * B - A + 0.00000 00012 38425
      B = F * A - B + 0.00000 00085 13091
      A = F * B - A + 0.00000 00568 16966
      B = F * A - B + 0.00000 05135 87727
      A = F * B - A + 0.00000 72475 91100
      B = F * A - B + 0.00017 27006 30778
      A = F * B - A + 0.00844 51226 24921
      A = F * A - B + 2.01655 84109 17480
      BESIK=0.199471140200717*(A-B)/SQRT(V)
      IF(E) RETURN
      BESIK=EXP(V)*BESIK
      RETURN
      ENTRY EBESI1(X)
      E=.TRUE.
      GO TO 2
      ENTRY BESI1(X)
      E=.FALSE.
    2 L=.TRUE.
      V=ABS(X)
      IF(V .GE. 8.0) GO TO 3
    7 F=0.0625*X**2-2.0
      A =           +   0.00000 00000 00001
      B = F * A     +   0.00000 00000 00031
      A = F * B - A +   0.00000 00000 01679
      B = F * A - B +   0.00000 00000 79291
      A = F * B - A +   0.00000 00032 27617
      B = F * A - B +   0.00000 01119 46285
      A = F * B - A +   0.00000 32641 38122
      B = F * A - B +   0.00007 87567 85754
      A = F * B - A +   0.00154 30190 15627
      B = F * A - B +   0.02399 30791 47841
      A = F * B - A +   0.28785 55118 04672
      B = F * A - B +   2.57145 99063 47755
      A = F * B - A +  16.33455 05525 22066
      B = F * A - B +  69.39591 76337 34448
      A = F * B - A + 181.31261 60405 70265
      A = F * A - B + 259.89023 78064 77292
      BESIK=0.0625*(A-B)*X
      IF(L .AND. E) BESIK=EXP(-V)*BESIK
      IF(L) RETURN
      A =           +   0.00000 00000 00001
      B = F * A     +   0.00000 00000 00042
      A = F * B - A +   0.00000 00000 02163
      B = F * A - B +   0.00000 00000 96660
      A = F * B - A +   0.00000 00036 96783
      B = F * A - B +   0.00000 01193 67971
      A = F * B - A +   0.00000 32025 10692
      B = F * A - B +   0.00007 00106 27855
      A = F * B - A +   0.00121 70569 94516
      B = F * A - B +   0.01630 00492 89816
      A = F * B - A +   0.16107 43016 56148
      B = F * A - B +   1.10146 19930 04852
      A = F * B - A +   4.66638 70268 62842
      B = F * A - B +   9.36161 78313 95389
      A = F * B - A -   1.83923 92242 86199
      A = F * A - B -  26.68809 54808 62668
      BESIK=LOG(0.125*X)*BESIK+1.0/X-0.0625*(A-B)*X
      IF(E) BESIK=EXP(X)*BESIK
      RETURN
    3 F=32.0/V-2.0
      B =           + 0.00000 00000 00001
      A = F * B     + 0.00000 00000 00001
      B = F * A - B - 0.00000 00000 00005
      A = F * B - A - 0.00000 00000 00010
      B = F * A - B + 0.00000 00000 00026
      A = F * B - A + 0.00000 00000 00107
      B = F * A - B - 0.00000 00000 00053
      A = F * B - A - 0.00000 00000 01024
      B = F * A - B - 0.00000 00000 01804
      A = F * B - A + 0.00000 00000 05103
      B = F * A - B + 0.00000 00000 35408
      A = F * B - A + 0.00000 00000 81531
      B = F * A - B - 0.00000 00000 47563
      A = F * B - A - 0.00000 00014 01141
      B = F * A - B - 0.00000 00096 13873
      A = F * B - A - 0.00000 00659 61142
      B = F * A - B - 0.00000 06297 24239
      A = F * B - A - 0.00000 97321 46728
      B = F * A - B - 0.00027 72053 60764
      A = F * B - A - 0.02446 74429 63276
      A = F * A - B + 1.95160 12046 52572
      BESIK=0.199471140200717*(A-B)/SQRT(V)
      IF(X .LT. 0.0) BESIK=-BESIK
      IF(E) RETURN
      BESIK=EXP(V)*BESIK
      RETURN
      ENTRY EBESK0(X)
      E=.TRUE.
      GO TO 11
      ENTRY BESK0(X)
      E=.FALSE.
   11 IF(X .LE. 0.0) GO TO 9
      L=.FALSE.
      V=X
      IF(X .LT. 5.0) GO TO 8
      F=20.0/X-2.0
      A =           - 0.00000 00000 00002
      B = F * A     + 0.00000 00000 00011
      A = F * B - A - 0.00000 00000 00079
      B = F * A - B + 0.00000 00000 00581
      A = F * B - A - 0.00000 00000 04580
      B = F * A - B + 0.00000 00000 39044
      A = F * B - A - 0.00000 00003 64547
      B = F * A - B + 0.00000 00037 92996
      A = F * B - A - 0.00000 00450 47338
      B = F * A - B + 0.00000 06325 75109
      A = F * B - A - 0.00001 11066 85197
      B = F * A - B + 0.00026 95326 12763
      A = F * B - A - 0.01131 05046 46928
      A = F * A - B + 1.97681 63484 61652
      BESIK=0.626657068657750*(A-B)/SQRT(X)
      IF(E) RETURN
      BESIK=EXP(-X)*BESIK
      RETURN
      ENTRY EBESK1(X)
      E=.TRUE.
      GO TO 12
      ENTRY BESK1(X)
      E=.FALSE.
   12 IF(X .LE. 0.0) GO TO 9
      L=.FALSE.
      V=X
      IF(X .LT. 5.0) GO TO 7
      F=20.0/X-2.0
      A =           + 0.00000 00000 00002
      B = F * A     - 0.00000 00000 00013
      A = F * B - A + 0.00000 00000 00089
      B = F * A - B - 0.00000 00000 00663
      A = F * B - A + 0.00000 00000 05288
      B = F * A - B - 0.00000 00000 45757
      A = F * B - A + 0.00000 00004 35417
      B = F * A - B - 0.00000 00046 45555
      A = F * B - A + 0.00000 00571 32218
      B = F * A - B - 0.00000 08451 72048
      A = F * B - A + 0.00001 61850 63810
      B = F * A - B - 0.00046 84750 28167
      A = F * B - A + 0.03546 52912 43331
      A = F * A - B + 2.07190 17175 44716
      BESIK=0.626657068657750*(A-B)/SQRT(X)
      IF(E) RETURN
      BESIK=EXP(-X)*BESIK
      RETURN
    9 BESIK=0.
      PRINT 200,X
200   FORMAT(1X,32HBESIK...NON-POSITIVE ARGUMENT X=,E12.5)
      RETURN
      END
C
       subroutine FAC(x,B1,b2)
*******************************************************************
c
C version of 27.08.98
c FAC.for calculation of the magnetic field from field-aligned
c current-I and II using dipole configuration
c with equatorial current.
c     
c         ami1 - total R-I current in MA, 
c         tm1 - colatitude of  R-I current in degrees
c         ami2 - total R-II current in MA, 
c         tm2 - colatitude of  R-II current in degrees
c x(1-3) - GSM coord. of point, B1(1-3) - mag. field GSM coord. (nT),
c                               B2(1-3) - mag. field GSM coord. (nT)
c      B2=0 in this version of the paraboloid model
* Written by V. Kalegaev         					 *     
********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)

      COMMON/COR3/R,CT,ST
      COMMON/COR4/CTE,STE,CFIE,SFIE
      COMMON/TFAC/STM,CTM,BFAC0,BFAC1,TETAM,AJ0
      COMMON/T2/PI,R1,R2,BETA0,AL0,C0,E5,AL1,BT,CPSI,SPSI,PSI,Z0,B0,BD
      COMMON/SM/SSCP,SSP,simf,smd,ssd,ssr,smr,sbt,ss1,ss2
      COMMON/fac12/ami1,ami2,tm1,tm2
      DIMENSION x(3),B2(3),xsm(3),Bsm(3),xr(3),Br(3),sm2gsm(3,3)
     *,zu(3,3),Bsm1(3),b1(3),Br1(3)
      do i=1,3
      br(i)=0.
      b1(i)=0.
      b2(i)=0.
      end do
      tetam=tm1
      tm0=tm2
      aj0=ami1*1.e6
      ami=ami2
      p=pi/180. 
      call SMtoGSM(SM2GSM)
      call PERE2(x,xsm,sm2gsm,-1)              
      r=sqrt(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
      cte=xsm(3)/r
      ste=sqrt(1-cte*cte)
       if (ste.eq.0.) then 
        cfie=1.
        sfie=0.
       else
        cfie=xsm(1)/r/ste
        sfie=xsm(2)/r/ste
       end if
      xr(1)=r
      xr(3)=acos(cfie)/p
      if (xsm(2).lt.0.) xr(3)=360-xr(3)
      xr(2)=acos(cte)/p
      ZU(3,1)=Cte                                                        
      ZU(3,2)=-STe                                                      
      ZU(3,3)=0.                                                        
      ZU(2,1)=STe*SFIe                                                    
      ZU(2,2)=CTe*SFIe                                                    
      ZU(2,3)=CFIe                                                       
      ZU(1,1)=STe*CFIe                                                    
      ZU(1,2)=CTe*CFIe                                                    
      ZU(1,3)=-SFIe  
      ss2=0.             ! FAC2 Currents are cancelled!                                                    
      if (ss2.eq.0.) goto1
c        call FAC2d(ami,tm0,xr,br)     
      call PERE2(br,bsm,zu,1)              
      call PERE2(bsm,b2,sm2gsm,1)
1     if (ss1.eq.0.) return
        call bFAC(br1)     
      call PERE2(br1,bsm1,zu,1)              
      call PERE2(bsm1,b1,sm2gsm,1)

      return          
      end


*********************************************************
*     AUXILLIARY SUBROUTINES
*********************************************************
C
      SUBROUTINE PERE2(A,B,T,K)
**********************************************************
*     Transition A into B vectors by T (K>0)
*                       B=T*A      
*     or T^{-1} matrices (K<=0)
*                       B=T^{-1}*A      
* Written by V. Kalegaev         					 *     
**********************************************************
      IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION A(3),B(3),T(3,3)
      if (k) 1,1,2
2     b(1)=t(1,1)*a(1)+t(1,2)*a(2)+t(1,3)*a(3)
      b(2)=t(2,1)*a(1)+t(2,2)*a(2)+t(2,3)*a(3)
      b(3)=t(3,1)*a(1)+t(3,2)*a(2)+t(3,3)*a(3)
      return
1     b(1)=t(1,1)*a(1)+t(2,1)*a(2)+t(3,1)*a(3)
      b(2)=t(1,2)*a(1)+t(2,2)*a(2)+t(3,2)*a(3)
      b(3)=t(1,3)*a(1)+t(2,3)*a(2)+t(3,3)*a(3)
      return
      end
C
      SUBROUTINE SMtoGSM(SM2GSM)
****************************************************************
*     Calculation of the transition matrix from SM coordinates to 
*     GSM ones:
*                  VectGSM=(SM2GSM)*VectSM
* Written by V. Kalegaev         					 *     
****************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)

      COMMON/T2/PI,R1,BETA0,AL0,C0,E5,AL1,BT,CPSI,SPSI,PSI,Z0,b0,bd,bd0
      DIMENSION SM2GSM(3,3)
      SM2GSM(1,1)=CPSI
      SM2GSM(1,2)=0.
      SM2GSM(1,3)=-SPSI
      SM2GSM(2,1)=0.
      SM2GSM(2,2)=1.
      SM2GSM(2,3)=0.
      SM2GSM(3,1)=SPSI
      SM2GSM(3,2)=0.
      SM2GSM(3,3)=CPSI
      return
      END
C
      SUBROUTINE PSTATUS(X1,X2,X3,X4,X5,X6,X7)
*****************************************************************
*     Determination of the parameters providing the model tuning
* Written by V. Kalegaev         					 *     
*****************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)

      COMMON/SM/SSCP,SSP,simf,smd,ssd,ssr,smr,sbt,ss1,ss2
      SSD =x1 ! dipole field on/off (1/0)
      SSR =x2 ! RC field on/off (1/0)
      SBT =x3 ! tail current field on/off (1/0)
      SMD =x4 ! dipole shielding field on/off (1/0)
      SMR =x5 ! RC shielding field IMF on/off (1/0)
      SS1 =x6 ! Region 1 FAC field on/off (1/0)
      SIMF=x7 ! IMF on/off (1/0)
      SS2 =0. ! Region 2 FAC field off (to be zero yet!)
      return
      END
C
      SUBROUTINE ALX00_SUN(IYR,IDAY,IHOUR,MIN,ISEC,GST,
     +                     SLONG,SRASN,SDEC)
C
C  CALCULATES FOUR QUANTITIES NECESSARY FOR COORDINATE TRANSFORMATIONS
C  WHICH DEPEND ON SUN POSITION (AND, HENCE, ON UNIVERSAL TIME AND SEASON)
C
C-------  INPUT PARAMETERS:
C  IYR,IDAY,IHOUR,MIN,ISEC -  YEAR, DAY, AND UNIVERSAL TIME IN HOURS, MINUTES,
C    AND SECONDS  (IDAY=1 CORRESPONDS TO JANUARY 1).
C
C-------  OUTPUT PARAMETERS:
C  GST - GREENWICH MEAN SIDEREAL TIME, SLONG - LONGITUDE ALONG ECLIPTIC
C  SRASN - RIGHT ASCENSION,  SDEC - DECLINATION  OF THE SUN (RADIANS)
C  THIS SUBROUTINE HAS BEEN COMPILED FROM: RUSSELL C.T., COSM.ELECTRO-
C  DYN., 1971, V.2,PP.184-196.
C
C
C                   AUTHOR: Gilbert D. Mead
C
C
      IMPLICIT REAL*8 (A-H,O-Z)

      REAL*8 GST,SLONG,SRASN,SDEC,RAD,T,VL,G,OBLIQ,SOB,SLP,SIND,
     1       COSD,SC

      INTEGER IYR,IDAY,IHOUR,MIN,ISEC
      COMMON/ddd/ sind,cosd
      DOUBLE PRECISION DJ,FDAY
      DATA RAD/57.295779513/
      IF(IYR.LT.1901.OR.IYR.GT.2099) RETURN
      FDAY=DFLOAT(IHOUR*3600+MIN*60+ISEC)/86400.D0
      DJ=365*(IYR-1900)+(IYR-1901)/4+IDAY-0.5D0+FDAY
      T=DJ/36525.
      VL=DMOD(279.696678+0.9856473354*DJ,360.D0)
      GST=DMOD(279.690983+.9856473354*DJ+360.*FDAY+180.,360.D0)/RAD
      G=DMOD(358.475845+0.985600267*DJ,360.D0)/RAD
      SLONG=(VL+(1.91946-0.004789*T)*SIN(G)+0.020094*SIN(2.*G))/RAD
      IF(SLONG.GT.6.2831853) SLONG=SLONG-6.2831853
      IF (SLONG.LT.0.) SLONG=SLONG+6.2831853
      OBLIQ=(23.45229-0.0130125*T)/RAD
      SOB=SIN(OBLIQ)
      SLP=SLONG-9.924E-5
C
C   THE LAST CONSTANT IS A CORRECTION FOR THE ANGULAR ABERRATION  DUE TO
C   THE ORBITAL MOTION OF THE EARTH
C
      SIND=SOB*SIN(SLP)
      COSD=SQRT(1.-SIND**2)
      SC=SIND/COSD
      SDEC=ATAN(SC)
      SRASN=3.141592654-ATAN2(COS(OBLIQ)/SOB*SC,-COS(SLP)/COSD)
      RETURN
      END
