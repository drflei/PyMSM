# 1 "um530.f"
      SUBROUTINE UM530
     :          (mpos, mb, ifail)
C
C!    Evaluate the magnetic field
C_    er25
C                       
      INCLUDE 'structure.h'
cDEC$ IF DEFINED (_x86_)
cDEC$ ATTRIBUTES DLLEXPORT :: UM530
cDEC$ ENDIF
C
      EXTERNAL UM531, UM532, UM533, UM535
C      
C     INTERFACE
C
        TYPE(zgeo) :: mpos
        TYPE(zvec) :: mb
        INTEGER*4     ifail
C
      COMMON /UC140/  mint, mext, msun 
C                                        
        TYPE(zimf) :: mint
        TYPE(zemf) :: mext
        TYPE(zsun) :: msun
C
      COMMON /UC160
     :               /pi, deg, re, gmagmo, eclipt, geoid, uma
        REAL*8        pi, deg, re, gmagmo, eclipt, geoid(3), uma(30)
C
      COMMON /UC192
     :               /xrmin, xbmin, xtmin, xbmax, epslon, epsfl,
     :                pvet, fvet, epsomeg, dltlat
        REAL*8        xrmin, xbmin, xtmin, xbmax, epslon, epsfl,
     :                pvet, fvet, epsomeg, dltlat
C 
C     VARIABLES                     
C     
        REAL*8        bir, bit, bip, ber, bet, bep
        TYPE(zgeo) :: mdpos
        REAL*8        rre1, sit1, cot1, sip, cop
        REAL*8        rre2, sit2, cot2, sip2, cop2
        REAL*8        slt, csgmla, rm
        integer*4     k531, k532
C
C     CODE
C
      mb%dnrm = 0.0d0
      ifail   = 0
      k531    = 0
      k532    = 0
C      
C  a/ Check the position
C
      if( mpos%colat .lt. 0.0d0   .or.
     :    mpos%colat .gt. 180.0d0 )then
        ifail = -53005
        return
      endif
      if( mpos%radius .lt. xrmin  )then
c        print*,mpos%radius
        ifail = -53001
        return
      endif
      CALL UM533( mpos, rm, ifail )
      if( ifail .lt. 0 ) return
      ifail = 0
C
      slt    = SIN( mint%colat * deg )
      sit1   = SIN( mpos%colat * deg )
      cot1   = COS( mpos%colat * deg )
      sip    = SIN( mpos%elong * deg )
      cop    = COS( mpos%elong * deg )     
C
      csgmla = sit1 * cop * slt * COS( mint%elong * deg ) +
     :         sit1 * sip * slt * SIN( mint%elong * deg ) +
     :         cot1 * COS( mint%colat * deg)         
C
C  b/ Compute the geomagnetic field
C
      sip2   = SIN( (mpos%elong+mint%saarot) * deg )
      cop2   = COS( (mpos%elong+mint%saarot) * deg )     
c
      rre1   = mpos%radius / re
      if( mint%kinner .eq. 1 )then
        CALL UM535 (mpos,mdpos)
        sit2 = SIN( mdpos%colat * deg )
        cot2 = COS( mdpos%colat * deg )
        rre2 = mdpos%radius / re
        CALL UM531 (rre2, sit2, cot2, sip2, cop2, bir, bit, bip, k531)
      else if( mint%kinner .lt. 4 )then
        CALL UM531 (rre1, sit1, cot1, sip2, cop2, bir, bit, bip, k531)
      else
        CALL UM537 (rre1, sit1, cot1, sip2, cop2, bir, bit, bip, k531)
      endif
C      
C  c/ Compute the external magnetic field
C
C
      ber  = 0.0d0
      bet  = 0.0d0
      bep  = 0.0d0
      if( mext%kouter .ne. 0 )then
        if( mext%kouter .eq. 5 .and. rre1 .gt. 15.0d0 )then
          k532 = -53201   
        else
          CALL UM532 (rre1, sit1, cot1, sip, cop, ber, bet, bep)
        endif
      endif
C
C  d/ Total
C
      mb%rho   = bir + ber
      mb%theta = bit + bet
      mb%phi   = bip + bep
      mb%dnrm  = SQRT( mb%rho**2 + mb%theta**2 +
     :                 mb%phi**2 )
C
      if( mb%dnrm .lt. xbmin ) then
        ifail = -53002
      elseif( mb%dnrm .gt. xbmax ) then
        ifail = -53004
      elseif( ABS(csgmla) .gt. xtmin ) then
        ifail = -53003
      endif
C
      if( ifail .ge. 0 )ifail=k531
      if( ifail .ge. 0 )ifail=k532
C
      END
C----------------------------------------------------------------------
      SUBROUTINE UM531
     :          (rkm, st, ct, sph, cph, br, bt, bp, ifail)
C
C!    Geomagnetic field evaluation
C                       
C     REFERENCES
C     Subroutine ALLMAG of bint.for in blxtra
C     E G STASSINOPOULOS AND G D MEAD (see below)
      IMPLICIT NONE
C
      INCLUDE 'structure.h'
C
C     INTERFACE
C
        REAL*8    rkm, st, ct, sph, cph, br, bt, bp
        INTEGER*4 ifail
C
C
      COMMON /UC140/ mint, mext, msun 
C
        TYPE(zemf) :: mext
        TYPE(zimf) :: mint
        TYPE(zsun) :: msun         
C                                              
C
C     HISTORY
C *  GEOCENTRIC VERSION OF GEOMAGNETIC FIELD ROUTINE
C *  DOUBLE PRECISION DECK FOR CII-HB / DPS05
C *  LONG DECK, THROUGH NMAX=13, FIXED INDICES WITHOUT DO LOOPS
C *  EXECUTION TIME PER CALL FACTOR OF THREE LESS THAN SHORT DECK
C *  PROGRAM DESIGNED AND TESTED BY E G STASSINOPOULOS AND G D MEAD,
C *  CODE 641, NASA GODDARD SPACE FLT CTR, GREENBELT, MD 20771
C *
C *  MODIFIED JUNE 82 BY L.JENSEN ESTEC/TMA: FIELD MODEL COEFFICIENTS
C *           ARE READ IN FROM DATA-FILE RATHER THAN STORED IN CORE.
C *  JUN 92: D. HEYNDERICKX (BIRA) : MOVED SCHMIDT NORMALISATION TO
C *                                   SUBROUTINE READC
C *                                  FIELD COEFFICIENTS AND NMAX IN 
C *                                   COMMON COEFF
C **          RKM        DISTANCE IN EARTH RADII TO EARTH'S CENTRE.
C **          ST,CT      SIN & COS OF GEOCENTRIC COLATITUDE
C **          SPH,CPH    SIN & COS OF EAST LONGITUDE
C **  OUTPUT: BR,BT,BP   GEOCENTRIC FIELD COMPONENTS IN GAUSS
C **    NOTE: FOR GREATEST EFFICIENCY, COMPLETE ALL CALCULATIONS WITH
C **            ONE MODEL AND ONE TIME BEFORE CHANGING MODELS OR TIME.
C
C
C     VARIABLES
C
        REAL*8 AR, AOR
        REAL*8 P21, P22, SP2, CP2, DP21, DP22, C2, SP3,
     :         CP3, P31, P32, P33
        REAL*8 DP31, DP32, DP33, C3, SP4, CP4, P41, P42,
     :         P43, P44, DP41
        REAL*8 DP42, DP43, DP44, C4, SP5, CP5, P51, P52,
     :         P53, P54, P55
        REAL*8 DP51, DP52, DP53, DP54, DP55, C5, SP6, CP6,
     :         P61, P62, P63
        REAL*8 P64, P65, P66, DP61, DP62, DP63, DP64, DP65,
     :         DP66, C6
        REAL*8 SP7, CP7, P71, P72, P73, P74, P75, P76, P77, 
     :         DP71, DP72
        REAL*8 DP73, DP74, DP75, DP76, DP77, C7, SP8, CP8,
     :         P81, P82, P83
        REAL*8 P84, P85, P86, P87, P88, DP81, DP82, DP83,
     :         DP84, DP85, DP86
        REAL*8 DP87, DP88, C8, SP9, CP9, P91, P92, P93, P94,
     :         P95, P96, P97
        REAL*8 P98, P99, DP91, DP92, DP93, DP94, DP95, DP96,
     :         DP97, DP98
        REAL*8 DP99, C9, SP10, CP10, P101, P102, P103, P104,
     :         P105, P106
        REAL*8 P107, P108, P109, P1010, DP101, DP102, DP103,
     :         DP104, DP105
        REAL*8 DP106, DP107, DP108, DP109, DP1010, C10, SP11,
     :         CP11, P111
        REAL*8 P112, P113, P114, P115, P116, P117, P118, P119,
     :         P1110
        REAL*8 P1111, DP111, DP112, DP113, DP114, DP115, DP116,
     :         DP117
        REAL*8 DP118, DP119, DP1110, DP1111, C11, SP12, CP12, 
     :         P121, P122
        REAL*8 P123, P124, P125, P126, P127, P128, P129, P1210, 
     :         P1211
        REAL*8 P1212, DP121, DP122, DP123, DP124, DP125, DP126,
     :         DP127
        REAL*8 DP128, DP129, DP1210, DP1211, DP1212, C12, SP13,
     :         CP13, P131
        REAL*8 P132, P133, P134, P135, P136, P137, P138, P139,
     :         P1310
        REAL*8 P1311, P1312, P1313, DP131, DP132, DP133, DP134,
     :         DP135
        REAL*8 DP136, DP137, DP138, DP139, DP1310, DP1311, DP1312,
     :         DP1313, C13
        REAL*8 SP14, CP14, P141, P142, P143, P144, P145, 
     :         P146, P147, P148, P149, P1410, P1411, P1412,
     :         P1413, P1414, DP141, DP142, DP143, DP144, 
     :         DP146, DP147, DP148, DP149, DP1410,DP145,
     :         DP1411, DP1412, DP1413, DP1414, C14
        REAL*8 G(nx140,nx140)
        INTEGER*4 i, j
C
C     CODE
C
      ifail=mint%norder
C
      do i=1,ifail
        do j=1,ifail
          G(i,j) = mint%coef(i,j)
        enddo
      enddo
C                                                           N= 2
      P21=CT
      P22=ST
      AR=1.0D0/RKM
      SP2=SPH
      CP2=CPH
      DP21=-P22
      DP22=P21
      AOR=AR*AR*AR
      C2=G(2,2)*CP2+G(1,2)*SP2
      BR=-(AOR+AOR)*(G(2,1)*P21+C2*P22)
      BT=AOR*(G(2,1)*DP21+C2*DP22)
      BP=AOR*(G(1,2)*CP2-G(2,2)*SP2)*P22
      IF (mint%norder .LE. 2) GO TO 1
C                                                           N= 3
      SP3=(SP2+SP2)*CP2
      CP3=(CP2+SP2)*(CP2-SP2)
      P31=P21*P21-0.333333333D0
      P32=P21*P22
      P33=P22*P22
      DP31=-P32-P32
      DP32=P21*P21-P33
      DP33=-DP31
      AOR=AOR*AR
      C2=G(3,2)*CP2+G(1,3)*SP2
      C3=G(3,3)*CP3+G(2,3)*SP3
      BR=BR-3.0D0*AOR*(G(3,1)*P31+C2*P32+C3*P33)
      BT=BT+AOR*(G(3,1)*DP31+C2*DP32+C3*DP33)
      BP=BP-AOR*((G(3,2)*SP2-G(1,3)*CP2)*P32+2.0D0*(G(3,3)*SP3-G(2,3)*
     :   CP3)*P33)
      IF (mint%norder .LE. 3) GO TO 1
C                                                           N= 4
      SP4=SP2*CP3+CP2*SP3
      CP4=CP2*CP3-SP2*SP3
      P41=P21*P31-0.26666666D0*P21
      DP41=P21*DP31+DP21*P31-0.26666666D0*DP21
      P42=P21*P32-0.2D0*P22
      DP42=P21*DP32+DP21*P32-0.2D0*DP22
      P43=P21*P33
      DP43=P21*DP33+DP21*P33
      P44=P22*P33
      DP44=3.0D0*P43
      AOR=AOR*AR
      C2=G(4,2)*CP2+G(1,4)*SP2
      C3=G(4,3)*CP3+G(2,4)*SP3
      C4=G(4,4)*CP4+G(3,4)*SP4
      BR=BR-4.0D0*AOR*(G(4,1)*P41+C2*P42+C3*P43+C4*P44)
      BT=BT+AOR*(G(4,1)*DP41+C2*DP42+C3*DP43+C4*DP44)
      BP=BP-AOR*((G(4,2)*SP2-G(1,4)*CP2)*P42+2.0D0*(G(4,3)*SP3-G(2,4)*
     :   CP3)*P43+3.0D0*(G(4,4)*SP4-G(3,4)*CP4)*P44)
      IF (mint%norder .LE. 4) GO TO 1
C                                                           N= 5
      SP5=(SP3+SP3)*CP3
      CP5=(CP3+SP3)*(CP3-SP3)
      P51=P21*P41-0.25714285D0*P31
      DP51=P21*DP41+DP21*P41-0.25714285D0*DP31
      P52=P21*P42-0.22857142D0*P32
      DP52=P21*DP42+DP21*P42-0.22857142D0*DP32
      P53=P21*P43-0.14285714D0*P33
      DP53=P21*DP43+DP21*P43-0.14285714D0*DP33
      P54=P21*P44
      DP54=P21*DP44+DP21*P44
      P55=P22*P44
      DP55=4.0D0*P54
      AOR=AOR*AR
      C2=G(5,2)*CP2+G(1,5)*SP2
      C3=G(5,3)*CP3+G(2,5)*SP3
      C4=G(5,4)*CP4+G(3,5)*SP4
      C5=G(5,5)*CP5+G(4,5)*SP5
      BR=BR-5.0D0*AOR*(G(5,1)*P51+C2*P52+C3*P53+C4*P54+C5*P55)
      BT=BT+AOR*(G(5,1)*DP51+C2*DP52+C3*DP53+C4*DP54+C5*DP55)
      BP=BP-AOR*((G(5,2)*SP2-G(1,5)*CP2)*P52+2.0D0*(G(5,3)*SP3-G(2,5)*
     :   CP3)*P53+3.0D0*(G(5,4)*SP4-G(3,5)*CP4)*P54+4.0D0*(G(5,5)*SP5-
     :   G(4,5)*CP5)*P55)
      IF (mint%norder .LE. 5) GO TO 1
C                                                           N= 6
      SP6=SP2*CP5+CP2*SP5
      CP6=CP2*CP5-SP2*SP5
      P61=P21*P51-0.25396825D0*P41
      DP61=P21*DP51+DP21*P51-0.25396825D0*DP41
      P62=P21*P52-0.23809523D0*P42
      DP62=P21*DP52+DP21*P52-0.23809523D0*DP42
      P63=P21*P53-0.19047619D0*P43
      DP63=P21*DP53+DP21*P53-0.19047619D0*DP43
      P64=P21*P54-0.11111111D0*P44
      DP64=P21*DP54+DP21*P54-0.11111111D0*DP44
      P65=P21*P55
      DP65=P21*DP55+DP21*P55
      P66=P22*P55
      DP66=5.0D0*P65
      AOR=AOR*AR
      C2=G(6,2)*CP2+G(1,6)*SP2
      C3=G(6,3)*CP3+G(2,6)*SP3
      C4=G(6,4)*CP4+G(3,6)*SP4
      C5=G(6,5)*CP5+G(4,6)*SP5
      C6=G(6,6)*CP6+G(5,6)*SP6
      BR=BR-6.0D0*AOR*(G(6,1)*P61+C2*P62+C3*P63+C4*P64+C5*P65+C6*P66)
      BT=BT+AOR*(G(6,1)*DP61+C2*DP62+C3*DP63+C4*DP64+C5*DP65+C6*DP66)
      BP=BP-AOR*((G(6,2)*SP2-G(1,6)*CP2)*P62+2.0D0*(G(6,3)*SP3-G(2,6)*
     :   CP3)*P63+3.0D0*(G(6,4)*SP4-G(3,6)*CP4)*P64+4.0D0*(G(6,5)*SP5-
     :   G(4,6)*CP5)*P65+5.0D0*(G(6,6)*SP6-G(5,6)*CP6)*P66)
      IF (mint%norder .LE. 6) GO TO 1
C                                                           N= 7
      SP7=(SP4+SP4)*CP4
      CP7=(CP4+SP4)*(CP4-SP4)
      P71=P21*P61-0.25252525D0*P51
      DP71=P21*DP61+DP21*P61-0.25252525D0*DP51
      P72=P21*P62-0.24242424D0*P52
      DP72=P21*DP62+DP21*P62-0.24242424D0*DP52
      P73=P21*P63-0.21212121D0*P53
      DP73=P21*DP63+DP21*P63-0.21212121D0*DP53
      P74=P21*P64-0.16161616D0*P54
      DP74=P21*DP64+DP21*P64-0.16161616D0*DP54
      P75=P21*P65-0.09090909D0*P55
      DP75=P21*DP65+DP21*P65-0.09090909D0*DP55
      P76=P21*P66
      DP76=P21*DP66+DP21*P66
      P77=P22*P66
      DP77=6.0D0*P76
      AOR=AOR*AR
      C2=G(7,2)*CP2+G(1,7)*SP2
      C3=G(7,3)*CP3+G(2,7)*SP3
      C4=G(7,4)*CP4+G(3,7)*SP4
      C5=G(7,5)*CP5+G(4,7)*SP5
      C6=G(7,6)*CP6+G(5,7)*SP6
      C7=G(7,7)*CP7+G(6,7)*SP7
      BR=BR-7.0D0*AOR*(G(7,1)*P71+C2*P72+C3*P73+C4*P74+C5*P75+C6*P76+
     :   C7*P77)
      BT=BT+AOR*(G(7,1)*DP71+C2*DP72+C3*DP73+C4*DP74+C5*DP75+C6*DP76+C7*
     :   DP77)
      BP=BP-AOR*((G(7,2)*SP2-G(1,7)*CP2)*P72+2.0D0*(G(7,3)*SP3-G(2,7)*
     :   CP3)*P73+3.0D0*(G(7,4)*SP4-G(3,7)*CP4)*P74+4.0D0*(G(7,5)*SP5-
     :   G(4,7)*CP5)*P75+5.0D0*(G(7,6)*SP6-G(5,7)*CP6)*P76+6.0D0*
     :   (G(7,7)*SP7-G(6,7)*CP7)*P77)
      IF (mint%norder .LE. 7) GO TO 1
C                                                           N= 8
      SP8=SP2*CP7+CP2*SP7
      CP8=CP2*CP7-SP2*SP7
      P81=P21*P71-0.25174825D0*P61
      DP81=P21*DP71+DP21*P71-0.25174825D0*DP61
      P82=P21*P72-0.24475524D0*P62
      DP82=P21*DP72+DP21*P72-0.24475524D0*DP62
      P83=P21*P73-0.22377622D0*P63
      DP83=P21*DP73+DP21*P73-0.22377622D0*DP63
      P84=P21*P74-0.18881118D0*P64
      DP84=P21*DP74+DP21*P74-0.18881118D0*DP64
      P85=P21*P75-0.13986013D0*P65
      DP85=P21*DP75+DP21*P75-0.13986013D0*DP65
      P86=P21*P76-0.07692307D0*P66
      DP86=P21*DP76+DP21*P76-0.07692307D0*DP66
      P87=P21*P77
      DP87=P21*DP77+DP21*P77
      P88=P22*P77
      DP88=7.0D0*P87
      AOR=AOR*AR
      C2=G(8,2)*CP2+G(1,8)*SP2
      C3=G(8,3)*CP3+G(2,8)*SP3
      C4=G(8,4)*CP4+G(3,8)*SP4
      C5=G(8,5)*CP5+G(4,8)*SP5
      C6=G(8,6)*CP6+G(5,8)*SP6
      C7=G(8,7)*CP7+G(6,8)*SP7
      C8=G(8,8)*CP8+G(7,8)*SP8
      BR=BR-8.0D0*AOR*(G(8,1)*P81+C2*P82+C3*P83+C4*P84+C5*P85+C6*P86+
     :   C7*P87+C8*P88)
      BT=BT+AOR*(G(8,1)*DP81+C2*DP82+C3*DP83+C4*DP84+C5*DP85+C6*DP86+C7*
     :   DP87+C8*DP88)
      BP=BP-AOR*((G(8,2)*SP2-G(1,8)*CP2)*P82+2.0D0*(G(8,3)*SP3-G(2,8)*
     :   CP3)*P83+3.0D0*(G(8,4)*SP4-G(3,8)*CP4)*P84+4.0D0*(G(8,5)*SP5-
     :   G(4,8)*CP5)*P85+5.0D0*(G(8,6)*SP6-G(5,8)*CP6)*P86+6.0D0*
     :   (G(8,7)*SP7-G(6,8)*CP7)*P87+7.0D0*(G(8,8)*SP8-G(7,8)*CP8)*P88)
      IF (mint%norder .LE. 8) GO TO 1
C                                                           N= 9
      SP9=(SP5+SP5)*CP5
      CP9=(CP5+SP5)*(CP5-SP5)
      P91=P21*P81-0.25128205D0*P71
      DP91=P21*DP81+DP21*P81-0.25128205D0*DP71
      P92=P21*P82-0.24615384D0*P72
      DP92=P21*DP82+DP21*P82-0.24615384D0*DP72
      P93=P21*P83-0.23076923D0*P73
      DP93=P21*DP83+DP21*P83-0.23076923D0*DP73
      P94=P21*P84-0.20512820D0*P74
      DP94=P21*DP84+DP21*P84-0.20512820D0*DP74
      P95=P21*P85-0.16923076D0*P75
      DP95=P21*DP85+DP21*P85-0.16923076D0*DP75
      P96=P21*P86-0.12307692D0*P76
      DP96=P21*DP86+DP21*P86-0.12307692D0*DP76
      P97=P21*P87-0.06666666D0*P77
      DP97=P21*DP87+DP21*P87-0.06666666D0*DP77
      P98=P21*P88
      DP98=P21*DP88+DP21*P88
      P99=P22*P88
      DP99=8.0D0*P98
      AOR=AOR*AR
      C2=G(9,2)*CP2+G(1,9)*SP2
      C3=G(9,3)*CP3+G(2,9)*SP3
      C4=G(9,4)*CP4+G(3,9)*SP4
      C5=G(9,5)*CP5+G(4,9)*SP5
      C6=G(9,6)*CP6+G(5,9)*SP6
      C7=G(9,7)*CP7+G(6,9)*SP7
      C8=G(9,8)*CP8+G(7,9)*SP8
      C9=G(9,9)*CP9+G(8,9)*SP9
      BR=BR-9.0D0*AOR*(G(9,1)*P91+C2*P92+C3*P93+C4*P94+C5*P95+C6*P96+
     :   C7*P97+C8*P98+C9*P99)
      BT=BT+AOR*(G(9,1)*DP91+C2*DP92+C3*DP93+C4*DP94+C5*DP95+C6*DP96+C7*
     :   DP97+C8*DP98+C9*DP99)
      BP=BP-AOR*((G(9,2)*SP2-G(1,9)*CP2)*P92+2.0D0*(G(9,3)*SP3-G(2,9)*
     :   CP3)*P93+3.0D0*(G(9,4)*SP4-G(3,9)*CP4)*P94+4.0d0*(G(9,5)*SP5-
     :   G(4,9)*CP5)*P95+5.0D0*(G(9,6)*SP6-G(5,9)*CP6)*P96+6.0D0*
     :   (G(9,7)*SP7-G(6,9)*CP7)*P97+7.0D0*(G(9,8)*SP8-G(7,9)*CP8)*P98+
     :   8.0D0*(G(9,9)*SP9-G(8,9)*CP9)*P99)
      IF (mint%norder .LE. 9) GO TO 1
C                                                           N=10
      SP10=SP2*CP9+CP2*SP9
      CP10=CP2*CP9-SP2*SP9
      P101=P21*P91-0.25098039D0*P81
      DP101=P21*DP91+DP21*P91-0.25098039D0*DP81
      P102=P21*P92-0.24705882D0*P82
      DP102=P21*DP92+DP21*P92-0.24705882D0*DP82
      P103=P21*P93-0.23529411D0*P83
      DP103=P21*DP93+DP21*P93-0.23529411D0*DP83
      P104=P21*P94-0.21568627D0*P84
      DP104=P21*DP94+DP21*P94-0.21568627D0*DP84
      P105=P21*P95-0.18823529D0*P85
      DP105=P21*DP95+DP21*P95-0.18823529D0*DP85
      P106=P21*P96-0.15294117D0*P86
      DP106=P21*DP96+DP21*P96-0.15294117D0*DP86
      P107=P21*P97-0.10980392D0*P87
      DP107=P21*DP97+DP21*P97-0.10980392D0*DP87
      P108=P21*P98-0.05882352D0*P88
      DP108=P21*DP98+DP21*P98-0.05882352D0*DP88
      P109=P21*P99
      DP109=P21*DP99+DP21*P99
      P1010=P22*P99
      DP1010=9.0D0*P109
      AOR=AOR*AR
      C2=G(10,2)*CP2+G(1,10)*SP2
      C3=G(10,3)*CP3+G(2,10)*SP3
      C4=G(10,4)*CP4+G(3,10)*SP4
      C5=G(10,5)*CP5+G(4,10)*SP5
      C6=G(10,6)*CP6+G(5,10)*SP6
      C7=G(10,7)*CP7+G(6,10)*SP7
      C8=G(10,8)*CP8+G(7,10)*SP8
      C9=G(10,9)*CP9+G(8,10)*SP9
      C10=G(10,10)*CP10+G(9,10)*SP10
      BR=BR-10.0D0*AOR*(G(10,1)*P101+C2*P102+C3*P103+C4*P104+C5*P105+
     :   C6*P106+C7*P107+C8*P108+C9*P109+C10*P1010)
      BT=BT+AOR*(G(10,1)*DP101+C2*DP102+C3*DP103+C4*DP104+C5*DP105+C6*
     :   DP106+C7*DP107+C8*DP108+C9*DP109+C10*DP1010)
      BP=BP-AOR*((G(10,2)*SP2-G(1,10)*CP2)*P102+2.0D0*(G(10,3)*SP3-
     :   G(2,10)*CP3)*P103+3.0D0*(G(10,4)*SP4-G(3,10)*CP4)*P104+4.0D0*
     :   (G(10,5)*SP5-G(4,10)*CP5)*P105+5.0D0*(G(10,6)*SP6-G(5,10)*CP6)*
     :   P106+6.0D0*(G(10,7)*SP7-G(6,10)*CP7)*P107+7.0D0*(G(10,8)*SP8-
     :   G(7,10)*CP8)*P108+8.0D0*(G(10,9)*SP9-G(8,10)*CP9)*P109+9.0D0*
     :   (G(10,10)*SP10-G(9,10)*CP10)*P1010)
      IF (mint%norder .LE. 10) GO TO 1
C                                                           N=11
      SP11=(SP6+SP6)*CP6
      CP11=(CP6+SP6)*(CP6-SP6)
      P111=P21*P101-0.25077399D0*P91
      DP111=P21*DP101+DP21*P101-0.25077399D0*DP91
      P112=P21*P102-0.24767801D0*P92
      DP112=P21*DP102+DP21*P102-0.24767801D0*DP92
      P113=P21*P103-0.23839009D0*P93
      DP113=P21*DP103+DP21*P103-0.23839009D0*DP93
      P114=P21*P104-0.22291021D0*P94
      DP114=P21*DP104+DP21*P104-0.22291021D0*DP94
      P115=P21*P105-0.20123839D0*P95
      DP115=P21*DP105+DP21*P105-0.20123839D0*DP95
      P116=P21*P106-0.17337461D0*P96
      DP116=P21*DP106+DP21*P106-0.17337461D0*DP96
      P117=P21*P107-0.13931888D0*P97
      DP117=P21*DP107+DP21*P107-0.13931888D0*DP97
      P118=P21*P108-0.09907120D0*P98
      DP118=P21*DP108+DP21*P108-0.09907120D0*DP98
      P119=P21*P109-0.05263157D0*P99
      DP119=P21*DP109+DP21*P109-0.05263157D0*DP99
      P1110=P21*P1010
      DP1110=P21*DP1010+DP21*P1010
      P1111=P22*P1010
      DP1111=10.0D0*P1110
      AOR=AOR*AR
      C2=G(11,2)*CP2+G(1,11)*SP2
      C3=G(11,3)*CP3+G(2,11)*SP3
      C4=G(11,4)*CP4+G(3,11)*SP4
      C5=G(11,5)*CP5+G(4,11)*SP5
      C6=G(11,6)*CP6+G(5,11)*SP6
      C7=G(11,7)*CP7+G(6,11)*SP7
      C8=G(11,8)*CP8+G(7,11)*SP8
      C9=G(11,9)*CP9+G(8,11)*SP9
      C10=G(11,10)*CP10+G(9,11)*SP10
      C11=G(11,11)*CP11+G(10,11)*SP11
      BR=BR-11.0D0*AOR*(G(11,1)*P111+C2*P112+C3*P113+C4*P114+C5*P115+
     :   C6*P116+C7*P117+C8*P118+C9*P119+C10*P1110+C11*P1111)
      BT=BT+AOR*(G(11,1)*DP111+C2*DP112+C3*DP113+C4*DP114+C5*DP115+C6*
     :   DP116+C7*DP117+C8*DP118+C9*DP119+C10*DP1110+C11*DP1111)
      BP=BP-AOR*((G(11,2)*SP2-G(1,11)*CP2)*P112+2.0D0*(G(11,3)*SP3-
     :   G(2,11)*CP3)*P113+3.0D0*(G(11,4)*SP4-G(3,11)*CP4)*P114+4.0D0*
     :   (G(11,5)*SP5-G(4,11)*CP5)*P115+5.0D0*(G(11,6)*SP6-G(5,11)*CP6)*
     :   P116+6.0D0*(G(11,7)*SP7-G(6,11)*CP7)*P117+7.0D0*(G(11,8)*SP8-
     :   G(7,11)*CP8)*P118+8.0D0*(G(11,9)*SP9-G(8,11)*CP9)*P119+9.0D0*
     :   (G(11,10)*SP10-G(9,11)*CP10)*P1110+10.0D0*(G(11,11)*SP11-
     :   G(10,11)*CP11)*P1111)
      IF (mint%norder .LE. 11) GO TO 1
C                                                           N=12
      SP12=SP2*CP11+CP2*SP11
      CP12=CP2*CP11-SP2*SP11
      P121=P21*P111-0.25062656D0*P101
      DP121=P21*DP111+DP21*P111-0.25062656D0*DP101
      P122=P21*P112-0.24812030D0*P102
      DP122=P21*DP112+DP21*P112-0.24812030D0*DP102
      P123=P21*P113-0.24060150D0*P103
      DP123=P21*DP113+DP21*P113-0.24060150D0*DP103
      P124=P21*P114-0.22807017D0*P104
      DP124=P21*DP114+DP21*P114-0.22807017D0*DP104
      P125=P21*P115-0.21052631D0*P105
      DP125=P21*DP115+DP21*P115-0.21052631D0*DP105
      P126=P21*P116-0.18796992D0*P106
      DP126=P21*DP116+DP21*P116-0.18796992D0*DP106
      P127=P21*P117-0.16040100D0*P107
      DP127=P21*DP117+DP21*P117-0.16040100D0*DP107
      P128=P21*P118-0.12781954D0*P108
      DP128=P21*DP118+DP21*P118-0.12781954D0*DP108
      P129=P21*P119-0.09022556D0*P109
      DP129=P21*DP119+DP21*P119-0.09022556D0*DP109
      P1210=P21*P1110-0.04761904D0*P1010
      DP1210=P21*DP1110+DP21*P1110-0.04761904D0*DP1010
      P1211=P21*P1111
      DP1211=P21*DP1111+DP21*P1111
      P1212=P22*P1111
      DP1212=11.0D0*P1211
      AOR=AOR*AR
      C2=G(12,2)*CP2+G(1,12)*SP2
      C3=G(12,3)*CP3+G(2,12)*SP3
      C4=G(12,4)*CP4+G(3,12)*SP4
      C5=G(12,5)*CP5+G(4,12)*SP5
      C6=G(12,6)*CP6+G(5,12)*SP6
      C7=G(12,7)*CP7+G(6,12)*SP7
      C8=G(12,8)*CP8+G(7,12)*SP8
      C9=G(12,9)*CP9+G(8,12)*SP9
      C10=G(12,10)*CP10+G(9,12)*SP10
      C11=G(12,11)*CP11+G(10,12)*SP11
      C12=G(12,12)*CP12+G(11,12)*SP12
      BR=BR-12.0D0*AOR*(G(12,1)*P121+C2*P122+C3*P123+C4*P124+C5*P125+
     :   C6*P126+C7*P127+C8*P128+C9*P129+C10*P1210+C11*P1211+C12*P1212)
      BT=BT+AOR*(G(12,1)*DP121+C2*DP122+C3*DP123+C4*DP124+C5*DP125+C6*
     :   DP126+C7*DP127+C8*DP128+C9*DP129+C10*DP1210+C11*DP1211+C12*
     :   DP1212)
      BP=BP-AOR*((G(12,2)*SP2-G(1,12)*CP2)*P122+2.0D0*(G(12,3)*SP3-
     :   G(2,12)*CP3)*P123+3.0D0*(G(12,4)*SP4-G(3,12)*CP4)*P124+4.0D0*
     :   (G(12,5)*SP5-G(4,12)*CP5)*P125+5.0D0*(G(12,6)*SP6-G(5,12)*CP6)*
     :   P126+6.0D0*(G(12,7)*SP7-G(6,12)*CP7)*P127+7.0D0*(G(12,8)*SP8-
     :   G(7,12)*CP8)*P128+8.0D0*(G(12,9)*SP9-G(8,12)*CP9)*P129+9.0D0*
     :   (G(12,10)*SP10-G(9,12)*CP10)*P1210+10.0D0*(G(12,11)*SP11-
     :   G(10,12)*CP11)*P1211+11.0D0*(G(12,12)*SP12-G(11,12)*CP12)*
     :   P1212)
      IF (mint%norder .LE. 12) GO TO 1
C                                                           N=13
      SP13=(SP7+SP7)*CP7
      CP13=(CP7+SP7)*(CP7-SP7)
      P131=P21*P121-0.25051759D0*P111
      DP131=P21*DP121+DP21*P121-0.25051759D0*DP111
      P132=P21*P122-0.24844720D0*P112
      DP132=P21*DP122+DP21*P122-0.24844720D0*DP112
      P133=P21*P123-0.24223602D0*P113
      DP133=P21*DP123+DP21*P123-0.24223602D0*DP113
      P134=P21*P124-0.23188405D0*P114
      DP134=P21*DP124+DP21*P124-0.23188405D0*DP114
      P135=P21*P125-0.21739130D0*P115
      DP135=P21*DP125+DP21*P125-0.21739130D0*DP115
      P136=P21*P126-0.19875776D0*P116
      DP136=P21*DP126+DP21*P126-0.19875776D0*DP116
      P137=P21*P127-0.17598343D0*P117
      DP137=P21*DP127+DP21*P127-0.17598343D0*DP117
      P138=P21*P128-0.14906832D0*P118
      DP138=P21*DP128+DP21*P128-0.14906832D0*DP118
      P139=P21*P129-0.11801242D0*P119
      DP139=P21*DP129+DP21*P129-0.11801242D0*DP119
      P1310=P21*P1210-0.08281573D0*P1110
      DP1310=P21*DP1210+DP21*P1210-0.08281573D0*DP1110
      P1311=P21*P1211-0.04347826D0*P1111
      DP1311=P21*DP1211+DP21*P1211-0.04347826D0*DP1111
      P1312=P21*P1212
      DP1312=P21*DP1212+DP21*P1212
      P1313=P22*P1212
      DP1313=12.0D0*P1312
      AOR=AOR*AR
      C2=G(13,2)*CP2+G(1,13)*SP2
      C3=G(13,3)*CP3+G(2,13)*SP3
      C4=G(13,4)*CP4+G(3,13)*SP4
      C5=G(13,5)*CP5+G(4,13)*SP5
      C6=G(13,6)*CP6+G(5,13)*SP6
      C7=G(13,7)*CP7+G(6,13)*SP7
      C8=G(13,8)*CP8+G(7,13)*SP8
      C9=G(13,9)*CP9+G(8,13)*SP9
      C10=G(13,10)*CP10+G(9,13)*SP10
      C11=G(13,11)*CP11+G(10,13)*SP11
      C12=G(13,12)*CP12+G(11,13)*SP12
      C13=G(13,13)*CP13+G(12,13)*SP13
      BR=BR-13.0D0*AOR*(G(13,1)*P131+C2*P132+C3*P133+C4*P134+C5*P135+
     :   C6*P136+C7*P137+C8*P138+C9*P139+C10*P1310+C11*P1311+C12*P1312+
     :   C13*P1313)
      BT=BT+AOR*(G(13,1)*DP131+C2*DP132+C3*DP133+C4*DP134+C5*DP135+C6*
     :   DP136+C7*DP137+C8*DP138+C9*DP139+C10*DP1310+C11*DP1311+C12*
     :   DP1312+C13*DP1313)
      BP=BP-AOR*((G(13,2)*SP2-G(1,13)*CP2)*P132+2.0D0*(G(13,3)*SP3-
     :   G(2,13)*CP3)*P133+3.0D0*(G(13,4)*SP4-G(3,13)*CP4)*P134+4.0D0*
     :   (G(13,5)*SP5-G(4,13)*CP5)*P135+5.0D0*(G(13,6)*SP6-G(5,13)*CP6)*
     :   P136+6.0D0*(G(13,7)*SP7-G(6,13)*CP7)*P137+7.0D0*(G(13,8)*SP8-
     :   G(7,13)*CP8)*P138+8.0D0*(G(13,9)*SP9-G(8,13)*CP9)*P139+9.0D0*
     :   (G(13,10)*SP10-G(9,13)*CP10)*P1310+10.0D0*(G(13,11)*SP11-
     :   G(10,13)*CP11)*P1311+11.0D0*(G(13,12)*SP12-G(11,13)*CP12)*
     :   P1312+12.0D0*(G(13,13)*SP13-G(12,13)*CP13)*P1313)
      IF (mint%norder .LE. 13) GO TO 1
CCC Code for order 14 added by D. Heynderickx (DH Consultancy)
CCC 25 May 2015
C                                                           N=14
      SP14=SP2*CP13+CP2*SP13
      CP14=CP2*CP13-SP2*SP13
      P141=P21*P131-0.25043478D0*P121
      DP141=P21*DP131+DP21*P131-0.25043478D0*DP121
      P142=P21*P132-0.24869565D0*P122
      DP142=P21*DP132+DP21*P132-0.24869565D0*DP122
      P143=P21*P133-0.24347826D0*P123
      DP143=P21*DP133+DP21*P133-0.24347826D0*DP123
      P144=P21*P134-0.23478261D0*P124
      DP144=P21*DP134+DP21*P134-0.23478261D0*DP124
      P145=P21*P135-0.22260870D0*P125
      DP145=P21*DP135+DP21*P135-0.22260870D0*DP125
      P146=P21*P136-0.20695652D0*P126
      DP146=P21*DP136+DP21*P136-0.20695652D0*DP126
      P147=P21*P137-0.18782609D0*P127
      DP147=P21*DP137+DP21*P137-0.18782609D0*DP127
      P148=P21*P138-0.16521739D0*P128
      DP148=P21*DP138+DP21*P138-0.16521739D0*DP128
      P149=P21*P139-0.13913044D0*P129
      DP149=P21*DP139+DP21*P139-0.13913044D0*DP129
      P1410=P21*P1310-0.10956522D0*P1210
      DP1410=P21*DP1310+DP21*P1310-0.10956522D0*DP1210
      P1411=P21*P1311-0.07652174D0*P1211
      DP1411=P21*DP1311+DP21*P1311-0.07652174D0*DP1211
      P1412=P21*P1312-0.04000000D0*P1212
      DP1412=P21*DP1312+DP21*P1312-0.04000000D0*DP1212
      P1413=P21*P1313
      DP1413=P21*DP1313+DP21*P1313
      P1414=P22*P1313
      DP1414=13.0D0*P1413
      AOR=AOR*AR
      C2=G(14,2)*CP2+G(1,14)*SP2
      C3=G(14,3)*CP3+G(2,14)*SP3
      C4=G(14,4)*CP4+G(3,14)*SP4
      C5=G(14,5)*CP5+G(4,14)*SP5
      C6=G(14,6)*CP6+G(5,14)*SP6
      C7=G(14,7)*CP7+G(6,14)*SP7
      C8=G(14,8)*CP8+G(7,14)*SP8
      C9=G(14,9)*CP9+G(8,14)*SP9
      C10=G(14,10)*CP10+G(9,14)*SP10
      C11=G(14,11)*CP11+G(10,14)*SP11
      C12=G(14,12)*CP12+G(11,14)*SP12
      C13=G(14,13)*CP13+G(12,14)*SP13
      C14=G(14,14)*CP14+G(13,14)*SP14
      BR=BR-14.0D0*AOR*(G(14,1)*P141+C2*P142+C3*P143+C4*P144+C5*P145+
     :   C6*P146+C7*P147+C8*P148+C9*P149+C10*P1410+C11*P1411+C12*P1412+
     :   C13*P1413+C14*P1414)
      BT=BT+AOR*(G(14,1)*DP141+C2*DP142+C3*DP143+C4*DP144+C5*DP145+C6*
     :   DP146+C7*DP147+C8*DP148+C9*DP149+C10*DP1410+C11*DP1411+C12*
     :   DP1412+C13*DP1413+C14*DP1414)
      BP=BP-AOR*((G(14,2)*SP2-G(1,14)*CP2)*P142+2.0D0*(G(14,3)*SP3-
     :   G(2,14)*CP3)*P143+3.0D0*(G(14,4)*SP4-G(3,14)*CP4)*P144+4.0D0*
     :   (G(14,5)*SP5-G(4,14)*CP5)*P145+5.0D0*(G(14,6)*SP6-G(5,14)*CP6)*
     :   P146+6.0D0*(G(14,7)*SP7-G(6,14)*CP7)*P147+7.0D0*(G(14,8)*SP8-
     :   G(7,14)*CP8)*P148+8.0D0*(G(14,9)*SP9-G(8,14)*CP9)*P149+9.0D0*
     :   (G(14,10)*SP10-G(9,14)*CP10)*P1410+10.0D0*(G(14,11)*SP11-
     :   G(10,14)*CP11)*P1411+11.0D0*(G(14,12)*SP12-G(11,14)*CP12)*
     :   P1412+12.0D0*(G(14,13)*SP13-G(12,14)*CP13)*P1413+13.0D0*
     :   (G(14,14)*SP14-G(13,14)*CP14)*P1414)
      IF (mint%norder .LE. 14) GO TO 1
      ifail=-53101
C
    1 BR = BR / 1.0D5
      BT = BT / 1.0D5
      if( abs(st).lt.1.e-6 )then
        ifail=-53102
        BP = 0
      else       
        BP = BP / ST / 1.0D5
      endif       
C
C
      END
C----------------------------------------------------------------------
      SUBROUTINE UM532
     :          (r, sla, cla, slo, clo, dbr, dbt, dbp)
C
C!    External magnetic field evaluation
C                       
C     REFERENCES
C     File bext.for in blxtra
C     Subroutines trans3 and trans4 in transfos.for
C
      INCLUDE 'structure.h'
C
      EXTERNAL  mead, tsy87s, tsy87l, t89c, t96_01, t01_01, t04_s, 
     &          bom97, a2000, bxyzmu, bdyn
C      
C     INTERFACE
C
        REAL*8        r, sla, cla, slo, clo, dbr, dbt, dbp     
C
      COMMON /UC140/  mint, mext, msun
C
        TYPE(zimf) :: mint
        TYPE(zemf) :: mext
        TYPE(zsun) :: msun
        TYPE(zdat) :: adate

      COMMON /ALEXEEV/ adate

C
      COMMON /UC160/  pi, deg, re, gmagmo, eclipt, geoid, uma
C
        REAL*8        pi, deg, re, gmagmo, eclipt, geoid(3), uma(30)
C                                                          
C     VARIABLES
C
        REAL*8        gx, gy, gz, ax, ay, az
        REAL*8        ps, abx, aby, abz
        REAL*8        bbx, bby, bbz
        REAL*8        t(3,3),ptsy(10)
        REAL*8        dvc1(3), dvc2(3)
        REAL*8        ut, xa(3), bma(3), bba(7,3)
        INTEGER*4     i,j
        INTEGER       year, month, day
C
C     CODE
C
      do i=1,3
        do j=1,3
          t(i,j) = mext%trans(i,j)
        enddo
      enddo
C        
      ps = mext%tilt * deg
C      
C  a/ CONVERTS GEO COORDINATES TO SM OR GSM
C
      gx = r * sla * clo
      gy = r * sla * slo
      gz = r * cla
      ax = t(1,1) * gx + t(1,2) * gy + t(1,3) * gz
      ay = t(2,1) * gx + t(2,2) * gy + t(2,3) * gz
      az = t(3,1) * gx + t(3,2) * gy + t(3,3) * gz
C      
c  b/ CALL THE MODEL
C
      if     ( mext%kouter .eq. 1 ) then
        call mead(ax, ay, az, mext%tilt, mext%ikp, abx, aby, abz)
      elseif ( mext%kouter .eq. 2 ) then
        call tsy87s(mext%ikp, ps, ax, ay, az, abx, aby, abz)
      elseif ( mext%kouter .eq. 3 ) then
        call tsy87l(mext%ikp, ps, ax, ay, az, abx, aby, abz)
      elseif ( mext%kouter .eq. 4 ) then
        call t89c(mext%ikp, ptsy, ps, ax, ay, az, abx, aby, abz)
      elseif ( mext%kouter .eq. 5 ) then
        call bxyzmu(ax, ay, az, mext%tilt, abx, aby, abz)
      elseif ( mext%kouter .eq. 6 ) then
        call bdyn(mext%wdens, mext%wvel, mext%vdst, 
     :            ax, ay, az, abx, aby, abz)
      elseif ( mext%kouter .eq. 7 ) then
        ptsy(1) = mext%pdyn
        ptsy(2) = mext%vdst
        ptsy(3) = mext%byimf
        ptsy(4) = mext%bzimf
        call T96_01(mext%ikp,ptsy,ps,ax,ay,az,abx,aby,abz)
      elseif ( mext%kouter .eq. 8 ) then
        dvc1(1) = ax
        dvc1(2) = ay
        dvc1(3) = az
        call bom97(dvc1,dvc2)
        abx = dvc2(1)
        aby = dvc2(2)
        abz = dvc2(3)
      elseif ( mext%kouter .eq. 9 ) then
        ptsy(1) = mext%pdyn
        ptsy(2) = mext%vdst
        ptsy(3) = mext%byimf
        ptsy(4) = mext%bzimf
        ptsy(5) = mext%g1
        ptsy(6) = mext%g2
        call T01_01(mext%ikp,ptsy,ps,ax,ay,az,abx,aby,abz)
      elseif ( mext%kouter .eq. 10 ) then
        ptsy(1) = mext%pdyn
        ptsy(2) = mext%vdst
        ptsy(3) = mext%byimf
        ptsy(4) = mext%bzimf
        ptsy(5) = mext%w1
        ptsy(6) = mext%w2
        ptsy(7) = mext%w3
        ptsy(8) = mext%w4
        ptsy(9) = mext%w5
        ptsy(10) = mext%w6
        call T04_s(mext%ikp,ptsy,ps,ax,ay,az,abx,aby,abz)
      elseif ( mext%kouter .eq. 11 ) then
        xa(1) = ax
        xa(2) = ay
        xa(3) = az
        ut = (adate%ihour+(adate%imin+adate%secs/60.0D0)/60.0D0) /
     &       24.0D0
        call pstatus(1.0D0,1.0D0,1.0D0,1.0D0,1.0D0,1.0D0,1.0D0)
        call a2000(ut,adate%iyear,adate%imonth,adate%iday,mext%wdens,
     &             mext%wvel,mext%bzimf,mext%vdst,mext%val,xa,bma,bba)
        abx = bma(1) - bba(1,1)
        aby = bma(2) - bba(1,2)
        abz = bma(3) - bba(1,3)
      else
        abx = 0.0d0
        aby = 0.0d0
        abz = 0.0d0
      endif
c
      abx = abx * 1.0D-5
      aby = aby * 1.0D-5
      abz = abz * 1.0D-5
C      
C  c/ CONVERTS SM OR GSM COORDINATES TO GEO
C
      bbx = t(1,1) * abx + t(2,1) * aby + t(3,1) * abz
      bby = t(1,2) * abx + t(2,2) * aby + t(3,2) * abz
      bbz = t(1,3) * abx + t(2,3) * aby + t(3,3) * abz
C      
C  d/ TRANSFORM CARTESIAN VECTOR COMPONENTS 
C     TO SPHERICAL VECTOR COMPONENTS
C
      dbr = sla * clo * bbx + sla * slo * bby + cla * bbz
      dbt = cla * clo * bbx + cla * slo * bby - sla * bbz
      dbp = -slo * bbx + clo * bby
C
C
      END
C----------------------------------------------------------------------
      SUBROUTINE UM533
     :                ( mgeo, rmagp, ifail )
C
C!    Distance to the magnetopause
C
C     REFERENCE
C     RMAGP in SHELL.for; Fairfield JGR 76.28; Carl McIlwain
C
      INCLUDE 'structure.h'
C
C     INTERFACE
C
        TYPE(zgeo) :: mgeo
        REAL*8        rmagp
        INTEGER*4     ifail
C
      COMMON /UC140/  mint, mext, msun 
C
        TYPE(zimf) :: mint
        TYPE(zemf) :: mext
        TYPE(zsun) :: msun
     
      COMMON /UC160/  pi, deg, re, gmagmo, eclipt, geoid, uma
C
        REAL*8        pi, deg, re, gmagmo, eclipt, geoid(3), uma(30)

      COMMON /UC190/ prop, stepx, stpmin, umsq, upsq, uk2, uk3, 
     :                epskm, epsrel, stplst, xclat, kmflg, kum533       

      REAL*8         prop, stepx, stpmin
      REAL*8         umsq, upsq, uk2, uk3
      REAL*8         epskm, epsrel, stplst, xclat
      INTEGER*4      kmflg, kum533      
C
C     VARIABLES
C
        REAL*8        rmp(7)
        REAL*8        gs, sn, aa, bb, dd, delta
C
      DATA rmp/
     :   0.028D0, 0.35D0, -0.58D0, 17.87D0, -233.7D0, 1.057D0, 0.0526D0/
C
C     CODE
C
      ifail = 0
C      
      gs    = - COS( ( msun%utdeg + mgeo%elong ) * deg )
      sn    = SQRT( 1.0D0 - gs**2 )
      dd    = sn**2 + rmp(1) * sn * gs + rmp(2) * gs**2
      aa    = 0.5d0 * ( rmp(3) * sn + rmp(4) * gs ) / dd
      bb    = aa**2 - rmp(5) / dd
      if( bb .ge. 0.0d0 )then
        rmagp = re * ( -aa + SQRT( bb ) ) * 
     :          ( rmp(6) - rmp(7) * mext%vkp 
     :                     / ( 1.0D0 + 0.1D0 * mext%vkp ) )
        delta = rmagp - mgeo%radius 
      else
        ifail = 10
        delta = 0.0d0
        rmagp = 0.0d0
      endif
C      
      if( delta .lt. 0 .and. kum533 .gt. 0)then
        ifail = -53301
        return
      endif
C
C
      END
C----------------------------------------------------------------------
      SUBROUTINE UM535
     :                ( mgeo, mgde )
C
C!    Convert from geocentric to geodetic
C
C     REFERENCE
C     IAU, 1964; ASTRON. J. VOL.66 P.15, 1961
C
      INCLUDE 'structure.h'
cDEC$ IF DEFINED (_x86_)
cDEC$ ATTRIBUTES DLLEXPORT :: UM535
cDEC$ ENDIF
C
C     INTERFACE
C
        TYPE(zgeo) :: mgeo
        TYPE(zgeo) :: mgde
C
      COMMON /UC160/  pi, deg, re, gmagmo, eclipt, geoid, uma
C
        REAL*8        pi, deg, re, gmagmo, eclipt, geoid(3), uma(30)
C
C     VARIABLES
C              
        REAL*8        A2, A4, A6, A8, rer
        REAL*8        scl, ccl, s2cl, c2cl, s4cl 
        REAL*8        c4cl, s8cl, s6cl, dltcl
C
C     CODE
C
      RER = mgeo%radius / geoid(1)
      A2  = ((-1.4127348D-8/RER+0.94339131D-8)/RER+0.33523288D-2) / RER
      A4  = (((-1.2545063D-10/RER+0.11760996D-9)/RER+0.11238084D-4)/RER
     :       -0.2814244D-5) / RER
      A6  = ((54.939685D-9/RER-28.301730D-9)/RER+3.5435979D-9) / RER
      A8  = (((320.0D0/RER-252.0D0)/RER+64.0D0)/RER-5.0D0) / RER *
     :       0.98008304D-12
C
      SCL   = COS( mgeo%colat * deg )
      CCL   = SQRT(1.0D0-SCL*SCL)
      S2CL  = 2.0D0 * SCL * CCL
      C2CL  = 2.0D0 * CCL * CCL - 1.0D0
      S4CL  = 2.0D0 * S2CL * C2CL
      C4CL  = 2.0D0 * C2CL * C2CL - 1.0D0
      S8CL  = 2.0D0 * S4CL * C4CL
      S6CL  = S2CL * C4CL + C2CL * S4CL
      DLTCL = S2CL * A2 + S4CL * A4 + S6CL * A6 + S8CL * A8
C
      mgde%elong  = mgeo%elong
      mgde%colat  = mgeo%colat - dltcl / deg
      mgde%radius = mgeo%radius + re - geoid(1) 
     :                           / SQRT(1.0D0+geoid(3)*SCL*SCL)
C
C
      END
C----------------------------------------------------------------------
      SUBROUTINE UM536
     :                ( mgde, mgeo )
C
C!    Convert from geodetic to geocentric
C
C     REFERENCE
C     IAU, 1964; ASTRON. J. VOL.66 P.15, 1961
C
      INCLUDE 'structure.h'
cDEC$ IF DEFINED (_x86_)
cDEC$ ATTRIBUTES DLLEXPORT :: UM536
cDEC$ ENDIF
C 
C     INTERFACE
C
        TYPE(zgeo) :: mgde
        TYPE(zgeo) :: mgeo
C
      COMMON /UC160/  pi, deg, re, gmagmo, eclipt, geoid, uma
C
        REAL*8        pi, deg, re, gmagmo, eclipt, geoid(3), uma(30)
C
C     VARIABLES
C
        REAL*8        sinlat, coslat, costh, sinth, rgeoid, x, y
C
C
        sinlat      = COS( mgde%colat * deg )
        coslat      = SQRT(1.0D0-sinlat**2)
        costh       = sinlat / SQRT((geoid(2)*coslat)**2+sinlat**2)
        sinth       = SQRT(1.0d0-costh**2)
        rgeoid      = geoid(1) / SQRT(1.0d0+geoid(3)*costh**2)
        x           = rgeoid * sinth + (mgde%radius-re) * coslat
        y           = rgeoid * costh + (mgde%radius-re) * sinlat
C
        mgeo%radius = SQRT(x*x+y*y)
        mgeo%colat  = 90.0d0 -  DATAN(y/x)/deg
        mgeo%elong  = mgde%elong
C
C
      END              
C----------------------------------------------------------------------
      SUBROUTINE UM537
     :          (rkm, st, ct, sph, cph, br, bt, bp, ifail)
C
C!    Kluge evaluation of the geomagnetic field
C                       
C     REFERENCES
C     G. Kluge, 1907, ESOC Internal Note 61, a generalised method for
C     the calculation of the geomagnetic field from multipole
C     expansions
C
      INCLUDE 'structure.h'
C
C     INTERFACE
C
        REAL*8    rkm, st, ct, sph, cph, br, bt, bp
        INTEGER*4 ifail
C
C
      COMMON /UC140/ mint, mext, msun 
C
        TYPE(zemf) :: mext
        TYPE(zimf) :: mint
        TYPE(zsun) :: msun         
C
C
C     VARIABLE
        REAL*8        gnew(nx140,nx140)
        REAL*8        S, T, RSQ
        TYPE(zxyz) :: xi, b
        INTEGER*4     n, kl, km
C
      rsq  = 1 / rkm**2
      t    = 2 * rsq / rkm
      xi%z = ct / rkm
      xi%x = st * cph / rkm
      xi%y = st * sph / rkm     
C
      kl   = mint%norder
      gnew(kl,1)      = mint%coef(kl,1)
      do km = 2, kl
        gnew(kl,km)   = mint%coef(kl,km)
        gnew(km-1,kl) = mint%coef(km-1,kl)
      enddo
C
      do n = 0 , 1
C
        do kl = mint%norder-1, n+1, -1
c km=1
          gnew(kl,1)      = mint%coef(kl,1) + ( 
     :                               xi%x * ( 2*gnew(kl+1,2) )
     :                         +     xi%y * ( 2*gnew(1,kl+1)   )
     :                         + 2 * xi%z *    gnew(kl+1,1)
     :                       ) / (kl-n)
c km=2
          if( kl.gt.1 )then
            gnew(kl,2)      = mint%coef(kl,2) + ( 
     :                        xi%x * ( gnew(kl+1,3)-gnew(kl+1,1) )
     :                  +     xi%y *         gnew(2,kl+1)         
     :                  + 2 * xi%z *         gnew(kl+1,2)
     :                ) / (kl-n)
            gnew(1,kl)      = mint%coef(1,kl) + (
     :                        xi%x *         gnew(2,kl+1)  
     :                  -     xi%y * ( gnew(kl+1,3)+gnew(kl+1,1) )
     :                  + 2 * xi%z *         gnew(1,kl+1)
     :                ) / (kl-n)
            do km = 3, kl
              gnew(kl,km)   = mint%coef(kl,km) + ( 
     :                      xi%x * ( gnew(kl+1,km+1)-gnew(kl+1,km-1) )
     :                +     xi%y * (  gnew(km,kl+1)+gnew(km-2,kl+1)  )
     :                + 2 * xi%z *          gnew(kl+1,km)
     :              ) / (kl-n)
              gnew(km-1,kl) = mint%coef(km-1,kl) + (
     :                      xi%x * (  gnew(km,kl+1)-gnew(km-2,kl+1)  )
     :                -     xi%y * ( gnew(kl+1,km+1)+gnew(kl+1,km-1) )
     :                + 2 * xi%z *          gnew(km-1,kl+1)
     :              ) / (kl-n)
            enddo
          endif
        enddo
C
      enddo 
C
      s    = gnew(1,1)*0.5d0 + 2.0d0 * 
     :        ( gnew(2,1)*xi%z+gnew(2,2)*xi%x+gnew(1,2)*xi%y )
      b%x  = t * ( gnew(2,2) - s * st * cph * rkm )
      b%y  = t * ( gnew(1,2) - s * st * sph * rkm )
      b%z  = t * ( gnew(2,1) - s * ct * rkm )
C
      br   = b%x*st*cph + b%y*st*sph + b%z*ct
      bt   = b%x*ct*cph + b%y*ct*sph - b%z*st
      bp   =  -b%x*sph  + b%y*cph
C
      ifail=0
C
      END                                              
C----------------------------------------------------------------------
      SUBROUTINE UM538 (mpos, amjd, xmlt, xlat, ifail) 
C
C     Evaluation of the magnetic local time and latitude
C     To be included in UNILIB 2.03
C
C     mpos:  Geographic position in geocentric coordinates (GEO) [in]
C     amjd:  Modified Julian day [in]
C     xmlt:  Magnetic Local Time [out]
C     xlat:  Magnetic Latitude [out]
C     ifail: Error flag [out] 
C
C     The subroutine UM538 evaluates the magnetic local time at a given 
C     geographic location and time. The time (argument amjd) passed to 
C     the subroutine is used to re-evaluate the position of the Sun and the
C     GSM or SM coordinate transformation used in magnetic field evaluations 
C     (only the external magnetic field model is affected by this
C     re-evaluation). 
C     Note that the magnetic field models have to be firstly initialized, 
C     e.g. by a call to the subroutine UM510 and UM520. 
C
C
C     Code:
C                       
      INCLUDE 'structure.h'
C
C     INTERFACE
C
        TYPE(zgeo) :: mpos
        REAL*8        amjd, xmlt, xlat
        INTEGER*4     ifail
C
C     COMMON
C
      COMMON /UC140/  mint, mext, msun
      TYPE(zimf) ::   mint
      TYPE(zsun) ::   msun
      TYPE(zemf) ::   mext          
      COMMON /UC160/  pi, deg, re, gmagmo, eclipt, geoid, uma
      REAL*8          pi, deg, re, gmagmo
      REAL*8          eclipt, geoid(3), uma(30)               
C
C     INTERNAL VARIABLES
C
      INTEGER*4       kunit, i, j
      TYPE(zxyz) ::   mdps, mdsm
      REAL*8          st, rot(3,3)
C
      DATA            kunit/-6/
C
C
C
      ifail    = 0
C
      mdps%z   = cos(mpos%colat*deg)
      mdps%x   = sin(mpos%colat*deg)*cos(mpos%elong*deg)
      mdps%y   = sin(mpos%colat*deg)*sin(mpos%elong*deg)
C
      CALL UM522 (amjd,kunit)
      CALL UM524 (kunit)
C                            
      do i=1,3
        do j=1,3
          rot(i,j)=mext%trans(i,j)
        enddo
      enddo
C
C-MK-2005july13: Replace the following lines
C
C#### if ( mext.kouter .ge. 2 .and. mext.kouter .le. 4 ) 
C####:   CALL UM523 (kunit)
C
C-MK-2005july13: by the lines
C                                                                 
      if (mext%kouter .eq.2 .or. mext%kouter .eq.3 
     :    .or. mext%kouter .eq.4 .or. mext%kouter .eq.7) then
        CALL UM523 (kunit)
C       Note-> for the other kouter, UM524 has already been called
      endif
C
C-MK-2005july13: end-of-correction 
C
      mdsm%x   = mdps%x*rot(1,1) + mdps%y*rot(1,2) + mdps%z*rot(1,3)
      mdsm%y   = mdps%x*rot(2,1) + mdps%y*rot(2,2) + mdps%z*rot(2,3)
      mdsm%z   = mdps%x*rot(3,1) + mdps%y*rot(3,2) + mdps%z*rot(3,3)
C
      st       = sqrt(mdsm%x**2+mdsm%y**2)
C
      if( st .lt. 1.d-10 )then
C       Problem:
C         Location alignes with dipole axis!
C
        ifail  = -53802
        xmlt   = -999.0D0
        xlat   = 90*mdsm%z
        return
      endif
C
      xlat     = asin(mdsm%z)/deg
      xmlt     = atan2(mdsm%y,mdsm%x)/deg/15.0D0+12.0D0
C                            
      end
C----------------------------------------------------------------------
      SUBROUTINE UM539
     :                (mpos, amjd, mb, ifail)
C
C!    Evaluate the magnetic field
C                       
      INCLUDE 'structure.h'
cDEC$ IF DEFINED (_x86_)
cDEC$ ATTRIBUTES DLLEXPORT :: UM539
cDEC$ ENDIF
C
C      
C     INTERFACE
C
        TYPE(zgeo) :: mpos
        TYPE(zvec) :: mb
        REAL*8        amjd
        INTEGER*4     ifail
*
* see UM530 except:
*     amjd = Modified Julian day
*
c The subroutine UM539 evaluates the magnetic field vector 
c at a given time and geographic location. The subroutine
c differs from subroutine UM530 by its argument amjd from
c which the time is specified. Note that the magnetic field
c models have to be firstly initialized, e.g. by a call to
c the subroutine UM510 and UM520.
c <p>
c The time (argument amjd) passed to the subroutine is used
c to re-evaluate the position of the Sun and the GSM or SM
c coordinate transformation. Note that only the external magnetic
c field model is affected by the re-evaluation.
C
      COMMON /UC140/ mint, mext, msun

      TYPE(zimf) :: mint
      TYPE(zsun) :: msun
      TYPE(zemf) :: mext          
C
      INTEGER*4 kunit
      DATA kunit/-6/
C                                        
      CALL UM522 (amjd,kunit)
C
C-MK-2006january18: Replace the following lines
C
C#### if ( mext%kouter .ge. 2 .and. mext%kouter .le. 4 ) 
C####:   CALL UM523 (kunit)
C
C-MK-2006january18: by the lines
C                                                                 
      if (mext%kouter .eq.2 .or. mext%kouter .eq.3 
     :    .or. mext%kouter .eq.4 .or. mext%kouter .eq.7) then
        CALL UM523 (kunit)
C
C-MK-2006january18: end-of-correction 
C
      else
        CALL UM524 (kunit)
      endif
      CALL UM530 (mpos, mb, ifail)
C
      END
C----------------------------------------------------------------------
