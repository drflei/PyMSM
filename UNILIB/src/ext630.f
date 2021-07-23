# 1 "ext630.f"
C FILE EXT630.F
      subroutine ext630
C          dummy subroutine
      end
C
C modif:   change gtd6    --> nwgtd6
C                 integer --> integer*4
C                 real    --> real*4
C                 logical --> logical*4
C          add    implicit none
C          remove write(6,*)
C
C........................................................................
C MSIS: 
C
C------------------------------------------------------------------------=
C
C    FILE ( MODULE )   :   asem49.f
C
C------------------------------------------------------------------------=
C     LAST CODE CHANGE  on  JAN 9 1996  ON acd.ucar.edu 
C------------------------------------------------------------------------=
C    NEW PROGRAM FOR GTD6 (SEMI-EMPIRIC ATMOSPHERIC MODEL FROM A.E. HEDIN)
C    BASED ON MODIFIED VERSION msis27.f
C    BY SIMON CHABRILLAT - I.A.S. - OCT 1994 - simonc@oma.be 
C
C------------------------------------------------------------------------=
C------------------------------------------------------------------------=
      SUBROUTINE NWGTD6(IYD,ALT,GLAT,GLONG,STL,F107A,F107,AP,MASS,D,T,
     &                  IW)    
      implicit none
C        Neutral Atmosphere Empirical Model from the surface to lower
C          exosphere  MSISE93 (JGR, 96, 1159-1172, 1991). Formulas (A.xx)=
C          refer to math appendix of JGR, 88, 10170-10188, 1983.
C         A.E.Hedin 4/24/90;6/3/91(add SAVE)
C           See subroutine GHP6 to specify a pressure rather than
C           altitude.
C
C     INPUT : 
C        IYD - DAY AS DDD (day of year from 1 to 365)
C        ALT - ALTITUDE(KM)
C        GLAT - GEODETIC LATITUDE(DEG)
C        GLONG - GEODETIC LONGITUDE(DEG)
C        F107A - 3 MONTH AVERAGE OF F10.7 FLUX
C        F107 - DAILY F10.7 FLUX FOR PREVIOUS DAY
C
C        AP - MAGNETIC INDEX(DAILY) OR WHEN    IW(9) = -1  :
C           - ARRAY CONTAINING:
C             (1) DAILY AP
C             (2) 3 HR AP INDEX FOR CURRENT TIME
C             (3) 3 HR AP INDEX FOR 3 HRS BEFORE CURRENT TIME
C             (4) 3 HR AP INDEX FOR 6 HRS BEFORE CURRENT TIME
C             (5) 3 HR AP INDEX FOR 9 HRS BEFORE CURRENT TIME
C             (6) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 12 TO 33 HRS
C                    PRIOR TO CURRENT TIME
C             (7) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 36 TO 59 HRS
C                    PRIOR TO CURRENT TIME
C
C        MASS - MASS NUMBER : ONLY DENSITY FOR SELECTED GAS IS
C               CALCULATED. SEE OUTPUT DENSITIES D(1..8).
C               MASS 0 FOR TEMPERATURE ONLY. =

C               MASS 48 TO GET ALL NUMBER DENSITIES AND TOTAL MASS DENSITY D(8)
C
C        IW - SWITCHES TO TURN ON AND OFF PARTICULAR VARIATIONS TERMS
C             IW IS A 26 ELEMENT ARRAY. CONCERNED TERMS & RECOMMENDED VALUES :
C
C            IW(1)  - F10.7 EFFECT ON MEAN                            :1
C            IW(2)  - TIME INDEPENDENT                                :1
C            IW(3)  - SYMMETRICAL ANNUAL                              :1
C            IW(4)  - SYMMETRICAL SEMIANNUAL                          :1
C            IW(5)  - ASYMMETRICAL ANNUAL                             :1
C            IW(6)  - ASYMMETRICAL SEMIANNUAL                         :1
C            IW(7)  - DIURNAL                                         :1
C            IW(8)  - SEMIDIURNAL                                     :1
C            IW(9)  - DAILY OR 3hrs AP (see below)                    :1
C            IW(10) - ALL UT/LONG EFFECTS                             :1
C            IW(11) - LONGITUDINAL                                    :1
C            IW(12) - UT AND MIXED UT/LONG                            :1
C            IW(13) - MIXED AP/UT/LONG                                :1
C            IW(14) - TERDIURNAL                                      :1
C            IW(15) - DIFFUSIVE,MIXED OR NET DENSITIES (see below)    :1
C            IW(16) - ALL TINF VAR                                    :1
C            IW(17) - ALL TLB VAR                                     :1
C            IW(18) - ALL TN1 VAR                                     :1
C            IW(19) - ALL S VAR                                       :1
C            IW(20) - ALL TN2 VAR                                     :1
C            IW(21) - ALL DdLB VAR                                    :1
C            IW(22) - ALL TN3 VAR                                     :1
C            IW(23) - TURBO SCALE HEIGHT VAR                          :1
C            IW(24) - LATITUDE AND EARTH RADIUS VARIATION OF GRAVITY  :1
C            IW(25) - for security : set to 0  to set all PRECEDING
C                     SWitches to 1  anyway. Set to 1  to transmit the
C                     SWitches you choose.         Recommended value  :0
C            IW(26) - 0  to jump calculations already done, for which
C                     results do not change. Set to 0  to spare
C                     computing time. Set to 1  if not sure of the
C                     SAVE statement efficiency.   Recommended value  :0
C
C             POSSIBLE VALUES FOR IW(1) TO IW(24), IW(9) AND IW(15) EXCEPTED :
C               IW(i)=0  : Main terms & Cross terms OFF
C               IW(i)=1  : Main terms & Cross terms ON
C               IW(i)=2  : Main terms OFF & Cross terms ON
C               IW(i)=3  : Main terms ON & Cross terms OFF
C
C             SET IW(9) = -1  TO GIVE 3 HRS AP : AP(1..7)
C                          0  TO NEGLECT AP EFFECT
C                          1  TO GIVE DAILY AP : AP(1) ONLY
C
C             SET IW(15) = -1  TO GET MIXING DENSITIES ( 0 ABOVE ALTL(i) )
C                               INCLUDING CORRECTION FACTORS CCOR 
C                           0  TO GET DIFFUSIVE DENSITIES ( 0 BELOW 72.5km)
C                           1  TO GET NET (REAL*4) DENSITIES
C           ---------------------------------------
C     NOTE :  F107, F107A, and AP effects are not large below 80 km 
C            and these can be set to 150., 150., and 4. respectively.
C
C     OUTPUT:   
C        D(1) - H (1)  NUMBER DENSITY(CM-3)   
C        D(2) - N (14) NUMBER DENSITY(CM-3)
C        D(3) - O (16) NUMBER DENSITY(CM-3)
C        D(4) - He (4) NUMBER DENSITY(CM-3)   
C        D(5) - O2(32) NUMBER DENSITY(CM-3)   
C        D(6) - Ar(40) NUMBER DENSITY(CM-3)
C        D(7) - N2(28) NUMBER DENSITY(CM-3)   
C        D(8) - TOTAL MASS DENSITY(GM/CM3)
C        T(1) - EXOSPHERIC TEMPERATURE
C        T(2) - TEMPERATURE AT ALT
C
C
C      O, H, and N set to zero below 72.5 km
C      
C----------------------------------------------------------------------
C       VERSION 0.0 : CALCULATES T(ALT) IF LEVEL=0
C               0.1 : CALCULATES DdLB(i), DIFFUSIVE DENSITY AT 120km 
C               0.2 : CALCULATES Ddiff(ALT) IF LEVEL=0, FCT Bates(Diff)DEns
C               1.0 : CALCULATES T(ALT)
C               1.1 : CALCULATES DIFFUSIVE DENSITIES IN LEVELS 0 & 1
C               1.2 : CALCULATES TOTAL DENSITIES IN LEVEL 0
C               1.3 : CALCULATES TOTAL DENSITIES IN LEVEL 0 & 1
C               2.0 : THESE DENSITIES COMPUTED IN BIG ROUTINE CalcDensUp
C               2.1 : CALCULATES TOTAL DENSITIES DOWNTO ZMIX = 62.5 km
C               2.2 : CALCULATES TOTAL DENSITIES IN ALL LEVELS 
C               3.0 : NEW SW(26), IF 1. -> COMPUTES ALL ( EVEN IN GLOBE6 )
C                     AT EACH GTD6 CALL ( USE IF NOT SURE OF *SAVE* ) 
C               3.1 : JUMPS CALCULATIONS DONE AT PRECEDING GTD6 CALL
C               3.2 : SW(15)=-1 -> MIX DENS ( 0 ABOVE ALTL(i) ) 
C                             0 -> DIFF DENS ( 0 BELOW 72.5 km)
C                             1 -> NET DENS
C               3.3 : ADDED GHP6 ROUTINE ( AT END )
C               3.4 : UNDERSTOOD DmixZH WAS TO RENAME DmLB, CHANGE ORDER
C                     AND NAMES IN CalcDensUp TO EXPLAIN Dmix CALCULATION
C               3.5 : CHANGED Dmix CALCULATION : in CalcDensUp, CALC DmixZA
C                     DIRECTLY INSTEAD OF DmLB
C               4.0 : REMOVED DmixZA FROM GTD6 ( UNUSED )
C               4.1 : CHANGED IYD, STL and DAY COMPUTING
C               4.2 : BUG WHEN MASS<>48 OR 28 AND SW(26)=1 FIXED
C               4.4 : line 50 of GLOBE6 re-modified 
C                       TINF = TINF+ABS(SWM(i))*T(i) following D. Heynderickx
C                     subroutine setswitch modified for SW(25)<>1
C                     Values of SW(9) and SWM(9) corrected
C               4.5 : adapted for PC (Microsoft fortran PowerStation)
C                     line 1123 : rout SPLINE, dimension U(NMAX=5)
C               4.6 : BUG WHEN MASS=0, SW(26)=0 AND ALT<122.807 FIXED
C                     Test on jumping complete calculation of TN1(2:5) completed
C
C  Modifications by Claude Camy-Peyret :
C
C               4.7 : bug corrected :
C                       calculation of DMC made AFTER line labelled 200
C                     SW  (REAL*4) changed to IW  (INTEGER*4)
C                     SWL (REAL*4) changed to IWL (INTEGER*4)
C                     SW9,SWM,SWC (REAL*4) changed to INTEGER*4 in COMMON /CSW/
C                     and 0./1. changed to 0/1 when necessary to be consistent
C               4.8 : EQUIVALENCE for more rigorous calls of GLOBE6 and GLOB6S
C               
C  by S. Chabrillat again :
C
C               4.9 : Bug between 62.5-72.5 km finally corrected :
C                       test on goto 200 more strict.
C------------------------------------------------------------------------=
      INTEGER*4 IYD,MASS,IW(26)
      REAL*4 ALT,GLAT,GLONG,STL,F107A,F107,AP(7),D(8),T(2)
C
      INTEGER*4 MN1,MN2,MN3
      parameter (MN1=5,MN2=4,MN3=5)
      REAL*4 ZN2(MN2),ZN3(MN3)
      REAL*4 TN1(MN1),TN2(MN2),TN3(MN3),TGN1(2),TGN2(2),TGN3(2)
      REAL*4 XS2(MN2),YS2(MN2),Y2OUT2(MN2)
      REAL*4 XS3(MN3),YS3(MN3),Y2OUT3(MN3)
      REAL*4 DdLB(7),Ddiff(7),Dmix(7),D12(8)
      character ver*40
      LOGICAL*4 jump
      INTEGER*4 Level
C
      INTEGER*4 SW9,SWM(26),SWC(24)
      COMMON /CSW/ SW9,SWM,SWC
      REAL*4 DAY,fixLAT
      common /FixInput/  DAY,fixLAT
      REAL*4 GSURF,RE,RGAS,DGTR,DR
      common /PhysCte/   GSURF,RE,RGAS,DGTR,DR
      REAL*4 XMM,ALTL(8),alfa(7)
      INTEGER*4 MassNo(7),NoPDM(7),NoPD(7)
      common /FixData/   XMM,MassNo,ALTL,alfa,NoPDM,NoPD
      REAL*4 ZN1(5),ZLB,ZA,S,ZGA,ZGDIF1
      common /FixAlt/    ZN1,ZLB,ZA,S,ZGA,ZGDIF1
      REAL*4 TA,TLB,TINF,DTA
      common /FixTemp/   TA,TLB,TINF,DTA
      REAL*4 XS1(5),YS1(5),Y2OUT1(5),YI1
      common /Splin1/    XS1,YS1,Y2OUT1,YI1
      REAL*4 PTM(10),PDM(10,8)
      COMMON /LOWER6/    PTM,PDM
      REAL*4 PT(150),PD(150,9),PS(150),PDL(25,2),PTL(100,4),PMA(100,10)
      COMMON /PARM6/     PT,PD,PS,PDL,PTL,PMA
      REAL*4 PAVGM(10)
      COMMON /MAVG6/     PAVGM
c
      REAL*4 pma1(100),pma2(100),pma3(100),pma4(100),pma5(100)
      REAL*4 pma6(100),pma7(100),pma8(100),pma9(100),pma10(100)
      EQUIVALENCE (PMA(1,1),pma1(1))
      EQUIVALENCE (PMA(1,2),pma2(1))
      EQUIVALENCE (PMA(1,3),pma3(1))
      EQUIVALENCE (PMA(1,4),pma4(1))
      EQUIVALENCE (PMA(1,5),pma5(1))
      EQUIVALENCE (PMA(1,6),pma6(1))
      EQUIVALENCE (PMA(1,7),pma7(1))
      EQUIVALENCE (PMA(1,8),pma8(1))
      EQUIVALENCE (PMA(1,9),pma9(1))
      EQUIVALENCE (PMA(1,10),pma10(1))
c
      REAL*4 ptl1(100),ptl2(100),ptl3(100),ptl4(100)
      EQUIVALENCE (ptl(1,1),ptl1)
      EQUIVALENCE (ptl(1,2),ptl2)
      EQUIVALENCE (ptl(1,3),ptl3)
      EQUIVALENCE (ptl(1,4),ptl4)
c
      REAL*4 pd4(150),pdi(150)
      EQUIVALENCE (pd(1,4),pd4)
C
      SAVE
      EXTERNAL GTD6BK
C
      REAL*4 GLOBE6,BatesDens,GLOB6S,SPLINI,SplinDens,TotMassDens
      REAL*4 ZMIX,ZETA,ZZ,ZL,ZG0,ALAST,ALT1,ALT2,SoTL,SEC,XLAT,G,XM
      REAL*4 X1,X2,X3,ZGDIF2,ZGDIF3,YI2,YI3,DMC,DMR,DMU
      REAL*4 TGLB,Tup,Tmid,D23N2
      INTEGER*4 i,j,LastLev,IYDa,IDAYinc
      REAL*4 v0, v24, v365, v3600, v86400
C
      DATA ZN1/122.807,110.,100.,90.,72.5/
      DATA ZN2/72.5,55.,45.,32.5/
      DATA ZN3/32.5,20.,15.,10.,0./
      DATA ZMIX/62.5/,LastLev/-9/,ALAST/-999./
      DATA ALTL/320.,450.,400.,200.,200.,240.,160.,450./
      data MassNo/1,14,16,4,32,40,28/, alfa/-.4,0.,0.,-.4,0.,0.,0./
      data NoPD/7,8,2,1,5,6,3/, NoPDM/6,7,2,1,4,5,3/
      DATA RGAS/831.4/,DGTR/1.74533E-2/,DR/1.72142E-2/
      data v0, v24, v365, v3600, v86400/0., 24., 365., 3600., 86400./
C                        Internal fct
      ZETA(ZZ,ZL)=(ZZ-ZL)*(RE+ZL)/(RE+ZZ)
C
      if (IW(9).ne.-1) then
        do j=2,7
          AP(j) = AP(1)
        enddo
      endif
C       TESTS IF NO OTHER INPUT CHANGE THAN ALT HAPPENED & IF IW(26).NE.1
C       IF jump is .TRUE., JUMPS CALCULATIONS ALREADY DONE TO SPARE TIME
      call TestInput (IYD,GLAT,GLONG,STL,F107A,F107,AP,MASS,IW,jump)
C
C       NEW TIME ADJUSTMENTS
      if (jump) goto 5
      IDAYinc = STL/v24
      IYDa = IYD + IDAYinc
      SoTL = MOD(STL,v24)
      DAY = IYDa
      DAY = MOD(DAY,v365) + SoTL/v24
      fixLAT = GLAT

      SEC = ( SoTL- GLONG/15.) * v3600
      if (SEC.lt.v0) SEC = SEC + v86400 
      if (SEC.gt.v86400) SEC = SEC - v86400
C
C      if (LastLev.eq.-9.and.IW(25).ne.1) write(6,*)
C     +  ' IW(25)<>1 : ALL VARIATIONS TERMS COMPUTED '
C      if (LastLev.eq.-9.and.IW(26).eq.1) then
C      write(6,*)
C     +  ' IW(26) = 1 : COMPUTING ALL PARAMETERS FOR EVERY POINT'
C      write(6,*) ' LOSS OF COMPUTING TIME'
C      endif
      call setswitch(IW)
c      if (SWM(15).eq.0.and.ALT.lt.ZN1(5))
c     +     write(6,*) 'Warning : IW(15)=0 - i.e. you want diffusive ',
c     +                'densities - & ALT = ',ALT,' km < 72.5 km !'    
C
C
C         Latitude variation of gravity (none for SWM(24)=0)
      XLAT=GLAT
      IF(SWM(24).EQ.0) XLAT=45.
      CALL GLATF(XLAT,GSURF,RE)
C
    5 continue
C       ZLB=PTM(6)=120km ; ZN1(1)=ZA=PDL(16,2)=PK1(41)=122.807km

      ZLB = PTM(6)
C     ZN1(1) = PDL(16,2)
      ZA = ZN1(1)
C       XMM=PDM(5,3)=28.95
      XMM = PDM(5,3)
      Level = 0
      if (ALT.lt.ZN1(1)) Level=1
      if (ALT.lt.ZN2(1)) Level=2
      if (ALT.lt.ZN3(1)) Level=3
      if (ALT.lt.0.) write(6,*) 'NOT ABLE TO COMPUTE ATMOSPHERIC ',
     + 'CONDITIONS UNDER THE LEVEL OF THE SEA !'
c     write(7,*) ' Level = ',Level,' LastLev = ',LastLev,' jump = ',
c    + jump,' ALAST = ', ALAST
C        

C-.......................................................................
C       LEVEL 0 : BATES PROFILE TEMPERATURE
C-.......................................................................
      do J=1,8
        D(J) = 0.
      enddo
C..(A.7),(A.8),(A.6),(A.5) to calculate S=SIGMA........................
      if (jump) goto 10
      TINF = PTM(1) * PT(1) *
     $  (1.+SWM(16)*GLOBE6(DAY,SEC,GLAT,GLONG,SoTL,F107A,F107,AP,PT))
      T(1) = TINF
C         TLB : Temperature at 120 km
      TLB = PTM(2)*(1.+SWM(17)*GLOBE6(DAY,SEC,GLAT,GLONG,SoTL,
     $                                F107A,F107,AP,PD4))*PD(1,4)
C         TGLB : Temperature gradient at 120 km      
      TGLB = PTM(4) * PS(1) *
     $   (1.+SWM(19)*GLOBE6(DAY,SEC,GLAT,GLONG,SoTL,F107A,F107,AP,PS))
      S = TGLB/(TINF-TLB)
C......Temperature & gradient of T at ZA from Bates profile..............=
      ZGA = ZETA(ZA,ZLB)
      TA = TINF - (TINF-TLB)*EXP(-S*ZGA)
      DTA = (TINF-TA) * S * ((RE+ZLB)/(RE+ZA))**2
   10 continue
C..(A.1a),(A.4a) : Bates temperature.....................................=
      if (Level.eq.0) then
        ZG0 = ZETA(ALT,ZLB)
        if( s .lt. 0.0 )then
          print*,' *** NWGTD6 ***'
          print*,'   IYD = ',iyd
          print*,'   ALT = ',alt
          print*,'  GLAT = ',glat
          print*,' GLONG = ',glong
          print*,'   STL = ',stl
          print*,' F107A = ',f107a
          print*,'  F107 = ',f107
          print*,'    AP = ',(ap(i),i=1,4)
          print*,'    AP = ',(ap(i),i=5,7)
          print*,'    IW = ',(iw(i),i=1,5)
          print*,'    IW = ',(iw(i),i=6,10)
          print*,'    IW = ',(iw(i),i=11,15)
          print*,'    IW = ',(iw(i),i=16,20)
          print*,'    IW = ',(iw(i),i=21,25)
          print*,'    IW = ',(iw(i),i=26,26)
          print*,' **************'
          print*,'  TINF = ',tinf
          print*,'   TLB = ',tlb
          print*,'  TGLB = ',tglb
          print*,'     S = ',s
          print*,' **************'
          stop
        endif
        T(2) = TINF - (TINF-TLB)*EXP(-S*ZG0)
        if (MASS.eq.0) goto 999
      endif
C..(A.18) : Calculates DdLB(i), Diffusive Density at ZLB=120km for each specy..
      if ( .not.jump .or. (ALAST.gt.ALTL(8).and.ALT.le.ALTL(8)) ) then
        do i=1,7
          if (MASS.eq.48 .or. MASS.eq.MassNo(i)
     +        .or. (MassNo(i).eq.28.and.ALT.le.ALTL(8)) ) then
            do J=1,150
             pdi(J)=PD(J,NoPD(i))
            enddo
            G = SWM(21) * GLOBE6(DAY,SEC,GLAT,GLONG,SoTL,F107A,F107,AP,
     +                           pdi)
            DdLB(i) = PDM(1,NoPDM(i)) * EXP(G) * PD(1,NoPD(i))
          endif
        enddo
      endif
C........................................................................=
      if (Level.eq.0) then
C            (A.13) : Calculates Ddiff(i), Diffusive Density at ALT
        do i=1,7
          if ( MASS.eq.48 .or. MASS.eq.MassNo(i)
     +        .or. (MassNo(i).eq.28.and.ALT.le.ALTL(8)) ) then
            XM =  ( MassNo(i) )
            Ddiff(i) = BatesDens (ZG0,XM,T(2),DdLB(i),alfa(i))
          endif
        enddo
      endif
      if ( ALT.ge.ALTL(8) .or. (SWM(15).eq.0.and.Level.eq.0) ) then
        do i=1,7
          D(i) = Ddiff(i)
          if (SWM(15).eq.-1) D(i) = v0
        enddo
        goto 998
      endif
C
C-.......................................................................=
C       COMPUTING T & grad T AT NODES
Ce       Lower thermosphere temp variations not significant for
Ce        density above ALTL(8) = 450 km
C-.......................................................................=
      TN1(1) = TA
      TGN1(1) = DTA
      IF (ALT.LT.ALTL(8)) THEN
        if ( .not.jump .or. ALAST.ge.ALTL(8)
     $                .or. (Lastlev.lt.1 .and. Level.ge.1)) then
          TN1(2) = PTM(7) * PTL(1,1) / (1.-SWM(18)*GLOB6S(PTL1))
          TN1(3) = PTM(3) * PTL(1,2) / (1.-SWM(18)*GLOB6S(PTL2))
          TN1(4) = PTM(8) * PTL(1,3) / (1.-SWM(18)*GLOB6S(PTL3))
          TN1(5) = PTM(5) * PTL(1,4)
     +               / (1.-SWM(18)*SWM(20)*GLOB6S(PTL4))
          TGN1(2) = (1.+SWM(18)*SWC(20)*GLOB6S(PMA9))
     $               * PTM(9) * PMA(1,9) * TN1(5)*TN1(5)
     $               / (PTM(5)*PTL(1,4))**2
        endif
       ELSE
         TN1(2) = PTM(7) * PTL(1,1)
         TN1(3) = PTM(3) * PTL(1,2)
         TN1(4) = PTM(8) * PTL(1,3)
         TN1(5) = PTM(5) * PTL(1,4)
         TGN1(2) = PTM(9)*PMA(1,9)*TN1(5)*TN1(5)/(PTM(5)*PTL(1,4))**2
      ENDIF
      if (Level.ge.2) then
        if ( .not.jump .or. LastLev.lt.2 ) then
          TN2(1) = TN1(5)
          TGN2(1) = TGN1(2)
          TN2(2) = PMA(1,1)*PAVGM(1) / (1.-SWM(20)*GLOB6S(PMA1))
          TN2(3) = PMA(1,2)*PAVGM(2) / (1.-SWM(20)*GLOB6S(PMA2))
          TN2(4) = PMA(1,3)*PAVGM(3)
     +             / (1.-SWM(20)*SWC(22)*GLOB6S(PMA3))
          TGN2(2) = (1.+SWM(20)*SWM(22)*GLOB6S(PMA10))
     $   * PAVGM(9)*PMA(1,10)* TN2(4)*TN2(4) / (PMA(1,3)*PAVGM(3))**2
        endif
      endif
      if (Level.eq.3) then
        if ( .not.jump .or. LastLev.lt.3 ) then
          TN3(1) = TN2(4)
          TGN3(1) = TGN2(2)
          TN3(2) = PMA(1,4) * PAVGM(4) / (1.-SWM(22)*GLOB6S(PMA4))
          TN3(3) = PMA(1,5) * PAVGM(5) / (1.-SWM(22)*GLOB6S(PMA5))
          TN3(4) = PMA(1,6) * PAVGM(6) / (1.-SWM(22)*GLOB6S(PMA6))
          TN3(5) = PMA(1,7) * PAVGM(7) / (1.-SWM(22)*GLOB6S(PMA7))
          TGN3(2) = PMA(1,8) *PAVGM(8) * (1.+SWM(22)*GLOB6S(PMA8))
     $               * TN3(5)*TN3(5) / (PMA(1,7)*PAVGM(7))**2
        endif
      endif
C-.......................................................................=
C        COMPUTING T(z) WITH SPLINES, integrating inverse T profile if needed
C-.......................................................................=
C  To compute densities, even if Level>1, YI1 needed at ZN1(5)=ZN2(1)=72.5km !
C                        even if Level=3, YI2 needed at ZN2(4)=ZN3(1)=32.5km !
C-.......................................................................=
      if ( (Level.eq.1.and.MASS.eq.0) .or. (MASS.ne.0) ) then
        if ( jump .and. Level.gt.1 .and. LastLev.gt.1 ) goto 100
        ALT1 = MAX(ALT,ZN1(5))
        call CalcTemp (ALT1,MN1,ZN1,TN1,TGN1,Tup
     +                                   ,XS1,YS1,Y2OUT1,X1,ZGDIF1)
        if (Level.eq.1) T(2) = Tup
        if (MASS.ne.0) YI1 = SPLINI (XS1,YS1,Y2OUT1,MN1,X1)
  100   continue
      endif
      if ( (Level.eq.2.and.MASS.eq.0) .or.
     +                        (Level.ge.2.and.MASS.ne.0) ) then
        if ( jump .and. Level.gt.2 .and. LastLev.gt.2 ) goto 110
        ALT2 = MAX(ALT,ZN2(4))
        call CalcTemp (ALT2,MN2,ZN2,TN2,TGN2,Tmid
     +                                    ,XS2,YS2,Y2OUT2,X2,ZGDIF2)
        if (Level.eq.2) T(2) = Tmid
        if (MASS.ne.0) YI2 = SPLINI (XS2,YS2,Y2OUT2,MN2,X2)
  110   continue
      endif
      if (Level.eq.3) then
        call CalcTemp (ALT,MN3,ZN3,TN3,TGN3,T(2)
     +                                    ,XS3,YS3,Y2OUT3,X3,ZGDIF3)
        if (MASS.ne.0) YI3 = SPLINI (XS3,YS3,Y2OUT3,MN3,X3)
      endif
      if (MASS.eq.0) goto 999
C-.......................................................................=
C        COMPUTING DENSITIES ( WITH CalcDensUp ABOVE OR AT ZN1(5) = 72.5km )
C-.......................................................................=
C............................  LEVELS 0 & 1  ............................=
      if (Level.le.1) then
        call CalcDensUp (Level,ALT,MASS,T(2),DdLB,Ddiff,Dmix,D)
        goto 998
      endif
C                   H, N & O set to zero below 72.5 km'
      do i = 1,3
        D(i) = 0.0
      enddo
      if (SWM(15).eq.0) then
        do i = 4,8
          D(i) = 0.0
        enddo
        goto 999
      endif
C......................  LEVEL 2 ABOVE ZMIX = 62.5km  .................=
      if ( Level.eq.2 .and. ALT.ge.ZMIX ) then
C          Linear transition to full mixing at ZMIX from almost
C          full mixing at ZN2(1)= ZN1(5) to improve efficiency
        if ( jump .and. LastLev.eq.2 .and. ALAST.ge.ZMIX) goto 200
        call CalcDensUp (1,ZN1(5),MASS,TN1(5),DdLB,Ddiff,Dmix,D12)
C          Computing N2 density
        DMR = D12(7)/Dmix(7) - 1.
  200   continue
        DMC = 1. - (ZN2(1)-ALT) / (ZN2(1)-ZMIX)
        D(7) =SplinDens (ZN2(1),ZGDIF2,XMM,TN2(1),T(2),Dmix(7),v0,YI2)
        if (SWM(15).eq.1) D(7) = D(7) * ( 1. + DMC*DMR )
        if (MASS.eq.28) goto 999
        if (SWM(15).eq.-1) goto 250
C          For He, O2, AR : densities = N2 dens * mix ratio * linear transition
        do i = 4,6
          if (MASS.eq.48.or. MASS.eq.MassNo(i)) then      
            DMU = ( D12(i) / (D12(7)*PDM(2,NoPDM(i))) )  - 1.
            D(i) = D(7) * PDM(2,NoPDM(i)) * (1.+DMU*DMC)
          endif
        enddo
        goto 998
      endif
C............N2 DENSITY AT LEVEL 2  BELOW  ZMIX = 62.5km  .............=
      if ( Level.eq.2 .and. ALT.lt.ZMIX ) then
        if ( jump .and. LastLev.eq.2 ) goto 210 
        call CalcDensUp (1,ZN1(5),28,TN1(5),DdLB,Ddiff,Dmix,D12)
  210   continue
        D(7) =SplinDens (ZN2(1),ZGDIF2,XMM,TN2(1),T(2),Dmix(7),v0,YI2)
        if (MASS.eq.28) goto 999
      endif
C...................... N2 DENSITY AT  LEVEL 3 ..........................=
      if ( Level.eq.3) then
        if ( jump .and. LastLev.eq.3 ) goto 220
        call CalcDensUp (1,ZN1(5),28,TN1(5),DdLB,Ddiff,Dmix,D12)
        D23N2 = SplinDens (ZN2(1),ZGDIF2,XMM,TN2(1),TN2(4),Dmix(7),
     +                                                       v0,YI2)
  220   continue
        D(7) = SplinDens (ZN3(1),ZGDIF3,XMM,TN3(1),T(2),D23N2,v0,YI3)
        if (MASS.eq.28) goto 999
      endif
C ..........For He, O2, AR : densities = N2 dens * mix ratio ..........=
  250 continue  
      do i = 4,6
        if (MASS.eq.48 .or. MASS.eq.MassNo(i)) 
     +                  D(i) = D(7) * PDM(2,NoPDM(i))
      enddo    

C-.......................................................................=
  998 if (MASS.eq.48) D(8) = TotMassDens(D,MassNo)
C-.......................................................................=
  999 LastLev = Level
      ALAST = ALT
c
c     Ordres d'ecriture utilises pour mettre en evidence ce qui n'a pas ete
c     bien initialise quand SW(26)=0 et que les altitudes croissent et passent
c     au dessus de ZMIX=62.5 km
c
c     write(7,*) ' DdLB '
c     write(7,*)   DdLB
c     write(7,*) ' Ddiff '
c     write(7,*)   Ddiff
c     write(7,*) ' Dmix '
c     write(7,*)   Dmix
c     write(7,*) ' D12 '
c     write(7,*)   D12
c     write(7,*) ' DMR = ',DMR,' DMC= ',DMC
      return
      entry namever(ver)
           ver = ' asem49.f ,january 1996 on acd.ucar.edu'
      end
C.................................................................
C------------------------------------------------------------------------
C------------------------------------------------------------------------
C
      subroutine CalcDensUp (Level,ALT,MASS,TZ,DdLB,Ddiff,Dmix,D)
      implicit none
C
      INTEGER*4 Level,MASS
      REAL*4 ALT,TZ,DdLB(7),Ddiff(7),Dmix(7),D(8)
C
      INTEGER*4 MN1
      parameter (MN1=5)
      REAL*4 DmixZA(7)
C
C   THIS ROUTINE - EQUIVALENT OF GTS6 DENSITY PART - CALCULATES DENSITIES IN
C   LEVEL 0 ( ALT > ZN1(1)=ZA=122.807km ) 
C   OR LEVEL 1 ( ZN1(5)=ZN2(1)=72.5km < ALT < 122.807km ) .
C   IT IS CALLED WITH ALT=72.5 BELOW 72.5 km 
C
      INTEGER*4 SW9,SWM(26),SWC(24)
      COMMON /CSW/ SW9,SWM,SWC
      REAL*4 DAY,fixLAT
      common /FixInput/  DAY,fixLAT
      REAL*4 GSURF,RE,RGAS,DGTR,DR
      common /PhysCte/   GSURF,RE,RGAS,DGTR,DR
      REAL*4 XMM,ALTL(8),alfa(7)
      INTEGER*4 MassNo(7),NoPDM(7),NoPD(7)
      common /FixData/   XMM,MassNo,ALTL,alfa,NoPDM,NoPD
      REAL*4 ZN1(5),ZLB,ZA,S,ZGA,ZGDIF1
      common /FixAlt/    ZN1,ZLB,ZA,S,ZGA,ZGDIF1
      REAL*4 TA,TLB,TINF,DTA
      common /FixTemp/   TA,TLB,TINF,DTA
      REAL*4 XS1(5),YS1(5),Y2OUT1(5),YI1
      common /Splin1/    XS1,YS1,Y2OUT1,YI1
      REAL*4 PTM(10),PDM(10,8)
      COMMON /LOWER6/    PTM,PDM
      REAL*4 PT(150),PD(150,9),PS(150),PDL(25,2),PTL(100,4),PMA(100,10)
      COMMON /PARM6/     PT,PD,PS,PDL,PTL,PMA
      SAVE
C
      INTEGER*4 i
      REAL*4 ZETA,ZZ,ZL,GA,SZA,GAMZA,ZG,EXPL,XM,ZH,Xzh,Yzh,Tzh,YIzhZa
      REAL*4 DdZH,DmZH,RL,HC01,ZC01,HCC01,ZCC01,RC01
      REAL*4 HC04,ZC04,HC32,ZC32,HC40,ZC40
      REAL*4 HC14,ZC14,HCC14,ZCC14,RC14
      REAL*4 HC16,ZC16,HCC16,ZCC16,RC16
      REAL*4 DdiffZA,BatesDens,SplinDens,SPLINT,SPLINI,DNET,CCOR
C     
      real*4 v0
C                          Internal fct
      ZETA(ZZ,ZL) = (ZZ-ZL)*(RE+ZL) / (RE+ZZ)
C
      data v0/0./
C
C        Preparing mixed density from Bates profile ( cf BatesDEns ) beginning
C        at ZA = 122.8 km instead of ZLB = 120 km
      GA = GSURF / (1.+ZA/RE)**2    
      SZA = DTA/(TINF-TA)
      GAMZA = XMM*GA / (SZA*RGAS*TINF)
      ZG = ZETA(ALT,ZA)
      EXPL = SZA*GAMZA*ZG
      if (EXPL.gt.50.) EXPL = 50.
C-.......................................................................
C        BIG LOOP ON 7 DIFFERENT DENSITIES TO COMPUTE       
C-.......................................................................=
      do i=1,7
        XM =  ( MassNo(i) )
        if ( MASS.eq.48 .or. MASS.eq.MassNo(i) .or.
     +             (ALT.le.ALTL(8) .and. MassNo(i).eq.28) ) then
          DdiffZA = BatesDens (ZGA,XM,TA,DdLB(i),alfa(i))
          if (Level.eq.1) 
     +       Ddiff(i) = SplinDens (ZA,ZGDIF1,XM,TA,TZ,DdiffZA,
     +                                   alfa(i),YI1)
C................ MIXED DENSITIES CALCULATIONS ..........................=
          if ( SWM(15).ne.0 .and.
     +             ( ALT.le.ALTL(i) .or. MassNo(i).eq.28 ) ) then
C              Turbopause height for each specy, variation for N2 turbopause
            ZH = PDM(3,NoPDM(i))
            if (MassNo(i).eq.28)
     +        ZH = ZH * PDL(25,2) * ( 1. + SWM(23)*SWC(5)*PDL(25,1)*
     +                    SIN(DGTR*fixLAT) * COS(DR*(DAY-PT(14))) )
C              Temperature at turbopause ZH
            Xzh = ZETA(ZH,ZN1(1)) / ZGDIF1
            Yzh = SPLINT (XS1,YS1,Y2OUT1,MN1,Xzh)
            Tzh = 1./Yzh
C              Diffusive = mixed density at ZH            
            YIzhZa = SPLINI (XS1,YS1,Y2OUT1,MN1,Xzh)
            DdZH = SplinDens (ZA,ZGDIF1,XM,TA,Tzh,DdiffZA,
     +                                   alfa(i),YIzhZA)
            DmZH = DdZH
C              Mixed density at ZA = 122.8 km
            DmixZA(i) =SplinDens(ZA,ZGDIF1,XMM,Tzh,TA,DmZH,v0,-YIzhZa)
C              Mixed density at ALT 
            if (Level.eq.0) then
C              Computing mixed density from Bates profile as in BatesDens,
C              but beginning at ZA = 122.8, not ZLB = 120 km   
              Dmix(i) = DmixZA(i) * (TA/TZ)**(1.+GAMZA) * EXP(-EXPL)
             else
              Dmix(i) =SplinDens(ZA,ZGDIF1,XMM,TA,TZ,DmixZA(i),v0,YI1)
            endif
            if (SWM(15).eq.-1) then
              D(i) = Dmix(i)
             else
              D(i) = DNET (ALT,Ddiff(i),Dmix(i),XMM,XM)      
            endif
           else
            D(i) = Ddiff(i)
            if ( SWM(15).eq.-1 .and. ALT.gt.ALTL(i) ) D(i) = 0.0
          endif
        endif
      enddo
C-.......................................................................=
C        CORRECTIONS TO SPECIFIED MIXING RATIOS AT GROUND (A20b, A20a) , 
C        CHEMISTRY CORRECTIONS FOR H, N & O (A21)
C-.......................................................................=
C          H - 01 - D(1)
      if ( ALT.le.ALTL(1) .and. SWM(15).ne.0 .and.
     +            (MASS.eq.48 .or. MASS.eq.MassNo(1)) ) then
        RL = LOG ( DmixZA(7)*PDM(2,6)*ABS(PDL(18,2)) / DmixZA(1) )
        HC01 = PDM(6,6)*PDL(12,2)
        ZC01 = PDM(5,6)*PDL(11,2)
        D(1) = D(1) * CCOR(ALT,RL,HC01,ZC01)
        HCC01 = PDM(8,6)*PDL(20,2)
        ZCC01 = PDM(7,6)*PDL(19,2)
        RC01 = PDM(4,6)*PDL(21,2)
        D(1) = D(1) * CCOR(ALT,RC01,HCC01,ZCC01)
      endif
C          N - 14 - D(2)
      if (ALT.le.ALTL(2) .and. SWM(15).ne.0 .and.
     +            (MASS.eq.48 .or. MASS.eq.MassNo(2)) ) then
        RL = LOG ( DmixZA(7)*PDM(2,7)*ABS(PDL(3,1)) / DmixZA(2) )
        HC14 = PDM(6,7)*PDL(2,1)
        ZC14 = PDM(5,7)*PDL(1,1)
        D(2) = D(2) * CCOR(ALT,RL,HC14,ZC14)
        HCC14 = PDM(8,7)*PDL(5,1)
        ZCC14 = PDM(7,7)*PDL(4,1)
        RC14 = PDM(4,7)*PDL(6,1)
        D(2) = D(2) * CCOR(ALT,RC14,HCC14,ZCC14)
      endif
C          O - 16 - D(3)
      if (ALT.le.ALTL(3) .and. SWM(15).ne.0 .and.
     +            (MASS.eq.48 .or. MASS.eq.MassNo(3)) ) then
        RL = LOG ( DmixZA(7)*PDM(2,2)*ABS(PDL(17,2)) / DmixZA(3) )
        HC16 = PDM(6,2)*PDL(4,2)
        ZC16 = PDM(5,2)*PDL(3,2)
        D(3) = D(3) * CCOR(ALT,RL,HC16,ZC16)
        HCC16 = PDM(8,2)*PDL(14,2)
        ZCC16 = PDM(7,2)*PDL(13,2)
        RC16 = PDM(4,2)*PDL(15,2)
        D(3) = D(3) * CCOR(ALT,RC16,HCC16,ZCC16)
      endif
C          He - 04 - D(4)
      if (ALT.le.ALTL(4) .and. SWM(15).ne.0 .and.
     +            (MASS.eq.48 .or. MASS.eq.MassNo(4)) ) then
        RL = LOG ( DmixZA(7)*PDM(2,1) / DmixZA(4) )
        ZC04 = PDM(5,1)*PDL(1,2)
        HC04 = PDM(6,1)*PDL(2,2)
        D(4) = D(4) * CCOR(ALT,RL,HC04,ZC04)
      endif
C          O2 - 32 - D(5)
      if (ALT.le.ALTL(5) .and. SWM(15).ne.0 .and.
     +            (MASS.eq.48 .or. MASS.eq.MassNo(5)) ) then
        RL = LOG ( DmixZA(7)*PDM(2,4)/ DmixZA(5) )
        HC32 = PDM(6,4)*PDL(8,2)
        ZC32 = PDM(5,4)*PDL(7,2)
        D(5) = D(5) * CCOR(ALT,RL,HC32,ZC32)
      endif
C          Ar - 40 - D(6)
      if (ALT.le.ALTL(6) .and. SWM(15).ne.0 .and.
     +            (MASS.eq.48 .or. MASS.eq.MassNo(6)) ) then
        RL = LOG ( DmixZA(7)*PDM(2,5) / DmixZA(6) )
        HC40 = PDM(6,5)*PDL(10,2)
        ZC40 = PDM(5,5)*PDL(9,2)
        D(6) = D(6) * CCOR(ALT,RL,HC40,ZC40)
      endif
C
      return   
      end
C-----------------------------------------------------------------------
      SUBROUTINE setswitch(IW)
      implicit none
C      Using switches IW given by the driver, this routine defines for i=1,24 :
C      SWM FOR MAIN TERMS ON (1) or OFF (0)
C      SWC FOR CROSS TERMS ON (1) or OFF (0)
C       IW(i)=0 : Main terms & Cross terms OFF
C       IW(i)=1 : Main terms & Cross terms ON
C       IW(i)=2 : Main terms OFF & Cross terms ON
C       IW(i)=3 : Main terms ON & Cross terms OFF
C           Replaces Hedin's subroutine TSELEC
      INTEGER*4 IW(26)
      INTEGER*4 SW9,SWM(26),SWC(24)
      COMMON /CSW/ SW9,SWM,SWC
      save
C
      INTEGER*4 i
C
      if (IW(25).ne.1) then
        do i=1,24
          IW(i) = 1
        enddo
      endif
      do i = 1,24
        if ((IW(i).ne.0).and.(IW(i).ne.1).and.(IW(i).ne.2)
     +            .and.(IW(i).ne.3).and.(IW(i).ne.-1)) then
          write(6,*) ' WARNING !! IW(',i,') GIVEN AS : ',IW(i),
     +     ' AND IW(25) WAS SET TO 1!! (WRONG IW SET TO 1) !'
          IW(i) = 1
        endif
        SWM(i)=1
        if (IW(i).eq.0 .or. IW(i).eq.2) SWM(i)=0
        SWC(i)=1
        if (IW(i).eq.0 .or. IW(i).eq.3) SWC(i)=0
      enddo
      SW9 = IW(9)
      if (SW9.eq.-1) SWM(9) = 1
      if (IW(15).eq.-1) SWM(15) = -1
      SWM(26) = 0
      if (IW(26).eq.1) SWM(26) = 1
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE GLATF(LAT,GV,REFF)
      implicit none
C      CALCULATE LATITUDE VARIABLE GRAVITY (GV) AND EFFECTIVE
C      RADIUS (REFF)
      REAL*4 LAT,LATL,GV,REFF
      INTEGER*4 SW9,SWM(26),SWC(24)
      COMMON /CSW/ SW9,SWM,SWC
      SAVE
C
      REAL*4 DGTR,a,b,C1,C2,S1
      DATA DGTR/1.74533E-2/,LATL/-999./
      data a/6378.137/,b/6356.752/
C
      IF( LAT.NE.LATL .or. SWM(26).eq.1 ) then
        C2 = COS(2.*DGTR*LAT)
        C1 = COS(DGTR*LAT)
        S1 = SIN(DGTR*LAT)
      endif
      LATL=LAT
      GV = 980.616*(1.-.0026373*C2)
      REFF = sqrt(a*a*C1*C1+b*b*S1*S1)
      RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      FUNCTION GLOBE6(DAY,SEC,LAT,LONG,TLOC,F107A,F107,AP,P)
      implicit none
C       CALCULATE G(L) FUNCTION =
C       Upper Thermosphere Parameters
C
      REAL*4 GLOBE6,DAY,SEC,LAT,LONG,TLOC,F107A,F107
      REAL*4 LONGL
      REAL*4 P(150),AP(7),T(14)
C      
      REAL*4 GSURF,RE,RGAS,DGTR,DR
      common /PhysCte/ GSURF,RE,RGAS,DGTR,DR
      INTEGER*4 SW9,SWM(26),SWC(24)
      COMMON /CSW/ SW9,SWM,SWC
      REAL*4   PLG(9,4),CTLOC,STLOC,C2TLOC,S2TLOC,C3TLOC,
     $                 S3TLOC,dDAY,DF,DFA,APD,APDF,APT(4),XLONG,
     +                 CLONG,SLONG
      COMMON /LPOLY/   PLG,CTLOC,STLOC,C2TLOC,S2TLOC,C3TLOC,
     $                 S3TLOC,dDAY,DF,DFA,APD,APDF,APT,XLONG,
     +                 CLONG,SLONG
      SAVE
C      
      REAL*4 XL,TLL
      REAL*4 DAYL,P14,P18,P32
      REAL*4 HR,SR,P39
C
      INTEGER*4 I,J
      REAL*4 C,S,C2,C4,S2,CD14,CD18,CD32,CD39,F1,F2
      REAL*4 T71,T72,T81,T82,P44,P45,EXP1,EXP2,TINF
      REAL*4 G0,A,SUMEX,EX,SG0
C
      DATA XL/1000./,TLL/1000./
      DATA DAYL/-1./,P14/-1000./,P18/-1000./,P32/-1000./
      DATA HR/.2618/,SR/7.2722E-5/,P39/-1000./
      DATA LONGL/-999./
C
C                       3hr Magnetic activity functions
C
      G0(A)=(A-4.+(P(26)-1.)*(A-4.+(EXP(-ABS(P(25))*(A-4.))-1.)/ABS(P(25
     *))))
       SUMEX(EX)=1.+(1.-EX**19)/(1.-EX)*EX**(.5)
      SG0(EX)=(G0(AP(2))+(G0(AP(3))*EX+G0(AP(4))*EX*EX+G0(AP(5))*EX**3
     $ +(G0(AP(6))*EX**4+G0(AP(7))*EX**12)*(1.-EX**8)/(1.-EX))
     $ )/SUMEX(EX)
C
      dDAY = DAY
      DO 10 J=1,14
       T(J)=0.
   10 CONTINUE
      XLONG=LONG
      IF (XL.EQ.LAT .and. SWM(26).ne.1)   GO TO 15
C          CALCULATE LEGENDRE POLYNOMIALS
      C = SIN(LAT*DGTR)
      S = COS(LAT*DGTR)
      C2 = C*C
      C4 = C2*C2
      S2 = S*S
      PLG(2,1) = C
      PLG(3,1) = 0.5*(3.*C2 -1.)
      PLG(4,1) = 0.5*(5.*C*C2-3.*C)
      PLG(5,1) = (35.*C4 - 30.*C2 + 3.)/8.
      PLG(6,1) = (63.*C2*C2*C - 70.*C2*C + 15.*C)/8.
      PLG(7,1) = (11.*C*PLG(6,1) - 5.*PLG(5,1))/6.
C     PLG(8,1) = (13.*C*PLG(7,1) - 6.*PLG(6,1))/7.
      PLG(2,2) = S
      PLG(3,2) = 3.*C*S
      PLG(4,2) = 1.5*(5.*C2-1.)*S
      PLG(5,2) = 2.5*(7.*C2*C-3.*C)*S
      PLG(6,2) = 1.875*(21.*C4 - 14.*C2 +1.)*S
      PLG(7,2) = (11.*C*PLG(6,2)-6.*PLG(5,2))/5.
C     PLG(8,2) = (13.*C*PLG(7,2)-7.*PLG(6,2))/6.
C     PLG(9,2) = (15.*C*PLG(8,2)-8.*PLG(7,2))/7.
      PLG(3,3) = 3.*S2
      PLG(4,3) = 15.*S2*C
      PLG(5,3) = 7.5*(7.*C2 -1.)*S2
      PLG(6,3) = 3.*C*PLG(5,3)-2.*PLG(4,3)
      PLG(7,3)=(11.*C*PLG(6,3)-7.*PLG(5,3))/4.
      PLG(8,3)=(13.*C*PLG(7,3)-8.*PLG(6,3))/5.
      PLG(4,4) = 15.*S2*S
      PLG(5,4) = 105.*S2*S*C 
      PLG(6,4)=(9.*C*PLG(5,4)-7.*PLG(4,4))/2.
      PLG(7,4)=(11.*C*PLG(6,4)-8.*PLG(5,4))/3.
      XL=LAT
   15 CONTINUE
      IF (TLL.EQ.TLOC .and. SWM(26).ne.1)   GO TO 16
      IF(SWM(7).EQ.0.AND.SWM(8).EQ.0.AND.SWM(14).EQ.0) GOTO 16
      STLOC = SIN(HR*TLOC)
      CTLOC = COS(HR*TLOC)
      S2TLOC = SIN(2.*HR*TLOC)
      C2TLOC = COS(2.*HR*TLOC)
      S3TLOC = SIN(3.*HR*TLOC)
      C3TLOC = COS(3.*HR*TLOC)
      TLL = TLOC
   16 CONTINUE
      IF (LONG.NE.LONGL .or. SWM(26).eq.1) THEN
        CLONG=COS(DGTR*LONG)
        SLONG=SIN(DGTR*LONG)
      ENDIF
      LONGL=LONG
      IF ( (DAY.NE.DAYL.OR.P(14).NE.P14) .or. SWM(26).eq.1 )
     +  CD14=COS(DR*(DAY-P(14)))
      IF ( (DAY.NE.DAYL.OR.P(18).NE.P18) .or. SWM(26).eq.1 )
     +  CD18=COS(2.*DR*(DAY-P(18)))
      IF ( (DAY.NE.DAYL.OR.P(32).NE.P32) .or. SWM(26).eq.1 )
     +  CD32=COS(DR*(DAY-P(32)))
      IF ( (DAY.NE.DAYL.OR.P(39).NE.P39) .or. SWM(26).eq.1 )
     + CD39=COS(2.*DR*(DAY-P(39)))
      DAYL = DAY
      P14 = P(14)
      P18 = P(18)
      P32 = P(32)
      P39 = P(39)
C         F10.7 EFFECT
      DF = F107 - F107A
      DFA=F107A-150.
      T(1) =  P(20)*DF + P(21)*DF*DF + P(22)*DFA
     & + P(30)*DFA**2
      F1 = 1. + (P(48)*DFA +P(20)*DF+P(21)*DF*DF)*SWC(1)
      F2 = 1. + (P(50)*DFA+P(20)*DF+P(21)*DF*DF)*SWC(1)
C        TIME INDEPENDENT
      T(2) =
     1  (P(2)*PLG(3,1) + P(3)*PLG(5,1)+P(23)*PLG(7,1))
     & +(P(15)*PLG(3,1))*DFA*SWC(1)
     2 +P(27)*PLG(2,1)
C        SYMMETRICAL ANNUAL
      T(3) =
     1 (P(19) )*CD32
C        SYMMETRICAL SEMIANNUAL
      T(4) =
     1 (P(16)+P(17)*PLG(3,1))*CD18
C        ASYMMETRICAL ANNUAL
      T(5) =  F1*
     1  (P(10)*PLG(2,1)+P(11)*PLG(4,1))*CD14
C         ASYMMETRICAL SEMIANNUAL
      T(6) =    P(38)*PLG(2,1)*CD39
C        DIURNAL
      IF(SWM(7).EQ.0) GOTO 200
      T71 = (P(12)*PLG(3,2))*CD14*SWC(5)
      T72 = (P(13)*PLG(3,2))*CD14*SWC(5)
      T(7) = F2*
     1 ((P(4)*PLG(2,2) + P(5)*PLG(4,2) + P(28)*PLG(6,2)
     2 + T71)*CTLOC
     4 + (P(7)*PLG(2,2) + P(8)*PLG(4,2) +P(29)*PLG(6,2)
     5 + T72)*STLOC)
  200 CONTINUE
C        SEMIDIURNAL
      IF(SWM(8).EQ.0) GOTO 210
      T81 = (P(24)*PLG(4,3)+P(36)*PLG(6,3))*CD14*SWC(5) 
      T82 = (P(34)*PLG(4,3)+P(37)*PLG(6,3))*CD14*SWC(5)
      T(8) = F2*
     1 ((P(6)*PLG(3,3) + P(42)*PLG(5,3) + T81)*C2TLOC
     3 +(P(9)*PLG(3,3) + P(43)*PLG(5,3) + T82)*S2TLOC)
  210 CONTINUE
C        TERDIURNAL
      IF(SWM(14).EQ.0) GOTO 220
      T(14) = F2*
     1 ((P(40)*PLG(4,4)+(P(94)*PLG(5,4)+P(47)*PLG(7,4))*CD14*SWC(5))*
     $ S3TLOC
     2 +(P(41)*PLG(4,4)+(P(95)*PLG(5,4)+P(49)*PLG(7,4))*CD14*SWC(5))*
     $ C3TLOC)
  220 CONTINUE
C          MAGNETIC ACTIVITY BASED ON DAILY AP
      IF(SW9.EQ.-1)  GO TO 30
      APD=(AP(1)-4.)
      P44=P(44)
      P45=P(45)
      IF(P44.LT.0.) P44=1.E-5
      APDF = (APD+(P45-1.)*(APD+(EXP(-P44  *APD)-1.)/P44  ))
      IF(SW9.EQ.0) GOTO 40
      T(9)=APDF*(P(33)+P(46)*PLG(3,1)+P(35)*PLG(5,1)+
     $ (P(101)*PLG(2,1)+P(102)*PLG(4,1)+P(103)*PLG(6,1))*CD14*SWC(5)+
     $ (P(122)*PLG(2,2)+P(123)*PLG(4,2)+P(124)*PLG(6,2))*SWC(7)*
     $ COS(HR*(TLOC-P(125))))
      GO TO 40
   30 CONTINUE
      IF(P(52).EQ.0.)  GO TO 40
      EXP1 = EXP(-10800.*ABS(P(52))/(1.+P(139)*(45.-ABS(LAT))))
      IF(EXP1.GT..99999) EXP1=.99999
      EXP2 = EXP(-10800.*ABS(P(54)))
      IF(EXP2.GT..99999) EXP2=.99999
      IF(P(25).LT.1.E-4) P(25)=1.E-4
      APT(1)=SG0(EXP1)
      APT(3)=SG0(EXP2)
      IF(SW9.EQ.0)  GOTO 40
      T(9) = APT(1)*(P(51)+P(97)*PLG(3,1)+P(55)*PLG(5,1)+
     $ (P(126)*PLG(2,1)+P(127)*PLG(4,1)+P(128)*PLG(6,1))*CD14*SWC(5)+
     $ (P(129)*PLG(2,2)+P(130)*PLG(4,2)+P(131)*PLG(6,2))*SWC(7)*
     $ COS(HR*(TLOC-P(132))))
  40  CONTINUE
      IF(SWM(10).EQ.0 .OR.LONG.LE.-1000.) GO TO 49
C        LONGITUDINAL
      IF(SWM(11).EQ.0) GOTO 230
      T(11)= (1.+P(81)*DFA*SWC(1))*
     $((P(65)*PLG(3,2)+P(66)*PLG(5,2)+P(67)*PLG(7,2)
     $ +P(104)*PLG(2,2)+P(105)*PLG(4,2)+P(106)*PLG(6,2)
     $ +SWC(5)*(P(110)*PLG(2,2)+P(111)*PLG(4,2)+P(112)*PLG(6,2))*CD14)*
     $     CLONG
     $ +(P(91)*PLG(3,2)+P(92)*PLG(5,2)+P(93)*PLG(7,2)
     $ +P(107)*PLG(2,2)+P(108)*PLG(4,2)+P(109)*PLG(6,2)
     $ +SWC(5)*(P(113)*PLG(2,2)+P(114)*PLG(4,2)+P(115)*PLG(6,2))*CD14)*
     $  SLONG)
  230 CONTINUE
C        UT AND MIXED UT,LONGITUDE
      IF(SWM(12).EQ.0) GOTO 240
      T(12)=(1.+P(96)*PLG(2,1))*(1.+P(82)*DFA*SWC(1))*
     $(1.+P(120)*PLG(2,1)*SWC(5)*CD14)*
     $((P(69)*PLG(2,1)+P(70)*PLG(4,1)+P(71)*PLG(6,1))*
     $     COS(SR*(SEC-P(72))))
      T(12)=T(12)+SWC(11)*
     $ (P(77)*PLG(4,3)+P(78)*PLG(6,3)+P(79)*PLG(8,3))*
     $     COS(SR*(SEC-P(80))+2.*DGTR*LONG)*(1.+P(138)*DFA*SWC(1))
  240 CONTINUE
C        UT,LONGITUDE MAGNETIC ACTIVITY
      IF(SWM(13).EQ.0) GOTO 48
      IF(SW9.EQ.-1)  GO TO 45
      T(13)= APDF*SWC(11)*(1.+P(121)*PLG(2,1))*
     $((P( 61)*PLG(3,2)+P( 62)*PLG(5,2)+P( 63)*PLG(7,2))*
     $     COS(DGTR*(LONG-P( 64))))
     $ +APDF*SWC(11)*SWC(5)*
     $ (P(116)*PLG(2,2)+P(117)*PLG(4,2)+P(118)*PLG(6,2))*
     $     CD14*COS(DGTR*(LONG-P(119)))
     $ + APDF*SWC(12)*
     $ (P( 84)*PLG(2,1)+P( 85)*PLG(4,1)+P( 86)*PLG(6,1))*
     $     COS(SR*(SEC-P( 76)))
      GOTO 48
   45 CONTINUE
      IF(P(52).EQ.0.) GOTO 48
      T(13)=APT(1)*SWC(11)*(1.+P(133)*PLG(2,1))*
     $((P(53)*PLG(3,2)+P(99)*PLG(5,2)+P(68)*PLG(7,2))*
     $     COS(DGTR*(LONG-P(98))))
     $ +APT(1)*SWC(11)*SWC(5)*
     $ (P(134)*PLG(2,2)+P(135)*PLG(4,2)+P(136)*PLG(6,2))*
     $     CD14*COS(DGTR*(LONG-P(137)))
     $ +APT(1)*SWC(12)*
     $ (P(56)*PLG(2,1)+P(57)*PLG(4,1)+P(58)*PLG(6,1))*
     $     COS(SR*(SEC-P(59)))
   48 CONTINUE
   49 CONTINUE
      TINF=P(31)
      DO 50 I = 1,14
   50 TINF = TINF + ABS(SWM(I))*T(I)
      GLOBE6 = TINF
      RETURN
      END
C-----------------------------------------------------------------------
      FUNCTION GLOB6S(P)
      implicit none
C             VERSION OF GLOBE FOR LOWER ATMOSPHERE 1/17/90
C
      REAL*4 GLOB6S
      REAL*4 P(100),T(14)
      REAL*4 GSURF,RE,RGAS,DGTR,DR
      common /PhysCte/ GSURF,RE,RGAS,DGTR,DR
      REAL*4   PLG(9,4),CTLOC,STLOC,C2TLOC,S2TLOC,C3TLOC,
     $                 S3TLOC,dDAY,DF,DFA,APD,APDF,APT(4),LONG,
     +                 CLONG,SLONG
      COMMON /LPOLY/   PLG,CTLOC,STLOC,C2TLOC,S2TLOC,C3TLOC,
     $                 S3TLOC,dDAY,DF,DFA,APD,APDF,APT,LONG,
     +                 CLONG,SLONG
      INTEGER*4 SW9,SWM(26),SWC(24)
      COMMON /CSW/ SW9,SWM,SWC
      SAVE
      REAL*4 DAYL,P32,P18,P14,P39
      REAL*4 DAY,CD14,CD18,CD32,CD39,T71,T72,T81,T82,TT
      INTEGER*4 I,J
C
      DATA DAYL/-1./,P32,P18,P14,P39/4*-1000./
C
      DAY = dDAY    
      DO 10 J=1,14
        T(J)=0.
   10 CONTINUE
      IF ( (DAY.NE.DAYL.OR.P32.NE.P(32)) .or. SWM(26).eq.1 )
     +  CD32=COS(DR*(DAY-P(32)))
      IF ( (DAY.NE.DAYL.OR.P18.NE.P(18)) .or. SWM(26).eq.1 )
     +  CD18=COS(2.*DR*(DAY-P(18)))
      IF ( (DAY.NE.DAYL.OR.P14.NE.P(14)) .or. SWM(26).eq.1 )
     +  CD14=COS(DR*(DAY-P(14)))
      IF ( (DAY.NE.DAYL.OR.P39.NE.P(39)) .or. SWM(26).eq.1 )
     +  CD39=COS(2.*DR*(DAY-P(39)))
      DAYL=DAY
      P32=P(32)
      P18=P(18)
      P14=P(14)
      P39=P(39)
C
C       F10.7
      T(1)=P(22)*DFA
C       TIME INDEPENDENT
      T(2)=P(2)*PLG(3,1)+P(3)*PLG(5,1)+P(23)*PLG(7,1)
     $     +P(27)*PLG(2,1)+P(28)*PLG(4,1)+P(29)*PLG(6,1)
C       SYMMETRICAL ANNUAL
      T(3)=(P(19)+P(48)*PLG(3,1)+P(30)*PLG(5,1))*CD32
C       SYMMETRICAL SEMIANNUAL
      T(4)=(P(16)+P(17)*PLG(3,1)+P(31)*PLG(5,1))*CD18
C       ASYMMETRICAL ANNUAL
      T(5)=(P(10)*PLG(2,1)+P(11)*PLG(4,1)+P(36)*PLG(6,1))*CD14
C       ASYMMETRICAL SEMIANNUAL
      T(6)=(P(38)*PLG(2,1))*CD39
C        DIURNAL
      IF(SWM(7).EQ.0) GOTO 200
      T71 = P(12)*PLG(3,2)*CD14*SWC(5)
      T72 = P(13)*PLG(3,2)*CD14*SWC(5)
      T(7) = 
     1 ((P(4)*PLG(2,2) + P(5)*PLG(4,2)
     2 + T71)*CTLOC
     4 + (P(7)*PLG(2,2) + P(8)*PLG(4,2)
     5 + T72)*STLOC)
  200 CONTINUE
C        SEMIDIURNAL
      IF(SWM(8).EQ.0) GOTO 210
      T81 = (P(24)*PLG(4,3)+P(47)*PLG(6,3))*CD14*SWC(5) 
      T82 = (P(34)*PLG(4,3)+P(49)*PLG(6,3))*CD14*SWC(5)
      T(8) = 
     1 ((P(6)*PLG(3,3) + P(42)*PLG(5,3) + T81)*C2TLOC
     3 +(P(9)*PLG(3,3) + P(43)*PLG(5,3) + T82)*S2TLOC)
  210 CONTINUE
C        TERDIURNAL
      IF(SWM(14).EQ.0) GOTO 220
      T(14) = P(40)*PLG(4,4)*S3TLOC
     $ +P(41)*PLG(4,4)*C3TLOC
  220 CONTINUE
C       MAGNETIC ACTIVITY
      IF(SW9.EQ.0) GOTO 40
      IF(SW9.EQ.1)
     $ T(9)=APDF*(P(33)+P(46)*PLG(3,1)*SWC(2))
      IF(SW9.EQ.-1)
     $ T(9)=(P(51)*APT(3)+P(97)*PLG(3,1)*APT(3)*SWC(2))
   40 CONTINUE
      IF(SWM(10).EQ.0.OR.SWM(11).EQ.0.OR.LONG.LE.-1000.) GO TO 49
C        LONGITUDINAL
      T(11)= (1.+PLG(2,1)*(P(81)*SWC(5)*COS(DR*(DAY-P(82)))
     $           +P(86)*SWC(6)*COS(2.*DR*(DAY-P(87))))
     $        +P(84)*SWC(3)*COS(DR*(DAY-P(85)))
     $           +P(88)*SWC(4)*COS(2.*DR*(DAY-P(89))))
     $ *((P(65)*PLG(3,2)+P(66)*PLG(5,2)+P(67)*PLG(7,2)
     $   +P(75)*PLG(2,2)+P(76)*PLG(4,2)+P(77)*PLG(6,2)
     $    )*CLONG
     $  +(P(91)*PLG(3,2)+P(92)*PLG(5,2)+P(93)*PLG(7,2)
     $   +P(78)*PLG(2,2)+P(79)*PLG(4,2)+P(80)*PLG(6,2)
     $    )*SLONG)
   49 CONTINUE
      TT=0.
      DO 50 I=1,14
   50 TT = TT + SWM(I)*T(I)
      GLOB6S=TT
      RETURN
      END
C-----------------------------------------------------------------------
      subroutine TestInput(IYD,GLAT,GLONG,STL,F107A,F107,AP,MASS,IW,
     +                        jump)
      implicit none
C       Based on Hedin's VTST . Tests if geophysical variables or switches
C       other than ALT changed and saves last values.
C       Output : LOGICAL*4 jump

      INTEGER*4 IYD,MASS,IW(26)
      REAL*4 GLAT,GLONG,STL,F107A,F107,AP(7)
      LOGICAL*4 jump
C
      SAVE
      INTEGER*4 IWL(26)
      REAL*4 APL(7)
      INTEGER*4 i,IYDL,MASSL
      REAL*4 GLATL,GLL,STLL,FAL,FL
      DATA IYDL/-999/,GLATL/-999./,GLL/-999./
      DATA STLL/-999./,FAL/-999./,FL/-999./,APL/7*-999./,MASSL/-9/
      DATA IWL/26*-999/
C 
      jump = .false.
      if (IYD.ne.IYDL) goto 10
      if (GLAT.ne.GLATL) goto 10
      if (GLONG.ne.GLL) goto 10
      if (STL.ne.STLL) goto 10
      if (F107A.ne.FAL) goto 10      
      if (F107.ne.FL) goto 10
      do i=1,7
        if (AP(i).ne.APL(i)) goto 10
      enddo
      if (MASS.ne.MASSL) goto 10
      do i=1,26
        if (IW(i).ne.IWL(i)) goto 10
      enddo
      if (IW(26).eq.1) goto 10
      jump = .true.
      goto 20
   10 continue
      IYDL = IYD
      GLATL = GLAT
      GLL = GLONG
      STLL = STL
      FAL = F107A
      FL = F107
      do i=1,7
        APL(i) = AP(i)
      enddo
      MASSL = MASS
      do i=1,26
        IWL(i) = IW(i)
      enddo
   20 continue
      return
      end
C-----------------------------------------------------------------------
      subroutine CalcTemp (Z,MN,ZN,TN,TGN,TZ,XS,YS,Y2OUT,X,ZGDIF)
      implicit none
C        Calculates TZ, temp at Z<ZA=122.8, and array of spline-temp nodes 
C        XS,YS,Y2OUT & spline coeff X,ZGDIF
C
      INTEGER*4 MN
      REAL*4 Z,ZN(MN),TN(MN),TGN(2),TZ,XS(MN),YS(MN),Y2OUT(MN),X,ZGDIF
C
      REAL*4 GSURF,RE,RGAS,DGTR,DR
      common /PhysCte/ GSURF,RE,RGAS,DGTR,DR
      SAVE
C
      INTEGER*4 K
      REAL*4 ZETA,ZZ,ZL,Z1,Z2,T1,T2,ZG,YD1,YD2,Y,SPLINT
C                         Internal fct

      ZETA(ZZ,ZL)=(ZZ-ZL)*(RE+ZL)/(RE+ZZ)
C
      Z1=ZN(1)
      Z2=ZN(MN)
      T1=TN(1)
      T2=TN(MN)
C      Geopotential difference from Z1
      ZG=ZETA(Z,Z1)
      ZGDIF=ZETA(Z2,Z1)
C       Set up spline nodes
      do K=1,MN
        XS(K)=ZETA(ZN(K),Z1)/ZGDIF
        YS(K)=1./TN(K)
      enddo
C        End node derivatives
      YD1=-TGN(1)/(T1*T1)*ZGDIF
      YD2=-TGN(2)/(T2*T2)*ZGDIF*((RE+Z2)/(RE+Z1))**2
C       Calculate spline coefficients : Y2OUT = 2d deriv of YS(k) on XS(k)
      CALL SPLINE(XS,YS,MN,YD1,YD2,Y2OUT)
      X=ZG/ZGDIF
      Y = SPLINT(XS,YS,Y2OUT,MN,X)
C       temperature at altitude Z
      TZ=1./Y
      return
      end
C-----------------------------------------------------------------------
      function BatesDens (ZG,XM,T,D,ALPHA)
      implicit none
C        Calculates diffusive density 
C        by Bates density profile (integration of Bates T profile) 
C        above ZA=122.8km (Level 0). Uses (A.17a),(A.16a),(A.15),A(14a) & 
C        first factors of (A.13)
C        
      REAL*4 BatesDens,ZG,XM,T,D,ALPHA
      REAL*4 GSURF,RE,RGAS,DGTR,DR
      common /PhysCte/ GSURF,RE,RGAS,DGTR,DR
      REAL*4 ZN1(5),ZLB,ZA,S,ZGA,ZGDIF1
      common /FixAlt/    ZN1,ZLB,ZA,S,ZGA,ZGDIF1
      REAL*4 TA,TLB,TINF,DTA
      common /FixTemp/ TA,TLB,TINF,DTA
C
      REAL*4 GLB,GAMMA,EXPL
C
      GLB = GSURF / (1.+ZLB/RE)**2
      GAMMA = XM*GLB / (S*RGAS*TINF)
      EXPL = S*GAMMA*ZG
      if (EXPL.gt.50.) EXPL = 50.
      BatesDens = D * (TLB/T)**(1.+ALPHA+GAMMA) * EXP(-EXPL)
      return
      end
C-----------------------------------------------------------------------
      function SplinDens (Zref,ZGDIF,XM,Tref,T,D,ALPHA,YI)
      implicit none
C        Calculates Diffusive or Mixed ( depends of XM & nature of D ) =

C        density using YI from SPLINI (integration of Spline-based Temperature 
C        profile between ZA & current ALT).
C        Works below ZA=122.8km (Level>0). (A.17b),(A.16b),(A.15)
C        last factors of -(A.13) if diffusive, (A.19) if mixing - profile
C
      REAL*4 SplinDens,Zref,ZGDIF,XM,Tref,T,D,ALPHA,YI
      REAL*4 GSURF,RE,RGAS,DGTR,DR
      common /PhysCte/  GSURF,RE,RGAS,DGTR,DR
      REAL*4 ZN1(5),ZLB,ZA,S,ZGA,ZGDIF1
      common /FixAlt/    ZN1,ZLB,ZA,S,ZGA,ZGDIF1
      REAL*4 TA,TLB,TINF,DTA
      common /FixTemp/  TA,TLB,TINF,DTA
C
      REAL*4 GZref,GAMM,EXPL
C
      GZref = GSURF / (1.+Zref/RE)**2
      GAMM = XM*GZref*ZGDIF / RGAS
      EXPL = GAMM*YI
      if (EXPL.gt.50.) EXPL = 50.
C       Density at altitude
      SplinDens = D * (Tref/T)**(1.+ALPHA) * EXP(-EXPL)
      RETURN
      END
C-----------------------------------------------------------------------
      function TotMassDens(D,MassNo)
      implicit none
C
C     Calculate total Mass Density (g/cm3) from the 7 number densities (cm-3)
C
      REAL*4 TotMassDens,D(8)
      INTEGER*4 MassNo(7)
c
      INTEGER*4 i
      REAL*4 TMD
c
      TMD = 0.
      do i = 1,7
        TMD = TMD + MassNo(i)*D(i)
      enddo   
      TotMassDens =  1.66E-24 * TMD
      return
      end
C-----------------------------------------------------------------------
      SUBROUTINE SPLINE(X,Y,N,YP1,YPN,Y2)
      implicit none
C        CALCULATE 2ND DERIVATIVES OF CUBIC SPLINE INTERP FUNCTION
C        ADAPTED FROM NUMERICAL RECIPES BY PRESS ET AL
C        X,Y: ARRAYS OF TABULATED FUNCTION IN ASCENDING ORDER BY X
C        N: SIZE OF ARRAYS X,Y
C        YP1,YPN: SPECIFIED DERIVATIVES AT X(1) AND X(N); VALUES
C                 >= 1E30 SIGNAL SECOND DERIVATIVE ZERO
C        Y2: OUTPUT ARRAY OF SECOND DERIVATIVES
C
      INTEGER*4 N
      REAL*4 X(N),Y(N),YP1,YPN,Y2(N)
C
      INTEGER*4 NMAX
      parameter (NMAX=5)
      REAL*4 U(NMAX)
      SAVE
C
      INTEGER*4 I,K
      REAL*4 SIG,P,QN,UN
C    
      IF(YP1.GT..99E30) THEN
        Y2(1)=0.
        U(1)=0.
       ELSE
        Y2(1)=-.5
        U(1)=(3./(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      ENDIF
      DO 11 I=2,N-1
        SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
        P=SIG*Y2(I-1)+2.
        Y2(I)=(SIG-1.)/P
        U(I)=(6.*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))
     $    /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
   11 CONTINUE
      IF(YPN.GT..99E30) THEN
        QN=0.
        UN=0.
      ELSE
        QN=.5
        UN=(3./(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      ENDIF
      Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.)
      DO 12 K=N-1,1,-1
        Y2(K)=Y2(K)*Y2(K+1)+U(K)
   12 CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------
      function SPLINT(XA,YA,Y2A,N,X)
      implicit none
C        CALCULATE CUBIC SPLINE INTERP VALUE
C        ADAPTED FROM NUMERICAL RECIPES BY PRESS ET AL.
C        XA,YA: ARRAYS OF TABULATED FUNCTION IN ASCENDING ORDER BY X
C        Y2A: ARRAY OF SECOND DERIVATIVES
C        N: SIZE OF ARRAYS XA,YA,Y2A
C        X: ABSCISSA FOR INTERPOLATION
C
      INTEGER*4 N
      REAL*4 SPLINT,XA(N),YA(N),Y2A(N),X
      SAVE
C
      INTEGER*4 K,KLO,KHI
      REAL*4 H,A,B
C
      KLO=1
      KHI=N
    1 CONTINUE
      IF(KHI-KLO.GT.1) THEN
        K=(KHI+KLO)/2
        IF(XA(K).GT.X) THEN
          KHI=K
        ELSE
          KLO=K
        ENDIF
        GOTO 1
      ENDIF
      H=XA(KHI)-XA(KLO)
      IF(H.EQ.0.) WRITE(6,*) 'BAD XA INPUT TO SPLINT'
      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      SPLINT = A*YA(KLO)+B*YA(KHI)+
     $  ((A*A*A-A)*Y2A(KLO)+(B*B*B-B)*Y2A(KHI))*H*H/6.
      RETURN
      END
C-----------------------------------------------------------------------
      function SPLINI(XA,YA,Y2A,N,X)
      implicit none
C       routine SPLINI similar to (A.14b) 
C       INTEGRATE CUBIC SPLINE FUNCTION FROM XA(1) TO X
C        XA,YA: ARRAYS OF TABULATED FUNCTION IN ASCENDING ORDER BY X
C        Y2A: ARRAY OF SECOND DERIVATIVES
C        N: SIZE OF ARRAYS XA,YA,Y2A
C        X: ABSCISSA ENDPOINT FOR INTEGRATION
C
      INTEGER*4 N
      REAL*4 SPLINI,XA(N),YA(N),Y2A(N),X
      SAVE
C
      INTEGER*4 KLO,KHI
      REAL*4 YI,XX,A,B,A2,B2,H
C
      YI=0.
      KLO=1
      KHI=2
    1 CONTINUE
      IF(X.GT.XA(KLO).AND.KHI.LE.N) THEN
        XX=X
        IF(KHI.LT.N) XX=MIN(X,XA(KHI))
        H=XA(KHI)-XA(KLO)
        A=(XA(KHI)-XX)/H
        B=(XX-XA(KLO))/H
        A2=A*A
        B2=B*B
        YI=YI+((1.-A2)*YA(KLO)/2.+B2*YA(KHI)/2.+
     $     ((-(1.+A2*A2)/4.+A2/2.)*Y2A(KLO)+
     $     (B2*B2/4.-B2/2.)*Y2A(KHI))*H*H/6.)*H
        KLO=KLO+1
        KHI=KHI+1
        GOTO 1
      ENDIF
      SPLINI = YI
      RETURN
      END
C-----------------------------------------------------------------------
      FUNCTION DNET(ALT,DD,DM,XMM,XM)
      implicit none
C       TURBOPAUSE CORRECTION - combining diff & mixed densities
C       Root mean density
C       8/20/80
C          DD - diffusive density
C          DM - full mixed density
C          ZMH has been replaced by 28.
C          XMM - full mixed molecular weight
C          XM  - species molecular weight
C          DNET - combined density
C
      REAL*4 DNET,ALT,DD,DM,XMM,XM
      SAVE
C
      REAL*4 A,YLOG
C
      A=28./(XMM-XM)
      IF(DM.GT.0.0.AND.DD.GT.0.0) GOTO 5
        WRITE(6,*) 'DNET LOG ERROR : DM : ',DM,' ; DD : ',DD,
     +   ' for XM : ',XM,' at ',ALT,' km'
        IF(DD.EQ.0.0.AND.DM.EQ.0.0) DD=1.
        IF(DM.EQ.0.0) GOTO 10
        IF(DD.EQ.0.0) GOTO 20
    5 CONTINUE
      YLOG=A*LOG(DM/DD)
      IF(YLOG.LT.-10.) GO TO 10
      IF(YLOG.GT.10.)  GO TO 20
        DNET=DD*(1.+EXP(YLOG))**(1/A)
        GO TO 50
   10 CONTINUE
        DNET=DD
        GO TO 50
   20 CONTINUE
        DNET=DM
        GO TO 50
   50 CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------
      FUNCTION  CCOR(ALT, R,H1,ZH)
      implicit none
C        CHEMISTRY/DISSOCIATION CORRECTION FOR MSIS MODELS
C        ALT - altitude
C        R - target ratio
C        H1 - transition scale length
C        ZH - altitude of 1/2 R
C
      REAL*4 CCOR,ALT,R,H1,ZH
      SAVE
C
      REAL*4 E,EX
C
      E=(ALT-ZH)/H1
      IF(E.GT.70.) GO TO 20
      IF(E.LT.-70.) GO TO 10
        EX=EXP(E)
        CCOR=R/(1.+EX)
        GO TO 50
   10   CCOR=R
        GO TO 50
   20   CCOR=0.
        GO TO 50
   50 CONTINUE
      CCOR=EXP(CCOR)
       RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE GHP6(IYD,ALT,GLAT,GLONG,STL,F107A,F107,AP,
     $  D,T,PRESS,IW)
      implicit none
C       FIND ALTITUDE OF PRESSURE SURFACE (PRESS) FROM GTD6
C        PRESS - PRESSURE LEVEL(MB) ; MASS fixed to 48
C     OUTPUT:
C        ALT - ALTITUDE(KM)
C
      INTEGER*4 IYD
      REAL*4 ALT,GLAT,GLONG,STL,F107A,F107,PRESS
      REAL*4 GSURF,RE,RGAS,DGTR,DR
      common /PhysCte/  GSURF,RE,RGAS,DGTR,DR
      REAL*4 D(8),T(2),AP(7)
      INTEGER*4 IW(26)
      SAVE
      INTEGER*4 IDAY,L
      REAL*4 BM,TEST,PL,Z,ZI,CA,CD,CL,CL2,XM,XN,P,DIFF,G,SH
      DATA BM/1.3806E-19/
      DATA TEST/.00043/
C
      PL=LOG10(PRESS)
C      Initial altitude estimate
      IF(PL.GE.-5.) THEN
         IF(PL.GT.2.5) ZI=18.06*(3.00-PL)
         IF(PL.GT..75.AND.PL.LE.2.5) ZI=14.98*(3.08-PL)
         IF(PL.GT.-1..AND.PL.LE..75) ZI=17.8*(2.72-PL)
         IF(PL.GT.-2..AND.PL.LE.-1.) ZI=14.28*(3.64-PL)
         IF(PL.GT.-4..AND.PL.LE.-2.) ZI=12.72*(4.32-PL)
         IF(PL.LE.-4.) ZI=25.3*(.11-PL)
         IDAY=MOD(IYD,1000)
         CL=GLAT/90.
         CL2=CL*CL
         IF(IDAY.LT.182) CD=1.-IDAY/91.25
         IF(IDAY.GE.182) CD=IDAY/91.25-3.
         CA=0.
         IF(PL.GT.-1.11.AND.PL.LE.-.23) CA=1.0
         IF(PL.GT.-.23) CA=(2.79-PL)/(2.79+.23)
         IF(PL.LE.-1.11.AND.PL.GT.-3.) CA=(-2.93-PL)/(-2.93+1.11)
         Z=ZI-4.87*CL*CD*CA-1.64*CL2*CA+.31*CA*CL
      ENDIF
      IF(PL.LT.-5.) Z=22.*(PL+4.)**2+110.
      L=0
C      ITERATION LOOP
   10 CONTINUE
        L=L+1
        CALL NWGTD6(IYD,Z,GLAT,GLONG,STL,F107A,F107,AP,48,D,T,IW)
        XN=D(1)+D(2)+D(3)+D(4)+D(5)+D(6)+D(7)
        P=BM*XN*T(2)
        DIFF=PL-LOG10(P)
        IF(ABS(DIFF).LT.TEST .OR. L.EQ.6) GOTO 20
        XM=D(8)/XN/1.66E-24
        G=GSURF/(1.+Z/RE)**2
        SH=RGAS*T(2)/(XM*G)
C         New altitude estimate using scale height
        Z=Z-SH*DIFF*2.302
        GOTO 10
   20 CONTINUE
      IF(L.EQ.6) WRITE(6,100) PRESS,DIFF
  100 FORMAT(1X,29HGHP6 NOT CONVERGING FOR PRESS,1PE12.2,E12.2)
      ALT=Z
      RETURN
      END
C-----------------------------------------------------------------------
      BLOCK DATA GTD6BK
C          MSISE 90 12-MAR-90   &   ASEM (SIMON CHABRILLAT) DEC 94 
      REAL*4 PT1(50),PT2(50),PT3(50),PA1(50),PA2(50),PA3(50),
     $       PB1(50),PB2(50),PB3(50),PC1(50),PC2(50),PC3(50),
     $       PD1(50),PD2(50),PD3(50),PE1(50),PE2(50),PE3(50),
     $       PF1(50),PF2(50),PF3(50),PG1(50),PG2(50),PG3(50),
     $       PH1(50),PH2(50),PH3(50),PI1(50),PI2(50),PI3(50),
     $       PJ1(50),PJ2(50),PJ3(50),PK1(50),PL1(50),PL2(50),
     $       PM1(50),PM2(50),PN1(50),PN2(50),PO1(50),PO2(50),
     $       PP1(50),PP2(50),PQ1(50),PQ2(50),PR1(50),PR2(50),
     $       PS1(50),PS2(50),PU1(50),PU2(50),PV1(50),PV2(50),
     $       PW1(50),PW2(50),PX1(50),PX2(50),PY1(50),PY2(50),
     $       PZ1(50),PZ2(50)
      COMMON/PARM6/PT1,PT2,PT3,PA1,PA2,PA3,
     $             PB1,PB2,PB3,PC1,PC2,PC3,
     $             PD1,PD2,PD3,PE1,PE2,PE3,
     $             PF1,PF2,PF3,PG1,PG2,PG3,
     $             PH1,PH2,PH3,PI1,PI2,PI3,
     $             PJ1,PJ2,PJ3,PK1,PL1,PL2,
     $             PM1,PM2,PN1,PN2,PO1,PO2,
     $             PP1,PP2,PQ1,PQ2,PR1,PR2,
     $             PS1,PS2,PU1,PU2,PV1,PV2,
     $             PW1,PW2,PX1,PX2,PY1,PY2,
     $             PZ1,PZ2
      REAL*4 PTM(10),PDM(10,8)
      COMMON/LOWER6/PTM,PDM
      REAL*4 PAVGM(10)
      COMMON/MAVG6/PAVGM
C         TEMPERATURE
      DATA PT1/
     *  9.96040E-01, 3.85528E-02, 3.03445E-03,-1.05531E-01,-6.07134E-03,
     * -5.16278E-04,-1.15622E-01, 2.02240E-03, 9.90156E-03,-1.27371E-01,
     * -3.02449E-02, 1.23512E-02,-5.26277E-03,-8.45398E+00, 0.00000E+00,
     *  1.42370E-02, 0.00000E+00, 1.25818E+02, 8.05486E-03, 1.64419E-03,
     * -6.21452E-06, 3.11701E-03, 0.00000E+00, 3.86578E-03, 1.32397E-01,
     *  2.13315E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00,-6.41110E-06,
     *  0.00000E+00, 3.00150E+01, 5.33297E-03, 3.89146E-03, 2.04725E-03,
     *  0.00000E+00, 0.00000E+00,-1.92645E-02, 2.75905E+00, 1.47284E-03,
     *  3.41345E-04,-1.17388E-03,-3.54589E-04, 1.13139E-01, 1.69134E-01,
     *  5.08295E-03, 3.65016E-05, 4.26385E-03, 1.15102E-04, 5.11819E-03/
      DATA PT2/
     *  6.09108E-03, 4.04995E-05, 1.53049E-03, 2.41470E-05, 2.30764E-03,
     *  1.55267E-03, 1.33722E-03,-1.82318E-03,-2.63007E+02, 0.00000E+00,
     *  1.37337E-03, 9.95774E-04, 0.00000E+00,-1.08983E+02, 5.62606E-03,
     *  5.94053E-03, 1.09358E-03, 0.00000E+00,-1.33410E-02,-2.43409E-02,
     * -1.35688E-02, 3.11370E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     * -2.83023E+03, 8.45583E-04, 5.38706E-04, 0.00000E+00, 2.47956E+02,
     *  2.92246E-03, 0.00000E+00, 0.00000E+00, 7.47703E-05, 8.87993E-04,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     * -1.16540E-02,-4.49173E-03,-3.53189E-04,-1.73933E-04,-1.53218E-04,
     * -5.65411E-01, 7.77272E-03,-9.11784E+01, 6.45187E-04, 0.00000E+00/
      DATA PT3/
     * -8.37685E-04, 2.42318E-03, 4.73796E-03,-3.01801E-03,-4.23564E-03,
     * -2.48289E-03, 9.19286E-04, 2.16372E-03, 8.63968E-04, 1.89689E-03,
     *  4.15654E-03, 0.00000E+00, 1.18068E-02, 3.31190E-03, 0.00000E+00,
     *  1.20222E-03, 0.00000E+00, 0.00000E+00,-3.07246E+00, 0.00000E+00,
     *  0.00000E+00, 6.72403E-04, 1.08930E-03, 9.72278E-04, 4.68242E+00,
     * -3.15034E-04, 4.00059E-03, 5.15036E-03, 1.62989E-03, 1.08824E-03,
     *  9.95261E-04, 4.18955E+00,-3.64059E-01, 1.70182E-03, 0.00000E+00,
     *  0.00000E+00,-3.20120E+00, 0.00000E+00, 5.80206E-03, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C         HE DENSITY
      DATA PA1/
     *  1.04934E+00,-2.88362E-02,-2.07095E-01,-1.03314E-01,-7.02373E-03,
     *  1.29664E-02, 4.08853E-01,-9.19895E-03,-1.88660E-02, 1.40927E+00,
     *  1.75033E-01, 1.87351E-02, 1.10979E-01,-7.42871E+00, 0.00000E+00,
     *  2.67143E-01,-5.95979E-02, 1.05038E+02,-8.40963E-02,-6.97632E-04,
     *  2.06521E-06, 7.65306E-04, 0.00000E+00, 0.00000E+00, 1.26762E-01,
     *  1.28876E-01,-5.04479E-02,-1.30735E-02,-2.24348E-02, 0.00000E+00,
     *  0.00000E+00,-1.50832E+02,-6.29928E-03, 0.00000E+00,-4.07760E-03,
     *  0.00000E+00, 0.00000E+00, 5.25725E-02,-3.11486E+01,-3.13351E-03,
     *  2.75838E-03, 0.00000E+00, 0.00000E+00, 1.11247E-01, 1.08815E-01,
     * -4.66713E-02, 0.00000E+00,-3.29329E-03, 0.00000E+00, 1.67838E-03/
      DATA PA2/
     * -9.16691E-03, 3.45044E-05,-9.71806E-03, 0.00000E+00,-2.04672E-03,
     * -7.86899E-03,-7.98285E-03, 5.36515E-03,-5.31172E+03, 0.00000E+00,
     * -6.42781E-03,-1.71690E-03, 0.00000E+00,-6.79131E+01,-1.79912E-02,
     * -1.58305E-02,-7.12313E-03, 0.00000E+00, 2.53477E-02, 8.52960E-02,
     *  1.02163E-01, 2.95009E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     * -6.84625E+03,-6.19098E-03,-2.69289E-03, 0.00000E+00,-5.20231E+02,
     * -6.33463E-03, 0.00000E+00, 0.00000E+00,-6.02428E-03,-4.07077E-03,
     *  5.42264E-03, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  4.07560E-02, 2.82288E-02, 9.08088E-03, 0.00000E+00, 0.00000E+00,
     * -4.05204E-01,-5.97931E-02,-7.31823E+01,-2.06620E-03, 0.00000E+00/
      DATA PA3/
     * -3.72723E-03,-1.88146E-02,-1.01794E-02, 8.04633E-03, 1.01090E-02,
     *  8.73253E-03, 2.38268E-02, 4.80444E-03, 1.71088E-03, 3.96369E-02,
     * -2.13809E-02, 0.00000E+00,-1.02588E-01,-5.91702E-03, 0.00000E+00,
     *  2.70923E-03, 0.00000E+00, 0.00000E+00,-1.75043E+02, 6.03489E-01,
     * -6.17589E-01, 8.38098E-03, 1.83871E-03,-7.05329E-04,-4.06644E+00,
     * -5.09347E-03,-2.84344E-02,-1.24160E-02, 1.33665E-02, 3.93410E-03,
     * -5.03723E-04,-4.57683E+00,-5.29542E-01,-4.25812E-03, 0.00000E+00,
     *  0.00000E+00, 1.91541E+01, 0.00000E+00, 3.23247E-03, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C         O DENSITY
      DATA PB1/
     *  9.31113E-01,-1.38721E-01,-1.33457E-01,-5.29542E-02,-4.44983E-03,
     *  1.35264E-02, 5.98075E-02,-3.62880E-02,-3.12798E-02, 3.72068E-01,
     *  2.95974E-02, 1.20509E-02, 5.21995E-02,-7.78888E+00, 0.00000E+00,
     *  1.18634E-01,-2.04495E-02, 1.03280E+02, 9.82432E-02, 4.77694E-04,
     *  0.00000E+00, 2.74372E-03, 0.00000E+00, 0.00000E+00, 7.57809E-02,
     *  1.71403E-01,-1.05205E-02, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-8.73348E+00,-5.81094E-03, 0.00000E+00,-8.14944E-03,
     *  0.00000E+00, 0.00000E+00, 5.17255E-02,-1.53028E+01,-3.48932E-03,
     *  9.61771E-04, 5.57732E-03,-4.54180E-04, 9.88213E-02, 9.40456E-02,
     * -3.18797E-02, 0.00000E+00, 0.00000E+00, 0.00000E+00, 2.32122E-03/
      DATA PB2/
     * -6.00220E-03, 2.77654E-05,-3.22019E-03, 0.00000E+00,-3.78551E-03,
     * -3.34809E-03,-1.70668E-03, 0.00000E+00, 6.36184E+03, 0.00000E+00,
     *  1.59986E-03,-3.88204E-03,-1.64825E-03,-7.47955E+01,-1.05360E-02,
     * -9.45723E-03,-1.59824E-03,-7.06730E-04,-1.68513E-02,-1.13023E-01,
     * -6.36637E-02,-1.37709E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     * -1.52368E+04,-5.86061E-03,-2.53108E-03, 0.00000E+00,-2.54837E+03,
     * -3.28988E-03, 0.00000E+00, 0.00000E+00,-2.76364E-03, 9.67923E-03,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  4.34255E-02, 1.14020E-02,-6.18447E-03, 0.00000E+00, 0.00000E+00,
     * -3.02568E-01,-3.27694E-02,-6.71589E+01,-2.28340E-03, 0.00000E+00/
      DATA PB3/
     *  3.06230E-03,-4.65113E-03,-9.73421E-03, 1.28326E-02, 7.88553E-03,
     *  7.97197E-03,-1.20760E-02,-7.67547E-03,-1.20755E-03,-2.98523E-02,
     * -1.26560E-02, 0.00000E+00,-5.68350E-02,-1.53039E-02, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 2.42911E-03,-4.01347E-03,-2.19074E-03, 3.11281E+00,
     *  3.23251E-03,-6.39523E-03,-6.63069E-03,-3.04403E-04,-4.01920E-03,
     * -1.18708E-03, 4.15211E+00,-2.01896E-01, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C         N2 DENSITY
      DATA PC1/
     *  1.06903E+00, 0.00000E+00, 0.00000E+00, 3.66210E-03, 0.00000E+00,
     *  1.90412E-02,-1.78929E-03, 0.00000E+00,-3.92257E-02,-1.19444E-01,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00,-8.45398E+00, 0.00000E+00,
     *  2.08180E-02, 0.00000E+00, 1.39638E+02, 8.98481E-02, 0.00000E+00,
     *  0.00000E+00, 3.77113E-04, 0.00000E+00, 0.00000E+00, 1.32397E-01,
     *  2.13315E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-2.36325E+01, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-2.43022E-03,
     * -3.99776E-06, 6.32343E-03, 5.48144E-03, 1.13139E-01, 1.69134E-01,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PC2/
     *  0.00000E+00, 2.41470E-05, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PC3/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C         TLB
      DATA PD1/
     *  9.76619E-01, 0.00000E+00, 0.00000E+00,-2.00200E-02, 0.00000E+00,
     * -9.38391E-03,-1.95833E-03, 0.00000E+00, 1.31480E-02,-1.92414E-02,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00,-8.45398E+00, 0.00000E+00,
     *  1.07674E-02, 0.00000E+00, 8.93820E+01, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 5.68478E-04, 0.00000E+00, 0.00000E+00, 1.32397E-01,
     *  2.13315E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 4.66814E-03, 0.00000E+00, 0.00000E+00,
     *  5.11651E-05, 2.55717E-03, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-2.60147E-03,-8.08556E-04, 1.13139E-01, 1.69134E-01,
     *  6.64196E-03, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PD2/
     *  5.82026E-03, 2.41470E-05, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 6.21998E-03, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PD3/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C         O2 DENSITY
      DATA PE1/
     *  9.31402E-01, 1.37976E-01, 0.00000E+00, 3.23736E-04, 0.00000E+00,
     * -9.10906E-03, 7.07506E-02, 0.00000E+00,-5.16650E-02, 6.89755E-02,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00,-8.45398E+00, 0.00000E+00,
     *  2.81140E-02, 0.00000E+00, 7.36009E+01, 5.96604E-02, 0.00000E+00,
     *  0.00000E+00,-1.51792E-03, 0.00000E+00, 0.00000E+00, 1.32397E-01,
     *  2.13315E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 9.48758E+00, 8.84541E-03, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 1.13139E-01, 1.69134E-01,
     *  1.45192E-02, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PE2/
     *  1.07906E-02, 2.99942E-05, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-1.48930E-02,
     * -7.87184E-03, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     * -6.83420E-02,-4.41778E-02, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 2.29730E-02, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PE3/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C         AR DENSITY
      DATA PF1/
     *  8.68053E-01, 2.36364E-01, 1.34306E-01, 1.03086E-02, 0.00000E+00,
     * -3.79164E-03,-1.57806E-01, 0.00000E+00,-5.87644E-02,-3.12508E-01,
     *  0.00000E+00, 4.37387E-02,-3.54091E-02,-2.23636E+01, 0.00000E+00,
     * -5.33976E-02, 0.00000E+00, 1.14091E+02, 5.17497E-02, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 1.32397E-01,
     *  2.13315E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 3.42702E+02, 1.57033E-02, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-3.66278E-03,
     * -1.16193E-03, 0.00000E+00, 0.00000E+00, 1.13139E-01, 1.69134E-01,
     *  1.78431E-02, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PF2/
     *  1.62864E-02, 3.16963E-05, 1.27968E-02, 0.00000E+00, 0.00000E+00,
     * -7.04599E-03, 2.07921E-03, 6.36660E-03, 2.29940E+04, 0.00000E+00,
     *  1.27833E-02,-2.08036E-03,-4.61820E-03,-6.29391E+01,-1.20745E-02,
     *  1.36675E-02, 1.36011E-02,-5.37162E-03, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  1.92509E+04, 8.35522E-03, 4.19439E-03, 0.00000E+00, 1.20366E+04,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00,-1.00034E-02,-2.33267E-03,
     *  9.72374E-03, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     * -2.65079E-02,-2.09125E-02,-1.09465E-02, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 2.17252E-02,-7.12385E+01,-1.89428E-03, 0.00000E+00/
      DATA PF3/
     * -6.02006E-03, 1.69058E-02, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 2.90646E-02,
     *  3.48971E-03, 0.00000E+00, 5.01174E-02, 5.50595E-02, 0.00000E+00,
     * -9.55897E-03, 0.00000E+00, 0.00000E+00,-1.51693E+03, 0.00000E+00,
     *  0.00000E+00, 1.29306E-02, 2.69567E-03, 0.00000E+00, 3.92243E+00,
     * -8.47690E-03, 1.16896E-02, 0.00000E+00, 1.48967E-02, 5.44521E-03,
     *  0.00000E+00, 5.64918E+00, 0.00000E+00,-7.72178E-03, 0.00000E+00,
     *  0.00000E+00,-7.34042E+01, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C          H DENSITY
      DATA PG1/
     *  1.27515E+00,-2.10472E-01,-1.77924E-01, 2.18900E-01, 2.88436E-02,
     *  1.90077E-02, 2.91001E-01, 2.17437E-02,-1.05186E-02, 4.36141E-01,
     *  1.07605E-01, 3.30755E-02, 4.00581E-02,-9.58051E+00, 0.00000E+00,
     *  1.54028E-02, 0.00000E+00, 7.34194E+01, 4.96540E-02,-5.95906E-03,
     *  3.84512E-05,-1.36000E-02, 0.00000E+00, 0.00000E+00, 1.32397E-01,
     *  2.13315E-01,-4.16610E-02, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 1.46276E+02,-1.98408E-02, 0.00000E+00, 1.32530E-02,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-1.04687E-04,
     * -1.47562E-03, 0.00000E+00, 0.00000E+00, 1.13139E-01, 1.69134E-01,
     * -1.26913E-02, 0.00000E+00, 0.00000E+00, 0.00000E+00,-6.08370E-03/
      DATA PG2/
     * -2.57587E-02, 3.19022E-05, 0.00000E+00, 0.00000E+00, 1.56644E-02,
     *  1.03640E-02, 1.05771E-03, 0.00000E+00, 3.57949E+03, 0.00000E+00,
     * -1.25672E-03, 1.52783E-03, 1.30518E-03, 7.55558E+00,-9.20341E-03,
     * -2.09142E-02,-1.34106E-02, 0.00000E+00,-4.83312E-02, 8.30900E-02,
     *  9.88009E-02,-1.41148E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     * -1.05513E+03, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 6.73442E-03, 2.01691E-03,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  5.98019E-02, 6.33298E-03,-1.12871E-03, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-1.28604E-02, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PG3/
     * -4.94960E-03,-1.36415E-02,-1.15039E-02, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-5.86860E-03,-1.41732E-03, 2.13697E-03, 2.63845E+00,
     * -8.34186E-03,-1.87336E-02,-1.90870E-02,-8.03810E-03,-2.84279E-03,
     *  2.56722E-03, 1.71429E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C          N DENSITY
      DATA PH1/
     *  5.73587E+01,-3.98747E-01, 0.00000E+00,-5.29554E-01,-5.82186E-03,
     *  7.14177E-02,-6.79279E-01,-1.67715E-01,-6.42434E-02,-2.11569E-01,
     * -1.59922E-01,-1.71024E-04,-1.15885E-01, 6.51603E+00, 0.00000E+00,
     * -1.76683E-01, 6.50395E-02, 1.43504E+00, 9.28208E-02, 5.11662E-03,
     *  0.00000E+00, 9.95121E-03, 0.00000E+00, 0.00000E+00, 1.32397E-01,
     *  2.13315E-01, 1.01451E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 5.67667E+01, 2.38192E-03, 0.00000E+00,-1.88240E-02,
     *  0.00000E+00, 0.00000E+00, 4.76218E-02, 2.35206E+01, 4.75901E-03,
     *  5.76162E-03, 1.51815E-02,-1.92730E-02, 1.13139E-01, 1.69134E-01,
     * -2.88771E-02, 0.00000E+00, 0.00000E+00, 0.00000E+00, 1.18418E-03/
      DATA PH2/
     * -3.68927E-03, 3.14704E-05, 8.82198E-03, 0.00000E+00,-1.92562E-02,
     * -2.58674E-03,-2.19913E-02, 0.00000E+00, 4.38655E+03, 0.00000E+00,
     *  7.60126E-03, 2.59438E-03, 1.72310E-03, 7.79204E+01, 7.97786E-04,
     * -7.70510E-03, 1.90982E-03, 2.72707E-03, 1.01016E-02, 1.16537E-01,
     * -3.12236E-03, 1.39783E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     * -1.30712E+03, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00,-3.20544E-03,-2.06970E-02,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  1.59010E-02,-1.91427E-03,-3.42829E-02, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-3.45379E-02, 8.94518E+01, 1.71556E-03, 0.00000E+00/
      DATA PH3/
     * -7.65278E-03,-2.08987E-04,-1.57393E-02, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-8.60673E-03,-1.19922E-02,-6.46356E-03,-3.00107E+00,
     * -9.32511E-03,-1.50205E-02,-8.67835E-03,-7.64801E-03,-1.31495E-02,
     * -6.76720E-03,-1.82396E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C        SPARE
      DATA PI1/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00,-8.45398E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 1.32397E-01,
     *  2.13315E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 1.13139E-01, 1.69134E-01,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PI2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PI3/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C          S PARAM  
      DATA PJ1/
     *  9.51363E-01,-4.67542E-02, 1.20260E-01, 0.00000E+00, 0.00000E+00,
     *  1.91357E-02, 0.00000E+00, 0.00000E+00, 1.25429E-03,-1.33240E-01,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00,-8.45398E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 2.52317E-03, 0.00000E+00,-9.73404E-03, 1.32397E-01,
     *  2.13315E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00,-7.18482E-04, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 7.87683E-03,-2.33698E-03, 1.13139E-01, 1.69134E-01,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PJ2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PJ3/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C          TURBO
      DATA PK1/
     *  9.33804E-01, 5.47446E+00, 1.53263E-01, 9.19303E-01, 1.64109E+01,
     *  4.27083E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 1.40925E-01,
     *  1.15897E+00, 4.71094E-01, 1.09459E+00, 5.25012E+00, 1.00000E+00,
     *  1.00000E+00, 1.03999E+00, 7.67132E-01, 1.10514E+00, 1.75636E+00,
     *  1.10845E+00, 2.33439E+00, 7.96532E-01, 4.31520E+00, 4.07300E+00,
     *  1.22807E+02, 2.39547E-01, 2.53791E-06, 8.42931E-01, 1.04192E+00,
     *  2.00202E+00, 1.00000E+00, 1.00000E+00, 1.00000E+00, 9.62736E-01/
C         LOWER BOUNDARY
      DATA PTM/
     L  1.04130E+03, 3.86000E+02, 1.95000E+02, 1.66728E+01, 2.13000E+02,
     L  1.20000E+02, 2.40000E+02, 1.87000E+02,-2.00000E+00, 0.00000E+00/
      DATA PDM/
     L  2.45600E+07, 6.71072E-06, 1.00000E+02, 0.00000E+00, 1.10000E+02,
     L  1.00000E+01, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
C
     L  8.59400E+10, 5.40000E-01, 1.05000E+02,-8.00000E+00, 1.10000E+02,
     L  1.00000E+01, 9.00000E+01, 2.00000E+00, 0.00000E+00, 0.00000E+00,
C
     L  2.81000E+11, 0.00000E+00, 1.05000E+02, 2.80000E+01, 2.89500E+01,
     L  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
C
     L  3.30000E+10, 2.68270E-01, 1.05000E+02, 0.00000E+00, 1.10000E+02,
     L  1.00000E+01, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
C
     L  1.33000E+09, 1.19615E-02, 1.05000E+02, 0.00000E+00, 1.10000E+02,
     L  1.00000E+01, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
C
     L  1.76100E+05, 1.00000E+00, 9.50000E+01,-8.00000E+00, 1.10000E+02,
     L  1.00000E+01, 9.00000E+01, 2.00000E+00, 0.00000E+00, 0.00000E+00,
C
     L  1.00000E+07, 1.00000E+00, 1.05000E+02,-8.00000E+00, 1.10000E+02,
     L  1.00000E+01, 9.00000E+01, 2.00000E+00, 0.00000E+00, 0.00000E+00,
C
     L  1.00000E+07, 1.00000E+00, 1.05000E+02,-8.00000E+00, 1.10000E+02,
     L  1.00000E+01, 9.00000E+01, 2.00000E+00, 0.00000E+00, 0.00000E+00/
C         TN1(2)
      DATA PL1/
     *  1.02083E+00, 4.08449E-02,-2.34582E-02, 4.38274E-04,-1.52380E-02,
     * -2.09089E-02, 4.46355E-03,-3.41250E-03,-1.12961E-02,-7.03277E-02,
     * -4.82724E-02, 0.00000E+00, 0.00000E+00,-6.20496E+00, 0.00000E+00,
     * -9.80197E-03,-1.45065E-02,-1.13226E+02, 2.28455E-02, 0.00000E+00,
     *  0.00000E+00, 4.93658E-04, 0.00000E+00, 3.79078E-03, 1.32397E-01,
     *  2.13315E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-8.89051E+03, 2.25900E-03, 1.76142E-03, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-2.55015E-04,
     *  2.21388E-03,-5.99073E-04,-3.52331E-03, 1.13139E-01, 1.69134E-01,
     *  7.79156E-03,-1.93458E-03,-1.08596E-02,-4.39285E-04, 0.00000E+00/
      DATA PL2/
     *  3.83994E-03, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 6.76608E-03, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C         TN1(3)
      DATA PM1/
     *  9.24880E-01, 7.41986E-02,-6.37629E-03, 6.00575E-03, 1.29382E-03,
     *  6.97550E-03,-1.70782E-03, 2.80584E-03,-8.87214E-03,-4.35703E-02,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 4.31515E+00, 0.00000E+00,
     * -1.81474E-02,-6.06627E-02,-8.43503E+01, 8.46944E-03, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00,-2.17081E-02,-2.19500E-03, 1.32397E-01,
     *  2.13315E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 2.47580E+02, 4.41585E-03, 7.80466E-03, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 6.44155E-04,
     * -2.49166E-03, 2.90482E-03,-3.40501E-04, 1.13139E-01, 1.69134E-01,
     * -6.01460E-03,-1.63368E-03, 0.00000E+00,-4.31340E-03, 0.00000E+00/
      DATA PM2/
     *  4.53979E-03, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-5.43660E-03, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C         TN1(4)
      DATA PN1/
     *  9.72669E-01,-4.26748E-02, 1.12876E-02,-8.44951E-03, 7.04114E-03,
     *  1.26036E-02,-3.88164E-03,-5.20509E-04,-6.09710E-04, 1.31603E-01,
     *  1.13804E-01, 0.00000E+00, 0.00000E+00,-6.15970E+00, 0.00000E+00,
     * -2.14214E-02,-6.62913E-02,-2.02884E-01, 2.35350E-02, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 1.13573E-02,-1.84905E-03, 1.32397E-01,
     *  2.13315E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 1.42645E+00,-2.64405E-03,-5.57771E-04, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00,-2.20621E+01,-1.10313E-03,
     *  3.97063E-05, 5.47632E-05, 3.57577E-03, 1.13139E-01, 1.69134E-01,
     *  0.00000E+00, 1.18897E-03, 0.00000E+00, 7.62305E-04, 0.00000E+00/
      DATA PN2/
     * -3.52015E-03, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-9.52550E-04,
     *  8.56253E-04, 4.33114E-04, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 1.21223E-03,
     *  2.38694E-04, 9.15245E-04, 1.28385E-03, 8.67668E-04,-5.61425E-06,
     *  1.04445E+00, 3.41112E+01, 0.00000E+00,-8.40704E-01,-2.39639E+02,
     *  7.06668E-01,-2.05873E+01,-3.63696E-01, 2.39245E+01, 1.00000E+01,
     * -1.06657E-03,-7.67292E-04, 1.54534E-04, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C         TN1(5) TN2(1)
      DATA PO1/
     *  9.99368E-01, 4.33893E-02,-2.07009E-03, 1.09617E-03, 1.05440E-03,
     *  4.83408E-04, 9.77040E-04, 9.24791E-04, 4.80247E-04, 4.94737E-02,
     *  1.05985E-03, 0.00000E+00, 0.00000E+00, 2.74409E+00, 0.00000E+00,
     * -4.96656E-03,-1.51684E-02, 4.65158E+01,-7.51133E-03, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 6.63808E-04, 1.32397E-01,
     *  2.13315E-01,-2.06652E-03,-6.32046E-03, 0.00000E+00, 0.00000E+00,
     *  5.94545E-03,-1.90958E+02, 0.00000E+00,-4.16892E-03, 0.00000E+00,
     * -1.67499E-02, 0.00000E+00, 2.58987E-03, 5.97781E+02, 0.00000E+00,
     *  0.00000E+00, 4.44890E-04, 4.66444E-04, 1.13139E-01, 1.69134E-01,
     *  0.00000E+00, 7.11360E-04, 1.32186E-02, 2.23948E-03, 0.00000E+00/
      DATA PO2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 1.60571E-03,
     *  6.28078E-04, 5.05469E-05, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-1.57829E-03,
     * -4.00855E-04, 5.04077E-05,-1.39001E-03,-2.33406E-03,-4.81197E-04,
     *  1.46758E+00, 6.20332E+00, 0.00000E+00, 3.66476E-01,-6.19760E+01,
     *  3.09198E-01,-1.98999E+01, 0.00000E+00,-3.29933E+02, 0.00000E+00,
     * -1.10080E-03,-9.39310E-05, 1.39638E-04, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C          TN2(2)
      DATA PP1/
     *  9.81637E-01,-1.41317E-03, 3.87323E-02, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-3.58707E-02,
     * -8.63658E-03, 0.00000E+00, 0.00000E+00,-2.02226E+00, 0.00000E+00,
     * -8.69424E-03,-1.91397E-02, 8.76779E+01, 4.52188E-03, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-7.07572E-03, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     * -4.11210E-03, 3.50060E+01, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  2.23760E-02, 0.00000E+00,-8.36657E-03, 1.61347E+01, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00,-1.45130E-02, 0.00000E+00, 0.00000E+00/
      DATA PP2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 1.24152E-03,
     *  6.43365E-04, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 1.33255E-03,
     *  2.42657E-03, 1.60666E-03,-1.85728E-03,-1.46874E-03,-4.79163E-06,
     *  1.22464E+00, 3.53510E+01, 0.00000E+00, 4.49223E-01,-4.77466E+01,
     *  4.70681E-01, 8.41861E+00,-2.88198E-01, 1.67854E+02, 0.00000E+00,
     *  7.11493E-04, 6.05601E-04, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C          TN2(3)
      DATA PQ1/
     *  1.00422E+00,-7.11212E-03, 5.24480E-03, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-5.28914E-02,
     * -2.41301E-02, 0.00000E+00, 0.00000E+00,-2.12219E+01, 0.00000E+00,
     * -3.28077E-03, 1.65727E-02, 1.68564E+00,-6.68154E-03, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 8.42365E-03, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-4.34645E-03,-1.03830E-02,-8.08279E-03, 2.16780E-02,
     *  0.00000E+00,-1.38459E+02, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  1.45155E-02, 0.00000E+00, 7.04573E-03,-4.73204E+01, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 1.08767E-02, 0.00000E+00, 0.00000E+00/
      DATA PQ2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 5.21769E-04,
     * -2.27387E-04, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 3.26769E-03,
     *  3.16901E-03, 4.60316E-04,-1.01431E-04, 1.02131E-03, 9.96601E-04,
     *  1.25707E+00, 2.50114E+01, 0.00000E+00, 4.24472E-01,-2.77655E+01,
     *  3.44625E-01, 2.75412E+01, 0.00000E+00, 7.94251E+02, 0.00000E+00,
     *  2.45835E-03, 1.38871E-03, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C          TN2(4) TN3(1)
      DATA PR1/
     *  1.01890E+00,-2.46603E-02, 1.00078E-02, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-6.70977E-02,
     * -4.02286E-02, 0.00000E+00, 0.00000E+00,-2.29466E+01, 0.00000E+00,
     *  2.26580E-03, 2.63931E-02, 3.72625E+01,-6.39041E-03, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-1.85291E-03,-7.47019E-03,-7.07265E-03, 0.00000E+00,
     *  0.00000E+00, 1.39717E+02, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  9.58383E-03, 0.00000E+00, 9.19771E-03,-3.69121E+02, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00,-1.57067E-02, 0.00000E+00, 0.00000E+00/
      DATA PR2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-2.92953E-03,
     * -2.77739E-03,-4.40092E-04, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 2.47280E-03,
     *  2.95035E-04,-1.81246E-03, 2.81945E-03, 4.27296E-03, 9.78863E-04,
     *  1.40545E+00,-6.19173E+00, 0.00000E+00, 0.00000E+00,-7.93632E+01,
     *  4.44643E-01,-4.03085E+02, 0.00000E+00, 1.15603E+01, 0.00000E+00,
     *  2.25068E-03, 8.48557E-04,-2.98493E-04, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C          TN3(2)
      DATA PS1/
     *  9.75801E-01, 3.80680E-02,-3.05198E-02, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 3.85575E-02,
     *  5.04057E-02, 0.00000E+00, 0.00000E+00,-1.76046E+02, 0.00000E+00,
     * -1.48297E-03,-3.68560E-03, 3.02185E+01,-3.23338E-03, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00,-1.15558E-02, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 4.89620E-03, 1.44594E-02, 9.91215E-03,-1.00616E-02,
     * -8.21324E-03,-1.57757E+02, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  1.53569E-02, 0.00000E+00, 6.63564E-03, 4.58410E+01, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00,-2.51280E-02, 0.00000E+00, 0.00000E+00/
      DATA PS2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-8.73148E-04,
     * -1.29648E-03,-7.32026E-05, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-4.68110E-03,
     * -4.66003E-03,-1.31567E-03,-7.39390E-04, 6.32499E-04,-4.65588E-04,
     * -1.29785E+00,-1.57139E+02, 0.00000E+00, 2.58350E-01,-3.69453E+01,
     *  4.10672E-01, 9.78196E+00,-1.52064E-01,-3.85084E+03, 0.00000E+00,
     * -8.52706E-04,-1.40945E-03,-7.26786E-04, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C          TN3(3)
      DATA PU1/
     *  9.60722E-01, 7.03757E-02,-3.00266E-02, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 2.22671E-02,
     *  4.10423E-02, 0.00000E+00, 0.00000E+00,-1.63070E+02, 0.00000E+00,
     *  5.40747E-04, 7.79481E-03, 1.44908E+02, 1.51484E-04, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00,-1.41844E-02, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 5.77884E-03, 1.06073E-02, 5.36685E-03, 9.74319E-03,
     *  0.00000E+00,-2.88015E+03, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  1.97547E-02, 0.00000E+00,-4.44902E-03,-2.92760E+01, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 2.34419E-02, 0.00000E+00, 0.00000E+00/
      DATA PU2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-4.65325E-04,
     * -5.50628E-04, 3.31465E-04, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-2.06179E-03,
     * -3.08575E-03,-7.93589E-04,-1.08629E-04, 5.95511E-04,-9.05050E-04,
     *  1.18997E+00, 4.15924E+01, 0.00000E+00,-4.72064E-01,-9.47150E+02,
     *  3.98723E-01, 1.98304E+01, 0.00000E+00, 3.73219E+03, 0.00000E+00,
     * -1.50040E-03,-1.14933E-03,-1.56769E-04, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C          TN3(4)
      DATA PV1/
     *  1.03123E+00,-7.05124E-02, 8.71615E-03, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-3.82621E-02,
     * -9.80975E-03, 0.00000E+00, 0.00000E+00, 2.89286E+01, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 8.66153E+01, 7.91938E-04, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 4.68917E-03, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 7.86638E-03, 9.57341E-03, 5.72268E-03, 9.90827E-03,
     *  0.00000E+00, 6.55573E+01, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00,-4.00200E+01, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 7.07457E-03, 0.00000E+00, 0.00000E+00/
      DATA PV2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-2.04970E-04,
     *  1.21560E-03,-8.05579E-06, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-2.49941E-03,
     * -4.57256E-04,-1.59311E-04, 2.96481E-04,-1.77318E-03,-6.37918E-04,
     *  1.02395E+00, 1.28172E+01, 0.00000E+00, 1.49903E-01,-2.63818E+01,
     *  0.00000E+00, 4.70628E+01,-2.22139E-01, 4.82292E-02, 0.00000E+00,
     * -8.67075E-04,-5.86479E-04, 5.32462E-04, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C          TN3(5) SURFACE TEMP TSL
      DATA PW1/
     *  1.00828E+00,-9.10404E-02,-2.26549E-02, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-2.32420E-02,
     * -9.08925E-03, 0.00000E+00, 0.00000E+00, 3.36105E+01, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00,-1.24957E+01,-5.87939E-03, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 2.79765E+01, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 2.01237E+03, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00,-1.75553E-02, 0.00000E+00, 0.00000E+00/
      DATA PW2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 3.29699E-03,
     *  1.26659E-03, 2.68402E-04, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 1.17894E-03,
     *  1.48746E-03, 1.06478E-04, 1.34743E-04,-2.20939E-03,-6.23523E-04,
     *  6.36539E-01, 1.13621E+01, 0.00000E+00,-3.93777E-01, 2.38687E+03,
     *  0.00000E+00, 6.61865E+02,-1.21434E-01, 9.27608E+00, 0.00000E+00,
     *  1.68478E-04, 1.24892E-03, 1.71345E-03, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C          TGN3(2) SURFACE GRAD TSLG
      DATA PX1/
     *  1.57293E+00,-6.78400E-01, 6.47500E-01, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-7.62974E-02,
     * -3.60423E-01, 0.00000E+00, 0.00000E+00, 1.28358E+02, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 4.68038E+01, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-1.67898E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 2.90994E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 3.15706E+01, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PX2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C          TGN2(1) TGN1(2)
      DATA PY1/
     *  8.66492E-01, 3.55807E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-1.12111E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 1.82458E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 1.01024E+02, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 6.54251E+02, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PY2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-1.56959E-02,
     *  1.91001E-02, 3.15971E-02, 1.00982E-02,-6.71565E-03, 2.57693E-03,
     *  1.38692E+00, 2.82132E-01, 0.00000E+00, 0.00000E+00, 3.81511E+02,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C          TGN3(1) TGN2(2)
      DATA PZ1/
     *  1.06029E+00,-5.25231E-02, 3.73034E-01, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 3.31072E-02,
     * -3.88409E-01, 0.00000E+00, 0.00000E+00,-1.65295E+02, 0.00000E+00,
     * -4.38916E-02,-3.22716E-01,-8.82393E+01, 1.18458E-01, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-1.19782E-01,-2.13801E-01, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 2.62229E+01, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     * -4.35863E-01, 0.00000E+00, 0.00000E+00,-5.37443E+01, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00,-4.55788E-01, 0.00000E+00, 0.00000E+00/
      DATA PZ2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 3.84009E-02,
     *  3.96733E-02, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 5.05494E-02,
     *  7.39617E-02, 1.92200E-02,-8.46151E-03,-1.34244E-02, 1.96338E-02,
     *  1.50421E+00, 1.88368E+01, 0.00000E+00, 0.00000E+00,-5.13114E+01,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  5.11923E-02, 3.61225E-02, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C         MIDDLE ATMOSPHERE AVERAGES
      DATA PAVGM/
     M  2.61000E+02, 2.64000E+02, 2.29000E+02, 2.17000E+02, 2.17000E+02,
     M  2.23000E+02, 2.86760E+02,-2.93940E+00, 2.50000E+00, 0.00000E+00/
      END
C-----------------------------------------------------------------------
