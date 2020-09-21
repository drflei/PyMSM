      PROGRAM MSM
c........+.........+.........+.........+.........+.........+.........+..
c
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-B)
      IMPLICIT REAL*8 (D-H)
      IMPLICIT REAL*8 (O-Z)

      INCLUDE 'interface.h'

      CHARACTER*4  CYEAR
c
c........+.........+.........+.........+.........+.........+.........+..

      COMMON /BANDL/   BBPT, FLPT, GLMDAR
      COMMON /CUTRIG/  RLGVPV,RUGVPV,RCGVPV, RLGVPE,RUGVPE,RCGVPE, 
     &                 RLGVPW,RUGVPW,RCGVPW, PTV
      COMMON /KYMDHM/  KYEAR,KMONTH,KDATE,KHOUR,KMIN,KSEC,KDOY,KUT

      REAL*8   BBPT1, FLPT1, GLMDAR1 
      REAL*8   RLGVPV1,RUGVPV1,RCGVPV1, RLGVPE1,RUGVPE1,RCGVPE1, 
     &                 RLGVPW1,RUGVPW1,RCGVPW1, PTV1
      REAL*8   BBPT2, FLPT2, GLMDAR2 
      REAL*8   RLGVPV2,RUGVPV2,RCGVPV2, RLGVPE2,RUGVPE2,RCGVPE2, 
     &                 RLGVPW2,RUGVPW2,RCGVPW2, PTV2

      REAL*8   TMSTART, TMEND, TMDUR, OSEC, TSTART, TEND
      INTEGER*4 OYEAR, OMON, ODAY, OHR, OMN, NTRAJ
      INTEGER*4 DATEVALS(8), KSKIP /0/, HEADER(2000), SAOLUN /3/
      CHARACTER*80 PRODEF, TITLE, ORBTIT, FOOTER
      CHARACTER*3 MONTH(12) /'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
     &                       'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'/

      CHARACTER*1000 filein, fileout, arg, dir
c
      INTEGER NRC
      PARAMETER (NRC=7)
      REAL*8 RC(NRC), F(NRC)
c
C Orbit header common block
      COMMON /OHEAD/ HEADER
      COMMON /MSMHOME/ dir

c Name list definition
      namelist/msmlist/filein, fileout, dir, fkpidx
c
      TYPE(META_VAR) :: MISVAR, ORBVAR

c........+.........+.........+.........+.........+.........+.........+..
c     \/ The following data statement initializes control !DS 2 Jun 2002
c            need for the first call to SETKP3UT          !DS 2 Nov 2004
c                  these values are modified in SETKP3UT  !DS 2 Nov 2004
c........+.........+.........+.........+.........+.........+.........+..

      DATA IYROLD,IMOLD,IHROLD,M5OLD,KPDXOLD,CYEAR 
     &        / 0,    0,    25,   -1,     11, '1995'/
      DATA RC/0.1, 0.15, 0.213, 0.315, 0.35, 0.5, 1.0/
c........+.........+.........+.........+.........+.........+.........+..
c     UNIT 1  World grid cutoff data input; used in subroutine SETKP3UT
c     UNIT 3  Input file of locations used for interpolation
c     UNIT 7  File of interpolated cutoff rigidities for locations on UNIT 3
c     UNIT 8  File of interpolated proton cutoff energies 
c                  for locations on UNIT 3
c     UNIT 9  File of interpolated proton cutoff energies , "L" & LAMBDA
c                  for locations on UNIT 3
c
c    take inputs from command line if supplied, otherwise from namelist file msm.nml
c
c........+.........+.........+.........+.........+.........+.........+..
      integer :: count
      real*8 :: xin(3,1),xout(3,1)
      character*256 :: SSYEAR, SSTIME

      INCLUDE 'putils.version'

      CALL DATE_AND_TIME(VALUES=DATEVALS)
      WRITE (SDATE,'(I2,''-'',A,''-'',I4)')
     &  DATEVALS(3), MONTH(DATEVALS(2)), DATEVALS(1)
      WRITE (STIME,'(2(I2.2,'':''),I2.2)') (DATEVALS(I),I=5,7)

      count = command_argument_count()
      if (count .eq. 4) then
         CALL get_command_argument(1, arg)
         filein = TRIM(ADJUSTL(arg))
         CALL get_command_argument(2, arg)
         read(arg,*) fkpidx
         CALL get_command_argument(3, arg)
         fileout = TRIM(ADJUSTL(arg))
         CALL get_command_argument(4, arg)
         dir = TRIM(ADJUSTL(arg))
      elseif (count .eq. 0) then
c  read in the name list definition
c
         nmlstat = 0
         open(unit=9,file='msm.nml',status='old',err=5)
         read(9,msmlist,err=5)
         close(9)
         nmlstat = 1
 5       if (nmlstat .eq. 0) then
            print *, ' Namelist file - msm.nml - error'
            call exit
         end if
      else
        print*, 'To Run: $ ./msm '
        print*, '   or   $ ./msm orbit-file kp output-trans-file MSMHOME
     &_DIR'
        call exit()
      end if
c........+.........+.........+.........+.........+.........+.........+..
c     \/ Load in default kp0 t00 world grid (UNIT1)
c........+.........+.........+.........+.........+.........+.........+..
      KPIDX = 0
      IUT   = 0
      CALL SETKP3UT (CYEAR,KPIDX,KPDXOLD, IUT,IUTOLD, FKPIDX)                  ! SUB02
C
      KPIDX = NINT(fkpidx)
c
      OPEN (SAOLUN,file=filein, status='old')
      read(saolun,*)
C Read orbit header
c      CALL READ_ORB(PRODEF, TITLE, NTRAJ, TMSTART, TMEND, TMDUR, ORBTIT,
c     &              OYEAR, OMON, ODAY, OHR, OMN, OSEC, TSTART, TEND,
c     &              SAOLUN)
c     CALL MV_READ('ORB', ORBVAR)
c      CALL MV_READ('MIS', MISVAR)

      OPEN (UNIT=9,FILE=fileout, BLANK='ZERO',STATUS='UNKNOWN')    ! AFRL
c      CALL MF_ORB(9, PRODEF, TITLE, MISVAR, ORBVAR, NRC, RC)
c
c........+.........+.........+.........+.........+.........+.........+..
c    !!! Beginning of logical loop, will continue until end of data
c     \/ Read data point
c        Get year month date hour of data point
c            A negative IYR will terminate input
c........+.........+.........+.........+.........+.........+.........+..

c      \/ Set date and time
c........+.........+.........+.........+.........+.........+.........+..

      KMONTH = 1
      KDATE  = 1
      KHOUR  = 0
      KMIN   = 0
c........+.........+.........+.........+.........+.........+.........+..
c     \/ Set up world grid cutoff data interpolation files
c        change setup every time Kp changes
c         or when UT hour (KHOUR or GH) advances into next 3-hr interval
c            iut =  0  !!!  22:30 - 00 - 01:30
c            iut =  3  !!!  01:30 - 03 - 04:30
c            iut =  6  !!!  04:30 - 06 - 07:30
c            iut =  9  !!!  07:30 - 09 - 10:30
c            iut = 12  !!!  10:30 - 12 - 13:30
c            iut = 15  !!!  13:30 - 15 - 1630
c            iut = 18  !!!  16:30 - 18 - 19:30
c            iut = 21  !!!  19:30 - 21 - 22:30
c........+.........+.........+.........+.........+.........+.........+..
      ntimes =1
c
c  100 CALL READ_EPH(FMJD, ALTKM, FLOND, FLATD, SPS, FOOTER, KSKIP)
c      IF (FOOTER .NE. '*CONTINUE*') GO TO 900
  100 READ (SAOLUN,*,IOSTAT=IOSTAT, ERR=110, END=900)
     &          SSYEAR, SSTIME, X,Y,Z
      READ(SSYEAR,'(I4,x,I2,x,I2)')KYEAR,KMONTH,KDATE
      READ(SSTIME,'(I2,x,I2,x,F6.3)')KHOUR,KMIN,SEC
c........+.........+.........+.........+.........+.........+.........+..
      KSEC = INT(SEC)
      print *, KYEAR, KMONTH, KDATE, KHOUR,KMIN,KSEC

c      Dec_y = KYEAR + (JULDAY(KYEAR,KMONTH,KDATE) - JULDAY(KYEAR,1,1)
c     &     + FMJD - DINT(FMJD) - 0.5D0)/365.D0  ! JULDAY starts at noon
c      KDOY = DINT(Dec_y)
      KDOY = GET_DOY(KYear,KMonth, Kdate)

c      CALL DECY2DATE_AND_TIME(Dec_y, KY, KM, KD, kdoy, KHOUR,KMIN,ks,UT)
      
      KHP = 0
      IF (KMIN.LE.30)  THEN
         KHP = 0
       ELSE
         KHP = 1
      ENDIF
      IUT = ((KHOUR + 1 +KHP)/3)*3
      IF (IUT.GE.24)  IUT = 0
      KUT = IUT*3600  ! this is used for L calculation, rather than the actual ut
      UT = KUT*1.D0
c
c      PRINT *,KY,KM, KD, KHOUR, KMIN, ks, ut, KPIDX, IUT
c
      xin(1,1) = x/6371.2
      xin(2,1) = y/6371.2
      xin(3,1) = z/6371.2
      utime = khour*3600. + kmin*60. + sec
      call coord_trans_vec1(ntimes,5,0,kyear,kdoy,utime,xIN,xOUT)

      ALTKM = xout(1,1)
      FLATD = xout(2,1)
      FLOND = xout(3,1)
      print *,xin
      print *,xout

c........+.........+.........+.........+.........+.........+.........+..
c        Remember, cutoff are calculated for > 20 KM,
c                  IF ALT <20, use 20 KM
c........+.........+.........+.........+.........+.........+.........+..
c
      IF (ALTKM.LT.20.0)  ALTKM = 20.0
  110 CONTINUE

c........+.........+.........+.........+.........+.........+.........+..

      IF (FLATD.GT. 89.98)  FLATD =  89.98
      IF (FLATD.LT.-89.98)  FLATD = -89.98
      IF (FLOND.LT.0.0)     FLOND = FLOND + 360.0
      IF (FLOND.GT.360.0)   FLOND = FLOND - 360.0
c
C
      IF (KYEAR .GE. 2015) THEN
         LYEAR = 2015
      ELSE IF (KYEAR .GE. 2010) THEN
         LYEAR = 2010
      ELSE IF (KYEAR .GE. 2005) THEN
         LYEAR = 2005
      ELSE IF (KYEAR .GE. 2000) THEN
         LYEAR = 2000
      ELSE IF (KYEAR .GE. 1995) THEN
         LYEAR = 1995
      ELSE IF (KYEAR .GE. 1990) THEN
         LYEAR = 1990
      ELSE IF (KYEAR .GE. 1985) THEN
         LYEAR = 1985
      ELSE IF (KYEAR .GE. 1980) THEN
         LYEAR = 1980
      ELSE IF (KYEAR .GE. 1975) THEN
         LYEAR = 1975
      ELSE IF (KYEAR .GE. 1970) THEN
         LYEAR = 1970
      ELSE IF (KYEAR .GE. 1965) THEN
         LYEAR = 1965
      ELSE IF (KYEAR .GE. 1960) THEN
         LYEAR = 1960
      ELSE 
         LYEAR = 1955
      END IF
      WRITE(CYEAR,'(I4)') LYEAR
c........+.........+.........+.........+.........+.........+.........+..
c     \/ Call vertical cutoff interpolation subroutine
c........+.........+.........+.........+.........+.........+.........+..
      
      CALL CUTUT3 (CYEAR,KPIDX,KPDXOLD,IUT,IUTOLD,
     & FLOND,FLATD,ALTKM,FKPIDX)          ! SUB03

c........+.........+.........+.........+.........+.........+.........+..
c     \/ Call east-west cutoff interpolation subroutine (arguments in common)
c        Get invariant latitude degrees
c........+.........+.........+.........+.........+.........+.........+..

      CALL CUTEW (FLOND,FLATD,ALTKM) ! MOD 27 AUG 02

C NOW DECIDE IF SEND INTERPOLATION IS REQUIRED
c
      IF ((LYEAR .GT. 1955) .AND. (LYEAR .LE. 2015)) THEN
         YEAR = KYEAR + 
     &        (KDOY+KUT/86400.D0)/365.D0
         DFY = YEAR - LYEAR
C         PRINT *, 'dy = ',DFY
         IF (DFY .GT. 0.1) THEN
            BBPT1 = BBPT; FLPT1 = FLPT ; GLMDAR1 = GLMDAR 
c            RLGVPV1 = RLGVPV; RUGVPV1 = RUGVPV; 
            RCGVPV1 = RCGVPV 
c            RLGVPE1 = RLGVPE; RUGVPE1 = RUGVPE; RCGVPE1 = RCGVPE 
c            RLGVPW1 = RLGVPW; RUGVPW1 = RUGVPW; RCGVPW1 = RCGVPW 
c            PTV1 = PTV
C
            LYEAR = LYEAR + 5
            IF (LYEAR .GT. 2015) LYEAR = 2015
            F1 = DFY/5.D0 ; F2 = 1. - F1
            WRITE(CYEAR,'(I4)') LYEAR
C            PRINT *,CYEAR, YEAR, F1, F2

            CALL CUTUT3(CYEAR,KPIDX,KPDXOLD,IUT,
     &           IUTOLD,FLOND,FLATD,ALTKM,FPKIDX)
            CALL CUTEW (FLOND,FLATD,ALTKM)
C 
            BBPT = F1*BBPT1+F2*BBPT; FLPT = F1*FLPT1 + F2*FLPT
            GLMDAR = F1*GLMDAR1 + F2*GLMDAR 
c            RLGVPV = F1*RLGVPV1 + F2*RLGVPV
c            RUGVPV = F1*RUGVPV1 + F2*RUGVPV
            RCGVPV = F1*RCGVPV1 + F2*RCGVPV 
c            RLGVPE = F1*RLGVPE1 + F2*RLGVPE 
c            RUGVPE = F1*RUGVPE1 + F2*RUGVPE 
c            RCGVPE = F1*RCGVPE1 + F2*RCGVPE 
c            RLGVPW = F1*RLGVPW1 + F2*RLGVPW 
c            RUGVPW = F1*RUGVPW1 + F2*RUGVPW 
c            RCGVPW = F1*RCGVPW1 + F2*RCGVPW 
c            PTV = F1*PTV1 + F2*PTV
         END IF
      END IF

      FLMDAD = GLMDAR*57.2957795

C     NOW CALCULATE THE TRANSMISSION FUNCTION
C
      CALL EXPFAC(RCGVPV, GLMDAR, RC, F, NRC)
   
c      WRITE (9,*) FLATD,FLOND,FLPT,RCGVPV,FACSHADOW(1.D0+ALTKM/6371.D0),
c      WRITE (9,*) FLATD,FLOND,FLPT,RCGVPV,EPNRCV
      WRITE(9, '(F8.4,1P,100('','',E13.4:),0P)') 
     &  FKPIDX, FLPT, FLMDAD, RCGVPV, FACSHADOW(1.D0+ALTKM/6371.D0), F

      GO TO 100

c........+.........+.........+.........+.........+.........+.........+..
c    /\ End of logical loop, end of data
c........+.........+.........+.........+.........+.........+.........+..

  900 WRITE(9, '(A)') '''End of File'''
      CLOSE(9, STATUS='KEEP')
      CLOSE(SAOLUN)

      STOP
      END

!**********************************************************************!
!     SUBROUTINE SETKP3UT                                      ! SUB02 !
!     LOAD IN 450 KM 5X5 WORLD GRID BY KP AND UT                       !
!**********************************************************************!

      SUBROUTINE SETKP3UT (CYEAR,KPIDX,KPDXOLD,IUT,IUTOLD,FKPIDX)
C........+.........+.........+.........+.........+.........+.........+..
C     Subroutine to load in proper 450 km 5x5 world grid in CYEAR and
c                Set up "L" interpolation
C     This is designed to be a operational version
C     Diagnostic outputs must be explicitly turned on to obtain tables
C........+.........+.........+.........+.........+.........+.........+..
c     Note: - The programming adheres to the conventional FORTRAN
c             default standard that variables beginning with
c             'I','J','K','L','M',or 'N' are integer variables
c             Variables beginning with "C" are character variables
C             Exception, the variable  str  is a character variable     ! AFRL
C             All other variables are REAL*8
C........+.........+.........+.........+.........+.........+.........+..
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-B)
      IMPLICIT REAL*8 (D-H)
      IMPLICIT REAL*8 (O-Z)
C
C........+.........+.........+.........+.........+.........+.........+..
C
      COMMON /SETRCL/ VKRC(37,75),  RCG(37,75), FLG(37,75)
c     *                ,VKRU(37,75), VKRL(37,75),VKRC(37,75), PWV(37,75)
      COMMON /DIAGC/  LTS(37),     LNS(75)
C
      COMMON /KYMDHM/  KYEAR,KMONTH,KDATE,KHOUR,KMIN,KSEC,KDOY,KUT
C........+.........+.........+.........+.........+.........+.........+..
C
c     Changes required for the use of the IRBEM library.  F.Lei
c     https://svn.code.sf.net/p/irbem/code/tags/IRBEM-4.4.0/manual/user_guide.html
c
      INTEGER*4    options(5)
      real*8       maginput(25,1), lm(1,1)
C
c
      CHARACTER*1  CKP
      CHARACTER*2  CUT

      CHARACTER*1000 dir                                                 ! AFRL
      CHARACTER*1000 str                                                 ! AFRL
      CHARACTER*4 CYEAR, CYEAROLD
      CHARACTER*21 CID
c
      COMMON /MSMHOME/ dir

      data options/0,30,0,0,0/

      SAVE                                                              ! AFRL
C
      IF ((KPIDX.EQ.KPDXOLD).AND.(IUT.EQ.IUTOLD) 
     &       .AND.(CYEAR.EQ.CYEAROLD)) RETURN
C
      CLOSE(1, STATUS='KEEP')
C
C........+.........+.........+.........+.........+.........+.........+..
C........+.........+.........+.........+.........+.........+.........+..
C     \/ Load data input files with -1.0
C        Used for checking for missing data
C........+.........+.........+.........+.........+.........+.........+..
C
      DO 90 J=1,37
         DO 80 K = 1,75
            FLG(J,K) = -1.0
            RCG(J,K) = -1.0
c            RLG(J,K) = -1.0
c            RUG(J,K) = -1.0
c            PWV(J,K) = -1.0
   80   CONTINUE
   90   CONTINUE
C
C........+.........+.........+.........+.........+.........+.........+..
C
C     \/ Setup, (1) Read in World Grid of cutoff rigidities
C        KPIDX Specifies Kp of world grid data file to be used
C........+.........+.........+.........+.........+.........+.........+..
  120 CONTINUE

c      call GETENV( "MSMHOME", dir )
      IF (KPIDX.LE.9) THEN
         WRITE(CKP,'(I1)')KPIDX
      ELSE
         CKP = 'X'
      ENDIF
      WRITE(CUT,'(I2.2)')IUT
      str = TRIM(ADJUSTL(dir))//'/MAPS/'//CYEAR//'/AVKP'//CKP//'T'//
     &      CUT//'.AVG'
c      write(*,*) TRIM(str)
c
      OPEN (1,FILE=str,BLANK='ZERO',STATUS='OLD')
C
C........+.........+.........+.........+.........+.........+.........+..

      PLAT = 0.0
  150 CONTINUE
      READ (1, 2010, IOSTAT=IOSTAT,ERR=150,END=170)
     *              RLATD,RLOND,IALT,RDFL,SSRU,SSRCK,SSRL
      PENW = SSRU - SSRL
c      WRITE (*, *) RLATD,RLOND,SSRU,SSRL,SSRCK  
  160 CONTINUE
 2010 FORMAT (1X, F7.2, F7.2,I4,F9.2,F6.2,F6.2,F6.2)
C........+.........+.........+.........+.........+.........+.........+..
C     \/ Get index of location LAT & LON data on UNIT 1
C        RLATD is Latitude  in Degrees
C        RLOND is Longitude in Degrees (East of Greenwich)
C        J is Latitude index
C        K is Longitude index
C          K(1)  =  -5 Deg Long,   is not used
C          K(2)  =   0 Deg Long
C          K(74) = 360 Deg Long
C          K(75)   365 Deg Long  =  5 Deg Long
C        +.........+.........+.........+.........+.........+.........+..
C        \/ Find Latitude index
C        +.........+.........+.........+.........+.........+.........+..
C
         J = (90.0-RLATD)/5.0+1.0001
         IF (RLATD.GE. 89.8)  J = 1
         IF (RLATD.LE.-89.8)  J = 37
C
C        +.........+.........+.........+.........+.........+.........+..
C        \/ Find Longitude index
C        +.........+.........+.........+.........+.........+.........+..
C
         IF (RLOND.LT.0.0)  RLOND = 360.0+RLOND
         K = RLOND/5.+2.0001
C
         FLG(J,K) = RDFL
         RCG(J,K) = SSRCK
c         RLG(J,K) = SSRL
c         RUG(J,K) = SSRU
c         PWV(J,K) = PENW
C
      IF (IOSTAT.EQ.0)  GO TO 150
C
  170 CONTINUE
C
      KPDXOLD = KPIDX
      IUTOLD  = IUT
      CYEAROLD = CYEAR
c
C
      KAL = 0
C
      ALT450 = 450.D0
      UT = KUT*1.D0
C
      PLAT = 0.0
      DO 210 JCK=1,37
         DO 200 KCK = 2,73
            DLATD = 90-((JCK-1)*5)
            IF (DLATD.GT. 89.98)  DLATD =  89.98
            IF (DLATD.LT.-89.98)  DLATD = -89.98
            DLOND = (KCK-2)*5.0
            IF (DLOND.LT.0.0)  DLOND = DLOND + 360.0
            IF (FLG(JCK,KCK) .LE. 0.0)  THEN
C
c               CALL INVARA (MODEL,TM,DLATD,DLOND,ALT450,ERR,BB,FL) 
c	F.Lei  - Removed the INVARA sub and changed to use equivalent ones in the IRBEM lib               
c                call make_lstar1(ntime,kext,options,sysaxes,iyear,idoy,ut, x1,x2,x3, maginput, lm,lstar,blocal,bmin,xj,mlt)
c
               maginput(1,1) = FKPIDX*10.
               if (FKPIDX.GE.9.) maginput(1,1) = 89.999D0
c      The IRBEM accepts KP 0-9 only, and YEAR 1965 - 2015
               READ(CYEAR,'(I4)')LYEAR
               IF (LYEAR .GE. 2015) LYEAR = 2014
               IF (LYEAR .LT. 1965) LYEAR = 1965
c      remember this is for the map grid and currently this already be precalculated!
	           call make_lstar1(1,4,options,0,LYEAR,1,UT,ALT450,
     *              DLATD,DLOND,maginput, lm,lstar,blocal,bmin,xj,mlt)
	       FL = DABS(lm(1,1))
c	end of F.Lei's changes
c
               IF (FL.GT.999.99)  FL = 999.99
               FLG(JCK,KCK) = FL
               KAL = KAL + 1
C
            ENDIF
  200    CONTINUE
  210 CONTINUE

C........+.........+.........+.........+.........+.........+.........+..
C     \/ Continue setup, fill in edges (-5, 365)
C        input data file now covers 0 to 360
C........+.........+.........+.........+.........+.........+.........+..
C
      DO 350 LT = 1, 37
         RCG(LT,1) = RCG(LT,73) 
c         RLG(LT,1) = RLG(LT,73) 
c         RUG(LT,1) = RUG(LT,73) 
         FLG(LT,1) = FLG(LT,73) 
c         PWV(LT,1) = PWV(LT,73) 

         RCG(LT,75) = RCG(LT,3) 
c         RLG(LT,75) = RLG(LT,3) 
c         RUG(LT,75) = RUG(LT,3) 
         FLG(LT,75) = FLG(LT,3) 
c         PWV(LT,75) = PWV(LT,3) 
  350 CONTINUE
C
C........+.........+.........+.........+.........+.........+.........+..
C     \/ Check for missing data
C        DLATD is Latitude  in degrees
C        DLOND is Longitude in degrees (East of Greenwich)
C........+.........+.........+.........+.........+.........+.........+..
C
      IFCNT = 0
c      KMIS = 0
c      DO 370 J=1,37
c         DO 360 K = 1,75
c            DLATD = 90-((J-1)*5)
c            IF (DLATD.GT. 89.98)  DLATD =  89.98
c            IF (DLATD.LT.-89.98)  DLATD = -89.98
c            DLOND = (K-2)*5.0
c            IF (RCG(J,K).LT.0)  THEN
c               KMIS = KMIS + 1
c               WRITE (20, 2210)  J, K, DLATD, DLOND
c            ENDIF
c  360    CONTINUE
c  370 CONTINUE
cC
c      IF (KMIS.GT.0)  THEN
c         WRITE (*, 2220) KMIS
c         STOP 2220
c      ENDIF
c 2210 FORMAT (' 2210 Missing at ',2I4, 2F7.2)
c 2220 FORMAT (' 2220 Fatal setup error, File ', a46/ 
c     &         I5,' missing input locations, PROGRAM STOPPED.'/
c     &        '  See  MISSING.DAT  file!!!')
C
C........+.........+.........+.........+.........+.........+.........+..
C     \/ Replace zero cutoff values with 1.0-6 to avoid possible "blowup" 
C........+.........+.........+.........+.........+.........+.........+..
C
      DO 410 J = 1, 37
         DO 400 K = 1, 75
            IF (RCG(J,K).EQ.0.0)  RCG(J,K) = 1.0E-06
c            IF (RLG(J,K).EQ.0.0)  RLG(J,K) = 1.0E-06
c            IF (RLG(J,K).EQ.0.0)  RLG(J,K) = 1.0E-06
  400    CONTINUE
  410 CONTINUE

C........+.........+.........+.........+.........+.........+.........+..
C     \/ Set up VKRC, VKRL and VKRU for each grid point
C        Assume form of  R = VK*L**2
C........+.........+.........+.........+.........+.........+.........+..
C
      DO 430 J = 1, 37
         DO 420 K = 1, 75
            VKRC(J,K) = RCG(J,K)*FLG(J,K)**2
c            VKRL(J,K) = RLG(J,K)*FLG(J,K)**2
c            VKRU(J,K) = RUG(J,K)*FLG(J,K)**2
  420    CONTINUE
  430 CONTINUE
C
C........+.........+.........+.........+.........+.........+.........+..
C     \/ Diagnostic printout of cutoff is turned off in operational version
C     ...+.........+.........+.........+.........+.........+.........+..
C     \/ Diagnostic printout of cutoff values RU, RC, RL
C        (1) Write table in landscape mode in the following sequence
C        (2) Write out line of RU cutoff values
C        (3) Write out line of RL cutoff values
C        (4) Write out line of RC cutoff values
C        (5) Write out line of L values
C        (6) Write out line of VKRU values
C        (7) Write out line of VKRL values
C        (8) Write out line of VKRC values
C     ...+.........+.........+.........+.........+.........+.........+..
C     Generate table, Table orientation is landscape mode, 15 CPI
C     \/ Generate longitude index for printing table 
C     ...+.........+.........+.........+.........+.........+.........+..
C
      DO 440 I = 1, 75
         IFCNT = IFCNT+1
         LNS(IFCNT) = (IFCNT-2)*5
  440 CONTINUE
C
C     ...+.........+.........+.........+.........+.........+.........+..
C     \/ Generate latitudes for table
C     ...+.........+.........+.........+.........+.........+.........+..
C
      DO 450 LT = 1, 37
         LTS(LT) = 90-((LT-1)*5)
  450 CONTINUE
C
      RETURN
      END

!**********************************************************************!
!     SUBROUTINE CUTUT3                                       ! SUB03 !
!                Control routine                                       !
!**********************************************************************!

      SUBROUTINE CUTUT3 (CYEAR,KPIDX,KPDXOLD,IUT,IUTOLD,
     &                   FLOND,FLATD,ALTKM,FPKIDX)
c........+.........+.........+.........+.........+.........+.........+..
c     Control routine for interpolating vertical cutoff rigidities
c     Version 4.03; 3-HR UT Version; Uses 3-hr avg RU, RC, RL & PW for all KP
c             INCOMPATIBLE WITH LOWER LEVEL VERSIONS 
c........+.........+.........+.........+.........+.........+.........+
c     Remember: particle cutoff rigidity in GV   IN  COMMON /CUTRIG/
c     Remember: proton   cutoff energy   in MEV  IN  COMMON /CUTERGY/
c     Remember: "B & L" related parameters       IN  COMMON /BANDL/
c........+.........+.........+.........+.........+.........+.........+..
cLast Mod 28 Oct 2004  Checking, Correction to CUTRIG common list
c     Mod 02 Jun 2003  Mod for Penumbral transparency vs Kp
c     Mod 20 Feb 2003  Identify AFRL GEOSPACE additions for SGI compatibility
c     Mod 25 Jun 2002  Cutoff values in COMMON CUTRIG AND CUTERGY 
c     Mod 25 May 2002  Current YR, MO, DT, HR, KMIN in COMMON /KYMDHM/
c     Mod 28 Nov 2001  Take care of 24 HR IUT = 0
c     Mod 23 Nov 2001  KHOUR & KMIN used to converted to (3 hour UT segments)
c     Mod 09 Oct 2000  REAL*8 for compatibility with NASA JSC SRAG
c........+.........+.........+.........+.........+.........+.........+..
c     Programmed by Don F. Smart  (SSSRC@MSN.COM)
c     Note: - The programming adheres to the conventional FORTRAN
c             default standard that variables beginning with
c             'I','J','K','L','M',or 'N' are integer variables
c             Variables beginning with "C" are character variables
C             Exception, the variable  str  is a character variable     ! AFRL
C             All other variables are REAL*8
c........+.........+.........+.........+.........+.........+.........+..
c
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-B)
      IMPLICIT REAL*8 (D-H)
      IMPLICIT REAL*8 (O-Z)
C
      CHARACTER*4 CYEAR, CYEAROLD
c
c........+.........+.........+.........+.........+.........+.........+..

      COMMON /BANDL/   BBPT, FLPT, GLMDAR
      COMMON /CUTRIG/  RLGVPV,RUGVPV,RCGVPV, RLGVPE,RUGVPE,RCGVPE, 
     &                 RLGVPW,RUGVPW,RCGVPW, PTV
c      COMMON /CUTERGY/ EPNRLV,EPNRUV,EPNRCV, EPNRLE,EPNRUE,EPNRCE,
c     &                 EPNRLW,EPNRUW,EPNRCW 
      COMMON /KYMDHM/  KYEAR,KMONTH,KDATE,KHOUR,KMIN,KSEC,KDOY,KUT
c      COMMON /PTV/     PTVDAT(40)

      SAVE                                                              ! AFRL

c........+.........+.........+.........+.........+.........+.........+..
c     \/ Entry values for diagnostic
c........+.........+.........+.........+.........+.........+.........+..

c     WRITE ( *,3010)  KPIDX,KPDXOLD,KHOUR,KMIN,IUTOLD,            !DIAG
c    &                 FLATD,FLOND,ALTKM                           !DIAG
C     WRITE (66,3010)  KPIDX,KPDXOLD,KHOUR,KMIN,IUTOLD,            !DIAG
C    *                 FLATD,FLOND,ALTKM                           !DIAG
C3010 FORMAT (' 3010 ',5I4,3F7.2)                                  !DIAG
C     WRITE (66, 3012)  KHOUR, IUT
C3012 FORMAT (' 3012 KHOUR, IUT',2I5)

      IF (KPIDX.NE.KPDXOLD .OR. IUT.NE.IUTOLD 
     &     .OR. CYEAR.NE.CYEAROLD)  THEN
C         WRITE ( *, 3020)  KPIDX,KPDXOLD, IUT,IUTOLD         ! DIAG
C         WRITE (66, 3020)  KPIDX,KPDXOLD, IUT,IUTOLD         ! DIAG

          CALL SETKP3UT (CYEAR,KPIDX,KPDXOLD, IUT,IUTOLD,FKPIDX)           ! SUB02
          KPDXOLD = KPIDX
          IUTOLD  = IUT
          CYEAROLD = CYEAR

C         WRITE ( *,3030)  KPIDX,KPDXOLD, IUT,IUTOLD          ! DIAG
C         WRITE (66,3030)  KPIDX,KPDXOLD, IUT,IUTOLD          ! DIAG
      ENDIF
C3020 FORMAT (' 3020 CALL  SETKP3UT '4I4)                     ! DIAG
C3030 FORMAT (' 3030 RET   SETKP3UT '4I4)                     ! DIAG
c
c........+.........+.........+.........+.........+.........+.........+..
c    \/ Get interpolated vertical cutoff rigidity for this position
c       Remember, cutoff values are now stored in COMMON 
c........+.........+.........+.........+.........+.........+.........+..
c
c     WRITE (66,3040)                                                  ! DIAG
c     WRITE (66,3050)  FLATD, FLOND, ALTKM                             ! DIAG

c      CALL PENTRAKP (KPIDX)                                            ! SUB07
      CALL LINT5X5 (FLATD,FLOND,ALTKM,FKPIDX)                           ! SUB04

c     WRITE (66,3060)                                                  ! DIAG
c     WRITE (66,3070)  FLATD,FLOND,ALTKM,BBPT,FLPT,                    ! DIAG
c    &                 PTV,RLGVPV,RUGVPV,RCGVPV                        ! DIAG
c
  490 CONTINUE
c3040 FORMAT (' 3040 CALL LINT5X5 ( FLATD, FLOND,    ALTKM')           ! DIAG
c3050 FORMAT (' 3050', 14X,2F7.2,F10.3)                                ! DIAG
c3060 FORMAT (' 3060 RET  LINT5X5 ( FLATD, FLOND,    ALTKM,   BBPT,',  ! DIAG
c    *        '  FLPT,   PTV, RLGVPV RUGVPV, RCGVPV)')                 ! DIAG
c3070 FORMAT (' 3070', 14X,2F7.2,F10.3,F8.5, 5F7.3)                    ! DIAG

      RETURN
      END

!**********************************************************************!
!     SUBROUTINE LINT5X5                                       ! SUB04 !
!     "L" INTERPOLATION PROGRAM (VERSION 4.03)                         !
!**********************************************************************!

      SUBROUTINE LINT5X5 (FLATD,FLOND,ALTKM,FKPIDX)
c........+.........+.........+.........+.........+.........+.........+..
c     "L" Interpolation program
c     VERSION 4.03; 3-HR UT Version; Uses 3-hr avg RU, RC, RL & PW for all KP
C             INCOMPATIBLE WITH LOWER LEVEL VERSIONS
c     Remember: particle cutoff rigidity in GV   IN  COMMON /CUTRIG/
c     Remember: proton   cutoff energy   in MEV  IN  COMMON /CUTERGY/
c     Remember: "B & L" related parameters       IN  COMMON /BANDL/
c........+.........+.........+.........+.........+.........+.........+..
cLast Mod 20 Jun 2003 @ 1600  Testing 
c     Mod 02 Jun 2003  Check IGH never .lt. 1
c     Mod 02 Jun 2003  Mod for Penumbral transparency vs Kp
c     Mod 02 Jun 2003  Add sub PENTRAKP (Penumbral transparency vs Kp)
c     Mod 20 Feb 2003  Identify AFRL GEOSPACE additions for SGI compatibility
c     Mod 28 Aug 2002  Radial distance adjust (reduce Geosynch)
c     Mod 25 JUN 2002  Position and outputs in common
c     Mod 09 OCT 2000  Check limits of FLATD and FLOND
c     Mod 09 Oct 2000  REAL*8 for compatibility with NASA JSC SRAG
c     Use CALL INVARA (MODEL,TM,FLATD,FLOND,ALT450, ERR,BBPT,FLPT) 
c........+.........+.........+.........+.........+.........+.........+..
c     Programmed by Don F. Smart  (SSSRC@MSN.COM)
c     Note: - The programming adheres to the conventional FORTRAN
c             default standard that variables beginning with
c             'I','J','K','L','M',or 'N' are integer variables
c             Variables beginning with "C" are character variables
C             Exception, the variable  str  is a character variable     ! AFRL
C             All other variables are REAL*8 
c........+.........+.........+.........+.........+.........+.........+..
c
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-B)
      IMPLICIT REAL*8 (D-H)
      IMPLICIT REAL*8 (O-Z)
      
      INTEGER*4    options(5)
      real*8       maginput(25,1), lm(1,1)
      data options/0,30,0,0,0/
c
c........+.........+.........+.........+.........+.........+.........+..
c
c     Diagnostic outputs
c        UNIT 10='T10RLBOX'
c        UNIT 11='T11RUBOX'
c        UNIT 12='T12RUBOX'
c        UNIT 14='T14VK'
c
c........+.........+.........+.........+.........+.........+.........+..
c
      COMMON /BANDL/   BBPT, FLPT, GLMDAR
      COMMON /CUTRIG/  RLGVPV,RUGVPV,RCGVPV, RLGVPE,RUGVPE,RCGVPE, 
     &                 RLGVPW,RUGVPW,RCGVPW, PTV
c      COMMON /CUTERGY/ EPNRLV,EPNRUV,EPNRCV, EPNRLE,EPNRUE,EPNRCE,
c     &                 EPNRLW,EPNRUW,EPNRCW 
      COMMON /KYMDHM/  KYEAR,KMONTH,KDATE,KHOUR,KMIN,KSEC,KDOY,KUT
c      COMMON /PTV/     PTVDAT(40)

      COMMON /SETRCL/ VKRC(37,75), RCG(37,75), FLG(37,75)
c     *                ,VKRU(37,75), VKRL(37,75),VKRC(37,75), PWV(37,75)
      COMMON /DIAGC/  LTS(37),     LNS(75)
c
c........+.........+.........+.........+.........+.........+.........+..
c
      SAVE                                                              ! AFRL

C      MODEL = 14
C      TM = 1995.0
c
C........+.........+.........+.........+.........+.........+.........+..
C     \/ Diagnostic printout of cutoff is turned off in operational version
C........+.........+.........+.........+.........+.........+.........+..
C
c      IDG10 = 0
c      IDG11 = 0
c      IDG12 = 0
c      IDG14 = 0

c      IF (IDG10.NE.0)  OPEN  (10,FILE='T10RLBOX',STATUS='UNKNOWN')
c      IF (IDG11.NE.0)  OPEN  (11,FILE='T11RUBOX',STATUS='UNKNOWN')
c      IF (IDG11.NE.0)  OPEN  (12,FILE='T12RCBOX',STATUS='UNKNOWN')
c      IF (IDG14.NE.0)  OPEN  (14,FILE='T14VK',   STATUS='UNKNOWN')
c
c........+.........+.........+.........+.........+.........+.........+..
c
c     \/ Get index of location
c        J is Latitude  index (1 is  N Pole)
c        K is Longitude index
c          K(1)  =  -5 Deg LONG, is not used
c          K(2)  =   0 Deg LONG
c          K(74) = 360 Deg LONG
c          K(75) = 365 Deg LONG  =  5 Deg LONG
c........+.........+.........+.........+.........+.........+.........+..
C     \/ First , check coordinates
c........+.........+.........+.........+.........+.........+.........+..
c
      IF (FLATD.GT.90.0 .OR. FLATD.LT.-90.0)  THEN
         WRITE (*, 4001)  FLATD
c         WRITE (66,4001)  FLATD
      ENDIF
 4001 FORMAT (' 4001 ERROR LAT > +- 90; LAT = ',F10.3)

c........+.........+.........+.........+.........+.........+.........+..
c     \/ Check to prevent "blowup"
c        The B & L routine will also "blowup at 90 degrees
c........+.........+.........+.........+.........+.........+.........+..

      J = (90.0-FLATD)/5.0+1.00001
      IF (FLATD.GT.89.98)  THEN
         J =  0
         FLATD = 89.98
      ENDIF
      IF (FLATD.LT.-89.98)  THEN
         J = 37
         FLATD = -89.98
      ENDIF
c
      IF (FLOND.GT. 360.0)  FLOND = MOD(FLOND,360.0)
      IF (FLOND.LT.-360.0)  FLOND = MOD(FLOND,360.0)
      IF (FLOND.LT.0.0)     FLOND = 360.0+FLOND
      K = (FLOND/5.0)+2.00001
c
c........+.........+.........+.........+.........+.........+.........+..
c     \/ Set up interpolation
c                     point will be in a box
c            J,K      BLT           BRT      J,K+1
c                     BL     P      BR
c            J+1,K    BLB           BLB      J+1,K+1
c        Interpolate in box
c        Get left and right cutoff at lat FLATD
c........+.........+.........+.........+.........+.........+.........+..
c
c      ERR = 0.01
      ALT450 = 450.0D0
      UT = KUT*1.D0
c
      KP = K + 1
      JP = J + 1
c
c     ...+.........+.........+.........+.........+.........+.........+..
c     \/ Calculate weighting factors in longitude and latitude
c     ...+.........+.........+.........+.........+.........+.........+..
c
      FLONDK  = LNS(K)
      FLONDDF = FLOND - FLONDK
      WFLNR = FLONDDF/5.0
      WFLNL = 1.0 - WFLNR
c
      IF (JP.LE.37)  THEN
         FLATDJP = LTS(JP)
         FLATDJ  = LTS(J)
         FLATDDF = FLATD  - FLATDJP
         WFLTT = FLATDDF/5.0
      ENDIF
      IF (JP.GT.37)  THEN
         JP = 37
         WFLTT = 1.0
      ENDIF

      WFLTB = 1.0 - WFLTT
c
c     ...+.........+.........+.........+.........+.........+.........+..
c     \/ Calculate "L" for this position at 450 KM
c     ...+.........+.........+.........+.........+.........+.........+..
c
c      CALL INVARA (MODEL,TM,FLATD,FLOND,ALT450, ERR,BBPT,FLPT)
c	F.Lei - change to use the IRBEM lib               
c       call make_lstar1(ntime,kext,options,sysaxes,iyear,idoy,ut, x1,x2,x3, maginput, lm,lstar,blocal,bmin,xj,mlt)
      maginput(1,1) = FKPIDX*10.
      if (FKPIDX.GE.9.) maginput(1,1) = 89.999D0
c      KDOY = GET_DOY(KYEAR,KMONTH,KDATE)
c	this function is not working (need to change '01' to '1' in source code) 
c      KDOY = JULDAY(KYEAR,KMONTH,KDATE)-JULDAY(KYEAR, 1, 1)
c      UT = KHOUR*3600. +KMIN*60. + KSEC
      LYEAR = KYEAR
      IF (LYEAR .GE. 2015) LYEAR = 2014
      IF (LYEAR .LT. 1965) LYEAR = 1965
c   we have to use the map date/time rather than the given one
      call make_lstar1(1,4,options,0,LYEAR,1,UT,ALT450,FLATD,
     &	FLOND,maginput, lm,lstar,blocal,bmin,xj,mlt)
      FLPT = DABS(lm(1,1))
c

c     WRITE (*, 4002)  FLATD,FLOND,ALT450, ERR,BBPT,FLPT                ! Diag
c     WRITE (66,4002)  FLATD,FLOND,ALT450, ERR,BBPT,FLPT                ! Diag
c4002 FORMAT(' 4002 ', 2F7.2, F7.1, 3F7.2)                              ! Diag

      IF (FLPT.GT.999.99 .or. FLPT .lt. 0.)  FLPT = 999.99
      FLPSQ = FLPT*FLPT
c
c........+.........+.........+.........+.........+.........+.........+..
c    \/ Interpolate RL in latitude within "box" at 450 KM
c       Get Left side, Right side and weighted average of Left and Right
c........+.........+.........+.........+.........+.........+.........+..
c
c      RLBLTI = VKRL(J, K)/FLPSQ
c      RLBLBI = VKRL(JP,K)/FLPSQ
c      RLBLI  = RLBLTI*WFLTT + RLBLBI*WFLTB
c
c      RLBRTI = VKRL(J, KP)/FLPSQ
c      RLBRBI = VKRL(JP,KP)/FLPSQ
c      RLBRI  = RLBRTI*WFLTT + RLBRBI*WFLTB
c
c      RLPT = RLBLI*WFLNL +  RLBRI*WFLNR
c
c     ...+.........+.........+.........+.........+.........+.........+..
c     \/ Diagnostic RL outputs; normally turned off in operational version
c     ...+.........+.........+.........+.........+.........+.........+..
c
c      IF (IDG10.NE.0)  THEN                                             ! Diag 
c         WRITE (10,4010)                                                ! Diag
c     1   LTS(J), LNS(K),RLG(J, K), VKRL(J,K), FLG(J,K), J,K,            ! Diag
c     *   J,KP,          RLG(J, KP),VKRL(J,KP),FLG(J,KP),LTS(J),LNS(KP)  ! Diag
c         WRITE (10,4020)                                                ! Diag
c     2   FLATD,   FLOND,  RLBLI, RLBLTI,WFLTT,                          ! Diag
c     *                   RLBRI, RLBRTI,WFLTB,  FLPT,FLATD,FLOND         ! Diag
c         WRITE (10,4030)                                                ! Diag
c     3   LTS(JP),LNS(K),RLG(JP,K), VKRL(JP,K),FLG(JP,K),                ! Diag
c     *   JP,K, JP,KP,  RLG(JP,KP),VKRL(JP,KP),FLG(JP,KP),LTS(JP),LNS(KP)! Diag
c         WRITE (10,4040)  FLATD, FLOND, RLPT,RLBLI, WFLNL, RLBRI, WFLNR ! Diag
c      ENDIF                                                             ! Diag
c 4010 FORMAT (' 4010 BLT',2I7,  F7.2,  F7.2,F8.3, 2I3,2I3,              ! Diag
c     *                   F7.2,  F7.2,F8.3,2I7,  ' BRT RL')              ! Diag
c 4020 FORMAT (' 4020 BL ',2F7.2,F7.2,  F7.2,  F7.2,  6X,                ! Diag
c     *                   F7.2,  F7.2,  F7.2, F8.3,2F7.2,' BR  RL')      ! Diag
c 4030 FORMAT (' 4030 BLB',2I7,  F7.2,  F7.2,F8.3, 2I3,2I3,              ! Diag
c     *                   F7.2,  F7.2,F8.3,2I7  ,' BRB RL')              ! Diag
c 4040 FORMAT (' 4040    ',3F7.2,22X,4F7.2,' RL'/)                       ! Diag

c........+.........+.........+.........+.........+.........+.........+..
c    \/ Interpolate RU in latitude within "box" at 450 KM
c       Get Left side, Right side and weighted average of Left and Right
c........+.........+.........+.........+.........+.........+.........+..
c
c      RUBLTI = VKRU(J, K)/FLPSQ
c      RUBLBI = VKRU(JP,K)/FLPSQ
c      RUBLI  = RUBLTI*WFLTT + RUBLBI*WFLTB

c      RUBRTI = VKRU(J ,KP)/FLPSQ
c      RUBRBI = VKRU(JP,KP)/FLPSQ
c      RUBRI  = RUBRTI*WFLTT + RUBRBI*WFLTB
c
c      RUPT = RUBLI*WFLNL +  RUBRI*WFLNR
c
c     ...+.........+.........+.........+.........+.........+.........+..
c     \/ Diagnostic RU outputs; normally turned off in operational version
c     ...+.........+.........+.........+.........+.........+.........+..

c      IF (IDG11.NE.0)  THEN                                             ! Diag 
c         WRITE (11,4050)                                                ! Diag 
c     1   LTS(J), LNS(K),RUG(J, K), VKRU(J,K), FLG(J,K), J,K,            ! Diag 
c     *   J,KP,          RUG(J, KP),VKRU(J,KP),FLG(J,KP),LTS(J),LNS(KP)  ! Diag 
c         WRITE (11,4060)                                                ! Diag 
c     2   FLATD,   FLOND,  RUBLI, RUBLTI,WFLTT,                          ! Diag 
c     *                    RUBRI, RUBRTI,WFLTB,FLPT,FLATD,FLOND          ! Diag 
c         WRITE (11,4070)                                                ! Diag 
c     3   LTS(JP),LNS(K),RUG(JP,K),VKRU(JP,K), FLG(JP,K),                ! Diag 
c     *   JP,K, JP,KP,  RUG(JP,KP),VKRU(JP,KP),FLG(JP,KP),LTS(JP),LNS(KP)! Diag 
c         WRITE (11,4080)  FLATD, FLOND, RUPT, RUBLI, WFLNL, RUBRI, WFLNR! Diag 
c      ENDIF                                                             ! Diag 
c 4050 FORMAT (' 4050 BLT',2I7,  F7.2,  F7.2,F8.3, 2I3,2I3,              ! Diag 
c     *                   F7.2,  F7.2,F8.3,2I7,  ' BRT RL')              ! Diag 
c 4060 FORMAT (' 4060 BL ',2F7.2,F7.2,  F7.2,  F7.2,  6X,                ! Diag 
c     *                   F7.2,  F7.2,  F7.2, F8.3,2F7.2,' BR  RL')      ! Diag 
c 4070 FORMAT (' 4070 BLB',2I7,  F7.2,  F7.2,F8.3, 2I3,2I3,              ! Diag 
c     *                   F7.2,  F7.2,F8.3,2I7  ,' BRB RL')              ! Diag 
c 4080 FORMAT (' 4080    ',3F7.2,22X,4F7.2,' RU'/)                       ! Diag 
c
c........+.........+.........+.........+.........+.........+.........+..
c    \/ Interpolate RC in latitude within "box" at 450 KM
c       Get Left side, Right side and weighted average of Left and Right
c........+.........+.........+.........+.........+.........+.........+..
c
      RCBLTI = VKRC(J, K)/FLPSQ
      RCBLBI = VKRC(JP,K)/FLPSQ
      RCBLI  = RCBLTI*WFLTT + RCBLBI*WFLTB

      RCBRTI = VKRC(J ,KP)/FLPSQ
      RCBRBI = VKRC(JP,KP)/FLPSQ
      RCBRI  = RCBRTI*WFLTT + RCBRBI*WFLTB

      RCPT = RCBLI*WFLNL +  RCBRI*WFLNR
c
c     ...+.........+.........+.........+.........+.........+.........+..
c     \/ Diagnostic RC outputs; normally turned off in operational version
c     ...+.........+.........+.........+.........+.........+.........+..

c      IF (IDG12.NE.0)  THEN                                             ! Diag 
c         WRITE (12,4150)                                                ! Diag 
c     1   LTS(J), LNS(K),RCG(J, K), VKRC(J,K), FLG(J,K), J,K,            ! Diag 
c     *   J,KP,          RCG(J, KP),VKRC(J,KP),FLG(J,KP),LTS(J),LNS(KP)  ! Diag 
c         WRITE (12,4160)                                                ! Diag 
c     2   FLATD,   FLOND,  RCBLI, RCBLTI,WFLTT,                          ! Diag 
c     *                  RCBRI, RCBRTI,WFLTB,  FLPT,FLATD,FLOND          ! Diag 
c         WRITE (12,4170)                                                ! Diag 
c     3   LTS(JP),LNS(K),RCG(JP,K), VKRC(JP,K),FLG(JP,K),                ! Diag 
c     *   JP,K, JP,KP,  RCG(JP,KP),VKRC(JP,KP),FLG(JP,KP),LTS(JP),LNS(KP)! Diag 
c         WRITE (12,4180)                                                ! Diag 
c     *         FLATD, FLOND, RCPT,    RCBLI, WFLNL, RCBRI, WFLNR        ! Diag 
c      ENDIF                                                             ! Diag 
c 4150 FORMAT (' 4150 BLT',2I7,  F7.2,  F7.2,F8.3, 2I3,2I3,              ! Diag 
c     *                   F7.2,  F7.2,F8.3,2I7,  ' BRT RC')              ! Diag 
c 4160 FORMAT (' 4160 BL ',2F7.2,F7.2,  F7.2,  F7.2,  6X,                ! Diag 
c     *                   F7.2,  F7.2,  F7.2, F8.3,2F7.2,' BR  RC')      ! Diag 
c 4170 FORMAT (' 4170 BLB',2I7,  F7.2,  F7.2,F8.3, 2I3,2I3,              ! Diag 
c     *                   F7.2,  F7.2,F8.3,2I7  ,' BRB RC')              ! Diag 
c 4180 FORMAT (' 4180    ',3F7.2,22X,4F7.2,' RC'/)                       ! Diag 
c
c........+.........+.........+.........+.........+.........+.........+..
c     \/ Get VK constant for this point for altitude interpolation
c........+.........+.........+.........+.........+.........+.........+..
c
c      VKRLP =  RLPT*FLPSQ
c      VKRUP =  RUPT*FLPSQ
      VKRCP =  RCPT*FLPSQ

c      IF (IDG14.NE.0)  THEN                                             ! Diag 
c         WRITE (14, 4200)  FLPSQ, VKRLP,RLPT, VKRUP,RUPT, VKRCP,RCPT,   ! Diag 
c     *                     FLATD, FLOND                                 ! Diag 
c      ENDIF                                                             ! Diag 
c 4200 FORMAT (' 4200 ',F9.1, 6F8.2, 2F8.2)                              ! Diag 
c
c........+.........+.........+.........+.........+.........+.........+..
c........+.........+.........+.........+.........+.........+.........+..
c     >>>>>>>>>>   Interpolated cutoff values at 450 KM   <<<<<<<<<<
c........+.........+.........+.........+.........+.........+.........+..
c........+.........+.........+.........+.........+.........+.........+..
c
c      RLGVPV = RLPT
c      RUGVPV = RUPT
      RCGVPV = RCPT
c
c........+.........+.........+.........+.........+.........+.........+..
c     \/ Altitude interpolation
c        Store vertical cutoff rigidity in COMMON CUTOFF ARRAYS
c........+.........+.........+.........+.........+.........+.........+..
c
      IF (ALTKM.NE.450.0)  THEN
c         CALL INVARA (MODEL,TM,FLATD,FLOND,ALTKM, ERR, BBPT,FLKM)
c	F.Lei - change to use the IRBEM lib               
c            call make_lstar1(ntime,kext,options,sysaxes,iyear,idoy,ut, x1,x2,x3, maginput, lm,lstar,blocal,bmin,xj,mlt)
        maginput(1,1) = KPIDX*10.
        if (KPIDX.GE.9) maginput(1,1) = 89.999D0
         LYEAR = KYEAR
         IF (LYEAR .GE. 2015) LYEAR = 2014
         IF (LYEAR .LT. 1965) LYEAR = 1965   
	     call make_lstar1(1,4,options,0,LYEAR,KDOY,UT,ALTKM,FLATD,
     &                 FLOND,maginput,lm,lstar,blocal,bmin,xj,mlt)
	    FLKM = DABS(lm(1,1))
c
	    IF (FLKM.GT.999.99 .or. FLKM .lt. 0.)  FLKM = 999.99
        FLPT = FLKM
        FLKMSQ = FLKM*FLKM
c         RLGVPV = VKRLP/FLKMSQ
c         RUGVPV = VKRUP/FLKMSQ
        RCGVPV = VKRCP/FLKMSQ
c        WRITE (66,4220) FLKMSQ,RLGVPV,VKRLP,RUGVPV,VKRUP,RCGVPV,VKRCP, ! Diag 
c    *                   FLATD,VFLOND,VALTKM,FLKM                       ! Diag
      ENDIF
c4220 FORMAT (' 4220 ', F9.1, 6F8.2, 2F8.2,F8.1,F8.2)                   ! Diag 
c
c........+.........+.........+.........+.........+.........+.........+..
c     \/ Specify penumbral transparency in vertical direction
c        Penumbral Transparency now specified for each half GV
c                  for each Kp values
c........+.........+.........+.........+.........+.........+.........+..
c
c      IGH = RLGVPV*2.0
c      IF (IGH.LT.1)  IGH = 1
c      PTV = PTVDAT(IGH)
c
c........+.........+.........+.........+.........+.........+.........+..
c     Radial Distance Adjustment                           ! DS 02 08 28
c........+.........+.........+.........+.........+.........+.........+..
c     email from Don on 09/12/2013
c     " If you examine the cutoff interpolation FORTRAN code in detail, 
c       you may notice that there is an adjustment in the "L" value altitude
c       interpolation process at the end of Subroutine LINT5X5 in a section
c       labeled "Radial Distance Adjustment".  While the "L" interpolation
c       equation has the basic form of L**-2, when this exact form is used to
c       extend to geosynchronous altitude, the cutoff values extrapolated from
c       the near Earth low altitudes are too high; approximately 0.3 GV at the
c       magnetic equator.  (See figure 7 of Shea & Smart, JGR, 72, 3447, 1967)
c       We incorporated a "patch" (actually an "ad-hoc" exponential function
c       that makes adjustments so at 6.6 earth radii the vertical cutoff rigidity
c       for local noon at the magnetic equator under extremely quiet magnetic
c       conditions is about zero (or extremely small).
c       This "ad-hoc" exponential function is not going to be reliable
c       beyond geosynchronous distances."
c
c       This is probally unneccessary and the magcosmics simulation shown pricisely
c       ~0.4 Gv cutoff at the magnetic equator. So this correction is taking out

c      RADIST = (ALTKM + 6371.2)/6371.2                    ! DS 02 08 28
c      RADISTSQ = RADIST*RADIST                            ! DS 02 08 28
c      RCORT = (LOG(RADISTSQ))/11.0                        ! DS 02 08 28

c      RCGVPV = RCGVPV - rcort                             ! DS 02 08 28

c      If (RCGVPV.le.0.0)  RCGVPV = 0.0                    ! DS 02 08 28

c........+.........+.........+.........+.........+.........+.........+..
c
      RETURN
c
      END

!**********************************************************************!
!     SUBROUTINE CUTEW                                         ! SUB05 !
!     "L" INTERPOLATION PROGRAM (VERSION 4.03)                         !
!**********************************************************************!

       SUBROUTINE CUTEW (FLOND,FLATD,ALTKM)
c........+.........+.........+.........+.........+.........+.........+..
c     Interpolate vertical cutoffs to 90 deg East and West cutoffs
c     Version 4.03; 3-HR UT Version; Uses 3-hr avg RU, RC, RL & PW for all KP
C             INCOMPATIBLE WITH LOWER LEVEL VERSIONS
c     Remember: particle cutoff rigidity in GV   IN  COMMON /CUTRIG/
c     Remember: proton   cutoff energy   in MEV  IN  COMMON /CUTERGY/
c     Remember: "B & L" related parameters       IN  COMMON /BANDL/
c........+.........+.........+.........+.........+.........+.........+..
CLast Mod 06 Apr 2003  Testing 
c     Mod 20 Feb 2003  Identify AFRL GEOSPACE additions for SGI compatibility
c     Mod 25 Jun 2002  Cutoff values in COMMON CUTRIG AND CUTERGY 
c     Mod 25 May 2002  Current YR, MO, DT, HR, KMIN in COMMON /KYMDHM/
c     Mod 09 Oct 2000  REAL*8 for compatibility with NASA JSC SRAG
c........+.........+.........+.........+.........+.........+.........+..
c     Programmed by Don F. Smart  (SSSRC@MSN.COM)
c     Note: - The programming adheres to the conventional FORTRAN
c             default standard that variables beginning with
c             'I','J','K','L','M',or 'N' are integer variables
c             Variables beginning with "C" are character variables
C             Exception, the variable  str  is a character variable     ! AFRL
C             All other variables are REAL*8 
c........+.........+.........+.........+.........+.........+.........+..
c
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-B)
      IMPLICIT REAL*8 (D-H)
      IMPLICIT REAL*8 (O-Z)
c      
c     Added by F.Lei      
      Real*8 xIn(3), xOUT(3)
c
c........+.........+.........+.........+.........+.........+.........+..

      COMMON /BANDL/   BBPT, FLPT, GLMDAR

      COMMON /CUTRIG/  RLGVPV,RUGVPV,RCGVPV, RLGVPE,RUGVPE,RCGVPE, 
     &                 RLGVPW,RUGVPW,RCGVPW, PTV
c     COMMON /CUTERGY/ EPNRLV,EPNRUV,EPNRCV, EPNRLE,EPNRUE,EPNRCE,
c     &                 EPNRLW,EPNRUW,EPNRCW 

      COMMON /KYMDHM/ KYEAR,KMONTH,KDATE,KHOUR,KMIN,KSEC,KDOY,KUT

      SAVE                                                              ! AFRL

c........+.........+.........+.........+.........+.........+.........+..
c     \/ Get equivalent magnetic latitude
C         Get Corrected geomagnetic latitude (GMLATC) at sub-satellite point
c         Then, get Invariant latitude at satellite position
c                   INVARIANT LAT = ACOS(1.0/SQRT(L))
c         Select the smaller value for magnetic latitude
c       Remember, keep degrees and radians in proper context
c       Also, we will always use absolute value of equivalent magnetic latitude
c........+.........+.........+.........+.........+.........+.........+..

c       CALL CGLALO95 (FLATD, FLOND, GMLATCD, GMLONCD)
c	F.Lei change to use IRBEM lib
c          call coord_trans_vec1(ntime,sysaxesIN,sysaxesOUT,iyr,idoy,secs,xIN,xOUT)
      xIN(1) = 0.      ! on the surface
      xIN(2) = FLATD
      xIN(3) = FLOND
c     GEO->MAG 
      LYEAR = KYEAR
      IF (LYEAR .GE. 2015) LYEAR = 2014
      IF (LYEAR .LT. 1965) LYEAR = 1965
      call coord_trans_vec1(1,0,6,LYEAR,1,
     &            UT,xIN,xOUT)
c     calculate the latitude     
      call CAR_SPH(xOUT,r,flati,flongi)
c     flati is the magnetic latitude in degrees. Note this is not the same as 
c     the corrected geomagnetic latitude, but the difference should be small
      GMLATCD = flati
      
      GMLATCR = GMLATCD/57.2957795

      GLMDAR = 0.0
      IF (FLPT.GT.1.0D0)  GLMDAR = DACOS((1.0D0/SQRT(FLPT)))
      GLMDAD = GLMDAR*57.2957795 

c      IF (DABS(GMLATCR).LT.(DABS(TCGLMDA)))  GLMDAR = DABS(GMLATCR)
c flei - should this line be corrected as
      IF (DABS(GMLATCR).LT.(DABS(GLMDAR)))  GLMDAR = DABS(GMLATCR)
c........+.........+.........+.........+.........+.........+.........+..
c     \/ Get cosine cubed equivalent magnetic invariant latitude
c........+.........+.........+.........+.........+.........+.........+..

c      TCGLMDA =DCOS(GLMDAR)
c      TCGL3   = TCGLMDA**3
c      IF (TCGL3.GT.10D0)  TCGL3   = 1.0

c........+.........+.........+.........+.........+.........+.........+..
c     \/ Get East and West factors
c........+.........+.........+.........+.........+.........+.........+..

c      EFACTOR = 1.0
c      WFACTOR = 1.0
c      EFACTOR = (SQRT(1.0 - TCGL3) + 1.0)**2
c      WFACTOR = (SQRT(1.0 + TCGL3) + 1.0)**2

c     WRITE (88,5110)  TCGLMDA, TCGL3, FLPT, RLGVPV, RUGVPV, RCGVPV     ! Diag
c5110 FORMAT (' 5110',    2F7.4,        F7.3, 3F7.3)                    ! Diag

c........+.........+.........+.........+.........+.........+.........+..
c     \/ Do rigidity to energy conversion for protons
C               need to define proton atomic number, charge and mass
c        Convert geomagnetic cutoff rigidity to proton cutoff energy
c        Rigidity in units of GV;
c        Rigidity to energy subroutine wants rigidity in units of MV
c........+.........+.........+.........+.........+.........+.........+..
c
c      NA = 1
c      NZ = c          call coord_trans_vec1(ntime,sysaxesIN,sysaxesOUT,iyr,idoy,secs,xIN,xOUT)1
c      PAMU = 1.0081451
c
      IF (RCGVPV.LE.0.0)  RCGVPV = 1.0E-06
c      IF (RLGVPV.LE.0.0)  RLGVPV = 1.0E-06
c      IF (RUGVPV.LE.0.0)  RUGVPV = 1.0E-06
c
c........+.........+.........+.........+.........+.........+.........+..
c     \/ Get vertical proton cutoff energy
c        First, get vertical cutoff rigidity in MV for this position
c........+.........+.........+.........+.........+.........+.........+..
c
c      RCMVPV = RCGVPV*1000.0
c      RLMVPV = RLGVPV*1000.0
c      RUMVPV = RUGVPV*1000.0
c
c     WRITE (88,5120)  RLMVPV, RUMVPV, RCMVPV
c5120 FORMAT (' VMV', 23X, 3F7.0)
c
c     WRITE (66,3130)                                                   ! DIAG
c     WRITE (66,3140)        NA,NZ,PAMU,  RLMVPV,  EPN, BETA            ! DIAG

c      CALL AZRGEG (NA,NZ,PAMU,RCMVPV,EPNRCV,BETACV)           !SUB06
c      CALL AZRGEG (NA,NZ,PAMU,RLMVPV,EPNRLV,BETALV)           !SUB06
c      CALL AZRGEG (NA,NZ,PAMU,RUMVPV,EPNRUV,BETAUV)           !SUB06

c     WRITE (66,5150)
c     WRITE (66,5160)        NA,NZ,PAMU,  RCMVPV, EPNRCV, BETACV        ! DIAG
c     WRITE (66,5160)        NA,NZ,PAMU,  RLMVPV, EPNRLV, BETALV        ! DIAG
c     WRITE (66,5160)        NA,NZ,PAMU,  RUMVPV, EPNRUV, BETAUV        ! DIAG
c5130 FORMAT (' CALL AZRGEG (NA,NZ,     PAMU,   RLVVPV,',               ! DIAG
c    *        '      EPN,     BETA)')                                   ! DIAG
c5140 FORMAT (' 5140',8X, 2I3,3F10.5,F10.5)                             ! DIAG
c5150 FORMAT (' RET  AZRGEG (NA,NZ,     PAMU,   RLMVPV,',               ! DIAG
c    *        '      EPN,     BETA)')                                   ! DIAG
c5160 FORMAT (' 5160',8X, 2I3,3F10.3,F10.5)                             ! DIAG
c
c      IF (EPNRCV.LE.0.0)  EPNRCV = 1.0E-06
c      IF (EPNRLV.LE.0.0)  EPNRLV = 1.0E-06
c      IF (EPNRUV.LE.0.0)  EPNRUV = 1.0E-06
c
c     WRITE (88,5170)  EPNRLV, EPNRUV, EPNRCV                           ! DIAG
c5170 FORMAT (' VPE', 23X, 3F7.0)                                       ! DIAG
c
c........+.........+.........+.........+.........+.........+.........+..
c     \/ Get interpolate  West proton cutoff energy
c........+.........+.........+.........+.........+.........+.........+..
c
c      RCGVPW = 4.0*RCGVPV/WFACTOR
c      RLGVPW = 4.0*RLGVPV/WFACTOR
c      RUGVPW = 4.0*RUGVPV/WFACTOR
c
c      IF (RCGVPW.LE.0.0)  RCGVPW = 1.0E-06
c      IF (RLGVPW.LE.0.0)  RLGVPW = 1.0E-06
c      IF (RUGVPW.LE.0.0)  RUGVPW = 1.0E-06
c
c      RCMVPW = RCGVPW*1000.0
c      RLMVPW = RLGVPW*1000.0
c      RUMVPW = RUGVPW*1000.0
c
c     WRITE (88,5180)  RLMVPW, RUMVPV, RCMVPV, WFACTOR                  ! DIAG
c5180 FORMAT (' WMV', 23X, 3F7.0,F7.2)                                  ! DIAG
c
c      CALL AZRGEG (NA,NZ,PAMU,RCMVPW,EPNRCW,BETACW)           !SUB06
c      CALL AZRGEG (NA,NZ,PAMU,RLMVPW,EPNRLW,BETALW)           !SUB06
c      CALL AZRGEG (NA,NZ,PAMU,RUMVPW,EPNRUW,BETAUW)           !SUB06
c
c      IF (EPNRCW.LE.0.0)  EPNRCW = 1.0E-06
c      IF (EPNRLW.LE.0.0)  EPNRLW = 1.0E-06
c      IF (EPNRUW.LE.0.0)  EPNRUW = 1.0E-06
c
c     WRITE (88,5190)  EPNRLW, EPNRUW, EPNRCW                           ! DIAG
c5190 FORMAT (' WPE', 23X, 3F7.0)                                       ! DIAG
c
c........+.........+.........+.........+.........+.........+.........+..
c     \/ Get interpolate  East proton cutoff energy
c........+.........+.........+.........+.........+.........+.........+..
c
c      RCGVPE = 4.0*RCGVPV/EFACTOR
c      RLGVPE = 4.0*RLGVPV/EFACTOR
c      RUGVPE = 4.0*RUGVPV/EFACTOR
c
c      IF (RCGVPE.LE.0.0)  RCGVPE = 1.0E-06
c      IF (RLGVPE.LE.0.0)  RLGVPE = 1.0E-06
c      IF (RUGVPE.LE.0.0)  RUGVPE = 1.0E-06
c
c      RCMVPE = RCGVPE*1000.0
c      RLMVPE = RLGVPE*1000.0
c      RUMVPE = RUGVPE*1000.0
c
c     WRITE (88,5200)  RLMVPE, RUMVPE, RCMVPE, EFACTOR                  ! DIAG
c5200 FORMAT (' EMV', 23X, 3F7.0,F7.2)                                  ! DIAG
c
c      CALL AZRGEG (NA,NZ,PAMU,RCMVPE,EPNRCE,BETACE)           !SUB06
c      CALL AZRGEG (NA,NZ,PAMU,RLMVPE,EPNRLE,BETALE)           !SUB06
c      CALL AZRGEG (NA,NZ,PAMU,RUMVPE,EPNRUE,BETAUE)           !SUB06
c
c      IF (EPNRCE.LE.0.0)  EPNRCE = 1.0E-06
c      IF (EPNRLE.LE.0.0)  EPNRLE = 1.0E-06
c      IF (EPNRUE.LE.0.0)  EPNRUE = 1.0E-06

c     WRITE (88,5210)  EPNRLE, EPNRUE, EPNRCE
c5210 FORMAT (' EPE', 23X, 3F7.0)
c
      RETURN
      END

!**********************************************************************!
!     SUBROUTINE AZRGEG                                        ! SUB06 !
!     RIGIDITY TO ENERGY AND VISA VERSA                                !
!**********************************************************************!

c      SUBROUTINE AZRGEG (NA,NZ,PAMU,RIGIN,EPN,BETA)
c
c........+.........+.........+.........+.........+.........+.........+..
c     SUBROUTINE TO CONVERT RIGIDITY TO ENERGY AND VISA VERSA
CLast Mod 20 Feb 2003  Testing 
c     Mod 20 Feb 2003  Identify AFRL GEOSPACE additions for SGI compatibility
c     Mod 09 Oct 2000  REAL*8 for compatibility with NASA JSC SRAG
c     Mod 17 Mar 94 
c........+.........+.........+.........+.........+.........+.........+..
c     Programmed by Don F. Smart  (SSSRC@MSN.COM)
c     Note: - The programming adheres to the conventional FORTRAN
c             default standard that variables beginning with
c             'I','J','K','L','M',or 'N' are integer variables
c             Variables beginning with "C" are character variables
C             All other variables are REAL*8 
c........+.........+.........+.........+.........+.........+.........+..
c
c      IMPLICIT INTEGER (I-N)
c      IMPLICIT DOUBLE PRECISION (A-B)
c      IMPLICIT DOUBLE PRECISION (D-H)
c      IMPLICIT DOUBLE PRECISION (O-Z)

c      SAVE                                                              ! AFRL
c
c........+.........+.........+.........+.........+.........+.........+..
c
c     WRITE (*,6010) NA,NZ,PAMU,RIGIN,EPN                               ! DIAG
c6010 FORMAT (' AZRGEG1 ',2I5,3F15.3)                                   ! DIAG
c
c........+.........+.........+.........+.........+.........+.........+..
c     EPAMU IS REST MASS ENERGY PER ATOMIC MASS UNIT
c     CHECK, IF NA, NZ, OR PAMU NOT SPECIFIED, PUT IN DEFAULT FOR PROTONS
c........+.........+.........+.........+.........+.........+.........+..
c
c      IF (PAMU.LE.0.0)  PAMU = 1.0081451
c      IF (NA.LE.0)      NA = 1
c      IF (NZ.LE.0)      NZ = 1
c
c      EPAMU = 931.141
c
c      ANUC = NA
c      ZCHARG = NZ
c      RSMSPN = (PAMU/ANUC)*EPAMU
c
c      TRIG = RIGIN
c
c........+.........+.........+.........+.........+.........+.........+..
c     \/ IF TRIG .LE. 0.0    DO ENERGY   TO RIGIDITY CONVERSION
c        IF TRIG .GT. 0.0    DO RIGIDITY TO ENERGY   CONVERSION
c........+.........+.........+.........+.........+.........+.........+..
c
c      IF (TRIG.LE.0.0)  THEN
c                                    ENERGY TO RIGIDITY CONVERSION
c         GMAEG = (ANUC*EPN+ANUC*RSMSPN)/(ANUC*RSMSPN)
c         GMAEGG = (EPN+RSMSPN)/RSMSPN
c         RIGIN = SQRT(GMAEG*GMAEG-1.0)*RSMSPN*ANUC/ZCHARG
c         RELGMA = GMAEG
c      ELSE
c                      RIGIDITY TO ENERGY CONVERSION
c         GMARG = SQRT(((RIGIN*ZCHARG)/(RSMSPN*ANUC))**2+1.0)
c         EPN = (GMARG-1.0)*RSMSPN
c         RELGMA = GMARG
c      ENDIF
c
c     WRITE (*,6020) ANUC,ZCHARG,RSMSPN,RIGIN,EPN,GMAEG,GMAEGG,GMARG    ! DIAG
c     WRITE (*,6030) NA,NZ,PAMU,RIGIN,EPN                               ! DIAG
c6020 FORMAT (' AZRGEG2 ',8F10.3)                                       ! DIAG
c6030 FORMAT (' AZRGEG3 ',2I5,3F15.5)                                   ! DIAG
c
c      BETA = SQRT(1.0-1.0/(RELGMA*RELGMA))
c
c      RETURN
c
c          BETA IS V/C (SPEED AS FRACTION OF LIGHT SPEED)
c          ENERGY IN MEV
c          EPAMU IS MASS-ENERGY CONVERSION = 931.141 MEV PER AMU
c          EPN IS KINETIC ENERGY PER NUCLEON
c          NA IS ATOMIC NUMBER
c          NZ IS CHARGE
c          PAMU IS REST MASS IN PHYSICAL ATOMIC MASS UNITS
c          RELGAM IS RELATIVISTIC FACTOR 'GAMMA'
c          RIGIDITY IN MV
c          RSMSPN IS REST MASS PER NUCLEON IN MEV
c
c      END
C
C----------------------------------------------------------------------
C
      SUBROUTINE MF_ORB(LUN, PRODEF, TITLE, MISVAR, ORBVAR, NRC, RC)

      IMPLICIT NONE

      INCLUDE 'interface.h'

      REAL*8 RC(*)
      INTEGER*4 LUN, NRC, MF(8), I
      CHARACTER*(*) PRODEF, TITLE
      CHARACTER*80 SDUM

      TYPE(META_VAR) :: MISVAR, ORBVAR

      MF(3) = MISVAR%N + ORBVAR%N + 2
      MF(4) = 0
      MF(5) = 6
      MF(6) = NRC + 5
      MF(7) = -1
      MF(8) = 0
      CALL MF_HEADER(LUN, MF, PRODEF, TITLE)
      WRITE (LUN,'(A)') '''MOD_ABB'', -1,''MSM'''
      CALL MF_ANNOT(LUN, MISVAR)
      CALL MF_ANNOT(LUN, ORBVAR)
      WRITE (SDUM,'(A,I2,A)') '(A,I3,1H,,', NRC, '(F7.2,1H,),A)'
      WRITE (LUN, TRIM(SDUM)) '''RIGIDITY'',', NRC, (RC(I),I=1,NRC),
     &  '''GV'''
      WRITE(LUN, '(A)') '''Kp'','' '', 1,''Kp'''
      WRITE(LUN, '(A)') '''Lm'',''Re '', 1,''McIlwain L'''
      WRITE(LUN, '(2A)') '''M_lati'',''degrees'', 1,',
     &  '''Magnetic latitude'''
      WRITE(LUN, '(2A)') '''RC'',''GV'', 1,',
     &  '''Vertical rigidity cut-off'''
      WRITE(LUN, '(2A)') '''SF'','' '', 1,',
     &  '''Earth trans/shadow factor'''
      WRITE(LUN, '(A,I3,A)') '''MS'','' '',', NRC,
     &  ',''Magnetic transmission factor for '',''RIGIDITY'''

      RETURN
      END

C----------------------------------------------------------------------
c      SUBROUTINE EXPFAC(RADIUS, MLAT, E, Z, A, F, N, ISTORM)
C computes an exposure factor at B,L for particles Z,A 
C arriving from interplanetary space with energies E(N)
C by averaging over arrival directions. Exposure factor
C returned in F(N).
C                    E.J. Daly ESA/ESTEC/WMA  9/87
c Updated by F. Lei 03/14
c 1) change to rigidity based in stead of energy, so no need of Z and A
c 2) the Rcv is given as input, no need for RADIUS and ISTORM
c 3) use omcao2
      SUBROUTINE EXPFAC(RCV, MLAT, RC, F, N)

      IMPLICIT NONE

      REAL*8 RC(*), F(*), RCV, MLAT
C  table of One-Minus-Cos-Angles-Over-2 :
      REAL*4 omcao2(9)/0., .067, .146, .25, .5, .75, .854, .933, 1./
c      REAL*8 opcao2(9)/1., .933, .854, .75, .5, .25, .146, .067, 0./
      REAL*8 ang(9)/0.01, .5236, .785, 1.047, 1.571, 2.094, 2.356,
     &              2.618, 3.1416/
      REAL*8 CUTFAC, cut(9),cutf
      INTEGER*4 N, IANGLE, NANGLE, IRC

      DATA nangle /9/

      IF (RCV .LE. 0.1) THEN
        DO IRC=1,N
           F(IRC) = 1.0
        END DO
      ELSE
        DO iangle=1,nangle
           cutf = CUTFAC(MLAT, ANG(IANGLE))
           CUT(IANGLE) = 4*RCV * cutf
        END DO
c        print *,RCV, mlat, cut(1), cut(5), cut(9)
c
        DO IRC = 1,N
           F(IRC) = 0.0
*  is there an value in the array above the lowest (Eastward, from West!) cutoff?:
           IF (RC(IRC) .GE. cut(1)) THEN
              F(IRC) = 1.0
              DO iangle=2,nangle
*     find the angular location where the cutoff goes over the rc:
                 IF (RC(IRC) .LE. cut(iangle)) THEN
                    F(IRC) = omcao2(iangle-1)+ (RC(IRC)-cut(iangle-1))
     &                      *(omcao2(iangle)-omcao2(iangle-1)) /
     &                      (cut(iangle)-cut(iangle-1))
                    GO TO 1
                 END IF
              END DO
 1            CONTINUE
           END IF
        END DO
      END IF

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      REAL*8 FUNCTION FACSHADOW(R)

C   This is a correction factor for the earth's shadow on 
C   the spacecraft according to simple geometrical optics.

      IMPLICIT NONE

      REAL*8 R

      FACSHADOW = 1.0D0 - 0.5D0 * (1.0D0-DSQRT(R**2-1.0D0)/R)

      RETURN

      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      REAL*8 FUNCTION CUTFAC(FMLAT,A)

      REAL*8 FMLAT, A, COSA, COSL

      COSA = DCOS(A)
      COSL = DCOS(FMLAT)
      CUTFAC = 1.D0/(1.D0+DSQRT(1.0D0+COSA*COSL**3))**2
      RETURN
      END
