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
      REAL*8   RCGVPV1

      INTEGER*4 NTIME_MAX
      PARAMETER (NTIME_MAX=100000)
      REAL*8 FMJD(NTIME_MAX), LON(NTIME_MAX), LAT(NTIME_MAX)
      REAL*8 ALT(NTIME_MAX), DT(NTIME_MAX)
      REAL*8 FKPIDX(100) /100*0.0D0/, A_ERR /5.0D0/, B_ERR /0.1D0/
      REAL*8 TMSTART, TMEND, TMDUR, OSEC, TSTART, TEND
      INTEGER*4 NKP /1/, OYEAR, OMON, ODAY, OHR, OMN, NTRAJ, IKP
      INTEGER*4 DATEVALS(8), KSKIP /0/, HEADER(2000), SAOLUN /3/
      CHARACTER*80 PRODEF, TITLE, ORBTIT, FOOTER
      CHARACTER*3 MONTH(12) /'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
     &                       'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'/

      CHARACTER*1000 filein, fileout, arg, dir
c
      INTEGER NRC
      PARAMETER (NRC=34)
      REAL*8 RC(NRC), F(NRC), FT1(NRC), FT2(NRC)
      REAL*8 FLUCUT(NRC,100), FLUCUTE(NRC,100)

c
C Orbit header common block
      COMMON /OHEAD/ HEADER
      COMMON /MSMHOME/ dir

c Name list definition
      namelist/msmlist/filein, fileout, dir, fkpidx, nkp, a_err, b_err
c
      TYPE(META_VAR) :: MISVAR, ORBVAR

c........+.........+.........+.........+.........+.........+.........+..
c     \/ The following data statement initializes control !DS 2 Jun 2002
c            need for the first call to SETKP3UT          !DS 2 Nov 2004
c                  these values are modified in SETKP3UT  !DS 2 Nov 2004
c........+.........+.........+.........+.........+.........+.........+..

      DATA IYROLD,IMOLD,IHROLD,M5OLD,KPDXOLD,CYEAR 
     &        / 0,    0,    25,   -1,     11, '1995'/
c
c     the rigidity values to be used for calculating the transmission function
c
      DATA RC/0.1, 0.2, 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5
     &, 6., 6.5, 7., 7.5, 8., 9., 10., 11., 12., 13, 14., 15., 16., 17.,
     &  20., 25., 20., 30., 40., 50., 60./
c
c    take inputs from command line if supplied, otherwise from namelist file msm.nml
c
c........+.........+.........+.........+.........+.........+.........+..
      integer :: count

      INCLUDE 'putils.version'

      CALL DATE_AND_TIME(VALUES=DATEVALS)
      WRITE (SDATE,'(I2,''-'',A,''-'',I4)')
     &  DATEVALS(3), MONTH(DATEVALS(2)), DATEVALS(1)
      WRITE (STIME,'(2(I2.2,'':''),I2.2)') (DATEVALS(I),I=5,7)

      count = command_argument_count()
      if (count .ge. 4) then
         CALL get_command_argument(1, arg)
         filein = TRIM(ADJUSTL(arg))
         CALL get_command_argument(2, arg)
         read(arg,*) fkpidx(1)
         nkp = 1
         CALL get_command_argument(3, arg)
         fileout = TRIM(ADJUSTL(arg))
         CALL get_command_argument(4, arg)
         dir = TRIM(ADJUSTL(arg))
         a_err = 5./100
         b_err = 0.1
      elseif (count .eq. 6) then
         call get_command_argument(5, arg)
         read(arg,*) a_err
         a_err = a_err/100.
         call get_command_argument(6, arg)
         read(arg,*) b_err
      elseif (count .eq. 0) then
c
c     read in the name list definition
c
         nmlstat = 0
         open(unit=9,file='msm.nml',status='old',err=5)
         read(9,msmlist,err=5)
         close(9)
         a_err = a_err/100.
         nmlstat = 1
 5       if (nmlstat .eq. 0) then
            print *, ' Namelist file - msm.nml - error'
            call exit
         end if
      else
        print*, 'To Run: $ ./msm '
        print*, '   or   $ ./msm orbit-file kp output-trans-file MSMHOME
     &_DIR'
        print*, '   or   $ ./msm orbit-file kp output-trans-file MSMHOME
     &_DIR a_err b_err'
        call exit()
      end if
c
C Read orbit header
c
      OPEN (SAOLUN, file=filein, status='old')
      CALL READ_ORB(PRODEF, TITLE, NTRAJ, TMSTART, TMEND, TMDUR, ORBTIT,
     &              OYEAR, OMON, ODAY, OHR, OMN, OSEC, TSTART, TEND,
     &              SAOLUN)
      CALL MV_READ('ORB', ORBVAR)
      CALL MV_READ('MIS', MISVAR)
      NT = 0
      FOOTER = '*CONTINUE*'
C Read in the orbit
      DO WHILE (FOOTER .EQ. '*CONTINUE*')
        NT = NT + 1
        IF (NT .GT. NTIME_MAX) THEN
          WRITE(6,*) 'Maximum number of trajectory points exceeded'
          STOP 1
        END IF
        CALL READ_EPH(FMJD(NT), ALT(NT), LON(NT), LAT(NT), SPS, FOOTER,
     &                KSKIP)
      END DO
      NT = NT - 1
      CLOSE(SAOLUN)
C Calculate time intervals
      DT(1) = (FMJD(2)-FMJD(1)) / 2.0D0
      DO I=2,NT-1
        DT(I) = (FMJD(I+1)-FMJD(I-1)) / 2.0D0
      END DO
      DT(NT) = (FMJD(NT)-FMJD(NT-1)) / 2.0D0

c open the output file
      OPEN (UNIT=9, FILE=fileout, STATUS='UNKNOWN')

C Loop over Kp values
      DO IKP=1,NKP
c........+.........+.........+.........+.........+.........+.........+..
c     \/ Load in default kp0 t00 world grid (UNIT1)
c........+.........+.........+.........+.........+.........+.........+..
        KPIDX = 0
        IUT   = 0
        CALL SETKP3UT (CYEAR,KPIDX,KPDXOLD, IUT,IUTOLD, FKPIDX(IKP))                  ! SUB02
        KPIDX = NINT(fkpidx(IKP))
        CALL MF_ORB(9, PRODEF, TITLE, MISVAR, ORBVAR, NRC, RC, NT, IKP,
     &            NKP)
        DO I=1,NRC
            FLUCUT(I,IKP) = 0.0D0
            FLUCUTE(I,IKP) = 0.0D0
        END DO
c
c........+.........+.........+.........+.........+.........+.........+..
c    !!! Beginning of logical loop, will continue until end of data
c     \/ Read data point
c        Get year month date hour of data point
c........+.........+.........+.........+.........+.........+.........+..

c      \/ Set the default date and time
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
c
C Loop over the trajectory
        DO IT=1,NT

c........+.........+.........+.........+.........+.........+.........+..
c        Remember, cutoff are calculated for > 20 KM,
c                  IF ALT <20, use 20 KM
c........+.........+.........+.........+.........+.........+.........+..
c
            ALTKM = MAX(ALT(IT), 20.0D0)
c........+.........+.........+.........+.........+.........+.........+..

            FLATD = LAT(IT)
            IF (FLATD.GT. 89.98)  FLATD =  89.98
            IF (FLATD.LT.-89.98)  FLATD = -89.98
            FLOND = LON(IT)
            IF (FLOND.LT.0.0)     FLOND = FLOND + 360.0
            IF (FLOND.GT.360.0)   FLOND = FLOND - 360.0

c........+.........+.........+.........+.........+.........+.........+..
            JDAY = INT(FMJD(IT) + JULDAY(1950,1,1) - 0.5D0)
            CALL CALDAT(JDAY, KYEAR, KMONTH, KDATE)
            KDOY = GET_DOY(KYEAR, KMONTH, KDATE)
            UT = (FMJD(IT)-INT(FMJD(IT))) * 86400.0D0
c      print *, KYEAR, KMONTH, KDATE, (FMJD(IT) - DINT(FMJD) )
c
            CALL DOY_AND_UT2DATE_AND_TIME(KYEAR, KDOY, UT, KMONTH,
     &                                    KDATE, KHOUR, KMIN, KSEC)
c     
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
     &                  FLOND,FLATD,ALTKM,FKPIDX(IKP))          ! SUB03

c........+.........+.........+.........+.........+.........+.........+..
c     \/ Call east-west cutoff interpolation subroutine (arguments in common)
c        Get invariant latitude degrees
c........+.........+.........+.........+.........+.........+.........+..

            CALL CUTEW (FLOND,FLATD,ALTKM) ! MOD 27 AUG 02

C NOW DECIDE IF INTERPOLATION IS REQUIRED IN BETWEEN THE YEARS
c
            IF ((LYEAR .GT. 1955) .AND. (LYEAR .LE. 2015)) THEN
                YEAR = KYEAR + (KDOY+KUT/86400.D0)/365.D0
                DFY = YEAR - LYEAR
                IF (DFY .GT. 0.1) THEN
                    BBPT1 = BBPT; FLPT1 = FLPT ; GLMDAR1 = GLMDAR
                    RCGVPV1 = RCGVPV
                    LYEAR = LYEAR + 5
                    IF (LYEAR .GT. 2015) LYEAR = 2015
                    F1 = DFY/5.D0 ; F2 = 1. - F1
                    WRITE(CYEAR,'(I4)') LYEAR

                    CALL CUTUT3(CYEAR,KPIDX,KPDXOLD,IUT,
     &                       IUTOLD,FLOND,FLATD,ALTKM,FKPIDX(IKP))
                    CALL CUTEW (FLOND,FLATD,ALTKM)
C 
                    BBPT = F1*BBPT1+F2*BBPT; FLPT = F1*FLPT1
     &                         + F2*FLPT
                    GLMDAR = F1*GLMDAR1 + F2*GLMDAR
                    RCGVPV = F1*RCGVPV1 + F2*RCGVPV
                END IF
            END IF

            FLMDAD = GLMDAR*57.2957795
C
C     NOW CALCULATE THE TRANSMISSION FUNCTION
C
            CALL EXPFAC(RCGVPV, GLMDAR, RC, F, NRC)

c       calculate the associated errors and store them in FT1
            dR = a_err*RCGVPV
            if (b_err .gt. dR) dR = b_err
            call EXPFAC(RCGVPV+dR, GLMDAR, RC, FT1, NRC)
            call EXPFAC(RCGVPV-dR, GLMDAR, RC, FT2, NRC)
            DO I=1,NRC
                FT1(I) = abs(FT1(I)-FT2(I))/2.
            END DO
c      the Earth showdow
            SHADOW = FACSHADOW(1.D0+ALTKM/6371.D0)
   
c      WRITE (9,*) FLATD,FLOND,FLPT,RCGVPV,FACSHADOW(1.D0+ALTKM/6371.D0),
c      WRITE (9,*) FLATD,FLOND,FLPT,RCGVPV,EPNRCV
            WRITE(9, '(F8.4,1P,100('','',E13.4:),0P)')
     &              FKPIDX(IKP), FLPT, FLMDAD, RCGVPV, SHADOW, F, FT1

c
            DO I=1,NRC
                FF = DT(IT) * SHADOW
                FLUCUT(I,IKP) = FLUCUT(I,IKP) + FF*F(I)
c                FLUCUTE(I,IKP) = FLUCUTE(I,IKP) + FF**2 * FT1(I)**2
c      the error cannot be treated as standard deviation and it is a function of its owen!
                FLUCUTE(I,IKP) = FLUCUTE(I,IKP) + FF*FT1(I)
            END DO

c........+.........+.........+.........+.........+.........+.........+..
c    /\ End of logical loop, end of data, i.e. NT loop
c........+.........+.........+.........+.........+.........+.........+..
        END DO
        WRITE(9, '(A)') '''End of Block'''

C End loop over Kp
      END DO

C Write trajectory average cut-off spectrum
      ORBDUR = FMJD(NT) - FMJD(1)
      DO I=1,NRC
        DO J=1,NKP
          FLUCUT(I,J) = FLUCUT(I,J) / ORBDUR
          FLUCUTE(I,J) = FLUCUTE(I,J) / ORBDUR
        END DO
      END DO
      CALL MF_TRI(9, PRODEF, TITLE, NRC, RC, FLUCUT, FLUCUTE, FKPIDX,
     &            MISVAR, ORBVAR, NKP)

      CLOSE(9)

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
      IF (IOSTAT.EQ.0)  GO TO 150
C
  170 CONTINUE
C
      KPDXOLD = KPIDX
      IUTOLD  = IUT
      CYEAROLD = CYEAR
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
               maginput(1,1) = FKPIDX*10.
               if (FKPIDX.GE.9.) maginput(1,1) = 89.9D0
C
c      The IRBEM accepts KP 0-9 only, and YEAR 1965 - 2015
c       No longer the case since 2015, now it is > 1900
C
               READ(CYEAR,'(I4)')LYEAR
c               IF (LYEAR .GE. 2015) LYEAR = 2014
               IF (LYEAR .LT. 1955) LYEAR = 1955
C
c      remember this is for the map grid and currently this already be precalculated!
C
	           call make_lstar1(1,4,options,0,LYEAR,1,UT,ALT450,
     *              DLATD,DLOND,maginput, lm,lstar,blocal,bmin,xj,mlt)
	           FL = DABS(lm(1,1))
c
               IF (FL.GT.999.99)  FL = 999.99
               FLG(JCK,KCK) = FL
               KAL = KAL + 1
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
         FLG(LT,1) = FLG(LT,73) 
         RCG(LT,75) = RCG(LT,3) 
         FLG(LT,75) = FLG(LT,3)
  350 CONTINUE
C
C........+.........+.........+.........+.........+.........+.........+..
C     \/ Check for missing data
C        DLATD is Latitude  in degrees
C        DLOND is Longitude in degrees (East of Greenwich)
C........+.........+.........+.........+.........+.........+.........+..
C
      IFCNT = 0
C........+.........+.........+.........+.........+.........+.........+..
C     \/ Replace zero cutoff values with 1.0-6 to avoid possible "blowup" 
C........+.........+.........+.........+.........+.........+.........+..
C
      DO 410 J = 1, 37
         DO 400 K = 1, 75
            IF (RCG(J,K).EQ.0.0)  RCG(J,K) = 1.0E-06
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
  420    CONTINUE
  430 CONTINUE
C
C........+.........+.........+.........+.........+.........+.........+..
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
     &                   FLOND,FLATD,ALTKM,FKPIDX)
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
      COMMON /KYMDHM/  KYEAR,KMONTH,KDATE,KHOUR,KMIN,KSEC,KDOY,KUT

      SAVE                                                              ! AFRL

      IF (KPIDX.NE.KPDXOLD .OR. IUT.NE.IUTOLD 
     &     .OR. CYEAR.NE.CYEAROLD)  THEN
          CALL SETKP3UT (CYEAR,KPIDX,KPDXOLD, IUT,IUTOLD,FKPIDX)           ! SUB02
          KPDXOLD = KPIDX
          IUTOLD  = IUT
          CYEAROLD = CYEAR

      ENDIF
      CALL LINT5X5 (FLATD,FLOND,ALTKM,FKPIDX)                           ! SUB04

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
C      
      INTEGER*4    options(5)
      real*8       maginput(25,1), lm(1,1)
      data options/0,30,0,0,0/
c
c........+.........+.........+.........+.........+.........+.........+..
c
      COMMON /BANDL/   BBPT, FLPT, GLMDAR
      COMMON /CUTRIG/  RLGVPV,RUGVPV,RCGVPV, RLGVPE,RUGVPE,RCGVPE, 
     &                 RLGVPW,RUGVPW,RCGVPW, PTV
      COMMON /KYMDHM/  KYEAR,KMONTH,KDATE,KHOUR,KMIN,KSEC,KDOY,KUT


      COMMON /SETRCL/ VKRC(37,75), RCG(37,75), FLG(37,75)

      COMMON /DIAGC/  LTS(37),     LNS(75)
c
c........+.........+.........+.........+.........+.........+.........+..
c
      SAVE                                                              ! AFRL

c........+.........+.........+.........+.........+.........+.........+..
C     \/ First , check coordinates
c........+.........+.........+.........+.........+.........+.........+..
c
      IF (FLATD.GT.90.0 .OR. FLATD.LT.-90.0)  THEN
         WRITE (*, 4001)  FLATD

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
c       use the IRBEM lib
c       call make_lstar1(ntime,kext,options,sysaxes,iyear,idoy,ut, x1,x2,x3, maginput, lm,lstar,blocal,bmin,xj,mlt)
      maginput(1,1) = FKPIDX*10.
      if (FKPIDX.GE.9.) maginput(1,1) = 89.9D0
      LYEAR = KYEAR
c      IF (LYEAR .GE. 2015) LYEAR = 2014
      IF (LYEAR .LT. 1955) LYEAR = 1955
C
c   we have to use the map date/time rather than the given one
c   however it is debatable if one should use the map DOY = 1, or the given one KDOY
c   ideally we should have multiple sets of maps corresponding to different DOY and interpolate
c   between the sets of maps. At the moment we only have maps for DOY=1.
c   we will use DOY=1 for the grid inetrpolation and then recalculate the L value using the 
c   actual DOY when rescaling the cut-off for altitude. Thus taking the DOY effect into consideration.
c
      call make_lstar1(1,4,options,0,LYEAR,1,UT,ALT450,FLATD,
     &  FLOND,maginput, lm,lstar,blocal,bmin,xj,mlt)
      FLPT = DABS(lm(1,1)) 
c
      IF (FLPT.GT.999.99 .or. FLPT .lt. 0.)  FLPT = 999.99
      FLPSQ = FLPT*FLPT
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
c........+.........+.........+.........+.........+.........+.........+..
c     \/ Get VK constant for this point for altitude interpolation
c........+.........+.........+.........+.........+.........+.........+..
c
      VKRCP =  RCPT*FLPSQ

c........+.........+.........+.........+.........+.........+.........+..
c     >>>>>>>>>>   Interpolated cutoff values at 450 KM   <<<<<<<<<<
c........+.........+.........+.........+.........+.........+.........+..
      RCGVPV = RCPT
c
c........+.........+.........+.........+.........+.........+.........+..
c     \/ Altitude interpolation
c        Store vertical cutoff rigidity in COMMON CUTOFF ARRAYS
c........+.........+.........+.........+.........+.........+.........+..
c
      IF (ALTKM.NE.450.0 .OR. KDOY.NE.1)  THEN
        maginput(1,1) = KPIDX*10.
        if (KPIDX.GE.9) maginput(1,1) = 89.9D0
         LYEAR = KYEAR
c         IF (LYEAR .GE. 2015) LYEAR = 2014
         IF (LYEAR .LT. 1955) LYEAR = 1955   
	     call make_lstar1(1,4,options,0,LYEAR,KDOY,UT,ALTKM,FLATD,
     &                 FLOND,maginput,lm,lstar,blocal,bmin,xj,mlt)
	    FLKM = DABS(lm(1,1))
c
	    IF (FLKM.GT.999.99 .or. FLKM .lt. 0.)  FLKM = 999.99
        FLPT = FLKM
        FLKMSQ = FLKM*FLKM
        RCGVPV = VKRCP/FLKMSQ
      ENDIF

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
c
c       Changed from /11. to /14. by f.lei
c
      RADIST = (ALTKM + 6371.2)/6371.2                    ! DS 02 08 28
      RADISTSQ = RADIST*RADIST                            ! DS 02 08 28
      RCORT = (LOG(RADISTSQ))/14.0                        ! DS 02 08 28
      RCGVPV = RCGVPV - rcort                             ! DS 02 08 28
      If (RCGVPV.le.0.0)  RCGVPV = 0.0                    ! DS 02 08 28

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

      xIN(1) = 0.      ! on the surface
      xIN(2) = FLATD
      xIN(3) = FLOND
C
c     GEO->MAG 
C
      LYEAR = KYEAR
c      IF (LYEAR .GE. 2015) LYEAR = 2014
      IF (LYEAR .LT. 1955) LYEAR = 1955
      call coord_trans_vec1(1,0,6,LYEAR,1,
     &            UT,xIN,xOUT)
c     calculate the latitude     
      call CAR_SPH(xOUT,r,flati,flongi)
C
c     flati is the magnetic latitude in degrees. Note this is not the same as 
c     the corrected geomagnetic latitude, but the difference should be small
C
      GMLATCD = flati
      GMLATCR = GMLATCD/57.2957795

      GLMDAR = 0.0
      IF (FLPT.GT.1.0D0)  GLMDAR = DACOS((1.0D0/SQRT(FLPT)))
      GLMDAD = GLMDAR*57.2957795 

      IF (DABS(GMLATCR).LT.(DABS(GLMDAR)))  GLMDAR = DABS(GMLATCR)

c
      IF (RCGVPV.LE.0.0)  RCGVPV = 1.0E-06
c
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE MF_ORB(LUN, PRODEF, TITLE, MISVAR, ORBVAR, NRC, RC, NT,
     &                  IKP, NKP)

      IMPLICIT NONE

      INCLUDE 'interface.h'

      REAL*8 RC(*)
      INTEGER*4 LUN, NRC, NT, IKP, NKP, MF(8), I
      CHARACTER*(*) PRODEF, TITLE
      CHARACTER*80 SDUM

      TYPE(META_VAR) :: MISVAR, ORBVAR

      MF(3) = MISVAR%N + ORBVAR%N + 2
      MF(4) = 2
      MF(5) = 7
      MF(6) = 2*NRC + 5
      MF(7) = NT
      MF(8) = NKP - IKP + 1
      CALL MF_HEADER(LUN, MF, PRODEF, TITLE)
      WRITE (LUN,'(A)') '''MOD_ABB'', -1,''MSM'''
      CALL MF_ANNOT(LUN, MISVAR)
      CALL MF_ANNOT(LUN, ORBVAR)
      WRITE (SDUM,'(A,I2,A)') '(A,I3,1H,,', NRC, '(F7.2,1H,),A)'
      WRITE (LUN, TRIM(SDUM)) '''RIGIDITY'',', NRC, (RC(I),I=1,NRC),
     &  '''GV'''
      WRITE (LUN,'(A)') '''PS Annotation'', 0, 1'
      WRITE (LUN,'(A)') '''PS Annotation'', 0, 0'
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
      WRITE(LUN, '(A,I3,A)') '''MS_Err'','' '',', NRC,
     &  ',''Error in magnetic transmission factor for '',''RIGIDITY'''

      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE MF_TRI(LUN, PRODEF, TITLE, NRC, RC, FLUCUT,FLUCUTE,
     &                  FKP, MISVAR, ORBVAR, NKP)

      IMPLICIT NONE

      INCLUDE 'interface.h'

      REAL*8 RC(*), FLUCUT(NRC,100), FLUCUTE(NRC,100), FKP(*)
      INTEGER*4 LUN, NRC, IKP, NKP, MF(8) /0,0,0,2,0,0,0,0/, I, J
      CHARACTER*(*) PRODEF, TITLE
      CHARACTER*80 SDUM

      TYPE(META_VAR) :: MISVAR, ORBVAR

      MF(3) = MISVAR%N + ORBVAR%N + 2
      MF(5) = 3
      MF(6) = 1 + 2 * NKP
      MF(7) = NRC
      CALL MF_HEADER(LUN, MF, PRODEF, TITLE)
      WRITE (LUN,'(A)') '''MOD_ABB'', -1,''MSM'''
      CALL MF_ANNOT(LUN, MISVAR)
      CALL MF_ANNOT(LUN, ORBVAR)
      WRITE (SDUM,'(A,I3,A)') '(A,I3,', NKP, '(1H,,F8.4),A)'
      WRITE (LUN,SDUM) '''KP'',', NKP, (FKP(I),I=1,NKP), ', '' '''
      WRITE (LUN,'(A)') '''PS Annotation'', 0, 1'
      WRITE (LUN,'(A)') '''PS Annotation'', 0, 0'
      WRITE (LUN,'(A)') '''Rigidity'',''GV'', 1,''Rigidity'''
      WRITE(LUN, '(A,I3,A)') '''AMS'','' '',', NKP,
     &  ',''Average magnetic transmission factor for Kp='',''KP'''
      WRITE(LUN, '(A,I3,A)') '''AMS_err'','' '',', NKP, ',''Error of the
     &  average magnetic transmission factor for Kp='',''KP'''

      DO I=1,NRC
        WRITE (LUN,'(F7.2,1P,100(1H,,E12.4:),0P)') RC(I),
     &    (FLUCUT(I,J),J=1,NKP), (FLUCUTE(I,J),J=1,NKP)
      END DO
      WRITE (LUN,'(A)') '''End of File'''

      RETURN
      END

C----------------------------------------------------------------------
c     SUBROUTINE EXPFAC(RADIUS, MLAT, E, Z, A, F, N, ISTORM)
C     computes an exposure factor at B,L for particles Z,A 
C     arriving from interplanetary space with energies E(N)
C     by averaging over arrival directions. 
C     Exposure factor returned in F(N).
C     First implementation E.J. Daly ESA/ESTEC/WMA  9/87
c     Updated by F. Lei 03/14
c     1) change to rigidity based in stead of energy, so no need of Z and A
c     2) the Rcv is given as input, no need for RADIUS and ISTORM
c     3) use omcao2
c
      SUBROUTINE EXPFAC(RCV, MLAT, RC, F, N)

      IMPLICIT NONE

      REAL*8 RC(*), F(*), RCV, MLAT
C     table of One-Minus-Cos-Angles-Over-2 :
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
c
c        print *,RCV, mlat, cut(1), cut(5), cut(9)
c
        DO IRC = 1,N
           F(IRC) = 0.0
c  is there an value in the array above the lowest (Eastward, from West!) cutoff?:
           IF (RC(IRC) .GE. cut(1)) THEN
              F(IRC) = 1.0
              DO iangle=2,nangle
c     find the angular location where the cutoff goes over the rc:
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
