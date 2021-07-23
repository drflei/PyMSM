      PROGRAM example5
C
      INTEGER*4    kunit
      LOGICAL*4    list
C
C     open the namelist file data.txt
C
      kunit = 1
      OPEN( UNIT=kunit, FILE='ex5.nml', STATUS='old') 
C
C     initialize the library
C
      list  = .TRUE.
      CALL initul (kunit, list)
C
C     close the namelist file
C
      CLOSE (kunit)
C
      END
C     -----------------------------------------------------------------
      SUBROUTINE initul (kunit, list)
C
C     initialisation of UNILIB with a namelist
C
C          kunit : unit number of the namelist file
C          list  : when set to true, prints the namelist contents
C
      INCLUDE 'structure.h'
C
      INTEGER*4     kunit
      LOGICAL*4     list
C
C     common blocks UC150 and UC190
C
      COMMON /UC150/matm, ntspec, nnspec, kspec, kflag
C 
      TYPE(zatm) matm
      INTEGER*4     ntspec, nnspec, kspec(30)
      INTEGER*4     kflag(50)
C
      COMMON /UC190/ prop, stepx, stpmin, umsq, upsq, uk2, uk3, 
     :               epskm, epsrel, stplst, xclat, kmflg, nxstp       
C
      REAL*8        prop, stepx, stpmin
      REAL*8        umsq, upsq, uk2, uk3
      REAL*8        epskm, epsrel, stplst, xclat
      INTEGER*4     kmflg, nxstp
C
C     definition of the namelist 'unilib'
C
C       note that datemagn and dateatmo are integer and
C       max_kp is a float
C
      INTEGER*4     kinner, kouter, datemagn(6), katmo, dateatmo(6)
      REAL*8        bltime, dst, val_kp, winddens, windvel
      REAL*8        rzss, f107, f107a, max_kp, ap(6) 
C
      NAMELIST /unilib/ kinner, bltime, kouter, datemagn, val_kp, dst,
     :                  winddens, windvel, kmflg, katmo, dateatmo,
     :                  rzss, f107, f107a, max_kp, ap, kflag
C
C     other variables
C
      INTEGER*4     kout, kini, ifail
      CHARACTER*32  label
      TYPE(zdat) mdm, mda
      REAL*8        param(10)
C
C     initialize UNILIB
C
      kout        =    6
      kini        =    1
C
      CALL UT990 (kout, kini, ifail)
      IF( ifail .LT. 0 )THEN
        WRITE(kout,1000) 'UT990',ifail
        STOP
      ENDIF
C
C     set the default values
C
C       note that kmflg and kflag values are set by
C       the subroutine UT990
C
      kinner      =    0
      bltime      = 1995.0
      kouter      =    0
      datemagn(1) = 1995
      datemagn(2) =    1
      datemagn(3) =    1
      datemagn(4) =    0
      datemagn(5) =    0
      datemagn(6) =    0
      val_kp      =    0.0
      dst         =  -30.0
      winddens    =   25.0
      windvel     =  300.0
      katmo       =    1
      dateatmo(1) = 1995
      dateatmo(2) =    1
      dateatmo(3) =    1
      dateatmo(4) =    0
      dateatmo(5) =    0
      dateatmo(6) =    0
      rzss        =   80.0
      f107        =   70.0
      f107a       =   50.0
      max_kp      =    3.0
      ap(1)       =   20.0
      ap(2)       =   20.0
      ap(3)       =   20.0
      ap(4)       =   20.0
      ap(5)       =   20.0
      ap(6)       =   20.0
C
C     read the namelist
C
      READ(kunit,unilib,IOSTAT=ifail)
      IF( ifail.ne.0 )THEN
        WRITE(kout,1010) ifail
        STOP
      ENDIF
C
C     transform the date
C
      mdm%iyear   = datemagn(1)
      mdm%imonth  = datemagn(2)
      mdm%iday    = datemagn(3)
      mdm%ihour   = datemagn(4)
      mdm%imin    = datemagn(5)
      mdm%secs    = datemagn(6)
      CALL UT540 (mdm)
C
      mda%iyear   = dateatmo(1)
      mda%imonth  = dateatmo(2)
      mda%iday    = dateatmo(3)
      mda%ihour   = dateatmo(4)
      mda%imin    = dateatmo(5)
      mda%secs    = dateatmo(6)
      CALL UT540 (mda)
C
C     set the different models
C
      CALL UM510 (kinner, bltime, label, kout, ifail)
      IF( ifail .LT. 0 )THEN
        WRITE(kout,1000) 'UM510',ifail
        STOP
      ENDIF
C
      param(1)    = val_kp
      param(2)    = dst
      param(3)    = 0.0d0
      param(4)    = winddens
      param(5)    = windvel
      param(6)    = 0.0d0
      param(7)    = 0.0d0
      param(8)    = 0.0d0
      param(9)    = 0.0d0
      param(10)   = 0.0d0
C
      CALL UM520 (kouter, mdm%amjd, param,
     :           label, kout, ifail)
      IF( ifail .LT. 0 )THEN
        WRITE(kout,1000) 'UM520',ifail
        STOP
      ENDIF
C
      CALL UA610 (katmo, mda, rzss, f107a, f107, max_kp,
     :            ap(1), ap(2), ap(3), ap(4), ap(5), ap(6),
     :            label, kout, ifail)
      IF( ifail .LT. 0 )THEN
        WRITE(kout,1000) 'UA610',ifail
        STOP
      ENDIF
C
C     print the namelist
C
      IF( list )then
        WRITE(kout,1020)
        WRITE(kout,unilib)
        WRITE(kout,1020)
      ENDIF
C
C     format
C
 1000 FORMAT(' *** Error in subroutine ',a5,', IFAIL =',i7,' ***')
 1010 FORMAT(' *** Error in the namelist, IOSTAT =',i7,' ***')
 1020 FORMAT(25('-'))
C
      END

