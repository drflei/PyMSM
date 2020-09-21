      SUBROUTINE FILPAR(PRONAM)      
C
C Parses the command line for a project name. If no command line
C arguments are given, it prompts for input.
C
      IMPLICIT NONE

      INCLUDE 'interface.h'

      INTEGER*4 VALUES(8), I, I1, I2
      INTEGER NARG /1/
      CHARACTER(*) PRONAM
      CHARACTER(80) ARG
      CHARACTER*3 MONTH(12) /'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
     &                       'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'/

      WRITE (ARG,'(80(1H ))')
      CALL GETARG(NARG, ARG)
      CALL STR_TRIM(ARG, I1, I2)
      DO WHILE (ARG(I1:I1) .EQ. ' ')
        WRITE (6,*) ' Project name ?'
        READ (5,'(A80)') ARG
        CALL STR_TRIM(ARG, I1, I2)
      END DO
      PROJECT = ARG(I1:I2)
      WRITE (NAMPRO,'(80(1H ))')
      NAMPRO = PRONAM

      CALL DATE_AND_TIME(VALUES=VALUES)
      WRITE (SDATE,'(I2,''-'',A,''-'',I4)') VALUES(3), MONTH(VALUES(2)),
     &                                      VALUES(1)
      WRITE (STIME,'(2(I2.2,'':''),I2.2)') (VALUES(I),I=5,7)

      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE FILE_OPEN(LUN, TYPE, EXT, FORM, RECL, STATUS, ACCESS,
     &                     IOSTAT)

      IMPLICIT NONE

      INTERFACE
        SUBROUTINE ERR_WRITE(ERRSTR, FCLOSE)
          CHARACTER(*) ERRSTR
          LOGICAL, OPTIONAL :: FCLOSE
        END SUBROUTINE
      END INTERFACE

      INTEGER*4 LUN, ISTAT, ISTATA, I1, I2, I3, I4
      INTEGER*4, OPTIONAL :: RECL, IOSTAT
      CHARACTER(*) TYPE, EXT
      CHARACTER(*), OPTIONAL :: FORM, STATUS, ACCESS
      CHARACTER*80 FMAT, FSTAT, FACC, FILNAM
      LOGICAL FILOPN

      CHARACTER*80 PROJECT, NAMPRO
      CHARACTER*20 VERSTR
      CHARACTER*11 SDATE
      CHARACTER*8 STIME

      COMMON /INFO/ PROJECT, NAMPRO, VERSTR, SDATE, STIME

      IF (PRESENT(FORM)) THEN
        FMAT = FORM
      ELSE
        FMAT = 'FORMATTED'
      END IF
      IF (PRESENT(STATUS)) THEN
        FSTAT = STATUS
      ELSE
        FSTAT = 'UNKNOWN'
      END IF
      IF (PRESENT(ACCESS)) THEN
        FACC = ACCESS
      ELSE
        FACC = 'SEQUENTIAL'
      END IF
      IF (PRESENT(IOSTAT)) THEN
        ISTATA = IOSTAT
      ELSE
        ISTATA = 0
      END IF
      ISTAT = 0
      CALL STR_TRIM(TYPE, I3, I4)
      CALL STR_TRIM(PROJECT, I1, I2)
      IF (TYPE(I3:I3) .EQ. ' ') THEN
        FILNAM = PROJECT(I1:I2)//'.'//EXT
      ELSE
        FILNAM = PROJECT(I1:I2)//'_'//TYPE(I3:I4)//'.'//EXT
      END IF
      CALL STR_TRIM(FILNAM, I1, I2)
      INQUIRE (FILE=FILNAM(I1:I2), OPENED=FILOPN)
      IF (.NOT. FILOPN) THEN
        IF (PRESENT(RECL)) THEN
          OPEN (LUN, FILE=FILNAM(I1:I2), FORM=FMAT, RECL=RECL,
     &          STATUS=FSTAT, ACCESS=FACC, IOSTAT=ISTAT)
        ELSE
          OPEN (LUN, FILE=FILNAM(I1:I2), FORM=FMAT,
     &          STATUS=FSTAT, ACCESS=FACC, IOSTAT=ISTAT)
        END IF
        IF ((ISTAT .NE. 0) .AND. (ISTATA .EQ. 0))
     &    CALL ERR_WRITE('Error opening file '//FILNAM(I1:I2))
      END IF
      IF (PRESENT(IOSTAT)) IOSTAT = ISTAT

      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ERR_WRITE(ERRSTR, FCLOSE)

      IMPLICIT NONE

      INTERFACE
        SUBROUTINE FILE_OPEN(LUN, TYPE, EXT, FORM, RECL, STATUS, ACCESS,
     &                       IOSTAT)
          INTEGER*4 LUN
          INTEGER*4, OPTIONAL :: RECL, IOSTAT
          CHARACTER(*) TYPE, EXT
          CHARACTER(*), OPTIONAL :: FORM, STATUS, ACCESS
        END SUBROUTINE
      END INTERFACE

      INTEGER*4 ISTAT, I1, I2
      CHARACTER(*) ERRSTR
      CHARACTER(80) FILNAM
      LOGICAL FILOPN, FILCLO
      LOGICAL, OPTIONAL :: FCLOSE

      CHARACTER*80 PROJECT, NAMPRO
      CHARACTER*20 VERSTR
      CHARACTER*11 SDATE
      CHARACTER*8 STIME

      COMMON /INFO/ PROJECT, NAMPRO, VERSTR, SDATE, STIME

      IF (PRESENT(FCLOSE)) THEN
        FILCLO = FCLOSE
      ELSE
        FILCLO = .TRUE.
      END IF
      ISTAT = 0
      CALL STR_TRIM(PROJECT, I1, I2)
      FILNAM = PROJECT(I1:I2)//'.err'
      CALL STR_TRIM(FILNAM, I1, I2)
      INQUIRE (FILE=FILNAM(I1:I2), OPENED=FILOPN, IOSTAT=ISTAT)
      IF ((.NOT. FILOPN) .OR. (ISTAT .NE. 0))
     &  OPEN (62, FILE=FILNAM(I1:I2), STATUS='UNKNOWN', IOSTAT=ISTAT)
      IF (ISTAT .NE. 0) THEN
        WRITE (*,*) 'Error opening error file '//FILNAM(I1:I2)
      ELSE
        WRITE (62,'(A)') ERRSTR
        IF (.NOT. FILCLO) RETURN
      END IF
      STOP ' '

      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE WAR_WRITE(WARSTR)

      IMPLICIT NONE

      INCLUDE 'interface.h'

      INTEGER*4 LUN /64/
      CHARACTER(*) WARSTR

      CALL FILE_OPEN(LUN, ' ', 'war')
      WRITE (LUN,'(A)') WARSTR

      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE UEXIT
  
      STOP

      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE MV_STORE(I, NT, LTXT, IFAIL)
C
C I: index of header line
C NT: total number of header lines
C LTXT: content of header line
C IFAIL: 1 (default) --> ok
C        0 --> stop calling MV_STORE
C       -1 --> error
C
      IMPLICIT NONE

      INCLUDE 'interface.h'

      INTEGER*4 I, NT, IFAIL
      CHARACTER(*) LTXT

      TYPE(META_VAR) :: MVARS

      COMMON /METAVAR/ MVARS

      IF (I .GT. NT) THEN
        IFAIL = -1
      ELSE
        MVARS%TXT(I) = LTXT
        MVARS%N = NT
        IFAIL = 1
      END IF

      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE MV_READ(NAME, VAR)

      IMPLICIT NONE

      INCLUDE 'interface.h'

      INTEGER*4 I, J, I1, I2
      CHARACTER(*) NAME

      TYPE(META_VAR) :: MVARS, VAR

      COMMON /METAVAR/ MVARS

      CALL STR_TRIM(NAME, I1, I2)
      
      IF (NAME(I1:I1) .EQ. ' ') THEN
        VAR%N = MVARS%N
        DO I=1,MVARS%N
          VAR%TXT(I) = MVARS%TXT(I)
        END DO
      ELSE
        VAR%N = 0
        DO I=1,MVARS%N
c          J = INDEX(''//MVARS%TXT(I), NAME)
          J = INDEX(MVARS%TXT(I), NAME)
          IF ((J .GT. 0) .AND. (J .LE. 2)) THEN
            VAR%N = VAR%N + 1
            VAR%TXT(VAR%N) = MVARS%TXT(I)
          END IF
        END DO
      END IF
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE PS_STORE(I, NT, LTXT, IFAIL)
C
C I: index of header line
C NT: total number of header lines
C LTXT: content of header line
C IFAIL: 1 (default) --> ok
C        0 --> stop calling PS_STORE
C       -1 --> error
C
      IMPLICIT NONE

      INTEGER*4 I, NT, IFAIL, N
      CHARACTER*80 TXT(1000)
      CHARACTER(*) LTXT

      COMMON /PSHEAD/ TXT, N

      IF (I .GT. NT) THEN
        IFAIL = -1
      ELSE IF ( LEN_TRIM(LTXT) .GT. 80 ) THEN
C MK 20100519:
C            : LEN_TRIM returns the length of a character string, ignoring any trailing blanks
C            :
        IFAIL = -1
      ELSE
        TXT(I) = LTXT
        N = NT
        IFAIL = 1
      END IF

      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE PS_READ(PS, NB)

      IMPLICIT NONE

      INCLUDE 'interface.h'

      INTEGER*4 PSN, NB, I, J, IL, IB, JB, NL
      CHARACTER*80 PSTXT(1000)

      TYPE(PSANN) :: PS

      COMMON /PSHEAD/ PSTXT, PSN

      IL = 1
      DO JB=1,NB
        I = INDEX(PSTXT(IL), '''PS Annotation''')
        IF (I .LT. 1)
     &    CALL ERR_WRITE(' Error in PS annotation: header label')
        I = INDEX(PSTXT(IL)(I+1:), '''')
        READ (PSTXT(IL)(I+3:),*) NL, IB
        IF (MOD(NL, 2) .NE. 0)
     &    CALL ERR_WRITE(' Error in PS annotation: odd number of lines')
        IF (IL+NL .GT. PSN)
     &    CALL ERR_WRITE(' Error in PS annotation: number of lines')
        IF (NB-JB .GT. IB) THEN
          CALL ERR_WRITE(' Error in PS annotation: not so many blocks')
        END IF
        IF (JB .EQ. NB) THEN
          PS%N = NL / 2
          DO I=1,PS%N
            PS%TXT(I) = PSTXT(2*I-1+IL)
            J = INDEX(PS%TXT(I), '''')
            IF (J .LT. 1)
     &        CALL ERR_WRITE(' Error in PS annotation: label text')
            PS%TXT(I) = PS%TXT(I)(J+1:)
            J = INDEX(PS%TXT(I), '''')
            IF (J .LT. 1)
     &        CALL ERR_WRITE(' Error in PS annotation: label text')
            PS%TXT(I) = PS%TXT(I)(1:J-1)
            READ (PSTXT(2*I+IL),*) (PS%POS(I,J),J=1,3)
          END DO
          RETURN
        END IF
        IL = IL + NL + 1
      END DO

      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE HEADER_HTML(LUN, PRODEF, TITLE)

      IMPLICIT NONE

      INCLUDE 'interface.h'

      INTEGER*4 LUN, I1, I2
      CHARACTER(*) PRODEF, TITLE

      CALL STR_TRIM(NAMPRO, I1, I2)
      WRITE (LUN,'(A,//,A,/,3A,/,A,//,A,/)') '<HTML>', '<HEAD>',
     &  '<TITLE>', NAMPRO(I1:I2), '</TITLE>', '</HEAD>',
     &  '<BODY bgcolor="white" link="blue" vlink="blue" alink="blue">'

      WRITE (LUN,'(4(A,/),2A,6(''&nbsp;''),A)')
     &  '<P>', '<CENTER>', '<TABLE border="0" cellspacing="2">', '<TR>',
     &  '<TD align="left"><B>SPENVIS&nbsp;', VERSTR, '</B></TD>'

      WRITE (LUN,'(2A,3(''&nbsp;''),2A,/,A)') 
     &  '<TD align="right"><B>', SDATE, STIME, '</B></TD>', '</TR>'

      WRITE (LUN,'(3A)') '<TR><TD colspan="2" align="center"><B>',
     &  NAMPRO(I1:I2), '</B></TD></TR>'

      CALL STR_TRIM(PRODEF, I1, I2)
      WRITE (LUN,'(3A)')
     &  '<TR><TD colspan="2" align="center"><B>Project: ',
     &  PRODEF(I1:I2), '</B></TD></TR>'

      CALL STR_TRIM(TITLE, I1, I2)
      WRITE (LUN,'(2A,4(A,/))')
     &  '<TR><TD colspan="2" align="center"><B>', TITLE(I1:I2),
     &  '</B></TD></TR>', '</TABLE>', '</CENTER>', '<P>'

      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE CLOSE_HTML(LUN)

      IMPLICIT NONE

      INTEGER*4 LUN

      WRITE (LUN,'(/,A,//,A)') '</BODY>', '</HTML>'

      CLOSE (LUN)

      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE MF_HEADER(LUN, MF, PRODEF, TITLE)
C
C  MF: INTEGER*4(8): header information for UNIRAD file format
C                    MF(1): output: set to MF(2) + ... + MF(5) + 1
C                    MF(2): output: set to 1
C                    MF(3): input: number of user metavariables
C                    MF(4-8): inputs
C
      IMPLICIT NONE

      INCLUDE 'interface.h'

      INTEGER*4 LUN, MF(8), I1, I2
      CHARACTER(*) PRODEF, TITLE

      MF(2) = 1
      MF(3) = MF(3) + 2
      MF(1) = MF(2) + MF(3) + MF(4) + MF(5) + 1
      WRITE (LUN,'(A,8('','',I6))') '''*''', MF
      WRITE (LUN,'(4A,1X,2A)') '''SPENVIS ', VERSTR, ' - ', SDATE,
     &  STIME, ''''
      CALL STR_TRIM(PRODEF, I1, I2)
      WRITE (LUN,'(3A)') '''PRJ_DEF'', -1,''', PRODEF(I1:I2), ''''
      CALL STR_TRIM(TITLE, I1, I2)
      WRITE (LUN,'(3A)') '''PRJ_HDR'', -1,''', TITLE(I1:I2), ''''
      MF(3) = MF(3) - 2

      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE MF_ANNOT(LUN, VAR)

      IMPLICIT NONE

      INCLUDE 'interface.h'

      INTEGER*4 LUN, I, I1, I2

      TYPE(META_VAR) :: VAR

      DO I=1,VAR%N
        CALL STR_TRIM(VAR%TXT(I), I1, I2)
        WRITE (LUN,'(A)') VAR%TXT(I)(I1:I2)
      END DO

      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE PS_ANNOT(LUN, PS, NB)

      IMPLICIT NONE

      INCLUDE 'interface.h'

      INTEGER*4 LUN, NB, I, J, I1, I2

      TYPE(PSANN) :: PS

      WRITE (LUN,'(A,2('','',I3))') '''PS Annotation''', 2*PS%N, NB
      DO I=1,PS%N
        CALL STR_TRIM(PS%TXT(I), I1, I2)
        WRITE (LUN,'(3A)') '''', PS%TXT(I)(I1:I2), ''''
        WRITE (LUN,'(2(F5.2,'',''),F5.2)') (PS%POS(I,J),J=1,3)
      END DO

      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE STR_TRIM(STR, I1, I2)

      IMPLICIT NONE

      CHARACTER(*) STR
      INTEGER*4 I1, I2

      DO I1=1,LEN(STR)
        IF (STR(I1:I1) .NE. ' ') GO TO 1
      END DO
      I1 = 1
   1  DO I2=LEN(STR),I1+1,-1
        IF (STR(I2:I2) .NE. ' ') GO TO 2
      END DO
      I2 = I1
   2  RETURN

      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE STR_REPL_QUOTE(STR, I1, I2)

      IMPLICIT NONE

      CHARACTER(*) STR
      INTEGER*4 I1, I2, I, IN, IX

      IN=I1
      IX=I2
      IF ( IN .LT. 1 ) IN = 1
      IF ( IX .GT. LEN(STR) ) IX = LEN(STR)
      IF ( IX .LT. 1 ) IX = LEN(STR) + IX
      IF ( IN .GT. IX ) RETURN
      DO I=IN,IX
       IF ( STR(I:I) .EQ. '''' .OR.
     &       STR(I:I) .EQ. '"'
     &     ) STR(I:I) = '+'
      END DO
      END
C
C-----------------------------------------------------------------------
C
      LOGICAL FUNCTION STR_EQ(STR1, STR2)
C Compares the first string with the second, ignoring
C any white space at the end of each string.

      IMPLICIT NONE

      INTEGER*4 I1, I2, J1, J2
      CHARACTER(*) STR1, STR2

      CALL STR_TRIM(STR1, I1, I2)
      CALL STR_TRIM(STR2, J1, J2)
      STR_EQ = (STR1(I1:I2) .EQ. STR2(J1:J2))

      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE SKIP_BLOCK(LUN, KSKIP)

      IMPLICIT NONE

      INCLUDE 'interface.h'

      INTEGER*4 LUN, KSKIP, MF(8), I, J
      CHARACTER*3 A

      DO I=1,KSKIP
        READ (LUN, *, ERR=9000) A, MF
        DO J=1,MF(1)-1
          READ (LUN, *, ERR=9000) A
        END DO

        A = '   '
        DO WHILE (A(1:1) .NE. '''')
          READ (LUN, '(A)', ERR=9000) A
        END DO
      END DO

      RETURN

 9000 CALL ERR_WRITE(' Error skipping data block')

      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE CALEND(DJ, IDIR, IY, IM, ID, IH, IMIN, SECS)

C  VERSION Aug 87 based on routines by H. Klinkrad/ESTEC/WMM
C                 coded E.J. Daly

C  IDIR=+1 CALCULATE DATE AND HOUR FROM MODIFIED JULIAN DAY
C  IDIR=-1 CALCULATE MODIFIED JULIAN DAY FROM DATE AND HOUR
C  DJ =MODIFIED JULIAN DAY   (JULIAN DAY - 2433282.5)
C  IY=YEAR, IM=MONTH,    ID=DAY
C  IH=HOUR, IMI=MINUTES, DS=SECONDS
C
C  modified Julian Day referenced to 0hrs 1-jan-1950 (1950 jan 1.0)
C  which has JD number 2433282.5.

      IMPLICIT NONE

      REAL*8 DJ, SECS
      INTEGER*4 IDIR, IY, IM, ID, IH, IMIN

      IF (IDIR .EQ. -1) THEN
        IF (IY .LT. 1900) IY = IY + 1900
        CALL DATOJD(IY, IM, ID, IH, IMIN, SECS, DJ)
      ELSE
        CALL JDTODA(DJ, IY, IM, ID, IH, IMIN, SECS)
        IF (IY .GT. 1900) IY = IY - 1900
      END IF

      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DATOJD(IY, IM, ID, IH, IMIN, SECS, DJ)

C Compute Modified Julian Day from Date
C Method from H. Klinkrad

      IMPLICIT NONE

      REAL*8 SECS, DJ
      INTEGER*4 IY, IM, ID, IH, IMIN, JJ, L, JD

      JJ = (14-IM) / 12
      L = IY - 1900 - JJ
      JD = ID - 18234 + (1461*L) / 4 + (367*(IM-2+JJ*12)) / 12
      DJ = JD + ((IH*60+IMIN)*60+SECS) / 86400.0D0

      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE JDTODA(DJ, IY, IM, ID, IH, IMIN, SECS)

C Compute Date from Modified Julian Day
C Method from H. Klinkrad

      IMPLICIT NONE

      REAL*8 DJ, SECS, DDAY
      INTEGER*4 IY, IM, ID, IH, IMIN, JDAY, L, M, N, JJ

      DDAY = DJ + 1.0D-10
      JDAY = DDAY
      L = (4000*(JDAY+18204)) / 1461001
      N = JDAY - (1461*L) / 4 + 18234
      M = (80*N) / 2447

      ID = N - (2447*M) / 80
      JJ = M / 11
      IM = M + 2 - 12 * JJ
      IY = 1900 + L + JJ
      SECS = (DDAY-JDAY) * 24.0D0
      IH = SECS
      SECS = (SECS-IH) * 60.0D0
      IMIN = SECS
      SECS = (SECS-IMIN) * 60.0D0

      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE READ_ORB(PRODEF, MISTIT, NTRAJ, TMSTART, TMEND, TMDUR,
     &                    ORBTIT, YEAR, MON, DAY, HR, MN, SEC, TSTART,
     &                    TSEND, LUN)

      IMPLICIT NONE

      INCLUDE 'interface.h'

      REAL*8 TMSTART, TMEND, TMDUR, TSTART, TSEND, SEC
      REAL*4 SSEC
      INTEGER*4 NTRAJ, YEAR, MON, DAY, HR, MN, NBLOCK, LUN
      INTEGER*4 HEAD(2000), HSIZE /1000/, KSKIP /0/, MSIZE, I
      CHARACTER*80 PRODEF, MISTIT, ORBTIT, SDUM

      TYPE(META_VAR) :: MISVAR, ORBVAR, METAVAR

      EXTERNAL MV_STORE, PS_STORE

      COMMON /OHEAD/ HEAD

      CALL READ_HEADER_PLUS(LUN, HEAD, HSIZE, KSKIP, NBLOCK, MV_STORE,
     &                      PS_STORE)
      MSIZE = -1
      CALL MV_READ('PRJ_DEF', METAVAR)
      IF ((METAVAR%N) .GT. 0) THEN
        CALL GET_META_STRING(HEAD, 'PRJ_DEF', PRODEF, MSIZE)
      ELSE
        WRITE (PRODEF,'(80(1H ))')
      END IF
      CALL GET_META_STRING(HEAD, 'PRJ_HDR', MISTIT, MSIZE)
      CALL GET_META_STRING(HEAD, 'ORB_HDR', ORBTIT, MSIZE)
      MSIZE = 1
      CALL GET_META_INT(HEAD, 'MIS_NTR', NTRAJ, MSIZE)
      CALL GET_META_INT(HEAD, 'ORB_YEA', YEAR, MSIZE)
      CALL GET_META_INT(HEAD, 'ORB_MON', MON, MSIZE)
      CALL GET_META_INT(HEAD, 'ORB_DAY', DAY, MSIZE)
      CALL GET_META_INT(HEAD, 'ORB_HOU', HR, MSIZE)
      CALL GET_META_INT(HEAD, 'ORB_MIN', MN, MSIZE)
      CALL GET_META(HEAD, 'ORB_SEC', SSEC, MSIZE)
      SEC = SSEC

      CALL MV_READ('MIS_STA', MISVAR)
      READ (MISVAR%TXT(1),*) SDUM, I, TMSTART
      CALL MV_READ('MIS_END', MISVAR)
      READ (MISVAR%TXT(1),*) SDUM, I, TMEND
      CALL MV_READ('MIS_DUR', MISVAR)
      READ (MISVAR%TXT(1),*) SDUM, I, TMDUR

      CALL MV_READ('ORB_MJD', ORBVAR)
      READ (ORBVAR%TXT(1),*) SDUM, I, TSTART
      CALL MV_READ('ORB_TSE', ORBVAR)
      READ (ORBVAR%TXT(1),*) SDUM, I, TSEND
      
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE READ_EPH(MJD, ALT, LON, LAT, ALPHA, FOOTER, KSKIP)

      IMPLICIT NONE

      INTEGER*4 HEAD(2000), KOLUMN(100), KSIZE, KSKIP
      INTEGER*4 NCOL, NUMBER
      PARAMETER (NCOL=5)
      REAL*8 COORD(NCOL), MJD, ALT, LON, LAT, ALPHA
      CHARACTER*80 TCOORD(NCOL), UNITS(NCOL), CTITLE(NCOL), FOOTER

      COMMON /OHEAD/ HEAD

      DATA TCOORD /'MJD', 'Longitude', 'Latitude', 'Altitude', 'SPS'/

      KSIZE = 100
      CALL SELECT_VAR(HEAD, TCOORD, NCOL, UNITS, CTITLE, KOLUMN, KSIZE)

      NUMBER = 1
      CALL GET_VAR(HEAD, KOLUMN, KSIZE, KSKIP, NUMBER, COORD, NCOL,
     &             FOOTER)
      MJD = COORD(1)
      LON = COORD(2)
      LAT = COORD(3)
      ALT = COORD(4)
      ALPHA = COORD(5)

      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE READ_POSH(PRODEF, TITLE, NP, MAP, LBL, LUN)

      IMPLICIT NONE

      INCLUDE 'interface.h'

      INTEGER*4 NP, MAP, LUN
      INTEGER*4 HEAD(1000), HSIZE /1000/, KSKIP /0/, MSIZE, NBLOCK /0/
      CHARACTER*80 PRODEF, TITLE, LBL

      TYPE(META_VAR) :: METAVAR

      EXTERNAL MV_STORE, PS_STORE

      COMMON /GHEAD/ HEAD

      CALL READ_HEADER_PLUS(LUN, HEAD, HSIZE, KSKIP, NBLOCK, MV_STORE,
     &                      PS_STORE)
      MSIZE = -1
      CALL MV_READ('PRJ_DEF', METAVAR)
      IF ((METAVAR%N) .GT. 0) THEN
        CALL GET_META_STRING(HEAD, 'PRJ_DEF', PRODEF, MSIZE)
      ELSE
        WRITE (PRODEF,'(80(1H ))')
      END IF
      CALL GET_META_STRING(HEAD, 'PRJ_HDR', TITLE, MSIZE)
      CALL GET_META_STRING(HEAD, 'GRD_LBL', LBL, MSIZE)
      MSIZE = 1
      CALL GET_META_INT(HEAD, 'GRD_TYP', MAP, MSIZE)
      NP = HEAD(2)

      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE READ_POSD(ALT, LAT, LON, LT, UT)

      IMPLICIT NONE

      INTEGER*4 HEAD(1000), KOLUMN(100), KSIZE, KSKIP /0/
      INTEGER*4 NCOL, NUMBER
      PARAMETER (NCOL=5)
      REAL*8 COORD(NCOL), ALT, LAT, LON, LT, UT
      CHARACTER*80 TCOORD(NCOL), UNITS(NCOL), CTITLE(NCOL), FOOTER

      COMMON /GHEAD/ HEAD

      DATA TCOORD /'Altitude', 'Latitude', 'Longitude', 'LT', 'UT'/

      KSIZE = 100
      CALL SELECT_VAR(HEAD, TCOORD, NCOL, UNITS, CTITLE, KOLUMN, KSIZE)

      NUMBER = 1
      CALL GET_VAR(HEAD, KOLUMN, KSIZE, KSKIP, NUMBER, COORD, NCOL,
     &             FOOTER)
      ALT = COORD(1)
      LAT = COORD(2)
      LON = COORD(3)
      LT = COORD(4)
      UT = COORD(5)

      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE READ_POSD2(ALT, LAT, LON)

      IMPLICIT NONE

      INTEGER*4 HEAD(1000), KOLUMN(100), KSIZE, KSKIP /0/
      INTEGER*4 NCOL, NUMBER
      PARAMETER (NCOL=3)
      REAL*8 COORD(NCOL), ALT, LAT, LON
      CHARACTER*80 TCOORD(NCOL), UNITS(NCOL), CTITLE(NCOL), FOOTER

      COMMON /GHEAD/ HEAD

      DATA TCOORD /'Altitude', 'Latitude', 'Longitude'/

      KSIZE = 100
      CALL SELECT_VAR(HEAD, TCOORD, NCOL, UNITS, CTITLE, KOLUMN, KSIZE)

      NUMBER = 1
      CALL GET_VAR(HEAD, KOLUMN, KSIZE, KSKIP, NUMBER, COORD, NCOL,
     &             FOOTER)
      ALT = COORD(1)
      LAT = COORD(2)
      LON = COORD(3)

      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE READ_SPEC(LUN, SPEC, PRODEF, TITLE, KPLANET, 
     &                     NTRAJ, SEGDUR, MISDUR)

      IMPLICIT NONE

      INCLUDE 'interface.h'

      INTEGER*4 LUN, NTRAJ, OMNI, I, NCOL, KOLUMN(10), NUMBER /-1/
      INTEGER*4 KPLANET
      INTEGER*4 HEAD(2000), HSIZE /1000/, KSKIP /0/, KSIZE /10/, MSIZE
      INTEGER*4 NBLOCK
      PARAMETER (NCOL=3)
      REAL*8 SEGDUR, MISDUR, VALUES(NCOL, N_Energy+1)
      REAL*4 SMISDUR, ORBSTART, ORBEND
      CHARACTER*80 TVALUE(NCOL), UNITS(NCOL), CTITLE(NCOL), FOOTER
      CHARACTER*80 PRODEF, TITLE
      CHARACTER*3 MODABB

      TYPE(SPECTRUM) :: SPEC

      TYPE(META_VAR) :: METAVAR

      EXTERNAL MV_STORE, PS_STORE

      DATA TVALUE /'Energy', 'IFlux', 'DFlux'/

      CALL READ_HEADER_PLUS(LUN, HEAD, HSIZE, KSKIP, NBLOCK, MV_STORE,
     &                      PS_STORE)
      MSIZE = 1
      CALL GET_META_INT(HEAD, 'MIS_PLA', KPLANET, MSIZE)
      MSIZE = -1
      CALL GET_META_STRING(HEAD, 'MOD_ABB', MODABB, MSIZE)
      IF (MODABB .EQ. 'LET')
     &  CALL GET_META_STRING(HEAD, MODABB//'_ABS', TVALUE(1), MSIZE)
      CALL MV_READ('PRJ_DEF', METAVAR)
      IF ((METAVAR%N) .GT. 0) THEN
        CALL GET_META_STRING(HEAD, 'PRJ_DEF', PRODEF, MSIZE)
      ELSE
        WRITE (PRODEF,'(80(1H ))')
      END IF
      CALL GET_META_STRING(HEAD, 'PLT_HDR', TITLE, MSIZE)
      WRITE (SPEC%LABEL,'(80(1H ))')
      SPEC%LABEL = TITLE
      WRITE (TITLE,'(80(1H ))')
      CALL GET_META_STRING(HEAD, 'PRJ_HDR', TITLE, MSIZE)
      MSIZE = 1
      CALL GET_META_INT(HEAD, 'MIS_NTR', NTRAJ, MSIZE)
      CALL GET_META_INT(HEAD, MODABB//'_OMN', OMNI, MSIZE)
      SPEC%OMNI = OMNI
      IF (MODABB .NE. 'SEP') THEN
        CALL GET_META(HEAD, 'ORB_MJD', ORBSTART, MSIZE)
        IF (MSIZE .GT. 0) THEN
          CALL GET_META(HEAD, 'ORB_TSE', ORBEND, MSIZE)
          SEGDUR = ORBEND - ORBSTART
        ELSE
          SEGDUR = 0.0D0
        END IF
      ELSE
        SEGDUR = 0.0D0
      END IF
      
      MSIZE = 1
      CALL GET_META(HEAD, 'MIS_DUR', SMISDUR, MSIZE)
      MISDUR = SMISDUR
      CALL SELECT_VAR(HEAD, TVALUE, NCOL, UNITS, CTITLE, KOLUMN, KSIZE)
      NUMBER = N_Energy + 1
      CALL GET_VAR(HEAD, KOLUMN, KSIZE, KSKIP, NUMBER, VALUES, NCOL,
     &             FOOTER)
      IF (FOOTER(1:3) .NE. 'End')
     &  CALL ERR_WRITE(' Error reading spectrum file')

      NUMBER = MIN(N_Energy, NUMBER)
      SPEC%N = 0
      DO I=1,NUMBER
c        IF (VALUES(2,I) .GT. 0.0D0) THEN
          SPEC%N = SPEC%N + 1
          SPEC%ENERGY(SPEC%N) = VALUES(1,I)
          SPEC%IFLUX(SPEC%N) = MAX(VALUES(2,I), 0.0D0)
          SPEC%DFLUX(SPEC%N) = MAX(VALUES(3,I), 0.0D0)
c        END IF
      END DO
        
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE READ_SPH(MODHD, OMNI, NENER, ENERGY, MVAR, PS, LUN)

      IMPLICIT NONE

      INCLUDE 'interface.h'

      REAL*8 ENERGY(N_Energy)
      INTEGER*4 OMNI, NENER, LUN, NB /3/, I
      INTEGER*4 HEAD(2000), HSIZE /1000/, KSKIP /0/, MSIZE, NBLOCK
      CHARACTER*80 SDUM
      CHARACTER(*) MODHD

      TYPE(META_VAR) :: MVAR

      TYPE(PSANN) :: PS

      EXTERNAL MV_STORE, PS_STORE

      COMMON /TPHEAD/ HEAD

      CALL READ_HEADER_PLUS(LUN, HEAD, HSIZE, KSKIP, NBLOCK, MV_STORE,
     &                      PS_STORE)

      MSIZE = -1
      CALL GET_META_STRING(HEAD, 'TRP_MOD', MODHD, MSIZE)
      MSIZE = 1
      CALL GET_META_INT(HEAD, 'TRP_OMN', OMNI, MSIZE)
      CALL MV_READ('ENERGY', MVAR)
      READ (MVAR%TXT(1),*) SDUM, NENER, (ENERGY(I),I=1,NENER)
      CALL MV_READ('TRP', MVAR)
      CALL PS_READ(PS, NB)

      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE READ_SPD(NENER, FLUX)

      IMPLICIT NONE

      INTEGER*4 NENER, HEAD(2000), KOLUMN(NENER), KSIZE, KSKIP /0/
      INTEGER*4 NUMBER, NDIM, INDEX
      REAL*8 FLUX(*)
      CHARACTER*80 UNITS, CTITLE, FOOTER

      COMMON /TPHEAD/ HEAD

      KSIZE = 0
      CALL SELECT_ONE(HEAD, 'Flux', UNITS, CTITLE, NDIM, INDEX, KSIZE,
     &                KOLUMN, NENER)

      NUMBER = 1
      CALL GET_VAR(HEAD, KOLUMN, KSIZE, KSKIP, NUMBER, FLUX, KSIZE,
     &             FOOTER)

      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE READ_SHIELD(LUN, ZDIST)

      IMPLICIT NONE

      INCLUDE 'interface.h'

      INTEGER*4 LUN, I, NCOL, KOLUMN(10), NUMBER /-1/, NBLOCK
      INTEGER*4 HEAD(2000), HSIZE /1000/, KSKIP /0/, KSIZE /10/, MSIZE
      PARAMETER (NCOL=5)
      REAL*8 VALUES(NCOL, 1000)
      CHARACTER*80 TVALUE(NCOL), UNITS(NCOL), CTITLE(NCOL), FOOTER

      TYPE(SHLDDIST) :: ZDIST

      EXTERNAL MV_STORE, PS_STORE

      DATA TVALUE /'LOW_LIM', 'SHI_THI', 'HI_LIM', 'SHI_PER', 'SHI_CUM'/

      CALL READ_HEADER_PLUS(LUN, HEAD, HSIZE, KSKIP, NBLOCK, MV_STORE,
     &                      PS_STORE)

      MSIZE = 1
      CALL GET_META_INT(HEAD, 'SEC_UNT', ZDIST%UNITS, MSIZE)
      CALL SELECT_VAR(HEAD, TVALUE, NCOL, UNITS, CTITLE, KOLUMN, KSIZE)
      NUMBER = N_Shield + 1
      CALL GET_VAR(HEAD, KOLUMN, KSIZE, KSKIP, NUMBER, VALUES, NCOL,
     &             FOOTER)
      IF (FOOTER(1:3) .NE. 'End')
     &  CALL ERR_WRITE(' Error reading shielding distribution file')

      IF (NUMBER .GT. N_Shield) CALL
     &  ERR_WRITE(' Too many entries in shielding distribution file')
      ZDIST%N = NUMBER
      DO I=1,ZDIST%N
        ZDIST%ZLow(I) = VALUES(1,I)
        ZDIST%Z(I) = VALUES(2,I)
        ZDIST%ZHigh(I) = VALUES(3,I)
        ZDIST%Fract(I) = VALUES(4,I)
        ZDIST%Cum(I) = VALUES(5,I)
      END DO

      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE READ_SHIELD_SSAT(LUN, ZDIST, DOSE)

C ******************************************************************************
C * MK-20090326:
C *     READ_SHIELD_SSAT() is called by niel.for, sd.for and sd2.for when 
C *     JSHLD=2 (shielding depths imported from SSAT).  The values read from
C *     spenvis_sso.txt are divided by 0.27 and marked as mm (UNITS=3).  This 
C *     behaviour is only correct when the values are provided in g/cm2.  But 
C *     SSAT could provide either g/cm2 (ok), cm or radL (?).  Units should
C *     thus be checked and an error issued if needed.
C *
C *     spenvis_sso.txt is produced by SSAT, it consists of:
C *                   - two blocks per material, 
C *                   - two extra blocks for all,
C *                   - possibly three extra blocks when SSAT in dose mode 
C *       even Blocks = Shielding distribution
C *         Metas:      SSAT, Target_Position, Shielding_Material, 
C *                     Shielding_Material_Density 
C *         Variables:  B_lo, cm, scalar, Lower edge of the bin
C *                     B_up, cm, scalar, Upper edge of the bin
C *                     B_mean, cm, scalar, Mean pos. of the recorded bin data
C *                     Value, (frac. of solid ang.), scalar, Bin value
C *                     Error, (frac. of solid ang.), scalar, Error in bin value
C *     READ_SHIELD_SSAT() reads the "Shielding distribution for all".
C *
C ******************************************************************************

      IMPLICIT NONE

      INCLUDE 'interface.h'

      INTEGER*4 LUN, DOSE, I, NCOL, KOLUMN(10), NUMBER /-1/, NBLOCK
      INTEGER*4 HEAD(2001), HSIZE /2000/, KSKIP /0/, KSIZE /10/
      PARAMETER (NCOL=5)
      REAL*8 VALUES(NCOL, 1000)
      CHARACTER*80 TVALUE(NCOL), UNITS(NCOL), CTITLE(NCOL), FOOTER

      TYPE(SHLDDIST) :: ZDIST

      EXTERNAL MV_STORE, PS_STORE

      DATA TVALUE /'B_lo', 'B_up', 'B_mean', 'Value', 'Error'/

      REWIND(LUN)
      CALL READ_HEADER_PLUS(LUN, HEAD, HSIZE, KSKIP, NBLOCK, MV_STORE,
     &                      PS_STORE)

      REWIND(LUN)

      IF (DOSE .LE. 0) THEN
        KSKIP = NBLOCK
      ELSE
        KSKIP = NBLOCK - 3
      END IF
      CALL READ_HEADER_PLUS(LUN, HEAD, HSIZE, KSKIP, NBLOCK, MV_STORE,
     &                      PS_STORE)

      CALL SELECT_VAR(HEAD, TVALUE, NCOL, UNITS, CTITLE, KOLUMN, KSIZE)

      IF ( UNITS(3) .NE. 'g/cm2' )
     &  CALL ERR_WRITE(' Shielding distribution file SHALL be in g/cm2')

      NUMBER = N_Shield + 1
      KSKIP = 0
      CALL GET_VAR(HEAD, KOLUMN, KSIZE, KSKIP, NUMBER, VALUES, NCOL,
     &             FOOTER)
      IF (FOOTER(1:3) .NE. 'End')
     &  CALL ERR_WRITE(' Error reading shielding distribution file')

      IF (NUMBER .GT. N_Shield) CALL
     &  ERR_WRITE(' Too many entries in shielding distribution file')
      ZDIST%N = NUMBER
      DO I=1,ZDIST%N
        ZDIST%ZLow(I) = VALUES(1,I) 
        ZDIST%Z(I) = VALUES(3,I) 
        ZDIST%ZHigh(I) = VALUES(2,I) 
        IF ((ZDIST%ZLow(I) .GT. 0.0D0) .AND. (ZDIST%Z(I) .LE. 0.0D0))
     &    ZDIST%Z(I) = (ZDIST%ZLow(I)+ZDIST%ZHigh(I)) / 2.0D0
        ZDIST%Fract(I) = VALUES(4,I) * 100.0D0
        IF (I .EQ. 1) THEN
          ZDIST%Cum(I) = ZDIST%Fract(I)
        ELSE
          ZDIST%Cum(I) = ZDIST%Cum(I-1) + ZDIST%Fract(I)
        END IF
      END DO

      ZDIST%UNITS = 2

      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE READ_ATH(LUN)

      IMPLICIT NONE

      INCLUDE 'interface.h'

      INTEGER*4 LUN, HEAD(2000), HSIZE /1000/, KSKIP /0/, NBLOCK

      COMMON /ATHEAD/ HEAD

      CALL READ_HEADER(LUN, HEAD, HSIZE, KSKIP, NBLOCK)

      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE READ_ATD(ATTANI, SATVEL, SATSUN)

      IMPLICIT NONE

      INTEGER*4 HEAD(2000), KOLUMN(15), KSIZE, KSKIP /0/
      INTEGER*4 NUMBER, NCOL /3/, NDIM /15/, I, J
      REAL*8 ATTANI(3,3), SATVEL(3), SATSUN(3), COORD(15)
      CHARACTER*80 TCOORD(3), UNITS(15), CTITLE(15), FOOTER

      COMMON /ATHEAD/ HEAD

      DATA TCOORD /'Attitude', 'SatVel', 'SatSun'/

      KSIZE = 100
      CALL SELECT_VAR(HEAD, TCOORD, NCOL, UNITS, CTITLE, KOLUMN, KSIZE)

      NUMBER = 1
      CALL GET_VAR(HEAD, KOLUMN, KSIZE, KSKIP, NUMBER, COORD, NDIM,
     &             FOOTER)

      DO I=1,3
        DO J=1,3
          ATTANI(J,I) = COORD(3*(I-1)+J)
        END DO
        SATVEL(I) = COORD(9+I)
        SATSUN(I) = COORD(12+I)
      END DO

      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DERIV_LAG(NP, X, Y, DY)

      IMPLICIT NONE

      INTEGER*4 NP, I
      REAL*8 X(*), Y(*), DY(*), X01, X02, X12

      DO I=1,NP
        IF (Y(I) .GE. 0.0D0) THEN
          DY(I) = 0.0D0
        ELSE
          DY(I) = Y(I)
        END IF
      END DO

      IF (NP .LT. 3) THEN
        DO I=1,NP
          DY(I) = 0.0D0
        END DO
      ELSE
        X01 = X(1) - X(2)
        X02 = X(1) - X(3)
        X12 = X(2) - X(3)
        DY(1) = -Y(1) * (X01+X02) / X01 / X02 +
     &          Y(2) * X02 / X01 / X12 -
     &          Y(3) * X01 / X02 / X12
        DY(1) = MAX(DY(1), 0.0D0)
        DO I=2,NP-1
          X01 = X(I-1) - X(I)
          X02 = X(I-1) - X(I+1)
          X12 = X(I) - X(I+1)
          DY(I) = -Y(I-1) * X12 / X01 / X02 -
     &            Y(I) * (1.0D0/X12-1.0D0/X01) +
     &            Y(I+1) * X01 / X02 / X12
          DY(I) = MAX(DY(I), 0.0D0)
        END DO
        X01 = X(NP-2) - X(NP-1)
        X02 = X(NP-2) - X(NP)
        X12 = X(NP-1) - X(NP)
        DY(NP) = Y(NP-2) * X12 / X01 / X02 -
     &           Y(NP-1) * X02 / X01 / X12 +
     &           Y(NP) * (X02+X12) / X02 / X12
        DY(NP) = MAX(DY(NP), 0.0D0)
      END IF

      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLAMDA(GCLAT, DLONG, B, FL, RINV, LAMBDA, COLAT, ELONG,
     &                  GMAGMO, MLAT, MLON)

      IMPLICIT NONE

      REAL*8 GCLAT, DLONG, B, FL, RINV, LAMBDA, COLAT, ELONG, GMAGMO
      REAL*8 DEGRAD, COSTD, SINTD, COSPD, SINPD, COST, SINT, COSP, SINP
      REAL*8 THETA, MLAT, MLON, X, Y, Z, AD, R1

      DEGRAD = 45.0D0 / ATAN(1.0D0)

      COSTD = COS(COLAT/DEGRAD)
      SINTD = SIN(COLAT/DEGRAD)
      COSPD = COS(ELONG/DEGRAD)
      SINPD = SIN(ELONG/DEGRAD)
      COST = SIN(GCLAT/DEGRAD)
      SINT = COS(GCLAT/DEGRAD)
      COSP = COS(DLONG/DEGRAD)
      SINP = SIN(DLONG/DEGRAD)

      Z = SINTD * COSPD * SINT * COSP + SINTD * SINPD * SINT * SINP +
     &    COSTD * COST
      THETA = ASIN(1.0D0) - ACOS(Z)
      MLAT = THETA * DEGRAD
      X = COSTD * COSPD * SINT * COSP + COSTD * SINPD * SINT * SINP -
     &    SINTD * COST
      Y = -SINPD * SINT * COSP + COSPD * SINT * SINP
      MLON = ATAN2(Y, X) * DEGRAD
      
c-HE
      CALL RLAMDA_BL( B, FL, GMAGMO, RINV, LAMBDA)
      IF (THETA .LT. 0.0D0 .AND. LAMBDA .GT. 0) LAMBDA = -LAMBDA

      RETURN
      END
C
C-----------------------------------------------------------------------
C
C     Calculate Radius and Lambda, using
C     method defined by Roberts, JGR vol 69, No.23, Dec 1964
C     pp 5089-5090, but using Newton Root finder to solve
C     for (COS(lambda)^2) = ps

      SUBROUTINE RLAMDA_BL( B, FL, GMAGMO, RINV, LAMBDA)
      IMPLICIT NONE
      REAL*8 B, FL, GMAGMO, RINV, LAMBDA
      REAL*8 AD, R1, DEGRAD
      DEGRAD = 45.0D0 / ATAN(1.0D0)

      IF (FL .GT. 0.0D0) THEN
        AD = (B/GMAGMO)**2 * FL**6
        R1 = 0.0D0
        RINV = 0.5D0
        DO WHILE (ABS((R1-RINV)/RINV) .GT. 0.000001D0)
          R1 = RINV
          RINV = R1 - (AD*R1**6.0D0+3.0D0*R1-4.0D0) /
     &           (6.0D0*AD*R1**5.0D0+3.0D0)
        END DO
        IF (RINV .LE. 1.0D0) THEN
          LAMBDA = ACOS(SQRT(RINV)) * DEGRAD
c-HE      IF (THETA .LT. 0.0D0) LAMBDA = -LAMBDA
        ELSE
          LAMBDA = -1.0D20
        END IF
        RINV = RINV * FL
      ELSE
        RINV = -1.0D20
        LAMBDA = -1.0D20
      END IF

      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE VRANGE(VARBLE, NVAL, VAL1, VAL2, STEP, IDIR, ILOG,
     &                  NMAX)

      IMPLICIT NONE

      REAL*8 VARBLE(*), VAL1, VAL2, STEP
      INTEGER*4 NVAL, IDIR, ILOG, NMAX, I

      IF (IDIR .LE. 0) THEN
        IF (NVAL .EQ. 1) THEN
          VAL2 = VAL1
          STEP = 0.0D0
        ELSE
          IF (ILOG .LE. 0) THEN
            STEP = (VAL2-VAL1) / (NVAL-1)
          ELSE
            STEP = LOG10(VAL2/VAL1) / (NVAL-1)
          END IF
        END IF
      ELSE
        IF (VAL1 .EQ. VAL2) THEN
          NVAL = 1
          VARBLE(1) = VAL1
        ELSE
          IF (ILOG .LE. 0) THEN
            NVAL = NINT((VAL2-VAL1)/STEP) + 1
          ELSE
            NVAL = NINT((VAL2-VAL1)/10.0D0**STEP) + 1
          END IF
        END IF
      END IF
      NVAL = MIN(NVAL, NMAX)
      DO I=1,NVAL
        IF (ILOG .LE. 0) THEN
          VARBLE(I) = VAL1 + (I-1) * STEP
        ELSE
          VARBLE(I) = VAL1 * 10.0D0**((I-1)*STEP)
        END IF
      END DO

      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE MAKEBDF(NVAL, VALUES, LENFMT, DATFMT, LUN)

      IMPLICIT NONE

      REAL*4 VALUES(*)
      INTEGER*4 NVAL, LENFMT, LUN, LENGTH, L1, L2, L3, I, J
      CHARACTER(*) DATFMT
      CHARACTER*80 FMT1, FMT2, FMT3, FMT4

      LENGTH = LENFMT + 1
      L1 = 66 / LENGTH
      L2 = NVAL / L1
      L3 = MOD(NVAL, L1)
      WRITE (FMT1,'(''(A8,I2,A1,A'',I2,'',A6)'')') LEN(DATFMT)
      WRITE (FMT2,FMT1) '(5X,''&'',', L1, '(', DATFMT, ','',''))'
      IF (L3 .GT. 0) THEN
        IF (L3 .GT. 1) THEN
          WRITE (FMT3,'(''(A8,I2,A1,A'',I2,'',A6,A'',I2,'',A5)'')')
     &          LEN(DATFMT), LEN(DATFMT)
          WRITE (FMT4,FMT3) '(5X,''&'',', L3-1, '(', DATFMT,
     &                      ','',''),', DATFMT, ',''/'')'
        ELSE
          WRITE (FMT3,'(''(A8,A'',I2,'',A5)'')') LEN(DATFMT)
          WRITE (FMT4,FMT3) '(5X,''&'',', DATFMT, ',''/'')'
        END IF
      ELSE
        WRITE (FMT3,'(''(A8,I2,A1,A'',I2,'',A6,A'',I2,'',A5)'')')
     &        LEN(DATFMT), LEN(DATFMT)
        WRITE (FMT4,FMT3) '(5X,''&'',', L1-1, '(', DATFMT,
     &                    ','',''),', DATFMT, ',''/'')'
      END IF

      IF (L3 .EQ. 0) THEN
        L2 = L2 - 1
        L3 = 1
      END IF
      DO I=1,L2
        WRITE (LUN,FMT2) (VALUES(J),J=(I-1)*L1+1,I*L1)
      END DO
      IF (L3 .GT. 0) WRITE (LUN,FMT4) (VALUES(J),J=L2*L1+1,NVAL)

      RETURN
      END
C-----------------------------------------------------------------------
C
      SUBROUTINE READ_MLET(LUN, MLET, PRODEF, TITLE)

      IMPLICIT NONE

      INCLUDE 'interface.h'

      INTEGER*4 LUN, NTRAJ, I, NCOL, KOLUMN(10), NUMBER /-1/
      INTEGER*4 HEAD(2000), HSIZE /1000/, KSKIP /0/, KSIZE /10/, MSIZE
      INTEGER*4 NBLOCK
      PARAMETER (NCOL=2)
      REAL*8 VALUES(NCOL, N_Energy+1)
      CHARACTER*80 TVALUE(NCOL), UNITS(NCOL), CTITLE(NCOL), FOOTER
      CHARACTER*80 PRODEF, TITLE
      CHARACTER*3 MODABB
      TYPE(STOPPOW) :: MLET

      TYPE(META_VAR) :: METAVAR

      EXTERNAL MV_STORE, PS_STORE

      DATA TVALUE /'Energy', 'LET'/

      CALL READ_HEADER_PLUS(LUN, HEAD, HSIZE, KSKIP, NBLOCK, MV_STORE,
     &                      PS_STORE)
c      WRITE(*,*) 'HEAD(1)=', HEAD(1), ', NBLOCK=', NBLOCK
      MSIZE = -1
      CALL GET_META_STRING(HEAD, 'MOD_ABB', MODABB, MSIZE)
c      write(*,*) 'MODABB=', MODABB
      IF (MODABB .EQ. 'MLET')
     &  CALL GET_META_STRING(HEAD, MODABB//'_ABS', TVALUE(1), MSIZE)
      CALL MV_READ('PRJ_DEF', METAVAR)
      IF ((METAVAR%N) .GT. 0) THEN
        CALL GET_META_STRING(HEAD, 'PRJ_DEF', PRODEF, MSIZE)
      ELSE
        WRITE (PRODEF,'(80(1H ))')
      END IF
      CALL GET_META_STRING(HEAD, 'PLT_HDR', TITLE, MSIZE)
      WRITE (MLET%LABEL,'(80(1H ))')
      MLET%LABEL = TITLE
      WRITE (TITLE,'(80(1H ))')
      CALL GET_META_STRING(HEAD, 'PRJ_HDR', TITLE, MSIZE)
      MSIZE = 1
      CALL SELECT_VAR(HEAD, TVALUE, NCOL, UNITS, CTITLE, KOLUMN, KSIZE)
      NUMBER = N_Energy + 1
      CALL GET_VAR(HEAD, KOLUMN, KSIZE, KSKIP, NUMBER, VALUES, NCOL,
     &             FOOTER)
c      write(*,*) 'footer=', FOOTER
      IF (FOOTER(1:3) .NE. 'End')
     &  CALL ERR_WRITE(' Error reading spectrum file')

      NUMBER = MIN(N_Energy, NUMBER)
      MLET%N = 0
      DO I=1,NUMBER
c        IF (VALUES(2,I) .GT. 0.0D0) THEN
          MLET%N = MLET%N + 1
          MLET%ENERGY(MLET%N) = VALUES(1,I)
          MLET%LET(MLET%N) = MAX(VALUES(2,I), 0.0D0)
c        END IF
      END DO
      
      RETURN
      END
C
C-----------------------------------------------------------------------
C Subroutine to read:
c - the GCR and ACR spectra from gcf.txt 
c - the solar fluence/flux spectra from sep.text/sefflare.txt
C - the trapped proton spectra from tri.text
C INPUT:
C - LUN : filename
C - SPEC : spectrumname
C - NTRAJ : segment number
C - TERM : short (=0) or long term SEU rate (>0)

      SUBROUTINE READ_SPECSPEC(LUN, SPEC2, PRODEF, TITLE, KPLANET,
     &                         NTRAJ, MISDUR, TERM)

      IMPLICIT NONE

      INCLUDE 'interface.h'
      INTEGER*4 N_VAL, HSIZE, TERM
      PARAMETER (N_VAL=N_Species*2+1, HSIZE=20000)

      INTEGER*4 LUN, NTRAJ, OMNI, I, NCOL, KOLUMN(N_VAL), NUMBER /-1/
      INTEGER*4 HEAD(HSIZE), KSKIP, KSIZE /10/, MSIZE
      INTEGER*4 NBLOCK, J, NELM, KPLANET, ELM1, GCRIEL, SEPIEL, FLAIEL
      INTEGER*4 KTRAJ
      REAL*8 MISDUR, VALUES(N_VAL, N_Energy+1)
      REAL*8 DUMMY
      REAL*4 SMISDUR, ORBSTART, ORBEND
      CHARACTER*80 TVALUE(3), UNITS(N_VAL), CTITLE(N_VAL)
      CHARACTER*80 PRODEF, TITLE, FOOTER
      CHARACTER*3 MODABB
      CHARACTER*4 PROTON 

      TYPE(SPECSPEC) :: SPEC2

      EXTERNAL MV_STORE, PS_STORE

      DATA TVALUE /'Energy', 'IFlux', 'DFlux'/
      
c
      WRITE(*,*) 'READ_SPECSPEC: lun,ntraj,term= ',LUN,NTRAJ,TERM
      KSKIP = 0
      IF ((TERM .GT. 0) .AND. (NTRAJ .GT. 1)) KSKIP = NTRAJ*TERM
   10 continue
      WRITE(*,*) 'CALL READ_HEADER_PLUS: kskip= ',KSKIP
      CALL READ_HEADER_PLUS(LUN, HEAD, HSIZE, KSKIP, NBLOCK, MV_STORE,
     &                      PS_STORE)
      KSKIP = 0
      WRITE(*,*) 'HEAD(1)=', HEAD(1), ', NBLOCK=', NBLOCK
      IF ( HEAD(1) .EQ. -1 ) 
     &     CALL ERR_WRITE(' Error reading spectrum file (EOF)')
      MSIZE = -1
      CALL GET_META_STRING(HEAD, 'MOD_ABB', MODABB, MSIZE)
      WRITE(*,*) 'CALLED GET_META_STRING: modabb= ',MODABB
c
      IF (MODABB .NE. 'TRP') THEN
        PROTON = 'proton'
        MSIZE = -N_Species
        CALL GET_META_STRING(HEAD, 'SPECIES', SPEC2%SPECIES, MSIZE)
        WRITE(*,*) 'CALLED GET_META_STRING(SPECIES) MSIZE =', MSIZE
        IF ( MSIZE .EQ. 0 ) THEN
          WRITE(*,*) 'Error caught: jump to next block'
          CALL GET_VAR(HEAD, 1, 1, -1, 0, DUMMY, 1, FOOTER)
          GOTO 10
        ELSE IF ( SPEC2%SPECIES(1) .eq. PROTON  ) THEN
          WRITE(*,*) 'Proton: jump to next block'
          CALL GET_VAR(HEAD, 1, 1, -1, 0, DUMMY, 1, FOOTER)
          GOTO 10
        END IF
        WRITE(*,*) 'SPEC2%SPECIES(1)= ', SPEC2%SPECIES(1)
        NELM = abs( MSIZE)
      ELSE
        NELM=1	  
      END IF
      
      SPEC2%NSPEC = NELM
      MSIZE = -1
      CALL GET_META_STRING(HEAD, 'PRJ_DEF', PRODEF, MSIZE)
      IF ( MSIZE .EQ. 0 ) WRITE (PRODEF,'(80(1H ))')
      MSIZE = -1
      CALL GET_META_STRING(HEAD, 'PLT_HDR', TITLE, MSIZE)
      IF ( MSIZE .EQ. 0 ) WRITE (TITLE,'(80(1H ))')
      SPEC2%LABEL = TITLE
      WRITE (TITLE,'(80(1H ))')
      CALL GET_META_STRING(HEAD, 'PRJ_HDR', TITLE, MSIZE)
   
      KPLANET = 0
      MSIZE = 1
      CALL GET_META_INT(HEAD, 'MIS_PLA', KPLANET, MSIZE)
      KTRAJ = 0
      MSIZE = 1
      CALL GET_META_INT(HEAD, 'MIS_NTR', KTRAJ, MSIZE)
      WRITE(*,*) 'CALLED GET_META_INT: ktraj,msize= ',KTRAJ, MSIZE

      OMNI = -1
      MSIZE = -1
      CALL GET_META_STRING(HEAD, 'MOD_ABB', MODABB, MSIZE)
      IF ( MSIZE .NE. 0 ) THEN
          MSIZE = 1
          CALL GET_META_INT(HEAD, MODABB//'_OMN', OMNI, MSIZE)
      ELSE
           WRITE(*,*) 'Error caught: MOD_ABB not defined'
      END IF
      WRITE(*,*) 'MODABB_OMN: omni= ',OMNI
      SPEC2%OMNI = OMNI
      SPEC2%TYPE = -1

c      SEGDUR = 0.0D0
c      MSIZE = 1
c      WRITE(*,*) 'GET_META_STRING(ORB_MJD) MSIZE =', MSIZE
c      CALL GET_META(HEAD, 'ORB_MJD', ORBSTART, MSIZE)
c      WRITE(*,*) 'MSIZE =', MSIZE
c      IF (MSIZE .EQ. 1) THEN
c          WRITE(*,*) 'GET_META_STRING(ORB_TSE) MSIZE =', MSIZE
c          CALL GET_META(HEAD, 'ORB_TSE', ORBEND, MSIZE)
c          WRITE(*,*) 'MSIZE =', MSIZE
c          IF (MSIZE .EQ. 1) SEGDUR = ORBEND - ORBSTART
c      ELSE
c        WRITE(*,*) 'Error caught: ORB_MJD not defined'
c      END IF
      
      SMISDUR = 0.0
      MSIZE = 1
      CALL GET_META(HEAD, 'MIS_DUR', SMISDUR, MSIZE)
      MISDUR = SMISDUR
      NCOL = 3
      KSIZE = N_VAL
c      write(*,*) 'NELM=', NELM, ', KSIZE =', KSIZE, ', NCOL =', NCOL
      WRITE(*,*) 'CALL SELECT_VAR: ncol,ksize=', ncol,ksize
      CALL SELECT_VAR(HEAD, TVALUE, NCOL, UNITS, CTITLE, KOLUMN, KSIZE)
      WRITE(*,*) 'CALLED SELECT_VAR: ncol,ksize=', ncol,ksize
c      write(*,*) 'NELM=', NELM, ', KSIZE =', KSIZE, ', NCOL =', NCOL
      MSIZE = 1
      IF ( KSIZE .NE. SPEC2%NSPEC*2+1 ) 
     &     CALL ERR_WRITE(' Error reading spectrum file (NSPEC)')
      IF (MODABB .EQ. 'GCR') THEN
        CALL GET_META_INT(HEAD, 'GCR_IEL', GCRIEL, MSIZE)
        ELM1 = GCRIEL
      ELSEIF (MODABB .EQ. 'FLA') THEN
        CALL GET_META_INT(HEAD, 'FLA_IEL', FLAIEL, MSIZE)
        ELM1 = FLAIEL
      ELSEIF (MODABB .EQ. 'SEP') THEN
        CALL GET_META_INT(HEAD, 'SEP_IEL', SEPIEL, MSIZE)
        ELM1 = SEPIEL
      ELSE
        ELM1 = 1
      END IF
	  
      NUMBER = N_Energy + 1
      NCOL = N_VAL
      KSKIP = 0
      WRITE(*,*) 'CALL GET_VAR: number,ncol= ',NUMBER,NCOL
      CALL GET_VAR(HEAD, KOLUMN, KSIZE, KSKIP, NUMBER, VALUES, NCOL,
     &             FOOTER)
      IF (FOOTER(1:3) .NE. 'End')
     &     CALL ERR_WRITE(' Error reading spectrum file (NENERG)')
     
      NUMBER = MIN(N_Energy, NUMBER)
      SPEC2%N = NUMBER
      DO I=1,NUMBER
          SPEC2%ENERGY(I) = VALUES(1,I)
          DO J=1,NELM
            SPEC2%IFLUX(I,ELM1+J-1) = MAX(VALUES(J+1,I), 0.0D0)
            SPEC2%DFLUX(I,ELM1+J-1) = MAX(VALUES(J+1+NELM,I), 0.0D0)
          END DO
      END DO
        
      RETURN
      END
C
C-----------------------------------------------------------------------
C Subroutine to read the view factor matrix from vfo.txt

      SUBROUTINE READ_VIEWF(LUN, VWF, PRODEF, TITLE)

      IMPLICIT NONE

      INCLUDE 'interface.h'
      INTEGER*4 N_VAL, HSIZE, NCOL
      PARAMETER (N_VAL=N_Surfaces+2, HSIZE=20000)
      INTEGER*4 LUN, NTRAJ, OMNI, I, KOLUMN(N_VAL), NUMBER /-1/
      INTEGER*4 HEAD(HSIZE), KSKIP/0/, KSIZE /10/, MSIZE
      INTEGER*4 NBLOCK, J, NRSRF
      REAL*8 VALUES(N_VAL, N_VAL)
      CHARACTER*80 TVALUE(3), UNITS(N_VAL), CTITLE(N_VAL)
      CHARACTER*80 PRODEF, TITLE, FOOTER

      TYPE(VWFACTOR) :: VWF

      EXTERNAL MV_STORE, PS_STORE

      DATA TVALUE /'Area', 'Vfactor', '1-sum'/
      
      CALL READ_HEADER_PLUS(LUN, HEAD, HSIZE, KSKIP, NBLOCK, MV_STORE,
     &                      PS_STORE)
     
c      IF ( HEAD(1) .EQ. -1 ) 
c     &     CALL ERR_WRITE(' Error reading viewfactor file (EOF)')
      MSIZE = -N_Surfaces
      CALL GET_META_STRING(HEAD, 'VOL_NAM', VWF%VOL_NAM, MSIZE)
      MSIZE = -N_Surfaces
      CALL GET_META_STRING(HEAD, 'SRF_NAM', VWF%SRF_NAM, MSIZE)
      MSIZE = -1
      CALL GET_META_STRING(HEAD, 'PRJ_DEF', PRODEF, MSIZE)
      IF ( MSIZE .EQ. 0 ) WRITE (PRODEF,'(80(1H ))')
      MSIZE = -1
      CALL GET_META_STRING(HEAD, 'PLT_HDR', TITLE, MSIZE)
      IF ( MSIZE .EQ. 0 ) WRITE (TITLE,'(80(1H ))')
      VWF%LABEL = TITLE
      WRITE (TITLE,'(80(1H ))')
      CALL GET_META_STRING(HEAD, 'PRJ_HDR', TITLE, MSIZE)
      MSIZE = 1
      CALL GET_META_INT(HEAD, 'VWF_GRP', VWF%NGRP, MSIZE)
      MSIZE = 1
      CALL GET_META_INT(HEAD, 'VWF_FCS', NRSRF, MSIZE)

      NCOL = 3
      KSIZE = N_VAL
      CALL SELECT_VAR(HEAD, TVALUE, NCOL, UNITS, CTITLE, KOLUMN, KSIZE)
      MSIZE = 1
      NUMBER = N_Surfaces +2 
      NCOL = N_VAL
      CALL GET_VAR(HEAD, KOLUMN, KSIZE, KSKIP, NUMBER, VALUES, NCOL,
     &             FOOTER)
      IF (FOOTER(1:3) .NE. 'End')
     &     CALL ERR_WRITE(' Error reading viewfactor file (EOF)')
     
      NUMBER = MIN(N_Surfaces, NRSRF)
      VWF%NSRF = NUMBER
      DO I=1,VWF%NSRF
          VWF%AREA(I) = VALUES(1,I)
          DO J=1,VWF%NSRF
            VWF%MATRIX(I,J) = MAX(VALUES(J+1,I), 0.0D0)
          END DO
          VWF%LOSSC(I) = VALUES(VWF%NSRF+2,I)
          VWF%LOSSR(I) = VALUES(I+1,VWF%NSRF+1)
      END DO
        
      RETURN
      END
C
C-----------------------------------------------------------------------

