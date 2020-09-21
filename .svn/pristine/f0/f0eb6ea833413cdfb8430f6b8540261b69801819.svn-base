      INTERFACE
        SUBROUTINE ERR_WRITE(ERRSTR, FCLOSE)
          CHARACTER(*) ERRSTR
          LOGICAL, OPTIONAL :: FCLOSE
        END SUBROUTINE
        SUBROUTINE FILE_OPEN(LUN, TYPE, EXT, FORM, RECL, STATUS, ACCESS,
     &                       IOSTAT)
          INTEGER*4 LUN
          INTEGER*4, OPTIONAL :: RECL, IOSTAT
          CHARACTER(*) TYPE, EXT
          CHARACTER(*), OPTIONAL :: FORM, STATUS, ACCESS
        END SUBROUTINE
      END INTERFACE

      TYPE PSANN
        SEQUENCE
        REAL*8 POS(100,3)
        INTEGER*4 N, IDUM
        CHARACTER*80 TXT(100)
      END TYPE PSANN

      TYPE META_VAR
        SEQUENCE
        INTEGER*4 N
        CHARACTER*500 TXT(100)
      END TYPE META_VAR

      INTEGER*4 N_Energy
      PARAMETER (N_Energy=300)
      TYPE SPECTRUM
        SEQUENCE
        REAL*8 Energy(N_Energy), IFlux(N_Energy), DFlux(N_Energy)
        INTEGER*4 N, OMNI, TYPE, IDUM
        CHARACTER*80 LABEL
      END TYPE SPECTRUM
      
      INTEGER*4 N_Species
      PARAMETER (N_Species=92)
      TYPE SPECSPEC
        SEQUENCE
        INTEGER*4 N, OMNI, TYPE, NSPEC
        REAL*8 Energy(N_Energy), IFlux(N_Energy,N_Species) 
        REAL*8 DFlux(N_Energy,N_Species)
        CHARACTER*80 LABEL
        CHARACTER*4 SPECIES(N_Species)
      END TYPE SPECSPEC

      TYPE STOPPOW
        SEQUENCE
        REAL*8 Energy(N_Energy), LET(N_Energy)
        INTEGER*4 N, TYPE
        CHARACTER*80 LABEL
      END TYPE STOPPOW

      INTEGER*4 N_Shield
      PARAMETER (N_Shield=200)
      TYPE SHLDDIST
        SEQUENCE
        REAL*8 ZLow(N_Shield), ZHigh(N_Shield), Z(N_Shield)
        REAL*8 Fract(N_Shield), Cum(N_Shield)
        INTEGER*4 N, UNITS
      END TYPE SHLDDIST

      INTEGER*4 NTRAJM
      PARAMETER (NTRAJM=100)

      TYPE G4_GPS
        SEQUENCE
        INTEGER*4 N
        CHARACTER*80 MACRO(1000)
      END TYPE G4_GPS
      
      INTEGER*4 N_Surfaces
      PARAMETER (N_Surfaces=60)
      TYPE VWFACTOR
        SEQUENCE
        INTEGER*4 NSRF, NGRP
        REAL*8 AREA(N_Surfaces), MATRIX(N_Surfaces,N_Surfaces)
        REAL*8 LOSSC(N_Surfaces), LOSSR(N_Surfaces)       
        CHARACTER*40 SRF_NAM(N_Surfaces), LABEL(N_Surfaces) 
        CHARACTER*40 VOL_NAM(N_Surfaces)
      END TYPE VWFACTOR

      CHARACTER*80 PROJECT, NAMPRO
      CHARACTER*20 VERSTR
      CHARACTER*11 SDATE
      CHARACTER*8 STIME

      COMMON /INFO/ PROJECT, NAMPRO, VERSTR, SDATE, STIME
