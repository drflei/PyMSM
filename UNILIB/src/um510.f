# 1 "um510.f"
      SUBROUTINE UM510
     :          (kint, year, lbint, kunit, ifail)
C
C!    Select a geomagnetic field model
C_    er25
C
C     REFERENCES
C     Subroutine readc, bint.for of blxtra
C
      INCLUDE 'structure.h'
cDEC$ IF DEFINED (_x86_)
cDEC$ ATTRIBUTES DLLEXPORT :: UM510
cDEC$ ENDIF
C
      EXTERNAL  UM511, UM512, UM513
C
C     INTERFACE               
C
        REAL*8        year
        INTEGER*4     kint, ifail, kunit
        CHARACTER*(*) lbint
C
      COMMON /UC140/  mint, mext, msun
C      
        TYPE(zimf) :: mint
        TYPE(zsun) :: msun
        TYPE(zemf) :: mext          
C
      COMMON /UC160 
     :               /pi, deg, re, gmagmo, eclipt, geoid, uma
C
        REAL*8        pi, deg, re, gmagmo, eclipt, geoid(3), uma(30)
C     
C     VARIABLES
C
        CHARACTER*32  lbl
        REAL*8        tm, gfc(2,2)
C
C     CODE
C
      ifail = 0
C
      mint%saarot = 0.0d0
      tm    = year
      IF     ( kint .EQ. 0 ) THEN
        CALL UM513 (tm, lbl)
      ELSEIF ( kint .EQ. -1) THEN
        CALL UM511 (tm, lbl)
        mint%saarot = 0.3d0 * (year-tm)
        lbl(22:32) = '+ SAA rot. '
      ELSEIF ( kint .EQ. 1 ) THEN
        CALL UM511 (tm, lbl)
      ELSEIF ( kint .EQ. -2) THEN
        tm = 1970.0d0
        CALL UM512 (tm, lbl)
        mint%saarot = 0.3d0 * (year-tm)
        lbl(22:32) = '+ SAA rot. '
      ELSEIF ( kint .EQ. 2 ) THEN
        CALL UM512 (tm, lbl)
      ELSEIF ( kint .EQ. 3 ) THEN
        CALL UM513 (tm, lbl)
        lbl         = 'Dipolar magnetic field          '
        mint%norder = 2         
      ELSEIF ( kint .EQ. 4 ) THEN
c
c 2015-07-01: disabled by Yves Christophe, according to suggestion of 
c   Daniel Heynderickx
c
        ifail       = -51001
        mint%norder =  0
        lbl         = '** DEPRECATED: IGRF/DGRF w Kluge'
        lbint       = lbl( 1:MIN( LEN(lbint),LEN(lbl) ) )
        return                                 
c
c        CALL UM513 (tm, lbl)
c        CALL UM517 (mint%norder,mint%coef)
c
      ELSE
        ifail       = -51001
        mint%norder =  0
        lbl         = 'ERROR***************************'
        lbint       = lbl( 1:MIN( LEN(lbint),LEN(lbl) ) )
        return                                 
      ENDIF
c
      mint%label  = lbl(1:32)
      mint%kinner = kint
      lbint       = lbl( 1:MIN( LEN(lbint),LEN(lbl) ) )
      gfc(1,2)    = mint%coef(1,2)
      gfc(2,1)    = mint%coef(2,1)
      gfc(2,2)    = mint%coef(2,2)
      mint%gmmo   = 1.0D-5 * SQRT( gfc(2,1)**2 + gfc(1,2)**2 + 
     :                            gfc(2,2)**2 )
      mint%epoch  = tm
      if( gfc(2,2).eq.0 .or. gfc(2,1).eq.0 )then
        mint%elong  = 0
        mint%colat  = 0
        ifail=-51002
        if( kunit.le.0 )return
        write(kunit,1200,err=2000)
        write(kunit,1140,err=2000) kint, lbl, mint%tzero, mint%norder, 
     :                    tm, 
     :                    mint%colat, mint%elong, mint%gmmo, 
     :                    mint%saarot
        return
      endif
      mint%elong  = ATAN( gfc(1,2)/gfc(2,2) ) / deg - mint%saarot
      mint%colat  = ATAN( SQRT(gfc(1,2)**2+gfc(2,2)**2)/gfc(2,1) )
     :                 / deg
c
      IF( kunit.gt.0 )THEN
        write(kunit,1003,err=2000)
        write(kunit,1140,err=2000) kint, lbl, mint%tzero, mint%norder,
     :                    tm, 
     :                    mint%colat, mint%elong, mint%gmmo, 
     :                    mint%saarot
      ENDIF
      ifail       = mint%norder
c
      return
 2000 ifail       = -51003
C
C
 1003 format(/' --- Geomagnetic field model ---')
 1140 format(/6x,'Model (',i2,'): ',a32,3x,'Epoch',f6.0,/29x,
     : 'Order + 1 =',i5,/21x,'Calculation epoch =',f7.1,5x,' year',/9x,
     :           'Colatitude of the dipole pole =',f8.2,4x,' deg',/10x,
     :            'Longitude of the dipole pole =',f8.2,4x,' deg',/19x,
     :                         'Earth dipole moment =',f12.6,' G/Re^3',
     :        /10x,'Correction for the SAA drift =',f8.2,4x,' deg')
 1200 format(/' --- Geomagnetic field uncorrectly loaded ---')
      END
C----------------------------------------------------------------------
C
      SUBROUTINE UM511
     :           (tm, lbl)
C
C!    Jensen & Cain 1960 model coefficients
C
      INCLUDE 'structure.h'
C
C     INTERFACE               
C
        REAL*8        tm
        CHARACTER*32  lbl
C
      COMMON /UC140/  mint, mext, msun
C      
        TYPE(zimf) :: mint
        TYPE(zsun) :: msun
        TYPE(zemf) :: mext          
C
C     VARIABLES
C
        REAL*8        gjc(27,2)
        INTEGER*4     i, j, k, n
C
      DATA ((gjc(i,j),j=1,2),i=1,27) 
     : / 30411.2d0,    0.0d0,  2147.4d0, -5798.9d0,  2403.5d0,    0.0d0,
     :   -5125.3d0, 3312.4d0, -1338.1d0,  -157.9d0, -3151.8d0,    0.0d0,
     :    6213.0d0, 1487.0d0, -2489.8d0,  -407.5d0,  -649.6d0,   21.0d0,
     :   -4179.4d0,    0.0d0, -4529.8d0, -1182.5d0, -2179.5d0, 1000.6d0,
     :     700.8d0,   43.0d0,  -204.4d0,   138.5d0,  1625.6d0,    0.0d0,
     :   -3440.7d0,  -79.6d0, -1944.7d0,  -200.0d0,   -60.8d0,  459.7d0,
     :     277.5d0,  242.1d0,    69.7d0,  -121.8d0, -1952.3d0,    0.0d0,
     :    -485.3d0, -575.8d0,   321.2d0,  -873.5d0,  2141.3d0, -340.6d0,
     :     105.1d0,  -11.8d0,    22.7d0,  -111.6d0,   111.5d0,  -32.5d0/
C
C     CODE
C
      mint%tzero     = 1960.0D0
      lbl            = 'Jensen & Cain 1960              '
      tm             = mint%tzero
      mint%norder    = 7
      mint%coef(1,1) = 0.0D0
      n = 1
      DO k=2,mint%norder
        mint%coef(k,1)     = gjc(n,1)
        DO i=2,k
          mint%coef(k,i)   = gjc(n+i-1,1)
          mint%coef(i-1,k) = gjc(n+i-1,2)
        END DO
        n = n+k
      END DO
C
C
      END
C----------------------------------------------------------------------
C
      SUBROUTINE UM512
     :           (tm, lbl)
C
C!    GSFC 12/66 model coefficients
C
      INCLUDE 'structure.h'
C
      EXTERNAL UM515
C
C     INTERFACE  
C             
        REAL*8        tm
        CHARACTER*32  lbl
C
      COMMON /UC140/  mint, mext, msun
C      
        TYPE(zimf) :: mint
        TYPE(zsun) :: msun
        TYPE(zemf) :: mext          
C
C     VARIABLES
C
        REAL*8        ggsfc(65,2), ggsfc1(65,2), ggsfc2(65,2)
        REAL*8        gg(nx140,nx140), ggt(nx140,nx140), 
     :                ggtt(nx140,nx140), t
        INTEGER*4     i, j, k, n
C
      DATA ((ggsfc(i,j),j=1,2),i=1,65)
     : / -30401.2d0,0.0d0,-2163.8d0,5778.2d0,-1540.1d0,0.0d0,
     : 2997.9d0,-1932.0d0,1590.3d0,202.9d0,1307.1d0,0.0d0,-1988.9d0,
     : -425.4d0,1276.8d0,227.8d0,881.2d0,-133.8d0,949.3d0,0.0d0,
     : 803.5d0,160.3d0,502.9d0,-274.3d0,-397.7d0,2.3d0,266.5d0,
     : -246.6d0,-233.5d0,0.0d0,355.7d0,5.1d0,228.4d0,117.8d0,-28.8d0,
     : -114.8d0,-157.9d0,-108.9d0,-62.2d0,82.4d0,49.2d0,0.0d0,57.5d0,
     : -12.1d0,-0.8d0,104.4d0,-238.3d0,56.6d0,-1.5d0,-23.4d0,-2.0d0,
     : -14.8d0,-108.9d0,-13.3d0,72.2d0,0.0d0,-53.7d0,-53.7d0,7.9d0,
     : -27.4d0,15.6d0,-8.1d0,-24.3d0,7.0d0,-3.6d0,24.3d0,15.5d0,
     : -22.5d0,3.6d0,-21.4d0,8.5d0,0.0d0,6.5d0,5.4d0,-9.3d0,-11.7d0,
     : -9.6d0,4.2d0,-6.1d0,-15.3d0,5.5d0,4.6d0,-8.1d0,21.9d0,13.0d0,
     : -0.7d0,7.4d0,-17.1d0,10.4d0,0.0d0,5.8d0,-22.4d0,7.5d0,13.8d0,
     : -15.1d0,6.3d0,12.1d0,-3.0d0,4.7d0,-1.9d0,0.2d0,9.0d0,1.6d0,
     : 11.5d0,0.9d0,0.1d0,0.2d0,-1.5d0,-2.9d0,0.0d0,-0.9d0,-0.1d0,
     : -2.2d0,4.5d0,0.8d0,-1.0d0,-2.8d0,2.6d0,6.4d0,-4.4d0,4.7d0,
     : -1.3d0,-0.2d0,-3.6d0,1.8d0,4.0d0,2.0d0,1.0d0,1.1d0,-2.0d0 /     
      DATA ((ggsfc1(i,j),j=1,2),i=1,65)
     : / 14.03d0,0.00d0,8.76d0,-3.71d0,-23.29d0,0.00d0,-0.09d0,
     : -14.31d0,-4.56d0,-16.62d0,-0.93d0,0.00d0,-10.62d0,5.20d0,2.31d0,
     : 2.53d0,-5.89d0,-6.98d0,1.45d0,0.00d0,0.90d0,-2.19d0,-1.75d0,
     : -0.14d0,0.66d0,1.88d0,-3.01d0,-6.52d0,1.61d0,0.00d0,0.60d0,
     : 2.24d0,3.34d0,1.59d0,-0.04d0,-2.61d0,-0.60d0,0.50d0,1.76d0,
     : -0.12d0,-0.42d0,0.00d0,0.82d0,0.05d0,0.82d0,0.09d0,2.35d0,
     : 2.55d0,0.83d0,-1.19d0,0.01d0,0.33d0,0.23d0,0.84d0,-0.57d0,
     : 0.00d0,-0.34d0,-0.96d0,-1.44d0,0.01d0,-0.90d0,0.43d0,0.03d0,
     : 0.75d0,-0.60d0,-0.33d0,-0.17d0,0.49d0,-0.64d0,0.90d0,0.35d0,
     : 0.00d0,0.50d0,-0.50d0,1.70d0,-0.21d0,-0.11d0,0.03d0,0.34d0,
     : -0.79d0,-0.07d0,0.05d0,0.43d0,0.10d0,-0.15d0,-0.36d0,-0.42d0,
     : -0.43d0,-0.10d0,0.00d0,-0.13d0,0.66d0,-1.20d0,0.54d0,0.08d0,
     : 0.03d0,-0.08d0,0.35d0,-0.39d0,-0.03d0,-0.36d0,-0.01d0,0.47d0,
     : 0.45d0,0.37d0,-0.05d0,-0.46d0,0.75d0,-0.01d0,0.00d0,-0.13d0,
     : -0.61d0,0.88d0,-0.64d0,-0.18d0,0.02d0,0.17d0,0.05d0,-0.02d0,
     : -0.63d0,0.05d0,-0.07d0,0.17d0,0.07d0,0.16d0,-0.03d0,0.31d0,
     : -0.02d0,-0.23d0,-0.45d0 /
      DATA ((ggsfc2(i,j),j=1,2),i=1,65)
     : / -0.062d0,0.000d0,0.114d0,-0.043d0,-0.154d0,0.000d0,-0.018d0,
     : 0.054d0,-0.253d0,-0.016d0,-0.123d0,0.000d0,-0.027d0,0.095d0,
     : 0.028d0,-0.007d0,-0.183d0,0.079d0,0.001d0,0.000d0,-0.044d0,
     : 0.004d0,0.017d0,0.056d0,0.007d0,-0.035d0,-0.097d0,-0.047d0,
     : 0.045d0,0.000d0,0.001d0,-0.046d0,0.075d0,0.007d0,0.008d0,
     : -0.007d0,0.015d0,0.001d0,0.056d0,-0.024d0,-0.006d0,0.000d0,
     : 0.015d0,0.020d0,0.010d0,-0.011d0,0.050d0,0.015d0,-0.011d0,
     : -0.029d0,0.026d0,0.029d0,0.023d0,-0.010d0,-0.014d0,0.000d0,
     : -0.006d0,-0.014d0,-0.034d0,0.016d0,-0.004d0,0.014d0,-0.006d0,
     : 0.005d0,-0.027d0,-0.008d0,-0.001d0,0.016d0,-0.004d0,0.011d0,
     : 0.006d0,0.000d0,0.008d0,-0.015d0,0.039d0,-0.012d0,-0.008d0,
     : 0.005d0,0.015d0,-0.011d0,-0.002d0,0.000d0,0.005d0,-0.003d0,
     : -0.008d0,-0.009d0,-0.007d0,-0.003d0,-0.005d0,0.000d0,-0.001d0,
     : 0.022d0,-0.027d0,0.007d0,0.005d0,-0.002d0,-0.007d0,0.009d0,
     : -0.006d0,0.006d0,-0.009d0,-0.001d0,0.006d0,0.009d0,0.005d0,
     : -0.004d0,-0.009d0,0.019d0,-0.003d0,0.000d0,-0.003d0,-0.012d0,
     : 0.020d0,-0.014d0,-0.008d0,0.001d0,0.007d0,0.001d0,0.001d0,
     : -0.011d0,0.001d0,-0.001d0,0.001d0,0.001d0,0.005d0,-0.001d0,
     : 0.004d0,0.001d0,-0.002d0,-0.006d0 /
C
C     CODE
C
      mint%tzero      = 1960.0D0
      lbl             = 'GSFC 12/66 120 Term             '
C
      mint%norder     = 11
      gg(1,1)         = 0.0D0
      ggt(1,1)        = 0.0D0
      ggtt(1,1)       = 0.0D0        
      n = 1
      DO k=2,mint%norder
        gg(k,1)       = ggsfc(n,1)
        ggt(k,1)      = ggsfc1(n,1)
        ggtt(k,1)     = ggsfc2(n,1)
        DO j=2,k
          gg(k,j)     = ggsfc(n+j-1,1)
          gg(j-1,k)   = ggsfc(n+j-1,2)
          ggt(k,j)    = ggsfc1(n+j-1,1)
          ggt(j-1,k)  = ggsfc1(n+j-1,2)
          ggtt(k,j)   = ggsfc2(n+j-1,1)
          ggtt(j-1,k) = ggsfc2(n+j-1,2)
        END DO
        n = n + k
      END DO
C
      t = tm - mint%tzero
      CALL UM515(t, mint%norder, gg, ggt, ggtt, mint%coef)
C
C
      END
C----------------------------------------------------------------------
C
      SUBROUTINE UM513
     :           (tm, lbl)
C
C!    DGRF/IGRF model coefficients
C
      INCLUDE 'structure.h'
C
C-v112      EXTERNAL dgrfbd
C      
C     INTERFACE
C                          
        REAL*8         tm
        CHARACTER*32   lbl
C
      COMMON /UC140/   mint, mext, msun
        TYPE(zimf) ::  mint
        TYPE(zsun) :: msun
        TYPE(zemf) :: mext          
C
C-v112      COMMON /DGRF
C-v112     :                /tzero, label, ggbd, ggtbd, nmax
        INTEGER*4      nx
        PARAMETER (nx=24)          
        REAL*8         tzero(nx)
        REAL*8         ggbd(nx,nx140,nx140), ggtbd(nx,nx140,nx140)
        INTEGER*4      nmax(nx)
        CHARACTER*32   label(nx)
C
C     VARIABLES
C
        INTEGER*4      i,j,k
        REAL*8         gg(nx140,nx140), ggt(nx140,nx140)
        REAL*8         ggtt(nx140,nx140), t
C
C-V112-
C
      INCLUDE 'um513.h'
C
C     CODE
C
      if( tzero(1).ne.1900.0d0 )then
        mint%tzero     = 0.0d0
        mint%norder    = 0
        lbl            = '!!!  IGRF/DGRF data missing  !!!'
        mint%coef(1,2) = 0.0d0
        mint%coef(2,2) = 0.0d0
        mint%coef(2,1) = 0.0d0
        return
      endif
C
      DO i=1,nx
        IF ( tzero(i) .GT. tm ) GO TO 100
      END DO
      i = nx + 1
  100 i = i - 1
      mint%tzero    = tzero(i)
      lbl           = label(i)
      mint%norder   = MIN( nmax(i) + 1 , nx140 )
      DO k=1,mint%norder
        DO j=1,mint%norder
          gg(k,j)   = ggbd(i,k,j)
          ggt(k,j)  = ggtbd(i,k,j)
          ggtt(k,j) = 0.0D0
        END DO
      END DO
      t = tm - mint%tzero
      CALL UM515(t, mint%norder, gg, ggt, ggtt, mint%coef)
C
C
      END
C----------------------------------------------------------------------
      SUBROUTINE UM515
     :           (t, nmax, gg, ggt, ggtt, g)
C
C!    Schmidt normalisation
C
      INCLUDE 'structure.h'
C
C     INTERFACE
C
        INTEGER*4 nmax          
        REAL*8    t, gg(nx140,nx140), ggt(nx140,nx140)
        REAL*8    ggtt(nx140,nx140), g(nx140,nx140)
C
C
C     VARIABLES
C                   
        INTEGER*4 n, m, jj
        REAL*8    shmit(nx140,nx140)
        SAVE      shmit
C
      INTEGER*4 nxsqr
      PARAMETER (nxsqr=nx140*nx140)
      real*8     dum
      DATA        shmit
     :            / nxsqr*0.0d0 /
C
C     CODE
C
C     INITIALIZATION
      IF ( shmit(1,1) .NE. -1.0D0 ) THEN
        shmit(1,1)       = -1.0D0
        DO n=2,nx140
          shmit(n,1)     = (2*n-3) * shmit(n-1,1) / (n-1)
          jj = 2
          DO m=2,n
            dum          = (n-m+1)*jj
            shmit(n,m)   = shmit(n,m-1)*SQRT(Dum/(n+m-2))
            shmit(m-1,n) = shmit(n,m)
            jj           = 1
          END DO
        END DO
      END IF
C
C
      DO n=1,nmax
        DO m=1,nmax
          g(n,m) = (gg(n,m)+t*ggt(n,m)+t*t*ggtt(n,m)) * shmit(n,m)
        END DO
      END DO
C
C
      END
C----------------------------------------------------------------------
      SUBROUTINE UM517
     :           (nmax, gcoef)
C
C!    Transform from "Schmidt" to "Kluge" normalisation
C
C see ESOC Internal Note no 61, G. Kluge, 1970
C
      INCLUDE 'structure.h'
C
C     INTERFACE
C
        INTEGER*4 nmax
        REAL*8    gcoef(nx140,nx140)
C
C     VARIABLES
C                   
        INTEGER*4 kl, km
        REAL*8    fctl, pwr, x
C
C
      fctl              = 1.0d0
      gcoef( 1, 1)      = - gcoef( 1, 1)
C
      do kl = 2, nmax
        fctl            = fctl * 0.5d0 * (kl-1)
        gcoef(kl, 1)    = -fctl * gcoef(kl, 1)
        pwr             = fctl * sqrt(0.5d0)
        do km = 2, kl
          x             = kl+km-2
          pwr           = pwr * sqrt(x/(kl-km+1))
          gcoef(kl,km)  = - pwr * gcoef(kl,km)
          gcoef(km-1,kl)= gcoef(kl,km)
        enddo
      enddo
C 
      END          
C----------------------------------------------------------------------
