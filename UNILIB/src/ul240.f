# 1 "ul240.f"
      SUBROUTINE UL240
     :          (mlab, ifail)
C
C!    Evaluate Hilton's function (er31)
C
C     REF.: TREND TN1 eqs. 2.30 and 2.31
C           HILTON,H.H.,JGR,76,6952,1971
C
      INCLUDE 'structure.h'
cDEC$ IF DEFINED (_x86_)
cDEC$ ATTRIBUTES DLLEXPORT :: UL240
cDEC$ ENDIF
C
C     INTERFACE
C
        TYPE(ZLBL) :: mlab
        INTEGER*4     ifail
C     
      COMMON /UC140/  mint, mext, msun
C      
        TYPE(zimf) :: mint
        TYPE(zsun) :: msun
        TYPE(zemf) :: mext          
C                                     
      COMMON /UC160/  pi, deg, re, gmagmo, eclipt, geoid, uma
C
        REAL*8        pi, deg, re, gmagmo, eclipt, geoid(3), uma(30)
C                                              
      COMMON /UC190/  prop, stepx, stpmin, umsq, upsq, uk2, uk3, epskm,
     :                epsrel, stplst, xclat, kmflg, nxstp       
C
        REAL*8        prop, stepx, stpmin, umsq, upsq, uk2, uk3, epskm,
     :                epsrel, stplst, xclat
        INTEGER*4     kmflg, nxstp       
C                
C     VARIABLES
C
        REAL*8        gm, eq230, g, tier, g3
        REAL*8        fim, fbm
C
C     CODE
C
      if( .not. mlab%linv )then
        ifail       = -24002
        mlab%llmi   = .false.
        mlab%lkauf  = .false.
        return
      endif
      if( .not. mlab%lbmp )then
        ifail       = -24003
        mlab%llmi   = .false.
        mlab%lkauf  = .false.
        return
      endif
C
      fim           = mlab%finv
      fbm           = mlab%fbmp
      tier          = 1.0d0 / 3.0d0
C
      if( mod(kmflg,10) .eq. 0 )then
        gm          = gmagmo
        ifail       = 0
      else
        gm          = mint%gmmo
        ifail       = 1
      endif
C
      g             = ( abs(fim) / re )**3 * fbm / gm
      g3            = g**tier
      eq230         = 1.0D0 + 1.350474D0 * g3 + 0.465380D0 * 
     :                g3**2 + 0.047546D0 * g
      mlab%flmi     = ( eq230 * gm / fbm )**tier
      mlab%llmi     = .true.
C
      mlab%fkauf    = abs(fim) * sqrt(fbm)
      mlab%lkauf    = .true.
C
      if( fim .lt. 0.0d0 ) 
     :        ifail = -24001 
C
C
      END
C----------------------------------------------------------------------
      SUBROUTINE UL242
     :                (mlab, ifail)
C
C     Inverse Hilton's function
C
      INCLUDE 'structure.h'
cDEC$ IF DEFINED (_x86_)
cDEC$ ATTRIBUTES DLLEXPORT :: UL242
cDEC$ ENDIF
C
C     INTERFACE
C
        TYPE(zlbl) ::  mlab
        INTEGER*4      ifail
C
      COMMON /UC140/  mint, mext, msun
C
        TYPE(zimf) :: mint
        TYPE(zsun) :: msun
        TYPE(zemf) :: mext          
C         
      COMMON /UC160/  pi, deg, re, gmagmo, eclipt, geoid, uma
C
        REAL*8        pi, deg, re
        REAL*8        gmagmo
        REAL*8        eclipt, geoid(3), uma(30)
C         
      COMMON /UC190/  prop, stepx, stpmin, umsq, upsq, uk2, uk3, 
     :                 epskm, epsrel, stplst, xclat, kmflg, nxstp       
C
        REAL*8        prop, stepx, stpmin
        REAL*8        umsq, upsq, uk2, uk3
        REAL*8        epskm, epsrel, stplst, xclat
        INTEGER*4     kmflg, nxstp       
C             
C     VARIABLES
C 
        REAL*8        gm, eq230, fbm, xz, xz3, g3, fim
        REAL*8        ahilt(4), coef(6), xq, sxq, xu
        logical       ldebug
C
        save ahilt,coef
        data ahilt/ 1.0D0, 1.350474d0, 0.465380d0, 0.0475546d0/
        data coef / 6*0.0d0/
C
C     CODE
C
      mlab%linv  = .FALSE.
      mlab%lkauf = .FALSE.                       
      ldebug     = ifail.lt.-900
c
      if (mod(kmflg,10) .eq. 0) then
        gm       = gmagmo
        ifail    = 0
      else
        gm       = mint%gmmo
        ifail    = 1
      endif                    
c
      if (mlab%flmi .le. 0. .or. .not. mlab%llmi ) then
        ifail    = -24202
        return
      endif
      if (mlab%fbmp .le. 0. .or. .not. mlab%lbmp ) then
        if (mlab%lalp0) then
          mlab%fbmp = gm/(sin(mlab%falp0*deg)**2*mlab%flmi**3)
          mlab%lbmp = .true.
          if (ldebug) print*,'Unilib-UL242: L, B, B0 =',mlab%fbmp,
     :                       gm/mlab%flmi**3
        else
          ifail    = -24201
          return
        endif
      endif  
c
      if( coef(1) .eq. 0.0D0 )then
        coef(1)  = ahilt(2)/ahilt(4)
        coef(2)  = ahilt(3)/ahilt(4)
        coef(3)  = 2*coef(2)**3-9*coef(1)*coef(2)
        coef(5)  = coef(2)**2-3*coef(1)
        coef(4)  = 4*coef(5)**3
        coef(6)  = 1.0d0/3.0d0
      endif
c
      if (ldebug) print*,'Unilib-UL242: L, B, B0 =',mlab%flmi,
     :                       mlab%fbmp,gm/mlab%flmi**3
c
      fbm        = mlab%fbmp
      eq230      = mlab%flmi**3 * fbm / gm
      xq         = (ahilt(1)-eq230)/ahilt(4)
      xu         = coef(3)+27*xq
      sxq        = xu**2-coef(4)
      if( sxq .le. 0 )then
         ifail   = -24203
         return
      endif
      xz3        = 0.5d0*(sqrt(sxq)-xu)
      xz         = xz3**coef(6)
      g3         = (xz-coef(2)+coef(5)/xz)/3
      fim        = g3*(gm/fbm)**coef(6) * re
      if( abs(fim) .lt. 1.e-7 ) fim=0
      if( fim .lt. 0 )then
c       B is less than B0
        ifail    = -24204
        if (ldebug) print*,'Unilib-UL242: I, g3, xz =',fim,
     :                       g3, xz
        return
      endif
      mlab%linv  = .TRUE.
      mlab%finv  = fim
      mlab%lkauf = .TRUE.
      mlab%fkauf = fim*SQRT(fbm)
C
      END         
C----------------------------------------------------------------------
      SUBROUTINE UL245
     :            ( mlin, mlab, ifail )
C
C     Evaluate the equatorial pitch angle
C
      INCLUDE 'structure.h'
C
C     INTERFACE
C
        TYPE(zfln) :: mlin
        TYPE(zlbl) :: mlab
        INTEGER*4     ifail
C
      COMMON /UC140/  mint, mext, msun
C      
        TYPE(zimf) :: mint
        TYPE(zsun) :: msun
        TYPE(zemf) :: mext          
C                                     
      COMMON /UC160/  pi, deg, re, gmagmo, eclipt, geoid, uma
C
        REAL*8        pi, deg, re
        REAL*8        gmagmo
        REAL*8        eclipt, geoid(3), uma(30) 
C                                                   
      COMMON /UC190/  prop, stepx, stpmin, umsq, upsq, uk2, uk3, 
     :                 epskm, epsrel, stplst, xclat, kmflg, nxstp       
C
        REAL*8        prop, stepx, stpmin
        REAL*8        umsq, upsq, uk2, uk3
        REAL*8        epskm, epsrel, stplst, xclat
        INTEGER*4     kmflg, nxstp 
C        
      COMMON /UC192/  xrmin, xbmin, xtmin, xbmax, epslon, epsfl,
     :                fvet, pvet, epsomeg, dltlat
C
        REAL*8        xrmin, xbmin, xtmin, xbmax, epslon, epsfl
        REAL*8        fvet, pvet, epsomeg, dltlat                   
C
C     VARIABLES
C
        REAL*8        gm, b0, sn2
C
C     CODE
C
      mlab%lalp0 = .FALSE.
      if (.not. mlab%lbmp) then
        ifail    = -24501
        return
      endif
C
      if ( mod(kmflg,10) .eq. 0 )then
        gm = gmagmo
      else
        gm = mint%gmmo
      endif
C
      if (kmflg .lt. 10) then 
        if (.not. mlab%llmi) then
          ifail  = -24502
          return
        endif
        b0       = gm/mlab%flmi**3
      else
        b0       = mlin%equat%b%dnrm
        if (b0 .lt. xbmin .or. b0 .gt. xbmin) then
          ifail  = -24503
          return
        endif
      endif
C
      sn2        = b0/mlab%fbmp
      if (sn2.gt.1) then
        if (sn2-1.0D0.gt.1.0D-9 ) ifail = -24504
        sn2      = 1.0d0
      endif
      mlab%lalp0 = .TRUE.
      mlab%falp0 = ASIN (SQRT(sn2))/deg
C
C      
      END
C----------------------------------------------------------------------