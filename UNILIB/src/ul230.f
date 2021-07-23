# 1 "ul230.f"
      SUBROUTINE UL230
     :          (finv, ifail)
C
C!    Evaluate the integral invariant coordinate I
C_    er09
C
      INCLUDE 'structure.h'
cDEC$ IF DEFINED (_x86_)
cDEC$ ATTRIBUTES DLLEXPORT :: UL230
cDEC$ ENDIF
C
C     INTERFACE
C
        REAL*8        finv
        INTEGER*4     ifail
C
      COMMON /UC170/ nsg, kgp, mlab, mlin, mele
C
        INTEGER*4     nsg, kgp
        TYPE(zlbl) :: mlab
        TYPE(zfln) :: mlin
        TYPE(zseg) :: mele(nx170)
C
      COMMON /UC190
     :               /prop, stepx, stpmin, umsq, upsq, uk2, uk3,
     :                epskm, epsrel, stplst, xclat, kmflg, nxstp       
C
        REAL*8        prop, stepx, stpmin, umsq, upsq, uk2, uk3, 
     :                epskm, epsrel, stplst, xclat
        INTEGER*4     kmflg, nxstp              
C                
C     VARIABLES
C
        REAL*8        fbm, tmp, ds, r2, r3, r4, res, tier
        INTEGER*4     j, jmpn, jmps
C
C     CODE
C
      mlab%linv         = .FALSE.
      mlab%lkauf        = .FALSE.
      jmpn              = mlin%ind%jmirpn
      jmps              = mlin%ind%jmirps
C      
      ifail             = 0
      if( jmpn .gt. jmps )then
        ifail           = -23001
        return
      endif
      fbm               = 0.5d0 * ( mele(jmpn)%beg%b%dnrm + 
     :                             mele(jmps)%beg%b%dnrm )
      if( ABS( mele(jmpn)%beg%b%dnrm - mele(jmps)%beg%b%dnrm )
     :                                 .gt. 5.d0 * fbm * epsrel )then
c        print*,mele(jmpn)%beg%b%dnrm,mele(jmps)%beg%b%dnrm,
c     :   ABS( mele(jmpn)%beg%b%dnrm - mele(jmps)%beg%b%dnrm )/fbm
        ifail           = -23002
        return
      endif
C
      if (.not. mlab%lbmp) then
        mlab%lbmp       = .TRUE.
        mlab%fbmp       = fbm
      else if( ABS(mlab%fbmp - fbm) .gt. 
     :         5.d0 *fbm*epsrel )then
        ifail           = -23004
        return
      endif
      fbm               = mlab%fbmp    
      mele(jmpn)%csalp  = 0.0d0
      res               = 0.0d0
      tier              = 1.0d0 / 3.0d0
C
      do j = jmpn, jmps-2
        if(   mele(j+1)%beg%b%dnrm  .gt. fbm )then
          ifail         = -23003
c          print*,'ul230',jmpn,jmps,j,mele(jmps)%beg%coord%radius
c          print*, mele(j+1)%beg%b%dnrm
          return
        endif
        if( mele(j)%rkstp(2) .lt. fbm ) then
          r2            = SQRT( 1.0d0 - mele(j)%rkstp(2)/fbm )
        else
          r2            = 0.0d0
          ifail         = ifail + 1
        endif
        if( mele(j)%rkstp(3) .lt. fbm ) then
          r3            = SQRT( 1.0d0 - mele(j)%rkstp(3)/fbm )
        else
          r3            = 0.0d0
          ifail         = ifail + 1
        endif
        if( mele(j)%rkstp(1) .lt. fbm ) then
          r4            = SQRT( 1.0d0 - mele(j)%rkstp(1)/fbm )
        else
          r4            = 0.0d0
          ifail = ifail + 1
        endif
        tmp             = ( mele(j)%csalp + r4 ) * 0.5d0
     :                           + umsq * r2 + upsq * r3
        ds              = mele(j+1)%arcl - mele(j)%arcl
        res             = res + tmp * ds * tier
        mele(j+1)%csalp = SQRT( 1.0d0 - mele(j+1)%beg%b%dnrm / fbm )
      enddo
C
      j                 = jmps - 1
      if( j .ge. jmpn )then
        if( mele(j)%rkstp(2) .lt. fbm ) then
          r2            = SQRT( 1.0d0 - mele(j)%rkstp(2)/fbm )
        else
          r2            = 0.0d0
          ifail         = ifail + 1
        endif
        if( mele(j)%rkstp(3) .lt. fbm ) then
          r3            = SQRT( 1.0d0 - mele(j)%rkstp(3)/fbm )
        else
          r3            = 0.0d0
          ifail         = ifail + 1
        endif
        if( mele(j)%rkstp(1) .lt. fbm ) then
          r4            = SQRT( 1.0d0 - mele(j)%rkstp(1)/fbm )
        else
          r4            = 0.0d0
          ifail = ifail + 1
        endif
        tmp             = (mele(j)%csalp + r4) * 0.5d0 
     :                    + umsq * r2 + upsq * r3
        ds              = mele(j+1)%arcl - mele(j)%arcl
        res             = res + tmp * ds * tier
      endif
C
      mele(jmps)%csalp  = 0.0d0
      if( mele(jmps)%arcl .gt. mele(jmpn)%arcl )then
        finv            = res
      else
        finv            = -res
      endif
C
      mlab%linv         = .TRUE.
      mlab%finv         = finv
      mlab%lkauf        = .TRUE.
      mlab%fkauf        = finv*SQRT(fbm)
C
C
      END
C----------------------------------------------------------------------