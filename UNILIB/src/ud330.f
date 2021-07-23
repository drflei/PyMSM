# 1 "ud330.f"
      subroutine ud330(phi,star,ifail)
c
c     evaluate the third invariant
c
      INCLUDE 'structure.h'
cDEC$ IFDEFINED (_x86_)
cDEC$ ATTRIBUTES DLLEXPORT :: UD330
cDEC$ ENDIF
C      
C     INTERFACE
C
        REAL*8        phi,star
        INTEGER*4     ifail
C
      COMMON /UC110/  mlbl, k1st, klmp
C
        INTEGER*4     k1st, klmp        
        TYPE(zlbl) :: mlbl
C  
      COMMON /UC140/ mint, mext, msun
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
     :                epskm, epsrel, stplst, xclat, kmflg, nxstp       
C
        REAL*8        prop, stepx, stpmin
        REAL*8        umsq, upsq, uk2, uk3
        REAL*8        epskm, epsrel, stplst, xclat
        INTEGER*4     kmflg, nxstp
C
C     VARIABLES
C
        REAL*8        colat,gm
c
C     CODE
C
        mlbl%lphi=.false.
        colat=0.0d0
        call ud331(colat,phi,ifail)
        if(ifail.lt.0)return
c
        if(mod(kmflg,10).eq.0)then
          gm=gmagmo
        else
          gm=mint%gmmo
        endif
        star=2*pi*gm*re**2/abs(phi)
c
        mlbl%lphi=.true.
        mlbl%fphi=phi
c
C
      END
C------------------------------------------------------------------
      subroutine ud331(colat,phi,ifail)
c
c     evaluate the magnetic flux through a spherical cap
c
      INCLUDE 'structure.h'
C
C     INTERFACE
C      
        REAL*8        colat,phi
        INTEGER*4     ifail
C        
      COMMON /UC110/  mlbl, k1st, klmp
C
        INTEGER*4     k1st, klmp        
        TYPE(zlbl) :: mlbl
C        
      COMMON /UC120/  nbrfl, kurfl, mfl
C
        INTEGER*4     nbrfl, kurfl
        TYPE(zfln) :: mfl(nx120)
C         
C     VARIABLES     
C
        REAL*8        phi0
        INTEGER*4     nc,knext
C
C     CODE
C
      if( k1st.le.0)then
        ifail=-33101
        return
      endif
      kurfl=k1st
c
      phi=0.0d0
      nc=0
      if( colat.eq.0.0d0 )then
  100   continue
         knext=mfl(kurfl)%keast
         call ud332(colat,mfl(kurfl)%footpn,
     :                    mfl(knext)%footpn,phi0,ifail)
         if(ifail.lt.0)return
         nc=nc+ifail
         phi=phi+phi0
         kurfl=knext
        if(kurfl.ne.k1st)goto 100
      elseif(colat.eq.180.0d0 )then
  200   continue
         knext=mfl(kurfl)%kwest
         call ud332(colat,mfl(kurfl)%footps,
     :                    mfl(knext)%footps,phi0,ifail)
         if(ifail.lt.0)return
         nc=nc+ifail
         phi=phi+phi0
         kurfl=knext
        if(kurfl.ne.k1st)goto 200
      else
        ifail=-33102
        return
      endif
      ifail=nc
C
C      
      END
C----------------------------------------------------------------
      subroutine ud332(colat,mpos1,mpos2,phi,ifail)
C
C     evaluate the magnetic flux through a spherical pie
C
      INCLUDE 'structure.h'
C
C     INTERFACE
C
        REAL*8        colat,phi
        TYPE(zgeo) :: mpos1, mpos2
        INTEGER*4     ifail
C        
      COMMON /UC160/  pi, deg, re, gmagmo, eclipt, geoid, uma
C
        REAL*8        pi, deg, re
        REAL*8        gmagmo
        REAL*8        eclipt, geoid(3), uma(30)
C
      COMMON /UC190/  prop, stepx, stpmin, umsq, upsq, uk2, uk3, 
     :                epskm, epsrel, stplst, xclat, kmflg, nxstp       
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
        REAL*8        cl, cls1, cls2, dscl, dlg, omeg,
     :                dslg, elg, els1, els2, spec, elg1, elg2
        INTEGER*4     nlat, nlon, nc
        TYPE(zgeo) ::  mgeo
        TYPE(zvec) ::  mb
C
C     CODE
C   
      elg1 = mpos1%elong
      elg2 = mpos2%elong
      dlg= elg2-elg1
      if( dlg .ge. 180.0d0 )then
        elg1=elg1+360.0d0
        dlg=dlg-360.0d0
      elseif( dlg .le. -180.0d0 )then
        elg2=elg2+360.0d0
        dlg=dlg+360.0d0
      endif
C
      if( colat .eq. 0.0d0 )then
        if( mpos1%colat .lt. mpos2%colat )then
           cls1=mpos1%colat
           cls2=mpos2%colat
           els1=elg1
           els2=elg2
        else
           cls1=mpos2%colat
           cls2=mpos1%colat
           els1=elg2
           els2=elg1
        endif
        nlat=cls1/dltlat+1
        dscl=cls1/nlat
      elseif( colat .eq. 180.0d0 )then
        if( mpos1%colat .lt. mpos2%colat )then
           cls1=mpos2%colat
           cls2=mpos1%colat
           els1=elg2
           els2=elg1
        else
           cls1=mpos1%colat
           cls2=mpos2%colat
           els1=elg1
           els2=elg2
        endif
        nlat=1-(cls1-180.0d0)/dltlat
        dscl=(cls1-180.0d0)/nlat
      else
        ifail=-33201
        return
      endif
C
c      print*,' UD332: ',mpos1%radius-re
      mgeo%radius=(mpos1%radius+mpos2%radius)*0.5d0
      if( abs(mpos1%radius-mgeo%radius).gt.2.0D0*epskm )then
        ifail=-33202
        print*,' UD332: FPs =',mpos1%radius-re,mpos2%radius-re
        print*,' UD332: del =',epskm
        return
      endif
      nc=0
c
      phi=0.0d0
      cl=dscl*0.5d0+colat
      do while ((cl-colat)*(cl-cls1) .lt. 0)
c
        mgeo%colat=cl
        omeg=( sin(cl*deg)*dscl*dlg*deg**2 )
        nlon=abs(omeg)/epsomeg+1
        dslg=dlg/nlon
        elg=elg1+0.5d0*dslg       
C
        omeg=omeg/nlon
        do while ((elg-elg1)*(elg-elg2).lt.0)
c
          mgeo%elong=elg
c-test-          write(4,*)mgeo%radius,mgeo%colat,mgeo%elong
          call um530(mgeo,mb,ifail)
          if(ifail .lt. 0 .and. 
     :       (ifail.ne.-53003 .and.ifail.ne.-53102) )return
          phi=phi+omeg*mb%rho
          nc=nc+1
c
          elg=elg+dslg
        enddo
c
        cl=cl+dscl
      enddo
c
      dscl=cls2-cls1
      nlat=abs(dscl/dltlat)+1
      dscl=dscl/nlat
      spec=(els2-els1)/nlat
      cl=cls1+dscl*0.5d0
      els1=els1+0.5d0*spec
      do while ((cl-cls1)*(cl-cls2) .lt. 0)
        mgeo%colat=cl
        dlg=els2-els1
        omeg=( sin(cl*deg)*dscl*dlg*deg**2 )*0.5d0
        nlon=abs(omeg)/epsomeg+1
        dslg=dlg/nlon
        elg=els1+0.5d0*dslg       
        omeg=omeg/nlon
        do while ((elg-els1)*(elg-els2).lt.0)
          mgeo%elong=elg
c-test-          write(4,*)mgeo%radius,mgeo%colat,mgeo%elong
          call um530(mgeo,mb,ifail)
          if(ifail .lt. 0 .and. 
     :        (ifail.ne.-53003.and.ifail.ne.-53102) )return
          phi=phi+omeg*mb%rho
          nc=nc+1
          elg=elg+dslg
        enddo       
        cl=cl+dscl
        els1=els1+spec
      enddo

c
      ifail=nc
      phi=phi*mgeo%radius**2
C
C      
      END
C-------------------------------------------------------------------------