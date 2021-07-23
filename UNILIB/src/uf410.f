# 1 "uf410.f"
      SUBROUTINE UF410
     :          (flm, fbm, glon, guess, mpos, falt, ifail)
C
C!    Search the geographic position of a magnetic field line segment
C_    er05
C
      INCLUDE 'structure.h'
cDEC$ IF DEFINED (_x86_)
cDEC$ ATTRIBUTES DLLEXPORT :: UF410
cDEC$ ENDIF
C
C
      EXTERNAL UF411, UF420
C      
C     INTERFACE
C
        REAL*8        flm, fbm, glon, guess, falt
        INTEGER*4     ifail
        TYPE(zgeo) :: mpos    
C
      COMMON /UC160
     :               /pi, deg, re, gmagmo, eclipt, geoid, uma
        REAL*8        pi, deg, re, gmagmo, eclipt, geoid(3), uma(30)
C           
      COMMON /UC190
     :               /prop, stepx, stpmin, umsq, upsq, uk2, uk3,
     :                epskm, epsrel, stplst, xclat, kmflg, nxstp
C
        REAL*8        prop, stepx, stpmin, umsq, upsq, uk2, uk3,
     :                epskm, epsrel, stplst, xclat
        INTEGER*4     kmflg, nxstp                                
C              
      COMMON /UC192
     :               /xrmin, xbmin, xtmin, xbmax, epslon, epsfl,
     :                fvet, pvet, epsomeg, dltlat
        REAL*8        xrmin, xbmin, xtmin, xbmax, epslon, epsfl, 
     :                fvet, pvet, epsomeg, dltlat
C                                                   
C     VARIABLES
C
        TYPE(zlbl) :: ml(5), mdummy
        TYPE(zpnt) :: mpeq(5)
        REAL*8        rgmin, rgmax, fld, flu, colat, fid, fiu, fim 
        INTEGER*4     k, k1, k2, k3, kcnt
c
c
        LOGICAL       debug(6)
        INTEGER*4     kdebug
        real*8        faltmx
C
C     CODE
C
      if( ifail.lt.0 .and. ifail .gt.-32)then
        kdebug=ifail
        ifail=-ifail
        do k=1,5
          debug(k)= mod(ifail,2).eq.1
          ifail=ifail/2
        enddo
      else
        debug(1)=.false.
        debug(2)=.false.
        debug(3)=.false.
        debug(4)=.false.
        debug(5)=.false.
      endif
      debug(6)=debug(1).or.debug(2)
      if(debug(6)) write(*,1100) 
      if(debug(6)) write(*,1120)flm, fbm, glon, guess, falt
C
      if( flm .le. 0.0d0 .or. fbm .gt. xbmax .or. fbm .lt. xbmin )then
        ifail  =-41005
        return
      endif
C
C
      fld      = flm * (1.0d0 - epsfl)
      flu      = flm * (1.0d0 + epsfl)
C
      mdummy%lbmp = .true.
      mdummy%llmi = .true.
      mdummy%fbmp = fbm
      mdummy%flmi = flu
      call ul242(mdummy,ifail)
      if(debug(5))write(*,1110)'UL242',ifail
      if( ifail.lt.0 )return
      fiu=mdummy%finv
      mdummy%flmi = fld
      call ul242(mdummy,ifail)
      if(debug(5))write(*,1110)'UL242',ifail
      if( ifail.eq.-24204 )then
        fid=0.0D0
      elseif( ifail.lt.0 ) then
        return
      else
        fid=mdummy%finv
      endif
      fim=0.5d0*(fiu+fid)
      if(debug(6))write(*,1130) fld, flu, fid, fiu
C
      kcnt     = 0
      ifail    = 0
      if( guess .le. xrmin ) guess = ABS( flm * 0.8d0 ) * re
      rgmin    = MIN( 2.0d0 * xrmin, re )
      rgmax    = 5.0d0 * flm * re
C
C The point at the equator, along the magnetic field line
C must have its longitude equal to glon and its radius
C between rgmin and rgmax
C
C In the first version of this subroutine, the strategy was based
C on the value of L, now it is based on the value of I. The algorithm
C is the same
C
C The subroutine uf411 search a magnetic field line such as
C the point with the lowest magnetic field intensity along the
C line has a longitude equal to glon and a magn. field intensity
C less than fbm
C
C
      k1       = 1
      k2       = 2
      k3       = 3
C
C  a/ First attempt
C     -------------
C
      ifail=kdebug
      CALL UF411( glon, guess, rgmin, rgmax, -1.0d0, fbm,
     :            mpeq(k1), ml(k1), ifail )
      if(debug(5))write(*,1110)'UF411',ifail
      if( ifail.lt.0 ) return
      kcnt      = kcnt + ifail
      colat     = mpeq(k1)%coord%colat
C
      call ul240(ml(k1),ifail)
      if(debug(5))write(*,1110)'UL240',ifail
      if( ifail.lt.0 ) return    
      if(debug(1))write(*,1010)' First',ml(k1)%finv,
     :               mpeq(k1)%coord%radius,
     :               mpeq(k1)%coord%colat,mpeq(k1)%coord%elong,kcnt     
C            
C  b/ Second attempt
C     --------------
C
      if( ml(k1)%finv .gt. fiu )then
C      
C  b.1/ Lm1 > Lm --> Alt2 < Alt1 & Lm2 < Lm
C
        rgmax   = guess
        guess   = MAX( rgmin, guess * 0.9d0 * fld/ml(k1)%flmi )
  210   continue
        ifail=kdebug
        CALL UF411( glon, guess, rgmin, rgmax, colat, fbm,
     :              mpeq(k2), ml(k2), ifail)
        if(debug(5))write(*,1110)'UF411',ifail
        if( ifail.lt.0 ) return
        kcnt    = kcnt + ifail
C
       if(debug(1))write(*,1010)' Lower',ml(k2)%finv,
     :               mpeq(k2)%coord%radius,
     :               mpeq(k2)%coord%colat,mpeq(k2)%coord%elong,kcnt     
C                 
        if( ml(k2)%finv .gt. fiu )then
          rgmax = guess
          guess = ( MAX( rgmin, guess * 0.8d0) + guess ) * 0.5d0
          if( ABS( guess - mpeq(k2)%coord%radius )
     :                               .lt. epskm )then
            ifail = -41001
            return
          endif
          goto 210
        else if( ml(k2)%finv .ge. fid )then
C     ... second attempt succeeded...
C
          mpos    = mpeq(k2)%coord
          ifail   = kcnt
          goto 500       
        endif
        k       = k1
        k1      = k2
        k2      = k          
      else if( ml(k1)%finv .lt. fid )then
C      
C  b.2/ Lm1 < Lm --> Alt2 > Alt1 & Lm2 > Lm
C
        rgmin   = guess
        guess   = MIN( rgmax,
     :                 guess * 1.1d0 * flu/ml(k1)%flmi )
  220   continue
        ifail=kdebug
        CALL UF411( glon, guess, rgmin, rgmax, colat, fbm,
     :              mpeq(k2), ml(k2), ifail)
        if(debug(5))write(*,1110)'UF411',ifail
        if( ifail.lt.0 ) return
        kcnt    = kcnt + ifail
C
       if(debug(1))write(*,1010)'Higher',ml(k2)%finv,
     :               mpeq(k2)%coord%radius,
     :               mpeq(k2)%coord%colat,mpeq(k2)%coord%elong,kcnt    
C                 
        if( ml(k2)%finv .lt. fid )then
          rgmin = guess
          guess = ( MIN( rgmax, guess * 1.2d0) + guess ) * 0.5d0
          if( ABS( guess - mpeq(k2)%coord%radius ) 
     :                              .lt. epskm )then
            ifail = -41002
            return
          endif
          goto 220
        else if( ml(k2)%finv .le. fiu )then
C     ... second attempt succeeded...
C
          mpos    = mpeq(k2)%coord
          ifail   = kcnt
          goto 500       
        endif
      else
C      
C  b.3/ First attempt succeeded...
C
        mpos    = mpeq(k1)%coord
        ifail   = kcnt
        goto 500  
      endif  
C              
C  c/ Iterative search
C     ----------------
C
C  ===>  Lm1 < Lm < Lm2
C
  300 continue
      rgmin     = mpeq(k1)%coord%radius
      rgmax     = mpeq(k2)%coord%radius
      colat     = ( mpeq(k1)%coord%colat 
     :            + mpeq(k2)%coord%colat ) * 0.5d0
      guess     = rgmin + ( rgmax - rgmin ) * 
     :        ( fim-ml(k1)%finv ) / ( ml(k2)%finv-ml(k1)%finv )
      ifail=kdebug
      CALL UF411( glon, guess, rgmin, rgmax, colat, fbm,
     :            mpeq(k3), ml(k3), ifail)
      if(debug(5))write(*,1110)'UF411',ifail
      if( ifail.lt.0 ) return
      kcnt      = kcnt + ifail
C
      if(debug(1))write(*,1010)' Iter.',ml(k3)%finv,
     :               mpeq(k3)%coord%radius,
     :               mpeq(k3)%coord%colat,mpeq(k3)%coord%elong,kcnt
C                 
      if( kcnt .gt. 100 )then
        kcnt    = -41003
      elseif( ml(k3)%finv .gt. fiu )then
C      
C  c.1/ too high
C
        k       = k3
        k3      = k2
        k2      = k
        goto 300
      elseif( ml(k3)%finv .lt. fid )then
C      
C  c.2/ too low
C
        k       = k3
        k3      = k1
        k1      = k
        goto 300
      endif
C
      mpos      = mpeq(k3)%coord
      ifail     = kcnt
      if( ifail .lt. 0 )return
C
 500  continue
c
      if( falt .ge. -500.0d0 )then
        CALL UF420(mpos, fbm, 1, falt, k)
        if( k .gt. nx170 ) k=-42699
        if(debug(5))write(*,1110)'UF420',k
        if( k .lt. 0 ) then
          if(debug(1))write(*,1030)k
          faltmx = -999.0d0
          CALL UF420(mpos, fbm, 1, faltmx, ifail)
          if( ifail .gt. nx170 ) ifail=-42699
          if( ifail .ge. 0 )then
            ifail = -41004
          else
c-29jan04-            print*,'UF410: UF420-> (UF428) ifail=',k
c-29jan04-            print*,'       -> UF420[-999] -> ifail=',ifail
          endif
          return
        endif
      endif    
C
      call uf415(k)
        if(debug(5))write(*,1110)'UF515',k
        if( k .lt. 0 ) ifail=k
      if(debug(1))write(*,1020)ifail
      if( ifail.lt.0 ) return    
C
 1000 format('  Input: L=',f6.3,' B=',f6.4,' Long=',f5.1,
     :       ' guess=',f7.1,' foot alt=',f6.1,
     :      /7x,': Limits on L=(',f6.3,',',f6.3,') I=(',f7.0,','
     :       f7.0,')')
 1010 format(7x,': ',a6,' attempt, I=',f7.0,' Eq=(',f7.1,';',f5.1,
     :       ';',f5.1,') cnt=',i2)
 1020 format(7x,': The last attempt was successfull, ifail =',i6)
 1030 format(7x,': Problem with the altitude boundary, ifail =',i6)
 1100 format(/' UNILIB: Subroutine UF410 switched to debug mode')
 1110 format(7x,': call to ',a5,' returned with IFAIL =',i7)
 1120 format(7x,': L=',f6.3,' B=',f6.4,' Long=',f5.1,
     :       ' guess=',f7.1,' foot alt=',f6.1)
 1130 format(7x,': Limits on L=(',f6.3,',',f6.3,') I=(',f7.0,','
     :       f7.0,')')
C
      END
C----------------------------------------------------------------------
      SUBROUTINE UF411
     :           (glon, radius, rmin, rmax, 
     :            gsclat, fbm, mpeq, ml, ifail)
C
C!    Search a local magnetic equator and compute the L value
C!    associate with Bm
C
      INCLUDE 'structure.h'
C
      EXTERNAL  UF420, UM530, UL230, UL240
C
C     INTERFACE
C                  
        TYPE(zpnt) :: mpeq
        TYPE(zlbl) :: ml
        REAL*8        fbm, glon, radius, rmin, rmax, gsclat   
        INTEGER*4     ifail          
C                                                             
      COMMON /UC140/  mint, mext, msun
C       
        TYPE(zimf) :: mint
        TYPE(zemf) :: mext
        TYPE(zsun) :: msun
C             
      COMMON /UC160
     :               /pi, deg, re, gmagmo, eclipt, geoid, uma
        REAL*8        pi, deg, re, gmagmo, eclipt, geoid(3), uma(30)
C        
      COMMON /UC170/  nsg, kgp, mlab, mlin, mele
C
        INTEGER*4     nsg, kgp
        TYPE(zlbl) :: mlab
        TYPE(zfln) :: mlin
        TYPE(zseg) :: mele(nx170)
C                              
      COMMON /UC190
     :               /prop, stepx, stpmin, umsq, upsq, uk2, uk3,
     :                epskm, epsrel, stplst, xclat, kmflg, nxstp       
        REAL*8        prop, stepx, stpmin, umsq, upsq, uk2, uk3,
     :                epskm, epsrel, stplst, xclat
        INTEGER*4     kmflg, nxstp                        
C              
      COMMON /UC192
     :               /xrmin, xbmin, xtmin, xbmax, epslon, epsfl, 
     :                fvet, pvet, epsomeg, dltlat
        REAL*8        xrmin, xbmin, xtmin, xbmax, epslon, epsfl,
     :                fvet, pvet, epsomeg, dltlat   
C                                               
C     VARIABLES
C
        REAL*8        xcl, ycl, finv, tlon, cmin, cmax
        INTEGER*4     kcount, ifl, kc1, kc2, kc3, kc4
        TYPE(zgeo) :: mpos
        TYPE(zvec) :: mb
C        
C * 
C * glon    in     = longitude of the local magnetic equator
C * radius  in/out = guessed altitude of the l.m.e.
C * rmin    in/out = lowest altitude
C * rmax    in/out = greast altitude
C * gsclat  in     = guessed co-latitude
C * fbm     in     = value of Bm
C * 
C
c
c
        LOGICAL       debug(6)
        INTEGER*4     kdebug,k
C
C     CODE
C
      if( ifail.lt.0 .and. ifail .gt.-32)then
        kdebug=ifail
        ifail=-ifail
        do k=1,5
          debug(k)= mod(ifail,2).eq.1
          ifail=ifail/2
        enddo
      else
        debug(1)=.false.
        debug(2)=.false.
        debug(3)=.false.
        debug(4)=.false.
        debug(5)=.false.
      endif
      debug(6)=debug(2)              
      if( debug(2) )write(*,1000)glon,radius,rmin,rmax,gsclat,fbm
C
C
      if( rmax-rmin .lt. epskm*0.1d0 )then
        ifail=-41107
        return
      endif
C
      cmin         = rmin
      cmax         = rmax
      kcount       = 0
      tlon         = glon
      kc1=0
      kc2=0
      kc3=0
      kc4=0
C
      if( gsclat .lt. 0 )then
        if( mint%colat .gt. 2.0d0 )then 
          xcl        = 1.0D0 / TAN( mint%colat * deg )
          ycl        = -COS(deg * ( mint%elong - glon ) )
          mpeq%coord%colat = ATAN2( xcl, ycl ) / deg
        else
          mpeq%coord%colat = 90.0d0
        endif
      else
        mpeq%coord%colat = gsclat
      endif
C
  100 continue
C  
      kcount       = kcount + 1
      if( kcount.gt.25 )then
C computation aborted after 10 iterations (case 1234)
        kcount=max(kc1,kc2,kc3,kc4)
        ifail      = -41114
        if( kc1.eq.kcount) ifail=-41111
        if( kc2.eq.kcount) ifail=-41112
        if( kc3.eq.kcount) ifail=-41113
        if(debug(4))write(*,1010)kc1,kc2,kc3,kc4
c        print*,' Remark: k1234=',kc1,kc2,kc3,kc4
        return
      endif
C      
      mpos%radius  = radius
      mpos%colat   = mpeq%coord%colat
      mpos%elong   = tlon
C-test-
      call um530( mpos, mb, ifail)
C-test-
      CALL UF420( mpos, fbm, 1, -999.0d0, ifail )
      if( ifail .gt. nx170 ) ifail=-42699
      if( debug(3) )write(*,1040)radius,mpeq%coord%colat,tlon,
     :              mb%dnrm,ifail,
     :              mlin%equat%coord%radius,mlin%equat%coord%colat,
     :              mlin%equat%coord%elong,mlin%equat%b%dnrm
C      
      mpeq         = mlin%equat
C      
      if(     ifail .eq. -53301 .or.
     :        ifail .eq. -53002 .or.
     :        ifail .eq. -53003 .or.
     :        ifail .eq. -53102 .or.
     :        ifail .eq. -53201 .or.
     :        ifail .eq. -42305     ) then
C MK, 12-feb-2004: trap also -53201 (15RE limit from O&P)
        kc1=kc1+1
        cmax       = MIN( radius, cmax)
        radius     = ( MAX( radius * 0.6d0, cmin) + radius )
     :                            * 0.5d0
        if( ABS(mpos%radius - radius) .ge. epskm 
     :      .and. radius .lt. cmax )then
          goto 100
        else
          ifail    = -41101
          if(debug(4))write(*,1020)ifail,cmin,cmax
          return
        endif
      elseif( ifail .eq. -53004 .or.
     :        ifail .eq. -53001 .or.
     :        ifail .eq. -42010 .or.
     :        ifail .eq. -42008 .or.
     :        ifail .eq. -42304     ) then
        kc2=kc2+1
        cmin       = MAX( radius, cmin)
        radius     = ( MIN( radius * 1.4d0, cmax) + radius )
     :                             * 0.5d0
        if( ABS(mpos%radius - radius) .ge. epskm 
     :      .and. radius .gt. cmin ) then
          goto 100
        else
          ifail    = -41102
          if(debug(4))write(*,1020)ifail,cmin,cmax
          return
        endif
      elseif( ifail .lt. 0  ) then
        ifl        = -41110
        if( ifail .eq. -42009 .or. ifail .eq. -42301 ) ifl = -41104
        if( ifail .eq. -42301 .or. ifail .eq. -42401 ) ifl = -41106
        if( ifail .eq. -42306 ) ifl = -41109
        if( ifail/100.ge.-428 .and. ifail/100.le.-426) ifl = -41108
        if(debug(6).or.debug(4))write(*,1030)'UF420',ifail
c-29jan04-        if( ifl.eq.-41110 )print*,'UF411 -> UF420',ifail
        ifail      = ifl
        return
      elseif( mpeq%b%dnrm .gt. fbm ) then
c        print*,'Remark: contact M. Kruglanski'
        kc3=kc3+1
        cmin       = MAX( radius, cmin )
        radius     = ( MIN( radius * 1.2d0, cmax) + radius )
     :                           * 0.5d0
        if(  ABS(mpos%radius - radius) .ge. epskm 
     :      .and. radius .gt. cmin )then
          mpeq%coord%colat = mpos%colat
          goto 100
        else
          ifail    = -41103
          return
        endif
      elseif( ABS(MOD((mpeq%coord%elong - glon)+540.0D0,
     :                               360.0D0)-180.0d0) .ge. epslon )then
         kc4=kc4+1
         tlon=glon+0.5d0*(MOD((tlon-mpeq%coord%elong)+540.0D0,
     :                               360.0D0)-180.0d0)
         radius=max(rmin,min(rmax,mpeq%coord%radius))
c         radius=max(rmin,min(rmax,mpeq%coord%radius*0.9d0+radius*0.1d0))
         cmin=rmin
         cmax=rmax
        goto 100
      endif
c
      radius = mpeq%coord%radius
C
      CALL UL230(finv, ifail)
      if( ifail .lt. 0 ) then
        if(debug(2))write(*,1030)'UL230',ifail
        ifail      = -41105
        return
      endif
      ml=mlab
      ifail=kcount
C
 1000 format(7x,': UF411 (',f6.1,',',f8.1,',',f8.1,',',f8.1,',',
     :          f6.1,',',f7.4,')')  
 1010 format(7x,': UF411 aborted after ',i2,'+',i2,'+',i2,'+',i2,
     :           ' iterations')    
 1020 format(7x,': in UF411, UF420 returned with IFAIL=',i6,' and',
     :       /7x,': Rmin =',f8.1,' Rmax =',f8.1)   
 1030 format(7x,': in UF411, ',a5,' returned with IFAIL=',i6)
 1040 format(7x,': in UF411, call to UF420 with mpos =',f9.1,f7.2,f7.2,
     :      ' where B =',1p,e12.5,0p,
     :      /19x,'returned with IFAIL =',i7,', equat =',f9.1,
     :      f7.2,f7.2,' and Beq =',f6.4)
C
      END
C----------------------------------------------------------------------
      SUBROUTINE UF415 (ifail)
C
C!    Rebuild the labels of the field line
C
      INCLUDE 'structure.h'
C
C
C     INTERFACE
C
      integer*4   ifail
C                      
      COMMON /UC170/  nsg, kgp, mlab, mlin, mele
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
C     variables
C
      real*8 bmpn,bmps,finv
C
C     Code    
C
      mlab%LINV=.false.
      mlab%LBMP=.false.
      mlab%LKAUF=.false.
      mlab%LLMI=.false.
      mlab%LALP0=.false.
      mlab%LPHI=.false.
      mlab%LTIM=.false.
C
      if(mlin%ind%jmirpn.le.0 .or. mlin%ind%jmirps.le.0)then
        ifail=-41502
c       internal error
        return
      endif
      bmpn=mele(mlin%ind%jmirpn)%beg%b%dnrm
      bmps=mele(mlin%ind%jmirps)%beg%b%dnrm
c MK, June 15, 1998: 1 line modified
c      mlab%fbmp=(bmpn+bmps)*0.5d0
      mlab%fbmp=max(bmpn,bmps)
c      print*,'uf415',mlin%ind%jmirpn,bmpn,mlin%ind%jmirps,bmps
      if( abs(bmpn-mlab%fbmp) .gt. epsrel*mlab%fbmp )then
        ifail=-41501
        return
      endif
      mlab%LBMP=.true.
C
      call ul230(finv,ifail)
      if( ifail.lt.0 )return   
C
      call ul240(mlab,ifail)
      if( ifail.lt.0 )return
C
      call ul245( mlin,mlab,ifail)
      if( ifail.lt.0 )return
C
      end         
C----------------------------------------------------------------------
      SUBROUTINE uf417
     :          (fim, fbm, glon, guess, falt, ifail)
C
C!    Search the geographic position of a magnetic field line segment
C
C     = UF410 of UNILIB where fim instead of flm as input, 
C                             without mpos as output,     
C
C_    er05
C
      INCLUDE 'structure.h'
C
C
      EXTERNAL UF411, UF420
C      
C     INTERFACE
C
        REAL*8        fbm, glon, guess, falt
        INTEGER*4     ifail
        TYPE(zgeo) :: mpos    
C
      COMMON /UC160
     :               /pi, deg, re, gmagmo, eclipt, geoid, uma
        REAL*8        pi, deg, re, gmagmo, eclipt, geoid(3), uma(30)
C
      COMMON /UC170/ nsg, kgp, mlab, mlin, mele
      INTEGER*4      nsg, kgp
      TYPE(zlbl) ::  mlab
      TYPE(zfln) ::  mlin
      TYPE(zseg) ::  mele(nx170) 
C           
      COMMON /UC190
     :               /prop, stepx, stpmin, umsq, upsq, uk2, uk3,
     :                epskm, epsrel, stplst, xclat, kmflg, nxstp
C
        REAL*8        prop, stepx, stpmin, umsq, upsq, uk2, uk3,
     :                epskm, epsrel, stplst, xclat
        INTEGER*4     kmflg, nxstp                                
C              
      COMMON /UC192
     :               /xrmin, xbmin, xtmin, xbmax, epslon, epsfl,
     :                fvet, pvet, epsomeg, dltlat
        REAL*8        xrmin, xbmin, xtmin, xbmax, epslon, epsfl, 
     :                fvet, pvet, epsomeg, dltlat
C                                                   
C     VARIABLES
C
        TYPE(zlbl) :: ml(5)
        TYPE(zpnt) :: mpeq(5)
        REAL*8        rgmin, rgmax, colat, fid, fiu, fim 
        INTEGER*4     k, k1, k2, k3, kcnt
c
c
        real*8        faltmx,fact
C
C     CODE
C
C
      if( fim .lt. 0.0d0 .or. fbm .gt. xbmax .or. fbm .lt. xbmin )then
        ifail  =-41705
        return
      endif
C
C
      fid      = (fim+re)*(1.0d0-epsfl)-re
      fiu      = (fim+re)*(1.0d0+epsfl)-re
 1000 format(' ****'/' UF417 fim=',f8.1,' fbm=',f6.4,' glon=',f8.1,
     :       ' guess=',f8.1,/'       fid=',f8.1,' fiu=',f8.1)
c
C
      kcnt     = 0
      ifail    = 0
      if( guess .le. xrmin ) guess = 2 * re
      rgmin    = 0.9D0* re
      rgmax    = 15.0D0 * re
C
C The point at the equator, along the magnetic field line
C must have its longitude equal to glon and its radius
C between rgmin and rgmax
C
C In the first version of this subroutine, the strategy was based
C on the value of L, now it is based on the value of I. The algorithm
C is the same
C
C The subroutine uf411 search a magnetic field line such as
C the point with the lowest magnetic field intensity along the
C line has a longitude equal to glon and a magn. field intensity
C less than fbm
C
C
      k1       = 1
      k2       = 2
      k3       = 3
C
C  a/ First attempt
C     -------------
C
      ifail=0
      CALL UF411( glon, guess, rgmin, rgmax, -1.0d0, fbm,
     :            mpeq(k1), ml(k1), ifail )
      if( ifail.lt.0 ) return
      kcnt      = kcnt + ifail
      colat     = mpeq(k1)%coord%colat
C
 1070 format(' UF417 first attempt in mpeq=',3f9.1,' inv=',f8.1)
C            
C  b/ Second attempt
C     --------------
C
      if( ml(k1)%finv .gt. fiu )then
C      
C  b.1/ Lm1 > Lm --> Alt2 < Alt1 & Lm2 < Lm
C
        rgmax   = guess
        guess   = (re+rgmax)/2
  210   continue
        ifail=0
        CALL UF411( glon, guess, rgmin, rgmax, colat, fbm,
     :              mpeq(k2), ml(k2), ifail)
        if( ifail.lt.0 ) then
c-29jan04-           print*,' B.1 UF411 mpeq',mpeq(k2)%coord%radius,
c-29jan04-     :          mpeq(k2)%coord%colat,
c-29jan04-     :          mpeq(k2)%coord%elong,' ifail=',ifail

           return
        endif
        kcnt    = kcnt + ifail
cccccc D. Heynderickx   Jun 2014
        IF (kcnt .GT. 250) THEN
          IFAIL = -41707
          RETURN
        END IF
cccccc
C                 
        if( ml(k2)%finv .gt. fiu )then
          rgmax = guess
          if( guess .gt. 6800. ) then
            guess=(guess+2.0D0*re)/3.0D0
          else
            guess =MAX( rgmin, guess * 0.65D0)
          endif
          if( ABS( guess - mpeq(k2)%coord%radius )
     :                               .lt. epskm )then
            ifail = -41701
            return
          endif
          goto 210
        else if( ml(k2)%finv .ge. fid )then
C     ... second attempt succeeded...
C
          mpos    = mpeq(k2)%coord
          ifail   = kcnt
 1010 format(' UF417 attempt (',a3,') succeeded after',i4)
          goto 500       
        endif
        k       = k1
        k1      = k2
        k2      = k
      else if( ml(k1)%finv .lt. fid )then
C      
C  b.2/ Lm1 < Lm --> Alt2 > Alt1 & Lm2 > Lm
C
        rgmin   = guess
cccccccc        guess   = MIN( rgmax, guess*2 )
        guess   = (rgmax+4*guess)/5.0D0
  220   continue
        ifail=0
        CALL UF411( glon, guess, rgmin, rgmax, colat, fbm,
     :              mpeq(k2), ml(k2), ifail)
        if( ifail.lt.0 ) return
        kcnt    = kcnt + ifail
cccccc D. Heynderickx   Mar 2015
        IF (kcnt .GT. 250) THEN
          IFAIL = -41708
          RETURN
        END IF
cccccc
C                 
        if( ml(k2)%finv .lt. fid )then
          rgmin = guess
          guess = ( MIN( rgmax, guess * 1.6d0) + guess ) * 0.5d0
          if( ABS( guess - mpeq(k2)%coord%radius ) 
     :                              .lt. epskm )then
            ifail = -41702
            return
          endif
          goto 220
        else if( ml(k2)%finv .le. fiu )then
C     ... second attempt succeeded...
C
          mpos    = mpeq(k2)%coord
          ifail   = kcnt
          goto 500       
        endif
      else
C      
C  b.3/ First attempt succeeded...
C
        mpos    = mpeq(k1)%coord
        ifail   = kcnt
        goto 500  
      endif  
 1040 format(' UF417 first pass in',i4,/'   LO=',i1,' mpeq=',3f9.1,
     :       ' inv=',f8.1,/'   HI=',i1,' mpeq=',3f9.1,
     :       ' inv=',f8.1)
C              
C  c/ Iterative search
C     ----------------
C
C  ===>  Lm1 < Lm < Lm2
C
  300 continue
      rgmin     = mpeq(k1)%coord%radius
      rgmax     = mpeq(k2)%coord%radius
      colat     = ( mpeq(k1)%coord%colat 
     :            + mpeq(k2)%coord%colat ) * 0.5d0
      fact      = ( fim-ml(k1)%finv ) / ( ml(k2)%finv-ml(k1)%finv )
      if( fact .lt. 0.1D0 )then
         fact=0.15D0
      endif
      if( fact .gt. 0.9D0 )then
         fact=0.85D0
      endif
 1060 format(7x,'modif fact from',f8.3,' to',f8.3)
      guess     = rgmin + ( rgmax - rgmin ) * fact

      ifail=0
 1020 format(' UF417 -> UF411 guess=',f8.1,' rgmin=',f8.1,
     : ' rgmax=',f9.1,' colat=',F7.1)
      CALL UF411( glon, guess, rgmin, rgmax, colat, fbm,
     :            mpeq(k3), ml(k3), ifail)
      if( ifail.lt.0 ) return
 1030 format('       <- UF411 mpeq=',3f9.1,' inv=',f8.1,' ifail=',i7)
      kcnt      = kcnt + ifail
C
C                 
      if( kcnt .gt. 100 )then
        kcnt    = -41703
      elseif( ml(k3)%finv .gt. fiu )then
C      
C  c.1/ too high
C
        k       = k3
        k3      = k2
        k2      = k
 1050 format(7x,'NW=',i1,a8,' LO=',i1,' HI=',i1)
        goto 300
      elseif( ml(k3)%finv .lt. fid )then
C      
C  c.2/ too low
C
        k       = k3
        k3      = k1
        k1      = k
        goto 300
      endif
C
      mpos      = mpeq(k3)%coord
      ifail     = kcnt
      if( ifail .lt. 0 )return
C
 500  continue
c
      if( falt .ge. -500.0d0 )then
        if( mlin%equat%coord%radius-re .lt. falt )then
          ifail = -41706
c          print*,'equat',mlin%equat%coord%radius-re, falt
        else
          CALL UF420(mpos, fbm, 1, falt, k)
          if( k .lt. 0 ) then
c-29jan04-            print*,'***',mlin%equat%coord%radius-re, falt
            CALL UF420(mpos, fbm, 1, faltmx, ifail)
            if( ifail .ge. 0 )then
              ifail = -41704
            else
c-29jan04-              print*,'UF410: UF420-> (UF428) ifail=',k
c-29jan04-              print*,'       -> UF420[-999] -> ifail=',ifail
              return
            endif
          endif
        endif
      endif    
C
      call uf415(k)
      if( k .lt. 0 )then
        if( k.ne.-24504)ifail=k
      endif
      if( ifail.lt.0 ) return    
C
C
      END
C----------------------------------------------------------------------
