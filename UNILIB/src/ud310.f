# 1 "ud310.f"
       SUBROUTINE UD310
     :          (fbm0, flm0, falt, knfl, ktyplus, ifail)
C
C!    trace a drift shell
C_    er13
C
      Include 'structure.h'
cDEC$ IF DEFINED (_x86_)
cDEC$ ATTRIBUTES DLLEXPORT :: UD310
cDEC$ ENDIF
C
C     INTERFACE
      REAL*8    fbm0, flm0, falt
      INTEGER*4 knfl, ifail, ktyplus
C
      EXTERNAL uf410, ud319
C
C     VARIABLES
      integer*4 nx315
      parameter( nx315=40)
      REAL*8    glon0, dlon, guess, gl, dlon2
      INTEGER*4 nfl1, nfl2, id(5), i, k, kind, j
      TYPE(ZGEO) :: mgpos, mlist(nx315)
      integer*4 klist(nx315),kdebug
      real*8    glist(nx315)
      integer*4 ktyp,k41004
C
C
      k41004=0
      kdebug=   0
      if( ifail.lt.0 ) kdebug=ifail
      ifail =   0
      glon0 = -90.0d0
      guess =  -1.0d0
      kind  =  -1
      ktyp  = abs(ktyplus)
C
C     1/ select nfl1 and nfl2
C
      if (knfl .le. 0) knfl=120
      do i=1,5
         id(i)=1
      enddo
      if ( mod(knfl,2)  .eq. 0 ) id(1) =  2
      if ( mod(knfl,3)  .eq. 0 ) id(2) =  3
      if ( mod(knfl,5)  .eq. 0 ) id(3) =  5
      if ( mod(knfl,7)  .eq. 0 ) id(4) =  7
      if ( mod(knfl,11) .eq. 0 ) id(5) = 11
      nfl1 = max( id(1)*id(2)*id(3), id(3)*id(4), id(2)*id(5),
     :            id(1)*id(5), id(4)*id(2), id(4)*id(1) )
      nfl2 = knfl / nfl1
C
C     2/ first loop
C
      if( ktyplus .lt. 0 ) then
          print 1000, fbm0, flm0, knfl,nfl1,nfl2
          print 1020
      endif
      dlon = 360.0d0 / nfl1
      do k = 1, nfl1
        gl     = glon0 + (k-1) * dlon
        if( ktyplus .lt. 0 ) print 1010,gl
        ifail=kdebug
        call UF410 (flm0, fbm0, gl, guess, mgpos, falt, ifail)
        if( ifail .eq. -41004 )then
          k41004 = ifail
          ifail  = 0
        elseif( ifail .lt. 0 ) then
          return
        endif                   
c        print*
c        print*,'ici',mgpos.elong,mgpos.colat,mgpos.radius,ifail
        call UD319(flm0, fbm0, kind, ktyp, ifail)
        if( ifail .lt. 0 ) return
        klist(k) = kind
        mlist(k) = mgpos
        glist(k) = guess
        if( ktyplus .lt. 0 .and. mod(k,10) .eq. 0 ) print 1020
      enddo
      if( ktyplus .lt. 0 .and. mod(nfl1,10) .ne. 0 ) print 1020
      mlist(nfl1+1) = mlist(1)
      glist(nfl1+1) = glist(1)
C
c        3/ second loop
c
      j = 1
      do k = 1, nfl1
        dlon2 = ( mlist(k+1)%elong - mlist(k)%elong ) / nfl2
        if( dlon2 .le. 0 )then
          dlon2 = ( 360.0d0+mlist(k+1)%elong - mlist(k)%elong ) / nfl2 
        endif
        guess = ( glist(k)+glist(k+1)*(nfl2-1) ) / nfl2
        kind  = klist(k)
        do i = 2, nfl2
          gl   = mlist(k)%elong + (i-1) * dlon2
          if( ktyplus .lt. 0 ) write(*,1010) gl
          ifail=kdebug
          call UF410 (flm0, fbm0, gl, guess, mgpos, falt, ifail)
        if( ifail .eq. -41004 )then
          k41004 = ifail
          ifail  = 0
        elseif( ifail .lt. 0 ) then
          return
        endif                   
c          print*
c          print*,'la',mgpos%elong,mgpos%colat,mgpos%radius,ifail
          call UD319(flm0, fbm0, kind, ktyp, ifail)
          if( ifail .lt. 0 ) return
          if( ktyplus .lt. 0 .and. mod(j,10) .eq. 0 ) print 1020
          j    = j + 1
        enddo
      enddo
      if( ktyplus .lt. 0 ) print 1030
C
      ifail=klist(1)
      if( k41004 .ne. 0 )ifail=k41004
C
 1000 format(' (UD310) Trace the drift shell (',f6.4,',',f5.2,
     :       ') using',i4,' =',i3,' x',i3,' field lines.')
 1010 format('+',f6.1,';',$)
 1020 format('        ',$)
 1030 format('+-done-')
      END
C----------------------------------------------------------------------
      subroutine ud315(mpn, mps, ifail)
C     Search the mirror points with the lowest altitude in both 
C     hemispheres
      INCLUDE 'structure.h'
      TYPE(zgeo) :: mpn, mps
      integer*4 ifail
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
      COMMON /UC130/  nbrsg, kursg, mseg
C
        INTEGER*4     nbrsg, kursg
        TYPE(zseg) :: mseg(nx130)
C        
c
      integer*4 k, kpn, kps, i, jn, js
      integer*4 jprev,jnext
      real*8    estradius
c
      ifail=0
c
      if( k1st .lt. 0 .or. k1st .gt. nbrfl )then
        ifail=-31501
        return
      endif
      i=nbrfl
      k=k1st
      kpn=mfl(k)%ind%jmirpn
      kps=mfl(k)%ind%jmirps
      mpn=mseg(kpn)%beg%coord
      mps=mseg(kps)%beg%coord
      jn = k
      js = k
  100 continue
        k=mfl(k)%keast
        i=i-1
        if( i .lt. 0 )then
          ifail=-31502
          return
        endif
        kpn=mfl(k)%ind%jmirpn
        kps=mfl(k)%ind%jmirps
        if( mpn%radius .gt. mseg(kpn)%beg%coord%radius )then
          mpn=mseg(kpn)%beg%coord
          jn=k
        endif   
        if( mps%radius .gt. mseg(kps)%beg%coord%radius )then
          mps=mseg(kps)%beg%coord
          js=k
        endif
      if( k .ne. k1st )goto 100
c
      if( mps%radius .lt. mpn%radius )then
        klmp = js
      else
        klmp = jn
      endif
c
      jnext=mfl(js)%keast
      jprev=mfl(js)%kwest
      call ud327(mseg(mfl(jprev)%ind%jmirps)%beg%coord,
     :           mseg(mfl(js)%ind%jmirps)%beg%coord,
     :           mseg(mfl(jnext)%ind%jmirps)%beg%coord,
     :           mseg(mfl(jprev)%ind%jmirps)%beg%coord%radius,
     :           mseg(mfl(js)%ind%jmirps)%beg%coord%radius,
     :           mseg(mfl(jnext)%ind%jmirps)%beg%coord%radius,
     :           mps,estradius,ifail)
      mps%radius=estradius
c
      jnext=mfl(jn)%keast
      jprev=mfl(jn)%kwest
      call ud327(mseg(mfl(jprev)%ind%jmirpn)%beg%coord,
     :           mseg(mfl(jn)%ind%jmirpn)%beg%coord,
     :           mseg(mfl(jnext)%ind%jmirpn)%beg%coord,
     :           mseg(mfl(jprev)%ind%jmirpn)%beg%coord%radius,
     :           mseg(mfl(jn)%ind%jmirpn)%beg%coord%radius,
     :           mseg(mfl(jnext)%ind%jmirpn)%beg%coord%radius,
     :           mpn,estradius,i)
      mpn%radius=estradius
      if(i.ne.0)ifail=i
c
c
c
      end
C----------------------------------------------------------------------
      subroutine ud316(mequat, ifail)
C     Search the equatorial point with the lowest magnetic field 
C     intensity
      INCLUDE 'structure.h'
      TYPE(zpnt) :: mequat
      integer*4 ifail
*
* mequat = equatorial point with the lowest B
* ifail  = error flag
*
C  The subroutine UD316 scans the field line segments of a drift shell and
c  determines the point with the lowest magnetic field intensity.
c  The drift shell has to be
c  defined as a set of field line segments stored in the common blocks UC110, 
c  UC120, and UC130. 
c diagnostics: see ud315
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
C        
c
      integer*4 k, i
c
      ifail=0
c
      if( k1st .lt. 0 .or. k1st .gt. nbrfl )then
        ifail=-31601
        return
      endif
      i=nbrfl
      k=k1st
      mequat=mfl(k)%equat
      klmp=k
  100 continue
        k=mfl(k)%keast
        i=i-1
        if( i .lt. 0 )then
          ifail=-31502
          return
        endif
        if( mfl(k)%equat%b%dnrm  .lt. mequat%b%dnrm )then
          mequat=mfl(k)%equat
          klmp=k
        endif
      if( k .ne. k1st )goto 100
c
      ifail=klmp
c
      end
C----------------------------------------------------------------------
       SUBROUTINE UD317
     :          (mlab0, falt, ktyplus, altmin, ifail)
C
C!    trace a drift shell
C     = UD310 of UNILIB where mlab0 instead of fbm0 and flm0 as input, 
C                             without knfl as input (fixed to 120 + ...),
C                             I is used instead of L,
C                             altmin as output (lowest MP alt.)    
C_    er13
C
      Include 'structure.h'
cDEC$ IF DEFINED (_x86_)
cDEC$ ATTRIBUTES DLLEXPORT :: UD317
cDEC$ ENDIF
C
C     INTERFACE
      TYPE(ZLBL) :: mlab0
      REAL*8    fbm0, fim0, falt, flm0, altmin
      INTEGER*4 ifail, ktyplus
C
      EXTERNAL uf410, ud319
C
      COMMON /UC120/ nbrfl, kurfl, mfl
        INTEGER*4     nbrfl, kurfl
        TYPE(zfln) :: mfl(nx120) 
      COMMON /UC160/ pi, deg, re, gmagmo, eclipt, geoid, uma
        REAL*8        pi, deg, re
        REAL*8        gmagmo
        REAL*8        eclipt, geoid(3), uma(30)
      COMMON /UC170/  nsg, kgp, mlab, mlin, mele
        INTEGER*4     nsg, kgp
        TYPE(zlbl) :: mlab
        TYPE(zfln) :: mlin
        TYPE(zseg) :: mele(nx170)
C
C
C     VARIABLES
      integer*4 nx315
      parameter( nx315=40)
      REAL*8    glon0, dlon, guess, gl
      INTEGER*4 nfl1, i
      integer*4 k, kind
      real*8    altmnn, altmns, x, altmx
      integer*4 kptmnn, kptmns, klist(6), kptmx, karousel 
      integer*4 ktyp,k41004
      logical   fdebug
C
C
      fdebug = ktyplus .lt. 0
      altmin=-999.0d0
      altmns=1000*re
      altmnn=1000*re
      altmx=-re
      ifail =   0
      if( .not. mlab0%linv )then 
        if( fdebug ) ifail = -999
        call ul242(mlab0,ifail)
        if( ifail.lt.0 ) return
      elseif( .not. mlab0%lbmp )then 
        if( fdebug ) print*,'Unilib-UD317: B missing'
        ifail=-31702
        return
      endif
      if( fdebug ) print*,'Unilib-UD317: B,I = ',mlab0%fbmp,mlab0%finv
c
      fbm0 = mlab0%fbmp
      fim0 = mlab0%finv
      if( .not. mlab0%llmi )then 
        call ul240(mlab0,ifail)
        if( ifail.lt.0 ) return
      endif
      flm0 = mlab0%flmi
c
      k41004=0
      ifail =   0
      glon0 = -60.0d0
cccccc
      guess =  re+flm0*re
      kind  =  -1
      ktyp  = abs(ktyplus)
C
C     2/ first loop
C
      nfl1 = 90
      dlon = 360.0d0 / nfl1
      do k = 1, nfl1
        gl     = glon0 + (k-1) * dlon
        ifail=0
        if( fdebug ) print*,'Unilib-UD317> UF417 ',gl,guess
        call UF417 (fim0, fbm0, gl, guess, falt, ifail)
        if( ifail .eq. -41704 .or. ifail .eq. -41706)then
c          print*,'UF417 -> -41704'
          k41004 = k41004+1
          ifail  = 0
        elseif( ifail .lt. 0 ) then
          return
        endif                   
        if(.not. mlab%linv)then
          print*,'oups',k,ifail
        endif
        if( fdebug ) print*,'Unilib-UD317> UD319 ',kind,ktyp
        call UD319(flm0, fbm0, kind, ktyp, ifail)
        if( ifail .lt. 0 ) return
        i = mlin%ind%jmirps
        x = mele(i)%beg%coord%radius-re
        if( altmns .gt. x ) then
           kptmns=kind
           altmns=x
        endif
        i = mlin%ind%jmirpn
        x = mele(i)%beg%coord%radius-re
        if( altmnn .gt. x ) then
           kptmnn=kind
           altmnn=x
        endif
        x = guess-re
        if( altmx.lt.x ) then
           kptmx=kind
           altmx=x
        endif
      enddo
C
C     3/
C
      do karousel=1,3
        if( fdebug )print*,'Unilib-UD317: Iteration ',karousel
        klist(1)=mfl(kptmnn)%kwest
        klist(2)=kptmnn
        klist(3)=mfl(kptmns)%kwest
        if( klist(3).eq.klist(1) .or. klist(3).eq.klist(2) )klist(3)=-1
        klist(4)=kptmns
        if( klist(4).eq.klist(1) .or. klist(4).eq.klist(2) )klist(4)=-1
        klist(5)=mfl(kptmx)%kwest
        if( klist(5).eq.klist(1) .or. klist(5).eq.klist(2) .or.
     :    klist(5).eq.klist(3) .or. klist(5).eq.klist(4) )klist(5)=-1
        klist(6)=kptmx
        if( klist(6).eq.klist(1) .or. klist(6).eq.klist(2) .or.
     :    klist(6).eq.klist(3) .or. klist(6).eq.klist(4) )klist(6)=-1
ccc        print*,karousel,'(karousel)',klist
C
        dlon=dlon/2
        do k=1,6
          kind=klist(k)
          if( kind.gt.0 )then
            gl=mfl(kind)%equat%coord%elong+dlon
            guess=mfl(kind)%equat%coord%radius
            ifail=0
            if( fdebug ) print*,'Unilib-UD317> UF417 ',gl,guess
            call UF417 (fim0, fbm0, gl, guess, falt, ifail)
            if( ifail .eq. -41704 .or. ifail.eq. -41706)then
              k41004 = k41004+1
              ifail  = 0
            elseif( ifail .lt. 0 ) then
              return
            endif                   
            if( fdebug ) print*,'Unilib-UD317> UD319 ',kind,ktyp
            call UD319(flm0, fbm0, kind, ktyp, ifail)
            if( ifail .lt. 0 ) return
            i = mlin%ind%jmirps
            x = mele(i)%beg%coord%radius-re
            if( altmns .gt. x ) then
               kptmns=kind
               altmns=x
            endif
            i = mlin%ind%jmirpn
            x = mele(i)%beg%coord%radius-re
            if( altmnn .gt. x ) then 
               kptmnn=kind
               altmnn=x
            endif
            x = guess-re
            if( altmx.lt.x ) then
               kptmx=kind
               altmx=x
            endif
          endif
        enddo
      enddo
C
C
      altmin=min(altmnn,altmns)
C
C
      if( k41004 .ne. 0 )then
        if( fdebug ) print*,'Unilib-UD317: -31701 stat = ',k41004
        ifail=-31701
      endif
      end
C----------------------------------------------------------------------
      SUBROUTINE UD319
     :           (flm, fbm, kind, ktyp, ifail)
C
C     Transfer a field line segment from UC170 to UC130 
C
      INCLUDE 'structure.h'
C                       
C     INTERFACE
C
        INTEGER*4     kind, ifail, ktyp
        real*8        flm, fbm
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
      COMMON /UC130/  nbrsg, kursg, mseg
C
        INTEGER*4     nbrsg, kursg
        TYPE(zseg) :: mseg(nx130)
C        
      COMMON /UC160/ pi, deg, re, gmagmo, eclipt, geoid, uma
C
        REAL*8        pi, deg, re
        REAL*8        gmagmo
        REAL*8        eclipt, geoid(3), uma(30)
C       
      COMMON /UC170/  nsg, kgp, mlab, mlin, mele
C
        INTEGER*4     nsg, kgp
        TYPE(zlbl) :: mlab
        TYPE(zfln) :: mlin
        TYPE(zseg) :: mele(nx170)
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
C    *
C    * KTYP   in   1 = transfer the whole field line segment
C    *             2 = transfer the field line segment comprised between
C    *                 the two mirror points (preferred)
C    *             3 = only transfer the segment description
C    *
C    * KIND   in  -1 = new drift shell
C    *        in  >0 = index of the field line to be on the west
C    *        out      index of the current field line in UC120
C    *
C
C     VARIABLES
C
        INTEGER*4     jbeg, jend, j, kwest, keast
C
C     CODE
C
      ifail                 = 0
C      
      if (ktyp.lt.1 .or. ktyp.gt.3) then
        ifail               = -31907
        return
      endif
C
      if ( nsg .le. 0  .or.
     :     .not. mlab%lbmp  .or.
     :     .not. mlab%linv  ) then
        ifail               = -31901
        return
      endif       
C
      if (kind .lt. 0 ) then
        k1st                = 1
        klmp                = -1
        mlbl                = mlab
        mlbl%fbmp           = fbm
        mlbl%flmi           = flm
        mlbl%lphi           = .FALSE.
        mlbl%ltim           = .FALSE.
        nbrfl               = 0
        kurfl               = 0
        nbrsg               = 0
        kursg               = 0
        kwest               = 1
        keast               = 1
      else if (kind .gt. nbrfl)then
        ifail               = -31908
        return
      else if (.not. mlbl%lbmp  .or.
     :         .not. mlbl%linv ) then
        ifail               = -31902
        return
      else if ( ABS(mlab%fbmp - mlbl%fbmp) .gt.
     :                  2*epsrel*mlbl%fbmp ) then
c        print*,'alors',mlab%fbmp,mlbl%fbmp,2*epsrel*mlbl%fbmp
        ifail               = -31903
        return
      else if ( 0.0d0*ABS(mlab%finv - mlbl%finv) .gt.               
     :  epsfl*(mlab%finv + mlbl%finv) ) then
c        print*,1.95*ABS(mlab%finv-mlbl%finv),epsfl*(mlab%finv+mlbl%finv)
        ifail               = -31904
        return
      else
        kwest               = kind
        keast               = mfl(kind)%keast
      endif
C        
      nbrfl                 = nbrfl + 1
      kurfl                 = nbrfl
      if ( kurfl .ge. nx120 ) then
        ifail               = -31905
        return
      endif
C
      mfl(kurfl)            = mlin
C
C
      if(     ktyp .eq. 1)then
C
        jbeg                  = mlin%ind%jbeg
        jend                  = mlin%ind%jend
        if ( nbrsg+jend-jbeg+1 .ge. nx130 ) then
          ifail               = -31906
          return
        endif            
C
        nbrsg                 = nbrsg + 1
        kursg                 = nbrsg
        mfl(kurfl)%ind%jbeg   = kursg
        mfl(kurfl)%ind%jmirpn = mlin%ind%jmirpn - mlin%ind%jbeg + kursg
        do j= jbeg, jend
          mseg(nbrsg)         = mele(j)
          nbrsg               = nbrsg + 1
        enddo
        nbrsg                 = nbrsg - 1
        mfl(kurfl)%ind%jend   = nbrsg
        mfl(kurfl)%ind%jmirps = mlin%ind%jmirps - mlin%ind%jend + nbrsg
C
      elseif( ktyp .eq. 2)then
C
        jbeg                  = mlin%ind%jmirpn
        jend                  = mlin%ind%jmirps
        if ( nbrsg+jend-jbeg+1 .ge. nx130 ) then
          ifail               = -31906
          return
        endif            
C
        nbrsg                 = nbrsg + 1
        kursg                 = nbrsg
        mfl(kurfl)%ind%jbeg   = kursg
        mfl(kurfl)%ind%jmirpn = kursg
        do j= jbeg, jend
          mseg(nbrsg)         = mele(j)
          nbrsg               = nbrsg + 1
        enddo
        nbrsg                 = nbrsg - 1
        do j=2,3
          mseg(nbrsg)%rkstp(j)= 0.0d0
        enddo
        mfl(kurfl)%ind%jend   = nbrsg
        mfl(kurfl)%ind%jmirps = nbrsg
C
      else
C
        jbeg                  = mlin%ind%jmirpn
        jend                  = mlin%ind%jmirps
        if ( nbrsg+2 .ge. nx130 ) then
          ifail               = -31906
          return
        endif            
C
        nbrsg                 = nbrsg + 1
        kursg                 = nbrsg
        mfl(kurfl)%ind%jbeg   = kursg
        mfl(kurfl)%ind%jmirpn = kursg
        mseg(nbrsg)           = mele(jbeg)
        nbrsg                 = nbrsg + 1
        mseg(nbrsg)           = mele(jend)
        mfl(kurfl)%ind%jend   = nbrsg
        mfl(kurfl)%ind%jmirps = nbrsg
        do j=2,3
          mseg(kursg)%rkstp(j)= 0.0d0
          mseg(nbrsg)%rkstp(j)= 0.0d0
        enddo
C
      endif
C
      mfl(kwest)%keast      = kurfl
      mfl(keast)%kwest      = kurfl
      mfl(kurfl)%keast      = keast
      mfl(kurfl)%kwest      = kwest
C
      kind                  = kurfl
      ifail                 = jend - jbeg + 1
C
C
      END
C----------------------------------------------------------------------