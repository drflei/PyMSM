# 1 "ud320.f"
      SUBROUTINE UD320
     :         (kpar1, rpar2, zsub, result, rnorm, rdrift, ifail)     
C
C!    evaluate a double time integral over a drift shell
C
      include 'structure.h'
c
c     INTERFACE
      REAL*8    rpar2, result(*), rnorm, rdrift
      INTEGER*4 kpar1, ifail, kdim, kdims
      EXTERNAL  zsub
C
*
*     kpar1   control parammeter passed to zsub
*     rpar2   control parameter passed to zsub
*     zsub    external subroutine
*     result  result of the integration
*     rnorm   time integral of the unity
*     rdrift  single time integral of the unity
*     ifail   error flag
*                                              
C
C     The synopsis of the subroutine zsub has to 
C     correspond to:
C           subroutine ZSUB (kpar1, rpar2, mpos, rval, ifail)
C           include 'structure.h'
C           record /zgeo/ mpos 
C           integer*4     kpar1, ifail
C           real*8        rpar2, rval
C     The input arguments kpar1, rpar2, mpos are two user-specific 
C     control parameters and the geocentric position where the
C     user function has to be evaluated. As result the argument
C     rval has to be return with the value of the user function.
C     If the argument ifail is set to a negative value by the
C     subroutine ZSUB, the integration will be automatically aborted.
C
C     The subroutine ua636 can be used to integrated
C     atmospheric densities over a drift shell. In that
C     case, the argument kpar1 and rpar2 corresponds to
C     the argument ktrl and eng of subroutine ua636.
C
c
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
      COMMON /UC130/ nbrsg, kursg, mseg
C
        INTEGER*4     nbrsg, kursg
        TYPE(zseg) :: mseg(nx130)
C
      COMMON /UC160
     :               /pi, deg, re, gmagmo, eclipt, geoid, uma
C
        REAL*8        pi, deg, re, gmagmo, eclipt, geoid(3), uma(30)
C
C         
C     VARIABLES     
C            
      REAL*8         rone,rline(300),bm,rn
      TYPE(zxyz) ::  mvlcr, mpsnx, mpscr, mdspr,mdscr
      integer*4      kprev,knext,i
      REAL*8         st,sp,ct,cp,fct,dnrm
C
C     CODE
C
      ifail=0
      rdrift=0
      rnorm=0.0d0
      kdim=1
ccccc      print *,'                      new_ud320'
      goto 17
c
      entry ud322 (kpar1, rpar2, zsub, result, kdims, 
     :             rnorm, rdrift, ifail)
      kdim=kdims
ccccc      print *,'                      entry ud322',kdim
c
   17 if(kdim.lt.1.or.kdim.gt.300)then
        ifail=-32201
        return
      endif
      do i=1,kdim
        result(i)=0.0d0
      enddo  
C
      if( k1st.le.0)then
        ifail=-32001
c       invalid drift shell description
        return                             
      endif
      kurfl=k1st
C
      bm = mlbl%fbmp
      if( bm.le.0 .or. .not.mlbl%lbmp)then
        ifail=-32002
c       invalid drift shell description
        return
      endif
C
      kprev=mfl(kurfl)%kwest
c
c     mpscr = GEO cartesian position in km at kurfl      
c     mpsnx = GEO cartesian position in km at knext
c     mdspr = cartesian distance between kprev and kurfl
c     mdscr = cartesian distance between kurfl and knext      
c     mvlcr = cartesian comp. of unit vector // to drift vel 
c
c     mfl%drift in km-1 Gauss-1 at kurfl
c
c    conversion to velocity: mass c^2/2/q beta^2 gamma (2-b/bm)
c    beta = v/c 
c    gamma = sqrt(1-beta^2)
c
c    note: here only (2-b/bm)
c
c
      st  = sin( mfl(kprev)%equat%coord%colat*deg )
      ct  = cos( mfl(kprev)%equat%coord%colat*deg )
      sp  = sin( mfl(kprev)%equat%coord%elong*deg )
      cp  = cos( mfl(kprev)%equat%coord%elong*deg )
      mdspr%z = ct    * mfl(kprev)%equat%coord%radius
      mdspr%x = st*cp * mfl(kprev)%equat%coord%radius
      mdspr%y = st*sp * mfl(kprev)%equat%coord%radius
c
      st  = sin( mfl(kurfl)%equat%coord%colat*deg )
      ct  = cos( mfl(kurfl)%equat%coord%colat*deg )
      sp  = sin( mfl(kurfl)%equat%coord%elong*deg )
      cp  = cos( mfl(kurfl)%equat%coord%elong*deg )
      mpscr%z = ct     * mfl(kurfl)%equat%coord%radius
      mpscr%x = st*cp * mfl(kurfl)%equat%coord%radius
      mpscr%y = st*sp * mfl(kurfl)%equat%coord%radius
c
      mdspr%z = mpscr%z - mdspr%z
      mdspr%x = mpscr%x - mdspr%x
      mdspr%y = mpscr%y - mdspr%y
c
  100 continue
ccccc      print*,'                      (100 continue)',bm
c
        fct = 2-mfl(kurfl)%equat%b%dnrm/bm
        mvlcr%z = mfl(kurfl)%drift%rho*ct - mfl(kurfl)%drift%theta*st
        mvlcr%z = mvlcr%z/mfl(kurfl)%drift%dnrm
        mvlcr%x =   mfl(kurfl)%drift%rho * st*cp
     :            + mfl(kurfl)%drift%theta * ct*cp
     :            - mfl(kurfl)%drift%phi * sp
        mvlcr%x = mvlcr%x/mfl(kurfl)%drift%dnrm
        mvlcr%y =   mfl(kurfl)%drift%rho * st*sp
     :            + mfl(kurfl)%drift%theta * ct*sp
     :            + mfl(kurfl)%drift%phi * cp
        mvlcr%y = mvlcr%y/mfl(kurfl)%drift%dnrm
        dnrm    = mfl(kurfl)%drift%dnrm * fct
ccccc        print*,'                      dnrm',mfl(kurfl)%drift%dnrm,fct
c
        knext=mfl(kurfl)%keast
c
        st  = sin( mfl(knext)%equat%coord%colat*deg )
        ct  = cos( mfl(knext)%equat%coord%colat*deg )
        sp  = sin( mfl(knext)%equat%coord%elong*deg )
        cp  = cos( mfl(knext)%equat%coord%elong*deg )
        mpsnx%z = ct    * mfl(knext)%equat%coord%radius
        mpsnx%x = st*cp * mfl(knext)%equat%coord%radius
        mpsnx%y = st*sp * mfl(knext)%equat%coord%radius
c
        mdscr%z = mpsnx%z - mpscr%z
        mdscr%x = mpsnx%x - mpscr%x
        mdscr%y = mpsnx%y - mpscr%y
c
        rn=(   (mdspr%z+mdscr%z)*mvlcr%z
     :       + (mdspr%x+mdscr%x)*mvlcr%x
     :       + (mdspr%y+mdscr%y)*mvlcr%y ) / 2 / dnrm
c
c       unit: Gauss
c
ccccc      print*,'                      --> UD323',kdim
        call ud323(mfl(kurfl)%ind,mseg, 
     :             kpar1, rpar2, zsub, rline, kdim, rone, ifail)
        if(ifail.lt.0)return
c
c  mfl%drift (/zvec/) set by uf424 (arg mxnb)
C  uf424 called by uf420 called by uf411 called by uf410
c        called by ud310: "trace a magnetic drift shell" 
c  mfl%drift = mxnb 
c            = Vector pointing to the binormal vector of the 
c              magnetic field line at EQUAT
c
c  see JS Lew (JGR 66,1961, p2681-2685) eq (3)
c
c  mxnb = e x n /rcurv /b
c
c  Drift velocity at EQUAT is prop. to mxnb (2-B/Bm)
c
c
c
        rdrift = rdrift+rn
        rnorm  = rnorm+rone*rn
ccccc        print*,'            ici?'
        do i=1,kdim
ccccc          print*,i,result(i),rline(i),rn
          result(i) = result(i)+rline(i)*rn
        enddo
ccccc        print*,'            ou la?'
c
c
        mpscr=mpsnx
        mdspr=mdscr
        kprev=kurfl
        kurfl=knext
      if(kurfl.ne.k1st)goto 100
C
      end       
C----------------------------------------------------------------------
      SUBROUTINE UD321
     :         (mind, mele, kpar1, rpar2, zsub, result, rnorm, ifail)   
C
C!    evaluate a time integral over a magnetic field line
C
c     ud321 has been checked by MK on October 8, 1998
c     the definition of zsub has been changed on the same date
c
      include 'structure.h'
c
c     INTERFACE
      REAL*8    rpar2, result(*), rnorm
      INTEGER*4 kpar1, ifail, kdim, kdims
      EXTERNAL  zsub
      TYPE(zseg) :: mele(*)
      TYPE(zind) :: mind
C
      COMMON /UC160/ pi, deg, re, gmagmo, eclipt, geoid, uma
      REAL*8         pi, deg, re, gmagmo, eclipt, geoid(3), uma(30)   
*
*     mind    pointers into mele
*     mele    Set of the elementary field line segments
*     kpar1   control parammeter passed to zsub
*     rpar2   control parameter passed to zsub
*     zsub    external subroutine
*     result  result of the integration
*     rnorm   time integral of the unity
*     ifail   error flag
*                                              
c
c      Intgr. of fct. times ds/sqrt(1-B/Bm) unit: Re times (unit of fct)
c      and       one  times ds/sqrt(1-B/Bm) unit: Re       
c
C         
C     VARIABLES     
C            
      integer*4     jstart,jend,jtmp,ja,jb,nsize,nlog
      parameter(    nlog=6,nsize=65 )
c                          nsize=1+2^nlog
      TYPE(zpnt) :: mgeo(nsize)
      REAL*8        sqval(nsize),atval(300,nsize),sderv
      real*8        bm,rs0(300),rs1(300),rn0,rn1,sba,frac,eps
      integer*4     istep, i0, ip, ikor, i
      logical*1     flag
c
c
c
      ifail=0
      kdim=1
      goto 17
      entry ud323 (mind, mele, kpar1, rpar2, zsub, result, kdims,
     :             rnorm, ifail)
      kdim=kdims
   17 if(kdim.lt.0 .or.kdim.gt.300)then
        ifail=-32301
        return
      endif
c
      jstart=mind%jmirpn
      jend=mind%jmirps
      if( jend .lt. jstart )then
        jstart=mind%jmirps
        jend=mind%jmirpn
      endif
      if( jstart.lt.0 )then
        ifail=-32101
c       invalid field line description
        return                             
      endif
c
c     B at mirror point
      bm=max(mele(jstart)%beg%b%dnrm,mele(jend)%beg%b%dnrm)
ccc      print 1000,jstart,jend,bm,kpar1,rpar2
 1000 format('# UD321: Start, End, Bm, Kpar1, Rpar2=',2I5,F8.5,i4,1P,
     :       E13.4)
 1010 format('# UD321 (',a1,'): Eval., Alt., Lat., Lon., Value, Store',
     :        /i10,f10.2,2f8.2,1p,E13.4,i4)
 1030 format('# UD321: from',i5,F8.5,F10.3,' to',i5,F8.5,F10.3,f8.1)
 1020 format('# UD321: istep, norm, integr=',i4,1P,3e13.5)      
c
c     nbr of iterations, integration of zsub, integration of 1
      ikor=0
      do i=1,kdim
        result(i)=0.0d0
      enddo
      rnorm=0.0d0
c
c     for future use
      mele(jstart)%dtbnd=0.0d0
c
c
c     *** the segment includes more than two points
      if ( jend-jstart .gt. 1 ) then
        jtmp=jstart+2
        if (jtmp.gt.jend) jtmp=-1
        ja=jstart
c 
c       eval sqrt(1-B/Bm) and function zsub at ja=jstart
        sqval(1)=0.0D0
        ikor=ikor+1
        call zsub(kpar1, rpar2, mele(ja)%beg%coord, mele(ja)%beg%b, 
     :            atval(1,1), ifail)
        if( ifail.lt.0)return
ccc        print1010,'a',ikor,mele(ja)%beg%coord%radius-re,90.0-
ccc     :      mele(ja)%beg%coord%colat,mele(ja)%beg%coord%elong,
ccc     :      atval(1),1
        do jb=jstart+1,jend
c
c         first element --> ja=jstart, jb=jstart+1, 
c                           jtmp=jstart+2 or -1 if only one element
c         after --> jtmp=jb-2, ja=jb-1, jb=jb
c
c         sqval() = sqrt( 1 - B/Bm )
c         sqval(1) refers to ja; sqval(nsize) refers to jb
c         when ja=jstart, sqval(1)=0
c         when jb=jend, sqval(nsize)=0
c
c         length of the segment
          sba=mele(jb)%arcl-mele(ja)%arcl
          if( sba.eq.0 )then
            ifail=-32102
c           wrong elementary segment
            return
          endif
ccc          print1030,ja,mele(ja)%beg%b%dnrm,mele(ja)%arcl,
ccc     :            jb,mele(jb)%beg%b%dnrm,mele(jb)%arcl,
ccc     :            mele(jb)%arcl-mele(ja)%arcl
cc          print*,' from pts',ja,' to',jb,sba,' km'
cc          print*,' Bm, Ba, Bb=',bm,mele(ja)%beg%b%dnrm,
cc     :                             mele(jb)%beg%b%dnrm
cc          print*,' Sa, Sb=',bm,mele(ja)%arcl, mele(jb)%arcl
c
c         evaluate sqval(nsize) and atval(nsize):
          if( jb .ne. jend )then
            sqval(nsize)=sqrt(1-mele(jb)%beg%b%dnrm/bm)
          else
            sqval(nsize)=0.0D0
          endif
          ikor=ikor+1
          call zsub(kpar1, rpar2, mele(jb)%beg%coord, mele(jb)%beg%b, 
     :              atval(1,nsize), ifail)
          if( ifail.lt.0)return
cc          print*,' ud321 b ',ikor,atval(nsize),nsize      
ccc          print1010,'b',ikor,mele(jb)%beg%coord%radius-re,90.0-
ccc     :      mele(jb)%beg%coord%colat,mele(jb)%beg%coord%elong,
ccc     :      atval(nsize),nsize
c
c         first evaluation of the integral (rn = integr of 0.5)
          rn0=1.0d0/(sqval(1)+sqval(nsize))
          do i=1,kdim
            rs0(i)=rn0*(atval(i,1)+atval(i,nsize))
          enddo
cc          print*,' ud321->',rn0,rs0
ccc          print1020,0,rn0*2,rs0
c
c         compute 15 intermediates points (nsize=17)
          call ud329( mele,ja,jb,jtmp,mgeo,nsize,ifail)
          if( ifail.lt.0)return
c
c         rn0 and rs0 are evaluated between pts 1 and 17
c         rn1 and rs1 will be evaluated in two steps 1-9 
c         and 9-17 (istep=8,i.e. (nsize-1)/2)
c         if needed then in four steps, eight or sixteen steps
c
          istep=(nsize-1)/2
   10     continue
            flag=.true.
            frac=istep
            frac=frac/(nsize-1)
            rn1=0.0d0
            do i=1,kdim
              rs1(i)=0.0d0
            enddo
            ip=1
            do i0=1+istep,nsize,istep
c
c             first time when istep=8, [ip,i0] are equal to [1,9] 
c                                                           and [9,17]
c             second time when istep=4, they are equal to [1,5], [5,9],
c                                                   [9,13] and [13,17]
c             third and fourth time, idem
c
c             note that only half the point have to be evaluated: 
c               9 for [1,9], nothing for [9,17], 5 for [1,5], nothing 
c               for [5,9], 13 for [9,13], etc
C             --> flag = true, false, true, false, ...
              if(flag)then
                ikor=ikor+1
                call zsub(kpar1,rpar2,mgeo(i0)%coord,
     :                    mgeo(i0)%b,atval(1,i0),ifail)
                if( ifail.lt.0)return
c                print*,' ud321 c ',ikor,atval(i0),i0      
cc                 print1010,'c',ikor,mele(i0)%beg%coord%radius-re,90.0-
cc     :              mele(i0)%beg%coord%colat,mele(i0)%beg%coord%elong,
cc     :              atval(i0),i0
                sqval(i0)=sqrt(1.0D0-mgeo(i0)%b%dnrm/bm)
c                print*,' ajout: #, Bi:',i0,mgeo(i0)%b%dnrm
              endif
c
c             evaluation of the integral with 2, 4, 8 or 16 steps
              rn1=rn1+frac/(sqval(i0)+sqval(ip))
              do i=1,kdim
                rs1(i)=rs1(i)+frac*(atval(i,i0)+atval(i,ip))
     :                    /(sqval(i0)+sqval(ip))
              enddo
c
c             next point
              ip=i0
              flag=.not.flag
            enddo
c
c           relative error in per thousand between current and 
c           previous evaluation
            eps=0
            do i=1,kdim
              eps=max(eps,
     :                2000.0d0*(rs1(i)-rs0(i)) /
     :                max(1.d-12,abs(rs1(i))+abs(rs0(i))) )
c            print*,' ud321=>',rn1,rs1,eps
ccc            print1020,istep,rn1*2,rs1,eps
c
              rs0(i)=rs1(i)
            enddo
            rn0=rn1
c
c           repeat until eps is greater than 1
            istep=istep/2
          if(istep.gt.0 .and. abs(eps).ge.5)goto 10
c
          rnorm=rnorm+2*rn1*sba
          do i=1,kdim
            result(i)=result(i)+rs1(i)*sba
          enddo
c
cc          stop
c
c         for future use 
          mele(ja)%dtbnd=mele(ja)%dtbnd+sba*rn1
          mele(jb)%dtbnd=sba*rn1
c
c         End of loop
c           --> step to the next segment
          jtmp=ja
          ja=jb
          sqval(1)=sqval(nsize)
          do i=1,kdim
            atval(i,1)=atval(i,nsize)
          enddo
        enddo
      else
c     ***  here jend=jstart (a single point) or jend=jstart+1 (a single
c                                                               element)
        do ja=jstart,jend
          ikor=ikor+1
          call zsub(kpar1, rpar2, mele(ja)%beg%coord,mele(ja)%beg%b, 
     :              atval(1,1),ifail)
c          print*,' ud321 d ',ikor,atval(1),ja      
ccc          print1010,'d',ikor,mele(ja)%beg%coord%radius-re,90.0-
ccc     :      mele(ja)%beg%coord%colat,mele(ja)%beg%coord%elong,
ccc     :      atval(1),1
          if(ifail.lt.0)return
c
ccc          print*,'# UD320 single point or single element'
ccc          print*,'#       jstart, jend=',jstart,jend
ccc          print*,'#       calling UD328' 
          call ud328(mele(ja)%beg%coord,sderv,ifail)
          sderv=-sderv
          if(ifail.lt.0)return
          if(sderv.le.0)then
            ifail=-32103
c           problem with an equatorial single point
            if (.false.) print*,'# UD321: UD328 ->',sderv,ifail
            return
          endif
c
          mele(ja)%dtbnd=pi*sqrt(2.0D0/sderv)/(jend-jstart+1)
          rnorm = rnorm + mele(ja)%dtbnd
          do i=1,kdim
            result(i) = result(i) + atval(i,1) * mele(ja)%dtbnd
          enddo
        enddo
      endif
c
      ifail=ikor
c
      end
C----------------------------------------------------------------------
      SUBROUTINE UD324
     :         (mind, mele, avdens, rnorm, ifail)   
C
C!    evaluate average number densities over a magnetic field line
C
      include 'structure.h'
c
c     INTERFACE
      REAL*8    avdens(30), rnorm
      INTEGER*4 ifail
      TYPE(zseg) :: mele(*)
      TYPE(zind) :: mind
C
c
c     = copy of ud321 but call to ua630 instead of zsub
c
*
*     mind    pointers into mele
*     mele    Set of the elementary field line segments
*     result  result of the integration
*     rnorm   time integral of the unity
*     ifail   error flag
*                                              
c
c      Intgr. of fct. times ds/sqrt(1-B/Bm) unit: Re times (unit of fct)
c      and       one  times ds/sqrt(1-B/Bm) unit: Re       
c
C         
C     VARIABLES     
C            
      integer*4     jstart,jend,jtmp,ja,jb,nsize,nlog
      parameter(    nlog=4,nsize=17 )
c                          nsize=1+2^nlog
      TYPE(zpnt) :: mgeo(nsize)
      REAL*8        sqval(nsize),atval(30,nsize),sderv
      real*8        bm,rs0(30),rs1(30),rn0,rn1,sba,frac,eps
      integer*4     istep, i0, ip, ikor,is
      logical*1     flag
      real*8        tmass,temp(3)
c
      COMMON /UC150/ matm, ntspec, nnspec, kspec, kswitch
      TYPE(zatm) :: matm
      INTEGER*4     ntspec, nnspec, kspec(30)
      INTEGER*4     kswitch(50)      
c
c
c
      ifail=0
      jstart=mind%jmirpn
      jend=mind%jmirps
      if( jend .lt. jstart )then
        jstart=mind%jmirps
        jend=mind%jmirpn
      endif
      if( jstart.lt.0 )then
        ifail=-32104
c       invalid field line description
        return                             
      endif
c
      bm=(mele(jstart)%beg%b%dnrm+mele(jend)%beg%b%dnrm)/2
      ikor=0
      do is=1,30
        avdens(is)=0.0d0
      enddo
      rnorm=0.0d0
c
c      print*,'ud324',jstart,jend,bm
c       
      mele(jstart)%dtbnd=0.0d0
      if (jend.gt.jstart) then
        jtmp=jstart+2
        if (jtmp.gt.jend) jtmp=-1
        ja=jstart
        do jb=jstart+1,jend
          sba=mele(jb)%arcl-mele(ja)%arcl
          if( sba.eq.0 )then
            ifail=-32104
c           wrong elementary segment
            return
          endif
c
c          print*,'ud329'
          call ud329( mele,ja,jb,jtmp,mgeo,nsize,ifail)
          if( ifail.lt.0)return
c
          sqval(1)=sqrt(max(0.0d0,1.0D0-mgeo(1)%b%dnrm/bm))
c          print*,'ua630a'
          call ua630(mgeo(1)%coord,tmass,atval(1,1),temp, ifail)
          if( ifail.lt.0)return
c
          sqval(nsize)=sqrt(max(0.0d0,1.0D0-mgeo(nsize)%b%dnrm/bm))
c          print*,'ua630b'
          call ua630(mgeo(nsize)%coord,tmass,atval(1,nsize),temp, ifail)
          if( ifail.lt.0)return
c
          rn0=1.0d0/(sqval(1)+sqval(nsize))
          do is=1,30
            rs0(is)=rn0*(atval(is,1)+atval(is,nsize))
          enddo
c
          istep=nsize-1
   10     continue
            flag=.true.
            frac=istep
            frac=frac/(nsize-1)
            rn1=0.0d0
            do is=1,30
              rs1(is)=0.0d0
            enddo
            ip=1
            do i0=ip+istep,nsize,istep
              if(flag)then
c                print*,'ua630c'
                call ua630(mgeo(i0)%coord,tmass,atval(1,i0),temp, ifail)
                if( ifail.lt.0)return
                sqval(i0)=sqrt(max(0.0d0,1.0D0-mgeo(i0)%b%dnrm/bm))
              endif
              rn1=rn1+frac/(sqval(i0)+sqval(ip))
              do is=1,30
                rs1(is)=rs1(is)+frac*(atval(is,i0)+atval(is,ip))
     :                         /(sqval(i0)+sqval(ip))
              enddo
              ip=i0
              flag=.not.flag
            enddo
            eps=0.0d0
            do is=1,ntspec
              if( rs1(kspec(is)) .ne.  rs0(kspec(is)) )
     :           eps=eps+(2*(rs1(kspec(is))-rs0(kspec(is))) /
     :                   (rs1(kspec(is))+rs0(kspec(is))))**2
              rs0(kspec(is))=rs1(kspec(is))
            enddo
            eps=eps/ntspec*100.0d0
            rn0=rn1
            istep=istep/2
          if(istep.gt.0 .and. eps.ge.1)goto 10
          if(eps.ge.1)ikor=ikor+nsize
          rnorm=rnorm+rn1*sba
          do is=1,30
            avdens(is)=avdens(is)+rs1(is)*sba
          enddo
          mele(ja)%dtbnd=mele(ja)%dtbnd+sba*rn1/2
          mele(jb)%dtbnd=sba*rn1/2
c
          jtmp=ja
          ja=jb
        enddo
      else
c       here jend=jstart:  a single point
c        print*,'ua630d'
        call ua630(mele(jstart)%beg%coord,tmass,atval(1,1),temp, ifail)
        if(ifail.lt.0)return
c        print*,'ud328'
        call ud328(mele(jstart)%beg%coord ,sderv,ifail)
        if(ifail.lt.0)return
        if(sderv.le.0)then
          ifail=-32403
c         problem with an equatorial single point
        endif
        rnorm=sqrt(8.0d0/sderv)
        do is=1,30
          avdens(is)=atval(is,1)*rnorm
        enddo
        mele(jstart)%dtbnd=rnorm
      endif
c
      ifail=ikor
c
      end
C----------------------------------------------------------------------
      SUBROUTINE UD325
     :         (avdens, rnorm, rdrift, ifail)   
C
C!    evaluate average atmospheric densities over a drift shell
C
      include 'structure.h'
c
c     INTERFACE
      REAL*8    avdens(30), rnorm, rdrift
      INTEGER*4 ifail
C
c
c     Mainly based on ud320 but call ua630 instead of zsub
c     Note that the subroutine returns in argument AVDENS
c     the average number densities, i.e. the equivalent of
c     RESULT/RNORM of subroutine ud320 
c
*
*     avdens  average number densities (see note in subroutine ua630)
*     rnorm   time integral of the unity
*     rdrift  single time integral of the unity
*     ifail   error flag
*                                              
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
      COMMON /UC130/ nbrsg, kursg, mseg
C
        INTEGER*4     nbrsg, kursg
        TYPE(zseg) :: mseg(nx130)
C
      COMMON /UC160
     :               /pi, deg, re, gmagmo, eclipt, geoid, uma
C
        REAL*8        pi, deg, re, gmagmo, eclipt, geoid(3), uma(30)
C
C         
C     VARIABLES     
C            
      REAL*8         rone,rline(30),bm,rn
      TYPE(zxyz) ::  mvlcr, mpsnx, mpscr, mdspr,mdscr
      integer*4      kprev,knext,is
      REAL*8         st,sp,ct,cp,fct,dnrm
C
C     CODE
C
      ifail=0
      rdrift=0
      rnorm=0.0d0
      do is=1,30
        avdens(is)=0.0d0  
      enddo
C
      if( k1st.le.0)then
        ifail=-32501
c       invalid drift shell description
        return                             
      endif
      kurfl=k1st
C
      bm = mlbl%fbmp
      if( bm.le.0 .or. .not.mlbl%lbmp)then
        ifail=-32502
c       invalid drift shell description
        return
      endif
C
      kprev=mfl(kurfl)%kwest
c
c     mpscr = GEO cartesian position in km at kurfl      
c     mpsnx = GEO cartesian position in km at knext
c     mdspr = cartesian distance between kprev and kurfl
c     mdscr = cartesian distance between kurfl and knext      
c     mvlcr = cartesian comp. of unit vector // to drift vel 
c
c     mfl%drift in km-1 Gauss-1 at kurfl
c
c    conversion to velocity: mass c^2/2/q beta^2 gamma (2-b/bm)
c    beta = v/c 
c    gamma = sqrt(1-beta^2)
c
c    note: here only (2-b/bm)
c
c
      st  = sin( mfl(kprev)%equat%coord%colat*deg )
      ct  = cos( mfl(kprev)%equat%coord%colat*deg )
      sp  = sin( mfl(kprev)%equat%coord%elong*deg )
      cp  = cos( mfl(kprev)%equat%coord%elong*deg )
      mdspr%z = ct    * mfl(kprev)%equat%coord%radius
      mdspr%x = st*cp * mfl(kprev)%equat%coord%radius
      mdspr%y = st*sp * mfl(kprev)%equat%coord%radius
c
      st  = sin( mfl(kurfl)%equat%coord%colat*deg )
      ct  = cos( mfl(kurfl)%equat%coord%colat*deg )
      sp  = sin( mfl(kurfl)%equat%coord%elong*deg )
      cp  = cos( mfl(kurfl)%equat%coord%elong*deg )
      mpscr%z = ct    * mfl(kurfl)%equat%coord%radius
      mpscr%x = st*cp * mfl(kurfl)%equat%coord%radius
      mpscr%y = st*sp * mfl(kurfl)%equat%coord%radius
c
      mdspr%z = mpscr%z - mdspr%z
      mdspr%x = mpscr%x - mdspr%x
      mdspr%y = mpscr%y - mdspr%y
c
  100 continue
c
        fct = 2.0D0-mfl(kurfl)%equat%b%dnrm/bm
        mvlcr%z = mfl(kurfl)%drift%rho*ct - mfl(kurfl)%drift%theta*st
        mvlcr%z = mvlcr%z/mfl(kurfl)%drift%dnrm
        mvlcr%x =   mfl(kurfl)%drift%rho * st*cp
     :            + mfl(kurfl)%drift%theta * ct*cp
     :            - mfl(kurfl)%drift%phi * sp
        mvlcr%x = mvlcr%x/mfl(kurfl)%drift%dnrm
        mvlcr%y =   mfl(kurfl)%drift%rho * st*sp
     :            + mfl(kurfl)%drift%theta * ct*sp
     :            + mfl(kurfl)%drift%phi * cp
        mvlcr%y = mvlcr%y/mfl(kurfl)%drift%dnrm
        dnrm    = mfl(kurfl)%drift%dnrm * fct
c
        knext=mfl(kurfl)%keast
c
        st  = sin( mfl(knext)%equat%coord%colat*deg )
        ct  = cos( mfl(knext)%equat%coord%colat*deg )
        sp  = sin( mfl(knext)%equat%coord%elong*deg )
        cp  = cos( mfl(knext)%equat%coord%elong*deg )
        mpsnx%z = ct    * mfl(knext)%equat%coord%radius
        mpsnx%x = st*cp * mfl(knext)%equat%coord%radius
        mpsnx%y = st*sp * mfl(knext)%equat%coord%radius
c
        mdscr%z = mpsnx%z - mpscr%z
        mdscr%x = mpsnx%x - mpscr%x
        mdscr%y = mpsnx%y - mpscr%y
c
        rn=(   (mdspr%z+mdscr%z)*mvlcr%z
     :       + (mdspr%x+mdscr%x)*mvlcr%x
     :       + (mdspr%y+mdscr%y)*mvlcr%y ) / 2.0D0 / dnrm
c
c       unit: Gauss
c
        call ud324(mfl(kurfl)%ind,mseg, 
     :             rline, rone, ifail)
        if(ifail.lt.0)return
c
c  mfl.drift (/zvec/) set by uf424 (arg mxnb)
C  uf424 called by uf420 called by uf411 called by uf410
c        called by ud310: "trace a magnetic drift shell" 
c  mfl%drift = mxnb 
c            = Vector pointing to the binormal vector of the 
c              magnetic field line at EQUAT
c
c  see JS Lew (JGR 66,1961, p2681-2685) eq (3)
c
c  mxnb = e x n /rcurv /b
c
c  Drift velocity at EQUAT is prop. to mxnb (2-B/Bm)
c
c
c
        rdrift = rdrift+rn
        rnorm  = rnorm+rone*rn
        do is=1,30
          avdens(is) = avdens(is)+rline(is)*rn
        enddo
c
c
        mpscr=mpsnx
        mdspr=mdscr
        kprev=kurfl
        kurfl=knext
      if(kurfl.ne.k1st)goto 100
C
      if( rnorm.le.0.0d0)then
        ifail=-32503
      else
        do is=1,30
          avdens(is) = avdens(is)/rnorm
        enddo
      endif      
C
      end       
C----------------------------------------------------------------------
      SUBROUTINE UD327
     :              ( mgp1, mgp2, mgp3, fct1, fct2, fct3, mextr, 
     :                est, ifail )
C
      include 'structure.h'
c
      TYPE(zgeo) :: mgp1, mgp2, mgp3, mextr
      REAL*8        fct1, fct2, fct3
      INTEGER*4     ifail
c      
      COMMON /UC160/ pi, deg, re, gmagmo, eclipt, geoid, uma
c
      REAL*8         pi, deg, re
      REAL*8         gmagmo
      REAL*8         eclipt, geoid(3), uma(30)    
c
c
      real*8   st,a,b,d12,d23,dx,est
      TYPE(zxyz) :: mx1,mx2,mx3,mres
c
c
      ifail=0
      mextr=mgp2
c
      st     = sin( mgp1%colat*deg ) * mgp1%radius
      mx1%z  = cos( mgp1%colat*deg ) * mgp1%radius
      mx1%x  = st * cos( mgp1%elong*deg )
      mx1%y  = st * sin( mgp1%elong*deg )
      st     = sin( mgp2%colat*deg ) * mgp2%radius
      mx2%z  = cos( mgp2%colat*deg ) * mgp2%radius
      mx2%x  = st * cos( mgp2%elong*deg )
      mx2%y  = st * sin( mgp2%elong*deg )
      st     = sin( mgp3%colat*deg ) * mgp3%radius
      mx3%z  = cos( mgp3%colat*deg ) * mgp3%radius
      mx3%x  = st * cos( mgp3%elong*deg )
      mx3%y  = st * sin( mgp3%elong*deg )
c
      d12    = sqrt((mx2%z-mx1%z)**2+(mx2%x-mx1%x)**2+(mx2%y-mx1%y)**2)
      d23    = sqrt((mx2%z-mx3%z)**2+(mx2%x-mx3%x)**2+(mx2%y-mx3%y)**2)
      if( d12*d23.eq.0 )then
        ifail=-32701
        return
      endif 
c
      a=(fct1-fct2)*d23**2-(fct3-fct2)*d12**2
      b=(fct2-fct3)*d12+(fct2-fct1)*d23
      if( b.eq.0 )then
        ifail=-32701
        return
      endif 
      dx     = -a/b/2
      est    = fct2+a**2/(b*4*d12*d23*(d12+d23))
c
      if( dx .lt. 0 )then
        mres%x = (-mx1%x*dx+mx2%x*(dx+d12))/d12
        mres%y = (-mx1%y*dx+mx2%y*(dx+d12))/d12
        mres%z = (-mx1%z*dx+mx2%z*(dx+d12))/d12
      else
        mres%x = (mx2%x*(d23-dx)+mx3%x*dx)/d23
        mres%y = (mx2%y*(d23-dx)+mx3%y*dx)/d23
        mres%z = (mx2%z*(d23-dx)+mx3%z*dx)/d23
      endif
c
      mextr%radius  = sqrt(mres%x**2+mres%y**2+mres%z**2)
      mextr%colat   = acos(mres%z/mextr%radius)/deg
      if( (mres%x).eq.0 .and. (mres%y).eq.0)then
        mextr%elong = 0
      else
        mextr%elong = atan2(mres%y,mres%x)/deg
      endif          
 1000 format(1x,a1,4f15.5)
c
      end
C----------------------------------------------------------------------
      SUBROUTINE UD328
     :                ( mgp, sderv, ifail )
C
C     Evaluate the quantity: - [ d2B / ds2 ] / B
C
      INCLUDE 'structure.h'
C
C     INTERFACE
C
        TYPE(zgeo) :: mgp
        REAL*8        sderv
        INTEGER*4     ifail           
C
      COMMON /UC160
     :               /pi, deg, re, gmagmo, eclipt, geoid, uma
C              
        REAL*8        pi, deg, re, gmagmo, eclipt, geoid(3), uma(30)
C                              
      COMMON /UC170/  nsg, kgp, mlab, mlin, mele
C
        INTEGER*4     nsg, kgp
        TYPE(zlbl) :: mlab
        TYPE(zfln) :: mlin
        TYPE(zseg) :: mele(nx170)
C
C     VARIABLES
C          
      integer*4       k0, kold, kretard
      real*8          rc, scale, pb1, palt, pb2, dst1, dst2, pb0
C
C     CODE          
C
      sderv                     = 0.0D0
      kretard                   = 0
      ifail                     = 0
      kold                      = kgp
C
      if( MAX(nsg,kgp) .gt. nx170-10 )then
        ifail                   = -32801
c          idem que -99801
        return
      endif
C
      k0                        = nx170-10
      mele(k0)%beg%coord        = mgp
      mele(k0)%arcl             = 0.0d0
      CALL UM530 (mgp, mele(k0)%beg%b, ifail)
      if( ifail .le. -53002 .and. ifail .ge. -53004)then
        kretard=ifail
        ifail=0
      endif
      if( ifail .lt. 0 ) return
      pb0                       = mele(k0)%beg%b%dnrm
      CALL UT999 (mgp, rc)
c
      scale                     = 0.2d0
      kgp                       = k0
      CALL UF423(scale, rc, pb1, palt, ifail)
      if( ifail .le. -42304 .and. ifail .ge. -42306)then
        kretard=ifail
        ifail=0
      endif
      if(ifail .lt. 0) return
      dst1                      = mele(kgp)%arcl - mele(k0)%arcl
c
      scale                     = -0.2d0
      kgp                       = k0
      CALL UF423(scale, rc, pb2, palt, ifail)
      if(ifail .lt. 0) return
      if( ifail .le. -42304 .and. ifail .ge. -42306)then
        kretard=ifail
        ifail=0
      endif
      dst2                      = mele(kgp)%arcl - mele(k0)%arcl
c
      kgp                       = kold
c
      sderv= -2*((pb2/pb0-1)/dst2-(pb1/pb0-1)/dst1)/(dst2-dst1)
      ifail=abs(kretard)
c
      end
C----------------------------------------------------------------------
      SUBROUTINE UD329
     :                ( mele, ja, jb, jtmp, mgeo, npt,ifail)
C
      include 'structure.h'
c
c!    Interpolate between two or three points
c
      TYPE(zseg) :: mele(*)
      INTEGER*4     ja, jb, jtmp, ifail, npt
      TYPE(zpnt) :: mgeo(npt)
c
*
*     mele     field line segment
*     ja       index of the first point
*     jb       index of the last point
*     jtmp     index of an other point
*     mgeo     interpolated points (include the points indiced by
*                                    ja and jb)
*     npt      number of points
*     ifail    error flag
C
C     The subroutine UD329 interpolates npt-2 points along
C     a magnetic field line segment.
C     The indices of the points are passed by the arguments
C     ja, jb and jtmp. If the argument jtmp is set to a
C     value less than zero, the subroutine makes a linear 
C     interpolation between the two points mele(ja) and
C     mele(jb). Otherwize, a parabolic interpolation is
C     applied.  The magnetic field vector is evaluated
C     at the interpolated points 
C
c
      TYPE(zgeo) :: mtmp
      TYPE(zxyz) :: ma, mba, mta, m0, mcoefa, mcoefb
      real*8        sba, sta, sbt, d, s0a, r
      integer*4     i, ifl
c
      if( npt .lt. 2 )then
        ifail    = -32901
c       wrong value of the argument npt
        return
      endif
      mgeo(1)    = mele(ja)%beg
      mgeo(npt)  = mele(jb)%beg
      sba        = mele(jb)%arcl - mele(ja)%arcl
c
      call ut541(mgeo(1)%coord,ma) 
      call ut541(mgeo(npt)%coord,mba) 
      mba%x      = mba%x - ma%x
      mba%y      = mba%y - ma%y
      mba%z      = mba%z - ma%z
      r          = sba / (npt-1)
      ifl        = 0
      if( jtmp .gt. 0 )then
        mtmp     = mele(jtmp)%beg%coord
        call ut541(mtmp,mta) 
        mta%x    = mta%x - ma%x
        mta%y    = mta%y - ma%y
        mta%z    = mta%z - ma%z
        sta      = mele(jtmp)%arcl - mele(ja)%arcl
        sbt      = sba - sta
        if( sba*sbt*sta .eq. 0)  then
          ifail  = -32902
c         not separated points
          return
        endif
        d        = 1 / sbt
        mcoefa%x = d * (mba%x/sba-mta%x/sta)
        mcoefa%y = d * (mba%y/sba-mta%y/sta)
        mcoefa%z = d * (mba%z/sba-mta%z/sta)
        mcoefb%x = d * (sba*mta%x/sta-sta*mba%x/sba)
        mcoefb%y = d * (sba*mta%y/sta-sta*mba%y/sba)
        mcoefb%z = d * (sba*mta%z/sta-sta*mba%z/sba)
        do i=2, npt-1
          s0a    = r * (i-1)
          m0%x   = ma%x + s0a * ( mcoefb%x + s0a * mcoefa%x )       
          m0%y   = ma%y + s0a * ( mcoefb%y + s0a * mcoefa%y )       
          m0%z   = ma%z + s0a * ( mcoefb%z + s0a * mcoefa%z )       
          call ut546(m0,mgeo(i)%coord)
          ifl    = 0
          call um530(mgeo(i)%coord,mgeo(i)%b,ifl)
          if( ifl.lt.0 )ifail=ifl
          call ut999(mgeo(i)%coord,mgeo(i)%rcurv)
        enddo
      else
        if( sba .eq. 0)  then
          ifail  = -32903
c         not separated points (same as -32903)
          return
        endif
        mcoefa%x = mba%x/sba
        mcoefa%y = mba%y/sba
        mcoefa%z = mba%z/sba
        do i=2, npt-1
          s0a    = r * (i-1)
          m0%x   = ma%x + s0a * mcoefa%x        
          m0%y   = ma%y + s0a * mcoefa%y        
          m0%z   = ma%z + s0a * mcoefa%z        
          call ut546(m0,mgeo(i)%coord)
          ifl    = 0
          call um530(mgeo(i)%coord,mgeo(i)%b,ifl)
          if( ifl.lt.0 )ifail=ifl
          call ut999(mgeo(i)%coord,mgeo(i)%rcurv)
        enddo
      endif
c
      end      
C------------------------------------------------------------------