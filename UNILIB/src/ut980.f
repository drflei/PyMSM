# 1 "ut980.f"
      SUBROUTINE UT980(lcall, kfail, kunit, ifail)
C
C     print the error messages
C
      implicit none
cDEC$ IF DEFINED (_x86_)
cDEC$ ATTRIBUTES DLLEXPORT :: UT980
cDEC$ ENDIF
C
C     INTERFACE
C
        INTEGER*4     ifail, kunit, kfail
        CHARACTER*(*) lcall
C  
C     VARIABLES
C      
        INTEGER*4     icall, isub, i, ierr
        CHARACTER*72  null
        CHARACTER*2   lXX(17),lyy(17)
        INTEGER*4     kxx(17)
        character*33  lmm(17)
C
        data kxx/  22,  23,  24,  31,  32,  33,  41,  42,
     :             51,  52,  53,  54,  55,  61,  63,  98,  99/
        data lxx/'UL','UL','UL','UD','UD','UD','UF','UF',
     :           'UM','UM','UM','UT','UT','UA','UA','UT','UT'/
        data lyy/'ul','ul','ul','ud','ud','ud','uf','uf',
     :           'um','um','um','ut','ut','ua','ua','ut','ut'/
        data lmm/'the magnetic coordinate treatment',
     :           'the integral invariant evaluation',
     :           'the Hilton''s transformation      ',
     :           'the magnetic drift shell tracing ',
     :           'the averaging over drift shell   ',
     :           'the third invariant computation  ',
     :           'the magn. field line segm. search',
     :           'the magnetic field line tracing  ',
     :           'the geomagnetic field model setup',
     :           'the extrn. mag. field model setup',
     :           'the magnetic field evaluation    ',
     :           'a coordinate or time conversion  ',
     :           'a coordinate or vector transform.',
     :           'the atmospheric model setup      ',
     :           'the atmospheric evaluation       ',
     :           'accessing a UNILIB tool          ',
     :           'accessing a UNILIB tool          '/
c
c     CODE
C
      ifail=0
      if( kunit .le. 0 )return
      write(kunit,1100,err=2000)
c
c
      null = ' '
      if( len(lcall) .lt. 5 )then
        ifail = -98001
C       Invalid value of the argument LCALL
        write(kunit,1140,err=2000)1,lcall,kfail
        write(kunit,1120,err=2000)
        return
      endif
      read(lcall,1000,err=20)icall
      i = icall/10
      do isub=1,17
        if( i .eq. kxx(isub) )goto 100
      enddo
   20 continue
      ifail = -98006
C       Incorrect name of subroutine 
        write(kunit,1140,err=2000)6,lcall,kfail
        write(kunit,1120,err=2000)
      return
  100 continue
      if( lcall(1:2).ne.lxx(isub) .and.
     :    lcall(1:2).ne.lyy(isub) )then
        ifail = -98005
C       Incorrect name of subroutine
        write(kunit,1140,err=2000)5,lcall,kfail
        write(kunit,1120,err=2000)
        return
      endif
      i=(-kfail)/1000
      do ierr=1,17
        if( i .eq. kxx(ierr) )goto 110
      enddo
        ifail = -98002
c       Invalid error code
        i=icall
        write(kunit,1130,err=2000) lxx(isub),icall,kfail,lyy(isub),i
        write(kunit,1120,err=2000)
        return
  110 continue
c
c
c
      i=(-kfail)/100
      write(kunit,1110,err=2000) lxx(isub),icall,lmm(ierr),
     :    lxx(ierr),i,kfail,lyy(ierr),i
      write(kunit,1120,err=2000)
      return
c
 2000 if( ifail .ge. 0 )ifail = -98004
c    
 1000 format(2x,i3)
 1100 format(1x,32('-'),' UNILIB error ',32('-'))
 1110 format(3x,'The subroutine ',a2,i3.3,' of the UNILIB library has',
     :       ' been called and an error has',
     :      /3x,' occured during ',a33,/3x,' in subroutine ',a2,i3.3,
     :       ' with the error code ',i6.5,'.',
     :      /3x,'Please, consult the documentation')
 1120 format(1x,78('-'))
 1130 format(3x,'The subroutine ',a2,i3.3,' of the UNILIB library has',
     :       ' been called and an error has',
     :      /3x,' occured with the error code ',i10,'. The source of ',
     :       'the error is not traceable.'  
     :      /3x,'Please, consult the documentation')
 1140 format(3x,'An error not referenced has occured in a subroutine',
     :       ' of the UNILIB library.',/50x,'<',i1,': ',a,1x,i6,'>',
     :      /3x,'Please, consult the documentation')
C
      END
C----------------------------------------------------------------------
      SUBROUTINE UT981(k1, k2, ifail)
        INTEGER*4 k1, k2, ifail
        ifail=0
      END
C----------------------------------------------------------------------
      SUBROUTINE UT982(k1, k2, k3, ifail)
        INTEGER*4 k1, k2, k3, ifail
        ifail=0
      END
C----------------------------------------------------------------------
      SUBROUTINE UT985 ( kindex, klength, mdata )  
C
C     Transfer a field line segment 
C
      INCLUDE 'structure.h'
C
C     INTERFACE
C
        INTEGER*4     kindex, klength
        TYPE(zseg) :: mdata(*)
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
      COMMON /UC170/  nsg, kgp, mlab, mlin, mele
C
        INTEGER*4     nsg, kgp
        TYPE(zlbl) :: mlab
        TYPE(zfln) :: mlin
        TYPE(zseg) :: mele(nx170)
C
C     VARIABLES
C
        INTEGER*4     j0, j1, jn, js, j
C
C     CODE
C
      if (kindex .eq. 0) then
        j0               = mlin%ind%jbeg
        j1               = mlin%ind%jend
        jn               = mlin%ind%jmirpn
        js               = mlin%ind%jmirps 
        if ( j0.gt.j1 .or. jn.lt.j0 .or. js.lt.j0 .or.
     :       jn.gt.j1 .or. js.gt.j1 .or. j0.lt.1 .or.
     :       j1.gt.nx170  .or. jn.gt.js ) then
          klength        = -999999999
          print*,'First pnt:',j0
          print*,'Last pnt: ',j1
          print*,'North MP: ',jn
          print*,'South MP: ',js                                          
          return
        endif
        if (klength .lt. j1-j0+1) then
          klength        = -(j1-j0+1) 
          return
        endif
        klength          = 0
        do j=j0,j1
          klength        = klength+1
          mdata(klength) = mele(j)
        enddo 
      else
        if(kindex .le. 0) kindex = k1st
        if(kindex.le.0 .or. kindex.gt.nbrfl)then
          klength        = -999999999
          print*,'k1st: ',k1st
          print*,'nbrfl:',nbrfl
          return
        endif
        j0               = mfl(kindex)%ind%jbeg
        j1               = mfl(kindex)%ind%jend
        jn               = mfl(kindex)%ind%jmirpn
        js               = mfl(kindex)%ind%jmirps
        if ( j0.gt.j1 .or. jn.lt.j0 .or. js.lt.j0 .or.
     :       jn.gt.j1 .or. js.gt.j1 .or. j0.lt.1 .or.
     :       j1.gt.nx130 .or. jn.gt.js ) then
          klength        = -999999999
          print*,'First pnt:',j0
          print*,'Last pnt: ',j1
          print*,'North MP: ',jn
          print*,'South MP: ',js                                          
          return
        endif    
        if (klength .lt. j1-j0+1) then
          klength        = -(j1-j0+1) 
          return
        endif  
        klength          = 0
        do j=j0,j1
          klength        = klength+1
          mdata(klength) = mseg(j)
        enddo
        kindex           = mfl(kindex)%keast
      endif
C
C
      END
C----------------------------------------------------------------------
      SUBROUTINE UT986 ( gm, colat, elong, msunidl, reidl, bsun )
c
      include 'structure.h'
cDEC$ IF DEFINED (_x86_)
cDEC$ ATTRIBUTES DLLEXPORT :: UT986
cDEC$ ENDIF
c
c    pass general variable to IDL
c
      real*8 gm,colat,elong,reidl,bsun
      TYPE(zxyz) :: msunidl
c
      COMMON /UC160/ pi, deg, re, gmagmo, eclipt, geoid, uma
c
      REAL*8         pi, deg, re
      REAL*8         gmagmo
      REAL*8         eclipt, geoid(3), uma(30)    
c
      COMMON /UC140/ mint, mext, msun
c
      TYPE(zimf) :: mint
      TYPE(zsun) :: msun
      TYPE(zemf) :: mext
c
      COMMON /UC190/ prop, stepx, stpmin, umsq, upsq, uk2, uk3, 
     :               epskm, epsrel, stplst, xclat, kmflg, nxstp       

      REAL*8         prop, stepx, stpmin
      REAL*8         umsq, upsq, uk2, uk3
      REAL*8         epskm, epsrel, stplst, xclat
      INTEGER*4      kmflg, nxstp      
c
c
c
      colat = mint%colat
      elong = mint%elong
      msunidl=msun%dir
      reidl=re
      bsun = mint%gmmo*sqrt(1.0D0+3.0D0*sin(mext%tilt*deg)**2)
c
      if( mod(kmflg,2) .eq. 0 )then
        gm=gmagmo
      else
        gm=mint%gmmo          
      endif
c
      end
C----------------------------------------------------------------------
