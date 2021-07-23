# 1 "ut990.f"
      SUBROUTINE UT990 ( kunit, kinit, ifail )
C
C!    Initialize and/or print the different constants
C
      INCLUDE 'structure.h'
cDEC$ IF DEFINED (_x86_)
cDEC$ ATTRIBUTES DLLEXPORT :: UT990
cDEC$ ENDIF
C
C     INTERFACE
C
        INTEGER*4     kunit, kinit, ifail
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
      COMMON /UC150/  matm, ntspec, nnspec, kspec, kflag
C 
        TYPE(zatm) :: matm
        INTEGER*4     ntspec, nnspec, kspec(30)
        INTEGER*4     kflag(50)
C
      COMMON /UC160/  pi, deg, re, gmagmo, eclipt, geoid, uma
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
      COMMON /UC190
     :               /prop, stepx, stpmin, umsq, upsq, uk2,
     :                uk3, epskm, epsrel, stplst, xclat,
     :                kmflg, kum533      
C
        REAL*8        prop, stepx, stpmin, umsq, upsq, uk2,
     :                uk3, epskm, epsrel, stplst, xclat
        INTEGER*4     kmflg, kum533       
C
      COMMON /UC192
     :               /xrmin, xbmin, xtmin, xbmax, epslon, epsfl,
     :                fvet, pvet, epsomeg, dltlat
        REAL*8        xrmin, xbmin, xtmin, xbmax, epslon, epsfl
        REAL*8        fvet, pvet, epsomeg, dltlat
C                
C     VARIABLES
C
!        character*33  message
        character*33  message, lver
!        INTEGER*4     iver,i
        INTEGER*4     i
        REAL*8        am(30)
C
      DATA AM / 1.d0,  4.d0, 20.d0, 40.d0, 84.d0, 131.d0,  2.d0, 28.d0,
     :  32.d0, 28.d0, 44.d0, 17.d0, 16.d0,  0.d0,   0.d0, 16.d0, 14.d0,
     :   0.d0,  1.d0,  1.d0,  4.d0, 12.d0, 14.d0,  16.d0, 20.d0, 28.d0,
     :  30.d0, 32.d0,  0.d0,  0.0/      
C
C        
C     CODE 
C
!      iver = 222

      message  = ' Common block content not altered'
C      
      if( kinit.ge.0 .or. kinit.eq.-160 )then
C      
C  a.1/ UC160
C
        message  = 'Common block UC160 re-initialized'
        pi       =  ATAN(1.0d0) * 4.0d0
        deg      =  pi / 180.0d0
        re       =  6371.2d0    
        gmagmo   =  0.311653d0
        eclipt   =  23.4415117563d0
        geoid(1) =  6378.1600d0
        geoid(2) =  ( geoid(1) / 6356.7746d0 )**2
        geoid(3) =  geoid(2) - 1.0d0
C
        do i=1,30
          uma(i) = am(i) * 1.6604D-24
        end do
        uma(19) = 9.1091D-28
c
      endif
      if( kinit.ge.0 .or. kinit.eq.-190 )then
C        
C  a.2/ UC190
C
        message  = 'Common block UC190 re-initialized'
        prop     =  0.2d0
C v1.10        stepx    =  0.4d0
C v1.10        stpmin   =  5.0d0
C V1.10        epskm    =  0.25d0
C v1.10        stplst   =  2.5d1
C v1.10        nxstp    =  400
        stepx    =  0.075d0
        stpmin   =  2.0d0
        umsq     =  1.0d0 - SQRT(0.5d0)
        upsq     =  1.0d0 + SQRT(0.5d0)
        uk2      =  -1.0d0
        uk3      =  5.0d0 - 3.0d0 * SQRT(2.0d0)
        epskm    =  0.2d0
        epsrel   =  6.0d-6
        stplst   =  1.5d1
        kmflg    =  0
        kum533   =  1
c
      endif
      if( kinit.ge.0 .or. kinit.eq.-192 )then
C        
C  a.3/ UC192
C
        message  = 'Common block UC192 re-initialized'
C V 1.10       dltlat   =  2.0d0
        xrmin    =  1.0d3
        xbmin    =  4.0d-5
        xtmin    =  COS( 5.0d0 * deg )
        xbmax    =  1.0d6
        epslon   =  0.5d-1
        epsfl    =  1.0d-3
        fvet     =  0.65d0
        pvet     =  3.452d0
        epsomeg  =  0.008d0
        dltlat   =  1.0d0
c
      endif
C        
C  b/ other
C
      if( kinit.ge.0 .or. kinit.eq.-110 )then
        message  = 'Common block UC110 re-initialized'
        k1st     = -1
        klmp     = -1
c
      endif
      if( kinit.ge.0 .or. kinit.eq.-120 )then
        message  = 'Common block UC120 re-initialized'
        nbrfl    = 0
        kurfl    = -1
c
      endif
      if( kinit.ge.0 .or. kinit.eq.-130 )then
        message  = 'Common block UC130 re-initialized'
        nbrsg    = 0
        kursg    = -1
c
      endif
      if( kinit.ge.0 .or. kinit.eq.-170 )then
        message  = 'Common block UC170 re-initialized'
        nsg      = 0
        kgp      = -1
c
      endif
      if( kinit.ge.0 .or. kinit.eq.-150 )then
        message  = 'Common block UC150 re-initialized'
        ntspec   = -1
        nnspec   = -1
        do i=1,24
          kflag(i)=1
        enddo
        kflag(25)=0
        kflag(26)=0
        do i=27,38
          kflag(i)=-1
        enddo
        kflag(39)=1
        do i=40,50
          kflag(i)=0
        enddo
C        
      endif
c
c
c
!      ifail    = iver
      ifail    = 0
C        
c
      if( stplst .lt. 4*stpmin )then
        ifail = -99002
        return
      endif
C      
C  c/ print the information
C
      if( kunit.gt.0 )then
!        write(kunit,1000,err=2000) iver/100, MOD(iver,100)
        call myversion( lver )
        write(kunit,1000,err=2000) trim(lver)
        if( kinit.ge.0 )then
          write(kunit,1001,err=2000)
        else
          write(kunit,1002,err=2000)
        endif
        write(kunit,1160,err=2000) geoid(1), geoid(1)/SQRT(geoid(2)), 
     :                    re, eclipt,
     :                    gmagmo
        if( kum533 .ne. 1) write(kunit,1195, err=2000) kum533
        write(kunit,1190,err=2000) stepx*re, stplst, stpmin, 
     :                    epskm*1.0d3,
     :                    epsrel*1.0d2
        if( kmflg .ne. 0 ) write(kunit,1191, err=2000) kmflg 
        write(kunit,1192,err=2000) ASIN(xtmin)/deg, xrmin/re
        write(kunit,1193,err=2000) epsfl*1.0d2, epslon
        if( kinit.lt.0 )write(kunit,1200,err=2000) message
      endif
c
c
!      if( kinit.lt. 0 ) ifail    = 10000 + iver
      if( kinit.lt. 0 ) ifail    = 1
c
      return
 2000 ifail=-99001
C
C
! 1000 FORMAT(
!     :    /36x,'SPENVIS',
!     :    /33x,'UNILIB Library',/34x,'Version',i3,'.',i2.2,
!     :    /20x,'Generated by SPENVIS team'/)
 1000 FORMAT(
     :    /36x,'SPENVIS',
     :    /33x,'UNILIB Library',/33x,'Version ',a,
     :    /27x,'Generated by SPENVIS team'/)
 1001 FORMAT(/' --- Main control parameters ---')
 1002 FORMAT(/' --- Main control parameters (modified) ---')
! 1003 FORMAT(/' --- Main control parameters (RERMM specific) ---')
 1160 FORMAT(/6x,'UC160 (general constants):',/22x,'Geoid major axis =',
     :     f8.2,4x,' km',/22x,'Geoid minor axis =',f8.2,4x,' km',/16x,
     :                  'Mean Earth radius (Re) =',f8.2,4x,' km',/18x,
     :                     'Ecliptic inclination =',f12.6,' deg',/17x,
     :                      'McIlwain Earth dipole =',f12.6,' G/Re^3')
 1190 FORMAT(/6x,'UC190 (field line tracing):',/21x,
     :      'Maximum step size =',f7.1,5x,' km at one Re',/14x,
     :                 'Maximum "last" step size =',f7.1,5x,' km',/21x,
     : 'Minimum step size =',f7.1,5x,' km',/20x,'Altitude precision =',
     :     f6.0,6x,' m',/14x,'Magnetic field precision =',f10.4,2x,' %')
 1191 FORMAT(24x,'Flag for UL240 =',i5)
 1192 FORMAT(/6x,'UC192 (magnetic field evaluation):',/13x,
     :               'Max. geomagnetic latitude =',f7.1,5x,' deg',/24x,
     :                                'Minimum radius =',f9.3,3x,' Re')
 1193 FORMAT(12x,'(drift shell tracing):',/24x,
     :               'Precision on L =',f11.5,1x,' %',/12x,
     :               'Precision on the longitude =',f8.2,4x,' deg')
 1195 FORMAT(24x,'Flag for UM533 =',i5)
 1200 FORMAT(/' --- Note ---',//6x,a33)

      END
C----------------------------------------------------------------------
      SUBROUTINE UT991
     :          (kunit, ifail)
C
C!    Print the magnetic field line stored in UC170
C
      INCLUDE 'structure.h'
C
C     INTERFACE
C
        INTEGER*4     ifail,kunit           
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
      CHARACTER*20 atxt(5),antdet
      integer*4    i,i1,i2
      data antdet/'-- Not determined --'/
C
C     CODE          
C
      ifail=0
      if( kunit.le.0 .or. nsg.le.0 ) then
        ifail=-99101
        return
      endif
C           
      do i=1,5
        atxt(i)=antdet
      enddo
      if( mlab%linv  )write(atxt(1),1001,err=2000) mlab%finv/re
      if( mlab%lbmp  )write(atxt(2),1002,err=2000) mlab%fbmp
      if( mlab%lkauf )write(atxt(3),1003,err=2000)
     :                                        mlab%fkauf/sqrt(re*gmagmo)
      if( mlab%llmi  )write(atxt(4),1004,err=2000) mlab%flmi
      if( mlab%lalp0 )write(atxt(5),1005,err=2000) mlab%falp0
C
      write(kunit,1000,err=2000) nsg,atxt
C
      if(mlin%ind%jend .gt. 0 .and. mlin%ind%jbeg .gt. 0) then
        i1=mlin%ind%jbeg
        i2=mlin%ind%jend
        ifail=mlin%ind%jend-mlin%ind%jbeg+1
        write(kunit,1020,err=2000) ifail,
     :      mlin%equat%coord%radius,mlin%equat%coord%colat,
     :      mlin%equat%coord%elong,mlin%equat%rcurv,mlin%equat%b%dnrm,
     :      mlin%equat%b%rho,mlin%equat%b%theta,mlin%equat%b%phi,
     :      mlin%footpn%radius,mlin%footpn%colat,mlin%footpn%elong,
     :      mlin%footps%radius,mlin%footps%colat,mlin%footps%elong,
     :      mlin%drift%rho,mlin%drift%theta,mlin%drift%phi,
     :      mlin%drift%dnrm,mlin%dtdft
      elseif( kgp .gt. 0) then
        i1=1
        i2=kgp
        ifail=i2-i1+1
        write(kunit,1050,err=2000) ifail
      else
        ifail=-99102
        return
      endif
C
      write(kunit,1025,err=2000)
      do i=i1,i2
          write(kunit,1030,err=2000) mele(i)%beg%coord%radius,
     :      mele(i)%beg%coord%colat,mele(i)%beg%coord%elong,
     :      mele(i)%beg%b%dnrm,mele(i)%beg%b%rho,mele(i)%beg%b%theta,
     :      mele(i)%beg%b%phi,mele(i)%beg%rcurv,mele(i)%arcl,
     :      mele(i)%csalp,mele(i)%dtbnd,mele(i)%rkstp
      enddo
C
      write(kunit,1040,err=2000)
C
      return
 2000 ifail=-99103
C
 1000 format('# Magnetic field line segment alone, stored in UC170,',
     :        i4,' pts',
     :      /'#      Integral invariant function =',a20,    
     :      /'# Magn. field intensity at the MPs =',a20,
     :      /'#   Kaufmann adiabatic invariant K =',a20, 
     :     /'#           McIlwain''s L parameter =',a20,         
     :      /'#           Equatorial pitch angle =',a20)
 1001 format(1x,f9.5,' Re       ')
 1002 format(1x,f9.5,' G        ')
 1003 format(1x,f9.5,' (M Re)1/2')
 1004 format(1x,f9.5,' Re       ')
 1005 format(1x,f9.5,' deg      ')
 1006 format(1x,f9.5,' M/Re ??? ')
 1007 format(1x,f9.5,' ?????????')           
 1020 format('##### The magnetic field line includes',i4,' pts ##### ',
     :/'# Equatorial pnt.: (',f8.1,';',f8.4,';',f8.4,') RC=',f8.1,' km',
     :/'#              B :',f7.4,' G, [',f7.4,';',f7.4,';',f7.4,']',
     :/'# North. footpnt.: (',f8.1,';',f8.4,';',f8.4,')',
     :/'# South. footpnt.: (',f8.1,';',f8.4,';',f8.4,')',
     :/'# Drift inform.  : [',f9.4,';',f9.4,';',f9.4,'] ',2f10.4)
 1025 format('#   radius    colat    elong   B       rho     theta ',
     :       '  phi     rcurv       arcl',/'#               ',
     :       '        csalp   dtbnd ',
     :       '   rkstp 4 rkstp 2 rkstp 3        ')
 1030 format(1x,f9.1,f9.4,f9.4,f7.4,3f8.4,f9.1,f11.1,/f29.3,f11.4,3f8.5)
 1040 format('#')
 1050 format('##### Temporary magnetic field line,',i4,' pts ##### ',
     :       5(/'#'))
C
      END          
C----------------------------------------------------------------------
      SUBROUTINE UT992
     :          (kunit, ifail)
C
C!    Print the magnetic field drift shell
C
      INCLUDE 'structure.h'
C
C     INTERFACE
C
        INTEGER*4     ifail,kunit           
C
      COMMON /UC110/ mlbl, k1st, klmp
C
             INTEGER*4     k1st, klmp
             TYPE(zlbl) :: mlbl
C
      COMMON /UC120/ nbrfl, kurfl, mfl
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
      CHARACTER*20 atxt(5),antdet
      integer*4    i,k
      TYPE(zfln) :: mlin
      data antdet/'-- Not determined --'/
C
C     CODE          
C
      ifail=0
      if( kunit.le.0 .or. nbrfl.le.0 ) then
        ifail=-99201
        return
      endif
C           
      do i=1,5
        atxt(i)=antdet
      enddo
      if( mlbl%linv  )write(atxt(1),1001,err=2000) mlbl%finv/re
      if( mlbl%lbmp  )write(atxt(2),1002,err=2000) mlbl%fbmp
      if( mlbl%lkauf )write(atxt(3),1003,err=2000) 
     :                                       mlbl%fkauf/sqrt(re*gmagmo)
      if( mlbl%llmi  )write(atxt(4),1004,err=2000) mlbl%flmi
      if( mlbl%lalp0 )write(atxt(5),1005,err=2000) mlbl%falp0
C
      write(kunit,1000,err=2000) nbrfl,atxt
C
      k=k1st
  100 continue
        if( k.le.0 )then
          ifail=-99202
          return
        endif
        mlin=mfl(k)
C
        i=mlin%ind%jend-mlin%ind%jbeg+1
        write(kunit,1020,err=2000) i,
     :      mlin%equat%coord%radius,mlin%equat%coord%colat,
     :      mlin%equat%coord%elong,mlin%equat%rcurv,mlin%equat%b%dnrm,
     :      mlin%equat%b%rho,mlin%equat%b%theta,mlin%equat%b%phi,
     :      mlin%footpn%radius,mlin%footpn%colat,mlin%footpn%elong,
     :      mlin%footps%radius,mlin%footps%colat,mlin%footps%elong,
     :      mlin%drift%rho,mlin%drift%theta,mlin%drift%phi,
     :      mlin%drift%dnrm,mlin%dtdft
C
        write(kunit,1025,err=2000)
        do i=mlin%ind%jbeg,mlin%ind%jend
          write(kunit,1030) mseg(i)%beg%coord%radius,
     :      mseg(i)%beg%coord%colat,mseg(i)%beg%coord%elong,
     :      mseg(i)%beg%b%dnrm,mseg(i)%beg%b%rho,mseg(i)%beg%b%theta,
     :      mseg(i)%beg%b%phi,mseg(i)%beg%rcurv,mseg(i)%arcl,
     :      mseg(i)%csalp,mseg(i)%dtbnd,mseg(i)%rkstp
        enddo
        k=mlin%keast
      if( k.ne.k1st )goto 100   
C
      write(kunit,1040,err=2000)
c
      return
 2000 ifail=-99203
C
 1000 format('# Magnetic drift shell stored in UC110, UC120 and UC130,',
     :       i4,' fld. lines'
     :      /'#      Integral invariant function =',a20,    
     :      /'# Magn. field intensity at the MPs =',a20,
     :      /'#   Kaufmann adiabatic invariant K =',a20, 
     :     /'#           McIlwain''s L parameter =',a20,         
     :      /'#           Equatorial pitch angle =',a20)
 1001 format(1x,f9.5,' Re       ')
 1002 format(1x,f9.5,' G        ')
 1003 format(1x,f9.5,' (M Re)1/2')
 1004 format(1x,f9.5,' Re       ')
 1005 format(1x,f9.5,' deg      ')
 1006 format(1x,f9.5,' M/Re ??? ')
 1007 format(1x,f9.5,' ?????????')           
 1020 format('##### The magnetic field line includes',i4,' pts ##### ',
     :/'# Equatorial pnt.: (',f8.1,';',f8.4,';',f8.4,') RC=',f8.1,' km',
     :/'#              B :',f7.4,' G, [',f7.4,';',f7.4,';',f7.4,']',
     :/'# North. footpnt.: (',f9.1,';',f9.4,';',f9.4,')',
     :/'# South. footpnt.: (',f9.1,';',f9.4,';',f9.4,')',
     :/'# Drift inform.  : [',f9.4,';',f9.4,';',f9.4,'] ',2f10.4)
 1025 format('#   radius    colat    elong   B       rho     theta ',
     :       '  phi     rcurv       arcl',/'#               ',
     :       '        csalp   dtbnd ',
     :       '   rkstp 4 rkstp 2 rkstp 3        ')
 1030 format(1x,f9.1,f9.4,f9.4,f7.4,3f8.4,f9.1,f11.1,/f29.3,f11.4,3f8.5)
 1040 format('#')

      end    
C----------------------------------------------------------------------
      SUBROUTINE UT993
     :          (kusumma, kushell, ifail)
C
C!    Store a magnetic field drift shell
C
c
c     rem: based on ut992
c
c    The subroutine UT993 prints two files with data relative to 
c    the magnetic drift
c    shell stored in the common blocks 
c    UC110, UC120 and UC130. The arguments kusumma and kushell
c    are used to specify the file units. Both files are printed in
c    ASCII format.
c
c    The body of the first file (kusumma) contains one text line
c    per magnetic field line which includes the location of
c    the foot and mirror points, the location where the magnetic
c    field intensity is minimum on the field line, the value of
c    minimum intensity, and the total length between the two
c    mirror points.  The header of the file contains the exact
c    definitions of each column present in the body part.
c
c    The second file (kushell) contains all the magnetic field line
c    segments of the drift shell.  The body of the file is divided
c    in blocks, each block containing a field line segment.
c    The first line of each block contains the number of points
c    listed in the block.
c    Each following text line is relative to a single geographic point along 
c    the field line. Only the points located between the two mirror
c    points are included. The text line includes the geographic location,
c    the intensity and spherical components of the magnetic field,
c    the local radius of curvature of the field line, and
c    the arc length along the field line.
c    See the header of the file for an exact
c    definitions of each column present in the body part.
c
c    The header of both files has a similar format that includes
c    the number of columns (first integer), 
c    the number of field line segments (third integer),
c    and the column definitions. Note that the format has been choosen
c    to be easily read by IDL or other graphical tools.
c
c
      INCLUDE 'structure.h'
C
C     INTERFACE
C
        INTEGER*4     ifail,kusumma,kushell           
C
      COMMON /UC110/ mlbl, k1st, klmp
C
             INTEGER*4     k1st, klmp
             TYPE(zlbl) :: mlbl
C
      COMMON /UC120/ nbrfl, kurfl, mfl
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
      CHARACTER*80 title
      INTEGER*4    k, i, n, i0, i1
      TYPE(zfln) :: mlin
      TYPE(zseg) :: mnmp,msmp
C
C     CODE          
C
      ifail=0
      if( kusumma.le.0 .or. kushell.le.0 ) then
        ifail=-99301
        return
      endif
      if( nbrfl.le.0 .or. 
     :    .not. (mlbl%lbmp .and. mlbl%llmi)) then
        ifail=-99302
        return
      endif
C
C
C           
      write (title,2000) mlbl%fbmp, mlbl%flmi
      write (kusumma,2020,err=4000) 15,0,nbrfl,0,0
      write (kusumma,2040,err=4000) 
      write (kushell,2020,err=4001)  9,0,nbrfl,1,0
      write (kushell,2030,err=4001)
C
C
C
      k=k1st
  100 continue
          if( k.le.0 )then
            ifail=-99303
            return
          endif
          mlin  = mfl(k)
          i0    = mlin%ind%jmirpn
          i1    = mlin%ind%jmirps
          if( k.le.0 .or. i0 .le. 0 .or. i1 .le. 0 )then
            ifail=-99304
            return
          endif
          if( i1 .lt. i0 )then
            i=i1
            i1=i0
            i0=i
          endif
C
          n     = i1-i0+1
          mnmp  = mseg(mlin%ind%jmirpn)
          msmp  = mseg(mlin%ind%jmirps) 
          write (kushell,2050,err=4001) n
          write (kusumma,2060,err=4000) mnmp%beg%coord%radius-re,
     :      90.0d0-mnmp%beg%coord%colat, mnmp%beg%coord%elong,
     :      msmp%beg%coord%radius-re, 90.0d0-msmp%beg%coord%colat,
     :      msmp%beg%coord%elong,mlin%equat%coord%radius-re,
     :      90.0d0-mlin%equat%coord%colat,mlin%equat%coord%elong,
     :      mlin%equat%b%dnrm,90.0d0-mlin%footpn%colat,
     :      mlin%footpn%elong,90.0d0-mlin%footps%colat,
     :      mlin%footps%elong,abs(mnmp%arcl-msmp%arcl)/re
c
c
c
          do i=i0,i1
            write (kushell,2070,err=4001) mseg(i)%beg%coord%radius/re,
     :        mseg(i)%beg%coord%colat*deg,mseg(i)%beg%coord%elong*deg,
     :        mseg(i)%beg%b%dnrm,mseg(i)%beg%b%rho,mseg(i)%beg%b%theta,
     :        mseg(i)%beg%b%phi,mseg(i)%arcl,mseg(i)%beg%rcurv
          enddo
          k=mlin%keast
      if( k.ne.k1st )goto 100   
C
      return
 4000 ifail=-99305
      return
 4001 ifail=-99306
C
C
 2000 format(' Visualisation of the magnetic drift shell Bm =',f8.5,
     :       ' Gauss, L =',f9.5,' Re. ')
 2020 format(5i8)
 2030 format(' Radial distance (Re)',
     :       /' Colatitude (rad)',
     :       /' Longitude (rad)',
     :       /' !8B!3 (Gauss)',
     :       /' !8B!dR!n!3 (Gauss)',
     :       /' !8B!d!7h!n!3 (Gauss)',
     :       /' !8B!d!7u!n!3 (Gauss)',
     :       /' !8s!3 (km)',
     :       /' Radius of curvature (km)')
 2040 format( ' Altitude (km) of the northern mirror point',
     :       /' Latitude (deg) of the northern mirror point', 
     :       /' Longitude (deg) of the northern mirror point', 
     :       /' Altitude (km) of the southern mirror point',
     :       /' Latitude (deg) of the southern mirror point', 
     :       /' Longitude (deg) of the southern mirror point', 
     :       /' Altitude (km) of the minimum !8B!3 location',
     :       /' Latitude (deg) of the minimum !8B!3 location', 
     :       /' Longitude (deg) of the minimum !8B!3 location', 
     :       /' !8B!3 (Gauss) at the minimum !8B!3 location', 
     :       /' Latitude (deg) of the northern foot point', 
     :       /' Longitude (deg) of the northern foot point', 
     :       /' Latitude (deg) of the southern foot point', 
     :       /' Longitude (deg) of the southern foot point',
     :       /' Length (Re) of the magnetic field line segment')
 2050 format(i10)
 2060 format(3(f11.3,f8.3,f9.3),f9.6,2(f8.3,f9.3),f9.4)
 2070 format(f9.4,2f9.5,4f10.6,1p,2e13.5,0p)
C
      end    
C----------------------------------------------------------------------
      SUBROUTINE UT998
     :          (mgp, rc, mnr, mb, ifail)
C
C!    Evaluate the magnetic field intensity and the radius of curvature
C
C
      INCLUDE 'structure.h'
cDEC$ IF DEFINED (_x86_)
cDEC$ ATTRIBUTES DLLEXPORT :: UT998
cDEC$ ENDIF
C
      EXTERNAL  UF423, UF425, UM530, UT999
C
C     INTERFACE
C
        TYPE(zgeo) :: mgp
        REAL*8        rc
        TYPE(zvec) :: mnr
        TYPE(zvec) :: mb
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
        INTEGER*4     k0, kold, k530
        REAL*8        dum(2), dstep, scale
        TYPE(zvec) :: mxrna, mxrnb
        TYPE(zvec) :: mdum
C
C     CODE          
C
      ifail                     = 0
      k530                      = 0
      kold                      = kgp
      mnr%rho                   = -999.0d0
      mnr%theta                 = -999.0d0
      mnr%phi                   = -999.0d0
      mnr%dnrm                  =    0.0d0
      mb%rho                    = -999.0d0
      mb%theta                  = -999.0d0
      mb%phi                    = -999.0d0
      mb%dnrm                   =    0.0d0
C
      CALL UT999 (mgp, rc)
C
      CALL UM530 (mgp, mb, k530)
      if( k530 .lt. 0 .and. k530.ne.-53002 .and.
     :    k530.ne.-53003 .and. k530.ne.-53004 )then
        ifail=-99803
        return
      endif
C
      if( MAX(nsg,kgp) .gt. nx170-10 )then
        ifail                   = -99801
        return
      endif
C
      k0                        = nx170-10
      mele(k0)%beg%coord        = mgp
      mele(k0)%arcl             = 0.0d0
      mele(k0)%beg%rcurv        = rc
      mele(k0)%beg%b            = mb
C 
      scale                     = 0.2d0
      kgp                       = k0
      CALL UF423(scale, rc, dum(1), dum(2), ifail)
      if(ifail .lt. 0 .and. ifail.ne.-42305 .and.
     :   ifail.ne.-42306 .and. ifail.ne.-42304 ) then
        scale                     = -0.2d0
        kgp                       = k0
        CALL UF423(scale, rc, dum(1), dum(2), ifail)
        if(ifail.lt.0 .and. ifail.ne.-42305 .and.
     :   ifail.ne.-42306 .and. ifail.ne.-42304  ) then
          CALL UT999 (mgp, rc)
          ifail = -99802
          return
        endif
        dstep          = mele(kgp)%arcl - mele(k0)%arcl
        CALL UF425 (mele(k0)%beg%coord,
     :            mele(k0)%beg%b, 
     :            mele(kgp)%beg%coord,
     :            mele(kgp)%beg%b,
     :            dstep, mnr, mdum)
        ifail = -99804
      else 
        ifail = 0
        dstep          = mele(kgp)%arcl - mele(k0)%arcl
        CALL UF425 (mele(k0)%beg%coord,
     :            mele(k0)%beg%b,
     :            mele(kgp)%beg%coord, 
     :            mele(kgp)%beg%b,
     :            dstep, mxrna, mdum)
C 
        scale                     = -0.2d0
        kgp                       = k0
        CALL UF423(scale, rc, dum(1), dum(2), ifail)
        if(ifail.lt.0 .and. ifail.ne.-42305 .and.
     :     ifail.ne.-42306 .and. ifail.ne.-42304 )then
          mnr = mxrna
          ifail = -99804
        else
        ifail = 0
          dstep          = mele(kgp)%arcl - mele(k0)%arcl
          CALL UF425 (mele(k0)%beg%coord,
     :            mele(k0)%beg%b, 
     :            mele(kgp)%beg%coord,
     :            mele(kgp)%beg%b,
     :            dstep, mxrnb, mdum)
C
          mnr%rho    = ( mxrna%rho + mxrnb%rho ) * 0.5d0     
          mnr%theta  = ( mxrna%theta + mxrnb%theta ) * 0.5d0     
          mnr%phi    = ( mxrna%phi + mxrnb%phi ) * 0.5d0
          mnr%dnrm   =  SQRT(mnr%rho**2 + mnr%theta**2 + mnr%phi**2)   
        endif 
      endif
C
      rc         = 1.d0 / mnr%dnrm
C
      kgp        = kold
C
C
      END
C----------------------------------------------------------------------
      SUBROUTINE UT999
     :          (mgp, rc)
C
C!    Radius of curvature in a dipolar magnetic field 
C
      INCLUDE 'structure.h'
C
C     INTERFACE
C
        TYPE(zgeo) :: mgp
        REAL*8        rc
C
      COMMON /UC140/  mint, mext, msun 
C
        TYPE(zimf) :: mint
        TYPE(zemf) :: mext
        TYPE(zsun) :: msun       
C
      COMMON /UC160
     :               /pi, deg, re, gmagmo, eclipt, geoid, uma
C
        REAL*8        pi, deg, re, gmagmo, eclipt, geoid(3), uma(30)
C
C     VARIABLES
C
        REAL*8        sclat, axgp, aygp, azgp, cml, sml
        REAL*8        cml2
C
C     CODE
C
C  a/ Evaluate the sinus & cosinus of the magnetic latitude
C
      sclat = SIN( mgp%colat*deg )
      axgp  = sclat * COS( mgp%elong*deg )
      aygp  = sclat * SIN( mgp%elong*deg )
      azgp  = COS( mgp%colat*deg )
      cml   = ( axgp * COS( mint%elong * deg )
     :         + aygp * SIN( mint%elong * deg ) )
     :        * SIN( mint%colat * deg ) 
     :        + azgp * COS( mint%colat * deg )
      sml   = SQRT( 1.0d0 - cml**2 )
C      
C  b/ Radius of curvature
C ****                     ****
C **** TO BE CHECKED AGAIN ****
C ****                     ****
      cml2  = cml**2
      rc    = mgp%radius * SQRT( 1.0d0 + 3.0d0 * cml2 )**3 /
     :         ( 3.0d0 * sml * ( 1.0d0 + cml2 ) )
C
C
      END
C----------------------------------------------------------------------
