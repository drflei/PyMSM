# 1 "ua630.f"
      subroutine ua630
     :           (mpos, tmass, dens, temp, ifail)
c
c
      include 'structure.h'
cDEC$ IF DEFINED (_x86_)
cDEC$ ATTRIBUTES DLLEXPORT :: UA630
cDEC$ ENDIF
c
      TYPE(zgeo) :: mpos
      real*8        tmass, dens(30), temp(3)
      integer*4     ifail
C
*
*   mpos   location
*   tmass  total mass density
*   dens   number density
*   temp   temperature neutral, elec, ions
*   ifail
*
C
      COMMON /UC150/matm, ntspec, nnspec, kspec, kflag
C
            TYPE(zatm) :: matm
            INTEGER*4     ntspec, nnspec, kspec(30)
            INTEGER*4     kflag(50)
C
      COMMON /UC160/ pi, deg, re, gmagmo, eclipt, geoid, uma
C
        REAL*8        pi, deg, re
        REAL*8        gmagmo
        REAL*8        eclipt, geoid(3), uma(30)
C           
C
      TYPE(zgeo) :: mgde   
      integer*4     i
      real*8        tmneut, tmion, xalt, parm
c
      ifail=0
      do i=1,30
        dens(i)=0.0d0
      enddo
      tmneut = 0.0d0
      tmion  = 0.0d0
      temp(3)= 300.0d0
      temp(1)= 300.0d0
      temp(2)= 300.0d0
c
c
c
      call um535( mpos, mgde)
      if( mgde%radius .lt. re )then
        ifail=-63001
        return
      endif
c
c
c
      if( matm%katm .eq. 1 )then
c       msis
        call ua632(mgde,nnspec,kspec(1),dens,temp(1),ifail)
        if(ifail.lt.0)return
        call ua634(dens,nnspec,kspec(1),0,dens, tmneut,ifail)
        if(ifail.lt.0)return
         elseif( matm%katm .eq. 2 )then
c       a&f neutre
        call ua631(mgde,nnspec,kspec(1),dens,ifail)
        if(ifail.lt.0)return
        call ua634(dens,nnspec,kspec(1),0,dens, tmneut,ifail)
        if(ifail.lt.0)return        
      elseif( matm%katm .eq. 3 )then
c       allen
C     ! SIMPLE DENSITY DISTRIBUTION FUNCTION
C       PROFILE IS A THREE-INTERVAL FIT TO A
C       TABLE OF DENSITIES IN "ASTROPHYSICAL QUANTITIES" (ALLEN, 1985).C 

        xalt=mgde%radius-re
C
        if( xalt.lt.120.0d0 )then
          tmneut=exp((296.3923182d0-xalt)/6.553140330D0)
        elseif( xalt.lt.1000.0d0)then
          tmneut=exp( 60.65982247d0-7.048075391D0*log(xalt) )
        else
          tmneut=exp( 30.67634731d0-2.707522787d0*log(xalt))
        endif          
      elseif( matm%katm .eq. 4 )then
c       pfitzer
        xalt=mgde%radius-re
        parm=0.99d0+0.518d0*sqrt((matm%f107+matm%f107a)/110.0d0)
c
        if( xalt.le.120.0d0 )then
          tmneut = 2.7d-11 * exp((120.0d0-xalt)/6.4d0)
        elseif( xalt.le.1000.0d0 )then
          tmneut = 2.7d-11 * exp((120.0d0-xalt)/parm/sqrt(xalt-103.0d0))
        else
          tmneut =2.7d-11 * exp(-880.0d0/parm/sqrt(897.0d0)
     :                          -1.0d3/150.0d0*log(xalt/1.0d3))
        endif
c        
      endif
c
c
c
      if( matm%kion .eq. 1)then
C       iri
        i=ntspec-nnspec
        call ua633(mgde,i,kspec(nnspec+1),dens, temp(3), temp(2), ifail)
        if(ifail.lt.0)return
        call ua634(dens,i,kspec(nnspec+1),0,dens, tmion,ifail)
        if(ifail.lt.0)return        
      elseif( matm%kion .eq. 2 )then
c       a&f ion
        i=ntspec-nnspec
        call ua631(mgde,i,kspec(nnspec+1),dens,ifail)
        if(ifail.lt.0)return
        call ua634(dens,i,kspec(nnspec+1),0,dens, tmion,ifail)
        if(ifail.lt.0)return        
      elseif( matm%kion .eq. 3 )then
c       pfitser
        xalt=mgde%radius-re
        parm=90.0d0+0.6d0*matm%f107
        tmion=1.67d-24*10.0d0**(-0.3145d0*mgde%radius/re+3.9043d0+
     :        (1.27d-3*parm-6.35d-2)*exp(-(xalt/re-1.0d0)/1.5d0))
      endif
c
c
c
      tmass = tmneut+tmion
c
c
c
      end
C----------------------------------------------------------------------
      SUBROUTINE UA631 
     :           (mgde, nsp, ksp, dens, ifail)
C
C!    ANDERSON AND FRANCIS ATMOSPHERE
C
      include 'structure.h'
C
      TYPE(zgeo) :: mgde
      real*8        dens(30)
      integer*4     nsp, ksp(nsp), ifail
C
*   mgde geodetic coordinate
*   n, ksp constituents 
*   dens density number
*
C
C
      integer*4     k1to1(30)
C
C
C OLD (blxtra):      REAL*8 FUNCTION AFATMO(N, X)
C X=ALTITUDE IN KM OR THE REDUCED HEIGHT, N=CONSTITUENT NUMBER
C E 1/O2 2/N2 3/O 4/N 5/HE 6/H 7/O2+ 8/NO+ 9/O+ 10/HE+ 11/H+ 12/
C FIT TO ANDERSON AND FRANCIS ATMOSPHERE- FROM NOTE BY WALT AND
C NEWKIRK DATED 12/16/63
C
C
      COMMON /UC160/ pi, deg, re, gmagmo, eclipt, geoid, uma
C
        REAL*8        pi, deg, re
        REAL*8        gmagmo
        REAL*8        eclipt, geoid(3), uma(30)
C           
      REAL*8 A(12), B(12), ALT, rap, REDHT
      INTEGER*4 N,k,inx
c
      data k1to1/ 7, 6, -1, -1, -1, -1, -1, 3, 2, -1,
     :           -1, -1, -1, -1, -1, 4, 5, -1, 1, 12,
     :           11, -1, -1, 10, -1, -1, 9, 8, -1, -1/
C
C
C
C            
C
      DATA A / 1.5917534D1,          ! n0=8.2e6
     &         2.569709D1, 2.714195D1, 2.428401D1, 1.901951D1,
     &         1.527661D1, 1.154241D1, 1.720580D1, 1.753069D1,
     &         1.980918D1, 6.964002D0, 2.769234D0 /
      DATA B / -5.5461649D1,         ! scale height 115km
     &         -2.403941D2, -2.103555D2, -1.199667D2, -1.050822D2,
     &         -2.996343D1, -7.586136D0, -2.038100D2, -1.904198D2,
     &         -1.010177D2, 2.018458D1, 4.143615D1 /
C
C
      ifail=0
      alt=mgde%radius-re
      rap = alt / (mgde%radius)  
c
c
c
      do k = 1, nsp
        inx=ksp(k)
        if( inx.le.0 .or. inx.gt.30 )then
          ifail=-63102
c       invalid value for inx
          return
        endif
        dens(inx)=0.0d0
c
        n=k1to1(inx)
        if( n.le.0 )then
          ifail=-63101
C         atmospheric component not included in the model
          return
        endif
C
        if( n .eq. 1 )then
          if( rap .lt. 0.111883D0) then
            redht = A(N) + B(N) * rap
          elseif( rap .le. 0.255834D0)then
            redht = 10.808814D0 - 9.800442D0 * rap
          elseif( rap .le. 0.574230D0)then
            redht = 9.370167D0 - 4.177076D0 * rap
          else
            redht=14.211403-12.607905*rap
          endif
        elseif( n .lt. 11)then
          redht = A(N) + B(N) * rap
        elseif( n .eq. 11)then
          if(rap .lt. 0.12410D0)then
            redht = A(N) + B(N) * rap
          elseif(rap .ge. 0.30384D0)then
            redht = 16.09293D0 - 32.89545D0 * rap
          else
            redht = 11.79650D0 - 18.75496D0 * rap
          endif
        elseif( rap .le. 0.120126D0)then
          redht = A(N) + B(N) * rap
        elseif( rap .lt. 0.22987D0)then 
          redht = 7.23661D0 + 4.24706D0 * rap
        elseif( rap .lt. 0.43971D0)then
          redht = 8.911084D0 - 3.037239D0 * rap
        elseif( rap .lt. 0.61083D0)then
          redht = 10.208135D0 - 5.987027D0 * rap
        else
          redht=14.415333d0 - 12.87470d0 * rap
        endif
c
        dens(inx) = dexp(redht)
c                 ! keep in reasonable range
c
      enddo
c
      END
C----------------------------------------------------------------------
      subroutine ua632
     :           (mgde, nsp, ksp, dens, tmpn, ifail)   
C
C!    MSISE-90 ATMOSPHERE
C
      include 'structure.h'
C
      TYPE(zgeo) :: mgde
      real*8        dens(30), tmpn
      integer*4     nsp, ksp(nsp), ifail
C
*   mgde geodetic coordinate
*   n, ksp constituents 
*   dens density number
*   tmpn temperature
*
C
C
      COMMON /UC150/matm, ntspec, nnspec, kspec, kflag
C
            TYPE(zatm) :: matm
            INTEGER*4     ntspec, nnspec, kspec(30)
            INTEGER*4     kflag(50)
C
      COMMON /UC160/ pi, deg, re, gmagmo, eclipt, geoid, uma
C
        REAL*8        pi, deg, re
        REAL*8        gmagmo
        REAL*8        eclipt, geoid(3), uma(30)           
C                        
      integer*4 k1to1(30)
c                                     
      real*4 ALT, GLAT, GLONG, STL, F107A, F107, AP(7), d(8), t(2)
      integer*4 mass
      integer*4 i,n
c
      data k1to1/  1,  4, -1,  6, -1, -1, -1,  7,  5, -1,
     :            -1, -1, -1, -1, -1,  3,  2, -1, -1, -1,
     :            -1, -1, -1, -1, -1, -1, -1, -1, -1, -1/
c
c
      data mass/ 48/
c
c
c
      ifail = 0
      alt=mgde%radius-re
      if( alt .lt. 0.0)then
        ifail=-63201
        return
      endif
c
      glat=90.0-mgde%colat
      glong=mgde%elong
      stl=(mgde%elong+matm%ut)/15.d0
      f107a=matm%f107a
      f107=matm%f107
      do i=1,7
        ap(i)=matm%apind(i)
      enddo
c
c      print*,'Unilib-UA632: ',matm%kday ,ALT, GLAT, GLONG
c      print*,'            : ', STL, F107A, F107
c      print*,'       (Ap) : ', ap
c      print '(15H    (switch) : ,13I4,/15x,13I4)', (kflag(i),i=1,26)
      call NWGTD6( matm%kday ,ALT, GLAT, GLONG, STL, F107A, F107, AP,
     :             MASS, D, T, kflag(1) )
c      print*,' -----------: ',mass,d,t 
c
      tmpn = t(2)
      do i =1, nsp
        if( ksp(i).le.0 .or. ksp(i).gt.30 ) then
          ifail=-63202
          return
        endif
        n=k1to1(ksp(i))
        if( n .le. 0 )then
          ifail=-63203
          return
        endif
        dens(ksp(i))=d(n)
      enddo                        
c
      end  
C----------------------------------------------------------------------
      subroutine ua633
     :           (mgde, nsp, ksp, dens, tmpi, tmpe, ifail)
      include 'structure.h'
C
C     IRI
C
C
      TYPE(zgeo) :: mgde
      real*8        dens(30), tmpi, tmpe
      integer*4     nsp, ksp(nsp), ifail
C
*   mgde geodetic coordinate
*   nsp, ksp constituents 
*   dens density number
*   tmpn temperature
*
C
C
      COMMON /UC140/  mint,mext,msun
C      
        TYPE(zimf) :: mint
        TYPE(zsun) :: msun
        TYPE(zemf) :: mext                   
C          
C
      COMMON /UC150/matm, ntspec, nnspec, kspec, kflag
C
            TYPE(zatm) :: matm
            INTEGER*4     ntspec, nnspec, kspec(30)
            INTEGER*4     kflag(50)
C
      COMMON /UC160/ pi, deg, re, gmagmo, eclipt, geoid, uma
C
        REAL*8        pi, deg, re
        REAL*8        gmagmo
        REAL*8        eclipt, geoid(3), uma(30)           
C                        
      integer*4 k1to1(30)
c
c     (for iris12)
      real*4 stl, rz, deni(11,50), oarr(30)
      real*4 t1,t1000
      integer*4 k0, kd
      REAL*4 SALT, SLAT, SLON, LATC, LONGC
c     (internal)
      REAL*8 DALT
      INTEGER*4 I, J, n
      REAL*8 FLMAX, LATM, LONGM, THETAM, XC, YC, FL, RHHE, SUMZ
      REAL*8 SUMZ2, C, DENC, DENCA, RCOLAT, RELONG, RLAT, RLON 
      REAL*8 DENT(6), SCH(5), AG, BK
      integer*4 kam(6)
      REAL*8 DENIEX(8), TEMPI(3)
c
      data k1to1 / -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
     :             -1, -1, -1, -1, -1, -1, -1, -1,  1,  3,
     :              4, -1,  8,  2, -1, -1,  6,  5, -1, -1 /
      DATA BK /1.3805D-16/
      DATA t1,t1000,k0/1.0,1000.0,0/
      data kam/19,24,20,21,28,27/
c 
      ifail=0
c
c     begin "iriext"
c
      STL             = (matm%ut+mgde%elong)/15.0d0
      rz              = matm%rzss
      if( kflag(39) .lt. 0 ) rz= - matm%f107A  
      ag              = 981.0D0 * (re/7371.2D0)**2
      dalt            = mgde%radius-re
      kd              = -MATM%KDAY
      if (dalt .lt. 120.0D0) then
        do i=1,8
          deniex(i)   = 0.0D0
        end do
        do i=1,3
          tempi(i)    = 200.0D0
        end do                       
      else if (dalt .le. 1000.0D0) then
        salt          = dalt
        slat          = 90.0d0 - mgde%colat
        slon          = mgde%elong
        CALL IRIS12(kflag(27), k0, SLAT, SLON, rz, kd, STL,
     :              SALT, SALT, t1, DENI, OARR)
c
        deniex(1)     = deni(1,1) * 1.0E-6
        do j=2,8
          deniex(j)   = deni(j+3,1) * deniex(1) / 100.0D0
        end do
        do j=1,3
          tempi(j)    = deni(j+1,1)
        end do
      else
        flmax  = 5.6D0 - 0.46D0 * matm%fkpx
        rlat   = (90.0d0-mgde%colat) * deg
        rlon   = mgde%elong * deg
        rcolat = mint%colat * deg
        relong = mint%elong * deg
c
C       Geomagnetic coordinates
        latm   = dsin(rcolat) * dcos(relong) * dcos(rlat) * dcos(rlon)
     :         + dsin(rcolat) * dsin(relong) * dcos(rlat) * dsin(rlon)
     :         + dcos(rcolat) * dsin(rlat)
        latm   = PI / 2.0D0 - dacos(latm)
        xc     = dcos(rcolat) * dcos(relong) * dcos(rlat) * dcos(rlon)
     :         + dcos(rcolat) * dsin(relong) * dcos(rlat) * dsin(rlon)
     :         - dsin(rcolat) * dsin(rlat)
        yc     = -dsin(relong) * dcos(rlat) * dcos(rlon) +
     :            dcos(relong) * dcos(rlat) * dsin(rlon)
        longm  = datan2(yc, xc)
        if (yc .lt. 0.0D0) longm = PI * 2.0D0 + longm
c
C       Calculate L of dipole field line passing through point.
        fl = (1.0D0+dalt/re) / dcos(latm)**2
        if (fl .gt. flmax) then
c
C         Beyond the plasmapause!
c
          deniex(1)     = 10.0D0
          do i=2,8
            deniex(i)   = 0.0D0
          end do
            deniex(3)   = 10.0D0     
          do i=1,3
            tempi(i)    = 4000.0D0
          end do                                 
c
        else
c
C         Electron and ion densities and scale heights at 1000 km.
          CALL IRIS12(kflag(27), k0, SLAT, SLON, rz, Kd, STL,
     :                t1000, t1000, t1, DENI, OARR)
          dent(1)   = deni(1,1) * 1.0d-6
          do j=2,6
            dent(j) = deni(j+3,1) * dent(1) / 100.0d0
          end do
          rhhe      = dent(4) / dent(3)
          sumz      = -dent(1) * uma(19) / deni(4,1)
          sumz2     = dent(1) / deni(4,1)
          do j=2,6
            sumz    = sumz + dent(j) * uma(kam(j)) / deni(3,1)
            sumz2   = sumz2 + dent(j) / deni(3,1)
          end do
          do j=1,5
            sch(j)  = bk * deni(3,1) / ag / (uma(kam(j+1))-sumz/sumz2)
          end do
C  
C         Latitude where field line reaches 1000 km.
          thetam    = dacos(dsqrt(7371.2D0/re/fl))
          if (latm .lt. 0.0D0) thetam = - thetam
c                                     ! Same hemisphere.
C         Geocentric coordinates.
          latc = -dsin(rcolat) * dsin(thetam) * dcos(longm) +
     :            dcos(rcolat) * dcos(thetam)
          latc = PI / 2.0D0 - acos(latc)
          xc   = dcos(rcolat) * dcos(relong) * dsin(thetam) *dcos(longm)
     :         - dsin(relong) * dsin(thetam) * dsin(longm)
     :         + dsin(rcolat) * dcos(relong) * dcos(thetam)
          yc   = dcos(rcolat) * dsin(relong) * dsin(thetam) *dcos(longm)
     :         + dcos(relong) * dsin(thetam) * dsin(longm)
     :         + dsin(rcolat) * dsin(relong) * dcos(thetam)
          longc       = datan2(yc, xc)
          if (yc .lt. 0.0D0) longc = 2.0D0 * PI + longc
          CALL IRIS12(kflag(27), k0, LATC, LONGC, rz, kd, STL,
     :                t1000, t1000, t1, DENI, OARR)
          denc        = deni(1,1) / 1.0d6
          call ua639(FL, MATM%KDAY, matm%rzss,DENCA)
          c           = -7371.2D0 * dlog(denc/denca)
          deniex(1)   = dmax1(denca*dexp(-c/(re+dalt)), 0.0D0)
          deniex(2)   = dent(2) * dexp((1000.0D0-dalt)*7371.2D0/(re
     :                +dalt)/sch(1)*1.0D5)
          deniex(2)   = dmax1(deniex(2), 0.0D0) 
          do j=5,6
            deniex(j) = dent(j) * exp((1000.0D0-dalt)*7371.2D0/
     :                  (re+dalt)/sch(j-1)*1.0D5)
            deniex(j) = dmax1(deniex(j), 0.0D0) 
          end do
          deniex(3)   = (deniex(1) - deniex(2) - deniex(5) - deniex(6))
     :                / (1.0D0+rhhe)
          deniex(3)   = dmax1(deniex(3), 0.0D0)
          deniex(4)   = dmax1(deniex(3)*rhhe, 0.0D0)
          deniex(7)   = 0.0D0
          deniex(8)   = 0.0D0
          do i=1,3
            tempi(i)    = 4000.0D0
          end do                                 
        end if
c
      endif
c
c     end "iriext"
c
      do i =1, nsp
        if( ksp(i).le.0 .or. ksp(i).gt.30 ) then
          ifail=-63301
          return
        endif
        n=k1to1(ksp(i))
        if( n .le. 0 )then
          ifail=-63302
          return
        endif
        dens(ksp(i))=deniex(n)
      enddo                        
c
      tmpi=tempi(2)
      tmpe=tempi(3)
c
      end    
C----------------------------------------------------------------------
      subroutine ua634
     :           (dens,nsp,ksp, kflg, cs, tmass,ifail)
C
C     compute mass density from number density
C
      include 'structure.h'
C
      REAL*8       dens(30), cs(30), tmass
      INTEGER*4    nsp, ksp(30), kflg, ifail
C       
      COMMON /UC160/ pi, deg, re, gmagmo, eclipt, geoid, uma
C
        REAL*8        pi, deg, re
        REAL*8        gmagmo
        REAL*8        eclipt, geoid(3), uma(30)           
C
      integer*4    i
C
      ifail=0
C
      tmass=0.0d0
      if( kflg.eq.0 )then
        do i=1,nsp
          if(ksp(i).le.0 .or. ksp(i).gt.30)then
            ifail=-63401
            return
          endif
          tmass=tmass+dens(ksp(i))*uma(ksp(i))
        enddo
      elseif( kflg.eq.1 )then
        do i=1,nsp
          if(ksp(i).le.0 .or. ksp(i).gt.30)then
            ifail=-63401
            return
          endif
          tmass=tmass+dens(ksp(i))*uma(ksp(i))*cs(ksp(i))
        enddo
      elseif( kflg.eq.2 )then
        do i=1,nsp
          if(ksp(i).le.0 .or. ksp(i).gt.30)then
            ifail=-63401
            return
          endif
          tmass=tmass+dens(ksp(i))*cs(ksp(i))
        enddo
      else
        ifail=-63402
      endif
C
      end
C----------------------------------------------------------------------
      subroutine ua635
     :                 (dens,nsp,ksp,temp,debye,ifail)    
c
c     evaluate the debye length
c
C
      include 'structure.h'
C
      REAL*8       dens(30), debye, temp(3)
      INTEGER*4    nsp, ksp(30), ifail
*   temp   temperature neutral, elec, ions
C       
      integer*4 i,j
c
      ifail=0
      debye=0.0d0
      do i=1,nsp
        j=ksp(i)
        if(j.eq.19)then
          if(temp(2).le.0)then
      	  ifail=-63501
c           invalid temperature
            return
          endif
          debye =debye + dens(19) / temp(2)
        elseif(j.gt.19)then
          if(temp(3).le.0)then
      	  ifail=-63502
c           invalid temperature
            return
          endif
          debye =debye + dens(j) / temp(3)
        endif
      enddo
      if( debye.gt.0 )then
        debye = 1.0d0/sqrt(debye * 2.0998D2)
      else
        debye=0.0d0
      endif
c
      end
C----------------------------------------------------------------------
      subroutine ua636
     :                (ktrl, eng, mpos, mb, rmass, ifail)
C
C!  Evaluate a weighted atmospheric mass
c
      include 'structure.h'
c
      TYPE(zgeo) :: mpos
      TYPE(zvec) :: mb 
      integer*4     ktrl, ifail
      real*8        eng, rmass
*
* ktrl  control parameter
* eng   energy
* mpos  geographic location
* rmass weighted atmospheric mass
* ifail
*
c
c ktrl =  0 --> sum (nbr dens) * uma 
c         1 --> sum (nbr dens) * uma * ua637
c         2 --> sum (nbr dens) * ua637
c        11 --> sum (nbr dens) * uma * ua638
c        12 --> sum (nbr dens) * ua638
c
C
      COMMON /UC150/matm, ntspec, nnspec, kspec, kflag
C
            TYPE(zatm) :: matm
            INTEGER*4     ntspec, nnspec, kspec(30)
            INTEGER*4     kflag(50)
C
      integer*4 kflg, kmdcs, nispec
      real*8    debye, cs(30), dens(30), temp(3), tmass
c
      ifail = 0
      kmdcs = ktrl/10
      kflg  = ktrl-10
      rmass = 0.0d0
c
      call UA630 (mpos, tmass, dens, temp, ifail)
      if (ifail .lt. 0) return
c
      if (kflg .eq. 0) then
        rmass = tmass
      elseif (ntspec .ge. 1) then
        if (kmdcs .eq. 0) then
          nispec=ntspec-nnspec
          call ua635(dens,nispec,kspec(nnspec+1),temp,debye,ifail)
          if (ifail .lt. 0)  return
          call ua637(ntspec,kspec,eng,debye,cs,ifail)
          if (ifail .lt. 0)  return
          call ua634(dens,ntspec,kspec, kflg, cs, rmass,ifail)
          if (ifail .lt. 0)  return
        elseif (kmdcs .eq. 1) then
          call ua638(ntspec, kspec, cs, ifail)
          if (ifail .lt. 0)  return
          call ua634(dens,ntspec,kspec, kflg, cs, rmass,ifail)
          if (ifail .lt. 0)  return
        else
          ifail=-63602
c         illegal ktrl parameter
        endif
      else
        ifail=-63601
c       uncompatible Atmospheric model and ktrl parameter
      endif
c
      end
C----------------------------------------------------------------------
      subroutine ua637
     :                (nsp, ksp, eng, debye, cs, ifail)
c
c!  Compute proton cross sections
c
      include 'structure.h'
c
      integer*4 nsp, ksp(30), ifail
      real*8    eng, debye, cs(30)
c
c_
c_  Compute the cross sections for collisions between protons and
c_  atmospheric particles.
c_  Origin: V. Pierrard (1994)
c_  Adaptation by D. Heynderickx (July 94):
c_    electron cross sections have been removed
c_  Adaptation by D. Heynderickx (September 95):
c_    added NO+: set equal to O2+
c_
c_	E	in	Energy in eV
c_	FLAMBD	in	Coulomb logarithm, only used when L is greater
c_			than 16
c_	N	in	Particle type, see below
c_	YN	out	Cross section in mbarn, 10-27cm2
c_
C
      COMMON /UC160/ pi, deg, re, gmagmo, eclipt, geoid, uma
C
        REAL*8        pi, deg, re
        REAL*8        gmagmo
        REAL*8        eclipt, geoid(3), uma(30)           
C                        
      integer*4 k1to1(30)
c
      integer*4 i
      real*8 ex(5,149), sx(5,149)
      real*8 aneff(5,9), pot(5,9)
      real*8 abcd(4,13)
c
      real*8 x(156), y(156), x1(156), x2(156), y1(156), y2(156)
      real*8 pm(7)
      integer*4 ns(9)
      real*8 e, flambd, yn, xx, yen, yen1, yen2, ssi, sio, seff
      integer*4 j, l, na, nb, nc, km, nm, im, lm, nz
c
c
      integer*4 k,ik
c
      data k1to1/  1,  2,  3,  4,  5,  6,  7,  8,  9, 10,
     :            11, 12, 13, -1, -1, 15, 16, -1, 24, 25,
     :            17, 18, 19, 20, 21, 22, 23, 23, -1, -1/
c
c
c
      data pm/ 2.0d0, 12.0d0, 14.0d0, 16.0d0, 6.0d0, 28.0d0, 32.0d0/
      data ns/ 1, 5, 6, 7, 9, 0, 0, 0, 0/
c
c Energies and charge exchange cross sections measured (Tawara,1980)
c
      data (ex(1,i),i=1,140) /0.17,0.52,1.02,1.28,1.99, 
     &  3.10, 3.37, 4.13, 5.00, 5.50,  5.68, 9.00, 9.51,
     &  9.68, 10.00, 15.00, 19.00, 20.00, 22.00,  22.10, 23.00,  
     & 24.80, 25.00, 30.00, 30.50, 36.00, 38.00, 40.00, 40.10,  
     & 46.00, 50.00, 50.10, 55.00, 60.00, 69.80, 76.00, 100.00,  
     & 107.00, 110.00, 140.00, 143.00, 176.00, 
     &  199.00, 200.00, 230.00,  
     & 240.00, 286.00, 300.00, 310.00, 400.00,
     &  470.00, 480.00, 600.00,  
     & 610.00, 660.00, 760.00, 820.00, 
     &  900.00, 1000.00, 1090.00, 1160.00,  
     & 1450.00, 1600.00, 1800.00, 1920.00, 
     &  1950.00, 2030.00, 2410.00, 2500.00,  
     & 3000.00, 3030.00, 3040.00, 3500.00, 
     &  3600.00, 3820.00, 4000.00, 4500.00,  
     & 4800.00, 5000.00, 5500.00, 5780.00,
     &  6000.00, 6050.00, 7000.00, 7060.00,  
     & 8000.00, 9000.00, 9030.00, 9600.00,
     &  10000.00, 10500.00, 12000.00, 12100.00,  
     & 13200.00, 14000.00, 15200.00, 15600.00,
     &  17000.00, 19200.00, 20000.00, 24000.00, 
     & 24100.00, 28000.00, 30000.00, 30400.00,
     &  32000.00, 37000.00, 38000.00, 38200.00, 
     & 40000.00, 42300.00, 45000.00, 48000.00,
     &  49000.00, 50000.00, 55000.00, 58900.00,  
     &  60000.00, 60500.00, 70000.00, 76200.00,  
     & 77800.00, 80000.00, 87100.00, 90000.00, 
     &  97700.00, 100000.00, 
     & 109000.00, 118000.00, 119000.00, 125000.00,  
     & 129000.00, 150000.00, 175000.00, 200000.00,  
     & 225000.00, 250000.00, 400000.00, 500000.00,  
     & 600000.00/
      data(ex(2,i),i=1,15) /1000.00, 6000.00,12500.00,  
     & 18000.00, 25000.00, 37500.00, 50000.00,   
     & 62500.00, 75000.00, 100000.00, 125000.00,  
     & 150000.00, 250000.00, 350000.00, 500000.00/  
      data(ex(3,i),i=1,149) /50.00, 65.00, 80.00,   
     & 100.00, 150.00, 200.00, 230.00, 
     & 250.00, 300.00, 390.00, 400.00,  
     & 500.00, 890.00, 1000.00, 1530.00, 
     & 1600.00, 2000.00, 2130.00, 2300.00,
     & 2700.00, 2720.00, 2800.00, 3000.00, 
     & 3050.00, 3380.00, 3500.00, 4000.00, 4040.00,  
     & 4400.00, 4900.00, 5000.00, 5780.00, 5800.00,
     & 6000.00, 6200.00, 6400.00, 6890.00, 
     & 7000.00, 7360.00, 7700.00, 7790.00, 
     & 8000.00, 8600.00, 8960.00, 9000.00, 
     & 9100.00, 9200.00, 9300.00, 9740.00,  
     & 10000.00, 10300.00, 11300.00, 12000.00, 
     & 12300.00, 12600.00, 13000.00, 13500.00,
     & 14000.00, 14700.00, 15000.00, 16000.00, 
     & 16600.00, 16700.00, 17000.00, 18000.00,
     & 19000.00, 19700.00, 20000.00, 20800.00, 
     & 23000.00, 23200.00, 25000.00, 25400.00, 
     & 26000.00, 26600.00, 27000.00, 28000.00,
     &  29700.00, 30000.00, 30400.00, 32000.00, 33000.00,
     & 34000.00, 34700.00, 35000.00, 36000.00,  
     & 36700.00, 38700.00, 40000.00, 42600.00, 
     & 44500.00, 45000.00, 46000.00, 46900.00,   
     & 48300.00, 48700.00, 50000.00, 51300.00,  
     & 53900.00, 55000.00, 56900.00, 60000.00,  
     & 62500.00, 70000.00, 72900.00, 80000.00,  
     & 81100.00, 88800.00, 90000.00, 100000.00, 
     & 103000.00, 109000.00, 110000.00, 120000.00,  
     & 128000.00, 130000.00, 140000.00, 149000.00,  
     & 150000.00, 160000.00, 180000.00, 200000.00,  
     & 250000.00, 300000.00, 350000.00, 400000.00,  
     & 425000.00, 440000.00, 500000.00, 550000.00,  
     & 600000.00, 654000.00, 700000.00, 800000.00, 851000.00, 
     & 880000.00, 900000.00, 1000000.00, 1040000.00,  
     & 1060000.00,1250000.00, 1430000.00, 1500000.00, 
     & 1750000.00,2000000.00, 2050000.00, 2450000.00,  
     & 2560000.00, 3280000.00/  
      data(ex(4,i),i=1,41) /1500.00,  2000.00, 3000.00, 
     & 3500.00, 4000.00, 5000.00, 5500.00, 5900.00, 
     & 6000.00, 7000.00, 8000.00, 8100.00, 9000.00,
     & 10000.00, 11500.00, 12000.00, 14000.00, 15000.00, 
     & 16000.00, 17500.00, 18000.0, 20000.00, 22000.00,
     & 24000.00, 25000.00, 26000.00, 30000.00, 33500.00,
     & 37000.00, 40000.00, 43000.00, 45000.00, 50000.00, 
     & 55000.00, 60000.00, 65000.00, 70000.00, 80000.00, 
     & 90000.00, 100000.00, 120000.00/ 
      data (ex(5,i),i=1,41) /750.00, 1000.00, 1100.00,   
     & 1250.00, 1340.00, 1500.00, 1700.00, 1750.00,   
     & 1850.00, 2000.00, 2250.00, 2300.00, 2410.00,  
     & 2500.00, 2750.00, 3000.00, 4000.00, 5000.00,  
     & 6000.00, 7000.00, 8000.00, 9000.00, 10000.00,  
     & 11000.00, 12000.00, 13000.00, 14000.00, 15000.00,  
     & 16000.00, 17000.00, 18000.00, 19000.00, 20000.00,  
     & 22000.00, 23000.00, 24000.00, 25000.00, 26000.00, 
     & 27000.00, 29000.00, 30000.00/  
      data (sx(1,i),i=1,140) /59.4, 48.7, 46.2,
     &  44.2, 40.9, 56.0, 37.1, 35.0,
     &  53.0, 55.0, 36.7, 47.0, 39.3,
     &  37.0, 50.0, 46.0, 43.0, 48.0,
     &  36.0, 32.8, 40.0, 29.1, 44.0,
     &  38.5, 33.4, 36.0, 32.0, 41.0,
     &  30.0, 33.0, 39.0, 28.3, 35.0,
     &  32.0, 26.9, 27.0, 32.5, 31.0,
     &  23.7, 24.8, 26.0, 19.0, 22.7,
     &  28.0, 25.0, 24.9, 21.2, 22.4,
     & 30.4, 27.1, 22.0, 25.4, 22.9,
     & 21.0, 21.0, 22.9, 21.0, 21.5,
     & 19.5, 20.3, 20.7, 20.0, 18.8,
     & 17.0, 13.3, 17.0, 17.3, 13.6,
     & 15.9, 12.0, 15.0, 12.1, 12.0,
     & 14.25, 11.11, 13.15, 12.50, 10.50,
     & 11.50, 10.00, 12.70, 10.50, 9.85,
     & 10.40, 12.20, 9.50, 8.30, 10.40,
     &  8.60, 9.20, 9.00, 7.90, 7.50,
     &  8.00, 6.90, 6.50, 5.60, 5.90,
     & 5.00, 5.70, 4.70, 4.10, 3.40,
     & 4.10, 2.97, 2.70, 2.10, 2.30,
     & 1.86, 2.20, 1.58, 1.35, 1.10,
     & 0.59, 1.03, 0.73, 0.44, 0.64,
     & 0.60, 0.33, 0.30, 0.23, 0.23,
     & 0.165, 0.163, 0.108, 0.1215, 0.078,
     & 0.064, 0.069, 0.050, 0.035, 0.021,
     & 0.0132, 0.0078, 0.0040, 0.00181, 0.00016,
     & 0.0000545, 0.0000187/
      data (sx(2,i),i=1,15) /0.88, 2.30,
     & 2.55, 2.50, 2.20, 1.32, 0.80,
     & 0.40, 0.21, 0.11, 0.0572, 0.0334,
     & 0.0059, 0.0016, 0.0003/
      data (sx(3,i),i=1,149) /0.22600000, 0.45200000, 0.39600000,
     &  0.43300000, 0.49800000, 0.56300000, 1.10000000, 0.64800000,
     &  0.76100000, 2.30000000, 1.18000000, 1.65000000, 4.90000000,
     & 4.67000000, 6.74000000, 7.10000000, 7.54750000, 7.24000000,
     & 8.80000000, 7.90000000, 7.48000000, 7.65500000, 8.40000000,
     & 9.14000000, 7.70000000, 7.97000000, 8.68667000, 7.88000000,
     &  8.15000000, 7.92000000, 8.54667000, 8.10000000, 8.27000000,
     &  8.91000000, 8.40000000, 8.70000000, 7.93000000, 8.49000000,
     & 10.0000000, 8.02000000, 7.94000000, 8.50500000, 7.89000000,
     & 7.86000000, 8.55000000, 8.10000000, 8.00000000, 9.15000000,
     & 9.69000000, 8.55857000, 7.62000000, 7.44000000, 7.50000000,
     & 8.10000000, 7.10000000, 9.66000000, 6.84000000, 6.55000000,
     & 6.24000000, 7.43143000, 6.79000000, 6.50000000, 6.19000000,
     & 9.41000000, 5.94000000, 5.60000000, 5.46000000, 5.92000000,
     & 5.20000000, 4.80000000, 4.78000000, 5.27286000, 4.50000000,
     & 4.62000000, 4.18000000, 3.90000000, 4.20000000, 3.60000000,
     & 4.08556000, 3.62000000, 2.70000000, 3.10000000, 3.04000000,
     & 3.16000000, 3.66200000, 3.09000000, 2.60000000, 2.78000000,
     & 2.87000000, 2.36000000, 1.83000000, 1.89333000, 2.14000000,
     & 2.09000000, 1.53000000, 1.41000000, 1.54862000, 1.71000000,
     & 1.15000000, 1.11550000, 1.47000000, 1.03400000, 0.75200000,
     & 0.68433000, 0.55000000, 0.56500000, 0.35400000, 0.27200000,
     & 0.47000000, 0.27425000, 0.16800000, 0.14700000, 0.17000000,
     & 0.11600000, 0.08200000, 0.11000000, 0.06730000, 0.04580000,
     & 0.0692000,0 0.04180000, 0.02740000, 0.01745000, 0.00365000,
     & 0.00242000, 0.00097200, 0.00040300, 0.00028200, 0.00038000,
     & 0.00015700, 0.00011100, 0.00006570, 0.00004000, 0.00002950,
     & 0.00001730, 0.00001100, 0.00000920, 0.00000732, 0.00000428,
     & 0.00000580, 0.00000350, 0.00000110, 0.00000075, 0.00000035,
     &  0.00000013, 0.00000002, 0.00000017, 0.00000006, 0.00000004,
     & 0.00000002/
      data (sx(4,i),i=1,41) /0.02100, 0.02700, 0.01600,
     & 0.01350, 0.03400, 0.05367, 0.05000, 0.14800,
     & 0.06700, 0.10567, 0.16450, 0.16300, 0.18367,
     & 0.22140, 0.24000, 0.27133, 0.32100, 0.28500,
     & 0.33200, 0.34000, 0.32700, 0.35400, 0.30267,
     & 0.29700, 0.33650, 0.31000, 0.32000, 0.27500,
     & 0.23000, 0.20000, 0.08500, 0.17500, 0.09550,
     & 0.09700, 0.05950, 0.05900, 0.03550, 0.02400,
     & 0.01720, 0.01100, 0.00540/
      data (sx(5,i),i=1,41) /0.15400, 0.17400,
     & 0.23400, 0.28500, 0.29500, 0.29300, 0.32200,
     & 0.31200, 0.31800, 0.29100, 0.27700, 0.28000,
     &  0.27000, 0.25000, 0.23400, 0.25725, 0.28950,
     &  0.31100, 0.27267, 0.29400, 0.28550, 0.29200,
     & 0.31800, 0.30150, 0.27650, 0.27100, 0.26450,
     & 0.24500, 0.21500, 0.19600, 0.1865,0 0.17600,
     & 0.21233, 0.14400, 0.13250, 0.1100,0 0.12200,
     & 0.10000, 0.14800, 0.10100, 0.1400/
c
c Table of data pot (ionization potential) and aneff (effective number)
c in Lotz (1966)
c
      data pot /54.4, 28.2, 32.6, 37.6, 45.2, 
     & 0.0, 61.3, 52.5, 58.8, 68.9,
     & 0.0, 102.0, 100.0, 83.7, 101.0,
     & 0.0, 392.0, 155.0, 148.0, 132.0,
     & 0.0, 490.0, 552.0, 219.0, 167.0,  
     & 0.0, 0.0, 667.0, 739.0, 270.0, 
     & 0.0, 0.0, 0.0, 871.0, 380.0,
     & 0.0, 0.0, 0.0, 0.0, 1196.0,
     & 0.0, 0.0, 0.0, 0.0, 1362.0/   
      data aneff /1.1, 3.0, 4.1, 5.0, 5.5,
     & 0.0, 2.9, 3.0, 4.0, 6.0,
     & 0.0, 2.2, 3.0, 3.0, 5.0,
     & 0.0, 2.0, 2.2, 3.1, 4.0,
     & 0.0, 1.0, 2.0, 2.2, 3.0,
     & 0.0, 0.0, 1.0, 2.0, 3.1,
     & 0.0, 0.0, 0.0, 1.0, 2.3,
     & 0.0, 0.0, 0.0, 0.0, 2.0,
     & 0.0, 0.0, 0.0, 0.0, 1.0/  
c
c Table of a, b, c and d (parameters to compute the electron production
c cross section for collision with a proton (Rudd 1985)).
c
      data abcd /0.28, 1.15, 0.44,0.907,
     & 0.49, 0.62, 0.13,1.520,
     & 1.63, 0.73, 0.31, 1.140,
     & 3.85, 1.98, 1.89, 0.890,
     & 5.67, 5.50, 2.42, 0.650,
     & 7.33, 11.10, 4.12, 0.410,
     & 0.71, 1.63, 0.51, 1.240,
     & 3.82, 2.78, 1.80, 0.700,
     & 4.77, 0.00, 1.76, 0.930,
     & 3.67, 2.79, 2.08, 1.050,
     & 6.55, 0.00, 3.74, 1.160,
     & 4.01, 0.00, 1.73, 1.020,
     & 4.55, 2.07, 2.54, 1.080/
c
c
      ifail=0
      do k=1,nsp
        ik=ksp(k)
        if(ik.lt.1.or.ik.gt.30)then
          ifail=-63701
          return
        endif
        l=k1to1(ik)
        if(l.le.0)then
          ifail=-63702
          return
        endif
c
        if( ik.ge.19 )then
          flambd = 1.3890D9 * eng * debye * uma(ik)/(1.6604D-24+uma(ik))
          if( flambd .le.0.0d0 )then
            ifail=-63703
            flambd = 0.0d0
          else
            flambd = dlog(flambd)
          endif
        else  
          flambd = 0.0d0
        endif
c
        na = 0
        nb = 0
        nc = 0
        e = eng
c
        if (l .lt. 14) then
c
c         Electron production cross sections (Rudd, 1985)
c         H=1, He=2, Ne=3, Ar=4, Kr=5, Xe=6, H2=7, N2=8, O2=9, CO=10, 
c         CO2=11, NH3=12, CH4=13, ion=14, O=15, N=16 
c
          xx = e / 24969.6D0
          yn = 1.0D16 * 4.0D0*pi*0.529D-8**2 / (1.0D0/abcd(3,l)/
     :       xx**abcd(4,l)+xx/(abcd(1,l)*dlog(1.0D0+xx)+abcd(2,l)))
c
c         Charge exchange cross sections.
c         Cross sections of ph1s, ph2s, ph2p and ph2 from Tawara (1985).
c         For the data of N2, McDaniel (1964) and for He, Brandsden 
c         (1954).
c
          if (l .eq. 1) then
            km = 1
            nm = 140
          else if (l .eq. 2) then
            km = 2
            nm = 15
          else if (l .eq. 7) then
            km = 3
            nm = 149
          else
            goto 100
          end if
          do i=1,nm
            x(i) = ex(km,i)
            y(i) = sx(km,i)
          end do
          do i=1,nm
            if (x(i) .le. e) na = na + 1
          end do
c         
          if (e .le. x(1)) then
            na = 1
          else if (e .gt. x(nm)) then
            na = nm - 1
          end if
          yen = y(na) + (y(na+1)-y(na)) * (e-x(na)) / (x(na+1)-x(na))
          if (yen .lt. 0.0D0) yen = 0.0D0
c     
          if (l .eq. 1) then
            do i=1,41
              x1(i) = ex(4,i)
              x2(i) = ex(5,i)
              y1(i) = sx(4,i)
              y2(i) = sx(5,i)	
            end do
            do j=1,41
              if (x1(j) .le. e) nb = nb + 1
              if (x2(j) .le. e) nc = nc + 1
            end do
            if ((e .le. x1(1)) .or. (e .ge. x1(41))) then
              yen1 = 0.0D0
            else
              yen1 = y1(nb) + (y1(nb+1)-y1(nb)) * (e-x1(nb)) / 
     :             (x1(nb+1)-x1(nb))
            end if
            if ((e .le. x2(1)) .or. (e .ge. x2(41))) then
              yen2 = 0.0D0
            else
              yen2 = y2(nc) + (y2(nc+1)-y2(nc)) * (e-x2(nc)) / 
     :             (x2(nc+1)-x2(nc))
            end if
            yen = yen + yen1 + yen2
          end if
          yn = yn + yen
        else if (l .eq. 15) then
          if (e .le. 10000.0D0) then
            yn = 2.0D0 * (3.47D0-0.20D0*dlog10(e))**2
          else if (e .ge. 500000.0D0) then
            yn = 400.0D0 * 5.7D0 * dlog(e/1836.0D0/15.7D0) / 
     :          (15.7D0*e/1836.0D0)
          else
            yn = 10.0D0**(1.15405248793D0-0.971758372626D0*(e-1.0D4)/
     :          4.9D5)
          end if
        else if (l .eq. 16) then
          if (e .ge. 100000.0D0) then
            yn = 4.9D0 * 400.0D0 * dlog(e/1836.0D0/16.2D0)/
     :         (16.2D0*e/1836.0D0)
          else if (e .ge. 5000.0D0) then
            yn = 2.14230630103D0 + 10.0D0**(-4.018D0+4.8989D0*
     :         dexp(-e/377683.0D0))
          else
            yn = 0.0D0
          end if
        else if (l .gt. 16) then
c
c         Collision cross section between a proton and an ion 
c         depending on its charge and its atomic mass.
c         High energy: Lotz (1969).
c         Low energy: Banks and Kockarts (1973).
c         He+=1, C+=2, N+=3, O+=4, Ne+=5, N2+=6, O2+=7, e-=8, 
c         H+=9 (Coulomb only)
c
          im = l - 16
          if ((im .eq. 6) .or. (im .eq. 7)) then
            lm = im - 3
          else
            lm = im
          end if
          nz = 1
          if ((lm .le. 5) .and. (nz .le. ns(lm))) then
            ssi = aneff(lm,nz) * 400.D0 * dlog(8.D5/1836.D0/pot(lm,nz))
     :          / (pot(lm,nz)*8.0D5/1836.D0)
            if (e .ge. 8.0D5) then
              sio = aneff(lm,nz) * 400.0D0 * dlog(e/1836.0D0/pot(lm,nz))
     :            / (pot(lm,nz)*e/1836.0D0)
            else if (e .gt. 2.0D5) then
              sio = ssi * (2.0D0-(e-2.0D5)/6.0D5)
            else
              sio = 2.0D0 * ssi
            end if
            if (im .ge. 6) sio = 2.0D0 * sio
            seff = 3.25D2 * ((1.0D0+pm(im))/pm(im)*nz/e)**2 * flambd
            yn = sio + seff
          else if (im .eq. 8) then
            yn = 3.25D2 * ((1.0D0+5.461D-4)/5.461D-4/e)**2 * flambd
          else if (im .eq. 9) then
            yn = 3.25D2 * 4.0D0 / e**2 * flambd
            if (yn .lt. 1.0D-11) yn = 1.0D-11
          else 
            yn = 0.0D0
          end if
        else
          yn = 0.0D0
        end if
c
  100   continue
        cs(ik) = yn * 1.0D11
c
      enddo
c
      end
C----------------------------------------------------------------------
      subroutine ua638
     :                (nsp, ksp, cs, ifail)
c
c     hassitt cross section
c
      include 'structure.h'
c
      integer*4 nsp, ksp(30), ifail
      real*8    cs(30)
c
      integer*4 k1to1(30)
      real*8 csh(12)
      integer*4 k,ik,l
c
      data k1to1/  7,  6, -1, -1, -1, -1, -1,  3,  2, -1,
     :            -1, -1, -1, -1, -1,  4,  5, -1,  1, 12,
     :            11, -1, -1, 10, -1, -1,  9,  8, -1, -1/    
c
      DATA CSH /25.0D0, 179.0D0, 159.0D0, 89.5D0, 80.0D0,  25.2D0,
     :          13.1D0, 168.0D0, 158.0D0, 78.5D0, 12.6D0,   0.0D0/
c  
c
      ifail=0
      do k=1,nsp
        ik=ksp(k)
        if(ik.lt.1.or.ik.gt.30)then
          ifail=-63801
          return
        endif
        l=k1to1(ik)
        if(l.le.0)then
          ifail=-63802
          return
        endif
c
        cs(ik)=csh(l) * 1.0d+10
c               en mbarn=10-27 cm2
      enddo
c
      end
C----------------------------------------------------------------------
      subroutine ua639
     :                (FL, kDAYNR, SSN, dca)
C
C EQUATORIAL ELECTRON DENSITY AT FL (L) FOR SSN = 13 MONTH AVERAGE
C SUNSPOT NUMBER.
C TAKEN FRON CARPENTER AND ANDERSON (JGR 1992,97,1097 & PREPRINT).
C
      IMPLICIT NONE

      REAL*8  FL, AD, ALNE
      REAL*8  SSN, dca
      INTEGER*4 kDAYNR
c
      AD = 0.01721D0 * (kDAYNR+9.0D0)
      ALNE = -0.3145D0 * FL + 3.9043D0 + (0.15D0*(COS(AD)-
     :       0.5D0*DCOS(2.0D0*AD))+0.00127D0*SSN-0.0635D0) /
     :       DEXP((FL-2.0D0)/1.5D0)
      DCA = 10.0D0**ALNE
c
      end
C----------------------------------------------------------------------