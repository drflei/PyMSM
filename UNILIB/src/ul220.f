# 1 "ul220.f"
      SUBROUTINE UL220
     :          (mpos, alpha, nfbm,  
     :           fbm, flm, fkm, fsm, fbeq, fs, ifail)
C
C!    Get information on a magnetic field line segment
C_    er06
C
      INCLUDE 'structure.h'
cDEC$ IF DEFINED (_x86_)
cDEC$ ATTRIBUTES DLLEXPORT :: UL220
cDEC$ ENDIF
C

      EXTERNAL UM530, UF420, UL230, UL240
C      
C     INTERFACE
C
        TYPE(zgeo) :: mpos
        INTEGER*4     ifail, nfbm
        REAL*8        alpha(nfbm), fbm(nfbm), flm(nfbm)
        REAL*8        fkm(nfbm), fsm(nfbm), fbeq, fs(nfbm)
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
      COMMON /UC190
     :               /prop, stepx, stpmin, umsq, upsq, uk2, 
     :                uk3, epskm, epsrel, stplst, xclat, kmflg, nxstp
        REAL*8        prop, stepx, stpmin, umsq, upsq, uk2, uk3, 
     :                epskm, epsrel, stplst, xclat
        INTEGER*4     kmflg, nxstp              
      COMMON /UC192/  xrmin, xbmin, xtmin, xbmax, epslon, epsfl,
     :                fvet, pvet, epsomeg, dltlat
        REAL*8        xrmin, xbmin, xtmin, xbmax, epslon, epsfl
        REAL*8        fvet, pvet, epsomeg, dltlat
C
C     VARIABLES
C
        integer*4     kl426
        REAL*8        erflg, fim, falt
        INTEGER*4     k1, k2, k, i, kjn, kjs, ks0
        TYPE(zvec) :: mb
        TYPE(zlbl) :: mlab2
        TYPE(zfln) :: mlin2
      TYPE(zgeo) :: mdummy
C
C     CODE
C
c      print*,'(ul220) mpos',mpos%radius,mpos%colat,mpos%elong
c      print*,'(ul220) alpha',nfbm,alpha
      ifail = 0
      falt  = -999.0d0
C      
C  1/ Evaluate the local value of the magnetic field intensity
C     --------------------------------------------------------
C
      fbeq  = 0.0d0
      k1    = 1
      k2    = nfbm
      CALL UM530(mpos, mb, ifail)
      if( ifail.lt.0 )then
        erflg   = -10.0d0
        if( ifail .eq. -53003 ) erflg = -1.0d0
        if( ifail .eq. -53102 ) erflg = -1.0d0
        if( ifail .eq. -53301 ) erflg = -2.0d0              
        if( ifail .eq. -53002 ) erflg = -3.0d0
        if( ifail .eq. -53004 ) erflg = -5.0d0
        if( ifail .eq. -53001 ) erflg = -5.0d0            
        goto 101
      endif 
c      print*,'(ul220) mb',mb.dnrm
C      
C  2/ Deduce the magnetic field intensity at the mirror points
C     --------------------------------------------------------
C
      if( alpha(1) .le. 0.0d0 )then
        ifail = -22001
        goto 100
      endif
C      
      fbm(1) = mb%dnrm / SIN(alpha(1)*deg)**2
      do i = 2, nfbm
        if( alpha(i) .le. alpha(i-1) .or. alpha(i) .gt. 90.0d0 )then
          ifail = -22002
          goto 100
        endif
        fbm(i) = mb%dnrm / SIN(alpha(i)*deg)**2
      enddo
C      
C  3/ Trace the field line
C     --------------------
C
  150 continue
      k = k2 - k1 + 1
c      print*,'(ul220) bm',nfbm,fbm
      kl426=0
      CALL UF420(mpos, fbm(k1), k, falt, ifail)
      if( ifail .gt. nx170 )then
        ifail = ifail - nx170
        kl426 = -42699
      endif
c      print*,'(ul220) ifail from uf420',ifail
      if( ifail .lt. 0 )then
C                     (1) field line reaches high geomagnetic latitude
C                     (2) field line crosses the magnetopause
C                     (3) integration is outside valid field limits
C                     (4) maximum step number has been reached
C                     (5) field line tracing goes deep inside the Earth
C                     (6) out of memory
C                     (7) the magnetic field has a local maximum
C                     (8) an error occurred in the interpolation routine
C                     (9) an error occurred in the evaluation of I
C                    (10) not documented
        erflg   = -10.0d0
        if( ifail .eq. -53003 ) erflg = -1.0d0
        if( ifail .eq. -53102 ) erflg = -1.0d0     
        if( ifail .eq. -53301 ) erflg = -2.0d0              
        if( ifail .eq. -53002 ) erflg = -3.0d0
        if( ifail .eq. -42009 ) erflg = -4.0d0
        if( ifail .eq. -53004 ) erflg = -5.0d0
        if( ifail .eq. -53001 ) erflg = -5.0d0            
        if( ifail .eq. -42304 ) erflg = -5.0d0            
        if( ifail .eq. -42301 ) erflg = -6.0d0
        if( ifail .eq. -42401 ) erflg = -6.0d0
        if( ifail .eq. -42305)then
          erflg = -3.0d0
          if( sin(deg*mele(kgp)%beg%coord%colat)
     :        .le. xtmin ) erflg = -1.0D0
        endif
        if( ifail/100 .ge. -428 .and. ifail/100 .le. -426)
     :          erflg = -8.0d0
        flm(k1) = erflg
        fkm(k1) = erflg
        fsm(k1) = erflg
        fs(k1)  = erflg
        k1      = k1 + 1
        if( k1 .le. k2 )goto 150
        if( erflg.LE.-8.0D0 ) then
          if (.false.) print*,'UL220',erflg,ifail
        endif
        ifail   = -22004
        return            
      end if
      ks0 = kgp
C
      fbeq = mlin%equat%b%dnrm 
      if( kl426 .lt. 0 ) fbeq=-8.0D0
C      
C  4/ Compute the set of informations
C     -------------------------------
C
      mlab2 = mlab
      mlin2 = mlin
      kjn = mlin%ind%jmirpn
      kjs = mlin%ind%jmirps      
  200 continue
C  
C  4.1/ for the first pair of mirror points
C        
        fsm(k1) = mele(kjs)%arcl - mele(kjn)%arcl
        fs(k1)  = mele(ks0)%arcl - mele(kjn)%arcl
        CALL UL230(fim,ifail)
        if( ifail .lt. 0 )then
          erflg    = -9.0d0
          flm(k1)  = erflg
          fkm(k1)  = erflg
          if (.false.) PRINT*,'UL220: UL230',IFAIL,' --> L=-9'
        else
          CALL UL240(mlab, ifail)
          flm(k1)   = mlab%flmi
          fkm(k1)   = mlab%fkauf
          if( ifail .lt. 0 )then
            erflg    = -9.0d0
            flm(k1)  = erflg
            fkm(k1)  = erflg
            if (.false.) PRINT*,'UL220: UL240',IFAIL,' --> L=-9'
          endif
        endif
C        
C  4.2/ go to the next pair of mirror points
C
      k1 = k1 + 1
      if( k1.le.k2 )then
        if( fbm(k1) .lt. fbeq  .and. .false.)
     :        print*,'UL220: A problem for Michel !'
  210   continue
        kjn = kjn + 1
        if( mele(kjn)%beg%b%dnrm .gt. fbm(k1) )goto 210
        if( ABS(mele(kjn)%beg%b%dnrm-fbm(k1)) .gt.
     :      ABS(mele(kjn-1)%beg%b%dnrm-fbm(k1))  ) kjn=kjn-1
  220   continue
        kjs = kjs - 1
        if( mele(kjs)%beg%b%dnrm .gt. fbm(k1) )goto 220
        if( ABS(mele(kjs)%beg%b%dnrm-fbm(k1)) .gt.
     :      ABS(mele(kjs+1)%beg%b%dnrm-fbm(k1))  ) kjs=kjs+1
        if( kjn .gt. kjs .and. .false.)
     :          print*,'UL220: An other problem for Michel !'
        if(      abs(mele(kjn)%beg%b%dnrm-fbm(k1)) .le. fbm(k1)*epsrel 
     :     .and. abs(mele(kjs)%beg%b%dnrm-fbm(k1)) .le. fbm(k1)*epsrel )
     :  then
          mlin%ind%jmirpn  = kjn
          mlin%ind%jmirps  = kjs
          mlab%fbmp        = fbm(k1)    
          goto 200
        endif
        ifail = -22003
      else
        ifail = 0
      endif
C
      mlab    = mlab2
      mlin    = mlin2               
C
C  5/ Errors
C     ------
C
      mdummy = mele(ks0)%beg%coord
      if(      mpos%radius .ne. mdummy%radius 
     :    .or. mpos%colat  .ne. mdummy%colat
     :    .or. mpos%elong  .ne. mdummy%elong )
     :    ifail = -22005
C
  100 continue
C
      erflg = -10.0d0
  101 continue
      do i = k1, k2
        flm(i) = erflg
        fkm(i) = erflg
        fsm(i) = erflg
        fs(i)  = erflg             
      enddo
C
C
      END
c----------------------------------------------------------------------
      subroutine ul225( fbm, flm, nfbm, frd, fla, ifail )
      INCLUDE 'structure.h'
      INTEGER*4     ifail,nfbm
      REAL*8        fbm(nfbm), flm(nfbm), frd(nfbm), fla(nfbm)
c
      COMMON /UC160/pi, deg, re, gmagmo, eclipt, geoid, uma
      REAL*8        pi, deg, re, gmagmo, eclipt, geoid(3), uma(30)
      COMMON /UC190/prop, stepx, stpmin, umsq, upsq, uk2, uk3, epskm,
     :              epsrel, stplst, xclat, kmflg, nxstp       
      REAL*8        prop, stepx, stpmin, umsq, upsq, uk2, uk3, epskm,
     :                epsrel, stplst, xclat
      INTEGER*4     kmflg, nxstp       
      COMMON /UC140/mint, mext, msun
      TYPE(zimf) :: mint
      TYPE(zsun) :: msun
      TYPE(zemf) :: mext          
c
      REAL*8        GM, AD, R1, rinv
      integer*4     i
C
      if( mod(kmflg,10) .eq. 0 )then
        gm          = gmagmo
        ifail       = 0
      else
        gm          = mint%gmmo
        ifail       = 1
      endif
c
      do i = 1, nfbm
        if( fbm(i).le.0 .or. flm(i).le.0 )then
          frd(i)=-99.0d0
          fla(i)=-99.0D0
          ifail = -22501
        else
          AD    = (fbm(i)/GM)**2 * flm(i)**6
          R1    = 0.0d0
          RINV  = 0.5d0
          DO WHILE (ABS((R1-RINV)/RINV) .GT. 0.00005d0)
             R1 = RINV
             RINV = R1 - (AD*R1**6+3*R1-4) /(6*AD*R1**5+3)
          END DO
          IF (RINV .LE. 1.0d0) THEN
             fla(i) = ACOS(SQRT(RINV)) / deg
             frd(i) = rinv * flm(i)
          elseIF (RINV .LE. 1.0001d0) THEN
             fla(i) = 0.0D0
             frd(i) = flm(i)
          else
             frd(i)=-99.0d0
             fla(i)=-99.0D0
             ifail = -22501
          endif
        endif
      enddo
c
      end
c----------------------------------------------------------------------
      subroutine ul229( altmin, ifail )
      INCLUDE 'structure.h'
      REAL*8        altmin
      INTEGER*4     ifail
c
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
      integer*4       kn,ks
      TYPE(zgeo) :: mpn, mps
C
      ifail =0
      altmin=-999.0d0
      kn=mlin%ind%jmirpn
      ks=mlin%ind%jmirps
      if( kn.le.0 .or. kn.gt.nsg .or. ks.le.0 .or. ks.gt.nsg )then
        ifail=-22901
        return
      endif
      mpn=mele(kn)%beg%coord
      mps=mele(ks)%beg%coord
      altmin=min(mpn%radius,mps%radius)-re

      end
C----------------------------------------------------------------------