# 1 "uf420.f"
      SUBROUTINE UF420
     :          (mpos, fbm, nfbm, falt, ifail)
C
C!    Trace a magnetic field line segment passing through a given 
C!    position
C_    er08
C
      INCLUDE 'structure.h'
cDEC$ IF DEFINED (_x86_)      
cDEC$ ATTRIBUTES DLLEXPORT :: UF420
cDEC$ ENDIF
C      
      EXTERNAL UF421, UF422, UF424, UF429
C      
C     INTERFACE
C
        INTEGER*4     nfbm, ifail
        REAL*8        fbm(nfbm), falt
        TYPE(zgeo) :: mpos
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
     :               /prop, stepx, stpmin, umsq, upsq, uk2,
     :                uk3, epskm, epsrel, stplst,xclat, kmflg, nxstpold 
C
        REAL*8        prop, stepx, stpmin, umsq, upsq, uk2, uk3, epskm,
     :                epsrel, stplst,xclat
        INTEGER*4     kmflg, nxstpold              
C
      COMMON /UC192/ xrmin, xbmin, xtmin, xbmax, epslon, epsfl,
     :               fvet, pvet, epsomeg, dltlat
C
      REAL*8         xrmin, xbmin, xtmin, xbmax, epslon, epsfl
      REAL*8         fvet, pvet, epsomeg, dltlat           
C
C     VARIABLES
C
        LOGICAL       lalt, lfbm, lfal, l426
        INTEGER*4     ji0, np1, np2, ji1, nm1, ksgn
        REAL*8        pb, palt, scale, xlen, bx
        REAL*8        xupe, xume
        INTEGER*4     kbeg, kend, kmpbeg, kmpend
        INTEGER*4     i, k, kd, ks0
        TYPE(zgeo) :: mpx, mbeg, mpbeg, mend, mpend
        INTEGER*4     nxstp
        PARAMETER ( nxstp=nx170-20 )
C
C *
C * fbm(-)     boundary conditions on the magnetic field intensity
C *            |    fbm(1) > fbm(2) >...
C * falt       boundary condition on the altitude     
C * lalt       no condition on altitude when false
C * ji0        starting index of the field line segment
C current point
C * pb         magnetic field intensity [G]
C * palt       altitude [km]
C * mpx        GEO position
C
C * scale      = +/-1, parameter used by uf422
C
C Doc on UF422:
C       #     ifail =  2&3  magn. field. intensity at the start point
C       #                     is equal to val_B 
C       #              1&3  altitude of the start point is equal to falt
C       #           +=  10  pass through an extremum in the magn. field
C       #                     intensity
C       #               30  stop on an maximum in the magn. field
C       #                     intensity
C       #               50  stop on an minimum in the magn. field
C       #                     intensity
C       #           += 100  val_B has been reached
C       #              200  falt has been reached
C       #              300  falt and val_B have been reached
C
C       #           -42201  invalid SCALE
C       #           -42202  find a minimum in B and none condition will be found
C       #           -42205  starting point is an extremum not compatible with KDIR
C       #           -42206  caught in an infinite loop
C       #
C
C     CODE
C
c
c initialisation
c
      l426                 = .FALSE.
      mlab%linv            = .FALSE.
      mlab%llmi            = .FALSE.
      mlab%lkauf           = .FALSE.
      mlab%fbmp            = fbm(1)
      mlab%lbmp            = .TRUE.      
      lalt                 = falt .ge. -500.0d0
      scale                = 1.0d0     
      mlin%ind%jbeg        = -1
      mlin%ind%jend        = -1
      mlin%ind%jmirpn      = -1
      mlin%ind%jmirps      = -1
      mlin%keast           = -1
      mlin%kwest           = -1
      ifail                =  0
c (indices of the foot points)
      kbeg                 = -1
      kend                 = -1
c (indices of the MP[1] )
      kmpbeg               = -1
      kmpend               = -1
c (precision on the magnetic field intensity )
      xupe             = 1.0d0 + 0.3d0*epsrel
      xume             = 1.0d0 - 0.3d0*epsrel
C
C set current point as starting point
C
      mpx                  = mpos
C      
      CALL UF421 (1, mpx, pb, palt, xlen, ji0, ifail)
      if( ifail .lt. 0 )return
C      
C check if the starting point is inside the boundary conditions
C
      if(    (  lalt      .and. ( palt .lt. falt ) 
     :                    .and. ( pb .ge. fbm(1)*xupe))
     : .or.  (  .not.lalt .and. ( pb .ge. fbm(1)*xupe)))then
C     
C if it is not the case, move in the direction of decreasing B until
C a boundary condition is reached
C
        CALL UF422 (-1, fbm(1), falt, scale, mpx, i, ifail)
        if( ifail .lt. 0 )then
          mlin%equat%coord = mpx
          if( ifail .eq. -42205 .or.
     :        ifail .eq. -42202 .or.
     :        ifail .eq. -42206       ) then
            ifail = -42008           
          endif
          return
        endif
c
c set the new point as starting point
c
        CALL UF421 (1, mpx, pb, palt, xlen, ji0, ifail)
        if( ifail .lt. 0 )return
        scale          = -scale
      endif
c
c determine the number of condition in the direction of
c increasing B
c
c * lfal = true when starting point is at the alt cond
c * np2  = 1 or 0.  When np2=1 there is an alt cond in the increasing dir
c * lfbm = true when starting point is at one of the B cond
C * np1  = number of MP in the increasing dir
c          when np1>0 = index in FBM of the first MP in the increasing dir
c * nm1  = index in FBM of the first MP in the decreasing dir
c * ksgn = 1 or -1. Step in the index of FBM for the next MP in the decreasing dir
c
c default value:
      lfal      = .false.
      np2       = 0
      lfbm      = .false.
      np1       = 0
      nm1       = -1
      ksgn      = 1
c
c count the number of boundary conditions on B
C
      if( pb .gt. fbm(1) * xupe )then
c       B > B1, B2, ...Bn   
        np1     = 0
        nm1     = 1
      elseif( pb .ge. fbm(1) * xume )then
c       B = B1
        lfbm    = .true.
        kmpbeg  = ji0
        mpbeg   = mpx
        if( nfbm .eq. 1 )then
          nm1   = 1
          ksgn  = -1
        else
          nm1   = 2
        endif
      elseif( pb .lt. fbm(nfbm) * xume )then
c       B1, B2, ...Bn > B
        np1     = nfbm
        nm1     = nfbm
        ksgn    = -1
      elseif( pb .le. fbm(nfbm) * xupe )then
c       B = Bn
        lfbm    = .true.
        np1     = nfbm-1
        nm1     = nfbm
        ksgn    = -1
      else
c       B1 > B > Bn
        do k = 2, nfbm-1
          if( fbm(k) .ge. fbm(k-1) )then
            ifail=-42001
            return                 
          elseif( pb .gt. fbm(k) * xupe )then
c           B(k-1) > B > Bk
            np1 = k-1
            nm1 = k
            goto 100	
          elseif( pb .ge. fbm(k) * xume )then
c           B = Bk
            lfbm= .true.
            np1 = k-1
            nm1 = k+1
            goto 100	
          endif
        enddo
c           B(n-1) > B > Bn
            np1 = nfbm-1
            nm1 = nfbm
  100   continue
        if( fbm(nfbm) .ge. fbm(nfbm-1) )then
          ifail=-42001
          return
        endif
      endif
c
c check the boundary conditions on alt
C
      if( lalt )then
        if( palt .gt. falt + epskm )then
          np2   = 1
        elseif( palt .ge. falt - epskm )then
          lfal  = .true.
          kbeg  = ji0
          mbeg  = mpx    
        endif
      endif
c
c check for a problem in the algorithm
c
      if( .not. lalt .and. .not. lfbm .and.
     :     np1 .le. 0 .and. np2 .le. 0      )then
          if (.false.) print*,'Programme aborted in UF420 (1)'
C         
C         Unilib 2.09 is stopping here, stop statement has
C         been replace by reseting UC170 and a negatif ifail
          nsg                  = -nsg
          mlab%linv            = .FALSE.
          mlab%llmi            = .FALSE.
          mlab%lkauf           = .FALSE.
          mlab%lbmp            = .FALSE.      
          mlab%LALP0           = .FALSE.      
          mlin%ind%jbeg        = -1
          mlin%ind%jend        = -1
          mlin%ind%jmirpn      = -1
          mlin%ind%jmirps      = -1
          mlin%keast           = -1
          mlin%kwest           = -1
          ifail = -42095
          return
C
      endif
c
c * kd = 1, 0 or -1. When kd=1, force to follow the field line in the increasing dir
c *                  When kd=-1, force to follow the field line in the decreasing dir
c *                  When kd=0, continue to the follow the field line in the same dir
      kd           = 1
      ji1          = ji0
c
c If needed, follow the m.f.l. in the increasing dir
c
      if( np1 .gt. 0 .or. np2 .gt. 0 )then
c
c REPEAT
  200 continue
c update the condition on B
        bx         = xbmax*2
        if( np1 .gt. 0 ) bx=fbm(np1)
c follow the m.f.l. until a condition is reached
        CALL UF422 (kd, bx, falt, scale, mpx, ji1, ifail)
        if( ifail .ge. 1000 )then
          l426=.true.
          ifail=ifail-1000
        endif              
c latter, if needed, continue in the same direction
        kd         = 0
        if( ifail .lt.   0 ) then         
          mlin%equat%coord = mpx
          return
        elseif( ji1 .ge. nxstp )then
          ifail    = -42009            
          return
        elseif( ifail .lt. 100 )then
          if (.false.) print*,'Programme aborted in UF420 (2)'
C         
C         Unilib 2.09 is stopping here, stop statement has
C         been replace by reseting UC170 and a negatif ifail
          nsg                  = -nsg
          mlab%linv            = .FALSE.
          mlab%llmi            = .FALSE.
          mlab%lkauf           = .FALSE.
          mlab%lbmp            = .FALSE.      
          mlab%LALP0           = .FALSE.      
          mlin%ind%jbeg        = -1
          mlin%ind%jend        = -1
          mlin%ind%jmirpn      = -1
          mlin%ind%jmirps      = -1
          mlin%keast           = -1
          mlin%kwest           = -1
          ifail = -42096
          return
C
        endif
C        
C  check if a boundary condition on the altitude has been reached
C
        if( ifail .ge. 200 )then
          if( np2 .gt. 0 .and. kbeg   .lt. 0 )then
            kbeg   = ji1
            mbeg   = mpx
            np2    = np2 - 1
          else
            ifail  = -42002
            return
          endif
c remove the flag concerning the altitude
          ifail    = ifail - 200
        endif
C        
C  check if a boundary condition on B has been reached
C
        if( ifail .ge. 100 )then
          if( np1 .gt. 1 )then
            np1    = np1 - 1
          else if( np1 .eq. 1 .and. kmpbeg .lt. 0 )then
            np1    = np1 - 1
            kmpbeg = ji1
            mpbeg  = mpx      
          else
            ifail  = -42003
            return
          endif
        endif
c
c UNTIL all the np1+np2 conditions are found
c
        if( np1+np2 .gt. 0 )GOTO 200
c
C Transpose the traced field line segment
c
        if( ji1 .gt. ji0 )then
          if( kbeg .gt. 0 )   kbeg   = ji1 + ji0 - kbeg 
          if( kmpbeg .gt. 0 ) kmpbeg = ji1 + ji0 - kmpbeg
          CALL UF429(ji0,ji1)
        else
          if (.false.) print*,'Programme aborted in UF420 (3)'
C         
C         Unilib 2.09 is stopping here, stop statement has
C         been replace by reseting UC170 and a negatif ifail
          nsg                  = -nsg
          mlab%linv            = .FALSE.
          mlab%llmi            = .FALSE.
          mlab%lkauf           = .FALSE.
          mlab%lbmp            = .FALSE.      
          mlab%LALP0           = .FALSE.      
          mlin%ind%jbeg        = -1
          mlin%ind%jend        = -1
          mlin%ind%jmirpn      = -1
          mlin%ind%jmirps      = -1
          mlin%keast           = -1
          mlin%kwest           = -1
          ifail = -42097
          return
C
        endif
      endif
c
c  Now, have to go back to the starting point and to follow the 
c  m.f.l. in the other direction
c
c  reverse the direction
      scale        = - scale
      kd           = - kd
c  now, kd=0 or kd=-1
      ks0          = ji1
C
c REPEAT
  300 continue
c
c check if all conditions have been reached
c
      if( (lalt .and. kend .lt. 0) .or. kmpend .lt. 0 )then
c update the condition on B
        bx         = xbmax*2
        if( nm1 .gt. 0 ) bx=fbm(nm1)
c follow the m.f.l. until a condition is reached
        CALL UF422 (kd, bx, falt, scale, mpx, ji1, ifail)
        if( ifail .ge. 1000 )then
          l426=.true.
          ifail=ifail-1000
        endif              
c latter, if needed, continue in the same direction
        kd         = 0
        if( ifail .eq. -42202 )then
c a minimum in B has been found before any condition
c NOTE: bx = fbm(nm1) is small than the magn field int. at the minimum
          mlin%equat%coord = mpx
          if( nm1 .eq. 1 )then
c           B1 is not part of the m.f.l.
            ifail  = -42010
            return
          elseif( ksgn .eq. 1)then
c           not all the Bk are part of the m.f.l.
            ksgn   = -1
            nm1    = nm1 + ksgn
            goto 300
          else
c           ksgn is equal to -1
c           there is an minimum in B which is greater than
c           the seeked value 
c            ->  uf422 missed the value or the magnetic field
c                is very strange
c           return with ifail = -42202
            return 
          endif
        elseif( ifail .lt.   0 ) then         
          mlin%equat%coord = mpx
          return
        elseif( ji1 .ge. nxstp )then
c         What is the problem ?
          ifail        = -42009
          return
        elseif( ifail .ge. 200 )then
C        
C boundary condition on the altitude has been reached
C
          ifail        = ifail - 200
          if( kbeg .lt. 0 )then
            kbeg       = ji1
            mbeg       = mpx       
          elseif( kend .lt. 0 )then
            kend       = ji1
            mend       = mpx              
          else
c           What is the problem ?
            ifail      = -42004
            return
          endif
          if( ifail.lt.100 )goto 300
        endif
        if( ifail .ge. 100 )then
C        
C boundary condition on the magnetic field intensity has been reached
C
          ifail        = ifail - 100
c
          if( nm1 .eq. 1 )then
            if(     ksgn .eq. 1    .and. 
     :              kmpbeg .eq. -1 .and. 
     :              ifail .lt. 50       )then
c first MP
              kmpbeg   = ji1
              mpbeg    = mpx
              if( nfbm .eq. 1 )then
                ksgn   = -1
              else
                nm1    =  2
              endif
              goto 300
            elseif( ksgn .eq. 1 .and.
     :              kmpbeg .eq. -1 .and. 
     :              ifail/10 .eq. 5     )then
c first and last MP
              kmpbeg   = ji1
              mpbeg    = mpx
              kmpend   = ji1
              mpend    = mpx
              nm1      = 0
              goto 300
            elseif( ksgn .eq. -1   .and.
     :              kmpbeg .gt. 0  .and. 
     :              kmpend .eq. -1 .and. 
     :              ifail  .lt. 50     )then
c last MP
              kmpend   = ji1
              mpend    = mpx
              nm1      = 0
              goto 300
c MK, June 15, 1998: changed kmpbeg.eq.1 into
            elseif( kmpbeg .gt.0 .and. 
     :              ifail .ge. 52     )then
c starting on the first and last MP
              if( ksgn .eq. -1 .and. nfbm .gt. 1 )then
c               What is the problem ?
                ifail    = -42012
                return
              endif
              kmpend   = ji1
              mpend    = mpx
              nm1      = 0
              goto 300
c MK, June 15, 1998: 4 lines added
            elseif( kmpbeg .eq. 1 .and. 
     :              ifail .ge. 52     )then
              if (.false.) print*,'UF420 to be checked (June 15, 1998)'
C             
C             Unilib 2.09 is stopping here, stop statement has
C             been replace by reseting UC170 and a negatif ifail
              nsg                  = -nsg
              mlab%linv            = .FALSE.
              mlab%llmi            = .FALSE.
              mlab%lkauf           = .FALSE.
              mlab%lbmp            = .FALSE.      
              mlab%LALP0           = .FALSE.      
              mlin%ind%jbeg        = -1
              mlin%ind%jend        = -1
              mlin%ind%jmirpn      = -1
              mlin%ind%jmirps      = -1
              mlin%keast           = -1
              mlin%kwest           = -1
              ifail = -42098
              return
C
            else
c             What is the problem ?
              ifail    = -42006
              return
            endif
          elseif( ifail/10 .eq. 5 .and. ksgn .eq. 1 )then
c MP at a minimum
            ksgn       = -1
            nm1        = nm1 + ksgn
            goto 300
          elseif( nm1 + ksgn .gt. nfbm )then
c MP at Bn
            ksgn       = -1
            goto 300
          else
c MP, continue
            nm1        = nm1 + ksgn
            goto 300
          endif
        endif
        if (.false.) print*,'Programme aborted in UF420 (4)'
C
C       Unilib 2.09 is stopping here, stop statement has
C       been replace by reseting UC170 and a negatif ifail
        nsg                  = -nsg
        mlab%linv            = .FALSE.
        mlab%llmi            = .FALSE.
        mlab%lkauf           = .FALSE.
        mlab%lbmp            = .FALSE.      
        mlab%LALP0           = .FALSE.      
        mlin%ind%jbeg        = -1
        mlin%ind%jend        = -1
        mlin%ind%jmirpn      = -1
        mlin%ind%jmirps      = -1
        mlin%keast           = -1
        mlin%kwest           = -1
        ifail = -42099
        return
C
      endif
C      
C the job is done
c
c
C close the segment
C
      CALL UF421(2, mpx, pb, palt, xlen, ji1, ifail)
      if( ifail .lt. 0 )return
      mlin%ind%jbeg        = ji0
      mlin%ind%jend        = ji1
C      
C when necessary, transpose the field line segment
C
      if( ( kmpbeg.lt.kmpend .and. 
     :      mpbeg%colat .gt. mpend%colat ) .or.
     :    ( kmpbeg.gt.kmpend .and.
     :      mpbeg%colat .lt. mpend%colat ) )then
        CALL UF429(ji0,ji1)
        kmpbeg             = ji1 + ji0 - kmpbeg
        kmpend             = ji1 + ji0 - kmpend
        if( lalt )then
          kbeg             = ji1 + ji0 - kbeg
          kend             = ji1 + ji0 - kend
        endif
        ks0                = ji1 + ji0 - ks0
      endif
C      
C store the data specific to the field line segment
C
C mirror points
      if( mpbeg%colat .le. mpend%colat )then
        mlin%ind%jmirpn    = kmpbeg
        mlin%ind%jmirps    = kmpend
      else
        mlin%ind%jmirpn    = kmpend
        mlin%ind%jmirps    = kmpbeg
      endif
C foot prints
      if( lalt )then
        if( abs( mbeg%radius-re-falt ) .gt. 2*epskm .or.
     :      abs( mend%radius-re-falt ) .gt. 2*epskm )then 
           ifail = -42011
           return
        endif
        if( mbeg%colat .le. mend%colat )then
          mlin%footpn      = mbeg
          mlin%footps      = mend
        else
          mlin%footpn      = mend
          mlin%footps      = mbeg
        endif
      else
        mlin%footpn%radius = -1.0d0
        mlin%footpn%colat  =  0.0d0
        mlin%footpn%elong  =  0.0d0
        mlin%footps%radius = -1.0d0
        mlin%footps%colat  =  0.0d0
        mlin%footps%elong  =  0.0d0
      endif
C lowest magnetic field intensity point
      mlin%dtdft           =  0.0d0
      call UF424( mlin%equat,
     :           pb, palt,
     :           mlin%drift, ifail )
      if( ifail .lt. 0 )  return
C
      ifail                = ji1 - ji0 + 1
      kgp                  = ks0
C     the current point is set to the initial point, when this point
C     lies inside the segment
      if( l426 )then
        ifail = ifail+nx170
c       test: uf420 -> uf422 -> uf426 -> error
c       test: try to correct the error
c       test: side effect = possibility of incorrect Beq 
      endif
C
C
      END
C----------------------------------------------------------------------

      SUBROUTINE UF421
     :          (kiniclo, mgeo, pb, palt, plen, jind, ifail)
C
C!    Initialize and close a field line segment 
C
      INCLUDE 'structure.h'
C
      EXTERNAL UM530, UT999
C
C     INTERFACE
C
        REAL*8        pb, palt, plen
        INTEGER*4     ifail, jind, kiniclo
        TYPE(zgeo) :: mgeo     
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
C     VARIABLE
C
        TYPE(zgeo) :: mgeop        
C                     
C *
C *  kiniclo  in            Select initialisation(1) or ending(2)
C *  mgeo     in(1)/out(2)  Start/end point location
C *  pb, palt out           Magnetic field intensity and altitude
C *  plen     out           Arc length of the field line segment
C *  jind     out           Index in UC170
C *
C
C     CODE
C
      ifail                         = 0
C      
      if( kiniclo .eq. 1 )then
C
C  Initialize the position and the arc length, evaluate the radius of 
C  curvature (from a dipolar magnetic field), the magnetic field 
C  vector and the magnetic field intensity
C
        nsg                         = 0
        jind                        = 1            
        mele(jind)%beg%coord        = mgeo
        CALL UM530 ( mgeo, mele(jind)%beg%b, ifail )
        pb                          = mele(jind)%beg%b%dnrm
        if( ifail .lt. 0 )return
        mele(jind)%arcl             = 0.0d0
        mele(jind)%beg%b%dnrm       = pb
        mgeop                       = mgeo        
        CALL UT999 ( mgeop, mele(jind)%beg%rcurv )
        plen                        = 0.0d0
C        
      else
C      
        jind                        = kgp        
        if( jind .ge. nx170 )then
          ifail                     = -42102
          return
        elseif( jind .lt. nsg )then
          ifail                     = -42103
          return
        elseif( jind .eq. nsg )then
          ifail                     = 999
c          print*,'UF421: single point!'
        endif
        mgeo                        = mele(jind)%beg%coord
        pb                          = mele(jind)%beg%b%dnrm
        mele(jind)%rkstp(2)         = 0.0d0
        mele(jind)%rkstp(3)         = 0.0d0
        mele(jind)%rkstp(1)         = 0.0d0
        plen                        = ABS( mele(nsg+1)%arcl
     :                                   - mele(jind)%arcl )
C        
      endif
C
C
      nsg                           = jind
      kgp                           = jind
      palt                          = mgeo%radius - re
C
C
      END
C----------------------------------------------------------------------
      SUBROUTINE UF422
     :          (kdir, fbm, falt, scale, mpos, jind, ifail)
C
C!    Follow a magnetic field line until a boundary condition is reached
C
      INCLUDE 'structure.h'
C
      EXTERNAL UF423, UF426, UF427, UF428       
C
C     INTERFACE
C
        REAL*8        fbm, falt, scale
        INTEGER*4     kdir, ifail, jind
        TYPE(zgeo) :: mpos        
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
C     VARIABLES
C
        REAL*8        dpbm, dalt, pb, palt
        REAL*8        rc, pbx(0:1), dpbx(0:1)
        INTEGER*4     ifl, kbeg, k0, k1, ktype, k426, kcnt
        LOGICAL       lalt, ltb, lta, lprev
        REAL*8        scal0, dstep0
C        
C *
C *    KDIR    in     direction = +1 : higher B
C *                   |           -1 : lower B
C *                   |            0 : no check
C *    FBM     in     boundary condition on the magnetic field intensity
C *    FALT    in     boundary condition on the altitude
C *    SCALE   in/out guess direction (+/- 1.0d0)
C *    MPOS     out    location where sub. stops
C *    JIND    out    indice of POS in UC170
C *   
C *    dpbm        initial delta on the magn. field intensity
C *    dalt        initial delta on the altitude
C *    k0, k1      =0/1
C *    pbx(k0/k1)  previous/actual value of the magnetic field
C *    dpbx(k0/k1) previous/actual delta
C *    ifl         temporary ifail
C *    ktype       when extremum, +/- 1: maximum/minimum
C *    epsrel      relative precision on the magnetic field 
C *                | intensity: 0.001 [1]
C *
C **
C **  The subroutine UF422 follow the field line until one of
C **  the boundary conditions is OVERreached. The fine tuning has
C **  to be done later
C **
C
C     CODE
C
c      print*,'uf422:',kdir, fbm, falt, scale

      kcnt = 0
      ifl     = 0
      if( ABS(scale).lt.0.01d0 .or. ABS(scale).gt.10.0d0 )then
        ifail = -42201
        return
      endif
      k0      = 0
      k1      = 1
c>
c>      write(2,*)'uf422:start=',kgp,mele(kgp)%beg%coord%radius,
c>     :     mele(kgp)%beg%coord%colat,mele(kgp)%beg%coord%elong,
c>     :     mele(kgp)%beg%b%dnrm   
C      
C  a/ Control the boundary conditions
C     -------------------------------
C
      pbx(k0)  = mele(kgp)%beg%b%dnrm
      palt     = mele(kgp)%beg%coord%radius - re
      dpbm     = pbx(k0) - fbm
      lalt     = falt .ge. -500.0d0
c MK, June 15, 1998: 5 lines added
      lprev    = .false.
      if( kgp.gt.1 )then
        lprev     = .true.
        dpbx(k0)  = pbx(k0)-mele(kgp-1)%beg%b%dnrm
      endif
c --
      if( lalt )then
        dalt   = palt - falt
        if( ABS(dalt) .lt. epskm )then
          ifl  = 1
        endif
        if( ABS(dpbm) .lt. 0.3d0*epsrel*fbm )then
          ifl  = ifl + 2
        endif
      else
        dalt   = 0.0d0
C        --> test on the altitude is always false
        if( ABS(dpbm) .lt. 0.3d0*epsrel*fbm )then
          ifl  = 2
        endif
      endif
      scal0    = scale
      rc       = mele(kgp)%beg%rcurv
      dstep0   = MAX( stpmin, MIN( prop*rc, stepx*
     :                mele(kgp)%beg%coord%radius ) )
      scale    = SIGN(0.9d0 * stplst / dstep0,scal0)
C               
C  b/ Check the direction
C     -------------------
C
      kbeg       = kgp
      rc         = mele(kgp)%beg%rcurv
      CALL UF423 (scale, rc, pbx(k1), palt, ifail)
      if( ifail .lt. 0 ) return
c>
c>      write(2,*)'uf422:uf423->',kgp,mele(kgp)%beg%coord%radius,
c>     :     mele(kgp)%beg%coord.colat,mele(kgp)%beg%coord%elong,
c>     :     mele(kgp)%beg%b%dnrm
      dpbx(k1)   = pbx(k1) - pbx(k0)
      if( dpbx(k1) * kdir .lt. 0.0d0 )then
c>
c>        write(2,*)'uf422: scale=-scale'
c MK, June 15, 1998: 2 lines added
        lprev    = .true.
        dpbx(k0) = -dpbx(k1)
c --
        scale    = -scale
        kgp      = kbeg
        rc       = mele(kgp)%beg%rcurv
        CALL UF423 (scale, rc, pbx(k1), palt, ifail)
      if( ifail .lt. 0 ) return
c>
c>        write(2,*)'uf422:uf423->',kgp,mele(kgp)%beg%coord%radius,
c>     :     mele(kgp)%beg%coord.colat,mele(kgp)%beg%coord%elong,
c>     :     mele(kgp)%beg%b%dnrm
        dpbx(k1) = pbx(k1) - pbx(k0)
        if( dpbx(k1) * kdir .lt. 0.0d0 ) then
C
C IMPROVED JAN 1998
C
          if( kdir.lt.0 .and. ifl.eq.2 )then
            IFAIL = 152
            kgp   = kbeg
            jind  = kgp
            mpos  = mele(jind)%beg%coord
            scale = SIGN( ABS(scal0), scale )
c            print*,'UF422 improved',kgp,ifail
            return
          else
            ifail  = -42205
c MK, JUNE 12, 1998: 4 lines added:
            kgp   = kbeg
            jind  = kgp
            mpos  = mele(jind)%beg%coord
            mlin%equat=mele(jind)%beg
            return
          endif
C
C
C
        endif
      endif
      scale      = SIGN( ABS(scal0), scale )
C      
C   boundary condition for the first step
C
      ltb = ( pbx(k1) - fbm ) * dpbm .lt. 0.0d0 
      lta = (  palt - falt  ) * dalt .lt. 0.0d0   
      if( ifl .gt. 0)then
        if( MOD(ifl,2) .eq. 1 ) then
          lta                 = .false.
          dalt                = palt - falt
        endif
        if( ifl/2 .eq. 1 ) then
c MK, June 15, 1998: 2 lines replaced by 20 lines
c                    ltb      = .false.
c                    dpbm     = pbx(k1) - fbm
          if( .not. lprev )then
            kbeg     = kgp
            scal0    = - sign(1.85D0*stplst/dstep0,scale)
            call UF423(scal0,rc,pb,palt,ifail)
            if( ifail.lt.0 )return
            kgp      = kbeg
            rc       = mele(kgp)%beg%rcurv
            lprev    = .true.
            dpbx(k0) = pbx(k0) - pb
cccccccccc            print*,'UF422 to be checked (June 15, 1998)'
          endif
          if(dpbx(k0)*dpbx(k1).gt.0.0D0)then
            ltb               = .false.
            dpbm              = pbx(k1) - fbm
          else
            kgp               = kgp-1
            ifl               = 152
            if( dpbx(k0).gt.0.0D0 ) ifl = 132
            goto 200
          endif
c --             
        endif
      endif
C  c/ Loop until a boundary condition is encountered
C     ----------------------------------------------
C
  100 continue
      kcnt = kcnt + 1
      IF (kcnt .GT. 200) THEN
        ifail = -42206
        RETURN
      END IF
C  
C  c.1/ check the boundary conditions
C
C
      if( ltb )then
C
c>
c>        write(2,*)'uf422:uf427 interpolate to b ',kgp
        CALL UF427(fbm, pbx(k1), palt, ifail)
c>
c>        write(2,*)'uf422:uf427->',kgp,mele(kgp)%beg%coord%radius,
c>     :     mele(kgp)%beg%coord%colat,mele(kgp)%beg%coord%elong,
c>     :     mele(kgp)%beg%b%dnrm
        if( ifail .lt. 0 ) return
        if( .not. lalt )then
          ifl                 = ifl + 100
          goto 200
        else if( ABS( palt - falt ) .le. epskm )then
          ifl                 = ifl + 300
          goto 200
        else if( ( palt - falt ) * dalt .gt. 0.0d0  )then
          ifl                 = ifl + 100
          goto 200
        else if(
     :    ( mele(kgp-1)%beg%coord%radius - re - falt ) * 
     :                               dalt .lt. 0.0d0 )then
C         go back and look the alt. boundary
          kgp                 = kgp - 1
        else
C         go back and look the alt. boundary
          continue
        endif
C
      endif
C
      if( lta )then
C
        CALL UF428(falt, pbx(k1), palt, ifail)
        if( ifail .lt. 0 ) return
        ifl                   = ifl + 200
        if( ABS( fbm - pbx(k1) ) .le. fbm * 0.3d0*epsrel )
     :    ifl                 = ifl + 100
        goto 200
C
      endif
C
C  c.2/ step once along the field line
C
      k1                      = k0
      k0                      = 1 - k0
      CALL UF423 (scale, rc, pbx(k1), palt, ifail)     
      if( (ifail .le. -53001 .and. ifail .ge. -53004) .or.
     :     ifail .eq. -42305 .or. ifail .eq. -42304 )then
        scal0    = scale
        dstep0   = MAX( stpmin, MIN( prop*rc, stepx*
     :                     mele(kgp)%beg%coord%radius ) )
        scale    = SIGN( stplst / dstep0,scal0)
        CALL UF423 (scale, rc, pbx(k1), palt, ifail)     
        if( ifail.lt.0 )then
          return
        endif
        scale = scal0
      elseif( ifail .lt. 0 ) then
        return
      endif
c>
c>      write(2,*)'uf422:uf423->',kgp,mele(kgp)%beg%coord%radius,
c>     :     mele(kgp)%beg%coord%colat,mele(kgp)%beg%coord%elong,
c>     :     mele(kgp)%beg%b%dnrm
      dpbx(k1)                = pbx(k1) - pbx(k0)
C      
C  c.3/ check for an extremum in the magnetic field intensity
C
      if( dpbx(k0) * dpbx(k1) .LE. 0.0d0 )then
        ifl                   = MOD(ifl,10)+10
        mele(kgp)%beg%rcurv = rc
C        
C  c.3.1/ move the previous point to the extremum location
C
c>
c>        write(2,*)'uf422:uf426 interpolate to extr. ',kgp
        k426 = 0
        CALL UF426 (pb, palt, ktype, ifail)
c>
c>        write(2,*)'uf422:uf426->',kgp,mele(kgp)%beg%coord%radius,
c>     :     mele(kgp)%beg%coord%colat,mele(kgp)%beg%coord%elong,
c>     :     mele(kgp)%beg%b%dnrm
        if( ifail .lt. 0 ) then
          k426 = ifail
c-test-test-test-test-          return
          ifail = kgp-1
          ifl = mod(ifl,1000)+1000
        endif
c>
c>        write(2,*)'       (extr)',ifail,mele(ifail)%beg%coord%radius,
c>     :     mele(ifail)%beg%coord%colat,mele(ifail)%beg%coord%elong,
c>     :     mele(ifail)%beg%b%dnrm
c
c  NOTE:  now ifail points to the extremum in mele()
c         and kgp= ifail+1 
c
C        
C  c.3.2/ test if the extremum corresponds to a boundary condition
C
          if( ABS( fbm - pb ) .le. fbm * 0.3d0*epsrel )then
            ifl                 = ifl + 30 - ktype * 10
            kgp                 = ifail
c                               ^^^^^^^ see note
            ifl                 = ifl + 100
            if( lalt .and. (ABS( palt - falt ) .le. epskm) )
     :          ifl             = ifl + 200
            if( k426 .LT. 0 )then
              ifail=K426
              return
            endif
            GOTO 200
          else
     :    if( (pb-fbm)*dpbm .lt. 0.0d0 )then
            kgp                 = ifail
c                               ^^^^^^^ see note
            ltb = .true.
            if(kgp.gt.1)then 
              if((mele(kgp-1)%beg%b%dnrm-fbm)*dpbm .lt. 0.0d0 )
     :          kgp=kgp-1
            endif
            palt=mele(kgp)%beg%coord%radius-re
            lta = (palt-falt)*dalt .lt. 0.0d0
c...rajout
            if( kgp.gt.1 )then
              pbx(k0) = mele(kgp-1)%beg%b%dnrm
            else
              pbx(k0) = fbm
            endif
            pbx(k1) = mele(kgp)%beg%b%dnrm
            palt    = mele(kgp)%beg%coord%radius-re
            dpbx(k1)= pbx(k1) - pbx(k0)
c...
            GOTO 100
          else
     :    if( lalt .and. ABS( palt - falt ) .le. epskm )then
            ifl                 = ifl + 30 - ktype * 10
            kgp                 = ifail
c                               ^^^^^^^ see note
            ifl                 = ifl + 200
            if( k426 .LT. 0 )then
              ifail=K426
              return
            endif
            GOTO 200
          else
     :    if( (palt-falt)*dalt .lt. 0.0d0 )then
            kgp                 = ifail
c                               ^^^^^^^ see note
            ltb = .false.
            lta = .true.
            if(kgp.gt.1)then 
              if((mele(kgp-1)%beg%coord%radius-re-falt)*dalt .lt. 0.0d0)
     :          kgp=kgp-1
            endif
c...rajout
            if( kgp.gt.1 )then
              pbx(k0) = mele(kgp-1)%beg%b%dnrm
            else
              pbx(k0) = fbm
            endif
            pbx(k1) = mele(kgp)%beg%b%dnrm
            palt    = mele(kgp)%beg%coord%radius-re
            dpbx(k1)= pbx(k1) - pbx(k0)
c...
            GOTO 100
          endif
C              
C  c.3.3/ when the extremum is a minimum, check that, at least, one
C         of the boundary conditions may be reached 
C 
        if( ktype .lt. 0 )then
          if( ( lalt      .and. ( falt .gt. palt ) 
     :                    .and. (  fbm .lt. pb   ) )
     :    .or.( .not.lalt .and. (  fbm .lt. pb   ) )   )then      
            kgp               = ifail
c                               ^^^^^^^ see note
            ifail             = -42202
            jind              = kgp
            mpos              = mele(jind)%beg%coord
c MK, JUNE 12, 1998 add 1 line:
            mlin%equat=mele(jind)%beg
            return
          endif
        endif
c
c
c
        pbx(k0) =mele(kgp-1)%beg%b%dnrm
        pbx(k1) =mele(kgp)%beg%b%dnrm
        palt    =mele(kgp)%beg%coord%radius-re
        dpbx(k1)= pbx(k1) - pbx(k0)
      endif
C
C    boundary conditions for the next step
C
      ltb = ( pbx(k1) - fbm ) * dpbm .lt. 0.0d0 
      lta = (  palt - falt  ) * dalt .lt. 0.0d0       
C>
c>      write(2,*)' pbx',pbx(k0),pbx(k1),dpbx(k1),fbm
C
      goto 100
C
C  d/ End of the loop
C     ---------------
C 
  200 continue 
      mele(kgp)%beg%rcurv     = rc
      jind                    = kgp
      mpos                    = mele(jind)%beg%coord
      ifail                   = ifl
C
C
      END
C----------------------------------------------------------------------
      SUBROUTINE UF423
     :          (scale, rc, pb, palt, ifail)
C
C!    Step along the magnetic field line 
C
      INCLUDE 'structure.h'
C     
      EXTERNAL UM530, UF425
C      
C     INTERFACE
C
        REAL*8        scale, rc, pb, palt   
        INTEGER*4     ifail           
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
     :               /prop, stepx, stpmin, umsq, upsq, uk2,
     :                uk3, epskm, epsrel, stplst, xclat,
     :                kmflg, nxstp       
C
        REAL*8        prop, stepx, stpmin, umsq, upsq, uk2,
     :                uk3, epskm, epsrel, stplst, xclat
        INTEGER*4     kmflg, nxstp              
C                                       
C     VARIABLES
C
        REAL*8        dstep0, dstep, p5ds, p6ds, p7ds, rkm, sth
        INTEGER*4     k,ifl, k530
        REAL*8        dy2(3)
        REAL*8        q3(3), dy3(3)
        REAL*8        q4(3), dy4(3)
        REAL*8        qf(3), dyf(3)
        TYPE(zgeo) :: mys1, mys2, mys3, mys4
        TYPE(zvec) :: mvb1, mvb2, mvb3, mvb4, md1tds, md2tds
C *
C *   rc              radius of curvature [km]
C *   dstep           integration step size [km]
C *   rc,             radius of curvature [km]
C *   prop            0.2 [1]
C *   stepx           # new 0.4 [1] # old 1.0 * re [km]
C *   stpmin          1.0 [km]
C *   umsq            1-sqrt(.5)
C *   upsq            1+sqrt(.5)
C *   dy2(1..3)       comp. rho, theta, phi  [km,rad,rad]
C *   mys2            position of the first R-K trial step
C *   mvb2            corresponding magnetic field vector
C *                   and intensity
C *   dy3(1..3), mys3, mvb3   idem for the 2nd try
C *   dy4(1..3), mys4, mvb4   idem for the 3rd try
C *   dyf(1..3)       the final R-K step
C *
C *   md1tds          d 1_t / d s = 1/Rc 1_n [km-1] at the initial pos.
C *   md2tds          d 1_t / d s = 1/Rc 1_n [km-1] at the final pos.
C *
C
C     CODE
C
      ifail = 0
      ifl=0
      k530 = 0
      if( ABS(scale) .lt. 1.0d-10 )then
        ifail = -42302
        if (.false.) print*,'UF423',ifail,scale
        return
      endif
      if( ABS(scale) .gt. 10.0d0 )then
        ifail = -42303
        return
      endif
      k = kgp + 1
      if( k .ge. nx170 )then
        ifail = -42301
        return
      endif
C
      mele(kgp)%rkstp(1)=0.d0
      mele(kgp)%rkstp(2)=0.d0
      mele(kgp)%rkstp(3)=0.d0
C      
C  a/ Determine the step size
C     -----------------------
C
 100  rkm         = mele(kgp)%beg%coord%radius
      dstep0      = MAX( stpmin, MIN( prop*rc, stepx*rkm ) )
      dstep       = dstep0 * scale
c      print*,'uf423',kgp,rc,rkm,dstep, stepx*rkm
      p5ds        = 0.5d0 * dstep
      p6ds        = umsq * dstep
      p7ds        = upsq * dstep
C      
C  b/ Gill's Method
C     -------------
C
C  b.1/ start point
C
      mys1        = mele(kgp)%beg%coord
      sth         = SIN( mele(kgp)%beg%coord%colat * deg )
C      
C  b.2/ first trial step "0.5 k1"
C
      dy2(1)      = p5ds * mele(kgp)%beg%b%rho
     :                   / mele(kgp)%beg%b%dnrm
      dy2(2)      = p5ds * mele(kgp)%beg%b%theta
     :                   / mele(kgp)%beg%b%dnrm /rkm
      dy2(3)      = p5ds * mele(kgp)%beg%b%phi
     :                   / mele(kgp)%beg%b%dnrm /rkm/sth 
      mys2%radius = mys1%radius + dy2(1) 
      mys2%colat  = mys1%colat + dy2(2) / deg
      mys2%elong  = mys1%elong + dy2(3) / deg
C      
      if( mys2%colat .lt. 0.0d0 )then
        mys2%colat=-mys2%colat
        mys2%elong=mod(mys2%elong+180.0d0,360.0d0)
      elseif( mys2%colat .gt. 180.0d0 )then
        mys2%colat=360.0D0-mys2%colat
        mys2%elong=mod(mys2%elong+180.0d0,360.0d0)
      endif
      CALL UM530 ( mys2, mvb2, ifail )
      if( ifail .lt. 0 ) then
        if( ifail.gt.-53002 .or. ifail.lt.-53004 .or.
     :      mvb2%dnrm.lt.1.e-6 .or. mys2%colat.lt.1.0 .or.
     :      mys2%colat.gt.179.0) return
        k530 = ifail
        ifail=0
      endif
      mele(kgp)%rkstp(2) = mvb2%dnrm
C
      rkm         = mys2%radius
      sth         = SIN( mys2%colat * deg )
C
C  b.3/ second trial step "(-.5+sqrt(.5)) k1 + (1-sqrt(.5)) k2"
C
      q3(1)       = p6ds * mvb2%rho
     :                   / mele(kgp)%rkstp(2)          
      q3(2)       = p6ds * mvb2%theta 
     :                   / mele(kgp)%rkstp(2) /rkm     
      q3(3)       = p6ds * mvb2%phi
     :                   / mele(kgp)%rkstp(2) /rkm/sth 
      dy3(1)      = q3(1) - 2*umsq * dy2(1)
      dy3(2)      = q3(2) - 2*umsq * dy2(2)
      dy3(3)      = q3(3) - 2*umsq * dy2(3)
      q3(1)       = 2 * dy2(1) + 3 * dy3(1) - q3(1)
      q3(2)       = 2 * dy2(2) + 3 * dy3(2) - q3(2)
      q3(3)       = 2 * dy2(3) + 3 * dy3(3) - q3(3)
      mys3%radius = mys2%radius + dy3(1) 
      mys3%colat  = mys2%colat  + dy3(2) / deg
      mys3%elong  = mys2%elong  + dy3(3) / deg
C
      if( mys3%colat .lt. 0.0d0 )then
        mys3%colat=-mys3%colat
        mys3%elong=mod(mys3%elong+180.0d0,360.0d0)
      elseif( mys3%colat .gt. 180.0d0 )then
        mys3%colat=360.0D0-mys3%colat
        mys3%elong=mod(mys3%elong+180.0d0,360.0d0)
      endif
      CALL UM530 ( mys3, mvb3, ifail )
      if( ifail .lt. 0 ) then
        if( ifail.gt.-53002 .or. ifail.lt.-53004 .or.
     :      mvb3%dnrm.lt.1.e-6 .or. mys3%colat.lt.1.0 .or.
     :      mys3%colat.gt.179.0) return
        k530 = ifail
        ifail=0
      endif                     
      mele(kgp)%rkstp(3) = mvb3%dnrm
C
      rkm         = mys3%radius
      sth         = SIN( mys3%colat * deg )
C
C  b.4/ third trial step "-sqrt(.5) k2 + (1+sqrt(.5)) k3"
C
      q4(1)       = p7ds * mvb3%rho 
     :                   / mele(kgp)%rkstp(3)          
      q4(2)       = p7ds * mvb3%theta 
     :                   / mele(kgp)%rkstp(3) /rkm     
      q4(3)       = p7ds * mvb3%phi
     :                   / mele(kgp)%rkstp(3) /rkm/sth 
      dy4(1)      = q4(1) - upsq * q3(1)
      dy4(2)      = q4(2) - upsq * q3(2)
      dy4(3)      = q4(3) - upsq * q3(3)
      q4(1)       = q3(1) + 3 * dy4(1) - q4(1)
      q4(2)       = q3(2) + 3 * dy4(2) - q4(2)
      q4(3)       = q3(3) + 3 * dy4(3) - q4(3)
      mys4%radius = mys3%radius + dy4(1) 
      mys4%colat  = mys3%colat + dy4(2) / deg
      mys4%elong  = mys3%elong + dy4(3) / deg
C      
      if( mys4%colat .lt. 0.0d0 )then
        mys4%colat=-mys4%colat
        mys4%elong=mod(mys4%elong+180.0d0,360.0d0)
      elseif( mys4%colat .gt. 180.0d0 )then
        mys4%colat=360.0D0-mys4%colat
        mys4%elong=mod(mys4%elong+180.0d0,360.0d0)
      endif
      CALL UM530 ( mys4, mvb4, ifail )
      if( ifail .lt. 0 ) then
        if( ifail.gt.-53002 .or. ifail.lt.-53004 .or.
     :      mvb4%dnrm.lt.1.e-6 .or. mys4%colat.lt.1.0 .or.
     :      mys4%colat.gt.179.0) return
        k530 = ifail
        ifail=0
      endif
      mele(kgp)%rkstp(1) = mvb4%dnrm
C
      rkm         = mys4%radius
      sth         = SIN( mys4%colat * deg )
C
C  b.5/ compute the real step
C
C     "( k1 + 2(1-sqrt(.5)) k2 + 2(1+sqrt(.5)) k3 + k4 )/6"
      qf(1)       = p5ds * mvb4%rho 
     :                   / mele(kgp)%rkstp(1)          
      qf(2)       = p5ds * mvb4%theta
     :                   / mele(kgp)%rkstp(1) /rkm     
      qf(3)       = p5ds * mvb4%phi
     :                   / mele(kgp)%rkstp(1) /rkm/sth 
      dyf(1)      = ( qf(1) - q4(1) ) / 3.0d0
      dyf(2)      = ( qf(2) - q4(2) ) / 3.0d0
      dyf(3)      = ( qf(3) - q4(3) ) / 3.0d0
      mele(k)%beg%coord%radius = mys4%radius + dyf(1) 
      mele(k)%beg%coord%colat  = mys4%colat + dyf(2) / deg
      mele(k)%beg%coord%elong  = mys4%elong + dyf(3) / deg
C
      if( mele(k)%beg%coord%colat .lt. 0.0d0 )then
        mele(k)%beg%coord%colat=-mele(k)%beg%coord%colat
        mele(k)%beg%coord%elong=
     :                       mod(mele(k)%beg%coord%elong+180.d0,360.d0)
      elseif( mele(k)%beg%coord%colat .gt. 180.0d0 )then
        mele(k)%beg%coord%colat=360.0D0-mele(k)%beg%coord%colat
        mele(k)%beg%coord%elong=
     :                       mod(mele(k)%beg%coord%elong+180.d0,360.d0)
      endif
      CALL UM530 ( mele(k)%beg%coord, mvb1, ifail )
      if( ifail .lt. 0 ) then
        if( ifail.gt.-53002 .or. ifail.lt.-53004 .or.
     :      mvb1%dnrm.lt.1.e-6 ) return
        k530 = ifail
        ifail=0
      endif                         
      mele(k)%beg%b    = mvb1
      mele(k)%beg%b%dnrm = mvb1%dnrm
C
C
      pb           = mele(k)%beg%b%dnrm
      palt         = mele(k)%beg%coord%radius - re
      mele(k)%arcl = mele(kgp)%arcl + dstep
C
C
C  c/ Evaluate the radius of curvature
C
C
      CALL UF425 (mele(kgp)%beg%coord, mele(kgp)%beg%b,
     :            mele(k)%beg%coord, mele(k)%beg%b, dstep,
     :            md1tds, md2tds)
C
      mele(kgp)%beg%rcurv    = 1.d0 / md1tds%dnrm
c      print*,'uf423 rc old',rc,mele(kgp)%beg%rcurv
      rc                     = 1.d0 / md2tds%dnrm
c      print*,'uf423 rc new',rc,dstep0
C
c
      if( k530 .lt. 0 ) then
        ifail = -42305
        if( k530 .eq. -53004) ifail = -42304
        goto 500
      endif
c
c
c
C  d/ Check the step a posteriori
C
      if( dstep0 .gt. 1.5d0 * MAX( prop * rc, stpmin ) )then 
        ifl = ifl + 1
c        print*,'uf423          goto 100'
        if( ifl .lt. 20 )GOTO 100
        ifail = -42306
        goto 500
      endif
ctmptmptmptmp
c      if( ifl .ge. 1 )continue
c      print*,'uf423:step=',dstep
C
      ifail=ifl
c      print*,'uf423',mele(k)%beg%coord%radius,mele(k)%beg%coord%colat,
c     :         mele(k)%beg%b%dnrm
c
c
c
  500 continue
c
c
      mele(kgp)%dtbnd = 0.0d0
      kgp             = k
C
C
      END
C----------------------------------------------------------------------
      SUBROUTINE UF424
     :          ( mgp, pb, palt, mxnb, ifail )
C
C!    Search the point of a field line segment with the lowest
C!    magnetic field intensity 
C
      INCLUDE 'structure.h'
C      
      EXTERNAL UF423, UF425
C      
C     INTERFACE
C
        INTEGER*4     ifail
        REAL*8        pb, palt
        TYPE(zpnt) :: mgp
        TYPE(zvec) :: mxnb
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
C     VARIABLES
C
        INTEGER*4     k, k0, kbeg, kend
        REAL*8        scale, rc, dstep
        TYPE(zvec) :: mnr, mne
        TYPE(zxyz) :: ma, mtot
        REAL*8        va, st, ct, sp, cp
        REAL*8        dum1(2), dum2(2)
        TYPE(zvec) :: mdxnra, mdxnrb, mdum3, mdum4
C
C     CODE
C
C  a/ Search the minimum along the field line segment
C     -----------------------------------------------
C
      kbeg           = 1
      kend           = nsg
      mgp            = mele(kbeg)%beg
      pb             = mgp%b%dnrm
      kgp            = kbeg
C
      do k=kbeg+1, kend
        if( mele(k)%beg%b%dnrm .lt. pb )then
          pb         = mele(k)%beg%b%dnrm
          kgp        = k
        endif
      enddo
      k              = kgp
C
      mgp            = mele(kgp)%beg
      palt           = mgp%coord%radius - re
C      
C  b/ Evaluate the radius of curvature
C     --------------------------------
C
      if( nsg.gt.nx170-4 )then
        ifail        = -42401
        return
      endif
      k0             = nx170-2
      mele(k0)%beg%b%dnrm     = pb
      mele(k0)%arcl              = 0.0d0
      mele(k0)%beg               = mgp
C 
      scale                      = 0.2d0
      rc                         = mele(k0)%beg%rcurv 
      kgp                        = k0 
      CALL UF423( scale, rc, dum1(1), dum1(2), ifail )
      if(ifail.eq.-42302 .and. .false.)print*,'UF424 -> UF423 (1)',ifail
      if(ifail.lt.0)return
C      
      dstep                      = mele(kgp)%arcl
     :                            - mele(k0)%arcl
      mele(k0)%beg%b%dnrm        = mele(k0)%beg%b%dnrm
      mele(kgp)%beg%b%dnrm       = mele(kgp)%beg%b%dnrm
      CALL UF425 ( mele(k0)%beg%coord, mele(k0)%beg%b,
     :             mele(kgp)%beg%coord, mele(kgp)%beg%b,
     :             dstep, mdxnra, mdum3 )
C     
      ma%x           = mele(kgp)%beg%coord%radius *
     :                 SIN(mele(kgp)%beg%coord%colat * deg) *
     :                 COS(mele(kgp)%beg%coord%elong * deg)     
      ma%y           = mele(kgp)%beg%coord%radius *
     :                 SIN(mele(kgp)%beg%coord%colat * deg) *
     :                 SIN(mele(kgp)%beg%coord%elong * deg)     
      ma%z           = mele(kgp)%beg%coord%radius *
     :                 COS(mele(kgp)%beg%coord%colat * deg)          
C 
      scale          = -0.2d0
      kgp            = k0
      CALL UF423(scale, rc, dum2(1), dum2(2), ifail)
      if(ifail.eq.-42302 .and. .false.)print*,'UF424 -> UF423 (2)',ifail
      if(ifail.lt.0)return
C      
      dstep          = mele(kgp)%arcl - mele(k0)%arcl
      mele(k0)%beg%b%dnrm  = mele(k0)%beg%b%dnrm
      mele(kgp)%beg%b%dnrm = mele(kgp)%beg%b%dnrm
      CALL UF425 ( mele(k0)%beg%coord, mele(k0)%beg%b,
     :             mele(kgp)%beg%coord, mele(kgp)%beg%b,
     :             dstep, mdxnrb, mdum4 )
C     
      mtot%x         = ma%x - mele(kgp)%beg%coord%radius *
     :                    SIN(mele(kgp)%beg%coord%colat * deg) *
     :                    COS(mele(kgp)%beg%coord%elong * deg)  
      mtot%y         = ma%y - mele(kgp)%beg%coord%radius *
     :                    SIN(mele(kgp)%beg%coord%colat * deg) *
     :                    SIN(mele(kgp)%beg%coord%elong * deg) 
      mtot%z         = ma%z - mele(kgp)%beg%coord%radius *
     :                    COS(mele(kgp)%beg%coord%colat * deg)    
C
      kgp          = k
C
      mnr%rho         = ( mdxnra%rho + mdxnrb%rho ) * 0.5d0     
      mnr%theta       = ( mdxnra%theta + mdxnrb%theta ) * 0.5d0     
      mnr%phi         = ( mdxnra%phi + mdxnrb%phi ) * 0.5d0     
      mnr%dnrm        = SQRT(mnr%rho**2 + mnr%theta**2 + mnr%phi**2)
C
      va             = 1.0d0 /
     :                 SQRT( mtot%x**2 + mtot%y**2 + mtot%z**2 )
      ct             = COS(mele(kgp)%beg%coord%colat * deg)
      st             = SIN(mele(kgp)%beg%coord%colat * deg)
      cp             = COS(mele(kgp)%beg%coord%elong * deg)
      sp             = SIN(mele(kgp)%beg%coord%elong * deg)
      mne%rho   = va * 
     :      (st * cp * mtot%x + st * sp * mtot%y + ct * mtot%z)
      mne%theta = va *
     :      (ct * cp * mtot%x + ct * sp * mtot%y - st * mtot%z)
      mne%phi   = va * (-sp * mtot%x + cp * mtot%y)
C
      mgp%rcurv            = 1.d0 / mnr%dnrm
      mele(kgp)%beg%rcurv  = mgp%rcurv
C
      mxnb%rho       = ( mne%theta * mnr%phi - mne%phi * mnr%theta )
     :                 / pb
      mxnb%theta     = ( mne%phi   * mnr%rho - mne%rho * mnr%phi   )
     :                 / pb
      mxnb%phi       = ( mne%rho * mnr%theta - mne%theta * mnr%rho )
     :                 / pb
C
      mxnb%dnrm      = mnr%dnrm/pb
C
      ifail          = kgp
C
C
      END
C----------------------------------------------------------------------
      SUBROUTINE UF425
     :          (mgp1, mvb1, mgp2, mvb2, ds, mxnr1, mxnr2)
C
C!    Evaluate the curvature of the field line 
C
      INCLUDE 'structure.h'
C
C     INTERFACE
C
        TYPE(ZGEO) :: mgp1, mgp2
        TYPE(ZVEC) :: mvb1, mvb2, mxnr1, mxnr2
        REAL*8        ds  
C
      COMMON /UC160
     :               /pi, deg, re, gmagmo, eclipt, geoid, uma 
        REAL*8        pi, deg, re, gmagmo, eclipt, geoid(3), uma(30)
C                              
C
C     VARIABLES
C
        REAL*8        brb1, brb2, btb1, btb2, bpb1, bpb2, t1, t2
        REAL*8        dbrbds, dbtbds, dbpbds
        REAL*8        ur1, ur2, rc
        REAL*8        dstep,dy1(3),dy2(3)
        REAL*8        d2rds2,d2tds2,d2pds2,drds,dtds,dpds
        INTEGER*4     ifail
        TYPE(zgeo) :: mpm
        TYPE(zvec) :: mbm
C
C     CODE
C
c      print*,' uf425',mgp1%radius,mgp1%colat,mgp1%elong
c      print*,'     |',mvb1%rho,mvb1%theta,mvb1%phi,mvb1%dnrm
c      print*,'     |',mgp2%radius,mgp2%colat,mgp2%elong
c      print*,'     |',mvb2%rho,mvb2%theta,mvb2%phi,mvb2%dnrm
c      print*,'     +',ds
c
      brb1        =   mvb1%rho / mvb1%dnrm
      btb1        =   mvb1%theta / mvb1%dnrm
      bpb1        =   mvb1%phi / mvb1%dnrm
      brb2        =   mvb2%rho / mvb2%dnrm
      btb2        =   mvb2%theta / mvb2%dnrm
      bpb2        =   mvb2%phi  / mvb2%dnrm
      t1          =   1.0d0 / TAN( mgp1%colat * deg )
      t2          =   1.0d0 / TAN( mgp2%colat * deg )
      ur1         =   1.0d0 / mgp1%radius       
      ur2         =   1.0d0 / mgp2%radius
C
      dbrbds      =   ( brb2 - brb1 ) / ds
      dbtbds      =   ( btb2 - btb1 ) / ds
      dbpbds      =   ( bpb2 - bpb1 ) / ds      
C
      mxnr1%rho   =   dbrbds - btb1**2 * ur1 - bpb1**2 * ur1
      mxnr1%theta =   brb1 * btb1 * ur1 + dbtbds
     :                - t1 * ur1 * bpb1**2
      mxnr1%phi   =   ur1 * brb1 * bpb1 + t1 * ur1 * btb1 * bpb1
     :                + dbpbds
      mxnr1%dnrm  = sqrt(mxnr1%rho**2+mxnr1%theta**2+mxnr1%phi**2)     
C
      mxnr2%rho   =   dbrbds - btb2**2 * ur2 - bpb2**2 * ur2
      mxnr2%theta =   brb2 * btb2 * ur2 + dbtbds
     :                - t2 * ur2 * bpb2**2
      mxnr2%phi   =   ur2 * brb2 * bpb2 + t2 * ur2 * btb2 * bpb2
     :                + dbpbds
      mxnr2%dnrm  = sqrt(mxnr2%rho**2+mxnr2%theta**2+mxnr2%phi**2)     
C
      if( mxnr1%dnrm .gt. 0.00025d0 .or. mxnr2%dnrm .gt. 0.00025d0 )then
        dstep=ds*0.5d0
        dy1(1)      = dstep * mvb1%rho/ mvb1%dnrm
        dy1(2)      = dstep * mvb1%theta/ mvb1%dnrm /mgp1%radius
        dy1(3)      = dstep * mvb1%phi/ mvb1%dnrm 
     :                /mgp1%radius/sin(deg*mgp1%colat) 
        mpm%radius  = mgp1%radius + dy1(1)
        mpm%colat   = mgp1%colat  + dy1(2)/deg
        mpm%elong   = mgp1%elong  + dy1(3)/deg
        if( mpm%colat .lt. 0.0d0 )then
          mpm%colat=-mpm%colat
          mpm%elong=mod(mpm%elong+180.0d0,360.0d0)
        elseif( mpm%colat .gt. 180.0d0 )then
          mpm%colat=360.0D0-mpm%colat
          mpm%elong=mod(mpm%elong+180.0d0,360.0d0)
        endif
        call um530(mpm,mbm,ifail)
        if(ifail.ge.0)then
          dy2(1)      = dstep * mbm%rho/ mbm%dnrm
          dy2(2)      = dstep * mbm%theta/ mbm%dnrm /mpm%radius
          dy2(3)      = dstep * mbm%phi/ mbm%dnrm 
     :                  /mpm%radius/sin(deg*mpm%colat) 
          mpm%radius  = mgp1%radius + 0.5d0*(dy1(1)+dy2(1))
          mpm%colat   = mgp1%colat  + 0.5d0*(dy1(2)+dy2(2))/deg
          mpm%elong   = mgp1%elong  + 0.5d0*(dy1(3)+dy2(3))/deg
c
          d2rds2      = (mgp1%radius-2*mpm%radius+mgp2%radius)/dstep**2
          d2tds2      = deg*(mgp1%colat-2*mpm%colat+mgp2%colat)/dstep**2
          d2pds2      = deg*(mgp1%elong-2*mpm%elong+mgp2%elong)/dstep**2
          if( mxnr1%dnrm .gt. 0.00025d0 )then
            rc=1.0d0/mxnr1%dnrm
            drds=(-1.5d0*mgp1%radius+2.d0*mpm%radius-0.5d0*mgp2%radius)
     :           /dstep
            dtds=deg*(-1.5d0*mgp1%colat+2.d0*mpm%colat
     :                -0.5d0*mgp2%colat)/dstep
            dpds=deg*(-1.5d0*mgp1%elong+2.d0*mpm%elong
     :                -0.5d0*mgp2%elong)/dstep
            dbrbds=d2rds2
            dbtbds=drds*dtds+mgp1%radius*d2tds2
            dbpbds=(drds*dpds+mgp1%radius*d2pds2)*sin(mgp1%colat*deg)+
     :            mgp1%radius*cos(mgp1%colat*deg)*dtds*dpds
              mxnr1%rho   =   dbrbds - btb1**2 * ur1 - bpb1**2 * ur1
              mxnr1%theta =   brb1 * btb1 * ur1 + dbtbds
     :                        - t1 * ur1 * bpb1**2
              mxnr1%phi   =   ur1 * brb1 * bpb1 + t1 * ur1 * btb1 * bpb1
     :                        + dbpbds
              mxnr1%dnrm  = sqrt(mxnr1%rho**2+mxnr1%theta**2+
     :                           mxnr1%phi**2)
c            if(abs(rc-1.0d0/mxnr1%dnrm).gt. 0.1d0*rc )then
c                print*,'uf425 spec',rc,1.0d0/mxnr1%dnrm
c            endif     
          endif
          if( mxnr2%dnrm .gt. 0.00025d0 )then
            rc=1.0d0/mxnr2%dnrm
            drds=(0.5d0*mgp1%radius-2.d0*mpm%radius+1.5d0*mgp2%radius)
     :           /dstep
            dtds=deg*(0.5d0*mgp1%colat-2.d0*mpm%colat
     :                +1.5d0*mgp2%colat)/dstep
            dpds=deg*(0.5d0*mgp1%elong-2.d0*mpm%elong
     :                +1.5d0*mgp2%elong)/dstep
            dbrbds=d2rds2
            dbtbds=drds*dtds+mgp2%radius*d2tds2
            dbpbds=(drds*dpds+mgp2%radius*d2pds2)*sin(mgp2%colat*deg)+
     :            mgp2%radius*cos(mgp2%colat*deg)*dtds*dpds
              mxnr2%rho   =   dbrbds - btb2**2 * ur2 - bpb2**2 * ur2
              mxnr2%theta =   brb2 * btb2 * ur2 + dbtbds
     :                        - t2 * ur2 * bpb2**2
              mxnr2%phi   =   ur2 * brb2 * bpb2 + t2 * ur2 * btb2 * bpb2
     :                        + dbpbds
              mxnr2%dnrm  = sqrt(mxnr2%rho**2+mxnr2%theta**2+
     :                           mxnr2%phi**2)     
c            if(abs(rc-1.0d0/mxnr2%dnrm).gt. 0.1d0*rc )then
c                print*,'uf425 sp_c',rc,1.0d0/mxnr2%dnrm
c            endif     
          endif
        endif
      endif                 
C
c      print*,'uf426 nrm',1/mxnr1%dnrm
C
      END
C----------------------------------------------------------------------
      SUBROUTINE UF426
     :          (bextr, aextr, ktype, ifail)
C
C!    Interpolating an extremum of the magnetic field intensity 
C
      INCLUDE 'structure.h'
C      
      EXTERNAL UF423
C      
C     INTERFACE
C
        REAL*8        bextr, aextr
        INTEGER*4     ifail, ktype        
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
C
C     VARIABLES
C
        INTEGER*4     kpa, kpb, ksave, kcount, kpc
        REAL*8        bab, bcb, sab, scb, cf1, cf2, delta
        REAL*8        smb, dstep0, scale, smb2
        REAL*8        rc, pb, palt, sma
        TYPE(zseg) :: msava, msavb
        LOGICAL       lfail
C
C *
C * epskm            precision to find an extremum in the magnetic
C *                  |   field intensity: 1.0 [km]
C * kpa, kpb         "a" = "c" - 2; "b" = "c" - 1
C * ksave            memorize the value of kgp
C * store            memorize "a"
C * cf1, cf2         B - B_b = cf1 * ( S - S_b ) + cf2 * ( S - S_b )^2
C * smb              value of the arc length from "b" to the extremum
C * dstep0           evaluation of the initial step size
C * sda              arc length from "a" to "old b"
C * sna              arc length from "a" to "new b"
C * kcount           number of iterations
C * bextr            magnetic field intensity at the extremum
C * aextr            altitude at the extremum
C * ktype            = +/- 1: maximum/minimum
C *
c
c     CODE
c
      cf2       = 0.0d0
      kcount    = 0
      lfail     =.false.
C      
C  a/ Save the index of the current point
C     -----------------------------------
C
      ksave     = kgp
      kpc       = kgp
      kpa       = kpc - 2
      kpb       = kpc - 1
      bextr     = mele(kpb)%beg%b%dnrm
      aextr     = mele(kpb)%beg%coord%radius - re        
      ktype     = 0
      msava = mele(kpa)
      msavb = mele(kpb)
C      
C  b/ Fit to a parabolic function
C     ---------------------------
C
  150 continue
      bab       = mele(kpa)%beg%b%dnrm - mele(kpb)%beg%b%dnrm
      bcb       = mele(kpc)%beg%b%dnrm - mele(kpb)%beg%b%dnrm
      sab       = mele(kpa)%arcl - mele(kpb)%arcl
      scb       = mele(kpc)%arcl - mele(kpb)%arcl
      delta     = sab * scb * ( scb - sab )
      if( abs(delta) .lt. 1.e-6 )then
c        print*,'uf426',ifail
        ifail   = -42601
        mele(kpa)=msava
        mele(kpb)=msavb   
        kgp     = ksave         
        bextr     = mele(kpb)%beg%b%dnrm
        aextr     = mele(kpb)%beg%coord%radius - re        
        return
      endif
      cf2       = ( bcb * sab - bab * scb ) / delta
      cf1       = ( bab / sab - cf2 * sab 
     :            + bcb /scb - cf2 * scb ) * 0.5d0
C     
C  c/ Estimate the position of the extremum
C     -------------------------------------
C
      if( cf2 .eq. 0 )then
c        print*,'uf426',ifail
        ifail   = -42602
        mele(kpa)=msava
        mele(kpb)=msavb   
        kgp     = ksave
        bextr     = mele(kpb)%beg%b%dnrm
        aextr     = mele(kpb)%beg%coord%radius - re        
        return
      endif           
      smb       = - cf1 / cf2 * 0.5d0
      ktype     = -SIGN(1.0d0,cf2)
      if( ABS(smb) .lt. epskm ) then
C      
C  c.1/ stop when the extremum is reached
C
        GOTO 100
      endif      
C     distance from "a" to the extremum
      sma       = smb - sab
C      
C  c.2/ estimate the scale
C
      rc        = mele(kpa)%beg%rcurv
      dstep0    = MAX( stpmin, MIN( prop*rc, stepx*
     :                mele(kpa)%beg%coord%radius ) )
      scale     = sma / dstep0
      if( ABS(scale) .gt. 4.0d0 )then
c        print*,'uf426',ifail
        ifail   = -42604
        mele(kpa)=msava
        mele(kpb)=msavb   
        kgp     = ksave
        bextr     = mele(kpb)%beg%b%dnrm
        aextr     = mele(kpb)%beg%coord%radius - re        
        return
      endif  
C       
C  d/ Perform a single RK integration and iterate from b/
C     ---------------------------------------------------
C
      kgp       = kpa
      CALL UF423 (scale, rc, pb, palt, ifail)
c      if(ifail.eq.-42302)print*,'UF426 -> UF423 (a)',ifail
      if( ifail.lt.0 ) then
        mele(kpa)=msava
        mele(kpb)=msavb   
        kgp     = ksave         
        bextr     = mele(kpb)%beg%b%dnrm
        aextr     = mele(kpb)%beg%coord%radius - re        
        return
      endif
      mele(kpb)%beg%rcurv = rc
      kcount    = kcount+1
      bextr     = mele(kpb)%beg%b%dnrm
      aextr     = mele(kpb)%beg%coord%radius - re
C
      if( kcount .lt. 50 )GOTO 150
      lfail     = .true.
C
C  e/ Go from the extremum to point "c" with a RK integr.
C     ---------------------------------------------------
C
  100 continue
      sma = mele(kpb)%arcl-mele(kpa)%arcl
      smb = sign(0.9d0*stplst,sma)
      smb2 = sign(0.95d0*stplst,sma)
      if( abs(sma).gt.abs(smb2) )then
        rc        = mele(kpa)%beg%rcurv
        dstep0    = MAX( stpmin, MIN( prop*rc, stepx*
     :                mele(kpa)%beg%coord%radius ) )
        scale     = (sma-smb) / dstep0
        kgp       = kpa
        CALL UF423 (scale, rc, pb, palt, ifail)
c        if(ifail.eq.-42302)print*,'UF426 -> UF423 (b)',ifail
c       from a to new_b
        if( ifail.lt.0 ) then
          mele(kpa)=msava
          mele(kpb)=msavb   
          kgp     = ksave
          bextr     = mele(kpb)%beg%b%dnrm
          aextr     = mele(kpb)%beg%coord%radius - re        
          return
        endif         
        dstep0    = MAX( stpmin, MIN( prop*rc, stepx*
     :                mele(kgp)%beg%coord%radius ) )
        scale     = smb / dstep0
        CALL UF423 (scale, rc, pb, palt, ifail)
c        if(ifail.eq.-42302)print*,'UF426 -> UF423 (c)',ifail
c       from new_b to m=new_c
        if( ifail.lt.0 ) then
          mele(kpa)=msava
          mele(kpb)=msavb   
          kgp     = ksave
          bextr     = mele(kpb)%beg%b%dnrm
          aextr     = mele(kpb)%beg%coord%radius - re        
          return
        endif         
        mele(kpa)%beg%rcurv = rc
c
        ksave     = kgp
      else
        ksave     = kpb
      endif
      kgp         = ksave
      rc          = mele(kgp)%beg%rcurv
      dstep0      = MAX( stpmin, MIN( prop*rc, stepx*
     :                mele(kgp)%beg%coord%radius ) )
      scale       = smb / dstep0
      CALL UF423 (scale, rc, pb, palt, ifail)
c      if(ifail.eq.-42302)print*,'UF426 -> UF423 (e)',ifail
c     from m to d or new_c
        if( ifail.lt.0 ) then
          mele(kpa)=msava
          mele(kpb)=msavb   
          kgp     = ksave
          bextr     = mele(kpb)%beg%b%dnrm
          aextr     = mele(kpb)%beg%coord%radius - re        
          return
        endif         
C
      ifail=ksave
      if(lfail)then
c        print*,'uf426',ifail
        ifail=-42603
        return
      endif
C
C
      END
C----------------------------------------------------------------------
        SUBROUTINE UF427
     :          (fbm, bi, ai, ifail)
C
C!    Interpolate to a  magnetic field intensity
C
      INCLUDE 'structure.h'
C      
      EXTERNAL UF423
C       
C     INTERFACE
C
        REAL*8        fbm, bi, ai
        INTEGER*4     ifail        
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
C     VARIABLES
C
        REAL*8        s0a, sac, sbc, sba, s0c, bba, b0a, bbc, b0c,
     :                s0b, bac, sxa, s0a1, s0a2
c-v
c-        REAL*8        s0an, bca, sca
c-^
        INTEGER*4     kpb, kpa, kcount
        REAL*8        dstep0, scale, rc, pb, palt
        TYPE(zseg) ::  msavb, msava
C
C     CODE
C
      kpb       = kgp
      kpa       = kpb - 1
c
c-v
c-      write(2,'(a8,f12.5)') 'uf427   ',fbm
c-      write(2,'(a8,4f12.5)')'      a:',mele(kpa)%beg%coord%radius-re,
c-     :                   90.0-mele(kpa)%beg%coord%colat,
c-     :                   mele(kpa)%beg%coord%elong,
c-     :                   mele(kpa)%beg%b%dnrm
c-      write(2,'(a8,4f12.5)')'      b:',mele(kpb)%beg%coord%radius-re,
c-     :                   90.0-mele(kpb)%beg%coord%colat,
c-     :                   mele(kpb)%beg%coord%elong,
c-     :                   mele(kpb)%beg%b%dnrm
c-^
c
      kcount    = 0
c note:  point a is before the seeked value
c        point b is after or at the value
c
c note: msava have to be always before the value
c       msavb have to be always after the value
C
C Estimation of the correct location with a linear fit
C 
      msavb     = mele(kpb)
      pb        = msavb%beg%b%dnrm
      msava     = mele(kpa)
      if( (msava%beg%b%dnrm-fbm)*(msavb%beg%b%dnrm-fbm) .gt. 0 )then
        ifail   = -42708
        return
      endif
      b0a       = fbm - msava%beg%b%dnrm
      bba       = pb -  msava%beg%b%dnrm
      sba       = msavb%arcl - msava%arcl
      if( abs(bba) .lt. 1.0e-9 )then
        ifail   = -42701
        return
      elseif( abs(sba) .lt. stpmin*0.1 )then
        ifail   = -42705
        return
      endif
      s0a       = ( b0a / bba ) * sba
      s0b       = s0a - sba
c
c check if the job has already been done
c
      if( ABS( bba - b0a ) .lt. 0.2d0 * epsrel * fbm .and.
     :   abs( sba ) .lt. stplst )then
        ifail   = kcount
        bi      = pb 
        ai      = msavb%beg%coord%radius - re
        return
      endif
c
c check if the last step is not too large
c
  200 continue
      if( abs(s0a) .gt. stplst*0.9 )then
c       add a point
cc        write(2,'(a8,2f12.5)')'s0a,dlt:',s0a,stplst*0.9
        rc      = mele(kpa)%beg%rcurv
        dstep0  = MAX( stpmin, MIN( prop*rc, stepx*
     :                mele(kpa)%beg%coord%radius ) )
        scale   = sign((abs(s0a)-stplst*0.5), s0a ) / dstep0
        kgp     = kpa
        kcount  = kcount + 1
        CALL UF423 (scale, rc, pb, palt, ifail)
        if( ifail.lt.0 )return
        mele(kpb)%beg%rcurv        = rc                             
c       note:  point b is lost but MSAVB can always be used
c              the new point (c) is before or after the value 
c-v
c-        write(2,'(a8,4f12.5)')'      ->',mele(kpb)%beg%coord%radius-re,
c-     :                   90.0-mele(kpb)%beg%coord%colat,
c-     :                   mele(kpb)%beg%coord%elong,
c-     :                   mele(kpb)%beg%b%dnrm
c-^
c
        b0c     = fbm - pb
c
        if( b0c*b0a .ge. 0 )then
c         point c is before the value
c
          msava=mele(kpb)
c         estimation of the correct location with a linear fit between c and msavb
          bbc     = msavb%beg%b%dnrm - pb
          if( abs(bbc) .lt.  0.001d0 * epsrel * fbm )then
c           reducing the step has no effect on B, abort
            s0a = s0a*1.05
            if (.false.) print*,' UF427 to be check'
            if( kcount.lt.15 )goto 200
            ifail = -42707
            return
          endif
          sbc     = msavb%arcl - msava%arcl
          s0c     = ( b0c / bbc ) * sbc
c
          if( s0a*s0c .lt. 0 )then
            ifail = -42705
            return   
          elseif( abs(s0c) .gt. stplst )then
c           but too far from the value
c-v
c-            write(2,*)' before but too far'
c-^
            s0a   = s0c * 1.05 + ( msava%arcl - mele(kpa)%arcl )
c
            if( kcount.lt.15 )goto 200
            ifail = -42702
            return
          elseif( abs(s0c) .lt. stpmin*2 )then
c           but too near!!!!
c-v
c-            write(2,*)' before but too near'
c-^
            s0a   = sign( abs(s0a) - 2.5 * stpmin, s0a)
            if( kcount.lt.15 )goto 200
            ifail = -42702
            return
          endif
        else
c         point c is after the value instead of before
          msavb=mele(kpb)
c         estimation of the correct location with a linear fit between c and msava
          bac     = msava%beg%b%dnrm - pb
          if( abs(bac) .lt.  0.001d0 * epsrel * fbm )then
c           reducing the step has no effect on B, abort
            s0a = s0a*0.95
            if (.false.) print*,' UF427 to be check'
            if( kcount.lt.15 )goto 200
            ifail = -42707
            return
          endif
          sac     = msava%arcl - msavb%arcl
          s0c     = ( b0c / bac ) * sac
c-v
c-          write(2,*)' after'
c-^
          if( s0a*s0c .ge. 0 )then
            ifail = -42706
            return   
          else
            s0a   = ( msavb%arcl - mele(kpa)%arcl ) + s0c * 1.05
            if( kcount.lt.15 )goto 200
            ifail = -42702
            return
          endif
        endif
c       set point c as new a
        kpa     = kgp
        msava   = mele(kpa)
        pb      = msavb%beg%b%dnrm
        b0a     = fbm - mele(kpa)%beg%b%dnrm
        bba     = pb -  mele(kpa)%beg%b%dnrm
        sba     = msavb%arcl - mele(kpa)%arcl
        kpb     = -1
        s0a     = ( b0a / bba ) * sba
      endif
c
c  now:  point a (kpa) is before and not too far
c        point b (msavb) is after
c  note that even the point b is at the value, the point has to be
c            re-evaluated
c
      s0a1      = 0
      s0a2      = 0
      if( (msava%beg%b%dnrm-fbm)*(msavb%beg%b%dnrm-fbm) .gt. 0 )then
        ifail   = -42799
        return
      endif
  300 continue
      s0a2      = s0a1
      s0a1      = s0a
c,      print'(a8,4f12.5)','uf427 a:',mele(kpa)%beg%coord%radius-re,
c,     :                   90.0-mele(kpa)%beg%coord%colat,
c,     :                   mele(kpa)%beg%coord%elong,
c,     :                   mele(kpa)%beg%b%dnrm
c,      print'(a8,5f12.5)','      b:',msavb%beg%coord%radius-re,
c,     :                   90.0-msavb%beg%coord%colat,
c,     :                   msavb%beg%coord%elong,
c,     :                   msavb%beg%b%dnrm, sba
c,      print'(a8,2f12.5)','     -> ',fbm,s0a
c     step to the new loc
      rc        = mele(kpa)%beg%rcurv
      dstep0    = MAX( stpmin, MIN( prop*rc, stepx*
     :                mele(kpa)%beg%coord%radius ) )
      scale     = s0a / dstep0
      kgp       = kpa
      kcount    = kcount + 1
      CALL UF423 (scale, rc, pb, palt, ifail)
      if( ifail.lt.0 )return
      mele(kgp)%beg%rcurv        = rc 
c
      b0c       = fbm - pb
      if( abs(b0c) .lt. 0.1 * epsrel * fbm )then
c       get it
        bi      = pb
        ai      = palt
        ifail   = kcount
        return
      elseif( b0c*b0a .lt. 0 )then
c       new point after
        msavb   = mele(kgp)
        sxa     = msava%arcl-mele(kpa)%arcl
      else
c       new point before
        msava   = mele(kgp)
        sxa     = msavb%arcl-mele(kpa)%arcl
      endif
        b0c     = fbm - msava%beg%b%dnrm
        bbc     = msavb%beg%b%dnrm - msava%beg%b%dnrm
        sbc     = msavb%arcl - msava%arcl
        s0c     = ( b0c / bbc ) * sbc
        s0a     = s0c + (msava%arcl-mele(kpa)%arcl)
        s0a     = 0.9 * s0a + 0.1 * sxa
c
      if( kcount.lt.30 )goto 300
      ifail = -42702
      return
C
      END
C----------------------------------------------------------------------
        SUBROUTINE UF428
     :          (falt, bi, ai, ifail)
C
C!    Interpolate to an altitude
C
      INCLUDE 'structure.h'
C      
      EXTERNAL UF423
C       
C     INTERFACE
C
        REAL*8        falt, bi, ai
        INTEGER*4     ifail        
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
C     VARIABLES
C
        REAL*8        s0a, sac, sbc, sba, s0c, bba, b0a, bbc, b0c,
     :                s0b, bac, sxa
        INTEGER*4     kpb, kpa, kcount
        REAL*8        dstep0, scale, rc, pb, palt
        TYPE(zseg) :: msavb, msava
C
C     CODE  based on UF427
C
      kpb       = kgp
      kpa       = kpb - 1
c
c
      kcount    = 0
C
C Estimation of the correct location with a linear fit
C 
      msavb     = mele(kpb)
      palt      = msavb%beg%coord%radius-re
      msava     = mele(kpa)
      b0a       = falt - (msava%beg%coord%radius-re)
      bba       = palt -  (msava%beg%coord%radius-re)
      sba       = msavb%arcl - msava%arcl
      if( abs(bba) .lt. 1.0e-3 )then
        ifail   = -42801
        return
      elseif( abs(sba) .lt. stpmin*0.1 )then
        ifail   = -42805
        return
      endif
      s0a       = ( b0a / bba ) * sba
      s0b       = s0a - sba
c
c check if the job has already been done
c
      if( ABS( bba - b0a ) .lt. 0.5d0 * epskm .and.
     :   abs( sba ) .lt. stplst )then
        ifail   = kcount
        bi      = msavb%beg%b%dnrm
        ai      = msavb%beg%coord%radius - re
        return
      endif
c
c check if the last step is not too large
c
  200 continue
      if( abs(s0a) .gt. stplst*0.9 )then
c       add a point
        rc      = mele(kpa)%beg%rcurv
        dstep0  = MAX( stpmin, MIN( prop*rc, stepx*
     :                mele(kpa)%beg%coord%radius ) )
        scale   = sign((abs(s0a)-stplst*0.5), s0a ) / dstep0
        kgp     = kpa
        kcount  = kcount + 1
        CALL UF423 (scale, rc, pb, palt, ifail)
        if( ifail.lt.0 )return
        mele(kpb)%beg%rcurv        = rc                             
c
        b0c     = falt - palt
c
        if( b0c*b0a .ge. 0 )then
c         point c is before the value
c
          msava=mele(kpb)
c         estimation of the correct location with a linear fit between c and msavb
          bbc     = msavb%beg%coord%radius-re - palt
          if( abs(bbc) .lt.  0.001d0 * epskm )then
c           reducing the step has no effect on B, abort
            s0a = s0a*1.05
            if( kcount.lt.15 )goto 200
            ifail = -42807
            return
          endif
          sbc     = msavb%arcl - msava%arcl
          s0c     = ( b0c / bbc ) * sbc
c
          if( s0a*s0c .lt. 0 )then
            ifail = -42805
            return   
          elseif( abs(s0c) .gt. stplst )then
c           but too far from the value
            s0a   = s0c * 1.05 + ( msava%arcl - mele(kpa)%arcl )
c
            if( kcount.lt.15 )goto 200
            ifail = -42802
            return
          elseif( abs(s0c) .lt. stpmin*2 )then
c           but too near!!!!
            s0a   = sign( abs(s0a) - 2.5 * stpmin, s0a)
            if( kcount.lt.15 )goto 200
            ifail = -42802
            return
          endif
        else
c         point c is after the value instead of before
          msavb=mele(kpb)
c         estimation of the correct location with a linear fit between c and msava
          bac     = msava%beg%coord%radius-re - palt
          if( abs(bac) .lt.  0.001d0 * epskm )then
c           reducing the step has no effect on alt, abort
            s0a = s0a*0.95
            if( kcount.lt.15 )goto 200
            ifail = -42807
            return
          endif
          sac     = msava%arcl - msavb%arcl
          s0c     = ( b0c / bac ) * sac
          if( s0a*s0c .ge. 0 )then
            ifail = -42806
            return   
          else
            s0a   = ( msavb%arcl - mele(kpa)%arcl ) + s0c * 1.05
            if( kcount.lt.15 )goto 200
            ifail = -42802
            return
          endif
        endif
c       set point c as new a
        kpa     = kgp
        palt      = msavb%beg%coord%radius-re
        b0a     = falt - (mele(kpa)%beg%coord%radius-re)
        bba     = palt -  (mele(kpa)%beg%coord%radius-re)
        sba     = msavb%arcl - mele(kpa)%arcl
        kpb     = -1
        s0a     = ( b0a / bba ) * sba
      endif
c
      msava     = mele(kpa)
c
  300 continue
c     step to the new loc
      rc        = mele(kpa)%beg%rcurv
      dstep0    = MAX( stpmin, MIN( prop*rc, stepx*
     :                mele(kpa)%beg%coord%radius ) )
      scale     = s0a / dstep0
      kgp       = kpa
      kcount    = kcount + 1
      CALL UF423 (scale, rc, pb, palt, ifail)
      if( ifail.lt.0 )return
      mele(kgp)%beg%rcurv        = rc 
c
      b0c       = falt - palt
      if( abs(b0c) .lt. 0.25 * epskm )then
c       get it
        bi      = palt
        ai      = palt
        ifail   = kcount
        return
      elseif( b0c*b0a .lt. 0 )then
c       new point after
        msavb   = mele(kgp)
        sxa     = msava%arcl-mele(kpa)%arcl
      else
c       new point before
        msava   = mele(kgp)
        sxa     = msavb%arcl-mele(kpa)%arcl
      endif
        b0c     = falt - (msava%beg%coord%radius-re)
        bbc     = msavb%beg%coord%radius - msava%beg%coord%radius
        sbc     = msavb%arcl - msava%arcl
        s0c     = ( b0c / bbc ) * sbc
        s0a     = s0c + (msava%arcl-mele(kpa)%arcl)
        s0a     = 0.9 * s0a + 0.1 * sxa
c
      if( kcount.lt.20 )goto 300
      ifail = -42802
      return
C
      END
C----------------------------------------------------------------------
      SUBROUTINE UF429
     :          (ji0, ji1)
C
C!    Transpose a field line segment 
C
      INCLUDE 'structure.h'
C
C     INTERFACE
C
        INTEGER*4     ji0, ji1    
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
C     VARIABLES
C
        TYPE(zseg) :: mstore
        INTEGER*4     i, k, j
C
C     CODE
C
      k = ( ji0 + ji1 - 1 ) / 2
      j = ji1
C      
      do i = ji0, k
C      
        mstore           = mele(i)       
        mele(i)          = mele(j)
        mele(j)          = mstore             
        j                = j - 1
C             
      enddo
C
      do i=ji0, ji1-1
        mele(i)%rkstp(2) = mele(i+1)%rkstp(2) + uk2 * 
     :              ( mele(i+1)%rkstp(1) - mele(i)%beg%b%dnrm )
        mele(i)%rkstp(3) = mele(i+1)%rkstp(3) + uk3 *
     :              ( mele(i+1)%rkstp(1) - mele(i)%beg%b%dnrm )
        mele(i)%rkstp(1) = mele(i+1)%beg%b%dnrm - 
     :              ( mele(i+1)%rkstp(1) - mele(i)%beg%b%dnrm )
        mele(i)%dtbnd    = mele(i+1)%dtbnd
      enddo
C      
      mele(ji1)%rkstp(2) = 0.0d0
      mele(ji1)%rkstp(3) = 0.0d0
      mele(ji1)%rkstp(1) = 0.0d0
      mele(ji1)%dtbnd    = 0.0d0
C
C
      END
C----------------------------------------------------------------------