# 1 "um520.f"
      SUBROUTINE UM520
     :          (kext, amjd, param, lbext, kunit, ifail)
C
C!    Select an external magnetic field model     
C_    er25
C
      INCLUDE 'structure.h'
cDEC$ IF DEFINED (_x86_)
cDEC$ ATTRIBUTES DLLEXPORT :: UM520
cDEC$ ENDIF
C
      EXTERNAL  UM521, UM522, UM523, UM524
C
C     INTERFACE
C
        REAL*8        amjd, param(20)
        INTEGER*4     kext, ifail, kunit
        CHARACTER*(*) lbext
C
      COMMON /UC140/  mint,mext,msun
C      
        TYPE(zimf) :: mint
        TYPE(zsun) :: msun
        TYPE(zemf) :: mext                   
        TYPE(zdat) :: adate
C          
      COMMON /UC160 
     :               /pi, deg, re, gmagmo, eclipt, geoid, uma
        REAL*8        pi, deg, re, gmagmo, eclipt, geoid(3), uma(30)
C
C     VARIABLES
C
        CHARACTER*32  lbl         
        REAL*8        dst, fkp, al, swd, swv, pdyn, sn
        REAL*8        bximf, byimf, bzimf, stdoff
        REAL*8        g1, g2, w1, w2, w3, w4, w5, w6

      COMMON /ALEXEEV/ adate

      adate%amjd = amjd
      CALL UT545(adate)
C
C     CODE
C 
      ifail = 0
C
C     Default value
C
      fkp   =   0.0d0
      dst   = -30.0d0
      al    =-100.0d0
      pdyn  =   3.0d0
      swd   =   6.0d0
      swv   = 300.0d0
      bximf =   6.0d0
      byimf =   0.0d0
      bzimf =   0.0d0
      stdoff=  10.0d0
      g1    =   6.0d0
      g2    =  10.0d0
      w1    =   4.0d0
      w2    =   3.0d0
      w3    =   4.0d0
      w4    =  10.0d0
      w5    =   7.0d0
      w6    =  20.0d0
C
C     Arguments change from version 1.06 to version 1.07
C
      if( param(1).gt.-999.0D0 ) fkp=param(1)
      if( param(2).gt.-999.0D0 ) dst=param(2)
      if( param(3).gt.-999.0D0 ) pdyn=param(3)
      if( param(4).gt.-999.0D0 ) swd=param(4)
      if( param(5).gt.-999.0D0 ) swv=param(5)
      if( param(6).gt.-999.0D0 ) bximf=param(6)
      if( param(7).gt.-999.0D0 ) byimf=param(7)
      if( param(8).gt.-999.0D0 ) bzimf=param(8)
      if( param(9).gt.-999.0D0 .and. kext.ne.6 )then
            stdoff=param(9)
      elseif( param(4).gt.0 .and. param(5).gt.0 )then 
            STDOFF = 98.0D0 / ((swd*swv**2)**(1.0D0/6.0D0))
      else
            stdoff = -10
      endif
      if( param(10).gt.-9999.0D0 ) al=param(10)
      if( param(11).gt.-999.0D0 ) g1=param(11)
      if( param(12).gt.-999.0D0 ) g2=param(12)
      if( param(13).gt.-999.0D0 ) w1=param(13)
      if( param(14).gt.-999.0D0 ) w2=param(14)
      if( param(15).gt.-999.0D0 ) w3=param(15)
      if( param(16).gt.-999.0D0 ) w4=param(16)
      if( param(17).gt.-999.0D0 ) w5=param(17)
      if( param(18).gt.-999.0D0 ) w6=param(18)
C
      if( kext.ge.0 .and. kext.le.11 )then
        mext%kouter  = kext
        CALL UM521 (lbl, fkp)
      else
        ifail  = -52001
        lbext  = 'ERROR'
        return                                 
      endif
C      
      lbext = lbl( 1 : MIN ( LEN(lbext), LEN(lbl) ) )
      mext%vdst   = dst
      mext%wdens  = swd
      mext%wvel   = swv
      mext%vkp    = fkp
      mext%pdyn   = pdyn
      mext%bximf  = bximf
      mext%byimf  = byimf
      mext%bzimf  = bzimf
      mext%stdoff = stdoff
      mext%val   = al
      mext%g1 = g1
      mext%g1 = g2
      mext%w1 = w1
      mext%w2 = w2
      mext%w3 = w3
      mext%w4 = w4
      mext%w5 = w5
      mext%w6 = w6
      mext%label  = lbl( 1 : MIN ( LEN(mext%label), LEN(lbl) ) )
C
      if( kunit.gt.0 )then
        write(kunit,1004,err=2000)
        if( kext.gt.0 .and. kext.lt.5 )then
          write(kunit,1142,err=2000) kext, lbl, fkp, mext%ikp
          if( fkp.lt.0 .or. fkp.gt.9 )write(kunit,1200,err=2000) 'Kp'
ccc          if( kext.eq.2 .or. kext.eq.3 )write(kunit,1205,err=2000)
        elseif( kext.eq.6 )then
          write(kunit,1143,err=2000) kext, lbl, dst, swd,
     :                    swv, abs(stdoff)
          if(dst.lt.-600.or.dst.gt.50 )write(kunit,1200,err=2000)  'Dst'
          if(swd.lt.1   .or.swd.gt.100 )write(kunit,1200,err=2000) 'SWD'
          if(swv.lt.100 .or.swv.gt.1200 )write(kunit,1200,err=2000)'SWV'
          if(stdoff.lt.2 )write(kunit,1201,err=2000)
          write(kunit,1250,err=2000) fkp
          if( fkp.lt.0 .or. fkp.gt.9 )write(kunit,1200,err=2000) 'Kp'
        elseif( kext.eq.7 .or. kext.eq.9 .or. kext.eq.10)then
          write(kunit,1144,err=2000) kext, lbl, pdyn, dst, byimf, bzimf
          if( kext.eq.9 ) write(kunit,1146,err=2000) g1, g2
          if( kext.eq.10 ) write(kunit,1147,err=2000) w1,w2,w3,w4,w5,w6
          write(kunit,1250,err=2000) fkp
          if(pdyn.lt.0.5.or.pdyn.gt.10 )write(kunit,1200,err=2000) 'SWP'
          if(dst.lt.-600.or.dst.gt.50 )write(kunit,1200,err=2000)  'Dst'
          if(byimf.lt.-50.or.byimf.gt.50)write(kunit,1200,err=2000)'By'
          if(bzimf.lt.-50.or.bzimf.gt.50)write(kunit,1200,err=2000)'Bz'
          if( kext.eq.9 )then
            if(g1.lt.0.0.or.g1.gt.10.0)write(kunit,1200,err=2000) 'G1'
            if(g2.lt.0.0.or.g2.gt.10.0)write(kunit,1200,err=2000) 'G2'
          endif
          if( kext.eq.10 )then
            if(w1.lt.0.0.or.w1.gt.20.0)write(kunit,1200,err=2000) 'W1'
            if(w2.lt.0.0.or.w2.gt.10.0)write(kunit,1200,err=2000) 'W2'
            if(w3.lt.0.0.or.w3.gt.20.0)write(kunit,1200,err=2000) 'W3'
            if(w4.lt.0.0.or.w4.gt.100.0)write(kunit,1200,err=2000) 'W4'
            if(w5.lt.0.0.or.w5.gt.100.0)write(kunit,1200,err=2000) 'W5'
            if(w6.lt.0.0.or.w6.gt.100.0)write(kunit,1200,err=2000) 'W6'
          endif
          if( fkp.lt.0 .or. fkp.gt.9 )write(kunit,1200,err=2000) 'Kp'
        elseif( kext.eq.8 )then
          write(kunit,1145,err=2000) kext, lbl, fkp, dst, pdyn, bzimf
          if( fkp.lt.0 .or. fkp.gt.9 )write(kunit,1200,err=2000) 'Kp'
          if(dst.lt.-600.or.dst.gt.50 )write(kunit,1200,err=2000)  'Dst'
          if(pdyn.lt.0.5.or.pdyn.gt.10 )write(kunit,1200,err=2000) 'SWP'
          if(bzimf.lt.-50.or.bzimf.gt.50)write(kunit,1200,err=2000)'Bz'
        elseif( kext.eq.11 )then
          write(kunit,1148,err=2000) kext, lbl, dst, al, swd, swv, bzimf
          write(kunit,1250,err=2000) fkp
          if(dst.lt.-600.or.dst.gt.50 )write(kunit,1200,err=2000)  'Dst'
          if(al.lt.-3000.or.al.gt.100 )write(kunit,1200,err=2000)  'AL'
          if(swd.lt.1   .or.swd.gt.100 )write(kunit,1200,err=2000) 'SWD'
          if(swv.lt.100 .or.swv.gt.1200 )write(kunit,1200,err=2000)'SWV'
          if(bzimf.lt.-50.or.bzimf.gt.50)write(kunit,1200,err=2000)'Bz'
          if( fkp.lt.0 .or. fkp.gt.9 )write(kunit,1200,err=2000) 'Kp'
        else
          write(kunit,1141,err=2000) kext, lbl
          write(kunit,1250,err=2000) fkp
          if( fkp.lt.0 .or. fkp.gt.9 )write(kunit,1200,err=2000) 'Kp'
        endif
      endif
C
      CALL UM522 (amjd,kunit)
C
      if (kext.eq.2 .or. kext.eq.3 .or. kext.eq.4 .or. kext.eq.7 .or.
     &    kext .eq. 9 .or. kext .eq. 10 .or. kext .eq. 11) then
        CALL UM523 (kunit)
      else
        CALL UM524 (kunit)
      endif
C
c
c
      if (kext.eq.8) then
c     om97
c
        sn = sin(mext%tilt*deg)
        CALL set_a (dst, pdyn, fkp, bzimf, sn)
c
      endif
c
c
c
      return
 2000 ifail=-52002
C
 1004 FORMAT (/' --- External magnetic field model ---')
 1141 FORMAT (/6x,'Model (',i2,'): ',a32)
 1142 FORMAT (/6x,'Model (',i1,'): ',a32,/11x,
     :            'Planetary activity index Kp =',f7.1,17x,'[',i1,']')
 1143 FORMAT (/6x,'Model (',i1,'): ',a32,/14x,
     :              'Storm activity index Dst =',f7.1,5x,' nT',/20x,
     :                'Solar wind density =',f7.1,5x,' 1/cm^3',/19x,
     :                      'Solar wind velocity =',f7.1,5x,' km/s',
     :                     /21x,'Standoff distance =',f7.1,5x,' Re')
 1144 FORMAT (/6x,'Model (',i2,'): ',a32,/19x,
     :            'Solar wind pressure =',f8.2,4x,' nP',/14x,
     :            'Storm activity index Dst =',f7.1,5x,' nT',/6x,
     :            'Interplanetary magnetic field, Y =',f10.4,2x,' nT',
     :                                      /37x,'Z =',f10.4,2x,' nT')
 1145 FORMAT (/6x,'Model (',i1,'): ',a32,/11x,
     :     'Planetary activity index Kp =',f7.1,/14x,
     :              'Storm activity index Dst =',f7.1,5x,' nT',
     :                   /19x,'Solar wind pressure =',f8.2,4x,' nP',
     :       /6x,'Interplanetary magnetic field, Z =',f10.4,2x,' nT')
 1146 FORMAT (36x,'G1 =',f7.1,5x,' km/s',/,
     :        36x,'G2 =',f7.1,5x,' km/s nT')
 1147 FORMAT (36x,'W1 =',f7.1,5x,' 1/(5 min)',/,
     :        36x,'W2 =',f7.1,5x,' 1/(5 min)',/,
     :        36x,'W3 =',f7.1,5x,' 1/(5 min)',/,
     :        36x,'W4 =',f7.1,5x,' 1/(5 min)',/,
     :        36x,'W5 =',f7.1,5x,' 1/(5 min)',/,
     :        36x,'W6 =',f7.1,5x,' 1/(5 min)')
 1148 FORMAT (/6x,'Model (',i2,'): ',a32,/14x,
     :            'Storm activity index Dst =',f7.1,5x,' nT',/30x,
     :            'AL index =',f7.1,5x,' nT',/20x,
     :            'Solar wind density =',f7.1,5x,' 1/cm^3',/19x,
     :            'Solar wind velocity =',f7.1,5x,' km/s',/6x,
     :            'Interplanetary magnetic field, Z =',f10.4,2x,' nT')
 1200 FORMAT (6x,'WARNING: the parameter ',a,' is out of range')
 1201 FORMAT (6x,'WARNING: the standoff distance is out of range')
 1205 FORMAT (6x,'WARNING: Dr. Tsyganenko recommends not using ',
     :                    'the obsolete 1987 versions')
 1250 FORMAT (11x,'Planetary activity index Kp =',f7.1,5x,
     :        ' (magnetopause check)')
      END
C----------------------------------------------------------------------
      SUBROUTINE UM521
     :          (lbl, fkp)
C
C!    Evalue the ground disturbances level based on Kp
C
C     REFERENCES
C     Subroutine blxtram of blxtra.for
C
      INCLUDE 'structure.h'
C
C     INTERFACE           
C    
        REAL*8        fkp
        CHARACTER*32  lbl
C
      COMMON /UC140/  mint, mext, msun
C
        TYPE(zimf) :: mint
        TYPE(zsun) :: msun
        TYPE(zemf) :: mext      
C
C     VARIABLES
C
        REAL*8        tab_outer(0:11,9)
        INTEGER*4     i, j, k
        CHARACTER*32  lextmod(12)
C
      DATA lextmod
     :       /'None                            ',
     :        'Mead-Fairfield 1975             ',
     :        'Tsyganenko 1987 short           ',
     :        'Tsyganenko 1987 long            ',
     :        'Tsyganenko 1989                 ',
     :        'Olson-Pfitzer 1977              ',
     :        'Olson-Pfitzer dynamic           ',
     :        'Tsyganenko 1996                 ',
     :        'Ostapenko-Maltsev 1997          ',
     :        'Tsyganenko 2001                 ',
     :        'Tsyganenko 2004                 ',
     :        'Paraboloid model                '/
C
      DATA   ((tab_outer(j,i),i=1,9),j=0,11) 
     : /0.0d0, 0.5d0, 2.0d0, 3.0d0, 9.999d2, 0.0d0, 0.0d0, 0.0d0, 0.0d0,
     :  0.0d0, 0.5d0, 2.0d0, 3.0d0, 9.999d2, 0.0d0, 0.0d0, 0.0d0, 0.0d0,
     :  0.0d0, 0.5d0,1.17d0,1.833d0,2.5d0, 3.5d0, 4.5d0, 5.17d0,9.999d2,
     :  0.0d0, 0.5d0, 1.5d0, 2.5d0, 3.5d0, 4.5d0, 9.999d2, 0.0d0, 0.0d0,
     :  0.0d0, 0.5d0, 1.5d0, 2.5d0, 3.5d0, 4.5d0, 5.5d0, 9.999d2, 0.0d0,
     :  0.0d0, 0.5d0, 1.5d0, 2.5d0, 3.5d0, 4.5d0, 9.999d2, 0.0d0, 0.0d0,
     :  0.0d0, 0.5d0, 1.5d0, 2.5d0, 3.5d0, 4.5d0, 9.999d2, 0.0d0, 0.0d0,
     :  0.0d0, 0.5d0, 1.5d0, 2.5d0, 3.5d0, 4.5d0, 9.999d2, 0.0d0, 0.0d0,
     :  0.0d0, 9.999d2, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0,
     :  0.0d0, 9.999d2, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0,
     :  0.0d0, 9.999d2, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0,
     :  0.0d0, 9.999d2, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/
C
C     CODE
C
      k   = mext%kouter
      if ( k .lt. 0 .or. k .gt. 11 ) k=0
      lbl = lextmod(k+1)
C     
      mext%ikp = 1
      do while ( fkp .ge. tab_outer(k,mext%ikp + 1) ) 
        mext%ikp = mext%ikp + 1
      end do
C
C
      END
C----------------------------------------------------------------------
      SUBROUTINE UM522
     :          (amjd,kunit)
C
C!    Compute the position of the Sun
C
C     REFERENCES
C     Subroutines of shell.for and transfos.for
C
      INCLUDE 'structure.h'
C
C     INTERFACE           
C    
        REAL*8        amjd
        INTEGER*4     kunit
C
      COMMON /UC140/  mint,mext,msun
C
        TYPE(zimf) :: mint
        TYPE(zsun) :: msun
        TYPE(zemf) :: mext   
C
      COMMON /UC160 
     :               /pi, deg, re, gmagmo, eclipt, geoid, uma
        REAL*8        pi, deg, re, gmagmo, eclipt, geoid(3), uma(30)
C                                                
C     VARIABLES
C
        TYPE(zxyz) :: msol
        REAL*8        ad2000, x7, x8, as
        REAL*8        resjd, t0, sth
        REAL*8        sst, cst
        INTEGER*4     intjd
        TYPE(zdat) :: mdate
        CHARACTER*3   month(12)
C
        DATA month/'Jan','Feb','Mar','Apr','May','Jun','Jul',
     :             'Aug','Sep','Oct','Nov','Dec'/
C
C     CODE
C
C  CALCULATE: the normalized direction of the Sun in GEO and
C             the Greenwich hour angle of equinox (cf GMST)
C  REF: VANFLANDERN, PULKKINEN: LOW-PRECISION FORMULAE FOR PLANETARY
C          POSITIONS. ASTROPH. J. SUPPL. 41 PP. 391-411, NOV. 1979
C  Remark: MJD(M.A.Hapgood) = AMJD + 33282
C
C  a/ Convert to julian day of anno 2000
C
      ad2000     = amjd - 18262.5d0
C      
C  b/ Sun position
C
      x7         = MOD(0.779072d0+0.273790931d-2*ad2000, 1.0d0)
     :             * 2.0d0 * pi  
      x8         = MOD(0.993126d0+0.273777850d-2*ad2000, 1.0d0)
     :             * 2.0d0 * pi
      as         = x7 + SIN(x8) * (335.0d-4+698.0d-6*COS(x8))
      msol%x     = COS(as)
      msol%y     = COS(eclipt*deg) * SIN(as)
      msol%z     = SIN(eclipt*deg) * SIN(as)
C
C  c/ Time in Julian centuries from 12h UT on 1 January 2000 to the
C     previous midnight
C
      intjd      = amjd
      resjd      = amjd - intjd
      msun%utdeg = resjd * 360.0d0
      t0         = (intjd - 18262.5d0) / 36525.0d0
C
C  d/ Greenwich hour angle of equinox
C
      sth        = ((-1.722d-9*t0+2.586222d-5)*t0+2.400051336907d3)
     :             * t0 + 6.69737455833d0
      if (sth .lt. 0.0d0) then
        sth = (IDINT(-sth/24.0d0)+1.0d0) * 24.0d0 + sth
      else
        sth = sth - IDINT(sth/24.0d0) * 24.0d0
      endif
      sth   = sth + resjd * 2.406570721d1
C
      msun%gha   = sth * 15.0d0
C
C  e/ Convert GEI coord. of the Sun to GEO coord.
C
      sst        = SIN( msun%gha * deg )
      cst        = COS( msun%gha * deg )
      msun%dir%x = cst * msol%x + sst * msol%y
      msun%dir%y = -sst * msol%x + cst * msol%y
      msun%dir%z = msol%z
C
      if ( kunit .gt. 0 ) then
        mdate%amjd = amjd
        call ut545 (mdate)
        write(kunit,1142,err=2000) mdate%iday, month(mdate%imonth),
     :                mdate%iyear,msun%utdeg/15.0d0, msun%gha,
     :                    msun%dir%x, msun%dir%y, msun%dir%z
      endif
C
      return
 2000 write(*,*)'%UNILIB: UM522, error on output device'
      stop        
C
 1142 FORMAT (/,34x,'Date = ',i3.2,'-',a3,'-',i4.4,
     :         /24x,'Universal Time =',f8.2,4x,' hrs',/12x,
     :             'Greenwich angle of equinox =',f8.2,4x,' deg',/11x,
     :               'GEO direction of the Sun, X =',f10.4,/37x,'Y =',
     :                                  f10.4,2x,/37x,'Z =',f10.4,2x) 
      END
C----------------------------------------------------------------------
      SUBROUTINE UM523
     :                (kunit)
C
C!    Transformation matrix from GEO to GSM
C
C     REFERENCES
C     Subroutine TRANS1 of transfos.for
C
      INCLUDE 'structure.h'
C
C     INTERFACE
C
        INTEGER*4     kunit               
C
      COMMON /UC140/  mint, mext, msun
C
        TYPE(zimf) :: mint
        TYPE(zsun) :: msun
        TYPE(zemf) :: mext     
C
      COMMON /UC160 
     :               /pi, deg, re, gmagmo, eclipt, geoid, uma
        REAL*8        pi, deg, re, gmagmo, eclipt, geoid(3), uma(30)
C                                                
C     VARIABLES
C
        REAL*8        slt, clt, sln, cln, dx, dy, dz
        REAL*8        aly, amy, any, amag
        REAL*8        dgsmx
        REAL*8        t(3,3)
        INTEGER*4     i, j
C
C     CODE
C
C  a/ Calculate cartesian GEO coordinates of boreal dipole pole
C
      slt     = SIN( mint%colat * deg )
      clt     = COS( mint%colat * deg )
      sln     = SIN( mint%elong * deg )
      cln     = COS( mint%elong * deg )
      dx      = slt * cln
      dy      = slt * sln
      dz      = clt       
C                
C  b/ X-Axis of GSM coordinate system coincides with sun direction
C
      t (1,1) = msun%dir%x
      t (1,2) = msun%dir%y
      t (1,3) = msun%dir%z
C      
C  c/ Calculate GEO coordinates of y-axis of GSM coordinate system
C
      aly     = dy * msun%dir%z - dz * msun%dir%y
      amy     = dz * msun%dir%x - dx * msun%dir%z
      any     = dx * msun%dir%y - dy * msun%dir%x
      amag    = 1.0d0 / SQRT(aly*aly+amy*amy+any*any)
      t (2,1) = aly * amag
      t (2,2) = amy * amag
      t (2,3) = any * amag
C      
C  d/ Calculate GEO coordinates of z-axis of GSM coordinate system
C
      t(3,1)  = t(1,2) * t(2,3) - t(1,3) * t(2,2)
      t(3,2)  = t(1,3) * t(2,1) - t(1,1) * t(2,3)
      t(3,3)  = t(1,1) * t(2,2) - t(1,2) * t(2,1)
C
C  e/ Dipole tilt angle is GSM colatitude of dipole axis
C     invert sign of tilt when dipole axis is tilted away from the sun,
C     i.e. when GSM x coordinate of dipole axis is negative
C
      mext%tilt = ACOS( t(3,1) * dx + t(3,2) * dy + t(3,3) * dz ) 
     :                              / deg
      dgsmx     = t(1,1) * dx + t(1,2) * dy + t(1,3) * dz
      if (dgsmx .lt. 0.0d0) mext%tilt = -mext%tilt
C
      do i=1,3
        do j=1,3
          mext%trans(j,i) = t(j,i)
        enddo
      enddo
C
      if( kunit.gt.0 )then
        write(kunit,1143,err=2000) mext%tilt
      endif
c
      return
 2000 write(*,*)'%UNILIB: UM522, error on output device'
      stop           
C
 1143 FORMAT(21x,'Dipole tilt angle =',f8.2,4x,' deg') 
      END
C----------------------------------------------------------------------
      SUBROUTINE UM524
     :                (kunit)
C
C!    Transformation matrix from GEO to SM
C
C     REFERENCES
C     Subroutine TRANS2 of transfos.for
C
      INCLUDE 'structure.h'
C
C     INTERFACE           
C    
        INTEGER*4     kunit               
C
      COMMON /UC140/  mint, mext, msun
C      
        TYPE(zemf) :: mext
        TYPE(zimf) :: mint
        TYPE(zsun) :: msun
C
      COMMON /UC160 
     :               /pi, deg, re, gmagmo, eclipt, geoid, uma
        REAL*8        pi, deg, re, gmagmo, eclipt, geoid(3), uma(30)
C                                                
C     VARIABLES
C
        REAL*8        slt, clt, sln, cln
        REAL*8        aly, amy, any, amag
        REAL*8        t(3,3)
        INTEGER*4     i,j
C
C     CODE
C
C  a/ boreal dipole pole is z-axis of SM coordinate system
C
      slt    = SIN(mint%colat * deg)
      clt    = COS(mint%colat * deg)
      sln    = SIN(mint%elong * deg)
      cln    = COS(mint%elong * deg)
      t(3,1) = slt * cln
      t(3,2) = slt * sln
      t(3,3) = clt
C
C  b/ calculate GEO coordinates of y-axis of SM coordinate system
C
      aly    = t(3,2) * msun%dir%z - t(3,3) * msun%dir%y
      amy    = t(3,3) * msun%dir%x - t(3,1) * msun%dir%z
      any    = t(3,1) * msun%dir%y - t(3,2) * msun%dir%x
      amag   = 1.0d0 / SQRT(aly*aly+amy*amy+any*any)  
      t(2,1) = aly * amag
      t(2,2) = amy * amag
      t(2,3) = any * amag
C
C  c/ calculate GEO coordinates of x-axis of SM coordinate system
C
      t(1,1) = t(2,2) * t(3,3) - t(2,3) * t(3,2)
      t(1,2) = t(2,3) * t(3,1) - t(2,1) * t(3,3)
      t(1,3) = t(2,1) * t(3,2) - t(2,2) * t(3,1)
C      
C  d/ dipole tilt angle is complement of colatitude of sun direction
C
      mext%tilt = ACOS( t(3,1) * msun%dir%x +
     :                  t(3,2) * msun%dir%y +
     :                  t(3,3) * msun%dir%z )
     :            /deg
      mext%tilt = 90.0d0 - mext%tilt
C
      do i=1,3
        do j=1,3
          mext%trans(j,i) = t(j,i)
        enddo
      enddo
C        
      if( kunit.gt.0 )then
        write(kunit,1144,err=2000) mext%tilt
      endif
c
      return
 2000 write(*,*)'%UNILIB: UM522, error on output device'
      stop           
C
C
 1144 FORMAT(21x,'Dipole tilt angle =',f8.2,4x,' deg') 

      END      
C----------------------------------------------------------------------