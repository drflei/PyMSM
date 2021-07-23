# 1 "ua610.f"
      subroutine ua610
     :           (kmodel, mdate, rzss, f107a, f107, fkpx,
     :            ap0, ap3, ap6, ap9, ap24, ap48, lbatm, kunit,ifail)

c
      include 'structure.h'
cDEC$ IF DEFINED (_x86_)
cDEC$ ATTRIBUTES DLLEXPORT :: UA610
cDEC$ ENDIF
C            
      integer*4     kmodel, kunit, ifail
      TYPE(zdat) :: mdate
      real*8        rzss, f107a, f107, fkpx,
     :              ap0, ap3, ap6, ap9, ap24, ap48
      character*(*) lbatm
C
C
      COMMON /UC150/matm, ntspec, nnspec, kspec, kflag
C 
            TYPE(zatm) :: matm
            INTEGER*4     ntspec, nnspec, kspec(30)
            INTEGER*4     kflag(50)
C                           
      TYPE(zdat) :: mjan1st
      character*32  label, lion, lneut
      REAL*8        x
C
C
      ifail=0
C
      mjan1st=mdate
      mjan1st%iday=1
      mjan1st%imonth=1
      call ut540(mjan1st)
      matm%kyear=mdate%iyear
      matm%kday=mdate%amjd-mjan1st%amjd+1.000005d0
      matm%ut=15.0d0*(mdate%ihour+(mdate%imin+mdate%secs/60.0d0)/60.0d0)
C
      if( matm%kyear .lt. 1900 .or. matm%kyear .gt. 2500 .or.
     :    matm%kday .le. 0 .or. matm%kday .gt. 366 )then
        ifail=-61001
c       incorrect date in argument mdate
        return
      endif
C
      matm%rzss=rzss
      matm%f107a=f107a
      matm%f107=f107
      matm%apind(2)=ap0
      matm%apind(3)=ap3    
      matm%apind(4)=ap6
      matm%apind(5)=ap9
      matm%apind(6)=ap24
      matm%apind(7)=ap48
      matm%apind(1)=((ap0+ap3+ap6+ap9)*0.25d0+ap24*0.5d0)/1.5d0
      matm%fkpx=fkpx
c
c
      ntspec=0
      nnspec=0
      matm%katm=0
      matm%kion=0
      if( kmodel.eq.0 )then
        label='No atmospheric model            '
      elseif( kmodel.eq.1 )then
        label='MSISE-90 neutral atmosphere     '
        matm%katm=1
        nnspec=7
        kspec(1)= 1
        kspec(2)= 2
        kspec(3)= 4
        kspec(4)= 8
        kspec(5)= 9
        kspec(6)=16
        kspec(7)=17
      elseif( kmodel.eq.2 )then
        label='MSISE-90 and IRI-90 models      '
        matm%katm=1
        matm%kion=1
        nnspec=7
        kspec(1)= 1
        kspec(2)= 2
        kspec(3)= 4
        kspec(4)= 8
        kspec(5)= 9
        kspec(6)=16
        kspec(7)=17
        ntspec=nnspec+7
        kspec(8)=  19
        kspec(9)=  20
        kspec(10)= 21
        kspec(11)= 23
        kspec(12)= 24
        kspec(13)= 27
        kspec(14)= 28
      elseif( kmodel.eq.3 )then
        label='Anderson & Francis model        '
        matm%katm=2
        matm%kion=2
        nnspec=6
        kspec(1)= 1
        kspec(2)= 2
        kspec(3)= 8
        kspec(4)= 9
        kspec(5)=16
        kspec(6)=17
        ntspec=nnspec+6
        kspec(7)=  19
        kspec(8)=  20
        kspec(9)=  21
        kspec(10)= 24
        kspec(11)= 27
        kspec(12)= 28
      elseif( kmodel.eq.4 )then
        label='simple atmospheric model        '
        matm%katm=3
      elseif( kmodel.eq.5 )then
        label='Pfitzer atmospheric model       '
        matm%katm=4
        matm%kion=3
      else
        ifail=-61002
c       invalid  model number
        return
      endif
      ntspec=max(ntspec,nnspec)
c
        lbatm = label( 1:MIN( LEN(lbatm),LEN(label) ) )
c
      if( kunit .gt. 0 )then
        write(kunit,1000,err=2000)kmodel,label
        call ua612(kspec(1),nnspec,lneut,ifail)
        if( ifail.lt.0 )return
        call ua612(kspec(1+nnspec),ntspec-nnspec,lion,ifail)
        if( ifail.lt.0 )return
        x=matm%kyear+(matm%kday-1)/365.0d0
        write(kunit,1010,err=2000)x,matm%ut,matm%rzss,matm%f107,
     :       matm%f107a,
     :       matm%fkpx,matm%apind,lneut,lion
      endif
c
      return
 2000 ifail=-61003
c
 1000 format(/' --- Atmospheric model ---',//6x,' Model (',i1,'): ',a32)
 1010 format(21x,'Calculation epoch =',f7.1,6x,'year',
     :      /24x,'Universal Time =',f7.2,5x,'deg',
     :      /22x,'Sun spots number =',f7.0,
     : /8x,'Daily radio flux F10.7 (day-1) =',f7.1,6x,'10-22 W m-2 Hz',
     : /9x,'81-Day avrg. radio flux F10.7 =',f7.1,6x,'10-22 W m-2 Hz',
     : /11x,'Max. value of the 3-hour Kp =',f7.1,
     : /18x,'Different Ap indices =',7f4.0,      
     :      /23x,'Neutral species =',a32,/23x,'Ionized species =',a32)
c
      end
C----------------------------------------------------------------------
      subroutine ua612(ksp,n,lbl,ifail)
c
      include 'structure.h'
C            
C
      integer*4 ksp(30),n,ifail
      character*(*) lbl
C
      character*4 ltbl(30)
      integer*4   i, klen
C
      data ltbl/ '   H','  He','  Ne','  Ar','  Kr','  Xe','  H2',
     :           '  N2','  O2','  CO',' CO2',' NH3',' CH4','####',
     :           '####','  O.','  N.','####','  e-','  p+',' He+',
     :           '  C+','  N+','  O+',' Ne+',' N2+',' NO+',' O2+',
     :           '####','####'/
C
      klen=len(lbl)/4
      klen=min(klen,n)
      ifail=0
C
      lbl=' '
C
      do i=1, klen
        if( ksp(i).le.0 .or. ksp(i).gt.30 )then
          ifail=-61201
          return
        endif
        lbl(i*4-3:i*4)=ltbl(ksp(i))
      enddo
C
      if( n .ne. klen )then
        lbl(len(lbl)-4:len(lbl))=' ...'
        ifail=n-klen
      endif
C
      end
C----------------------------------------------------------------------
