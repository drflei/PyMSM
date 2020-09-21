cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc          Trapped Radiation ENvironment Development (TREND) - 4             cc
cc      ESTEC/Contract No. 11711/95/NL/JG - CCN 1 to Work Order No. 3         cc
cc                                                                            cc
cc FORTRAN interface routines for UNIRAD tabular outputs                      cc
cc                                                                            cc
cc Author: Kruglanski M., BIRA/IASB                                           cc
cc                                                                            cc
cc Content:                                                                   cc
cc   subroutine read_header( lun, header, nsize, kskip, nmulti)               cc
cc   subroutine read_header_plus ( lun, header, nsize, kskip, nmulti, metavar,cc
cc                                 futur)                                     cc
cc   subroutine metavar    ( knt, nmeta, adummy, flag)  --- user supplied --- cc
cc   subroutine futur      ( knt, nfutur, adummy, flag) --- user supplied --- cc
cc   subroutine to_string  ( kval, string)                                    cc
cc   subroutine get_listof ( header, ktype, table, ksize)                     cc
cc   subroutine get_text   ( header, number, string)                          cc
cc   subroutine select_var ( header, name, number, unit, title, kolumn,       cc
cc  :                        ksize)                                           cc
cc   subroutine select_one ( header, name, unit, title, ndim, index,          cc
cc  :                        ksize, kolumn, nkol)                             cc
cc   subroutine info_var   ( header, name, unit, title, kolumn, ndim)         cc
cc   subroutine get_meta   ( header, meta, rvalue, ksize)                     cc
cc   subroutine get_meta_int( header, meta, kvalue, ksize)                    cc
cc   subroutine get_meta_string( header, meta, svalue, ksize)                 cc
cc   subroutine get_var    ( header, kolumn, ksize, kskip, number, value,     cc
cc  :                        ndim, footer)                                    cc
cc                                                                            cc
cc Documentation:                                                             cc
cc   http://www.magnet.oma.be/trend4/doc/soft/ssd/format/index.html           cc
cc                                                                            cc
cc History:                                                                   cc
cc   August 10, 1998    First release                                         cc
cc   December 17, 1998  MK: addition of select_one                            cc
cc   June 21, 1999      MK: get_meta can not be used on PC to retrieve        cc
cc                          a string metavariable. Addition of a new entry    cc
cc                          (get_meta_string) with the same syntax            cc
cc   September 2001     MK: Addition of the entry read_header_plus to access  cc
cc                          the 'future use' header block                     cc
cc   November 8, 2001   MK: Addition of the entry get_meta_int to retrieve    cc
cc                          integer metadata                                  cc
cc   March 21, 2002     DH: Addition of the metavar function in               cc
cc                          read_header_plus to store the metavariable lines  cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine read_header(lun,header,nsize,kskip,nmulti)
c
      implicit     none

      INCLUDE 'interface.h'

      integer*4    lun, nsize, kskip
      integer*4    header(*), nmulti
c
      integer*4    nlim
      parameter  ( nlim= 5000 )
      character    achaine*20000
      real*4       rchaine(nlim)
      integer*4    ichaine(nlim)
      equivalence  (achaine,rchaine)
      equivalence  (achaine,ichaine)
      integer*4    ispace
c
      character    adummy*20000
      integer*4    i, knt, nl, index, ici, ityp, kcol, k
      integer*4    nheader,ntext,nmeta,nfutur,nvar, j, flag
      character*80 atask
      external     metavar,futur
c
      flag      = 0
      goto 10
c
      entry read_header_plus(lun,header,nsize,kskip,nmulti,metavar,
     :                       futur)
      flag      = 1
c
   10 continue
      achaine   = '    '
      ispace    = ichaine(1)
c
      header(1) = lun
      header(2) = 0
      header(3) = 0
      header(4) = 0
      header(5) = 0
      header(6) = 0
      index     = 7
c
      nl  = 0
      knt = kskip
  100 continue
        atask = 'reading the first header line'
        read (lun,*,err=999,end=999) 
     &          achaine,nheader,(ichaine(i),i=2,7),nmulti
        nl = nl + 1
        if( achaine(1:4) .ne. '*   ' )goto 999
        if( knt .le. 0 )goto 120
        atask = 'skipping a whole block'
        do i=2,nheader
          read (lun,1010,err=999,end=999) adummy
          nl = nl + 1
        enddo
  110   continue
          read (lun,1010,err=999,end=999) adummy
          nl = nl + 1
        if( adummy(1:1) .ne. '''' )goto 110
C read 'End of block line' was not taken into account  (DH, June 2007)
c        if( adummy(1:3) .ne. '''*''' )goto 999
        knt = knt - 1
      goto 100
  120 continue
c
      header(2) = ichaine(7)
      header(3) = ichaine(6)
      ntext     = ichaine(2)
      nmeta     = ichaine(3)
      nfutur    = ichaine(4)
      nvar      = ichaine(5)
c
      atask = 'reading the text lines'
      ici = 4
      do knt=1,ntext
        read (lun,*,err=999,end=999) achaine
        nl = nl + 1
        do k=nlim,1,-1
          if( ichaine(k) .ne. ispace )goto 150
        enddo
        k = 1
  150   header(ici)   = index
        ici           = index
        header(ici)   = 0
        header(ici+1) = k
        do i=1,k
          header(ici+1+i) = ichaine(i)
        enddo
        index = ici+2+k
        if( index .gt. nsize )goto 998
      enddo
c
      atask = 'reading the metavariables'
      ici = 5
      do knt=1,nmeta
        read (lun,1010,err=999,end=999) adummy
        if( flag.ne.0 )then
          atask = 'calling user subroutine'
          call metavar(knt,nmeta,adummy,flag)
          if( flag.lt.0 )goto 999
          atask = 'reading the metavariables'
        endif
        nl = nl + 1
        read (adummy,*,err=999,end=999) achaine, ityp
        do k=nlim,1,-1
          if( ichaine(k) .ne. ispace )goto 200
        enddo
        k = 1
  200   header(ici)   = index
        ici           = index
        header(ici)   = 0
        header(ici+1) = k
        do i=1,k
          header(ici+1+i) = ichaine(i)
        enddo
        index = ici+2+k
        header(index) = ityp
        index = index+1
c
        if( ityp .lt. 0 )then
          do j=1,-ityp
            read (adummy,*,err=999,end=999) 
     &                 achaine,ityp,(achaine,k=1,j)
            do k=nlim,1,-1
              if( ichaine(k) .ne. ispace )goto 220
            enddo
            k = 1
  220       header(index)   = k
            do i=1,k
              header(index+i) = ichaine(i)
            enddo
            index = index+1+k
            if( index .gt. nsize )goto 998
          enddo
        elseif( ityp .gt. 0 )then
          read (adummy,*,err=999,end=999) 
     &               achaine,ityp,(rchaine(i),i=1,ityp)
          do i=1,ityp
            header(index+i-1) = ichaine(i)
          enddo
          index = index+ityp
        else
          atask = 'storing metavariable ' // achaine
          goto 999
        endif
        if( index .gt. nsize )goto 998
      enddo
c
      atask = 'skipping not used lines'
      do knt=1,nfutur
        read (lun,1010,err=999,end=999) adummy
        if( flag.ne.0 )then
          atask = 'calling user subroutine'
          call futur(knt,nfutur,adummy,flag)
          if( flag.lt.0 )goto 999
          atask = 'skipping not used lines'
        endif
        nl = nl + 1
      enddo
c
      atask = 'reading the variable description'
      ici = 6
      kcol = 1
      do knt=1,nvar
        read (lun,1010,err=999,end=999) adummy
        nl = nl + 1
        read (adummy,*,err=999,end=999) achaine
        do k=nlim,1,-1
          if( ichaine(k) .ne. ispace )goto 250
        enddo
        k = 1
  250   header(ici)   = index
        ici           = index
        header(ici)   = 0
        header(ici+1) = k
        do i=1,k
          header(ici+1+i) = ichaine(i)
        enddo
        index = ici+2+k
        read (adummy,*,err=999,end=999) achaine,achaine,ityp
        do k=nlim,1,-1
          if( ichaine(k) .ne. ispace )goto 260
        enddo
        k = 1
  260   header(index) = ityp
        header(index+1) = kcol
        kcol = kcol + ityp
        header(index+2) = k
        do i=1,k
          header(index+2+i) = ichaine(i)
        enddo
        index = index+3+k
        read (adummy,*,err=999,end=999) achaine,achaine,ityp,achaine
        do k=nlim,1,-1
          if( ichaine(k) .ne. ispace )goto 270
        enddo
        k = 1
  270   header(index) = k
        do i=1,k
          header(index+i) = ichaine(i)
        enddo
        index = index+1+k
        if( index .gt. nsize )goto 998
      enddo
      if( kcol .ne. header(3)+1 )goto 997
c
      return
c
  999 print 1000, nl, lun, atask
      header(1)=-1
      return
  998 print 1020, atask, lun
      header(1)=-1
      return
  997 print 1030, lun
      header(1)=-1
      return
 1000 format(' An error occured in routine read_header',/' after',i6,
     :       ' lines of unit',i3,/' when ',a80)
 1010 format(A)
 1020 format(' Out of memory when ',a80,/' in unit',i3)
 1030 format(' Wrong number of columns in unit',i3)
      end
C
C-----------------------------------------------------------------------
C
      subroutine to_string(kval,string)
c
      implicit      none
      integer*4     kval(*)
      character(*) string
c
      integer*4     klength, ichaine, i, k
      character*4   achaine
      equivalence   (ichaine,achaine)
c
      klength = min( kval(1), (len(string)-1)/4+1 )
      k = 2
      ichaine=kval(k)
      string=achaine
      k = k + 1
      do i=1,klength-1
        ichaine=kval(k)
        string=string(1:i*4)//achaine
        k=k+1
      enddo
c
      end
C
C-----------------------------------------------------------------------
C
      subroutine get_listof(header,ktype,table,ksize)
c
      implicit      none
      integer*4     header(*), ktype, ksize
      character(*) table(ksize)
c
      integer*4     ici, k
c
      if( header(1) .le. 0 )goto 999
c
      ici = ktype + 3
      k   = 0
  100 continue
        if( header(ici) .le. 0 )then
          ksize = k
          return
        endif
        ici = header(ici)
        k = k+1
        if( k .gt. ksize )return
        call to_string(header(ici+1),table(k))
      goto 100
  999 ksize=-1
      print 1000
 1000 format(' An error occured in routine get_listof')
       end
C
C-----------------------------------------------------------------------
C
      subroutine get_text(header,number,string)
c
      implicit      none
      integer*4     header(*), number
      character(*) string
c
      integer*4     ici, knt
c
      if( header(1) .le. 0 )goto 999
c
      ici = 4
      do knt=1,number
        if( header(ici) .le. 0 )goto 999
        ici = header(ici)
      enddo
c
      ici = ici+1
      call to_string(header(ici),string)
      return
  999 string = '*ERROR*'
      print 1000
 1000 format(' An error occured in routine get_text')
      end
C
C-----------------------------------------------------------------------
C
      subroutine select_var(header,name,number,unit,title,kolumn,ksize)
c
      implicit      none
      integer*4     header(*), kolumn(*), number, ksize
      character(*)  name(*), unit(*), title(*)
c
      integer*4     kvar, kol, idim, i
      logical       flag
c
      kol = 1
      flag = .FALSE.
      do kvar=1, number
        call info_var(header,name(kvar),unit(kol),title(kol),
     :                kolumn(kol),idim)
        if( kolumn(kol) .le. 0 .or. idim .le. 0)then
          print 1000,kvar,name(kvar)
          flag = .TRUE.
        elseif( kol+idim-1 .gt. ksize )then
          print 1010,kvar,name(kvar)
          flag = .TRUE.
        else
          do i=1,idim-1
            unit(kol+i) = unit(kol)
            title(kol+i) = title(kol)
            kolumn(kol+i) = kolumn(kol)+i
          enddo
          kol = kol + idim
        endif
      enddo
c  
      ksize = kol-1
      if( flag )then
        ksize = - ksize
        print 1020
      endif
      return
c
 1000 format(' Consequently variable #',i2,' not found : ',A12)
 1010 format(' Out of memory at variable #',i2,' : ',A12)
 1020 format(' One or several problems occured in routine select_var')
c
      end
C
C-----------------------------------------------------------------------
C
      subroutine select_one(header,name,unit,title,ndim,index,
     :                       ksize,kolumn,nkol)
c
      implicit      none
      integer*4     nkol, header(*), kolumn(nkol), ndim, index, ksize
      character(*) name, unit, title
c
      integer*4     kol,i
c
      index = -1
      if( ksize.lt.0 )then
          print 1020
          return
      endif
      kol=ksize+1
      call info_var(header,name,unit,title,kolumn(kol),ndim)
      if( kolumn(kol) .le. 0 .or. ndim .le. 0)then
          print 1000,name
          ksize = 1-kol
      elseif( kol+ndim-1 .gt. nkol )then
          print 1010,name
          ksize = 1-kol
      else
          index = kol
          do i=1,ndim-1
            kolumn(kol+i) = kolumn(kol)+i
          enddo
          kol = kol + ndim
          ksize = kol-1
      endif
      return
c
 1000 format(' Consequently variable not found : ',A12)
 1010 format(' Out of memory for variable : ',A12)
 1020 format(' Illegal value of argument ''ksize''')
      end
C
C-----------------------------------------------------------------------
C
      subroutine info_var(header,name,unit,title,kolumn,ndim)
c
      implicit      none
      integer*4     header(*), kolumn, ndim
      character(*) name, unit, title
c
      integer*4     ici, nlen, k
      character*80  string
c
      if( header(1) .le. 0 )goto 999
      nlen = min(80, len(name))
c
      ici = 6
  100 continue
        if( header(ici) .le. 0 )goto 999
        ici = header(ici)
        k = ici+1
        call to_string(header(k),string)
      if( string(1:nlen) .ne. name(1:nlen) )goto 100
c
      ici    = ici + header(k) + 4
      ndim   = header(ici-2)
      kolumn = header(ici-1)
      call to_string(header(ici),unit)
      ici    = ici + header(ici) + 1
      call to_string(header(ici),title)
c
      return
  999 kolumn = -1
      ndim   = 0
      unit   = '*ERROR*'
      title  = '*ERROR*'
      print 1000
 1000 format(' An error occured in routine info_var')
      end
C
C-----------------------------------------------------------------------
C
      subroutine get_meta(header,meta,value,ksize)
c
      implicit      none
      integer*4     header(*), value(*), ksize
c-HE
      real*8        Rvalue(*)
c-HE
      character(*) meta, svalue(*)
c
      integer*4     ici, k, ityp, nlen, ktemp
      real*4        rtemp
      character*80  string
c     lenterReal added by HE
      logical       lenterstring,lenterinteg,lenterReal
      equivalence   (ktemp,rtemp)
c
      character*4 dstring
      integer*4 dinteg
      equivalence (dstring,dinteg)
c
      lenterstring= .false.
      lenterinteg= .false.
      goto 10
      entry get_meta_string(header,meta,svalue,ksize)
      lenterstring= .true.
      lenterinteg= .false.
c-HE
      lenterreal = .false.
c-HE
      ksize=-abs(ksize)
      goto 10
      entry get_meta_int(header,meta,value,ksize)
      lenterstring= .false.
      lenterinteg= .true.
c-HE
      lenterreal = .false.
      goto 10
      entry get_meta_real(header,meta,rvalue,ksize)
      lenterstring= .false.
      lenterinteg= .false.
      lenterreal = .true.
c-HE

   10 continue
c
      if( header(1) .le. 0 )goto 999
      nlen = min(80, len(meta))
c
      ici = 5
  100 continue
        if( header(ici) .le. 0 )goto 999
        ici = header(ici)
        k = ici+1
        call to_string(header(k),string)
      if( string(1:nlen) .ne. meta(1:nlen) )goto 100
      ici = ici+header(k)+3
      ityp = header(ici-1)
      if( ityp*ksize .le.0 )goto 998
      if( lenterstring )then
        ityp = max(ityp,ksize)
        ksize = ityp
        do k=1,-ityp
          call to_string(header(ici),svalue(k))
          ici=ici+header(ici)+1
        enddo
      elseif( ityp .lt. 0 )then
        ityp = max(ityp,ksize)
        ksize = ityp
        do k=1,-ityp
          call to_string(header(ici),dstring)
          value(k) = dinteg
          ici=ici+header(ici)+1
        enddo
      elseif( lenterinteg )then
        ityp = min(ityp,ksize)
        ksize = ityp
        do k=1,ityp
          ktemp=header(ici-1+k)
          value(k)=rtemp
        enddo
c-HE  
      elseif( lenterreal )then
        ityp = min(ityp,ksize)
        ksize = ityp
        do k=1,ityp
          ktemp=header(ici-1+k)
          rvalue(k)=rtemp
        enddo
c-HE
      else
        ityp = min(ityp,ksize)
        ksize = ityp
        do k=1,ityp
          value(k)=header(ici-1+k)
        enddo
      endif
      return
  998 ksize = 0
c     print 1001
c-HE
      print 1001, meta
c-HE
      stop
  999 ksize = 0
c      print 1000
c 1000 format(' An error occured in routine get_meta: invalid name')
c 1001 format(' An error occured in routine get_meta: type mismatch')
c-HE
      print 1000, meta
 1000 format(' An error occured in routine get_meta: invalid name;',
     &       ' META=[',a,']')
 1001 format(' An error occured in routine get_meta: type mismatch;',
     &       ' META=[',a,']')
c-HE
      end
C
C-----------------------------------------------------------------------
C
      subroutine get_var(header,kolumn,ksize,kskip,number,
     :                   value,ndim,footer)
c
      implicit      none
      integer*4     header(*), kolumn(*), ksize, kskip, number, ndim
      character(*) footer
      real*8        value(ndim,*)
c
      character     achaine*32000, atask*80
      integer*4     lun, k, knumber, i, ncol, kcount
      real*8        dummy(2000)
c
      lun = header(1)
      ncol = header(3)
      footer = '*ERROR*'
      kcount = 0
      atask = 'checking header information'
      if( lun .le. 0 .or. ncol .le. 0)goto 999
c
      if( kskip .lt. 0 )then
        atask = 'skipping records to the end of the body section'
  100   read(lun,1000,err=999,end=999) achaine
        kcount = kcount + 1
        if( achaine(1:1) .ne. '''' )goto 100
        read(achaine,*,err=999,end=999) footer
        if( number.ne.0 )number = 0
        return
      endif
c
      atask = 'skipping records'
      do k=1, kskip
        read(lun,1000,err=999,end=999) achaine
        kcount = kcount + 1
        if( achaine(1:1) .eq. '''' )then
          read(achaine,*,err=999,end=999) footer
          number = 0
          return
        endif
      enddo
c
      atask = 'reading the first record'
      knumber = max( 0, number)
      do k=1, knumber
        read(lun,1000,err=999,end=999) achaine
        kcount = kcount + 1
        if( achaine(1:1) .eq. '''' )then
          read(achaine,*,err=999,end=999) footer
          number = k-1
          return
        endif
        atask = 'selecting data from a record'
        read(achaine,*,err=999,end=999) (dummy(i),i=1,ncol)
        do i=1, ksize
          value(i,k) = dummy( kolumn(i) )
        enddo
        atask = 'reading a record'
      enddo
      footer = '*CONTINUE*'
      number = knumber
      return
c
  999 number = -1
      footer = '*ERROR*'
      print 1050, kcount, lun, atask
c
 1000 format(A32000)
 1050 format(' An error occured in routine get_var',/' after',i6,
     :       ' lines of unit',i3,/' when ',a80)
c
      end