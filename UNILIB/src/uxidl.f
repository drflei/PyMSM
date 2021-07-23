c
      integer*4 function uxidl(argc, argv)
c
      implicit none         
c
c IDL interface to:
c
      external ut990,um510,um520,ua610,um530,um539
      external ua630,ul220,ud310,ud317,ud330,uf420
      external ut980,ut998,uf410,um535,um536,ut541
      external ut542,ut546,ut547,ut555,ut556,ut545
      external um522,um538,uf411,ul240,ul242,ut986
      external setkmflg
c
c --- select the proper integer size for your architecture
c       IDL on 32bit machines -> integer*4; on 64bit machines -> integer*8
c      integer*4      argc
c      integer*4      argv(*)      
      integer*8      argc
      integer*8      argv(*)      
c
      integer*4      arg_cnt
c
      character*5    subname(33), request
      character*64   dumstr
      integer*4      ksub, npara(33), idl2ftn, ftn2idl
      integer*4      lnreq, lndum
c
      parameter (idl2ftn=1,ftn2idl=-1)
c
      data subname/'UL220','UD310','UD317','UD330','UF410','UF420',
     :             'UM510','UM520','UM530','UM535','UM536','UM539',  
     :     'UT540','UT541','UT542','UT545','UT546','UT547','UT555', 
     :     'UT556','UA610','UA630','UT980','UT990','UT998','UT550',
     :     'UM522','UM538','UF411','UL240','UL242','UT986','KMFLG'/
      data npara/      10,      6,      5,      3,      7,      5,
     :                  5,      6,      3,      2,      2,      4,
     :          1,      2,      3,      1,      2,      3,      3,
     :          4,     15,      5,      4,      3,      5,      4,
     :          2,      5,      9,      2,      2,      6,      1/
c
      request = 'XXXXX'
      dumstr  = '- Dummy string -- Dummy string -- Dummy string -'
     :          //'- Dummy string -'
      lnreq   = len(request)
      lndum   = len(dumstr)
c
      arg_cnt = loc(argc)
      if(arg_cnt .lt. 1)then
         write(*,1000) subname
         write(*,1010) 
         uxidl = -1.0D0
         return
      endif
      call str_conv(%val(argv(1)),request,lnreq,idl2ftn)
c
      do ksub=1,33
        if( subname(ksub).eq.request )then
          if( arg_cnt .ne. npara(ksub)+1 )then
            write(*,1020) request, npara(ksub), arg_cnt-1
            write(*,1010)
            uxidl = -1.0D0
            return
          endif
c
c
          if(     ksub .eq.  1 )then
            call UL220( %val(argv(2)), %val(argv(3)), %val(argv(4)),
     :                  %val(argv(5)), %val(argv(6)), %val(argv(7)),
     :                  %val(argv(8)), %val(argv(9)), %val(argv(10)),
     :                  %val(argv(11)) )   
          elseif( ksub .eq.  2 )then
            call UD310( %val(argv(2)), %val(argv(3)), %val(argv(4)),
     :                  %val(argv(5)), %val(argv(6)), %val(argv(7)) )
          elseif( ksub .eq.  3 )then
            call UD317( %val(argv(2)), %val(argv(3)), %val(argv(4)),
     :                  %val(argv(5)), %val(argv(6)) )
          elseif( ksub .eq.  4 )then
            call UD330( %val(argv(2)), %val(argv(3)), %val(argv(4)) )
          elseif( ksub .eq.  5 )then
            call UF410( %val(argv(2)), %val(argv(3)), %val(argv(4)),
     :                  %val(argv(5)), %val(argv(6)), %val(argv(7)),
     :                  %val(argv(8)) )
          elseif( ksub .eq.  6 )then
            call UF420( %val(argv(2)), %val(argv(3)), %val(argv(4)),
     :                  %val(argv(5)), %val(argv(6)) )
          elseif( ksub .eq.  7 )then
            call str_conv(%val(argv(4)),dumstr,lndum,idl2ftn)
            call UM510( %val(argv(2)), %val(argv(3)), dumstr,
     :                  %val(argv(5)), %val(argv(6)) )
            call str_conv(%val(argv(4)),dumstr,lndum,ftn2idl)
          elseif( ksub .eq.  8 )then
            call str_conv(%val(argv(5)),dumstr,lndum,idl2ftn)
            call UM520( %val(argv(2)), %val(argv(3)), %val(argv(4)),
     :                  dumstr, %val(argv(6)), %val(argv(7)) )
            call str_conv(%val(argv(5)),dumstr,lndum,ftn2idl)
          elseif( ksub .eq.  9 )then
            call UM530( %val(argv(2)), %val(argv(3)), %val(argv(4)) )
          elseif( ksub .eq. 10 )then
            call UM535( %val(argv(2)), %val(argv(3)) )
          elseif( ksub .eq. 11 )then
            call UM535( %val(argv(2)), %val(argv(3)) )
          elseif( ksub .eq. 12 )then
            call UM539( %val(argv(2)), %val(argv(3)), %val(argv(4)),
     :                  %val(argv(5)) )
          elseif( ksub .eq. 13 )then
            call UT540( %val(argv(2)) )
          elseif( ksub .eq. 14 )then
            call UT541( %val(argv(2)), %val(argv(3)) )
          elseif( ksub .eq. 15 )then
            call UT542( %val(argv(2)), %val(argv(3)), %val(argv(4)) )
          elseif( ksub .eq. 16 )then
            call UT545( %val(argv(2)) )
          elseif( ksub .eq. 17 )then
            call UT546( %val(argv(2)), %val(argv(3)) )
          elseif( ksub .eq. 18 )then
            call UT547( %val(argv(2)), %val(argv(3)), %val(argv(4)) )
          elseif( ksub .eq. 19 )then
            call UT555( %val(argv(2)), %val(argv(3)), %val(argv(4)) )
          elseif( ksub .eq. 20 )then
            call UT556( %val(argv(2)), %val(argv(3)), %val(argv(4)),
     :                  %val(argv(5)) )
          elseif( ksub .eq. 21 )then
            call str_conv(%val(argv(14)),dumstr,lndum,idl2ftn)
            call UA610( %val(argv(2)), %val(argv(3)), %val(argv(4)),
     :                  %val(argv(5)), %val(argv(6)), %val(argv(7)),
     :                  %val(argv(8)), %val(argv(9)), %val(argv(10)),
     :                  %val(argv(11)), %val(argv(12)), %val(argv(13)),
     :                  dumstr, %val(argv(15)), %val(argv(16)) )
            call str_conv(%val(argv(14)),dumstr,lndum,ftn2idl)
          elseif( ksub .eq. 22 )then
            call UA630( %val(argv(2)), %val(argv(3)), %val(argv(4)),
     :                  %val(argv(5)), %val(argv(6)) )
          elseif( ksub .eq. 23 )then
            call str_conv(%val(argv(14)),dumstr,lndum,idl2ftn)
            call UT980( dumstr, %val(argv(3)), %val(argv(4)),
     :                  %val(argv(5)) )
            call str_conv(%val(argv(14)),dumstr,lndum,ftn2idl)
          elseif( ksub .eq. 24 )then
            call UT990( %val(argv(2)), %val(argv(3)), %val(argv(4)) )
          elseif( ksub .eq. 25 )then
            call UT998( %val(argv(2)), %val(argv(3)), %val(argv(4)),
     :                  %val(argv(5)), %val(argv(6)) )
          elseif( ksub .eq. 26 )then
            call UT550( %val(argv(2)), %val(argv(3)), %val(argv(4)),
     :                  %val(argv(5)) )
          elseif( ksub .eq. 27 )then
            call UM522( %val(argv(2)), %val(argv(3)) )
          elseif( ksub .eq. 28 )then
            call UM538( %val(argv(2)), %val(argv(3)), %val(argv(4)),
     :                  %val(argv(5)), %val(argv(6)) )
          elseif( ksub .eq. 29 )then
            call UF411( %val(argv(2)), %val(argv(3)), %val(argv(4)),
     :                  %val(argv(5)), %val(argv(6)), %val(argv(7)),
     :                  %val(argv(8)), %val(argv(9)), %val(argv(10)) )
          elseif( ksub .eq. 30 )then
            call UL240( %val(argv(2)), %val(argv(3)) )
          elseif( ksub .eq. 31 )then
            call UL242( %val(argv(2)), %val(argv(3)) )
          elseif( ksub .eq. 32 )then
            call UT986( %val(argv(2)), %val(argv(3)), %val(argv(4)),
     :                  %val(argv(5)), %val(argv(6)), %val(argv(7)) )
          elseif( ksub .eq. 33 )then
            call SETKMFLG( %val(argv(2)) )
          else
            write(*,1040) ksub
            write(*,1010)
            uxidl = -1
            return
          endif
c
c
          uxidl = 0
          return
        endif
      enddo
c
c
      write(*,1030) request, subname
      write(*,1010)
      uxidl = -1
      return        
c
c
 1000 format(' UXIDL:',/3x,'* Usage:',/5x,'status = CALL_EXTERNAL(',
     :       'Image, ''UXIDL'', Entry [, P1, ..., Pn])',/5x,'where',
     :       /7X,'- Image is the pathname of the sharable object file;',
     :       /7x,'- Entry is the name of the UNILIB subroutine;',
     :       /7x,'- P1, ..., Pn are the parameter to be passed.',
     :       /5x,'The supported entries are:', 5(/3x,6(4x,a5)))
 1010 format( 3X,'* Call_external aborted')
 1020 format(' UXIDL:',/3x,'* Subroutine ',a5,' needs',i3,' arguments',
     :      /5x,'(',i3,' arguments provided through the CALL_EXTERNAL)')
 1030 format(' UXIDL:',/3x,'* First argument of CALL_EXTERNAL (',a5,')',
     :       /5x,'not recognized as an supported entry.',
     :       /5x,'The supported entries are:', 5(/3x,6(4x,a5)))
 1040 format(' UXIDL:',/3x,'* Internal error (',i7.7,')')
      end
c
c
c
      subroutine str_conv(idlrec,ftnstr,lenstr,kdir)
         implicit       none
c
c kdir =  1........idl to fortran
c        -1........fortran to idl
c
c
         type string
           sequence
           integer*4    slen
           integer*2    stype
           integer      s
         end type string
c
         type(string) :: idlrec
         integer*2      slen2
         integer*4      kdir,lenstr
         character*(*)  ftnstr
c
         slen2 = idlrec%slen
c

         call str_cnv2(%val(idlrec%s),slen2,ftnstr,lenstr,kdir,
     :                 %val(slen2))
c
      end
c
      subroutine str_cnv2(idlstr,strlen,ftnstr,lenstr,kdir)
         implicit       none
c
         integer*4      kdir, lenstr
         character*(*)  idlstr
         character*(*)  ftnstr
         integer*2      strlen, kmin
c
         kmin = min(strlen,lenstr)
c
         if( kdir .eq. 1 )then
           ftnstr(1:kmin)=idlstr(1:kmin)
         elseif( kdir .eq. -1 )then
           idlstr(1:kmin)=ftnstr(1:kmin)
         else
           write(*,1000) kdir
         endif
 1000    format(' UXIDL:',/3x,'* Internal error in str_cnv2 (',i3.3,')',
     :         /3x,'* Execution not aborted !')
      end

      SUBROUTINE SETKMFLG(KMFLAG)

cDEC$ IF DEFINED (_x86_)
cDEC$ ATTRIBUTES DLLEXPORT :: KMFLG
cDEC$ ENDIF

      COMMON /UC190
     :               /prop, stepx, stpmin, umsq, upsq, uk2,
     :                uk3, epskm, epsrel, stplst, xclat,
     :                kmflg, kum533      
C
        REAL*8        prop, stepx, stpmin, umsq, upsq, uk2,
     :                uk3, epskm, epsrel, stplst, xclat
        INTEGER*4     kmflg, kum533       
C
      INTEGER*4 KMFLAG

      KMFLG = KMFLAG

      RETURN
      END
