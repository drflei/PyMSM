# 1 "ut540.f"
      SUBROUTINE UT540
     :          (mdate)
C
C!    Compute modified Julian Day from date
C_    er24
C
C     REFERENCES
C     Subroutine datojd, sapre5.for of sapre
C     Method from H. Klinkrad 
C
      INCLUDE 'structure.h'
cDEC$ IF DEFINED (_x86_)
cDEC$ ATTRIBUTES DLLEXPORT :: UT540
cDEC$ ENDIF
C
C     INTERFACE            
C   
        TYPE(zdat) :: mdate
C
C     VARIABLES
C
        INTEGER*4     jj, l, jd
C
C     CODE
C
      jj         = (14-mdate%imonth) / 12 
      l          = mdate%iyear - 1900 - jj 
      jd         = mdate%iday - 18234 + (1461*l) / 4
     :             + (367*(mdate%imonth-2+jj*12)) / 12 
      mdate%amjd = jd + ((mdate%ihour*60+mdate%imin)*60+mdate%secs)
     :                      / 86400.0d0
C
C
      END
C----------------------------------------------------------------------
      SUBROUTINE UT541
     :          (mrtp,mxyz)
C
C     Convert spherical coordinates to cartesian coordinates
C
      INCLUDE 'structure.h'
cDEC$ IF DEFINED (_x86_)
cDEC$ ATTRIBUTES DLLEXPORT :: UT541
cDEC$ ENDIF
C
C     INTERFACE 
C              
        TYPE(zgeo) :: mrtp
        TYPE(zxyz) :: mxyz
*
*  mrtp spherical coordinates
*  mxyz cartesian coordinates
*
C
C    Subroutine UT541 converts the spherical coordinates of
C    a geographic position into cartesian coordinates. Be
C    aware that this subroutine is not relevent for transformation
C    between cartesian and spherical components of a vector.
C
C    "History" The subroutine does not exist in version 1.05
C                    and earlier        
C    "See also" ut546
C                         
C
      COMMON /UC160
     :               /pi, deg, re, gmagmo, eclipt, geoid, uma
        REAL*8        pi, deg, re, gmagmo, eclipt, geoid(3), uma(30)
C
C
C     VARIABLES
C
      REAL*8 st
C
      st     = sin( mrtp%colat*deg ) * mrtp%radius
      mxyz%z = cos( mrtp%colat*deg ) * mrtp%radius
      mxyz%x = st * cos( mrtp%elong*deg )
      mxyz%y = st * sin( mrtp%elong*deg )
C
      end
C----------------------------------------------------------------------
      SUBROUTINE UT542
     :          (mpos,mvec,mxyz)
C
C     Convert spherical vector components to cartesian components
C
      INCLUDE 'structure.h'
cDEC$ IF DEFINED (_x86_)
cDEC$ ATTRIBUTES DLLEXPORT :: UT542
cDEC$ ENDIF
C
C     INTERFACE 
C              
        TYPE(zgeo) :: mpos
        TYPE(zvec) :: mvec
        TYPE(zxyz) :: mxyz
*
*  mrtp spherical coordinates
*  mxyz cartesian coordinates
*
C
C    Subroutine UT542 converts the spherical components of
C    a vector into cartesian components. The transformation
C    is function of the geographic position where the
C    is evaluated.
C
C    "History" The subroutine does not exist in version 1.05
C                    and earlier        
C    "See also" ut547
C                         
C
      COMMON /UC160
     :               /pi, deg, re, gmagmo, eclipt, geoid, uma
        REAL*8        pi, deg, re, gmagmo, eclipt, geoid(3), uma(30)
C
C
C     VARIABLES
C
      REAL*8 st,sp,ct,cp
C
      st     = sin( mpos%colat*deg )
      ct     = cos( mpos%colat*deg )
      sp     = sin( mpos%elong*deg )
      cp     = cos( mpos%elong*deg )
      mxyz%z = mvec%rho *  ct   - mvec%theta *  st
      mxyz%x = mvec%rho * st*cp + mvec%theta * ct*cp - mvec%phi * sp
      mxyz%y = mvec%rho * st*sp + mvec%theta * ct*sp + mvec%phi * cp
C
      end
C----------------------------------------------------------------------
      SUBROUTINE UT545
     :          (mdate)
C
C!    Compute date from modified Julian Day
C_    er24
C
C     REFERENCES
C     Subroutine jdtoda, sapre5.for of sapre
C     Method from H. Klinkrad 
C
      INCLUDE 'structure.h'
cDEC$ IF DEFINED (_x86_)
cDEC$ ATTRIBUTES DLLEXPORT :: UT545
cDEC$ ENDIF
C
C     INTERFACE 
C              
        TYPE(zdat) :: mdate
C
C     VARIABLES
C  
        REAL*8        dday, rhour, rmin
        INTEGER*4     jday, l, m, n, jj
C
C     CODE
C  
      dday         = mdate%amjd + 1.0d-11
C
      jday         = dday 
      l            = (4000*(jday+18204)) / 1461001 
      n            = jday - (1461*l) / 4 + 18234 
      m            = (80*n) / 2447 
C           
      mdate%iday   = n - (2447*m) / 80 
      jj           = m / 11 
      mdate%imonth = m + 2 - 12 * jj 
      mdate%iyear  = 1900 + l + jj 
C
      rhour        = (dday-jday) * 24.0d0 
      mdate%ihour  = rhour 
      rmin         = (rhour-mdate%ihour) * 60.0d0
      mdate%imin   = rmin 
      mdate%secs   = (rmin-mdate%imin) * 60.0d0
C
C
      END
C----------------------------------------------------------------------
      SUBROUTINE UT546
     :          (mxyz,mrtp)
C
C     Convert cartesian coordinates to spherical coordinates
C
      INCLUDE 'structure.h'
cDEC$ IF DEFINED (_x86_)
cDEC$ ATTRIBUTES DLLEXPORT :: UT546
cDEC$ ENDIF
C
C     INTERFACE 
C              
        TYPE(zgeo) :: mrtp
        TYPE(zxyz) :: mxyz
*
*  mxyz cartesian coordinates
*  mrtp spherical coordinates
*
C
C    Subroutine UT546 converts the cartesian coordinates of
C    a geographic position into spherical coordinates. Be
C    aware that this subroutine is not relevent for transformation
C    between cartesian and spherical components of a vector.
C                         
C    "history" The subroutine does not exist in version 1.05
C                    and earlier        
C    "see also" ut541
C
C
      COMMON /UC160
     :               /pi, deg, re, gmagmo, eclipt, geoid, uma
        REAL*8        pi, deg, re, gmagmo, eclipt, geoid(3), uma(30)
C
C
C
      mrtp%radius  = sqrt(mxyz%x**2+mxyz%y**2+mxyz%z**2)
      mrtp%colat   = acos(mxyz%z/mrtp%radius)/deg
      if( (mxyz%x).eq.0 .and. (mxyz%y).eq.0)then
        mrtp%elong = 0
      else
        mrtp%elong = atan2(mxyz%y,mxyz%x)/deg
      endif
C
      end
C----------------------------------------------------------------------
      SUBROUTINE UT547
     :          (mpos,mxyz,mvec)
C
C     Convert cartesian vector components to spherical components
C
      INCLUDE 'structure.h'
cDEC$ IF DEFINED (_x86_)
cDEC$ ATTRIBUTES DLLEXPORT :: UT547
cDEC$ ENDIF
C
C     INTERFACE 
C              
        TYPE(zgeo) :: mpos
        TYPE(zvec) :: mvec
        TYPE(zxyz) :: mxyz
*
*  mrtp spherical coordinates
*  mxyz cartesian coordinates
*
C
C    Subroutine UT547 converts the cartesian components of
C    a vector into spherical components. The transformation
C    is function of the geographic position where the
C    is evaluated.
C
C    "history" The subroutine does not exist in version 1.05
C                    and earlier        
C    "see also" ut542
C                         
C
      COMMON /UC160
     :               /pi, deg, re, gmagmo, eclipt, geoid, uma
        REAL*8        pi, deg, re, gmagmo, eclipt, geoid(3), uma(30)
C
C
C     VARIABLES
C
      REAL*8 st,sp,ct,cp
C
      st     = sin( mpos%colat*deg )
      ct     = cos( mpos%colat*deg )
      sp     = sin( mpos%elong*deg )
      cp     = cos( mpos%elong*deg )
      mvec%rho   = mxyz%z *  ct + mxyz%x * st*cp + mxyz%y * st*sp   
      mvec%theta = -st * mxyz%z + mxyz%x * ct*cp + mxyz%y * ct*sp
      mvec%phi   =                - sp * mxyz%x  + mxyz%y *  cp
      mvec%dnrm  = sqrt(mxyz%x**2+mxyz%y**2+mxyz%z**2)
C
      end
C----------------------------------------------------------------------