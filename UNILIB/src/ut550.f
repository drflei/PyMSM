# 1 "ut550.f"
      subroutine ut550
     :                 (kfrom,kto,trans,ifail)
c
c     select a coordinate transformation
c
      include 'structure.h'

cDEC$ IF DEFINED (_x86_)
cDEC$ ATTRIBUTES DLLEXPORT :: UT550
cDEC$ ENDIF
c
*     kfrom input coordinate system
*     kto   output coordinate system
*     trans matrix transformation
*     ifail
C
C     The sub. ut550 initializes the coordinate transformation.
C     The computed transformation can be applied then with the
C     help of sub ut555 or ut556  The arguments kfrom and kto are
C     used to specify the input and output coordinate systems,
C     respectively.  
c     The resulting matrix trans allows to transform cartesian
c     coordinate of type kfrom to type kto.
c
c     The correspondance between indices and
C     coordinate systems is shown in the table below
C
C     kfrom, kto     coordinate system
C       1              Geographic coordinate (GEO)
C       2              Geocentric equatorial inertial coordinate (GEI)
C       3              Geomagnetic coordinate (MAG)
C       4              Solar magnetic coordinate (SM)
C       5              Geocentric solar magnetospheric coordinate (GSM)
c
c     The characteristics of the five coordinate systems are
c     sumarized in the next table.  The origin of all these coordinate 
c     systems are located at the centre of the Earth.
c     
c
c                X axis                Y axis              Z axis
c      GEO    included in the Earth's                 parallel to the Earth's
c             equatorial plane and     Z x X          rotation axis and pointing
c             passing through the                     to the North
c             Greenwich meridian
c
c      GEI    pointing to the first    Z x X          parallel to the Earth's
c             point of Aries, i.e.                    rotation axis and pointing
c             at the intersection                     to the North
c             between the Earth's 
c             equatorial plane and 
c             the ecliptic plane
c
c      MAG       Y x Z          perpendicular to      parallel to the geomagnetic
c                               both the geomagnetic  dipole axis and pointing
c                               dipole axis and the   to the North
c                               Earth's rotation axis
c                               in the East direction
c
c      SM        Y x Z          perpendicular to the  parallel to the geomagnetic
c                               Earth-Sun line        dipole axis and pointing
c                               towards dusk          to the North
c
c      GSM    pointing towards  perpendicular to the       X x Y
c             the Sun           geomagnetic dipole
c                               axis such that the
c                               Z axis is pointing to
c                               the North
c
c     Note that the coordinate transformation generally depends
c     on the Sun position and/or the geomagnetic dipole field.
C     For this information, the subroutine ut550 makes use of the
c     data stored in the argument mint and msun of the common 
c     block uc140.  The common block uc140 can be initialized by the
c     subroutines um510 and um520. The Sun position can be modified
c     without affecting the external magnetic field model by a call
c     to subroutine um522
c
c     history....
c
c     reference (demander a Bart)
c              Olson 1970
c              Russell 1971
c
c     see also  ut555 ut556 um510 um522
c
c
      integer*4 kfrom, kto, ifail
      real*8    trans(3,3)
c
      COMMON /UC140/ mint, mext, msun

      TYPE(zimf) :: mint
      TYPE(zsun) :: msun
      TYPE(zemf) :: mext          
c
      COMMON /UC160/ pi, deg, re, gmagmo, eclipt, geoid, uma

      REAL*8         pi, deg, re
      REAL*8         gmagmo
      REAL*8         eclipt, geoid(3), uma(30)               
c
      real*8 cg,sg,ce,se,cc,sc,rot(3,3),st,dummy(3,3),tmp(3)
      TYPE(zxyz) :: maxe
c
      ifail      = 0
c
      dummy(1,1) = 1      
      dummy(1,2) = 0      
      dummy(1,3) = 0      
      dummy(2,1) = 0      
      dummy(2,2) = 1      
      dummy(2,3) = 0      
      dummy(3,1) = 0      
      dummy(3,2) = 0      
      dummy(3,3) = 1
c
      if(kfrom.lt.1.or.kfrom.gt.6.or.kto.lt.1.or.kto.gt.6)then
        ifail=-55001
c       invalid coordinate system
        return
      endif
c
      if( kfrom.eq.kto ) then
        trans(1,1) = 1      
        trans(1,2) = 0      
        trans(1,3) = 0      
        trans(2,1) = 0      
        trans(2,2) = 1      
        trans(2,3) = 0      
        trans(3,1) = 0      
        trans(3,2) = 0      
        trans(3,3) = 1                   
        return
      endif
c
      if( kfrom.eq.2) then
c         gei-->geo
c
        cg = cos(msun%gha*deg)
        sg = sin(msun%gha*deg)
        dummy(1,1) = cg
        dummy(2,2) = cg
        dummy(1,2) = sg
        dummy(2,1) = -sg
c
      elseif( kfrom.eq.3) then
c         mag-->geo
c
c         y-rot
        cc = cos(mint%colat*deg)
        sc = sin(mint%colat*deg)
        dummy(3,1) = -sc
        dummy(3,3) = cc
c         z-rot
        ce = cos(mint%elong*deg)
        se = sin(mint%elong*deg)         
        dummy(1,1) = ce*cc
        dummy(1,2) = -se
        dummy(1,3) = ce*sc
        dummy(2,1) = se*cc
        dummy(2,2) = ce
        dummy(2,3) = se*sc
      elseif( kfrom.eq.4) then
c         sm-->geo
c
        ce = cos(mint%elong*deg)
        se = sin(mint%elong*deg)   
        cc = cos(mint%colat*deg)
        sc = sin(mint%colat*deg)
c         z:
        dummy(1,3) = sc*ce
        dummy(2,3) = sc*se
        dummy(3,3) = cc              
c         y:
        maxe%x = sc*se * msun%dir%z - cc    * msun%dir%y
        maxe%y = cc    * msun%dir%x - sc*ce * msun%dir%z
        maxe%z = sc*ce * msun%dir%y - sc*se * msun%dir%x
        st = sqrt(maxe%x**2+maxe%y**2+maxe%z**2)
        if( st.eq.0 )then
          ifail=-55003
c         problem with the Sun direction
          return
        endif
c-old-        if(maxe%z.lt.0)st=-st
        if( ((maxe%y*dummy(3,3)-maxe%z*dummy(2,3))*msun%dir%x + 
     :     (maxe%z*dummy(1,3)-maxe%x*dummy(3,3))*msun%dir%y +
     :     (maxe%x*dummy(2,3)-maxe%y*dummy(1,3))*msun%dir%z   
     :     ) .lt. 0 ) st=-st
        dummy(1,2)   = maxe%x/st
        dummy(2,2)   = maxe%y/st
        dummy(3,2)   = maxe%z/st
c         x:
        dummy(1,1) = dummy(2,2) * dummy(3,3) - dummy(3,2) * dummy(2,3)
        dummy(2,1) = dummy(3,2) * dummy(1,3) - dummy(1,2) * dummy(3,3)
        dummy(3,1) = dummy(1,2) * dummy(2,3) - dummy(2,2) * dummy(1,3)
      elseif( kfrom.eq.5) then
c         gsm-->geo
c
        ce = cos(mint%elong*deg)
        se = sin(mint%elong*deg)   
        cc = cos(mint%colat*deg)
        sc = sin(mint%colat*deg)
c         x:sun
        dummy(1,1)   = msun%dir%x
        dummy(2,1)   = msun%dir%y
        dummy(3,1)   = msun%dir%z
c         y:
        maxe%x = sc*se * msun%dir%z - cc    * msun%dir%y
        maxe%y = cc    * msun%dir%x - sc*ce * msun%dir%z
        maxe%z = sc*ce * msun%dir%y - sc*se * msun%dir%x
        st = sqrt(maxe%x**2+maxe%y**2+maxe%z**2)
        if( st.eq.0 )then
          ifail=-55003
c         problem with the Sun direction
          return
        endif
        if( (dummy(1,1)*maxe%y-dummy(2,1)* maxe%x) .lt. 0)st=-st
        dummy(1,2)   = maxe%x/st
        dummy(2,2)   = maxe%y/st
        dummy(3,2)   = maxe%z/st
c         z:
        dummy(1,3) = dummy(2,1) * dummy(3,2) - dummy(3,1) * dummy(2,2)
        dummy(2,3) = dummy(3,1) * dummy(1,2) - dummy(1,1) * dummy(3,2)
        dummy(3,3) = dummy(1,1) * dummy(2,2) - dummy(2,1) * dummy(1,2)
      elseif( kfrom.eq.6) then
c         gse-->geo
c
        cg = cos(msun%gha*deg)
        sg = sin(msun%gha*deg)
c         x:sun
        dummy(1,1)   = msun%dir%x
        dummy(2,1)   = msun%dir%y
        dummy(3,1)   = msun%dir%z
c         z: xgse X x gei
        maxe%x       = msun%dir%z*sg
        maxe%y       = msun%dir%z*cg
        maxe%z       = -msun%dir%x*sg-msun%dir%y*cg
        st = sqrt(maxe%x**2+maxe%y**2+maxe%z**2)
        if( st.eq.0 )then
          ifail=-55003
c         problem with the Sun direction
          return
        endif
        if(maxe%z.lt.0)st=-st
        dummy(1,3)   = maxe%x/st
        dummy(2,3)   = maxe%y/st
        dummy(3,3)   = maxe%z/st
c         y:
        dummy(1,2) = dummy(3,1) * dummy(2,3) - dummy(2,1) * dummy(3,3)
        dummy(2,2) = dummy(1,1) * dummy(3,3) - dummy(3,1) * dummy(1,3)
        dummy(3,2) = dummy(2,1) * dummy(1,3) - dummy(1,1) * dummy(2,3)
c                                                                 
      endif
c
c
c
      if( kto.eq.1)then
c
        trans(1,1) = dummy(1,1)      
        trans(1,2) = dummy(1,2)       
        trans(1,3) = dummy(1,3)       
        trans(2,1) = dummy(2,1)       
        trans(2,2) = dummy(2,2)       
        trans(2,3) = dummy(2,3)       
        trans(3,1) = dummy(3,1)       
        trans(3,2) = dummy(3,2)       
        trans(3,3) = dummy(3,3)                       
      elseif( kto.eq.2) then
c         geo-->gei
c
        cg = cos(msun%gha*deg)
        sg = sin(msun%gha*deg)
        trans(1,1) = cg*dummy(1,1) - sg*dummy(2,1)
        trans(1,2) = cg*dummy(1,2) - sg*dummy(2,2)
        trans(1,3) = cg*dummy(1,3) - sg*dummy(2,3)
        trans(2,1) = sg*dummy(1,1) + cg*dummy(2,1)
        trans(2,2) = sg*dummy(1,2) + cg*dummy(2,2) 
        trans(2,3) = sg*dummy(1,3) + cg*dummy(2,3) 
        trans(3,1) = dummy(3,1) 
        trans(3,2) = dummy(3,2) 
        trans(3,3) = dummy(3,3) 
c
      elseif( kto.eq.3) then
c         geo-->mag
c
c         z-rot
        ce = cos(mint%elong*deg)
        se = sin(mint%elong*deg)   
        tmp(1)     = ce*dummy(1,1) + se*dummy(2,1)
        tmp(2)     = ce*dummy(1,2) + se*dummy(2,2)
        tmp(3)     = ce*dummy(1,3) + se*dummy(2,3)
        trans(2,1) = - se*dummy(1,1) + ce*dummy(2,1)
        trans(2,2) = - se*dummy(1,2) + ce*dummy(2,2)
        trans(2,3) = - se*dummy(1,3) + ce*dummy(2,3)
c         y-rot
        cc = cos(mint%colat*deg)
        sc = sin(mint%colat*deg)
        trans(1,1) = cc*tmp(1) - sc*dummy(3,1)
        trans(1,2) = cc*tmp(2) - sc*dummy(3,2)
        trans(1,3) = cc*tmp(3) - sc*dummy(3,3)
        trans(3,1) = sc*tmp(1) + cc*dummy(3,1)
        trans(3,2) = sc*tmp(2) + cc*dummy(3,2)
        trans(3,3) = sc*tmp(3) + cc*dummy(3,3)
c
      elseif( kto.eq.4 ) then
c         geo-->sm
c
        ce = cos(mint%elong*deg)
        se = sin(mint%elong*deg)   
        cc = cos(mint%colat*deg)
        sc = sin(mint%colat*deg)
c         z:
        rot(3,1) = sc*ce
        rot(3,2) = sc*se
        rot(3,3) = cc              
c         y:
        maxe%x = sc*se * msun%dir%z - cc    * msun%dir%y
        maxe%y = cc    * msun%dir%x - sc*ce * msun%dir%z
        maxe%z = sc*ce * msun%dir%y - sc*se * msun%dir%x
        st = sqrt(maxe%x**2+maxe%y**2+maxe%z**2)
        if( st.eq.0 )then
          ifail=-55004
c         problem with the Sun direction
          return
        endif
        if( ( (maxe%y*rot(3,3)-maxe%z*rot(3,2))*msun%dir%x +
     :      (maxe%z*rot(3,1)-maxe%x*rot(3,3))*msun%dir%y +
     :      (maxe%x*rot(3,2)-maxe%y*rot(3,1))*msun%dir%z
     :    ).lt.0)st=-st
        rot(2,1)   = maxe%x/st
        rot(2,2)   = maxe%y/st
        rot(2,3)   = maxe%z/st
c         x:
        rot(1,1)   = rot(2,2) * rot(3,3) - rot(2,3) * rot(3,2)
        rot(1,2)   = rot(2,3) * rot(3,1) - rot(2,1) * rot(3,3)
        rot(1,3)   = rot(2,1) * rot(3,2) - rot(2,2) * rot(3,1)
c                                                                 
        trans(1,1) = rot(1,1)*dummy(1,1)+rot(1,2)*dummy(2,1)+
     :               rot(1,3)*dummy(3,1)    
        trans(1,2) = rot(1,1)*dummy(1,2)+rot(1,2)*dummy(2,2)+
     :               rot(1,3)*dummy(3,2)     
        trans(1,3) = rot(1,1)*dummy(1,3)+rot(1,2)*dummy(2,3)+
     :               rot(1,3)*dummy(3,3)     
        trans(2,1) = rot(2,1)*dummy(1,1)+rot(2,2)*dummy(2,1)+
     :               rot(2,3)*dummy(3,1)        
        trans(2,2) = rot(2,1)*dummy(1,2)+rot(2,2)*dummy(2,2)+
     :               rot(2,3)*dummy(3,2)        
        trans(2,3) = rot(2,1)*dummy(1,3)+rot(2,2)*dummy(2,3)+
     :               rot(2,3)*dummy(3,3)        
        trans(3,1) = rot(3,1)*dummy(1,1)+rot(3,2)*dummy(2,1)+
     :               rot(3,3)*dummy(3,1)        
        trans(3,2) = rot(3,1)*dummy(1,2)+rot(3,2)*dummy(2,2)+
     :               rot(3,3)*dummy(3,2)        
        trans(3,3) = rot(3,1)*dummy(1,3)+rot(3,2)*dummy(2,3)+
     :               rot(3,3)*dummy(3,3)              
      elseif( kto.eq.5) then
c         geo-->gsm
c
        ce = cos(mint%elong*deg)
        se = sin(mint%elong*deg)   
        cc = cos(mint%colat*deg)
        sc = sin(mint%colat*deg)
c         x:sun
        rot(1,1)   = msun%dir%x
        rot(1,2)   = msun%dir%y
        rot(1,3)   = msun%dir%z
c         y:
        maxe%x = sc*se * msun%dir%z - cc    * msun%dir%y
        maxe%y = cc    * msun%dir%x - sc*ce * msun%dir%z
        maxe%z = sc*ce * msun%dir%y - sc*se * msun%dir%x
        st = sqrt(maxe%x**2+maxe%y**2+maxe%z**2)
        if( st.eq.0 )then
          ifail=-55002
c         problem with the Sun direction
          return
        endif
c-old-        if(maxe%z.lt.0)st=-st
        if( (rot(1,1)*maxe%y - rot(1,2)*maxe%x) .lt. 0) st=-st
        rot(2,1)   = maxe%x/st
        rot(2,2)   = maxe%y/st
        rot(2,3)   = maxe%z/st
c         z:
        rot(3,1)   = rot(1,2) * rot(2,3) - rot(1,3) * rot(2,2)
        rot(3,2)   = rot(1,3) * rot(2,1) - rot(1,1) * rot(2,3)
        rot(3,3)   = rot(1,1) * rot(2,2) - rot(1,2) * rot(2,1)
c                                                                 
        trans(1,1) = rot(1,1)*dummy(1,1)+rot(1,2)*dummy(2,1)+
     :               rot(1,3)*dummy(3,1)    
        trans(1,2) = rot(1,1)*dummy(1,2)+rot(1,2)*dummy(2,2)+
     :               rot(1,3)*dummy(3,2)     
        trans(1,3) = rot(1,1)*dummy(1,3)+rot(1,2)*dummy(2,3)+
     :               rot(1,3)*dummy(3,3)     
        trans(2,1) = rot(2,1)*dummy(1,1)+rot(2,2)*dummy(2,1)+
     :               rot(2,3)*dummy(3,1)        
        trans(2,2) = rot(2,1)*dummy(1,2)+rot(2,2)*dummy(2,2)+
     :               rot(2,3)*dummy(3,2)        
        trans(2,3) = rot(2,1)*dummy(1,3)+rot(2,2)*dummy(2,3)+
     :               rot(2,3)*dummy(3,3)        
        trans(3,1) = rot(3,1)*dummy(1,1)+rot(3,2)*dummy(2,1)+
     :               rot(3,3)*dummy(3,1)        
        trans(3,2) = rot(3,1)*dummy(1,2)+rot(3,2)*dummy(2,2)+
     :               rot(3,3)*dummy(3,2)        
        trans(3,3) = rot(3,1)*dummy(1,3)+rot(3,2)*dummy(2,3)+
     :               rot(3,3)*dummy(3,3)              
c
      elseif( kto.eq.6) then
c         geo-->gse
c
        cg = cos(msun%gha*deg)
        sg = sin(msun%gha*deg)
c         x:sun
        rot(1,1)   = msun%dir%x
        rot(1,2)   = msun%dir%y
        rot(1,3)   = msun%dir%z
c         z: xgse X x gei
        maxe%x       = msun%dir%z*sg
        maxe%y       = msun%dir%z*cg
        maxe%z       = -msun%dir%x*sg-msun%dir%y*cg
        st = sqrt(maxe%x**2+maxe%y**2+maxe%z**2)
        if( st.eq.0 )then
          ifail=-55003
c         problem with the Sun direction
          return
        endif
        if(maxe%z.lt.0)st=-st
        rot(3,1)   = maxe%x/st
        rot(3,2)   = maxe%y/st
        rot(3,3)   = maxe%z/st
c         y:
        rot(2,1) = rot(1,3) * rot(3,2) - rot(1,2) * rot(3,3)
        rot(2,2) = rot(1,1) * rot(3,3) - rot(1,3) * rot(3,1)
        rot(2,3) = rot(1,2) * rot(3,1) - rot(1,1) * rot(3,2)
c
        trans(1,1) = rot(1,1)*dummy(1,1)+rot(1,2)*dummy(2,1)+
     :               rot(1,3)*dummy(3,1)    
        trans(1,2) = rot(1,1)*dummy(1,2)+rot(1,2)*dummy(2,2)+
     :               rot(1,3)*dummy(3,2)     
        trans(1,3) = rot(1,1)*dummy(1,3)+rot(1,2)*dummy(2,3)+
     :               rot(1,3)*dummy(3,3)     
        trans(2,1) = rot(2,1)*dummy(1,1)+rot(2,2)*dummy(2,1)+
     :               rot(2,3)*dummy(3,1)        
        trans(2,2) = rot(2,1)*dummy(1,2)+rot(2,2)*dummy(2,2)+
     :               rot(2,3)*dummy(3,2)        
        trans(2,3) = rot(2,1)*dummy(1,3)+rot(2,2)*dummy(2,3)+
     :               rot(2,3)*dummy(3,3)        
        trans(3,1) = rot(3,1)*dummy(1,1)+rot(3,2)*dummy(2,1)+
     :               rot(3,3)*dummy(3,1)        
        trans(3,2) = rot(3,1)*dummy(1,2)+rot(3,2)*dummy(2,2)+
     :               rot(3,3)*dummy(3,2)        
        trans(3,3) = rot(3,1)*dummy(1,3)+rot(3,2)*dummy(2,3)+
     :               rot(3,3)*dummy(3,3)              
c                                                                 
      endif
c
c
c
      return
c
      end
C----------------------------------------------------------------------
      subroutine ut551
     :                (alpha, beta, gamma, ktrl, trans)
c
c     Initialize an Euler rotation matrix
c
      include 'structure.h'
c
*     alpha first Euler angle expressed in degrees             
*     beta  second Euler angle expressed in degrees             
*     gamma third Euler angle expressed in degrees
*     ktrl  control parameter
*     trans matrix transformation
c
c     The subr ut551 initializes or modifies the coordinate
c     transformation with the help of Euler angles. The 
c     rotation matrix is defined by three consecutive rotations:
c       1. a rotation around the z-axis by the angle alpha
c       2. a rotation around the y-axis by the angle beta
c       3. a rotation around the z-axis by the angle gamma
c     Note that the rotations are done every time around an axis
c     formulated from the previous rotation.
c
c     The argument ktrl controls the way how the rotation
c     matrix trans id produced.  When ktrl is equal to zero the
c     Euler rotation matrix is copied into the argument
c     trans; otherwise, the Euler rotation matrix is combined
c     with the rotation matrix already stored in argument
c     trans.  When the value of arguemnt ktrl is negative,
c     the Euler rotation is applied before the rotation
c     stored in argument trans. When the value is positive,
c     the Euler rotation is applied after the stored rotation.
c     The argument ktrl allows thus to easily combine the
c     functionalities of subroutines ut550, ut551 and ut552.
c
c     History....
c     See also ut550, ut552
c              !!!! dans ut550 rajouter see also vers ut551, ut552
c
      real*8    alpha, beta, gamma, trans(3,3)
      integer*4 ktrl
c
c
      COMMON /UC160/ pi, deg, re, gmagmo, eclipt, geoid, uma

      REAL*8         pi, deg, re
      REAL*8         gmagmo
      REAL*8         eclipt, geoid(3), uma(30)               
c                                       
c
      real*8    temp(2,3,3)
      integer*4 ik,i,j
      real*8    ca,sa,cg,sg,cb,sb           
c
c
      if( ktrl.eq.0 )then
        temp(1,1,1) =1.0d0
        temp(1,1,2) =0.0d0
        temp(1,1,3) =0.0d0
        temp(1,2,1) =0.0d0
        temp(1,2,2) =1.0d0
        temp(1,2,3) =0.0d0
        temp(1,3,1) =0.0d0
        temp(1,3,2) =0.0d0
        temp(1,3,3) =1.0d0
        ik          = 2
      else
        if( ktrl.gt.0 )then
          ik        = 2
        else
          ik        = 1
        endif
        do i=1,3
          do j=1,3
            temp(3-ik,i,j)=trans(i,j)
          enddo
        enddo
      endif
c                       
c
      ca  = cos(alpha*deg)
      sa  = sin(alpha*deg)
      cb  = cos(beta*deg)
      sb  = sin(beta*deg)
      cg  = cos(gamma*deg)
      sg  = sin(gamma*deg)
c
      temp(ik,1,1) = cg*cb*ca - sg*sa
      temp(ik,1,2) = cg*cb*sa + sg*ca
      temp(ik,1,3) = - cg*sb
      temp(ik,2,1) = - sg*cb*ca - cg*sa
      temp(ik,2,2) = - sg*cb*sa + cg*ca
      temp(ik,2,3) = sg*sb
      temp(ik,3,1) = sb*ca
      temp(ik,3,2) = sb*sa
      temp(ik,3,3) = cb                 
c
      do i=1,3
        do j=1,3
c
          trans(i,j)=temp(2,i,1)*temp(1,1,j)+
     :               temp(2,i,2)*temp(1,2,j)+
     :               temp(2,i,3)*temp(1,3,j)

        enddo
      enddo
c                    
      end
C----------------------------------------------------------------------
      subroutine ut552
     :                (quat, ktrl, trans)
c
c     Initialize a quaternion rotation matrix
c
      include 'structure.h'
c 
*     quat  quaternion
*     ktrl  control parameter
*     trans matrix transformation
c
c     The subr ut551 initializes or modifies the coordinate
c     transformation with the help of one quaternion. A quaternion
c     is a four-element vector that generalizes the complex
c     number and which can be used to represent a rotation in
c     a three dimensional space. Note that 
c     q_1^2+q_2^2+q_3^2+q_4^2 has to be equal to one.
c
c     The argument ktrl controls the way how the rotation
c     matrix trans id produced.  When ktrl is equal to zero the
c     Euler rotation matrix is copied into the argument
c     trans; otherwise, the Euler rotation matrix is combined
c     with the rotation matrix already stored in argument
c     trans.  When the value of arguemnt ktrl is negative,
c     the Euler rotation is applied before the rotation
c     stored in argument trans. When the value is positive,
c     the Euler rotation is applied after the stored rotation.
c     The argument ktrl allows thus to easily combine the
c     functionalities of subroutines ut550, ut551 and ut552.
c
c     History....
c     See also ut550, ut552
c              [x]-> ...?
c
      real*8    quat(4), trans(3,3)
      integer*4 ktrl
c
c
      real*8    temp(2,3,3)
      integer*4 ik,i,j           
c                                       
c
      if( ktrl.eq.0 )then
        temp(1,1,1) =1.0d0
        temp(1,1,2) =0.0d0
        temp(1,1,3) =0.0d0
        temp(1,2,1) =0.0d0
        temp(1,2,2) =1.0d0
        temp(1,2,3) =0.0d0
        temp(1,3,1) =0.0d0
        temp(1,3,2) =0.0d0
        temp(1,3,3) =1.0d0
        ik          = 2
      else
        if( ktrl.gt.0 )then
          ik        = 2
        else
          ik        = 1
        endif
        do i=1,3
          do j=1,3
            temp(3-ik,i,j)=trans(i,j)
          enddo
        enddo
      endif
c                       
      temp(ik,1,1) = quat(1)**2-quat(2)**2-quat(3)**2+quat(4)**2
      temp(ik,1,2) = 2.0d0 * (quat(1)*quat(2)+quat(3)*quat(4))
      temp(ik,1,3) = 2.0d0 * (quat(1)*quat(3)-quat(2)*quat(4))
      temp(ik,2,1) = 2.0d0 * (quat(1)*quat(2)-quat(3)*quat(4))
      temp(ik,2,2) = -quat(1)**2+quat(2)**2-quat(3)**2+quat(4)**2
      temp(ik,2,3) = 2.0d0 * (quat(2)*quat(3)+quat(1)*quat(4))
      temp(ik,3,1) = 2.0d0 * (quat(1)*quat(3)+quat(2)*quat(4))
      temp(ik,3,2) = 2.0d0 * (quat(2)*quat(3)-quat(1)*quat(4))
      temp(ik,3,3) = -quat(1)**2-quat(2)**2+quat(3)**2+quat(4)**2
c
      do i=1,3
        do j=1,3
c
          trans(i,j)=temp(2,i,1)*temp(1,1,j)+
     :               temp(2,i,2)*temp(1,2,j)+
     :               temp(2,i,3)*temp(1,3,j)

        enddo
      enddo
c                    
      end
C----------------------------------------------------------------------
      subroutine ut555
     :                (mfrom,mto,trans)
c
c     coordinate conversion
c
      include 'structure.h'             
cDEC$ IF DEFINED (_x86_)
cDEC$ ATTRIBUTES DLLEXPORT :: UT555
cDEC$ ENDIF
c
c
c     The subr ut555 converts coordinates from a coordinate system
c     to another according to the transformation matrix trans passed
c     as argument. The matrix trans corresponds generally to a rotation
c     and can be evaluated by the help of subroutine ut550.
c     The subr converts the location mfrom expressed in the 'old'
c     coordinate system to the same location but expressed in the
c     'new' coordinate system
c
c     history
c     see also ut556
c
*
*     mfrom location in the 'old' system
*     mto   location in the 'new' system
*     trans transformation matrix
*
      TYPE(zgeo) :: mfrom, mto
      real*8        trans(3,3)
c
c
c
      COMMON /UC160/ pi, deg, re, gmagmo, eclipt, geoid, uma

      REAL*8         pi, deg, re
      REAL*8         gmagmo
      REAL*8         eclipt, geoid(3), uma(30)               
c                                       
c
c
      TYPE(zxyz) :: mcafr,mcato
      real*8       st
c                        
c
c
      st     = sin( mfrom%colat*deg ) * mfrom%radius
      mcafr%z = cos( mfrom%colat*deg ) * mfrom%radius
      mcafr%x = st * cos( mfrom%elong*deg )
      mcafr%y = st * sin( mfrom%elong*deg )
c
      mcato%x = mcafr%x*trans(1,1) + mcafr%y*trans(1,2) + 
     :          mcafr%z*trans(1,3)
      mcato%y = mcafr%x*trans(2,1) + mcafr%y*trans(2,2) + 
     :          mcafr%z*trans(2,3)
      mcato%z = mcafr%x*trans(3,1) + mcafr%y*trans(3,2) + 
     :          mcafr%z*trans(3,3)             
c
      mto%radius  = sqrt(mcato%x**2+mcato%y**2+mcato%z**2)
      if( abs(mcato%z) .lt. mto%radius )then
        mto%colat   = acos(mcato%z/mto%radius)/deg
      else
        mto%colat   = 0
      endif
      if( (mcato%x).eq.0 .and. (mcato%y).eq.0)then
        mto%elong = 0
      else
        mto%elong = atan2(mcato%y,mcato%x)/deg
      endif          
c
      end
C----------------------------------------------------------------------
      subroutine ut556
     :                (moldpos,mvfrom,mvto,trans)
c
c     vector conversion
c
      include 'structure.h'             
cDEC$ IF DEFINED (_x86_)
cDEC$ ATTRIBUTES DLLEXPORT :: UT556
cDEC$ ENDIF
c
c
c     The subr ut556 converts vector components from a coordinate 
c     system to another according to the transformation matrix trans 
c     and the location moldpos in the 'old' coordinate system, both 
c     passed as argument. Since the vector is expressed by its spherical
c     components, the transformation is function of the position where 
c     it is evaluated. This position is passed to the subroutine
c     in the 'old' coordinate system
c
*
*     moldpos  position in the 'old' coordinate system
*     mvfrom   vector components in the 'old' coordinate system
*     mvto     vector components in the 'new' coordinate system
*     trans   
*
C
      TYPE(zgeo) :: moldpos
      TYPE(zvec) :: mvfrom,mvto
      real*8 trans(3,3)
c
c
c
      COMMON /UC160/ pi, deg, re, gmagmo, eclipt, geoid, uma

      REAL*8         pi, deg, re
      REAL*8         gmagmo
      REAL*8         eclipt, geoid(3), uma(30)               
c
c
c
      TYPE(zxyz) :: mvfr,mvtc,mcto
      real*8       colat,elong
      real*8       st,sp,ct,cp
c
c
c
      st     = sin( moldpos%colat*deg )
      ct     = cos( moldpos%colat*deg )
      sp     = sin( moldpos%elong*deg )
      cp     = cos( moldpos%elong*deg )
      mvfr%z = mvfrom%rho * ct   - mvfrom%theta * st
      mvfr%x = mvfrom%rho *st*cp + mvfrom%theta *ct*cp - mvfrom%phi *sp
      mvfr%y = mvfrom%rho *st*sp + mvfrom%theta *ct*sp + mvfrom%phi *cp
c
c
      mcto%x = st*cp*trans(1,1) + st*sp*trans(1,2) + ct*trans(1,3)
      mcto%y = st*cp*trans(2,1) + st*sp*trans(2,2) + ct*trans(2,3)
      mcto%z = st*cp*trans(3,1) + st*sp*trans(3,2) + ct*trans(3,3)
c
      if(mcto%z .gt. 1)mcto%z=1
      if(mcto%z .lt. -1)mcto%z=-1
      colat  = acos(mcto%z)
      if( mcto%x.eq.0 .and. mcto%y.eq.0)then
        elong = 0
      else
        elong = atan2(mcto%y,mcto%x)
      endif          
c
c
      mvtc%x = mvfr%x*trans(1,1) + mvfr%y*trans(1,2) + mvfr%z*trans(1,3)
      mvtc%y = mvfr%x*trans(2,1) + mvfr%y*trans(2,2) + mvfr%z*trans(2,3)
      mvtc%z = mvfr%x*trans(3,1) + mvfr%y*trans(3,2) + mvfr%z*trans(3,3)
c
c
      st     = sin( colat )
      ct     = cos( colat )
      sp     = sin( elong )
      cp     = cos( elong )
      mvto%rho   = mvtc%z *  ct + mvtc%x * st*cp + mvtc%y * st*sp
      mvto%theta = -st * mvtc%z + mvtc%x * ct*cp + mvtc%y * ct*sp
      mvto%phi   =                - sp * mvtc%x  + mvtc%y *  cp
      mvto%dnrm  = sqrt(mvtc%x**2+mvtc%y**2+mvtc%z**2)
c
      end
c----------------------------------------------------------------------