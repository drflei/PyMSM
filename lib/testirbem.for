      PROGRAM testirbem
c........+.........+.........+.........+.........+.........+.........+..
c
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-B)
      IMPLICIT REAL*8 (D-H)
      IMPLICIT REAL*8 (O-Z)

C
      INTEGER*4    options(5), DOY
      real*8       maginput(25,1), lm(1,1), lstar(1,1) ,xIN(3), xOUT(3)
C
      data options/1,30,0,0,0/

      FKP = 0.; LYEAR=2000
c     calculation of L and L*
      print *, "Tests L & L* calculations:"
      DOY = 1
      UT = 0.
      FKP = 0.
      maginput(1,1) = FKP
      call make_lstar1(1,4,options,0,LYEAR,DOY,UT, 450.D0,
     *              60.D0,100.D0,maginput, lm,lstar,blocal,bmin,xj,mlt)
      print *,'KP = ',FKP,'DOY = ',DOY,'UT= ',UT,'L = ',
     &     lm(1,1),'L*=',lstar(1,1)
      FKP = 3.
      maginput(1,1) = FKP
      call make_lstar1(1,4,options,0,LYEAR,DOY,UT, 450.D0,
     *              60.D0,100.D0,maginput, lm,lstar,blocal,bmin,xj,mlt)
      print *,'KP = ',FKP,'DOY = ',DOY,'UT= ',UT,'L = ',
     &     lm(1,1),'L*=',lstar(1,1)
      FKP = 9.
      maginput(1,1) = FKP
      call make_lstar1(1,4,options,0,LYEAR,DOY,UT, 450.D0,
     *              60.D0,100.D0,maginput, lm,lstar,blocal,bmin,xj,mlt)
      print *,'KP = ',FKP,'DOY = ',DOY,'UT= ',UT,'L = ',
     &     lm(1,1),'L*=',lstar(1,1)
      FKP = 10.
      maginput(1,1) = FKP
      call make_lstar1(1,4,options,0,LYEAR,DOY,UT, 450.D0,
     *              60.D0,100.D0,maginput, lm,lstar,blocal,bmin,xj,mlt)
      print *,'KP = ',FKP,'DOY = ',DOY,'UT= ',UT,'L = ',
     &     lm(1,1),'L*=',lstar(1,1)
      FKP = 40.
      maginput(1,1) = FKP
      call make_lstar1(1,4,options,0,LYEAR,DOY,UT, 450.D0,
     *              60.D0,100.D0,maginput, lm,lstar,blocal,bmin,xj,mlt)
      print *,'KP = ',FKP,'DOY = ',DOY,'UT= ',UT,'L = ',
     &     lm(1,1),'L*=',lstar(1,1)
      FKP = 50.
      maginput(1,1) = FKP
      call make_lstar1(1,4,options,0,LYEAR,DOY,UT, 450.D0,
     *              60.D0,100.D0,maginput, lm,lstar,blocal,bmin,xj,mlt)
      print *,'KP = ',FKP,'DOY = ',DOY,'UT= ',UT,'L = ',
     &     lm(1,1),'L*=',lstar(1,1)
      FKP = 60.
      maginput(1,1) = FKP
      call make_lstar1(1,4,options,0,LYEAR,DOY,UT, 450.D0,
     *              60.D0,100.D0,maginput, lm,lstar,blocal,bmin,xj,mlt)
      print *,'KP = ',FKP,'DOY = ',DOY,'UT= ',UT,'L = ',
     &     lm(1,1),'L*=',lstar(1,1)
      FKP = 70.
      maginput(1,1) = FKP
      call make_lstar1(1,4,options,0,LYEAR,DOY,UT, 450.D0,
     *              60.D0,100.D0,maginput, lm,lstar,blocal,bmin,xj,mlt)
      print *,'KP = ',FKP,'DOY = ',DOY,'UT= ',UT,'L = ',
     &     lm(1,1),'L*=',lstar(1,1)
      FKP = 80.
      maginput(1,1) = FKP
      call make_lstar1(1,4,options,0,LYEAR,DOY,UT, 450.D0,
     *              60.D0,100.D0,maginput, lm,lstar,blocal,bmin,xj,mlt)
      print *,'KP = ',FKP,'DOY = ',DOY,'UT= ',UT,'L = ',
     &     lm(1,1),'L*=',lstar(1,1)

      DOY = 1
      UT = 0.
      FKP = 90.
      maginput(1,1) = FKP
      call make_lstar1(1,4,options,0,LYEAR,DOY,UT, 450.D0,
     *              60.D0,100.D0,maginput, lm,lstar,blocal,bmin,xj,mlt)
      print *,'KP = ',FKP,'DOY = ',DOY,'UT= ',UT,'L = ',
     &     lm(1,1),'L*=',lstar(1,1)
      DOY = 100
      UT = 0.
      FKP = 0
      maginput(1,1) = FKP
      call make_lstar1(1,4,options,0,LYEAR,DOY,UT, 450.D0,
     *              60.D0,100.D0,maginput, lm,lstar,blocal,bmin,xj,mlt)
      print *,'KP = ',FKP,'DOY = ',DOY,'UT= ',UT,'L = ',
     &     lm(1,1),'L*=',lstar(1,1)
      DOY = 1
      UT = 12.*3600.
      FKP = 0
      maginput(1,1) = FKP
      call make_lstar1(1,4,options,0,LYEAR,DOY,UT, 450.D0,
     *              60.D0,100.D0,maginput, lm,lstar,blocal,bmin,xj,mlt)
      print *,'KP = ',FKP,'DOY = ',DOY,'UT= ',UT,'L = ',
     &     lm(1,1),'L*=',lstar(1,1)
      print *, " new IGRF calculated per day"
      options(2) = 0
      DOY = 1
      UT = 0.
      call make_lstar1(1,4,options,0,LYEAR,DOY,UT, 450.D0,
     *              60.D0,260.D0,maginput, lm,lstar,blocal,bmin,xj,mlt)
      print *,'KP = ',FKP,'DOY = ',DOY,'UT= ',UT,'L = ',
     &     lm(1,1),'L*=',lstar(1,1)

      call make_lstar1(1,4,options,0,LYEAR,DOY,UT, 450.D0,
     *              60.D0,270.D0,maginput, lm,lstar,blocal,bmin,xj,mlt)
      print *,'KP = ',FKP,'DOY = ',DOY,'UT= ',UT,'L = ',
     &     lm(1,1),'L*=',lstar(1,1)

      call make_lstar1(1,4,options,0,LYEAR,DOY,UT, 450.D0,
     *              60.D0,280.D0,maginput, lm,lstar,blocal,bmin,xj,mlt)
      print *,'KP = ',FKP,'DOY = ',DOY,'UT= ',UT,'L = ',
     &     lm(1,1),'L*=',lstar(1,1)

c
c          call coord_trans_vec1(ntime,sysaxesIN,sysaxesOUT,iyr,idoy,secs,xIN,xOUT)
      xIN(1) = 0.D0      ! on the surface
      xIN(2) = 0.D0
      xIN(3) = 100.D0
c     GEO->MAG 
      LYEAR = 1995
      call coord_trans_vec1(1,0,6,LYEAR,1,0.D0,xIN,xOUT)
      call CAR_SPH(xOUT,r,flati,flongi)
c     flati is the magnetic latitude in degrees. Note this is not the same as 
c     the corrected geomagnetic latitude, but the difference should be small
      print *, r, flati, flongi
      xIN(1)= 450.D0  ! 450 km altitude
      call coord_trans_vec1(1,0,6,LYEAR,1,0.D0,xIN,xOUT)
      call CAR_SPH(xOUT,r,flati,flongi)
      print *, r, flati, flongi
     
      end
