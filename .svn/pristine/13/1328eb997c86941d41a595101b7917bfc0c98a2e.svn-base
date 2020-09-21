      PROGRAM testirbem
c........+.........+.........+.........+.........+.........+.........+..
c
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-B)
      IMPLICIT REAL*8 (D-H)
      IMPLICIT REAL*8 (O-Z)

C
      INTEGER*4    options(5)
      real*8       maginput(25,1), lm(1,1), lstar(1,1) ,xIN(3), xOUT(3)
      character (len=56) :: arg
C
      data options/1,30,0,0,0/

      CALL get_command_argument(1, arg)
      read(arg,*) ky
      CALL get_command_argument(2, arg)
      read(arg,*) kd
      CALL get_command_argument(3, arg)
      read(arg,*) ut
      ut = ut * 3600.D0
      CALL get_command_argument(4, arg)
      read(arg,*) kp
      CALL get_command_argument(5, arg)
      read(arg,*) rlat
      CALL get_command_argument(6, arg)
      read(arg,*) rlon
      CALL get_command_argument(7, arg)
      read(arg,*) alt
C

c     calculation of L and L*
      print *, "Tests L & L* calculations:"
      maginput(1,1) = kp
      call make_lstar1(1,4,options,0,ky,kd,ut, alt,
     *              rlat,rlon,maginput, lm,lstar,blocal,bmin,xj,mlt)
      print *,'KP = ',KP,'DOY = ',kd,'L = ',
     &     lm(1,1),'L*=',lstar(1,1)
      end
