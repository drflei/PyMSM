      program PREMAP
c     program to convert the magcos out file to the AVG file required by RCUT3
c     and add the L parameter for each grid position
c
C........+.........+.........+.........+.........+.........+.........+..
c
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-B)
      IMPLICIT REAL*8 (D-H)
      IMPLICIT REAL*8 (O-Z)

      INTEGER*4    options(5), year, kp, iut
      real*8       maginput(25,1), lm(1,1), lstar(1,1)
c
      real*8 rlat(37,73), rlon(37,73), rc(37,73,3)

      character*256 filein, fileout, fileout1,ckp

      data options/0,30,0,0,0/

c      INTEGER :: i
      CHARACTER(len=56) :: arg
          
      CALL get_command_argument(1, arg)
      read(arg,*) year
      CALL get_command_argument(2, arg)
      ckp = TRIM(arg)
      read(arg,*) kp
      if (kp.eq.10) ckp = "X"
      fileout1 = "AVKP"//ckp
      CALL get_command_argument(3, arg)
      read (arg,*) iut
      fileout = TRIM(fileout1)//"T"//TRIM(arg)//".AVG"
      CALL get_command_argument(4, arg)
      filein = TRIM(arg)
C
      open(1,file=filein, status='old')
c
      read(1,*)
      read(1,*)
      read(1,*)
      read(1,*)
      read(1,*)
c
      do i = 37, 1, -1
         do j = 1,72
            read(1,*) rlat(i,j),rlon(i,j),rc(i,j,1),rc(i,j,2),rc(i,j,3)
         enddo
         rlat(i,73) = rlat(i,72)
         rlon(i,73) = 360.
         rc(i,73,1) = rc(i,1,1)
         rc(i,73,2) = rc(i,1,2)
         rc(i,73,3) = rc(i,1,3)
      enddo
      close(1)
c
      maginput(1,1) = kp*10.D0
      if ( maginput(1,1) .ge. 90.D0) maginput(1,1) = 89.9D0
c      The IRBEM acceptes KP 0-9 (0 - 90) only
      ut = iut*3600.D0 
      write(*,*) year, kp, iut, ut, trim(fileout)
     & , trim(filein)
c
      open(2,file=fileout,status='unknown')
      do i = 1, 37
         do j = 1,73            
            call make_lstar1(1,4,options,0,year,1,ut,450.D0,
     *      rlat(i,j),rlon(i,j),maginput,lm,lstar,blocal,bmin,xj,mlt)
            FL = DABS(lm(1,1))
            if (FL .gt. 100.D0) FL=99.99D0
            write (2,2000) rlat(i,j),rlon(i,j),"450",FL,
     *           rc(i,j,1),rc(i,j,2),rc(i,j,3)
         enddo
      enddo
 2000 FORMAT (1X, F7.2, F7.2,A4,F9.2,F6.2,F6.2,F6.2)
      close(2)
      end
