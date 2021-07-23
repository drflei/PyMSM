      SUBROUTINE myversion ( version )
!
!!    Returns version of unilib based on SVN path
!!       expects svnver = '$'+'URL: '+'<path>/<version>/src/<filename>'+' $'
!
        implicit none
!
        CHARACTER(len=*)    version
!
        CHARACTER(len=*), PARAMETER :: svnver='$URL: https://www.mag-unilib.eu/svn/branches/VALIRENE/src/version.f90 $'
        CHARACTER(len=*), PARAMETER :: undef='*undefined*'
        INTEGER                     :: iversion,isrc,iflnm
!
        version = undef
        iflnm = SCAN( svnver, '/' , back=.TRUE. )
        IF (iflnm .LE. 2) RETURN
        isrc = SCAN( svnver(:iflnm-1), '/' , back=.TRUE. )
        IF (isrc .LT. 8 .OR. svnver(isrc+1:iflnm-1) .NE. 'src' .OR. svnver(:6) .NE. '$URL: ') RETURN
        iversion = SCAN( svnver(:isrc-1), '/' , back=.TRUE. )
        IF (isrc-iversion-1 .GT. len(version) ) RETURN
        version = svnver(iversion+1:isrc-1)

      END SUBROUTINE myversion
