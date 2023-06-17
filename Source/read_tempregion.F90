SUBROUTINE read_tempregion(nout,nx,ny,nz,lfile,filename,fileformat)

USE crunchtype
USE medium
USE CrunchFunctions
USE params
USE temperature
USE transport
USE flow
USE strings
USE io

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                   :: nout
INTEGER(I4B), INTENT(IN)                                   :: nx
INTEGER(I4B), INTENT(IN)                                   :: ny
INTEGER(I4B), INTENT(IN)                                   :: nz
CHARACTER (LEN=mls), INTENT(IN)                            :: filename
CHARACTER (LEN=mls), INTENT(IN)                            :: fileformat
INTEGER(I4B), INTENT(IN)                                   :: lfile

!  Internal variables and arrays

INTEGER(I4B)                                                :: jx
INTEGER(I4B)                                                :: jy
INTEGER(I4B)                                                :: jz
LOGICAL(LGT)                                               :: ext
REAL(DP)     :: dummy1
CHARACTER (LEN=mls)                                           :: FileTemp
INTEGER(I4B)                                                  :: FileNameLength
REWIND nout

!!!!!!!!!!!!!!!!!!!!
WRITE(iunit2,*)
WRITE(iunit2,*) '  Reading temp zone from file: ',filename(1:lfile)
WRITE(iunit2,*)





  IF (ALLOCATED(temp_region)) THEN
    DEALLOCATE(temp_region)
    ALLOCATE(temp_region(1:nx,1:ny,1:nz))
  ELSE
    ALLOCATE(temp_region(1:nx,1:ny,1:nz))
  END IF

  !! READ PUMP LOCATIONS:

IF (fileformat == 'SingleColumn') THEN

  INQUIRE(FILE=filename,EXIST=ext)
  IF (.NOT. ext) THEN
    WRITE(*,*)
    WRITE(*,*) ' temperature zone file not found: ',filename(1:lfile)
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF

  OPEN(UNIT=23,FILE=filename,STATUS='old',ERR=8005)
  FileTemp = filename
    CALL stringlen(FileTemp,FileNameLength)
    DO jz = 1,nz
      DO jy = 1,ny
    DO jx = 1,nx
  READ(23,*,END=1020) temp_region(jx,jy,jz)
  !READ(23,*,END=1020) dummy1
      END DO
    END DO
END DO
  CLOSE(UNIT=23,STATUS='keep')

ELSE
  WRITE(*,*) ''
  WRITE(*,*) ' Can only read Singlecolumnfile for temperature region'
  WRITE(*,*) ''
  STOP
END IF



RETURN

1020  WRITE(*,*) ' End of file during read'
WRITE(*,*) ' Trying to read the file: ', filename(1:lfile)
READ(*,*)
STOP

8005   WRITE(*,*) ' Error opening temp region file: STOP'
        READ(*,*)
        STOP
END SUBROUTINE read_tempregion