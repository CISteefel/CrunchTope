SUBROUTINE read_watertable_timeseries(nout,nx,ny,nz,lfile,watertablefile,WatertableFileFormat)

USE crunchtype
USE CrunchFunctions
USE params
USE concentration
USE transport
USE flow
USE strings
USE RunTime
USE io

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                   :: nout
INTEGER(I4B), INTENT(IN)                                   :: nx
INTEGER(I4B), INTENT(IN)                                   :: ny
INTEGER(I4B), INTENT(IN)                                   :: nz
INTEGER(I4B), INTENT(IN)                                   :: lfile
CHARACTER (LEN=mls), INTENT(IN)                            :: watertablefile
CHARACTER (LEN=mls), INTENT(IN)                            :: WatertableFileFormat

!  Internal variables and arrays

INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs
INTEGER(I4B)                                                :: nlen1
INTEGER(I4B)                                                :: nco
INTEGER(I4B)                                                :: intbnd_tmp
INTEGER(I4B)                                                :: jx
INTEGER(I4B)                                                :: jy
INTEGER(I4B)                                                :: jz
REAL(DP)                                                    :: tpdis
CHARACTER (LEN=mls)                                           :: FileTemp
REAL(DP)                                                      :: xdum
REAL(DP)                                                      :: ydum
REAL(DP)                                                      :: zdum
REAL(DP)                                                      :: tottime
REAL(DP)                                                      :: check
INTEGER(I4B)                                                   :: tp
INTEGER(I4B)                                                  :: FileNameLength
LOGICAL(LGT)                                               :: ext

REAL(DP), DIMENSION(:), ALLOCATABLE             :: check1
REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE      :: check2
REAL(DP), DIMENSION(:,:,:), ALLOCATABLE      :: check3
REAL(DP), DIMENSION(:,:,:), ALLOCATABLE      :: check4
REAL(DP), DIMENSION(:,:,:), ALLOCATABLE      :: check5
REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE      :: qgdummy
REWIND nout

!! Time discretization (in years) and geochemical conditions of the wells

10 READ(nout,'(a)',END=500) zone
nlen1 = LEN(zone)
CALL majuscules(zone,nlen1)
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)

IF(ls /= 0) THEN
  lzs=ls
  CALL convan(ssch,lzs,res)
  IF (ssch == 'watertabletimeseries') THEN
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF(ls /= 0) THEN
      lzs=ls
      CALL convan(ssch,lzs,res)
      IF (res == 'n') THEN
        tpdis = DNUM(ssch)
      ELSE                !  An ascii string--so bag it.
        WRITE(*,*)
        WRITE(*,*) ' Cant interpret string following "watertableseries"'
        WRITE(*,*) ' Looking for numerical value' !should provide the length of the time series
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF




    ELSE
      WRITE(*,*)
      WRITE(*,*) ' No time series length given'
      WRITE(*,*) ' Time series ignored'
      WRITE(*,*)
      GO TO 10
    END IF
  ELSE
    GO TO 10
  END IF
ELSE
  GO TO 10
END IF

GO TO 10
500 CONTINUE
!!!!!!!!!!!!!!!!!!!!
WRITE(iunit2,*)
WRITE(iunit2,*) '  Reading watertabletimeseries from file: ',watertablefile(1:lfile)
WRITE(iunit2,*)



  IF (ALLOCATED(twatertable)) THEN
    DEALLOCATE(twatertable)
    ALLOCATE(twatertable(1:int(tpdis)))
  ELSE
    ALLOCATE(twatertable(1:int(tpdis)))
  END IF

  IF (ALLOCATED(pressurebct)) THEN
    DEALLOCATE(pressurebct)
    ALLOCATE(pressurebct(1:int(tpdis),1,1:ny,1)) !! only work for a 2D model XY
  ELSE
    ALLOCATE(pressurebct(1:int(tpdis),1,1:ny,1))
  END IF

 !! IF (ALLOCATED(pres)) THEN
 !!   DEALLOCATE(pres)
  !!  ALLOCATE(pres(0,ny,1))
 !! ELSE
  !!  ALLOCATE(pres(0,ny,1))
 !! END IF

IF (WatertableFileFormat == 'SingleFile3D') THEN

    INQUIRE(FILE=watertablefile,EXIST=ext)
    IF (.NOT. ext) THEN
      WRITE(*,*)
      WRITE(*,*) ' 3D watertable file not found: ',watertablefile(1:lfile)
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF

    
    check=tpdis

    OPEN(UNIT=23,FILE=watertablefile,STATUS='old',ERR=8005)
    FileTemp = watertablefile
    CALL stringlen(FileTemp,FileNameLength)
    DO jz = 1,nz
      DO jy = 1,ny
        DO jx = 0,0
          DO tp = 1,int(tpdis)
          READ(23,*,END=1020) xdum,ydum,zdum,twatertable(tp),pressurebct(tp,jx,jy,jz)
          if (pressurebct(tp,jx,jy,jz)/=0) then
          pressurebct(tp,jx,jy,jz)=(10**pressurebct(tp,jx,jy,jz))
          else
          permx(jx,jy,jz)=0  
          end if
            END DO
            pres(jx,jy,jz)=pressurebct(1,jx,jy,jz)
            activecellPressure(jx,jy,jz) = 0
        END DO
      END DO
    END DO
    check1=twatertable
  check2=pressurebct
check3=pres
check4=activecellPressure
check5=permx
    CLOSE(UNIT=23,STATUS='keep')
END IF

RETURN

1020  WRITE(*,*) ' End of file during read'
WRITE(*,*) ' Trying to read the file: ', FileTemp(1:FileNameLength)
READ(*,*)
STOP

8005   WRITE(*,*) ' Error opening 3D permeability file: STOP'
        READ(*,*)

END SUBROUTINE read_watertable_timeseries