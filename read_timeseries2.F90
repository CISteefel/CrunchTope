SUBROUTINE read_timeseries2(nout,nx,ny,nz,tout,qtout,lfile,tsfile,tslength)

USE crunchtype
USE medium
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
INTEGER(I4B), INTENT(IN)                                   :: tslength
CHARACTER (LEN=mls), INTENT(IN)                            :: tsfile
REAL(DP),  INTENT(IN OUT)            :: tout(tslength)
REAL(DP) ,  INTENT(IN OUT)            :: qtout(tslength)
!  Internal variables and arrays

INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs
INTEGER(I4B)                                                :: nlen1
INTEGER(I4B)                                                :: jx
INTEGER(I4B)                                                :: jy
INTEGER(I4B)                                                :: jz
REAL(DP)                                               :: tpdis
CHARACTER (LEN=mls)                                           :: FileTemp
REAL(DP)                                                      :: check
INTEGER(I4B)                                                   :: tp
INTEGER(I4B)                                                  :: FileNameLength
LOGICAL(LGT)                                               :: ext
INTEGER(I4B)                                                     :: ierr
INTEGER(I4B)                                                  :: lenarray
REAL(DP)                                                      :: t_dum
REAL(DP)                                                      :: q_dum
REAL(DP), ALLOCATABLE, DIMENSION(:)           :: tout_dum
REAL(DP) , ALLOCATABLE, DIMENSION(:)           :: qtout_dum
REWIND nout


WRITE(iunit2,*)
WRITE(iunit2,*) '  Reading time series from file: ',tsfile(1:lfile)
WRITE(iunit2,*)



  ! IF (ALLOCATED(tout)) THEN
  !   DEALLOCATE(tout)
  !   ALLOCATE(tout(1:8761))
  ! ELSE
  !   ALLOCATE(tout(1:8761))
  ! END IF

  ! IF (ALLOCATED(qtout)) THEN
  !   DEALLOCATE(qtout)
  !   ALLOCATE(qtout(1:8761))
  ! ELSE
  !   ALLOCATE(qtout(1:8761))
  ! END IF
IF (ALLOCATED(tout_dum)) THEN
  DEALLOCATE(tout_dum)
END IF
IF (ALLOCATED(qtout_dum)) THEN
  DEALLOCATE(qtout_dum)
END IF
ALLOCATE(tout_dum(tslength))
ALLOCATE(qtout_dum(tslength))

!! READ PUMP TIMESERIES:

    INQUIRE(FILE=tsfile,EXIST=ext)
    IF (.NOT. ext) THEN
      WRITE(*,*)
      WRITE(*,*) ' time series file not found: ',timeseriesfile(1:lfile)
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF

    
    !check=tpdis

    OPEN(UNIT=23,FILE=tsfile,STATUS='old',ERR=8005)
    FileTemp = tsfile
    CALL stringlen(FileTemp,FileNameLength)
          
      ! DO WHILE (ierr == 0)
      !   READ(23,*,iostat=IERR) t_dum,q_dum
      !     IF (ierr /= 0) THEN
      !     lenarray = lenarray + 1
      !     ENDIF
      !       tout(lenarray) = t_dum
      !       qtout(lenarray) = ((q_dum)/1000.0d0)*dxx(nx)*dzz(nx,ny,nz) !! Converting from mm/year to m3/year
      !     END DO

      DO tp = 1,tslength
        READ(23,*,iostat=IERR) t_dum,q_dum
        tout_dum(tp) = t_dum
        qtout_dum(tp) = q_dum
        !if (ny==1 .AND. nz==1) then !! 1D case
        !qtout_dum(tp) = ((q_dum)/1000.0d0) !! Converting from mm/year to m3/year
        !else
        !qtout_dum(tp) = ((q_dum)/1000.0d0)*dxx(nx)*dzz(nx,ny,nz) !! Converting from mm/year to m3/year for 2D case
        !endif
      END DO
      
    CLOSE(UNIT=23,STATUS='keep')

tout = tout_dum(1:tslength)
qtout = qtout_dum(1:tslength)

RETURN

8005   WRITE(*,*) ' Error opening timeseries file', timeseriesfile(1:lfile)
        READ(*,*)

END SUBROUTINE read_timeseries2