SUBROUTINE read_pump_timeseries(nout,nx,ny,nz,nchem2,lfile,pumpfile,PumpFileFormat)

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
INTEGER(I4B), INTENT(IN)                                   :: nchem2
INTEGER(I4B), INTENT(IN)                                   :: lfile
CHARACTER (LEN=mls), INTENT(IN)                            :: pumpfile
CHARACTER (LEN=mls), INTENT(IN)                            :: PumpFileFormat

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
REAL(DP)                                               :: tpdis
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
REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE      :: check3
REAL(DP), DIMENSION(:,:,:), ALLOCATABLE      :: check4
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
  IF (ssch == 'pumptimeseries') THEN
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF(ls /= 0) THEN
      lzs=ls
      CALL convan(ssch,lzs,res)
      IF (res == 'n') THEN
        tpdis = DNUM(ssch)
      ELSE                !  An ascii string--so bag it.
        WRITE(*,*)
        WRITE(*,*) ' Cant interpret string following "pump"'
        WRITE(*,*) ' Looking for numerical value'
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF

!  Now, look for geochemical condition following pumping rate (only used if rate is positive)

      id = ids + ls
      CALL sschaine(zone,id,iff,ssch,ids,ls)
      IF(ls /= 0) THEN
        lzs=ls
        CALL convan(ssch,lzs,res)



        DO nco = 1,nchem2
          IF (ssch == condlabel(nco)) THEN
            GO TO 50
          END IF
        END DO
        WRITE(*,*)
        WRITE(*,*) ' Geochemical condition for pumping well not found'
        WRITE(*,*) ' Label = ',ssch
        WRITE(*,*)
        READ(*,*)
        STOP
        50         continue
        intbnd_tmp = nco
      ELSE         !  Blank string
        WRITE(*,*)
        WRITE(*,*) ' No geochemical condition for pumping well provided'
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF




    ELSE
      WRITE(*,*)
      WRITE(*,*) ' No pumping rate given'
      WRITE(*,*) ' Pumping zone ignored'
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
WRITE(iunit2,*) '  Reading pump from file: ',pumpfile(1:lfile)
WRITE(iunit2,*)



  IF (ALLOCATED(tpump)) THEN
    DEALLOCATE(tpump)
    ALLOCATE(tpump(1:int(tpdis)))
  ELSE
    ALLOCATE(tpump(1:int(tpdis)))
  END IF

  IF (ALLOCATED(qgt)) THEN
    DEALLOCATE(qgt)
    ALLOCATE(qgt(1:int(tpdis),1:nx,1:ny,1:nz))
  ELSE
    ALLOCATE(qgt(1:int(tpdis),1:nx,1:ny,1:nz))
  END IF

  IF (ALLOCATED(intbnd)) THEN
    DEALLOCATE(intbnd)
    ALLOCATE(intbnd(1,1:nx,1:ny,1:nz))
  ELSE
    ALLOCATE(intbnd(1,1:nx,1:ny,1:nz))
  END IF

  IF (ALLOCATED(npump)) THEN
    DEALLOCATE(npump)
    ALLOCATE(npump(1:nx,1:ny,1:nz))
  ELSE
    ALLOCATE(npump(1:nx,1:ny,1:nz))
  END IF

  IF (ALLOCATED(qg)) THEN
    DEALLOCATE(qg)
    ALLOCATE(qg(1,nx,ny,nz))
  ELSE
    ALLOCATE(qg(1,nx,ny,nz))
  END IF

IF (PumpFileFormat == 'SingleFile3D') THEN

    INQUIRE(FILE=pumpfile,EXIST=ext)
    IF (.NOT. ext) THEN
      WRITE(*,*)
      WRITE(*,*) ' 3D pump file not found: ',pumpfile(1:lfile)
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF

    
    check=tpdis

    OPEN(UNIT=23,FILE=pumpfile,STATUS='old',ERR=8005)
    FileTemp = pumpfile
    CALL stringlen(FileTemp,FileNameLength)
    DO jz = 1,nz
      DO jy = 1,ny
        DO jx = 1,nx
          DO tp = 1,int(tpdis)
          READ(23,*,END=1020) xdum,ydum,zdum,tpump(tp),qgt(tp,jx,jy,jz)
          if (qgt(tp,jx,jy,jz)>-30 .AND. qgt(tp,jx,jy,jz)/=500) then
          qgt(tp,jx,jy,jz)=10**qgt(tp,jx,jy,jz)
          npump(jx,jy,jz)=1
          elseif (qgt(tp,jx,jy,jz)==-30) then
          qgt(tp,jx,jy,jz)=0
          npump(jx,jy,jz)=1
          elseif (qgt(tp,jx,jy,jz)==500) then
          qgt(tp,jx,jy,jz)=0
          npump(jx,jy,jz)=0
          end if
            END DO
            intbnd(1,jx,jy,jz)=intbnd_tmp
            qg(1,nx,ny,nz)=qgt(1,jx,jy,jz)
        END DO
      END DO
    END DO

    check1=tpump
  check2=qgt
check3=intbnd
check4=npump

    CLOSE(UNIT=23,STATUS='keep')
END IF

RETURN

1020  WRITE(*,*) ' End of file during read'
WRITE(*,*) ' Trying to read the file: ', FileTemp(1:FileNameLength)
READ(*,*)
STOP

8005   WRITE(*,*) ' Error opening 3D permeability file: STOP'
        READ(*,*)

END SUBROUTINE read_pump_timeseries