SUBROUTINE read_pump_timeseries2(nout,nx,ny,nz,nchem2,lfile,pumptimeseriesfile,PumptimeseriesFileFormat,lfile2,pumplocationsfile,PumplocationsFileFormat)

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
INTEGER(I4B), INTENT(IN)                                   :: nchem2
INTEGER(I4B), INTENT(IN)                                   :: lfile
CHARACTER (LEN=mls), INTENT(IN)                            :: pumptimeseriesfile
CHARACTER (LEN=mls), INTENT(IN)                            :: PumptimeseriesFileFormat
INTEGER(I4B), INTENT(IN)                                   :: lfile2
CHARACTER (LEN=mls), INTENT(IN)                            :: pumplocationsfile
CHARACTER (LEN=mls), INTENT(IN)                            :: PumplocationsFileFormat

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
REAL(DP)                                               :: locpdis
CHARACTER (LEN=mls)                                           :: FileTemp
REAL(DP)                                                      :: xdum
REAL(DP)                                                      :: ydum
REAL(DP)                                                      :: zdum
REAL(DP)                                                      :: tottime
REAL(DP)                                                      :: check
REAL(DP)                                                      :: dummy
INTEGER(I4B)                                                   :: tp
INTEGER(I4B)                                                   :: locp
INTEGER(I4B)                                                  :: FileNameLength
LOGICAL(LGT)                                               :: ext

REWIND nout

!! Time discretization (in years), nb of pumps and geochemical conditions of the wells

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
        WRITE(*,*) ' Cant interpret string following "pumptimeseries"'
        WRITE(*,*) ' Looking for numerical value'
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF

      id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
      IF(ls /= 0) THEN
        lzs=ls
        CALL convan(ssch,lzs,res)
        IF (res == 'n') THEN
          locpdis = DNUM(ssch)
        ELSE                !  An ascii string--so bag it.
          WRITE(*,*)
          WRITE(*,*) ' Cant interpret string following length pump time series'
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
      WRITE(*,*) ' No length time series given'
      WRITE(*,*) ' Pumping zone ignored'
      WRITE(*,*)
      GO TO 10
    END IF
  ELSE
    WRITE(*,*)
      WRITE(*,*) ' No nb pump given'
      WRITE(*,*)
      STOP
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
WRITE(iunit2,*) '  Reading pump from file: ',pumptimeseriesfile(1:lfile)
WRITE(iunit2,*)



  IF (ALLOCATED(tpump)) THEN
    DEALLOCATE(tpump)
    ALLOCATE(tpump(1:int(tpdis)))
  ELSE
    ALLOCATE(tpump(1:int(tpdis)))
  END IF

  IF (ALLOCATED(qgt)) THEN
    DEALLOCATE(qgt)
    ALLOCATE(qgt(1:int(tpdis)))
  ELSE
    ALLOCATE(qgt(1:int(tpdis)))
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



!! READ PUMP TIMESERIES:


IF (PumptimeseriesFileFormat == 'SingleFile3D') THEN

    INQUIRE(FILE=pumptimeseriesfile,EXIST=ext)
    IF (.NOT. ext) THEN
      WRITE(*,*)
      WRITE(*,*) ' 3D pump file not found: ',pumptimeseriesfile(1:lfile)
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF

    
    check=tpdis

    OPEN(UNIT=23,FILE=pumptimeseriesfile,STATUS='old',ERR=8005)
    FileTemp = pumptimeseriesfile
    CALL stringlen(FileTemp,FileNameLength)
          DO tp = 1,int(tpdis)
          READ(23,*,END=1020) tpump(tp),qgt(tp)
            qgt(tp)=((qgt(tp))/1000.0d0)*dxx(nx)*dzz(jx,jy,jz) !! Converting from mm/year to m3/year
            END DO
            qg(1,nx,ny,nz)=qgt(1)

    CLOSE(UNIT=23,STATUS='keep')
END IF



!! READ PUMP LOCATIONS:


IF (PumplocationsFileFormat == 'SingleFile3D') THEN

  INQUIRE(FILE=pumplocationsfile,EXIST=ext)
  IF (.NOT. ext) THEN
    WRITE(*,*)
    WRITE(*,*) ' 3D pump file not found: ',pumplocationsfile(1:lfile2)
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF

  OPEN(UNIT=23,FILE=pumplocationsfile,STATUS='old',ERR=8005)
  FileTemp = pumplocationsfile
  CALL stringlen(FileTemp,FileNameLength)
  DO jz = 1,nz
    DO jy = 1,ny
      DO jx = 1,nx
        DO locp = 1,int(locpdis)
        READ(23,*,END=1020) xdum,ydum,zdum
        if (jx==xdum .and. jy==ydum .and. jz==zdum) then
          npump(jx,jy,jz)=1
          else
          npump(jx,jy,jz)=0
        end if
          END DO
          intbnd(1,jx,jy,jz)=intbnd_tmp
      END DO
    END DO
  END DO

  CLOSE(UNIT=23,STATUS='keep')
END IF


RETURN

1020  WRITE(*,*) ' End of file during read'
WRITE(*,*) ' Trying to read the file: ', FileTemp(1:FileNameLength)
READ(*,*)
STOP

8005   WRITE(*,*) ' Error opening 3D permeability file: STOP'
        READ(*,*)

END SUBROUTINE read_pump_timeseries2