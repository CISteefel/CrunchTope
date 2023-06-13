
SUBROUTINE read_AqueousFluxSeries(nout,nAqueousFluxSeries,nx,ny,nz)
USE crunchtype
USE CrunchFunctions
USE params
USE concentration
USE strings
USE runtime

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nout
INTEGER(I4B), INTENT(OUT)                                   :: nAqueousFluxSeries
INTEGER(I4B), INTENT(IN)                                    :: nx
INTEGER(I4B), INTENT(IN)                                    :: ny
INTEGER(I4B), INTENT(IN)                                    :: nz

!  Internal variables and arrays

INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs
INTEGER(I4B)                                                :: nxyz
INTEGER(I4B)                                                :: nlen1
INTEGER(I4B)                                                :: i
CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE              :: workchar1

nAqueousFluxSeries = 0

IF (ALLOCATED(FluxSeriesFile)) THEN
  DEALLOCATE(FluxSeriesFile)
END IF
ALLOCATE(FluxSeriesFile(1))

nxyz = nx*ny*nz

REWIND nout

10  READ(nout,'(a)',END=1000) zone
nlen1 = LEN(zone)
!!CALL majuscules(zone,nlen1)
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  CALL convan(ssch,lzs,res)
  
!  Check to see if initial substring is "FluxSeries" (could be more than 1)
  
  IF (ssch == 'writeaqueousflux') THEN
    nAqueousFluxSeries = nAqueousFluxSeries + 1
  ELSE
    GO TO 10
  END IF
ELSE         ! No string found
  GO TO 10
END IF

id = ids + ls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  i = size(FluxSeriesFile,1)
  ALLOCATE(workchar1(i))
  workchar1 = FluxSeriesFile
  DEALLOCATE(FluxSeriesFile)
  ALLOCATE(FluxSeriesFile(nAqueousFluxSeries))
  IF(nseries /= 0) FluxSeriesFile(1:FluxSeries-1) = workchar1(1:FluxSeries-1)
  DEALLOCATE(workchar1)
  FluxSeriesFile(FluxSeries) = ssch
ELSE
  IF (FluxSeries == 1) THEN
    FluxSeriesFile(1) = 'Flux.out'
  ELSE
    WRITE(*,*) 
    WRITE(*,*) ' File name for Fluxes required when more than one is used'
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
END IF

IF (nxyz == 1) THEN
  jxAqueousFluxSeries(nAqueousFluxSeries) = 1
  jyAqueousFluxSeries(nAqueousFluxSeries) = 1
  jzAqueousFluxSeries(nAqueousFluxSeries) = 1
  GO TO 500
END IF

id = ids + ls
CALL sschaine_hyph(zone,id,iff,ssch_a,ssch_b,ids,ls_a,ls_b,ls)
IF (ls /= 0) THEN
  lzs=ls_a
  CALL convan(ssch_a,lzs,res)
  IF (res == 'n') THEN
    jxAqueousFluxSeries_lo(nAqueousFluxSeries) = JNUM(ssch_a)
  ELSE                !  An ascii string--so bag it.
    WRITE(*,*)
    WRITE(*,*) ' A grid location should follow zone specification'
    WRITE(*,*) ' Dont know what to do with this string'
    WRITE(*,*)
    STOP
  END IF
  IF (ls_b /= 0) THEN
    lzs=ls_b
    CALL convan(ssch_b,lzs,res)
    IF (res == 'n') THEN
      jxAqueousFluxSeries_hi(nAqueousFluxSeries) = JNUM(ssch_b)
    ELSE                !  An ascii string--so bag it.
      WRITE(*,*)
      WRITE(*,*) ' A grid location should follow zone specification'
      WRITE(*,*) ' Dont know what to do with this string after "permeability"'
      WRITE(*,*)
      STOP
    END IF
  ELSE
    jxAqueousFluxSeries_hi(nAqueousFluxSeries) = jxAqueousFluxSeries_lo(nAqueousFluxSeries)   !  Assume jxAqueousFluxSeries_hi=jxAqueousFluxSeries_lo
  END IF
ELSE                  ! Zero length trailing string
  WRITE(*,*)
  WRITE(*,*) ' No X or Y grid location given for FluxSeries'
  WRITE(*,*) ' FluxSeries zone ',nAqueousFluxSeries
  WRITE(*,*)
  STOP
END IF

!!!  Got to here

IF (ny == 1) THEN
  jyAqueousFluxSeries(nAqueousFluxSeries) = 1
  jzAqueousFluxSeries(nAqueousFluxSeries) = 1
  GO TO 500
END IF

id = ids + ls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  CALL convan(ssch,lzs,res)
  IF (res == 'n') THEN
    jyAqueousFluxSeries(nAqueousFluxSeries) = JNUM(ssch)
  ELSE                !  An ascii string--so bag it.
    WRITE(*,*)
    WRITE(*,*) ' A grid location (Y) should follow "FluxSeries" label'
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
ELSE
  WRITE(*,*)
  WRITE(*,*) ' Y location for FluxSeries must be specified'
  WRITE(*,*)
  READ(*,*)
  STOP
END IF

IF (nz == 1) THEN
  jzAqueousFluxSeries(nAqueousFluxSeries) = 1
  GO TO 500
END IF

id = ids + ls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  CALL convan(ssch,lzs,res)
  IF (res == 'n') THEN
    jzAqueousFluxSeries(nAqueousFluxSeries) = JNUM(ssch)
  ELSE                !  An ascii string--so bag it.
    WRITE(*,*)
    WRITE(*,*) ' A grid location (Z) should follow "FluxSeries" label'
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
ELSE
  WRITE(*,*)
  WRITE(*,*) ' Z location for FluxSeries must be specified'
  WRITE(*,*)
  READ(*,*)
  STOP
END IF 

500 GO TO 10

1000  RETURN
END SUBROUTINE read_AqueousFluxSeries
