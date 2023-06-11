SUBROUTINE read_infiltration_2(nout,nx,ny,nz,infiltration_file,rate,lfile,boolfix,booltimeseries,tslength)
USE crunchtype
USE CrunchFunctions
USE params
USE flow
USE strings
USE medium

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nout
INTEGER(I4B), INTENT(IN)                                    :: nx
INTEGER(I4B), INTENT(IN)                                    :: ny
INTEGER(I4B), INTENT(IN)                                    :: nz
INTEGER(I4B), INTENT(OUT)                                   :: lfile
INTEGER(I4B), INTENT(OUT)                                   :: tslength
CHARACTER (LEN=mls), INTENT(OUT)                            :: infiltration_file
REAL(DP), INTENT(OUT)                                       :: rate
LOGICAL(LGT), INTENT(IN OUT)                                :: boolfix
LOGICAL(LGT), INTENT(IN OUT)                                :: booltimeseries

!  Internal variables and arrays

INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs
INTEGER(I4B)                                                :: nxyz
INTEGER(I4B)                                                :: nlen1
INTEGER(I4B)                                                :: lenformat

nxyz = nx*ny*nz
infiltration_file = ' '
rate = 0.0d0

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
  
  IF (ssch == 'infiltration_forcing') THEN
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF(ls /= 0) THEN
        lzs=ls
      !!  *****************************************************************     
      !!check for numerical values
      CALL convan(ssch,lzs,res)
      IF (res == 'n') THEN
      boolfix = .true.
      rate = DNUM(ssch)
      !!  *****************************************************************
      !!check for the name of timeseries file
      ELSEIF (res /= 'n') THEN
        CALL stringtype(ssch,lzs,res)
        infiltration_file = ssch
        lfile = ls
        booltimeseries = .true.
        id = ids + ls
      CALL sschaine(zone,id,iff,ssch,ids,ls)
        IF(ls /= 0) THEN
        lzs=ls
        CALL convan(ssch,lzs,res)
        IF (res == 'n') THEN
        tslength = int(DNUM(ssch))
        ELSE
          WRITE(*,*)
          WRITE(*,*) ' Must provide a numerical value for length time series infiltration '
          WRITE(*,*)
          READ(*,*)
          STOP  
        ENDIF
      ELSE
        WRITE(*,*)
        WRITE(*,*) ' Must provide a length for time series infiltration '
        WRITE(*,*)
        READ(*,*)
        STOP  
        ENDIF

      ENDIF
  !!  *****************************************************************

      ELSE
        WRITE(*,*)
        WRITE(*,*) ' No info provided on infiltration '
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
      
  ELSE
    GO TO 10
  END IF
  GO TO 10
  
ELSE         ! No string found
  GO TO 10
END IF

1000 RETURN
END SUBROUTINE read_infiltration_2
