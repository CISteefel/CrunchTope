SUBROUTINE read_boundary_condition_Richards(nout, Richards_steady, BC_location, BC_type, BC_file, value, lfile, constant_BC, time_length)
! This subroutine reads the boundary condition for the Richards equation
USE crunchtype
USE CrunchFunctions
USE params
USE strings

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nout
LOGICAL(LGT), INTENT(IN)                                    :: Richards_steady ! logical to define if the boundary condition is for the steady or time-dependent problem 
INTEGER(I4B), INTENT(IN)                                    :: BC_location ! ingeger to define the location of the boundary condition (0: lower boundary condition; 1: upper boundary condition)
INTEGER(I4B), INTENT(OUT)                                   :: lfile
INTEGER(I4B), INTENT(OUT)                                   :: time_length ! length of the time series to define the size of the array for the time-series
CHARACTER (LEN=mls), INTENT(OUT)                            :: BC_type ! the type of the boundary condition
CHARACTER (LEN=mls), INTENT(OUT)                            :: BC_file ! filename for time-dependent boundary condition
REAL(DP), INTENT(OUT)                                       :: value ! value used in the boundary condition
LOGICAL(LGT), INTENT(IN OUT)                                :: constant_BC ! logical to define if the boundary condition is constant or time-dependent

!  Internal variables and arrays

INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs
INTEGER(I4B)                                                :: nlen1
CHARACTER (LEN=mls)                                         :: BC_string ! the string to match each boundary condition for steady-state and time-dependent problem
BC_type = ' '
BC_file = ' '
value = 0.0d0

! ****************************************************************************************************
! determine BC_string depending on the location of the boundary and the type of the problem (steady or time-dependent)
IF (BC_location == 0) THEN
  IF (Richards_steady) THEN
    BC_string = 'x_begin_bc_type_steady'
  ELSE
    BC_string = 'x_begin_bc_type'
  END IF
ELSE IF (BC_location == 1) THEN
  IF (Richards_steady) THEN
    BC_string = 'x_end_bc_type_steady'
  ELSE
    BC_string = 'x_end_bc_type'
  END IF
ELSE
  WRITE(*,*)
  WRITE(*,*) ' Error in the input value for the BC_location. '
  STOP
END IF
! ****************************************************************************************************

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
  
  IF (ssch == BC_string) THEN
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF(ls /= 0) THEN
      lzs=ls
      CALL convan(ssch,lzs,res)
      IF (res == 'n') THEN
        WRITE(*,*) ' Must provide a character to specify the type of the boudnary condition. '
      ELSE
        ! select the type of boundary condition
        CALL stringtype(ssch,lzs,res)
        BC_type = ssch
        id = ids + ls
        
        ! read further information on the boundary condition
        bc_case: SELECT CASE (BC_type)
          
        CASE ('constant_dirichlet', 'constant_neumann', 'constant_flux') bc_case
          constant_BC = .TRUE.
          CALL sschaine(zone,id,iff,ssch,ids,ls)
          ! obtain the value of the constant boundary condition
          IF(ls /= 0) THEN
            lzs=ls
            CALL convan(ssch,lzs,res)
            IF (res == 'n') THEN
              value = REAL(DNUM(ssch), DP)
            ELSE
            WRITE(*,*)
            WRITE(*,*) ' Must provide a numerical value for the constant boundary condition for ', BC_string
            WRITE(*,*)
            READ(*,*)
            STOP  
            ENDIF
          ELSE
            WRITE(*,*)
            WRITE(*,*) ' No information provided on the constant boundary condition for ', BC_string
            WRITE(*,*)
            READ(*,*)
            STOP  
          ENDIF
        CASE ('variable_dirichlet', 'variable_neumann', 'variable_flux') bc_case
          constant_BC = .FALSE.
          ! variable boundary condition
          CALL sschaine(zone,id,iff,ssch,ids,ls)
          ! obtain the file name for the variable boundary condition
          IF(ls /= 0) THEN
            lzs=ls
            CALL convan(ssch,lzs,res)
            IF (res == 'n') THEN
              WRITE(*,*)
              WRITE(*,*) ' Must provide a filename for the variable boundary condition for ', BC_string
              WRITE(*,*)
              READ(*,*)
            ELSE
              ! obtain the file name for the variable boundary condition
              CALL stringtype(ssch,lzs,res)
              BC_file = ssch
              lfile = ls
              id = ids + ls
              
              ! read further information on the boundary condition
              CALL sschaine(zone,id,iff,ssch,ids,ls)
              IF(ls /= 0) THEN
                lzs=ls
                CALL convan(ssch,lzs,res)
                IF (res == 'n') THEN
                time_length = int(DNUM(ssch))
                ELSE
                  WRITE(*,*)
                  WRITE(*,*) ' Must provide a numerical value for length time series for the boundary condition for ', BC_string
                  WRITE(*,*)
                  READ(*,*)
                  STOP  
                ENDIF
              ELSE
                WRITE(*,*)
                WRITE(*,*) ' Must provide a length for time series for the boundary condition for ', BC_string
                WRITE(*,*)
                READ(*,*)
                STOP  
              ENDIF
            ENDIF
          ELSE
            WRITE(*,*)
            WRITE(*,*) ' No information provided on the filename for the variable boundary condition for ', BC_string
            WRITE(*,*)
            READ(*,*)
            STOP  
          ENDIF
        CASE ('environmental_forcing') bc_case
          constant_BC = .FALSE.
          
        CASE DEFAULT bc_case
          WRITE(*,*)
          WRITE(*,*) ' The boundary condition type ', BC_type, ' is not supported for ', BC_string
          WRITE(*,*)
        END SELECT bc_case  
      END IF
    ELSE
      WRITE(*,*)
      WRITE(*,*) ' No information provided on the boundary condition for ', BC_string
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
END SUBROUTINE read_boundary_condition_Richards