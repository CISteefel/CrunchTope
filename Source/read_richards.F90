MODULE read_richards_module

USE crunchtype
USE Richards_module, ONLY: Richards_Base, RichardsBC, RichardsVariableBC, Richards_Variable_BC, Richards_Options

IMPLICIT NONE

PRIVATE

PUBLIC RichardsReadBoundaryCondition!, &

! Interface block for the external subroutine
INTERFACE
  SUBROUTINE convan(ssch,ls,res)
    USE crunchtype
    USE params
    CHARACTER (LEN=mls), INTENT(IN OUT) :: ssch
    CHARACTER (LEN=1), INTENT(OUT) :: res
    INTEGER(I4B), INTENT(IN OUT) :: ls
  END SUBROUTINE convan

  SUBROUTINE majuscules(zone,nbc)
    USE crunchtype
    USE params
    INTEGER(I4B), INTENT(IN) :: nbc
    CHARACTER (LEN=mls), INTENT(IN OUT) :: zone
  END SUBROUTINE majuscules

  SUBROUTINE sschaine(zone,id,iff,ssch,ids,ls)
    USE crunchtype
    USE params
    CHARACTER (LEN=mls), INTENT(OUT) :: ssch
    CHARACTER (LEN=mls), INTENT(IN) :: zone
    INTEGER(I4B), INTENT(IN) :: id
    INTEGER(I4B), INTENT(IN) :: iff
    INTEGER(I4B), INTENT(OUT) :: ids
    INTEGER(I4B), INTENT(OUT) :: ls
  END SUBROUTINE sschaine

  SUBROUTINE sschaine_hyph(zone,id,iff,ssch_a,ssch_b,ids,ls_a,ls_b,ls)
    USE crunchtype
    USE params
    CHARACTER (LEN=mls), INTENT(OUT) :: ssch_a
    CHARACTER (LEN=mls), INTENT(OUT) :: ssch_b
    CHARACTER (LEN=mls), INTENT(IN)  :: zone
    INTEGER(I4B), INTENT(IN) :: id
    INTEGER(I4B), INTENT(IN) :: iff
    INTEGER(I4B), INTENT(OUT) :: ids
    INTEGER(I4B), INTENT(OUT) :: ls
    INTEGER(I4B), INTENT(OUT) :: ls_a
    INTEGER(I4B), INTENT(OUT) :: ls_b
  END SUBROUTINE sschaine_hyph

  SUBROUTINE stringtype(ssch,ls,res)
    USE crunchtype
    USE params
    CHARACTER (LEN=mls), INTENT(IN OUT) :: ssch
    CHARACTER (LEN=1), INTENT(OUT) :: res
    INTEGER(I4B), INTENT(IN OUT) :: ls
  END SUBROUTINE stringtype

END INTERFACE

  CONTAINS

! ************************************************************************** !
SUBROUTINE InsertBoundaryCondition(nout,nx,ny,nz,BCs, jx_low, jx_high, jy_low, jy_high, jz_low, jz_high, &
                                   BC_type, value, is_variable, variable_BC_index, error)
USE crunchtype
USE params, ONLY: mls, pressure_air, rho_water

IMPLICIT NONE

!  External variables and arrays
INTEGER(I4B), INTENT(IN) :: nout, nx, ny, nz
INTEGER(I4B), INTENT(IN) :: jx_low, jx_high, jy_low, jy_high, jz_low, jz_high
INTEGER(I4B), INTENT(IN) :: BC_type
INTEGER(I4B), INTENT(IN) :: variable_BC_index
LOGICAL(LGT), INTENT(IN) :: is_variable
REAL(DP), INTENT(IN) :: value ! value used in the boundary condition
INTEGER(I4B), INTENT(OUT) :: error
INTEGER(I4B) :: i, jx, jy
TYPE(RichardsBC), POINTER :: BCs(:)

error = 0

! check the bounds
IF (jx_high > nx+1) THEN
  WRITE(*,*)
  WRITE(*,*) 'You have specified a BoundaryCondition at JX > NX+1'
  WRITE(*,*)
  error = 1
  RETURN
END IF
IF (jy_high > ny+1) THEN
  WRITE(*,*)
  WRITE(*,*) 'You have specified a BoundaryCondition at JY > NY'
  WRITE(*,*)
  error = 1
  RETURN
END IF
IF (jz_high > nz+1) THEN
  WRITE(*,*)
  WRITE(*,*) 'You have specified a BoundaryCondition at JZ > NZ'
  WRITE(*,*)
  error = 1
  RETURN
END IF
IF (jx_low < 0) THEN
  WRITE(*,*)
  WRITE(*,*) 'You have specified a BoundaryCondition at JX < 0'
  WRITE(*,*)
  error = 1
  RETURN
END IF
IF (jy_low < 0) THEN
  WRITE(*,*)
  WRITE(*,*) 'You have specified a BoundaryCondition at JY < 1'
  WRITE(*,*)
  error = 1
  RETURN
END IF
IF (jz_low < 0) THEN
  WRITE(*,*)
  WRITE(*,*) 'You have specified a BoundaryCondition at JZ < 1'
  error = 1
  RETURN
END IF

IF (nx > 1 .AND. ny == 1 .AND. nz == 1) THEN ! one-dimensional problem
  IF (jx_low == 0) THEN ! left boundary
    i = 1 
  ELSE IF (jx_low == nx + 1) THEN ! right boundary
    i = 2
  ELSE
    WRITE(*,*)
    WRITE(*,*) ' Invlid boundary condition location for 1D Richards solver '
    WRITE(*,*)
    error = 1
    RETURN
  END IF

  BCs(i)%BC_type = BC_type
  BCs(i)%is_variable = is_variable
  IF (is_variable) THEN
    BCs(i)%variable_BC_index = variable_BC_index
  ELSE
    BCs(i)%BC_value = value
  END IF

ELSE IF (nx > 1 .AND. ny > 1 .AND. nz == 1) THEN ! two-dimensional problem

  IF (Richards_Base%spatial_domain == 'regular') THEN  
    DO jy = jy_low, jy_high
      DO jx = jx_low, jx_high

        IF (jx == 0) THEN ! left boundary
          i = 2*nx+2*ny-jy+1
        ELSE IF (jx == nx+1) THEN ! right boundary
          i = nx+jy
        ELSE IF (jy == 0) THEN ! bottom boundary
          i = jx
        ELSE IF (jy == ny+1) THEN ! top boundary
          i = 2*nx+ny - jx + 1
        ELSE
          WRITE(*,*)
          WRITE(*,*) ' Invlid boundary condition location for 2D Richards solver for regular spatial domain '
          WRITE(*,*)
          error = 1
          RETURN
        END IF

        BCs(i)%BC_type = BC_type
        BCs(i)%is_variable = is_variable
        IF (is_variable) THEN
          BCs(i)%variable_BC_index = variable_BC_index
        ELSE
          BCs(i)%BC_value = value
        END IF

      END DO
    END DO

  ELSE
    WRITE(*,*)
    WRITE(*,*) ' Currently, two-dimensional Richards solver does not support the shape ', Richards_Base%spatial_domain
    WRITE(*,*)
    error = 1
    RETURN
  END IF

ELSE IF (nx > 1 .AND. ny > 1 .AND. nz > 1) THEN
  WRITE(*,*)
  WRITE(*,*) ' Currently, three-dimensional Richards solver is supported.'
  WRITE(*,*)
  error = 1
  RETURN
END IF

END SUBROUTINE InsertBoundaryCondition

! ************************************************************************** !
SUBROUTINE RichardsReadBoundaryCondition(nout, nx, ny, nz, dist_scale,time_scale,BCs, is_steady, error)
USE crunchtype
IMPLICIT NONE

INTEGER(I4B), INTENT(IN) :: nout, nx, ny, nz
INTEGER(I4B), INTENT(OUT) :: error
REAL(DP), INTENT(IN) :: dist_scale, time_scale
INTEGER(I4B) :: nBoundaryConditionZone_Richards = 0
TYPE(RichardsBC), POINTER :: BCs(:)
LOGICAL(LGT), INTENT(IN) :: is_steady

IF (.NOT. ASSOCIATED(BCs)) THEN
  error = 1
ELSE
  CALL ReadBoundaryConditionByZone(nout,nx,ny,nz,dist_scale,time_scale,BCs,is_steady,nBoundaryConditionZone_Richards, error)
  IF (error == 2) THEN ! error in ReadBoundaryConditionByZone
    RETURN
  ELSE
    error = 0
  END IF

END IF

END SUBROUTINE RichardsReadBoundaryCondition
! ************************************************************************** !
SUBROUTINE ReadBoundaryConditionByZone(nout,nx,ny,nz,dist_scale,time_scale,BCs,is_steady,nBoundaryConditionZone_Richards, error)
USE crunchtype
USE params, ONLY: mls, pressure_air, rho_water

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN) :: nout, nx, ny, nz
REAL(DP), INTENT(IN) :: dist_scale,time_scale
TYPE(RichardsBC), POINTER :: BCs(:)
LOGICAL(LGT), INTENT(IN) :: is_steady
INTEGER(I4B), INTENT(INOUT) :: nBoundaryConditionZone_Richards
INTEGER(I4B), INTENT(OUT) :: error

!  Internal variables and arrays

INTEGER(I4B) :: i, id, iff, ids, ls, lzs, nlen1,ls_a,ls_b
INTEGER(I4B) :: jx_low, jx_high, jy_low, jy_high, jz_low, jz_high
CHARACTER (LEN=1) :: res
CHARACTER (LEN=mls) :: zone
CHARACTER (LEN=mls) :: ssch,ssch_a,ssch_b
CHARACTER (LEN=mls) :: BC_name
CHARACTER (LEN=mls) :: keyword
INTEGER(I4B) :: mBoundaryConditionZone = 1000
INTEGER(I4B) :: BC_type
LOGICAL(LGT) :: is_variable
REAL(DP) :: value = 0.0 ! value used in the boundary condition
INTEGER(I4B) :: lfile
INTEGER(I4B) :: variable_BC_index = 0
INTEGER(I4B) :: time_length ! length of the time series to define the size of the array for the time-series
CHARACTER (LEN=mls) :: BC_file ! filename for time-dependent boundary condition
TYPE(RichardsVariableBC), POINTER :: ptr ! temporary pointer
TYPE(RichardsVariableBC), POINTER :: tail ! pointer to tail of the list
REAL(DP), ALLOCATABLE :: BC_time(:)
REAL(DP), ALLOCATABLE :: BC_values(:)
INTEGER(I4B) :: error2

REWIND nout

! change the keyword name for steady state Richards solver
IF (is_steady) THEN
  keyword = 'boundary_condition_steady'
ELSE
  keyword = 'boundary_condition'
END IF

!search keyword
DO i = 1,mBoundaryConditionZone

  READ(nout,'(a)',END=500) zone

  nlen1 = LEN(zone)
  CALL majuscules(zone,nlen1)
  id = 1
  iff = mls
  CALL sschaine(zone,id,iff,ssch,ids,ls)
  lzs=ls
  CALL convan(ssch,lzs,res)

  IF (ssch == keyword) THEN
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF(ls == 0) THEN
      WRITE(*,*)
      WRITE(*,*) ' No boundary condtiion information is provided after ', keyword
      error = 2
      RETURN
    ELSE
      lzs=ls
      CALL convan(ssch,lzs,res)
      IF (ssch /= 'zone') THEN
        WRITE(*,*)
        WRITE(*,*) ' Dont understand string following BoundaryCondition specification'
        WRITE(*,*)   ssch(1:ls)
        WRITE(*,*)
        error = 2
        RETURN
      ELSE
! "Zone" specified, so look for locations
        nBoundaryConditionZone_Richards = nBoundaryConditionZone_Richards + 1

        IF (nBoundaryConditionZone_Richards > mBoundaryConditionZone) THEN

          WRITE(*,*)
          WRITE(*,*)  ' Number of BoundaryCondition zones dimensioned too small'
          WRITE(*,*)  ' Number of BoundaryCondition zones = ', nBoundaryConditionZone_Richards
          WRITE(*,*)  ' Dimension of BoundaryCondition zones = ', mBoundaryConditionZone
          WRITE(*,*)  ' Contact C.I. Steefel at Berkeley Lab: "CISteefel@lbl.gov"'
          WRITE(*,*)
          error = 2
          RETURN

        END IF

        id = ids + ls

        CALL ReadHyphenBlcok(zone, id, jx_low, jx_high, jy_low, jy_high, jz_low, jz_high, error)

! Look for the type of boundary condition (Dirichlet, flux, etc.)

        CALL sschaine(zone,id,iff,BC_name,ids,ls)
        IF (ls /= 0) THEN
          lzs=ls
          CALL convan(BC_name,lzs,res)
          bc_case: SELECT CASE (BC_name)
          CASE ('dirichlet', 'neumann', 'flux', 'atomosphere') bc_case
            is_variable = .FALSE.
            id = ids + ls
            CALL sschaine(zone,id,iff,ssch,ids,ls)
! obtain the value of the constant boundary condition
            IF(ls /= 0) THEN
              lzs=ls
              CALL convan(ssch,lzs,res)
              IF (res == 'n') THEN
                value = REAL(DNUM(ssch), DP)
              ELSE
                WRITE(*,*)
                WRITE(*,*) ' Must provide a numerical value for the constant boundary condition for ', nBoundaryConditionZone_Richards
                WRITE(*,*)
                error = 2
                RETURN  
              END IF
            ELSE
              WRITE(*,*)
              WRITE(*,*) ' No information provided on the constant boundary condition for ', nBoundaryConditionZone_Richards
              WRITE(*,*)
              error = 2
              RETURN 
            ENDIF

            IF (BC_name == 'dirichlet') THEN
              BC_type = 1
              value = value/dist_scale
              IF (.NOT. Richards_Options%psi_is_head) THEN
                  value = (value - pressure_air)/(rho_water*9.80665d0)
              END IF
            ELSE IF (BC_name == 'neumann') THEN
              BC_type = 2
            ELSE IF (BC_name == 'flux') THEN
              BC_type = 3
              value = value/(dist_scale * time_scale)
            ELSE IF (BC_name == 'atomosphere') THEN
              BC_type = 4
              value = value/(dist_scale * time_scale)
            END IF

          CASE ('variable_dirichlet', 'variable_neumann', 'variable_flux', 'variable_atomosphere') bc_case
            is_variable = .TRUE.

! variable boundary condition
            CALL sschaine(zone,id,iff,ssch,ids,ls)
! obtain the file name for the variable boundary condition
            IF(ls /= 0) THEN
              lzs=ls
              CALL convan(ssch,lzs,res)
              IF (res == 'n') THEN
                WRITE(*,*)
                WRITE(*,*) ' Must provide a filename for the variable boundary condition for ', nBoundaryConditionZone_Richards
                error = 2
                RETURN
              ELSE

! obtain the file name for the variable boundary condition
                CALL stringtype(ssch,lzs,res)
                BC_file = ssch
                lfile = ls
                id = ids + ls
                variable_BC_index = variable_BC_index + 1

! read further information on the boundary condition
                CALL sschaine(zone,id,iff,ssch,ids,ls)
                IF(ls /= 0) THEN
                  lzs=ls
                  CALL convan(ssch,lzs,res)
                  IF (res == 'n') THEN
                  time_length = int(DNUM(ssch))
                  ELSE
                    WRITE(*,*)
                    WRITE(*,*) ' Must provide a numerical value for length time series for the boundary condition for ', nBoundaryConditionZone_Richards
                    WRITE(*,*)
                    error = 2
                    RETURN
                  END IF
                ELSE
                  WRITE(*,*)
                  WRITE(*,*) ' Must provide a length for time series for the boundary condition for ', nBoundaryConditionZone_Richards
                  WRITE(*,*)
                  error = 2
                  RETURN
                END IF
              END IF
            ELSE
              WRITE(*,*)
              WRITE(*,*) ' No information provided on the filename for the variable boundary condition for ', nBoundaryConditionZone_Richards
              WRITE(*,*)
              error = 2
              RETURN
            END IF

! read the file and store the values in an array

            IF (ALLOCATED(BC_time)) THEN
              DEALLOCATE(BC_time)
            END IF

            IF (ALLOCATED(BC_values)) THEN
              DEALLOCATE(BC_values)
            END IF

            ALLOCATE(BC_time(time_length))
            ALLOCATE(BC_values(time_length))

            CALL read_timeseries(nout, nx, ny, nz, BC_time, BC_values, lfile, BC_file, time_length)

            IF (BC_name == 'dirichlet') THEN
              BC_type = 1
              BC_values = BC_values/dist_scale
              IF (.NOT. Richards_Options%psi_is_head) THEN
                  BC_values = (BC_values - pressure_air)/(rho_water*9.80665d0)
              END IF
            ELSE IF (BC_name == 'neumann') THEN
              BC_type = 2
            ELSE IF (BC_name == 'flux') THEN
              BC_type = 3
              BC_values = BC_values/(dist_scale * time_scale)
            ELSE IF (BC_name == 'atomosphere') THEN
              BC_type = 4
              BC_values = BC_values/(dist_scale * time_scale)
            END IF

! store the data
            IF (.NOT. ASSOCIATED(Richards_Variable_BC)) THEN
              ALLOCATE(Richards_Variable_BC)
              tail => Richards_Variable_BC
            ELSE
              ALLOCATE(tail%p)
              tail => tail%p
            END IF

            NULLIFY(tail%p)
            tail%BC_time = BC_time
            tail%BC_values = BC_values

          CASE DEFAULT bc_case
            WRITE(*,*)
            WRITE(*,*) ' Dont recognize this type of boundary condition'
            WRITE(*,*) ' Trying to read ByGrid', nBoundaryConditionZone_Richards
            WRITE(*,*)
            error = 2
            RETURN
          END SELECT bc_case

        ELSE
            WRITE(*,*)
            WRITE(*,*) ' No boundary condtiion type is provided for the Richards solver '
            error = 2
            RETURN
        END IF

      END IF
    END IF

! store the extracted values to boundary condition container for each face
    CALL InsertBoundaryCondition(nout,nx,ny,nz,BCs, jx_low, jx_high, jy_low, jy_high, jz_low, jz_high, &
                                     BC_type, value, is_variable, variable_BC_index, error2)

    IF (error2 /= 0) THEN
      WRITE(*,*)
      WRITE(*,*) ' Error when inserting boundary condition values for ', nBoundaryConditionZone_Richards
      error = 2
      RETURN
    END IF
  END IF
END DO    

500 Continue

END SUBROUTINE ReadBoundaryConditionByZone

! ************************************************************************** !

SUBROUTINE ReadHyphenBlcok(zone, id, jx_low, jx_high, jy_low, jy_high, jz_low, jz_high, error)
USE crunchtype
USE params, ONLY: mls, pressure_air, rho_water

IMPLICIT NONE

!  External variables and arrays
INTEGER(I4B), INTENT(OUT) :: jx_low, jx_high, jy_low, jy_high, jz_low, jz_high
INTEGER(I4B), INTENT(OUT) :: error
INTEGER(I4B), INTENT(INOUT) :: id
CHARACTER (LEN=mls), INTENT(IN) :: zone
!  Internal variables and arrays

INTEGER(I4B) :: iff, ids, ls, lzs, ls_a,ls_b
CHARACTER (LEN=1) :: res
CHARACTER (LEN=mls) :: ssch_a,ssch_b

iff = mls
CALL sschaine_hyph(zone,id,iff,ssch_a,ssch_b,ids,ls_a,ls_b,ls)
IF (ls /= 0) THEN
  lzs=ls_a
  CALL convan(ssch_a,lzs,res)
  IF (res == 'n') THEN
    jx_low = JNUM(ssch_a)
  ELSE !  An ascii string--so bag it.
    WRITE(*,*)
    WRITE(*,*) ' A grid location should follow zone specification'
    WRITE(*,*) ' Dont know what to do with this string'
    WRITE(*,*)
    error = 2
    RETURN
  END IF
  IF (ls_b /= 0) THEN
    lzs=ls_b
    CALL convan(ssch_b,lzs,res)
    IF (res == 'n') THEN
      jx_high = JNUM(ssch_b)
    ELSE !  An ascii string--so bag it.
      WRITE(*,*)
      WRITE(*,*) ' A grid location should follow zone specification'
      WRITE(*,*) ' Dont know what to do with this string after "boundary_condition"'
      WRITE(*,*)
      error = 2
      RETURN
    END IF
  ELSE
    jx_high = jx_low
  END IF
ELSE ! Zero length trailing string
  WRITE(*,*)
  WRITE(*,*) ' No X location given for BoundaryCondition'
  WRITE(*,*) ' BoundaryCondition zone '
  WRITE(*,*)
  error = 2
  RETURN
END IF

! y coordinate
id = ids + ls
CALL sschaine_hyph(zone,id,iff,ssch_a,ssch_b,ids,ls_a,ls_b,ls)
IF (ls /= 0) THEN
  lzs=ls_a
  CALL convan(ssch_a,lzs,res)
  IF (res == 'n') THEN
    jy_low = JNUM(ssch_a)
  ELSE !  An ascii string--so bag it.
    WRITE(*,*)
    WRITE(*,*) ' A grid location should follow zone specification'
    WRITE(*,*) ' Dont know what to do with this string'
    WRITE(*,*)
    error = 2
    RETURN
  END IF
  IF (ls_b /= 0) THEN
    lzs=ls_b
    CALL convan(ssch_b,lzs,res)
    IF (res == 'n') THEN
      jy_high = JNUM(ssch_b)
    ELSE !  An ascii string--so bag it.
      WRITE(*,*)
      WRITE(*,*) ' A grid location should follow zone specification'
      WRITE(*,*) ' Dont know what to do with this string after "boundary_condition"'
      WRITE(*,*)
      error = 2
      RETURN
    END IF
  ELSE
    jy_high = jy_low
  END IF
ELSE ! Zero length trailing string
  WRITE(*,*)
  WRITE(*,*) ' No Y grid location given for BoundaryCondition'
  WRITE(*,*) ' BoundaryCondition zone '
  WRITE(*,*)
  error = 2
  RETURN
END IF

! z coordinate
id = ids + ls
CALL sschaine_hyph(zone,id,iff,ssch_a,ssch_b,ids,ls_a,ls_b,ls)
IF (ls /= 0) THEN
  lzs=ls_a
  CALL convan(ssch_a,lzs,res)
  IF (res == 'n') THEN
    jz_low = JNUM(ssch_a)
  ELSE !  An ascii string--so bag it.
    WRITE(*,*)
    WRITE(*,*) ' A grid location should follow zone specification'
    WRITE(*,*) ' Dont know what to do with this string'
    WRITE(*,*)
    error = 2
    RETURN
  END IF
  IF (ls_b /= 0) THEN
    lzs=ls_b
    CALL convan(ssch_b,lzs,res)
    IF (res == 'n') THEN
      jz_high = JNUM(ssch_b)
    ELSE !  An ascii string--so bag it.
      WRITE(*,*)
      WRITE(*,*) ' A grid location should follow zone specification'
      WRITE(*,*) ' Dont know what to do with this string after "boundary_condition"'
      WRITE(*,*)
      error = 2
      RETURN
    END IF
  ELSE
    jz_high = jz_low
  END IF
ELSE ! Zero length trailing string
  WRITE(*,*)
  WRITE(*,*) ' No Z grid location given for BoundaryCondition'
  WRITE(*,*) ' BoundaryCondition zone '
  WRITE(*,*)
  error = 2
  RETURN
END IF

id = ids + ls

END SUBROUTINE ReadHyphenBlcok

END MODULE read_richards_module