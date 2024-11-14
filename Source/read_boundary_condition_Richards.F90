SUBROUTINE read_boundary_condition_Richards_1D(nout, Richards_steady, BC_location, BC_type, BC_file, value, lfile, constant_BC, time_length)
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
          
        CASE ('constant_dirichlet', 'constant_neumann', 'constant_flux', 'constant_atomosphere') bc_case
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
        CASE ('variable_dirichlet', 'variable_neumann', 'variable_flux', 'variable_atomosphere') bc_case
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
END SUBROUTINE read_boundary_condition_Richards_1D
  
SUBROUTINE read_boundary_condition_Richards(nout,nx,ny,nz)
USE crunchtype
USE params
USE strings
USE Richards_module, ONLY: Richards_BCs
IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nout
INTEGER(I4B), INTENT(IN)                                    :: nx
INTEGER(I4B), INTENT(IN)                                    :: ny
INTEGER(I4B), INTENT(IN)                                    :: nz


!  Internal variables and arrays

INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs
INTEGER(I4B)                                                :: nlen1
INTEGER(I4B)                                                :: nco
INTEGER(I4B)                                                :: ik
INTEGER(I4B)                                                :: kk
INTEGER(I4B)                                                :: k
INTEGER(I4B)                                                :: nex
INTEGER(I4B)                                                :: ns
REAL(DP)                                      :: value ! value used in the boundary condition


REWIND nout

10  READ(nout,'(a)',END=800) zone
nlen1 = LEN(zone)
CALL majuscules(zone,nlen1)
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch_a,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  CALL convan(ssch_a,lzs,res)
  
  SELECT CASE (ssch_a)
    CASE ('x_begin', 'x_end', 'y_begin', 'y_end', 'z_begin', 'z_end')
      !   Look for condition label following the boundary
      id = ids + ls
      CALL sschaine(zone,id,iff,ssch,ids,ls)
      IF(ls /= 0) THEN
        lzs=ls
        CALL convan(ssch,lzs,res)
        IF (res == 'n') THEN
          WRITE(*,*) ' Must provide a character to specify the type of the boudnary condition. '
        ELSE
          ! read further information on the boundary condition
          bc_case: SELECT CASE (ssch)
          
          CASE ('dirichlet', 'neumann', 'flux', 'atomosphere') bc_case
            CALL sschaine(zone,id,iff,ssch,ids,ls)
            ! obtain the value of the constant boundary condition
            IF(ls /= 0) THEN
              lzs=ls
              CALL convan(ssch,lzs,res)
              IF (res == 'n') THEN
                value = REAL(DNUM(ssch), DP)
              ELSE
                WRITE(*,*)
                WRITE(*,*) ' Must provide a numerical value for the constant boundary condition for ', ssch_a
                WRITE(*,*)
                READ(*,*)
                STOP  
              END IF
            ELSE
              WRITE(*,*)
              WRITE(*,*) ' No information provided on the constant boundary condition for ', ssch_a
              WRITE(*,*)
              READ(*,*)
              STOP  
            END IF
                  
          CASE DEFAULT bc_case
            WRITE(*,*)
            WRITE(*,*) ' The boundary condition type ', ssch, ' is not supported for ', ssch_a
            WRITE(*,*)
          END SELECT bc_case  
        END IF
      ELSE
        WRITE(*,*)
        WRITE(*,*) ' Blank string following boundary condition for Richards solver'
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
    
    CASE DEFAULT
      GO TO 10
  END SELECT
  
ELSE
  GO TO 10
END IF

800 CONTINUE
!IF (nx > 1 .AND. jc(1) == -1) THEN
!  WRITE(*,*)
!  WRITE(*,*) ' No boundary conditions found for JX = 1 '
!  WRITE(*,*) 
!  STOP
!END IF
!IF (nx > 1 .AND. jc(2) == -1) THEN
!  WRITE(*,*)
!  WRITE(*,*) ' No boundary conditions found for JX = NX '
!  WRITE(*,*) 
!  STOP
!END IF
!IF (ny > 1 .AND. jc(3) == -1) THEN
!  WRITE(*,*)
!  WRITE(*,*) ' No boundary conditions found for JY = 1 '
!  WRITE(*,*) 
!  STOP
!END IF
!IF (ny > 1 .AND. jc(4) == -1) THEN
!  WRITE(*,*)
!  WRITE(*,*) ' No boundary conditions found for JY = NY '
!  WRITE(*,*) 
!  STOP
!END IF
!IF (nz > 1 .AND. jc(5) == -1) THEN
!  WRITE(*,*)
!  WRITE(*,*) ' No boundary conditions found for JZ = 1 '
!  WRITE(*,*) 
!  STOP
!END IF
!IF (nz > 1 .AND. jc(6) == -1) THEN
!  WRITE(*,*)
!  WRITE(*,*) ' No boundary conditions found for JZ = NZ '
!  WRITE(*,*) 
!  STOP
!END IF
       
END SUBROUTINE read_boundary_condition_Richards
  
  
SUBROUTINE read_boundary_condition_RichardsByZone(nout,nx,ny,nz,nBoundaryConditionZone_Richards)
USE crunchtype
USE CrunchFunctions
USE params
USE strings
USE flow

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nout
INTEGER(I4B), INTENT(IN)                                    :: nx
INTEGER(I4B), INTENT(IN)                                    :: ny
INTEGER(I4B), INTENT(IN)                                    :: nz
INTEGER(I4B), INTENT(OUT)                                   :: nBoundaryConditionZone_Richards

!  Internal variables and arrays

INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs
INTEGER(I4B)                                                :: nxyz
INTEGER(I4B)                                                :: nlen1
INTEGER(I4B)                                                :: ls_a
INTEGER(I4B)                                                :: ls_b
INTEGER(I4B)                                                :: l
INTEGER(I4B)                                                :: ncond

CHARACTER (LEN=mls)                                         :: BC_ConditionName

INTEGER(I4B), PARAMETER                                     :: mBoundaryConditionZone=500

IF (ALLOCATED(BoundaryZone_Richards)) THEN
  DEALLOCATE(BoundaryZone_Richards)
END IF
ALLOCATE(BoundaryZone_Richards(1000))

IF (ALLOCATED(BoundaryValue_Richards)) THEN
  DEALLOCATE(BoundaryValue_Richards)
END IF
ALLOCATE(BoundaryValue_Richards(mBoundaryConditionZone))


IF (ALLOCATED(jxxBC_Richards_lo)) THEN
  DEALLOCATE(jxxBC_Richards_lo)
END IF
ALLOCATE(jxxBC_Richards_lo(mBoundaryConditionZone))

IF (ALLOCATED(jyyBC_Richards_lo)) THEN
  DEALLOCATE(jyyBC_Richards_lo)
END IF
ALLOCATE(jyyBC_Richards_lo(mBoundaryConditionZone))

IF (ALLOCATED(jzzBC_Richards_lo)) THEN
  DEALLOCATE(jzzBC_Richards_lo)
END IF
ALLOCATE(jzzBC_Richards_lo(mBoundaryConditionZone))

IF (ALLOCATED(jxxBC_Richards_hi)) THEN
  DEALLOCATE(jxxBC_Richards_hi)
END IF
ALLOCATE(jxxBC_Richards_hi(mBoundaryConditionZone))

IF (ALLOCATED(jyyBC_Richards_hi)) THEN
  DEALLOCATE(jyyBC_Richards_hi)
END IF
ALLOCATE(jyyBC_Richards_hi(mBoundaryConditionZone))

IF (ALLOCATED(jzzBC_Richards_hi)) THEN
  DEALLOCATE(jzzBC_Richards_hi)
END IF
ALLOCATE(jzzBC_Richards_hi(mBoundaryConditionZone))

nxyz = nx*ny*nz

REWIND nout

nBoundaryConditionZone_Richards = 0

DO nCond = 1,mBoundaryConditionZone

  READ(nout,'(a)',END=500) zone

  nlen1 = LEN(zone)
  CALL majuscules(zone,nlen1)
  id = 1
  iff = mls
  CALL sschaine(zone,id,iff,ssch,ids,ls)

  lzs=ls
  CALL convan(ssch,lzs,res)

  IF (ssch == 'boundarycondition' .OR. ssch == 'BoundaryCondition' .OR. ssch == 'BOUNDARYCONDITION') THEN
    
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF(ls /= 0) THEN
      lzs=ls
      CALL convan(ssch,lzs,res)
      IF (res == 'a') THEN
        IF (ssch == 'zone') THEN
            
!  "Zone" specified, so look for locations
            
          nBoundaryConditionZone_Richards = nBoundaryConditionZone_Richards + 1

          IF (nBoundaryConditionZone_Richards > mBoundaryConditionZone) THEN

            WRITE(*,*)
            WRITE(*,*)  ' Number of BoundaryCondition zones dimensioned too small'
            WRITE(*,*)  ' Number of BoundaryCondition zones = ', nBoundaryConditionZone_Richards
            WRITE(*,*)  ' Dimension of BoundaryCondition zones = ', mBoundaryConditionZone
            WRITE(*,*)  ' Contact C.I. Steefel at Berkeley Lab: "CISteefel@lbl.gov"'
            WRITE(*,*)
            READ(*,*)
            STOP

          END IF
                        
          id = ids + ls
          CALL sschaine_hyph(zone,id,iff,ssch_a,ssch_b,ids,ls_a,ls_b,ls)
          IF(ls /= 0) THEN
            lzs=ls_a
            CALL convan(ssch_a,lzs,res)
            IF (res == 'n') THEN
              jxxBC_Richards_lo(nBoundaryConditionZone_Richards) = JNUM(ssch_a)
            ELSE                !  An ascii string--so bag it.
              WRITE(*,*)
              WRITE(*,*) ' A grid location should follow zone specification'
              WRITE(*,*) ' Dont know what to do with this string'
              WRITE(*,*)
              READ(*,*)
              STOP
            END IF
            IF (ls_b /= 0) THEN
              lzs=ls_b
              CALL convan(ssch_b,lzs,res)
              IF (res == 'n') THEN
                jxxBC_Richards_hi(nBoundaryConditionZone_Richards) = JNUM(ssch_b)
              ELSE                !  An ascii string--so bag it.
                WRITE(*,*)
                WRITE(*,*) ' A grid location should follow zone specification'
                WRITE(*,*) ' Dont know what to do with this string after "boundarycondition"'
                WRITE(*,*)
                READ(*,*)
                STOP
              END IF
            ELSE
              jxxBC_Richards_hi(nBoundaryConditionZone_Richards) = jxxBC_Richards_lo(nBoundaryConditionZone_Richards) 
            END IF
          ELSE                  ! Zero length trailing string
            WRITE(*,*)
            WRITE(*,*) ' No X or Y grid location given for BoundaryCondition'
            WRITE(*,*) ' BoundaryCondition zone ',nBoundaryConditionZone_Richards
            WRITE(*,*)
            READ(*,*)
            STOP
          END IF
            
          id = ids + ls
          CALL sschaine_hyph(zone,id,iff,ssch_a,ssch_b,ids,ls_a,ls_b,ls)
          IF(ls /= 0) THEN
            lzs=ls_a
            CALL convan(ssch_a,lzs,res)
            IF (res == 'n') THEN
              jyyBC_Richards_lo(nBoundaryConditionZone_Richards) = JNUM(ssch_a)
            ELSE                !  An ascii string--so bag it.
              WRITE(*,*)
              WRITE(*,*) ' No Y location for BoundaryCondition '
              WRITE(*,*)
              READ(*,*)
              STOP
            END IF
            IF (ls_b /= 0) THEN
              lzs=ls_b
              CALL convan(ssch_b,lzs,res)
              IF (res == 'n') THEN
                jyyBC_Richards_hi(nBoundaryConditionZone_Richards) = JNUM(ssch_b)
              ELSE                !  An ascii string--so bag it.
                WRITE(*,*)
                WRITE(*,*) ' A grid location should follow zone specification'
                WRITE(*,*) ' Dont know what to do with this string after "BoundaryCondition"'
                WRITE(*,*)
                READ(*,*)
                STOP
              END IF
            ELSE
              jyyBC_Richards_hi(nBoundaryConditionZone_Richards) = jyyBC_Richards_lo(nBoundaryConditionZone_Richards)  
            END IF
          ELSE                  ! Zero length trailing string
            WRITE(*,*)
            WRITE(*,*) ' No Y location given for BoundaryCondition zone'
            WRITE(*,*) ' BoundaryCondition zone number ',nBoundaryConditionZone_Richards
            WRITE(*,*)
            READ(*,*)
            STOP
          END IF    
            
          id = ids + ls
          CALL sschaine_hyph(zone,id,iff,ssch_a,ssch_b,ids,ls_a,ls_b,ls)
          IF(ls /= 0) THEN
            lzs=ls_a
            CALL convan(ssch_a,lzs,res)
            IF (res == 'n') THEN
              jzzBC_Richards_lo(nBoundaryConditionZone_Richards) = JNUM(ssch_a)
            ELSE                !  An ascii string--so bag it.
              WRITE(*,*)
              WRITE(*,*) ' No Z location for BoundaryCondition '
              WRITE(*,*)
              READ(*,*)
              STOP
            END IF
            IF (ls_b /= 0) THEN
              lzs=ls_b
              CALL convan(ssch_b,lzs,res)
              IF (res == 'n') THEN
                jzzBC_Richards_hi(nBoundaryConditionZone_Richards) = JNUM(ssch_b)
              ELSE                !  An ascii string--so bag it.
                WRITE(*,*)
                WRITE(*,*) ' A grid location should follow zone specification'
                WRITE(*,*) ' Dont know what to do with this string after "BoundaryCondition"'
                WRITE(*,*)
                READ(*,*)
                STOP
              END IF
            ELSE
              jzzBC_Richards_hi(nBoundaryConditionZone_Richards) = jzzBC_Richards_lo(nBoundaryConditionZone_Richards)   !  Assume jxxpermx_hi=jxxpermx_lo
            END IF
          ELSE                  ! Zero length trailing string
            WRITE(*,*)
            WRITE(*,*) ' No Z location given for BoundaryCondition zone'
            WRITE(*,*) ' BoundaryCondition zone number ',nBoundaryConditionZone_Richards
            WRITE(*,*)
            READ(*,*)
            STOP
          END IF    

          ELSE
            WRITE(*,*)
            WRITE(*,*) ' Dont understand string following BoundaryCondition specification'
            WRITE(*,*)   ssch(1:ls)
            WRITE(*,*)
            READ(*,*)
            STOP
          END IF
          
        ELSE                !  A number--so bag it.
          WRITE(*,*)
          WRITE(*,*) ' Cant interpret string following BoundaryCondition value'
          WRITE(*,*) ' Looking for an ASCII string'
          WRITE(*,*)
          READ(*,*)
          STOP
        END IF

      ELSE

        WRITE(*,*)
        WRITE(*,*) ' Boundary condition is not provided for the 2D Richards solver '
        WRITE(*,*)
        READ(*,*)
        STOP

      END IF
      
      !   Look for the type of boundary condition (Dirichlet, flux, etc.)
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF (ls /= 0) THEN
      lzs=ls
      CALL convan(ssch,lzs,res)
      IF (ssch == 'dirichlet') THEN
        BoundaryZone_Richards(nBoundaryConditionZone_Richards) = 1
      ELSE IF (ssch == 'neumann') THEN
        BoundaryZone_Richards(nBoundaryConditionZone_Richards) = 2
      ELSE IF (ssch == 'flux') THEN
        BoundaryZone_Richards(nBoundaryConditionZone_Richards) = 3
      ELSE
        WRITE(*,*)
        WRITE(*,*) ' Dont recognize this type of boundary condition'
        WRITE(*,*) ' Trying to read ByGrid', nBoundaryConditionZone_Richards
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
    ELSE
      WRITE(*,*)
      WRITE(*,*) ' No boundary condtiion type is provided for the Richards solver '
      READ(*,*)
      STOP
    END IF
    
    ! look for numerical value for the boundary condition
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF(ls /= 0) THEN
      lzs=ls
      CALL convan(ssch,lzs,res)
      IF (res == 'n') THEN
        BoundaryValue_Richards(nBoundaryConditionZone_Richards) = DNUM(ssch)
      ELSE                !  An ascii string--so bag it.
        WRITE(*,*)
        WRITE(*,*) ' Cant interpret string following the zone boundary condition'
        WRITE(*,*) ' Looking for numerical value'
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
    ELSE
      WRITE(*,*)
      WRITE(*,*) ' No value given for a boundary condition for the Richards solver '
      READ(*,*)
      STOP
    END IF
    

  END IF
    

  

END DO

500 DO l = 1,nBoundaryConditionZone_Richards
  IF (jxxBC_Richards_hi(l) > nx+1) THEN
    WRITE(*,*)
    WRITE(*,*) 'You have specified a BoundaryCondition at JX > NX+1'
    WRITE(*,*)
    STOP
  END IF
  IF (jyyBC_Richards_hi(l) > ny+1) THEN
    WRITE(*,*)
    WRITE(*,*) 'You have specified a BoundaryCondition at JY > NY'
    WRITE(*,*)
    STOP
  END IF
  IF (jzzBC_Richards_hi(l) > nz+1) THEN
    WRITE(*,*)
    WRITE(*,*) 'You have specified a BoundaryCondition at JZ > NZ'
    WRITE(*,*)
    STOP
  END IF
  IF (jxxBC_Richards_lo(l) < 0) THEN
    WRITE(*,*)
    WRITE(*,*) 'You have specified a BoundaryCondition at JX < 0'
    WRITE(*,*)
    STOP
  END IF
  IF (jyyBC_Richards_lo(l) < 0) THEN
    WRITE(*,*)
    WRITE(*,*) 'You have specified a BoundaryCondition at JY < 1'
    WRITE(*,*)
    STOP
  END IF
  IF (jzzBC_Richards_lo(l) < 0) THEN
    WRITE(*,*)
    WRITE(*,*) 'You have specified a BoundaryCondition at JZ < 1'
    STOP
  END IF
END DO

RETURN
END SUBROUTINE read_boundary_condition_RichardsByZone