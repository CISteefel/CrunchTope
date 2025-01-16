MODULE Crunch_Modules
  
  USE crunchtype
  USE hydraulic_function_module

  IMPLICIT NONE

  PUBLIC
  
  CONTAINS

! ************************************************************************** !
  SUBROUTINE AllocateArray(array, lbound_x, ubound_x, lbound_y, ubound_y, lbound_z, ubound_z)
    USE crunchtype
    IMPLICIT NONE

    REAL(DP), ALLOCATABLE, INTENT(INOUT) :: array(:, :, :)
    INTEGER(I4B), INTENT(IN) :: lbound_x, ubound_x, lbound_y, ubound_y, lbound_z, ubound_z

!   Check if the array is allocated, deallocate if necessary, then allocate with the new dimensions
    IF (ALLOCATED(array)) THEN
      DEALLOCATE (array)
      ALLOCATE (array(lbound_x:ubound_x, lbound_y:ubound_y, lbound_z:ubound_z))
    ELSE
      ALLOCATE (array(lbound_x:ubound_x, lbound_y:ubound_y, lbound_z:ubound_z))
    END IF

  END SUBROUTINE AllocateArray

! ************************************************************************** !
    SUBROUTINE AllocateArrayPlus(array, neqn, lbound_x, ubound_x, lbound_y, ubound_y, lbound_z, ubound_z)
    USE crunchtype
    IMPLICIT NONE

    REAL(DP), ALLOCATABLE, INTENT(INOUT) :: array(:, :, :, :)
    INTEGER(I4B), INTENT(IN) :: neqn, lbound_x, ubound_x, lbound_y, ubound_y, lbound_z, ubound_z

!   Check if the array is allocated, deallocate if necessary, then allocate with the new dimensions
    IF (ALLOCATED(array)) THEN
      DEALLOCATE (array)
      ALLOCATE (array(neqn, lbound_x:ubound_x, lbound_y:ubound_y, lbound_z:ubound_z))
    ELSE
      ALLOCATE (array(neqn, lbound_x:ubound_x, lbound_y:ubound_y, lbound_z:ubound_z))
    END IF

    END SUBROUTINE AllocateArrayPlus
    
! ************************************************************************** !
    SUBROUTINE AllocateArrayPlusPlus(array, neqn1, neqn2, lbound_x, ubound_x, lbound_y, ubound_y, lbound_z, ubound_z)
    USE crunchtype
    IMPLICIT NONE

    REAL(DP), ALLOCATABLE, INTENT(INOUT) :: array(:, :, :, :, :)
    INTEGER(I4B), INTENT(IN) :: neqn1, neqn2, lbound_x, ubound_x, lbound_y, ubound_y, lbound_z, ubound_z

!   Check if the array is allocated, deallocate if necessary, then allocate with the new dimensions
    IF (ALLOCATED(array)) THEN
      DEALLOCATE (array)
      ALLOCATE (array(neqn1, neqn2, lbound_x:ubound_x, lbound_y:ubound_y, lbound_z:ubound_z) )
    ELSE
      ALLOCATE (array(neqn1, neqn2, lbound_x:ubound_x, lbound_y:ubound_y, lbound_z:ubound_z) )
    END IF

  END SUBROUTINE AllocateArrayPlusPlus

! ************************************************************************** !
    
  SUBROUTINE AllocateArray_1D(array, lbound_x, ubound_x)
    USE crunchtype
    IMPLICIT NONE

    REAL(DP), ALLOCATABLE, INTENT(INOUT) :: array(:)
    INTEGER(I4B), INTENT(IN) :: lbound_x, ubound_x

!   Check if the array is allocated, deallocate if necessary, then allocate with the new dimensions
    IF (ALLOCATED(array)) THEN
      DEALLOCATE (array)
      ALLOCATE (array(lbound_x:ubound_x))
    ELSE
      ALLOCATE (array(lbound_x:ubound_x))
    END IF

  END SUBROUTINE AllocateArray_1D

! ************************************************************************** !
  SUBROUTINE AllocateArray_2D(array, lbound_x, ubound_x, lbound_y, ubound_y)
    USE crunchtype
    IMPLICIT NONE

    REAL(DP), ALLOCATABLE, INTENT(INOUT) :: array(:, :)
    INTEGER(I4B), INTENT(IN) :: lbound_x, ubound_x, lbound_y, ubound_y

!   Check if the array is allocated, deallocate if necessary, then allocate with the new dimensions
    IF (ALLOCATED(array)) THEN
      DEALLOCATE (array)
      ALLOCATE (array(lbound_x:ubound_x, lbound_y:ubound_y))
    ELSE
      ALLOCATE (array(lbound_x:ubound_x, lbound_y:ubound_y))
    END IF

  END SUBROUTINE AllocateArray_2D
! ************************************************************************** 
  
  END MODULE Crunch_Modules