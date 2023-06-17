SUBROUTINE read_vanGenuchten_parameters(nout, lchar, parchar, parfind, section, nx, ny, nz, VG_error)
! This subroutine reads the van Genuchten parameters from the input file
USE crunchtype
USE params
USE flow

IMPLICIT NONE

INTERFACE
  SUBROUTINE read_multpar(nout,lchar,parchar,parfind,realmult,lenarray,section)
  USE crunchtype
  USE params
  USE strings
  IMPLICIT NONE
  INTEGER(I4B), INTENT(IN)                                    :: nout
  INTEGER(I4B), INTENT(OUT)                                   :: lchar
  CHARACTER (LEN=mls), INTENT(IN)                             :: parchar
  CHARACTER (LEN=mls), INTENT(IN OUT)                         :: parfind
  INTEGER(I4B), INTENT(OUT)                                   :: lenarray
  REAl(DP), DIMENSION(:), INTENT(IN OUT)                      :: realmult
  CHARACTER (LEN=mls), INTENT(IN)                             :: section
  END SUBROUTINE read_multpar
END INTERFACE


INTEGER(I4B)                                                 :: i
INTEGER(I4B)                                                 :: ii
INTEGER(I4B)                                                 :: jx
INTEGER(I4B)                                                 :: jxx


INTEGER(I4B), INTENT(INOUT)                                  :: lchar ! length of character array
INTEGER(I4B)                                                 :: lenarray ! length of the array in the input file
INTEGER(I4B), INTENT(IN)                                     :: nx
INTEGER(I4B), INTENT(IN)                                     :: ny
INTEGER(I4B), INTENT(IN)                                     :: nz
INTEGER(I4B), INTENT(IN)                                     :: nout ! unit number for output file
INTEGER(I4B)                                                 :: nzone ! number of zones for each van Genuchten parameter
INTEGER(I4B), INTENT(INOUT)                                  :: VG_error ! error flag for reading van Genuchten parameters

CHARACTER (LEN=mls), INTENT(INOUT)                           :: parchar ! character for each van Genuchten parameter in the input file
CHARACTER (LEN=mls), INTENT(INOUT)                           :: parfind ! character found in the input file

CHARACTER (LEN=mls), INTENT(IN)                              :: section ! the name of the section
REAL(DP), DIMENSION(:), ALLOCATABLE                          :: realmult
REAL(DP), DIMENSION(:), ALLOCATABLE                          :: ncells_VG ! arrays for the number of cells for each van Genuchten parameter
REAL(DP), DIMENSION(:), ALLOCATABLE                          :: value_VG ! array for value of van Genuchten parameter for each cell

parfind = ' '
ALLOCATE(realmult(1000))
realmult = 0.0
CALL read_multpar(nout,lchar,parchar,parfind,realmult,lenarray,section)

IF (parfind == ' ') THEN

ELSE
  IF (realmult(1) <= 0.0) THEN
    WRITE(*,*) ' The number of cells must be positive. '
    VG_error = 1
  ELSE
!  Check to see if there are an even number values given
!  One for the number of cells, the second for parameter value
    IF (MOD(lenarray, 2) == 0) THEN
      nzone = lenarray/2
      
      ALLOCATE(ncells_VG(nzone))
      ALLOCATE(value_VG(nzone))
      
      DO i = 1, lenarray, 2
        ii = (i+1)/2
        ncells_VG(ii) = realmult(i)
        value_VG(ii) = realmult(i+1)
      END DO
      
      SELECT CASE (parchar)
      CASE ('vg_theta_r')
        IF (ALLOCATED(theta_r)) THEN
          DEALLOCATE(theta_r)
          ALLOCATE(theta_r(nx, ny, nz))
        ELSE
          ALLOCATE(theta_r(nx, ny, nz))
        END IF
        
        jxx = 0
        DO i = 1, nzone
          DO jx = 1, ncells_VG(i)
            jxx = jxx + 1
            theta_r(jxx, ny, nz) = value_VG(i)
          END DO
        END DO
        
        
      CASE ('vg_theta_s')
        IF (ALLOCATED(theta_s)) THEN
          DEALLOCATE(theta_s)
          ALLOCATE(theta_s(nx, ny, nz))
        ELSE
          ALLOCATE(theta_s(nx, ny, nz))
        END IF
        
        jxx = 0
        DO i = 1, nzone
          DO jx = 1, ncells_VG(i)
            jxx = jxx + 1
            theta_s(jxx, ny, nz) = value_VG(i)
          END DO
        END DO
        
        
      CASE ('vg_alpha')
        IF (ALLOCATED(VG_alpha)) THEN
          DEALLOCATE(VG_alpha)
          ALLOCATE(VG_alpha(nx, ny, nz))
        ELSE
          ALLOCATE(VG_alpha(nx, ny, nz))
        END IF
        
        jxx = 0
        DO i = 1, nzone
          DO jx = 1, ncells_VG(i)
            jxx = jxx + 1
            VG_alpha(jxx, ny, nz) = value_VG(i)
          END DO
        END DO
      
      CASE ('vg_n')
        IF (ALLOCATED(VG_n)) THEN
          DEALLOCATE(VG_n)
          ALLOCATE(VG_n(nx, ny, nz))
        ELSE
          ALLOCATE(VG_n(nx, ny, nz))
        END IF
        
        jxx = 0
        DO i = 1, nzone
          DO jx = 1, ncells_VG(i)
            jxx = jxx + 1
            VG_n(jxx, ny, nz) = value_VG(i)
          END DO
        END DO
      
      ! this is for modified van Genuchten model        
      !CASE ('vg_psi_s')
      !  IF (ALLOCATED(psi_s)) THEN
      !    DEALLOCATE(psi_s)
      !    ALLOCATE(psi_s(nx, ny, nz))
      !  ELSE
      !    ALLOCATE(psi_s(nx, ny, nz))
      !  END IF
      !  
      !  jxx = 0
      !  DO i = 1, nzone
      !    DO jx = 1, ncells_VG(i)
      !      jxx = jxx + 1
      !      psi_s(jxx, ny, nz) = value_VG(i)
      !    END DO
      !  END DO
        
      CASE DEFAULT
        WRITE(*,*) ' The parameter could not be allocated and inserted. '
        VG_error = 1
      END SELECT

    ELSE
      WRITE(*, *) ' Both cell number and parameter value must be specified. '
      VG_error = 1
    END IF
  END IF
END IF
      
END SUBROUTINE read_vanGenuchten_parameters