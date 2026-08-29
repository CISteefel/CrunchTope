SUBROUTINE OpenInputFile(numproc,InputFileCounter,InputFile,InputFileName,scalar)

  USE TYPES_AND_PARAMETERS
  #include "petsc/finclude/petsc.h"
  USE petscmat  
  USE mpi

  IMPLICIT NONE
  
  INTEGER(I4B)                                 :: numproc
  INTEGER(I4B)                                 :: InputFileCounter
  CHARACTER(LEN=max_len_str), DIMENSION(1000)  :: InputFile
  CHARACTER (LEN=max_len_str)                  :: InputFileName
  TYPE(SCALAR_V)                               :: scalar 

  LOGICAL(LGT)                                 :: file1_exists
  LOGICAL(LGT)                                 :: file2_exists
  INTEGER(I4B)                                 :: errcode
  INTEGER(I4B)                                 :: ierr
  
  numproc = 0   !!! No MPI for CrunchTope yet
  
  IF (numproc == 0) THEN
  
    IF (InputFileCounter == 1) THEN

      file1_exists = .FALSE.
      file2_exists = .FALSE.
      INQUIRE(FILE='CrunchControl.in',EXIST=file1_exists)
      INQUIRE(FILE='PestControl.ant',EXIST=file2_exists)
    
      IF (file1_exists) THEN          !!  CrunchControl.in file exists, so read input filename from it rather than prompting user
        OPEN(scalar%iunit1,FILE='CrunchControl.in',STATUS='old',ERR=43)
        READ(scalar%iunit1,'(a)') InputFileName
        CLOSE(scalar%iunit1,STATUS='keep')
    
      ELSE IF (file2_exists) THEN          !!  Pestcontrol.ant file exists, so read input filename from it rather than prompting user
        OPEN(scalar%iunit1,FILE='PestControl.ant',STATUS='old',ERR=42)
        READ(scalar%iunit1,'(a)') InputFileName
        CLOSE(scalar%iunit1,STATUS='keep')

      ELSE !!  No Pestcontrol.ant and no CrunchControl.in file, so prompt user for the file name
        WRITE(*,*)
        WRITE(*,*) ' Type in your input file name'
        READ(*,'(a)') InputFileName
      END IF
  
    ELSE
      InputFileName = InputFile(InputFileCounter)
    END IF 
  
    file1_exists = .FALSE.
  
    INQUIRE(FILE=InputFileName,EXIST=file1_exists)
    IF (.NOT. file1_exists) THEN
      WRITE(*,*) 
      WRITE(*,*) ' Cannot find input file: ', TRIM(ADJUSTL(InputFileName))
      WRITE(*,*)
      READ(*,*)
      !!! CALL MPI_Abort(PETSC_COMM_WORLD, errcode, ierr)
    END IF
    
  END IF  
  
  RETURN

  ! * Error Messages
  42 WRITE(*,*)
  WRITE(*,*) ' Error opening CrunchControl.in'
  WRITE(*,*)
  READ(*,*)
  CALL MPI_Abort(PETSC_COMM_WORLD, errcode, ierr)
  43 WRITE(*,*)
  WRITE(*,*) ' Error opening PestControl.ant'
  WRITE(*,*)
  READ(*,*)
  CALL MPI_Abort(PETSC_COMM_WORLD, errcode, ierr)
  
END SUBROUTINE OpenInputFile