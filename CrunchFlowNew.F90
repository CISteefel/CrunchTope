!  CrunchFlowNew.f90 
!
PROGRAM CrunchFlowNew
!!USE crunchtype
!!USE params
!!USE RunTime, ONLY: InputFile

IMPLICIT NONE

interface
  subroutine CrunchFlow(NumInputfiles,InputFileCounter,NewInput)
    LOGICAL(KIND(.TRUE.))                                          :: NewInput
    INTEGER(SELECTED_INT_KIND(9))                                  :: NumInputFiles
    INTEGER(SELECTED_INT_KIND(9))                                  :: InputFileCounter
  end subroutine CrunchFlow
end interface

LOGICAL(KIND(.TRUE.))                                               :: NewInput
INTEGER(SELECTED_INT_KIND(9))                                       :: NumInputFiles
INTEGER(SELECTED_INT_KIND(9))                                       :: InputFileCounter

NumInputFiles = 1   
InputFileCounter = 1  
NewInput = .TRUE.


DO WHILE (NewInput)
  CALL CrunchFlow(NumInputFiles,InputFileCounter,NewInput)
END DO

STOP
END 