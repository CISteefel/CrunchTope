
    
MODULE Easy_Read_and_Store_M    

CONTAINS
! the following subroutines and functions
! * 1 * Is_Char_OK: test if a character is not too fancy, i.e. not eligible in input files.
! * 2 * IsCommentChar: test if a character is a comment character (! or #).
! * 3 * ContainsMinusChar: test if a string contains a minus (-) char. Returns 0 if not, otherwise returns the position of the first occurence in the string
! * 4 * IsHyphenatedRange: test if a string is an hyphenated range of integer, e.g. 1-23 or 45-56
! * 5 * LowerCase: change a string into its lower case version
! * 6 * FindTheEnd: Find the END keyword corresponding to an entry block
! * 7 * IsDigit: test if a character is a digit character (from 0 to 9)
! * 8 * IsUpperCase: test if a character is an uppercase letter (from A to Z)
! * 9 * IsLowerCase: test if a character is an lowercase letter (from a to z)

! * ReadALine
! * CopyALine
! * ReadAFile
! * WriteACleanFile
! * ExtractBlocksFromInputFile
! * ExtractBlocksFromDatabase


!* 1 *************************************************************************************

PURE FUNCTION IsCharOK(dumchar)
  USE TYPES_AND_PARAMETERS  
  IMPLICIT NONE
  ! see ANSI Character Code Chart at https://www.intel.com/content/www/us/en/docs/fortran-compiler/developer-guide-reference/2024-1/ansi-character-codes-chart-windows.html 

  CHARACTER(len=1),INTENT(IN) :: dumchar

  LOGICAL(LGT)                :: IsCharOK
  
  SELECT CASE(ICHAR(dumchar))  
  CASE(32:126)  
    IsCharOK = .TRUE.  
  CASE DEFAULT 
    IsCharOK = .FALSE.  
  END SELECT 
  RETURN
END FUNCTION


!* 2 *************************************************************************************

PURE FUNCTION IsCommentChar(dumchar)
  USE TYPES_AND_PARAMETERS  
  IMPLICIT NONE

  CHARACTER(len=1),INTENT(IN) :: dumchar

  LOGICAL(LGT)                :: IsCommentChar
  
  IsCommentChar = .FALSE. 
  IF ( dumchar == "!" .OR. dumchar == "#" ) IsCommentChar = .TRUE.  

  RETURN
END FUNCTION

!* 3 *************************************************************************************

PURE FUNCTION ContainsMinusChar(dumstr)
  USE TYPES_AND_PARAMETERS  
  IMPLICIT NONE

  CHARACTER(LEN=max_len_str)   ,INTENT(IN)  :: dumstr

  INTEGER(I4B)                         :: ContainsMinusChar  
  INTEGER(I4B)                         :: i
  
  ContainsMinusChar = 0
  DO i = 1, LEN(dumstr)
    IF (dumstr(i:i) == "-") THEN
      ContainsMinusChar = i
      EXIT
    END IF  
  END DO

  RETURN
END FUNCTION

!* 4 *************************************************************************************

PURE FUNCTION IsHyphenatedRange(dumstr)
  USE TYPES_AND_PARAMETERS  
  IMPLICIT NONE

  CHARACTER(LEN=max_len_str)   ,INTENT(IN)  :: dumstr

  LOGICAL(LGT)                         :: IsHyphenatedRange
                                       
  INTEGER(I4B)                         :: i  
  INTEGER(I4B)                         :: posminus
  
  IsHyphenatedRange = .TRUE.
  posminus = ContainsMinusChar(dumstr)
  
  IF ( (posminus > 1) .AND. (posminus < LEN(dumstr)) ) THEN
  
    DO i = 1, posminus - 1
      IF ( ( ICHAR(dumstr(i:i)) < 48 ) .OR. ( ICHAR(dumstr(i:i)) > 57 ) ) IsHyphenatedRange = .FALSE.
    END DO
  
    DO i = posminus + 1, LEN(TRIM(ADJUSTL(dumstr)))
      IF ( ( ICHAR(dumstr(i:i)) < 48 ) .OR. ( ICHAR(dumstr(i:i)) > 57 ) ) IsHyphenatedRange = .FALSE.
    END DO
    
  ELSE
  
    IsHyphenatedRange = .FALSE.
  
  END IF
  
  RETURN
END FUNCTION

!* 5 *************************************************************************************

PURE FUNCTION LowerCase(dumstr)

  USE TYPES_AND_PARAMETERS  
  IMPLICIT NONE

  CHARACTER(LEN=max_len_str),INTENT(IN)    :: dumstr
  CHARACTER(LEN=max_len_str)               :: LowerCase

  CHARACTER(LEN=max_len_str)     :: dumstrcopy

  INTEGER(I4B)                                               :: ic
  INTEGER(I4B)                                               :: i

  
  dumstrcopy = dumstr
  DO i = 1,LEN(dumstr)
    ic = ICHAR(dumstr(i:i))
    IF (ic > 64 .AND. ic < 91) dumstrcopy(i:i) = CHAR(ic+32)
  END DO
  
  LowerCase = dumstrcopy
  
  RETURN
   
END FUNCTION LowerCase


!* 6 *************************************************************************************

SUBROUTINE FindTheEnd(i,j,FileContent,dumstr)

  USE TYPES_AND_PARAMETERS
  #include "petsc/finclude/petsc.h"
  USE petscmat

  IMPLICIT NONE

  INTEGER(I4B)                                 :: i
  INTEGER(I4B)                                 :: j
  TYPE(FILE_CONTENT)                           :: FileContent 
  CHARACTER(LEN=max_len_str)   ,INTENT(IN)          :: dumstr
  INTEGER(I4B)                                 :: ierr
  INTEGER(I4B)                                 :: errcode  
  INTEGER(I4B)                                 :: refline

  refline = i
  j = 0
  DO 
    i = i + 1
    j = j + 1
    IF ( i > FileContent%NbLinesInFileContent ) THEN
      WRITE(*,'(3a,I,a)') " ERROR: could not find the END statement corresponding to the ", dumstr, " block on line ", refline, " in input file."
      WRITE(*,'(a)') " STOP "
      CALL MPI_Abort(PETSC_COMM_WORLD, errcode, ierr)
      READ(*,*)
    END IF
    
    IF ( (.NOT. FileContent%FileLine(i)%Empty) .AND. (.NOT. FileContent%FileLine(i)%Comment) ) THEN
      IF ( LowerCase(FileContent%FileLine(i)%LineWordsAndNumbers(1)%WordOrNumber) == "end" ) THEN
        j = j - 1
        EXIT
      END IF
    END IF
    
  END DO

END SUBROUTINE FindTheEnd

!* 7 *************************************************************************************

PURE FUNCTION IsDigit(dumchar)
  USE TYPES_AND_PARAMETERS  
  IMPLICIT NONE
  ! see ANSI Character Code Chart at https://www.intel.com/content/www/us/en/docs/fortran-compiler/developer-guide-reference/2024-1/ansi-character-codes-chart-windows.html 

  CHARACTER(len=1),INTENT(IN) :: dumchar

  LOGICAL(LGT)                :: IsDigit
  
  SELECT CASE(ICHAR(dumchar))  
  CASE(48:57)  
    IsDigit = .TRUE.  
  CASE DEFAULT 
    IsDigit = .FALSE.  
  END SELECT 
  RETURN
END FUNCTION IsDigit

!* 8 *************************************************************************************

PURE FUNCTION IsUpperCase(dumchar)
  USE TYPES_AND_PARAMETERS  
  IMPLICIT NONE
  ! see ANSI Character Code Chart at https://www.intel.com/content/www/us/en/docs/fortran-compiler/developer-guide-reference/2024-1/ansi-character-codes-chart-windows.html 

  CHARACTER(len=1),INTENT(IN) :: dumchar

  LOGICAL(LGT)                :: IsUpperCase
  
  SELECT CASE(ICHAR(dumchar))  
  CASE(65:90)  
    IsUpperCase = .TRUE.  
  CASE DEFAULT 
    IsUpperCase = .FALSE.  
  END SELECT 
  RETURN
END FUNCTION

!* 9 *************************************************************************************

PURE FUNCTION IsLowerCase(dumchar)
  USE TYPES_AND_PARAMETERS  
  IMPLICIT NONE
  ! see ANSI Character Code Chart at https://www.intel.com/content/www/us/en/docs/fortran-compiler/developer-guide-reference/2024-1/ansi-character-codes-chart-windows.html 

  CHARACTER(len=1),INTENT(IN) :: dumchar

  LOGICAL(LGT)                :: IsLowerCase
  
  SELECT CASE(ICHAR(dumchar))  
  CASE(97:122)  
    IsLowerCase = .TRUE.  
  CASE DEFAULT 
    IsLowerCase = .FALSE.  
  END SELECT 
  RETURN
END FUNCTION


!***************************************************************************************************************************************!
!********************************************             ******************************************************************************!
!******************************************** Read a Line ******************************************************************************!
!********************************************             ******************************************************************************!
!***************************************************************************************************************************************!



SUBROUTINE ReadALine(FileName,numfile,LineToRead,LineJump,LineNumber)

  USE TYPES_AND_PARAMETERS
  
  IMPLICIT NONE
  
  CHARACTER(LEN=max_len_str)     :: FileName
  INTEGER(I4B), INTENT(IN)       :: numfile
  TYPE(LINE_IN_FILE)             :: LineToRead
  INTEGER(I4B)                   :: LineJump
  INTEGER(I4B)                   :: LineNumber
                                 
  INTEGER(I4B), PARAMETER        :: buflen=1024
  CHARACTER(LEN=buflen)          :: buffer
  CHARACTER(LEN=max_len_str*10)  :: line
  CHARACTER(LEN=max_len_str*10)  :: linecopy
  INTEGER(I4B)                   :: ier
  INTEGER(I4B)                   :: last
  INTEGER(I4B)                   :: isize
  CHARACTER(LEN=max_len_str)     :: WordOrNumber
  INTEGER(I4B)                   :: i
  INTEGER(I4B)                   :: j 
  CHARACTER(LEN=1)               :: dumchar
  CHARACTER(LEN=1)               :: dumchar2
  INTEGER(I4B)                   :: dumint
  REAL(DP)                       :: dumfloat
  LOGICAL(LGT)                   :: Comment

  
  IF ( ALLOCATED(LineToRead%LineWordsAndNumbers)) THEN
    DEALLOCATE(LineToRead%LineWordsAndNumbers)
  END IF
  
  ! Read a complete line, which can be continued with \ and store it in a temporary Line variable
  line = ""
  ier = 0
  

  ! read characters from line and append to result
  INFINITE : DO
     READ(numfile,IOSTAT=ier,FMT='(a)',ADVANCE='no',SIZE=isize) buffer
     ! append what was read to result
     IF (isize > 0) line = TRIM(ADJUSTL(line//buffer(:isize)))
     ! if hit EOR, reading of the line is complete unless backslash ends the line
     IF (IS_IOSTAT_EOR(ier)) THEN   
        last = LEN(TRIM(ADJUSTL(line)))
        IF (last /= 0) THEN
           ! if line ends in backslash it is assumed a continued line 
           IF (line(last:last) == '\') THEN
              LineJump = LineJump + 1 ! keep track of the number of lines to be substracted from the total
              ! remove backslash
              linecopy = " "
              linecopy = line(:last-1)
              line = linecopy
              ! continue on and read next line and append to result
              ier = 0
              CYCLE INFINITE
           END IF 
        END IF
        ! hitting end of record is not an error 
        ier = 0
        ! end of reading line
        EXIT INFINITE
     ! end of file 
     ELSE IF (ier /= 0) THEN
        EXIT INFINITE
     END IF
  END DO INFINITE
  
  LineToRead%FullLine = TRIM(ADJUSTL(line))

  ! replace not eligible characters with space characters
  DO i = 1, LEN(line) 
    dumchar = line(i:i)
    IF ( IsCommentChar(dumchar) ) EXIT 
    IF ( .NOT. IsCharOK(dumchar) ) THEN
      line(i:i) = " "
      WRITE(*,'(3a,I,a,I,3a)') " In file ", TRIM(ADJUSTL(filename)), ", on line ", LineNumber, " at position ", i, ", the character """, dumchar, """ is not authorized (e.g. a tab). It has been replaced with a space character"
    END IF 
  END DO
  
  ! remove leading spaces and trailing blank characters
  line = TRIM(ADJUSTL(line))
  
  ! The "clean" line is stored
  LineToRead%CleanLine = line

  Comment = IsCommentChar(line(1:1))
 
  ! Now split the "line" in words and numbers
  IF (line == "") THEN ! The line is empty
  
    LineToRead%Empty = .TRUE.
    LineToRead%Comment = .FALSE.
    LineToRead%NumberOfItems = 0
  ELSE IF ( Comment ) THEN ! The full line is a comment

    LineToRead%Empty = .FALSE.
    LineToRead%Comment = .TRUE.
    LineToRead%NumberOfItems = 1

    ALLOCATE(LineToRead%LineWordsAndNumbers(1))
    
    LineToRead%LineWordsAndNumbers(1)%StrLength = LEN(TRIM(ADJUSTL(line)))
    LineToRead%LineWordsAndNumbers(1)%WordOrNumber = " "
    LineToRead%LineWordsAndNumbers(1)%IsWord = .TRUE.
    LineToRead%LineWordsAndNumbers(1)%IsNumber = .FALSE.
    LineToRead%LineWordsAndNumbers(1)%IsInteger = .FALSE.
    LineToRead%LineWordsAndNumbers(1)%IsComment = .TRUE. 
    
  ELSE ! Something interesting here
  
    LineToRead%Empty = .FALSE.
    LineToRead%Comment = .FALSE.
    
    ! Look for the number of items in the line
    LineToRead%NumberOfItems = 1
    DO i = 2, LEN(TRIM(ADJUSTL(line))) ! "2" because first character is not a space char and not a comment char. There is no trailing space character either.
      dumchar = line(i:i)
      dumchar2 = line(i-1:i-1)
      Comment = IsCommentChar(dumchar)
      IF ( Comment ) THEN
        IF ( dumchar2 /= " " ) LineToRead%NumberOfItems = LineToRead%NumberOfItems + 1 
        EXIT
      ELSE IF ( dumchar  == " " .AND. dumchar2 /= " " ) THEN
        LineToRead%NumberOfItems = LineToRead%NumberOfItems + 1    
      END IF
    END DO
    
    ! Extract each word, number and comment from the line

    ALLOCATE(LineToRead%LineWordsAndNumbers(LineToRead%NumberOfItems))
        
    j = 1
    WordOrNumber = " "
    WordOrNumber(1:1) = line(1:1)
    DO i = 2, LEN(line) ! "2" because first character is not a space char and not a comment char. There is no trailing space character either.
      dumchar = line(i:i)
      Comment = IsCommentChar(dumchar)
      IF ( line(i:i) == " ") THEN
        IF ( line(i-1:i-1) /= " ") THEN
          LineToRead%LineWordsAndNumbers(j)%WordOrNumber = WordOrNumber
          LineToRead%LineWordsAndNumbers(j)%IsComment = .FALSE.
          j = j + 1
          WordOrNumber = " "
        ELSE 
          CONTINUE ! just multiple space characters
        END IF  
      ELSE IF ( Comment ) THEN ! store the item if not empty, get the rest of the line as a comment and exit
        IF ( WordOrNumber /= " ") THEN 
          LineToRead%LineWordsAndNumbers(j)%WordOrNumber = WordOrNumber
          j = j + 1
        END IF  
        LineToRead%LineWordsAndNumbers(j)%WordOrNumber = " "
        LineToRead%LineWordsAndNumbers(j)%WordOrNumber = line(i:LEN(TRIM(ADJUSTL(line))))
        LineToRead%LineWordsAndNumbers(j)%IsComment = .TRUE.
        EXIT       
      ELSE IF (i == LEN(TRIM(ADJUSTL(line)))) THEN ! we are at the end of the line on the last word or number without hitting a space character
        WordOrNumber = TRIM(ADJUSTL(WordOrNumber)) // line(i:i)
        LineToRead%LineWordsAndNumbers(j)%WordOrNumber = " "
        LineToRead%LineWordsAndNumbers(j)%WordOrNumber = WordOrNumber
      ELSE ! Add the character to the word or number  
        WordOrNumber = TRIM(ADJUSTL(WordOrNumber)) // line(i:i)
      END IF
    END DO
    
    ! Extract information concerning the items
    
    DO i = 1, LineToRead%NumberOfItems
      
      LineToRead%LineWordsAndNumbers(i)%StrLength = LEN(TRIM(ADJUSTL(LineToRead%LineWordsAndNumbers(i)%WordOrNumber)))

      LineToRead%LineWordsAndNumbers(i)%IsWord = .TRUE.
      LineToRead%LineWordsAndNumbers(i)%IsNumber = .FALSE.
      LineToRead%LineWordsAndNumbers(i)%IsInteger = .FALSE.     
      LineToRead%LineWordsAndNumbers(i)%IsComment = .FALSE.
      LineToRead%LineWordsAndNumbers(i)%IsHyphRange = .FALSE.
      LineToRead%LineWordsAndNumbers(i)%IntValue = 0
      LineToRead%LineWordsAndNumbers(i)%FloatValue = 0D0
     
      dumchar = LineToRead%LineWordsAndNumbers(i)%WordOrNumber(1:1) 
      Comment = IsCommentChar( dumchar )
       
      IF ( Comment ) THEN
      
        LineToRead%LineWordsAndNumbers(i)%IsComment = .TRUE. 
        
      ELSE IF ( IsHyphenatedRange(LineToRead%LineWordsAndNumbers(i)%WordOrNumber) ) THEN
      
        LineToRead%LineWordsAndNumbers(i)%IsWord = .FALSE.
        LineToRead%LineWordsAndNumbers(i)%IsHyphRange = .TRUE.      
      
      ELSE IF ( LowerCase(LineToRead%LineWordsAndNumbers(i)%WordOrNumber) == "e-" ) THEN
      
        CONTINUE ! e- is recognized as a float otherwise   
        
      ELSE
      
        READ(LineToRead%LineWordsAndNumbers(i)%WordOrNumber,*,IOSTAT=ier) dumfloat
        IF (ier == 0) THEN
          LineToRead%LineWordsAndNumbers(i)%IsWord = .FALSE.
          LineToRead%LineWordsAndNumbers(i)%IsNumber = .TRUE.
          LineToRead%LineWordsAndNumbers(i)%IsInteger = .FALSE.     
          LineToRead%LineWordsAndNumbers(i)%IntValue = 0
          LineToRead%LineWordsAndNumbers(i)%FloatValue = dumfloat
        END IF
        
        READ(LineToRead%LineWordsAndNumbers(i)%WordOrNumber,'(I)',IOSTAT=ier) dumint
        IF (ier == 0) THEN
          LineToRead%LineWordsAndNumbers(i)%IsWord = .FALSE.
          LineToRead%LineWordsAndNumbers(i)%IsNumber = .TRUE.
          LineToRead%LineWordsAndNumbers(i)%IsInteger = .TRUE.     
          LineToRead%LineWordsAndNumbers(i)%IntValue = dumint
          LineToRead%LineWordsAndNumbers(i)%FloatValue = REAL(dumint)
        END IF
        CONTINUE
      END IF  
      
    END DO
    
  END IF
  
  RETURN

END SUBROUTINE ReadALine


!***************************************************************************************************************************************!
!********************************************             ******************************************************************************!
!******************************************** Copy a Line ******************************************************************************!
!********************************************             ******************************************************************************!
!***************************************************************************************************************************************!

SUBROUTINE CopyALine(LineToRead,LineToWrite)

  USE TYPES_AND_PARAMETERS
  
  IMPLICIT NONE
  

  TYPE(LINE_IN_FILE)             :: LineToRead
  TYPE(LINE_IN_FILE)             :: LineToWrite
                                 
  INTEGER(I4B)                   :: i
  INTEGER(I4B)                   :: NbItems
  
  IF ( ALLOCATED(LineToWrite%LineWordsAndNumbers)) THEN
    DEALLOCATE(LineToWrite%LineWordsAndNumbers)
  END IF
  
  LineToWrite%FullLine             = LineToRead%FullLine           
  LineToWrite%CleanLine            = LineToRead%CleanLine          
  LineToWrite%PosLineInFile        = LineToRead%PosLineInFile      
  LineToWrite%NameOfBlock          = LineToRead%NameOfBlock        
  LineToWrite%PosLineInBlock       = LineToRead%PosLineInBlock     
  LineToWrite%NumberOfItems        = LineToRead%NumberOfItems      
  LineToWrite%Empty                = LineToRead%Empty              
  LineToWrite%Comment              = LineToRead%Comment            
  
  NbItems = SIZE(LineToRead%LineWordsAndNumbers)  
  ALLOCATE(LineToWrite%LineWordsAndNumbers(NbItems))

  DO i = 1, NbItems
    LineToWrite%LineWordsAndNumbers(i)%StrLength    = LineToRead%LineWordsAndNumbers(i)%StrLength   
    LineToWrite%LineWordsAndNumbers(i)%WordOrNumber = LineToRead%LineWordsAndNumbers(i)%WordOrNumber
    LineToWrite%LineWordsAndNumbers(i)%IsWord       = LineToRead%LineWordsAndNumbers(i)%IsWord      
    LineToWrite%LineWordsAndNumbers(i)%IsNumber     = LineToRead%LineWordsAndNumbers(i)%IsNumber    
    LineToWrite%LineWordsAndNumbers(i)%IsInteger    = LineToRead%LineWordsAndNumbers(i)%IsInteger   
    LineToWrite%LineWordsAndNumbers(i)%IsComment    = LineToRead%LineWordsAndNumbers(i)%IsComment   
    LineToWrite%LineWordsAndNumbers(i)%IsHyphRange  = LineToRead%LineWordsAndNumbers(i)%IsHyphRange 
    LineToWrite%LineWordsAndNumbers(i)%IntValue     = LineToRead%LineWordsAndNumbers(i)%IntValue    
    LineToWrite%LineWordsAndNumbers(i)%FloatValue   = LineToRead%LineWordsAndNumbers(i)%FloatValue  
  END DO
  
  RETURN

END SUBROUTINE CopyALine













!***************************************************************************************************************************************!
!********************************************             ******************************************************************************!
!******************************************** Read a File ******************************************************************************!
!********************************************             ******************************************************************************!
!***************************************************************************************************************************************!

SUBROUTINE ReadAFile(FileName,FileContent)

  USE TYPES_AND_PARAMETERS
  #include "petsc/finclude/petsc.h"
  USE petscmat
  
  IMPLICIT NONE
  
  CHARACTER(LEN=max_len_str)     :: FileName
  TYPE(FILE_CONTENT)             :: FileContent 
  
  INTEGER(I4B)                   :: numfile
  INTEGER(I4B)                   :: NumberOfLines
  LOGICAL(LGT)                   :: FileExists
  INTEGER(I4B)                   :: i
  INTEGER(I4B)                   :: LineJump
  INTEGER(I4B)                   :: TotalLineJump
  TYPE(LINE_IN_FILE)             :: LineToRead                                 


  INTEGER(I4B)                   :: ier
  CHARACTER(LEN=max_len_str)      :: WordOrNumber
  INTEGER(I4B)                   :: j 
  INTEGER(I4B)                   :: ierr
  INTEGER(I4B)                   :: errcode  
  

  numfile = 42
  
  INQUIRE(FILE=TRIM(ADJUSTL(FileName)), EXIST=FileExists)

  IF ( FileExists ) THEN
     OPEN(numfile,FILE=TRIM(ADJUSTL(FileName)),STATUS='UNKNOWN',FORM='FORMATTED',ACTION="READ") 
  ELSE 
    ! Error
    WRITE(*,'(3a)') " Error: File ", TRIM(ADJUSTL(FileName)), " does not exist. STOP"
    STOP
    
  END IF
 
  ! Get the number of lines in the file
  NumberOfLines = 0
  DO
    NumberOfLines = NumberOfLines + 1
    READ(numfile,'(a)', IOSTAT=ier) WordOrNumber
    IF ( ier < 0 ) EXIT
  END DO

  REWIND(numfile) 
    
  IF ( NumberOfLines == 1 .AND. WordOrNumber == "") NumberOfLines = 0 ! File is empty 
  
  IF ( NumberOfLines /= 0 ) THEN
  
    FileContent%NbLinesInFile = NumberOfLines
    FileContent%NbLinesInFileContent = NumberOfLines
    
    ALLOCATE(FileContent%FileLine(NumberOfLines))
    
    i = 1
    TotalLineJump = 0
    DO WHILE ( i <= NumberOfLines )
      LineJump = 0
      
      LineToRead%FullLine = "" 
      LineToRead%CleanLine = "" 
      LineToRead%NameOfBlock = ""
      LineToRead%PosLineInBlock = 0
      LineToRead%NumberOfItems = 0
      LineToRead%Empty = .FALSE.
      LineToRead%Comment = .FALSE.
            
      CALL ReadALine(FileName,numfile,LineToRead,LineJump,i)
      FileContent%FileLine(i-TotalLineJump)%FullLine = LineToRead%FullLine
      FileContent%FileLine(i-TotalLineJump)%CleanLine = LineToRead%CleanLine
      FileContent%FileLine(i-TotalLineJump)%NameOfBlock = ""
      FileContent%FileLine(i-TotalLineJump)%PosLineInBlock = 0
      FileContent%FileLine(i-TotalLineJump)%NumberOfItems = LineToRead%NumberOfItems
      FileContent%FileLine(i-TotalLineJump)%Empty = LineToRead%Empty
      FileContent%FileLine(i-TotalLineJump)%Comment = LineToRead%Comment
      FileContent%FileLine(i-TotalLineJump)%PosLineInFile = i
      
      !WRITE(*,'(a)') LineToRead%CleanLine
      
      ALLOCATE( FileContent%FileLine(i-TotalLineJump)%LineWordsAndNumbers(LineToRead%NumberOfItems) )
      
      DO j = 1, LineToRead%NumberOfItems

        FileContent%FileLine(i-TotalLineJump)%LineWordsAndNumbers(j)%StrLength    = LineToRead%LineWordsAndNumbers(j)%StrLength   
        FileContent%FileLine(i-TotalLineJump)%LineWordsAndNumbers(j)%WordOrNumber = LineToRead%LineWordsAndNumbers(j)%WordOrNumber
        FileContent%FileLine(i-TotalLineJump)%LineWordsAndNumbers(j)%IsWord       = LineToRead%LineWordsAndNumbers(j)%IsWord      
        FileContent%FileLine(i-TotalLineJump)%LineWordsAndNumbers(j)%IsNumber     = LineToRead%LineWordsAndNumbers(j)%IsNumber    
        FileContent%FileLine(i-TotalLineJump)%LineWordsAndNumbers(j)%IsInteger    = LineToRead%LineWordsAndNumbers(j)%IsInteger   
        FileContent%FileLine(i-TotalLineJump)%LineWordsAndNumbers(j)%IsComment    = LineToRead%LineWordsAndNumbers(j)%IsComment   
        FileContent%FileLine(i-TotalLineJump)%LineWordsAndNumbers(j)%IsHyphRange  = LineToRead%LineWordsAndNumbers(j)%IsHyphRange 
        FileContent%FileLine(i-TotalLineJump)%LineWordsAndNumbers(j)%IntValue     = LineToRead%LineWordsAndNumbers(j)%IntValue    
        FileContent%FileLine(i-TotalLineJump)%LineWordsAndNumbers(j)%FloatValue   = LineToRead%LineWordsAndNumbers(j)%FloatValue  
      
      END DO 
      
      FileContent%NbLinesInFileContent = FileContent%NbLinesInFileContent - LineJump
      i = i + 1 + LineJump
      TotalLineJump = TotalLineJump + LineJump
      
    END DO
  
  ELSE 
    ! Error
    WRITE(*,'(3a)') " Error: File ", TRIM(ADJUSTL(FileName)), " is empty. STOP"
    READ(*,*)
    CALL MPI_Abort(PETSC_COMM_WORLD, errcode, ierr)
    
  END IF

  IF ( ALLOCATED(LineToRead%LineWordsAndNumbers) ) THEN
    DEALLOCATE(LineToRead%LineWordsAndNumbers)   
  END IF  

  CLOSE(numfile)
  
  RETURN
  
END SUBROUTINE ReadAFile


!***************************************************************************************************************************************!
!********************************************                    ***********************************************************************!
!******************************************** Write a Clean File ***********************************************************************!
!********************************************                    ***********************************************************************!
!***************************************************************************************************************************************!

SUBROUTINE WriteACleanFile(FileName,FileContent)

  USE TYPES_AND_PARAMETERS
  
  IMPLICIT NONE
  
  CHARACTER(LEN=max_len_str)     :: FileName
  TYPE(FILE_CONTENT)             :: FileContent 
  
  INTEGER(I4B)                   :: numfile
  INTEGER(I4B)                   :: i
  INTEGER(I4B)                   :: j


  numfile = 42  
  
  OPEN(numfile,FILE=TRIM(ADJUSTL("Clean_" // FileName)),STATUS='UNKNOWN',FORM='FORMATTED',ACTION="WRITE") 
  
  DO i = 1, FileContent%NbLinesInFileContent
  
    IF ( (.NOT. FileContent%FileLine(i)%Empty) .AND. (.NOT. FileContent%FileLine(i)%Comment) ) THEN
    
      DO j = 1, FileContent%FileLine(i)%NumberOfItems
        
        IF ( .NOT. FileContent%FileLine(i)%LineWordsAndNumbers(j)%IsComment ) THEN
          
          WRITE(numfile,'(a)',ADVANCE="NO") TRIM(ADJUSTL(FileContent%FileLine(i)%LineWordsAndNumbers(j)%WordOrNumber))
          WRITE(numfile,'(a)',ADVANCE="NO") " "
          
        END IF
      END DO  
      WRITE(numfile, *) ! go tot next line
    END IF    
             
  END DO  
  
  CLOSE(numfile)
  
END SUBROUTINE WriteACleanFile
  
  

!*****************************************************************************************************************************************!
!********************************************                                  ***********************************************************!
!******************************************** Extract blocks from Input files  ***********************************************************!
!********************************************                                  ***********************************************************!
!*****************************************************************************************************************************************!

SUBROUTINE ExtractBlocksFromInputFile(scalar,FileContent,RunTimeBlock ,DiscretizationBlock,OutputBlock,TemperatureBlock,InitialConditionsBlock,RestartConditionsBlock,BoundaryConditionsBlock,FlowBlock,PorosityBlock,PrimarySpeciesBlock,SecondarySpeciesBlock,IonExchangeBlock,MineralsBlock,GasesBlock,SurfaceComplexationBlock,SurfacesAvailableForSurfaceComplexationBlock,ElectrodesBlock,KdBlock,TransportBlock,ConditionsBlock)

  USE TYPES_AND_PARAMETERS
  #include "petsc/finclude/petsc.h"
  USE petscmat



  IMPLICIT NONE

  TYPE(SCALAR_V)                               :: scalar 
  TYPE(FILE_CONTENT)                           :: FileContent 
  TYPE(BLOCK_CONTENT)                          :: RunTimeBlock 
  TYPE(BLOCK_CONTENT)                          :: DiscretizationBlock
  TYPE(BLOCK_CONTENT)                          :: OutputBlock
  TYPE(BLOCK_CONTENT)                          :: TemperatureBlock
  TYPE(BLOCK_CONTENT)                          :: InitialConditionsBlock
  TYPE(BLOCK_CONTENT)                          :: RestartConditionsBlock
  TYPE(BLOCK_CONTENT)                          :: BoundaryConditionsBlock
  TYPE(BLOCK_CONTENT)                          :: FlowBlock
  TYPE(BLOCK_CONTENT)                          :: PorosityBlock
  TYPE(BLOCK_CONTENT)                          :: PrimarySpeciesBlock
  TYPE(BLOCK_CONTENT)                          :: SecondarySpeciesBlock
  TYPE(BLOCK_CONTENT)                          :: IonExchangeBlock
  TYPE(BLOCK_CONTENT)                          :: MineralsBlock
  TYPE(BLOCK_CONTENT)                          :: GasesBlock
  TYPE(BLOCK_CONTENT)                          :: SurfaceComplexationBlock
  TYPE(BLOCK_CONTENT)                          :: SurfacesAvailableForSurfaceComplexationBlock
  TYPE(BLOCK_CONTENT)                          :: ElectrodesBlock
  TYPE(BLOCK_CONTENT)                          :: KdBlock
  TYPE(BLOCK_CONTENT)                          :: TransportBlock
  TYPE(BLOCK_CONTENT),DIMENSION(:),ALLOCATABLE :: ConditionsBlock

  INTEGER(I4B)                                 :: i
  INTEGER(I4B)                                 :: j
  INTEGER(I4B)                                 :: k
  INTEGER(I4B)                                 :: m
  CHARACTER(:),ALLOCATABLE                     :: dumstr
  CHARACTER(:),ALLOCATABLE                     :: dumstr2
  INTEGER(I4B)                                 :: nbconditions

  LOGICAL(LGT)                                 :: ier
  INTEGER(I4B)                                 :: ierr
  INTEGER(I4B)                                 :: errcode  
  
  RunTimeBlock%Exists                                 = .FALSE. 
  DiscretizationBlock%Exists                          = .FALSE.
  OutputBlock%Exists                                  = .FALSE.
  TemperatureBlock%Exists                             = .FALSE.
  InitialConditionsBlock%Exists                       = .FALSE.
  RestartConditionsBlock%Exists                       = .FALSE.
  BoundaryConditionsBlock%Exists                      = .FALSE.
  FlowBlock%Exists                                    = .FALSE.
  PorosityBlock%Exists                                = .FALSE.
  PrimarySpeciesBlock%Exists                          = .FALSE.
  SecondarySpeciesBlock%Exists                        = .FALSE.
  IonExchangeBlock%Exists                             = .FALSE.
  MineralsBlock%Exists                                = .FALSE.
  GasesBlock%Exists                                   = .FALSE.
  SurfaceComplexationBlock%Exists                     = .FALSE.
  SurfacesAvailableForSurfaceComplexationBlock%Exists = .FALSE.
  ElectrodesBlock%Exists                              = .FALSE.
  KdBlock%Exists                                      = .FALSE.
  TransportBlock%Exists                               = .FALSE.

  RunTimeBlock%NbLinesInBlock                                 = 0 
  DiscretizationBlock%NbLinesInBlock                          = 0
  OutputBlock%NbLinesInBlock                                  = 0
  TemperatureBlock%NbLinesInBlock                             = 0
  InitialConditionsBlock%NbLinesInBlock                       = 0
  RestartConditionsBlock%NbLinesInBlock                       = 0
  BoundaryConditionsBlock%NbLinesInBlock                      = 0
  FlowBlock%NbLinesInBlock                                    = 0
  PorosityBlock%NbLinesInBlock                                = 0
  PrimarySpeciesBlock%NbLinesInBlock                          = 0
  SecondarySpeciesBlock%NbLinesInBlock                        = 0
  IonExchangeBlock%NbLinesInBlock                             = 0
  MineralsBlock%NbLinesInBlock                                = 0
  GasesBlock%NbLinesInBlock                                   = 0
  SurfaceComplexationBlock%NbLinesInBlock                     = 0
  SurfacesAvailableForSurfaceComplexationBlock%NbLinesInBlock = 0
  ElectrodesBlock%NbLinesInBlock                              = 0
  KdBlock%NbLinesInBlock                                      = 0
  TransportBlock%NbLinesInBlock                               = 0

  
  RunTimeBlock%BlockDescription                                 = "" 
  DiscretizationBlock%BlockDescription                          = ""
  OutputBlock%BlockDescription                                  = ""
  TemperatureBlock%BlockDescription                             = ""
  InitialConditionsBlock%BlockDescription                       = ""
  RestartConditionsBlock%BlockDescription                       = ""
  BoundaryConditionsBlock%BlockDescription                      = ""
  FlowBlock%BlockDescription                                    = ""
  PorosityBlock%BlockDescription                                = ""
  PrimarySpeciesBlock%BlockDescription                          = ""
  SecondarySpeciesBlock%BlockDescription                        = ""
  IonExchangeBlock%BlockDescription                             = ""
  MineralsBlock%BlockDescription                                = ""
  GasesBlock%BlockDescription                                   = ""
  SurfaceComplexationBlock%BlockDescription                     = ""
  SurfacesAvailableForSurfaceComplexationBlock%BlockDescription = ""
  ElectrodesBlock%BlockDescription                              = ""
  KdBlock%BlockDescription                                      = ""
  TransportBlock%BlockDescription                               = ""
  
  
  ! Get the number of conditions
  i = 1
  nbconditions = 0
  DO WHILE (i <= FileContent%NbLinesInFileContent)  
    IF ( (.NOT. FileContent%FileLine(i)%Empty) .AND. (.NOT. FileContent%FileLine(i)%Comment) ) THEN
      dumstr2 = LowerCase(FileContent%FileLine(i)%LineWordsAndNumbers(1)%WordOrNumber)
      IF (dumstr2 == "condition") nbconditions = nbconditions + 1
    END IF
    i = i + 1
  END DO 
  
  ALLOCATE( ConditionsBlock(nbconditions) )
  ConditionsBlock(:)%Exists                              = .FALSE.
  ConditionsBlock(:)%NbLinesInBlock                      = 0
  
  
  i = 1
  m = 0
  DO WHILE (i <= FileContent%NbLinesInFileContent)
  
    IF ( (.NOT. FileContent%FileLine(i)%Empty) .AND. (.NOT. FileContent%FileLine(i)%Comment) ) THEN
    
      dumstr2 = LowerCase(FileContent%FileLine(i)%LineWordsAndNumbers(1)%WordOrNumber)
      SELECT CASE( dumstr2 )
      
      CASE("title")
      
        IF ( FileContent%FileLine(i)%NumberOfItems > 1 ) THEN
          DO j = 2, FileContent%FileLine(i)%NumberOfItems
            scalar%ltitle = TRIM(ADJUSTL(scalar%ltitle)) // " " // FileContent%FileLine(i)%LineWordsAndNumbers(j)%WordOrNumber
          END DO       
        END IF          
     
      CASE("runtime")

        dumstr = "RUNTIME"
        RunTimeBlock%Exists = .TRUE.
        RunTimeBlock%PosBlockInFile = FileContent%FileLine(i)%PosLineInFile
        j = 0
        CALL FindTheEnd(i,j,FileContent,dumstr)
        
        ! j is the number of lines of information (or comment/empty lines) in the block
        RunTimeBlock%NbLinesInBlock = j
        ALLOCATE(RunTimeBlock%BlockLine(j))       
        DO k = 1, j        
          CALL CopyALine( FileContent%FileLine(RunTimeBlock%PosBlockInFile + k), RunTimeBlock%BlockLine(k) )
        END DO
      
      CASE("discretization")

        dumstr = "DISCRETIZATION"
        DiscretizationBlock%Exists = .TRUE.
        DiscretizationBlock%PosBlockInFile = FileContent%FileLine(i)%PosLineInFile
        j = 0
        CALL FindTheEnd(i,j,FileContent,dumstr)
        
        ! j is the number of lines of information (or comment/empty lines) in the block
        DiscretizationBlock%NbLinesInBlock = j
        ALLOCATE(DiscretizationBlock%BlockLine(j))       
        DO k = 1, j        
          CALL CopyALine( FileContent%FileLine(DiscretizationBlock%PosBlockInFile + k), DiscretizationBlock%BlockLine(k) )
        END DO
            
      CASE("output")

        dumstr = "OUTPUT"
        OutputBlock%Exists = .TRUE.
        OutputBlock%PosBlockInFile = FileContent%FileLine(i)%PosLineInFile
        j = 0
        CALL FindTheEnd(i,j,FileContent,dumstr)
        
        ! j is the number of lines of information (or comment/empty lines) in the block
        OutputBlock%NbLinesInBlock = j
        ALLOCATE(OutputBlock%BlockLine(j))       
        DO k = 1, j     
          CALL CopyALine( FileContent%FileLine(OutputBlock%PosBlockInFile + k), OutputBlock%BlockLine(k) )
        END DO
        
      CASE("temperature")
      
        dumstr = "TEMPERATURE"
        TemperatureBlock%Exists = .TRUE.
        TemperatureBlock%PosBlockInFile = FileContent%FileLine(i)%PosLineInFile
        j = 0
        CALL FindTheEnd(i,j,FileContent,dumstr)
        
        ! j is the number of lines of information (or comment/empty lines) in the block
        TemperatureBlock%NbLinesInBlock = j
        ALLOCATE(TemperatureBlock%BlockLine(j))       
        DO k = 1, j        
          CALL CopyALine( FileContent%FileLine(TemperatureBlock%PosBlockInFile + k), TemperatureBlock%BlockLine(k) )
        END DO
      
      CASE("initial_conditions")
      
        dumstr = "INITIAL_CONDITIONS"
        InitialConditionsBlock%Exists = .TRUE.
        InitialConditionsBlock%PosBlockInFile = FileContent%FileLine(i)%PosLineInFile
        j = 0
        CALL FindTheEnd(i,j,FileContent,dumstr)
        
        ! j is the number of lines of information (or comment/empty lines) in the block
        InitialConditionsBlock%NbLinesInBlock = j
        ALLOCATE(InitialConditionsBlock%BlockLine(j))       
        DO k = 1, j        
          CALL CopyALine( FileContent%FileLine(InitialConditionsBlock%PosBlockInFile + k), InitialConditionsBlock%BlockLine(k) )
        END DO

      CASE("restart_conditions")
      
        dumstr = "RESTART_CONDITIONS"
        RestartConditionsBlock%Exists = .TRUE.
        RestartConditionsBlock%PosBlockInFile = FileContent%FileLine(i)%PosLineInFile
        j = 0
        CALL FindTheEnd(i,j,FileContent,dumstr)
        
        ! j is the number of lines of information (or comment/empty lines) in the block
        RestartConditionsBlock%NbLinesInBlock = j
        ALLOCATE(RestartConditionsBlock%BlockLine(j))       
        DO k = 1, j        
          CALL CopyALine( FileContent%FileLine(RestartConditionsBlock%PosBlockInFile + k), RestartConditionsBlock%BlockLine(k) )
        END DO
        
      CASE("boundary_conditions")
      
        dumstr = "BOUNDARY_CONDITIONS"
        BoundaryConditionsBlock%Exists = .TRUE.
        BoundaryConditionsBlock%PosBlockInFile = FileContent%FileLine(i)%PosLineInFile
        j = 0
        CALL FindTheEnd(i,j,FileContent,dumstr)
        
        ! j is the number of lines of information (or comment/empty lines) in the block
        BoundaryConditionsBlock%NbLinesInBlock = j
        ALLOCATE(BoundaryConditionsBlock%BlockLine(j))       
        DO k = 1, j        
          CALL CopyALine( FileContent%FileLine(BoundaryConditionsBlock%PosBlockInFile + k), BoundaryConditionsBlock%BlockLine(k) )
        END DO

      CASE("flow")
      
        dumstr = "FLOW"
        FlowBlock%Exists = .TRUE.
        FlowBlock%PosBlockInFile = FileContent%FileLine(i)%PosLineInFile
        j = 0
        CALL FindTheEnd(i,j,FileContent,dumstr)
        
        ! j is the number of lines of information (or comment/empty lines) in the block
        FlowBlock%NbLinesInBlock = j
        ALLOCATE(FlowBlock%BlockLine(j))       
        DO k = 1, j        
          CALL CopyALine( FileContent%FileLine(FlowBlock%PosBlockInFile + k), FlowBlock%BlockLine(k) )
        END DO

      CASE("porosity")
      
        dumstr = "POROSITY"
        PorosityBlock%Exists = .TRUE.
        PorosityBlock%PosBlockInFile = FileContent%FileLine(i)%PosLineInFile
        j = 0
        CALL FindTheEnd(i,j,FileContent,dumstr)
        
        ! j is the number of lines of information (or comment/empty lines) in the block
        PorosityBlock%NbLinesInBlock = j
        ALLOCATE(PorosityBlock%BlockLine(j))       
        DO k = 1, j        
          CALL CopyALine( FileContent%FileLine(PorosityBlock%PosBlockInFile + k), PorosityBlock%BlockLine(k) )
        END DO

      CASE("primary_species")
      
        dumstr = "PRIMARY_SPECIES"
        PrimarySpeciesBlock%Exists = .TRUE.
        PrimarySpeciesBlock%PosBlockInFile = FileContent%FileLine(i)%PosLineInFile
        j = 0
        CALL FindTheEnd(i,j,FileContent,dumstr)
        
        ! j is the number of lines of information (or comment/empty lines) in the block
        PrimarySpeciesBlock%NbLinesInBlock = j
        ALLOCATE(PrimarySpeciesBlock%BlockLine(j))       
        DO k = 1, j       
          CALL CopyALine( FileContent%FileLine(PrimarySpeciesBlock%PosBlockInFile + k), PrimarySpeciesBlock%BlockLine(k) ) 
        END DO
      
      CASE("secondary_species")
      
        dumstr = "SECONDARY_SPECIES"
        SecondarySpeciesBlock%Exists = .TRUE.
        SecondarySpeciesBlock%PosBlockInFile = FileContent%FileLine(i)%PosLineInFile
        j = 0
        CALL FindTheEnd(i,j,FileContent,dumstr)
        
        ! j is the number of lines of information (or comment/empty lines) in the block
        SecondarySpeciesBlock%NbLinesInBlock = j
        ALLOCATE(SecondarySpeciesBlock%BlockLine(j))       
        DO k = 1, j        
          CALL CopyALine( FileContent%FileLine(SecondarySpeciesBlock%PosBlockInFile + k), SecondarySpeciesBlock%BlockLine(k) ) 
        END DO
      
      CASE("ion_exchange")
      
        dumstr = "ION_EXCHANGE"
        IonExchangeBlock%Exists = .TRUE.
        IonExchangeBlock%PosBlockInFile = FileContent%FileLine(i)%PosLineInFile
        j = 0
        CALL FindTheEnd(i,j,FileContent,dumstr)
        
        ! j is the number of lines of information (or comment/empty lines) in the block
        IonExchangeBlock%NbLinesInBlock = j
        ALLOCATE(IonExchangeBlock%BlockLine(j))       
        DO k = 1, j        
          CALL CopyALine( FileContent%FileLine(IonExchangeBlock%PosBlockInFile + k), IonExchangeBlock%BlockLine(k) ) 
        END DO
      
      CASE("minerals")
      
        dumstr = "MINERALS"
        MineralsBlock%Exists = .TRUE.
        MineralsBlock%PosBlockInFile = FileContent%FileLine(i)%PosLineInFile
        j = 0
        CALL FindTheEnd(i,j,FileContent,dumstr)
        
        ! j is the number of lines of information (or comment/empty lines) in the block
        MineralsBlock%NbLinesInBlock = j
        ALLOCATE(MineralsBlock%BlockLine(j))       
        DO k = 1, j        
          CALL CopyALine( FileContent%FileLine(MineralsBlock%PosBlockInFile + k), MineralsBlock%BlockLine(k) ) 
        END DO
      
      CASE("gases")
      
        dumstr = "GASES"
        GasesBlock%Exists = .TRUE.
        GasesBlock%PosBlockInFile = FileContent%FileLine(i)%PosLineInFile
        j = 0
        CALL FindTheEnd(i,j,FileContent,dumstr)
        
        ! j is the number of lines of information (or comment/empty lines) in the block
        GasesBlock%NbLinesInBlock = j
        ALLOCATE(GasesBlock%BlockLine(j))       
        DO k = 1, j        
          CALL CopyALine( FileContent%FileLine(GasesBlock%PosBlockInFile + k), GasesBlock%BlockLine(k) ) 
        END DO
      
      CASE("surface_complexation")
      
        dumstr = "SURFACE_COMPLEXATION"
        SurfaceComplexationBlock%Exists = .TRUE.
        SurfaceComplexationBlock%PosBlockInFile = FileContent%FileLine(i)%PosLineInFile
        j = 0
        CALL FindTheEnd(i,j,FileContent,dumstr)
        
        ! j is the number of lines of information (or comment/empty lines) in the block
        SurfaceComplexationBlock%NbLinesInBlock = j
        ALLOCATE(SurfaceComplexationBlock%BlockLine(j))       
        DO k = 1, j        
          CALL CopyALine( FileContent%FileLine(SurfaceComplexationBlock%PosBlockInFile + k), SurfaceComplexationBlock%BlockLine(k) ) 
        END DO
      
      CASE("surfaces_available_for_surface_complexation")
      
        dumstr = "SURFACES_AVAILABLE_FOR_SURFACE_COMPLEXATION"
        SurfacesAvailableForSurfaceComplexationBlock%Exists = .TRUE.
        SurfacesAvailableForSurfaceComplexationBlock%PosBlockInFile = FileContent%FileLine(i)%PosLineInFile
        j = 0
        CALL FindTheEnd(i,j,FileContent,dumstr)
        
        ! j is the number of lines of information (or comment/empty lines) in the block
        SurfacesAvailableForSurfaceComplexationBlock%NbLinesInBlock = j
        ALLOCATE(SurfacesAvailableForSurfaceComplexationBlock%BlockLine(j))       
        DO k = 1, j 
          CALL CopyALine( FileContent%FileLine(SurfacesAvailableForSurfaceComplexationBlock%PosBlockInFile + k), SurfacesAvailableForSurfaceComplexationBlock%BlockLine(k) ) 
        END DO

      CASE("electrodes")
      
        dumstr = "ELECTRODES"
        ElectrodesBlock%Exists = .TRUE.
        ElectrodesBlock%PosBlockInFile = FileContent%FileLine(i)%PosLineInFile
        j = 0
        CALL FindTheEnd(i,j,FileContent,dumstr)
        
        ! j is the number of lines of information (or comment/empty lines) in the block
        ElectrodesBlock%NbLinesInBlock = j
        ALLOCATE(ElectrodesBlock%BlockLine(j))       
        DO k = 1, j        
          CALL CopyALine( FileContent%FileLine(ElectrodesBlock%PosBlockInFile + k), ElectrodesBlock%BlockLine(k) ) 
        END DO
      
      CASE("kd")
      
        dumstr = "KD"
        KdBlock%Exists = .TRUE.
        KdBlock%PosBlockInFile = FileContent%FileLine(i)%PosLineInFile
        j = 0
        CALL FindTheEnd(i,j,FileContent,dumstr)
        
        ! j is the number of lines of information (or comment/empty lines) in the block
        KdBlock%NbLinesInBlock = j
        ALLOCATE(KdBlock%BlockLine(j))       
        DO k = 1, j        
          CALL CopyALine( FileContent%FileLine(KdBlock%PosBlockInFile + k), KdBlock%BlockLine(k) ) 
        END DO
      
      CASE("transport")
          
        dumstr = "TRANSPORT"
        TransportBlock%Exists = .TRUE.
        TransportBlock%PosBlockInFile = FileContent%FileLine(i)%PosLineInFile
        j = 0
        CALL FindTheEnd(i,j,FileContent,dumstr)
        
        ! j is the number of lines of information (or comment/empty lines) in the block
        TransportBlock%NbLinesInBlock = j
        ALLOCATE(TransportBlock%BlockLine(j))       
        DO k = 1, j        
          CALL CopyALine( FileContent%FileLine(TransportBlock%PosBlockInFile + k), TransportBlock%BlockLine(k) ) 
        END DO
        
      CASE("condition")
        m = m + 1 
      
        dumstr = "CONDITION"
        ConditionsBlock(m)%Exists = .TRUE.
        ConditionsBlock(m)%PosBlockInFile = FileContent%FileLine(i)%PosLineInFile
        ConditionsBlock(m)%BlockDescription = ""
        ier = .FALSE.
        
        IF ( FileContent%FileLine(i)%NumberOfItems > 1 ) THEN        
          ConditionsBlock(m)%BlockDescription = FileContent%FileLine(i)%LineWordsAndNumbers(2)%WordOrNumber          
        ELSE
          ier = .TRUE.
        END IF  

        IF (ConditionsBlock(m)%BlockDescription == "") ier = .TRUE.
        
        IF (ier == .TRUE.) THEN
        
          WRITE(*,'(a,I,a)') " ERROR: CONDITION block on line ", i, " in input file must have a name."
          WRITE(*,'(a)') " STOP "
          CALL MPI_Abort(PETSC_COMM_WORLD, errcode, ierr)
          READ(*,*)
          
        END IF
                  
        j = 0
        CALL FindTheEnd(i,j,FileContent,dumstr)
        
        ! j is the number of lines of information (or comment/empty lines) in the block
        ConditionsBlock(m)%NbLinesInBlock = j
        ALLOCATE(ConditionsBlock(m)%BlockLine(j))       
        DO k = 1, j        
          CALL CopyALine( FileContent%FileLine(ConditionsBlock(m)%PosBlockInFile + k), ConditionsBlock(m)%BlockLine(k) ) 
        END DO
      

      CASE DEFAULT
        
        CONTINUE
        
      END SELECT
       
    END IF
  
    i = i + 1 
  
  END DO
  
  ! Check that all conditions have a different name
  
  IF ( nbconditions > 1 ) THEN
  
    DO i = 1, nbconditions - 1
      DO j = i + 1, nbconditions
      
        IF ( ConditionsBlock(i)%BlockDescription == ConditionsBlock(j)%BlockDescription ) THEN
        
          WRITE(*,'(3a)') " ERROR: CONDITION ", TRIM(ADJUSTL(ConditionsBlock(i)%BlockDescription)), " is defined at least twice in the input file."
          WRITE(*,'(a)') " STOP "
          CALL MPI_Abort(PETSC_COMM_WORLD, errcode, ierr)
          READ(*,*)

        END IF
      
      END DO
    END DO
    
  END IF  
  
  
  RETURN
  
END SUBROUTINE ExtractBlocksFromInputFile
  











END MODULE Easy_Read_and_Store_M    