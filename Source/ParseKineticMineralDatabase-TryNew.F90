SUBROUTINE read_minkin(nout,ncomp,nspec,nkin,ngas,data1,namrl,  &
                       ispeciate,igenericrates,nammin_tmp)
USE crunchtype
USE CrunchFunctions
USE params
USE concentration
USE mineral
USE strings
USE NanoCrystal

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nout
INTEGER(I4B), INTENT(IN)                                    :: ncomp
INTEGER(I4B), INTENT(IN)                                    :: nspec
INTEGER(I4B), INTENT(OUT)                                   :: nkin
INTEGER(I4B), INTENT(IN)                                    :: ngas

INTEGER(I4B), INTENT(IN)                                    :: ispeciate
INTEGER(I4B), INTENT(IN)                                    :: igenericrates

CHARACTER (LEN=mls), DIMENSION(:), INTENT(IN OUT)           :: nammin_tmp
CHARACTER (LEN=mls), INTENT(IN)                             :: data1
CHARACTER (LEN=mls), DIMENSION(:), INTENT(IN OUT)           :: namrl

CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE                :: workchar1

!  Internal variables and arrays


INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs
INTEGER(I4B)                                                :: n0
INTEGER(I4B)                                                :: mnrl
INTEGER(I4B)                                                :: npar
INTEGER(I4B)                                                :: nreac
INTEGER(I4B)                                                :: i
INTEGER(I4B)                                                :: k
INTEGER(I4B)                                                :: ksave
INTEGER(I4B)                                                :: lsave
INTEGER(I4B)                                                :: ik
INTEGER(I4B)                                                :: ndep
INTEGER(I4B)                                                :: monod
INTEGER(I4B)                                                :: ls_spec
INTEGER(I4B)                                                :: lss
INTEGER(I4B)                                                :: itot_tmp
INTEGER(I4B)                                                :: ls_monod
INTEGER(I4B)                                                :: ksp
INTEGER(I4B)                                                :: ls_inhibit
INTEGER(I4B)                                                :: np
INTEGER(I4B)                                                :: lslabel
INTEGER(I4B)                                                :: lenlabel
INTEGER(I4B)                                                :: kk
INTEGER(I4B)                                                :: llen
INTEGER(I4B)                                                :: ls_mineral

INTEGER(I4B)                                                :: iunitDatabase

REAL(DP)                                                    :: depend_tmp

LOGICAL(LGT)                                                :: SpeciesFound
LOGICAL(LGT)                                                :: MineralFound
LOGICAL(LGT)                                                ::  itsopen 

CHARACTER (LEN=mls)                                         :: dummy1
CHARACTER (LEN=mls)                                         :: dumstring
CHARACTER (LEN=mls)                                         :: tempmin
CHARACTER (LEN=mls)                                         :: namdum
CHARACTER (LEN=mls)                                         :: nam_depend
CHARACTER (LEN=mls)                                         :: nam_inhibit
CHARACTER (LEN=mls)                                         :: nam_monod
CHARACTER (LEN=mls)                                         :: templabel

CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE              :: namspecies
CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE              :: NAME

REAL(DP), DIMENSION(:), ALLOCATABLE                         :: dep_tmp
REAL(DP), DIMENSION(:), ALLOCATABLE                         :: iuser
REAL(DP), DIMENSION(:), ALLOCATABLE                         :: rinhibit_tmp
REAL(DP), DIMENSION(:), ALLOCATABLE                         :: halfsat_tmp

!REAL(DP), DIMENSION(:), ALLOCATABLE                         :: sto
!REAL(DP), DIMENSION(:,:), ALLOCATABLE                       :: ax
!REAL(DP), DIMENSION(:,:), ALLOCATABLE                       :: ainv
!REAL(DP), DIMENSION(:,:), ALLOCATABLE                       :: bx
!REAL(DP), DIMENSION(:,:), ALLOCATABLE                       :: yx
!REAL(DP), DIMENSION(:,:), ALLOCATABLE                       :: indx

character (len=mls)                                         :: biomass
integer(i4b)                                                :: ios, im, ib, count_ib
integer(i4b),dimension(:),allocatable                       :: workint

namelist /BiomassDecay/                                        biomass

CHARACTER (LEN=mls)                                         :: DatabaseString

type :: RowData
  integer, allocatable :: cols(:)
end type

type(RowData), allocatable :: rows(:)

integer :: number_rows
integer :: i

nlines = 0
OPEN (1, file = 'file.txt')
DO
    READ (1,*, END=10)
    nlines = nlines + 1
END DO
10 CLOSE (1)

number_rows = count_lines(11)
allocate(rows(number_rows))

do i=1,number_rows
  rows(i)%cols = read_row(11)
enddo

!!!  ******************************************************************

!   First, find end of the minerals section

iunitDatabase = 8

inquire(unit=iunitDatabase, opened=itsopen) 
IF ( itsopen ) THEN
  ClOSE(UNIT=iunitDatabase,status='KEEP')
END IF
OPEN(UNIT=iunitDatabase,FILE=data1,STATUS='old',err=334)
  
  
200 READ(iunitDatabase,'(a)') DatabaseString
    
 
MineralKineticsFound = .FALSE
IF (DatabaseString == 'Begin mineral kinetics') THEN
  write(*,*) ' Mineral kinetics section found in database'
  MineralKineticsFound = .TRUE.
END IF

475 READ(iunitDatabase,'(a)',END=300) dummy1

IF (dummy1 == 'End of mineral kinetics') THEN
  dumstring = nammin_tmp(nkin)
  namdum = rlabel(np,nkin)
  CALL stringlen(namdum,lslabel)
  WRITE(*,*)
  WRITE(*,*) ' Kinetic mineral reaction not found in database'
  WRITE(*,*) ' Looking for ', dumstring(1:lsave)
  WRITE(*,*) ' Rate label = ',namdum(1:lslabel)
  WRITE(*,*) ' STOP'
  WRITE(*,*)
  READ(*,*)
  STOP
END IF

tempmin = nammin_tmp(nkin)
id = 1
iff = mls
CALL sschaine(dummy1,id,iff,ssch,ids,ls)
lzs=ls
!      call convan(ssch,lzs,res)
CALL stringtype(ssch,lzs,res)
!      write(*,*) ls
!      write(*,*) ssch(1:ls)
!      write(*,*) tempmin(1:ls)
IF (ssch == tempmin) THEN
  READ(iunitDatabase,'(a)',END=300) dummy1
  id = 1
  iff = mls
  CALL sschaine(dummy1,id,iff,ssch,ids,ls)
  lzs=ls
  CALL convan(ssch,lzs,res)
  IF (ssch == 'label') THEN
    id = ids + ls
    CALL sschaine(dummy1,id,iff,ssch,ids,ls)
    IF (ssch /= '=') THEN
      WRITE(*,*)
      WRITE(*,*) ' Error in database reading "label" '
      WRITE(*,*) ' Mineral = ',tempmin(1:lsave)
      WRITE(*,*) ' Label = ',templabel(1:lenlabel)
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
    id = ids + ls
    CALL sschaine(dummy1,id,iff,ssch,ids,ls)
!          write(*,*) ' Mineral and label '
!          write(*,*) tempmin
!          write(*,*) ssch
!          write(*,*)
!          write(*,*)
    IF (ssch /= rlabel(np,nkin)) THEN
      CALL breakfind
      GO TO 475
    END IF
    templabel = rlabel(np,nkin)
    CALL stringlen(templabel,lenlabel)
    
    READ(iunitDatabase,'(a)',END=300) dummy1
    IF (iuser(4) == 0) THEN  ! No user input, so read
      id = 1
      iff = mls
      CALL sschaine(dummy1,id,iff,ssch,ids,ls)
      lzs=ls
!!!      CALL convan(ssch,lzs,res)
      IF (ssch == 'type') THEN                              !! Read TYPE
        id = ids + ls
        CALL sschaine(dummy1,id,iff,ssch,ids,ls)
        IF (ssch == '=') THEN
          id = ids + ls
          CALL sschaine(dummy1,id,iff,ssch,ids,ls)
          IF (ssch == 'tst') THEN
            imintype(npar,nkin) = 1
          ELSE IF (ssch == 'monod' .OR. ssch == 'Monod' .OR. ssch == 'MONOD') THEN
            imintype(npar,nkin) = 2
          ELSE IF (ssch == 'nucleation') THEN
            imintype(npar,nkin) = 10
          ELSE IF (ssch == 'irreversible' .OR. ssch == 'Irreversible') THEN
            imintype(npar,nkin) = 3
          ELSE IF (ssch == 'PrecipitationOnly' .OR. ssch == 'Precipitationonly' .OR. ssch == 'precipitationonly') THEN
            imintype(npar,nkin) = 4
          ELSE IF (ssch == 'DissolutionOnly' .OR. ssch == 'Dissolutiononly' .OR. ssch == 'dissolutiononly') THEN
            imintype(npar,nkin) = 5
          ELSE IF (ssch == 'Forward' .OR. ssch == 'ForwardOnly' .OR. ssch == 'forwardonly') THEN
            imintype(npar,nkin) = 6
          ELSE IF (ssch == 'Reverse' .OR. ssch == 'ReverseOnly' .OR. ssch == 'reverseonly') THEN
            imintype(npar,nkin) = 7

! biomass
          ELSE IF (ssch == 'MonodBiomass' .OR. ssch == 'monodbiomass') THEN
            imintype(npar,nkin) = 8
! biomass end 

! sergi: biomass decay
         ELSE IF (ssch == 'BiomassDecay' .OR. ssch == 'biomassdecay') THEN
           imintype(npar,nkin) = 9


! sergi: biomass decay end
          ELSE
            WRITE(*,*)
            WRITE(*,*) 'Rate type not recognized in database'
            WRITE(*,*) ' Mineral = ',tempmin(1:lsave)
            WRITE(*,*) ' Label = ',templabel(1:lenlabel)
            WRITE(*,*)
            READ(*,*)
            STOP
          END IF
        ELSE
          WRITE(*,*)
          WRITE(*,*) ' Error in database: Rate "type" '
          WRITE(*,*) ' Mineral = ',tempmin(1:lsave)
          WRITE(*,*) ' Label = ',templabel(1:lenlabel)
          WRITE(*,*)
          READ(*,*)
          STOP
        END IF
      ELSE
        WRITE(*,*)
        WRITE(*,*) '  Error in database: Rate "type"'
        WRITE(*,*) ' Mineral = ',tempmin(1:lsave)
        WRITE(*,*) ' Label = ',templabel(1:lenlabel)
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
      
    END IF                                                      ! If input by user, skip block above
    
    READ(iunitDatabase,'(a)',END=300) dummy1
    
    IF (iuser(1) == 0) THEN ! No user input, so read
      id = 1
      iff = mls
      CALL sschaine(dummy1,id,iff,ssch,ids,ls)
      lzs=ls
      CALL convan(ssch,lzs,res)
      IF (ssch == 'rate(85c)') THEN                              !!  Read RATE
        id = ids + ls
        CALL sschaine(dummy1,id,iff,ssch,ids,ls)
        IF (ssch == '=') THEN
          id = ids + ls
          CALL sschaine(dummy1,id,iff,ssch,ids,ls)
          lzs = ls
          CALL stringtype(ssch,lzs,res)
          IF (res == 'a') THEN
            WRITE(*,*)
            WRITE(*,*) ' "Rate" should be followed by a number'
            WRITE(*,*) ' Mineral = ',namdum(1:lsave)
            WRITE(*,*) ' String = ',ssch(1:ls)
            WRITE(*,*)
            READ(*,*)
            STOP
          ELSE
            rate0(npar,nkin) = DNUM(ssch)
            continue
          END IF
        ELSE
          WRITE(*,*)
          WRITE(*,*) ' Error in database: Reaction "rate" '
          WRITE(*,*) ' Mineral = ',tempmin(1:lsave)
          WRITE(*,*) ' Label = ',templabel(1:lenlabel)
          WRITE(*,*)
          READ(*,*)
          STOP
        END IF
      ELSE
        WRITE(*,*)
        WRITE(*,*) '  Error in database: Reaction "rate" '
        WRITE(*,*) ' Mineral = ',tempmin(1:lsave)
        WRITE(*,*) ' Label = ',templabel(1:lenlabel)
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
      
    END IF                                                  ! If input by user, skip block above
    
    READ(iunitDatabase,'(a)',END=300) dummy1
    IF (iuser(8) == 0) THEN ! No user input, so read
      id = 1
      iff = mls
      CALL sschaine(dummy1,id,iff,ssch,ids,ls)
      lzs=ls
      CALL convan(ssch,lzs,res)
      IF (ssch == 'activation') THEN                         !! Read ACTIVATION ENERGY
        id = ids + ls
        CALL sschaine(dummy1,id,iff,ssch,ids,ls)
        IF (ssch == '=') THEN
          id = ids + ls
          CALL sschaine(dummy1,id,iff,ssch,ids,ls)
          lzs = ls
          CALL stringtype(ssch,lzs,res)
          IF (res == 'a') THEN
            WRITE(*,*)
            WRITE(*,*) ' "Activation" should be followed by a number'
            WRITE(*,*) ' Mineral = ',namdum(1:lsave)
            WRITE(*,*) ' String = ',ssch(1:ls)
            WRITE(*,*)
            READ(*,*)
            STOP
          ELSE
            ea(np,nkin) = DNUM(ssch)
          END IF
        ELSE
          WRITE(*,*)
          WRITE(*,*) ' Error in database: Reaction "activation" '
          WRITE(*,*) ' Mineral = ',tempmin(1:lsave)
          WRITE(*,*) ' Label = ',templabel(1:lenlabel)
          WRITE(*,*)
          READ(*,*)
          STOP
        END IF
      ELSE
        WRITE(*,*)
        WRITE(*,*) '  Error in database: Reaction "activation" '
        WRITE(*,*) ' Mineral = ',tempmin(1:lsave)
        WRITE(*,*) ' Label = ',templabel(1:lenlabel)
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
      
    END IF                                                             !! If input by user, skip block above
    
    IF (imintype(np,nkin) == 2) THEN                                   !! Read Monod parameters
      
!!    Look for Monod terms (species specification with half-saturation constants)
      
      READ(iunitDatabase,'(a)',END=300) dummy1
      IF (iuser(5) == 0) THEN ! No user input, so read
        id = 1
        iff = mls
        CALL sschaine(dummy1,id,iff,ssch,ids,ls)
        lzs=ls
        CALL convan(ssch,lzs,res)
        IF (ssch == 'monod_term' .OR. ssch == 'monod_terms') THEN
          nmonod(np,nkin) = 0
          id = ids + ls
          CALL sschaine(dummy1,id,iff,ssch,ids,ls)
          IF (ssch == ':') THEN
            333             id = ids + ls
            CALL sschaine(dummy1,id,iff,ssch,ids,ls)
            IF (ls /= 0) THEN
              nmonod(np,nkin) = nmonod(np,nkin) + 1
              IF (nmonod(np,nkin) > mpre) THEN
                WRITE(*,*)
                WRITE(*,*) ' Redimension parameter "mpre" in params.inc'
                WRITE(*,*) ' Mpre = ',mpre
                WRITE(*,*) ' Number of Monod terms = ', nmonod(np,nkin)
                WRITE(*,*)
                READ(*,*)
                STOP
              END IF
              lzs = ls
!                    call convan(ssch,lzs,res)
              CALL stringtype(ssch,lzs,res)
              NAME(nmonod(np,nkin)) = ssch
              id = ids + ls
              CALL sschaine(dummy1,id,iff,ssch,ids,ls)
              lzs = ls
              CALL stringtype (ssch,lzs,res)
              IF (res == 'a') THEN
                WRITE(*,*)
                WRITE(*,*) '  Error in database: Monod terms'
                WRITE(*,*) ' Mineral = ',tempmin(1:lsave)
                WRITE(*,*) ' Label = ',templabel(1:lenlabel)
                WRITE(*,*)
                READ(*,*)
                STOP
              ELSE
                halfsat_tmp(nmonod(np,nkin)) = DNUM(ssch)
              END IF
              GO TO 333
            END IF
            
            DO kk = 1,nmonod(np,nkin)
              dumstring = NAME(kk)
              IF (dumstring(1:4) == 'tot_') THEN
                CALL stringlen(dumstring,lss)
                NAME(kk) = dumstring(5:lss)
                itot_tmp = 1
              ELSE
                itot_tmp = 0
              END IF
              IF (NAME(kk) == 'self' .OR. NAME(kk) == 'Self' .OR. NAME(kk) == 'SELF') THEN           !!  Rate depends on mineral concentration, mol/L

                speciesfound = .TRUE.
                kmonod(kk,np,k) = 1                                  !!  Flag to indicate that Monod rate depends on the mineral concentration itself
                imonod(kk,np,k) = 0                                  !!  If cycling through species, this should tell the code (in Reaction.F90) to check kmonod(kk,np,k)
                halfsat(kk,np,nkin) = halfsat_tmp(kk)            

              ELSE
!!              Check to see if name in Monod term is recognized as a species
                speciesfound = .false.
                DO i = 1,ncomp                                  !! Loop over primary species
                  IF (ulab(i) == NAME(kk)) THEN
                    imonod(kk,np,nkin) = i
                    speciesfound = .true.
                    halfsat(kk,np,nkin) = halfsat_tmp(kk)
                    IF (itot_tmp == 1) THEN
                      itot_monod(kk,np,nkin) = 1
                    END IF
                  END IF
                END DO
              END IF
              IF (speciesfound) THEN
                CONTINUE
              ELSE
                DO ksp = 1,nspec                              !! Loop over secondary species
                  ik = ncomp + ksp
                  IF (ulab(ik) == NAME(kk)) THEN
                    imonod(kk,np,nkin) = ik
                    speciesfound = .true.
                    halfsat(kk,np,nkin) = halfsat_tmp(kk)
                    IF (itot_tmp == 0) THEN
                      CONTINUE
                    ELSE
                      WRITE(*,*)
                      WRITE(*,*) ' Cannot specify total concentration for secondary species'
                      WRITE(*,*) ' in Monod expression: ',dumstring(1:lss)
                      WRITE(*,*)
                      READ(*,*)
                      STOP
                    END IF
                  END IF
                END DO
              END IF
              IF (speciesfound) THEN
                CONTINUE
              ELSE  
                WRITE(*,*)
                WRITE(*,*) ' Cannot find species in Monod expression: ',dumstring(1:lss)
                dumstring = NAME(kk)
                CALL stringlen(dumstring,llen)
                WRITE(*,*) ' Looking for species: ',dumstring(1:llen)
                WRITE(*,*)
                READ(*,*)
                STOP
              END IF
            END DO    !  End of loop through "monod" species
            
          ELSE
            WRITE(*,*)
            WRITE(*,*) ' Error in database: Monod terms '
            WRITE(*,*) ' Mineral = ',tempmin(1:lsave)
            WRITE(*,*) ' Label = ',templabel(1:lenlabel)
            WRITE(*,*)
            READ(*,*)
            STOP
          END IF
        ELSE
          WRITE(*,*)
          WRITE(*,*) '  Error in database: Monod terms '
          WRITE(*,*) ' Mineral = ',tempmin(1:lsave)
          WRITE(*,*) ' Label = ',templabel(1:lenlabel)
          WRITE(*,*)
          READ(*,*)
          STOP
        END IF
        
      END IF        ! If input by user, skip block above
      
!  Look for inhibition terms
      
      READ(iunitDatabase,'(a)',END=300) dummy1
      IF (iuser(6) == 0) THEN ! No user input, so read
        id = 1
        iff = mls
        CALL sschaine(dummy1,id,iff,ssch,ids,ls)
        lzs=ls
        CALL convan(ssch,lzs,res)
        IF (ssch == 'inhibition') THEN
          ninhibit(np,nkin) = 0
          id = ids + ls
          CALL sschaine(dummy1,id,iff,ssch,ids,ls)
          IF (ssch == ':') THEN
            444             id = ids + ls
            CALL sschaine(dummy1,id,iff,ssch,ids,ls)
            IF (ls /= 0) THEN
              ninhibit(np,nkin) = ninhibit(np,nkin) + 1
              IF (ninhibit(np,nkin) > mpre) THEN
                WRITE(*,*)
                WRITE(*,*) ' Redimension parameter "mpre" in params.inc'
                WRITE(*,*) ' Mpre = ',mpre
                WRITE(*,*) ' Number of inhibition terms', ninhibit(np,nkin)
                WRITE(*,*)
                READ(*,*)
                STOP
              END IF
              lzs = ls
!                    call convan(ssch,lzs,res)
              CALL stringtype(ssch,lzs,res)
              NAME(ninhibit(np,nkin)) = ssch
              id = ids + ls
              CALL sschaine(dummy1,id,iff,ssch,ids,ls)
              lzs = ls
              CALL stringtype (ssch,lzs,res)
              IF (res == 'a') THEN
                WRITE(*,*)
                WRITE(*,*) '  Error in database: Monod inhibition terms'
                WRITE(*,*) ' Mineral = ',tempmin(1:lsave)
                WRITE(*,*) ' Label = ',templabel(1:lenlabel)
                WRITE(*,*)
                READ(*,*)
                STOP
              ELSE
                rinhibit_tmp(ninhibit(np,nkin)) = DNUM(ssch)
              END IF
              GO TO 444
            END IF
            
            DO kk = 1,ninhibit(np,nkin)
              dumstring = NAME(kk)
              IF (dumstring(1:4) == 'tot_') THEN
                CALL stringlen(dumstring,lss)
                NAME(kk) = dumstring(5:lss)
                itot_tmp = 1
              ELSE
                itot_tmp = 0
              END IF

              speciesfound = .false.          
              DO i = 1,ncomp
                IF (ulab(i) == NAME(kk)) THEN
                  inhibit(kk,np,nkin) = i
                  speciesfound = .true.
                  rinhibit(kk,np,nkin) = rinhibit_tmp(kk)
                  IF (itot_tmp == 1) THEN
                    itot_inhibit(kk,np,nkin) = 1
                  END IF
                END IF
              END DO
              IF (speciesfound) THEN
                CONTINUE
              ELSE
                DO ksp = 1,nspec
                  ik = ncomp + ksp
                  IF (ulab(ik) == NAME(kk)) THEN
                    inhibit(kk,np,nkin) = ik
                    speciesfound = .true.
                    rinhibit(kk,np,nkin) = rinhibit_tmp(kk)
                    IF (itot_tmp == 0) THEN
                      CONTINUE
                    ELSE
                      WRITE(*,*)
                      WRITE(*,*) ' Cannot specify total concentration for secondary species'
                      WRITE(*,*) ' in Monod expression: ',dumstring(1:lss)
                      WRITE(*,*)
                      READ(*,*)
                      STOP
                    END IF
                  END IF
                END DO
                IF (speciesfound) THEN
                  CONTINUE
                ELSE
                  WRITE(*,*)
                  WRITE(*,*) ' Cannot find species in Monod expression: ',dumstring(1:lss)
                  dumstring = NAME(kk)
                  CALL stringlen(dumstring,llen)
                  WRITE(*,*) ' Looking for species: ',dumstring(1:llen)
                  WRITE(*,*)
                  READ(*,*)
                  STOP
                END IF
              END IF
            END DO    !  End of loop through "inhibiting" species
            
          ELSE
            WRITE(*,*)
            WRITE(*,*) ' Error in database: Monod inhibition terms '
            WRITE(*,*) ' Mineral = ',tempmin(1:lsave)
            WRITE(*,*) ' Label = ',templabel(1:lenlabel)
            WRITE(*,*)
            READ(*,*)
            STOP
          END IF
        ELSE
          WRITE(*,*)
          WRITE(*,*) '  Error in database: Monod inhibition terms '
          WRITE(*,*) ' Mineral = ',tempmin(1:lsave)
          WRITE(*,*) ' Label = ',templabel(1:lenlabel)
          WRITE(*,*)
          READ(*,*)
          STOP
        END IF
        
      END IF                                                !! If input by user, skip block above

! biomass
    ELSE IF (imintype(np,nkin) == 8) THEN                                   !! Read Monod + thermodynamic parameters
      
!!    Look for Monod terms (species specification with half-saturation constants)
      
      READ(iunitDatabase,'(a)',END=300) dummy1
      IF (iuser(5) == 0) THEN ! No user input, so read
        id = 1
        iff = mls
        CALL sschaine(dummy1,id,iff,ssch,ids,ls)
        lzs=ls
        CALL convan(ssch,lzs,res)
        IF (ssch == 'monod_term' .OR. ssch == 'monod_terms') THEN
          nmonod(np,nkin) = 0
          id = ids + ls
          CALL sschaine(dummy1,id,iff,ssch,ids,ls)
          IF (ssch == ':') THEN
            666             id = ids + ls
            CALL sschaine(dummy1,id,iff,ssch,ids,ls)
            IF (ls /= 0) THEN
              nmonod(np,nkin) = nmonod(np,nkin) + 1
              IF (nmonod(np,nkin) > mpre) THEN
                WRITE(*,*)
                WRITE(*,*) ' Redimension parameter "mpre" in params.inc'
                WRITE(*,*) ' Mpre = ',mpre
                WRITE(*,*) ' Number of Monod terms = ', nmonod(np,nkin)
                WRITE(*,*)
                READ(*,*)
                STOP
              END IF
              lzs = ls
!                    call convan(ssch,lzs,res)
              CALL stringtype(ssch,lzs,res)
              NAME(nmonod(np,nkin)) = ssch
              id = ids + ls
              CALL sschaine(dummy1,id,iff,ssch,ids,ls)
              lzs = ls
              CALL stringtype (ssch,lzs,res)
              IF (res == 'a') THEN
                WRITE(*,*)
                WRITE(*,*) '  Error in database: Monod terms'
                WRITE(*,*) ' Mineral = ',tempmin(1:lsave)
                WRITE(*,*) ' Label = ',templabel(1:lenlabel)
                WRITE(*,*)
                READ(*,*)
                STOP
              ELSE
                halfsat_tmp(nmonod(np,nkin)) = DNUM(ssch)
              END IF
              GO TO 666
            END IF
            
            DO kk = 1,nmonod(np,nkin)
              dumstring = NAME(kk)
              IF (dumstring(1:4) == 'tot_') THEN
                CALL stringlen(dumstring,lss)
                NAME(kk) = dumstring(5:lss)
                itot_tmp = 1
              ELSE
                itot_tmp = 0
              END IF
              IF (NAME(kk) == 'self' .OR. NAME(kk) == 'Self' .OR. NAME(kk) == 'SELF') THEN           !!  Rate depends on mineral concentration, mol/L

                speciesfound = .TRUE.
                kmonod(kk,np,k) = 1                                  !!  Flag to indicate that Monod rate depends on the mineral concentration itself
                imonod(kk,np,k) = 0                                  !!  If cycling through species, this should tell the code (in Reaction.F90) to check kmonod(kk,np,k)
                halfsat(kk,np,nkin) = halfsat_tmp(kk)            

              ELSE
!!              Check to see if name in Monod term is recognized as a species
                speciesfound = .false.
                DO i = 1,ncomp                                  !! Loop over primary species
                  IF (ulab(i) == NAME(kk)) THEN
                    imonod(kk,np,nkin) = i
                    speciesfound = .true.
                    halfsat(kk,np,nkin) = halfsat_tmp(kk)
                    IF (itot_tmp == 1) THEN
                      itot_monod(kk,np,nkin) = 1
                    END IF
                  END IF
                END DO
              END IF
              IF (speciesfound) THEN
                CONTINUE
              ELSE
                DO ksp = 1,nspec                              !! Loop over secondary species
                  ik = ncomp + ksp
                  IF (ulab(ik) == NAME(kk)) THEN
                    imonod(kk,np,nkin) = ik
                    speciesfound = .true.
                    halfsat(kk,np,nkin) = halfsat_tmp(kk)
                    IF (itot_tmp == 0) THEN
                      CONTINUE
                    ELSE
                      WRITE(*,*)
                      WRITE(*,*) ' Cannot specify total concentration for secondary species'
                      WRITE(*,*) ' in Monod expression: ',dumstring(1:lss)
                      WRITE(*,*)
                      READ(*,*)
                      STOP
                    END IF
                  END IF
                END DO
              END IF
              IF (speciesfound) THEN
                CONTINUE
              ELSE  
                WRITE(*,*)
                WRITE(*,*) ' Cannot find species in Monod expression: ',dumstring(1:lss)
                dumstring = NAME(kk)
                CALL stringlen(dumstring,llen)
                WRITE(*,*) ' Looking for species: ',dumstring(1:llen)
                WRITE(*,*)
                READ(*,*)
                STOP
              END IF
            END DO    !  End of loop through "monod" species
            
          ELSE
            WRITE(*,*)
            WRITE(*,*) ' Error in database: Monod terms '
            WRITE(*,*) ' Mineral = ',tempmin(1:lsave)
            WRITE(*,*) ' Label = ',templabel(1:lenlabel)
            WRITE(*,*)
            READ(*,*)
            STOP
          END IF
        ELSE
          WRITE(*,*)
          WRITE(*,*) '  Error in database: Monod terms '
          WRITE(*,*) ' Mineral = ',tempmin(1:lsave)
          WRITE(*,*) ' Label = ',templabel(1:lenlabel)
          WRITE(*,*)
          READ(*,*)
          STOP
        END IF

!       commented out, now read at read_CatabolicPath        
!!        READ(iunitDatabase,*) chi(np,nkin),bq(np,nkin),direction(np,nkin)

      END IF        ! If input by user, skip block above
! biomass end

! sergi: biomass decay
  ELSE IF (imintype(np,nkin) == 9) THEN
    
!     first time around allocation
      if (.not.allocated(biomass_decay)) then

        allocate(biomass_decay(mreact,mrct))

      end if 

!     read 'biomass' namelist from file
      READ(iunitDatabase,nml=BiomassDecay,iostat=ios)     

!     result from read (IOS)
      if (ios == 0) then

!       successful read

!       find the biomass for this reaction in the mineral list
        ib = 0
        count_ib = 0
        do_biomass: do im=1,len(namrl(im))

          if (im == 1) then
            count_ib = count_ib + 1
          elseif (namrl(im-1) /= namrl(im)) then
              count_ib = count_ib + 1
          endif
            
          if (biomass == namrl(im)) then
               
            ib = count_ib
            exit do_biomass
               
          end if
            
        end do do_biomass
            
        if (ib == 0) then
          write(*,*)'biomass ',biomass,' for reaction: ',name,' is not in the list of minerals'
          write(*,*)'biomass decay needs to be after biomass has been defined!!'
          stop
        end if
        
        biomass_decay(np,nkin)=ib
      
      else if (ios < 0) then

!       no namelist to read
        write(*,*)'namelist not found. bye'
        stop

      else if (ios > 0) then

        write(*,*)'error reading namelist: goodbye'
        write(*,nml=BiomassDecay,iostat=ios)
        stop

      end if

! end sergi: biomass decay

      
    ELSE                      !! If not a Monod rate law, look for species (exponential) dependences
      
      READ(iunitDatabase,'(a)',END=300) dummy1
      IF (iuser(3) == 0) THEN ! No user input, so read
        id = 1
        iff = mls
        CALL sschaine(dummy1,id,iff,ssch,ids,ls)
        lzs=ls
        CALL convan(ssch,lzs,res)
        IF (ssch == 'dependence') THEN
          ndepend(np,nkin) = 0
          id = ids + ls
          CALL sschaine(dummy1,id,iff,ssch,ids,ls)
          
          IF (ssch == ':') THEN
              
            222             id = ids + ls
                     
!!   ***********************************************
!!    Now read string right after ":"

            CALL sschaine(dummy1,id,iff,ssch,ids,ls)
            IF (ls /= 0) THEN
              IF (ssch == 'HyperbolicInhibition') THEN
                  
                HyperbolicInhibition(np,nkin) = .TRUE.
!!                id = ids + ls
!!                CALL sschaine(dummy1,id,iff,ssch,ids,ls)
!!                lzs = ls
!!                CALL stringtype(ssch,lzs,res)     

                id = ids + ls
                CALL sschaine(dummy1,id,iff,ssch,ids,ls)
                lzs = ls
                Kformation(np,nkin) = DNUM(ssch)
                
                id = ids + ls
                CALL sschaine(dummy1,id,iff,ssch,ids,ls)
                lzs = ls
                HyperbolicInhibitionName(np,nkin) = ssch
                
                id = ids + ls
                CALL sschaine(dummy1,id,iff,ssch,ids,ls)
                lzs = ls
                HyperbolicInhibitionDepend(np,nkin) = DNUM(ssch)                
                
!! Now read first the K-formation term, then the species name, then the exponential dependence
!!                READ(iunitDatabase,*) Kformation(np,nkin),HyperbolicInhibitionName(np,nkin),HyperbolicInhibitionDepend(np,nkin)

                speciesfound = .false.
                DO i = 1,ncomp
                  IF (ulab(i) == HyperbolicInhibitionName(np,nkin)) THEN
                    speciesfound = .true.
                    HyperbolicInhibitionPointer(np,nkin) = i
                  END IF
                END DO
                IF (speciesfound) THEN
                  CONTINUE
                ELSE
                  WRITE(*,*)
                  WRITE(*,*) ' Primary species not found for HyperbolicInhibition term in Dependence'
                  WRITE(*,*) ' Looking for: ', HyperbolicInhibitionName(np,nkin)
                  WRITE(*,*)
                  READ(*,*)
                  STOP
                END IF
                
                GO TO 222            !!  Go back, even if this is obsolete syntax
                
              ELSE   
                ndepend(np,nkin) = ndepend(np,nkin) + 1
                IF (ndepend(np,nkin) > mpre) THEN
                  WRITE(*,*)
                  WRITE(*,*) ' Redimension parameter "mpre" in params.inc'
                  WRITE(*,*) ' Mpre = ',mpre
                  WRITE(*,*) ' Ndepend = ',ndepend(np,nkin)
                  WRITE(*,*)
                  READ(*,*)
                  STOP
                END IF
                lzs = ls
                CALL stringtype(ssch,lzs,res)
                NAME(ndepend(np,nkin)) = ssch
                id = ids + ls
                CALL sschaine(dummy1,id,iff,ssch,ids,ls)
                lzs = ls
                CALL stringtype (ssch,lzs,res)
                IF (res == 'a') THEN
                  WRITE(*,*)
                  WRITE(*,*) '  Error in database: Reaction "dependence" '
                  WRITE(*,*) ' Mineral = ',tempmin(1:lsave)
                  WRITE(*,*) ' Label = ',templabel(1:lenlabel)
                  WRITE(*,*)
                  READ(*,*)
                  STOP
                ELSE
                  dep_tmp(ndepend(np,nkin)) = DNUM(ssch)
                END IF
                
                GO TO 222            !!  Go back, even if this is obsolete syntax
                
              END IF 
              
            END IF           !! End of ls /= 0 block
            
            DO kk = 1,ndepend(np,nkin)
              dumstring = NAME(kk)
              IF (dumstring(1:4) == 'tot_') THEN
                CALL stringlen(dumstring,lss)
                NAME(kk) = dumstring(5:lss)
                itot_tmp = 1
              ELSE
                itot_tmp = 0
              END IF
              speciesfound = .false.
              DO i = 1,ncomp
                IF (ulab(i) == NAME(kk)) THEN
                  idepend(kk,np,nkin) = i
                  speciesfound = .true.
                  namdep_nyf(kk,np,nkin) = 'found'
                  depend(kk,np,nkin) = dep_tmp(kk)
                  IF (itot_tmp == 1) THEN
                    itot_min(kk,np,nkin) = 1
                  END IF
                END IF
              END DO
              IF (speciesfound) THEN
                CONTINUE
              ELSE
                IF (itot_tmp == 1) THEN
                  WRITE(*,*)
                  WRITE(*,*) ' Cannot specify total concentration for other than'
                  WRITE(*,*) ' primary species in rate law: ',dumstring(1:lss)
                  WRITE(*,*)
                  READ(*,*)
                  STOP
                END IF
                
!  If species is not among primary species, search later from list of secondary aqueous,
!     exchange, and surface complex species (exchange and surface complexes not yet read)
                
                namdep_nyf(kk,np,nkin) = NAME(kk)
                depend(kk,np,nkin) = dep_tmp(kk)
                
                
!                      do ksp = 1,nspec
!                        ik = ncomp + ksp
!                        if (ulab(ik).eq.name(kk)) then
!                          idepend(kk,np,nkin) = ik
!                          speciesfound = .true.
!                          depend(kk,np,nkin) = dep_tmp(kk)
!                          if (itot_tmp.eq.0) then
!                            continue
!                          else
!                            write(*,*)
!                            write(*,*) ' Cannot specify total concentration for secondary species'
!                            write(*,*) ' in rate law: ',dumstring(1:lss)
!                            write(*,*)
!                            stop
!                          endif
!                        endif
!                      end do
                
!                      if (speciesfound) then
!                        continue
!                      else
!                        write(*,*)
!                        write(*,*) ' Cannot find species in rate law: ',dumstring(1:lss)
!                        dumstring = name(kk)
!                        call stringlen(dumstring,llen)
!                        write(*,*) ' Looking for species: ',dumstring(1:llen)
!                        write(*,*)
!                        stop
!                      endif
                
              END IF
              
            END DO    !  End of loop through "ndepend" species
            
          ELSE
            WRITE(*,*)
            WRITE(*,*) ' Error in database: Reaction "dependence" '
            WRITE(*,*) ' Mineral = ',tempmin(1:lsave)
            WRITE(*,*) ' Label = ',templabel(1:lenlabel)
            WRITE(*,*)
            READ(*,*)
            STOP
          END IF
        ELSE
          WRITE(*,*)
          WRITE(*,*) '  Error in database: Reaction "dependence" '
          WRITE(*,*) ' Mineral = ',tempmin(1:lsave)
          WRITE(*,*) ' Label = ',templabel(1:lenlabel)
          WRITE(*,*)
          READ(*,*)
          STOP
        END IF
        
      END IF                                  !!  User input, so skip block above
      
      READ(iunitDatabase,'(a)',END=300) dummy1
      IF (dummy1(1:1) == '+') THEN
        CONTINUE
      ELSE
        id = 1
        iff = mls
        CALL sschaine(dummy1,id,iff,ssch,ids,ls)
        lzs=ls
        CALL convan(ssch,lzs,res)
        IF (ssch == 'affinitydependence') THEN
          id = ids + ls
          CALL sschaine(dummy1,id,iff,ssch,ids,ls)
          IF (ssch == '=') THEN                               !! OK, read 3 numbers for affinity dependence

            id = ids + ls
            CALL sschaine(dummy1,id,iff,ssch,ids,ls)
            lzs = ls
            CALL stringtype(ssch,lzs,res)
            IF (res == 'a') THEN
              WRITE(*,*)
              WRITE(*,*) ' "AffinityDependence = " should be followed by a number'
              WRITE(*,*) ' Mineral = ',namdum(1:lsave)
              WRITE(*,*) ' String = ',ssch(1:ls)
              WRITE(*,*)
              READ(*,*)
              STOP
            ELSE
              AffinityDepend1(np,nkin) = DNUM(ssch)                      !! Corresponds to Hellmann's "m2"
            END IF

            id = ids + ls
            CALL sschaine(dummy1,id,iff,ssch,ids,ls)
            lzs = ls
            CALL stringtype(ssch,lzs,res)
            IF (res == 'a') THEN
              WRITE(*,*)
              WRITE(*,*) ' "AffinityDependence = " should be followed by two numbers--only first one was read'
              WRITE(*,*) ' Mineral = ',namdum(1:lsave)
              WRITE(*,*) ' String = ',ssch(1:ls)
              WRITE(*,*)
              READ(*,*)
              STOP
            ELSE
              AffinityDepend2(np,nkin) = DNUM(ssch)                    !! Corresponds to Hellmann's "n"
            END IF

            id = ids + ls
            CALL sschaine(dummy1,id,iff,ssch,ids,ls)
            lzs = ls
            CALL stringtype(ssch,lzs,res)
            IF (res == 'a') THEN
              WRITE(*,*)
              WRITE(*,*) ' "AffinityDependence = " should be followed by three numbers--only first two were read'
              WRITE(*,*) ' Mineral = ',namdum(1:lsave)
              WRITE(*,*) ' String = ',ssch(1:ls)
              WRITE(*,*)
              READ(*,*)
              STOP
            ELSE
              AffinityDepend3(np,nkin) = DNUM(ssch)                    !! !! Corresponds to Hellmann's m1
            END IF
          
          ELSE
            WRITE(*,*)
            WRITE(*,*) ' Error in database: "AffinityDependence" '
            WRITE(*,*) ' Mineral = ',tempmin(1:lsave)
            WRITE(*,*) ' Label = ',templabel(1:lenlabel)
            WRITE(*,*)
            READ(*,*)
            STOP
          END IF
        ELSE
          CONTINUE
        END IF
      END IF                                                           !! End of read of AffinityDependence

!! Check to see if all of the exponents /= 1, in which case we have the Burch-Hellmann rate law
!!!      IF (AffinityDepend1(np,nkin) /= 1.0d0 .AND. AffinityDepend2(np,nkin) /= 1.0d0 .AND. AffinityDepend3(np,nkin) /= 1.0d0) THEN
      IF (AffinityDepend3(np,nkin) /= 1.0d0) THEN
        JacobianNumerical = .TRUE.
      END IF

    END IF                                                             !! End of choice between "dependence" and "monod_term"
    
!  Read stoichiometry here
    
!          READ(iunitDatabase,'(a)',end=300) dummy1
!          namsearch = nammin_tmp(nkin)
!          call stringlen(namsearch,ls_min)
!          templabel = rlabel(np,nkin)
!          call stringlen(templabel,ls_label)
!          call reaction_interp(dummy1,npar,nkin,namsearch,
!     &       templabel,namspecies,nsp,sto)
    
! Check to see if all of the species are present, skipping the mineral
    
!          kkh2o = 0
!          do kk = 2,nsp
!            do ik = 1,ncomp+nspec
!              if (namspecies(kk).eq.ulab(ik)) then
!                goto 550
!              endif
!            end do
!            if (namspecies(kk).eq.'H2O'.or.namspecies(kk).eq.'h2o') then
!                ikh2O = ik
!                goto 550
!              endif
!            namtemp = namspecies(kk)
!            call stringlen(namtemp,ls)
!            write(*,*)
!            write(*,*) ' Could not find species: ',namtemp(1:ls)
!            write(*,*) ' In mineral reaction:    ',namsearch(1:ls_min)
!            write(*,*) ' Label:                  ',templabel(1:ls_label)
!            write(*,*)
!            stop
!  550       continue
!          end do
!          nreac = nreac + 1
    
    
!          ax(nreac,nreac) = 1.0
    
!     do kk = 2,nsp
    
!       do i = 1, ncomp
!         if (namspecies(kk).eq.ulab(i)) then
!           bx(nreac,i) = sto(kk)
!           goto 1001
!              endif
!            end do
    
!       do ksp = 1, nspec
!         if (namspecies(kk).eq.ulab(ksp+ncomp)) then
!           ax(nreac,ksp) = -sto(kk)
!           goto 1001
!              endif
!            end do
    
!       do ngg = 1, ngas
!         if (namspecies(kk).eq.namg(ngg)) then
!           ax(nreac,ngg+nspec) = -sto(kk)
!           goto 1001
!              endif
!            end do
    
! 1001       continue
    
!          end do
    
  ELSE            ! Mineral name not followed by "label"
    WRITE(*,*)
    WRITE(*,*) ' Error in database in "label" statement'
    WRITE(*,*) ' Mineral = ',tempmin(1:lsave)
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
  
ELSE
  CALL breakfind
  GO TO 475
END IF

GO TO 100








2001 FORMAT(20(1X,f7.2))

DEALLOCATE(name)
DEALLOCATE(namspecies)
DEALLOCATE(dep_tmp)
DEALLOCATE(iuseR)
DEALLOCATE(rinhibit_tmp)
DEALLOCATE(halfsat_tmp)

i = size(namAssociate,1)
ALLOCATE(workchar1(i))
workchar1 = namAssociate
DEALLOCATE(namAssociate)
ALLOCATE(namAssociate(nkin))
IF(nkin /= 0) namAssociate(1:nkin) = workchar1(1:nkin)
DEALLOCATE(workchar1)

! mineral type
i = size(mintype,1)
ALLOCATE(workint(i))
workint = mintype
DEALLOCATE(mintype)
ALLOCATE(mintype(nkin+50))
IF(nkin /= 0) mintype(1:nkin+50) = workint(1:nkin+50)
DEALLOCATE(workint)

!!  Check that the "associated" minerals are in the list

IF (ALLOCATED(MineralID)) THEN
  DEALLOCATE(MineralID)
  ALLOCATE(MineralID(nkin))
ELSE
  ALLOCATE(MineralID(nkin))
END IF


MineralFound = .FALSE.
DO k = 1,nkin
  MineralID(k) = k
  IF (namAssociate(k) /= ' ') THEN         !! Check to make sure the mineral is in the list
    DO kk = 1,nkin
      IF (nammin_tmp(kk) == namAssociate(k)) THEN
        MineralFound = .TRUE.
        MineralID(k) = kk
      END IF
    END DO
    IF (MineralFound) THEN
      CONTINUE
    ELSE
      dumstring = namAssociate(k)
      CALL stringlen(dumstring,lss)
      tempmin = nammin_tmp(k)
      CALL stringlen(tempmin,lsave)
      WRITE(*,*)
      WRITE(*,*) ' Cannot find associated mineral for update of volume fraction '
      WRITE(*,*) ' Searching for mineral: ',dumstring(1:lss)
      WRITE(*,*) ' Associated with mineral: ', tempmin(1:lsave)
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
  ELSE
    CONTINUE
  END IF
END DO

IF (ALLOCATED(MineralAssociate)) THEN
  DEALLOCATE(MineralAssociate)
  ALLOCATE(MineralAssociate(nkin))
ELSE
  ALLOCATE(MineralAssociate(nkin))
END IF

DO k = 1,nkin
  IF (namAssociate(k) /= ' ') THEN
    MineralAssociate(k) = .TRUE.
  ELSE
    MineralAssociate(k) = .FALSE.
  END IF 
END DO


RETURN

556 WRITE(*,*) ' Error in reading rate dependence in database'
WRITE(*,*) ' Mineral = ',tempmin(1:lsave)
WRITE(*,*) ' Label = ',templabel(1:lenlabel)
WRITE(*,*)
READ(*,*)
STOP
334 WRITE(*,*) ' Error in opening database file'
WRITE(*,*) ' Looking for: ',data1
WRITE(*,*)
READ(*,*)
STOP

300 CLOSE(UNIT=iunitDatabase,STATUS='keep')