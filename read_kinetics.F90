!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 10:05:09
 
!************** (C) COPYRIGHT 1995,1998,1999 ******************
!*******************     C.I. Steefel      *******************
!                    All Rights Reserved

!  GIMRT98 IS PROVIDED "AS IS" AND WITHOUT ANY WARRANTY EXPRESS OR IMPLIED.
!  THE USER ASSUMES ALL RISKS OF USING GIMRT98. THERE IS NO CLAIM OF THE
!  MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.

!  YOU MAY MODIFY THE SOURCE CODE FOR YOUR OWN USE, BUT YOU MAY NOT
!  DISTRIBUTE EITHER THE ORIGINAL OR THE MODIFIED CODE TO ANY OTHER
!  WORKSTATIONS
!**********************************************************************

SUBROUTINE read_kinetics(nout,ncomp,nspec,nrct,ikin,data1)
USE crunchtype
USE params
USE concentration
USE mineral
USE strings

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nout
INTEGER(I4B), INTENT(IN)                                    :: ncomp
INTEGER(I4B), INTENT(IN)                                    :: nspec
INTEGER(I4B), INTENT(IN)                                    :: nrct
INTEGER(I4B), INTENT(OUT)                                   :: ikin
CHARACTER (LEN=mls), INTENT(IN)                             :: data1

!  Internal variables and arrays

LOGICAL(LGT)                                                :: speciesfound
CHARACTER (LEN=mls)                                         :: dummy1
CHARACTER (LEN=mls)                                         :: dummy2
CHARACTER (LEN=mls)                                         :: dumstring
CHARACTER (LEN=mls)                                         :: dumlabel
CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE              :: name

REAL(DP), DIMENSION(:), ALLOCATABLE                         :: skin
REAL(DP), DIMENSION(:), ALLOCATABLE                         :: dep_tmp
REAL(DP), DIMENSION(:), ALLOCATABLE                         :: rinhibit_tmp
REAL(DP), DIMENSION(:), ALLOCATABLE                         :: halfsat_tmp

CHARACTER (LEN=23)                                          :: str_endkin
CHARACTER (LEN=mls)                                         :: reactiontype

INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs
INTEGER(I4B)                                                :: i
INTEGER(I4B)                                                :: lsave
INTEGER(I4B)                                                :: nspecies
INTEGER(I4B)                                                :: j
INTEGER(I4B)                                                :: ii
INTEGER(I4B)                                                :: ll
INTEGER(I4B)                                                :: nk
INTEGER(I4B)                                                :: kk
INTEGER(I4B)                                                :: lss
INTEGER(I4B)                                                :: itot_tmp
INTEGER(I4B)                                                :: ksp
INTEGER(I4B)                                                :: llen
INTEGER(I4B)                                                :: ik
INTEGER(I4B)                                                :: k

REAL(DP), PARAMETER                                         :: tiny=1.E-12
REAL(DP)                                                    :: charge

ALLOCATE(rinhibit_tmp(mcomp+mspec))
ALLOCATE(halfsat_tmp(mcomp+mspec))
ALLOCATE(name(mpre))
ALLOCATE(skin(50))
ALLOCATE(dep_tmp(50))

str_endkin = 'End of aqueous kinetics'
OPEN(UNIT=18,FILE=data1,STATUS='old')
REWIND nout
ikin = 0

100 READ(nout,'(a)',END=300) zone
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  CALL stringtype(ssch,lzs,res)
  dummy2(1:ls) = ssch(1:ls)
  lsave = ls
ELSE
  GO TO 100
END IF

200 READ(18,'(a)') dummy1
IF (dummy1 == 'Begin aqueous kinetics') THEN
  GO TO 400
ELSE
  GO TO 200
END IF

400 READ(18,'(a)',END=300) dummy1
IF (dummy1 == str_endkin) THEN
  WRITE(*,*)
  WRITE(*,*) ' Aqueous kinetic reaction not found in database'
  WRITE(*,*) ' Looking for ',dummy2(1:lsave)
  WRITE(*,*) ' STOP'
  WRITE(*,*)
  READ(*,*)
  STOP
END IF

id = 1
iff = mls
CALL sschaine(dummy1,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  CALL stringtype(ssch,lzs,res)
  IF (res == 'a') THEN
    IF (ssch(1:ls) == dummy2(1:ls)) THEN
      ikin = ikin + 1
      IF (ikin > maqkin) THEN
        WRITE(*,*)
        WRITE(*,*) ' Parameter "maqkin" dimensioned too small'
        WRITE(*,*) ' Number of aqueous kinetic reactions = ',ikin
        WRITE(*,*) ' Dimension for aqueous kinetic reactions = ', maqkin
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
!!      WRITE(*,*)
!!      WRITE(*,*) ' Kinetic reaction found'
!!      WRITE(*,*) dummy2(1:ls)
!!      WRITE(*,*)
      BACKSPACE 18
      READ(18,*) namkin(ikin),nreactkin(ikin),reactiontype,nspecies,  &
          (skin(j),NAME(j),j=1,nspecies),keqkin(ikin)
      IF (nreactkin(ikin) <= 0) THEN
        WRITE(*,*)
        WRITE(*,*) ' Rate law label given, but no'
        WRITE(*,*) '   rate laws included in database'
        WRITE(*,*) dummy2(1:ls)
        WRITE(*,*)
        STOP
      END IF
      
      IF (nreactkin(ikin) > mreact) THEN
        WRITE(*,*)
        WRITE(*,*) ' Redimension parameter mreact in params.F90'
        WRITE(*,*) ' Mreact =  ',mreact
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
      
!!      WRITE(*,*) ' namkin = ', namkin(ikin)
      IF (reactiontype == 'tst') THEN
        WRITE(*,505) namkin(ikin)
        iaqtype(ikin) = 1
      ELSE IF (reactiontype == 'monod') THEN
        WRITE(*,506) namkin(ikin)
        IF (nreactkin(ikin) > 1) THEN
          WRITE(*,*)
          WRITE(*,*) ' Only one parallel Monod reaction allowed for the moment'
          WRITE(*,*)
          READ(*,*)
          STOP
        END IF
        iaqtype(ikin) = 2
      ELSE IF (reactiontype == 'irreversible') THEN
        WRITE(*,507) namkin(ikin)
        iaqtype(ikin) = 3
      ELSE IF (reactiontype == 'radioactive_decay') THEN
        WRITE(*,508) namkin(ikin)
        iaqtype(ikin) = 4

! biomass
      ELSE IF (reactiontype == 'monodbiomass') THEN
        WRITE(*,506) namkin(ikin)
        IF (nreactkin(ikin) > 1) THEN
          WRITE(*,*)
          WRITE(*,*) ' Only one parallel Monod reaction allowed for the moment'
          WRITE(*,*)
          READ(*,*)
          STOP
        END IF
        iaqtype(ikin) = 8
! biomass end  

      ELSE
        WRITE(*,*) ' Aqueous kinetic reaction type not recognized'
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF

505 FORMAT(1x,'Using TST rate formulation for reaction: ',a25)
506 FORMAT(1x,'Using Monod rate formulation for reaction: ',a25)
507 FORMAT(1x,'Using irreversible rate formulation for reaction: ',a25)
508 FORMAT(1x,'Using radioactive decay formulation for reaction: ',a25)


!!      WRITE(*,*) ' nreactkin = ',nreactkin(ikin)
!!      WRITE(*,*) ' nspecies = ',nspecies
!!      DO ii = 1,nspecies
!!        WRITE(*,*) skin(ii)
!!      END DO
!!      WRITE(*,*) ' keqkin = ',keqkin(ikin)
      
!  Construct stoichiometric coefficient matrix for this particular reaction
!  NOTE: We assume at this stage that the reaction is written in terms of
!    primary species only.  If not, give an error message and get out.
      
      DO i = 1,ncomp
        mukin(ikin,i) = 0.0
      END DO
      DO j = 1,nspecies
        IF (NAME(j) == 'h2o' .OR. NAME(j) == 'H2O') THEN
          GO TO 500    ! Don't require that H2O be in the primary species
        END IF
        speciesfound = .false.
!   Search for species among primary species list
        DO i = 1,ncomp
          IF (ulab(i) == NAME(j)) THEN
            mukin(ikin,i) = skin(j)
            speciesfound = .true.
          END IF
        END DO
        IF (speciesfound) THEN
          CONTINUE
        ELSE
          WRITE(*,*)
          WRITE(*,*) ' Species in kinetic reaction not found in primary species list'
          WRITE(*,*) ' Reaction = ',dummy2
          WRITE(*,*) ' Species = ',NAME(j)
          WRITE(*,*)
          READ(*,*)
          STOP
        END IF
        500          CONTINUE
      END DO
      
!!    Check electrical balance in reactions
      charge = 0.00
      DO i = 1,ncomp
        charge = charge + mukin(ikin,i)*chg(i)
      END DO
      IF (ABS(charge) > tiny) THEN
        WRITE(*,*)
        WRITE(*,501) 
        WRITE(*,502) charge
        WRITE(*,503) namkin(ikin)
        WRITE(*,*)
        STOP
      ENDIF

501 FORMAT(1X,'Aqueous kinetic reaction not electrically balanced')
502 FORMAT(1x,'Charge = ',1pe10.2)
503 FORMAT(1X,'Reaction label: ',a15)

      IF (iaqtype(ikin) == 1 .OR. iaqtype(ikin) == 3 .OR. iaqtype(ikin) == 4) THEN
        DO ll = 1,nreactkin(ikin)
          DO i = 1,ncomp
            itot(i,ll,ikin) = 0
          END DO
          READ(18,*) ratek(ll,ikin),nk,(NAME(kk), dep_tmp(kk),kk=1,nk)
!!          WRITE(*,*) ratek(ll,ikin)
!!          DO kk = 1,nk
!!            WRITE(*,*) dep_tmp(kk)
!!          END DO
          DO i = 1,ncomp
            dependk(i,ll,ikin) = 0.0
          END DO
          DO kk = 1,nk
            IF (kk > mpre) THEN
              WRITE(*,*)
              WRITE(*,*) ' Redimension parameter mpre in params.inc'
              WRITE(*,*) ' Mpre = ',mpre
              WRITE(*,*)
              READ(*,*)
              STOP
            END IF
!  First, check to see if we are dealing with a total concentration
            dumstring = NAME(kk)
            IF (dumstring(1:4) == 'tot_') THEN
              CALL stringlen(dumstring,lss)
              NAME(kk) = dumstring(5:lss)
              itot_tmp = 1
            ELSE
              itot_tmp = 0
            END IF
            speciesfound = .false.
            DO i = 1,ncomp                          !!  Searching through primary species list
              IF (ulab(i) == NAME(kk)) THEN
                speciesfound = .true.
                dependk(i,ll,ikin) = dep_tmp(kk)
                IF (itot_tmp == 1) THEN
                  itot(i,ll,ikin) = 1
                ELSE
                  itot(i,ll,ikin) = 0
                END IF
              END IF
            END DO
            IF (speciesfound) THEN
              CONTINUE
            ELSE
              WRITE(*,*)
              WRITE(*,*) ' Kinetic aqueous reaction includes a species not in the primary species list'
              WRITE(*,*) ' Reaction label = ',dummy2(1:ls)
              WRITE(*,*) ' Species not found = ', NAME(kk)
              WRITE(*,*)
              READ(*,*)
              STOP
            END IF
          END DO
          
        END DO   ! End of llth parallel rxn for ikinth overall kinetic rxn
        
      ELSE IF (iaqtype(ikin) == 2) THEN       !  Monod
        
!  Assume only one parallel reaction for the moment
        
        ll = 1
        DO i = 1,ncomp
          itot_monodaq(i,ikin) = 0
        END DO
        READ(18,*) ratek(ll,ikin),nmonodaq(ikin),(NAME(kk),  &
            halfsat_tmp(kk),kk=1,nmonodaq(ikin))
        DO kk = 1,nmonodaq(ikin)
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
              imonodaq(kk,ikin) = i
              speciesfound = .true.
              halfsataq(kk,ikin) = halfsat_tmp(kk)
              IF (itot_tmp == 1) THEN
                itot_monodaq(kk,ikin) = 1
              END IF
            END IF
          END DO
          IF (speciesfound) THEN
            CONTINUE
          ELSE
            DO ksp = 1,nspec
              ik = ncomp + ksp
              IF (ulab(ik) == NAME(kk)) THEN
                imonodaq(kk,ikin) = ik
                speciesfound = .true.
                halfsataq(kk,ikin) = halfsat_tmp(kk)
                IF (itot_tmp == 0) THEN
                  CONTINUE
                ELSE
                  WRITE(*,*)
                  WRITE(*,*) ' Cannot specify total concentration for secondary species'
                  WRITE(*,*) '   in Monod expression: ',dumstring(1:lss)
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
        END DO    !  End of loop through "monod" species
        
!  Now look for inhibition information
        
        READ(18,*) dumlabel
        IF (dumlabel /= 'Inhibition') THEN
          WRITE(*,*) ' For Monod law, there should be an entry for inhibition terms'
          WRITE(*,*)
          READ(*,*)
          STOP
        END IF
        BACKSPACE 18
        READ(18,*) dumlabel,ninhibitaq(ikin),  &
            (NAME(kk),rinhibit_tmp(kk),kk=1,ninhibitaq(ikin))
        DO kk = 1,ninhibitaq(ikin)
          dumstring = NAME(kk)
          IF (dumstring(1:4) == 'tot_') THEN
            CALL stringlen(dumstring,lss)
            NAME(kk) = dumstring(5:lss)
            itot_tmp = 1
          ELSE
            itot_tmp = 0
          END IF
          speciesfound = .false.
          DO i = 1,ncomp                         !! Search through primary species for inhibitor
            IF (ulab(i) == NAME(kk)) THEN
              inhibitaq(kk,ikin) = i
              speciesfound = .true.
              rinhibitaq(kk,ikin) = rinhibit_tmp(kk)
              IF (itot_tmp == 1) THEN
                itot_inhibitaq(kk,ikin) = 1
              END IF
            END IF
          END DO
          IF (speciesfound) THEN
            CONTINUE
          ELSE
            DO ksp = 1,nspec                     !! Search through secondary species for inhibitor
              ik = ncomp + ksp
              IF (ulab(ik) == NAME(kk)) THEN
                inhibitaq(kk,ikin) = ik
                speciesfound = .true.
                rinhibitaq(kk,ikin) = rinhibit_tmp(kk)
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
              DO k = 1,nrct                     !! Search through minerals for inhibitor
                IF (umin(k) == NAME(kk)) THEN           !! kk is the number of the inhibitor, ikin the reaction, pointing to the mineral
                  Inhibitaq(kk,ikin) = -k
                  speciesfound = .true.
                  rinhibitaq(kk,ikin) = rinhibit_tmp(kk)
                END IF
              END DO
            END IF
            IF (speciesfound) THEN
              CONTINUE
            ELSE       
              WRITE(*,*)
              WRITE(*,*) ' Cannot find species or mineral in Monod inhibition: ',dumstring(1:lss)
              dumstring = NAME(kk)
              CALL stringlen(dumstring,llen)
              WRITE(*,*) ' Looking for species: ',dumstring(1:llen)
              WRITE(*,*)
              READ(*,*)
              STOP
            END IF
          END IF
        END DO    !  End of loop through "inhibiting" species

! biomass
      ELSE IF (iaqtype(ikin) == 8) THEN       !  Monod, with thermodynamic factor, F_T
        
!  Assume only one parallel reaction for the moment
        
        ll = 1
        DO i = 1,ncomp
          itot_monodaq(i,ikin) = 0
        END DO
        READ(18,*) ratek(ll,ikin),nmonodaq(ikin),(NAME(kk),  &
            halfsat_tmp(kk),kk=1,nmonodaq(ikin))
        DO kk = 1,nmonodaq(ikin)
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
              imonodaq(kk,ikin) = i
              speciesfound = .true.
              halfsataq(kk,ikin) = halfsat_tmp(kk)
              IF (itot_tmp == 1) THEN
                itot_monodaq(kk,ikin) = 1
              END IF
            END IF
          END DO
          IF (speciesfound) THEN
            CONTINUE
          ELSE
            DO ksp = 1,nspec
              ik = ncomp + ksp
              IF (ulab(ik) == NAME(kk)) THEN
                imonodaq(kk,ikin) = ik
                speciesfound = .true.
                halfsataq(kk,ikin) = halfsat_tmp(kk)
                IF (itot_tmp == 0) THEN
                  CONTINUE
                ELSE
                  WRITE(*,*)
                  WRITE(*,*) ' Cannot specify total concentration for secondary species'
                  WRITE(*,*) '   in Monod expression: ',dumstring(1:lss)
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
        END DO    !  End of loop through "monod" species
! biomass end

! add inhibition 2011-08-15
!!  Now look for inhibition information
!        
        READ(18,*) dumlabel
        IF (dumlabel /= 'Inhibition') THEN
          WRITE(*,*) ' For Monod Biomass law, there should be an entry for inhibition terms NOW'
          WRITE(*,*)
          READ(*,*)
          STOP
        END IF
        BACKSPACE 18
        READ(18,*) dumlabel,ninhibitaq(ikin),  &
            (NAME(kk),rinhibit_tmp(kk),kk=1,ninhibitaq(ikin))
        DO kk = 1,ninhibitaq(ikin)
          dumstring = NAME(kk)
          IF (dumstring(1:4) == 'tot_') THEN
            CALL stringlen(dumstring,lss)
            NAME(kk) = dumstring(5:lss)
            itot_tmp = 1
          ELSE
            itot_tmp = 0
          END IF
          speciesfound = .false.
          DO i = 1,ncomp                         !! Search through primary species for inhibitor
            IF (ulab(i) == NAME(kk)) THEN
              inhibitaq(kk,ikin) = i
              speciesfound = .true.
              rinhibitaq(kk,ikin) = rinhibit_tmp(kk)
              IF (itot_tmp == 1) THEN
                itot_inhibitaq(kk,ikin) = 1
              END IF
            END IF
          END DO
          IF (speciesfound) THEN
            CONTINUE
          ELSE
            DO ksp = 1,nspec                     !! Search through secondary species for inhibitor
              ik = ncomp + ksp
              IF (ulab(ik) == NAME(kk)) THEN
                inhibitaq(kk,ikin) = ik
                speciesfound = .true.
                rinhibitaq(kk,ikin) = rinhibit_tmp(kk)
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
              DO k = 1,nrct                     !! Search through minerals for inhibitor
                IF (umin(k) == NAME(kk)) THEN           !! kk is the number of the inhibitor, ikin the reaction, pointing to the mineral
                  Inhibitaq(kk,ikin) = -k
                  speciesfound = .true.
                  rinhibitaq(kk,ikin) = rinhibit_tmp(kk)
                END IF
              END DO
            END IF
            IF (speciesfound) THEN
              CONTINUE
            ELSE       
              WRITE(*,*)
              WRITE(*,*) ' Cannot find species or mineral in Monod inhibition: ',dumstring(1:lss)
              dumstring = NAME(kk)
              CALL stringlen(dumstring,llen)
              WRITE(*,*) ' Looking for species: ',dumstring(1:llen)
              WRITE(*,*)
              READ(*,*)
              STOP
            END IF
          END IF
        END DO    !  End of loop through "inhibiting" species
        
        continue
!end add   
        
      ELSE
        WRITE(*,*) ' No rate formulation specified for aqueous reaction '
        READ(*,*)
        STOP
      END IF
      
    ELSE
!!      WRITE(*,*) ' No match, going to next one'
      GO TO 400  ! No match
    END IF                  ! End of a single kinetic reaction
  ELSE          ! First substring is a number, so skip it
    GO TO 400
  END IF
  
  REWIND 18
!   close(unit=18,status='keep')
  GO TO 100      ! Check for the next aqueous kinetic reaction label
  
ELSE
  GO TO 400
END IF

600 WRITE(*,*)
WRITE(*,*) ' Aqueous kinetic reaction not found in database'
WRITE(*,*) ' STOP'
WRITE(*,*)
READ(*,*)
STOP


300 DEALLOCATE(rinhibit_tmp)
DEALLOCATE(halfsat_tmp)
DEALLOCATE(name)
DEALLOCATE(skin)
DEALLOCATE(dep_tmp)

RETURN
END SUBROUTINE read_kinetics
