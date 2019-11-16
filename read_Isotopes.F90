!!! *** Copyright Notice ***
!!! “CrunchFlow”, Copyright (c) 2016, The Regents of the University of California, through Lawrence Berkeley National Laboratory 
!!! (subject to receipt of any required approvals from the U.S. Dept. of Energy).  All rights reserved.
!!! 
!!! If you have questions about your rights to use or distribute this software, please contact 
!!! Berkeley Lab's Innovation & Partnerships Office at  IPO@lbl.gov.
!!! 
!!! NOTICE.  This Software was developed under funding from the U.S. Department of Energy and the U.S. Government 
!!! consequently retains certain rights. As such, the U.S. Government has been granted for itself and others acting 
!!! on its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the Software to reproduce, distribute copies to the public, 
!!! prepare derivative works, and perform publicly and display publicly, and to permit other to do so.
!!!
!!! *** License Agreement ***
!!! “CrunchFlow”, Copyright (c) 2016, The Regents of the University of California, through Lawrence Berkeley National Laboratory)
!!! subject to receipt of any required approvals from the U.S. Dept. of Energy).  All rights reserved."
!!! 
!!! Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
!!! 
!!! (1) Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
!!!
!!! (2) Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer 
!!! in the documentation and/or other materials provided with the distribution.
!!!
!!! (3) Neither the name of the University of California, Lawrence Berkeley National Laboratory, U.S. Dept. of Energy nor the names of 
!!! its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
!!!
!!! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, 
!!! BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT 
!!! SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
!!! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
!!! OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
!!! LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF 
!!! THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!!!
!!! You are under no obligation whatsoever to provide any bug fixes, patches, or upgrades to the features, functionality or 
!!! performance of the source code ("Enhancements") to anyone; however, if you choose to make your
!!! Enhancements available either publicly, or directly to Lawrence Berkeley National Laboratory, without 
!!! imposing a separate written license agreement for such 
!!! Enhancements, then you hereby grant the following license: a  non-exclusive, royalty-free perpetual license to install, use, 
!!! modify, prepare derivative works, incorporate into other computer software, distribute, and sublicense such enhancements or 
!!! derivative works thereof, in binary and source code form.

!!!      ****************************************
    
SUBROUTINE read_Isotopes(nout,ncomp,nspec,nrct,nsurf,nexchange,ngas)
USE crunchtype
USE CrunchFunctions
USE params
USE concentration
USE mineral
USE strings
USE isotope

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                     :: nout
INTEGER(I4B), INTENT(IN)                                     :: ncomp
INTEGER(I4B), INTENT(IN)                                     :: nspec
INTEGER(I4B), INTENT(IN)                                     :: nrct
INTEGER(I4B), INTENT(IN)                                     :: nsurf
INTEGER(I4B), INTENT(IN)                                     :: nexchange
INTEGER(I4B), INTENT(IN)                                     :: ngas
!!LOGICAL(LGT), INTENT(IN OUT)                                 :: speciesfound

!  Internal variables and arrays

CHARACTER (LEN=mls)                                          :: tempstring
!!!LOGICAL(LGT)                                                 :: constraingas

INTEGER(I4B)                                                 :: id
INTEGER(I4B)                                                 :: iff
INTEGER(I4B)                                                 :: ids
INTEGER(I4B)                                                 :: ls
INTEGER(I4B)                                                 :: lzs
INTEGER(I4B)                                                 :: nlen1
INTEGER(I4B)                                                 :: k
INTEGER(I4B)                                                 :: ids_save
INTEGER(I4B)                                                 :: ll
INTEGER(I4B)                                                 :: idsave
INTEGER(I4B)                                                 :: ls_save
INTEGER(I4B)                                                 :: lcond

CHARACTER (LEN=mls)                                          :: dumstring
CHARACTER (LEN=mls)                                          :: stringspecies

INTEGER(I4B)                                                 :: lenstring
INTEGER(I4B)                                                 :: iPrimaryRare
INTEGER(I4B)                                                 :: lenPrimaryRare
CHARACTER (LEN=mls)                                          :: namePrimaryRare
INTEGER(I4B)                                                 :: iPrimaryCommon
INTEGER(I4B)                                                 :: lenPrimaryCommon
CHARACTER (LEN=mls)                                          :: namePrimaryCommon

INTEGER(I4B)                                                 :: kMineralRare
INTEGER(I4B)                                                 :: lenMineralRare
CHARACTER (LEN=mls)                                          :: nameMineralRare
INTEGER(I4B)                                                 :: kMineralCommon
INTEGER(I4B)                                                 :: lenMineralCommon
CHARACTER (LEN=mls)                                          :: nameMineralCommon

INTEGER(I4B)                                                 :: nnIsotope
INTEGER(I4B)                                                 :: isotopologue
INTEGER(I4B)                                                 :: kIsotopologue
INTEGER(I4B)                                                 :: IsotopologueRareSave
INTEGER(I4B)                                                 :: IsotopologueCommonSave
INTEGER(I4B)                                                 :: i

LOGICAL(LGT)                                                 :: RareIsotopeFound
LOGICAL(LGT)                                                 :: CommonIsotopeFound

IF (ALLOCATED(isotopeRare)) THEN
  DEALLOCATE(isotopeRare)
END IF
ALLOCATE(isotopeRare(25))
IF (ALLOCATED(isotopeCommon)) THEN
  DEALLOCATE(isotopeCommon)
END IF
ALLOCATE(isotopeCommon(25))
IF (ALLOCATED(nameIsotopeRare)) THEN
  DEALLOCATE(nameIsotopeRare)
END IF
ALLOCATE(nameIsotopeRare(25))
IF (ALLOCATED(nameIsotopeCommon)) THEN
  DEALLOCATE(nameIsotopeCommon)
END IF
ALLOCATE(nameIsotopeCommon(25))
IF (ALLOCATED(IsotopeReference)) THEN
  DEALLOCATE(IsotopeReference)
END IF
ALLOCATE(IsotopeReference(25))

IF (ALLOCATED(MoleFractionAqueousRare)) THEN
  DEALLOCATE(MoleFractionAqueousRare)
END IF
ALLOCATE(MoleFractionAqueousRare(25))
IF (ALLOCATED(MoleFractionAqueousCommon)) THEN
  DEALLOCATE(MoleFractionAqueousCommon)
END IF
ALLOCATE(MoleFractionAqueousCommon(25))

IF (ALLOCATED(iPointerIsotope)) THEN
  DEALLOCATE(iPointerIsotope)
END IF
ALLOCATE(iPointerIsotope(ncomp))

!!!IF (ALLOCATED(IsotopePrimaryRare)) THEN
!!!  DEALLOCATE(IsotopePrimaryRare)
!!!END IF
!!!ALLOCATE(IsotopePrimaryRare(ncomp))
!!!IF (ALLOCATED(IsotopePrimaryCommon)) THEN
!!!  DEALLOCATE(IsotopePrimaryCommon)
!!!END IF
!!!ALLOCATE(IsotopePrimaryCommon(ncomp))

IF (ALLOCATED(kIsotopeRare)) THEN
  DEALLOCATE(kIsotopeRare)
END IF
ALLOCATE(kIsotopeRare(25))

IF (ALLOCATED(kIsotopeCommon)) THEN
  DEALLOCATE(kIsotopeCommon)
END IF
ALLOCATE(kIsotopeCommon(25))

IF (ALLOCATED(nameIsotopeMineralRare)) THEN
  DEALLOCATE(nameIsotopeMineralRare)
END IF
ALLOCATE(nameIsotopeMineralRare(25))

IF (ALLOCATED(nameIsotopeMineralCommon)) THEN
  DEALLOCATE(nameIsotopeMineralCommon)
END IF
ALLOCATE(nameIsotopeMineralCommon(25))

IF (ALLOCATED(isotopeBackReactionOption)) THEN
  DEALLOCATE(isotopeBackReactionOption)
END IF
ALLOCATE(isotopeBackReactionOption(25))

IF (ALLOCATED(PointerToPrimaryIsotope)) THEN
  DEALLOCATE(PointerToPrimaryIsotope)
END IF
ALLOCATE(PointerToPrimaryIsotope(25))

IF (ALLOCATED(UseAqueousMoleFraction)) THEN
  DEALLOCATE(UseAqueousMoleFraction)
END IF
ALLOCATE(UseAqueousMoleFraction(25))

IF (ALLOCATED(kPointerIsotope)) THEN
  DEALLOCATE(kPointerIsotope)
END IF
ALLOCATE(kPointerIsotope(nrct))

IF (ALLOCATED(IsotopeMineralRare)) THEN
  DEALLOCATE(IsotopeMineralRare)
END IF
ALLOCATE(IsotopeMineralRare(nrct))
IF (ALLOCATED(IsotopeMineralCommon)) THEN
  DEALLOCATE(IsotopeMineralCommon)
END IF
ALLOCATE(IsotopeMineralCommon(nrct))

IF (ALLOCATED(MoleFractionMineralRare)) THEN
  DEALLOCATE(MoleFractionMineralRare)
END IF
ALLOCATE(MoleFractionMineralRare(25))
IF (ALLOCATED(MoleFractionMineralCommon)) THEN
  DEALLOCATE(MoleFractionMineralCommon)
END IF
ALLOCATE(MoleFractionMineralCommon(25))

IF (ALLOCATED(dMoleFractionAqueousCommon)) THEN
  DEALLOCATE(dMoleFractionAqueousCommon)
END IF
ALLOCATE(dMoleFractionAqueousCommon(ncomp,25))
IF (ALLOCATED(dMoleFractionAqueousRare)) THEN
  DEALLOCATE(dMoleFractionAqueousRare)
END IF
ALLOCATE(dMoleFractionAqueousRare(ncomp,25))

IF (ALLOCATED(lambda)) THEN
  DEALLOCATE(lambda)
END IF
ALLOCATE(lambda(nrct))
lambda = 0.0d0

IF (ALLOCATED(decay)) THEN
  DEALLOCATE(decay)
END IF
ALLOCATE(decay(nrct))

iPointerIsotope = 0
IsotopePrimaryRare   = .FALSE.
IsotopePrimaryCommon = .FALSE.

IsotopeMineralRare   = .FALSE.
IsotopeMineralCommon = .FALSE.

nIsotopePrimary = 0
nIsotopeMineral = 0

kPointerIsotope = 0

!!REWIND nout

100 READ(nout,'(a)',END=300) zone
nlen1 = LEN(zone)
!      call majuscules(zone,nlen1)
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
!        call convan(ssch,lzs,res)
  CALL stringtype(ssch,lzs,res)
  IF (res /= 'a') THEN
    WRITE(*,*)
    WRITE(*,*) ' Isotope input should start with an ASCII string'
    WRITE(*,*) '   String found', ssch(1:30)
    WRITE(*,*) '       ABORTING RUN  '
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
END IF

IF (ssch == 'primary') THEN
    
!!  If "primary", then read isotope pair (rare then common) followed by reference standard 
  nIsotopePrimary = nIsotopePrimary + 1

  id = ids + ls
  CALL sschaine(zone,id,iff,ssch,ids,ls)
  IF(ls /= 0) THEN
    lzs=ls
    CALL stringtype(ssch,lzs,res)
    
    IF (res /= 'a') THEN
      WRITE(*,*)
      WRITE(*,*) ' Isotope input following "primary" or "mineral" should start with an ASCII string'
      WRITE(*,*) '   String found', ssch(1:30)
      WRITE(*,*) '       ABORTING RUN  '
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
    
    namePrimaryRare = ssch
    CALL stringlen(namePrimaryRare,lenPrimaryRare)
!!    write(*,*) ssch(1:lenPrimaryRare)
    iprimaryRare = GetPrimarySpeciesNumber (ncomp,namePrimaryRare)
    isotopeRare(nIsotopePrimary) = iPrimaryRare
!!  Array is length ncomp:  iPointerIsotope
    iPointerIsotope(iPrimaryRare) = nIsotopePrimary
    nameIsotopeRare(nIsotopePrimary) = namePrimaryRare
    
  ELSE
    WRITE(*,*)
    WRITE(*,*) ' No input following rare isotope in Isotope input'
    WRITE(*,*) '       ABORTING RUN  '
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
    
  id = ids + ls
  CALL sschaine(zone,id,iff,ssch,ids,ls)
  IF(ls /= 0) THEN
    lzs=ls
    CALL stringtype(ssch,lzs,res)
      
    IF (res /= 'a') THEN
      WRITE(*,*)
      WRITE(*,*) ' Isotope input following "primary" or "mineral" should start with an ASCII string'
      WRITE(*,*) '   String found', ssch(1:30)
      WRITE(*,*) '       ABORTING RUN  '
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
    
    namePrimaryCommon = ssch
    CALL stringlen(namePrimaryCommon,lenPrimaryCommon)
    iPrimaryCommon = GetPrimarySpeciesNumber (ncomp,namePrimaryCommon)
    isotopeCommon(nIsotopePrimary) = iPrimaryCommon
!!  Array is length ncomp:  iPointerIsotope
    iPointerIsotope(iPrimaryCommon) = nIsotopePrimary
    nameIsotopeCommon(nIsotopePrimary) = namePrimaryCommon
      
  ELSE
    WRITE(*,*)
    WRITE(*,*) ' No input following "primary" in Isotope input'
    WRITE(*,*) '       ABORTING RUN  '
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
      
  id = ids + ls
  CALL sschaine(zone,id,iff,ssch,ids,ls)
  IF(ls /= 0) THEN
    lzs=ls
    CALL stringtype(ssch,lzs,res)
      
    IF (res /= 'n') THEN
      WRITE(*,*)
      WRITE(*,*) ' Isotope pair should be followed by the reference standard value (rare/common)' 
      WRITE(*,*) '   String found', ssch(1:30)
      WRITE(*,*) '       ABORTING RUN  '
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
        
    IsotopeReference(nIsotopePrimary) = DNUM(ssch)
      
  ELSE
    WRITE(*,*)
    WRITE(*,*) ' No input following rare isotope in Isotope input'
    WRITE(*,*) '       ABORTING RUN  '
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
      
      
ELSE IF (ssch == 'mineral') THEN

!!  If "mineral", then read mineral pair (rare then common) followed back reaction approach (bulk, surface, none)
  nIsotopeMineral = nIsotopeMineral + 1
  
  id = ids + ls
  CALL sschaine(zone,id,iff,ssch,ids,ls)
  IF(ls /= 0) THEN
    lzs=ls
    CALL stringtype(ssch,lzs,res)
      
    IF (res /= 'a') THEN
      WRITE(*,*)
      WRITE(*,*) ' Isotope input following "primary" or "mineral" should start with an ASCII string'
      WRITE(*,*) '   String found', ssch(1:30)
      WRITE(*,*) '       ABORTING RUN  '
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
    
    nameMineralRare = ssch
    CALL stringlen(nameMineralRare,lenMineralRare)
    CALL GetMineralNumber (nrct,nameMineralRare,kMineralRare)
    kIsotopeRare(nIsotopeMineral) = kMineralRare
!!  Array is length nrct:  kPointerIsotope
    kPointerIsotope(kMineralRare) = nIsotopeMineral
    nameIsotopeMineralRare(nIsotopeMineral) = nameMineralRare
      
  ELSE
    WRITE(*,*)
    WRITE(*,*) ' No input following "mineral" in Isotope input'
    WRITE(*,*) '       ABORTING RUN  '
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
  
  id = ids + ls
  CALL sschaine(zone,id,iff,ssch,ids,ls)
  IF(ls /= 0) THEN
    lzs=ls
    CALL stringtype(ssch,lzs,res)
      
    IF (res /= 'a') THEN
      WRITE(*,*)
      WRITE(*,*) ' Isotope input following rare isotope mineral should start with an ASCII string'
      WRITE(*,*) '   String found', ssch(1:30)
      WRITE(*,*) '       ABORTING RUN  '
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
    
    nameMineralCommon = ssch
    CALL stringlen(nameMineralCommon,lenMineralCommon)
    CALL GetMineralNumber (nrct,nameMineralCommon,kMineralCommon)
    kIsotopeCommon(nIsotopeMineral) = kMineralCommon
!!  Array is length nrct:  kPointerIsotope
    kPointerIsotope(kMineralCommon) = nIsotopeMineral
    nameIsotopeMineralCommon(nIsotopeMineral) = nameMineralCommon
      
  ELSE
    WRITE(*,*)
    WRITE(*,*) ' No input following rare isotope mineral in Isotope input'
    WRITE(*,*) '       ABORTING RUN  '
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF

  id = ids + ls
  CALL sschaine(zone,id,iff,ssch,ids,ls)
  IF(ls /= 0) THEN
    lzs=ls
    CALL convan(ssch,lzs,res)
    CALL stringtype(ssch,lzs,res)
      
    IF (res /= 'a') THEN
      WRITE(*,*)
      WRITE(*,*) ' Following mineral isotope pair should be option for calculating back reaction'
      WRITE(*,*) ' Options are:  Bulk  Surface  None'
      WRITE(*,*) '   String found: ', ssch(1:30)
      WRITE(*,*) '       ABORTING RUN  '
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF

    isotopeBackReactionOption(nIsotopeMineral) = ssch
    dumstring = ssch
    CALL stringlen(dumstring,lenstring)
    IF (dumstring(1:lenstring) == 'bulk' .OR. dumstring(1:lenstring) == 'surface' .OR. dumstring(1:lenstring) == 'none') THEN 
       CONTINUE
    ELSE
       WRITE(*,*)  
       WRITE(*,*) ' Isotope back reaction option for minerals not recognized'
       WRITE(*,*) ' Options are:  Bulk  Surface  None'
       WRITE(*,*) '   String found: ', ssch(1:30)
       WRITE(*,*) '       ABORTING RUN  '
       WRITE(*,*)
       READ(*,*)
       STOP 
    END IF 
    
  id = ids + ls
  CALL sschaine(zone,id,iff,ssch,ids,ls)
  IF(ls /= 0) THEN
    lzs=ls
    CALL convan(ssch,lzs,res)
    CALL stringtype(ssch,lzs,res)
    
    IF (res == 'a') THEN
      WRITE(*,*)
      WRITE(*,*) ' Following mineral isotope option, reading number for decay constant lambda' 
      WRITE(*,*) ' Should be a number, not a character string'
      WRITE(*,*) '   String found: ', ssch(1:30)
      WRITE(*,*) '       ABORTING RUN  '
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
    
    lambda(nIsotopeMineral) = DNUM(ssch)
    
    END IF
  
  ELSE

    WRITE(*,*)
    WRITE(*,*) ' No input following isotope mineral pair'
    WRITE(*,*) '       ABORTING RUN  '
    WRITE(*,*)
    READ(*,*)
    STOP

  END IF  
  

    
ELSE

!!!  WRITE(*,*)
!!!  WRITE(*,*) ' Isotope input should start with either "primary" or "mineral" '
!!!  WRITE(*,*) '   Input found', ssch(1:30)
!!!  WRITE(*,*) '   ABORTING RUN  '
!!!  WRITE(*,*)
!!!  READ(*,*)
!!!  STOP
    
END IF

GO TO 100

300 CONTINUE



DO nnIsotope = 1,nIsotopePrimary
  IsotopePrimaryRare(isotopeRare(nnIsotope))     = .TRUE.
  IsotopePrimaryCommon(isotopeCommon(nnIsotope)) = .TRUE.
END DO

!! Find which primary isotope system goes with the mineral isotope system
DO kIsotopologue = 1,nIsotopeMineral

  CommonIsotopeFound = .FALSE.
  RareIsotopeFound   = .FALSE.

!! Look in MineralRare for Rare Primary
  kMineralRare     = kIsotopeRare(kIsotopologue)
  kMineralCommon   = kIsotopeCommon(kIsotopologue)

  DO isotopologue = 1,nIsotopePrimary
    iPrimaryRare   = isotopeRare(isotopologue)
    iPrimaryCommon = isotopeCommon(isotopologue)

    DO i = 1,ncomp
      IF (mumin(1,kMineralRare,i) /= 0.0) THEN
        IF ( ulab(iPrimaryRare) == ulab(i) ) THEN
          RareIsotopeFound = .TRUE.
          isotopologueRareSave = isotopologue
        END IF
      END IF
      IF (mumin(1,kMineralCommon,i) /= 0.0) THEN
        IF ( ulab(iPrimaryCommon) == ulab(i) ) THEN
          CommonIsotopeFound = .TRUE.
          isotopologueCommonSave = isotopologue
        END IF
      END IF
    END DO

    IF (CommonIsotopeFound .AND. RareIsotopeFound) THEN
!!  Check that we got the same answer from Rare and Common
      IF (isotopologueRareSave /= isotopologueCommonSave) THEN
         WRITE(*,*)
         WRITE(*,*) ' Rare and common isotope pointer should be to the same system'
         WRITE(*,*)
         READ(*,*)
         STOP
      END IF
      PointerToPrimaryIsotope(kIsotopologue) = isotopologueRareSave
    ELSE
      CONTINUE
    END IF
  END DO

  IF (CommonIsotopeFound .AND. RareIsotopeFound) THEN
    CONTINUE
  ELSE
    WRITE(*,*)
    WRITE(*,*) ' Primary isotope system not found in mineral'
    WRITE(*,*) ' Looking for rare isotope system in mineral:   ', umin(kMineralRare)
    WRITE(*,*) ' Looking for common isotope system in mineral: ', umin(kMineralCommon)
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF

  IsotopeMineralRare(kIsotopeRare(kIsotopologue))     = .TRUE.
  IsotopeMineralCommon(kIsotopeCommon(kIsotopologue)) = .TRUE.
END DO


RETURN
END SUBROUTINE read_Isotopes

