!!! *** Copyright Notice ***
!!! �CrunchFlow�, Copyright (c) 2016, The Regents of the University of California, through Lawrence Berkeley National Laboratory 
!!! (subject to receipt of any required approvals from the U.S. Dept. of Energy).� All rights reserved.
!!!�
!!! If you have questions about your rights to use or distribute this software, please contact 
!!! Berkeley Lab's Innovation & Partnerships Office at��IPO@lbl.gov.
!!!�
!!! NOTICE.� This Software was developed under funding from the U.S. Department of Energy and the U.S. Government 
!!! consequently retains certain rights. As such, the U.S. Government has been granted for itself and others acting 
!!! on its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the Software to reproduce, distribute copies to the public, 
!!! prepare derivative works, and perform publicly and display publicly, and to permit other to do so.
!!!
!!! *** License Agreement ***
!!! �CrunchFlow�, Copyright (c) 2016, The Regents of the University of California, through Lawrence Berkeley National Laboratory)
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

REAL(DP)                                                    :: depend_tmp

LOGICAL(LGT)                                                :: SpeciesFound
LOGICAL(LGT)                                                :: MineralFound

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
integer(i4b)                                                :: ios, im, ib
integer(i4b),dimension(:),allocatable                       :: workint

namelist /BiomassDecay/                                        biomass


imonod = 0
kmonod = 0

JacobianNumerical = .FALSE.

ALLOCATE(name(mpre))
ALLOCATE(namspecies(50))
ALLOCATE(dep_tmp(mpre))
!!ALLOCATE(iuser(8))
ALLOCATE(iuser(7))
ALLOCATE(rinhibit_tmp(mcomp+mspec))
ALLOCATE(halfsat_tmp(mcomp+mspec))

IF (ALLOCATED(namAssociate)) THEN
  DEALLOCATE(namAssociate)
  ALLOCATE(namAssociate(mrct))
ELSE
  ALLOCATE(namAssociate(mrct))
END IF

IF (ALLOCATED(mintype)) THEN
  DEALLOCATE(mintype)
  ALLOCATE(mintype(mrct))
ELSE
  ALLOCATE(mintype(mrct))
END IF
! initialize to 0 (not biomass)
mintype(:)=0

!ALLOCATE(sto(50))
!ALLOCATE(ax(mreact*mrct,mreact*mrct))
!ALLOCATE(ainv(mreact*mrct,mreact*mrct))
!ALLOCATE(bx(mreact*mrct,mreact*mrct))
!ALLOCATE(yx(mreact*mrct))
!ALLOCATE(indx(mreact*mrct))

n0 = mreact*mrct

OPEN(UNIT=18,FILE=data1,STATUS='old',err=334)
REWIND nout

mnrl = 0
nkin = 0
npar = 0
nreac = 0

LocalEquilibrium = .FALSE.

100 READ(nout,'(a)',END=300) zone

iuser = 0

id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
!is        lzs=ls
!is        call convan(ssch,lzs,res)
  DO k = 1,nkin
    IF (ssch == nammin_tmp(k)) THEN
!!  Mineral already loaded -- treat as parallel reaction
      nreactmin(k) = nreactmin(k) + 1
      npar = npar + 1
      ksave = k
      IF (nreactmin(k) > mreact) THEN
        WRITE(*,*)
        WRITE(*,*) ' Parameter "mreact" dimensioned too small'
        WRITE(*,*) ' Number of parallel reactions = ',nreactmin(k)
        WRITE(*,*) ' Dimension for parallel reactions = ', mreact
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
      GO TO 110
    END IF
  END DO
  nkin = nkin + 1
  namAssociate(nkin) = ' '
  ksave = nkin
  nammin_tmp(nkin) = ssch
  IF (nkin > nm) THEN
    WRITE(*,*)
    WRITE(*,*) ' Parameter "NM" dimensioned too small'
    WRITE(*,*) ' Number of kinetic mineral reactions = ',nkin
    WRITE(*,*) ' Dimension for minerals              = ',nm
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
  nreactmin(k) = 1
  npar = 1
  lsave = ls
  
  IF (ispeciate == 1 .OR. igenericrates == 1) GO TO 555
  
!        write(*,*)
!        write(*,*) ' Nkin = ',nkin
!        write(*,*) ' Mineral = ',ssch(1:lsave)
!        write(*,*)
  
!  Look for various options
!  Current options include:
!        -label        if not present, assume default -- points to a database entry
!        -rate         log rate, in mol/m**2/s
!        -activation   temperature dependence, kcal/mole
!        -dependence   far from equilibrium dependence
!        -type         type of rate law, e.g. "tst", "monod", "irreversible"
!        -monod_terms  dependence on species and its half-saturation constant
!        -inhibition   inhibition by species and inhibition constant
!!CIS(6/08/09)    
!!       -associate    link a particular mineral/mineral rate law to the volume fraction of another
  
  110   DO i = 1,ncomp
    itot_min(i,npar,nkin) = 0
    itot_monod(i,npar,nkin) = 0
    itot_inhibit(i,npar,nkin) = 0
  END DO
  DO ik = 1,ncomp+nspec
    depend(ik,npar,nkin) = 0.0
  END DO
  
  rlabel(npar,nkin) = 'default'
  crossaff(npar,nkin) = 'none'
  kcrossaff(npar,nkin) = 0
  ssa(npar,nkin) = 0.0
  
  ndep = 0
  monod = 0
  
  111   id = ids + ls
  CALL sschaine(zone,id,iff,ssch,ids,ls)
  
  IF (ls == 0) THEN      !  End of line, now search database for mineral
    GO TO 425
  END IF
  
  IF (ssch == '-label') THEN
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    rlabel(npar,nkin) = ssch
  ELSE IF (ssch == '-local_equilibrium' .OR. ssch == '-equilibrium' .OR. ssch == '-localequilibrium') THEN
    LocalEquilibrium(nkin) = .TRUE.
  
  ELSE IF (ssch == '-SSA' .OR. ssch == '-ssa') THEN
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    lzs = ls
    CALL stringtype(ssch,lzs,res)
    IF (res == 'a') THEN
      WRITE(*,*)
      WRITE(*,*) ' "-SSA" should be followed by a number'
      WRITE(*,*) ' Mineral = ',namdum(1:lsave)
      WRITE(*,*) ' String = ',ssch(1:ls)
      WRITE(*,*)
      READ(*,*)
      STOP
    ELSE
      ssa(npar,nkin) = DNUM(ssch)
    END IF
    
  ELSE IF (ssch == '-cross_affinity') THEN
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    lzs = ls
    crossaff(npar,nkin) = ssch
  ELSE IF (ssch == '-rate') THEN
    iuser(1) = 1
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    lzs = ls
    CALL stringtype(ssch,lzs,res)
    IF (res == 'a') THEN
      WRITE(*,*)
      WRITE(*,*) ' "-rate" should be followed by a number'
      WRITE(*,*) ' Mineral = ',namdum(1:lsave)
      WRITE(*,*) ' String = ',ssch(1:ls)
      WRITE(*,*)
      READ(*,*)
      STOP
    ELSE
      rate0(npar,nkin) = DNUM(ssch)
      continue
    END IF
!!  ELSE IF (ssch == '-sigma') THEN
!!    iuser(8) = 1
!!    id = ids + ls
!!    CALL sschaine(zone,id,iff,ssch,ids,ls)
!!    lzs = ls
!!    CALL stringtype(ssch,lzs,res)
!!    IF (res == 'a') THEN
!!      WRITE(*,*)
!!      WRITE(*,*) ' "-sigma" should be followed by a number'
!!      WRITE(*,*) ' Mineral = ',namdum(1:lsave)
!!      WRITE(*,*) ' String = ',ssch(1:ls)
!!      WRITE(*,*)
!!      READ(*,*)
!!      STOP
!!    ELSE
!!      sigma(nkin) = DNUM(ssch)
!!      IF (sigma(nkin) /= 0.0) THEN
!!        CrystalSizeDistribution(nkin) = .TRUE.
!!      END IF
!!    END IF

  ELSE IF (ssch == '-activation') THEN
    iuser(2) = 1
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
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
      ea(npar,nkin) = DNUM(ssch)
    END IF
  ELSE IF (ssch == '-dependence') THEN
    iuser(3) = 1
    ndep = ndep + 1
    namdep_nyf(ndep,npar,nkin) = ' '
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    ls_spec = ls
    lzs=ls
!          call convan(ssch,lzs,res)
    CALL stringtype(ssch,lzs,res)
    nam_depend = ssch
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    lzs=ls
    CALL stringtype(ssch,lzs,res)
    IF (res == 'a') THEN
      WRITE(*,*)
      WRITE(*,*) ' Looking for an exponent following species'
      WRITE(*,*) ' Dependence on species ',nam_depend(1:ls_spec)
      WRITE(*,*) ' Mineral = ',namdum(1:lsave)
      WRITE(*,*)
      READ(*,*)
      STOP
    ELSE
      depend_tmp = DNUM(ssch)
    END IF
!  First, check to see if we are dealing with a total concentration
    dumstring = nam_depend
    IF (dumstring(1:4) == 'tot_') THEN
      CALL stringlen(dumstring,lss)
      nam_depend = dumstring(5:lss)
      itot_tmp = 1
    ELSE
      itot_tmp = 0
    END IF
    speciesfound = .false.
    DO i = 1,ncomp
      IF (ulab(i) == nam_depend) THEN
        idepend(ndep,npar,nkin) = i
        speciesfound = .true.
        namdep_nyf(ndep,npar,nkin) = 'found'
        depend(ndep,npar,nkin) = depend_tmp
        IF (itot_tmp == 1) THEN
          itot_min(ndep,npar,nkin) = 1
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
      
      namdep_nyf(ndep,npar,nkin) = nam_depend
      depend(ndep,npar,nkin) = depend_tmp
      
      
!            do ksp = 1,nspec
!              ik = ncomp + ksp
!              if (ulab(ik).eq.nam_depend) then
!                idepend(ndep,npar,nkin) = ik
!                speciesfound = .true.
!                depend(ndep,npar,nkin) = depend_tmp
!                if (itot_tmp.eq.0) then
!                  continue
!                else
!                  write(*,*)
!                  write(*,*) ' Cannot specify total concentration for secondary species'
!                  write(*,*) ' in rate law: ',dumstring(1:lss)
!                  write(*,*)
!                  stop
!                endif
!              endif
!            end do
!            if (speciesfound) then
!              continue
!            else
!              write(*,*)
!              write(*,*) ' Cannot find species in rate law: ',dumstring(1:lss)
!              write(*,*)
!              stop
!            endif
   
    END IF
  ELSE IF (ssch == '-type') THEN
    iuser(4) = 1
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    lzs = ls
    CALL convan(ssch,lzs,res)
    IF (ssch == 'tst') THEN
      imintype(npar,nkin) = 1
    ELSE IF (ssch == 'monod') THEN
      imintype(npar,nkin) = 2
    ELSE IF (ssch == 'nucleation') THEN
      imintype(npar,nkin) = 10
    ELSE IF (ssch == 'irreversible') THEN
      imintype(npar,nkin) = 3
    ELSE IF (ssch == 'ripening') THEN
      imintype(npar,nkin) = 11
    ELSE
      imintype(npar,nkin) = 1
    END IF
  ELSE IF (ssch == '-associate') THEN
    iuser(7) = 1
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    ls_mineral = ls
    lzs=ls
    CALL stringtype(ssch,lzs,res)
    namAssociate(nkin) = ssch
  ELSE IF (ssch == '-monod_term' .OR. ssch == '-monod_terms') THEN
    iuser(5) = 1
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    ls_monod = ls
    lzs=ls
!          call convan(ssch,lzs,res)
    CALL stringtype(ssch,lzs,res)
    nam_monod = ssch
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    CALL stringtype(ssch,lzs,res)
    IF (res == 'a') THEN
      WRITE(*,*)
      WRITE(*,*) ' Looking for half-saturation constant for Monod term'
      WRITE(*,*) ' Monod species',nam_monod(1:ls_monod)
      WRITE(*,*) ' Mineral = ',namdum(1:lsave)
      WRITE(*,*)
      READ(*,*)
      STOP
    ELSE
      halfsat_tmp(1) = DNUM(ssch)
    END IF
!  First, check to see if we are dealing with a total concentration
    dumstring = nam_monod
    IF (dumstring(1:4) == 'tot_') THEN
      CALL stringlen(dumstring,lss)
      nam_monod = dumstring(5:lss)
      itot_tmp = 1
    ELSE
      itot_tmp = 0
    END IF
    speciesfound = .false.
    DO i = 1,ncomp
      IF (ulab(i) == nam_monod) THEN
        speciesfound = .true.
        halfsat(i,npar,nkin) = halfsat_tmp(1)
        IF (itot_tmp == 1) THEN
          itot_min(i,npar,nkin) = 1
        END IF
      END IF
    END DO
    IF (speciesfound) THEN
      CONTINUE
    ELSE
      DO ksp = 1,nspec
        ik = ncomp + ksp
        IF (ulab(ik) == nam_monod) THEN
          speciesfound = .true.
          halfsat(ik,npar,nkin) = halfsat_tmp(1)
          IF (itot_tmp == 0) THEN
            CONTINUE
          ELSE
            WRITE(*,*)
            WRITE(*,*) ' Cannot specify total concentration for secondary species'
            WRITE(*,*) ' in Monod rate law: ',dumstring(1:lss)
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
        WRITE(*,*) ' Cannot find species in Monod rate law: ',dumstring(1:lss)
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
    END IF
    
  ELSE IF (ssch == '-inhibition') THEN
    iuser(6) = 1
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    ls_inhibit = ls
    lzs=ls
!          call convan(ssch,lzs,res)
    CALL stringtype(ssch,lzs,res)
    nam_inhibit = ssch
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    CALL stringtype(ssch,lzs,res)
    IF (res == 'a') THEN
      WRITE(*,*)
      WRITE(*,*) ' Looking for inhibition constant for Monod term'
      WRITE(*,*) ' Inhibiting species',nam_inhibit(1:ls_inhibit)
      WRITE(*,*) ' Mineral = ',namdum(1:lsave)
      WRITE(*,*)
      READ(*,*)
      STOP
    ELSE
      rinhibit_tmp(1) = DNUM(ssch)
    END IF
!  First, check to see if we are dealing with a total concentration
    dumstring = nam_inhibit
    IF (dumstring(1:4) == 'tot_') THEN
      CALL stringlen(dumstring,lss)
      nam_inhibit = dumstring(5:lss)
      itot_tmp = 1
    ELSE
      itot_tmp = 0
    END IF
    speciesfound = .false.
    DO i = 1,ncomp
      IF (ulab(i) == nam_inhibit) THEN
        speciesfound = .true.
        rinhibit(i,npar,nkin) = rinhibit_tmp(1)
        IF (itot_tmp == 1) THEN
          itot_min(i,npar,nkin) = 1
        END IF
      END IF
    END DO
    IF (speciesfound) THEN
      CONTINUE
    ELSE
      DO ksp = 1,nspec
        ik = ncomp + ksp
        IF (ulab(ik) == nam_inhibit) THEN
          speciesfound = .true.
          rinhibit(ik,npar,nkin) = rinhibit_tmp(1)
          IF (itot_tmp == 0) THEN
            CONTINUE
          ELSE
            WRITE(*,*)
            WRITE(*,*) ' Cannot specify total concentration for secondary species'
            WRITE(*,*) ' in Monod inhibition term: ',dumstring(1:lss)
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
        WRITE(*,*) ' Cannot find inhibiting species in Monod rate law: ',dumstring(1:lss)
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
    END IF
  ELSE IF (ssch == '-biomass') THEN
! label this mineral as biomass
    id = ids + ls
    lzs=ls
    WRITE(*,*) ' Mineral = ',namdum(1:lsave),' is biomass.'
    WRITE(*,*) ' In the database, specify number of cells per mol of biomass in the molar volume field'
    mintype(nkin) = 1
  ELSE
    CONTINUE
  END IF
  
  GO TO 111
  
!  NOTE: "nkin" is number of reacting minerals (not number of rate laws)
  
ELSE
  GO TO 100
END IF

!    *********SEARCH FOR KINETIC MINERAL REACTION IN DATABASE**********

425 np = npar

ndepend(np,nkin) = ndep
!      nmonod(np,nkin) = ndep


555 mnrl = mnrl + 1
namrl(mnrl) = nammin_tmp(nkin)
!!WRITE(*,*) nammin_tmp(nkin)

IF (ispeciate == 1 .OR. igenericrates == 1) GO TO 100

REWIND 18

!          First, find end of the minerals section

200 READ(18,'(a)') dummy1
IF (dummy1 == 'Begin mineral kinetics') THEN
  READ(18,'(a)') dummy1
  GO TO 475
ELSE
  GO TO 200
END IF

475 READ(18,'(a)',END=300) dummy1

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
  READ(18,'(a)',END=300) dummy1
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
    
    READ(18,'(a)',END=300) dummy1
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
    
    READ(18,'(a)',END=300) dummy1
    IF (iuser(1) == 0) THEN ! No user input, so read
      id = 1
      iff = mls
      CALL sschaine(dummy1,id,iff,ssch,ids,ls)
      lzs=ls
      CALL convan(ssch,lzs,res)
      IF (ssch == 'rate(25c)') THEN                              !!  Read RATE
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
    
    READ(18,'(a)',END=300) dummy1
    IF (iuser(2) == 0) THEN ! No user input, so read
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
      
      READ(18,'(a)',END=300) dummy1
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
      
      READ(18,'(a)',END=300) dummy1
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
      
      READ(18,'(a)',END=300) dummy1
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
!!        READ(18,*) chi(np,nkin),bq(np,nkin),direction(np,nkin)

      END IF        ! If input by user, skip block above
! biomass end

! sergi: biomass decay
  ELSE IF (imintype(np,nkin) == 9) THEN
    
!     first time around allocation
      if (.not.allocated(biomass_decay)) then

        allocate(biomass_decay(mreact,mrct))

      end if 

!     read 'biomass' namelist from file
      read(18,nml=BiomassDecay,iostat=ios)     

!     result from read (IOS)
      if (ios == 0) then

!       successful read

!       find the biomass for this reaction in the mineral list
        ib = 0
        do_biomass: do im=1,nkin-1
            
          if (biomass == namrl(im)) then
               
            ib = im
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
      
      READ(18,'(a)',END=300) dummy1
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
!!                READ(18,*) Kformation(np,nkin),HyperbolicInhibitionName(np,nkin),HyperbolicInhibitionDepend(np,nkin)

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
      
      READ(18,'(a)',END=300) dummy1
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
    
!          read(18,'(a)',end=300) dummy1
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

300 CLOSE(UNIT=18,STATUS='keep')

!      do k = 1, nreac
!        do l = 1, nreac
!          wx(l,k) = ax(l,k)
!        end do
!      end do

! call ludcmp(ax,nreac,n0,indx,dd)

!      write(*,*) ' nreac = ',nreac
!      pause

!      do i = 1, nreac
!        do k = 1, nreac
!     yx(k) = 0.
!      if (k .eq. i) then
!            yx(k) = 1.0
!          end if
!        end do

!        call lubksb(ax,nreac,n0,indx,yx)

!        do l = 1, nreac
!          ainv(l,i) = yx(l)
!        end do
!      end do

!      nr = 0
!      do k = 1, nkin
!        do np = 1,nreactmin(k)
!      do k = 1,nreac
!          nr = nr + 1
!     do i = 1, ncomp
!       sum = 0.0
!       do kk = 1, nreac
!         sum = sum + ainv(k,kk)*bx(kk,i)
!c            end do
!       mumin(1,k,i) = sum
!          end do

!  Note:  "Ainv" is the transformation matrix for log Ks and coefficients also

!        end do
!      end do


!      do k = 1,nkin
!        do np = 1,nreactmin(k)
!         do k = 1,nreac
!          write(*,2001) (mumin(1,k,i),i=1,ncomp)
!          write(*,*)
!        end do
!      end do
!      pause

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


!DEALLOCATE(sto)
!DEALLOCATE(ax)
!DEALLOCATE(ainv)
!DEALLOCATE(bx)
!DEALLOCATE(yx)
!DEALLOCATE(indx)

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

END SUBROUTINE read_minkin
