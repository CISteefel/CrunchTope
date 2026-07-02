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

integer :: no_rows
integer :: irow

!!!  ******************************************************************

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

n0 = mreact*mrct

!!!REWIND nout

mnrl = 0
nkin = 0
npar = 0
nreac = 0
LocalEquilibrium = .FALSE.

!!! NOTE: "nout" is unit number for input file
!!! NOTE: "iunitDatabase" is the unit number for the database file

100 READ(nout,'(a)',END=300) zone

iuser = 0

id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN

  DO k = 1,nkin

    IF (ssch == nammin_tmp(k)) THEN
!!    Mineral already loaded -- treat as parallel reaction
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
      
      CYCLE       !! Exit loop so that a new mineral is not added (only a parallel reaction is added)
      
    ELSE
      
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
      
    END IF
  END DO
  
  IF (ispeciate == 1 .OR. igenericrates == 1) GO TO 555
  
!  Look for various options
!  Current options include:
!        -label        if not present, assume default -- points to a database entry
!        -rate         log rate, in mol/m**2/s
!        -activation   temperature dependence, kcal/mole
!        -dependence   far from equilibrium dependence
!        -type         type of rate law, e.g. "tst", "monod", "irreversible"
!        -monod_terms  dependence on species and its half-saturation constant
!        -inhibition   inhibition by species and inhibition constant  
!!       -associate    link a particular mineral/mineral rate law to the volume fraction of another
  
  DO i = 1,ncomp
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
  
  111 id = ids + ls      !! Check for more rate law options on same line
  CALL sschaine(zone,id,iff,ssch,ids,ls)
  
  IF (ls == 0) THEN      !  End of line, now search database for mineral
    
    CONTINUE    !! Going to line 425 to start database search
    
  ELSE
   
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
    
!     First, check to see if we are dealing with a total concentration
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
      
!       If species is not among primary species, search later from list of secondary aqueous,
!       exchange, and surface complex species (exchange and surface complexes not yet read)   
        namdep_nyf(ndep,npar,nkin) = nam_depend
        depend(ndep,npar,nkin) = depend_tmp
   
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
!     First, check to see if we are dealing with a total concentration
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
  !   label this mineral as biomass
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
    GO TO 100      !!! Go to next "zone" read
  END IF
    
END IF             
  

!    *********SEARCH FOR KINETIC MINERAL REACTION IN DATABASE**********

425 np = npar

ndepend(np,nkin) = ndep

555 mnrl = mnrl + 1
namrl(mnrl) = nammin_tmp(nkin)
!!WRITE(*,*) nammin_tmp(nkin)

IF (ispeciate == 1 .OR. igenericrates == 1) GO TO 100
!!!IF (ispeciate /= 1 .AND. igenericrates /= 1) THEN

!   First, find end of the minerals section

iunitDatabase = 8

inquire(unit=iunitDatabase, opened=itsopen) 
IF ( itsopen ) THEN
  ClOSE(UNIT=iunitDatabase,status='KEEP')
END IF
OPEN(UNIT=iunitDatabase,FILE=data1,STATUS='old',err=334)

type :: RowData
  integer, allocatable :: cols(:)
end type

type(RowData), allocatable :: rows(:)

integer :: no_rows
integer :: i

open(unit=11, file='text.txt')

no_rows = count_lines(11)
allocate(rows(no_rows))

do i=1,no_rows
  rows(i)%cols = read_row(11)
enddo










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
