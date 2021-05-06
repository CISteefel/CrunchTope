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
!
! New version of read_kinetics - September 2011
! 
! Written by S. Molins based on an earlier read_kinetics version 
! and using subroutines contained in CrunchFlow as of Sept 2011.
! For these: 
!
!
! the new version of the code 
!
! 1. includes new option for line entries in input file under
!    the keyword 'AQUEOUS_KINETICS':
!
!    -pathway 'pathway' 0.5
!
! 2. reads new database format in separate database file 
!    named 'aqueous.dbs' formatted as a set of namelists
!    &Aqueous and &AqueousKinetics
!
! 3. no new global variables have been added to the code
!
! 4. backward compatibility: TBD?
!
!******************************************************************

SUBROUTINE read_kinetics_Bio(nout,ncomp,nspec,nrct,ikin,nkin,data1)
USE crunchtype
USE params
USE concentration
use mineral, only:     umin,      &
                       mukinTMP,  & 
                       keqkinTMP, &
                       p_cat_kin, &
                       ibiomass_kin, &
                       bq_kin, chi_kin, direction_kin, &
                       UseMetabolicLagAqueous,LagTimeAqueous, &
                       MetabolicLagAqueous,                   &
                       RampTimeAqueous,ThresholdConcentrationAqueous,  &
                       SubstrateForLagAqueous,  &
                       nMonodBiomassAqueous
USE strings

use io
use CrunchFunctions

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nout
INTEGER(I4B), INTENT(IN)                                    :: ncomp
INTEGER(I4B), INTENT(IN)                                    :: nspec
INTEGER(I4B), INTENT(IN)                                    :: nrct
INTEGER(I4B), INTENT(OUT)                                   :: ikin
integer(i4b),intent(in)                                     :: nkin !mineral reactions
CHARACTER (LEN=mls), INTENT(IN)                             :: data1

!  Internal variables and arrays

LOGICAL(LGT)                                                :: speciesfound
CHARACTER (LEN=mls)                                         :: dummy1
CHARACTER (LEN=mls)                                         :: dummy2
CHARACTER (LEN=mls)                                         :: dumstring
CHARACTER (LEN=mls)                                         :: dumlabel
CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE              :: name_
CHARACTER (LEN=mls)                                         :: filename2

LOGICAL(LGT)                                                :: ext

REAL(DP), DIMENSION(:), ALLOCATABLE                         :: skin
REAL(DP), DIMENSION(:), ALLOCATABLE                         :: dep_tmp
REAL(DP), DIMENSION(:), ALLOCATABLE                         :: rinhibit_tmp
REAL(DP), DIMENSION(:), ALLOCATABLE                         :: halfsat_tmp

!!CHARACTER (LEN=23)                                          :: str_endkin
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

REAL(DP), PARAMETER                                         :: tiny=1.d-10 !!1.E-12
REAL(DP)                                                    :: charge

integer(i4b),dimension(:),allocatable                         :: p_kin_cat



! namelists -----------------------------------------------------------------
!
!   type
    type reaction_t
      real(dp)                                              :: mu
      character(len=mls)                                    :: name
    end type reaction_t

    type terms_t
      character(len=mls)                                    :: name
      real(dp)                                              :: coeff
    end type terms_t

!   input data structure     
    character(len=mls)                                      :: name
    character(len=mls)                                      :: label
    character(len=mls)                                      :: type
    type(reaction_t),dimension(mcomp)                       :: stoichiometry
    real(dp)                                                :: keq

    real(dp)                                                :: rate25C
    real(dp)                                                :: activation
    type(terms_t),dimension(mcomp)                          :: monod_terms
    type(terms_t),dimension(mcomp)                          :: dependence
    type(terms_t),dimension(mcomp)                          :: inhibition


    character(len=mls)                                            :: biomass
    integer(i4b)                                                  :: chi
    real(dp)                                                      :: bq
    integer(i4b)                                                  :: direction
    logical(LGT)                                                  :: UseMetabolicLag
    REAL(DP)                                                      :: LagTime
    REAL(DP)                                                      :: RampTime
    REAL(DP)                                                      :: ThresholdConcentration
    character(len=mls)                                            :: SubstrateForLag
    character(len=mls)                                            :: planktonic

    namelist /Aqueous/                                         name,          &
                                                               label,         &
                                                               type,          &
                                                               stoichiometry, &
                                                               keq

    namelist /AqueousKinetics/                                 name,                   &
                                                               label,                  &
                                                               type,                   &
                                                               rate25C,                &
                                                               activation,             &                                                               
                                                               dependence,             &
                                                               monod_terms,            &
                                                               inhibition,             &
                                                               biomass,                &
                                                               chi,                    &
                                                               bq,                     &
                                                               direction,              &
                                                               UseMetabolicLag,        &
                                                               LagTime,                &
                                                               RampTime,               &
                                                               ThresholdConcentration, &
                                                               SubstrateForLag,        &
                                                               planktonic
!
! end namelists -------------------------------------------------------------

!
CHARACTER (LEN=mls), dimension(:), allocatable                :: name_pathway
CHARACTER (LEN=mls)                                           :: rlabel_
integer                                                       :: ios
integer                                                       :: ic, jj, kpath, &
                                                                 im, ib, ifound, &
                                                                 jkin, jreac, ksave, &
                                                                 npath,is
integer                                                       :: cat_pathway, &
                                                                 ana_pathway
integer(i4b),parameter                                        :: jmax = mls


integer(i4b),dimension(:),allocatable                         :: chi_
integer(i4b),dimension(:),allocatable                         :: direction_
!!CIS added Nov. 29, 2011
logical(lgt),dimension(:),allocatable                         :: UseMetabolicLag_
real(dp),dimension(:),allocatable                             :: LagTime_
real(dp),dimension(:),allocatable                             :: RampTime_
real(dp),dimension(:),allocatable                             :: ThresholdConcentration_
integer(i4b),dimension(:),allocatable                         :: SubstrateForLag_


integer(i4b),dimension(:),allocatable                         :: ibiomass_kin_temp
integer, DIMENSION(:), ALLOCATABLE                            :: iuser

real(dp),dimension(:),allocatable                             :: bq_
real(dp),dimension(:,:),allocatable                           :: mukin_
real(dp),dimension(:,:),allocatable                           :: mukin_cat
real(dp),dimension(:),allocatable                             :: mu_
real(dp),dimension(:),allocatable                             :: keq_
real(dp),dimension(:),allocatable                             :: keq_cat
real(dp)                                                      :: rate_0
real(dp),dimension(:),allocatable                             :: multiplier
real(dp)                                                      :: addup
!

ALLOCATE(rinhibit_tmp(mcomp+mspec))
ALLOCATE(halfsat_tmp(mcomp+mspec))
ALLOCATE(name_(mpre))
ALLOCATE(skin(50))
ALLOCATE(dep_tmp(50))

! initialize number of parallel reactions 
nreactkin(:) = 0 !allocated in FirstAllocation.F90, will be resized in reallocate.F90

! initialize temporary variables
allocate(mukin_(mpre,ncomp)); mukin_ = 0.0d0
allocate(mukin_cat(mpre,ncomp)); mukin_cat = 0.0d0
allocate(mu_(ncomp)); mu_ = 0.0d0
allocate(keq_(ncomp)); keq_ = 0.0d0
allocate(keq_cat(ncomp)); keq_cat = 0.0d0
allocate(ibiomass_kin_temp(mpre)); ibiomass_kin_temp = 0
allocate(chi_(mpre)); chi_ = 1
allocate(bq_(mpre)); bq_ = 0.0d0
allocate(direction_(mpre)); direction_ = -1
!!CIS added Nov. 29, 2011
allocate(UseMetabolicLag_(mpre)); UseMetabolicLag_ = .false.
allocate(LagTime_(mpre)); LagTime_ = 0.0d0
allocate(RampTime_(mpre)); RampTime_ = 0.0d0
allocate(ThresholdConcentration_(mpre)); ThresholdConcentration_ = 0.0d0
allocate(SubstrateForLag_(mpre)); SubstrateForLag_ = 0

allocate(name_pathway(mpre)); name_pathway(:) = ''
ALLOCATE(iuser(7)); iuser = 0
allocate(multiplier(mpre)); multiplier = 0.0d0

IF (ALLOCATED(LagTimeAqueous)) THEN
  DEALLOCATE(LagTimeAqueous)
END IF
allocate(LagTimeAqueous(mpre))  

IF (ALLOCATED(RampTimeAqueous)) THEN
  DEALLOCATE(RampTimeAqueous)
END IF       
allocate(RampTimeAqueous(mpre))  

IF (ALLOCATED(ThresholdConcentrationAqueous)) THEN
  DEALLOCATE(ThresholdConcentrationAqueous)
END IF     
allocate(ThresholdConcentrationAqueous(mpre))

IF (ALLOCATED(UseMetabolicLagAqueous)) THEN
  DEALLOCATE(UseMetabolicLagAqueous)
END IF   
allocate(UseMetabolicLagAqueous(mpre))

IF (ALLOCATED(SubstrateForLagAqueous)) THEN
  DEALLOCATE(SubstrateForLagAqueous)
END IF  
allocate(SubstrateForLagAqueous(mpre))

!!str_endkin = 'End of aqueous kinetics'

if (data1 == ' ') then

  INQUIRE(FILE='AqueousControl.ant',EXIST=ext)
  IF (EXT) THEN          !!  Aqueous Control file exists, so read input filename from it rather than prompting user
    OPEN(113,FILE='AqueousControl.ant',STATUS='old',ERR=708)
    READ(113,'(a)') filename2
    CLOSE(113,STATUS='keep')
  !! ELSE
  !! INQUIRE(FILE=trim(adjustl(data1))//'x',EXIST=ext)
  !! IF (EXT) then          !! try with the name of the database + and 'x' at the end
  !!  filename2 = trim(adjustl(data1))//'x'
  ELSE                   !!  No AqueousControl.ant file, so just use "aqueous.dbs"
    filename2 = 'aqueous.dbs'
  END IF
  !! END IF
else 
 
  filename2 = data1
  
end if 

OPEN(UNIT=112,FILE=filename2,STATUS='old')
REWIND nout
ikin = 0

100 READ(nout,'(a)',END=300) zone
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  
  ksave = 0
  DO k = 1,ikin
    IF (ssch == namkin(k)) THEN
!!  Reaction already loaded -- treat as parallel reaction
      nreactkin(k) = nreactkin(k) + 1
      !!npar = npar + 1
      ksave = k
      IF (nreactkin(k) > mreact) THEN
        WRITE(*,*)
        WRITE(*,*) ' Parameter "mreact" dimensioned too small'
        WRITE(*,*) ' Number of parallel reactions = ',nreactkin(k)
        WRITE(*,*) ' Dimension for parallel reactions = ', mreact
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
      GO TO 110
    END IF
  END DO

  lzs=ls
  CALL stringtype(ssch,lzs,res)
  dummy2(1:ls) = ssch(1:ls)
  lsave = ls
  npath = 0 ! no paths so far
  cat_pathway = 0 ! no catabolic pathway
  ana_pathway = 0 ! no anabolic pathway
  
! sergi added 0001

!  Look for various options
!  Current options include:
!        -pathway 'pathway' multiplier  name of pathways involved in building the overall
!                                       reaction stoichiometry weighted by the multiplier
!                                       if not present, the reaction name is the 'pathway'
!                                       this allows for standard entry of reactions 
!                                       but also a new entry format based on a number of 
!                                       different pathways e.g. catabolic and anabolic 
!                                       pathways in microbially mediated reactions
!        -label 'label'                 points to a database entry
!                                       if not present, label = default
!        -rate rate                     rate, if not present rate = 1.0

  110 continue
  
  rlabel_  = 'default'
  iuser(:) = 0
  
  111   id = ids + ls
  CALL sschaine(zone,id,iff,ssch,ids,ls)
  
  IF (ls == 0) THEN      !  End of line, now search database for pathways
    GO TO 200
  END IF
  
  IF (ssch == '-pathway') THEN

    npath = npath + 1

!   read label for current pathway    
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    name_pathway(npath) = ssch

!   read portion of final stoichiometry that goes to this pathway    
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    
    lzs = ls
    CALL stringtype(ssch,lzs,res)
    IF (res == 'a') THEN
      WRITE(*,*)
      WRITE(*,*) ' "-pathway pathway" should be followed by a number'
      WRITE(*,*) ' Aqueous kinetic reaction = ',dummy2(1:lsave)
      WRITE(*,*) ' String = ',ssch(1:ls)
      WRITE(*,*)
      READ(*,*)
      STOP
    ELSE
      multiplier(npath) = DNUM(ssch)
    END IF

  ELSE IF (ssch == '-label') THEN
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    rlabel_ = ssch

  ELSE IF (ssch == '-rate') THEN
    iuser(1) = 1
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    lzs = ls
    CALL stringtype(ssch,lzs,res)
    IF (res == 'a') THEN
      WRITE(*,*)
      WRITE(*,*) ' "-rate" should be followed by a number'
      WRITE(*,*) ' Aqueous kinetic reaction = ',dummy2(1:lsave)
      WRITE(*,*) ' String = ',ssch(1:ls)
      WRITE(*,*)
      READ(*,*)
      STOP
    ELSE
      rate_0 = DNUM(ssch)
    END IF
  ELSE
    CONTINUE
  END IF
  
  GO TO 111
  
! sergi added 0001  

ELSE
  GO TO 100
END IF


!================================= read database =============================================
200 continue
! done reading line in input file for this aqueous kinetic reaction

! check whether this is a parallel reaction 
! skip stoichiometry and eq. const.
if (ksave /= 0) then
  go to 333
end if

! initialize stoichiometry, and equilibrium constant
!!mukin_ = 0.d0
mu_    = 0.d0
!!keq_   = 0.d0

! first, check whether pathways were provided

if (npath == 0) then
! no specific pathway provided, reaction name is the only pathway with weight 1.0

  npath               = 1
  name_pathway(npath) = dummy2(1:lsave)
  multiplier(npath)   = 1.0d0
  
end if

addup = 0.0d0

!loop over pathways provided in input file for this reaction
do_input_pathways: do kpath=1,npath

! look for reaction in database that matches this pathway
  ifound = 0

! read all necessary pathways (reactions) from file
  do_pathways: do 

!   initialize namelist variables before reading it in from file
    name    = 'default'
    label   = 'default'
    type    = 'default'

!   stoichiometry    
    stoichiometry(:)%mu  = 0.0
    stoichiometry(:)%name = 'default'
    
!   equilibrium constant
    keq = 0.0d0

!   read 'Aqueous' namelist from file
    read(112,nml=Aqueous,iostat=ios)     

!   result from read (IOS)
    if (ios == 0) then

!     successful read, compare to current reaction
      if (name_pathway(kpath) == name) then

!       found: point to it, rewind and let me know
        ifound = kpath
        rewind(112)
!!        write(*,nml=Aqueous)
                    
        if (type == 'catabolic') then

          cat_pathway = cat_pathway + 1
            
          if (cat_pathway > 1) then 
            write(*,*)'only one catabolic pathway permitted'
            stop
          end if
            
        else if (type == 'anabolic') then
           
          ana_pathway = ana_pathway + 1
            
          if (ana_pathway > 1) then 
            write(*,*)'only one anabolic pathway permitted'
            stop
          end if

        end if          
          
        exit do_pathways

      else
      
        cycle do_pathways

      end if

    else if (ios < 0) then

!     no more pathways to read
      write(*,*)'End of file'
      exit do_pathways

    else if (ios > 0) then

!!      write(*,nml=Aqueous)
      write(*,*)' Error reading namelist: goodbye'
      stop

    end if

  end do do_pathways

!================================= check if found or not =====================================
  if_found: if (ifound == 0) then        
! input file reaction not found in database: quit

    write(*,*)'Could not find the pathway ',name_pathway(kpath), &
              ' of reaction ',dummy2(1:lsave),' in aqueous.dbs database file'
    stop

  else
! pathway found in database

!   add this reaction to the list of aqueous kinetic reaction
    if (kpath == 1) then
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
      namkin(ikin) = dummy2(1:lsave)
    end if
    
!   STOICHIOMETRY ----------------------------------------------------------------------------
!   build stoichiometry for overall reaction for mass balance 
!   for catabolic pathway, save stoichiometry for thermodynamic term

!   assemble the stoichiometric matrix in any pathway
    mu_    = 0.d0

!   go over all reactants read from namelist
    do_reactants: do jj=1,jmax

      if (stoichiometry(jj)%name == 'default') then
!       end of list, get out              
        jreac = jj - 1
        exit do_reactants
      end if

      do_components: do ic=1,ncomp

        if (stoichiometry(jj)%name == ulab(ic)) then

        ! ic = component
          mu_(ic) = stoichiometry(jj)%mu
             
          cycle do_reactants

        end if

      end do do_components
      
      if (stoichiometry(jj)%name == 'h2o' .or. &
          stoichiometry(jj)%name == 'H2O') then
        cycle do_reactants  
      end if
            
!     if here, no match has been found for a reactant
      WRITE(*,*)
      WRITE(*,*) ' Species in kinetic reaction not found in primary species list'
      WRITE(*,*) ' Reaction = ',name_pathway(kpath)
      WRITE(*,*) ' Species = ',stoichiometry(jj)%name
      WRITE(*,*)
      READ(*,*)
      STOP

    end do do_reactants 

!   Check electrical balance in reactions
    charge = 0.00
    DO i = 1,ncomp
      charge = charge + mu_(i)*chg(i)
    END DO
    IF (ABS(charge) > tiny) THEN
      WRITE(*,*)
      WRITE(*,501) 
      WRITE(*,502) charge
      WRITE(*,503) name_pathway(kpath)
      WRITE(*,*)
      STOP
    ENDIF

501 FORMAT(1X,'Aqueous kinetic reaction not electrically balanced')
502 FORMAT(1x,'Charge = ',1pe10.2)
503 FORMAT(1X,'Reaction: ',a25)

!   store the stoichiometry in the right variable
    do ic=1,ncomp

      mukin_(ikin,ic) = mukin_(ikin,ic) + multiplier(kpath) * mu_(ic)

    end do

!   equilibrium constant
    keq_(ikin)  = keq_(ikin) + multiplier(kpath) * keq ! it is actually converted in reactkin * clg         ! convert to ln

!   keep track of the multiplier
    addup = addup + multiplier(kpath)

!   store the stoichiometry in catabolic variable
    if (type == 'catabolic') then

      keq_cat(ikin) = keq !!* clg ! it is actually converted in reactkin * clg
      
      do ic=1,ncomp

        mukin_cat(ikin,ic) = mu_(ic)

      end do

    end if

  end if if_found

end do do_input_pathways

write(*,*)'multipliers for reaction ',namkin(ikin),' add up to: ',addup

IF ( addup < (1.0d0-tiny) .OR. addup > (1.0d0+tiny) ) THEN
  dumstring = namkin(ikin)
  WRITE(*,*)
  WRITE(*,*) ' Pathway fractions should sum to 1.0'
  WRITE(*,*) ' For mineral:  ', dumstring(1:20)
  WRITE(*,*)
  STOP
END IF
  

! here we have the overall stoichiometry and equilibrium constant
! and if a catabolic pathway has been entered also its stoichiometry 

 333 continue ! parallel reaction, nreactkin > 1, skipped stoichiometry step

! look for the rate expression in database file as given in input file 

! initialize ifound
ifound = 0
  
do_aqueouskinetics: do

! initialize namelist parameters to be read
  name    = 'default'
  label   = 'default'
  type    = 'default'
  rate25C = 1.0d0
  activation = 0.0d0

! monod_terms    
  monod_terms(:)%name  = 'default'
  monod_terms(:)%coeff = 0.0d0

! dependence
  dependence(:)%name  = 'default'
  dependence(:)%coeff = 0.0d0

! inhibition    
  inhibition(:)%name  = 'default'
  inhibition(:)%coeff = 0.0d0
  
! biomass related (MonodBiomass)
  biomass = 'default'
  chi = 1
  bq = 0.0d0
  direction = -1

  UseMetabolicLag= .false.
  LagTime = 0.0d0 
  RampTime = 0.0d0   
  ThresholdConcentration = 0.0d0
  SubstrateForLag = ' '
  
! read 'AqueousKinetics' namelist from file
  read(112,nml=AqueousKinetics,iostat=ios)     

! result from read (IOS)
  if (ios == 0) then

!   successful read, compare to current reaction
    if (namkin(ikin) == name .and. rlabel_ == label) then

!     found: point to it, rewind and let me know
      ifound = 1
      rewind(112)
!!      write(*,nml=AqueousKinetics)
        
      exit do_aqueouskinetics

    else
    
      cycle do_aqueouskinetics

    end if

  else if (ios < 0) then

!   no more pathways to read
    write(*,*)'end of file'
    exit do_aqueouskinetics

  else if (ios > 0) then

    write(*,*)'error reading namelist: goodbye'
    stop

  end if

end do do_aqueouskinetics
    
if_found_rate: if (ifound == 0) then
  
  write(*,*)'could not find the kinetic rate ',rlabel_, &
            ' of reaction ',namkin(ikin),' in aqueous.dbs database file'
  stop

else
! the rate expression has been found, check what has been read 

! get the appropriate index to this reaction > jkin
  if (ksave /= 0) then
    jkin = ksave ! parallel reaction
  else
    jkin = ikin  ! first time this reaction is read
    nreactkin(jkin) = 1
  end if
  
  IF (nreactkin(jkin) > mreact) THEN
    WRITE(*,*)
    WRITE(*,*) ' Redimension parameter mreact in params.F90'
    WRITE(*,*) ' Mreact =  ',mreact
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF

! Reaction rate type: TST, Monod, Irreversible, Radioactive decay, and MonodBiomass
  select case (type)

!   TST - transition state theory 'tst' type
    case ("tst")

      if (nreactkin(jkin) == 1) then
        WRITE(*,505) namkin(jkin)
        iaqtype(jkin) = 1
      end if

!   MONOD - Monod terms
    case("Monod")

      if (nreactkin(jkin) == 1) then
        WRITE(*,506) namkin(jkin)
        iaqtype(jkin) = 2
      else IF (nreactkin(jkin) > 1) THEN
          WRITE(*,*)
          WRITE(*,*) ' Only one parallel Monod reaction allowed for the moment'
          WRITE(*,*)
          READ(*,*)
          STOP
      END IF
!     save direction of reaction (-1 for case of going from left (negative stoichiometric coefficients) to right (positive stoich)
      direction_(jkin) = direction

!   IRREVERSIBLE
    case("irreversible")
    
      if (nreactkin(jkin) == 1) then
        WRITE(*,507) namkin(jkin)
        iaqtype(jkin) = 3
      end if  
        
!   RADIOACTIVE DECAY
    case("radioactive_decay")
    
      if (nreactkin(jkin) == 1) then
        WRITE(*,508) namkin(jkin)
        iaqtype(jkin) = 4
      end if

!   MONOD-BIOMASS - microbially mediated reactions
    case("MonodBiomass")

      if (nreactkin(jkin) == 1) then

        WRITE(*,506) namkin(jkin)
        IF (nreactkin(jkin) > 1) THEN
          WRITE(*,*)
          WRITE(*,*) ' Only one parallel Monod reaction allowed for the moment'
          WRITE(*,*)
          READ(*,*)
          STOP
        END IF
        iaqtype(jkin) = 8

      end if

!     save thermodynamic limitation paramters: chi, bq, and direction
      chi_(jkin)       = chi
      bq_(jkin)        = bq
      direction_(jkin) = direction
      UseMetabolicLag_(jkin) = UseMetabolicLag
      LagTime_(jkin) = LagTime
      Ramptime_(jkin) = RampTime
      ThresholdConcentration_(jkin) = ThresholdConcentration


!!    Find the substrate to trigger the metabolic lag stage (associated with ThresholdConcentrationAqueous)
      is = 0
      do_substrateAqueous: DO I = 1,ncomp
        IF (SubstrateForLag == ulab(i)) THEN
          is = i
          EXIT do_substrateAqueous
        END IF
      END DO do_substrateAqueous

      IF (is == 0 .and. SubstrateForLag /= ' ') THEN
        WRITE(*,*)
        WRITE(*,*) ' Primary species to be used a metabolic lag substrate not found: ',SubstrateForLag(1:20)
        STOP
      ELSE
        SubstrateForLag_(jkin) = is
      END IF

!     find the biomass for this reaction in the mineral list        
      ib = 0
      do_biomass: do im=1,nkin
        
        if (biomass == umin(im)) then

          ib = im
          exit do_biomass
           
        end if
        
      end do do_biomass
        
      if (ib == 0) then
        write(*,*)'biomass ',biomass,' for reaction: ',name,' is not in the list of minerals'
        stop
      end if

!     biomass pointer    
      ibiomass_kin_temp(jkin)   = ib
      
    case default

      WRITE(*,*) ' Aqueous kinetic reaction type not recognized'
      WRITE(*,*)
      READ(*,*)
      STOP

  end select

505 FORMAT(1x,'Using TST rate formulation for reaction: ',a15)
506 FORMAT(1x,'Using Monod rate formulation for reaction: ',a15)
507 FORMAT(1x,'Using irreversible rate formulation for reaction: ',a15)
508 FORMAT(1x,'Using radioactive decay formulation for reaction: ',a15)  

! DEPENDENCE TERMS
  IF (iaqtype(jkin) == 1 .OR. iaqtype(jkin) == 3 .OR. iaqtype(jkin) == 4) THEN

    !sergi    !!DO ll = 1,nreactkin(ikin)
    
    ll = nreactkin(jkin) !sergi
    
    DO i = 1,ncomp
      itot(i,ll,jkin) = 0
    END DO

   !sergi       READ(18,*) ratek(ll,ikin),nk,(NAME_(kk), dep_tmp(kk),kk=1,nk

!------------------------------------------------
!   go over all dependence terms from namelist
    do_dependence: do jj=1,jmax

      if (dependence(jj)%name == 'default') then
!       end of list, get out
        nk = jj - 1
        exit do_dependence            
      end if

      NAME_(jj)   = dependence(jj)%name
      dep_tmp(jj) = dependence(jj)%coeff

    end do do_dependence

    if (iuser(1)==1) then
      ratek(ll,jkin) = rate_0
    else 
      ratek(ll,jkin) = rate25C
    end if
!------------------------------------------------

!!          WRITE(*,*) ratek(ll,ikin)
!!          DO kk = 1,nk
!!            WRITE(*,*) dep_tmp(kk)
!!          END DO

    DO i = 1,ncomp
      dependk(i,ll,jkin) = 0.0
    END DO
    DO kk = 1,nk
      IF (kk > mpre) THEN
        WRITE(*,*)
        WRITE(*,*) ' Redimension parameter "mpre" in params.F90'
        WRITE(*,*) ' Mpre = ',mpre
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
!  First, check to see if we are dealing with a total concentration
      dumstring = NAME_(kk)
      IF (dumstring(1:4) == 'tot_') THEN
        CALL stringlen(dumstring,lss)
        NAME_(kk) = dumstring(5:lss)
        itot_tmp = 1
      ELSE
        itot_tmp = 0
      END IF
      speciesfound = .false.
      DO i = 1,ncomp                          !!  Searching through primary species list
        IF (ulab(i) == NAME_(kk)) THEN
          speciesfound = .true.
          dependk(i,ll,jkin) = dep_tmp(kk)
          IF (itot_tmp == 1) THEN
            itot(i,ll,jkin) = 1
          ELSE
            itot(i,ll,jkin) = 0
          END IF
        END IF
      END DO
      IF (speciesfound) THEN
        CONTINUE
      ELSE
        WRITE(*,*)
        WRITE(*,*) ' Kinetic aqueous reaction includes a species not in the primary species list'
        WRITE(*,*) ' Reaction label = ',namkin(jkin)
        WRITE(*,*) ' Species not found = ', NAME_(kk)
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
    END DO

      !sergi  END DO   ! End of llth parallel rxn for ikinth overall kinetic rxn
        
! MONOD TERMS
  ELSE IF (iaqtype(jkin) == 2) THEN       !  Monod
        
!  Assume only one parallel reaction for the moment
        
    ll = 1
    DO i = 1,ncomp
      itot_monodaq(i,jkin) = 0
    END DO
 !sergi   READ(18,*) ratek(ll,ikin),nmonodaq(ikin),(NAME_(kk),  &
 !sergi       halfsat_tmp(kk),kk=1,nmonodaq(ikin))

!------------------------------------------------
!   go over all monod terms from namelist
    do_monod: do jj=1,jmax

      if (monod_terms(jj)%name == 'default') then
!       end of list, get out
        nmonodaq(jkin) = jj - 1
        exit do_monod           
      end if

      NAME_(jj)       = monod_terms(jj)%name
      halfsat_tmp(jj) = monod_terms(jj)%coeff

    end do do_monod

    if (nmonodaq(jkin) == 0) then
      write(*,*)'no monod terms were entered'
      stop
    end if

    if (iuser(1)==1) then
      ratek(ll,jkin) = rate_0
    else 
      ratek(ll,jkin) = rate25C
    end if
!------------------------------------------------

    DO kk = 1,nmonodaq(jkin)
      dumstring = NAME_(kk)
      IF (dumstring(1:4) == 'tot_') THEN
        CALL stringlen(dumstring,lss)
        NAME_(kk) = dumstring(5:lss)
        itot_tmp = 1
      ELSE
        itot_tmp = 0
      END IF
      speciesfound = .false.
      DO i = 1,ncomp
        IF (ulab(i) == NAME_(kk)) THEN
          imonodaq(kk,jkin) = i
          speciesfound = .true.
          halfsataq(kk,jkin) = halfsat_tmp(kk)
          IF (itot_tmp == 1) THEN
            itot_monodaq(kk,jkin) = 1
          END IF
        END IF
      END DO
      IF (speciesfound) THEN
        CONTINUE
      ELSE
        DO ksp = 1,nspec
          ik = ncomp + ksp
          IF (ulab(ik) == NAME_(kk)) THEN
            imonodaq(kk,jkin) = ik
            speciesfound = .true.
            halfsataq(kk,jkin) = halfsat_tmp(kk)
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
          dumstring = NAME_(kk)
          CALL stringlen(dumstring,llen)
          WRITE(*,*) ' Looking for species: ',dumstring(1:llen)
          WRITE(*,*)
          READ(*,*)
          STOP
        END IF
      END IF
    END DO    !  End of loop through "monod" species
 
 
 !  Now look for inhibition information
        
!    READ(18,*) dumlabel
!    IF (dumlabel /= 'Inhibition') THEN
!      WRITE(*,*) ' For Monod law, there should be an entry for inhibition terms'
!      WRITE(*,*)
!      READ(*,*)
!      STOP
!    END IF
!    BACKSPACE 18
!    READ(18,*) dumlabel,ninhibitaq(jkin),  &
!        (NAME_(kk),rinhibit_tmp(kk),kk=1,ninhibitaq(jkin))

!------------------------------------------------
!   go over all monod terms from namelist
    do_inhib: do jj=1,jmax

      if (inhibition(jj)%name == 'default') then
!       end of list, get out
        ninhibitaq(jkin) = jj - 1
        exit do_inhib           
      end if

      NAME_(jj)        = inhibition(jj)%name
      rinhibit_tmp(jj) = inhibition(jj)%coeff

    end do do_inhib
!------------------------------------------------

    DO kk = 1,ninhibitaq(jkin)
      dumstring = NAME_(kk)
      IF (dumstring(1:4) == 'tot_') THEN
        CALL stringlen(dumstring,lss)
        NAME_(kk) = dumstring(5:lss)
        itot_tmp = 1
      ELSE
        itot_tmp = 0
      END IF
      speciesfound = .false.
      DO i = 1,ncomp                         !! Search through primary species for inhibitor
        IF (ulab(i) == NAME_(kk)) THEN
          inhibitaq(kk,jkin) = i
          speciesfound = .true.
          rinhibitaq(kk,jkin) = rinhibit_tmp(kk)
          IF (itot_tmp == 1) THEN
            itot_inhibitaq(kk,jkin) = 1
          END IF
        END IF
      END DO
      IF (speciesfound) THEN
        CONTINUE
      ELSE
        DO ksp = 1,nspec                     !! Search through secondary species for inhibitor
          ik = ncomp + ksp
          IF (ulab(ik) == NAME_(kk)) THEN
            inhibitaq(kk,jkin) = ik
            speciesfound = .true.
            rinhibitaq(kk,jkin) = rinhibit_tmp(kk)
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
            IF (umin(k) == NAME_(kk)) THEN           !! kk is the number of the inhibitor, jkin the reaction, pointing to the mineral
              Inhibitaq(kk,jkin) = -k
              speciesfound = .true.
              rinhibitaq(kk,jkin) = rinhibit_tmp(kk)
            END IF
          END DO
        END IF
        IF (speciesfound) THEN
          CONTINUE
        ELSE       
          WRITE(*,*)
          WRITE(*,*) ' Cannot find species or mineral in Monod inhibition: ',dumstring(1:lss)
          dumstring = NAME_(kk)
          CALL stringlen(dumstring,llen)
          WRITE(*,*) ' Looking for species: ',dumstring(1:llen)
          WRITE(*,*)
          READ(*,*)
          STOP
        END IF
      END IF
    END DO    !  End of loop through "inhibiting" species

! biomass
  ELSE IF (iaqtype(jkin) == 8) THEN       !  Monod, with thermodynamic factor, F_T
        
!  Assume only one parallel reaction for the moment
    
    ll = 1
    DO i = 1,ncomp
      itot_monodaq(i,jkin) = 0
    END DO
    !!READ(18,*) ratek(ll,jkin),nmonodaq(jkin),(NAME_(kk),  &
    !!    halfsat_tmp(kk),kk=1,nmonodaq(jkin))

!------------------------------------------------
!   go over all monod terms from namelist
    do_monod_biomass: do jj=1,jmax

      if (monod_terms(jj)%name == 'default') then
!       end of list, get out
        nmonodaq(jkin) = jj - 1
        exit do_monod_biomass           
      end if

      NAME_(jj)       = monod_terms(jj)%name
      halfsat_tmp(jj) = monod_terms(jj)%coeff

    end do do_monod_biomass

    if (nmonodaq(jkin) == 0) then
      write(*,*)'no monod biomass terms were entered'
      stop
    end if

    if (iuser(1)==1) then
      ratek(ll,jkin) = rate_0
    else 
      ratek(ll,jkin) = rate25C
    end if
!------------------------------------------------    
    
    DO kk = 1,nmonodaq(jkin)
      dumstring = NAME_(kk)
      IF (dumstring(1:4) == 'tot_') THEN
        CALL stringlen(dumstring,lss)
        NAME_(kk) = dumstring(5:lss)
        itot_tmp = 1
      ELSE
        itot_tmp = 0
      END IF
      speciesfound = .false.
      DO i = 1,ncomp
        IF (ulab(i) == NAME_(kk)) THEN
          imonodaq(kk,jkin) = i
          speciesfound = .true.
          halfsataq(kk,jkin) = halfsat_tmp(kk)
          IF (itot_tmp == 1) THEN
            itot_monodaq(kk,jkin) = 1
          END IF
        END IF
      END DO
      IF (speciesfound) THEN
        CONTINUE
      ELSE
        DO ksp = 1,nspec
          ik = ncomp + ksp
          IF (ulab(ik) == NAME_(kk)) THEN
            imonodaq(kk,jkin) = ik
            speciesfound = .true.
            halfsataq(kk,jkin) = halfsat_tmp(kk)
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
          dumstring = NAME_(kk)
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
!    READ(18,*) dumlabel
!    IF (dumlabel /= 'Inhibition') THEN
!      WRITE(*,*) ' For Monod Biomass law, there should be an entry for inhibition terms NOW'
!      WRITE(*,*)
!      READ(*,*)
!      STOP
!    END IF
!    BACKSPACE 18
!    READ(18,*) dumlabel,ninhibitaq(jkin),  &
!        (NAME_(kk),rinhibit_tmp(kk),kk=1,ninhibitaq(jkin))

!------------------------------------------------
!   go over all monod terms from namelist
    do_inhib_biomass: do jj=1,jmax

      if (inhibition(jj)%name == 'default') then
!       end of list, get out
        ninhibitaq(jkin) = jj - 1
        exit do_inhib_biomass           
      end if

      NAME_(jj)        = inhibition(jj)%name
      rinhibit_tmp(jj) = inhibition(jj)%coeff

    end do do_inhib_biomass
!------------------------------------------------

    DO kk = 1,ninhibitaq(jkin)
      dumstring = NAME_(kk)
      IF (dumstring(1:4) == 'tot_') THEN
        CALL stringlen(dumstring,lss)
        NAME_(kk) = dumstring(5:lss)
        itot_tmp = 1
      ELSE
        itot_tmp = 0
      END IF
      speciesfound = .false.
      DO i = 1,ncomp                         !! Search through primary species for inhibitor
        IF (ulab(i) == NAME_(kk)) THEN
          inhibitaq(kk,jkin) = i
          speciesfound = .true.
          rinhibitaq(kk,jkin) = rinhibit_tmp(kk)
          IF (itot_tmp == 1) THEN
            itot_inhibitaq(kk,jkin) = 1
          END IF
        END IF
      END DO
      IF (speciesfound) THEN
        CONTINUE
      ELSE
        DO ksp = 1,nspec                     !! Search through secondary species for inhibitor
          ik = ncomp + ksp
          IF (ulab(ik) == NAME_(kk)) THEN
            inhibitaq(kk,jkin) = ik
            speciesfound = .true.
            rinhibitaq(kk,jkin) = rinhibit_tmp(kk)
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
            IF (umin(k) == NAME_(kk)) THEN           !! kk is the number of the inhibitor, jkin the reaction, pointing to the mineral
              Inhibitaq(kk,jkin) = -k
              speciesfound = .true.
              rinhibitaq(kk,jkin) = rinhibit_tmp(kk)
            END IF
          END DO
        END IF
        IF (speciesfound) THEN
          CONTINUE
        ELSE       
          WRITE(*,*)
          WRITE(*,*) ' Cannot find species or mineral in Monod inhibition: ',dumstring(1:lss)
          dumstring = NAME_(kk)
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

! check if there are rates provided again
  IF (nreactkin(jkin) <= 0) THEN
    WRITE(*,*)
    WRITE(*,*) ' Rate law label given, but no'
    WRITE(*,*) '   rate laws included in database'
    WRITE(*,*) type
    WRITE(*,*)
    STOP
  END IF

! read next line in input file  
  GO TO 100
  
end if if_found_rate

300 continue

! everything read and found, now store in permanent variables

! allocate and intialize variables for catabolic pathway 
IF (ALLOCATED(mukinTMP)) THEN
  DEALLOCATE(mukinTMP)
END IF  
allocate(mukinTMP(ikin,ncomp))

IF (ALLOCATED(keqkinTMP)) THEN
  DEALLOCATE(keqkinTMP)
END IF  
allocate(keqkinTMP(ikin))

IF (ALLOCATED(ibiomass_kin)) THEN
  DEALLOCATE(ibiomass_kin)
END IF  
allocate(ibiomass_kin(ikin))
 
mukinTMP(:,:)   = 0.0d0
keqkin(:)       = 0.0d0
keqkinTMP(:)    = 0.0d0
ibiomass_kin(:) = 0

nMonodBiomassAqueous = 0
! loop over all reactions
do jj=1,ikin

! equilibrium constant
  keqkin(jj) = keq_(jj)

! copy stoichiometric coefficients
  do ic=1,ncomp
  
    if (dabs(mukin_(jj,ic)) > tiny) then
      mukin(jj,ic) = mukin_(jj,ic)
    end if
 
    if (iaqtype(jj) == 8) then
   
      if (dabs(mukin_(jj,ic)) > tiny) then
        mukinTMP(jj,ic) = mukin_cat(jj,ic)
      end if

    end if
    
  end do

! print in iunit2 file with the other reactions    
  if (jj == 1) then
    WRITE(iunit2,*)
    WRITE(iunit2,*) '    AQUEOUS KINETIC REACTIONS     '
    WRITE(iunit2,*)
    WRITE(iunit2,*) 'Aqueous kinetics  log K       Stoichiometric Coefficients'
    WRITE(iunit2,599) (ulab(j),j=1,ncomp)
  end if
  
  !!WRITE(iunit2,600) namkin(jj),keqkin(jj)/clg,(mukin(jj,i),i=1,ncomp)
  WRITE(iunit2,600) namkin(jj),keqkin(jj),(mukin(jj,i),i=1,ncomp)

  if (iaqtype(jj) == 8) then
  
  ! equilibrium constant
    keqkinTMP(jj) = keq_cat(jj)

    WRITE(iunit2,600) 'catabolic',keqkinTMP(jj)/clg,(mukinTMP(jj,i),i=1,ncomp)

!   copy pointer to biomass 'mineral'
    ibiomass_kin(jj) = ibiomass_kin_temp(jj)
    
!   chi, bq, and direction
    chi_kin(jj)       = chi_(jj)
    bq_kin(jj)        = bq_(jj)
    direction_kin(jj) = direction_(jj)
    UseMetabolicLagAqueous(jj) = UseMetabolicLag_(jj)
    LagTimeAqueous(jj) = LagTime_(jj)
    RampTimeAqueous(jj) = RampTime_(jj)
    ThresholdConcentrationAqueous(jj) = ThresholdConcentration_(jj)
    SubstrateForLagAqueous(jj) = SubstrateForLag_(jj)

    nMonodBiomassAqueous = nMonodBiomassAqueous + 1

  else

    direction_kin(jj) = direction_(jj)

  end if

end do

599 FORMAT(' ',27X,100A7)
600 FORMAT(' ',a12,1PE12.4,100(0PF7.2))


!deallocate temporary variables 
deallocate(mukin_)
deallocate(mukin_cat)
deallocate(mu_)
deallocate(keq_)
deallocate(keq_cat)
deallocate(ibiomass_kin_temp)
deallocate(chi_)
deallocate(bq_)
deallocate(direction_)
deallocate(name_pathway)
deallocate(multiplier)
deallocate(iuser)

DEALLOCATE(rinhibit_tmp)
DEALLOCATE(halfsat_tmp)
DEALLOCATE(name_)
DEALLOCATE(skin)
DEALLOCATE(dep_tmp)

RETURN

708 WRITE(*,*)
WRITE(*,*) ' Error opening Aqueous control file (AqueousControl.ant)'
WRITE(*,*)
READ(*,*)
STOP

END SUBROUTINE read_kinetics_Bio
