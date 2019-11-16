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

!! Written by Sergi Molins, 2014
    
    subroutine read_CatabolicPath(ncomp,nkin,ikin)

    use params
    use concentration, only: ulab, namkin, iaqtype
    use mineral, only: umin, imintype, &
                       muminTMP,   & 
                       keqminTMP,  &
                       p_cat_min,  &
                       ibiomass_min,  & 
                       bq_min, chi_min, direction_min,                &
                       UseMetabolicLagMineral,LagTimeMineral,         &
                       RampTimeMineral,ThresholdConcentrationMineral, &
                       SubstrateForLagMineral,                        &
                       nMonodBiomassMineral

    

    implicit none

!   arguments
    integer(i4b),intent(in)                                       :: ncomp
    integer(i4b),intent(in)                                       :: nkin
    integer(i4b),intent(in)                                       :: ikin

!   internal
    integer                                                       :: ios

    

          
!   type
    type reaction_t
      real(dp)                                                    :: mu
      character(len=mls)                                          :: name
    end type reaction_t

!   input data structure     
    character(len=mls)                                            :: name
    type(reaction_t),dimension(mcomp)                             :: reaction     
    real(dp)                                                      :: keq 
    character(len=mls)                                            :: biomass
    integer(i4b)                                                  :: chi
    real(dp)                                                      :: bq
    integer(i4b)                                                  :: direction
    logical(LGT)                                                  :: UseMetabolicLag
    REAL(DP)                                                      :: LagTime
    REAL(DP)                                                      :: RampTime
    REAL(DP)                                                      :: ThresholdConcentration
    character(len=mls)                                            :: SubstrateForLag

    namelist /CatabolicPathway/                                      name,                   &
                                                                     reaction,               &
                                                                     keq,                    & 
                                                                     biomass,                &
                                                                     chi,                    &
                                                                     bq,                     &
                                                                     direction,              &
                                                                     UseMetabolicLag,        &
                                                                     LagTime,                &
                                                                     RampTime,               &
                                                                     ThresholdConcentration, &
                                                                     SubstrateForLag
    
    
    integer(i4b)                                                  :: ii,jj,ic,ifound,jreac,np,is,i
    integer(i4b)                                                  :: im,ib
    integer(i4b)                                                  :: kPath,nx
    integer(i4b),parameter                                        :: jmax = mls
    integer(i4b),dimension(:),allocatable                         :: pp
    integer(i4b),dimension(:),allocatable                         :: p_min_cat

CHARACTER (LEN=mls)                                         :: filename

LOGICAL(LGT)                                                :: ext

! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

! -------------------------------------------------------------------------------------------------
!   Mineral kinetic reactions
! -------------------------------------------------------------------------------------------------
!   hardwired
    np = 1 ! number of parallel pathways
    
!   count catabolic pathways and create pointer to mineral reactions
    IF (ALLOCATED(pp)) THEN
      DEALLOCATE(pp)
    END IF  
    allocate(pp(mls)); pp(:) = 0

    IF (ALLOCATED(p_cat_min)) THEN
      DEALLOCATE(p_cat_min)
    END IF  
    allocate(p_cat_min(nkin)); p_cat_min(:) = 0
    
    nx = 0
    do_count_min: do ii=1,nkin

      if_MonodBiomass: if (imintype(np,ii) == 8) then

        nx = nx + 1
        pp(nx) = ii           !! Maps mineral number to catabolic pathway number
        p_cat_min(ii) = nx    !! Maps catabolic pathway number to mineral number

      end if if_MonodBiomass

    end do do_count_min

    if_read_catabolic_min: if (nx > 0) then

!     open file

      INQUIRE(FILE='CatabolicControl.ant',EXIST=ext)
      IF (EXT) THEN          !!  Aqueous Control file exists, so read input filename from it rather than prompting user
        OPEN(113,FILE='CatabolicControl.ant',STATUS='old',ERR=708)
        READ(113,'(a)') filename
        CLOSE(113,STATUS='keep')
      ELSE                   !!  No AqueousControl.ant file, so just use "aqueous.dbs" 
        filename = 'CatabolicPathways.in'
      END IF

      OPEN(UNIT=112,FILE=filename,STATUS='old')
 !!     open(unit=112,file='CatabolicPathways.in',status='unknown')
   
!     allocate stoichiometric coefficients, equilibrium constants, pointers
      IF (ALLOCATED(muminTMP)) THEN
        DEALLOCATE(muminTMP)
      END IF  
      allocate(muminTMP(np,nx,ncomp))           !! stoichiometric coefficient for catabolic reaction

      IF (ALLOCATED(keqminTMP)) THEN
        DEALLOCATE(keqminTMP)
      END IF  
      allocate(keqminTMP(np,nx))                !! equilibrium constant for catabolic reaction

      IF (ALLOCATED(ibiomass_min)) THEN
        DEALLOCATE(ibiomass_min)
      END IF  
      allocate(ibiomass_min(nx))                !!

      IF (ALLOCATED(p_min_cat)) THEN
        DEALLOCATE(p_min_cat)
      END IF  
      allocate(p_min_cat(nx))                   !! ??  

      IF (ALLOCATED(LagTimeMineral)) THEN
        DEALLOCATE(LagTimeMineral)
      END IF  
      allocate(LagTimeMineral(np,nx))   

      IF (ALLOCATED(RampTimeMineral)) THEN
        DEALLOCATE(RampTimeMineral)
      END IF       
      allocate(RampTimeMineral(np,nx))   

      IF (ALLOCATED(ThresholdConcentrationMineral)) THEN
        DEALLOCATE(ThresholdConcentrationMineral)
      END IF  
      allocate(ThresholdConcentrationMineral(np,nx))

      IF (ALLOCATED(UseMetabolicLagMineral)) THEN
        DEALLOCATE(UseMetabolicLagMineral)
      END IF  
      allocate(UseMetabolicLagMineral(np,nx))

      IF (ALLOCATED(SubstrateForLagMineral)) THEN
        DEALLOCATE(SubstrateForLagMineral)
      END IF  
      allocate(SubstrateForLagMineral(np,nx))

      nMonodBiomassMineral = nx

!     initialize
      muminTMP(:,:,:) = 0.0d0
      keqminTMP(:,:)  = 0.0d0
      ibiomass_min(:) = 0
      p_min_cat(:)    = 0

      LagTimeMineral(:,:) = 0.0d0
      RampTimeMineral(:,:) = 0.0d0
      ThresholdConcentrationMineral(:,:) = 0.0d0
      UseMetabolicLagMineral(:,:) = .false.

!     transfer to permanent pointer     
      p_min_cat(1:nx) = pp(1:nx)
      deallocate(pp)

!!    NOTE:  "nx" here is the number of catabolic pathways for minerals

!     loop again over catabolic pathways for minerals
      do_minerals: do kPath = 1,nx    !! Number of pathways

!       point to reaction
        ii = p_min_cat(kPath) 
        
!       look for reaction in input file and database that may require this catabolic pathway
        ifound = 0

!       read all necessary catabolic pathways from file
        do_pathways: do 

!         initialize namelist variables before reading it in from file
          name = 'default'
          reaction(:)%mu  = 0.0
          reaction(:)%name = 'default'      
          keq = 0.0d0
          biomass = 'default'
          chi = 1
          bq = 0.0d0
          direction = -1
          UseMetabolicLag= .false.
          LagTime = 0.0d0 
          RampTime = 0.0d0   
          ThresholdConcentration = 0.0d0
          SubstrateForLag = ' '

!         read 'CatabolicPathway' namelist from file
          read(112,nml=CatabolicPathway,iostat=ios)     

!         result from read (IOS)
          if (ios == 0) then

!           successful read, compare to current mineral reaction
            if (name == umin(ii)) then

!             found: point to it, rewind and let me know
              ifound = ii
              rewind(112)
!!              write(*,nml=CatabolicPathway)
              exit do_pathways

            else 

!             not found, keep trying
              cycle do_pathways             

            end if
      
          else if (ios < 0) then

!           no more pathways to read
            write(*,*)' End of file during read of catabolic pathways'
            exit do_pathways

          else if (ios > 0) then

            write(*,*)'error reading namelist for catabolic pathways: Adios'
            stop

          end if

        end do do_pathways


!       found or not found?

        if_found: if (ifound == 0) then        
!       not found

          write(*,*)'could not find the catabolic pathway for this reaction: ',umin(ii)
          stop          

        else
!       found, check and store

        ! np = parallel pathway hardwired
        ! kPath = pointer to this catabolic pathway
        
!         chi, bq, and direction
          chi_min(np,kPath)       = chi
          bq_min(np,kPath)        = bq
          direction_min(np,kPath) = direction

          UseMetabolicLagMineral(np,kPath) = UseMetabolicLag
          LagTimeMineral(np,kPath) = LagTime
          RamptimeMineral(np,kPath) = RampTime
          ThresholdConcentrationMineral(np,kPath) = ThresholdConcentration

!!  Find the substrate to trigger the metabolic lag stage (associated with ThresholdConcentrationMineral)
          is = 0
          do_substrateMineral: DO I = 1,ncomp
            IF (SubstrateForLag == ulab(i)) THEN
              is = i
              EXIT do_substrateMineral
            END IF
          END DO do_substrateMineral

          IF (is == 0 .and. SubstrateForLag /= ' ') THEN
            WRITE(*,*)
            WRITE(*,*) ' Primary species to be used a metabolic lag substrate not found: ',SubstrateForLag(1:20)
            STOP
          ELSE
            SubstrateForLagMineral(np,kPath) = is
          END IF


!         find the biomass for this reaction in the mineral list
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

!         equilibrium constant and biomass pointer          
          keqminTMP(np,kPath)   = keq * clg         ! convert to ln
          ibiomass_min(kPath)   = ib

!         go over all reactants read from namelist
          do_reactants: do jj=1,jmax

            if (reaction(jj)%name == 'default') then
!             end of list, get out              
              jreac = jj - 1
              exit do_reactants              
            end if

!! it was here
          
            do_components: do ic=1,ncomp

              if (reaction(jj)%name == ulab(ic)) then

              ! np = parallel pathway hardwired
              ! kPath = pointer to this catabolic pathway
              ! ic = component
                muminTMP(np,kPath,ic) = reaction(jj)%mu
 !! moved up too
 !!               keqminTMP(np,kPath)   = keq * clg         ! convert to ln
 !!               ibiomass_min(kPath)   = ib
                
                cycle do_reactants

              end if

            end do do_components
            
!           if here, no match has been found for a reactant
            write(*,*)'no match has been found for catabolic pathway reactant: ',reaction(jj)%name
            write(*,*)'in catabolic pathway: ',name,'. au revoir'          
            stop            

          end do do_reactants        
         
        end if if_found

      end do do_minerals
      
      close(112)
      
    else ! if_read_catabolic_min

      deallocate(pp)
      deallocate(p_cat_min)

    end if if_read_catabolic_min

!! -------------------------------------------------------------------------------------------------
!!   Aqueous kinetic reactions
!! -------------------------------------------------------------------------------------------------
!
!!   count catabolic pathways and create pointer to aqueous kinetic reactions
!    allocate(pp(mls)); pp(:) = 0 
!    allocate(p_cat_kin(ikin)); p_cat_kin(:) = 0
!    nx = 0
!    do_count_kin: do ii=1,ikin
!
!      if_MonodBiomass_: if (iaqtype(ii) == 8) then
!
!        nx = nx + 1
!        pp(nx) = ii
!        p_cat_kin(ii) = nx
!
!      end if if_MonodBiomass_
!
!    end do do_count_kin
!    
!    if_read_catabolic_kin: if (nx > 0) then
!
!!     open file
!      open(unit=112,file='CatabolicPathways.in',status='unknown')
!
!!     allocate stoichiometric coefficients, equilibrium constants, pointers
!      allocate(mukinTMP(nx,ncomp))
!      allocate(keqkinTMP(nx))
!      allocate(ibiomass_kin(nx))
!      allocate(p_kin_cat(nx))    
!
!!     initialize
!      mukinTMP(:,:)   = 0.0d0
!      keqkinTMP(:)    = 0.0d0
!      ibiomass_kin(:) = 0
!      p_kin_cat(:)    = 0
!
!!     transfer to permanent pointer     
!      p_kin_cat(1:nx) = pp(1:nx)
!      deallocate(pp)
!
!!     loop again over catabolic pathways for minerals
!
!      do_aqueous_kin: do k=1,nx
!
!!       point to reaction
!        ii = p_kin_cat(k)
!      
!!       look for reaction in input file and database that may require this catabolic pathway
!        ifound = 0
!
!!       read all necessary catabolic pathways from file
!        do_pathways_: do 
!
!!         initialize namelist variables before reading it in from file
!          name = 'default'
!          reaction(:)%mu  = 0.0
!          reaction(:)%name = 'default'      
!          keq = 0.0d0
!          biomass = 'default'
!          chi = 1
!          bq = 0.0d0
!          direction = -1
!
!!         read 'CatabolicPathway' namelist from file
!          read(112,nml=CatabolicPathway,iostat=ios)     
!
!!         result from read (IOS)
!          if (ios == 0) then
!
!!           successful read, compare to current aqueous kinetic reaction
!            if (name == namkin(ii)) then
!
!!             found: point to it, rewind and let me know
!              ifound = ii
!              rewind(112)
!              write(*,nml=CatabolicPathway)
!              exit do_pathways_
!
!            else 
!
!!             not found, keep trying
!              cycle do_pathways_
!
!            end if
!      
!          else if (ios < 0) then
!
!!           no more pathways to read
!            write(*,*)'end of file'
!            exit do_pathways_
!
!          else if (ios > 0) then
!
!            write(*,*)'error reading namelist: goodbye'
!            stop
!
!          end if
!
!        end do do_pathways_
!
!
!!       found or not found?
!
!        if_found_: if (ifound == 0) then
!!       not found
!
!          write(*,*)'could not find this reaction s catabolic pathway : ',namkin(ii)
!          stop
!
!        else
!!       found, check and store
!
!        ! k = pointer to this catabolic pathway
!        
!!         chi, bq, and direction
!          chi_kin(k)       = chi
!          bq_kin(k)        = bq
!          direction_kin(k) = direction        
!
!!         find the biomass for this reaction in the mineral list
!          ib = 0
!          do_biomass_: do im=1,nkin
!            
!            if (biomass == umin(im)) then
!                
!              ib = im
!              exit do_biomass_
!               
!            end if
!            
!          end do do_biomass_
!            
!          if (ib == 0) then
!            write(*,*)'biomass ',biomass,' for reaction: ',name,' is not in the list of minerals'
!            stop
!          end if
!
!!         equilibrium constant and biomass pointer
!          keqkinTMP(k)    = keq             ! do NOT convert to ln
!          ibiomass_kin(k) = ib
!
!!         go over all reactants read from namelist
!          do_reactants_: do jj=1,jmax
!
!            if (reaction(jj)%name == 'default') then
!!             end of list, get out              
!              jreac = jj - 1
!              exit do_reactants_              
!            end if
!
!!! biomass was here
!
!            do_components_: do ic=1,ncomp
!
!              if (reaction(jj)%name == ulab(ic)) then
!               
!               ! k = pointer to this catabolic pathway
!               ! ic = component
!                mukinTMP(k,ic)  = reaction(jj)%mu                
!!! moved up
!!!                keqkinTMP(k)    = keq             ! do NOT convert to ln
!!!                ibiomass_kin(k) = ib
!                
!                cycle do_reactants_
!
!              end if
!
!            end do do_components_
!            
!!           if here, no match has been found for a reactant
!            write(*,*)'no match has been found for catabolic pathway reactant: ',reaction(jj)%name
!            write(*,*)'in catabolic pathway: ',name,'. au revoir'
!            stop
!
!          end do do_reactants_
!
!        end if if_found_
!      
!      end do do_aqueous_kin
!  
!      close(112)
!    
!    else ! if_read_catabolic_kin
!
!      deallocate(pp)
!      deallocate(p_cat_kin)
!
!    end if if_read_catabolic_kin

    return

708 WRITE(*,*)
WRITE(*,*) ' Error opening CatabolicControl file (CatabolicControl.ant)'
WRITE(*,*)
READ(*,*)
STOP

  end subroutine read_CatabolicPath