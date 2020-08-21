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

SUBROUTINE find_condition(nin,nout,found,phfound,ncomp,  &
    nspec,nrct,nkin,ngas,nexchange,nsurf,ndecay,         & 
    ph,guessph,constraint,nchem,unitsflag,jpor,          &
    DensityModule,RunningPest)
USE crunchtype
USE params
USE concentration
USE mineral
USE medium
USE temperature
USE Runtime, ONLY: SkipAdjust

IMPLICIT NONE

!  ****************  Beginning of interface blocks  ***********************************

INTERFACE
  SUBROUTINE rocalc(tc,rotemp,nchem,DensityModule,unitsflag)
  USE crunchtype
  USE params
  IMPLICIT NONE
  REAL(DP), INTENT(IN)                                  :: tc
  REAL(DP), INTENT(OUT)                                 :: rotemp
  INTEGER(I4B), INTENT(IN)                              :: nchem
  INTEGER(I4B), DIMENSION(:), INTENT(IN)                :: unitsflag
  CHARACTER (LEN=mls)                                   :: DensityModule
  END SUBROUTINE rocalc
END INTERFACE

INTERFACE
  SUBROUTINE read_ph(nout,ph,guessph,i,isolution,constraint,nrct,phfound)
  USE crunchtype
  USE params
  IMPLICIT NONE
  INTEGER(I4B), INTENT(IN)                                    :: nout
  REAL(DP), DIMENSION(:), INTENT(IN OUT)                      :: ph
  REAL(DP), DIMENSION(:), INTENT(IN OUT)                      :: guessph
  INTEGER(I4B), INTENT(IN)                                    :: i
  INTEGER(I4B), INTENT(IN)                                    :: isolution
  CHARACTER (LEN=mls), DIMENSION(:,:), INTENT(IN OUT)         :: constraint
  INTEGER(I4B), INTENT(IN)                                    :: nrct
  LOGICAL(LGT), INTENT(IN OUT)                                :: phfound   
  END SUBROUTINE read_ph
END INTERFACE

INTERFACE
  SUBROUTINE read_concentration(nout,i,isolution,constraint,  &
    ncomp,nspec,nrct,ngas,speciesfound)
  USE crunchtype
  USE params
  IMPLICIT NONE
  INTEGER(I4B), INTENT(IN)                                     :: nout
  INTEGER(I4B), INTENT(IN)                                     :: i
  INTEGER(I4B), INTENT(IN)                                     :: isolution
  CHARACTER (LEN=mls), DIMENSION(:,:), INTENT(IN OUT)          :: constraint
  INTEGER(I4B), INTENT(IN)                                     :: ncomp
  INTEGER(I4B), INTENT(IN)                                     :: nspec
  INTEGER(I4B), INTENT(IN)                                     :: nrct
  INTEGER(I4B), INTENT(IN)                                     :: ngas
  LOGICAL(LGT), INTENT(IN OUT)                                 :: speciesfound
  END SUBROUTINE read_concentration
END INTERFACE

INTERFACE
  SUBROUTINE units_concentration(nout,unitsflag,nchem,labeltemp,lcond)
  USE crunchtype
  USE params
  IMPLICIT NONE
  INTEGER(I4B), INTENT(IN)                                        :: nout
  INTEGER(I4B), DIMENSION(:), INTENT(INOUT)                       :: unitsflag
  INTEGER(I4B), INTENT(IN)                                        :: nchem
  CHARACTER (LEN=mls), INTENT(IN)                                 :: labeltemp
  INTEGER(I4B), INTENT(IN)                                        :: lcond 
  END SUBROUTINE units_concentration
END INTERFACE

!  ****************  End of interface blocks  ***********************************

!  External variables

INTEGER(I4B), INTENT(IN)                                     :: nin
INTEGER(I4B), INTENT(IN)                                     :: nout
INTEGER(I4B), INTENT(IN)                                     :: ncomp
INTEGER(I4B), INTENT(IN)                                     :: nspec
INTEGER(I4B), INTENT(IN)                                     :: nrct
INTEGER(I4B), INTENT(IN)                                     :: nkin
INTEGER(I4B), INTENT(IN)                                     :: ngas
INTEGER(I4B), INTENT(IN)                                     :: nexchange
INTEGER(I4B), INTENT(IN)                                     :: nsurf
INTEGER(I4B), INTENT(IN)                                     :: ndecay
INTEGER(I4B), INTENT(OUT)                                    :: nchem
INTEGER(I4B), INTENT(IN)                                     :: jpor

INTEGER(I4B), DIMENSION(:), INTENT(INOUT)                    :: unitsflag

LOGICAL(LGT), INTENT(IN OUT)                                 :: found
LOGICAL(LGT), INTENT(IN OUT)                                 :: phfound

REAL(DP), DIMENSION(:), INTENT(OUT)                          :: ph
REAL(DP), DIMENSION(:), INTENT(OUT)                          :: guessph

CHARACTER (LEN=mls), DIMENSION(:,:), INTENT(IN OUT)          :: constraint
CHARACTER (LEN=mls),  INTENT(IN)                             :: DensityModule    
LOGICAL(LGT), INTENT(IN)                                     :: RunningPest 

!  *****  Internal variables  *****************************************

LOGICAL(LGT)                                                 :: endoffile
LOGICAL(LGT)                                                 :: speciesfound
LOGICAL(LGT)                                                 :: mineralfound
LOGICAL(LGT)                                                 :: tempfound
LOGICAL(LGT)                                                 :: PorFound
LOGICAL(LGT)                                                 :: SaturationFound
LOGICAL(LGT)                                                 :: equilfound
LOGICAL(LGT)                                                 :: complexfound
LOGICAL(LGT)                                                 :: topefound
LOGICAL(LGT)                                                 :: SkipSaturationAdjust

INTEGER(I4B)                                                 :: ncount
INTEGER(I4B)                                                 :: i
INTEGER(I4B)                                                 :: iunivar
INTEGER(I4B)                                                 :: ix
INTEGER(I4B)                                                 :: k
INTEGER(I4B)                                                 :: ks
INTEGER(I4B)                                                 :: id
INTEGER(I4B)                                                 :: idecayall
INTEGER(I4B)                                                 :: kd
INTEGER(I4B)                                                 :: lcond
INTEGER(I4B)                                                 :: ls

REAL(DP)                                                     :: tempc
REAL(DP)                                                     :: RoTemp
REAL(DP)                                                     :: permole
REAL(DP)                                                     :: PorTemp
REAL(DP)                                                     :: SaturationTemp
REAL(DP)                                                     :: sum
REAL(DP)                                                     :: MeanSaltConcentration
REAL(DP)                                                     :: RoSolution
REAL(DP)                                                     :: WritePorosity
REAL(DP)                                                     :: MineralMolality
REAL(DP)                                                     :: TotalVolumeMinerals

CHARACTER (LEN=mls)                                          :: dumstring
CHARACTER (LEN=mls)                                          :: labeltemp

CHARACTER (LEN=mls)                                           :: parchar
CHARACTER (LEN=mls)                                           :: parfind

INTEGER(I4B)                                                  :: lchar

REWIND nin

nchem = 0


endoffile = .false.
200 nchem = nchem + 1

CALL read_condition(nin,nout,found,ncount,nchem,endoffile)

IF (endoffile) THEN
  nchem = nchem - 1
  RETURN
END IF

IF (found) THEN
!  WRITE(*,*)
!  WRITE(*,5022) nchem
!  WRITE(*,*) condlabel(nchem),condtitle(nchem)
!  WRITE(*,*)

  labeltemp = condlabel(nchem)
  CALL stringlen(labeltemp,lcond)

  IF (nchem > mchem) THEN
    WRITE(*,*)
    WRITE(*,*) ' Too many geochemical conditions'
    WRITE(*,*) ' Dimension mchem larger in params.f90'
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
  
!!  Look for logical telling the code to skip the mineral saturation adjustment

  parchar = 'SkipSaturationAdjust'
  parfind = ' '
  SkipSaturationAdjust = .FALSE.
  CALL read_logical(nout,lchar,parchar,parfind,SkipSaturationAdjust)
  
  SkipAdjust(nchem) = SkipSaturationAdjust

! Search for a temperature for this geochemical condition
  
  tempfound = .false.
  CALL read_temp(nout,nchem,tempc,tempfound,labeltemp,lcond)
  IF (tempfound) THEN
    tempcond(nchem) = tempc
    IF (tempc /= tinit .AND. RunIsothermal) THEN
      WRITE(*,*) 
      WRITE(*,*) ' Isothermal run specified, but temperature points for conditions not all the same'
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
  ELSE
    WRITE(*,*)
    WRITE(*,*) ' No temperature given for condition:  ', labeltemp(1:lcond)
    WRITE(*,*) ' Using value given in TEMPERATURE keyword block: ',tinit
    WRITE(*,*)
    tempcond(nchem) = tinit
  END IF

  tempc = tempcond(nchem) 

! Search for a local porosity for this geochemical condition (only if jpor = -1)

  IF (jpor == -1 .OR. jpor == 2) THEN    !  This is a global parameter setting a constant porosity--can be overridden in a condition
    porfound = .false.
    CALL read_porosity(nout,nchem,portemp,porfound)
    IF (porfound) THEN
      porcond(nchem) = portemp
    ELSE
      IF (.NOT. RunningPest) THEN
        WRITE(*,*)
        WRITE(*,*) ' No porosity given for condition:  ', labeltemp(1:lcond)
        WRITE(*,*) ' Using value given in POROSITY keyword block: ',constantpor
        WRITE(*,*)
      END IF
      porcond(nchem) = constantpor
    END IF
  ELSE
    CONTINUE    !  Calculate local porosity below based on (1 - sum mineral volume fractions)
  END IF 

!!  Search for liquid saturation within geochemical condition

  SaturationFound = .FALSE. 
  CALL ReadSaturation(nout,nchem,SaturationTemp,SaturationFound)
  IF (SaturationFound) THEN
    SaturationCond(nchem) = SaturationTemp
    IF (SaturationTemp < 1.0d0) THEN
      isaturate = 1
    END IF
  ELSE
    IF (.NOT. RunningPest) THEN
      WRITE(*,*)
      WRITE(*,*) ' No saturation given for condition:  ', labeltemp(1:lcond)
      WRITE(*,*) ' Using value given in FixSaturation keyword block: ', FixSaturation
      WRITE(*,*)
    END IF
    SaturationCond(nchem) = FixSaturation
  END IF

!  Look for concentration units

   CALL units_concentration(nout,unitsflag,nchem,labeltemp,lcond)
  
  DO i = 1,ncomp     ! Cycle through primary species list
    iunivar = 0
    phfound = .false.
    speciesfound = .false.
    IF (ulab(i) == 'h+' .OR. ulab(i) == 'H+') THEN
      guess(i,nchem) = 0.0
      ph(nchem) = -500.0
      guessph(nchem) = -500.0
      constraint(i,nchem) = ' '
!  Check to see if pH is given, rather than H+
      CALL read_ph(nout,ph,guessph,i,nchem,constraint,nrct,phfound)
      IF (phfound) THEN                !              itype(i,nchem) = 7
        IF (ph(nchem) /= -500.0) THEN
!          WRITE(*,*)
!          WRITE(*,*) ' Specified pH = ',ph(nchem)
        END IF
        IF (constraint(i,nchem) /= ' ') THEN
!          WRITE(*,*)
!          WRITE(*,*) ' Constraint mineral for pH provided',  &
!              constraint(i,nchem)
!          WRITE(*,*)
        END IF
        IF (ph(nchem) == -500.0 .AND. guessph(nchem) == -500.0) THEN
          guess(i,nchem) = 1.e-07
        ELSE IF (guessph(nchem) == -500.0 .AND.  &
              ph(nchem) /= -500.0) THEN  ! no guess needed
          guess(i,nchem) = 10**(-ph(nchem))
        ELSE IF (guessph(nchem) /= -500.0 .AND.  &
              ph(nchem) == -500.0) THEN  ! case of mineral constraint
          guess(i,nchem) = 10**(-guessph(nchem))
        ELSE IF (guessph(nchem) /= -500.0 .AND.  &
              ph(nchem) /= -500.0) THEN         ! redundant case
          guessph(nchem) = ph(nchem)
          guess(i,nchem) = 10**(-ph(nchem))
        END IF
        
      ELSE
        
!  Case where concentration of H+ is given
        
        CALL read_concentration(nout,i,nchem,constraint,  &
            ncomp,nspec,nrct,ngas,speciesfound)

        IF (unitsflag(nchem) == 2) THEN                   !  Units in PPM
          IF (wtaq(i) /= 0.0) THEN
            ctot(i,nchem) = ctot(i,nchem)*0.001/wtaq(i)
            guess(i,nchem) = guess(i,nchem)*0.001/wtaq(i)
          ELSE
            dumstring = ulab(i)
            CALL stringlen(dumstring,ls)
            WRITE(*,*)
            WRITE(*,*) ' Database has molecular weight for species = 0'
            WRITE(*,*) ' Cant convert from PPM'
            WRITE(*,*) ' Offending species: ', dumstring(1:ls)
            WRITE(*,*) ' In condition: ', labeltemp(1:lcond)
            WRITE(*,*) 
            READ(*,*)
            STOP
          END IF
        ELSE IF (unitsflag(nchem) == 3) THEN              !  Units in mmol/kg
          ctot(i,nchem) = ctot(i,nchem)*0.001
          guess(i,nchem) = guess(i,nchem)*0.001
        ELSE IF (unitsflag(nchem) == 4) THEN              !  Units in umol/kg
          ctot(i,nchem) = ctot(i,nchem)*1.E-06
          guess(i,nchem) = guess(i,nchem)*1.E-06
        ELSE                                              !  Units already in mol/kg
          CONTINUE
        END IF
        
        IF (guess(i,nchem) < 1.e-15) THEN
          IF (ctot(i,nchem) > 0.0) THEN
            guess(i,nchem) = 1.e-05
          ELSE
            guess(i,nchem) = 1.e-08
          END IF
        END IF
        
      END IF
    ELSE                ! For other concentrations
      
      guess(i,nchem) = 1.E-06
      
!     write(*,*)
!     write(*,*) ' Starting to read concentrations'
!     write(*,*) ' Looking for ',ulab(i)
      
      CALL read_concentration(nout,i,nchem,constraint,  &
          ncomp,nspec,nrct,ngas,speciesfound)

      IF (unitsflag(nchem) == 2) THEN                   !  Units in PPM
        IF (wtaq(i) /= 0.0) THEN
          ctot(i,nchem) = ctot(i,nchem)*0.001/wtaq(i)
          guess(i,nchem) = guess(i,nchem)*0.001/wtaq(i)
        ELSE
          dumstring = ulab(i)
          CALL stringlen(dumstring,ls)
          WRITE(*,*)
          WRITE(*,*) ' Database has molecular weight for species = 0'
          WRITE(*,*) ' Cant convert from PPM'
          WRITE(*,*) ' Offending species: ', dumstring(1:ls)
          WRITE(*,*) ' In condition: ', labeltemp(1:lcond)
          WRITE(*,*) 
          READ(*,*)
          STOP
        END IF
      ELSE IF (unitsflag(nchem) == 3) THEN              !  Units in mmol/kg
        ctot(i,nchem) = ctot(i,nchem)*0.001
        guess(i,nchem) = guess(i,nchem)*0.001
      ELSE IF (unitsflag(nchem) == 4) THEN              !  Units in umol/kg
        ctot(i,nchem) = ctot(i,nchem)*1.E-06
        guess(i,nchem) = guess(i,nchem)*1.E-06

      ELSE                                              !  Units already in mol/kg
        CONTINUE
      END IF
      
!     write(*,*) ctot(i,nchem),guess(i,nchem)
!            write(*,*) constraint(i,nchem)
!            write(*,*) ' Gas pp = ',gaspp(i,nchem)
!     write(*,*) ' Itype = ',itype(i,nchem)
      IF (itype(i,nchem) == 4) THEN
        ctot(i,nchem) = gaspp(i,nchem)
      END IF
      IF (guess(i,nchem) == 0.0) THEN
        IF (itype(i,nchem) == 1) THEN
          guess(i,nchem) = ctot(i,nchem)
        END IF
        IF (itype(i,nchem) == 7) THEN
          guess(i,nchem) = ctot(i,nchem)
        END IF
        IF (itype(i,nchem) == 8) THEN
          guess(i,nchem) = ctot(i,nchem)
        END IF
        IF (itype(i,nchem) == 4) THEN
          IF (ulab(i) == 'o2(aq)' .OR. ulab(i) == 'O2(aq)') THEN
            guess(i,nchem) = 0.0013*ctot(i,nchem)
          ELSE
            guess(i,nchem) = 1.e-10
          END IF
        END IF
        IF (itype(i,nchem) == 3) THEN
          guess(i,nchem) = 1.e-15
        END IF
        IF (itype(i,nchem) == 2) THEN
          guess(i,nchem) = 1.e-06
        END IF
      END IF
      
    END IF
    IF (speciesfound .OR. phfound) THEN
      CONTINUE
    ELSE
      WRITE(*,*)
      WRITE(*,*) ' Primary species not found'
      WRITE(*,*) ' Species = ',ulab(i)
      WRITE(*,*) ' In geochemical condition: ',labeltemp(1:lcond)
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
    
  END DO
  
  IF (DensityModule /= 'temperature') THEN
!   Check that total concentration constraint has been used (otherwise ctot has not been calculated)
    IF (itype(MeanSalt(1),nchem) /= 1 .OR. itype(MeanSalt(2),nchem) /= 1) THEN
      CALL stringlen(DensityModule,ls)
      WRITE(*,*) 
      WRITE(*,*) ' Density calculation depending on concentration requires that relevant '
	  WRITE(*,*) '     species be specified with a total concentration constraint'
      WRITE(*,*) ' Density module used: ', DensityModule(1:ls)
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
  END IF

  tempc = tempcond(nchem) 
  CALL rocalc(tempc,rotemp,nchem,DensityModule,unitsflag)
  rocond(nchem) = rotemp

  IF (unitsflag(nchem) == 5) THEN              !  Units in molarity 

    IF (DensityModule /= 'temperature') THEN
!     Calculate the correction for the mass fraction of water:  kg_solution/kg_water
      MeanSaltConcentration = 0.001*(wtaq(MeanSalt(1))*ctot(MeanSalt(1),nchem) +     &
      wtaq(MeanSalt(2))*ctot(MeanSalt(2),nchem))                                       !! Mean salt in kg salt/L solution (g/cm**3)
      RoSolution = rocond(nchem)/1000.0d0                                              !! Convert fluid density from kg/m**3 to kg/L (or g/cm**3), same as salt
      OneOverMassFraction(nchem) = RoSolution/(RoSolution - MeanSaltConcentration)
      conversion(nchem) = OneOverMassFraction(nchem)/RoSolution
    ELSE
      OneOverMassFraction(nchem) = 1.0d0
      conversion(nchem) = 1000.0d0/rocond(nchem)
    END IF

    DO i = 1,ncomp
      ctot(i,nchem) = ctot(i,nchem)*conversion(nchem)               
      guess(i,nchem) = guess(i,nchem)*conversion(nchem)
    END DO

!    write(*,*) 
!    write(*,*) ' Conversion factor from molarity to molality: ', conversion
!    read(*,*)
  ELSE 
    conversion(nchem) = 1.0d0
    OneOverMassFraction(nchem) = 1.0d0
  END IF

!  Change "unitsflag" back to molality now and recompute the density
  IF (unitsflag(nchem) == 5) THEN
    unitsflag(nchem) = 1
  END IF
  
!  CALL rocalc(tempc,rotemp,nchem,DensityModule,unitsflag)
!  rocond(nchem) = rotemp
 
!  Look for mineral volume fractions and surface areas
  
  DO k = 1,nkin
    dumstring = umin(k)
    CALL stringlen(dumstring,ls)
    mineralfound = .false.
!   write(*,*)
!   write(*,*) ' Starting to read minerals'
!   write(*,*) ' Looking for ',umin(k)
    CALL read_mineral(nout,k,nchem, ncomp,nspec,nrct,ngas,mineralfound)
!          write(*,*) umin(k),volin(k,nchem),areain(k,nchem)
    IF (mineralfound) THEN
      
!!    Specific surface area specified
      IF (iarea(k,nchem) == 1) THEN    !!  Specific surface area specified
        IF (volmol(k) == 0.0D0 .OR. volmol(k) == 500.0d0) THEN
          WRITE(*,*) 
          WRITE(*,*) ' Zero molar volume or "500.00" for mineral: ', dumstring(1:ls)
          WRITE(*,*) ' Cannot calculate bulk surface area (divide by 0)'
          WRITE(*,*) ' In geochemical condition: ', labeltemp(1:lcond)
          WRITE(*,*)
          READ(*,*)
          STOP
        END IF
        IF (volin(k,nchem) > 0.0d0) THEN
          areain(k,nchem) = volin(k,nchem)*specific(k,nchem)*wtmin(k)/volmol(k)
        ELSE
          areain(k,nchem) = voltemp(k,nchem)*specific(k,nchem)*wtmin(k)/volmol(k)
        END IF
        
!!    Bulk surface area specified      
      ELSE                             !!  Bulk surface area specified
        IF (volin(k,nchem) /= 0.0) THEN
          specific(k,nchem) = areain(k,nchem)*volmol(k)/(volin(k,nchem)*wtmin(k))
        ELSE
          specific(k,nchem) = 0.0
!!          WRITE(*,*)
!!          WRITE(*,*) ' Use specific surface area option with a zero initial volume fraction'
!!          WRITE(*,*) ' For mineral: ',dumstring(1:ls)
!!          WRITE(*,*) ' In geochemical condition: ', labeltemp(1:lcond)
!!          WRITE(*,*)
!!          STOP
        END IF
        
      END IF
      
    ELSE
      IF (.NOT. RunningPest) THEN
        WRITE(*,*)
        WRITE(*,*) ' WARNING:  Mineral not found '
        WRITE(*,*) umin(k)
        WRITE(*,*) ' In geochemical condition: ', labeltemp(1:lcond)
        WRITE(*,*) ' Assuming no mineral initially present and'
        WRITE(*,*) ' surface area = 100 m^2/m^3'
        WRITE(*,*)
      END IF
!            if (ispeciate.eq.0 .and. igenericrates.eq.0) then
!              pause
!            endif
      iarea(k,nchem) = 0
      volin(k,nchem) = 0.0d0
      areain(k,nchem) = 1.0d0
      MineralMoles(k,nchem) = 0.0d0
    END IF
  END DO

! Calculate the porosity for this geochemical condition if Jpor = 0 or Jpor = 1 
!   (set in POROSITY keyword block)

  IF (jpor == 1 .OR. jpor == 0 .OR. jpor== 3) THEN
    sum = 0.0
    DO k = 1,nkin
      sum = sum + volin(k,nchem)
    END DO
    porcond(nchem) = 1.0 - sum
  END IF


!! Look for specification of solid density
!! No input:  Calculate solid density from sum of mineral volume fractions and their molar volumes
!! Solid_density:  Overrides mineral volume fractions
!! SolidDensityFrom(nchem) = 0, Calculate from mineral volume fractions
!! SolidDensityFrom(nchem) = 1, from specification of solid density
!! SolidDensityFrom(nchem) = 2, from specification of solid:solution ratio

  CALL read_SolidDensity(nout,nchem)

  DO k = 1,nkin
    IF (MineralMoles(k,nchem) > 0.0d0 .AND. volin(k,nchem) == 0.0d0) THEN
      volin(k,nchem) = MineralMoles(k,nchem)*porcond(nchem)*rocond(nchem)*SaturationCond(nchem)*volmol(k)
    END IF
  END DO

!!  SolidDensity(nchem) = 0.0d0

  IF (SolidDensityFrom(nchem) == 3) THEN                     !!  Calculate solid density and solid:solution ratio from mineral volume fractions

    sum = 0.0d0
    DO k = 1,nkin
      sum = sum + volin(k,nchem)
    END DO
    IF (sum == 0.0d0) THEN 
      SolidDensity(nchem) = 0.0d0
    ELSE
      TotalVolumeMinerals = 0.0d0
      DO k = 1,nkin
        TotalVolumeMinerals = TotalVolumeMinerals + volin(k,nchem)
        IF (volmol(k) /= 0.0d0) THEN
          SolidDensity(nchem) = SolidDensity(nchem) + (volin(k,nchem)/sum)*0.001d0*wtmin(k)/volmol(k) 
        END IF
      END DO
!!      porcond(nchem) = 1.0d0 - TotalVolumeMinerals
    END IF
    SolidSolutionRatio(nchem) = OneOverMassFraction(nchem)*1000.d0*SolidDensity(nchem)*(1.0-porcond(nchem))/(SaturationCond(nchem)*porcond(nchem)*rocond(nchem))
    IF (porcond(nchem) == 1.0d0 .AND. SolidDensity(nchem) /= 0.0d0) THEN
      WRITE(*,*)
      WRITE(*,*) ' Solid density is non-zero, but porosity = 100% (no solids)' 
      WRITE(*,*) ' In Geochemical Condition: ', nchem
      WRITE(*,*) ' Enter to continue'
      READ(*,*)
   END IF

  ELSE IF (SolidDensityFrom(nchem) == 1) THEN        !! Solid density specified, calculate solid solution ratio

    SolidSolutionRatio(nchem) = OneOverMassFraction(nchem)*1000.d0*SolidDensity(nchem)*(1.0-porcond(nchem))/(SaturationCond(nchem)*porcond(nchem)*rocond(nchem))   
    IF (porcond(nchem) == 1.0d0) THEN
      WRITE(*,*)
      WRITE(*,*) ' Solid density is non-zero, but porosity = 100% (no solids)' 
      WRITE(*,*) ' In Geochemical Condition: ', nchem
      WRITE(*,*) ' Enter to continue'
      READ(*,*)
    END IF

  ELSE IF (SolidDensityFrom(nchem) == 2) THEN          !! Solid solution ratio specified

    IF (porcond(nchem) == 1.0d0) THEN                  !! About to divide by zero
      SolidDensity(nchem) = 2700.0d0                   !! Ignoring porosity specification

!!      WRITE(*,*) 
!!      WRITE(*,*) ' --> 0% solids (100% porosity) incompatible with solid:solution ratio'
!!      WRITE(*,*) ' ----> Solid density computed is infinite'
!!      WritePorosity = 1.0d0/(1.0d0 + SolidSolutionRatio(nchem)*rocond(nchem)/(2700.0d0*1000.0d0) )
!!      WRITE(*,*) ' For this solid:solution ratio and '
!!      WRITE(*,*) '   a solid density = 2700 kg/m^3, Porosity = ', WritePorosity
!!      WRITE(*,*)
!!      READ(*,*)
!!      STOP

    ELSE
      SolidDensity(nchem) = SolidSolutionRatio(nchem)*(SaturationCond(nchem)*porcond(nchem)*rocond(nchem))/(OneOverMassFraction(nchem)*1000.d0*(1.0-porcond(nchem)))
    END IF

  ELSE
    CONTINUE
  END IF

!  Note:  Density calculation may depend on solute concentrations.  Here, the value of CTOT is used since
!         the exact concentration has not been calculated yet

!!NOTE:  Need to change so that CEC is in units of g/L and then do the change to g/kgw

! Look for exchange species
  
  DO ix = 1,nexchange
    IF (iexchange(ix) == 0) THEN                                          !!  Bulk exchange on the sediment
      CALL read_ionexchangeBS(nout,ix,nchem,ncomp,nexchange,speciesfound)
      IF (icec(ix) == 1) THEN
        totexch(ix,nchem) = cec(ix,nchem)*SolidSolutionRatio(nchem)
!!          SolidSolutionRatio(nchem) = OneOverMassFraction(nchem)*1000.d0*SolidDensity(nchem)*(1.0-porcond(nchem))/(porcond(nchem)*rocond(nchem))
      ELSE
        IF (SolidSolutionRatio(nchem) == 0.0d0) THEN
          cec(ix,nchem) = 0.0d0
        ELSE
          cec(ix,nchem) = totexch(ix,nchem)/SolidSolutionRatio(nchem)
        END IF
      END IF
    ELSE IF (iexchange(ix) == 1) THEN   !  Exchange on a specific mineral
      k = kexch(ix)
      CALL read_ionexchangeMIN(nout,ix,nchem,ncomp,nexchange,speciesfound)
      IF (icec(ix) == 1) THEN                                             !!  Calculate from equivalents/g mineral and the mineral volume fraction
!!        IF (volin(k,nchem) <= 0.0d0) THEN
!!          dumstring = umin(k)
!!          CALL stringlen(dumstring,ls)
!!          WRITE(*,*)
!!          WRITE(*,*) ' Initial volume fraction exchanger mineral must be > 0'
!!          WRITE(*,*) ' Exchange on mineral: ', dumstring(1:ls)
!!          WRITE(*,*) ' In geochemical condition: ', nchem
!!          WRITE(*,*)
!!          READ(*,*)
!!          STOP
!!        END IF
        totexch(ix,nchem) = OneOverMassFraction(nchem)*cec(ix,nchem)*wtmin(k)*volin(k,nchem)/(volmol(k)*SaturationCond(nchem)*porcond(nchem)*rocond(nchem))  
      ELSE IF (icec(ix) == 0) THEN                                         !!  Direct specification of totexch (equivalents/kgw)--need to calculate CEC for later use if mineral fraction changes
        cec(ix,nchem) = totexch(ix,nchem)*volmol(k)*SaturationCond(nchem)*porcond(nchem)*rocond(nchem)/(wtmin(k)*volin(k,nchem)*OneOverMassFraction(nchem))
      ELSE
        WRITE(*,*)
        WRITE(*,*) ' CEC option not recognized'
        WRITE(*,*)
        STOP
      END IF
      
    ELSE
      WRITE(*,*) ' Dont recognize type of ion exchange'
      WRITE(*,*) ' Should be on a mineral or bulk sediment'
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
        
    IF (speciesfound) THEN
      CONTINUE
    ELSE
      WRITE(*,*)
      WRITE(*,*) ' Exchange species not found'
      WRITE(*,*) ' Species = ',namexc(ix)
      WRITE(*,*) ' In geochemical condition ',labeltemp(1:lcond)
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
    

  END DO

! Read surface complexes
  
DO ks = 1,nsurf
    complexfound = .false.
!    WRITE(*,*)
!    WRITE(*,*) ' Starting to read surface complexes'
!    WRITE(*,*) ' Looking for ',namsurf(ks)
!    WRITE(*,*)
    CALL read_surfacecomplex(nout,ks,nchem,ncomp,nkin,ngas,complexfound)
    IF (complexfound) THEN
      CONTINUE
    ELSE
      dumstring = namsurf(ks)
      WRITE(*,*)
      WRITE(*,*) ' Surface complex not found: ',dumstring(1:15)
      WRITE(*,*) ' In Condition: ',nchem
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
  END DO
  
!  Calculate concentration of surface hydroxyls per kg H2O

  DO ks = 1,nsurf
    k = ksurf(ks)
    IF (specific(k,nchem) == 0.0) THEN
      WRITE(*,*) 
      WRITE(*,*) ' Specific surface area for mineral = 0 serving as sorbate for surface complex'
      dumstring = umin(k)
      CALL stringlen(dumstring,ls)
      WRITE(*,*) ' --->Mineral: ', dumstring(1:ls)
      WRITE(*,*) ' --->In geochemical condition: ', labeltemp(1:lcond)
      WRITE(*,*) ' Use "specific_surface_area" followed by non-zero number for mineral'
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
      
!   Site_density: moles site/m^2 mineral
!   Specific:     m^2/g mineral
!   Wtmin:        g/mole mineral

    permole = site_density(ks,nchem)*specific(k,nchem)*wtmin(k)    !  Mole sites/Mole mineral

!   Now convert to moles sites per kg solution
!!  volin(m^3 mineral/m^3 porous medium) /( volmol[m^3/mol] * rocond[kg/m^3 fluid] * porcond[m^3 pore/m^3 PM] * SaturationCond[m^3 fluid/m^3 pore)

!!  Moles mineral/kgw
    IF (MineralMoles(k,nchem) > 0.0d0 .AND. volin(k,nchem) == 0.0d0) THEN
      c_surf(ks,nchem) = permole*MineralMoles(k,nchem)
    ELSE IF (volin(k,nchem) == 0.0d0) THEN
      MineralMolality = voltemp(k,nchem)/( volmol(k)*rocond(nchem)*porcond(nchem)*SaturationCond(nchem) )
      c_surf(ks,nchem) = permole*MineralMolality
    ELSE
      MineralMolality = volin(k,nchem)/( volmol(k)*rocond(nchem)*porcond(nchem)*SaturationCond(nchem) )
      c_surf(ks,nchem) = permole*MineralMolality
    END IF

    IF (c_surf(ks,nchem) < 1.D-30) THEN
      c_surf(ks,nchem) = 1.D-30
    END IF

!    WRITE(*,*) ' k = ',k
!    WRITE(*,*) ' Site density ',site_density(ks,nchem)
!    WRITE(*,*) ' Specific area ',specific(ks,nchem)
!    WRITE(*,*) ' Wtmin = ',wtmin(k)
!    WRITE(*,*) ' Permole = ',permole
!    WRITE(*,*) ' Volin = ',volin(k,nchem)
!    WRITE(*,*) ' volmol = ',volmol(k)
  END DO
  
  IF (nsurf > msurf) THEN
    WRITE(*,*)
    WRITE(*,*) ' Parameter "msurf" dimensioned too small '
    WRITE(*,*) ' Msurf = ',msurf
    WRITE(*,*) ' Number of surface complexes = ',nsurf
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
  
  DO ks = 1,nsurf
    guess_surf(ks,nchem) = c_surf(ks,nchem)
  END DO
  
! ***********  Look for the topes  ***************
  
!  DO id = 1,ndecay
!    i = idecay(id)
    
!    topefound = .false.
    
!  First, check to see if the decay element is followed by "all"; in this case, use ratio
!      for all of the minerals
    
!    idecayall = 1
    
!    CALL read_toperatio(nout,ncomp,ndecay,idecayall,id,kd,nchem,topefound)
    
!    IF (topefound) THEN
!      CONTINUE
!    ELSE
!      idecayall = 0
!      DO kd = 1,nmindecay(id)
!        k = kdecay(kd,id)
!        CALL read_toperatio(nout,ncomp,ndecay,idecayall,id,kd,nchem,topefound)
!        IF (topefound) THEN
!          CONTINUE
!        ELSE
!          WRITE(*,*)
!          WRITE(*,*)  ' Isotope ratio not found'
!          WRITE(*,*)
!          STOP
!        END IF
!      END DO
!    END IF
!  END DO
  
  equilfound = .FALSE.
  CALL read_equil(nout,nchem,equilfound)
  IF (equilfound) THEN
    equilibrate(:,nchem) = .TRUE.
  END IF
  
  
!  Check for more geochemical conditions
  
  GO TO 200
  
ELSE
  WRITE(*,*)
  WRITE(*,*) ' Failed to find any geochemical conditions'
  WRITE(*,*) ' This is going to be a boring simulation'
  WRITE(*,*) '             BAG IT!'
  READ(*,*)
  STOP
END IF

5022 FORMAT(1X,'Geochemical condition number ',i2,' found')

RETURN
END SUBROUTINE find_condition
