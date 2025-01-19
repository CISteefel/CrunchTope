
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
SUBROUTINE equilib_co2(ncomp,nspec,nrct,ngas,nsurf,igamma,ikph,     &
    nco,nexchange,nexch_sec,nsurf_sec,npot,neqn,tempc,portemp,      &
    DensityModule,ipest,PestUnit,pest)
USE crunchtype
USE params
USE concentration
USE mineral
USE solver
USE io
USE temperature
USE medium
USE runtime, ONLY: HanfordStrontium, Duan, Duan2006,SkipAdjust

IMPLICIT NONE

EXTERNAL dgetrf
EXTERNAL dgetrs

!  ********************  Beginning of interface blocks  ******************
interface
  SUBROUTINE ludcmp90(a,indx,d,n)
  USE crunchtype
  REAL(DP), DIMENSION(:,:), intent(in out)                   :: a
  INTEGER(I4B), DIMENSION(:), intent(out)                    :: indx
  INTEGER(I4B), INTENT(IN)                                   :: n
  REAL(DP), intent(out)                                      :: d
  END SUBROUTINE ludcmp90
END interface

interface
  SUBROUTINE lubksb90(a,indx,b,n)
  USE crunchtype
  REAL(DP), DIMENSION(:,:), intent(in)                       :: a
  INTEGER(I4B), DIMENSION(:), intent(in)                     :: indx
  REAL(DP), DIMENSION(:), intent(inout)                      :: b
  INTEGER(I4B), INTENT(IN)                                   :: n
  END SUBROUTINE lubksb90
END interface

interface
  SUBROUTINE gaussj(a,b,n)
  USE crunchtype
  IMPLICIT NONE
  REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: a,b
  INTEGER(I4B), INTENT(IN)                :: n
  END SUBROUTINE gaussj
END interface

!  ********************  End of interface blocks  ******************

!  External variables

INTEGER(I4B), INTENT(IN)                                   :: ncomp
INTEGER(I4B), INTENT(IN)                                   :: nspec
INTEGER(I4B), INTENT(IN)                                   :: nrct
INTEGER(I4B), INTENT(IN)                                   :: ngas
INTEGER(I4B), INTENT(IN)                                   :: nsurf
INTEGER(I4B), INTENT(IN)                                   :: igamma
INTEGER(I4B), INTENT(IN)                                   :: ikph
INTEGER(I4B), INTENT(IN)                                   :: nco
INTEGER(I4B), INTENT(IN)                                   :: nexchange
INTEGER(I4B), INTENT(IN)                                   :: nexch_sec
INTEGER(I4B), INTENT(IN)                                   :: nsurf_sec
INTEGER(I4B), INTENT(IN)                                   :: npot
INTEGER(I4B), INTENT(IN)                                   :: neqn

REAL(DP), INTENT(IN)                                       :: tempc
REAL(DP), INTENT(IN)                                       :: portemp


CHARACTER (LEN=mls)                                        :: DensityModule

INTEGER(I4B), INTENT(IN)                                   :: ipest
INTEGER(I4B), INTENT(IN)                                   :: PestUnit
LOGICAL(LGT), INTENT(IN)                                   :: pest

!  Internal variables and arrays

REAL(DP)                                                   :: pg

REAL(DP)                                                   :: PressureTemp

REAL(DP), DIMENSION(:), ALLOCATABLE                        :: u
REAL(DP), DIMENSION(:), ALLOCATABLE                        :: feq
REAL(DP), DIMENSION(:), ALLOCATABLE                        :: beta
REAL(DP), DIMENSION(:,:), ALLOCATABLE                      :: betagauss
REAL(DP), DIMENSION(:), ALLOCATABLE                        :: spbase
REAL(DP), DIMENSION(:), ALLOCATABLE                        :: ffscale
REAL(DP), DIMENSION(:), ALLOCATABLE                        :: sind
REAL(DP), DIMENSION(:), ALLOCATABLE                        :: tol
REAL(DP), DIMENSION(:), ALLOCATABLE                        :: totex_bas
REAL(DP), DIMENSION(:), ALLOCATABLE                        :: totequiv
REAL(DP), DIMENSION(:), ALLOCATABLE                        :: totSurf_bas

REAL(DP), DIMENSION(:,:), ALLOCATABLE                      :: fj
REAL(DP), DIMENSION(:,:), ALLOCATABLE                      :: fj_num

INTEGER(I4B), DIMENSION(:), ALLOCATABLE                    :: indx
INTEGER(I4B), DIMENSION(:), ALLOCATABLE                    :: iTotNegative

CHARACTER (LEN=mls)                                        :: namtemp
CHARACTER (LEN=mls)                                        :: dumstring

REAL(DP), PARAMETER                                        :: tolf=2.e-14
REAL(DP)                                                   :: tk
REAL(DP)                                                   :: sumchg
REAL(DP)                                                   :: sum
REAL(DP)                                                   :: check
REAL(DP)                                                   :: sumiap
REAL(DP)                                                   :: fxmaxx
REAL(DP)                                                   :: det
REAL(DP)                                                   :: errx
REAL(DP)                                                   :: betamax
REAL(DP)                                                   :: atol
REAL(DP)                                                   :: sion_tmp
REAL(DP)                                                   :: pHprint
REAL(DP)                                                   :: peprint
REAL(DP)                                                   :: Ehprint
REAL(DP)                                                   :: FaradayKJ
REAL(DP)                                                   :: FaradayKcal
REAL(DP)                                                   :: prt_bas
REAL(DP)                                                   :: spprint
REAL(DP)                                                   :: totcharge
REAL(DP)                                                   :: silnTMP
REAL(DP)                                                   :: siprnt
REAL(DP)                                                   :: sqrt_sion
REAL(DP)                                                   :: delta_z
REAL(DP)                                                   :: feq_tmp
REAL(DP)                                                   :: temp1
REAL(DP)                                                   :: actprint
REAL(DP)                                                   :: actprint10
REAL(DP)                                                   :: faraday
REAL(DP)                                                   :: gramsperL
REAL(DP)                                                   :: ChargeCheck
REAL(DP)                                                   :: AqueousToBulk
REAL(DP)                                                   :: actcoeffprint
REAL(DP)                                                   :: MolePerGSolids
REAL(DP)                                                   :: EquivPerGSolids
REAL(DP),DIMENSION(nrct)                                                   :: SurfaceChargeOriginal
REAL(DP)                                                   :: deltaSurfaceCharge

REAL(DP)                                                   :: CheckValue
REAL(DP)                                                   :: rone
REAL(DP)                                                   :: eps = 1.0D-11

REAL(DP), PARAMETER                                        :: perturb=1.e-07
REAL(DP), PARAMETER                                        :: corrmax=2.0
REAL(DP), PARAMETER                                        :: SmallLogNumber=-35.0

INTEGER(I4B), PARAMETER                                    :: ntrial=5000
INTEGER(I4B), PARAMETER                                    :: ione=1
INTEGER(I4B)                                               :: i
INTEGER(I4B)                                               :: i2
INTEGER(I4B)                                               :: nmass
INTEGER(I4B)                                               :: nmineq
INTEGER(I4B)                                               :: ngaseq
INTEGER(I4B)                                               :: nfixact
INTEGER(I4B)                                               :: nfixcon
INTEGER(I4B)                                               :: ncharge
INTEGER(I4B)                                               :: ik
INTEGER(I4B)                                               :: ihalf
INTEGER(I4B)                                               :: ktrial
INTEGER(I4B)                                               :: ix
INTEGER(I4B)                                               :: is
INTEGER(I4B)                                               :: k
INTEGER(I4B)                                               :: kk
INTEGER(I4B)                                               :: kg
INTEGER(I4B)                                               :: ichg
INTEGER(I4B)                                               :: i3
INTEGER(I4B)                                               :: ix2
INTEGER(I4B)                                               :: is2
INTEGER(I4B)                                               :: nex
INTEGER(I4B)                                               :: ns
INTEGER(I4B)                                               :: lscond
INTEGER(I4B)                                               :: npt
INTEGER(I4B)                                               :: npt2
INTEGER(I4B)                                               :: nrow
INTEGER(I4B)                                               :: ncol
INTEGER(I4B)                                               :: ksp
INTEGER(I4B)                                               :: ks
INTEGER(I4B)                                               :: info
INTEGER(I4B)                                               :: nl
INTEGER(I4B)                                               :: ico2
INTEGER(I4B)                                               :: iH2O

LOGICAL(LGT)                                               :: bagit
LOGICAL(LGT)                                               :: ChargeBalance
LOGICAL(LGT)                                               :: NeedChargeBalance
LOGICAL(LGT)                                               :: ChargeOK
LOGICAL(LGT)                                               :: InitializeMineralEquilibrium
LOGICAL(LGT)                                               :: ikH2Ofound


REAL(DP)                                                   :: MeanSaltConcentration
REAL(DP)                                                   :: MassFraction
REAL(DP)                                                   :: Kd
REAL(DP)                                                   :: denmol
REAL(DP)                                                   :: KvNaSr
REAL(DP)                                                   :: KvNaCa
REAL(DP)                                                   :: LogKvNaSr
REAL(DP)                                                   :: LogKvNaCa
REAL(DP)                                                   :: TotalExchangeEquivalents
REAL(DP)                                                   :: TotalExchangeMoles

REAL(DP)                                                   :: vrInOut
REAL(DP)                                                   :: fxMaxPotential
REAL(DP)                                                   :: ChargeSum
REAL(DP)                                                   :: lnActivity
REAL(DP)                                                   :: gammawatertmp

CHARACTER (LEN=3)                                          :: ulabPrint


CHARACTER (LEN=1)                                          :: trans

REAL(DP)                                                   :: ph2o, ln_fco2, fco2, ln_fco2_check

!!real(dp), dimension(26) :: press
!!real(dp), dimension(10) :: tempe, fugac
!!integer :: iii,jjj

vrInOut = 1.0

IF (nco==1) THEN
  IF (ALLOCATED(fweight)) THEN
    DEALLOCATE(fweight)
  END IF
  ALLOCATE(fweight(neqn,neqn))
END IF

trans = 'N'

bagit = .FALSE.

ALLOCATE(u(neqn))
ALLOCATE(feq(neqn))
ALLOCATE(beta(neqn))
ALLOCATE(betagauss(neqn,1))
ALLOCATE(spbase(ncomp+nspec))
ALLOCATE(ffscale(neqn))
ALLOCATE(sind(nrct))
ALLOCATE(tol(neqn))
ALLOCATE(totex_bas(ncomp))
ALLOCATE(totequiv(ncomp))
ALLOCATE(totSurf_bas(ncomp))
ALLOCATE(fj(neqn,neqn))
ALLOCATE(fj_num(neqn,neqn))
ALLOCATE(indx(neqn))
ALLOCATE(iTotNegative(ncomp))
iTotNegative = 0

tk = tempc + 273.15
denmol = 1.e05/(8.314*tk)   ! P/RT = n/V, with pressure converted from bars to Pascals
faraday = 96485.0 
rone = 1.0 
NeedChargeBalance = .FALSE.
ChargeOK = .TRUE.
ChargeBalance = .TRUE.

IF (DensityModule /= 'temperature') THEN
! Calculate the correction for the mass fraction of water:  kg_solution/kg_water
  MeanSaltConcentration = 0.001*(wtaq(MeanSalt(1))*ctot(MeanSalt(1),nco) +   &
        wtaq(MeanSalt(2))*ctot(MeanSalt(2),nco)) 
  MassFraction = 1.0/(1.0 + MeanSaltConcentration)
ELSE
  MassFraction = 1.0
END IF
AqueousToBulk = MassFraction*SaturationCond(nco)*porcond(nco)*rocond(nco)

sqrt_sion = 0.0

tol = 1.0
fexch = 0.0

nmass = 0
nmineq = 0
ngaseq = 0
nfixact = 0
nfixcon = 0
ncharge = 0

dumstring = condlabel(nco)
CALL stringlen(dumstring,lscond)

!  Add up the number of constraints and compare to the number of
!  unknowns.

NeedChargeBalance = .FALSE.

pg = 0.0d0                         !!!co2

DO i = 1,ncomp
  IF (itype(i,nco) == 1) THEN
    nmass = nmass + 1
    IF (ctot(i,nco) > 0.0) THEN
      tol(i) = ctot(i,nco)
    ELSE
      tol(i) = 1.0
    END IF
  ELSE IF (itype(i,nco) == 3) THEN
    nmineq = nmineq + 1
  ELSE IF (itype(i,nco) == 7) THEN
    nfixact = nfixact + 1
  ELSE IF (itype(i,nco) == 8) THEN
    nfixcon = nfixcon + 1
  ELSE IF (itype(i,nco) == 4) THEN
    ngaseq = ngaseq + 1
    
    pg = pg + ctot(i,nco)   !! Not just CO2, but all of the gases.  These are augmented below by the partial pressures of gases not specified as constraints
    
  ELSE IF (itype(i,nco) == 2) THEN
    
    ncharge = ncharge + 1
    NeedChargeBalance = .TRUE.
    
  ELSE
    WRITE(*,*)
    WRITE(*,*) ' Incorrect option for initialization specified'
    WRITE(*,*) ulab(i),itype(i,nco)
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
END DO

!!! Find species number for H2O
ikH2Ofound = .FALSE.
DO I = 1,ncomp
  IF ( ulab(i) == 'H2O' .or. ulab(i) == 'h2o' .or. ulab(i) == 'hho' .or. ulab(i) == 'HHO' ) THEN
    ih2o = i
    ikH2Ofound = .TRUE.
  END IF
END DO
  
IF (ikH2Ofound) THEN
  CONTINUE
ELSE
  write(*,*)
  write(*,*) ' Water species number not found'
  write(*,*)
  read(*,*)
  stop
END IF

IF (Duan .OR. Duan2006) THEN            !! Duan CO2 routine
 
  ph2o = 0.0d0
  call calc_ph2o(tk,ph2o)
  pg = pg + ph2o

! need to re-calculate Keq for CO2 for this condition (pg) 
  CALL keqcalc2_init_co2(ncomp,nrct,nspec,ngas,nsurf_sec,tempc,pg)
 
END IF

!!  Check to see if there are any negative coefficients appearing in stoichiometric coefficient matrix

DO i = 1,ncomp
  DO ksp = 1,nspec
    IF (muaq(ksp,i) < 0.0) THEN
       iTotNegative(i) = 1
    END IF
  END DO
  DO ks = 1,nsurf_sec
    IF (musurf(ks,i) < 0.0) THEN
       iTotNegative(i) = 1
    END IF
  END DO
  DO nex = 1,nexch_sec
    IF (muexc(nex,i) < 0.0) THEN
       iTotNegative(i) = 1
    END IF
  END DO
END DO

sumchg = 0.0
DO ik = 1,ncomp+nspec
  sumchg = sumchg + sptmp10(ik)*chg(ik)
END DO

DO i = 1,ncomp
  
  IF (itype(i,nco) == 2) THEN
    IF (chg(i) < 0.0) THEN
      IF (sumchg > 0.0) THEN
        sptmp10(i) = sumchg
        sptmp(i) = LOG(sptmp10(i))
      END IF
    ELSE IF (chg(i) > 0.0) THEN
      IF (sumchg < 0.0) THEN
        sptmp10(i) = -sumchg
        sptmp(i) = LOG(sptmp10(i))
      END IF
    ELSE
      WRITE(*,*)
      WRITE(*,*) ' Cant charge balance on a neutral species'
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
  END IF
END DO

IF (ncomp /= nmass+nmineq+ngaseq+ncharge+nfixact+nfixcon) THEN
  WRITE(*,*) 'There must be a condition for each component'
  READ(*,*)
  STOP
END IF

IF (ncharge > 1) THEN
  WRITE(*,*) 'You can only balance the charge with one equation'
  READ(*,*)
  STOP
END IF

!  Get some better guesses for pH and O2(aq)

DO i = 1,ncomp
  IF (itype(i,nco) == 1 .AND. guess(i,nco) > (1.E-06+eps) .AND. guess(i,nco) < (1.0E-06+eps)) THEN
    IF (ctot(i,nco) > 0.0) THEN
      sptmp(i) = LOG(ctot(i,nco))
      sptmp10(i) = EXP(sptmp(i))
    ELSE
      IF (ulab(i) == 'o2(aq)' .OR. ulab(i) == 'O2(aq)') THEN
        IF (ctot(i,nco) <= 0.0) THEN
          sptmp(i) = clg*(-70.0)
          sptmp10(i) = EXP(sptmp(i))
        ELSE
          sptmp(i) = clg*(-7.0)
          sptmp10(i) = EXP(sptmp(i))
        END IF
      ELSE IF (ulab(i) == 'h+' .OR. ulab(i) == 'H+') THEN
        IF (ctot(i,nco) <= 0.0) THEN
          sptmp(i) = clg*(-10.0)
          sptmp10(i) = EXP(sptmp(i))
        ELSE
          sptmp(i) = clg*(-3.0)
          sptmp10(i) = EXP(sptmp(i))
        END IF
      END IF
    END IF
  END IF
  
!!!  IF (ulab(i) == 'H2O' .AND. itype(i,nco) == 1) THEN
  IF (ulab(i) == 'H2O') THEN
    sptmp10(i) = 55.50843506
    sptmp(i) = DLOG(sptmp10(i))
  END IF
END DO

ihalf = 0
212 bagit = .FALSE.

!! Find out component number for carbonate species so as to include CO2(g) fugacity coefficient

  ico2 = 0
  DO i = 1,ncomp
    IF (ulab(i) == 'CO2(aq)' .OR. ulab(i) == 'HCO3-' .OR. ulab(i) == 'CO3--') THEN
      ico2 = i
    END IF
  END DO

!  **********Start of NEWTON-RAPHSON routine***************
 
DO  ktrial = 1,ntrial

  IF (Duan .OR. Duan2006) THEN    !!  Duan CO2 routine
 
    IF (pg == 0.0d0) pg = 1.0
    CALL gases_init_co2(ncomp,ngas,tempc,pg,vrInOut)
    sum = 0.0d0
    DO kk = 1,ngas
      sum = sum + spgastmp10(kk)
    END DO
    pg = sum

    ph2o = 0.0d0
    CALL calc_ph2o(tk,ph2o)
    pg = pg + ph2o

!!  need to re-calculate Keq for CO2 for this condition (pg) 
    CALL keqcalc2_init_co2(ncomp,nrct,nspec,ngas,nsurf_sec,tempc,pg)
 
  end if

  IF (ihalf > 2000) THEN

    WRITE(*,*)
    WRITE(*,*) ' Halved the concentrations too many times in condition: ', dumstring(1:lscond)
    WRITE(*,*) ' Exiting EQUILIBRATION routine'
    WRITE(*,*)
    GOTO 3001
!!    STOP
  END IF

!!!  IF (igamma /= 0 .AND. ktrial > 5) THEN
  IF (igamma /= 0) THEN
    
    if (Duan .OR. Duan2006) then
      call gamma_init_co2(ncomp,nspec,tempc,sqrt_sion,pg)
    else
      CALL gamma_init(ncomp,nspec,igamma,tempc,sqrt_sion,sion_tmp)
    end if
  ELSE

    ChargeSum = 0.0d0
    DO ik = 1,ncomp+nspec
      ChargeSum = ChargeSum + sptmp10(ik)*chg(ik)*chg(ik)
    END DO
    sion_tmp = 0.50D0*ChargeSum
    
    IF (sion_tmp < 25.0d0) THEN
      sqrt_sion = DSQRT(sion_tmp)
    ELSE
      sqrt_sion = 0.0d0
    END IF

  END IF
  
  IF (SaltCreep) THEN
    PressureTemp = PressureCond(nco)
    call keqcalc2_init(ncomp,nrct,nspec,ngas,nsurf_sec,tempc,PressureTemp)
  END IF

  CALL species_init(ncomp,nspec)
  
  if (Duan .OR. Duan2006) then
    CALL gases_init_co2(ncomp,ngas,tempc,pg,vrInOut)
  else
    CALL gases_init(ncomp,ngas,tempc,nco)
  end if
  
  CALL surf_init(ncomp,nspec,nsurf,nsurf_sec,nco)
  CALL exchange_init(ncomp,nspec,nexchange,nexch_sec,nco)
  CALL totconc_init(ncomp,nspec,nexchange,nexch_sec,nsurf,nsurf_sec,nco)
  CALL totgas_init(ncomp,nspec,ngas)
  CALL totsurf_init(ncomp,nsurf,nsurf_sec)
!        call exactivity_init(ncomp,nexchange,nexch_sec,nco)
  CALL totexchange_init(ncomp,nexchange,nexch_sec,nco)
  CALL jacobian_init(ncomp,nspec)
  CALL jacexchange_init(ncomp,nexchange,nexch_sec,neqn,nco)
  CALL jacsurf_init(ncomp,nsurf,nsurf_sec,nco)
  
  IF (ktrial == 1) THEN
    DO i = 1,ncomp
      u(i) = sptmp(i)
    END DO
    DO ix = 1,nexchange
      u(ix+ncomp) = spextmp(ix)
    END DO
    DO is = 1,nsurf
      u(is+ncomp+nexchange) = spsurftmp(is)
    END DO
    DO npt = 1,npot
      u(npt+ncomp+nexchange+nsurf) = 0.0
    END DO
  END IF
  
!  Calculate residuals
  
  ChargeBalance = .TRUE.

  DO  i = 1,ncomp
    
    bagit = .FALSE.
    IF (itype(i,nco) == 1) THEN
      
      IF (ctot(i,nco) /= 0.0) THEN
        check = stmp(i)/ctot(i,nco)  !! Ratio of calculated total concentration to imposed
        continue
      END IF
      feq(i) = stmp(i) - ctot(i,nco)
       
      IF (check > 1000000.0 .AND. iTotNegative(i) == 0) THEN
        
        sptmp(i) = sptmp(i) - clg    !! If calculated total concentration too high, reduce primary species
        sptmp10(i) = DEXP(sptmp(i))
        
        if (ihalf < 2000) then
          ihalf = ihalf + 1
          bagit = .TRUE.
        end if
        
      ELSE IF (check < 1.e-07  .AND. iTotNegative(i) == 0) THEN
        
        sptmp(i) = sptmp(i) + clg
        sptmp10(i) = DEXP(sptmp(i))
        if (ihalf < 2000) then
          ihalf = ihalf + 1
          bagit = .TRUE.
        end if
      END IF

      IF (bagit) THEN
         EXIT
      END IF

    ELSE IF (itype(i,nco) == 3) THEN

      DO k = 1,nrct
        IF (umin(k) == ncon(i,nco)) THEN
          kk = k + nspec
          GO TO 400
        END IF
      END DO
      
      WRITE(*,*)
      WRITE(*,*)
      WRITE(*,*) 'Constraint mineral not listed in input file'
      WRITE(*,542) ulab(i)
      WRITE(*,544) nco
      WRITE(*,543) ncon(i,nco)
      WRITE(*,*)
      READ(*,*)
      STOP
      
      400       CONTINUE
      sumiap = 0.0D0
      lnActivity = 0.0d0
      DO  i2 = 1,ncomp
        
        ulabPrint = ulab(i2)
        IF (ulabPrint(1:3) == 'H2O' .or. ulabPrint(1:3) == 'HHO') THEN
          lnActivity = gamtmp(i2) 
         ELSE
          lnActivity = (sptmp(i2)+gamtmp(i2))
        endif

        sumiap = sumiap + mumin(1,k,i2) * lnActivity

      END DO

      feq(i) = sumiap - keqmin_tmp(1,k)

      IF (SkipAdjust(nco)) THEN
        CONTINUE
      ELSE

        IF (feq(i) > clg) THEN
           sptmp(i) = sptmp(i) - clg
           sptmp10(i) = DEXP(sptmp(i))
          if (ihalf < 2000) then
            ihalf = ihalf + 1
            bagit = .TRUE.
          end if
        END IF

        IF (bagit) THEN
          EXIT
        END IF

!!!         sind(k) = feq(i)/clg

       END IF

    ELSE IF (itype(i,nco) == 4) THEN             !!! Gas constraint

      DO kg = 1,ngas
        IF (namg(kg) == ncon(i,nco)) THEN
          kk = kg + nspec + nrct
          GOTO 500
        END IF
      END DO
      
      WRITE(*,*)
      WRITE(*,*)
      WRITE(*,*) 'Constraint gas not listed in input file'
      WRITE(*,542) ulab(i)
      WRITE(*,544) nco
      WRITE(*,643) ncon(i,nco)
      WRITE(*,*)
      WRITE(*,*)
      READ(*,*)
      STOP
      
      500       CONTINUE
      sumiap = 0.0D0
      lnActivity = 0.0d0
      DO i2 = 1,ncomp
        
        ulabPrint = ulab(i2)
        IF (ulabPrint(1:3) == 'H2O' .or. ulabPrint(1:3) == 'HHO') THEN
          lnActivity = gamtmp(i2) 
!!!         lnActivity = 0.0d0
        ELSE
          lnActivity = (sptmp(i2)+gamtmp(i2))
        endif

        sumiap = sumiap + mugas(kg,i2) * lnActivity

      END DO

      IF (Duan) THEN
        ln_fco2 = 0.0d0
        IF (i == ico2) THEN
          ln_fco2 = 0.0d0   ! fugacity coefficient for CO2(g)
          CALL fugacity_co2(pg,tk,ln_fco2,vrInOut)  
        END IF
      
        check = sumiap + keqgas_tmp(kg) - ln_fco2 - LOG(ctot(i,nco))
        
      ELSE IF (Duan2006) THEN
        ln_fco2 = 0.0d0
        IF (i == ico2) THEN
          ln_fco2 = 0.0d0   ! fugacity coefficient for CO2(g)
          CALL fugacity_co24(pg,tk,ln_fco2,vrInOut)  
        END IF
      
        check = sumiap + keqgas_tmp(kg) - ln_fco2 - LOG(ctot(i,nco))
      ELSE
        check = sumiap + keqgas_tmp(kg) - LOG(ctot(i,nco))
      END IF

      feq(i) = check

    ELSE IF (itype(i,nco) == 7) THEN
      
      lnActivity = 0.0d0
      ulabPrint = ulab(i)
      IF (ulabPrint(1:3) == 'H2O' .or. ulabPrint(1:3) == 'HHO') THEN
        
!!!        lnActivity = gamtmp(i)
        lnActivity = 0.0d0
      else 
        lnActivity = sptmp(i)+gamtmp(i)
      end if
      
      feq(i) = lnActivity - DLOG(guess(i,nco))
      
    ELSE IF (itype(i,nco) == 8) THEN
      
      feq(i) = sptmp(i) -  DLOG(guess(i,nco))
      
    ELSE IF (itype(i,nco) == 2) THEN
      
      IF (chg(i) < 0.0) THEN
        ichg = -1
      ELSE IF (chg(i) > 0.0) THEN
        ichg = 1
      ELSE
        WRITE(*,*) ' Cant charge balance on a neutral species'
        READ(*,*)
        STOP
      END IF

      sumchg = 0.0
      DO i2 = 1,ncomp
        IF (i2 /= i) THEN      !  Calculate charge before considering charge balancing species
          sumchg = sumchg + chg(i2)*stmp(i2)   
        END IF
      END DO
      IF (sumchg < 0.0 .AND. ichg == -1 .AND. itotNegative(i)==0) THEN
        ChargeBalance = .FALSE.
        feq(i) = 0.0
      ELSE IF (sumchg > 0.0 .AND. ichg == 1 .AND. itotNegative(i)==0) THEN
        ChargeBalance = .FALSE.
        feq(i) = 0.0
      ELSE
        sumchg = 0.0
        DO i2 = 1,ncomp+nspec
          sumchg = sumchg + chg(i2)*sptmp10(i2)
        END DO
        feq(i) = sumchg
      END IF

    ELSE
      WRITE(*,*) 'No other conditions possible'
      READ(*,*)
      STOP
    END IF
    
  END DO

  IF (bagit) then
    EXIT
  ENDIF
  
  bagit = .FALSE.

  DO ix = 1,nexchange

    check = sumactivity(ix)
    IF (bagit) THEN
      EXIT
    END IF
    feq(ix+ncomp) =  sumactivity(ix) - 1.0
    
  END DO

  IF (bagit) then
    EXIT
  ENDIF

  bagit = .FALSE.
  
  DO is = 1,nsurf
    IF (c_surf(is,nco) /= 0.0) THEN
      check = ssurftmp(is)/c_surf(is,nco)
    ELSE
      check = 1.0
    END IF
    feq(is+ncomp+nexchange) = ssurftmp(is) - c_surf(is,nco)

    IF (bagit) THEN
      EXIT
    END IF
  END DO

  IF (bagit) then
    EXIT
  ENDIF
  
!          ********** End of residual calculation ***************

  fj = 0.0d0

!  Form the Jacobian
  
  DO i = 1,ncomp
     
    IF (itype(i,nco) == 1) THEN

      DO i2 = 1,ncomp
        fj(i,i2) = dpsi(i,i2) + fexch(i,i2) + fsurftmp(i,i2)
      END DO

    ELSE IF (itype(i,nco) == 3) THEN
      
      DO k = 1,nrct
        IF (umin(k) == ncon(i,nco)) THEN
          kk = k + nspec
          GO TO 401
        END IF
      END DO
      
      WRITE(*,*) 'Constraint mineral not listed in input file'
      READ(*,*)
      STOP
401   CONTINUE
      
      DO i2 = 1,ncomp
        fj(i,i2) = mumin(1,k,i2)
      END DO

    ELSE IF (itype(i,nco) == 4) THEN
      
      DO kg = 1,ngas
        IF (namg(kg) == ncon(i,nco)) THEN
          kk = kg + nspec + nrct
          GO TO 501
        END IF
      END DO
      
      WRITE(*,*) 'Constraint gas not listed in input file'
      READ(*,*)
      STOP
501   CONTINUE
      
      DO i2 = 1,ncomp
        fj(i,i2) = mugas(kg,i2)
      END DO
      
    ELSE IF (itype(i,nco) == 7) THEN

      DO i2 = 1,ncomp
        IF (i2 == i) THEN
!!!          fj(i,i2) = DEXP(sptmp(i)+gamtmp(i))
!!!          fj(i,i2) = DEXP(gamtmp(i))
          fj(i,i2) = 1.0d0
        ELSE
          fj(i,i2) = 0.0
        END IF
      END DO
      
    ELSE IF (itype(i,nco) == 8) THEN

      DO i2 = 1,ncomp
        IF (i2 == i) THEN
!!          fj(i,i2) = sptmp10(i)
          fj(i,i2) = 1.0d0
        ELSE
          fj(i,i2) = 0.0
        END IF
      END DO
      
    ELSE IF (itype(i,nco) == 2) THEN

      IF (ChargeBalance) THEN
        DO i2 = 1,ncomp
          fj(i,i2) = 0.0
          DO ksp = 1,nspec
            fj(i,i2) = fj(i,i2) + chg(ksp+ncomp)*muaq(ksp,i2)*sptmp10(ksp+ncomp)
          END DO
          fj(i,i2) = fj(i,i2) + chg(i2)*sptmp10(i2)
        END DO
      ELSE
        DO i2 = 1,ncomp
!          fj(i,i2) = dpsi(i,i2) + fexch(i,i2) + fsurftmp(i,i2)
          fj(i,i2) = dpsi(i,i2)
        END DO    
      END IF    
    ELSE
      WRITE(*,*) 'No other conditions possible'
      READ(*,*)
      STOP
    END IF
    DO ix2 = 1,nexchange
!!      fj(i,ix2+ncomp) = fexch(i,ix2+ncomp)
      fj(i,ix2+ncomp) = 0.0
    END DO
    DO is2 = 1,nsurf
!!      fj(i,is2+ncomp+nexchange) = fsurftmp(i,is2+ncomp)
      fj(i,is2+ncomp+nexchange) = 0.0
    END DO
    
  END DO
  
  DO ix = 1,nexchange
    DO i2 = 1,ncomp+nexchange
      fj(ix+ncomp,i2) = fexch(ix+ncomp,i2)
    END DO
  END DO
  DO is = 1,nsurf
    DO i2 = 1,ncomp
      fj(is+ncomp+nexchange,i2) = fsurftmp(is+ncomp,i2)
    END DO
    DO is2 = 1,nsurf
      fj(is+ncomp+nexchange,is2+ncomp+nexchange) = fsurftmp(is+ncomp,is2+ncomp)
    END DO
  END DO

! Start coulombic correction section for surface complexation

  IF (sqrt_sion > 0.0) THEN

    CALL SurfaceCharge_init(ncomp,nspec,nsurf,nsurf_sec,npot,portemp,nco)

    DO npt = 1,npot
      nrow = npt + ncomp + nexchange + nsurf
      temp1 = SINH(LogPotential_tmp(npt))
      feq(nrow) = 0.1174d0*sqrt_sion*temp1 - surfcharge_init(kpot(npt))
    END DO

!!  Dependence of electrostatic potential (with mineral surface charge) on other variables (ncomp, nsurf, npot)
    DO npt = 1,npot

      nrow = npt + ncomp + nexchange + nsurf

      k = kpot(npt)

      IF (volin(k,nco) == 0.0d0) THEN
        gramsperL = wtmin(k)*voltemp(k,nco)/(rocond(nco)*volmol(k)*portemp)  ! units of g/L    
      ELSE
        gramsperL = wtmin(k)*volin(k,nco)/(rocond(nco)*volmol(k)*portemp)  ! units of g/L
      END IF
 
      DO i2 = 1,ncomp
        ncol = i2
        sum = 0.0
        DO ns = 1,nsurf_sec
          IF (ksurf(islink(ns)) == kpot(npt)) THEN
            sum = sum - musurf(ns,i2)*zsurf(ns+nsurf)*spsurftmp10(ns+nsurf)
          END IF
        END DO
        fj(nrow,ncol) = sum*faraday/(specific(k,nco)*gramsperL)
      END DO

      DO is2 = 1,nsurf
        ncol = is2 + ncomp + nexchange
        sum = 0.0
        DO ns = 1,nsurf_sec
          IF (ksurf(islink(ns)) == kpot(npt)) THEN
              sum = sum - musurf(ns,is2+ncomp)*zsurf(ns+nsurf)*spsurftmp10(ns+nsurf)
          END IF
        END DO
        fj(nrow,ncol) = sum*faraday/(specific(k,nco)*gramsperL)

        IF (ksurf(is2) == kpot(npt)) THEN
          fj(nrow,ncol) = fj(nrow,ncol)  &
            - zsurf(is2)*spsurftmp10(is2)*faraday/(specific(k,nco)*gramsperL)
        END IF

      END DO

      DO npt2 = 1,npot
        sum = 0.0
        ncol = npt2 + ncomp + nexchange + nsurf
        DO ns = 1,nsurf_sec
          delta_z = zsurf(ns+nsurf) - zsurf(islink(ns))   !! Charge associated with changing from primary to secondary surface complex
          IF ( ksurf(islink(ns)) == kpot(npt2) .AND. kpot(npt) == kpot(npt2) ) THEN
            sum = sum - zsurf(ns+nsurf)*spsurftmp10(ns+nsurf)*delta_z*2.0
          END IF
        END DO
        
        IF (nrow == ncol) THEN
          CheckValue = 0.1174d0*sqrt_sion*COSH(LogPotential_tmp(npt))
          fj(nrow,ncol) = 0.1174d0*sqrt_sion*COSH(LogPotential_tmp(npt)) - sum*faraday/(specific(k,nco)*gramsperL)
        ELSE
          fj(nrow,ncol) = 0.0d0
!!          fj(nrow,ncol) =  -sum*faraday/(specific(k,nco)*gramsperL)
        END IF

      END DO

    END DO

    DO i = 1,ncomp

      IF (equilibrate(i,nco) .AND. itype(i,nco) == 1) THEN
        nrow = i
        DO npt2 = 1,npot
          ncol = npt2+ncomp+nexchange+nsurf
          sum = 0.0
          DO ns = 1,nsurf_sec
            delta_z = zsurf(ns+nsurf) - zsurf(islink(ns))
            IF (ksurf(islink(ns)) == kpot(npt)) THEN
              sum = sum - 2.0d0*delta_z*musurf(ns,i)*spsurftmp10(ns+nsurf)
            END IF
          END DO
          fj(nrow,ncol) = sum 
          fj(nrow,ncol) = 0.0d0    
        END DO
      END IF
    END DO

    DO is = 1,nsurf
      nrow = is + ncomp + nexchange
      DO npt2 = 1,npot
        ncol = npt2+ncomp+nexchange+nsurf
          sum = 0.0
          DO ns = 1,nsurf_sec
            delta_z = zsurf(ns+nsurf) - zsurf(islink(ns))
             IF (ksurf(islink(ns)) == kpot(npt2)) THEN
              sum = sum - 2.0d0*delta_z*musurf(ns,is+ncomp)*spsurftmp10(ns+nsurf)
            END IF
          END DO
          fj(nrow,ncol) = sum  
      END DO
    END DO

!!!    SurfacechargeOriginal = surfcharge_init
!!!    fj_num = 0.0d0
!!!    DO i2 = 1,ncomp
!!!      ncol = i2
!!!      sptmp(i2) = sptmp(i2) + perturb
!!!      CALL surf_init(ncomp,nspec,nsurf,nsurf_sec,nco)
!!!      CALL SurfaceCharge_init(ncomp,nspec,nsurf,nsurf_sec,npot,portemp,nco)
!!!      DO npt = 1,npot
!!!        nrow = npt + ncomp + nexchange + nsurf
!!!        k = kpot(npt)
!!!        deltaSurfaceCharge = surfcharge_init(k)-SurfaceChargeOriginal(k)
!!!        feq_tmp = 0.1174*sqrt_sion*SINH(LogPotential_tmp(npt)) - surfcharge_init(k)
!!!        fj_num(nrow,ncol) = (feq_tmp - feq(nrow))/perturb
!!!        IF (DABS( fj_num(nrow,ncol) ) < 1.0D-25) THEN
!!!           fj_num(nrow,ncol) = 0.0d0
!!!        END IF
!!!      END DO
!!!      sptmp(i2) = sptmp(i2) - perturb

!!!      write(*,*)
!!!      write(*,*) ulab(i2)
!!!      do nrow=42,47
!!!         write(*,*) nrow, fj(nrow,ncol),fj_num(nrow,ncol)
!!!      end do
!!!      write(*,*)

!!!    END DO

!!!    DO is2 = 1,nsurf
!!!      ncol = is2+ncomp+nexchange
!!!      spsurftmp(is2) = spsurftmp(is2) + perturb
!!!      CALL surf_init(ncomp,nspec,nsurf,nsurf_sec,nco)
!!!      CALL SurfaceCharge_init(ncomp,nspec,nsurf,nsurf_sec,npot,portemp,nco)
!!!      DO npt = 1,npot
!!!        nrow = npt + ncomp + nexchange + nsurf
!!!        k = kpot(npt)
!!!        if (k==ksurf(is2)) then
!!!        deltaSurfaceCharge = surfcharge_init(k)-SurfaceChargeOriginal(k)
!!!        feq_tmp = 0.1174*sqrt_sion*SINH(LogPotential_tmp(npt)) - surfcharge_init(k)
!!!        fj_num(nrow,ncol) = (feq_tmp - feq(nrow))/perturb
!!!        else
!!!          fj_num(nrow,ncol) = 0.0d0
!!!        end if
!!!        IF (DABS( fj_num(nrow,ncol) ) < 1.0D-25) THEN
!!!           fj_num(nrow,ncol) = 0.0d0
!!!        END IF
!!!      END DO
!!!
!!!      spsurftmp(is2) = spsurftmp(is2) - perturb

!!!      write(*,*)
!!!     write(*,*) namsurf(is2)
!!!      do nrow=42,47
!!!         write(*,*) nrow,fj(nrow,ncol),fj_num(nrow,ncol)
!!!      end do
!!!      write(*,*)

!!!    END DO

!!!    DO npt2 = 1,npot
!!!      ncol = npt2+ncomp+nexchange+nsurf
!!!
!!!      LogPotential_tmp(npt2) = LogPotential_tmp(npt2) + perturb
!!!
!!!      CALL surf_init(ncomp,nspec,nsurf,nsurf_sec,nco)
!!!      CALL SurfaceCharge_init(ncomp,nspec,nsurf,nsurf_sec,npot,portemp,nco)
!!!      DO npt = 1,npot
!!!        nrow = npt + ncomp + nexchange + nsurf
!!!        if (npt2==npt) then
!!!          k = kpot(npt)
!!!          feq_tmp = 0.1174*sqrt_sion*SINH(LogPotential_tmp(npt)) - surfcharge_init(k)
!!!          fj_num(nrow,ncol) = (feq_tmp - feq(nrow))/perturb
!!!        else
!!!          fj_num(nrow,ncol) = 0.0d0
!!!        end if
!!!      END DO 

!!!      write(*,*)
!!!     write(*,*) ' Potential'
!!!      do nrow=42,47
!!!         write(*,*) nrow,fj(nrow,ncol),fj_num(nrow,ncol)
!!!      end do
!!      read(*,*)
  
!!!      write(*,*)
!!!      write(*,*) ' Dependence on potential ',npt2
!!!      do nrow=38,43
!!!          write(*,*) ' on Potential ',nrow-37, fj(nrow,ncol),fj_num(nrow,ncol)
!!!       end do
!!!       read(*,*)
   
!!!      IF (equilibrate(i,nco) .AND. itype(i,nco) == 1) THEN
!!!      CALL totconc_init(ncomp,nspec,nexchange,nexch_sec,nsurf,nsurf_sec,nco)
!!!      DO i = 1,ncomp
!!!       IF (itype(i,nco) == 1) THEN
!!!          nrow = i
!!!          feq_tmp = stmp(i) - ctot(i,nco)
!!!          fj_num(nrow,ncol) = (feq_tmp - feq(nrow))/perturb
!!!        END IF
!!!!!       fj(nrow,ncol) = fj_num(nrow,ncol)
!!!      END DO
!!!      end if

!!!     CALL surf_init(ncomp,nspec,nsurf,nsurf_sec,nco)
!!!     CALL totsurf_init(ncomp,nsurf,nsurf_sec)
!!!      DO is = 1,nsurf
!!!        nrow = is + ncomp + nexchange
!!!        feq_tmp = ssurftmp(is) - c_surf(is,nco)
!!!        if (kpot(npt2) == ksurf(is)) then
!!!        fj_num(nrow,ncol) = (feq_tmp - feq(nrow))/perturb
!!!        else
!!!         fj_num(nrow,ncol) = 0.0d0
!!!        end if
!!       fj(nrow,ncol) = fj_num(nrow,ncol)
!!!      END DO

!!!      write(*,*)
!!!      write(*,*) ' Dependence on potential ',npt2
!!!      do nrow=1,ncomp
!!!          write(*,*) ' of primary ',nrow, fj(nrow,ncol),fj_num(nrow,ncol)
!!!       end do
!!!      write(*,*)
!!!      write(*,*) ' Dependence on potential ',npt2
!!!      do nrow=ncomp+1,ncomp+nexchange
!!!          write(*,*) ' of exchange ',nrow, fj(nrow,ncol),fj_num(nrow,ncol)
!!!       end do
!!!      write(*,*)
!!!      write(*,*) ' Dependence on potential ',npt2
!!!      do nrow=ncomp+nexchange+1,ncomp+nexchange+nsurf
!!!          write(*,*) ' of surface ',nrow, fj(nrow,ncol),fj_num(nrow,ncol)
!!!       end do
!!!       read(*,*)

!!!      LogPotential_tmp(npt2) = LogPotential_tmp(npt2) - perturb

!!!    END DO

  ELSE

    DO npt = 1,npot
      nrow = npt + ncomp + nexchange + nsurf
      feq(nrow) = LogPotential_tmp(npt)
      fj(nrow,nrow) = 1.0
    END DO 

  END IF

  fxmaxx = 0.0
  DO i = 1,ncomp+nexchange+nsurf
    IF (DABS(feq(i)) > fxmaxx .AND.  ulab(i) /= 'H2O') THEN
      fxmaxx = DABS(feq(i))
    END IF
  END DO
  fxmaxPotential = 0.0
  DO i = ncomp+nexchange+nsurf,ncomp+nexchange+nsurf+npot
    IF (DABS(feq(i)) > fxmaxPotential) THEN
      fxmaxPotential = DABS(feq(i))
    END IF
  END DO

  beta = -feq

!  Solve the set of equations using LU decomposition. 

  CALL dgetrf(neqn,neqn,fj,neqn,indx,info)
  CALL dgetrs(trans,neqn,ione,fj,neqn,indx,beta,neqn,info)

!  Check the error in the concentration corrections.
  
  errx = 0.0D0
  betamax = 0.0D0
  DO i = 1,ncomp+nexchange+nsurf+npot
    IF (DABS(beta(i)) > betamax) THEN
      betamax = DABS(beta(i))
    END IF
    errx = errx + DABS(beta(i))
  END DO
  
!  Damp the Newton step if necessary.
  
  DO i = 1,neqn-npot
    IF (DABS(beta(i)) > corrmax) THEN
      beta(i) = SIGN(corrmax,beta(i))
    ELSE
      CONTINUE
    END IF
    u(i) = u(i) + beta(i)
  END DO

  DO i = 1,npot
    IF (DABS(beta(i+ncomp+nexchange+nsurf)) > 0.9d0) THEN
       beta(i+ncomp+nexchange+nsurf) = SIGN(0.9d0,beta(i+ncomp+nexchange+nsurf))
     ELSE
       CONTINUE
     END IF
    u(i+ncomp+nexchange+nsurf) = u(i+ncomp+nexchange+nsurf) + beta(i+ncomp+nexchange+nsurf)
  END DO

  DO i = 1,ncomp
    sptmp(i) = u(i)
    sptmp10(i) = DEXP(sptmp(i))
  END DO
  DO ix = 1,nexchange
    spextmp(ix) = u(ix+ncomp)
    spextmp10(ix) = DEXP(spextmp(ix))
  END DO
  DO is = 1,nsurf
    spsurftmp(is) = u(is+ncomp+nexchange)
    spsurftmp10(is) = DEXP(spsurftmp(is))
  END DO
  DO npt = 1,npot
    LogPotential_tmp(npt) = u(npt+ncomp+nexchange+nsurf)
  END DO
  
  DO ik = 1,ncomp+nspec
    spbase(ik) = sptmp(ik)/clg
  END DO
  
  atol = tolf

  ChargeOK = .TRUE.

  IF (NeedChargeBalance) THEN
    IF (ChargeBalance) THEN
      ChargeOK = .TRUE.
    ELSE
      ChargeOK = .FALSE.
    END IF
  END IF

   if (ktrial == 4999) THEN
     continue
   end if
   
!!!   write(*,*) ' fxmaxx =',fxmaxx
  
  IF (fxmaxx < atol .AND. DABS(fxmaxPotential) < 1.0E-08 .AND. ktrial > 20 .AND. ChargeOK) THEN
    
    IF (igamma /= 0) THEN
      if (Duan .OR. Duan2006) then
        call gamma_init_co2(ncomp,nspec,tempc,sqrt_sion,pg)
      else
        CALL gamma_init(ncomp,nspec,igamma,tempc,sqrt_sion,sion_tmp)
      end if
    END IF
    
    CALL species_init(ncomp,nspec)
    
    if (Duan .OR. Duan2006) then
      CALL gases_init_co2(ncomp,ngas,tempc,pg,vrInOut)
    else
      CALL gases_init(ncomp,ngas,tempc,nco)
    end if
    
    CALL surf_init(ncomp,nspec,nsurf,nsurf_sec,nco)
    CALL exchange_init(ncomp,nspec,nexchange,nexch_sec,nco)
    CALL totconc_init(ncomp,nspec,nexchange,nexch_sec,nsurf, nsurf_sec,nco)
    CALL totgas_init(ncomp,nspec,ngas)
    CALL totsurf_init(ncomp,nsurf,nsurf_sec)
    
    IF (Duan .OR. Duan2006) THEN
      vrInitial(nco) = vrInOut
    END IF
     
    WRITE(iunit2,*)
    WRITE(iunit2,*) '*********************************************'
    WRITE(iunit2,*) ' ---> GEOCHEMICAL CONDITION: ', dumstring(1:lscond)
    WRITE(iunit2,*)

555 FORMAT(2X, 'Temperature (C)         = ',f10.3)
558 FORMAT(2X, 'Porosity                = ',f10.3)
559 FORMAT(2X, 'Conversion (M->m)       = ',f10.3)
556 FORMAT(2X, 'Liquid Saturation       = ',f10.3)
557 FORMAT(2X, 'Liquid Density (kg/m^3) = ',f10.3)
560 FORMAT(2X, 'Solid Density (kg/m^3)  = ',f10.3)
562 FORMAT(2X, 'Solid Density           = ',f10.3, ' (Assumed--cannot be computed with porosity = 100%)')
561 FORMAT(2X, 'Solid:Solution Ratio    = ',f10.3)
563 FORMAT(2X, 'Activity of water       = ',f10.3)
411 FORMAT(2X, 'Ionic Strength          = ',f10.3)
5022 FORMAT(2X,'Solution pH             = ',f10.3)
5023 FORMAT(2X,'Solution pe             = ',f10.3)
5024 FORMAT(2X,'Solution Eh (volts)     = ',f10.3)

205  FORMAT(2X,'Total Charge            = ',1pe10.3)
207 FORMAT('                 ','          Log',1X,'       Log',1x,  &
'              ',1x,'          ', '       Activity')
204 FORMAT(' Species         ','     Molality',1X,'  Activity',1x,  &
'     Molality ',1x,'    Activity ', '  Coefficient','    Type')
515 FORMAT(a14,3X,1PE12.4,5X,1PE12.4,3x,1PE12.4)

    sum = 0.0D0
    DO ik = 1,ncomp+nspec

        sum = sum + sptmp10(ik)*chg(ik)*chg(ik)

    END DO
    sion_tmp = 0.50D0*sum
    
    gammawatertmp = EXP(gamtmp(ih2o))
    
    WRITE(iunit2,555) tempc
    WRITE(iunit2,558) porcond(nco)
    WRITE(iunit2,556) SaturationCond(nco)
    WRITE(iunit2,557) rocond(nco)
    IF (porcond(nco) == 1.0d0 .AND. SolidDensityFrom(nco) == 2) THEN
      WRITE(iunit2,562) SolidDensity(nco)
    ELSE
      WRITE(iunit2,560) SolidDensity(nco)
    END IF
    WRITE(iunit2,561) SolidSolutionRatio(nco) 
    WRITE(iunit2,563) gammawatertmp
    
    WRITE(iunit2,411) sion_tmp
    
    IF (conversion(nco) /= 1.0) THEN
      WRITE(iunit2,559) conversion(nco)
    END IF
    
    IF (ikph /= 0) THEN
      pHprint = -(sptmp(ikph)+gamtmp(ikph))/clg
      WRITE(iunit2,5022) pHprint
    END IF
    IF (O2Found) THEN     !  Go through the possible redox couples, calculating pe
      IF (npointH2gas == 0) THEN
        WRITE(iunit2,5025)
        5025 FORMAT(2x,'Solution pe and Eh require H2(g) be included')
      ELSE
        peprint = (-0.5*keqgas_tmp(npointH2gas) + 0.25*(sptmp(npointO2aq)+gamtmp(npointO2aq)) -     &  
           pHprint*clg + 0.25*keqgas_tmp(npointO2gas) )/clg
        WRITE(iunit2,5023) peprint
        FaradayKJ = 96.4935             !!  kJ/volt-eq
        FaradayKcal = 23.06               !!  kcal/volt-eq
!!      Ehprint = peprint*clg*rgas*tk/23.06
        Ehprint = peprint*clg*rgas*tk/FaradayKJ
        WRITE(iunit2,5024) Ehprint  
      END IF  
    END IF
    totcharge = 0.0
    DO i = 1,ncomp
      totcharge = totcharge + chg(i)*stmp(i)
    END DO
    WRITE(iunit2,205) totcharge
    TotChargeSave(nco) = totcharge
   
    equilibrate(:,nco) = .FALSE.

    CALL totconc_init(ncomp,nspec,nexchange,nexch_sec,nsurf,nsurf_sec,nco)
    WRITE(iunit2,*)
    WRITE(iunit2,*) 'Total Aqueous Concentrations of Primary Species  '
    WRITE(iunit2,*) '--------------------------------------- '
    WRITE(iunit2,*)
    WRITE(iunit2,*) 'Species                Molality     Constraint    Constraint Phase '

    DO i = 1,ncomp
      IF (itype(i,nco) == 1) THEN
        namtemp = 'Total        '
      ELSE IF (itype(i,nco) == 2) THEN
        namtemp = 'Charge       '
      ELSE IF (itype(i,nco) == 3) THEN
        namtemp = 'Mineral      '
      ELSE IF (itype(i,nco) == 4) THEN
        namtemp = 'Gas          '
      ELSE IF (itype(i,nco) == 7) THEN
        namtemp = 'Activity      '
      ELSE IF (itype(i,nco) == 8) THEN
        namtemp = 'Concentration '
      END IF
      WRITE(iunit2,511) ulab(i),stmp(i),namtemp,ncon(i,nco)
    END DO
    
!!!    write(123,*) PressureTemp,stmp(1)
    write(123,*) tempc,stmp(1)

    namtemp = 'Exchange'
    IF (nexchange > 0) THEN
      WRITE(iunit2,*)
      WRITE(iunit2,*) 'Exchangers        Equiv/kgw        Equiv/g solid  Equiv/m^3 bulk'
      TotalExchangeEquivalents = 0.0d0
      DO ix = 1,nexchange
        IF (SolidSolutionRatio(nco) == 0.0d0) THEN
          WRITE(iunit2,515) namexc(ix),totextmp(ix),0.0d0
        ELSE
          WRITE(iunit2,515) namexc(ix),totextmp(ix),totextmp(ix)/SolidSolutionRatio(nco),totextmp(ix)*porcond(nco)*rocond(nco)
        END IF
        TotalExchangeEquivalents = TotalExchangeEquivalents + totextmp(ix)
      END DO
    END IF

    namtemp = 'Surface Complex'
    IF (nsurf > 0) THEN
      WRITE(iunit2,*)
      WRITE(iunit2,*) 'Surface complex   Sites/kgw        Moles/g solid  Moles/m^3 bulk'
      DO is = 1,nsurf
        IF (SolidSolutionRatio(nco) == 0.0d0) THEN
          WRITE(iunit2,515) namsurf(is),ssurftmp(is),0.0d0
        ELSE
          WRITE(iunit2,515) namsurf(is),ssurftmp(is),ssurftmp(is)/SolidSolutionRatio(nco),ssurftmp(is)*porcond(nco)*rocond(nco)
        END IF
      END DO
    END IF

    IF (nsurf > 0) THEN
      namtemp = 'Surface'
      WRITE(iunit2,*)
      WRITE(iunit2,*) 'Total Concentrations on Surface Hydroxyl Sites '
      WRITE(iunit2,*) '---------------------------------------------- '
      WRITE(iunit2,*)
      WRITE(iunit2,*) 'Primary Species   Moles/kgw       Moles/g solid'
      totSurf_bas = 0.0
      DO i = 1,ncomp  
        DO ns = 1,nsurf_sec
          totSurf_bas(i) = totSurf_bas(i) + musurf(ns,i)*spsurftmp10(ns+nsurf)
        END DO
        IF (totSurf_bas(i) /= 0.0) THEN
          IF (SolidSolutionRatio(nco) == 0.0d0) THEN
            WRITE(iunit2,515) ulab(i),totSurf_bas(i),0.0D0 
          ELSE
            WRITE(iunit2,515) ulab(i),totSurf_bas(i),totSurf_bas(i)/SolidSolutionRatio(nco)
          END IF
        END IF
      END DO

      IF (pest) THEN

!! **********  PEST output  ****************************
!!  Convert to mol/g_solids
!!        mol/g = g/kg * mol/kgw * kgw/m^3 water * m^3 water/m^3 bulk *m^3 bulk/m^3 solid * m^3 solid/kg solid 
!!                = 1000.0* totex_bas * rocond        * porcond            /( (1-porcond)      * solid_density )

        IF (NPestSurface > 0) THEN

          IF (PestSurfaceUnits == 'mol/g') THEN
            WRITE(UnitPestSurface,*) 'Primary Species   Moles/kgw       Moles/g solid'
          ELSE IF (PestSurfaceUnits == 'mmol/g') THEN
            WRITE(UnitPestSurface,*) 'Primary Species   mMoles/kgw      mMoles/g solid '
          ELSE IF (PestSurfaceUnits == 'umol/g') THEN
            WRITE(UnitPestSurface,*) 'Primary Species   uMoles/kgw      uMoles/g solid '
          ELSE IF (PestSurfaceUnits == 'log') THEN
            WRITE(UnitPestSurface,*) 'Primary Species   Log Moles/kgw   Log Moles/g solid'
          ELSE
            WRITE(*,*)
            WRITE(*,*) ' PestSurfaceUnits not recognized in subroutine EQUILIB'
            WRITE(*,*)
            READ(*,*)
            STOP
          END IF
          DO nl = 1,NPestSurface
            CALL GetPrimarySpeciesNumber(ncomp,PestSurfaceList(nl),i)
            IF (totSurf_bas(i) /= 0.0) THEN  
              IF (SolidSolutionRatio(nco) == 0.0d0) THEN 
                MolePerGSolids = 0.0d0
              ELSE
                MolePerGSolids = totSurf_bas(i)/SolidSolutionRatio(nco)
                IF (PestSurfaceUnits == 'mol/g') THEN
                  WRITE(UnitPestSurface,701) ulab(i),totSurf_bas(i),MolePerGSolids
                ELSE IF (PestSurfaceUnits == 'mmol/g') THEN
                  WRITE(UnitPestSurface,701) ulab(i),totSurf_bas(i)*1.0d03,MolePerGSolids*1.0d03
                ELSE IF (PestSurfaceUnits == 'umol/g') THEN
                  WRITE(UnitPestSurface,701) ulab(i),totSurf_bas(i)*1.0d06,MolePerGSolids*1.0d06
                ELSE IF (PestSurfaceUnits == 'log') THEN
                  IF (totSurf_bas(i) <= 0.0d0 .OR. MolePerGSolids <= 0.0d0) THEN
                    WRITE(UnitPestSurface,701) ulab(i),SmallLogNumber,SmallLogNumber
                  ELSE
                    WRITE(UnitPestSurface,701) ulab(i),Log10(totSurf_bas(i)),Log10(MolePerGSolids) 
                  END IF
                ELSE
                  WRITE(*,*)
                  WRITE(*,*) ' Units for PestSurface not recognized in subroutine EQUILIB'
                  WRITE(*,*)
                  READ(*,*)
                  STOP
                END IF
              END IF
            END IF
          END DO
        END IF

      END IF

    END IF

    IF (nexchange > 0) THEN
      namtemp = 'Exchange'
      WRITE(iunit2,*)
      WRITE(iunit2,*) 'Total Concentrations in Exchange Sites '
      WRITE(iunit2,*) '-------------------------------------- '
      WRITE(iunit2,*)
      WRITE(iunit2,*) 'Primary Species   Moles/kgw       Moles/g solid   Equiv/g solid'
      totex_bas = 0.0d0
      totequiv = 0.0d0
      TotalExchangeMoles = 0.0d0
      DO i = 1,ncomp  
        DO nex = 1,nexch_sec
          ix = ixlink(nex)
          totex_bas(i) = totex_bas(i) + muexc(nex,i)*spextmp10(nex+nexchange)
          totequiv(i) = totequiv(i) + muexc(nex,i)*spextmp10(nex+nexchange)*muexc(nex,ix+ncomp)
        END DO
        IF (SolidSolutionRatio(nco) == 0.0d0) THEN
          MolePerGSolids = 0.0D0
          EquivPerGSolids = 0.0D0
        ELSE
          MolePerGSolids = totex_bas(i)/SolidSolutionRatio(nco)
          EquivPerGSolids = totequiv(i)/SolidSolutionRatio(nco)
        END IF
        IF (totex_bas(i) /= 0.0) THEN
          WRITE(iunit2,515) ulab(i),totex_bas(i),MolePerGSolids, EquivPerGSolids
        END IF
        TotalExchangeMoles = TotalExchangeMoles + totex_bas(i)
      END DO

      IF (HanfordStrontium) THEN

        KvNaSr = ( (DEXP( gamtmp(4) )*sptmp10(4) )**2.0d0 * totex_bas(6)/TotalExchangeMoles )/  &
            ( DEXP( gamtmp(6) )*sptmp10(6) * (totex_bas(4)/TotalExchangeMoles)**2.0d0 )

        KvNaCa = ( (DEXP( gamtmp(4) )*sptmp10(4) )**2.0d0 * totex_bas(3)/TotalExchangeMoles )/  &
            (DEXP( gamtmp(3) )*sptmp10(3) * (totex_bas(4)/TotalExchangeMoles)**2.0d0 )

        LogKvNaSr = DLOG10(KvNaSr)
        LogKvNaCa = DLOG10(KvNaCa)
  
        WRITE(88,*) totequiv(4)/TotalExchangeEquivalents, LogKvNaSr , LogKvNaCa
      END IF

      IF (pest) THEN

!! **********  PEST output  ****************************
!!  Convert to mol/g_solids
!!        mol/g = g/kg * mol/kgw * kgw/m^3 water * m^3 water/m^3 bulk *m^3 bulk/m^3 solid * m^3 solid/kg solid 
!!                = 1000.0* totex_bas * rocond        * porcond            /( (1-porcond)      * solid_density )

        IF (NPestExchange > 0) THEN
          IF (PestExchangeUnits == 'mol/g') THEN
            WRITE(UnitPestExchange,*) 'Primary Species   Moles/kgw     Moles/g solid'
          ELSE IF (PestExchangeUnits == 'mmol/g') THEN
            WRITE(UnitPestExchange,*) 'Primary Species   mMoles/kgw    mMoles/g solid '
          ELSE IF (PestExchangeUnits == 'umol/g') THEN
            WRITE(UnitPestExchange,*) 'Primary Species   uMoles/kgw    uMoles/g solid '
          ELSE IF (PestExchangeUnits == 'equiv/g') THEN
            WRITE(UnitPestExchange,*) 'Primary Species   Equiv/kgw     Equiv/g solid'
          ELSE IF (PestExchangeUnits == 'mequiv/g') THEN
            WRITE(UnitPestExchange,*) 'Primary Species   mEquiv/kgw    mEquiv/g solid'
          ELSE IF (PestExchangeUnits == 'uequiv/g') THEN
            WRITE(UnitPestExchange,*) 'Primary Species   uEquiv/kgw    uEquiv/g solid'
          ELSE IF (PestExchangeUnits == 'logequivalents') THEN
            WRITE(UnitPestExchange,*) 'Primary Species   Log Equiv/kgw Log Equiv/g solid'
          ELSE IF (PestExchangeUnits == 'logmoles') THEN
            WRITE(UnitPestExchange,*) 'Primary Species   Log Moles/kgw Log Moles/g solid'
          ELSE
            WRITE(*,*)
            WRITE(*,*) ' PestExchangeUnits not recognized in subroutine EQUILIB'
            WRITE(*,*)
            READ(*,*)
            STOP
          END IF

          DO nl = 1,NPestExchange
            CALL GetPrimarySpeciesNumber(ncomp,PestExchangeList(nl),i)
            IF (totex_bas(i) /= 0.0) THEN  
              IF (SolidSolutionRatio(nco) == 0.0d0) THEN 
                MolePerGSolids = 0.0d0
                EquivPerGSolids = 0.0d0
              ELSE
                MolePerGSolids = totex_bas(i)/SolidSolutionRatio(nco)
                EquivPerGSolids = totequiv(i)/SolidSolutionRatio(nco)
                IF (PestExchangeUnits == 'mol/g') THEN
                  WRITE(UnitPestExchange,701) ulab(i),totEx_bas(i),MolePerGSolids
                ELSE IF (PestExchangeUnits == 'mmol/g') THEN
                  WRITE(UnitPestExchange,701) ulab(i),totEx_bas(i)*1.0D03,MolePerGSolids*1.0d03
                ELSE IF (PestExchangeUnits == 'umol/g') THEN
                  WRITE(UnitPestExchange,701) ulab(i),totEx_bas(i)*1.0D06,MolePerGSolids*1.0d06
                ELSE IF (PestExchangeUnits == 'equiv/g') THEN
                  WRITE(UnitPestExchange,701) ulab(i),totequiv(i),EquivPerGSolids
                ELSE IF (PestExchangeUnits == 'mequiv/g') THEN
                  WRITE(UnitPestExchange,701) ulab(i),totequiv(i)*1.0d03,EquivPerGSolids*1.0d03
                ELSE IF (PestExchangeUnits == 'uequiv/g') THEN
                  WRITE(UnitPestExchange,701) ulab(i),totequiv(i)*1.0D06,EquivPerGSolids*1.0d06
                ELSE IF (PestExchangeUnits == 'logequivalents') THEN
                  IF (totequiv(i) <= 0.0d0 .OR. EquivPerGSolids <= 0.0d0) THEN
                    WRITE(UnitPestExchange,701) ulab(i),SmallLogNumber,SmallLogNumber
                  ELSE
                    WRITE(UnitPestExchange,701) ulab(i),Log10(totequiv(i)),Log10(EquivPerGSolids) 
                  END IF
                ELSE IF (PestExchangeUnits == 'logmoles') THEN
                  IF (totEx_bas(i) <= 0.0d0 .OR. MolePerGSolids <= 0.0d0) THEN
                    WRITE(UnitPestExchange,701) ulab(i),SmallLogNumber,SmallLogNumber
                  ELSE
                    WRITE(UnitPestExchange,701) ulab(i),Log10(totEx_bas(i)),Log10(MolePerGSolids) 
                  END IF
                ELSE
                  WRITE(*,*)
                  WRITE(*,*) ' PestExchangeUnits not recognized in subroutine EQUILIB'
                  WRITE(*,*)
                  READ(*,*)
                  STOP
                END IF
              END IF
            END IF
          END DO
        END IF

      END IF
    END IF

!! **********  End of PEST output  ****************************
   
    111 FORMAT(1pe11.4)

!!   Kd = porcond(nco)*totex_bas(1)/stmp(1)/1.500

!!   write(101,*)  log10(stmp(2)),log10(stmp(1)),Kd


    WRITE(iunit2,*)
    WRITE(iunit2,*) 'Concentrations of Individual Species, Exchangers, and Surface Complexes '
    WRITE(iunit2,*) '----------------------------------------------------------------------- '
    WRITE(iunit2,*)

    WRITE(iunit2,207)
    WRITE(iunit2,204) 
    namtemp = 'Exchange'
    DO nex = 1,nexch_sec
      ix = nex + nexchange
      IF (spextmp10(ix) == 0.0d0) THEN
        spprint = -100.0
        actprint = -100.0
        actprint10 = 0.0
      ELSE
        spprint = DLOG10(spextmp10(ix))
        actprint10 = aexch(nex) 
        actprint = DLOG10(actprint10)
      END IF

      WRITE(iunit2,211) nam_exchsec(nex),spprint,actprint,spextmp10(ix),actprint10,namtemp
    END DO

    namtemp = 'Surface'
    DO is = 1,nsurf
      spprint = DLOG10(spsurftmp10(is))
      IF (c_surf(is,nco) /= 0.0d0) THEN
        actprint = spprint - DLOG10(c_surf(is,nco))
        actprint10 = spsurftmp10(is)/c_surf(is,nco) 
      ELSE
        actprint = 0.0d0
        actprint10 = 0.0d0 
      END IF
      WRITE(iunit2,211) namsurf(is),spprint,actprint,spsurftmp10(is),actprint10,namtemp
    END DO
    DO ns = 1,nsurf_sec
      is = ns + nsurf
      IF (c_surf(islink(ns),nco) /= 0.0d0) THEN
        spprint = DLOG10(spsurftmp10(is))
        actprint = spprint - DLOG10(c_surf(islink(ns),nco))
        actprint10 = musurf(ns,islink(ns)+ncomp)*spsurftmp10(is)/c_surf(islink(ns),nco)
      ELSE
        actprint = -100.0d0
        actprint10 = 0.0d0 
        spprint = -100.0
      END IF
      WRITE(iunit2,211) namsurf_sec(ns),spprint,actprint,spsurftmp10(is),actprint10,namtemp
    END DO

    namtemp = 'Aqueous'
    DO ik = 1,ncomp+nspec
      spbase(ik) = DLOG10(sptmp10(ik))

      actprint = (sptmp(ik)+gamtmp(ik))/clg

      actprint10 = 10**(actprint)
      actcoeffprint = DEXP(gamtmp(ik))
      WRITE(iunit2,202) ulab(ik),spbase(ik),actprint,sptmp10(ik),actprint10,actcoeffprint,namtemp
    END DO  

202 FORMAT(2X,a18,2X,f8.3,3X,f8.3,2X,1PE12.3,2X,1PE12.3,2X,1PE12.3,2x,a8)
211 FORMAT(2X,a18,2X,f8.3,3X,f8.3,2X,1PE12.3,2X,1PE12.3,2X,'            ',2x,a8)

    WRITE(iunit2,*)

    if (Duan .OR. Duan2006) then
      CALL gases_init_co2(ncomp,ngas,tempc,pg,vrInOut)
      GasPressureTotalInit(nco) = pg
      write(iunit2,514)'Total gas pressure:             ',pg
      ln_fco2 = 0.0d0
      IF (Duan) THEN
          CALL fugacity_co2(pg,tk,ln_fco2,vrInOut)  
      ELSE IF (Duan2006) THEN
          CALL fugacity_co24(pg,tk,ln_fco2,vrInOut)  
          CALL fugacity_co2(pg,tk,ln_fco2_check,vrInOut)  
      END IF
      fco2 = DEXP(ln_fco2)
      write(iunit2,514)'Ln CO2 fugacity coefficient: ',ln_fco2
      write(iunit2,514)'CO2 fugacity coefficient:    ',fco2
    END IF

    WRITE(iunit2,*)
    WRITE(iunit2,*) ' ****** Partial pressure of gases (bars) *****'
    WRITE(iunit2,*)
    
    CALL GasPartialPressure_Init(ncomp,ngas,tempc,pg)

    DO i = 1,ngas
      WRITE(iunit2,513)  namg(i),1.0d0*spgastmp10(i)
    END DO

    WRITE(iunit2,*)
    WRITE(iunit2,*) ' ****** Saturation state of minerals (log[Q/K] *****'
    WRITE(iunit2,*)
    
    DO k = 1,nrct
      kk = k + nspec
      sumiap = 0.0D0
      DO i = 1,ncomp

          sumiap = sumiap + mumin(1,k,i)* (sptmp(i)+gamtmp(i))

      END DO
      silnTMP = sumiap - keqmin_tmp(1,k)
      siprnt = silnTMP/clg
      WRITE(iunit2,509) umin(k),siprnt
    END DO
    
    DEALLOCATE(u)
    DEALLOCATE(feq)
    DEALLOCATE(beta)
    DEALLOCATE(spbase)
    DEALLOCATE(ffscale)
    DEALLOCATE(tol)
    DEALLOCATE(totex_bas)
    DEALLOCATE(totSurf_bas)
    DEALLOCATE(totequiv)
    DEALLOCATE(fj)
    DEALLOCATE(indx)
    DEALLOCATE(sind)
    DEALLOCATE(iTotNegative)

    RETURN
    
  END IF
  
  
!*******************************************
  
END DO

IF (bagit) then
  GOTO 212
ENDIF

10    CONTINUE
!***************End of iteration loop*********************

DEALLOCATE(u)
DEALLOCATE(feq)
DEALLOCATE(beta)
DEALLOCATE(spbase)
DEALLOCATE(ffscale)
DEALLOCATE(tol)
DEALLOCATE(totex_bas)
DEALLOCATE(totequiv)
DEALLOCATE(fj)
DEALLOCATE(sind)

3001 IF (NeedChargeBalance) THEN
  IF (.not. ChargeBalance) THEN
    IF (sumchg < 0.0 .AND. ichg == -1) THEN
      WRITE(*,*) '  ----> Trying to charge balance a negatively charged solution with an anion'
      WRITE(*,2050) '  ----> Total charge without balancing anion: ',sumchg
      WRITE(*,*)
      WRITE(*,*) '       ******  TRY A CATION  ********'
      WRITE(*,*)
      READ(*,*)
      STOP
    ELSE IF (sumchg > 0.0 .AND. ichg == 1) THEN
      WRITE(*,*) '  ----> Trying to charge balance a positively charged solution with a cation'
      WRITE(*,2050) '  ----> Total charge without balancing cation: ',sumchg
      WRITE(*,*)
      WRITE(*,*) '       ******  TRY AN ANION  ********'
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
  END IF
END IF

2050 FORMAT(a47,1x,1PE12.4)

WRITE(*,*)
WRITE(*,*) ' -->EXCEEDED MAXIMUM ITERATIONS IN CONDITION: ', dumstring(1:lscond)
WRITE(*,*)
READ(*,*)
STOP

600 FORMAT(2X,f10.2,2X,a15)
201 FORMAT(2X,a18,2X,f8.2)
203 FORMAT(2X,a18)
701 FORMAT(a15,2x,1PE12.4,2x,1PE12.4)


509 FORMAT(2X,a18,2X,f12.4)
513 FORMAT(1X,a18,1X,1PE16.8) !FORMAT(1X,a18,1X,1PE12.4)
514 FORMAT(1X,a28,1X,1PE16.8) 
510 FORMAT(2X,'GEOCHEMICAL CONDITION NUMBER',i3)


512 FORMAT(' Basis species    ','     Molality  ', '  Constraint type')
511 FORMAT(1X,a18,1X,1PE12.5,5X,a13,1X,a18)

542  FORMAT(5X,'Primary species ',a18,1X,'at grid pt. ',i5)
544  FORMAT(5X,'Check geochem. condition # ',i2)
543  FORMAT(5X,'Constraint mineral ',2X,a18)
643  FORMAT(5X,'Constraint gas',2X,a18)

END SUBROUTINE equilib_co2
!****************************************************************
