
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

SUBROUTINE Flash(ncomp,nspec,nrct,ngas,nsurf,igamma,ikph,nco,nexchange,nexch_sec,nsurf_sec,npot,neqn,DensityModule,jx,jy,jz)
USE crunchtype
USE params
USE concentration
USE mineral
USE solver
USE io
USE temperature
USE medium
USE runtime, ONLY: Duan, Duan2006

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


CHARACTER (LEN=mls)                                        :: DensityModule

INTEGER(I4B), INTENT(IN)                                   :: jx
INTEGER(I4B), INTENT(IN)                                   :: jy
INTEGER(I4B), INTENT(IN)                                   :: jz


!  Internal variables and arrays

REAL(DP)                                                   :: pg

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

REAL(DP), PARAMETER                                        :: tolf=1.e-14
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

LOGICAL(LGT)                                               :: bagit
LOGICAL(LGT)                                               :: ChargeBalance
LOGICAL(LGT)                                               :: NeedChargeBalance
LOGICAL(LGT)                                               :: ChargeOK
LOGICAL(LGT)                                               :: InitializeMineralEquilibrium

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

REAL(DP)                                                   :: tempc
REAL(DP)                                                   :: portemp


CHARACTER (LEN=1)                                          :: trans

REAL(DP)                                                   :: ph2o, ln_fco2, fco2, ln_fco2_check


!!  Do we need these (they should be calculated)
!!!spgas10(kk,jx,jy,jz) = spcondgas10(kk,ConditionNumber)
!!!spgas(kk,jx,jy,jz)   = spcondgas(kk,ConditionNumber)
!!!spex10(ix+nexchange,jx,jy,jz) = convert*spcondex10(ix+nexchange,ConditionNumber)  ! Now in eq/m3 por. med.
!!!spsurf10(is,jx,jy,jz) = convert*spcondsurf10(is,ConditionNumber)


sptmp(:) = sp(:,jx,jy,jz)
sptmp10(:) = sp10(:,jx,jy,jz)
spextmp(:) = spex(:,jx,jy,jz)
spsurftmp(:) = spsurf(:,jx,jy,jz)
LogPotential_tmp(:) = LogPotential(:,jx,jy,jz)
portemp = por(jx,jy,jz)
tempc = t(jx,jy,jz)
tk = tempc + 273.15
denmol = 1.e05/(8.314*tk)   ! P/RT = n/V, with pressure converted from bars to Pascals
faraday = 96485.0 
rone = 1.0 

vrInOut = 1.0
trans = 'N'
iTotNegative = 0

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
 

!!!pg = 0.0d0                         !!!co2
pg = GasPressureTotal(jx,jy,jz)

DO i = 1,ncomp

  IF (itype(i,nco) == 1) THEN
    IF (ctot(i,nco) > 0.0) THEN
      tol(i) = ctot(i,nco)
    ELSE
      tol(i) = 1.0
    END IF

  ELSE IF (itype(i,nco) == 4) THEN
    
    pg = pg + ctot(i,nco)   !! Not just CO2, but all of the gases.  These are augmented below by the partial pressures of gases not specified as constraints
    
  ELSE IF (itype(i,nco) == 2) THEN

    NeedChargeBalance = .TRUE.
  ELSE

    CONTINUE
  END IF

END DO

IF (Duan .OR. Duan2006) THEN            !! Duan CO2 routine
 
  ph2o = 0.0d0
  call calc_ph2o(tk,ph2o)
  pg = pg + ph2o

! need to re-calculate Keq for CO2 for this condition (pg) 
  CALL keqcalc2_init_co2(ncomp,nrct,nspec,ngas,nsurf_sec,tempc,pg)
 
END IF

sumchg = 0.0
DO ik = 1,ncomp+nspec
  sumchg = sumchg + sptmp10(ik)*chg(ik)
END DO

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
 
  END IF

  IF (igamma /= 0) THEN
    if (Duan .OR. Duan2006) then
      call gamma_init_co2(ncomp,nspec,tempc,sqrt_sion,pg)
    else
      CALL gamma_init(ncomp,nspec,tempc,sqrt_sion,sion_tmp)
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

  CALL species_init(ncomp,nspec)
  
  if (Duan .OR. Duan2006) then
    CALL gases_init_co2(ncomp,ngas,tempc,pg,vrInOut)
  else
    CALL gases_init(ncomp,ngas,tempc)
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
      u(npt+ncomp+nexchange+nsurf) = LogPotential_tmp(npt)
    END DO

  END IF
  
!  Calculate residuals
  
  ChargeBalance = .TRUE.

  DO  i = 1,ncomp

    IF (itype(i,nco) == 1) THEN

      feq(i) = stmp(i) - ctot(i,nco)

    ELSE IF (itype(i,nco) == 3) THEN

      sumiap = 0.0D0
      DO  i2 = 1,ncomp

          sumiap = sumiap + mumin(1,k,i2)* (sptmp(i2)+gamtmp(i2))

      END DO

      feq(i) = sumiap - keqmin_tmp(1,k)

    ELSE IF (itype(i,nco) == 4) THEN             !!! Gas constraint

      sumiap = 0.0D0
      DO i2 = 1,ncomp

          sumiap = sumiap + mugas(kg,i2)* (sptmp(i2)+gamtmp(i2))
 
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
      feq(i) = (sptmp(i)+gamtmp(i)) - DLOG(guess(i,nco))
    ELSE IF (itype(i,nco) == 8) THEN
      feq(i) = sptmp(i) -  DLOG(guess(i,nco))
    ELSE IF (itype(i,nco) == 2) THEN
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
      CONTINUE
    END IF
    
  END DO

  DO ix = 1,nexchange
    feq(ix+ncomp) =  sumactivity(ix) - 1.0
  END DO
  
  DO is = 1,nsurf

    feq(is+ncomp+nexchange) = ssurftmp(is) - c_surf(is,nco)

  END DO

  
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
          EXIT
        END IF
      END DO

      DO i2 = 1,ncomp
        fj(i,i2) = mumin(1,k,i2)
      END DO

    ELSE IF (itype(i,nco) == 4) THEN

      DO kg = 1,ngas
        IF (namg(kg) == ncon(i,nco)) THEN
          kk = kg + nspec + nrct
          EXIT
        END IF
      END DO

      DO i2 = 1,ncomp
        fj(i,i2) = mugas(kg,i2)
      END DO
      
    ELSE IF (itype(i,nco) == 7) THEN

      DO i2 = 1,ncomp
        IF (i2 == i) THEN
          fj(i,i2) = 1.0d0
        ELSE
          fj(i,i2) = 0.0
        END IF
      END DO
      
    ELSE IF (itype(i,nco) == 8) THEN

      DO i2 = 1,ncomp
        IF (i2 == i) THEN
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
          fj(i,i2) = dpsi(i,i2)
        END DO    
      END IF    
    ELSE
      CONTINUE
    END IF

    DO ix2 = 1,nexchange
      fj(i,ix2+ncomp) = 0.0
    END DO

    DO is2 = 1,nsurf
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

  ELSE

    DO npt = 1,npot
      nrow = npt + ncomp + nexchange + nsurf
      feq(nrow) = LogPotential_tmp(npt)
      fj(nrow,nrow) = 1.0
    END DO 

  END IF

  fxmaxx = 0.0
  DO i = 1,ncomp+nexchange+nsurf
    IF (DABS(feq(i)) > fxmaxx) THEN
      fxmaxx = DABS(feq(i))
    END IF
  END DO
  fxmaxPotential = 0.0
  DO i = ncomp+nexchange+nsurf,ncomp+nexchange+nsurf+npot
    IF (DABS(feq(i)) > fxmaxPotential) THEN
      fxmaxPotential = DABS(feq(i))
    END IF
  END DO

!  Solve the set of equations using LU decomposition. 

  beta = -feq

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

  
  IF (fxmaxx < atol .AND. DABS(fxmaxPotential) < 1.0E-07 .AND. ktrial > 1 .AND. ChargeOK) THEN

    write(*,*) 
    write(*,*) ' Convergence in FLASH'
    write(*,*)

    sp(:,jx,jy,jz) = sptmp(:)
    sp10(:,jx,jy,jz) = sptmp10(:)
    spex(:,jx,jy,jz) = spextmp(:) 
    spsurf(:,jx,jy,jz) = spsurftmp(:)
    LogPotential(:,jx,jy,jz) = LogPotential_tmp(:)


  END IF

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

  RETURN
     
END DO


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


WRITE(*,*)
WRITE(*,*) ' -->EXCEEDED MAXIMUM ITERATIONS IN FLASH ', jx,jy,jz
WRITE(*,*)
READ(*,*)
STOP
 
END SUBROUTINE Flash

!!!    if (Duan .OR. Duan2006) then
!!!      call gamma_init_co2(ncomp,nspec,tempc,sqrt_sion,pg)
!!!    else
!!!      CALL gamma_init(ncomp,nspec,tempc,sqrt_sion)
!!!    end if
!!!    CALL species_init(ncomp,nspec)
    
!!!    if (Duan .OR. Duan2006) then
!!!      CALL gases_init_co2(ncomp,ngas,tempc,pg,vrInOut)
!!!    else
!!!      CALL gases_init(ncomp,ngas,tempc)
!!!    end if
    
!!!    CALL surf_init(ncomp,nspec,nsurf,nsurf_sec,nco)
!!!    CALL exchange_init(ncomp,nspec,nexchange,nexch_sec,nco)
!!!    CALL totconc_init(ncomp,nspec,nexchange,nexch_sec,nsurf, nsurf_sec,nco)
!!!    CALL totgas_init(ncomp,nspec,ngas)
!!!    CALL totsurf_init(ncomp,nsurf,nsurf_sec)
    
!!!    IF (Duan .OR. Duan2006) THEN
!!!      vrInitial(nco) = vrInOut
!!!    END IF
!****************************************************************
