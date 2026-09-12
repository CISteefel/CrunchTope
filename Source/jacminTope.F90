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
    
SUBROUTINE jacmin(ncomp,nspec,nexchange,nsurf,nkin,nrct,jx,jy,jz,time,round)
USE crunchtype
USE runtime
USE params
USE concentration
USE mineral
USE solver
USE medium
USE temperature
USE isotope
USE transport

IMPLICIT NONE

!********  INTERFACE BLOCKS  *************************

interface
  SUBROUTINE ReactionNumerical(ncomp,nkin,nrct,nspec,nexchange,nsurf,jx,jy,jz,rminTMP)
  USE crunchtype
  USE params
  USE concentration
  USE mineral
  USE medium
  USE temperature
  USE runtime
  IMPLICIT NONE
  INTEGER(I4B), INTENT(IN)                                        :: ncomp
  INTEGER(I4B), INTENT(IN)                                        :: nkin
  INTEGER(I4B), INTENT(IN)                                        :: nrct
  INTEGER(I4B), INTENT(IN)                                        :: nspec
  INTEGER(I4B), INTENT(IN)                                        :: nexchange
  INTEGER(I4B), INTENT(IN)                                        :: nsurf
  INTEGER(I4B), INTENT(IN)                                        :: jx
  INTEGER(I4B), INTENT(IN)                                        :: jy
  INTEGER(I4B), INTENT(IN)                                        :: jz
  REAL(DP),DIMENSION(:,:),INTENT(OUT)                             :: rminTMP
  END SUBROUTINE ReactionNumerical
END interface
interface
  SUBROUTINE AffinityNumerical(ncomp,nrct,jx,jy,jz,np,k,sppTMP,AffinityTerm,time)
    USE crunchtype
    USE concentration, ONLY: gamma,sn,sppTMP10,ulab,ikh2o
    IMPLICIT NONE
    INTEGER(I4B), INTENT(IN)                                        :: ncomp
    INTEGER(I4B), INTENT(IN)                                        :: nrct
    INTEGER(I4B), INTENT(IN)                                        :: jx
    INTEGER(I4B), INTENT(IN)                                        :: jy
    INTEGER(I4B), INTENT(IN)                                        :: jz
    INTEGER(I4B), INTENT(IN)                                        :: np
    INTEGER(I4B), INTENT(IN)                                        :: k
    REAL(DP),DIMENSION(:),INTENT(IN)                                :: sppTMP
    REAL(DP), INTENT(OUT)                                           :: AffinityTerm
    REAL(DP), INTENT(IN)                                            :: time
  END SUBROUTINE AffinityNumerical
END interface


!********  END INTERFACE BLOCKS  *************************

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                             :: ncomp
INTEGER(I4B), INTENT(IN)                             :: nspec
INTEGER(I4B), INTENT(IN)                             :: nexchange
INTEGER(I4B), INTENT(IN)                             :: nsurf
INTEGER(I4B), INTENT(IN)                             :: nkin
INTEGER(I4B), INTENT(IN)                             :: nrct
INTEGER(I4B), INTENT(IN)                             :: jx
INTEGER(I4B), INTENT(IN)                             :: jy
INTEGER(I4B), INTENT(IN)                             :: jz
REAL(DP), INTENT(IN)                                 :: time
INTEGER(I4B), INTENT(IN)                             :: round

!  Internal variables and arrays

REAL(DP)                                             :: tk
REAL(DP)                                             :: tkinv
REAL(DP)                                             :: reft
REAL(DP)                                             :: porfactor

REAL(DP)                                             :: AffinityTerm
REAL(DP)                                             :: AffinityJacobianTerm
REAL(DP)                                             :: perturb
REAL(DP)                                             :: termTMP
REAL(DP)                                             :: sign
REAL(DP)                                             :: term1
REAL(DP)                                             :: term2
REAL(DP)                                             :: term3
REAL(DP)                                             :: term2New
REAL(DP)                                             :: snormTMP

INTEGER(I4B)                                         :: k
INTEGER(I4B)                                         :: np
INTEGER(I4B)                                         :: i
INTEGER(I4B)                                         :: i2
INTEGER(I4B)                                         :: is
INTEGER(I4B)                                         :: iss
INTEGER(I4B)                                         :: ix
INTEGER(I4B)                                         :: kk
INTEGER(I4B)                                         :: ksp

REAL(DP)                                             :: RateFactor
REAL(DP)                                             :: RateFactorTMP
REAL(DP)                                             :: CubicTerm
REAL(DP)                                             :: SilogCubic

REAL(DP)                                             :: CorrectedTerm
REAL(DP)                                             :: CheckDiff
!!REAL(DP)                                             :: Kformation
REAL(DP)                                             :: Al
REAL(DP)                                             :: PartialSnorm

REAL(DP)                                             :: tmp1
REAL(DP)                                             :: tmp2
REAL(DP)                                             :: tmp3
REAL(DP)                                             :: tmp4

REAL(DP)                                             :: factor1
REAL(DP)                                             :: factor2
REAL(DP)                                             :: SaturationTerm
REAL(DP), PARAMETER                                  :: tiny=1.D-12               
!!  For numerical derivatives
!!REAL(DP),DIMENSION(nreactmax,nkin)                   :: rminTMP
!!REAL(DP),DIMENSION(nreactmax,nkin)                   :: rminORIG
!!REAL(DP),DIMENSION(ncomp)                            :: jac_RateCheck
!!REAL(DP), DIMENSION(ncomp,nreactmax,nkin)           :: jac_check
REAL(DP), DIMENSION(ncomp)                           :: jac_check

CHARACTER (LEN=5)                                    :: Dumstring

REAL(DP)                                             :: dependTMP
REAL(DP)                                             :: surftmp

INTEGER(I4B)                                                    :: kIsotopologue
INTEGER(I4B)                                                    :: Isotopologue
INTEGER(I4B)                                                    :: kMineralRare
INTEGER(I4B)                                                    :: kMineralCommon
INTEGER(I4B)                                                    :: iPrimaryRare
INTEGER(I4B)                                                    :: iPrimaryCommon

INTEGER(I4B)                                                    :: kIsotopologuePoint

!!!REAL(DP), DIMENSION(NCOMP)                           :: dMoleFraction40
!!!REAL(DP), DIMENSION(NCOMP)                           :: dMoleFraction44
!!!REAL(DP)                                             :: MoleFraction40Tmp
!!!REAL(DP)                                             :: MoleFraction44Tmp
!!!REAL(DP)                                             :: MoleFraction40
!!!REAL(DP)                                             :: MoleFraction44
!!!REAL(DP), DIMENSION(NCOMP)                           :: dMoleFraction32
!!!REAL(DP), DIMENSION(NCOMP)                           :: dMoleFraction34
!!!REAL(DP), DIMENSION(NCOMP)                           :: dMoleFraction32S
!!!REAL(DP), DIMENSION(NCOMP)                           :: dMoleFraction34S
!!!REAL(DP)                                             :: MoleFraction32Tmp
!!!REAL(DP)                                             :: MoleFraction34Tmp
!!!REAL(DP)                                             :: MoleFraction32
!!!REAL(DP)                                             :: MoleFraction34
!!!REAL(DP)                                             :: MoleFraction32S
!!!REAL(DP)                                             :: MoleFraction34S

REAL(DP)                                             :: DenomSquared
REAL(DP)                                             :: Denom
REAL(DP)                                             :: AffinityInitial

!!!REAL(DP)                                             :: MoleFraction44Mineral
!!!REAL(DP)                                             :: MoleFraction40Mineral
!!!!!!REAL(DP)                                             :: CurrentMoleFraction40
!!!REAL(DP)                                             :: CurrentMoleFraction44
REAL(DP)                                             :: denominator
LOGICAL(LGT)                                         :: SetToAqueousMoleFraction

!!!REAL(DP)                                             :: CalciumCarbonateRatioEffect
!!!REAL(DP)                                             :: RatioToPower
!!!REAL(DP), DIMENSION(NCOMP)                           :: DerivativeActivityRatio

REAL(DP)                                             ::  termHyperbolic

!!!REAL(DP), DIMENSION(ncomp,15)                           :: dMoleFractionAqueousCommon
!!!REAL(DP), DIMENSION(ncomp,15)                           :: dMoleFractionAqueousRare

REAL(DP)                                                        :: MoleFractionMineral

! biomass
INTEGER(I4B)                                                    :: ib, jj
! biomass end

!! Nucleation
!!!REAL(DP)                                                        :: SigmaNucleation
!!!REAL(DP)                                                        :: Bnucleation
!!!REAL(DP)                                                        :: Anucleation
REAL(DP)                                                        :: AffinityTMP
REAL(DP)                                                        :: sumiap
REAL(DP)                                                        :: silnTMP

REAL(DP)                                                        :: BtimesSigma

REAL(DP)                                                        :: NucleationExponentialTerm
REAL(DP)                                                        :: v_min
REAL(DP)                                                        :: checkNucleation
REAL(DP), PARAMETER                                             :: BoltzmannTerm=1.3806E-20
REAL(DP)                                                        :: NucleationTerm
REAL(DP)                                                        :: testSigma

REAL(DP)                                                        :: liqsat_fac, sat

REAL(DP)                                                        :: lnActivity
CHARACTER (LEN=3)                                               :: ulabPrint

!!!NoFractionationDissolution = .false.

SetToAqueousMoleFraction = .FALSE.

IF (JacobianNumerical) THEN
  
  write(*,*)
  write(*,*) ' Numerical Jacobian associated with Hellmann rate law no longer supported'
  write(*,*)
  stop

END IF

tk = t(jx,jy,jz) + 273.15D0
tkinv = 1.0d0/tk
reft = 1.0d0/298.15d0
porfactor = (por(jx,jy,jz)/porin(jx,jy,jz))**(0.666666D0)
!!porfactor = 1.0d0
!!porfactor = por(jx,jy,jz)
IF (DampRateInLowPorosity .AND. por(jx,jy,jz) < 0.001) THEN
  porfactor = porfactor*PorosityDamp
END IF

porfactor = 1.0d0

IF (nIsotopeMineral > 0) THEN
  dMoleFractionAqueousCommon = 0.0d0
  dMoleFractionAqueousRare   = 0.0d0
END IF

DO Isotopologue = 1,nIsotopePrimary
  iPrimaryRare = isotopeRare(Isotopologue)
  iPrimaryCommon = isotopeCommon(Isotopologue)
  Denom = ( sp10(iPrimaryRare,jx,jy,jz) + sp10(iPrimaryCommon,jx,jy,jz) )
  DenomSquared = Denom*Denom
  dMoleFractionAqueousRare(iPrimaryRare,isotopologue)   =        &
         sp10(iPrimaryRare,jx,jy,jz)*   (-sp10(iPrimaryRare,jx,jy,jz)/DenomSquared + 1.0d0/Denom)
  dMoleFractionAqueousCommon(iPrimaryCommon,isotopologue) =    &
         sp10(iPrimaryCommon,jx,jy,jz)* (-sp10(iPrimaryCommon,jx,jy,jz)/DenomSquared + 1.0d0/Denom)
  dMoleFractionAqueousRare(iPrimaryCommon,isotopologue)   =  sp10(iPrimaryCommon,jx,jy,jz)*   (-sp10(iPrimaryRare,jx,jy,jz)/DenomSquared)
  dMoleFractionAqueousCommon(iPrimaryRare,isotopologue)   =  sp10(iPrimaryRare,jx,jy,jz)*   (-sp10(iPrimaryCommon,jx,jy,jz)/DenomSquared)

END DO

jac_rmin = 0.0d0

DO k = 1,nkin
  
  IF (nisotopemineral > 0) THEN
    IF (kPointerIsotope(k)/= 0) THEN
      kIsotopologue  = kPointerIsotope(k)
      kMineralRare   = kIsotopeRare(kIsotopologue)
      KMineralCommon = kIsotopeCommon(kIsotopologue)
    END IF
  END IF

  checkSaturation = .FALSE.
  UseDissolutionOnly = .FALSE.

!!  IF (imintype(1,k) == 5 .AND. imintype(2,k) == 4 .AND. iterat >= 1) THEN
!!    checkSaturation = .TRUE.
!!  END IF
!!  IF (DABS(silog(1,k)) < 0.00001 .AND. checkSaturation) THEN
!!    UseDissolutionOnly = .TRUE.
!!  END IF


!!  Ivolume(k) == 1 means rate is set to just deplete current mineral volume fraction (so no Jacobian needed)
!!  Imintype(1,k) ==6 .OR. 7 means forwward and reverse irreversible, with surface area set constant to 1.0

  IF (ivolume(k) == 1 .OR. imintype(1,k) == 6 .OR. imintype(1,k) == 7) THEN

    jac_rmin(:,:,k) = 0.0d0
    
  ELSE   !! Case where Jacobian of mineral reactions are calculated (block above skips all computation of Jacobian)
    
    jac_pre = 0.0d0
    
    DO np = 1,nreactmin(k)
    
      jac_sat = 0.0d0

!***********  What surface area to use ************************
      
!!  Taken from reaction.F90
      
! ************* Saturation state dependence of rate ********
   
      IF (si(np,k) >= 1.0d0) THEN
        sign = 1.0d0
      ELSE
        sign = -1.0d0
      END IF

!!    If .NOT. biomass
      IF (imintype(np,k) /= 8) THEN

        IF (AffinityDepend1(np,k) == 1.0D0) THEN
          term1 = sign*DABS(snorm(np,k) - 1.0D0)
        ELSE
!!!       Should be the general case, with nonlinear exponents on "snorm" term possible (calculated in reactionTope)
          term1 = sign*DABS(snorm(np,k) - 1.0D0)**(AffinityDepend1(np,k))
        END IF

        IF (imintype(np,k) == 5 .OR. UseDissolutionOnly) THEN                          !! Dissolution only

          AffinityTerm = MIN(0.0d0,term1)
          IF (AffinityTerm == 0.0d0) THEN
              AffinityTerm = 1.0d-12
          END IF

!!          END IF
        ELSE IF (imintype(np,k) == 4) THEN                               !! Precipitation only

          AffinityTerm = MAX(0.0d0,term1)
          IF (AffinityTerm == 0.0d0) THEN
              AffinityTerm = 1.0d-12
          END IF
          
        ELSE IF (imintype(np,k) == 10) THEN                            !! Nucleation case
          
          v_min = 4.6629E-29
        
         IF (si(np,k) > 1.0) THEN
            
           checkNucleation = 16.0/3.0 * 3.1415 * v_min**2/(BoltzmannTerm**3 * Tk**3) *   &
               SigmaNucleation(np,k)**3 / siln(np,k)**2
          
            BtimesSigma = Bnucleation(np,k)*SigmaNucleation(np,k)*SigmaNucleation(np,k)*SigmaNucleation(np,k)
            NucleationTerm = checkNucleation
            NucleationExponentialTerm = DEXP(-checkNucleation)
          
            AffinityTerm = Azero25C(np,k)*NucleationExponentialTerm 
          
          ELSE
            AffinityTerm = 0.0d0
          END IF
        ELSE
          AffinityTerm = term1
        END IF

        IF (imintype(np,k) == 3 .OR. imintype(np,k) == 2) THEN           !! Monod or irreversible
          AffinityTerm = sign
        END IF

        IF (imintype(np,k) /= 3 .AND. imintype(np,k) /= 2) THEN         !! Everything but Monod and irreversible
            
!!        First do nucleation case (imintype = 10)
            
          IF (imintype(np,k) == 10) THEN
              
            IF (si(np,k) > 1.0d0) THEN
!!!  **********  Numerical derivative **********************
            perturb = 1.0D-09
            sppTMP(:)   = sp(:,jx,jy,jz)
            sppTMP10(:) = sp10(:,jx,jy,jz)

            AffinityInitial = AffinityTerm
            DO i = 1,ncomp
              sppTMP(i) = sppTMP(i) + perturb
              sppTMP10(i) = DEXP(sppTMP(i))
              
              sumiap = 0.0D0
              DO i2 = 1,ncomp
                
                ulabPrint = ulab(i)
                IF (ulabPrint(1:3) == 'H2O' .or. ulabPrint(1:3) == 'HHO') THEN
                  lnActivity = lngamma(i,jx,jy,jz)
                ELSE
                  lnActivity = sp(i,jx,jy,jz) + lngamma(i,jx,jy,jz)
                END IF

                sumiap = sumiap + mumin(1,k,i2)*lnActivity

              END DO
              
              silnTMP = (sumiap - keqmin(1,k,jx,jy,jz))
              
              checkNucleation = 16.0/3.0 * 3.1415 * v_min**2/(BoltzmannTerm**3 * Tk**3) *   &
                  SigmaNucleation(np,k)**3 / silnTMP**2
          
              NucleationTerm = checkNucleation
              NucleationExponentialTerm = DEXP(-checkNucleation)
              AffinityTMP = Azero25C(np,k)*NucleationExponentialTerm
            
              jac_sat(i) = (AffinityTMP - AffinityInitial)/perturb
!!!              jac_sat(i) = jac_check(i)
              sppTMP(i) = sp(i,jx,jy,jz) 
              sppTMP10(i) = sp10(i,jx,jy,jz)
            END DO
!!!  *******************************************************
            ELSE
                jac_sat = 0.0d0
            END IF
            
!!        Then strictly first order TST type rate laws:
          ELSE IF (AffinityDepend1(np,k) == 1.0d0 .AND. AffinityDepend2(np,k) == 1.0d0 .AND. AffinityDepend3(np,k) == 1.0d0) THEN

              IF (IsotopeMineralRare(k)) THEN

                jac_sat = 0.0d0

                DO i = 1,ncomp

                  jac_sat(i) = mumin(np,k,i)*si(np,k)  
            
                END DO

                kIsotopologue  = kPointerIsotope(k)
                kMineralRare   = kIsotopeRare(kIsotopologue)
                KMineralCommon = kIsotopeCommon(kIsotopologue)
                isotopologue   = PointerToPrimaryIsotope(kIsotopologue)
                iPrimaryCommon = isotopeCommon(Isotopologue)
                iPrimaryRare   = isotopeRare(Isotopologue)

                IF (MineralAssociate(kMineralCommon)) THEN
                  kIsotopologuePoint = kPointerIsotope( MineralID(kMineralCommon) )
                ELSE
                  kIsotopologuePoint = kIsotopologue
                END IF

                IF (isotopeBackReactionOption(kIsotopologuePoint) == 'none' .OR. UseAqueousMoleFraction(kIsotopologue)) THEN


                    jac_sat(iPrimaryRare)   = jac_sat(iPrimaryRare) + si(np,k) *                                                       &
                                              (-mumin(1,kMineralCommon,iPrimaryCommon))/MoleFractionAqueousRare(Isotopologue)          &
                                              * dMoleFractionAqueousRare(iPrimaryRare,isotopologue)

                    jac_sat(iPrimaryCommon) = si(np,k) *                                                                               &
                                              (-mumin(1,kMineralCommon,iPrimaryCommon))/MoleFractionAqueousRare(Isotopologue)          &
                                              * dMoleFractionAqueousRare(iPrimaryCommon,isotopologue)

                ELSE 

                  CONTINUE

                END IF

                IF (imintype(np,k) == 4 .AND. silog(np,k) < 0.0d0) THEN    !! Precipitation only
                  jac_sat = 0.0d0
                END IF

                IF (imintype(np,k) == 5 .AND. silog(np,k) > 0.0d0) THEN   !! Dissolution only
                  jac_sat = 0.0d0
                END IF

              ELSE IF (IsotopeMineralCommon(k)) THEN

                jac_sat = 0.0d0

                DO i = 1,ncomp

                  jac_sat(i) = mumin(np,k,i)*si(np,k)  
            
                END DO

                kIsotopologue  = kPointerIsotope(k)
                kMineralRare   = kIsotopeRare(kIsotopologue)
                KMineralCommon = kIsotopeCommon(kIsotopologue)
                isotopologue   = PointerToPrimaryIsotope(kIsotopologue)
                iPrimaryCommon = isotopeCommon(Isotopologue)
                iPrimaryRare   = isotopeRare(Isotopologue)

                IF (MineralAssociate(kMineralCommon)) THEN
                  kIsotopologuePoint = kPointerIsotope( MineralID(kMineralCommon) )
                ELSE
                  kIsotopologuePoint = kIsotopologue
                END IF

                IF (isotopeBackReactionOption(kIsotopologuePoint) == 'none' .OR. UseAqueousMoleFraction(kIsotopologue)) THEN


                    jac_sat(iPrimaryCommon) = jac_sat(iPrimaryCommon) + si(np,k) *                                                &
                                              (-mumin(1,kMineralCommon,iPrimaryCommon))/MoleFractionAqueousCommon(Isotopologue)   &
                                                * dMoleFractionAqueousCommon(iPrimaryCommon,isotopologue)
                    jac_sat(iPrimaryRare)   = si(np,k) *                                                                          &
                                              (-mumin(1,kMineralCommon,iPrimaryCommon))/MoleFractionAqueousCommon(Isotopologue)   &
                                                * dMoleFractionAqueousCommon(iPrimaryRare,isotopologue)
                ELSE

                  CONTINUE

                END IF

                IF (imintype(np,k) == 4 .AND. silog(np,k) < 0.0d0) THEN    !! Precipitation only
                  jac_sat = 0.0d0
                END IF

                IF (imintype(np,k) == 5 .AND. silog(np,k) > 0.0d0) THEN   !! Dissolution only
                  jac_sat = 0.0d0
                END IF

              ELSE
  
                DO i = 1,ncomp

                  jac_sat(i) = mumin(np,k,i)*si(np,k)  

                END DO

                IF (imintype(np,k) == 4 .AND. silog(np,k) < 0.0d0) THEN    !! Precipitation only
                  jac_sat = 0.0d0
                END IF

                IF (imintype(np,k) == 5 .AND. silog(np,k) > 0.0d0) THEN   !! Dissolution only
                  jac_sat = 0.0d0
                END IF

            END IF

          ELSE IF ( AffinityDepend1(np,k) /= 1.0d0 .AND. AffinityDepend2(np,k) == 1.0d0 .AND. AffinityDepend3(np,k) == 1.0d0) THEN
              
            SaturationTerm = (si(np,k) - 1.0d0)
            IF (saturationterm == 0.0d0) THEN
              term3 = 0.0d0
            ELSE
              term3 = (DABS(SaturationTerm))**(AffinityDepend1(np,k)-1.0d0)
            END IF
                        
            DO i = 1,ncomp      
  
              jac_sat(i) = mumin(np,k,i)*si(np,k)* AffinityDepend1(np,k)*AffinityDepend2(np,k) * term3

              IF (imintype(np,k) == 4 .AND. silog(np,k) < 0.0d0) THEN   !! Precipitation only
                jac_sat(i) = 0.0d0
              END IF

              IF (imintype(np,k) == 5 .AND. silog(np,k) > 0.0d0) THEN   !! Dissolution only
                  jac_sat(i) = 0.0d0
              END IF

            END DO         

          ELSE IF ( AffinityDepend1(np,k) == 1.0d0 .AND. AffinityDepend2(np,k) /= 1.0d0 .AND. AffinityDepend3(np,k) == 1.0d0) THEN   

            SaturationTerm = ( si(np,k)**AffinityDepend2(np,k) - 1.0d0 )
            IF (saturationterm == 0.0d0) THEN
              term3 = 0.0d0
            ELSE
              term3 = DABS(SaturationTerm)
            END IF
                        
            DO i = 1,ncomp      
  
              jac_sat(i) = mumin(np,k,i)*si(np,k)* AffinityDepend1(np,k)*AffinityDepend2(np,k)*       &
                           ( si(np,k)**(AffinityDepend2(np,k)-1.0d0) ) * term3

              IF (imintype(np,k) == 4 .AND. silog(np,k) < 0.0d0) THEN   !! Precipitation only
                jac_sat(i) = 0.0d0
              END IF

              IF (imintype(np,k) == 5 .AND. silog(np,k) > 0.0d0) THEN   !! Dissolution only
                jac_sat(i) = 0.0d0
              END IF

            END DO     
  
          ELSE IF ( (AffinityDepend1(np,k) /= 1.0d0 .OR. AffinityDepend2(np,k) /= 1.0d0) .AND. AffinityDepend3(np,k) == 1.0d0) THEN

            SaturationTerm = ( si(np,k)**AffinityDepend2(np,k) - 1.0d0 )
            IF (saturationterm == 0.0d0) THEN
              term3 = 0.0d0
            ELSE
              term3 = (DABS(SaturationTerm))**(AffinityDepend1(np,k)-1.0d0)

            END IF
                        
            DO i = 1,ncomp      
  
              jac_sat(i) = mumin(np,k,i)*si(np,k)* AffinityDepend1(np,k)*AffinityDepend2(np,k)*       &
                           ( si(np,k)**(AffinityDepend2(np,k)-1.0d0) ) * term3

              IF (imintype(np,k) == 4 .AND. silog(np,k) < 0.0d0) THEN   !! Precipitation only
                jac_sat(i) = 0.0d0
              END IF

              IF (imintype(np,k) == 5 .AND. silog(np,k) > 0.0d0) THEN   !! Dissolution only
                  jac_sat(i) = 0.0d0
              END IF

            END DO    
              
          ELSE          !!  Full numerical derivative
           
!!  Numerical derivative for full Burch-Hellmann rate law (all exponents /= 1.0)
            IF (.NOT. JacobianNumerical) THEN
              write(*,*) ' Logical JacobianNumerical not turned to .TRUE., despite non-unitary exponents in mineral rate law'
              write(*,*)
              read(*,*)
              STOP
            END IF
            write(*,*) 
            write(*,*) ' Hellmann rate law no longer supported'
            write(*,*)
            stop
            
          ENDIF

        END IF      !!! END of imintype(np,k) /= 2 or 3

      END IF        !!! END of imintype(np,k) /= 8

      IF (imintype(np,k) == 8) THEN                   !! Monod, with thermodynamic factor, F_T

        jj = p_cat_min(k) !! moved here from below

        IF( si(np,k) > 1.0d0) then
          IF (direction_min(np,jj) == 1) THEN
            IF (chi_min(np,jj) == 1.0d0) THEN
              snorm(np,k) = 1.0d0/si(np,k)
            ELSE
!!NOTE:  This is not going to work below (need Chain Rule derivative)
              snorm(np,k) = (1.0d0/si(np,k))**(1.0d0/chi_min(np,jj))
              WRITE(*,*) ' Not set up for chi .NE. 1.0 yet'
              STOP
            END IF
          ELSE
            snorm(np,k) = 1.0d0
          END IF
        ELSE
          IF (direction_min(np,jj) == -1) THEN
            IF (chi_min(np,jj) == 1.0d0) THEN
              snorm(np,k) = si(np,k)
            ELSE
!!NOTE:  This is not going to work below (need Chain Rule derivative)
              snorm(np,k) = si(np,k)**(1.0d0/chi_min(np,jj))
              WRITE(*,*) ' Not set up for chi .NE. 1.0 yet'
              STOP
            END IF
          ELSE
            snorm(np,k) = 1.0d0
          END IF
        END IF

!! NOTE:  This is assuming linear dependence on Delta G for now in the case of Monod/Biomass (with Ft)
        term1 = sign*DABS(snorm(np,k) - 1.0D0)

        IF (direction_min(np,jj) == 1) THEN                               !! Precipitation only
          AffinityTerm = MAX(0.0d0,term1)
        ELSE                         !! Dissolution only
          AffinityTerm = MIN(0.0d0,term1)
        END IF

!! NOTE: What follows assumes a linear dependence on saturation state for MonodBiomass reaction
        DO i = 1,ncomp        
          jac_sat(i) = muminTMP(np,jj,i)*snorm(np,k)
        END DO

      END IF
!   ! end biomass

     
!******** Far from equilibrium dependence of rate *********
 
      IF (imintype(np,k) == 1 .OR. imintype(np,k) == 3 .OR. imintype(np,k) == 4 .OR. imintype(np,k) == 10) THEN  !! TST, irreversible, or ppt only
!!      +++++++++++++++++++++++++++++++++++++++++IF (imintype(np,k) == 1 .OR. imintype(np,k) == 3 .OR. imintype(np,k) == 4) THEN  !! TST, irreversible, or ppt only

          DO kk = 1,ndepend(np,k)
            i = idepend(kk,np,k)
            IF (itot_min(kk,np,k) == 1) THEN    ! Dependence on total concentration
              IF (gimrt) THEN
                DO i2 = 1,ncomp
                  jac_pre(i2,np) = jac_pre(i,np) + fjac(i2,i,jx,jy,jz)*depend(kk,np,k)* pre_rmin(np,k)/s(i,jx,jy,jz)
                END DO
              ELSE
                DO i2 = 1,ncomp
                  jac_pre(i2,np) = jac_pre(i2,np) + fjac_loc(i2,i)*depend(kk,np,k)* pre_rmin(np,k)/s(i,jx,jy,jz)
                END DO
              END IF
            ELSE
              IF (i <= ncomp) THEN   !  Primary species dependence
                jac_pre(i,np) = jac_pre(i,np) + pre_rmin(np,k)*depend(kk,np,k)
              ELSE                   !  Secondary species dependence
                ksp = i - ncomp
                DO i2 = 1,ncomp
                  jac_pre(i2,np) = jac_pre(i2,np) + pre_rmin(np,k)*depend(kk,np,k)*  &
                    muaq(ksp,i2)
                END DO
              END IF
            END IF
          END DO     !  End of of loop through aqueous species dependencies

          DO kk = 1,ndependex(np,k)
            ix = ixdepend(kk,np,k)
            DO i2 = 1,ncomp+nexchange
              jac_pre(i2,np) = jac_pre(i2,np) + pre_rmin(np,k)*dependex(kk,np,k)*  &
                muexc(ix,i2)
            END DO
          END DO     !  End of of loop through exchange species dependencies

          DO kk = 1,ndependsurf(np,k)
            is = isdepend(kk,np,k)
            dependTMP = dependsurf(kk,np,k)
            IF (is <= nsurf) THEN    !!  Primary surface complex dependence
              jac_pre(is+ncomp+nexchange,np) = jac_pre(is+ncomp+nexchange,np) + pre_rmin(np,k)*dependTMP
            ELSE 
              iss = is - nsurf
              DO i2 = 1,ncomp+nsurf
                jac_pre(i2,np) = jac_pre(i2,np) + pre_rmin(np,k)*dependTMP*musurf(iss,i2)    
              END DO
            END IF
          END DO    
          
         IF (HyperbolicInhibition(np,k)) THEN
           Al = sp10(HyperbolicInhibitionPointer(np,k),jx,jy,jz)
  
           denominator = Kformation(np,k) + Al**HyperbolicInhibitionDepend(np,k)
           termHyperbolic = -Al*Kformation(np,k)*HyperbolicInhibitionDepend(np,k)* Al**(HyperbolicInhibitionDepend(np,k)-1.0d0)   &
               /(denominator*denominator)

!!           termHyperbolic = -Kformation(np,k)*HyperbolicInhibitionDepend(np,k)* Al*Al**(1.0-HyperbolicInhibitionDepend(np,k))   &
!!               /(denominator*denominator)
           
           
!!            jac_pre(ikAl,np) = -0.33333333d0*Al*Kformation/(Kformation + Al**0.333333333333)**2.0d0 *Al**0.666666666666
          
           jac_pre(HyperbolicInhibitionPointer(np,k),np) = jac_pre(HyperbolicInhibitionPointer(np,k),np) + termHyperbolic
        END IF
        
! REDIMENSION JAC_PRE(1:ncomp+nexchange+nsurf,np)
        
      ELSE IF (imintype(np,k) == 2) THEN   !  Monod rate law
        
        DO kk = 1,nmonod(np,k)
          IF (imonod(kk,np,k) == 0 .AND. kmonod(kk,np,k) == 1) THEN                !! Mineral Monod reaction depends on its own concentration
            CONTINUE                                                       !! No dependence on this term, since mineral concentration is not a primary unknown
          ELSE
            i = imonod(kk,np,k)
            IF (itot_monod(kk,np,k) == 1) THEN    ! Dependence on total concentration
              IF (gimrt) THEN
                DO i2 = 1,ncomp
                  jac_pre(i2,np) = jac_pre(i2,np) +  &
                    fjac(i2,i,jx,jy,jz)*pre_rmin(np,k)*   &
                       ( 1.0d0/s(i,jx,jy,jz) - 1.0d0/(s(i,jx,jy,jz)+halfsat(kk,np,k))  )
                END DO
              ELSE
                DO i2 = 1,ncomp
                  jac_pre(i2,np) = jac_pre(i2,np) +  &
                    fjac_loc(i2,i)*pre_rmin(np,k)*   &
                       ( 1.0d0/s(i,jx,jy,jz) - 1.0d0/(s(i,jx,jy,jz)+halfsat(kk,np,k))  )
                END DO
              END IF
            ELSE
              IF (i <= ncomp) THEN
                jac_pre(i,np) = jac_pre(i,np) +  &
                  pre_rmin(np,k)*( 1.0d0 - sp10(i,jx,jy,jz)/(sp10(i,jx,jy,jz)+halfsat(kk,np,k)) )
              ELSE
                ksp = i - ncomp
                DO i2 = 1,ncomp
                  jac_pre(i2,np) = jac_pre(i2,np) +  &
                    pre_rmin(np,k)* muaq(ksp,i2)*(1.0d0-sp10(i,jx,jy,jz)/(sp10(i,jx,jy,jz)+halfsat(kk,np,k)) )
                END DO
              END IF
            END IF
          END IF
        END DO

! biomass
      ELSE IF (imintype(np,k) == 8) THEN   !  Monod rate law with thermodynamic factor F_T
        
        DO kk = 1,nmonod(np,k)
          IF (imonod(kk,np,k) == 0 .AND. kmonod(kk,np,k) == 1) THEN     !! Mineral Monod reaction depends on its own concentration
            CONTINUE                                                    !! No dependence on this term, since mineral concentration is not a primary unknown
          ELSE
            i = imonod(kk,np,k)
            IF (itot_monod(kk,np,k) == 1) THEN    ! Dependence on total concentration
              IF (gimrt) then
                DO i2 = 1,ncomp
                  jac_pre(i2,np) = jac_pre(i2,np) +  &
                    fjac(i2,i,jx,jy,jz)*pre_rmin(np,k)*   &
                       ( 1.0d0/s(i,jx,jy,jz) - 1.0d0/(s(i,jx,jy,jz)+halfsat(kk,np,k))  )
                END DO
              ELSE
                DO i2 = 1,ncomp
                  jac_pre(i2,np) = jac_pre(i2,np) +  &
                    fjac_loc(i2,i)*pre_rmin(np,k)*   &
                       ( 1.0d0/s(i,jx,jy,jz) - 1.0d0/(s(i,jx,jy,jz)+halfsat(kk,np,k))  )
                END DO
              END IF
            ELSE
              IF (i <= ncomp) THEN
                jac_pre(i,np) = jac_pre(i,np) +  &
                  pre_rmin(np,k)*( 1.0d0 - sp10(i,jx,jy,jz)/(sp10(i,jx,jy,jz)+halfsat(kk,np,k)) )
              ELSE
                ksp = i - ncomp
                DO i2 = 1,ncomp
                  jac_pre(i2,np) = jac_pre(i2,np) +  &
                    pre_rmin(np,k)* muaq(ksp,i2)*(1.0d0-sp10(i,jx,jy,jz)/(sp10(i,jx,jy,jz)+halfsat(kk,np,k)) )
                END DO
              END IF
            END IF
          END IF
          
        END DO
!!      biomass end

      END IF
     
!!    Now the chain rule ---->
  
!!     DO i = 1,ncomp+nexchange+nsurf
     
      IF (imintype(np,k) == 8) THEN
          
!!      pointer to biomass for current reaction
        jj = p_cat_min(k)
        ib = ibiomass_min(jj)
        IF (UseMetabolicLagMineral(np,jj)) THEN
            
          DO i = 1,ncomp
            jac_rmin(i,np,k) =  MetabolicLagMineral(np,jj,jx,jy,jz)*surf(np,k)*actenergy(np,k)*rate0(np,k)*volfx(ib,jx,jy,jz)* &
                          ( pre_rmin(np,k)*jac_sat(i) + jac_pre(i,np)*AffinityTerm  )
          END DO
          
        ELSE
            
          DO i = 1,ncomp
            jac_rmin(i,np,k) =  surf(np,k)*actenergy(np,k)*rate0(np,k)*volfx(ib,jx,jy,jz)* &
                          ( pre_rmin(np,k)*jac_sat(i) + jac_pre(i,np)*AffinityTerm  )
          END DO
          
        END IF
!!    biomass end     
        
      ELSE

        IF (nIsotopeMineral > 0) THEN

          IF (IsotopeMineralCommon(k)) THEN
            kIsotopologue = kPointerIsotope(k)

            IF (MineralAssociate(kMineralCommon)) THEN
              kIsotopologuePoint = kPointerIsotope( MineralID(kMineralCommon) )
            ELSE
              kIsotopologuePoint = kIsotopologue
            END IF

            IF (isotopeBackReactionOption(kIsotopologuePoint) == 'none' .OR. UseAqueousMoleFraction(kIsotopologue) ) THEN
              isotopologue = PointerToPrimaryIsotope(kIsotopologue)
              MoleFractionMineral = MoleFractionAqueousCommon(isotopologue)
              DO i = 1,ncomp
                jac_rmin(i,np,k) =  surf(np,k)*actenergy(np,k)*rate0(np,k)* &
                     ( MoleFractionMineral*pre_rmin(np,k)*jac_sat(i) + MoleFractionMineral*jac_pre(i,np)*AffinityTerm  +  &
                      pre_rmin(np,k)*AffinityTerm*dMoleFractionAqueousCommon(i,isotopologue)  )
              END DO
            ELSE
              MoleFractionMineral = MoleFractionMineralCommon(kPointerIsotope(k))
              DO i = 1,ncomp
                jac_rmin(i,np,k) =  MoleFractionMineral*surf(np,k)*actenergy(np,k)*rate0(np,k)* &
                     ( pre_rmin(np,k)*jac_sat(i) + jac_pre(i,np)*AffinityTerm )
              END DO
            END IF

          ELSE IF (IsotopeMineralRare(k)) THEN

            kIsotopologue = kPointerIsotope(k)

            IF (MineralAssociate(kMineralCommon)) THEN
              kIsotopologuePoint = kPointerIsotope( MineralID(kMineralCommon) )
            ELSE
              kIsotopologuePoint = kIsotopologue
            END IF

            IF (isotopeBackReactionOption(kIsotopologuePoint) == 'none' .OR. UseAqueousMoleFraction(kIsotopologue)) THEN

              isotopologue = PointerToPrimaryIsotope(kIsotopologue)
              MoleFractionMineral = MoleFractionAqueousRare(isotopologue)
!!            Use the bulk surface area (common)
              surf(np,k) = surf(np,kIsotopeCommon(kIsotopologue))
              DO i = 1,ncomp
                jac_rmin(i,np,k) =  surf(np,k)*actenergy(np,k)*rate0(np,k)* &
                     ( MoleFractionMineral*pre_rmin(np,k)*jac_sat(i) +   &
                      MoleFractionMineral*jac_pre(i,np)*AffinityTerm  +  &
                      pre_rmin(np,k)*AffinityTerm*dMoleFractionAqueousRare(i,isotopologue)  )
              END DO
              
            ELSE
              
              MoleFractionMineral = MoleFractionMineralRare(kPointerIsotope(k))
!!            Use the bulk surface area (common)
              surf(np,k) = surf(np,kIsotopeCommon(kIsotopologue))
              KMineralCommon = kIsotopeCommon(kIsotopologue)
          
              IF (NoFractionationDissolution .and. IsotopeMineralRare(k) .and. si(np,k) < 1.0d0 ) THEN

                DO i = 1,ncomp
                  jac_rmin(i,np,k) =  MoleFractionMineral*surf(np,k)*actenergy(np,k)*rate0(np,kMineralCommon)* &
                     ( pre_rmin(np,k)*jac_sat(i) + jac_pre(i,np)*AffinityTerm )
                END DO
              ELSE
                DO i = 1,ncomp
                  jac_rmin(i,np,k) =  MoleFractionMineral*surf(np,k)*actenergy(np,k)*rate0(np,k)* &
                     ( pre_rmin(np,k)*jac_sat(i) + jac_pre(i,np)*AffinityTerm )     
                END DO
              END IF
              
            END IF

          ELSE
              
            DO i = 1,ncomp
              jac_rmin(i,np,k) =  surf(np,k)*actenergy(np,k)*rate0(np,k)* &
                     ( pre_rmin(np,k)*jac_sat(i) + jac_pre(i,np)*AffinityTerm )
            END DO
            
          END IF

        ELSE    !!  Non-isotope case
          
          !!! Steefel Checked
          DO i = 1,ncomp
            
            IF (umin(k)=='TOC_soil' .OR. umin(k)=='TOCsoil') THEN
              liqsat_fac = 1
              IF (satliq(jx,jy,jz) > thres_OM2) THEN
                liqsat_fac = 1/(1 + (thres_OM1/satliq(jx,jy,jz))**exp_OM)
              ELSE
                liqsat_fac = 1/(1 + (thres_OM1/thres_OM2)**exp_OM)
              ENDIF
              jac_rmin(i,np,k) =  liqsat_fac*surf(np,k)*actenergy(np,k)*rate0(np,k)* &
                     ( pre_rmin(np,k)*jac_sat(i) + jac_pre(i,np)*AffinityTerm )
            ELSEIF (umin(k)=='Root_respiration' .or. umin(k)=='Root_exudates') THEN
              liqsat_fac = 1/(1 + (thres_root/satliq(jx,jy,jz))**exp_root)
              jac_rmin(i,np,k) =  liqsat_fac*surf(np,k)*actenergy(np,k)*rate0(np,k)* &
                     ( pre_rmin(np,k)*jac_sat(i) + jac_pre(i,np)*AffinityTerm )
            ELSE
              jac_rmin(i,np,k) =  surf(np,k)*actenergy(np,k)*rate0(np,k)* &
                     ( pre_rmin(np,k)*jac_sat(i) + jac_pre(i,np)*AffinityTerm )
            ENDIF
            
          END DO
          
        END IF

      END IF    

      
    END DO    !! End of np loop through parallel reactions
    
  END IF      !! End of IF block for calculation of Jacobian
  
END DO        !! End of loop through minerals


RETURN
END SUBROUTINE jacmin
!  ***************************************************************
