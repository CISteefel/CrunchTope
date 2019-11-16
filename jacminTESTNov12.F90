!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 09:59:42
 
!************** (C) COPYRIGHT 1995,1998,1999 ******************
!*******************     C.I. Steefel      *******************
!                    All Rights Reserved

!  GIMRT98 IS PROVkkED "AS IS" AND WITHOUT ANY WARRANTY EXPRESS OR IMPLIED.
!  THE USER ASSUMES ALL RISKS OF USING GIMRT98. THERE IS NO CLAIM OF THE
!  MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.

!  YOU MAY MODIFY THE SOURCE CODE FOR YOUR OWN USE, BUT YOU MAY NOT
!  DISTRIBUTE EITHER THE ORIGINAL OR THE MODIFIED CODE TO ANY OTHER
!  WORKSTATIONS
!**********************************************************************

SUBROUTINE jacmin(ncomp,nspec,nexchange,nsurf,nkin,nrct,jx,jy,jz,time,round)
USE crunchtype
USE runtime
USE params
USE concentration
USE mineral
USE solver
USE medium
USE temperature

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
    USE concentration, ONLY: gam
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
interface
  SUBROUTINE RateFactorNumerical(ncomp,nrct,jx,jy,jz,np,k,sppTMP,RateFactor)
    USE crunchtype
    USE concentration, ONLY: gam
    IMPLICIT NONE
    INTEGER(I4B), INTENT(IN)                                        :: ncomp
    INTEGER(I4B), INTENT(IN)                                        :: nrct
    INTEGER(I4B), INTENT(IN)                                        :: jx
    INTEGER(I4B), INTENT(IN)                                        :: jy
    INTEGER(I4B), INTENT(IN)                                        :: jz
    INTEGER(I4B), INTENT(IN)                                        :: np
    INTEGER(I4B), INTENT(IN)                                        :: k
    REAL(DP),DIMENSION(:),INTENT(IN)                                :: sppTMP
    REAL(DP), INTENT(OUT)                                           :: RateFactor
  END SUBROUTINE RateFactorNumerical
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

REAL(DP), DIMENSION(NCOMP)                           :: dMoleFraction40
REAL(DP), DIMENSION(NCOMP)                           :: dMoleFraction44
REAL(DP)                                             :: MoleFraction40Tmp
REAL(DP)                                             :: MoleFraction44Tmp
REAL(DP)                                             :: MoleFraction40
REAL(DP)                                             :: MoleFraction44
REAL(DP), DIMENSION(NCOMP)                           :: dMoleFraction32
REAL(DP), DIMENSION(NCOMP)                           :: dMoleFraction34
REAL(DP), DIMENSION(NCOMP)                           :: dMoleFraction32S
REAL(DP), DIMENSION(NCOMP)                           :: dMoleFraction34S
REAL(DP)                                             :: MoleFraction32Tmp
REAL(DP)                                             :: MoleFraction34Tmp
REAL(DP)                                             :: MoleFraction32
REAL(DP)                                             :: MoleFraction34
REAL(DP)                                             :: MoleFraction32S
REAL(DP)                                             :: MoleFraction34S
REAL(DP)                                             :: DenomSquared
REAL(DP)                                             :: Denom

REAL(DP)                                             :: MoleFraction44Mineral
REAL(DP)                                             :: MoleFraction40Mineral
REAL(DP)                                             :: CurrentMoleFraction40
REAL(DP)                                             :: CurrentMoleFraction44
REAL(DP)                                             :: denominator
LOGICAL(LGT)                                         :: SetToAqueousMoleFraction

REAL(DP)                                             :: CalciumCarbonateRatioEffect
REAL(DP)                                             :: RatioToPower
REAL(DP), DIMENSION(NCOMP)                           :: DerivativeActivityRatio


! biomass
INTEGER(I4B)                                         :: ib, jj
! biomass end

SetToAqueousMoleFraction = .FALSE.
sppTMP(:) = sp(:,jx,jy,jz)
sppTMP10(:) = sp10(:,jx,jy,jz)

!! JennyDruhan block  !!**************************************************************
IF (JennyDruhan) THEN       
!!CIS Hardwired for 44Ca40CaCO3
!!  MoleFraction44 = sppTMP10(7)/( sppTMP10(6) + sppTMP10(7) )
!!  MoleFraction40 = 1.0d0 - MoleFraction44
  MoleFraction44 = sn(7,jx,jy,jz)/( sn(6,jx,jy,jz) + sn(7,jx,jy,jz) )
  MoleFraction40 = 1.0d0 - MoleFraction44
  dMoleFraction40 = 0.0d0
  dMoleFraction44 = 0.0d0

!! For dissolution, use the mole fraction of 44Ca and 40Ca in the mineral from previous time step

!!  IF (volfx(6,jx,jy,jz) == 0.0d0 .AND. volfx(5,jx,jy,jz) == 0.0d0 .OR. UseAqueousMoleFraction .OR. time==0.0d0) THEN
  IF (volfx(6,jx,jy,jz) == 0.0d0 .AND. volfx(5,jx,jy,jz) == 0.0d0 .OR. UseAqueousMoleFraction) THEN
    MoleFraction44Mineral = MoleFraction44
    MoleFraction40Mineral = MoleFraction40

  ELSE

    IF (time < LagSurface/365.0 .OR. UseBulkMineral) THEN
      MoleFraction44Mineral = volfx(6,jx,jy,jz)/( volfx(6,jx,jy,jz) + volfx(5,jx,jy,jz) )
      MoleFraction40Mineral = volfx(5,jx,jy,jz)/( volfx(6,jx,jy,jz) + volfx(5,jx,jy,jz) )
    ELSE
      MoleFraction44Mineral = VolSave(6,jx,jy,jz)/( VolSave(6,jx,jy,jz) + VolSave(5,jx,jy,jz) )
      MoleFraction40Mineral = VolSave(5,jx,jy,jz)/( VolSave(6,jx,jy,jz) + VolSave(5,jx,jy,jz) )
    END IF

  END IF

!!CIS Hardwired for 44Ca40CaCO3

!!CIS Hardwired for Fe34S32S
  MoleFraction34 = sppTMP10(13)/( sppTMP10(12) + sppTMP10(13) )
  MoleFraction32 = 1.0d0 - MoleFraction34
  dMoleFraction32 = 0.0d0
  dMoleFraction34 = 0.0d0
!!CIS Hardwired for Fe34S32S

!!CIS Hardwired for 34S32S (elemental sulfur)

  MoleFraction34S = sp10(13,jx,jy,jz)/( sp10(12,jx,jy,jz) + sp10(13,jx,jy,jz) )
  MoleFraction32S = 1.0d0 - MoleFraction34S
  dMoleFraction32S = 0.0d0
  dMoleFraction34S = 0.0d0

!!CIS Hardwired for 34S32S (elemental sulfur)

END IF    
!!  End of JennyDruhan block  ************************************************************

tk = t(jx,jy,jz) + 273.15D0
tkinv = 1.0d0/tk
reft = 1.0d0/298.15d0
porfactor = (por(jx,jy,jz)/porin(jx,jy,jz))**(0.666666)

!! JennyDruhan block   ********************************************************
IF (JennyDruhan) THEN      
!!CIS Hardwired for 44Ca40CaCO3
!!  NOTE:  Species 6 is 40Ca, Species 7 is 44Ca
  Denom = ( sppTMP10(6) + sppTMP10(7) )
  DenomSquared = ( sppTMP10(6) + sppTMP10(7) ) * ( sppTMP10(6) + sppTMP10(7) )
  tmp1 = -sppTMP10(6)/DenomSquared + 1.0d0/Denom
  tmp2 = -sppTMP10(6)/DenomSquared
  tmp3 = -sppTMP10(7)/DenomSquared + 1.0d0/Denom
  tmp4 = -sppTMP10(7)/DenomSquared

  dMoleFraction40(6) = spptmp10(6)*tmp1
  dMoleFraction40(7) = spptmp10(7)*tmp2
  dMoleFraction44(6) = spptmp10(6)*tmp4
  dMoleFraction44(7) = spptmp10(7)*tmp3
!!CIS Hardwired for 44Ca40CaCO3

!!CIS Hardwired for Fe34S32S
!!  NOTE:  Species 12 is 32S, Species 13 is 34S
  Denom = ( sppTMP10(12) + sppTMP10(13) )
  DenomSquared = ( sppTMP10(12) + sppTMP10(13) ) * ( sppTMP10(12) + sppTMP10(13) )
  tmp1 = -sppTMP10(12)/DenomSquared + 1.0d0/Denom
  tmp2 = -sppTMP10(12)/DenomSquared
  tmp3 = -sppTMP10(13)/DenomSquared + 1.0d0/Denom
  tmp4 = -sppTMP10(13)/DenomSquared

  dMoleFraction32(12) = spptmp10(12)*tmp1
  dMoleFraction32(13) = spptmp10(13)*tmp2
  dMoleFraction34(12) = spptmp10(12)*tmp4
  dMoleFraction34(13) = spptmp10(13)*tmp3
!!CIS Hardwired for Fe34S32S

!!CIS Hardwired for 34S32S (elemental sulfur)
  dMoleFraction32S(12) = spptmp10(12)*tmp1
  dMoleFraction32S(13) = spptmp10(13)*tmp2
  dMoleFraction34S(12) = spptmp10(12)*tmp4
  dMoleFraction34S(13) = spptmp10(13)*tmp3
!!CIS Hardwired for 34S32S (elemental sulfur)

END IF    
!!  End of JennyDruhan block  ********************************************************

!! ********** Numerical derivatives for 44Ca and 40Ca mole fractions *****************
!!      perturb = 1.d-10 
!!      i = 6
!!      sppTMP(i) = sppTMP(i) + perturb
!!      MoleFraction44Tmp = DEXP(sppTMP(7))/( DEXP(sppTMP(6)) + DEXP(sppTMP(7)) )
!!      MoleFraction40Tmp = 1.0d0 - MoleFraction44Tmp
!!      dMoleFraction40(6) = (MoleFraction40Tmp - MoleFraction40)/perturb
!!      dMoleFraction44(6) = (MoleFraction44Tmp - MoleFraction44)/perturb
!!      sppTMP(i) = sppTMP(i) - perturb 
!!      write(*,*) dMoleFraction40(6)
!!      write(*,*) dMoleFraction44(6)
!!      dMoleFraction44(7) = dMoleFraction40(6)
!!      dMoleFraction40(7) = dMoleFraction44(6)
!!      i = 7
!!      sppTMP(i) = sppTMP(i) + perturb
!!      sppTMP10(7) = DEXP(sppTMP(7))
!!      MoleFraction44Tmp = DEXP(sppTMP(7))/( DEXP(sppTMP(6)) + DEXP(sppTMP(7)) )
!!      MoleFraction40Tmp = 1.0d0 - MoleFraction44Tmp
!!      dMoleFraction40(7) = (MoleFraction40Tmp - MoleFraction40)/perturb
!!      dMoleFraction44(7) = (MoleFraction44Tmp - MoleFraction44)/perturb
!!      sppTMP(i) = sppTMP(i) - perturb 
!!      write(*,*) dMoleFraction40(7)
!!      write(*,*) dMoleFraction44(7)
!!      write(*,*)
!!      write(*,*) dMoleFraction40(6),spptmp10(6)*tmp1
!!      write(*,*) dMoleFraction40(7),spptmp10(7)*tmp2
!!      write(*,*) dMoleFraction44(6),spptmp10(6)*tmp4
!!      write(*,*) dMoleFraction44(7),spptmp10(7)*tmp3
!!      read(*,*)
!! ********** Numerical derivatives for 44Ca and 40Ca mole fractions *****************


jac_rmin = 0.0d0

DO k = 1,nkin

!!  Ivolume(k) == 1 means rate is set to just deplete current mineral volume fraction (so no Jacobian needed)
!!  Imintype(1,k) ==6 .OR. 7 means forwward and reverse irreversible, with surface area set constant to 1.0

  IF (ivolume(k) == 1 .OR. imintype(1,k) == 6 .OR. imintype(1,k) == 7) THEN

    jac_rmin(:,:,k) = 0.0d0

  ELSE   !! Case where Jacobian of mineral reactions are calculated (block above skips all computation of Jacobian)
    
    jac_pre = 0.0d0
    
    DO np = 1,nreactmin(k)
!!    DO np = nreactmin(k),1,-1
      
      jac_sat = 0.0d0

!***********  What surface area to use ************************
      
!!  Taken from reaction.F90
      
! ************* Saturation state dependence of rate ********
   
      IF (si(np,k) >= 1.0d0) THEN
        sign = 1.0d0
      ELSE
        sign = -1.0d0
      END IF  

      IF (KateMaher) THEN
        IF (k == KUCalcite .AND. rlabel(np,k) == 'recrystallize') THEN
          sign = 1.0d0
        END IF
      END IF

!!    If .NOT. biomass
      IF (imintype(np,k) /= 8) THEN

        IF (AffinityDepend1(np,k) == 1.0D0) THEN
          term1 = sign*DABS(snorm(np,k) - 1.0D0)
        ELSE
          term1 = sign*DABS(snorm(np,k) - 1.0D0)**(AffinityDepend1(np,k))
        END IF

        IF (imintype(np,k) == 4) THEN                               !! Precipitation only
          AffinityTerm = MAX(0.0d0,term1)
        ELSE IF (imintype(np,k) == 5) THEN                          !! Dissolution only
          AffinityTerm = MIN(0.0d0,term1)
        ELSE
          AffinityTerm = term1
        END IF

        IF (imintype(np,k) == 3 .OR. imintype(np,k) == 2) THEN           !! Monod or irreversible
          AffinityTerm = sign
        END IF

        IF (imintype(np,k) /= 3 .AND. imintype(np,k) /= 2) THEN         !! Everything but Monod and irreversible
            
          IF (AffinityDepend1(np,k) == 1.0d0 .AND. AffinityDepend2(np,k) == 1.0d0 .AND. AffinityDepend3(np,k) == 1.0d0) THEN
            DO i = 1,ncomp
              jac_sat(i) = mumin(np,k,i)*si(np,k)
              if (imintype(np,k) == 4 .AND. si(np,k) < 1.0d0) then   !! Precipitation only
                jac_sat(i) = 0.0d0
              end if
              if (imintype(np,k) == 5 .AND. si(np,k) > 1.0d0) then   !! Dissolution only
                jac_sat(i) = 0.0d0
              end if
!!              write(*,*) i,mumin(np,k,i),jac_check(i)
            END DO
!!            read
          ELSE IF (AffinityDepend1(np,k) /= 1.0d0 .OR. AffinityDepend2(np,k) /= 1.0d0 .AND. AffinityDepend3(np,k) == 1.0d0) THEN
              
!!            SaturationTerm = (si(np,k)**AffinityDepend2(np,k) - 1.0d0)
            SaturationTerm = snorm(np,k) - 1.0

            IF (SaturationTerm == 0.0d0) THEN  !! Trap divide by zero for pathological case in which exponent is less than 1
              SaturationTerm = tiny
            END IF
                
            DO i = 1,ncomp

              IF (AffinityDepend1(np,k) == 1.0d0) THEN
                factor1 = 1.0d0
              ELSE
                factor1 = AffinityDepend1(np,k)*                  &
                  DABS(SaturationTerm)**(AffinityDepend1(np,k)- 1.0d0) 
              END IF
              IF (AffinityDepend2(np,k) == 1.0d0) THEN
                factor2 = 1.0d0
              ELSE
                factor2 = AffinityDepend2(np,k)*(si(np,k))**(AffinityDepend2(np,k)-1.0d0)
              END IF      

!!              jac_check(i) = sign*mumin(np,k,i)*si(np,k)*factor1*factor2
!!                termA = 
                jac_sat(i) = mumin(np,k,i)*si(np,k)* AffinityDepend1(np,k)*AffinityDepend2(np,k)*       &
                           ( si(np,k)**(AffinityDepend2(np,k)-1.0d0) ) * (DABS(SaturationTerm))**(AffinityDepend1(np,k)-1.0d0)

!!!              if (imintype(np,k) == 4 .AND. si(np,k) < 1.0d0) then   !! Precipitation only
!!!                jac_sat(i) = 0.0d0
!!!              end if
!!!              if (imintype(np,k) == 5 .AND. si(np,k) > 1.0d0) then   !! Dissolution only
!!!                jac_sat(i) = 0.0d0
!!!              end if

            END DO              
              
          ELSE
           
!!  Numerical derivative for full Burch-Hellmann rate law (all exponents /= 1.0)
            perturb = 1.d-08 
            DO i = 1,ncomp
!!              IF (mumin(np,k,i) /= 0.0d0) THEN
                sppTMP(i) = sppTMP(i) + perturb
                sppTMP10(i) = DEXP(sppTMP(i))
                CALL AffinityNumerical(ncomp,nrct,jx,jy,jz,np,k,sppTMP,termTMP,time)
                jac_check(i) = (termTMP - AffinityTerm)/perturb
                sppTMP(i) = sppTMP(i) - perturb 
                sppTMP10(i) = sp10(i,jx,jy,jz)
!!              END IF
            END DO
            
          ENDIF
          
          perturb = 1.d-09 
         DO i = 1,ncomp
        !!    IF (mumin(np,k,i) /= 0.0d0) THEN
             sppTMP(i) = sppTMP(i) + perturb
             sppTMP10(i) = DEXP(sppTMP(i))
             CALL AffinityNumerical(ncomp,nrct,jx,jy,jz,np,k,sppTMP,termTMP,time)
             jac_check(i) = (termTMP - AffinityTerm)/perturb
             sppTMP(i) = sppTMP(i) - perturb 
            sppTMP10(i) = sp10(i,jx,jy,jz)
!!            END IF
         END DO
            
          jac_sat = jac_check
        
 !!         write(*,*) umin(k),si(np,k)
 !!         do i = 1,ncomp
 !!           write(*,*) jac_sat(i),jac_check(i)
 !!         end do
 !!         read(*,*)

        END IF


      END IF !!! END of imintype(np,k) /= 8

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
 
      IF (imintype(np,k) == 1 .OR. imintype(np,k) == 3 .OR. imintype(np,k) == 4) THEN  !! TST, irreversible, or ppt only
!!      +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++IF (imintype(np,k) == 1 .OR. imintype(np,k) == 3 .OR. imintype(np,k) == 4) THEN  !! TST, irreversible, or ppt only

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
          
!!        IF (k == kUPlag .AND. OelkersRateLaw) THEN                                     !! Oelkers' Rate Law
         IF (HyperbolicInhibition(np,k)) THEN
!!          Kformation = 0.0000149823D0
           Al = sp10(HyperbolicInhibitionPointer(np,k),jx,jy,jz)
  
           denominator = Kformation(np,k) + Al**HyperbolicInhibitionDepend(np,k)
!!           term2 = Kformation(np,k)/denominator
          
           jac_pre(HyperbolicInhibitionPointer(np,k),np) = jac_pre(HyperbolicInhibitionPointer(np,k),np) -   &
               Kformation(np,k)*HyperbolicInhibitionDepend(np,k)*Al**(HyperbolicInhibitionDepend(np,k)-1.0d0)/(denominator*denominator)
           
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
          IF (imonod(kk,np,k) == 0 .AND. kmonod(kk,np,k) == 1) THEN                !! Mineral Monod reaction depends on its own concentration
            CONTINUE                                                       !! No dependence on this term, since mineral concentration is not a primary unknown
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
          do i = 1,ncomp
            jac_rmin(i,np,k) =  MetabolicLagMineral(np,jj,jx,jy,jz)*surf(k)*actenergy(np,k)*rate0(np,k)*volfx(ib,jx,jy,jz)* &
                          ( pre_rmin(np,k)*jac_sat(i) + jac_pre(i,np)*AffinityTerm  )
          end do
        ELSE
          do i = 1,ncomp
            jac_rmin(i,np,k) =  surf(k)*actenergy(np,k)*rate0(np,k)*volfx(ib,jx,jy,jz)* &
                          ( pre_rmin(np,k)*jac_sat(i) + jac_pre(i,np)*AffinityTerm  )
          end do
        END IF

      ELSE


        IF (JennyDruhan) THEN       !!  ********************************************************
!!      Start of JennyDruhan block
       

!!CIS     Hardwired for 44Ca40CaCO3
          if (k == 5) then

            IF (DePaolo) THEN          !! DePaolo = .TRUE. 

              if (np == 1) then        !! Dissolution
!! Uni-directional dissolution of calcite

                IF (SetToAqueousMoleFraction) THEN
                  do i = 1,ncomp

                    jac_rmin(i,np,k) =  surf(k)*actenergy(np,k)*rate0(np,k)* &
                          ( MoleFraction40*pre_rmin(np,k)*jac_sat(i) + MoleFraction40*jac_pre(i,np)*AffinityTerm +  &
                           pre_rmin(np,k)*AffinityTerm*dMoleFraction40(i)  )
                 
                  end do

                ELSE

                  do i = 1,ncomp

                    jac_rmin(i,np,k) =  MoleFraction40Mineral*surf(k)*actenergy(np,k)*rate0(np,k)* &
                          ( pre_rmin(np,k)*jac_sat(i) + jac_pre(i,np)*AffinityTerm )

!!                    denominator = rmin(2,5) + rmin(2,6)
!!                    CurrentMoleFraction40 = rmin(2,5)/denominator
!!
!!                    jac_rmin(i,np,k) =  surf(k)*actenergy(np,k)*rate0(np,k)*                                                          &
!!                          (    CurrentMoleFraction40*pre_rmin(np,k)*jac_sat(i) + CurrentMoleFraction40*jac_pre(i,np)*AffinityTerm     &                                                                                                        &
!!                                              +                                                                                       &
!!                               pre_rmin(np,k)*AffinityTerm*                                                                           &
!!                          ( jac_rmin(i,2,5)/(denominator) * (-1.0)/(denominator*denominator)*( jac_rmin(i,2,5) + jac_rmin(i,2,6)  ))   ) 


!!                    denominator = rmin(2,5) + rmin(2,6) + rmin(1,5) + rmin(1,6)
!!                    CurrentMoleFraction40 = (rmin(2,5) + rmin(1,5))/denominator
!!                    jac_rmin(i,np,k) =  surf(k)*actenergy(np,k)*rate0(np,k)*                                                          &
!!                          (                                                                                                           &
!!                              rmin(2,5)/denominator*pre_rmin(np,k)*jac_sat(i) + rmin(2,5)/denominator*jac_pre(i,np)*AffinityTerm  +   &
!!                              pre_rmin(np,k)*AffinityTerm*                                                                            &
!!                              jac_rmin(i,2,5)/denominator * (-1.0)/(denominator*denominator) *                                         &
!!                              ( jac_rmin(i,2,5) + jac_rmin(i,2,6) + jac_rmin(i,1,5) + jac_rmin(i,1,6)  )                              &
!!                                       +                                                                                              &
!!                              rmin(1,5)/denominator*pre_rmin(np,k)*jac_sat(i) + rmin(1,5)/denominator*jac_pre(i,np)*AffinityTerm  +   &
!!                              pre_rmin(np,k)*AffinityTerm*                                                                            &
!!                              jac_rmin(i,1,5)/denominator * (-1.0)/(denominator*denominator) *                                         &
!!                              ( jac_rmin(i,2,5) + jac_rmin(i,2,6) + jac_rmin(i,1,5) + jac_rmin(i,1,6) )                         )

                  end do

!!                  do i = 1,ncomp
!!
!!                    denominator = rmin(2,5) + rmin(2,6) + rmin(1,5) + rmin(1,6)
!!                    CurrentMoleFraction40 = (rmin(2,5)+rmin(1,5))/denominator
!!                    rmin(np,k) = CurrentMoleFraction40*surf(k)*rate0(np,k)*actenergy(np,k)*pre_rmin(np,k)*AffinityTerm
!!
!!                  end do

                END IF

              else if (np == 2) then    !! Precipitation
!! Uni-directional precipitation of calcite

                do i = 1,ncomp
                  jac_rmin(i,np,k) =  surf(k)*actenergy(np,k)*rate0(np,k)* &
                          ( pre_rmin(np,k)*jac_sat(i) + jac_pre(i,np)*(-1.0d0*AffinityTerm)  )
                end do

              else
                write(*,*)
                write(*,*) ' There should be only two parallel rate laws for calcite in DePaolo model'
                stop
              endif

            ELSE                  !! DePaolo = .FALSE.

              DerivativeActivityRatio = 1.0d0
              CalciumCarbonateRatioEffect = 1.0
              RatioToPower = DEXP( 1.228*( sp(6,jx,jy,jz) - sp(27,jx,jy,jz) ) )
              DerivativeActivityRatio(6) = sp10(6,jx,jy,jz)*1.41082/( RatioToPower*sp10(27,jx,jy,jz) )
              DerivativeActivityRatio(14) = muaq(27,14)*sp10(27,jx,jy,jz)*1.41082*sp10(6,jx,jy,jz)/( RatioToPower*sp10(27,jx,jy,jz)*sp10(27,jx,jy,jz) )
              DerivativeActivityRatio(1) = muaq(27,1)*sp10(27,jx,jy,jz)*1.41082*sp10(6,jx,jy,jz)/( RatioToPower*sp10(27,jx,jy,jz)*sp10(27,jx,jy,jz) )
              
              IF (SetToAqueousMoleFraction) THEN
                do i = 1,ncomp

                  jac_rmin(i,np,k) =  surf(k)*actenergy(np,k)*rate0(np,k)* &
                     ( CalciumCarbonateRatioEffect*MoleFraction40*pre_rmin(np,k)*jac_sat(i) + CalciumCarbonateRatioEffect*MoleFraction40*jac_pre(i,np)*AffinityTerm  +  &
                           CalciumCarbonateRatioEffect*pre_rmin(np,k)*AffinityTerm*dMoleFraction40(i) )
                end do
                do i = 1,ncomp
                  jac_rmin(i,np,k) =  MoleFraction40Mineral*(surf(5))*actenergy(np,k)*rate0(np,k)* &
                          ( pre_rmin(np,k)*jac_sat(i) + jac_pre(i,np)*AffinityTerm  )
                end do
              END IF

            END IF

          else if (k == 6) then     !! Hardwired for 44CaCO3

            IF (DePaolo) THEN          !! DePaolo = .TRUE. 

              if (np == 1) then        !! Dissolution
!! Uni-directional dissolution of calcite

                IF (SetToAqueousMoleFraction) THEN
                  do i = 1,ncomp

                    jac_rmin(i,np,k) =  surf(5)*actenergy(np,k)*rate0(np,k)* &
                          ( MoleFraction44*pre_rmin(np,k)*jac_sat(i) + MoleFraction44*jac_pre(i,np)*AffinityTerm +  &
                           pre_rmin(np,k)*AffinityTerm*dMoleFraction44(i) )

                  end do

                ELSE

                  DO i = 1,ncomp

                    jac_rmin(i,np,k) =  MoleFraction44Mineral*surf(5)*actenergy(np,k)*rate0(np,k)* &
                          ( pre_rmin(np,k)*jac_sat(i) + jac_pre(i,np)*AffinityTerm )

!!                    denominator = rmin(2,5) + rmin(2,6)
!!                    CurrentMoleFraction44 = rmin(2,6)/denominator
!!                    jac_rmin(i,np,k) =  surf(5)*actenergy(np,k)*rate0(np,k)*                                                      &
!!                          ( CurrentMoleFraction44*pre_rmin(np,k)*jac_sat(i) + CurrentMoleFraction44*jac_pre(i,np)*AffinityTerm    &
!!                                              +                                                                                   &
!!                           pre_rmin(np,k)*AffinityTerm*                                                                           &
!!                          ( jac_rmin(i,2,6)/(denominator) * (-1.0)/(denominator*denominator)*( jac_rmin(i,2,5) + jac_rmin(i,2,6) ) )   ) 
 

!!                    denominator = rmin(2,5) + rmin(2,6) + rmin(1,5) + rmin(1,6)
!!                    CurrentMoleFraction44 = (rmin(2,6)+rmin(1,6))/denominator
!!                    jac_rmin(i,np,k) =  surf(5)*actenergy(np,k)*rate0(np,k)*                                                          &
!!                          (                                                                                                           &
!!                              rmin(2,6)/denominator*pre_rmin(np,k)*jac_sat(i) + rmin(2,6)/denominator*jac_pre(i,np)*AffinityTerm  +   &
!!                              pre_rmin(np,k)*AffinityTerm*                                                                            &
!!                              jac_rmin(i,2,6)/denominator * (-1.0)/(denominator*denominator)*                                         &
!!                              ( jac_rmin(i,2,5) + jac_rmin(i,2,6) + jac_rmin(i,1,5) + jac_rmin(i,1,6)  )                              &
!!                                       +                                                                                              &
!!                              rmin(1,6)/denominator*pre_rmin(np,k)*jac_sat(i) + rmin(1,6)/denominator*jac_pre(i,np)*AffinityTerm  +   &
!!                              pre_rmin(np,k)*AffinityTerm*                                                                            &
!!                              jac_rmin(i,1,6)/denominator * (-1.0)/(denominator*denominator)*                                         &
!!                              ( jac_rmin(i,2,5) + jac_rmin(i,2,6) + jac_rmin(i,1,5) + jac_rmin(i,1,6) )                         )

                  END DO

                END IF

              else if (np == 2) then    !! Precipitation
!! Uni-directional precipitation of calcite

                do i = 1,ncomp
                  jac_rmin(i,np,k) =  surf(5)*actenergy(np,k)*rate0(np,k)* &
                          ( pre_rmin(np,k)*jac_sat(i) + jac_pre(i,np)*(-1.0d0*AffinityTerm)  )
                end do

              else
                write(*,*)
                write(*,*) ' There should be only two parallel rate laws for calcite in DePaolo model'
                stop
              endif

            ELSE                       !! DePaolo = .FALSE.
                
              DerivativeActivityRatio = 1.0d0
              CalciumCarbonateRatioEffect = 1.0
              RatioToPower = DEXP( 1.228*( sp(6,jx,jy,jz) - sp(27,jx,jy,jz) ) )
              DerivativeActivityRatio(6) = sp10(6,jx,jy,jz)*1.41082/( RatioToPower*sp10(27,jx,jy,jz) )
              DerivativeActivityRatio(14) = muaq(27,14)*sp10(27,jx,jy,jz)*1.41082*sp10(6,jx,jy,jz)/( RatioToPower*sp10(27,jx,jy,jz)*sp10(27,jx,jy,jz) )
              DerivativeActivityRatio(1) = muaq(27,1)*sp10(27,jx,jy,jz)*1.41082*sp10(6,jx,jy,jz)/( RatioToPower*sp10(27,jx,jy,jz)*sp10(27,jx,jy,jz) )
              
              do i = 1,ncomp
                IF (SetToAqueousMoleFraction) THEN
                  jac_rmin(i,np,k) =  surf(5)*actenergy(np,k)*rate0(np,k)* &
                     ( CalciumCarbonateRatioEffect*MoleFraction44*pre_rmin(np,k)*jac_sat(i) + CalciumCarbonateRatioEffect*MoleFraction44*jac_pre(i,np)*AffinityTerm  +  &
                      CalciumCarbonateRatioEffect*pre_rmin(np,k)*AffinityTerm*dMoleFraction44(i)  )
                ELSE
                  jac_rmin(i,np,k) =  MoleFraction44Mineral*(surf(5))*actenergy(np,k)*rate0(np,k)* &
                          ( pre_rmin(np,k)*jac_sat(i) + jac_pre(i,np)*AffinityTerm  )
                END IF
              end do
            END IF

!!CIS Hardwired for 44Ca40CaCO3

!!CIS Hardwired for Fe34S32S
          else if (k == 8) then
            do i = 1,ncomp
              jac_rmin(i,np,k) =  surf(k)*actenergy(np,k)*rate0(np,k)* &
                          ( MoleFraction32*pre_rmin(np,k)*jac_sat(i) + MoleFraction32*jac_pre(i,np)*AffinityTerm +  &
                           pre_rmin(np,k)*AffinityTerm*dMoleFraction32(i) )
            end do

          else if (k == 9) then
            do i = 1,ncomp
             jac_rmin(i,np,k) =  surf(9)*actenergy(np,k)*rate0(np,k)* &
                          ( MoleFraction34*pre_rmin(np,k)*jac_sat(i) + MoleFraction34*jac_pre(i,np)*AffinityTerm +  &
                           pre_rmin(np,k)*AffinityTerm*dMoleFraction34(i) )
            end do
!!CIS Hardwired for Fe34S32S

!!CIS Hardwired for 34S32S (Elemental Sulfur)
          else if (k == 10) then
            do i = 1,ncomp
              jac_rmin(i,np,k) =  surf(k)*actenergy(np,k)*rate0(np,k)* &
                          ( MoleFraction32S*pre_rmin(np,k)*jac_sat(i) + MoleFraction32S*jac_pre(i,np)*AffinityTerm +  &
                           pre_rmin(np,k)*AffinityTerm*dMoleFraction32S(i) )
            end do

          else if (k == 11) then
            do i = 1,ncomp
              jac_rmin(i,np,k) =  surf(10)*actenergy(np,k)*rate0(np,k)* &
                          ( MoleFraction34S*pre_rmin(np,k)*jac_sat(i) + MoleFraction34S*jac_pre(i,np)*AffinityTerm +  &
                           pre_rmin(np,k)*AffinityTerm*dMoleFraction34S(i) )
            end do
!!CIS Hardwired for 34S32S (Elemental Sulfur)

          else
            do i = 1,ncomp
              jac_rmin(i,np,k) =  surf(k)*actenergy(np,k)*rate0(np,k)* &
                          ( pre_rmin(np,k)*jac_sat(i) + jac_pre(i,np)*AffinityTerm  )
            end do

          endif   
!!      End of JennyDruhan block  ******************************************************************

        ELSE

          do i = 1,ncomp
            jac_rmin(i,np,k) =  surf(k)*actenergy(np,k)*rate0(np,k)* &
                          ( pre_rmin(np,k)*jac_sat(i) + jac_pre(i,np)*AffinityTerm  )
          IF (si(np,k) >= IntervalBelowEquilibrium .AND. si(np,k) <= IntervalAboveEquilibrium .AND. rate0Prime(np,k) /= 0.0d0) THEN
!!!              jac_rmin(i,np,k) =  surf(k)*actenergy(np,k)*rate0Prime(np,k)* &
!!!                          ( pre_rmin(np,k)*jac_sat(i) + jac_pre(i,np)*AffinityTerm  )
            END IF
          end do

        END IF

      END IF    
!!    biomass end     
      
    END DO    !! End of np loop through parallel reactions
    
  END IF      !! End of IF block for calculation of Jacobian
  
END DO        !! End of loop through minerals


!!perturb = 1.d-06

!!do i = 1,ncomp
!!  sppTMP(i) = sppTMP(i) + perturb
!!  CALL ReactionNumerical(ncomp,nkin,nrct,nspec,nexchange,nsurf,jx,jy,jz,rminTMP)
!!  do k = 1,nkin
!!    do np = 1,nreactmin(k)
!!      jac_check(i,np,k) = (rminTMP(np,k) - rmin(np,k))/perturb
!!    end do
!!  end do
!!  sppTMP(i) = sppTMP(i) - perturb
!!end do

!!jac_rmin = jac_check

!!   do k = 1,nkin
!!     do np = 1,nreactmin(k)
!!     do i = 1,ncomp
!!       CheckDiff = (jac_check(i,np,k) - jac_rmin(i,np,k))/jac_rmin(i,np,k)
!!       IF (DABS(CheckDiff) > 1.E-08) THEN
!!         write(*,*) umin(k)
!!         write(*,*) CheckDiff,jac_check(i,np,k),jac_rmin(i,np,k)
!!         write(*,*) ' jx = ',jx
!!         read(*,*)
!!        WRITE(*,*)
!!       END IF
!!     end do
!!   end do
!!  end do

!!  continue

RETURN
END SUBROUTINE jacmin
!  ***************************************************************
