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

SUBROUTINE jacmin(ncomp,nspec,nexchange,nsurf,nkin,nrct,jx,jy,jz)
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
  SUBROUTINE AffinityNumerical(ncomp,nrct,jx,jy,jz,np,k,sppTMP,snormTMP)
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
    REAL(DP), INTENT(OUT)                                           :: snormTMP
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
REAL(DP)                                             :: Kformation
REAL(DP)                                             :: Al
REAL(DP)                                             :: PartialSnorm

!!  For numerical derivatives
!!REAL(DP),DIMENSION(nreactmax,nkin)                   :: rminTMP
!!REAL(DP),DIMENSION(nreactmax,nkin)                   :: rminORIG
!!REAL(DP),DIMENSION(ncomp)                            :: jac_RateCheck
!!REAL(DP), DIMENSION(ncomp,nreactmax,nkin)           :: jac_check

CHARACTER (LEN=5)                                    :: Dumstring

REAL(DP)                                             :: dependTMP
REAL(DP)                                             :: surftmp


!!  sppTMP(:) = sp(:,jx,jy,jz)


tk = t(jx,jy,jz) + 273.15D0
tkinv = 1.0d0/tk
reft = 1.0d0/298.15d0
porfactor = (por(jx,jy,jz)/porin(jx,jy,jz))**(0.666666)

jac_rmin = 0.0d0

DO k = 1,nkin

  IF (ivolume(k) == 1) THEN
    CONTINUE
  ELSE
    
    jac_pre = 0.0d0
    
    DO np = 1,nreactmin(k)
      
      jac_sat = 0.0d0
      jac_RateFactor = 0.0d0

!***********  What surface area to use ************************
      
!!  Taken from reaction.F90
      
! ************* Saturation state dependence of rate ********
   
      IF (si(np,k) > 1.0d0) THEN
        sign = 1.0d0
      ELSE
        sign = -1.0d0
      END IF  

      IF (KateMaher) THEN
        IF (k == KUCalcite .AND. rlabel(np,k) == 'recrystallize') THEN
          sign = 1.0d0
        END IF
      END IF

      IF (imintype(np,k) == 4) THEN                        !! Precipitation only rate law
        term1 = snorm(np,k)
      ELSE                                                 !!  TST forms
        IF (AffinityDepend1(np,k) == 1.0D0) THEN
          term1 = sign*DABS(snorm(np,k) - 1.0D0)
        ELSE
          term1 = sign*(DABS(snorm(np,k) - 1.0D0))**(AffinityDepend1(np,k))
        END IF
      END IF

      IF (imintype(np,k) == 3 .OR. imintype(np,k) == 2) THEN           !! Monod or irreversible
!!        AffinityTerm = 1.0d0
        AffinityTerm = sign
      ELSE                                                            !! TST forms
        AffinityTerm = term1
      END IF

!!      IF (AffinityDepend2(np,k) /= 1.0D0 .OR. AffinityDepend3(np,k) /= 1.0d0) THEN
      IF (AffinityDepend3(np,k) /= 1.0d0) THEN
        PartialSnorm = -sign*snorm(np,k)*AffinityDepend2(np,k)*AffinityDepend3(np,k)*DABS(siln(np,k))**(AffinityDepend3(np,k)-1.0d0)
        DO i = 1,ncomp
          jac_sat(i) =  PartialSnorm * decay_correct(i,k)*mumin(np,k,i)*si(np,k)/si(np,k)
        END DO 
!!        perturb = 1.d-08
!!        write(*,*) ' Mineral: ',k
!!        write(*,*) ' SI = ',si(np,k)
!!        do i = 1,ncomp
!!          IF (mumin(np,k,i) /= 0.0d0) THEN
!!            sppTMP(i) = sppTMP(i) + perturb
!!            CALL AffinityNumerical(ncomp,nrct,jx,jy,jz,np,k,sppTMP,snormTMP)
!!            jac_check(i) = (snormTMP - snorm(np,k))/perturb
!!            sppTMP(i) = sppTMP(i) - perturb
!!            write(*,*) i,jac_sat(i),jac_check(i)
!!            read(*,*)
!!          END IF
!!        end do
!!        write(*,*)
      ELSE IF (AffinityDepend2(np,k) /= 1.0d0) THEN
        IF (thresh(np,k) > 0.0d0 .AND. silog(np,k) > 0.0d0) THEN
          jac_sat(i) = 0.0d0
        ELSE
          PartialSnorm = AffinityDepend2(np,k)*si(np,k)**(AffinityDepend2(np,k)-1.0d0)
          DO i = 1,ncomp
            jac_sat(i) =  PartialSnorm * decay_correct(i,k)*mumin(np,k,i)*si(np,k)
          END DO 
        END IF
      ELSE
        IF (kcrossaff(np,k) /= 0) THEN
          kk = kcrossaff(np,k)
          DO i = 1,ncomp
            jac_sat(i) = decay_correct(i,k)*mumin(1,kk,i)*si(1,kk)
          END DO
        ELSE
          DO i = 1,ncomp
            jac_sat(i) = decay_correct(i,k)*mumin(np,k,i)*si(np,k)
          END DO
        END IF
      END IF
     
!******** Far from equilibrium dependence of rate calculation *********
      
      IF (imintype(np,k) == 1 .OR. imintype(np,k) == 3 .OR. imintype(np,k) == 4) THEN  !! TST, irreversible, or ppt only

        IF (k == kUPlag .AND. OelkersRateLaw) THEN                                     !! Oelkers' Rate Law
 
          Al = sp10(ikAl,jx,jy,jz)
          Kformation = 0.0000149823d0
          jac_pre(ikAl,np) = -0.33333333d0*Al*Kformation/(Kformation + Al**0.333333333333)**2.0d0*Al**0.666666666666

        ELSE

          DO kk = 1,ndepend(np,k)
            i = idepend(kk,np,k)
            IF (itot_min(kk,np,k) == 1) THEN    ! Dependence on total concentration
              DO i2 = 1,ncomp
                jac_pre(i2,np) = jac_pre(i2,np) +  &
                  fjac(i2,i,jx,jy,jz)*depend(kk,np,k)* pre_rmin(np,k)/s(i,jx,jy,jz)
              END DO
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
          END DO     !  End of of loop through aqueous species dependences
          DO kk = 1,ndependex(np,k)
            ix = ixdepend(kk,np,k)
            DO i2 = 1,ncomp+nexchange
              jac_pre(i2,np) = jac_pre(i2,np) + pre_rmin(np,k)*dependex(kk,np,k)*  &
                muexc(ix,i2)
            END DO
          END DO     !  End of of loop through exchange species dependences

          DO kk = 1,ndependsurf(np,k)
            is = isdepend(kk,np,k)
            dependTMP = dependsurf(kk,np,k)
            IF (is <= nsurf) THEN   !  Primary surface complex dependence
              jac_pre(is+ncomp+nexchange,np) = jac_pre(is+ncomp+nexchange,np)  &
                + pre_rmin(np,k)*dependTMP
            ELSE                   !  Secondary species dependence
              iss = is - nsurf
              DO i2 = 1,ncomp+nsurf
                jac_pre(i2,np) = jac_pre(i2,np) + pre_rmin(np,k)*dependTMP*  &
                  musurf(iss,i2)
              END DO
            END IF
          END DO     !  End of of loop through aqueous species dependences

        END IF
        
! REDIMENSION JAC_PRE(1:ncomp+nexchange+nsurf,np)
        
      ELSE IF (imintype(np,k) == 2) THEN   !  Monod rate law
        
        DO kk = 1,nmonod(np,k)
          IF (imonod(kk,np,k) == 0 .AND. kmonod(kk,np,k) == 1) THEN                !! Mineral Monod reaction depends on its own concentration
            CONTINUE                                                       !! No dependence on this term, since mineral concentration is not a primary unknown
          ELSE
            i = imonod(kk,np,k)
            IF (itot_monod(kk,np,k) == 1) THEN    ! Dependence on total concentration
              DO i2 = 1,ncomp
                jac_pre(i2,np) = jac_pre(i2,np) +  &
                  fjac(i2,i,jx,jy,jz)*pre_rmin(np,k)*   &
                      ( 1.0d0/s(i,jx,jy,jz) - 1.0d0/(s(i,jx,jy,jz)+halfsat(kk,np,k))  )
              END DO
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
        
      END IF
      
      IF (snorm(np,k) == 1.0d0 .OR. imintype(np,k) == 4) THEN   !!  If a ppt only rate law, don't use the TST (1 - snorm) form
        AffinityJacobianTerm = 1.0d0
      ELSE
        IF (AffinityDepend1(np,k) /= 1.0d0) THEN
          AffinityJacobianTerm = sign*AffinityDepend1(np,k)*(DABS(snorm(np,k) - 1.0d0))**(AffinityDepend1(np,k)-1.0d0)
        ELSE
!!          AffinityJacobianTerm = 1.0d0
          AffinityJacobianTerm = sign
        END IF
      END IF

! --> Following is for linear TST only, with threshold set greater than zero

      IF (thresh(np,k) > 0.0d0 .AND. AffinityDepend1(np,k) == 1.0d0 .AND. AffinityDepend2(np,k) == 1.0d0) THEN

        IF (silog(np,k) < -2.0d0*thresh(np,k)) THEN
          RateFactor = 1.0d0
        ELSE IF (silog(np,k) > 0.0d0) THEN
          RateFactor = 0.0d0
        ELSE 
          CorrectedTerm = silog(np,k) + thresh(np,k)
          CubicTerm = CorrectedTerm*CorrectedTerm*CorrectedTerm/    &
                       ( thresh(np,k)*thresh(np,k)*thresh(np,k) )
          RateFactor = 0.5d0 - 0.75d0*CorrectedTerm/thresh(np,k) + 0.25d0*CubicTerm
        END IF
   
        perturb = 1.d-08
        DO i = 1,ncomp
          IF (mumin(np,k,i) /= 0.0) THEN
            sppTMP(i) = sppTMP(i) + perturb
            CALL RateFactorNumerical(ncomp,nrct,jx,jy,jz,np,k,sppTMP,RateFactorTMP)
            jac_RateFactor(i) = (RateFactorTMP - RateFactor)/perturb
            sppTMP(i) = sppTMP(i) - perturb
          END IF
        END DO

        IF (silog(np,k) < -2.0d0*thresh(np,k)) THEN
          jac_RateFactor = 0.0
        ELSE IF (silog(np,k) >= 0.0d0) THEN
          jac_RateFactor = 0.0d0
        ELSE
          DO i = 1,ncomp
            jac_RateFactor(i) = ( -mumin(np,k,i)*0.325721d0/thresh(np,k) +   &
                             0.325721d0*mumin(np,k,i)*(thresh(np,k)+silog(np,k))**2.0d0/(thresh(np,k)**3.0d0) )
          END DO
        END IF               
          
        CONTINUE

!!  Or for case with threshold set, but nonlinear TST

      ELSE                               

        IF (thresh(np,k) > 0.0d0 .AND. silog(np,k) > 0.0d0) THEN
          RateFactor = 0.0d0
        ELSE
          RateFactor = 1.0d0
        END IF
        jac_RateFactor = 0.0d0
      END IF

!! Precipitation-only rate law, undersaturated

      IF (imintype(np,k) == 4 .AND. silog(np,k) < 0.0d0) THEN   
        RateFactor = 0.0d0
        jac_RateFactor = 0.0d0
      END IF

!! Monod or irreversible rate law (no saturation dependence)

      IF (imintype(np,k) == 3 .OR. imintype(np,k) == 2) THEN      
        jac_RateFactor = 0.0d0
      END IF

!! Or linear TST case with threshold set

      IF (thresh(np,k) > 0.0d0 .AND. AffinityDepend1(np,k) == 1.0d0 .AND. AffinityDepend2(np,k) == 1.0d0) THEN 
        jac_sat = 0
        AffinityJacobianTerm = 1.0d0
        AffinityTerm = 1.0d0
      END IF
     
!! Now the chain rule ---->
  
     DO i = 1,ncomp+nexchange+nsurf
       jac_rmin(i,np,k) =  surf(k)*actenergy(np,k)*rate0(np,k)*     &
            (sign*RateFactor*AffinityJacobianTerm*pre_rmin(np,k)*jac_sat(i) +  &
             RateFactor*jac_pre(i,np)*AffinityTerm +                  &  
             pre_rmin(np,k)*AffinityTerm*jac_RateFactor(i) )
     END DO
      
    END DO
    
  END IF
  
END DO


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
