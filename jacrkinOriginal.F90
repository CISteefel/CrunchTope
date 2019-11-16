SUBROUTINE jacrkin(ncomp,nspec,nrct,ikin,jx,jy,jz,AqueousToBulk)
USE crunchtype
USE params
USE runtime
USE concentration
USE mineral
USE solver
USE medium
USE isotope

! biomass
USE temperature, ONLY: ro,T
! biomass end

IMPLICIT NONE

!  External variables

INTEGER(I4B), INTENT(IN)                             :: ncomp
INTEGER(I4B), INTENT(IN)                             :: nspec
INTEGER(I4B), INTENT(IN)                             :: nrct
INTEGER(I4B), INTENT(IN)                             :: ikin
INTEGER(I4B), INTENT(IN)                             :: jx
INTEGER(I4B), INTENT(IN)                             :: jy
INTEGER(I4B), INTENT(IN)                             :: jz
REAL(DP), INTENT(IN)                                 :: AqueousToBulk

!  Internal variables

REAL(DP)                                             :: affinity
REAL(DP)                                             :: denom
REAL(DP)                                             :: denomSquared
REAL(DP)                                             :: derivative
REAL(DP)                                             :: jac_inhibition
REAL(DP)                                             :: MonodTerm
REAL(DP)                                             :: InhibitTerm

INTEGER(I4B)                                         :: ir
INTEGER(I4B)                                         :: i
INTEGER(I4B)                                         :: i2
INTEGER(I4B)                                         :: ll
INTEGER(I4B)                                         :: id
INTEGER(I4B)                                         :: ksp
INTEGER(I4B)                                         :: k

! biomass 
INTEGER(I4B)                                         :: ib, jj

INTEGER(I4B)                                         :: IsotopologueOther
INTEGER(I4B)                                         :: nnisotope

REAL(DP)                                             :: bqTMP
REAL(DP)                                             :: tk
REAL(DP)                                             :: sign
REAL(DP)                                             :: term1
REAL(DP)                                             :: perturb
REAL(DP)                                             :: termTMP
REAL(DP)                                             :: snormAqueous

REAL(DP) , DIMENSION(ncomp)                          :: check

REAL(DP)                                             :: MonodTermTmp1
REAL(DP)                                             :: MonodTermTmp2
REAL(DP)                                             :: MonodTerm1
REAL(DP)                                             :: MonodTerm2
REAL(DP)                                             :: term1Perturb
REAL(DP)                                             :: term2Perturb
!!EAL(DP) , DIMENSION(ncomp)                          :: stmp
REAL(DP) , DIMENSION(ncomp)                          :: Checkjac_prekin
REAL(DP) , DIMENSION(ncomp,ncomp)                    :: sTMPperturb

REAL(DP)                                             :: sum
INTEGER(I4B)                                         :: kk
INTEGER(I4B)                                         :: ii
INTEGER(I4B)                                         :: idPoint

REAL(DP)                                             :: tmp1
REAL(DP)                                             :: tmp2
REAL(DP)                                             :: tmp3
REAL(DP)                                             :: tmp4

!!REAL(DP), DIMENSION(ikin)                            :: MoleFraction
!!REAL(DP), DIMENSION(ncomp,ikin)                      :: dMoleFraction

ALLOCATE(stmp(ncomp))

sppTMP = sp(:,jx,jy,jz)
sppTMP10 = sp10(:,jx,jy,jz)

tk = t(jx,jy,jz) + 273.15D0


! biomass end

!  First, calculate the derivatives of the affinity terms w/respect
!    to the primary species

!! Numerical derivatives of total concentration for use in Monod terms

sppTMP10(:) = sp10(:,jx,jy,jz)
perturb = 1.d-09 

DO i2 = 1,ncomp

  sppTMP(i2) = sppTMP(i2) + perturb
  sppTMP10(i2) = DEXP(sppTMP(i2))

!!  New speciation
  DO ksp = 1,nspec
    IF (muaq(ksp,i2) /= 0.0) THEN  !! If the secondary species is not affectd by primary i2, then skip it
      sum = 0.0D0
      DO ii = 1,ncomp
        sum = sum + muaq(ksp,ii)*(sppTMP(ii) + gam(ii,jx,jy,jz))
      END DO
      sppTMP(ksp+ncomp) = keqaq(ksp,jx,jy,jz) - gam(ksp+ncomp,jx,jy,jz) + sum
      sppTMP10(ksp+ncomp) = EXP(sppTMP(ksp+ncomp))  
    END IF
  END DO

  DO id = 1,nmonodaq(ir)
    i = imonodaq(id,ir)        !! Pointer to the primary species in the Monod expression (e.g., electron donor or acceptor)
    sum=0.0D0
    DO ksp = 1,nspec
      kk = ksp + ncomp
      sum = sum + muaq(ksp,i)*sppTMP10(kk)
    END DO
    sTMPperturb(i2,i) = sum + sppTMP10(i)
  END DO

  sppTMP(i2) = sppTMP(i2) - perturb
  sppTMP10(i2) = DEXP(sppTMP(i2))

END DO

end if   !  End of JennyDruhan block

rdkin = 0.0

DO ir = 1,ikin
  
  jac_sat = 0.0
  jac_prekin = 0.0
  Checkjac_prekin = 0.0
  
  IF (iaqtype(ir) == 3 .OR. iaqtype(ir) == 2 .OR. iaqtype(ir) == 4) THEN   !  Monod or irreversible
!!  IF (iaqtype(ir) == 3 .OR. iaqtype(ir) == 4) THEN   !  irreversible
    affinity = 1.0
  ELSE
    affinity = 1.0 - satkin(ir)
    DO i = 1,ncomp
      jac_sat(i) = -mukin(ir,i)*satkin(ir)
    END DO
  END IF
  
!  Calculate "pre-affinity" terms (far from equilibrium)
  
  IF (iaqtype(ir) == 1 .OR. iaqtype(ir) == 3 .OR. iaqtype(ir) == 4) THEN     !  TST or irreversible
    
    DO ll = 1,nreactkin(ir)
      DO i = 1,ncomp
        IF (ierode == 1) THEN
!!          denom = s(i,jx,jy,jz)                           !  Include aqueous 
          denom = s(i,jx,jy,jz) + sch(i,jx,jy,jz)/AqueousToBulk       !  Include aqueous and sorbed
        ELSE
!!          denom = s(i,jx,jy,jz)                           !  Include aqueous 
          denom = s(i,jx,jy,jz) + (sNCexch_local(i)+sNCsurf_local(i))/AqueousToBulk          !  Include aqueous and sorbed
        END IF
        IF (dependk(i,ll,ir) /= 0.0) THEN
          IF (itot(i,ll,ir) == 1) THEN
            IF (os3d) THEN
              DO i2 = 1,ncomp
!!                derivative = fjac_loc(i2,i) 
                derivative = fjac_loc(i2,i) + fch_local(i2,i)/AqueousToBulk
                 jac_prekin(i2,ll) =  jac_prekin(i2,ll) +  &
                   derivative*dependk(i,ll,ir)* pre_raq(ll,ir)/denom
              END DO
            ELSE
              IF (ierode == 1) THEN
                DO i2 = 1,ncomp
!!                  derivative = fjac(i2,i,jx,jy,jz) 
                  derivative = fjac(i2,i,jx,jy,jz) + fch(i2,i,jx,jy,jz)/AqueousToBulk
                   jac_prekin(i2,ll) =  jac_prekin(i2,ll) +  &
                     derivative*dependk(i,ll,ir)* pre_raq(ll,ir)/denom
                END DO
              ELSE
                DO i2 = 1,ncomp
!!                  derivative = fjac(i2,i,jx,jy,jz) 
                  derivative = fjac(i2,i,jx,jy,jz) + fch_local(i2,i)/AqueousToBulk
                   jac_prekin(i2,ll) =  jac_prekin(i2,ll) +  &
                     derivative*dependk(i,ll,ir)* pre_raq(ll,ir)/denom
                END DO
              END IF
            END IF
          ELSE
             jac_prekin(i,ll) =  jac_prekin(i,ll) + pre_raq(ll,ir)*dependk(i,ll,ir)
          END IF
        ELSE
          CONTINUE
        END IF
      END DO
    END DO
    
  ELSE IF (iaqtype(ir) == 2) THEN       ! Monod
      
!! Normal Monod terms

    DO id = 1,nmonodaq(ir)
      i = imonodaq(id,ir)
      IF (itot_monodaq(id,ir) == 1) THEN                 ! Dependence on total concentration
        denom = (s(i,jx,jy,jz)+halfsataq(id,ir)) * (s(i,jx,jy,jz)+halfsataq(id,ir))
        MonodTerm = s(i,jx,jy,jz)/(halfsataq(id,ir)+s(i,jx,jy,jz))
        DO i2 = 1,ncomp
          IF (os3d) THEN
             jac_prekin(i2,1) =  jac_prekin(i2,1) +  &
                pre_raq(1,ir)/MonodTerm * fjac_loc(i2,i)*halfsataq(id,ir)/denom

!!              fjac_loc(i2,i)*pre_raq(1,ir)* ( 1.0/s(i,jx,jy,jz) - 1.0/(s(i,jx,jy,jz)+halfsataq(id,ir))  )
!!               fjac_loc(i2,i)*halfsataq(id,ir)/denom
          ELSE
             jac_prekin(i2,1) =  jac_prekin(i2,1) +  &
                pre_raq(1,ir)/MonodTerm * fjac(i2,i,jx,jy,jz)*halfsataq(id,ir)/denom

!!              fjac(i2,i,jx,jy,jz)*pre_raq(1,ir)*( 1.0/s(i,jx,jy,jz) - 1.0/(s(i,jx,jy,jz)+halfsataq(id,ir))  )
!!              fjac(i2,i,jx,jy,jz)*halfsataq(id,ir)/denom
          END IF
        END DO
      ELSE
        denom = (sp10(i,jx,jy,jz)+halfsataq(id,ir))*(sp10(i,jx,jy,jz)+halfsataq(id,ir))
        MonodTerm = sp10(i,jx,jy,jz)/(halfsataq(id,ir)+sp10(i,jx,jy,jz))
        IF (i <= ncomp) THEN
!!           jac_prekin(i,1) =  jac_prekin(i,1) +  &
!!              pre_raq(1,ir)*( 1.0 - sp10(i,jx,jy,jz)/(sp10(i,jx,jy,jz)+halfsataq(id,ir)) )
!! Or
           jac_prekin(i,1) =  jac_prekin(i,1) +  &
              pre_raq(1,ir)/MonodTerm * sp10(i,jx,jy,jz)*halfsataq(id,ir)/denom

!!              sp10(i,jx,jy,jz)*halfsataq(id,ir)/denom
        ELSE
          ksp = i - ncomp
          DO i2 = 1,ncomp
             jac_prekin(i2,1) =  jac_prekin(i2,1) +  &
              pre_raq(1,ir)/MonodTerm * muaq(ksp,i2)*sp10(i,jx,jy,jz)*halfsataq(id,ir)/denom

!!                pre_raq(1,ir)* muaq(ksp,i2)* (1.0 - sp10(i,jx,jy,jz)/(sp10(i,jx,jy,jz)+halfsataq(id,ir)) )
          END DO
        END IF
      END IF
    END DO

!!  Inhibition terms

    DO id = 1,ninhibitaq(ir)
      i = inhibitaq(id,ir)
      IF (inhibitaq(id,ir) < 0) THEN
        CONTINUE                                             !!  No dependence on mineral volume fraction
      ELSE
        IF (itot_inhibitaq(id,ir) == 1) THEN                 !!  Dependence on total concentration
          denom = ( rinhibitaq(id,ir)+s(i,jx,jy,jz) )* ( rinhibitaq(id,ir)+s(i,jx,jy,jz) )
          InhibitTerm = rinhibitaq(id,ir)/(rinhibitaq(id,ir)+s(i,jx,jy,jz))
          DO i2 = 1,ncomp
            IF (os3d) THEN
              jac_inhibition = -fjac_loc(i2,i)*rinhibitaq(id,ir)/denom
              jac_prekin(i2,1) =  jac_prekin(i2,1) + pre_raq(1,ir)/InhibitTerm * jac_inhibition     
            ELSE
              jac_inhibition = -fjac(i2,i,jx,jy,jz)*rinhibitaq(id,ir)/denom
              jac_prekin(i2,1) =  jac_prekin(i2,1) + pre_raq(1,ir)/InhibitTerm * jac_inhibition  
            END IF
          END DO
        ELSE
          denom = ( rinhibitaq(id,ir)+sp10(i,jx,jy,jz) )*( rinhibitaq(id,ir)+sp10(i,jx,jy,jz) )
          InhibitTerm = rinhibitaq(id,ir)/(rinhibitaq(id,ir)+sp10(i,jx,jy,jz))
          IF (i <= ncomp) THEN
            jac_inhibition = -sp10(i,jx,jy,jz)*rinhibitaq(id,ir)/denom
            jac_prekin(i,1) =  jac_prekin(i,1) + pre_raq(1,ir)/InhibitTerm * jac_inhibition 
          ELSE
            ksp = i - ncomp
            DO i2 = 1,ncomp
              jac_inhibition = -muaq(ksp,i2)*sp10(i,jx,jy,jz)*rinhibitaq(id,ir)/denom
              jac_prekin(i2,1) =  jac_prekin(i2,1) + pre_raq(1,ir)/InhibitTerm * jac_inhibition               
            END DO
          END IF
        END IF
      END IF

     END DO

! biomass
    ELSE IF (iaqtype(ir) == 8) THEN       ! Monod, but with thermodynamic factor, F_T
      
!! Normal Monod terms

    DO id = 1,nmonodaq(ir)
      i = imonodaq(id,ir)        !! Pointer to the primary species in the Monod expression (e.g., electron donor or acceptor)

      IF (IsotopePrimaryCommon(i)) THEN

        if (os3d) then
          write(*,*)
          write(*,*) ' Maggi Formulation not set up yet for OS3D option'
          write(*,*)
          read(*,*)
          stop
        end if

        IsotopologueOther = isotopeRare(iPointerIsotope(i))
        IF (itot_monodaq(id,ir) == 1) THEN                 ! Dependence on total concentration
          MonodTerm(id) = s(i,jx,jy,jz)/( s(i,jx,jy,jz)+ halfsataq(id,ir)*(1.0d0+s(IsotopologueOther,jx,jy,jz)/halfsataq(id,ir)) )
        ELSE                                               ! Dependence on individual species
          MonodTerm(id) = sp10(i,jx,jy,jz)/( sp10(i,jx,jy,jz)+ halfsataq(id,ir)*(1.0d0+sp10(IsotopologueOther,jx,jy,jz)/halfsataq(id,ir)) )
        END IF

      ELSE IF (IsotopePrimaryRare(i)) THEN

        IsotopologueOther = isotopeCommon(iPointerIsotope(i))
        IF (itot_monodaq(id,ir) == 1) THEN                 ! Dependence on total concentration
          MonodTerm(id) = s(i,jx,jy,jz)/( s(i,jx,jy,jz)+ halfsataq(id,ir)*(1.0d0+s(IsotopologueOther,jx,jy,jz)/halfsataq(id,ir)) )
        ELSE                                               ! Dependence on individual species
          MonodTerm(id) = sp10(i,jx,jy,jz)/( sp10(i,jx,jy,jz)+ halfsataq(id,ir)*(1.0d0+sp10(IsotopologueOther,jx,jy,jz)/halfsataq(id,ir)) )
        END IF

      ELSE    !general case - no isotopes
        IF (itot_monodaq(id,ir) == 1) THEN                 ! Dependence on total concentration
          denom = (s(i,jx,jy,jz)+halfsataq(id,ir)) * (s(i,jx,jy,jz)+halfsataq(id,ir))
          MonodTerm(id) = s(i,jx,jy,jz)/(halfsataq(id,ir)+s(i,jx,jy,jz))
        ELSE
          denom = (sp10(i,jx,jy,jz)+halfsataq(id,ir)) * (sp10(i,jx,jy,jz)+halfsataq(id,ir))
          MonodTerm(id) = sp10(i,jx,jy,jz)/(halfsataq(id,ir)+sp10(i,jx,jy,jz))
        END IF
      END IF

      DO i2 = 1,ncomp

        IF (os3d) THEN

          jac_prekin(i2,1) =  jac_prekin(i2,1) +  &
                pre_raq(1,ir)/MonodTerm * fjac_loc(i2,i)*halfsataq(id,ir)/denom

        ELSE


          IF (IsotopePrimaryCommon(i)) THEN


            IsotopologueOther = isotopeRare(iPointerIsotope(i))
            IF (itot_monodaq(id,ir) == 1) THEN                 ! Dependence on total concentration
              MonodTerm(id) = s(i,jx,jy,jz)/( s(i,jx,jy,jz)+ halfsataq(id,ir)*(1.0d0+s(IsotopologueOther,jx,jy,jz)/halfsataq(id,ir)) )
            ELSE                                               ! Dependence on individual species
              MonodTerm(id) = sp10(i,jx,jy,jz)/( sp10(i,jx,jy,jz)+ halfsataq(id,ir)*(1.0d0+sp10(IsotopologueOther,jx,jy,jz)/halfsataq(id,ir)) )
            END IF

          ELSE IF (IsotopePrimaryRare(i)) THEN

            IsotopologueOther = isotopeCommon(iPointerIsotope(i))
            IF (itot_monodaq(id,ir) == 1) THEN                 ! Dependence on total concentration
              MonodTerm(id) = s(i,jx,jy,jz)/( s(i,jx,jy,jz)+ halfsataq(id,ir)*(1.0d0+s(IsotopologueOther,jx,jy,jz)/halfsataq(id,ir)) )
            ELSE                                               ! Dependence on individual species
              MonodTerm(id) = sp10(i,jx,jy,jz)/( sp10(i,jx,jy,jz)+ halfsataq(id,ir)*(1.0d0+sp10(IsotopologueOther,jx,jy,jz)/halfsataq(id,ir)) )
            END IF

          ELSE    !general case - no isotopes
            IF (itot_monodaq(id,ir) == 1) THEN                 ! Dependence on total concentration
              denom = (s(i,jx,jy,jz)+halfsataq(id,ir)) * (s(i,jx,jy,jz)+halfsataq(id,ir))
              MonodTerm(id) = s(i,jx,jy,jz)/(halfsataq(id,ir)+s(i,jx,jy,jz))
            ELSE
              denom = (sp10(i,jx,jy,jz)+halfsataq(id,ir)) * (sp10(i,jx,jy,jz)+halfsataq(id,ir))
              MonodTerm(id) = sp10(i,jx,jy,jz)/(halfsataq(id,ir)+sp10(i,jx,jy,jz))
            END IF
          END IF


!!
              if (i == 10) then    !! First, do the 32Sulfate (i=10)
                if (maggi) then

                  MonodTermTmp2 = sTMPperturb(i2,10)/( sTMPperturb(i2,10) + halfsataq(2,1)*(1.0d0+sTMPperturb(i2,11)/halfsataq(2,2)) )
                  Term2Perturb = (MonodTermTmp2 - MonodTerm2)/perturb

                  Checkjac_prekin(i2) = Checkjac_prekin(i2) + termMonod(1,1)*Term2Perturb
                else
                  MonodTermTmp1 = sTMPperturb(i2,10)/( sTMPperturb(i2,10) + halfsataq(2,1)  )
                  Term1Perturb = (MonodTermTmp1 - MonodTerm1)/perturb

                  Checkjac_prekin(i2) = Checkjac_prekin(i2) + termMonod(1,1)*Term1Perturb
                end if

              else if (i==11) then   !! Then, do the 34Sulfate (i=11)

                if (maggi) then
                  MonodTermTmp2 = sTMPperturb(i2,11)/( sTMPperturb(i2,11) + halfsataq(2,2)*(1.0d0+sTMPperturb(i2,10)/halfsataq(2,1)) )
                  Term2Perturb = (MonodTermTmp2 - MonodTerm2)/perturb

                  Checkjac_prekin(i2) = Checkjac_prekin(i2) + termMonod(1,2)*Term2Perturb
                else
                  MonodTermTmp1 = sTMPperturb(i2,11)/( sTMPperturb(i2,11) + halfsataq(2,2) )
                  Term1Perturb = (MonodTermTmp1 - MonodTerm1)/perturb

                  Checkjac_prekin(i2) = Checkjac_prekin(i2) + termMonod(1,2)*Term1Perturb
                end if

              else if (i == 15) then  !! Then, do acetate

                MonodTermTmp1 = sTMPperturb(i2,15)/(halfsataq(id,ir) + sTMPperturb(i2,15))
                Term1Perturb = (MonodTermTmp1 - MonodTerm1)/perturb

                Checkjac_prekin(i2) = Checkjac_prekin(i2) + termMonod(2,ir)*Term1Perturb

              else
                write(*,*)
                write(*,*) ' There should be no other options in JACRKIN for the aqueous reactions--Check code'
                read(*,*)
                stop
              endif


!!               write(*,*) ' Kinetic reaction: ',ir
!!               write(*,*) ' Dependence on: ',ulab(i2)
!!               write(*,*) jac_prekin(i2,1),Checkjac_prekin(i2)
!!               read(*,*)

              jac_prekin(i2,1) = Checkjac_prekin(i2)

          END IF
        END DO

      ELSE
        denom = (sp10(i,jx,jy,jz)+halfsataq(id,ir))*(sp10(i,jx,jy,jz)+halfsataq(id,ir))
        MonodTerm = sp10(i,jx,jy,jz)/(halfsataq(id,ir)+sp10(i,jx,jy,jz))
        IF (i <= ncomp) THEN
!!           jac_prekin(i,1) =  jac_prekin(i,1) +  &
!!              pre_raq(1,ir)*( 1.0 - sp10(i,jx,jy,jz)/(sp10(i,jx,jy,jz)+halfsataq(id,ir)) )
!! Or
           jac_prekin(i,1) =  jac_prekin(i,1) +  &
              pre_raq(1,ir)/MonodTerm * sp10(i,jx,jy,jz)*halfsataq(id,ir)/denom

!!              sp10(i,jx,jy,jz)*halfsataq(id,ir)/denom
        ELSE
          ksp = i - ncomp
          DO i2 = 1,ncomp
             jac_prekin(i2,1) =  jac_prekin(i2,1) +  &
              pre_raq(1,ir)/MonodTerm * muaq(ksp,i2)*sp10(i,jx,jy,jz)*halfsataq(id,ir)/denom

!!                pre_raq(1,ir)* muaq(ksp,i2)* (1.0 - sp10(i,jx,jy,jz)/(sp10(i,jx,jy,jz)+halfsataq(id,ir)) )
          END DO
        END IF
      END IF
    END DO

!! add inhibition 2011-08-15
!!  Inhibition terms

    DO id = 1,ninhibitaq(ir)
      i = inhibitaq(id,ir)
      IF (inhibitaq(id,ir) < 0) THEN
        CONTINUE                                             !!  No dependence on mineral volume fraction
      ELSE
        IF (itot_inhibitaq(id,ir) == 1) THEN                 !!  Dependence on total concentration
          denom = ( rinhibitaq(id,ir)+s(i,jx,jy,jz) )* ( rinhibitaq(id,ir)+s(i,jx,jy,jz) )
          InhibitTerm = rinhibitaq(id,ir)/(rinhibitaq(id,ir)+s(i,jx,jy,jz))
          DO i2 = 1,ncomp
            IF (os3d) THEN
              jac_inhibition = -fjac_loc(i2,i)*rinhibitaq(id,ir)/denom
              jac_prekin(i2,1) =  jac_prekin(i2,1) + pre_raq(1,ir)/InhibitTerm * jac_inhibition     
            ELSE
              jac_inhibition = -fjac(i2,i,jx,jy,jz)*rinhibitaq(id,ir)/denom
              jac_prekin(i2,1) =  jac_prekin(i2,1) + pre_raq(1,ir)/InhibitTerm * jac_inhibition  
            END IF
          END DO
        ELSE
          denom = ( rinhibitaq(id,ir)+sp10(i,jx,jy,jz) )*( rinhibitaq(id,ir)+sp10(i,jx,jy,jz) )
          InhibitTerm = rinhibitaq(id,ir)/(rinhibitaq(id,ir)+sp10(i,jx,jy,jz))
          IF (i <= ncomp) THEN
            jac_inhibition = -sp10(i,jx,jy,jz)*rinhibitaq(id,ir)/denom
            jac_prekin(i,1) =  jac_prekin(i,1) + pre_raq(1,ir)/InhibitTerm * jac_inhibition 
          ELSE
            ksp = i - ncomp
            DO i2 = 1,ncomp
              jac_inhibition = -muaq(ksp,i2)*sp10(i,jx,jy,jz)*rinhibitaq(id,ir)/denom
              jac_prekin(i2,1) =  jac_prekin(i2,1) + pre_raq(1,ir)/InhibitTerm * jac_inhibition               
            END DO
          END IF
        END IF
      END IF

     END DO
!end add

!! No inhibition terms here, so calculate thermodynamic factor, F_T

    IF (direction_kin(ir) < 0) THEN
      sign = 1.0d0
    ELSE
      sign = -1.0d0
    END IF

    IF(satkin(ir) > 1.0d0) THEN
      snormAqueous = 1.0d0
    ELSE 
      snormAqueous = satkin(ir)
    ENDIF

    term1 = sign*DABS(snormAqueous - 1.0D0)

!!  Reaction assumed to be irreversible, so do not let it go in reverse

    affinity = term1
    
!!    jj = p_cat_kin(ir)
    jj = ir
    check(:) = -mukinTMP(jj,:)*snormAqueous
    jac_sat = check
! biomass end

  ELSE
    WRITE(*,*)
    WRITE(*,*) ' No other rate laws recognized in JACRKIN.F90 '
    READ(*,*)
    STOP
  END IF

! biomass
  if (iaqtype(ir) == 8) then

!!    jj = p_cat_kin(ir)
    jj = ir
!   pointer to biomass for current reaction
!!    ib = ibiomass_kin(p_cat_kin(ir))
    ib = ibiomass_kin(ir)
    
    IF (UseMetabolicLagAqueous(jj)) THEN
      DO i = 1,ncomp
        DO ll = 1,nreactkin(ir)
          rdkin(ir,i) = rdkin(ir,i) + MetabolicLagAqueous(jj,jx,jy,jz)*volfx(ib,jx,jy,jz) * ratek(ll,ir)*  &
            (MoleFraction(ir)*pre_raq(ll,ir)*jac_sat(i) +  MoleFraction(ir)*jac_prekin(i,ll)*affinity + pre_raq(ll,ir)*affinity*dMoleFraction(i,ir) )
        END DO
      END DO
    ELSE
      DO i = 1,ncomp
        DO ll = 1,nreactkin(ir)
          rdkin(ir,i) = rdkin(ir,i) + volfx(ib,jx,jy,jz) * ratek(ll,ir)*  &
            (MoleFraction(ir)*pre_raq(ll,ir)*jac_sat(i) +  MoleFraction(ir)*jac_prekin(i,ll)*affinity + pre_raq(ll,ir)*affinity*dMoleFraction(i,ir) )
        END DO
      END DO
    END IF
    
  else 

    DO i = 1,ncomp
      DO ll = 1,nreactkin(ir)
        rdkin(ir,i) = rdkin(ir,i) + ratek(ll,ir)*  &
            (pre_raq(ll,ir)*jac_sat(i) +  jac_prekin(i,ll)*affinity )
      END DO
    END DO

  end if
! biomass end
  
END DO

DEALLOCATE(stmp)

RETURN
END SUBROUTINE jacrkin

