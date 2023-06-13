!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 10:02:15
 
!************** (C) COPYRIGHT 1995,1998,1999 ******************
!*******************     C.I. Steefel      *******************
!                    All Rights Reserved

!  GIMRT98 IS PROVIDED "AS IS" AND WITHOUT ANY WARRANTY EXPRESS OR IMPLIED.
!  THE USER ASSUMES ALL RISKS OF USING GIMRT98. THERE IS NO CLAIM OF THE
!  MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.

!  YOU MAY MODIFY THE SOURCE CODE FOR YOUR OWN USE, BUT YOU MAY NOT
!  DISTRIBUTE EITHER THE ORIGINAL OR THE MODIFIED CODE TO ANY OTHER
!  WORKSTATIONS
!**********************************************************************

SUBROUTINE reaction(ncomp,nkin,nrct,nspec,nexchange,nsurf,ndecay,jx,jy,jz,delt)
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
INTEGER(I4B), INTENT(IN)                                        :: ndecay
INTEGER(I4B), INTENT(IN)                                        :: jx
INTEGER(I4B), INTENT(IN)                                        :: jy
INTEGER(I4B), INTENT(IN)                                        :: jz

REAL(DP), INTENT(IN)                                            :: delt

!  Internal variables and arrays

REAL(DP)                                                        :: tk
REAL(DP)                                                        :: tkinv
REAL(DP)                                                        :: reft
REAL(DP)                                                        :: porfactor
REAL(DP)                                                        :: sumiap
REAL(DP)                                                        :: term2
REAL(DP)                                                        :: term1
REAL(DP)                                                        :: vcheck 
REAL(DP)                                                        :: check 
REAL(DP)                                                        :: checkmonod
REAL(DP)                                                        :: sign

INTEGER(I4B)                                                    :: k
INTEGER(I4B)                                                    :: np
INTEGER(I4B)                                                    :: i
INTEGER(I4B)                                                    :: is
INTEGER(I4B)                                                    :: kk
INTEGER(I4B)                                                    :: nex
INTEGER(I4B)                                                    :: ixx

REAL(DP)                                                        :: RateFactor
REAL(DP)                                                        :: CubicTerm
REAL(DP)                                                        :: CorrectedTerm
REAL(DP)                                                        :: Kformation
REAL(DP)                                                        :: Al

REAL(DP)                                                        :: dependTMP
REAL(DP)                                                        :: surftmp
REAL(DP)                                                        :: MinConvert
REAL(DP)                                                        :: CheckRgas

!  This routine calculates the reaction rates of the individual
!  minerals.

tk = t(jx,jy,jz) + 273.15D0
tkinv = 1.0D0/tk
reft = 1.0D0/298.15D0
porfactor = (por(jx,jy,jz)/porin(jx,jy,jz))**(0.666666D0)

decay_correct = 1.0D0
snorm = 0.0d0

!DO id = 1,ndecay
!  i = idecay(id)
!  DO kd = 1,nmindecay(id)
!    k = kdecay(kd,id)
!    sum = 0.0
!    DO isotope = 1,nisotope(id)
!      sum = sum + ratio_isotope(isotope,kd,id,jx,jy)
!    END DO
!    decay_correct(i,k) = sum
!  END DO
!END DO

ivolume = 0

DO k = 1,nkin
  dppt(k,jx,jy,jz) = 0.0D0
  
  DO np = 1,nreactmin(k)
    
    IF (t(jx,jy,jz) == 25.0D0) THEN
      actenergy(np,k) = 1.0D0
    ELSE
      actenergy(np,k) = DEXP( (ea(np,k)/rgasKCAL)*(reft-tkinv) )
      CheckRgas = rgasKCAL
    END IF
    
    IF (KateMaher) THEN
      mumin(np,kUCalcite,ikCa) = muCalcium(jx,jy,jz)
      mumin(np,kUCalcite,ik234U) = muUranium234(jx,jy,jz)
      mumin(np,kUCalcite,ik238U) = muUranium238(jx,jy,jz)
    END IF

!!  Imintype = 1 --> TST
!!  Imintype = 2 --> Monod kinetics (NOTE: Should add the affinity term here)
!!  Imintype = 3 --> Irreversible kinetics
!!  Imintype = 4 --> Precipitation only

    IF (imintype(np,k) == 3 .OR. imintype(np,k) == 2) THEN
      silog(np,k) = -500.0D0
      si(np,k) = 0.0D0
    ELSE

      IF (kcrossaff(np,k) /= 0) THEN
        kk = kcrossaff(np,k)
        sumiap = 0.0D0
        DO i = 1,ncomp
          sumiap = sumiap + decay_correct(i,k)*mumin(1,kk,i)*(sp(i,jx,jy,jz)+gam(i,jx,jy,jz))
        END DO
        silog(np,k) = (sumiap - keqmin(1,kk,jx,jy,jz))/clg
        si(np,k) = 10**(silog(np,k))
      ELSE
        sumiap = 0.0D0
        DO i = 1,ncomp
          sumiap = sumiap + decay_correct(i,k)*mumin(np,k,i)*(sp(i,jx,jy,jz)+gam(i,jx,jy,jz))
        END DO
        silog(np,k) = (sumiap - keqmin(np,k,jx,jy,jz))/clg
        siln(np,k) = clg*silog(np,k)
        si(np,k) = 10**(silog(np,k))
      END IF
    END IF
    
!  Decide on the value of the reactive surface area depending on
!  whether the mineral is under or super-saturated and whether
!  it has a zero volume fraction.
    
    IF (silog(np,k) >= 0.0D0 .AND. iarea(k,jinit(jx,jy,jz)) == 0) THEN
      surf(k) = areain(k,jinit(jx,jy,jz))*porfactor
      IF (kcrossaff(np,k) /= 0) THEN
        surf(k) = 0.0D0
      END IF
    ELSE
      IF (porfactor < 0.01d0) THEN
        IF (volfx(k,jx,jy,jz) > 1.d-14) THEN
          surf(k) = area(k,jx,jy,jz)*porfactor
        ELSE
          surf(k) = 0.0D0
        END IF
      ELSE
        IF (volfx(k,jx,jy,jz) > 1.d-14) THEN
          surf(k) = area(k,jx,jy,jz)
        ELSE
          surf(k) = 0.0D0
        END IF
      END IF
    END IF

    IF (imintype(np,k) == 3 .OR. imintype(np,k) == 2) THEN
      surf(k) = area(k,jx,jy,jz)
    END IF

    IF (KateMaher) THEN
      IF (k==KUcalcite .AND. rlabel(np,k) == 'recrystallize') THEN
        surf(k) = area(kMarineCalcite,jx,jy,jz)
      END IF
    END IF
    
    IF (surf(k) == 0.0D0) THEN
      rmin(np,k) = 0.0D0
      dppt(k,jx,jy,jz) = 0.0D0
      pre_rmin(np,k) = 1.0D0
    ELSE
          
      IF (imintype(np,k) == 1 .OR. imintype(np,k) == 3 .OR. imintype(np,k) == 4) THEN         !! TST, irreversible, or ppt only
        
        term2 = 0.0D0
        DO kk = 1,ndepend(np,k)
          i = idepend(kk,np,k)
          IF (depend(kk,np,k) < 0.000) THEN
            IF (itot_min(kk,np,k) == 1) THEN
              term2 = term2 + depend(kk,np,k)*DLOG(s(i,jx,jy,jz))
            ELSE
              term2 = term2 + depend(kk,np,k)*(gam(i,jx,jy,jz)+sp(i,jx,jy,jz))
            END IF
          ELSE
            IF (itot_min(kk,np,k) == 1) THEN
              term2 = term2 + depend(kk,np,k)*DLOG(s(i,jx,jy,jz))
            ELSE
              term2 = term2 + depend(kk,np,k)*(gam(i,jx,jy,jz)+sp(i,jx,jy,jz))
            END IF
          END IF
        END DO
        DO kk = 1,ndependex(np,k)
          nex = ixdepend(kk,np,k)
          ixx = ixlink(nex)
!          term2 = term2 + dependex(kk,np,k)*spex(nex+nexchange,jx,jy,jz)/exchangesites(ixx,jx,jy,jz)
          term2 = term2 + dependex(kk,np,k)*spex(nex+nexchange,jx,jy,jz)
        END DO
        DO kk = 1,ndependsurf(np,k)
          is = isdepend(kk,np,k)
          dependTMP = dependsurf(kk,np,k)
          IF (is <= nsurf) THEN
            term2 = term2 + dependsurf(kk,np,k)*spsurf(is,jx,jy,jz)
          ELSE
            term2 = term2 + dependTMP*DLOG(spsurf10(is,jx,jy,jz))        
          END IF
        END DO

        IF (k == kUPlag .AND. OelkersRateLaw) THEN     !  Oelkerian rate law

          Kformation = 0.0000149823D0
          Al = sp10(ikAl,jx,jy,jz)
          term2 = Kformation/(Kformation + Al**0.333333)
        END IF
        
!        ***** * Monod kinetics **********
        
      ELSE IF (imintype(np,k) == 2) THEN     !  Monod rate law

        IF (volfx(k,jx,jy,jz) > 0.0d0) THEN
          surf(k) = 1.0d0                                                        !! Resets surface area to a constant 1
        ELSE
          surf(k) = 0.0d0
        END IF

        term2 = 1.0D0
        DO kk = 1,nmonod(np,k)
          IF (imonod(kk,np,k) == 0 .AND. kmonod(kk,np,k) == 1) THEN                !! Mineral Monod reaction depends on its own concentration
            MinConvert = volfx(k,jx,jy,jz)/(volmol(k)*por(jx,jy,jz)*ro(jx,jy,jz))  !! Converts mineral volume fraction to moles mineral per kg fluid (molality)                                  
            checkmonod =  MinConvert/(MinConvert+halfsat(kk,np,k))
            term2 = term2 * checkmonod 
          ELSE                                                                     !! Monod term depends on aqueous species
            i = imonod(kk,np,k)
            IF (itot_monod(kk,np,k) == 1) THEN
              checkmonod =  s(i,jx,jy,jz)/(s(i,jx,jy,jz)+halfsat(kk,np,k))
              term2 = term2 * checkmonod             
            ELSE
              checkmonod =  sp10(i,jx,jy,jz)/(sp10(i,jx,jy,jz)+halfsat(kk,np,k))
              term2 = term2 * checkmonod
            END IF
          END IF
        END DO
        DO kk = 1,ninhibit(np,k)
          i = inhibit(kk,np,k)
          IF (itot_inhibit(kk,np,k) == 1) THEN
            term2 = term2 * rinhibit(kk,np,k)/(rinhibit(kk,np,k)+s(i,jx,jy,jz))
          ELSE
            term2 = term2 * rinhibit(kk,np,k)/(rinhibit(kk,np,k)+sp10(i,jx,jy,jz))
          END IF
        END DO

      ELSE
        WRITE(*,*)
        WRITE(*,*) ' Rate formulation not recognized'
        WRITE(*,*) ' Imintype = ', imintype(np,k)
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
      
      IF (imintype(np,k) == 2) THEN     ! Monod kinetics
        pre_rmin(np,k) = term2
      ELSE IF (k == kUPlag .AND. OelkersRateLaw) THEN
        pre_rmin(np,k) = term2
      ELSE
        IF (term2 == 0.0D0) THEN
          pre_rmin(np,k) = 1.0D0
        ELSE
          pre_rmin(np,k) = DEXP(term2)
        END IF
      END IF
      
!***************** Dependence on mineral saturation state *************
     
      IF (si(np,k) > 1.0D0) THEN
        sign = 1.0D0
      ELSE
        sign = -1.0D0
      END IF

      IF (KateMaher) THEN
        IF (k == KUCalcite .AND. rlabel(np,k) == 'recrystallize') THEN
          sign = 1.0D0
        END IF
      END IF

      IF (imintype(np,k) == 3 .OR. imintype(np,k) == 2) THEN     !!  Monod or irreversible
        snorm(np,k) = 0.0D0
      ELSE
        IF (AffinityDepend2(np,k) == 1.0D0 .AND. AffinityDepend3(np,k) == 1.0d0) THEN
          snorm(np,k) = si(np,k)
        ELSE IF (AffinityDepend2(np,k) /= 1.0d0 .AND. AffinityDepend3(np,k) == 1.0d0) THEN
          snorm(np,k) = si(np,k)**AffinityDepend2(np,k)
        ELSE
          snorm(np,k) = DEXP(-AffinityDepend2(np,k)*(DABS(siln(np,k)))**AffinityDepend3(np,k) )
        END IF
      END IF 

      IF (imintype(np,k) == 4) THEN                        !! Precipitation only rate law
        term1 = snorm(np,k)
      ELSE                                                 !!  TST forms
        IF (AffinityDepend1(np,k) == 1.0D0) THEN
          term1 = sign*ABS(snorm(np,k) - 1.0D0)
        ELSE
          term1 = sign*(ABS(snorm(np,k) - 1.0D0))**(AffinityDepend1(np,k))
        END IF
      END IF

! --> Following is for linear TST kinetics where precipitation is inhibited by setting a threshold

      IF (thresh(np,k) > 0.0D0 .AND. AffinityDepend1(np,k) == 1.0d0 .AND. AffinityDepend2(np,k) == 1.0d0) THEN
        term1 = sign
      END IF

      IF (thresh(np,k) > 0.0d0 .AND. silog(np,k) > 0.0d0) THEN
        snorm(np,k) = 0.0d0
      END IF

      IF (thresh(np,k) > 0.0D0 .AND. AffinityDepend1(np,k) == 1.0d0 .AND. AffinityDepend2(np,k) == 1.0d0) THEN
        IF (silog(np,k) < -2.0D0*thresh(np,k)) THEN
           RateFactor = 1.0D0
        ELSE IF (silog(np,k) >= 0.0D0) THEN
           RateFactor = 0.0D0
        ELSE 
          CorrectedTerm = silog(np,k)+thresh(np,k)
          CubicTerm = CorrectedTerm*CorrectedTerm*CorrectedTerm/   &
                        ( thresh(np,k)*thresh(np,k)*thresh(np,k) )
          RateFactor = 0.5D0 - 0.75D0*CorrectedTerm/thresh(np,k) + 0.25D0*CubicTerm
        END IF
      ELSE
        IF (thresh(np,k) > 0.0d0 .AND. silog(np,k) > 0.0d0) THEN
          RateFactor = 0.0d0
        ELSE
          RateFactor = 1.0d0
        END IF
      END IF

      IF (imintype(np,k) == 4 .AND. silog(np,k) < 0.0d0) THEN
        RateFactor = 0.0d0
      END IF

!!  *****************  End of linear TST suppression of precipitation  ************************

!************** End of dependence on mineral saturation state *************
      
      rmin(np,k) = RateFactor*surf(k)*term1*rate0(np,k)*pre_rmin(np,k)*actenergy(np,k)
      
      dppt(k,jx,jy,jz) = dppt(k,jx,jy,jz) + rmin(np,k)
      
    END IF
    
  END DO   !  End of npth parallel reaction
  
!!  vcheck = volfx(k,jx,jy,jz) + dppt(k,jx,jy,jz)*volmol(k)*delt
  
!!  IF (vcheck < 0.0) THEN
!!    dppt(k,jx,jy,jz) = -volfx(k,jx,jy,jz)/(volmol(k)*delt)
!!    rmin(1,k) = dppt(k,jx,jy,jz)
!!    IF (nreactmin(k) > 1) THEN
!!      DO np = 2,nreactmin(k)
!!        rmin(np,k) = 0.0
!!      END DO
!!    END IF
!!    ivolume(k) = 1

!    WRITE(*,*)
!    WRITE(*,*) ' Mineral reaction rate exceeds available mineral mass'
!    WRITE(*,*) ' Setting rate to match mineral mass and moving rate to R.H.S.'
!    WRITE(*,*)
    
!!  END IF
  
END DO     !  End of kth mineral

RETURN
END SUBROUTINE reaction
!********************************************************

!!  For Colleen Hansel modeling
!!        if (k == 5) then
!!           surf(k) = 100.0
!!        end if
!!        IF (k == 2) THEN              !  Magnetite
!!          surf(k) = 0.001*area(5,jx,jy,jz) + area(2,jx,jy,jz) + 0.0001*area(4,jx,jy,jz)    !  Use surface area of green rust
!!        END IF
!!        if (k == 3) then
!!          surf(k) = area(3,jx,jy,jz) !! + 1.0*area(1,jx,jy,jz) + 0.0001*area(4,jx,jy,jz)
!!        end if
!!        if (k == 4) then
!!          surf(k) = area(4,jx,jy,jz) !!+ area(1,jx,jy,jz)
!!        end if
