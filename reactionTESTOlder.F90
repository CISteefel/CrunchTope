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

SUBROUTINE reaction(ncomp,nkin,nrct,nspec,nexchange,nsurf,ndecay,jx,jy,jz,delt,time)
USE crunchtype
USE params
USE concentration
USE mineral
USE medium
USE temperature
USE runtime
USE NanoCrystal

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
REAL(DP), INTENT(IN)                                            :: time

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
INTEGER(I4B)                                                    :: l
INTEGER(I4B)                                                    :: kJapaneseGlassNoSilica

REAL(DP)                                                        :: RateFactor
REAL(DP)                                                        :: CubicTerm
REAL(DP)                                                        :: CorrectedTerm
!!REAL(DP)                                                        :: Kformation
REAL(DP)                                                        :: Al

REAL(DP)                                                        :: dependTMP
REAL(DP)                                                        :: surftmp
REAL(DP)                                                        :: MinConvert
REAL(DP)                                                        :: CheckRgas
REAL(DP)                                                        :: AffinityTerm

REAL(DP)                                                        :: Astar
REAL(DP)                                                        :: Bstar
REAL(DP)                                                        :: Sstar
REAL(DP)                                                        :: denominator
REAL(DP)                                                        :: MoleFraction40
REAL(DP)                                                        :: MoleFraction44
REAL(DP)                                                        :: MoleFraction32
REAL(DP)                                                        :: MoleFraction34
REAL(DP)                                                        :: MoleFraction32S
REAL(DP)                                                        :: MoleFraction34S

REAL(DP)                                                        :: RateFraction
REAL(DP)                                                        :: MoleFraction44Mineral
REAL(DP)                                                        :: MoleFraction40Mineral
REAL(DP)                                                        :: CurrentMoleFraction40
REAL(DP)                                                        :: CurrentMoleFraction44
REAL(DP)                                                        :: CheckMoleFraction
REAL(DP)                                                        :: CheckLogQ
REAL(DP)                                                        :: CheckLogQ10
REAL(DP)                                                        :: kBackward
REAL(DP)                                                        :: kForward
REAL(DP)                                                        :: rminCheck

REAL(DP)                                                        :: CalciumCarbonateRatioEffect

! biomass
INTEGER(I4B)                                                    :: ib, jj
! end biomass

!! CrystalSizeDistribution
REAL(DP)                                                        :: SizeDepKeq
REAL(DP)                                                        :: AdjustKeq

rmin = 0.0d0

!!CSD    sigma(3) = 200.0/1000.0
!!CSD    sigma(4) = 90.0/1000.0

!  This routine calculates the reaction rates of the individual
!  minerals.

tk = t(jx,jy,jz) + 273.15D0
tkinv = 1.0D0/tk
reft = 1.0D0/298.15D0
porfactor = (por(jx,jy,jz)/porin(jx,jy,jz))**(0.666666D0)


If (JennyDruhan) THEN       !! **********************************************************

!!CIS Hardwired for 44Ca40CaCO3
  MoleFraction44 = sn(7,jx,jy,jz)/( sn(6,jx,jy,jz) + sn(7,jx,jy,jz) )   !! NOTE: Using old time step total concentrations (lagged)
  MoleFraction40 = 1.0d0 - MoleFraction44

  if (ulab(7) /= 'Ca44++') then
    write(*,*)
    write(*,*) ' Primary species 7 should be Ca44++'
    write(*,*)
    read(*,*)
    stop
  end if    
  if (ulab(6) /= 'Ca++') then
    write(*,*)
    write(*,*) ' Primary species 6 should be Ca++'
    write(*,*)
    read(*,*)
    stop
  end if 

  IF (MoleFraction40 == 0.0d0) THEN
    WRITE(*,*)
    WRITE(*,*) ' MoleFraction40 = 0'
    write(*,*)
    stop
  END IF
  IF (MoleFraction44 == 0.0d0) THEN
    WRITE(*,*)
    WRITE(*,*) ' MoleFraction44 = 0'
    write(*,*)
    stop
  END IF

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

  MoleFraction34 = sp10(13,jx,jy,jz)/( sp10(12,jx,jy,jz) + sp10(13,jx,jy,jz) )
  MoleFraction32 = 1.0d0 - MoleFraction34
  if (ulab(13) /= 'H2S34(aq)') then
    write(*,*)
    write(*,*) ' Primary species 13 should be H2S34(aq)+'
    write(*,*)
    read(*,*)
    stop
  end if    
  if (ulab(12) /= 'H2S(aq)') then
    write(*,*)
    write(*,*) ' Primary species 12 should be H2S(aq)'
    write(*,*)
  end if 

!!CIS Hardwired for Fe34S32S

!!CIS Hardwired for 34S32S (elemental sulfur)

  MoleFraction34S = sp10(13,jx,jy,jz)/( sp10(12,jx,jy,jz) + sp10(13,jx,jy,jz) )
  MoleFraction32S = 1.0d0 - MoleFraction34S

!!CIS Hardwired for 34S32S (elemental sulfur)

END IF    !!  End of JennyDruhan block  ********************************************************

decay_correct = 1.0D0
snorm = 0.0d0
ivolume = 0

DO k = 1,nkin
  dppt(k,jx,jy,jz) = 0.0D0

!!  **********  CrystalSizeDistribution  *********************************************
!!CSD  IF (CrystalSizeDistribution(k)) THEN  
!!CSD    DO l = 1,nCSD
!!CSD      AdjustKeq = 2.0d0*sigma(k)*volmol(k)/(rgas*Tk)* (1.0d0/radius(l))
!!CSD      AdjustKeqLog(l,k) = DLOG(AdjustKeq)
!!CSD    END DO
!!CSD  END IF
!!  **********  End of CrystalSizeDistribution  **************************************  

!!  DO np = 1,nreactmin(k)

   DO np = nreactmin(k),1,-1
    
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
!!  Imintype = 4 --> Precipitation only, thermodynamic term
!!  Imintype = 5 --> Dissolution only, thermodynamic term
!!  Imintype = 6 --> Forward one-way, surface area = 1.0 while mineral is present
!!  Imintype = 7 --> Reverse one-way, surface area = 1.0 while mineral is present
!!  Imintype = 8 --> Microbially mediated monod with theormodynamic factor, and biomass
!!  Imintype = 9 --> biomass first order decay !!smr-2011-04-29

!!  Calculate thermodynamic term
    IF (imintype(np,k) == 3 .OR. imintype(np,k) == 2) THEN  !! Simple Monod (Option 2 with no Ft) or Irreversible (Option 3)

      silog(np,k) = -500.0D0
      si(np,k) = 0.0D0

    ELSE

!!  Cross-Affinity Case
      IF (kcrossaff(np,k) /= 0) THEN

        kk = kcrossaff(np,k)
        sumiap = 0.0D0
        DO i = 1,ncomp
          sumiap = sumiap + decay_correct(i,k)*mumin(1,kk,i)*(sp(i,jx,jy,jz)+gam(i,jx,jy,jz))
        END DO
        silog(np,k) = (sumiap - keqmin(1,kk,jx,jy,jz))/clg
        si(np,k) = 10**(silog(np,k))

!!  Biomass case
      ELSE IF (imintype(np,k) == 8) THEN             !! Microbially mediated, with thermodynamic factor Ft

        jj = p_cat_min(k)
        sumiap = 0.0D0
        DO i = 1,ncomp
          sumiap = sumiap + muminTMP(np,jj,i)*(sp(i,jx,jy,jz)+gam(i,jx,jy,jz))
        END DO
        silog(np,k) = (sumiap - keqminTMP(np,jj) - BQ_min(np,jj)/(rgas*Tk))/clg    !!  BQ in kJ/mole
        siln(np,k) = clg*silog(np,k)
        si(np,k) = 10**(silog(np,k))

!!  Base case
      ELSE    
                                      
        sumiap = 0.0D0
        DO i = 1,ncomp
          sumiap = sumiap + decay_correct(i,k)*mumin(np,k,i)*(sp(i,jx,jy,jz)+gam(i,jx,jy,jz))
        END DO

!!  JennyDruhan case for isotopes, with modifications to IAP using solid solution activity (to be generalized)
        IF (JennyDruhan) THEN   !! ***********************************************************
!!CIS Hard coding 44Ca40CaCO3 
          IF (k == 5) THEN
            IF (MoleFraction40Mineral == 0.0d0) THEN
              CONTINUE
            ELSE
              sumiap = sumiap - DLOG(MoleFraction40Mineral)
            END IF
          END IF
          IF (k == 6) THEN
            IF (MoleFraction44Mineral == 0.0d0) THEN
              CONTINUE
            ELSE
              sumiap = sumiap - DLOG(MoleFraction44Mineral)
            END IF
          END IF
!!CIS End of hard coding 44Ca40CaCO3

!!CIS Hard coding Fe34S32S
          IF (k == 8) THEN
            sumiap = sumiap - DLOG(MoleFraction32)
          END IF
          IF (k == 9) THEN
            sumiap = sumiap - DLOG(MoleFraction34)
          END IF
!!CIS End of hard coding Fe34S32S

!!CIS Hard coding 34S32S (Elemental Sulfur)
          IF (k == 10) THEN
            sumiap = sumiap - DLOG(MoleFraction32S)
          END IF
          IF (k == 11) THEN
            sumiap = sumiap - DLOG(MoleFraction34S)
          END IF
!!CIS End of hard coding 34S32S (Elemental Sulfur)

        END IF    
!!  End of JennyDruhan block ***************************************************

        silog(np,k) = (sumiap - keqmin(np,k,jx,jy,jz))/clg
        
!!CSD!!  **********  CrystalSizeDistribution  **********************************
!!CSD        IF (CrystalSizeDistribution(k)) THEN
!!CSD
!!CSD!!        Cycle over crystal size distribution
!!CSD
!!CSD          DO l = 1,nCSD
!!CSD            SizeDepKeq = keqmin(np,k,jx,jy,jz) + AdjustKeqLog(l,k)
!!CSD            silogCrystalSize(l,np,k) =  (sumiap - SizeDepKeq)/clg
!!CSD!!          Assumes half sphere--> 0.5 * 4*pi*r^2
!!CSD            AreaCrystalSize(l,k) = NucleationScaleFactor*DFLOAT(nCrystal(l,k,jx,jy,jz)) * 2.0d0*pi*radius(l)*radius(l)
!!CSD            siCrystalSize(l,np,k) = 10**(silogCrystalSize(l,np,k)) 
!!CSD          END DO
!!CSD   
!!CSD        END IF
!!CSD!!  **********  CrystalSizeDistribution  *********************************************

        siln(np,k) = clg*silog(np,k)
        si(np,k) = 10**(silog(np,k))

      END IF

    END IF

!!  Decide on the value of the reactive surface area depending on
!!  whether the mineral is under or super-saturated and whether
!!  it has a zero volume fraction.

!!  Reset the surface areas for minerals that are associated with other minerals
!!  "MineralID" is the pointer to the mineral number whose volume fraction is being tracked
    
    IF (silog(np,k) >= 0.0D0 .AND. iarea(k,jinit(jx,jy,jz)) == 0) THEN       !! Supersaturated AND bulk_surface_area option

!!  Associate mineral with another mineral (surface area and volume fraction)
      IF (MineralAssociate(k)) THEN

        if (MineralID(k) < k) then                  !!  NOTE: This requires that the mineral that is associated with is earlier in list
          surf(k) = surf(MineralID(k))
        else
          surf(k) = areain(MineralID(k),jinit(jx,jy,jz))*porfactor
        end if

      ELSE

        if (JennyDruhan) then
          if (SetSurfaceAreaConstant) then
            surf(k) = areain(k,jinit(jx,jy,jz))*porfactor
          else
            surf(k) = area(k,jx,jy,jz)*porfactor
          endif
        else

          surf(k) = areain(k,jinit(jx,jy,jz))*porfactor

        end if

      END IF

!!  Case where either undersaturated OR using the specific_surface_area option  
    ELSE                                                

      IF (MineralAssociate(k)) THEN

        if (MineralID(k) < k) then
          surf(k) = surf(MineralID(k))
        else

          IF (porfactor < 0.01d0) THEN
            surf(k) = area(MineralID(k),jx,jy,jz)*porfactor
          ELSE
            surf(k) = area(MineralID(k),jx,jy,jz)
          END IF

        end if

      ELSE

        IF (porfactor < 0.01d0) THEN
          surf(k) = area(k,jx,jy,jz)*porfactor
        ELSE
          surf(k) = area(k,jx,jy,jz)
        END IF

        if (JennyDruhan) then
          if (SetSurfaceAreaConstant) then
            surf(k) = areain(k,jinit(jx,jy,jz))*porfactor
          else
            surf(k) = area(k,jx,jy,jz)*porfactor
          endif
        end if

      END IF

    END IF

!! Reset surface area for irreversible or Monod reaction
    IF (imintype(np,k) == 3 .OR. imintype(np,k) == 2) THEN

      IF (MineralAssociate(k)) THEN
        surf(k) = area(MineralID(k),jx,jy,jz)
      ELSE
        surf(k) = area(k,jx,jy,jz)
      END IF

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
          
!!    TST, irreversible, or PrecipitationOnly
      IF (imintype(np,k) == 1 .OR. imintype(np,k) == 3 .OR. imintype(np,k) == 4 .OR. imintype(np,k) == 5) THEN 
        
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

!!        IF (k == kUPlag .AND. OelkersRateLaw) THEN     !  Oelkerian rate law
!!          Kformation = 0.0000149823D0
!!          Al = sp10(ikAl,jx,jy,jz)
!!          term2 = Kformation/(Kformation + Al**0.333333)
!!        END IF
        
!!        ***** * Monod kinetics (without thermodynamic term)  **********
        
      ELSE IF (imintype(np,k) == 2) THEN     !  Monod rate law

!!  If in Monod rate law there is a Monod term for the mineral itself, do not use the mineral surface area (reset to 1.0) and use the hyperbolic expression  

        term2 = 1.0D0
        DO kk = 1,nmonod(np,k)
          IF (imonod(kk,np,k) == 0 .AND. kmonod(kk,np,k) == 1) THEN                !! Mineral Monod reaction depends on its own concentration
            surf(k) = 1.0d0                                                        !! Resets surface area to a constant 1
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
!!        IF (umin(k)=='JapaneseGlassSilicaOnly') THEN 
!!          kJapaneseGlassNoSilica = 2                                                     
!!          checkmonod =  0.01d0/(volfx(kJapaneseGlassNoSilica,jx,jy,jz) + 0.01d0)
!!          term2 = term2 * checkmonod 
!!        END IF
      ELSE IF (imintype(np,k) == 6 .OR. imintype(np,k) == 7) THEN  !! Forward or reverse, with surface area set constant at 1.0

        CONTINUE

!!    Monod rate law with biomass and thermodynamic factor, F_T     
      ELSE IF (imintype(np,k) == 8) THEN     

        jj = p_cat_min(k)

!!      Metabolic lag function (see Wood et al, 1995)
        IF (UseMetabolicLagMineral(np,jj)) THEN

!!        S* == Critical substrate concentration
!!        A* == Metabolic lag (in years)
!!        B* == Metabolic lag + ramp up period (so B* - A* = ramp up period (yrs)
          Astar = LagTimeMineral(np,jj)/365.0d0         !! Days converted to years
          Bstar = RampTimeMineral(np,jj)/365.0d0        !! Days converted to years
          Sstar = ThresholdConcentrationMineral(np,jj)  !! Critical concentration of substrate (acetate in this case)

          IF (sn(SubstrateForLagMineral(np,jj),jx,jy,jz) > Sstar .AND. tauZeroMineral(np,jj,jx,jy,jz) == 0.0d0) THEN
            tauZeroMineral(np,jj,jx,jy,jz) = time
          END IF

          IF (tauZeroMineral(np,jj,jx,jy,jz) == 0.0d0) THEN
            MetabolicLagMineral(np,jj,jx,jy,jz) = 0.0d0
          ELSE IF (tauZeroMineral(np,jj,jx,jy,jz) > 0.0d0 .AND. time < (tauZeroMineral(np,jj,jx,jy,jz)+Astar) ) THEN
            MetabolicLagMineral(np,jj,jx,jy,jz) = 0.0D0
          ELSE    
            denominator = (tauZeroMineral(np,jj,jx,jy,jz)+Astar+Bstar) - (tauZeroMineral(np,jj,jx,jy,jz)+Astar)
            IF (denominator == 0.0d0) THEN
              MetabolicLagMineral(np,jj,jx,jy,jz) = 1.0d0
            ELSE
              MetabolicLagMineral(np,jj,jx,jy,jz) =  (time -(tauZeroMineral(np,jj,jx,jy,jz)+Astar) )/denominator
            END IF
          END IF
          IF (MetabolicLagMineral(np,jj,jx,jy,jz) >= 1.0d0) THEN
            MetabolicLagMineral(np,jj,jx,jy,jz) = 1.0D0
          ELSE IF (MetabolicLagMineral(np,jj,jx,jy,jz) <= 0.0d0) THEN
            MetabolicLagMineral(np,jj,jx,jy,jz) = 0.0D0
          ELSE
            CONTINUE
          END IF

        END IF

!!  If in Monod rate law there is a Monod term for the mineral itself, do not use the mineral surface area (reset to 1.0) and use the hyperbolic expression  
             
        term2 = 1.0D0
        DO kk = 1,nmonod(np,k)
          IF (imonod(kk,np,k) == 0 .AND. kmonod(kk,np,k) == 1) THEN                !! Mineral Monod reaction depends on its own concentration
            surf(k) = 1.0d0                                                        !! Resets surface area to a constant 1
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

! frankfurt - sergi - add inhibition
        !DO kk = 1,ninhibit(np,k)
        !  i = inhibit(kk,np,k)
        !  IF (itot_inhibit(kk,np,k) == 1) THEN
        !    term2 = term2 * rinhibit(kk,np,k)/(rinhibit(kk,np,k)+s(i,jx,jy,jz))
        !  ELSE
        !    term2 = term2 * rinhibit(kk,np,k)/(rinhibit(kk,np,k)+sp10(i,jx,jy,jz))
        !  END IF
        ! END DO        
! biomass end

! sergi: first order biomass decay - may 2011
      
      ELSE IF (imintype(np,k) == 9) THEN

        term2 = volfx(biomass_decay(np,nkin),jx,jy,jz)
        !!write(*,*)'sergi: volume of biomass: ',volfx(biomass_decay(np,nkin),jx,jy,jz)
        pre_rmin(np,k) = term2
        continue
        
! end sergi: first order biomass decay - may 2011


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

! biomass
      ELSE IF (imintype(np,k) == 8) THEN   ! MonodBiomass kinetics

        pre_rmin(np,k) = term2
! biomass end        

! sergi: first order biomass decay - may 2011
      ELSE IF (imintype(np,k) == 9) THEN 

        pre_rmin(np,k) = term2
! end sergi: first order biomass decay - may 2011


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

! biomass
      ELSE IF (imintype(np,k) == 8) THEN    !!  Monod, but with thermodynamic factor (Ft) a la Jin and Bethke

! added
        jj = p_cat_min(k)

!!  Add in here the thermodynamic factor, F_T

!!  Formulation assumes that Delta G is negative and that it goes forward, but we have written the ppt reaction with 
!!  the mineral on the left.  So invert the Q/Keq if supersaturated (sign will capture the direction)

        IF( si(np,k) > 1.0d0) then

          if (direction_min(np,jj) == 1) THEN
            snorm(np,k) = (1.0d0/si(np,k))**(1.0d0/chi_min(np,jj))
          ELSE
            snorm(np,k) = 1.0d0
          END IF

        ELSE
 
          IF( direction_min(np,jj) == -1) THEN
            snorm(np,k) = si(np,k)**(1.0d0/chi_min(np,jj))
          ELSE
            snorm(np,k) = 1.0d0
          END IF

        END IF

 !!       snorm(np,k) = si(np,k)

        term1 = sign*DABS(snorm(np,k) - 1.0D0)

        IF (direction_min(np,jj) == 1) THEN               !! Precipitation only
          AffinityTerm = MAX(0.0d0,term1)
        ELSE                                         !! Dissolution only
          AffinityTerm = MIN(0.0d0,term1)
        ENDIF
! end Monodbiomass

! sergi: first order biomass decay - may 2011
      
      ELSE IF (imintype(np,k) == 9) THEN
      
        term1 = 1.d0
        AffinityTerm = 1.0d0
      
! end sergi: first order biomass decay - may 2011


      ELSE

        IF (AffinityDepend2(np,k) == 1.0D0 .AND. AffinityDepend3(np,k) == 1.0d0) THEN
          snorm(np,k) = si(np,k)
        ELSE IF (AffinityDepend2(np,k) /= 1.0d0 .AND. AffinityDepend3(np,k) == 1.0d0) THEN
          snorm(np,k) = si(np,k)**AffinityDepend2(np,k)
        ELSE
          snorm(np,k) = DEXP(-AffinityDepend2(np,k)*(DABS(siln(np,k)))**AffinityDepend3(np,k) )
        END IF

      END IF    !! Sergi:  Need to move this down?

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

!!    Sergi:  Need to move the END IF above here??


!************** End of dependence on mineral saturation state *************

! commented out        rmin(np,k) = surf(k)*rate0(np,k)*actenergy(np,k)*pre_rmin(np,k)*AffinityTerm 

!!    Biomass, with thermodynamic term, Ft         
      IF (imintype(np,k) == 8) then   !! MonodBiomass

!       pointer to biomass for current reaction
        jj = p_cat_min(k)
!       pointer to biomass for current reaction
        ib = ibiomass_min(jj)

        IF (UseMetabolicLagMineral(np,jj)) THEN
          rmin(np,k) = MetabolicLagMineral(np,jj,jx,jy,jz)*surf(k)*rate0(np,k)*           &
              actenergy(np,k)*pre_rmin(np,k)*volfx(ib,jx,jy,jz)*AffinityTerm
        ELSE
          rmin(np,k) = surf(k)*rate0(np,k)*           &
             actenergy(np,k)*pre_rmin(np,k)*volfx(ib,jx,jy,jz)*AffinityTerm
        END IF
        

      ELSE

        IF (JennyDruhan) THEN       !!  *******************************************************************
!!CIS Hardwired for 44Ca40CaCO3
          IF (k == 5) THEN

            if (umin(k) /= 'CalciteRifle') then
              write(*,*)
              write(*,*) ' Mineral 5 should be CalciteRifle'
              write(*,*)
              read(*,*)
              stop
            end if 

!!  DePaolo = .TRUE. case 

            IF (DePaolo) THEN

!!            Irreversible dissolution
              IF (np == 1) THEN

                rmin(np,k) = MoleFraction40Mineral*surf(k)*rate0(np,k)*actenergy(np,k)*pre_rmin(np,k)*AffinityTerm

!!                CurrentMoleFraction40 = rmin(2,5)/(rmin(2,6)+rmin(2,5))

!!                denominator = rmin(2,5) + rmin(2,6) + rmin(1,5) + rmin(1,6)
!!                if (time==0.0 .or. denominator==0.0) then
!!                  CurrentMoleFraction40 = MoleFraction40
!!                else
!!                  CurrentMoleFraction40 = (rmin(2,5)+rmin(1,5))/denominator
!!                endif
!!                if (CurrentMoleFraction40 < 0.0) then
!!                   write(*,*) ' CurrentMoleFraction40: ', CurrentMoleFraction40
!!                end if
!!                rmin(np,k) = CurrentMoleFraction40*surf(k)*rate0(np,k)*actenergy(np,k)*pre_rmin(np,k)*AffinityTerm
!!                rmin(np,k) = CurrentMoleFraction40*surf(k)*rate0(np,k)*pre_rmin(np,k)

                rminSaveForDepaolo(1,jx,jy,jz) = rmin(np,k)
              ELSE IF (np == 2) THEN
!!            Irreversible precipitation
                rmin(np,k) = surf(k)*rate0(np,k)*actenergy(np,k)*pre_rmin(np,k)*(-1.0d0*AffinityTerm)
                rminSaveForDepaolo(2,jx,jy,jz) = rmin(np,k)
              ELSE
                write(*,*)' No other option'
                read(*,*)
                stop
              END IF
            ELSE

!!  DePaolo = .FALSE.

              CalciumCarbonateRatioEffect = 1.0
              
              rmin(np,k) = CalciumCarbonateRatioEffect*MoleFraction40Mineral*(surf(5))*rate0(np,k)*actenergy(np,k)*pre_rmin(np,k)*AffinityTerm
              
              kBackward =  CalciumCarbonateRatioEffect*rate0(np,k)
              kForward = kBackward/10**( -8.4801)
              IF (JennyFirstOrder) THEN
                rminSaveForDepaolo(1,jx,jy,jz) = -kBackward*surf(5)*MoleFraction40Mineral*( 3.0*si(np,k)-(si(np,k)*si(np,k)) - 1.0)  !! Dissolution
                rminSaveForDepaolo(2,jx,jy,jz) = kForward* surf(5)*DEXP(gam(6,jx,jy,jz)+sp(6,jx,jy,jz))*DEXP(gam(27,jx,jy,jz)+sp(27,jx,jy,jz))   !! Precipitation
                rminCheck = rminSaveForDepaolo(1,jx,jy,jz) + rminSaveForDepaolo(2,jx,jy,jz)
              ELSE
                rminSaveForDepaolo(2,jx,jy,jz) = kForward* surf(5)*DEXP(2.0*(gam(6,jx,jy,jz)+sp(6,jx,jy,jz)))*DEXP(2.0*(gam(27,jx,jy,jz)+sp(27,jx,jy,jz)))   !! Precipitation
                rminSaveForDepaolo(1,jx,jy,jz) = rmin(np,k) - rminSaveForDepaolo(2,jx,jy,jz)  !! Dissolution

              END IF
            END IF


          ELSE IF (k == 6) THEN
 
            if (umin(k) /= 'Calcite44Rifle') then
              write(*,*)
              write(*,*) ' Mineral 6 should be Calcite44Rifle'
              write(*,*)
              read(*,*)
              stop
            end if
            IF (DePaolo) THEN
!!            Irreversible dissolution
              IF (np == 1) THEN

                rmin(np,k) = MoleFraction44Mineral*surf(5)*rate0(np,k)*actenergy(np,k)*pre_rmin(np,k)*AffinityTerm

!!                CurrentMoleFraction44 = rmin(2,6)/(rmin(2,6)+rmin(2,5))

!!                denominator = rmin(2,5) + rmin(2,6) + rmin(1,5) + rmin(1,6)
!!              if (time==0.0 .or. denominator==0.0) then
!!                  CurrentMoleFraction44 = MoleFraction44
!!              else
!!                CurrentMoleFraction44 = (rmin(2,6)+rmin(1,6))/denominator
!!              endif
!!                if (CurrentMoleFraction44 < 0.0) then
!!                   write(*,*) ' MoleFraction44: ', CurrentMoleFraction44,MoleFraction44
!!                end if

!!                rmin(np,k) = CurrentMoleFraction44*surf(5)*rate0(np,k)*actenergy(np,k)*pre_rmin(np,k)*AffinityTerm

              ELSE IF (np == 2) THEN
!!            Irreversible precipitation
                rmin(np,k) = surf(5)*rate0(np,k)*actenergy(np,k)*pre_rmin(np,k)*(-1.0d0*AffinityTerm)
              ELSE
                write(*,*)' No other option'
                read(*,*)
                stop
              END IF
            ELSE                !! DePaolo = .FALSE.
                
              CalciumCarbonateRatioEffect = 6.1878*DEXP( -0.228*( sp(6,jx,jy,jz) - sp(27,jx,jy,jz) ) )
              CalciumCarbonateRatioEffect = 1.0
              rmin(np,k) = CalciumCarbonateRatioEffect*MoleFraction44Mineral*(surf(5))*rate0(np,k)*actenergy(np,k)*pre_rmin(np,k)*AffinityTerm
            END IF
!!CIS Hardwired for 44Ca40CaCO3
!!CIS Hardwired for Fe34S32S
          ELSE IF (k == 8) then
            if (umin(k) /= 'FeS(am)') then
              write(*,*)
              write(*,*) ' Mineral 8 should be FeS(am)'
              write(*,*)
              read(*,*)
              stop
            end if
            rmin(np,k) = MoleFraction32*surf(k)*rate0(np,k)*actenergy(np,k)*pre_rmin(np,k)*AffinityTerm
          ELSE IF (k == 9) then
            if (umin(k) /= 'FeS34(am)') then
              write(*,*)
              write(*,*) ' Mineral 9 should be FeS34(am)'
              write(*,*)
              read(*,*)
              stop
            end if
            rmin(np,k) = MoleFraction34*surf(8)*rate0(np,k)*actenergy(np,k)*pre_rmin(np,k)*AffinityTerm
!!CIS Hardwired for Fe34S32S
!!CIS Hardwired for 34S32S (Elemental Sulfur)
          ELSE IF (k == 10) then
            if (umin(k) /= 'S32') then
              write(*,*)
              write(*,*) ' Mineral 10 should be S32'
              write(*,*)
              read(*,*)
              stop
            end if
            rmin(np,k) = MoleFraction32S*surf(k)*rate0(np,k)*actenergy(np,k)*pre_rmin(np,k)*AffinityTerm
          ELSE IF (k == 11) then
            if (umin(k) /= 'S34') then
              write(*,*)
              write(*,*) ' Mineral 11 should be S34'
              write(*,*)
              read(*,*)
              stop
            end if
            rmin(np,k) = MoleFraction34S*surf(10)*rate0(np,k)*actenergy(np,k)*pre_rmin(np,k)*AffinityTerm

!!CIS Hardwired for 34S32S (Elemental Sulfur)
          ELSE 
            rmin(np,k) = surf(k)*rate0(np,k)*actenergy(np,k)*pre_rmin(np,k)*AffinityTerm
          END IF 
!!      End of JennyDruhan conditional

        ELSE

          rmin(np,k) = surf(k)*rate0(np,k)*actenergy(np,k)*pre_rmin(np,k)*AffinityTerm


        END IF   

      END IF
! biomass end

    END IF


!!CSD        ELSE


    dppt(k,jx,jy,jz) = dppt(k,jx,jy,jz) + rmin(np,k)

!!CSD        END IF
!!      END IF
      
!!    END IF
    
  END DO   !  End of npth parallel reaction
  
  IF (MineralAssociate(k)) THEN
    vcheck = volfx(MineralId(k),jx,jy,jz) + dppt(k,jx,jy,jz)*volmol(MineralId(k))*delt
  ELSE
    vcheck = volfx(k,jx,jy,jz) + dppt(k,jx,jy,jz)*volmol(k)*delt
  END IF
  
  IF (vcheck < 0.0) THEN

    IF (MineralAssociate(k)) THEN
      dppt(k,jx,jy,jz) = -volfx(MineralID(k),jx,jy,jz)/(volmol(MineralId(k))*delt)
      rmin(1,k) = dppt(k,jx,jy,jz)
      IF (nreactmin(k) > 1) THEN
        DO np = 2,nreactmin(k)
          rmin(np,k) = 0.0
        END DO
      END IF
      ivolume(k) = 1

    ELSE

      dppt(k,jx,jy,jz) = -volfx(k,jx,jy,jz)/(volmol(k)*delt)
      rmin(1,k) = dppt(k,jx,jy,jz)
      IF (nreactmin(k) > 1) THEN
        DO np = 2,nreactmin(k)
          rmin(np,k) = 0.0
        END DO
      END IF
      ivolume(k) = 1

    END IF
    
  END IF
  
END DO     !  End of kth mineral

RETURN
END SUBROUTINE reaction
