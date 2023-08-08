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


SUBROUTINE reactkin(ncomp,nspec,nrct,ikin,jx,jy,jz,AqueousToBulk,time)
USE crunchtype
USE params
USE concentration
use mineral, only:     volfx,volmol,  &
                       mukinTMP,  & 
                       keqkinTMP, &
                       p_cat_kin, &
                       ibiomass_kin, &
                       bq_kin, chi_kin, direction_kin, &
                       UseMetabolicLagAqueous,LagTimeAqueous, &
                       MetabolicLagAqueous,                   &
                       RampTimeAqueous,ThresholdConcentrationAqueous,  &
                       SubstrateForLagMineral,SubstrateForLagAqueous,  &
                       nMonodBiomassMineral,nMonodBiomassAqueous,tauZeroAqueous, &
                       satlog
USE medium
USE temperature, ONLY: ro,T
USE strings
USE runtime, ONLY: JennyDruhan,Maggi
USE transport, ONLY: satliq,satliqold,dorm_a,dorm_b,dorm_Se, thres_OM1, thres_OM2, exp_OM
USE isotope

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                       :: ncomp
INTEGER(I4B), INTENT(IN)                                       :: nspec
INTEGER(I4B), INTENT(IN)                                       :: nrct
INTEGER(I4B), INTENT(IN)                                       :: ikin
INTEGER(I4B), INTENT(IN)                                       :: jx
INTEGER(I4B), INTENT(IN)                                       :: jy
INTEGER(I4B), INTENT(IN)                                       :: jz
REAL(DP), INTENT(IN)                                           :: AqueousToBulk
REAL(DP), INTENT(IN)                                           :: time

!  Internal variables and arrays

INTEGER(I4B)                                                   :: ir
INTEGER(I4B)                                                   :: i
INTEGER(I4B)                                                   :: ll
INTEGER(I4B)                                                   :: id
INTEGER(I4B)                                                   :: k
! biomass
INTEGER(I4B)                                                   :: jj, ib
! biomass end

INTEGER(I4B)                                                   :: IsotopologueOther
INTEGER(I4B)                                                   :: nnisotope

REAL(DP)                                                       :: affinity
REAL(DP)                                                       :: sum
!!REAL(DP)                                                       :: satlog
REAL(DP)                                                       :: term2
REAL(DP)                                                       :: sumkin
REAL(DP)                                                       :: term_inhibit
REAL(DP)                                                       :: MinConvert

REAL(DP)                                                       :: Astar
REAL(DP)                                                       :: Bstar
REAL(DP)                                                       :: Sstar
REAL(DP)                                                       :: denominator

! biomass
REAL(DP)                                                       :: bqTMP
REAL(DP)                                                       :: tk
REAL(DP)                                                       :: tkinv
REAL(DP)                                                       :: reft
REAL(DP)                                                       :: sign
REAL(DP)                                                       :: term1
REAL(DP)                                                       :: termTMP
REAL(DP)                                                       :: snormAqueous

real(dp)                                                       :: vol_temp
real(dp)                                                       :: satL
REAL(DP)                                                       :: liqsat_fac
!!REAL(DP), DIMENSION(ikin)                                      :: MoleFraction

tk = t(jx,jy,jz) + 273.15D0
satL = 0.5*(satliq(jx,jy,jz) + satliqold(jx,jy,jz) )
tkinv = 1.0D0/tk
reft = 1.0D0/(298.15D0)

!!MoleFractionCommon = 1.0d0
!!MoleFractionRare = 1.0d0

!!IF (nIsotopePrimary > 0) THEN
!!  DO nnisotope = 1,nIsotopePrimary
!!    denominator = ( sp10(isotopeRare(nnisotope),jx,jy,jz)+sp10(isotopeCommon(nnisotope),jx,jy,jz)  )
!!    MoleFractionCommon(nnisotope) = sp10(isotopeCommon(nnisotope),jx,jy,jz)/denominator
!!    MoleFractionRare(nnisotope) = sp10(isotopeRare(nnisotope),jx,jy,jz)/denominator
!!  END DO
!!END IF

!!CIS IF (.NOT. Maggi .AND. JennyDruhan) THEN 
!!CIS   MoleFraction(2) = sp10(11,jx,jy,jz)/( sp10(10,jx,jy,jz) + sp10(11,jx,jy,jz) )
!!CIS   MoleFraction(1) = 1.0d0 - MoleFraction(2)
!!CIS END IF

!! IAQTYPE options
!! iaqtype = 1 --> TST
!! iaqtype = 2 --> Simple Monod
!! iaqtype = 3 --> irreversible
!! iaqtype = 4 --> radioactive decay
!! iaqtype = 8 --> MonodBiomass

DO ir = 1,ikin
  
  IF (iaqtype(ir) == 3 .OR. iaqtype(ir) == 2 .OR. iaqtype(ir) == 4) THEN
    satkin(ir) = 1.0d0
    affinity = 1.0d0
  ELSE
    sum = 0.0d0
    DO i = 1,ncomp
      sum = sum + mukin(ir,i)*sp(i,jx,jy,jz)
    END DO
    satlog(ir,jx,jy,jz) = sum - clg*keqkin(ir)
    satkin(ir) = DEXP(satlog(ir,jx,jy,jz))
    affinity = 1.0 - satkin(ir)
  END IF
  
!! TST or irreversible or radioactive decay (non-hyperbolic expressions)
  IF (iaqtype(ir) == 1 .OR. iaqtype(ir) == 3 .OR. iaqtype(ir) == 4) THEN  

    DO ll = 1,nreactkin(ir)
      sum = 0.0d0
      DO i = 1,ncomp
        IF (itot(i,ll,ir) == 1) THEN
          IF (ierode == 1) THEN
!!            sum = sum + dependk(i,ll,ir)*LOG(s(i,jx,jy,jz) )                   !  Include only aqueous
            sum = sum + dependk(i,ll,ir)*LOG(s(i,jx,jy,jz)+sch(i,jx,jy,jz)/AqueousToBulk)    !  Include sorbed mass
          ELSE
!!            sum = sum + dependk(i,ll,ir)*LOG(s(i,jx,jy,jz))
            sum = sum + dependk(i,ll,ir)*LOG(s(i,jx,jy,jz) + (sNCexch_local(i)+sNCsurf_local(i))/AqueousToBulk)
          END IF
        ELSE
          sum = sum + dependk(i,ll,ir)*sp(i,jx,jy,jz)
        END IF
      END DO
      IF (sum == 0.0d0) THEN
        pre_raq(ll,ir) = 1.0
      ELSE
        pre_raq(ll,ir) = DEXP(sum)
      END IF

      IF (namkin(ir) == 'Biomass_dormancy') THEN
      satL = 0.5*(satliq(jx,jy,jz) + satliqold(jx,jy,jz) )
      pre_raq(ll,ir) = pre_raq(ll,ir) * (1/(1 + (satL/dorm_Se)**dorm_b))
      ELSEIF (namkin(ir) == 'Biomass_activation') THEN
      satL = 0.5*(satliq(jx,jy,jz) + satliqold(jx,jy,jz) )
      pre_raq(ll,ir) = pre_raq(ll,ir )* (1/(1 + dorm_a*(dorm_Se/satL)**dorm_b))
      ELSEIF (namkin(ir) == 'TOCsoil_release' .OR. namkin(ir) == 'TOCsoil_release') THEN
      satL = 0.5*(satliq(jx,jy,jz) + satliqold(jx,jy,jz) )
      liqsat_fac = 1
      IF (satL > thres_OM2) THEN
        liqsat_fac = 1/(1 + (thres_OM1/satL)**exp_OM)
      ELSE
        liqsat_fac = 1/(1 + (thres_OM1/thres_OM2)**exp_OM)
      ENDIF
      pre_raq(ll,ir) = pre_raq(ll,ir ) * liqsat_fac
      ENDIF

    END DO

  ELSE IF (iaqtype(ir) == 2) THEN    ! Simple Monod kinetics

!! NOTE: ---> This assumes only ONE parallel reaction for Monod
   
!!  Normal Monod terms

    term2 = 1.0d0

    DO id = 1,nmonodaq(ir)
      i = imonodaq(id,ir)
      IF (itot_monodaq(id,ir) == 1) THEN
        term2 = term2 * s(i,jx,jy,jz)/(s(i,jx,jy,jz) + halfsataq(id,ir))
      ELSE
        termTMP = sp10(i,jx,jy,jz)/( sp10(i,jx,jy,jz) + halfsataq(id,ir) )
        term2 = term2 * sp10(i,jx,jy,jz)/(sp10(i,jx,jy,jz) + halfsataq(id,ir))
      END IF
    END DO


!!  Inhibition terms

    DO id = 1,ninhibitaq(ir)
      i = inhibitaq(id,ir)
      IF (inhibitaq(id,ir) < 0) THEN                 !! Dependence on mineral volume fraction
        k = -inhibitaq(id,ir)
        MinConvert = volfx(k,jx,jy,jz)/(volmol(k)*por(jx,jy,jz)*ro(jx,jy,jz)*satliq(jx,jy,jz))  !! Converts mineral volume fraction to moles mineral per kg fluid (molality)                                  
        term_inhibit =  rinhibitaq(id,ir)/(MinConvert + rinhibitaq(id,ir))
        term2 = term2 * term_inhibit
      ELSE
        IF (itot_inhibitaq(id,ir) == 1) THEN
          term_inhibit = rinhibitaq(id,ir)/(rinhibitaq(id,ir)+s(i,jx,jy,jz))
          term2 = term2 * term_inhibit
        ELSE
          term_inhibit = rinhibitaq(id,ir)/(rinhibitaq(id,ir)+ sp10(i,jx,jy,jz))
          term2 = term2 * term_inhibit
        END IF
      END IF
    END DO
    
!! NOTE: ---> This assumes only ONE parallel reaction for Monod

    IF (direction_kin(ir) < 0) THEN
      sign = 1.0d0
    ELSE
      sign = -1.0d0
    END IF

    if (sign*term2 < 0.0) then
      continue
    endif
    pre_raq(1,ir) = sign*term2

!!    affinity = 1.0

!! Biomass option
  ELSE IF (iaqtype(ir) == 8 .OR. iaqtype(ir) == 9) THEN    ! Monod kinetics, but with thermodynamic fact, F_T

!! NOTE: ---> This assumes only ONE parallel reaction for Monod

!!  Normal Monod terms

    term2 = 1.0d0

    DO id = 1,nmonodaq(ir)
      i = imonodaq(id,ir)

!!!      IF (nIsotopePrimary > 0) THEN

        IF (IsotopePrimaryCommon(i)) THEN
         
          IsotopologueOther = isotopeRare(iPointerIsotope(i))
          IF (itot_monodaq(id,ir) == 1) THEN

            termMonod(id,ir) = s(i,jx,jy,jz)/( s(i,jx,jy,jz) + halfsataq(id,ir)*(1.0d0+s(IsotopologueOther,jx,jy,jz)/halfsataq(id,ir) ) )
            term2 = term2 * termMonod(id,ir)


!! For case where isotopologue half sats are NOT the same (hard wired, so needs additional pointers to find out which 
!!     aqueous reaction to use
!!CIS            termMonod(id,ir) = s(i,jx,jy,jz)/( s(i,jx,jy,jz) + halfsataq(id,ir)*(1.0d0+s(IsotopologueOther,jx,jy,jz)/halfsataq(2,2)) )

          ELSE      !! Case where individual species concentration is used in hyperbolic term

            termMonod(id,ir) = sp10(i,jx,jy,jz)/( sp10(i,jx,jy,jz) + halfsataq(id,ir)*(1.0d0+sp10(IsotopologueOther,jx,jy,jz)/halfsataq(id,ir) ) )
            term2 = term2 * termMonod(id,ir)

          END IF

        ELSE IF (IsotopePrimaryRare(i)) THEN
         
          IsotopologueOther = isotopeCommon(iPointerIsotope(i))
          IF (itot_monodaq(id,ir) == 1) THEN

            termMonod(id,ir) = s(i,jx,jy,jz)/( s(i,jx,jy,jz) + halfsataq(id,ir)*(1.0d0+s(IsotopologueOther,jx,jy,jz)/halfsataq(id,ir) ) )
            term2 = term2 * termMonod(id,ir)

          ELSE      !! Case where individual species concentration is used in hyperbolic term

            termMonod(id,ir) = sp10(i,jx,jy,jz)/( sp10(i,jx,jy,jz) + halfsataq(id,ir)*(1.0d0+sp10(IsotopologueOther,jx,jy,jz)/halfsataq(id,ir) ) )
            term2 = term2 * termMonod(id,ir)

          END IF
          
!!!        ELSE
          
!!!          WRITE(*,*)
!!!          WRITE(*,*) ' Should be either the Common or Rare isotope in reactkin.F90'
!!!          WRITE(*,*)
!!!          READ(*,*)
!!!          STOP
        
!!!        END IF

      ELSE   !! General case where no isotopes involved

        IF (itot_monodaq(id,ir) == 1) THEN

         termMonod(id,ir) = s(i,jx,jy,jz)/(s(i,jx,jy,jz)+halfsataq(id,ir))
         term2 = term2 * termMonod(id,ir)

        ELSE      !! Case where individual species concentration is used in hyperbolic term

         termMonod(id,ir) = sp10(i,jx,jy,jz)/(sp10(i,jx,jy,jz)+halfsataq(id,ir))
         term2 = term2 * termMonod(id,ir)
        END IF

      END IF

    END DO 

! prague - add inhibition - 

!!!  Inhibition terms
!
    DO id = 1,ninhibitaq(ir)
      i = inhibitaq(id,ir)
      IF (inhibitaq(id,ir) < 0) THEN                 !! Dependence on mineral volume fraction
        k = -inhibitaq(id,ir)
        MinConvert = volfx(k,jx,jy,jz)/(volmol(k)*por(jx,jy,jz)*ro(jx,jy,jz)*satliq(jx,jy,jz))  !! Converts mineral volume fraction to moles mineral per kg fluid (molality)                                  
        term_inhibit =  rinhibitaq(id,ir)/(MinConvert + rinhibitaq(id,ir))
        term2 = term2 * term_inhibit
        write(*,*) k,term_inhibit
      ELSE
        IF (itot_inhibitaq(id,ir) == 1) THEN
          term_inhibit = rinhibitaq(id,ir)/(rinhibitaq(id,ir)+s(i,jx,jy,jz))
          term2 = term2 * term_inhibit
        ELSE
          term_inhibit = rinhibitaq(id,ir)/(rinhibitaq(id,ir)+ sp10(i,jx,jy,jz))
          term2 = term2 * term_inhibit
        END IF
      END IF
    END DO
    
! end add inhibition

!!  No inhibition terms, but add thermodynamic factor, F_T


    !! jj = p_cat_kin(ir)
    jj = ir
!!    write(*,*) ' Pointer here?'
!!    write(*,*) p_cat_kin(ir),ir
!!    read(*,*)

!! NOTE:  Sergi seems to make jj equivalent to ir here, but the Lag arrays may require a pointer.

!!  Metabolic lag function (see Wood et al, 1995)
    IF (UseMetabolicLagAqueous(jj)) THEN

!!    S* == Critical substrate concentration
!!    A* == Metabolic lag (in years)
!!    B* == Metabolic lag + ramp up period (so B* - A* = ramp up period (yrs)
      Astar = LagTimeAqueous(jj)/365.0d0         !! Days converted to years
      Bstar = RampTimeAqueous(jj)/365.0d0        !! Days converted to years
      Sstar = ThresholdConcentrationAqueous(jj)  !! Critical concentration of substrate (acetate in this case)

      IF (sn(SubstrateForLagAqueous(jj),jx,jy,jz) > Sstar .AND. tauZeroAqueous(jj,jx,jy,jz) == 0.0d0) THEN
            tauZeroAqueous(jj,jx,jy,jz) = time
      END IF

      IF (tauZeroAqueous(jj,jx,jy,jz) == 0.0d0) THEN
        MetabolicLagAqueous(jj,jx,jy,jz) = 0.0d0
      ELSE IF (tauZeroAqueous(jj,jx,jy,jz) > 0.0d0 .AND. time < (tauZeroAqueous(jj,jx,jy,jz)+Astar) ) THEN
        MetabolicLagAqueous(jj,jx,jy,jz) = 0.0D0
      ELSE    
        denominator = (tauZeroAqueous(jj,jx,jy,jz)+Astar+Bstar) - (tauZeroAqueous(jj,jx,jy,jz)+Astar)
        IF (denominator == 0.0d0) THEN
           MetabolicLagAqueous(jj,jx,jy,jz) = 1.0d0
        ELSE
          MetabolicLagAqueous(jj,jx,jy,jz) =  (time -(tauZeroAqueous(jj,jx,jy,jz)+Astar) )/denominator
        END IF
      END IF
      IF (MetabolicLagAqueous(jj,jx,jy,jz) >= 1.0d0) THEN
        MetabolicLagAqueous(jj,jx,jy,jz) = 1.0D0
      ELSE IF (MetabolicLagAqueous(jj,jx,jy,jz) <= 0.0d0) THEN
        MetabolicLagAqueous(jj,jx,jy,jz) = 0.0D0
      ELSE
        CONTINUE
      END IF

    END IF
      
    sum = 0.0d0
    DO i = 1,ncomp
      sum = sum + mukinTMP(jj,i)*sp(i,jx,jy,jz)
    END DO
    i = chi_kin(jj)
    i = direction_kin(jj)
    bqTMP = bq_kin(jj) ! bq is negative according to carl's convention
    satlog(jj,jx,jy,jz) = sum - clg*keqkinTMP(jj) - bqTMP/(rgas*Tk)
    satkin(ir) = DEXP(satlog(jj,jx,jy,jz))

!! Assume chiAqueous = 1 for now, which should be the case for an anaerobic reaction written in terms of one electron

!!    IF( satkin(ir) > 1.0d0) THEN
!!      snormAqueous(ir) = (1.0d0/satkin(ir))**(1.0d0/chiAqueous(ir))
!!    ELSE 
!!      snormAqueous(ir) = (satkin(ir))**(1.0d0/chiAqueous(ir))
!!    ENDIF

    IF (direction_kin(jj) < 0) THEN
      sign = 1.0d0
    ELSE
      sign = -1.0d0
    END IF

    ! IF (satkin(ir) > 1.0d0) THEN
    !   snormAqueous = 1.0d0
    ! ELSE 
    !   snormAqueous = satkin(ir)
    ! END IF
    
    ! !!  Reaction assumed to be irreversible, so do not let it go in reverse

    ! term1 = sign*DABS(snormAqueous - 1.0D0)
    ! affinity = MAX(0.0d0,term1)

    !Lucien (regular affinity calculation):
    snormAqueous = satkin(ir)
    affinity = 1 - satkin(ir)

    ! if (satkin(ir) >= 1) then
    ! affinity = 0
    ! sign = 0
    ! endif
    
!! NOTE: ---> This assumes only ONE parallel reaction for Monod

    pre_raq(1,ir) = sign*term2
! biomass end

    
  ELSE

    WRITE(*,*)
    WRITE(*,*) ' Rate formulation not recognized'
    WRITE(*,*) ' Iaqtype = ', iaqtype(ir)
    WRITE(*,*)
    READ(*,*)
    STOP

  END IF

! biomass
  if (iaqtype(ir) == 8) then

!   pointer to biomass for current reaction
  !!  ib = ibiomass_kin(p_cat_kin(ir))
    ib = ibiomass_kin(ir)

    ! note on units:
    ! 
    ! Convert units of rate constant for micr.-mediated react. to mol/Kg-H2O/yr:
    !   raq_tot will later be multiplied by por * satL * ro (in assemble_local and assmeblePS02)
    !   There the mass balance is expressed per bulk volume (mol/m3-bulk/yr).
    !
    ! Rate constants for non-biomass reactions are defined as mol/Kg-H2O/yr. 
    ! For  microbially-mediated reactions with a biomass term, rate constants in input
    !    are in units of mol-reaction/mol-biomass/yr
    !  The volume fraction of the biomass is divided by the molar volume, the liquid saturation, the porosity,
    !    and the fluid density
    

    vol_temp = volfx(ib,jx,jy,jz) / (por(jx,jy,jz) * ro(jx,jy,jz) * satL) ! [mol biomass / Kg-H2O]
    
!!!    m3_min/m3_pm* mol/ m3_min * m3_pm/m3_fluid * m3_fluid/kgw = mol/kgw * mol/kgw/yr   -->  mol / kgw / yr

    IF (UseMetabolicLagAqueous(jj)) THEN
      sumkin = 0.0
      DO ll = 1,nreactkin(ir)
        raq(ll,ir) = MetabolicLagAqueous(jj,jx,jy,jz)*vol_temp*ratek(ll,ir)*pre_raq(ll,ir)*affinity
        sumkin = sumkin + raq(ll,ir)
      END DO
    ELSE
      sumkin = 0.0
      DO ll = 1,nreactkin(ir)
        ! ************************************
        ! Edit by Lucien Stolze, June 2023
        actenergyaq(ll,ir) = DEXP( (actk(ll,ir)/rgasKCAL)*(reft - tkinv) )
        raq(ll,ir) = ratek(ll,ir)*vol_temp*pre_raq(ll,ir)*affinity*actenergyaq(ll,ir) ![mol/mol-biomass/yr]*[mol-biomass/kgw/yr]
        ! ************************************
        sumkin = sumkin + raq(ll,ir)
      END DO
    END IF

  elseif (iaqtype(ir) == 9) then

!   pointer to biomass for current reaction
  !!  ib = ibiomass_kin(p_cat_kin(ir))
    ib = ibiomass_kin(ir)

    ! note on units:
    ! 
    ! Convert units of rate constant for micr.-mediated react. to mol/Kg-H2O/yr:
    !   raq_tot will later be multiplied by por * satL * ro (in assemble_local and assmeblePS02)
    !   There the mass balance is expressed per bulk volume (mol/m3-bulk/yr).
    !
    ! Rate constants for non-biomass reactions are defined as mol/Kg-H2O/yr. 
    ! For  microbially-mediated reactions with a biomass term, rate constants in input
    !    are in units of mol-reaction/mol-biomass/yr
    !  The volume fraction of the biomass is divided by the molar volume, the liquid saturation, the porosity,
    !    and the fluid density
    

    vol_temp = sp10(ib,jx,jy,jz) ! [mol biomass / Kg-H2O]
    
!!!    m3_min/m3_pm* mol/ m3_min * m3_pm/m3_fluid * m3_fluid/kgw = mol/kgw * mol/kgw/yr   -->  mol / kgw / yr

    IF (UseMetabolicLagAqueous(jj)) THEN
      sumkin = 0.0
      DO ll = 1,nreactkin(ir)
        raq(ll,ir) = MetabolicLagAqueous(jj,jx,jy,jz)*vol_temp*ratek(ll,ir)*pre_raq(ll,ir)*affinity
        sumkin = sumkin + raq(ll,ir)
      END DO
    ELSE
      sumkin = 0.0
      DO ll = 1,nreactkin(ir)
        ! ************************************
        ! Edit by Lucien Stolze, June 2023
        actenergyaq(ll,ir) = DEXP( (actk(ll,ir)/rgasKCAL)*(reft - tkinv) )
        raq(ll,ir) = ratek(ll,ir)*vol_temp*pre_raq(ll,ir)*affinity*actenergyaq(ll,ir) ![mol/mol-biomass/yr]*[mol-biomass/kgw]
        ! ************************************
        sumkin = sumkin + raq(ll,ir)
      END DO
    END IF

  else

    sumkin = 0.0
    DO ll = 1,nreactkin(ir)
      ! ************************************
      ! Edit by Lucien Stolze, June 2023
      actenergyaq(ll,ir) = DEXP( (actk(ll,ir)/rgasKCAL)*(reft - tkinv) )
      raq(ll,ir) = ratek(ll,ir)*pre_raq(ll,ir)*affinity*actenergyaq(ll,ir)
      ! ************************************
      sumkin = sumkin + raq(ll,ir)
    END DO

  end if
! biomass end
    
  raq_tot(ir,jx,jy,jz) = sumkin
  
END DO

RETURN
END SUBROUTINE reactkin

