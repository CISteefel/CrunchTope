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
    
SUBROUTINE gamma_init(ncomp,nspec,tempc,sqrt_sion,sion_tmp)
USE crunchtype
USE params
USE runtime, ONLY: Benchmark
USE concentration
USE temperature
USE mineral

IMPLICIT NONE

!  External variables

INTEGER(I4B), INTENT(IN)                                   :: ncomp
INTEGER(I4B), INTENT(IN)                                   :: nspec

REAL(DP), INTENT(IN)                                       :: tempc
REAL(DP), INTENT(OUT)                                      :: sqrt_sion

!  Internal variables

REAL(DP)                                                   :: TotalMoles
REAL(DP)                                                   :: ah
REAL(DP)                                                   :: bh
REAL(DP)                                                   :: bdt
REAL(DP)                                                   :: Chargesum
REAL(DP)                                                   :: sion_tmp
REAL(DP)                                                   :: aa1
REAL(DP)                                                   :: GamWaterCheck

REAL(DP)                                                   :: dhad,dhbd,tempk,tconv

INTEGER(I4B)                                               :: ik
INTEGER(I4B)                                               :: it
INTEGER(I4B)                                               :: ItPoint

CHARACTER (LEN=3)                                          :: ulabPrint

REAL(DP),PARAMETER                                    :: Avogadro=6.02252E23
REAL(DP),PARAMETER                                    :: BoltzmannConstant=1.36E-23
REAL(DP),PARAMETER                                    :: pi=3.14159265359
REAL(DP),PARAMETER                                    :: log_10=2.30258509299405
REAL(DP),PARAMETER                                    :: ConvertBarsToAtm=0.986923
REAL(DP)                                              :: U1
REAL(DP)                                              :: U2
REAL(DP)                                              :: U3
REAL(DP)                                              :: U4
REAL(DP)                                              :: U5
REAL(DP)                                              :: U6
REAL(DP)                                              :: U7
REAL(DP)                                              :: U8
REAL(DP)                                              :: U9
REAL(DP)                                              :: C_BradleyPitz
REAL(DP)                                              :: D1000_bradleyPitz
REAL(DP)                                              :: B_BradleyPitz
REAL(DP)                                              :: RgasAppelo

REAL(DP)                                              :: P_bars

REAL(DP)                                              :: eps_r

REAL(DP)                                              :: b1
REAL(DP)                                              :: b2
REAL(DP)                                              :: b3
REAL(DP)                                              :: b4
REAL(DP)                                              :: b5
REAL(DP)                                              :: b6

REAL(DP)                                              :: th
REAL(DP)                                              :: rho_0_sat
REAL(DP)                                              :: rho_0
REAL(DP)                                              :: p0
REAL(DP)                                              :: p1
REAL(DP)                                              :: p2
REAL(DP)                                              :: p3
REAL(DP)                                              :: p_sat
REAL(DP)                                              :: pa
REAL(DP)                                              :: e2_Dkt
REAL(DP)                                              :: DH_B
REAL(DP)                                              :: DH_A
REAL(DP)                                              :: tCritical
REAL(DP)                                              :: Av
REAL(DP)                                              :: PartialLogEps
REAL(DP)                                              :: act_H2O


ChargeSum = 0.0d0
TotalMoles = 0.0d0
tconv = 273.15d0
tempk = tempc + tconv

DO ik = 1,ncomp+nspec
  ulabPrint = ulab(ik)
!!!  IF (ulabPrint(1:3) /= 'H2O') THEN
    TotalMoles = TotalMoles + sptmp10(ik)
    ChargeSum = ChargeSum + sptmp10(ik)*chg(ik)*chg(ik)
!!!  ELSE
!!!    CONTINUE
!!!   END IF
END DO

sion_tmp = 0.50D0*ChargeSum

IF (sion_tmp < 25.0d0) THEN
  sqrt_sion = DSQRT(sion_tmp)
ELSE
  sion_tmp = 0.0d0
  sqrt_sion = 0.0d0
END IF

IF (ntemp == 1) THEN
  ah = adh(1)
  bh = bdh(1)
  bdt = bdot(1)
ELSE

  ItPoint = 0
  DO it = 1,ntemp
    IF (tempc == DatabaseTemperature(it)) THEN
      ItPoint = it
      ah = adh(ItPoint)
      bh = bdh(ItPoint)
      bdt = bdot(ItPoint)
    END IF
  END DO

  IF (ItPoint == 0) THEN
    ah = adhcoeff(1) + adhcoeff(2)*tempc  &
       + adhcoeff(3)*tempc*tempc + adhcoeff(4)*tempc*tempc*tempc  &
       + adhcoeff(5)*tempc*tempc*tempc*tempc
    bh = bdhcoeff(1) + bdhcoeff(2)*tempc  &
       + bdhcoeff(3)*tempc*tempc + bdhcoeff(4)*tempc*tempc*tempc  &
       + bdhcoeff(5)*tempc*tempc*tempc*tempc
    bdt = bdtcoeff(1) + bdtcoeff(2)*tempc  &
       + bdtcoeff(3)*tempc*tempc + bdtcoeff(4)*tempc*tempc*tempc  &
       + bdtcoeff(5)*tempc*tempc*tempc*tempc

  ELSE
    CONTINUE
  END IF

END IF

IF (Benchmark) THEN
  tconv = 273.15d0
  tempk = tempc + tconv
  call dhconst(ah,bh,tempk,tconv)
END IF

DO ik = 1,ncomp+nspec

  IF (chg(ik) == 0.0d0) THEN

    ulabPrint = ulab(ik)
    IF (ulabPrint(1:3) == 'H2O' .or. ulabPrint(1:3) == 'HHO') THEN
      
      gamWaterCheck = 1.0d0 - 0.017d0*TotalMoles
      gamtmp(IK) = clg*gamWaterCheck
      
!!!   Assumes molecular weight of H2O of 18.01528

    ELSE
        aa1 = 0.10d0*sion_tmp
        gamtmp(IK) = clg*aa1
    END IF

  ELSE

    IF (IncludeBdot) THEN  

      IF (acmp(ik) == 0.0d0) THEN      !! Davies equation
        aa1 = -0.5115d0*chg(IK)*chg(IK)* (sqrt_sion/(1.0d0 + sqrt_sion)  - 0.24d0*sion_tmp  )
      ELSE                  !! WATEQ extended Debye-Huckel

        aa1 = -(ah*chg(IK)*chg(IK)*sqrt_sion)/            &
              (1.0d0 + acmp(IK)*bh*sqrt_sion)                   &         
              + bdotParameter(ik)*sion_tmp
      END IF
    ELSE                    !!  Helgesonian-LLNL bdot expression based on extended Debye-Huckel

      aa1 = -(ah*chg(IK)*chg(IK)*sqrt_sion)/              &
              (1.0d0 + acmp(IK)*bh*sqrt_sion)                  &         
              + bdt*sion_tmp
    END IF

    gamtmp(IK) = clg*aa1

  END IF

END DO

IF (SaltCreep) THEN
  
  U1 = 3.4279E2
  U2 = -5.0866E-3
  U3 = 9.4690E-7
  U4 = -2.0525
  U5 = 3.1159E3
  U6 = -1.8289E2
  U7 = -8.0325E3
  U8 = 4.2142E6
  U9 = 2.1417
  
  C_BradleyPitz = U4 + U5/(U6 + tempk)
  D1000_bradleyPitz = U1 *EXP(U2*tempk + U3*tempk*tempk)
  B_BradleyPitz = U7 + U8/tempk + U9*tempk
  
  pa = 1.0
  P_bars = pa/ConvertBarsToAtm                        !! Conversion to bars
  RgasAppelo = 82.0574587      !! (atm cm^3 mol^-1 K^-1)

  eps_r = D1000_bradleyPitz + C_BradleyPitz * Log( (B_BradleyPitz + P_bars)/(B_BradleyPitz + 1000.0) )  
    
  tCritical = 647.096                      
  th = 1 - tempK / tCritical
  
  b1 = 1.99274064
  b2 = 1.09965342
  b3 = -0.510839303
  b4 = -1.75493479
  b5 = -45.5170352
  b6 = -6.7469445e5
  
  rho_0_sat = 322.0 * (1.0 + b1 * th**0.3333 + b2 * th**0.66666 + b3 * th**(5.0/3.0) + b4 * th**(16.0/3.0) +   &
                b5 * th**(43.0/3.0) + b6 * th**(110.0/3.0) )
  
  p0 =  5.1880000E-02 + tempc * (-4.1885519E-04 + tempc * ( 6.6780748E-06 + tempc * (-3.6648699E-08 + tempc *  8.3501912E-11)))
  p1 = -6.0251348E-06 + tempc * ( 3.6696407E-07 + tempc * (-9.2056269E-09 + tempc * ( 6.7024182E-11 + tempc * -1.5947241E-13)))
  p2 = -2.2983596E-09 + tempc * (-4.0133819E-10 + tempc * ( 1.2619821E-11 + tempc * (-9.8952363E-14 + tempc *  2.3363281E-16)))
  p3 =  7.0517647E-11 + tempc * ( 6.8566831E-12 + tempc * (-2.2829750E-13 + tempc * ( 1.8113313E-15 + tempc * -4.2475324E-18)))
  
  gamWaterCheck = 1.0d0 - 0.017d0*TotalMoles
      
  if (gamwatercheck <= 1.0) then
    p_sat = exp(11.6702 - 3816.44 / (tempk - 46.13)) * gamwatercheck
  else
    p_sat = exp(11.6702 - 3816.44 / (tempk - 46.13))
  end if  

  if (pa < p_sat) then
    pa = p_sat
  else
    pa = pa - (p_sat - 1e-6)
  end if

  rho_0 = rho_0_sat + pa * (p0 + pa * (p1 + pa * (p2 + dsqrt(pa) * p3)))

  e2_DkT = 1.671008d-3 / (eps_r * tempk)
  

  DH_B = dsqrt( 8.0 * pi * AVOGADRO * e2_DkT * rho_0 / 1.0E06 )     
     !!! // Debye length parameter, 1/cm(mol/kg)^-0.5

  DH_A = DH_B * e2_DkT / (2.0d0 * clg)
  
  DH_B = DH_B * 1.0E-08
  
!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
  
  aa1 = -(DH_A * chg(ikNa)*chg(ikNa)*sqrt_sion)/                      &
              (1.0d0 + 4.08 * DH_B * sqrt_sion)                       &     
              + 0.082*sion_tmp
  gamtmp(ikNa) = clg*aa1
  
  aa1 = -(DH_A  * chg(ikCl)*chg(ikCl)*sqrt_sion)/                      &
              (1.0d0 + 3.63 * DH_B * sqrt_sion)                        &  
              + 0.017*sion_tmp
  gamtmp(ikCl) = clg*aa1
  
END IF

RETURN
END SUBROUTINE gamma_init
!*****************************************************************
