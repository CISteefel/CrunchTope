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

SUBROUTINE keqcalc2(ncomp,nrct,nspec,ngas,nsurf_sec,jx,jy,jz)
USE crunchtype
USE params
USE concentration
USE mineral
USE temperature
USE medium, ONLY: stress

IMPLICIT NONE

!  External variables

INTEGER(I4B), INTENT(IN)                             :: ncomp
INTEGER(I4B), INTENT(IN)                             :: nrct
INTEGER(I4B), INTENT(IN)                             :: nspec
INTEGER(I4B), INTENT(IN)                             :: ngas
INTEGER(I4B), INTENT(IN)                             :: nsurf_sec
INTEGER(I4B), INTENT(IN)                             :: jx
INTEGER(I4B), INTENT(IN)                             :: jy
INTEGER(I4B), INTENT(IN)                             :: jz

!  Internal variables

REAL(DP)                                             :: temp
REAL(DP)                                             :: temp2
REAL(DP)                                             :: x1
REAL(DP)                                             :: x2
REAL(DP)                                             :: x3
REAL(DP)                                             :: x4
REAL(DP)                                             :: x5

REAL(DP)                                              :: DeltaV_r
REAL(DP),PARAMETER                                    :: AvogadroNumber=6.02214E23
REAL(DP),PARAMETER                                    :: BoltzmannConstant=1.36E-23
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
REAL(DP)                                              :: partialEpsInverse
REAL(DP)                                              :: partialLogEps
REAL(DP)                                              :: Vm_Na
REAL(DP)                                              :: Vm_Cl
REAL(DP)                                              :: Vm_NaCl
REAL(DP)                                              :: Vm_Na_zero
REAL(DP)                                              :: Vm_Cl_zero
REAL(DP)                                              :: bh
REAL(DP)                                              :: Chargesum
REAL(DP)                                              :: sqrt_sion
REAL(DP)                                              :: sion_tmp
REAL(DP)                                              :: P_appelo
REAL(DP)                                              :: P_bars
REAL(DP)                                              :: Av
REAL(DP)                                              :: eps_r
REAL(DP)                                              :: RgasAppelo
REAL(DP)                                              :: ConvertBarsToAtm
REAL(DP)                                              :: ConvertPaToBars

INTEGER(I4B)                                         :: ksp
INTEGER(I4B)                                         :: kg
INTEGER(I4B)                                         :: k
INTEGER(I4B)                                         :: ik
INTEGER(I4B)                                         :: msub
INTEGER(I4B)                                         :: ns
INTEGER(I4B)                                         :: np

temp = t(jx,jy,jz) + 273.15
temp2 = temp*temp

ChargeSum = 0.0d0

DO ik = 1,ncomp+nspec
  ChargeSum = ChargeSum + sp10(ik,jx,jy,jz)*chg(ik)*chg(ik)
END DO

sion_tmp = 0.50D0*ChargeSum
sqrt_sion = DSQRT(sion_tmp)

DO ksp = 1,nspec
  IF (ntemp == 1 .OR. RunIsothermal) THEN
    keqaq(ksp,jx,jy,jz) = -clg*eqhom(ksp)
  ELSE
    x1 = as1(ksp,1)
    x2 = as1(ksp,2)
    x3 = as1(ksp,3)
    x4 = as1(ksp,4)
    x5 = as1(ksp,5) 
    keqaq(ksp,jx,jy,jz) = -clg*(x1*DLOG(temp) + x2 + x3*temp + x4/temp + x5/(temp2))
  END IF
END DO

DO kg = 1,ngas
  ksp = kg + nspec
  IF (ntemp == 1 .OR. RunIsothermal) THEN
    keqgas(kg,jx,jy,jz) = -clg*eqgas(kg)
  ELSE
    x1 = as1(ksp,1)
    x2 = as1(ksp,2)
    x3 = as1(ksp,3)
    x4 = as1(ksp,4)
    x5 = as1(ksp,5)
    keqgas(kg,jx,jy,jz) = -clg*(x1*DLOG(temp) + x2 + x3*temp + x4/temp + x5/(temp2))
  END IF
END DO

msub = 0
DO k = 1,nrct
  DO np = 1,nreactmin(k)
    msub = msub + 1
    ksp = msub + ngas + nspec
    IF (ntemp == 1 .OR. RunIsothermal) THEN
      keqmin(np,k,jx,jy,jz) = clg*alnk(msub)
    ELSE
      x1 = as1(ksp,1)
      x2 = as1(ksp,2)
      x3 = as1(ksp,3)
      x4 = as1(ksp,4)
      x5 = as1(ksp,5)
      keqmin(np,k,jx,jy,jz) = clg*(x1*DLOG(temp) + x2 +  &
          x3*temp + x4/temp + x5/(temp2))
    END IF
  END DO
END DO

!!!  ****************************************************************************
!!!  For SaltCreep

IF (SaltCreep) THEN
  
  ConvertBarsToAtm = 0.986923
  ConvertPaToBars = 1.0E-05
  RgasAppelo = 82.0574587      !! (atm cm^3 mol^-1 K^-1)

  
  P_appelo = ABS(stress(jx,jy,jz)) * ConvertBarsToAtm * ConvertPaToBars   !! Conversion to ATM
  P_bars = ABS(stress(jx,jy,jz)) * ConvertPaToBars                        !! Conversion to bars
  
  IF (P_bars > 2000.0) THEN
    P_bars = 2000.0
  END IF
  
  bh = 0.3288
  
  U1 = 3.4279E2
  U2 = -5.0866E-3
  U3 = 9.4690E-7
  U4 = -2.0525
  U5 = 3.1159E3
  U6 = -1.8289E2
  U7 = -8.0325E3
  U8 = 4.2142E6
  U9 = 2.1417
  
  C_BradleyPitz = U4 + U5/(U6 + temp)
  D1000_bradleyPitz = U1 *EXP(U2*temp + U3*temp*temp)
  B_BradleyPitz = U7 + U8/temp + U9*temp
  
  eps_r = D1000_bradleyPitz + C_BradleyPitz * Log( (B_BradleyPitz + P_bars)/(B_BradleyPitz+1000.0) )  

  partialLogEps = C_BradleyPitz/       &
      (  (B_BradleyPitz + P_bars) * ( C_BradleyPitz * Log( (B_BradleyPitz + P_bars)/(B_BradleyPitz+1000.0) ) +   &
          D1000_bradleyPitz )  )
  
  partialEpsInverse = C_BradleyPitz/(  (B_BradleyPitz + P_bars) *    &
              ( C_BradleyPitz * Log( (B_BradleyPitz + P_bars)/(B_BradleyPitz+1000.0) ) + D1000_bradleyPitz )**2  )
  
  
  Av = (RgasAppelo * temp * 0.5114 * 2.0/3.0 * 2.303 * (3.0 * partialLogEps - 4.52E-05))
  
  Vm_Na_zero = 41.84 * ( 0.1 * 2.28 + 100.0*(-4.38)/(2600.0 + P_appelo) + (-4.10)/(temp-228.0) +  &
                 10000.0*(-0.586)/ ( (2600.0 + P_appelo)*(temp-228.0) )  - 0.09*1.0E+05 * partialEpsInverse )
  Vm_Cl_zero = 41.84 * (0.1 * 4.465 + 100.0*(4.801)/(2600.0 + P_appelo) + (4.325)/(temp-228.0) +   &
                 10000.0*(-2.847)/( (2600.0 + P_appelo)*(temp-228.0) ) - 1.748*1.0E+05 * partialEpsInverse)
  
  Vm_Na = Vm_Na_zero + Av*0.5*chg(ikNa)*chg(ikNa) * sqrt_sion/(1.0d0 + 4.0*0.3288*sqrt_sion)     &
                 + ( 0.30 + 52.0/(temp - 228.0)   + -3.33E-03 * (temp - 228.0) )*sion_tmp**0.566
  Vm_Cl = Vm_Cl_zero + Av*0.5*chg(ikCl)*chg(ikCl) * sqrt_sion/(1.0d0 + 0.0*sqrt_sion)     &
                + ( -0.331 + 20.16/(temp - 228.0) + 0.0 * (temp - 228.0) )*sion_tmp**1.00
  
  Vm_NaCl = 27.1

  DeltaV_r = Vm_Na + Vm_Cl - Vm_NaCl
  
  keqmin(1,1,jx,jy,jz) = ( 1.570 - DeltaV_r * (P_appelo - 1.0)/(2.303*RgasAppelo*temp) )
  
!!!  write(*,*) 'Pressure (atm) = ', P_appelo
  
!!!    write(*,*) ' DeltaV_r =   ', DeltaV_r
!!!    write(*,*) ' keqmin_tmp = ', keqmin_tmp(1,1)
!!!    write(*,*)
!!!    write(*,*)
!!!    read(*,*)
  
  keqmin(1,1,jx,jy,jz) = clg*keqmin(1,1,jx,jy,jz) 
  
END IF
!!!  *********************************************************

DO ns = 1,nsurf_sec
  ksp = msub + ngas + nspec + ns
  IF (ntemp == 1 .OR. RunIsothermal) THEN
    keqsurf(ns,jx,jy,jz) = -clg*eqsurf(ns)
  ELSE
    x1 = as1(ksp,1)
    x2 = as1(ksp,2)
    x3 = as1(ksp,3)
    x4 = as1(ksp,4)
    x5 = as1(ksp,5)
    keqsurf(ns,jx,jy,jz) = -clg*(x1*DLOG(temp) + x2 + x3*temp + x4/temp + x5/(temp2))
  END IF
END DO

RETURN
END SUBROUTINE keqcalc2
!*********************************************************************
