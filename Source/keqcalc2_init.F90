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

SUBROUTINE keqcalc2_init(ncomp,nrct,nspec,ngas,nsurf_sec,tempc,PressureTemp)
USE crunchtype
USE params
USE concentration
USE mineral
USE temperature
USE medium, ONLY: PressureFound
USE CrunchFunctions

IMPLICIT NONE

!  External variables

INTEGER(I4B), INTENT(IN)                             :: ncomp
INTEGER(I4B), INTENT(IN)                             :: nrct
INTEGER(I4B), INTENT(IN)                             :: nspec
INTEGER(I4B), INTENT(IN)                             :: ngas
INTEGER(I4B), INTENT(IN)                             :: nsurf_sec

REAL(DP), INTENT(IN)                                 :: tempc
REAL(DP), INTENT(IN)                                 :: PressureTemp


!  Internal variables

REAL(DP)                                             :: temp
REAL(DP)                                             :: temp2
REAL(DP)                                             :: x1
REAL(DP)                                             :: x2
REAL(DP)                                             :: x3
REAL(DP)                                             :: x4
REAL(DP)                                             :: x5

REAL(DP)                                              :: DeltaV_r1
REAL(DP)                                              :: DeltaV_r2
REAL(DP)                                              :: DeltaV_r3
REAL(DP)                                              :: DeltaV_r4
REAL(DP)                                              :: DeltaV_r5

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
REAL(DP)                                              :: Vm_CaCO3
REAL(DP)                                              :: Vm_Na_zero
REAL(DP)                                              :: Vm_Cl_zero

REAL(DP)                                              :: Vm_Ca
REAL(DP)                                              :: Vm_CO3
REAL(DP)                                              :: Vm_HCO3
REAL(DP)                                              :: Vm_CO2
REAL(DP)                                              :: Vm_OH

REAL(DP)                                              :: Vm_Ca_zero
REAL(DP)                                              :: Vm_CO3_zero
REAL(DP)                                              :: Vm_HCO3_zero
REAL(DP)                                              :: Vm_CO2_zero
REAL(DP)                                              :: Vm_OH_zero


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
REAL(DP)                                              :: kReal
REAL(DP)                                              :: tempT

REAL(DP)                                              :: A_temp1
REAL(DP)                                              :: A_temp2
REAL(DP)                                              :: A_temp3
REAL(DP)                                              :: A_temp4
REAL(DP)                                              :: A_temp5
REAL(DP)                                              :: A_temp6

REAL(DP)                                              :: PressureCorrection
REAL(DP)                                              :: PressureCorrection1
REAL(DP)                                              :: PressureCorrection2
REAL(DP)                                              :: PressureCorrection3
REAL(DP)                                              :: PressureCorrection4
REAL(DP)                                              :: PressureCorrection5

REAL(DP)                                              :: MolarEnergy

CHARACTER*1          :: string1
CHARACTER*2          :: string2
CHARACTER*3          :: string3
CHARACTER*4          :: string4

CHARACTER*5          :: StringReal


INTEGER(I4B)                                         :: ksp
INTEGER(I4B)                                         :: ik
INTEGER(I4B)                                         :: kg
INTEGER(I4B)                                         :: k
INTEGER(I4B)                                         :: msub
INTEGER(I4B)                                         :: ns
INTEGER(I4B)                                         :: np
INTEGER(I4B)                                         :: kkk
INTEGER(I4B)                                         :: ikOH
INTEGER(I4B)                                         :: ikHCO3

LOGICAL(LGT)                                         :: WriteInputP
LOGICAL(LGT)                                         :: WriteInputT

temp = tempc + 273.15
temp2 = temp*temp

WriteInputP = .FALSE.
WriteInputT = .FALSE.

ChargeSum = 0.0d0

DO ik = 1,ncomp+nspec
  ChargeSum = ChargeSum + sptmp10(ik)*chg(ik)*chg(ik)
END DO

sion_tmp = 0.50D0*ChargeSum
sqrt_sion = DSQRT(sion_tmp)

DO ksp = 1,nspec
  IF (ntemp == 1 .OR. RunIsothermal) THEN
    keqaq_tmp(ksp) = -clg*eqhom(ksp)
  ELSE
    x1 = as1(ksp,1)
    x2 = as1(ksp,2)
    x3 = as1(ksp,3)
    x4 = as1(ksp,4)
    x5 = as1(ksp,5)
    keqaq_tmp(ksp) = -clg*(x1*DLOG(temp) + x2 +  &
        x3*temp + x4/temp + x5/(temp2))
  END IF
END DO

DO kg = 1,ngas
  ksp = kg + nspec
  IF (ntemp == 1 .OR. RunIsothermal) THEN
    keqgas_tmp(kg) = -clg*eqgas(kg)
  ELSE
    x1 = as1(ksp,1)
    x2 = as1(ksp,2)
    x3 = as1(ksp,3)
    x4 = as1(ksp,4)
    x5 = as1(ksp,5)
    keqgas_tmp(kg) = -clg*(x1*DLOG(temp) + x2 +  &
        x3*temp + x4/temp + x5/(temp2))
  END IF
END DO

msub = 0
DO k = 1,nrct
  DO np = 1,nreactmin(k)
    msub = msub + 1
    ksp = msub + ngas + nspec
    IF (ntemp == 1 .OR. RunIsothermal) THEN
      keqmin_tmp(np,k) = clg*alnk(msub)
    ELSE
      x1 = as1(ksp,1)
      x2 = as1(ksp,2)
      x3 = as1(ksp,3)
      x4 = as1(ksp,4)
      x5 = as1(ksp,5)
      keqmin_tmp(np,k) = clg*(x1*DLOG(temp) + x2 +  &
          x3*temp + x4/temp + x5/(temp2))
    END IF
  END DO
END DO

!!  For SaltCreep
IF (SaltCreep) THEN
    
  ConvertBarsToAtm = 0.986923
  ConvertPaToBars = 1.0E-05
  RgasAppelo = 82.0574587      !! (atm cm^3 mol^-1 K^-1)
  
  if (WriteInputP) THEN
    
  do kkk = 1,3500,2
    kReal = FLOAT(kkk)
    
    IF (kkk < 10) THEN
      
      string1 = IntegerToCharacter(kkk)
       
      write(122,111) string1(1:1)
      write(122,112) string1(1:1)
      
111   format('Condition Halite',a1)
112   format('Pressure    ',a1)
      
    ELSE IF (kkk >= 10 .and. kkk<100) THEN
      
      string2 = IntegerToCharacter(kkk)
      write(122,113) string2
      write(122,114) string2
113   format('Condition Halite',a2)
114   format('Pressure    ',a2)
      
    ELSE IF (kkk >= 100 .and. kkk<1000) THEN
      string3 = IntegerToCharacter(kkk)
      write(122,115) string3
      write(122,116) string3
115   format('Condition Halite',a3)
116   format('Pressure    ',a3) 

    ELSE
      string4 = IntegerToCharacter(kkk)
      write(122,117) string4
      write(122,118) string4
117   format('Condition Halite',a4)
118   format('Pressure    ',a4) 
    ENDIF

    write(122,*) 'temperature     25.0'
    write(122,*) 'set_porosity    0.001'
    write(122,*) 'units           mol/kg'
    write(122,*) 'Na+             Halite'
    write(122,*) 'Cl-             charge'
    write(122,*) 'H2O             55.5'
    write(122,*) 'Tracer          1.0E-06'
    write(122,*) 'Halite          0.10  bsa 1.0'
    write(122,*) 'END'
    write(122,*)

  END DO
  
  ENDIF
  
  if (WriteInputT) THEN
    
  do kkk = 25,200
    
    IF (kkk < 10) THEN
      
      string1 = IntegerToCharacter(kkk)
      tempT = Real(kkk) 
!!!      tempT = tempT*0.2 + 22.0
      temp = tempT + 273.15
      StringReal = RealToCharacter(tempT)
      write(122,111) string1(1:1)
      write(122,121) stringReal
      
    ELSE IF (kkk >= 10 .and. kkk<100) THEN
      
      string2 = IntegerToCharacter(kkk)
      tempT = Real(kkk)
!!!      tempT = tempT*0.2 + 22.0
      temp = tempT + 273.15
      StringReal = RealToCharacter(tempT)
      write(122,113) string2
      write(122,121) stringReal
      
    ELSE IF (kkk >= 100 .and. kkk<1000) THEN
      
      string3 = IntegerToCharacter(kkk)
      tempT = Real(kkk)
!!!      tempT = tempT*0.2 + 22.0
      temp = tempT + 273.15
      StringReal = RealToCharacter(tempT)
      write(122,115) string3
      write(122,121) stringReal
      
    ELSE
      
      write(*,*) 'Temperatures should not be above 100C'
      write(*,*)
      read(*,*)
      stop
      
    ENDIF
    
119   format('Temperature    ',a1)
120   format('Temperature    ',a2)
121   format('Temperature    ',a5) 
      
    write(122,*) 'pressure        1.0'
    write(122,*) 'set_porosity    0.001'
    write(122,*) 'units           mol/kg'
    write(122,*) 'Na+             Halite'
    write(122,*) 'Cl-             charge'
    write(122,*) 'Tracer          1.0E-06'
    write(122,*) 'H2O             55.55'
    write(122,*) 'Halite          0.10  bsa 1.0'
    write(122,*) 'END'
    write(122,*)

  END DO
  
  ENDIF
  
  
  IF (PressureFound) THEN
    P_bars = PressureTemp
  ELSE
    P_bars = 1.0
  END IF
  
  IF (P_bars > 2000.0) THEN
    P_bars = 2000.0
  END IF
  
  write(*,*) ' Pressure (bars) = ',P_bars
  
  
  P_appelo = P_bars * ConvertBarsToAtm                        !! Conversion to bars
  
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
  
  Av = ( RgasAppelo * temp * 0.5114 * 2.0/3.0 * 2.303 * (3.0 * partialLogEps - 4.52E-05) )
  
 !!! Vm(T, pb, I) = 41.84 * (a1 * 0.1 + a2 * 100 / (2600 + pb)  + a3 / (T - 228) +  &
 !!!                       a4 * 1e4 / (2600 + pb) / (T - 228) - W * QBrn ) +         &
 !!!                       z^2 / 2 * Av * f(I^0.5) + (i1 + i2 / (T - 228) + i3 * (T - 228)) * I^i4
  
!!          -Vm   a1     a2     a3      a4     W      a0  i1     i2    i3         i4
!!! Na+     -Vm   2.28   -4.38  -4.1    -0.586  0.09   4   0.3    52     -3.33e-3   0.45
!!! Cl-     -Vm   4.465   4.801  4.325  -2.847  1.748  0  -0.331  20.16   0         1 
!!! Ca      -Vm  -0.3456 -7.252  6.149  -2.479  2.139  5   1.60  -57.1   -6.12E-03    1
!!! CO3     -Vm   4.91    0.00   0.00   -5.41   4.76   0   0.386  89.7   -1.57e-2   1
!!! HCO3    -Vm   8.54    0.00  -11.7    0.00   1.60   0   0.000  116.0   0.000     1
!!! CO2(aq) -Vm   20.85  -46.93 -46.93   27.9  -0.193  0   0.000  0       0.000     0
  
  Vm_Na_zero = 41.84 * ( 0.1 * 2.28 + 100.0*(-4.38)/(2600.0 + P_appelo) + (-4.10)/(temp-228.0) +  &
                 10000.0*(-0.586)/ ( (2600.0 + P_appelo)*(temp-228.0) )  - 0.09*1.0E+05 * partialEpsInverse )
  
  Vm_Cl_zero = 41.84 * (0.1 * 4.465 + 100.0*(4.801)/(2600.0 + P_appelo) + (4.325)/(temp-228.0) +   &
                 10000.0*(-2.847)/( (2600.0 + P_appelo)*(temp-228.0) ) - 1.748*1.0E+05 * partialEpsInverse)
  
  !!! +++++++++++++++++++++++++++++++++++++++
  
  Vm_Na = Vm_Na_zero + Av*0.5*chg(ikNa)*chg(ikNa) * sqrt_sion/(1.0d0 + 4.0*0.3288*sqrt_sion)     &
                 + ( 0.30 + 52.0/(temp - 228.0)      + (-3.33E-03) * (temp - 228.0) )*sion_tmp**0.566
  
  Vm_Cl = Vm_Cl_zero + Av*0.5*chg(ikCl)*chg(ikCl) * sqrt_sion/(1.0d0 + 0.0*sqrt_sion)     &
                + ( -0.331 + 20.16/(temp - 228.0)    + 0.0 * (temp - 228.0) )*sion_tmp**1.00
  
  Vm_NaCl = 27.1

  DeltaV_r1 = Vm_Na + Vm_Cl - Vm_NaCl
  
  
  A_temp1 = -713.4616
  A_temp2 = -.1201241
  A_temp3 = 37302.21
  A_temp4 = 262.4583
  A_temp5 = -2106915.
  A_temp6 = 0.0
  
  MolarEnergy = temp * 0.00831470
  
  !!! PHREEQc van't Hoff equation = logK_25 - Delta H * (298.15 - temp) /(2.303 * 298.15 * MolarEnergy)
  keqmin_tmp(1,1) = 1.570 - 1.37 *(298.15 - temp)/(2.303 * 298.15 * MolarEnergy)
  
  !!!PHREEQc log Ks
  
 !!! keqmin_tmp(1,1) = (A_temp1 + A_temp2*temp + A_temp3/temp + A_temp4*LOG10(temp) +    &
 !!!                     A_temp5/(temp*temp) + A_temp6*temp*temp )
  
  PressureCorrection = DeltaV_r1 * (P_appelo - 1.0)/(2.303*RgasAppelo*temp)
  
!!!  keqmin_tmp(1,1) = ( keqmin_tmp(1,1)/clg  - PressureCorrection )
  keqmin_tmp(1,1) = ( keqmin_tmp(1,1)  - PressureCorrection )
  
!!!    write(*,*) ' DeltaV_r =   ', DeltaV_r
!!!    write(*,*) ' keqmin_tmp = ', keqmin_tmp(1,1)
!!!    write(*,*)
!!!    write(*,*)
!!!    read(*,*)
  
  keqmin_tmp(1,1) = clg*keqmin_tmp(1,1)
  
END IF

!!  For CalciteCreep
IF (CalciteCreep) THEN
    
  ConvertBarsToAtm = 0.986923
  ConvertPaToBars = 1.0E-05
  RgasAppelo = 82.0574587      !! (atm cm^3 mol^-1 K^-1)
  
  
  IF (PressureFound) THEN
    P_bars = PressureTemp
  ELSE
    P_bars = 1.0
  END IF
  
  IF (P_bars > 2000.0) THEN
    P_bars = 2000.0
  END IF
  
  P_appelo = P_bars * ConvertBarsToAtm                        !! Conversion to bars
  
  write(*,*) ' Pressure (bars) = ',P_bars, P_appelo
  write(*,*)
  
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
  
  Av = ( RgasAppelo * temp * 0.5114 * 2.0/3.0 * 2.303 * (3.0 * partialLogEps - 4.52E-05) )
  
 !!! Vm(T, pb, I) = 41.84 * (a1 * 0.1 + a2 * 100 / (2600 + pb)  + a3 / (T - 228) +  &
 !!!                       a4 * 1e4 / (2600 + pb) / (T - 228) - W * QBrn ) +         &
 !!!                       z^2 / 2 * Av * f(I^0.5) + (i1 + i2 / (T - 228) + i3 * (T - 228)) * I^i4
  
!!          -Vm   a1     a2     a3      a4     W      a0  i1     i2    i3         i4
!!! Na+     -Vm   2.28   -4.38  -4.1    -0.586  0.09   4   0.3    52     -3.33e-3   0.45
!!! Cl-     -Vm   4.465   4.801  4.325  -2.847  1.748  0  -0.331  20.16   0         1 
!!! Ca      -Vm  -0.3456 -7.252  6.149  -2.479  1.239  5   1.60  -57.1   -6.12E-03  1
!!! CO3     -Vm   4.91    0.00   0.00   -5.41   4.76   0   0.386  89.7   -1.57e-2   1
!!! HCO3    -Vm   8.54    0.00  -11.7    0.00   1.60   0   0.000  116.0   0.000     1
!!! CO2(aq) -Vm   20.85  -46.93 -79.0   27.9  -0.193  0   0.000  0.00     0.000    0
!!! OH-     -Vm   -9.66   28.5   80.0    -22.9   1.89  0   1.09   0.00     0.00     1
  
  Vm_Na_zero = 41.84 * ( 0.1 * 2.28 + 100.0*(-4.38)/(2600.0 + P_appelo) + (-4.10)/(temp-228.0) +  &
                 10000.0*(-0.586)/ ( (2600.0 + P_appelo)*(temp-228.0) )  - 0.09*1.0E+05 * partialEpsInverse )
  
  Vm_Cl_zero = 41.84 * (0.1 * 4.465 + 100.0*(4.801)/(2600.0 + P_appelo) + (4.325)/(temp-228.0) +   &
                 10000.0*(-2.847)/( (2600.0 + P_appelo)*(temp-228.0) ) - 1.748*1.0E+05 * partialEpsInverse)
  
  Vm_Ca_zero = 41.84 * ( 0.1 * (-0.3456) + 100.0*(-7.252)/(2600.0 + P_appelo) + (6.149)/(temp-228.0) +  &
                 10000.0*(-2.479)/ ( (2600.0 + P_appelo)*(temp-228.0) )  - 1.239*1.0E+05 * partialEpsInverse )
  
  Vm_CO3_zero = 41.84 * ( 0.1 * (4.91) + 100.0*(0.00)/(2600.0 + P_appelo) + (0.00)/(temp-228.0) +  &
                 10000.0*(-5.41)/ ( (2600.0 + P_appelo)*(temp-228.0) )  - 4.76*1.0E+05 * partialEpsInverse )
  
  Vm_HCO3_zero = 41.84 * ( 0.1 * (8.54) + 100.0*(0.00)/(2600.0 + P_appelo) + (-11.7)/(temp-228.0) +  &
                 10000.0*(0.00)/ ( (2600.0 + P_appelo)*(temp-228.0) )  - 1.60*1.0E+05 * partialEpsInverse )
  
  Vm_CO2_zero = 41.84 * ( 0.1 * (20.85) + 100.0*(-46.93)/(2600.0 + P_appelo) + (-79.0)/(temp-228.0) +  &
                 10000.0*(27.9)/ ( (2600.0 + P_appelo)*(temp-228.0) )  - (-0.193)*1.0E+05 * partialEpsInverse )
  
  Vm_OH_zero  = 41.84 * ( 0.1 * (-9.66) + 100.0*( 28.50)/(2600.0 + P_appelo) + ( 80.0)/(temp-228.0) +  &
                 10000.0*(-22.9)/ ( (2600.0 + P_appelo)*(temp-228.0) )  - (1.89)*1.0E+05 * partialEpsInverse )
  
  write(*,*)
  write(*,*) ' Intrinsic molar volumes at 25C'
  write(*,*) ' Vm_Na_zero =', Vm_Na_zero
  write(*,*) ' Vm_Cl_zero =', Vm_Cl_zero
  write(*,*) ' Vm_Ca_zero =', Vm_Ca_zero
  write(*,*) ' Vm_CO3_zero =', Vm_CO3_zero
  write(*,*) ' Vm_HCO3_zero =', Vm_HCO3_zero
  write(*,*) ' Vm_CO2_zero =', Vm_CO2_zero
  write(*,*) ' Vm_OH_zero =', Vm_OH_zero
  write(*,*)
!!!  read(*,*)
  
  !!! +++++++++++++++++++++++++++++++++++++++

  ikCa = 2
  ikCO3 = -2
  ikHCO3 = -1.0
  ikOH = -1.0

!!! PRIMARY_SPECIES
!!! H+
!!! CO2(aq)
!!! Ca++
!!! Na+
!!! Cl-
!!! Tracer
!!! Bogus
!!! END

!!! SECONDARY_SPECIES
!!! OH-
!!! HCO3-
!!! CO3--
!!! END
  
  Vm_Na = Vm_Na_zero + Av*0.5*chg(ikNa)*chg(ikNa) * sqrt_sion/(1.0d0 + 4.0*0.3288*sqrt_sion)     &
                 + ( 0.30 + 52.0/(temp - 228.0)      + (-3.33E-03) * (temp - 228.0) )*sion_tmp**0.566
  
  Vm_Cl = Vm_Cl_zero + Av*0.5*chg(ikCl)*chg(ikCl) * sqrt_sion/(1.0d0 + 0.0*sqrt_sion)     &
                + ( -0.331 + 20.16/(temp - 228.0)    + 0.0 * (temp - 228.0) )*sion_tmp**1.00
  
  Vm_Ca = Vm_Ca_zero + Av*0.5*(2.0)*(2.0) * sqrt_sion/(1.0d0 + 5.0*0.3288*sqrt_sion)     &
                 + ( 1.60 + (-57.1)/(temp - 228.0)   + (-6.12E-03) * (temp - 228.0) )*sion_tmp**1.00
  
  Vm_CO3 = Vm_CO3_zero + Av*0.5*(-2.0)*(-2.0) * sqrt_sion/(1.0d0 + 0.0*0.3288*sqrt_sion)     &
                 + ( 0.386 + (89.7)/(temp - 228.0)   + (-1.57-02) * (temp - 228.0) )*sion_tmp**1.00
  
  Vm_HCO3 = Vm_HCO3_zero + Av*0.5*(-1.0)*(-1.0) * sqrt_sion/(1.0d0 + 0.0*0.3288*sqrt_sion)     &
                 + ( 0.0 + (116.0)/(temp - 228.0)    + (0.00) * (temp - 228.0) )*sion_tmp**1.00
  
  VM_CO2 = Vm_CO2_zero 
  
  Vm_OH = Vm_OH_zero + Av*0.5*(-1.0)*(-1.0) * sqrt_sion/(1.0d0 + 0.0*0.3288*sqrt_sion)     &
                 + ( 1.09 + (0.0)/(temp - 228.0)    + (0.00) * (temp - 228.0) )*sion_tmp**1.00
  
  Vm_NaCl = 27.015
  Vm_CaCO3 = 36.9340

  DeltaV_r1 = Vm_Na + Vm_Cl  - Vm_NaCl
  DeltaV_r2 = Vm_Ca + Vm_CO3 - Vm_CaCO3
  DeltaV_r3 = Vm_CO2 - Vm_HCO3  !! Apparent molal volume of H+ assumed = 0
  DeltaV_r4 = Vm_CO2 - Vm_CO3   !! Apparent molal volume of H+ assumed = 0
  DeltaV_r5 = -Vm_OH            !! Apparent molal volume of H+ assumed = 0
  
  !! +++++++++++++++++++++++++
  
  A_temp1 = -713.4616
  A_temp2 = -.1201241
  A_temp3 = 37302.21
  A_temp4 = 262.4583
  A_temp5 = -2106915.
  A_temp6 = 0.0
  
  MolarEnergy = temp * 0.00831470
  
  !!! PHREEQc van't Hoff equation = logK_25 - Delta H * (298.15 - temp) /(2.303 * 298.15 * MolarEnergy)
  !!! keqmin_tmp(1,1) = 1.570 - 1.37 *(298.15 - temp)/(2.303 * 298.15 * MolarEnergy)
  
  !!!PHREEQc log Ks
  
 !!! keqmin_tmp(1,1) = (A_temp1 + A_temp2*temp + A_temp3/temp + A_temp4*LOG10(temp) +    &
 !!!                     A_temp5/(temp*temp) + A_temp6*temp*temp )
  
 !!! ++++++++++++++++++++++++++++++++++++++++++++
  
  PressureCorrection1 = DeltaV_r1 * (P_appelo - 1.0)/(2.303*RgasAppelo*temp)   !! Halite
  PressureCorrection2 = DeltaV_r2 * (P_appelo - 1.0)/(2.303*RgasAppelo*temp)   !! Calcite
  PressureCorrection3 = DeltaV_r3 * (P_appelo - 1.0)/(2.303*RgasAppelo*temp)   !! HCO3-
  PressureCorrection4 = DeltaV_r4 * (P_appelo - 1.0)/(2.303*RgasAppelo*temp)   !! CO3--
  PressureCorrection5 = DeltaV_r5 * (P_appelo - 1.0)/(2.303*RgasAppelo*temp)   !! OH-
  
  keqmin_tmp(1,2) = keqmin_tmp(1,2)/clg
  keqmin_tmp(2,2) = keqmin_tmp(2,2)/clg
  keqmin_tmp(3,2) = keqmin_tmp(3,2)/clg
      
  keqaq_tmp(1)    = keqaq_tmp(1)/clg
  keqaq_tmp(2)    = keqaq_tmp(2)/clg
  keqaq_tmp(3)    = keqaq_tmp(3)/clg
  
  keqmin_tmp(1,2) = ( keqmin_tmp(1,2)  - PressureCorrection2 )
  keqmin_tmp(2,2) = ( keqmin_tmp(2,2)  - PressureCorrection2 )
  keqmin_tmp(3,2) = ( keqmin_tmp(3,2)  - PressureCorrection2 )
  
  keqaq_tmp(1) = ( keqaq_tmp(1)  - PressureCorrection5 )   !! OH-
  keqaq_tmp(2) = ( keqaq_tmp(2)  - PressureCorrection3 )   !! HCO3-
  keqaq_tmp(3) = ( keqaq_tmp(3)  - PressureCorrection4 )   !! CO3--
  
  keqmin_tmp(1,2) = clg*keqmin_tmp(1,2)
  keqmin_tmp(2,2) = clg*keqmin_tmp(2,2)
  keqmin_tmp(3,2) = clg*keqmin_tmp(3,2)
      
  keqaq_tmp(1) = clg*keqaq_tmp(1)
  keqaq_tmp(2) = clg*keqaq_tmp(2)
  keqaq_tmp(3) = clg*keqaq_tmp(3)

  
END IF   !!! End of CalciteCreep

DO ns = 1,nsurf_sec
  ksp = msub + ngas + nspec + ns
  IF (ntemp == 1 .OR. RunIsothermal) THEN
    keqsurf_tmp(ns) = -clg*eqsurf(ns)
  ELSE
    x1 = as1(ksp,1)
    x2 = as1(ksp,2)
    x3 = as1(ksp,3)
    x4 = as1(ksp,4)
    x5 = as1(ksp,5)
    keqsurf_tmp(ns) = -clg*(x1*DLOG(temp) + x2 +  &
        x3*temp + x4/temp + x5/(temp2))
  END IF
END DO


RETURN
END SUBROUTINE keqcalc2_init
!*********************************************************************
