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
    
SUBROUTINE gamma_init(ncomp,nspec,tempc,sqrt_sion)
USE crunchtype
USE params
USE runtime, ONLY: Benchmark
USE concentration
USE temperature

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

ChargeSum = 0.0d0
TotalMoles = 0.0d0
DO ik = 1,ncomp+nspec
  IF (ulab(ik) /= 'H2O') THEN
    TotalMoles = TotalMoles + sptmp10(ik)
    ChargeSum = ChargeSum + sptmp10(ik)*chg(ik)*chg(ik)
  ELSE
    CONTINUE
  END IF
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

tconv = 273.15d0
tempk = tempc + tconv
call dhconst(ah,bh,tempk,tconv)

DO ik = 1,ncomp+nspec

  IF (chg(ik) == 0.0d0) THEN

    IF (ulab(ik) == 'H2O') THEN

      gamWaterCheck = 1.0d0 - 0.017d0*TotalMoles
      gamtmp(ik) = DLOG(gamWaterCheck)
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

RETURN
END SUBROUTINE gamma_init
!*****************************************************************
