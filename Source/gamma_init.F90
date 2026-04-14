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
    
SUBROUTINE gamma_init(ncomp,nspec,igamma,tempc,sqrt_IonS,sion_tmp)
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

INTEGER(I4B), INTENT(IN)                                   :: igamma

REAL(DP), INTENT(IN)                                       :: tempc
REAL(DP), INTENT(OUT)                                      :: sqrt_IonS
REAL(DP), INTENT(OUT)                                      :: sion_tmp

!  Internal variables

REAL(DP)                                                   :: TotalMoles
REAL(DP)                                                   :: ah
REAL(DP)                                                   :: bh
REAL(DP)                                                   :: bdt
REAL(DP)                                                   :: Chargesum

REAL(DP)                                                   :: aa1
REAL(DP)                                                   :: GamWaterCheck

REAL(DP)                                                   :: dhad,dhbd,tempk,tconv

INTEGER(I4B)                                               :: ik
INTEGER(I4B)                                               :: it
INTEGER(I4B)                                               :: ItPoint

REAL(DP)                                                   :: gammawaterTMP
REAL(DP)                                                   :: bdotpar

CHARACTER (LEN=3)                                          :: ulabPrint

REAL(DP)                                                   :: tmp1
REAL(DP)                                                   :: tmp2
REAL(DP)                                                   :: tmp3

LOGICAL(LGT)                                              :: Davies
LOGICAL(LGT)                                              :: Wateq_Extended_DH
LOGICAL(LGT)                                              :: Helgeson
LOGICAL(LGT)                                              :: Unity

Davies = .FALSE.
Helgeson = .FALSE.
Unity = .FALSE.
Wateq_Extended_DH = .FALSE.

ChargeSum = 0.0d0
DO ik = 1,ncomp+nspec
  ulabPrint = ulab(ik)
  IF (ulabPrint(1:3) /= 'H2O' .AND. ulabPrint(1:3) /= 'HHO') THEN
    ChargeSum = ChargeSum + sptmp10(ik)*chg(ik)*chg(ik)
  END IF
END DO
sion_tmp = 0.50D0*ChargeSum
sqrt_IonS = SQRT(sion_tmp)

TotalMoles = 0.0d0
DO ik = 1,ncomp+nspec
  ulabPrint = ulab(ik)
  IF (ulabPrint(1:3) /= 'H2O' .and. ulabPrint(1:3) /= 'HHO') THEN
    TotalMoles = TotalMoles + sptmp10(ik)
  END IF
END DO
gammawaterTmp = 1.0d0 - 0.017d0*TotalMoles
!!!IF (gammawaterTmp < 0.01) THEN
!!!  gammawaterTmp = 0.01
!!!END IF

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

!*** Now calculate the activity coefficients for: water, other uncharged species, charged species

gamtmp  = 0D0

IF (igamma == 0) THEN ! unit activity coefficient except for water

  DO ik = 1,ncomp+nspec
    ulabPrint = ulab(ik)
    IF (ulabPrint(1:3) == 'H2O' .or. ulabPrint(1:3) == 'HHO') THEN
      
      gamtmp(ik) = LOG(gammawaterTMP)

    END IF
    
  END DO

ELSE

  DO ik = 1,ncomp+nspec
  
    IF (IncludeBdot) THEN
      IF (azero(ik) == 0.0d0) THEN
        Davies = .TRUE.
      ELSE
        Wateq_Extended_DH = .TRUE.
      END IF
    ELSE
      Helgeson = .TRUE. !!  Helgesonian-LLNL bdot expression based on extended Debye-Huckel
    END IF
    
    ulabPrint = ulab(ik)   
    IF (ulabPrint(1:3) == 'H2O' .or. ulabPrint(1:3) == 'HHO') THEN

      gamtmp(ik) = LOG(gammawaterTMP)
      CONTINUE

  
    ELSE IF (chg(ik) == 0D0) THEN ! neutral species
      
      gamtmp(ik) = clg * 0.1D0 * sion_tmp

  
    ELSE ! charged species
  
      IF (Davies) THEN
        aa1 = - clg * ah * chg(ik) * chg(ik)
        gamtmp(ik) = aa1 * (sqrt_IonS/(1.0d0 + sqrt_IonS)  - 0.3d0 * sion_tmp  )

      END IF 
  
      IF (Helgeson) then
        bdotpar = bdt
      ELSE IF (Wateq_Extended_DH) THEN
        bdotpar = bdotparameter(ik)
      ElSE
        bdotpar = 0.0d0
      END IF
  
      IF (Helgeson .OR. Wateq_Extended_DH) THEN
  
        tmp1 = 1D0 + bh * azero(IK) * sqrt_IonS 
        tmp2 = ah * chg(ik) * chg(ik)
        tmp3 = 0.5D0 * bh * azero(IK) * sion_tmp

        gamtmp(ik) = clg * ( - (tmp2 * sqrt_IonS) / tmp1 + bdotpar * sion_tmp)

      END IF
  
    END IF
  END DO  
END IF

IF (SaltCreep) THEN
  
  aa1 = -(ah*chg(ikNa)*chg(ikNa)*sqrt_IonS)/                &
              (1.0d0 + 4.08*bh*sqrt_IonS)                   &         
              + 0.082*sion_tmp
  gamtmp(ikNa) = clg*aa1
  
  aa1 = -(ah*chg(ikCl)*chg(ikCl)*sqrt_IonS)/                &
              (1.0d0 + 3.63*bh*sqrt_IonS)                   &         
              + 0.017*sion_tmp
  gamtmp(ikCl) = clg*aa1
  
END IF

RETURN
END SUBROUTINE gamma_init
!*****************************************************************
