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
    
SUBROUTINE gammaUpdated(ncomp,nspec,nsurf,nexchange,npot,jx,jy,jz,igamma)
USE crunchtype
USE params
USE concentration
USE temperature
USE runtime, ONLY: Benchmark
USE mineral

IMPLICIT NONE

!  External variables

INTEGER(I4B), INTENT(IN)                                   :: ncomp
INTEGER(I4B), INTENT(IN)                                   :: nspec
INTEGER(I4B), INTENT(IN)                                   :: nsurf
INTEGER(I4B), INTENT(IN)                                   :: nexchange
INTEGER(I4B), INTENT(IN)                                   :: npot
INTEGER(I4B), INTENT(IN)                                   :: jx 
INTEGER(I4B), INTENT(IN)                                   :: jy
INTEGER(I4B), INTENT(IN)                                   :: jz
INTEGER(I4B), INTENT(IN)                                   :: igamma

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

!!!INTEGER(I4B)                                               :: i
INTEGER(I4B)                                               :: ik
INTEGER(I4B)                                               :: it
INTEGER(I4B)                                               :: ItPoint
INTEGER(I4B)                                               :: pos_der

REAL(DP)                                                   :: tempc
REAL(DP)                                                   :: sqrt_IonS
REAL(DP)                                                   :: tmp1
REAL(DP)                                                   :: tmp2
REAL(DP)                                                   :: tmp3
REAL(DP)                                                   :: gammawaterTMP
REAL(DP)                                                   :: bdotpar

CHARACTER (LEN=3)                                          :: ulabPrint

LOGICAL(LGT)                                              :: Davies
LOGICAL(LGT)                                              :: Wateq_Extended_DH
LOGICAL(LGT)                                              :: Helgeson
LOGICAL(LGT)                                              :: Unity
  
Davies = .FALSE.
Helgeson = .FALSE.
Unity = .FALSE.
Wateq_Extended_DH = .FALSE.

gammawater(jx,jy,jz) = EXP(lngammawater(jx,jy,jz))
gammawaterTMP = gammawater(jx,jy,jz)
tempc = t(jx,jy,jz)
sion_tmp = sion(jx,jy,jz)
sqrt_IonS = SQRT(sion_tmp)

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

!*** Now calculate the activity coefficients for : water, other uncharged species, charged species

gamma(:,jx,jy,jz) = 1D0
lngamma(:,jx,jy,jz)  = 0D0
deriv_gamma(:,:,jx,jy,jz) = 0D0 

IF (igamma == 0) THEN ! unit activity coefficient except for water

  DO ik = 1,ncomp+nspec
    ulabPrint = ulab(ik)
    IF (ulabPrint(1:3) == 'H2O' .or. ulabPrint(1:3) == 'HHO') THEN
      gamma(ik,jx,jy,jz) = gammawaterTMP
      lngamma(ik,jx,jy,jz) = LOG(gammawaterTMP)
      pos_der = ncomp + nexchange + nsurf + npot + 1 + 1       ! gamma water
      deriv_gamma(ik,pos_der,jx,jy,jz) = gammawaterTMP                  ! gammawater is a primary variable  
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
      gamma(ik,jx,jy,jz) = gammawaterTMP
      lngamma(ik,jx,jy,jz) = LOG(gammawaterTMP)
      pos_der = ncomp + nexchange + nsurf + npot + 1 + 1      ! gamma water--note that activity of water is 1 + 1 
      deriv_gamma(ik,pos_der,jx,jy,jz) = gammawaterTMP            ! gammawater is a primary variable  
  
    ELSE IF (chg(ik) == 0D0) THEN ! neutral species
      
      lngamma(ik,jx,jy,jz) = clg * 0.1D0 * sion_tmp
      gamma(ik,jx,jy,jz) = EXP(lngamma(ik,jx,jy,jz))
      pos_der = ncomp + nexchange + nsurf + npot + 1          ! ionic strength
      deriv_gamma(ik,pos_der,jx,jy,jz) = gamma(ik,jx,jy,jz) * lngamma(ik,jx,jy,jz) 
  
    ELSE ! charged species
  
      IF (Davies) THEN
        aa1 = - clg * ah * chg(ik) * chg(ik)
        lngamma(ik,jx,jy,jz) = aa1 * (sqrt_IonS/(1.0d0 + sqrt_IonS)  - 0.3d0 * sion_tmp  )
        gamma(ik,jx,jy,jz) = EXP(lngamma(ik,jx,jy,jz))
  
        !*** Only one non-zero derivative = relative to ionic strength
        pos_der = ncomp + nexchange + nsurf + npot + 1
        deriv_gamma(ik,pos_der,jx,jy,jz) = aa1 * gamma(ik,jx,jy,jz) * sion_tmp *         &
                  ( ( 0.5D0 / (sion_tmp + sqrt_IonS) - 0.5D0 / ((1 + sqrt_IonS)*(1 + sqrt_IonS)) ) - 0.3D0 )
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
        
  
        lngamma(ik,jx,jy,jz) = clg * ( - (tmp2 * sqrt_IonS) / tmp1 + bdotpar * sion_tmp)
        gamma(ik,jx,jy,jz) = EXP(lngamma(ik,jx,jy,jz))
  
        !*** Only one non-zero derivative = relative to ionic strength
        pos_der = ncomp + nexchange + nsurf + npot + 1
        deriv_gamma(ik,pos_der,jx,jy,jz) = clg * gamma(ik,jx,jy,jz) *    &
                      ( (- 0.5D0 * tmp2 * sqrt_IonS) /tmp1 + tmp2 * tmp3 / (tmp1 * tmp1) + bdotpar * sion_tmp ) 
      END IF
  
    END IF
  END DO  
END IF

RETURN
END SUBROUTINE gammaUpdated
!*****************************************************************
