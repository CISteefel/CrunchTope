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

SUBROUTINE gamma_local(ncomp,nspec,jx,jy,jz)
USE crunchtype
USE params
USE concentration
USE temperature

IMPLICIT NONE
!fp! auto_par_loops=0;

!  External variables

INTEGER(I4B), INTENT(IN)                                   :: ncomp
INTEGER(I4B), INTENT(IN)                                   :: nspec
INTEGER(I4B), INTENT(IN)                                   :: jx 
INTEGER(I4B), INTENT(IN)                                   :: jy
INTEGER(I4B), INTENT(IN)                                   :: jz

!  Internal variables

REAL(DP)                                                   :: ctotal
REAL(DP)                                                   :: ah
REAL(DP)                                                   :: bh
REAL(DP)                                                   :: bdt
REAL(DP)                                                   :: sum
REAL(DP)                                                   :: sion_tmp
REAL(DP)                                                   :: tempc

REAL(DP)                                                   :: TotalMoles

REAL(DP)                                                   :: Chargesum
REAL(DP)                                                   :: aa1
REAL(DP)                                                   :: GamWaterCheck

INTEGER(I4B)                                               :: ik

CHARACTER (LEN=3)                                          :: ulabPrint

ctotal = 0.0
DO ik = 1,ncomp+nspec
  ulabPrint = ulab(ik)
  IF (ulabPrint(1:3) /= 'H2O') THEN
    ctotal = ctotal + sp10(ik,jx,jy,jz)
  END IF
END DO

tempc = t(jx,jy,jz)

IF (ntemp == 1) THEN
  ah = adh(1)
  bh = bdh(1)
  bdt = bdot(1)
ELSE
  ah = adhcoeff(1) + adhcoeff(2)*tempc  &
      + adhcoeff(3)*tempc*tempc + adhcoeff(4)*tempc*tempc*tempc  &
      + adhcoeff(5)*tempc*tempc*tempc*tempc
  bh = bdhcoeff(1) + bdhcoeff(2)*tempc  &
      + bdhcoeff(3)*tempc*tempc + bdhcoeff(4)*tempc*tempc*tempc  &
      + bdhcoeff(5)*tempc*tempc*tempc*tempc
  bdt = bdtcoeff(1) + bdtcoeff(2)*tempc  &
      + bdtcoeff(3)*tempc*tempc + bdtcoeff(4)*tempc*tempc*tempc  &
      + bdtcoeff(5)*tempc*tempc*tempc*tempc
END IF

sum = 0.0D0

DO ik = 1,ncomp+nspec
  sum = sum + sp10(ik,jx,jy,jz)*chg(ik)*chg(ik)
END DO
sion_tmp = 0.50D0*sum

DO ik = 1,ncomp+nspec
  IF (chg(ik) == 0.0D0) THEN
    
    ulabPrint = ulab(ik)
    IF (ulabPrint(1:3) /= 'H2O') THEN

      gamWaterCheck = 1.0d0 - 0.017d0*TotalMoles
!!!   Assumes molecular weight of H2O of 18.01528
      IF (gamWaterCheck < 0.0d0) THEN
        lngamma(ik,jx,jy,jz) = DLOG(1.0d0/55.50843506)
      ELSE
        lngamma(ik,jx,jy,jz) = DLOG(gamWaterCheck/55.50843506)
      END IF

    ELSE
      
      lngamma(ik,jx,jy,jz) = clg*0.10d0*sion_tmp
      
    END IF
  ELSE
    gam_local(ik) = -clg*( (ah*chg(ik)*chg(ik)*SQRT(sion_tmp))/  &
        (1.0+ acmp(ik)*bh*SQRT(sion_tmp)) + bdt*sion_tmp )
  END IF
END DO

RETURN
END SUBROUTINE gamma_local
!*****************************************************************
