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
    
SUBROUTINE GammaResidual(ncomp,nspec,tempc,jx,jy,jz)
USE crunchtype
USE params
USE concentration, ONLY: muaq,s,sn,sp,sp10,chg,adh,bdh,bdot,adhcoeff,bdhcoeff,bdtcoeff,acmp,  &
     NumericalJac,sppTMP,sppTMP10,gamTMP,keqaqTMP,stmp
USE solver, only: fjac

IMPLICIT NONE

!  External variables

INTEGER(I4B), INTENT(IN)                                  :: ncomp
INTEGER(I4B), INTENT(IN)                                  :: nspec
REAL(DP), INTENT(IN)                                      :: tempc
INTEGER(I4B), INTENT(IN)                                  :: jx
INTEGER(I4B), INTENT(IN)                                  :: jy
INTEGER(I4B), INTENT(IN)                                  :: jz


!  Internal variables

REAL(DP)                                                   :: TotalMoles
REAL(DP)                                                   :: ah
REAL(DP)                                                   :: bh
REAL(DP)                                                   :: bdt
REAL(DP)                                                   :: Chargesum
REAL(DP)                                                   :: sionTMP
REAL(DP)                                                   :: aa1
REAL(DP)                                                   :: GamWaterCheck

REAL(DP)                                                   :: dhad,dhbd,tempk,tconv

INTEGER(I4B)                                               :: ik
INTEGER(I4B)                                               :: kk
INTEGER(I4B)                                               :: it
INTEGER(I4B)                                               :: ItPoint
INTEGER(I4B)                                               :: ksp
INTEGER(I4B)                                               :: i

REAL(DP)                                                   :: sqrt_sion
REAL(DP)                                                   :: sum

CHARACTER (LEN=3)                                          :: ulabPrint

REAL(DP)                                                   :: perturb=1.0E-09

sTMP = 0.0d0

ah = adhcoeff(1) + adhcoeff(2)*tempc  &
       + adhcoeff(3)*tempc*tempc + adhcoeff(4)*tempc*tempc*tempc  &
       + adhcoeff(5)*tempc*tempc*tempc*tempc
bh = bdhcoeff(1) + bdhcoeff(2)*tempc  &
       + bdhcoeff(3)*tempc*tempc + bdhcoeff(4)*tempc*tempc*tempc  &
       + bdhcoeff(5)*tempc*tempc*tempc*tempc
bdt = bdtcoeff(1) + bdtcoeff(2)*tempc  &
       + bdtcoeff(3)*tempc*tempc + bdtcoeff(4)*tempc*tempc*tempc  &
       + bdtcoeff(5)*tempc*tempc*tempc*tempc

NumericalJac = 0.0d0
do ik = 1,ncomp

   sppTMP(ik) = sppTMP(ik) + perturb
   sppTMP10(ik) = DEXP(sppTMP(ik))
   
   ChargeSum = 0.0d0
   DO i = 1,ncomp+nspec    
     ChargeSum  = ChargeSum + sppTMP10(i)*chg(i)*chg(i)  
   END DO
   sionTMP = 0.50D0*ChargeSum
   
   if (sionTMP > 5.0d0) THEN
     sionTMP = 1.0d0
   end if
   
   sqrt_sion = DSQRT(sionTMP)
   
   DO i = 1,ncomp+nspec
     IF (chg(i) == 0.0D0) THEN  
       gamTMP(i) = clg*0.10d0*sionTMP    
     ELSE    
       aa1 = -(ah*chg(I)*chg(I)*sqrt_sion)/                  &
             (1.0d0 + acmp(I)*bh*sqrt_sion)                  &         
             + bdt*sionTMP
       gamTMP(i) = clg*aa1
     END IF
   END DO
   
   DO ksp = 1,nspec   
      sum = 0.0D0
      DO i = 1,ncomp
        sum = sum + muaq(ksp,i)*(sppTMP(i) + gamTMP(i))
      END DO
      sppTMP(ksp+ncomp) = keqaqTMP(ksp) - gamTMP(ksp+ncomp) + sum
      sppTMP10(ksp+ncomp) = DEXP(sppTMP(ksp+ncomp))           
   END DO
      
    DO i = 1,ncomp
      sum=0.0D0
      DO ksp = 1,nspec
        kk = ksp + ncomp
          sum = sum + muaq(ksp,i)*sppTMP10(kk)
      END DO
      sTMP(i) = sum + sppTMP10(i)
    END DO
       
    do i = 1,ncomp
      NumericalJac(ik,i) = ( sTMP(i) - s(i,jx,jy,jz) )/perturb
    end do
        
    sppTMP(ik) = sppTMP(ik) - perturb
    sppTMP10(ik) = DEXP(sppTMP(ik))
    
END DO

RETURN

END SUBROUTINE GammaResidual
!*****************************************************************
