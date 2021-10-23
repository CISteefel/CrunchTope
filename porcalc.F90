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


SUBROUTINE porcalc(nx,ny,nz,nkin,jpor)
USE crunchtype
USE params
USE concentration
USE mineral
USE medium

IMPLICIT NONE
!fp! auto_par_loops=0;

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                        :: nx
INTEGER(I4B), INTENT(IN)                                        :: ny
INTEGER(I4B), INTENT(IN)                                        :: nz
INTEGER(I4B), INTENT(IN)                                        :: nkin
INTEGER(I4B), INTENT(IN)                                        :: jpor

!  Internal variables and arrays

REAL(DP)                                                        :: sum
REAL(DP)                                                        :: vinit
REAL(DP)                                                        :: porfactor

INTEGER(I4B)                                                    :: jx
INTEGER(I4B)                                                    :: jy
INTEGER(I4B)                                                    :: jz
INTEGER(I4B)                                                    :: k

!  NOTE:  Initial volume fractions should be advected as well
!         Assumed immobile for the moment

DO jz = 1,nz
  DO jy = 1,ny
    DO jx = 1,nx
      
      sum = 0.0
      DO k = 1,nkin

        vinit = volinByGrid(k,jx,jy,jz)
        
        IF (LocalEquilibrium(k)) THEN                      !! Local equilibrium fantasy, so don't change the surface area
            
          IF (volfx(k,jx,jy,jz) > 0.0d0) THEN
            area(k,jx,jy,jz) = areainByGrid(k,jx,jy,jz)
          ELSE
            area(k,jx,jy,jz) = 0.0d0
          END IF
          if (mintype(k) == 0) sum = sum + volfx(k,jx,jy,jz) !! exclude biomass
          
        ELSE                                               !! Update reactive surface area

          IF (iarea(k,jinit(jx,jy,jz)) == 0) THEN                     !!  Bulk surface area
              
            IF (vinit == 0.0d0) THEN
              area(k,jx,jy,jz) = areainByGrid(k,jx,jy,jz)* (volfx(k,jx,jy,jz)/0.01)**0.6666
              sum = sum + volfx(k,jx,jy,jz)
            ELSE
              area(k,jx,jy,jz) = areainByGrid(k,jx,jy,jz)* (volfx(k,jx,jy,jz)/vinit)**0.6666
              if (mintype(k) == 0) sum = sum + volfx(k,jx,jy,jz) !! exclude biomass
            END IF
            
          ELSE                                                        !!  Specific surface area
              
            IF ( volinByGrid(k,jx,jy,jz) == 0.0d0 .AND. volfx(k,jx,jy,jz) < voltemp(k,jinit(jx,jy,jz)) ) THEN   !!  Initially a zero volume fraction, so use "threshold" volume fraction
              area(k,jx,jy,jz) = voltemp(k,jinit(jx,jy,jz))*specificByGrid(k,jx,jy,jz)*wtmin(k)/volmol(k)
            ELSE
              area(k,jx,jy,jz) = volfx(k,jx,jy,jz)*specificByGrid(k,jx,jy,jz)*wtmin(k)/volmol(k)
            END IF
            if (mintype(k) == 0) sum = sum + volfx(k,jx,jy,jz) !! exclude biomass
            
          END IF
          

        END IF

      END DO   !!!  End of mineral loop
      
    END DO
  END DO
END DO

RETURN
END SUBROUTINE porcalc
!  *******************************************************
