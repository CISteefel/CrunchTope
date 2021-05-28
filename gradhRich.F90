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

SUBROUTINE gradhRich(jx,jy,jz,rsend,rrecv)
USE crunchtype
USE params
USE medium
USE flow
USE CrunchFunctions

IMPLICIT NONE

INTEGER(I4B), INTENT(IN)                                       :: jx
INTEGER(I4B), INTENT(IN)                                       :: jy
INTEGER(I4B), INTENT(IN)                                       :: jz
! REAL(DP), INTENT(INOUT)                                          :: rsend_zm
! REAL(DP), INTENT(INOUT)                                          :: rsend_zp
! REAL(DP), INTENT(INOUT)                                          :: rrecv_zm
! REAL(DP), INTENT(INOUT)                                          :: rrecv_zp
REAL(DP), INTENT(INOUT), DIMENSION(-3:3)                         :: rsend
REAL(DP), INTENT(INOUT), DIMENSION(-3:3)                         :: rrecv

!  ****** PARAMETERS  ****************************
REAL(DP)                                       :: dhxm
REAL(DP)                                       :: dhxp
REAL(DP)                                       :: dhym
REAL(DP)                                       :: dhyp
REAL(DP)                                       :: dhzm
REAL(DP)                                       :: dhzp
REAL(DP)                                       :: dh_tot

! Calculate head gradient in each direction

! x direction
IF (Kfacx(jx,jy,jz) > 0.0) THEN
    dhxp = 2.0d0*(head(jx+1,jy,jz) - head(jx,jy,jz)) / (dxx(jx+1) + dxx(jx))
ELSE
    dhxp = 0.0d0
END IF
IF (Kfacx(jx-1,jy,jz) > 0.0) THEN
    dhxm = 2.0d0*(head(jx,jy,jz) - head(jx-1,jy,jz)) / (dxx(jx) + dxx(jx-1))
ELSE
    dhxm = 0.0d0
END IF

! ! y direction
IF (Kfacy(jx,jy,jz) > 0.0) THEN
    dhyp = 2.0d0*(head(jx,jy+1,jz) - head(jx,jy,jz)) / (dyy(jy+1) + dyy(jy))
ELSE
    dhyp = 0.0d0
END IF
IF (Kfacy(jx,jy-1,jz) > 0.0) THEN
    dhym = 2.0d0*(head(jx,jy,jz) - head(jx,jy-1,jz)) / (dyy(jy) + dyy(jy-1))
ELSE
    dhym = 0.0d0
END IF

! z direction
IF (Kfacz(jx,jy,jz) > 0.0) THEN
    dhzp = 2.0d0*(head(jx,jy,jz+1) - head(jx,jy,jz)) / (dzz(jx,jy,jz+1) + dzz(jx,jy,jz)) - 1.0d0
ELSE
    dhzp = 0.0d0
END IF
IF (Kfacz(jx,jy,jz-1) > 0.0) THEN
    dhzm = 2.0d0*(head(jx,jy,jz) - head(jx,jy,jz-1)) / (dzz(jx,jy,jz) + dzz(jx,jy,jz-1)) - 1.0d0
ELSE
    dhzm = 0.0d0
END IF

! ! ignore horizontal redistribution for now, Zhi Li 20200711
! dhxp = 0.0d0
! dhxm = 0.0d0
! dhyp = 0.0d0
! dhym = 0.0d0

! calculate split ratio
dh_tot = abs(dhxp) + abs(dhxm) + abs(dhyp) + abs(dhym) + abs(dhzp) + abs(dhzm)

! split ratio in x
IF (dhxp*dhxm >= 0.0) THEN
    IF (dhxp > 0.0 .OR. dhxm > 0.0) THEN
        rsend(-1) = (abs(dhxp) + abs(dhxm)) / dh_tot
        rrecv(1) = (abs(dhxp) + abs(dhxm)) / dh_tot
    ELSE
        rrecv(-1) = (abs(dhxp) + abs(dhxm)) / dh_tot
        rsend(1) = (abs(dhxp) + abs(dhxm)) / dh_tot
    END IF
ELSE
    IF (dhxp > 0.0) THEN
        rrecv(-1) = abs(dhxm) / dh_tot
        rrecv(1) = abs(dhxp) / dh_tot
    ELSE
        rsend(-1) = abs(dhxm) / dh_tot
        rsend(1) = abs(dhxp) / dh_tot
    END IF
END IF

! split ratio in y
IF (dhyp*dhym >= 0.0) THEN
    IF (dhyp > 0.0 .OR. dhym > 0.0) THEN
        rsend(-2) = (abs(dhyp) + abs(dhym)) / dh_tot
        rrecv(2) = (abs(dhyp) + abs(dhym)) / dh_tot
    ELSE
        rrecv(-2) = (abs(dhyp) + abs(dhym)) / dh_tot
        rsend(2) = (abs(dhyp) + abs(dhym)) / dh_tot
    END IF
ELSE
    IF (dhyp > 0.0) THEN
        rrecv(-2) = abs(dhym) / dh_tot
        rrecv(2) = abs(dhyp) / dh_tot
    ELSE
        rsend(-2) = abs(dhym) / dh_tot
        rsend(2) = abs(dhyp) / dh_tot
    END IF
END IF

! split ratio in z
IF (dhzp*dhzm >= 0.0) THEN
    IF (dhzp > 0.0 .OR. dhzm > 0.0) THEN
        rsend(-3) = (abs(dhzp) + abs(dhzm)) / dh_tot
        rrecv(3) = (abs(dhzp) + abs(dhzm)) / dh_tot
    ELSE
        rrecv(-3) = (abs(dhzp) + abs(dhzm)) / dh_tot
        rsend(3) = (abs(dhzp) + abs(dhzm)) / dh_tot
    END IF
ELSE
    IF (dhzp > 0.0) THEN
        rrecv(-3) = abs(dhzm) / dh_tot
        rrecv(3) = abs(dhzp) / dh_tot
    ELSE
        rsend(-3) = abs(dhzm) / dh_tot
        rsend(3) = abs(dhzp) / dh_tot
    END IF
END IF


RETURN
END SUBROUTINE gradhRich
