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

SUBROUTINE redist_sendRich(nx,ny,nz,jx,jy,jz,delV,rsend_zm,rsend_zp)
USE crunchtype
USE params
USE medium
USE flow
USE CrunchFunctions

IMPLICIT NONE

INTEGER(I4B), INTENT(IN)                                       :: nx
INTEGER(I4B), INTENT(IN)                                       :: ny
INTEGER(I4B), INTENT(IN)                                       :: nz
INTEGER(I4B), INTENT(IN)                                       :: jx
INTEGER(I4B), INTENT(IN)                                       :: jy
INTEGER(I4B), INTENT(IN)                                       :: jz
REAL(DP), INTENT(INOUT)                                           :: delV
REAL(DP), INTENT(INOUT)                                          :: rsend_zm
REAL(DP), INTENT(INOUT)                                          :: rsend_zp


!  ****** PARAMETERS  ****************************

REAL(DP)                                                              :: temp
REAL(DP)                                                              :: delVzm
REAL(DP)                                                              :: delVzp
INTEGER(I4B)                                                          :: iz

delVzm = 0.0d0
delVzp = 0.0d0
! send moisture up
IF (rsend_zm > 0.0) THEN
    delVzm = delV * rsend_zm
    iz = jz
    DO WHILE (iz .NE. 1)
        iz = iz - 1
        IF (room(jx,jy,iz) > 0.0) THEN
            ! if excess moisture < available space
            IF (room(jx,jy,iz) > delVzm) THEN
                wc(jx,jy,iz) = wc(jx,jy,iz) + delVzm/(dxx(jx)*dyy(jy)*dzz(jx,jy,iz))
                room(jx,jy,iz) = room(jx,jy,iz) - delVzm
                delVzm = 0.0d0
                EXIT
            ! if excess moisture > available space
            ELSE
                wc(jx,jy,iz) = wc(jx,jy,iz) + room(jx,jy,iz)/(dxx(jx)*dyy(jy)*dzz(jx,jy,iz))
                delVzm = delVzm - room(jx,jy,iz)
                room(jx,jy,iz) = 0.0d0
            END IF
        END IF
    END DO
    ! allow excess moisture to leave domain if permeablity > 0
    IF (Kfacz(jx,jy,iz-1) > 0.0) THEN
        delVzm = 0.0d0
    END IF
END IF

! send moisture down
IF (rsend_zp > 0.0) THEN
    delVzp = delV * rsend_zp
    iz = jz
    DO WHILE (iz .NE. nz)
        iz = iz + 1
        IF (room(jx,jy,iz) > 0.0) THEN
            ! if excess moisture < available space
            IF (room(jx,jy,iz) > delVzp) THEN
                wc(jx,jy,iz) = wc(jx,jy,iz) + delVzp/(dxx(jx)*dyy(jy)*dzz(jx,jy,iz))
                room(jx,jy,iz) = room(jx,jy,iz) - delVzp
                delVzp = 0.0d0
                EXIT
            ! if excess moisture > available space
            ELSE
                wc(jx,jy,iz) = wc(jx,jy,iz) + room(jx,jy,iz)/(dxx(jx)*dyy(jy)*dzz(jx,jy,iz))
                delVzp = delVzp - room(jx,jy,iz)
                room(jx,jy,iz) = 0.0d0
            END IF
        END IF
    END DO
    ! allow excess moisture to leave domain if permeablity > 0
    IF (Kfacz(jx,jy,iz) > 0.0) THEN
        delVzp = 0.0d0
    END IF
END IF


! If excess moisture remains, reverse sending directions
IF (delVzm + delVzp > 0.0) THEN
    delV = delVzm + delVzp
    temp = rsend_zm
    rsend_zm = rsend_zp
    rsend_zp = temp
ELSE
    delV = 0.0d0
END IF


RETURN
END SUBROUTINE redist_sendRich
