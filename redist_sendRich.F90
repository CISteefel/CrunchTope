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

SUBROUTINE redist_sendRich(nx,ny,nz,jx,jy,jz,delV,rsend,dtyr)
USE crunchtype
USE params
USE medium
USE flow
USE transport
USE CrunchFunctions

IMPLICIT NONE

INTEGER(I4B), INTENT(IN)                                       :: nx
INTEGER(I4B), INTENT(IN)                                       :: ny
INTEGER(I4B), INTENT(IN)                                       :: nz
INTEGER(I4B), INTENT(IN)                                       :: jx
INTEGER(I4B), INTENT(IN)                                       :: jy
INTEGER(I4B), INTENT(IN)                                       :: jz
REAL(DP), INTENT(IN)                                           :: dtyr
REAL(DP), INTENT(INOUT)                                           :: delV
! REAL(DP), INTENT(INOUT)                                          :: rsend_zm
! REAL(DP), INTENT(INOUT)                                          :: rsend_zp
REAL(DP), DIMENSION(-3:3)                                               :: rsend
REAL(DP)                                                       :: dt


!  ****** PARAMETERS  ****************************

REAL(DP)                                                              :: temp
REAL(DP)                                                              :: delVxm
REAL(DP)                                                              :: delVxp
INTEGER(I4B)                                                          :: ix
REAL(DP)                                                              :: delVym
REAL(DP)                                                              :: delVyp
INTEGER(I4B)                                                          :: iy
REAL(DP)                                                              :: delVzm
REAL(DP)                                                              :: delVzp
INTEGER(I4B)                                                          :: iz

delVxm = 0.0d0
delVxp = 0.0d0
delVym = 0.0d0
delVyp = 0.0d0
delVzm = 0.0d0
delVzp = 0.0d0

dt = dtyr * 365 * 86400

qx = qx/secyr
qy = qy/secyr
qz = qz/secyr

! ***************************************************************
!                   Send in x direction
! ***************************************************************

! ******************    send moisture in negative x direction
IF (rsend(-1) > 0.0) THEN
    delVxm = delV * rsend(-1)
    ix = jx - 1
    DO WHILE (ix .NE. 0)
        IF (room(ix,jy,jz) > 0.0) THEN
            ! if excess moisture < available space
            IF (room(ix,jy,jz) > delVxm) THEN
                wc(ix,jy,jz) = wc(ix,jy,jz) + delVxm/(dxx(ix)*dyy(jy)*dzz(ix,jy,jz))
                qx(ix,jy,jz) = qx(ix,jy,jz) - delVxm/(dyy(jy)*dzz(ix,jy,jz))/dt
                room(ix,jy,jz) = room(ix,jy,jz) - delVxm
                delVxm = 0.0d0
                EXIT
            ! if excess moisture > available space
            ELSE
                wc(ix,jy,jz) = wc(ix,jy,jz) + room(ix,jy,jz)/(dxx(ix)*dyy(jy)*dzz(ix,jy,jz))
                qx(ix,jy,jz) = qx(ix,jy,jz) - room(ix,jy,jz)/(dyy(jy)*dzz(ix,jy,jz))/dt
                delVxm = delVxm - room(ix,jy,jz)
                room(ix,jy,jz) = 0.0d0
            END IF
        END IF
        ix = ix - 1
    END DO

    IF (delVxm > 0.0 .AND. activecellPressure(0,jy,jz) == 1) THEN
        delVxm = 0.0d0
    END IF

    IF (y_is_vertical) THEN
        ! if extra water remains, send upward
        IF (delVxm > 0.0) THEN
            iy = jy - 1
            ix = ix + 1
            DO WHILE (activecellPressure(ix,iy,jz) == 1 .AND. iy > 0)
                IF (room(ix,iy,jz) > 0.0) THEN
                    ! if excess moisture < available space
                    IF (room(ix,iy,jz) > delVxm) THEN
                        wc(ix,iy,jz) = wc(ix,iy,jz) + delVxm/(dxx(ix)*dyy(iy)*dzz(ix,iy,jz))
                        qy(ix,iy,jz) = qy(ix,iy,jz) - delVxm/(dxx(ix)*dzz(ix,iy,jz))/dt
                        room(ix,iy,jz) = room(ix,iy,jz) - delVxm
                        delVxm = 0.0d0
                        EXIT
                    ! if excess moisture > available space
                    ELSE
                        wc(ix,iy,jz) = wc(ix,iy,jz) + room(ix,iy,jz)/(dxx(ix)*dyy(iy)*dzz(ix,iy,jz))
                        qy(ix,iy,jz) = qy(ix,iy,jz) - room(ix,iy,jz)/(dxx(ix)*dzz(ix,iy,jz))/dt
                        delVxm = delVxm - room(ix,iy,jz)
                        room(ix,iy,jz) = 0.0d0
                    END IF
                END IF
                iy = iy - 1
            END DO
        END IF
        ! if extra water still exists, send downward
        IF (delVxm > 0.0) THEN
            DO WHILE (iy .NE. ny)
                iy = iy + 1
                IF (room(ix,iy,jz) > 0.0) THEN
                    ! if excess moisture < available space
                    IF (room(ix,iy,jz) > delVxm) THEN
                        wc(ix,iy,jz) = wc(ix,iy,jz) + delVxm/(dxx(ix)*dyy(iy)*dzz(ix,iy,jz))
                        qy(ix,iy-1,jz) = qy(ix,iy-1,jz) + delVxm/(dxx(ix)*dzz(ix,iy,jz))/dt
                        room(ix,iy,jz) = room(ix,iy,jz) - delVxm
                        delVxm = 0.0d0
                        EXIT
                    ! if excess moisture > available space
                    ELSE
                        wc(ix,iy,jz) = wc(ix,iy,jz) + room(ix,iy,jz)/(dxx(ix)*dyy(iy)*dzz(ix,iy,jz))
                        qy(ix,iy-1,jz) = qy(ix,iy-1,jz) + room(ix,iy,jz)/(dxx(ix)*dzz(ix,iy,jz))/dt
                        delVxm = delVxm - room(ix,iy,jz)
                        room(ix,iy,jz) = 0.0d0
                    END IF
                END IF
            END DO
        END IF
    END IF

END IF



! ******************    send moisture in positive x direction
IF (rsend(1) > 0.0) THEN
    delVxp = delV * rsend(1)
    ix = jx + 1
    DO WHILE (ix .NE. nx+1)
        IF (room(ix,jy,jz) > 0.0) THEN
            ! if excess moisture < available space
            IF (room(ix,jy,jz) > delVxp) THEN
                wc(ix,jy,jz) = wc(ix,jy,jz) + delVxp/(dxx(ix)*dyy(jy)*dzz(ix,jy,jz))
                qx(ix-1,jy,jz) = qx(ix-1,jy,jz) + delVxp/(dyy(jy)*dzz(ix,jy,jz))/dt
                room(ix,jy,jz) = room(ix,jy,jz) - delVxp
                delVxp = 0.0d0
                EXIT
            ! if excess moisture > available space
            ELSE
                wc(ix,jy,jz) = wc(ix,jy,jz) + room(ix,jy,jz)/(dxx(ix)*dyy(jy)*dzz(ix,jy,jz))
                qx(ix-1,jy,jz) = qx(ix-1,jy,jz) + room(ix,jy,jz)/(dyy(jy)*dzz(ix,jy,jz))/dt
                delVxp = delVxp - room(ix,jy,jz)
                room(ix,jy,jz) = 0.0d0
            END IF
        END IF
        ix = ix + 1
    END DO

    IF (delVxp > 0.0 .AND. activecellPressure(nx+1,jy,jz) == 1) THEN
        delVxp = 0.0d0
    END IF

    IF (y_is_vertical) THEN
        ! if extra water remains, send upward
        IF (delVxp > 0.0) THEN
            iy = jy - 1
            ix = ix - 1
            DO WHILE (activecellPressure(ix,iy,jz) == 1 .AND. iy > 0)
                IF (room(ix,iy,jz) > 0.0) THEN
                    ! if excess moisture < available space
                    IF (room(ix,iy,jz) > delVxp) THEN
                        wc(ix,iy,jz) = wc(ix,iy,jz) + delVxp/(dxx(ix)*dyy(iy)*dzz(ix,iy,jz))
                        qy(ix,iy,jz) = qy(ix,iy,jz) - delVxp/(dxx(ix)*dzz(ix,iy,jz))/dt
                        room(ix,iy,jz) = room(ix,iy,jz) - delVxp
                        delVxp = 0.0d0
                        EXIT
                    ! if excess moisture > available space
                    ELSE
                        wc(ix,iy,jz) = wc(ix,iy,jz) + room(ix,iy,jz)/(dxx(ix)*dyy(iy)*dzz(ix,iy,jz))
                        qy(ix,iy,jz) = qy(ix,iy,jz) - room(ix,iy,jz)/(dxx(ix)*dzz(ix,iy,jz))/dt
                        delVxp = delVxp - room(ix,iy,jz)
                        room(ix,iy,jz) = 0.0d0
                    END IF
                END IF
                iy = iy - 1
            END DO
        END IF
        ! if extra water still exists, send downward
        IF (delVxp > 0.0) THEN
            DO WHILE (iy .NE. ny)
                iy = iy + 1
                IF (room(ix,iy,jz) > 0.0) THEN
                    ! if excess moisture < available space
                    IF (room(ix,iy,jz) > delVxp) THEN
                        wc(ix,iy,jz) = wc(ix,iy,jz) + delVxp/(dxx(ix)*dyy(iy)*dzz(ix,iy,jz))
                        qy(ix,iy-1,jz) = qy(ix,iy-1,jz) + delVxp/(dxx(ix)*dzz(ix,iy,jz))/dt
                        room(ix,iy,jz) = room(ix,iy,jz) - delVxp
                        delVxp = 0.0d0
                        EXIT
                    ! if excess moisture > available space
                    ELSE
                        wc(ix,iy,jz) = wc(ix,iy,jz) + room(ix,iy,jz)/(dxx(ix)*dyy(iy)*dzz(ix,iy,jz))
                        qy(ix,iy-1,jz) = qy(ix,iy-1,jz) + room(ix,iy,jz)/(dxx(ix)*dzz(ix,iy,jz))/dt
                        delVxp = delVxp - room(ix,iy,jz)
                        room(ix,iy,jz) = 0.0d0
                    END IF
                END IF
            END DO
        END IF
    END IF

END IF

! ***************************************************************
!                   Send in y direction
! ***************************************************************

! 2D x-y
IF (y_is_vertical) THEN
    ! send moisture up
    IF (rsend(-2) > 0.0) THEN
        delVym = delV * rsend(-2)
        iy = jy - 1
        DO WHILE (activecellPressure(jx,iy,jz) == 1 .AND. iy > 0)
            IF (room(jx,iy,jz) > 0.0) THEN
                ! if excess moisture < available space
                IF (room(jx,iy,jz) > delVym) THEN
                    wc(jx,iy,jz) = wc(jx,iy,jz) + delVym/(dxx(jx)*dyy(iy)*dzz(jx,iy,jz))
                    qy(jx,iy,jz) = qy(jx,iy,jz) - delVym/(dxx(jx)*dzz(jx,iy,jz))/dt
                    room(jx,iy,jz) = room(jx,iy,jz) - delVym
                    delVym = 0.0d0
                    EXIT
                ! if excess moisture > available space
                ELSE
                    wc(jx,iy,jz) = wc(jx,iy,jz) + room(jx,iy,jz)/(dxx(jx)*dyy(iy)*dzz(jx,iy,jz))
                    qy(jx,iy,jz) = qy(jx,iy,jz) - room(jx,iy,jz)/(dxx(jx)*dzz(jx,iy,jz))/dt
                    delVym = delVym - room(jx,iy,jz)
                    room(jx,iy,jz) = 0.0d0
                END IF
            END IF
            iy = iy - 1
        END DO
        ! allow excess moisture to leave domain if permeablity > 0
        ! IF (Kfacy(jx,iy,jz) > 0.0) THEN
        !     delVym = 0.0d0
        ! END IF
    END IF

    ! send moisture down
    IF (rsend(2) > 0.0) THEN
        delVyp = delV * rsend(2)
        iy = jy
        DO WHILE (iy .NE. ny)
            iy = iy + 1
            IF (room(jx,iy,jz) > 0.0) THEN
                ! if excess moisture < available space
                IF (room(jx,iy,jz) > delVyp) THEN
                    wc(jx,iy,jz) = wc(jx,iy,jz) + delVyp/(dxx(jx)*dyy(iy)*dzz(jx,iy,jz))
                    qy(jx,iy-1,jz) = qy(jx,iy-1,jz) + delVyp/(dxx(jx)*dzz(jx,iy,jz))/dt
                    room(jx,iy,jz) = room(jx,iy,jz) - delVyp
                    delVyp = 0.0d0
                    EXIT
                ! if excess moisture > available space
                ELSE
                    wc(jx,iy,jz) = wc(jx,iy,jz) + room(jx,iy,jz)/(dxx(jx)*dyy(iy)*dzz(jx,iy,jz))
                    qy(jx,iy-1,jz) = qy(jx,iy-1,jz) + room(jx,iy,jz)/(dxx(jx)*dzz(jx,iy,jz))/dt
                    delVyp = delVyp - room(jx,iy,jz)
                    room(jx,iy,jz) = 0.0d0
                END IF
            END IF
        END DO
        ! allow excess moisture to leave domain if permeablity > 0
        ! IF (Kfacz(jx,iy,jz) > 0.0) THEN
        !     delVyp = 0.0d0
        ! END IF
    END IF

ELSE
    ! 2D x-z
    ! send moisture up
    IF (rsend(-3) > 0.0) THEN
        delVzm = delV * rsend(-3)
        iz = jz - 1
        DO WHILE (activecellPressure(jx,jy,iz) == 1 .AND. iz > 0)
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
            iz = iz - 1
        END DO
        ! allow excess moisture to leave domain if permeablity > 0
        IF (Kfacz(jx,jy,iz) > 0.0) THEN
            delVzm = 0.0d0
        END IF
    END IF

    ! send moisture down
    IF (rsend(3) > 0.0) THEN
        delVzp = delV * rsend(3)
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
END IF


! If excess moisture remains, reverse sending directions
IF (y_is_vertical) THEN
    IF (delVym + delVyp > 0.0) THEN
        temp = rsend(-2)
        rsend(-2) = rsend(2)
        rsend(2) = temp
    END IF
    IF (delVxm + delVxp > 0.0) THEN
        temp = rsend(-1)
        rsend(-1) = rsend(1)
        rsend(1) = temp
    END IF

ELSE
    IF (delVzm + delVzp > 0.0) THEN
        delV = delVzm + delVzp
        temp = rsend(-3)
        rsend(-3) = rsend(3)
        rsend(3) = temp
    ELSE
        delV = 0.0d0
    END IF
END IF




IF (y_is_vertical) THEN
    delV = delVxm + delVxp + delVym + delVyp
ELSE
    delV = delVxm + delVxp + delVzm + delVzp
END IF

qx = qx*secyr
qy = qy*secyr
qz = qz*secyr


RETURN
END SUBROUTINE redist_sendRich
