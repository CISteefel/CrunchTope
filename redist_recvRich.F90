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

SUBROUTINE redist_recvRich(nx,ny,nz,jx,jy,jz,delV,rrecv)
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
! REAL(DP), INTENT(INOUT)                                          :: rrecv_zm
! REAL(DP), INTENT(INOUT)                                          :: rrecv_zp
REAL(DP), DIMENSION(-3:3)                                               :: rrecv

!  ****** PARAMETERS  ****************************
INTEGER(I4B)                                                          :: ix
INTEGER(I4B)                                                          :: iy
INTEGER(I4B)                                                          :: iz
REAL(DP)                                                              :: temp
REAL(DP)                                                              :: delVxm
REAL(DP)                                                              :: delVxp
REAL(DP)                                                              :: delVym
REAL(DP)                                                              :: delVyp
REAL(DP)                                                              :: delVzm
REAL(DP)                                                              :: delVzp

delVxm = 0.0d0
delVxp = 0.0d0
delVym = 0.0d0
delVyp = 0.0d0
delVzm = 0.0d0
delVzp = 0.0d0


! extract moisture in x direction
IF (rrecv(-1) > 0.0) THEN
    delVxm = delV * rrecv(1)
    ix = jx
    DO WHILE (ix .NE. 1)
        ix = ix - 1
        IF (wc(ix,jy,jz) > wcr) THEN
            IF ((wc(ix,jy,jz)-wcr) > delVxm/(dxx(ix)*dyy(jy)*dzz(ix,jy,jz))) THEN
                wc(ix,jy,jz) = wc(ix,jy,jz) - delVxm/(dxx(ix)*dyy(jy)*dzz(ix,jy,jz))
                room(ix,jy,jz) = room(ix,jy,jz) + delVxm
                delVxm = 0.0d0
                EXIT
            ELSE
                delVxm = delVxm - (wc(ix,jy,jz)-wcr) * (dxx(ix)*dyy(jy)*dzz(ix,jy,jz))
                room(ix,jy,jz) = room(ix,jy,jz) + (wc(ix,jy,jz)-wcr) * (dxx(ix)*dyy(jy)*dzz(ix,jy,jz))
                wc(ix,jy,jz) = wcr
            END IF
        END IF
        ! limit recv within 1 cell
        EXIT
    END DO
END IF

IF (rrecv(1) > 0.0) THEN
    delVxp = delV * rrecv(1)
    ix = jx
    DO WHILE (ix .NE. nx-1)
        ix = ix + 1
        IF (wc(ix,jy,jz) > wcr) THEN
            IF ((wc(ix,jy,jz)-wcr) > delVxp/(dxx(ix)*dyy(jy)*dzz(ix,jy,jz))) THEN
                wc(ix,jy,jz) = wc(ix,jy,jz) - delVxp/(dxx(ix)*dyy(jy)*dzz(ix,jy,jz))
                room(ix,jy,jz) = room(ix,jy,jz) + delVxp
                delVxp = 0.0d0
                EXIT
            ELSE
                delVxp = delVxp - (wc(ix,jy,jz)-wcr) * (dxx(ix)*dyy(jy)*dzz(ix,jy,jz))
                room(ix,jy,jz) = room(ix,jy,jz) + (wc(ix,jy,jz)-wcr) * (dxx(ix)*dyy(jy)*dzz(ix,jy,jz))
                wc(ix,jy,jz) = wcr
            END IF
        END IF
        ! limit recv within 1 cell
        EXIT
    END DO
END IF

! 2D x-y
IF (y_is_vertical) THEN
    ! extract moisture from up
    IF (rrecv(-2) > 0.0) THEN
        delVym = delV * rrecv(-2)
        iy = jy
        DO WHILE (activecellPressure(jx,iy,jz) == 1)
            iy = iy - 1
            IF (wc(jx,iy,jz) > wcr) THEN
                IF ((wc(jx,iy,jz)-wcr) > delVym/(dxx(jx)*dyy(iy)*dzz(jx,iy,jz))) THEN
                    IF (activecellPressure(jx,iy,jz) == 1) THEN
                        wc(jx,iy,jz) = wc(jx,iy,jz) - delVym/(dxx(jx)*dyy(iy)*dzz(jx,iy,jz))
                    END IF
                    room(jx,iy,jz) = room(jx,iy,jz) + delVym
                    delVym = 0.0d0
                    EXIT
                ELSE
                    delVym = delVym - (wc(jx,iy,jz)-wcr) * (dxx(jx)*dyy(iy)*dzz(jx,iy,jz))
                    room(jx,iy,jz) = room(jx,iy,jz) + (wc(jx,iy,jz)-wcr) * (dxx(jx)*dyy(iy)*dzz(jx,iy,jz))
                    IF (activecellPressure(jx,jy,jz) == 1) THEN
                        wc(jx,iy,jz) = wcr
                    END IF
                END IF
            END IF
            ! limit recv within 1 cell
            EXIT
        END DO
    END IF

    ! extract moisture from down
    IF (rrecv(2) > 0.0) THEN
        delVyp = delV * rrecv(2)
        iy = jy
        DO WHILE (iy .NE. ny-1)
            iy = iy + 1
            IF (wc(jx,iy,jz) > wcr) THEN
                IF ((wc(jx,iy,jz)-wcr) > delVyp/(dxx(jx)*dyy(iy)*dzz(jx,iy,jz))) THEN
                    wc(jx,iy,jz) = wc(jx,iy,jz) - delVyp/(dxx(jx)*dyy(iy)*dzz(jx,iy,jz))
                    room(jx,iy,jz) = room(jx,iy,jz) + delVyp
                    delVyp = 0.0d0
                    EXIT
                ELSE
                    delVyp = delVyp - (wc(jx,iy,jz)-wcr) * (dxx(jx)*dyy(iy)*dzz(jx,iy,jz))
                    room(jx,iy,jz) = room(jx,iy,jz) + (wc(jx,iy,jz)-wcr) * (dxx(jx)*dyy(iy)*dzz(jx,iy,jz))
                    wc(jx,iy,jz) = wcr
                END IF
            END IF
            ! limit recv within 1 cell
            EXIT
        END DO
    END IF


ELSE
    ! 2D x-z
    ! extract moisture from up
    IF (rrecv(-3) > 0.0) THEN
        delVzm = delV * rrecv(-3)
        iz = jz
        DO WHILE (iz .NE. 1)
            iz = iz - 1
            IF (wc(jx,jy,iz) > wcr) THEN
                IF ((wc(jx,jy,iz)-wcr) > delVzm/(dxx(jx)*dyy(jy)*dzz(jx,jy,iz))) THEN
                    IF (activecellPressure(jx,jy,iz) == 1) THEN
                        wc(jx,jy,iz) = wc(jx,jy,iz) - delVzm/(dxx(jx)*dyy(jy)*dzz(jx,jy,iz))
                    END IF
                    room(jx,jy,iz) = room(jx,jy,iz) + delVzm
                    delVzm = 0.0d0
                    EXIT
                ELSE
                    delVzm = delVzm - (wc(jx,jy,iz)-wcr) * (dxx(jx)*dyy(jy)*dzz(jx,jy,iz))
                    room(jx,jy,iz) = room(jx,jy,iz) + (wc(jx,jy,iz)-wcr) * (dxx(jx)*dyy(jy)*dzz(jx,jy,iz))
                    IF (activecellPressure(jx,jy,jz) == 1) THEN
                        wc(jx,jy,iz) = wcr
                    END IF
                END IF
            END IF
            ! limit recv within 1 cell
            EXIT
        END DO
        ! IF (Kfacz(jx,jy,iz-1) > 0.0) THEN
        !     delVzm = 0.0d0
        ! END IF
    END IF

    ! extract moisture from down
    IF (rrecv(3) > 0.0) THEN
        delVzp = delV * rrecv(3)
        iz = jz
        DO WHILE (iz .NE. nz-1)
            iz = iz + 1
            IF (wc(jx,jy,iz) > wcr) THEN
                IF ((wc(jx,jy,iz)-wcr) > delVzp/(dxx(jx)*dyy(jy)*dzz(jx,jy,iz))) THEN
                    wc(jx,jy,iz) = wc(jx,jy,iz) - delVzp/(dxx(jx)*dyy(jy)*dzz(jx,jy,iz))
                    room(jx,jy,iz) = room(jx,jy,iz) + delVzp
                    delVzp = 0.0d0
                    EXIT
                ELSE
                    delVzp = delVzp - (wc(jx,jy,iz)-wcr) * (dxx(jx)*dyy(jy)*dzz(jx,jy,iz))
                    room(jx,jy,iz) = room(jx,jy,iz) + (wc(jx,jy,iz)-wcr) * (dxx(jx)*dyy(jy)*dzz(jx,jy,iz))
                    wc(jx,jy,iz) = wcr
                END IF
            END IF
            ! limit recv within 1 cell
            EXIT
        END DO
        ! IF (Kfacz(jx,jy,iz) > 0.0) THEN
        !     delVzp = 0.0d0
        ! END IF
    END IF
END IF


! If excess moisture remains, reverse extracting directions
! IF (delVzm + delVzp > 0.0) THEN
!     delV = delVzm + delVzp
!     temp = rrecv(-3)
!     rrecv(-3) = rrecv(3)
!     rrecv(3) = temp
! ELSE
!     delV = 0.0d0
! END IF

IF (y_is_vertical) THEN
    delV = delVxm + delVxp + delVym + delVyp
ELSE
    delV = delVxm + delVxp + delVzm + delVzp
END IF


RETURN
END SUBROUTINE redist_recvRich
