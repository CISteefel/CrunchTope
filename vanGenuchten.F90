!!! *** Copyright Notice ***
!!! ìCrunchFlowî, Copyright (c) 2016, The Regents of the University of California, through Lawrence Berkeley National Laboratory
!!! (subject to receipt of any required approvals from the U.S. Dept. of Energy).† All rights reserved.
!!!†
!!! If you have questions about your rights to use or distribute this software, please contact
!!! Berkeley Lab's Innovation & Partnerships Office at††IPO@lbl.gov.
!!!†
!!! NOTICE.† This Software was developed under funding from the U.S. Department of Energy and the U.S. Government
!!! consequently retains certain rights. As such, the U.S. Government has been granted for itself and others acting
!!! on its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the Software to reproduce, distribute copies to the public,
!!! prepare derivative works, and perform publicly and display publicly, and to permit other to do so.
!!!
!!! *** License Agreement ***
!!! ìCrunchFlowî, Copyright (c) 2016, The Regents of the University of California, through Lawrence Berkeley National Laboratory)
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

!!! >>>>> Get the water retention parameters following van Genuchten model, Zhi Li20200629


SUBROUTINE vanGenuchten(nx,ny,nz)

USE crunchtype
USE params
USE solver
USE medium
USE transport
USE temperature
USE flow
USE CrunchFunctions

REAL(DP)                                                              :: m
REAL(DP)                                                              :: satu
REAL(DP), PARAMETER                                                   :: visc=0.001d0
REAL(DP), PARAMETER                                                   :: Ss=1.0D-05
REAL(DP), PARAMETER                                                   :: grav=9.8d0
REAL(DP), PARAMETER                                                   :: rho=1.0d3


DO jz = 0,nz+1
    DO jy = 0,ny+1
        DO jx = 0,nx+1
            m = 1.0d0 - 1.0d0/vgn(jx,jy,jz)
            h = head(jx,jy,jz)
            !!! Get saturation from head
            satu = (1 + abs(vga(jx,jy,jz)*h)**vgn(jx,jy,jz)) ** (-m)
            ! satu = (wc(jx,jy,jz)- wcr(jx,jy,jz)) / (wcs(jx,jy,jz)- wcr(jx,jy,jz))
            IF (satu > 1) THEN
                satu = 1.0d0
            ELSE IF (satu < 0) THEN
                satu = 0.0d0
            END IF

            !!! Get wc from head
            IF (h > 0.0d0) THEN
                wch(jx,jy,jz) = wcs(jx,jy,jz)
            ELSE
                wch(jx,jy,jz) = wcr(jx,jy,jz) + (wcs(jx,jy,jz) - wcr(jx,jy,jz)) * satu
            END IF
            !!! Get K from head
            IF (h > 0.0d0) THEN
                Kr(jx,jy,jz) = 1.0d0
            ELSE
                Kr(jx,jy,jz) = satu**(0.5) * (1.0d0 - (1.0d0 - satu**(1.0/m))**m)**2.0
                ! top boundary
                IF (activecellPressure(jx,jy,jz) == 0 .AND. jy <= ny) THEN
                    IF (activecellPressure(jx,jy+1,jz) == 1) THEN
                        Kr(jx,jy,jz) = 1.0d0
                    END IF
                END IF
                IF (Kr(jx,jy,jz) > 1.0d0) THEN
                    Kr(jx,jy,jz) = 1.0d0
                END IF
            END IF
            ! IF (activecellPressure(jx,jy,jz) == 1) THEN
            !     IF (activecellPressure(jx,jy-1,jz) == 0) THEN
            !         WRITE(*,*) 'head, wc, satu, m, Kr, a, n = ',h,wc(jx,jy,jz),satu,m,Kr(jx,jy,jz),vga(jx,jy,jz),vgn(jx,jy,jz)
            !     END IF
            ! END IF
            !!! Get C from head
            Ch(jx,jy,jz) = (vga(jx,jy,jz)*m*vgn(jx,jy,jz)*(wcs(jx,jy,jz)-wcr(jx,jy,jz))*abs(vga(jx,jy,jz)*h)**(vgn(jx,jy,jz)-1)) / (1.0d0 + abs(vga(jx,jy,jz)*h)**vgn(jx,jy,jz))**(m+1)
            IF (h > 0.0d0) THEN
                Ch(jx,jy,jz) = 0.0d0
            END IF
            !!! get h from wc
            IF (wc(jx,jy,jz) - wcr(jx,jy,jz) < 0.001) THEN
                wc(jx,jy,jz) = wcr(jx,jy,jz) + 0.001
            END IF
            IF (wc(jx,jy,jz) > wcs(jx,jy,jz)) THEN
                hwc(jx,jy,jz) = 0.0d0
            ELSE
                hwc(jx,jy,jz) = -(1.0/vga(jx,jy,jz)) * (((wcs(jx,jy,jz) - wcr(jx,jy,jz))/(wc(jx,jy,jz) - wcr(jx,jy,jz)))**(1.0/m) - 1.0) ** (1.0/vgn(jx,jy,jz))
            END IF
        END DO
    END DO
END DO

!!! Get face conductivity
DO jz = 1,nz
    DO jy = 1,ny
        DO jx = 0,nx
            IF (nx <= 3) THEN
                Kfacx(jx,jy,jz) = 0.0d0
            ELSE
                IF (upstream_weighting) THEN
                    IF (qx(jx,jy,jz) >= 0.0) THEN
                        Kfacx(jx,jy,jz) = permx(jx,jy,jz)*Kr(jx,jy,jz) * (ro(jx,jy,jz)*grav/visc)
                    ELSE
                        Kfacx(jx,jy,jz) = permx(jx+1,jy,jz)*Kr(jx+1,jy,jz) * (ro(jx+1,jy,jz)*grav/visc)
                    END IF
                ELSE
                    Kfacx(jx,jy,jz) = 0.5 * (permx(jx,jy,jz)*Kr(jx,jy,jz) + permx(jx+1,jy,jz)*Kr(jx+1,jy,jz)) * (ro(jx,jy,jz)*grav/visc)
                END IF
                ! zero conductivity if either cell is impermeable
                IF (permx(jx,jy,jz) == 0.0d0 .OR. permx(jx+1,jy,jz) == 0.0d0) THEN
                    Kfacx(jx,jy,jz) = 0.0d0
                END IF
                ! zero conductivity for inactive cells
                IF (activecellPressure(jx,jy,jz) == 0) THEN
                    Kfacx(jx,jy,jz) = 0.0d0
                    IF (jx > 0) THEN
                        Kfacx(jx-1,jy,jz) = 0.0d0
                    END IF
                END IF
            END IF
        END DO
    END DO
END DO
DO jz = 1,nz
    DO jy = 0,ny
        DO jx = 1,nx
            IF (ny == 1) THEN
                Kfacy(jx,jy,jz) = 0.0d0
            ELSE
                IF (upstream_weighting) THEN
                    IF (qy(jx,jy,jz) >= 0.0) THEN
                        Kfacy(jx,jy,jz) = permy(jx,jy,jz)*Kr(jx,jy,jz) * (ro(jx,jy,jz)*grav/visc)
                    ELSE
                        Kfacy(jx,jy,jz) = permy(jx,jy+1,jz)*Kr(jx,jy+1,jz) * (ro(jx,jy+1,jz)*grav/visc)
                    END IF
                ELSE
                    Kfacy(jx,jy,jz) = 0.5 * (permy(jx,jy,jz)*Kr(jx,jy,jz) + permy(jx,jy+1,jz)*Kr(jx,jy+1,jz)) * (ro(jx,jy,jz)*grav/visc)
                END IF
                IF (permy(jx,jy,jz) == 0.0d0 .OR. permy(jx,jy+1,jz) == 0.0d0) THEN
                    Kfacy(jx,jy,jz) = 0.0d0
                END IF
            END IF
            ! if 2D x-y
            IF (y_is_vertical) THEN
                IF (activecellPressure(jx,jy,jz) == 0) THEN
                    IF (activecellPressure(jx,jy+1,jz) == 1) THEN
                        ! saturated if pressure is prescribed
                        IF (headOld(jx,jy,jz) > 0.0) THEN
                            Kfacy(jx,jy,jz) = permy(jx,jy,jz) * (ro(jx,jy,jz)*grav/visc)
                        ELSE
                            Kfacy(jx,jy,jz) = 0.5 * (permy(jx,jy,jz)*Kr(jx,jy,jz) + permy(jx,jy+1,jz)*Kr(jx,jy+1,jz)) * (ro(jx,jy,jz)*grav/visc)
                            IF (upstream_weighting) THEN
                                Kfacy(jx,jy,jz) = permy(jx,jy,jz)*Kr(jx,jy,jz) * (ro(jx,jy,jz)*grav/visc)
                            END IF
                        END IF
                    ELSE
                        Kfacy(jx,jy,jz) = 0.0d0
                    END IF
                END IF
            END IF
        END DO
    END DO
END DO
DO jz = 0,nz
    DO jy = 1,ny
        DO jx = 1,nx
            IF (nz == 1) THEN
                Kfacz(jx,jy,jz) = 0.0d0
            ELSE
                Kfacz(jx,jy,jz) = 0.5 * (permz(jx,jy,jz)*Kr(jx,jy,jz) + permz(jx,jy,jz+1)*Kr(jx,jy,jz+1)) * (ro(jx,jy,jz)*grav/visc)
                IF (permz(jx,jy,jz) == 0.0d0 .OR. permz(jx,jy,jz+1) == 0.0d0) THEN
                    Kfacz(jx,jy,jz) = 0.0d0
                END IF
                IF (jz == 0 .AND. pres(jx,jy,jz) > 0) THEN
                    Kfacz(jx,jy,jz) = permz(jx,jy,jz) * (ro(jx,jy,jz)*grav/visc)
                END IF
                IF (activecellPressure(jx,jy,jz) == 0) THEN
                    IF (activecellPressure(jx,jy,jz+1) == 1) THEN
                        Kfacz(jx,jy,jz) = 0.5 * (permz(jx,jy,jz)*Kr(jx,jy,jz) + permz(jx,jy,jz+1)*Kr(jx,jy,jz+1)) * (ro(jx,jy,jz)*grav/visc)
                    ELSE
                        Kfacz(jx,jy,jz) = 0.0d0
                    END IF
                END IF
                ! impermeable bottom
                IF (jz == nz) THEN
                    Kfacz(jx,jy,jz) = 0.0d0
                END IF
            END IF
        END DO
    END DO
END DO




RETURN
END SUBROUTINE vanGenuchten
