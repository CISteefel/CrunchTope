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

SUBROUTINE redistributeRich(nx,ny,nz)
USE crunchtype
USE params
USE medium
USE flow
USE CrunchFunctions

IMPLICIT NONE

INTEGER(I4B), INTENT(IN)                                       :: nx
INTEGER(I4B), INTENT(IN)                                       :: ny
INTEGER(I4B), INTENT(IN)                                       :: nz


!  Internal variables and arrays
INTEGER(I4B)                                                          :: jx
INTEGER(I4B)                                                          :: jy
INTEGER(I4B)                                                          :: jz

!  ****** PARAMETERS  ****************************
REAL(DP)                                                              :: delV
REAL(DP)                                                              :: rsend_zp
REAL(DP)                                                              :: rsend_zm
REAL(DP)                                                              :: rrecv_zp
REAL(DP)                                                              :: rrecv_zm
INTEGER(I4B)                                                              :: dist_path1
INTEGER(I4B)                                                              :: dist_path2
INTEGER(I4B)                                                              :: dist_path3
INTEGER(I4B)                                                              :: dist_path4
LOGICAL(LGT)                                                          :: nextto_sat

dist_path1 = 0
dist_path2 = 0
dist_path3 = 0
dist_path4 = 0

DO jz = 1,nz
  DO jy = 1,ny
    DO jx = 1,nx
        rsend_zm = 0.0d0
        rsend_zp = 0.0d0
        rrecv_zm = 0.0d0
        rrecv_zp = 0.0d0
        ! For over-saturated cells
        IF (wc(jx,jy,jz) >= wcs(jx,jy,jz)) THEN
            IF (redistribute_wc) THEN
                delV = (wc(jx,jy,jz) - wcs(jx,jy,jz)) * dzz(jx,jy,jz) * dxx(jx) * dyy(jy)
                IF (wc(jx,jy,jz)-wcs(jx,jy,jz) > wcs(jx,jy,jz)) THEN
                    delV = wcs(jx,jy,jz) * dzz(jx,jy,jz) * dxx(jx) * dyy(jy)
                    WRITE(*,*) ' WARNING : Too much water needs to be redistributed. Should consider reducing dt!'
                END IF
                CALL gradhRich(nx,ny,nz,jx,jy,jz,rsend_zm,rsend_zp,rrecv_zm,rrecv_zp)
                CALL redist_sendRich(nx,ny,nz,jx,jy,jz,delV,rsend_zm,rsend_zp)
                IF (delV > 1e-20) THEN
                    CALL redist_sendRich(nx,ny,nz,jx,jy,jz,delV,rsend_zm,rsend_zp)
                END IF
            END IF
            wc(jx,jy,jz) = wcs(jx,jy,jz)
            dist_path1 = dist_path1 + 1
        ! For unsaturated cells
        ELSE
            nextto_sat = .FALSE.
            CALL nexttoSatRich(nx,ny,nz,jx,jy,jz,nextto_sat)
            ! For isolated cells
            IF (.NOT. nextto_sat) THEN
                IF (wc(jx,jy,jz) < 0.9999*wcs(jx,jy,jz)) THEN
                    head(jx,jy,jz) = hwc(jx,jy,jz)
                END IF
                dist_path2 = dist_path2 + 1
            ! For non-isolated cells
            ELSE
                IF (redistribute_wc) THEN
                    ! Receive moisture
                    IF (wch(jx,jy,jz) > wc(jx,jy,jz)) THEN
                        delV = (wch(jx,jy,jz) - wc(jx,jy,jz)) * dzz(jx,jy,jz) * dxx(jx) * dyy(jy)
                        IF (wch(jx,jy,jz)-wc(jx,jy,jz) > wcs(jx,jy,jz)) THEN
                            delV = wcs(jx,jy,jz) * dzz(jx,jy,jz) * dxx(jx) * dyy(jy)
                            WRITE(*,*) ' WARNINGg : Too much water needs to be redistributed. Should consider reducing dt!'
                        END IF
                        CALL gradhRich(nx,ny,nz,jx,jy,jz,rsend_zm,rsend_zp,rrecv_zm,rrecv_zp)
                        CALL redist_recvRich(nx,ny,nz,jx,jy,jz,delV,rrecv_zm,rrecv_zp)
                        IF (delV > 1e-20) THEN
                            CALL redist_recvRich(nx,ny,nz,jx,jy,jz,delV,rrecv_zm,rrecv_zp)
                        END IF
                        dist_path3 = dist_path3 + 1
                    ! Send moisture
                    ELSE
                        delV = (wc(jx,jy,jz) - wch(jx,jy,jz)) * dzz(jx,jy,jz) * dxx(jx) * dyy(jy)
                        IF (wc(jx,jy,jz)-wch(jx,jy,jz) > wcs(jx,jy,jz)) THEN
                            delV = wcs(jx,jy,jz) * dzz(jx,jy,jz) * dxx(jx) * dyy(jy)
                            WRITE(*,*) ' WARNINGg : Too much water needs to be redistributed. Should consider reducing dt!'
                        END IF
                        CALL gradhRich(nx,ny,nz,jx,jy,jz,rsend_zm,rsend_zp,rrecv_zm,rrecv_zp)
                        CALL redist_sendRich(nx,ny,nz,jx,jy,jz,delV,rsend_zm,rsend_zp)
                        IF (delV > 1e-20) THEN
                            CALL redist_sendRich(nx,ny,nz,jx,jy,jz,delV,rsend_zm,rsend_zp)
                        END IF
                        dist_path4 = dist_path4 + 1
                    END IF
                END IF
                wc(jx,jy,jz) = wch(jx,jy,jz)
                room(jx,jy,jz) = (wcs(jx,jy,jz) - wc(jx,jy,jz)) * dxx(jx) * dyy(jy) * dzz(jx,jy,jz)
            END IF
        END IF


    END DO
  END DO
END DO

IF (redistribute_wc) THEN
    WRITE(*,*) '  Redistribution via 4 mechanisms: ',dist_path1,dist_path2,dist_path3,dist_path4
END IF 



RETURN
END SUBROUTINE redistributeRich
