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

!!! >>>>> Solves the head-form of the Richards equation for variably-saturated subsurface flow

SUBROUTINE pressureRich (nx,ny,nz,dtyr,amatP,SteadyFlow)
USE crunchtype
USE params
USE solver
USE medium
USE transport
USE temperature
USE flow
USE CrunchFunctions

#include "petsc/finclude/petscmat.h"
USE petscmat

IMPLICIT NONE

INTEGER(I4B), INTENT(IN)                                       :: nx
INTEGER(I4B), INTENT(IN)                                       :: ny
INTEGER(I4B), INTENT(IN)                                       :: nz

REAL(DP)                                                       :: dt
REAL(DP), INTENT(IN)                                           :: dtyr

LOGICAL(LGT), INTENT(IN)                                        :: SteadyFlow

!  Internal variables and arrays

INTEGER(I4B)                                                          :: jx
INTEGER(I4B)                                                          :: jy
INTEGER(I4B)                                                          :: jz
INTEGER(I4B)                                                          :: j
INTEGER(I4B)                                                          :: i
INTEGER(I4B)                                                          :: ierr
INTEGER(I4B)                                                          :: itsiterate
INTEGER(I4B)                                                          :: nxyz
INTEGER(I4B)                                                          :: rows
INTEGER(I4B)                                                          :: cols

REAL(DP)                                                              :: AccumulationTerm
REAL(DP)                                                              :: RightHandSide
REAL(DP)                                                              :: DiagonalTerm
REAL(DP)                                                              :: ct
REAL(DP)                                                              :: Check
REAL(DP)                                                              :: AddPressureX
REAL(DP)                                                              :: AddPressureY
REAL(DP)                                                              :: AddPressureZ

!  ****** PARAMETERS  ****************************

REAL(DP), PARAMETER                                                   :: visc=0.000001d0
!REAL(DP), PARAMETER                                                   :: visc=31.54d0
REAL(DP), PARAMETER                                                   :: ctTransient=4.8D-07
REAL(DP), PARAMETER                                                   :: ctSteady=0.00
REAL(DP), PARAMETER                                                   :: big=100000.0D0
REAL(DP), PARAMETER                                                   :: zero=0.0d0
REAL(DP), PARAMETER                                                   :: grav=9.8d0
REAL(DP), PARAMETER                                                   :: alphaUnsat=0.001d0
REAL(DP), PARAMETER                                                   :: alphaSat=0.000d0
REAL(DP), PARAMETER                                                   :: Ss=1.0d-5

REAL(DP), DIMENSION(-3:3)                                             :: coef

REAL(DP)                                                              :: alphaBear

CHARACTER (LEN=1)                                                     :: Coordinate

! *******************begin PETSc declarations of f90 variables***********
INTEGER(I4B)             ::numprocs
INTEGER(I4B)             ::irank
INTEGER(I4B)             ::linefil
INTEGER(I4B), PARAMETER  ::maxitsksp=1000

!*********************end PETSc declarations ******************************

! ******************** PETSC declarations ********************************
PetscFortranAddr     userP(6)
Mat                  amatP
! ************************end PETSc declarations of PETSc variables ******

IF (SteadyFlow) THEN
  ct = ctSteady
ELSE
  ct = ctTransient
END IF

ct = ctTransient

alphaBear = alphaSat

CALL MatZeroEntries(amatP,ierr)

dt = dtyr * 365 * 86400

!! *********************** Calculate head from pressure ***************************

DO jz = 0,nz+1
  DO jy = 0,ny+1
    DO jx = 0,nx+1
      head(jx,jy,jz) = pres(jx,jy,jz) / (ro(jx,jy,jz) * grav)
    END DO
  END DO
END DO


!*********************begin Matrix Coefficients ******************************
!! --- Matrix coefficients for the interior of domain
DO jz = 2,nz-1
    DO jy = 2,ny-1
        DO jx = 2,nx-1
            j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1
            IF (activecellPressure(jx,jy,jz) == 0) THEN
                CALL MatSetValues(amatP,1,j,1,j-1,zero,INSERT_VALUES,ierr)
                CALL MatSetValues(amatP,1,j,1,j+1,zero,INSERT_VALUES,ierr)
                CALL MatSetValues(amatP,1,j,1,j-nx,zero,INSERT_VALUES,ierr)
                CALL MatSetValues(amatP,1,j,1,j+nx,zero,INSERT_VALUES,ierr)
                CALL MatSetValues(amatP,1,j,1,j-nx*ny,zero,INSERT_VALUES,ierr)
                CALL MatSetValues(amatP,1,j,1,j+nx*ny,zero,INSERT_VALUES,ierr)
                CALL MatSetValues(amatP,1,j,1,j,big,INSERT_VALUES,ierr)
                BvecCrunchP(j) = head(jx,jy,jz)*big
                XvecCrunchP(j) = head(jx,jy,jz)
            ELSE
                COORDINATE = 'X'
                coef(-1) = -2.0d0*dt*Kfacx(jx-1,jy,jz) / (dxx(jx)*(dxx(jx)+dxx(jx-1)))
                coef(1) = -2.0d0*dt*Kfacx(jx,jy,jz) / (dxx(jx)*(dxx(jx+1)+dxx(jx)))

                COORDINATE = 'Y'
                coef(-2) = -2.0d0*dt*Kfacy(jx,jy-1,jz) / (dyy(jy)*(dyy(jy)+dyy(jy-1)))
                coef(2)  = -2.0d0*dt*Kfacy(jx,jy,jz) / (dyy(jy)*(dyy(jy+1)+dyy(jy) ))

                COORDINATE = 'Z'
                coef(-3) = -2.0d0*dt*Kfacz(jx,jy,jz-1) / (dzz(jx,jy,jz)*(dzz(jx,jy,jz)+dzz(jx,jy,jz-1)))
                coef(3)  = -2.0d0*dt*Kfacz(jx,jy,jz) / (dzz(jx,jy,jz)*(dzz(jx,jy,jz+1)+dzz(jx,jy,jz) ))

                coef(0) = (Ch(jx,jy,jz)+Ss*wc(jx,jy,jz)/wcs(jx,jy,jz)) - coef(-3) - coef(-2) - coef(-1) - coef(1) - coef(2) - coef(3)

                CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
                CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
                CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
                CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
                CALL MatSetValues(amatP,1,j,1,j-nx*ny,coef(-3),INSERT_VALUES,ierr)
                CALL MatSetValues(amatP,1,j,1,j+nx*ny,coef(3),INSERT_VALUES,ierr)
                CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)

                BvecCrunchP(j) = (Ch(jx,jy,jz)+Ss*wc(jx,jy,jz)/wcs(jx,jy,jz))*headOld(jx,jy,jz) - (dt/dzz(jx,jy,jz))*(Kfacz(jx,jy,jz)-Kfacz(jx,jy,jz-1))
                XvecCrunchP(j) = head(jx,jy,jz)

                ! IF (jx == 2 .AND. jy == 2 .AND. jz == 2) THEN
                !     WRITE(*,*) '>>> MAT: jz, coef, rhs = ',jz,coef,BvecCrunchP(j)
                !     WRITE(*,*) '>>> K, dt, dz, C = ',Kfacz(jx,jy,jz)*1e6,dt,dzz(jx,jy,jz),Ch(jx,jy,jz)
                ! END IF
            END IF
            continue
        END DO
    END DO
END DO


!!! --- Matrix coefficients for the faces

! xm face
jx = 1
DO jz = 2,nz-1
  DO jy = 2,ny-1
    j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1
    IF (activecellPressure(jx,jy,jz) == 0) THEN
!!            CALL MatSetValues(amatP,1,j,1,j-1,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+1,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx*ny,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx*ny,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j,big,INSERT_VALUES,ierr)
        BvecCrunchP(j) = head(jx,jy,jz)*big
        XvecCrunchP(j) = head(jx,jy,jz)
    ELSE

        COORDINATE = 'X'
        coef(1) = -2.0d0*dt*Kfacx(jx,jy,jz) / (dxx(jx)*(dxx(jx+1)+dxx(jx)))
        IF (activecellPressure(0,jy,jz) == 0) THEN
            coef(-1) = -2.0d0*dt*Kfacx(jx-1,jy,jz) / (dxx(jx)*(dxx(jx)+dxx(jx-1)))
        ELSE
            coef(-1) = 0.0d0
        END IF
        COORDINATE = 'Y'
        coef(-2) = -2.0d0*dt*Kfacy(jx,jy-1,jz) / (dyy(jy)*(dyy(jy)+dyy(jy-1)))
        coef(2)  = -2.0d0*dt*Kfacy(jx,jy,jz) / (dyy(jy)*(dyy(jy+1)+dyy(jy) ))

        COORDINATE = 'Z'
        coef(-3) = -2.0d0*dt*Kfacz(jx,jy,jz-1) / (dzz(jx,jy,jz)*(dzz(jx,jy,jz)+dzz(jx,jy,jz-1)))
        coef(3)  = -2.0d0*dt*Kfacz(jx,jy,jz) / (dzz(jx,jy,jz)*(dzz(jx,jy,jz+1)+dzz(jx,jy,jz) ))

        coef(0) = (Ch(jx,jy,jz)+Ss*wc(jx,jy,jz)/wcs(jx,jy,jz)) - coef(-3) - coef(-2) - coef(-1) - coef(1) - coef(2) - coef(3)

    !!            CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx*ny,coef(-3),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx*ny,coef(3),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)

        IF (activecellPressure(0,jy,jz) == 0) THEN
            AddPressureX =  -coef(-1) * headOld(0,jy,jz)
        ELSE
            AddPressureX = 0.0d0
        END IF
        BvecCrunchP(j) = (Ch(jx,jy,jz)+Ss*wc(jx,jy,jz)/wcs(jx,jy,jz))*headOld(jx,jy,jz) - (dt/dzz(jx,jy,jz))*(Kfacz(jx,jy,jz)-Kfacz(jx,jy,jz-1)) + AddPressureX
        XvecCrunchP(j) = head(jx,jy,jz)
    END IF
  END DO
END DO

! xp face
jx = nx
DO jz = 2,nz-1
  DO jy = 2,ny-1
    j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1
    IF (activecellPressure(jx,jy,jz) == 0) THEN
        CALL MatSetValues(amatP,1,j,1,j-1,zero,INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j+1,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx*ny,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx*ny,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j,big,INSERT_VALUES,ierr)
        BvecCrunchP(j) = head(jx,jy,jz)*big
        XvecCrunchP(j) = head(jx,jy,jz)
    ELSE

        COORDINATE = 'X'
        coef(-1) = -2.0d0*dt*Kfacx(jx-1,jy,jz) / (dxx(jx)*(dxx(jx)+dxx(jx-1)))
        IF (activecellPressure(nx+1,jy,jz) == 0) THEN
            coef(1) = -2.0d0*dt*Kfacx(jx,jy,jz) / (dxx(jx)*(dxx(jx+1)+dxx(jx)))
        ELSE
            coef(1) = 0.0d0
        END IF
        COORDINATE = 'Y'
        coef(-2) = -2.0d0*dt*Kfacy(jx,jy-1,jz) / (dyy(jy)*(dyy(jy)+dyy(jy-1)))
        coef(2)  = -2.0d0*dt*Kfacy(jx,jy,jz) / (dyy(jy)*(dyy(jy+1)+dyy(jy) ))

        COORDINATE = 'Z'
        coef(-3) = -2.0d0*dt*Kfacz(jx,jy,jz-1) / (dzz(jx,jy,jz)*(dzz(jx,jy,jz)+dzz(jx,jy,jz-1)))
        coef(3)  = -2.0d0*dt*Kfacz(jx,jy,jz) / (dzz(jx,jy,jz)*(dzz(jx,jy,jz+1)+dzz(jx,jy,jz) ))

        coef(0) = (Ch(jx,jy,jz)+Ss*wc(jx,jy,jz)/wcs(jx,jy,jz)) - coef(-3) - coef(-2) - coef(-1) - coef(1) - coef(2) - coef(3)

        CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx*ny,coef(-3),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx*ny,coef(3),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)

        IF (activecellPressure(nx+1,jy,jz) == 0) THEN
            AddPressureX =  -coef(1) * headOld(nx+1,jy,jz)
        ELSE
            AddPressureX = 0.0d0
        END IF
        BvecCrunchP(j) = (Ch(jx,jy,jz)+Ss*wc(jx,jy,jz)/wcs(jx,jy,jz))*headOld(jx,jy,jz) - (dt/dzz(jx,jy,jz))*(Kfacz(jx,jy,jz)-Kfacz(jx,jy,jz-1)) + AddPressureX
        XvecCrunchP(j) = head(jx,jy,jz)
    END IF
  END DO
END DO

! ym face
jy = 1
DO jz = 2,nz-1
  DO jx = 2,nx-1
    j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1
    IF (activecellPressure(jx,jy,jz) == 0) THEN
        CALL MatSetValues(amatP,1,j,1,j-1,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+1,zero,INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j-nx,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx*ny,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx*ny,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j,big,INSERT_VALUES,ierr)
        BvecCrunchP(j) = head(jx,jy,jz)*big
        XvecCrunchP(j) = head(jx,jy,jz)
    ELSE

        COORDINATE = 'X'
        coef(-1) = -2.0d0*dt*Kfacx(jx-1,jy,jz) / (dxx(jx)*(dxx(jx)+dxx(jx-1)))
        coef(1) = -2.0d0*dt*Kfacx(jx,jy,jz) / (dxx(jx)*(dxx(jx+1)+dxx(jx)))

        COORDINATE = 'Y'
        coef(2)  = -2.0d0*dt*Kfacy(jx,jy,jz) / (dyy(jy)*(dyy(jy+1)+dyy(jy) ))
        IF (activecellPressure(jx,0,jz) == 0) THEN
            coef(-2) = -2.0d0*dt*Kfacy(jx,jy-1,jz) / (dyy(jy)*(dyy(jy)+dyy(jy-1)))
        ELSE
            coef(-2) = 0.0d0
        END IF

        COORDINATE = 'Z'
        coef(-3) = -2.0d0*dt*Kfacz(jx,jy,jz-1) / (dzz(jx,jy,jz)*(dzz(jx,jy,jz)+dzz(jx,jy,jz-1)))
        coef(3)  = -2.0d0*dt*Kfacz(jx,jy,jz) / (dzz(jx,jy,jz)*(dzz(jx,jy,jz+1)+dzz(jx,jy,jz) ))

        coef(0) = (Ch(jx,jy,jz)+Ss*wc(jx,jy,jz)/wcs(jx,jy,jz)) - coef(-3) - coef(-2) - coef(-1) - coef(1) - coef(2) - coef(3)

        CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx*ny,coef(-3),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx*ny,coef(3),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)

        IF (activecellPressure(jx,0,jz) == 0) THEN
            AddPressureY =  -coef(-2) * headOld(jx,0,jz)
        ELSE
            AddPressureY = 0.0d0
        END IF

        BvecCrunchP(j) = (Ch(jx,jy,jz)+Ss*wc(jx,jy,jz)/wcs(jx,jy,jz))*headOld(jx,jy,jz) - (dt/dzz(jx,jy,jz))*(Kfacz(jx,jy,jz)-Kfacz(jx,jy,jz-1)) + AddPressureY
        XvecCrunchP(j) = head(jx,jy,jz)
    END IF
  END DO
END DO

! yp face
jy = ny
DO jz = 2,nz-1
  DO jx = 2,nx-1
    j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1
    IF (activecellPressure(jx,jy,jz) == 0) THEN
        CALL MatSetValues(amatP,1,j,1,j-1,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+1,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx,zero,INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j+nx,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx*ny,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx*ny,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j,big,INSERT_VALUES,ierr)
        BvecCrunchP(j) = head(jx,jy,jz)*big
        XvecCrunchP(j) = head(jx,jy,jz)
    ELSE

        COORDINATE = 'X'
        coef(-1) = -2.0d0*dt*Kfacx(jx-1,jy,jz) / (dxx(jx)*(dxx(jx)+dxx(jx-1)))
        coef(1) = -2.0d0*dt*Kfacx(jx,jy,jz) / (dxx(jx)*(dxx(jx+1)+dxx(jx)))

        COORDINATE = 'Y'
        coef(-2) = -2.0d0*dt*Kfacy(jx,jy-1,jz) / (dyy(jy)*(dyy(jy)+dyy(jy-1)))
        IF (activecellPressure(jx,ny+1,jz) == 0) THEN
            coef(2)  = -2.0d0*dt*Kfacy(jx,jy,jz) / (dyy(jy)*(dyy(jy+1)+dyy(jy) ))
        ELSE
            coef(2) = 0.0d0
        END IF

        COORDINATE = 'Z'
        coef(-3) = -2.0d0*dt*Kfacz(jx,jy,jz-1) / (dzz(jx,jy,jz)*(dzz(jx,jy,jz)+dzz(jx,jy,jz-1)))
        coef(3)  = -2.0d0*dt*Kfacz(jx,jy,jz) / (dzz(jx,jy,jz)*(dzz(jx,jy,jz+1)+dzz(jx,jy,jz) ))

        coef(0) = (Ch(jx,jy,jz)+Ss*wc(jx,jy,jz)/wcs(jx,jy,jz)) - coef(-3) - coef(-2) - coef(-1) - coef(1) - coef(2) - coef(3)

        CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx*ny,coef(-3),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx*ny,coef(3),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)

        IF (activecellPressure(jx,ny+1,jz) == 0) THEN
            AddPressureY =  -coef(2) * headOld(jx,ny+1,jz)
        ELSE
            AddPressureY = 0.0d0
        END IF

        BvecCrunchP(j) = (Ch(jx,jy,jz)+Ss*wc(jx,jy,jz)/wcs(jx,jy,jz))*headOld(jx,jy,jz) - (dt/dzz(jx,jy,jz))*(Kfacz(jx,jy,jz)-Kfacz(jx,jy,jz-1)) + AddPressureY
        XvecCrunchP(j) = head(jx,jy,jz)

    END IF
  END DO
END DO

! zm face
jz = 1
DO jy = 2,ny-1
  DO jx = 2,nx-1
    j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1
    IF (activecellPressure(jx,jy,jz) == 0) THEN
        CALL MatSetValues(amatP,1,j,1,j-1,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+1,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx,zero,INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j-nx*ny,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx*ny,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j,big,INSERT_VALUES,ierr)
        BvecCrunchP(j) = head(jx,jy,jz)*big
        XvecCrunchP(j) = head(jx,jy,jz)
    ELSE

        COORDINATE = 'X'
        coef(-1) = -2.0d0*dt*Kfacx(jx-1,jy,jz) / (dxx(jx)*(dxx(jx)+dxx(jx-1)))
        coef(1) = -2.0d0*dt*Kfacx(jx,jy,jz) / (dxx(jx)*(dxx(jx+1)+dxx(jx)))

        COORDINATE = 'Y'
        coef(2)  = -2.0d0*dt*Kfacy(jx,jy,jz) / (dyy(jy)*(dyy(jy+1)+dyy(jy) ))
        coef(-2) = -2.0d0*dt*Kfacy(jx,jy-1,jz) / (dyy(jy)*(dyy(jy)+dyy(jy-1)))

        COORDINATE = 'Z'
        coef(3)  = -2.0d0*dt*Kfacz(jx,jy,jz) / (dzz(jx,jy,jz)*(dzz(jx,jy,jz+1)+dzz(jx,jy,jz) ))
        IF (activecellPressure(jx,jy,0) == 0) THEN
            coef(-3) = -2.0d0*dt*Kfacz(jx,jy,jz-1) / (dzz(jx,jy,jz)*(dzz(jx,jy,jz)+dzz(jx,jy,jz-1)))
        ELSE
            coef(-3) = 0.0d0
        END IF

        coef(0) = (Ch(jx,jy,jz)+Ss*wc(jx,jy,jz)/wcs(jx,jy,jz)) - coef(-3) - coef(-2) - coef(-1) - coef(1) - coef(2) - coef(3)

        CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j-nx*ny,coef(-3),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx*ny,coef(3),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)

        IF (activecellPressure(jx,jy,0) == 0) THEN
            AddPressureZ =  -coef(-3) * headOld(jx,jy,0)
        ELSE
            AddPressureZ = 0.0d0
        END IF

        BvecCrunchP(j) = (Ch(jx,jy,jz)+Ss*wc(jx,jy,jz)/wcs(jx,jy,jz))*headOld(jx,jy,jz) - (dt/dzz(jx,jy,jz))*(Kfacz(jx,jy,jz)-Kfacz(jx,jy,jz-1)) + AddPressureZ
        XvecCrunchP(j) = head(jx,jy,jz)

        ! IF (jx==2 .AND. jy==26) THEN
        !     WRITE(*,*) ' ---------- '
        !     WRITE(*,*) '    Ch, hOld, Km, Kp, BC = ',Ch(jx,jy,jz),headOld(jx,jy,jz),Kfacz(jx,jy,jz-1),Kfacz(jx,jy,jz),AddPressureZ
        !     WRITE(*,*) 'ym, zm, ct, zp, yp, rhs = ',coef(-2),coef(-3),coef(0),coef(3),coef(2),BvecCrunchP(j)
        !     WRITE(*,*) ' ---------- '
        ! END IF

    END IF
  END DO
END DO

! zp face
jz = nz
DO jy = 2,ny-1
  DO jx = 2,nx-1
    j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1
    IF (activecellPressure(jx,jy,jz) == 0) THEN
        CALL MatSetValues(amatP,1,j,1,j-1,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+1,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx*ny,zero,INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j+nx*ny,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j,big,INSERT_VALUES,ierr)
        BvecCrunchP(j) = head(jx,jy,jz)*big
        XvecCrunchP(j) = head(jx,jy,jz)
    ELSE

        COORDINATE = 'X'
        coef(-1) = -2.0d0*dt*Kfacx(jx-1,jy,jz) / (dxx(jx)*(dxx(jx)+dxx(jx-1)))
        coef(1) = -2.0d0*dt*Kfacx(jx,jy,jz) / (dxx(jx)*(dxx(jx+1)+dxx(jx)))

        COORDINATE = 'Y'
        coef(2)  = -2.0d0*dt*Kfacy(jx,jy,jz) / (dyy(jy)*(dyy(jy+1)+dyy(jy) ))
        coef(-2) = -2.0d0*dt*Kfacy(jx,jy-1,jz) / (dyy(jy)*(dyy(jy)+dyy(jy-1)))

        COORDINATE = 'Z'
        coef(-3) = -2.0d0*dt*Kfacz(jx,jy,jz-1) / (dzz(jx,jy,jz)*(dzz(jx,jy,jz)+dzz(jx,jy,jz-1)))
        IF (activecellPressure(jx,jy,nz+1) == 0) THEN
            coef(3)  = -2.0d0*dt*Kfacz(jx,jy,jz) / (dzz(jx,jy,jz)*(dzz(jx,jy,jz+1)+dzz(jx,jy,jz) ))
        ELSE
            coef(3) = 0.0d0
        END IF

        coef(0) = (Ch(jx,jy,jz)+Ss*wc(jx,jy,jz)/wcs(jx,jy,jz)) - coef(-3) - coef(-2) - coef(-1) - coef(1) - coef(2) - coef(3)

        CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx*ny,coef(-3),INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j+nx*ny,coef(3),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)

        IF (activecellPressure(jx,jy,nz+1) == 0) THEN
            AddPressureZ =  -coef(3) * headOld(jx,jy,nz+1)
        ELSE
            AddPressureZ = 0.0d0
        END IF

        BvecCrunchP(j) = (Ch(jx,jy,jz)+Ss*wc(jx,jy,jz)/wcs(jx,jy,jz))*headOld(jx,jy,jz) - (dt/dzz(jx,jy,jz))*(Kfacz(jx,jy,jz)-Kfacz(jx,jy,jz-1)) + AddPressureZ
        XvecCrunchP(j) = head(jx,jy,jz)

    END IF
  END DO
END DO
!!  ************************************************************************************************************************************************
!!  Now the edges

! im-jm edge
  jx = 1
  jy = 1
  DO jz = 2,nz-1
    j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1
    IF (activecellPressure(jx,jy,jz) == 0) THEN
        !!            CALL MatSetValues(amatP,1,j,1,j-1,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+1,zero,INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j-nx,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx*ny,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx*ny,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j,big,INSERT_VALUES,ierr)
        BvecCrunchP(j) = head(jx,jy,jz)*big
        XvecCrunchP(j) = head(jx,jy,jz)
    ELSE

        COORDINATE = 'X'
        coef(1) = -2.0d0*dt*Kfacx(jx,jy,jz) / (dxx(jx)*(dxx(jx+1)+dxx(jx)))
        IF (activecellPressure(0,jy,jz) == 0) THEN
            coef(-1) = -2.0d0*dt*Kfacx(jx-1,jy,jz) / (dxx(jx)*(dxx(jx)+dxx(jx-1)))
        ELSE
            coef(-1) = 0.0d0
        END IF

        COORDINATE = 'Y'
        coef(2)  = -2.0d0*dt*Kfacy(jx,jy,jz) / (dyy(jy)*(dyy(jy+1)+dyy(jy) ))

        IF (activecellPressure(jx,0,jz) == 0) THEN
            coef(-2) = -2.0d0*dt*Kfacy(jx,jy-1,jz) / (dyy(jy)*(dyy(jy)+dyy(jy-1)))
        ELSE
            coef(-2) = 0.0d0
        END IF

        COORDINATE = 'Z'
        coef(-3) = -2.0d0*dt*Kfacz(jx,jy,jz-1) / (dzz(jx,jy,jz)*(dzz(jx,jy,jz)+dzz(jx,jy,jz-1)))
        coef(3)  = -2.0d0*dt*Kfacz(jx,jy,jz) / (dzz(jx,jy,jz)*(dzz(jx,jy,jz+1)+dzz(jx,jy,jz) ))

        coef(0) = (Ch(jx,jy,jz)+Ss*wc(jx,jy,jz)/wcs(jx,jy,jz)) - coef(-3) - coef(-2) - coef(-1) - coef(1) - coef(2) - coef(3)

        !!            CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx*ny,coef(-3),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx*ny,coef(3),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)

        IF (activecellPressure(0,jy,jz) == 0) THEN
            AddPressureX =  -coef(-1) * headOld(0,jy,jz)
        ELSE
            AddPressureX = 0.0d0
        END IF
        IF (activecellPressure(jx,0,jz) == 0) THEN
            AddPressureY =  -coef(-2) * headOld(jx,0,jz)
        ELSE
            AddPressureY = 0.0d0
        END IF

        BvecCrunchP(j) = (Ch(jx,jy,jz)+Ss*wc(jx,jy,jz)/wcs(jx,jy,jz))*headOld(jx,jy,jz) - (dt/dzz(jx,jy,jz))*(Kfacz(jx,jy,jz)-Kfacz(jx,jy,jz-1)) + AddPressureX + AddPressureY
        XvecCrunchP(j) = head(jx,jy,jz)
    END IF
  END DO

!! ip-jm edge
  jx = nx
  jy = 1
  DO jz = 2,nz-1
    j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1
    IF (activecellPressure(jx,jy,jz) == 0) THEN
        CALL MatSetValues(amatP,1,j,1,j-1,zero,INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j+1,zero,INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j-nx,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx*ny,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx*ny,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j,big,INSERT_VALUES,ierr)
        BvecCrunchP(j) = head(jx,jy,jz)*big
        XvecCrunchP(j) = head(jx,jy,jz)
    ELSE
        COORDINATE = 'X'
        coef(-1) = -2.0d0*dt*Kfacx(jx-1,jy,jz) / (dxx(jx)*(dxx(jx)+dxx(jx-1)))
        IF (activecellPressure(nx+1,jy,jz) == 0) THEN
            coef(1) = -2.0d0*dt*Kfacx(jx,jy,jz) / (dxx(jx)*(dxx(jx+1)+dxx(jx)))
        ELSE
            coef(1) = 0.0d0
        END IF

        COORDINATE = 'Y'
        coef(2)  = -2.0d0*dt*Kfacy(jx,jy,jz) / (dyy(jy)*(dyy(jy+1)+dyy(jy) ))
        IF (activecellPressure(jx,0,jz) == 0) THEN
            coef(-2) = -2.0d0*dt*Kfacy(jx,jy-1,jz) / (dyy(jy)*(dyy(jy)+dyy(jy-1)))
        ELSE
            coef(-2) = 0.0d0
        END IF

        COORDINATE = 'Z'
        coef(-3) = -2.0d0*dt*Kfacz(jx,jy,jz-1) / (dzz(jx,jy,jz)*(dzz(jx,jy,jz)+dzz(jx,jy,jz-1)))
        coef(3)  = -2.0d0*dt*Kfacz(jx,jy,jz) / (dzz(jx,jy,jz)*(dzz(jx,jy,jz+1)+dzz(jx,jy,jz) ))

        coef(0) = (Ch(jx,jy,jz)+Ss*wc(jx,jy,jz)/wcs(jx,jy,jz)) - coef(-3) - coef(-2) - coef(-1) - coef(1) - coef(2) - coef(3)

        CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx*ny,coef(-3),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx*ny,coef(3),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)

        IF (activecellPressure(nx+1,jy,jz) == 0) THEN
          AddPressureX =  -coef(1) * headOld(nx+1,jy,jz)
        ELSE
          AddPressureX = 0.0d0
        END IF
        IF (activecellPressure(jx,0,jz) == 0) THEN
          AddPressureY =  -coef(-2) * headOld(jx,0,jz)
        ELSE
          AddPressureY = 0.0d0
        END IF

        BvecCrunchP(j) = (Ch(jx,jy,jz)+Ss*wc(jx,jy,jz)/wcs(jx,jy,jz))*headOld(jx,jy,jz) - (dt/dzz(jx,jy,jz))*(Kfacz(jx,jy,jz)-Kfacz(jx,jy,jz-1)) + AddPressureX + AddPressureY
        XvecCrunchP(j) = head(jx,jy,jz)

    END IF
  END DO

!! ip-jp edge
  jx = nx
  jy = ny
  DO jz = 2,nz-1
    j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1
    IF (activecellPressure(jx,jy,jz) == 0) THEN
        CALL MatSetValues(amatP,1,j,1,j-1,zero,INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j+1,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx,zero,INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j+nx,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx*ny,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx*ny,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j,big,INSERT_VALUES,ierr)
        BvecCrunchP(j) = head(jx,jy,jz)*big
        XvecCrunchP(j) = head(jx,jy,jz)
    ELSE
        COORDINATE = 'X'
        coef(-1) = -2.0d0*dt*Kfacx(jx-1,jy,jz) / (dxx(jx)*(dxx(jx)+dxx(jx-1)))
        IF (activecellPressure(nx+1,jy,jz) == 0) THEN
            coef(1) = -2.0d0*dt*Kfacx(jx,jy,jz) / (dxx(jx)*(dxx(jx+1)+dxx(jx)))
        ELSE
            coef(1) = 0.0d0
        END IF

        COORDINATE = 'Y'
        coef(-2) = -2.0d0*dt*Kfacy(jx,jy-1,jz) / (dyy(jy)*(dyy(jy)+dyy(jy-1)))
        IF (activecellPressure(jx,ny+1,jz) == 0) THEN
            coef(2)  = -2.0d0*dt*Kfacy(jx,jy,jz) / (dyy(jy)*(dyy(jy+1)+dyy(jy) ))
        ELSE
            coef(2) = 0.0d0
        END IF

        COORDINATE = 'Z'
        coef(-3) = -2.0d0*dt*Kfacz(jx,jy,jz-1) / (dzz(jx,jy,jz)*(dzz(jx,jy,jz)+dzz(jx,jy,jz-1)))
        coef(3)  = -2.0d0*dt*Kfacz(jx,jy,jz) / (dzz(jx,jy,jz)*(dzz(jx,jy,jz+1)+dzz(jx,jy,jz) ))

        coef(0) = (Ch(jx,jy,jz)+Ss*wc(jx,jy,jz)/wcs(jx,jy,jz)) - coef(-3) - coef(-2) - coef(-1) - coef(1) - coef(2) - coef(3)

        CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx*ny,coef(-3),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx*ny,coef(3),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)
        IF (activecellPressure(nx+1,jy,jz) == 0) THEN
          AddPressureX =  -coef(1) * headOld(nx+1,jy,jz)
        ELSE
          AddPressureX = 0.0d0
        END IF
        IF (activecellPressure(jx,ny+1,jz) == 0) THEN
          AddPressureY =  -coef(2) * headOld(jx,ny+1,jz)
        ELSE
          AddPressureY = 0.0d0
        END IF
        BvecCrunchP(j) = (Ch(jx,jy,jz)+Ss*wc(jx,jy,jz)/wcs(jx,jy,jz))*headOld(jx,jy,jz) - (dt/dzz(jx,jy,jz))*(Kfacz(jx,jy,jz)-Kfacz(jx,jy,jz-1)) + AddPressureX + AddPressureY
        XvecCrunchP(j) = head(jx,jy,jz)

    END IF
  END DO

!! im-jp edge
  jx = 1
  jy = ny
  DO jz = 2,nz-1
    j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1
    IF (activecellPressure(jx,jy,jz) == 0) THEN
        !!            CALL MatSetValues(amatP,1,j,1,j-1,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+1,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx,zero,INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j+nx,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx*ny,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx*ny,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j,big,INSERT_VALUES,ierr)
        BvecCrunchP(j) = head(jx,jy,jz)*big
        XvecCrunchP(j) = head(jx,jy,jz)
    ELSE
        COORDINATE = 'X'
        coef(1) = -2.0d0*dt*Kfacx(jx,jy,jz) / (dxx(jx)*(dxx(jx+1)+dxx(jx)))
        IF (activecellPressure(0,jy,jz) == 0) THEN
            coef(-1) = -2.0d0*dt*Kfacx(jx-1,jy,jz) / (dxx(jx)*(dxx(jx)+dxx(jx-1)))
        ELSE
            coef(-1) = 0.0d0
        END IF

        COORDINATE = 'Y'
        coef(-2) = -2.0d0*dt*Kfacy(jx,jy-1,jz) / (dyy(jy)*(dyy(jy)+dyy(jy-1)))
        IF (activecellPressure(jx,ny+1,jz) == 0) THEN
            coef(2)  = -2.0d0*dt*Kfacy(jx,jy,jz) / (dyy(jy)*(dyy(jy+1)+dyy(jy) ))
        ELSE
            coef(2) = 0.0d0
        END IF

        COORDINATE = 'Z'
        coef(-3) = -2.0d0*dt*Kfacz(jx,jy,jz-1) / (dzz(jx,jy,jz)*(dzz(jx,jy,jz)+dzz(jx,jy,jz-1)))
        coef(3)  = -2.0d0*dt*Kfacz(jx,jy,jz) / (dzz(jx,jy,jz)*(dzz(jx,jy,jz+1)+dzz(jx,jy,jz) ))

        coef(0) = (Ch(jx,jy,jz)+Ss*wc(jx,jy,jz)/wcs(jx,jy,jz)) - coef(-3) - coef(-2) - coef(-1) - coef(1) - coef(2) - coef(3)

        !!            CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx*ny,coef(-3),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx*ny,coef(3),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)

        IF (activecellPressure(0,jy,jz) == 0) THEN
          AddPressureX =  -coef(-1) * headOld(0,jy,jz)
        ELSE
          AddPressureX = 0.0d0
        END IF
        IF (activecellPressure(jx,ny+1,jz) == 0) THEN
          AddPressureY =  -coef(2) * headOld(jx,ny+1,jz)
        ELSE
          AddPressureY = 0.0d0
        END IF
        BvecCrunchP(j) = (Ch(jx,jy,jz)+Ss*wc(jx,jy,jz)/wcs(jx,jy,jz))*headOld(jx,jy,jz) - (dt/dzz(jx,jy,jz))*(Kfacz(jx,jy,jz)-Kfacz(jx,jy,jz-1)) + AddPressureX + AddPressureY
        XvecCrunchP(j) = head(jx,jy,jz)

    END IF
  END DO

!! im-km edge
  jz = 1
  jx = 1
  DO jy = 2,ny-1
    j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1
    IF (activecellPressure(jx,jy,jz) == 0) THEN
        !!            CALL MatSetValues(amatP,1,j,1,j-1,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+1,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx,zero,INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j-nx*ny,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx*ny,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j,big,INSERT_VALUES,ierr)
        BvecCrunchP(j) = head(jx,jy,jz)*big
        XvecCrunchP(j) = head(jx,jy,jz)
    ELSE
        COORDINATE = 'X'
        coef(1) = -2.0d0*dt*Kfacx(jx,jy,jz) / (dxx(jx)*(dxx(jx+1)+dxx(jx)))
        IF (activecellPressure(0,jy,jz) == 0) THEN
            coef(-1) = -2.0d0*dt*Kfacx(jx-1,jy,jz) / (dxx(jx)*(dxx(jx)+dxx(jx-1)))
        ELSE
            coef(-1) = 0.0d0
        END IF

        COORDINATE = 'Y'
        coef(2)  = -2.0d0*dt*Kfacy(jx,jy,jz) / (dyy(jy)*(dyy(jy+1)+dyy(jy) ))
        coef(-2) = -2.0d0*dt*Kfacy(jx,jy-1,jz) / (dyy(jy)*(dyy(jy)+dyy(jy-1)))

        COORDINATE = 'Z'
        coef(3)  = -2.0d0*dt*Kfacz(jx,jy,jz) / (dzz(jx,jy,jz)*(dzz(jx,jy,jz+1)+dzz(jx,jy,jz) ))
        IF (activecellPressure(jx,jy,0) == 0) THEN
            coef(-3) = -2.0d0*dt*Kfacz(jx,jy,jz-1) / (dzz(jx,jy,jz)*(dzz(jx,jy,jz)+dzz(jx,jy,jz-1)))
        ELSE
            coef(-3) = 0.0d0
        END IF

        coef(0) = (Ch(jx,jy,jz)+Ss*wc(jx,jy,jz)/wcs(jx,jy,jz)) - coef(-3) - coef(-2) - coef(-1) - coef(1) - coef(2) - coef(3)

        !!            CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j-nx*ny,coef(-3),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx*ny,coef(3),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)

        IF (activecellPressure(0,jy,jz) == 0) THEN
          AddPressureX =  -coef(-1) * headOld(0,jy,jz)
        ELSE
          AddPressureX = 0.0d0
        END IF
        IF (activecellPressure(jx,jy,0) == 0) THEN
          AddPressureZ =  -coef(-3) * headOld(jx,jy,0)
        ELSE
          AddPressureZ = 0.0d0
        END IF
        BvecCrunchP(j) = (Ch(jx,jy,jz)+Ss*wc(jx,jy,jz)/wcs(jx,jy,jz))*headOld(jx,jy,jz) - (dt/dzz(jx,jy,jz))*(Kfacz(jx,jy,jz)-Kfacz(jx,jy,jz-1)) + AddPressureX + AddPressureZ
        XvecCrunchP(j) = head(jx,jy,jz)

    END IF
  END DO

!! ip-km edge
  jz = 1
  jx = nx
  DO jy = 2,ny-1
    j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1
    IF (activecellPressure(jx,jy,jz) == 0) THEN
        CALL MatSetValues(amatP,1,j,1,j-1,zero,INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j+1,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx,zero,INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j-nx*ny,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx*ny,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j,big,INSERT_VALUES,ierr)
        BvecCrunchP(j) = head(jx,jy,jz)*big
        XvecCrunchP(j) = head(jx,jy,jz)
    ELSE
        COORDINATE = 'X'
        coef(-1) = -2.0d0*dt*Kfacx(jx-1,jy,jz) / (dxx(jx)*(dxx(jx)+dxx(jx-1)))
        IF (activecellPressure(nx+1,jy,jz) == 0) THEN
            coef(1) = -2.0d0*dt*Kfacx(jx,jy,jz) / (dxx(jx)*(dxx(jx+1)+dxx(jx)))
        ELSE
            coef(1) = 0.0d0
        END IF

        COORDINATE = 'Y'
        coef(2)  = -2.0d0*dt*Kfacy(jx,jy,jz) / (dyy(jy)*(dyy(jy+1)+dyy(jy) ))
        coef(-2) = -2.0d0*dt*Kfacy(jx,jy-1,jz) / (dyy(jy)*(dyy(jy)+dyy(jy-1)))

        COORDINATE = 'Z'
        coef(3)  = -2.0d0*dt*Kfacz(jx,jy,jz) / (dzz(jx,jy,jz)*(dzz(jx,jy,jz+1)+dzz(jx,jy,jz) ))
        IF (activecellPressure(jx,jy,0) == 0) THEN
            coef(-3) = -2.0d0*dt*Kfacz(jx,jy,jz-1) / (dzz(jx,jy,jz)*(dzz(jx,jy,jz)+dzz(jx,jy,jz-1)))
        ELSE
            coef(-3) = 0.0d0
        END IF

        coef(0) = (Ch(jx,jy,jz)+Ss*wc(jx,jy,jz)/wcs(jx,jy,jz)) - coef(-3) - coef(-2) - coef(-1) - coef(1) - coef(2) - coef(3)

        IF (activecellPressure(nx+1,jy,jz) == 0) THEN
          AddPressureX =  -coef(1) * headOld(nx+1,jy,jz)
        ELSE
          AddPressureX = 0.0d0
        END IF
        IF (activecellPressure(jx,jy,0) == 0) THEN
          AddPressureZ =  -coef(-3) * headOld(jx,jy,0)
        ELSE
          AddPressureZ = 0.0d0
        END IF
        BvecCrunchP(j) = (Ch(jx,jy,jz)+Ss*wc(jx,jy,jz)/wcs(jx,jy,jz))*headOld(jx,jy,jz) - (dt/dzz(jx,jy,jz))*(Kfacz(jx,jy,jz)-Kfacz(jx,jy,jz-1)) + AddPressureX + AddPressureZ
        XvecCrunchP(j) = head(jx,jy,jz)

        CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j-nx*ny,coef(-3),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx*ny,coef(3),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)
    END IF
  END DO

!! im-kp edges
  jz = nz
  jx = 1
  DO jy = 2,ny-1
    j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1
    IF (activecellPressure(jx,jy,jz) == 0) THEN
        !!            CALL MatSetValues(amatP,1,j,1,j-1,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+1,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx*ny,zero,INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j+nx*ny,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j,big,INSERT_VALUES,ierr)
        BvecCrunchP(j) = head(jx,jy,jz)*big
        XvecCrunchP(j) = head(jx,jy,jz)
    ELSE
        COORDINATE = 'X'
        coef(1) = -2.0d0*dt*Kfacx(jx,jy,jz) / (dxx(jx)*(dxx(jx+1)+dxx(jx)))
        IF (activecellPressure(0,jy,jz) == 0) THEN
            coef(-1) = -2.0d0*dt*Kfacx(jx-1,jy,jz) / (dxx(jx)*(dxx(jx)+dxx(jx-1)))
        ELSE
            coef(-1) = 0.0d0
        END IF

        COORDINATE = 'Y'
        coef(2)  = -2.0d0*dt*Kfacy(jx,jy,jz) / (dyy(jy)*(dyy(jy+1)+dyy(jy) ))
        coef(-2) = -2.0d0*dt*Kfacy(jx,jy-1,jz) / (dyy(jy)*(dyy(jy)+dyy(jy-1)))

        COORDINATE = 'Z'
        coef(-3) = -2.0d0*dt*Kfacz(jx,jy,jz-1) / (dzz(jx,jy,jz)*(dzz(jx,jy,jz)+dzz(jx,jy,jz-1)))
        IF (activecellPressure(jx,jy,nz+1) == 0) THEN
            coef(3)  = -2.0d0*dt*Kfacz(jx,jy,jz) / (dzz(jx,jy,jz)*(dzz(jx,jy,jz+1)+dzz(jx,jy,jz) ))
        ELSE
            coef(3) = 0.0d0
        END IF

        coef(0) = (Ch(jx,jy,jz)+Ss*wc(jx,jy,jz)/wcs(jx,jy,jz)) - coef(-3) - coef(-2) - coef(-1) - coef(1) - coef(2) - coef(3)

        IF (activecellPressure(0,jy,jz) == 0) THEN
          AddPressureX =  -coef(-1) * headOld(0,jy,jz)
        ELSE
          AddPressureX = 0.0d0
        END IF
        IF (activecellPressure(jx,jy,nz+1) == 0) THEN
          AddPressureZ =  -coef(3) * headOld(jx,jy,nz+1)
        ELSE
          AddPressureZ = 0.0d0
        END IF
        BvecCrunchP(j) = (Ch(jx,jy,jz)+Ss*wc(jx,jy,jz)/wcs(jx,jy,jz))*headOld(jx,jy,jz) - (dt/dzz(jx,jy,jz))*(Kfacz(jx,jy,jz)-Kfacz(jx,jy,jz-1)) + AddPressureX + AddPressureZ
        XvecCrunchP(j) = head(jx,jy,jz)

        !!            CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx*ny,coef(-3),INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j+nx*ny,coef(3),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)
    END IF
  END DO

!! ip-kp edge
  jz = nz
  jx = nx
  DO jy = 2,ny-1
    j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1
    IF (activecellPressure(jx,jy,jz) == 0) THEN
        CALL MatSetValues(amatP,1,j,1,j-1,zero,INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j+1,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx*ny,zero,INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j+nx*ny,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j,big,INSERT_VALUES,ierr)
        BvecCrunchP(j) = head(jx,jy,jz)*big
        XvecCrunchP(j) = head(jx,jy,jz)
    ELSE
        COORDINATE = 'X'
        coef(-1) = -2.0d0*dt*Kfacx(jx-1,jy,jz) / (dxx(jx)*(dxx(jx)+dxx(jx-1)))
        IF (activecellPressure(nx+1,jy,jz) == 0) THEN
            coef(1) = -2.0d0*dt*Kfacx(jx,jy,jz) / (dxx(jx)*(dxx(jx+1)+dxx(jx)))
        ELSE
            coef(1) = 0.0d0
        END IF

        COORDINATE = 'Y'
        coef(2)  = -2.0d0*dt*Kfacy(jx,jy,jz) / (dyy(jy)*(dyy(jy+1)+dyy(jy) ))
        coef(-2) = -2.0d0*dt*Kfacy(jx,jy-1,jz) / (dyy(jy)*(dyy(jy)+dyy(jy-1)))

        COORDINATE = 'Z'
        coef(-3) = -2.0d0*dt*Kfacz(jx,jy,jz-1) / (dzz(jx,jy,jz)*(dzz(jx,jy,jz)+dzz(jx,jy,jz-1)))
        IF (activecellPressure(jx,jy,nz+1) == 0) THEN
            coef(3)  = -2.0d0*dt*Kfacz(jx,jy,jz) / (dzz(jx,jy,jz)*(dzz(jx,jy,jz+1)+dzz(jx,jy,jz) ))
        ELSE
            coef(3) = 0.0d0
        END IF

        coef(0) = (Ch(jx,jy,jz)+Ss*wc(jx,jy,jz)/wcs(jx,jy,jz)) - coef(-3) - coef(-2) - coef(-1) - coef(1) - coef(2) - coef(3)

        IF (activecellPressure(nx+1,jy,jz) == 0) THEN
          AddPressureX =  -coef(1) * headOld(nx+1,jy,jz)
        ELSE
          AddPressureX = 0.0d0
        END IF
        IF (activecellPressure(jx,jy,nz+1) == 0) THEN
          AddPressureZ =  -coef(3) * headOld(jx,jy,nz+1)
        ELSE
          AddPressureZ = 0.0d0
        END IF
        BvecCrunchP(j) = (Ch(jx,jy,jz)+Ss*wc(jx,jy,jz)/wcs(jx,jy,jz))*headOld(jx,jy,jz) - (dt/dzz(jx,jy,jz))*(Kfacz(jx,jy,jz)-Kfacz(jx,jy,jz-1)) + AddPressureX + AddPressureZ
        XvecCrunchP(j) = head(jx,jy,jz)

        CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx*ny,coef(-3),INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j+nx*ny,coef(3),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)
    END IF
  END DO

!! jm-km edge
  jz = 1
  jy = 1
  DO jx = 2,nx-1
    j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1
    IF (activecellPressure(jx,jy,jz) == 0) THEN
        CALL MatSetValues(amatP,1,j,1,j-1,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+1,zero,INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j-nx,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx,zero,INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j-nx*ny,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx*ny,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j,big,INSERT_VALUES,ierr)
        BvecCrunchP(j) = head(jx,jy,jz)*big
        XvecCrunchP(j) = head(jx,jy,jz)
    ELSE
        COORDINATE = 'X'
        coef(1) = -2.0d0*dt*Kfacx(jx,jy,jz) / (dxx(jx)*(dxx(jx+1)+dxx(jx)))
        coef(-1) = -2.0d0*dt*Kfacx(jx-1,jy,jz) / (dxx(jx)*(dxx(jx)+dxx(jx-1)))

        COORDINATE = 'Y'
        coef(2)  = -2.0d0*dt*Kfacy(jx,jy,jz) / (dyy(jy)*(dyy(jy+1)+dyy(jy) ))
        IF (activecellPressure(jx,0,jz) == 0) THEN
            coef(-2) = -2.0d0*dt*Kfacy(jx,jy-1,jz) / (dyy(jy)*(dyy(jy)+dyy(jy-1)))
        ELSE
            coef(-2) = 0.0d0
        END IF

        COORDINATE = 'Z'
        coef(3)  = -2.0d0*dt*Kfacz(jx,jy,jz) / (dzz(jx,jy,jz)*(dzz(jx,jy,jz+1)+dzz(jx,jy,jz) ))
        IF (activecellPressure(jx,jy,0) == 0) THEN
            coef(-3) = -2.0d0*dt*Kfacz(jx,jy,jz-1) / (dzz(jx,jy,jz)*(dzz(jx,jy,jz)+dzz(jx,jy,jz-1)))
        ELSE
            coef(-3) = 0.0d0
        END IF

        coef(0) = (Ch(jx,jy,jz)+Ss*wc(jx,jy,jz)/wcs(jx,jy,jz)) - coef(-3) - coef(-2) - coef(-1) - coef(1) - coef(2) - coef(3)

        IF (activecellPressure(jx,0,jz) == 0) THEN
          AddPressureY =  -coef(-2) * headOld(jx,0,jz)
        ELSE
          AddPressureY = 0.0d0
        END IF
        IF (activecellPressure(jx,jy,0) == 0) THEN
          AddPressureZ =  -coef(-3) * headOld(jx,jy,0)
        ELSE
          AddPressureZ = 0.0d0
        END IF
        BvecCrunchP(j) = (Ch(jx,jy,jz)+Ss*wc(jx,jy,jz)/wcs(jx,jy,jz))*headOld(jx,jy,jz) - (dt/dzz(jx,jy,jz))*(Kfacz(jx,jy,jz)-Kfacz(jx,jy,jz-1)) + AddPressureY + AddPressureZ
        XvecCrunchP(j) = head(jx,jy,jz)

        CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j-nx*ny,coef(-3),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx*ny,coef(3),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)
    END IF
  END DO

!! jp-km
  jz = 1
  jy = ny
  DO jx = 2,nx-1
    j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1
    IF (activecellPressure(jx,jy,jz) == 0) THEN
        CALL MatSetValues(amatP,1,j,1,j-1,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+1,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx,zero,INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j+nx,zero,INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j-nx*ny,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx*ny,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j,big,INSERT_VALUES,ierr)
        BvecCrunchP(j) = head(jx,jy,jz)*big
        XvecCrunchP(j) = head(jx,jy,jz)
    ELSE
        COORDINATE = 'X'
        coef(1) = -2.0d0*dt*Kfacx(jx,jy,jz) / (dxx(jx)*(dxx(jx+1)+dxx(jx)))
        coef(-1) = -2.0d0*dt*Kfacx(jx-1,jy,jz) / (dxx(jx)*(dxx(jx)+dxx(jx-1)))

        COORDINATE = 'Y'
        coef(-2) = -2.0d0*dt*Kfacy(jx,jy-1,jz) / (dyy(jy)*(dyy(jy)+dyy(jy-1)))
        IF (activecellPressure(jx,ny+1,jz) == 0) THEN
            coef(2)  = -2.0d0*dt*Kfacy(jx,jy,jz) / (dyy(jy)*(dyy(jy+1)+dyy(jy) ))
        ELSE
            coef(2) = 0.0d0
        END IF

        COORDINATE = 'Z'
        coef(3)  = -2.0d0*dt*Kfacz(jx,jy,jz) / (dzz(jx,jy,jz)*(dzz(jx,jy,jz+1)+dzz(jx,jy,jz) ))
        IF (activecellPressure(jx,jy,0) == 0) THEN
            coef(-3) = -2.0d0*dt*Kfacz(jx,jy,jz-1) / (dzz(jx,jy,jz)*(dzz(jx,jy,jz)+dzz(jx,jy,jz-1)))
        ELSE
            coef(-3) = 0.0d0
        END IF

        coef(0) = (Ch(jx,jy,jz)+Ss*wc(jx,jy,jz)/wcs(jx,jy,jz)) - coef(-3) - coef(-2) - coef(-1) - coef(1) - coef(2) - coef(3)

        IF (activecellPressure(jx,ny+1,jz) == 0) THEN
          AddPressureY =  -coef(2) * headOld(jx,ny+1,jz)
        ELSE
          AddPressureY = 0.0d0
        END IF
        IF (activecellPressure(jx,jy,0) == 0) THEN
          AddPressureZ =  -coef(-3) * headOld(jx,jy,0)
        ELSE
          AddPressureZ = 0.0d0
        END IF
        BvecCrunchP(j) = (Ch(jx,jy,jz)+Ss*wc(jx,jy,jz)/wcs(jx,jy,jz))*headOld(jx,jy,jz) - (dt/dzz(jx,jy,jz))*(Kfacz(jx,jy,jz)-Kfacz(jx,jy,jz-1)) + AddPressureY + AddPressureZ
        XvecCrunchP(j) = head(jx,jy,jz)

        CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j-nx*ny,coef(-3),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx*ny,coef(3),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)
    END IF
  END DO

!! jm-kp
  jz = nz
  jy = 1
  DO jx = 2,nx-1
    j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1
    IF (activecellPressure(jx,jy,jz) == 0) THEN
        CALL MatSetValues(amatP,1,j,1,j-1,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+1,zero,INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j-nx,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx*ny,zero,INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j+nx*ny,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j,big,INSERT_VALUES,ierr)
        BvecCrunchP(j) = head(jx,jy,jz)*big
        XvecCrunchP(j) = head(jx,jy,jz)
    ELSE
        COORDINATE = 'X'
        coef(1) = -2.0d0*dt*Kfacx(jx,jy,jz) / (dxx(jx)*(dxx(jx+1)+dxx(jx)))
        coef(-1) = -2.0d0*dt*Kfacx(jx-1,jy,jz) / (dxx(jx)*(dxx(jx)+dxx(jx-1)))

        COORDINATE = 'Y'
        coef(2)  = -2.0d0*dt*Kfacy(jx,jy,jz) / (dyy(jy)*(dyy(jy+1)+dyy(jy) ))
        IF (activecellPressure(jx,0,jz) == 0) THEN
            coef(-2) = -2.0d0*dt*Kfacy(jx,jy-1,jz) / (dyy(jy)*(dyy(jy)+dyy(jy-1)))
        ELSE
            coef(-2) = 0.0d0
        END IF

        COORDINATE = 'Z'
        coef(-3) = -2.0d0*dt*Kfacz(jx,jy,jz-1) / (dzz(jx,jy,jz)*(dzz(jx,jy,jz)+dzz(jx,jy,jz-1)))
        IF (activecellPressure(jx,jy,nz+1) == 0) THEN
            coef(3)  = -2.0d0*dt*Kfacz(jx,jy,jz) / (dzz(jx,jy,jz)*(dzz(jx,jy,jz+1)+dzz(jx,jy,jz) ))
        ELSE
            coef(3) = 0.0d0
        END IF

        coef(0) = (Ch(jx,jy,jz)+Ss*wc(jx,jy,jz)/wcs(jx,jy,jz)) - coef(-3) - coef(-2) - coef(-1) - coef(1) - coef(2) - coef(3)

        IF (activecellPressure(jx,0,jz) == 0) THEN
          AddPressureY =  -coef(-2) * headOld(jx,0,jz)
        ELSE
          AddPressureY = 0.0d0
        END IF
        IF (activecellPressure(jx,jy,nz+1) == 0) THEN
          AddPressureZ =  -coef(3) * headOld(jx,jy,nz+1)
        ELSE
          AddPressureZ = 0.0d0
        END IF
        BvecCrunchP(j) = (Ch(jx,jy,jz)+Ss*wc(jx,jy,jz)/wcs(jx,jy,jz))*headOld(jx,jy,jz) - (dt/dzz(jx,jy,jz))*(Kfacz(jx,jy,jz)-Kfacz(jx,jy,jz-1)) + AddPressureY + AddPressureZ
        XvecCrunchP(j) = head(jx,jy,jz)

        CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx*ny,coef(-3),INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j+nx*ny,coef(3),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)

    END IF
  END DO

!! jp-kp
  jz = nz
  jy = ny
  DO jx = 2,nx-1
    j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1
    IF (activecellPressure(jx,jy,jz) == 0) THEN
        CALL MatSetValues(amatP,1,j,1,j-1,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+1,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx,zero,INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j+nx,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx*ny,zero,INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j+nx*ny,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j,big,INSERT_VALUES,ierr)
        BvecCrunchP(j) = head(jx,jy,jz)*big
        XvecCrunchP(j) = head(jx,jy,jz)
    ELSE
        COORDINATE = 'X'
        coef(1) = -2.0d0*dt*Kfacx(jx,jy,jz) / (dxx(jx)*(dxx(jx+1)+dxx(jx)))
        coef(-1) = -2.0d0*dt*Kfacx(jx-1,jy,jz) / (dxx(jx)*(dxx(jx)+dxx(jx-1)))

        COORDINATE = 'Y'
        coef(-2) = -2.0d0*dt*Kfacy(jx,jy-1,jz) / (dyy(jy)*(dyy(jy)+dyy(jy-1)))
        IF (activecellPressure(jx,ny+1,jz) == 0) THEN
            coef(2)  = -2.0d0*dt*Kfacy(jx,jy,jz) / (dyy(jy)*(dyy(jy+1)+dyy(jy) ))
        ELSE
            coef(2) = 0.0d0
        END IF

        COORDINATE = 'Z'
        coef(-3) = -2.0d0*dt*Kfacz(jx,jy,jz-1) / (dzz(jx,jy,jz)*(dzz(jx,jy,jz)+dzz(jx,jy,jz-1)))
        IF (activecellPressure(jx,jy,nz+1) == 0) THEN
            coef(3)  = -2.0d0*dt*Kfacz(jx,jy,jz) / (dzz(jx,jy,jz)*(dzz(jx,jy,jz+1)+dzz(jx,jy,jz) ))
        ELSE
            coef(3) = 0.0d0
        END IF

        coef(0) = (Ch(jx,jy,jz)+Ss*wc(jx,jy,jz)/wcs(jx,jy,jz)) - coef(-3) - coef(-2) - coef(-1) - coef(1) - coef(2) - coef(3)

        IF (activecellPressure(jx,ny+1,jz) == 0) THEN
          AddPressureY =  -coef(2) * headOld(jx,ny+1,jz)
        ELSE
          AddPressureY = 0.0d0
        END IF
        IF (activecellPressure(jx,jy,nz+1) == 0) THEN
          AddPressureZ =  -coef(3) * headOld(jx,jy,nz+1)
        ELSE
          AddPressureZ = 0.0d0
        END IF
        BvecCrunchP(j) = (Ch(jx,jy,jz)+Ss*wc(jx,jy,jz)/wcs(jx,jy,jz))*headOld(jx,jy,jz) - (dt/dzz(jx,jy,jz))*(Kfacz(jx,jy,jz)-Kfacz(jx,jy,jz-1)) + AddPressureY + AddPressureZ
        XvecCrunchP(j) = head(jx,jy,jz)

        CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx*ny,coef(-3),INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j+nx*ny,coef(3),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)

    END IF
  END DO


!!  **************************************************************************************************************************************************
!!  Now the corners

!! im-jm-km corner
    jx = 1
    jy = 1
    jz = 1
    j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1
    IF (activecellPressure(jx,jy,jz) == 0) THEN
        !!            CALL MatSetValues(amatP,1,j,1,j-1,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+1,zero,INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j-nx,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx,zero,INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j-nx*ny,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx*ny,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j,big,INSERT_VALUES,ierr)
        BvecCrunchP(j) = head(jx,jy,jz)*big
        XvecCrunchP(j) = head(jx,jy,jz)
    ELSE
        COORDINATE = 'X'
        coef(1) = -2.0d0*dt*Kfacx(jx,jy,jz) / (dxx(jx)*(dxx(jx+1)+dxx(jx)))
        IF (activecellPressure(0,jy,jz) == 0) THEN
            coef(-1) = -2.0d0*dt*Kfacx(jx-1,jy,jz) / (dxx(jx)*(dxx(jx)+dxx(jx-1)))
        ELSE
            coef(-1) = 0.0d0
        END IF

        COORDINATE = 'Y'
        coef(2)  = -2.0d0*dt*Kfacy(jx,jy,jz) / (dyy(jy)*(dyy(jy+1)+dyy(jy) ))
        IF (activecellPressure(jx,0,jz) == 0) THEN
            coef(-2) = -2.0d0*dt*Kfacy(jx,jy-1,jz) / (dyy(jy)*(dyy(jy)+dyy(jy-1)))
        ELSE
            coef(-2) = 0.0d0
        END IF

        COORDINATE = 'Z'
        coef(3)  = -2.0d0*dt*Kfacz(jx,jy,jz) / (dzz(jx,jy,jz)*(dzz(jx,jy,jz+1)+dzz(jx,jy,jz) ))
        IF (activecellPressure(jx,jy,0) == 0) THEN
            coef(-3) = -2.0d0*dt*Kfacz(jx,jy,jz-1) / (dzz(jx,jy,jz)*(dzz(jx,jy,jz)+dzz(jx,jy,jz-1)))
        ELSE
            coef(-3) = 0.0d0
        END IF

        coef(0) = (Ch(jx,jy,jz)+Ss*wc(jx,jy,jz)/wcs(jx,jy,jz)) - coef(-3) - coef(-2) - coef(-1) - coef(1) - coef(2) - coef(3)

        IF (activecellPressure(0,jy,jz) == 0) THEN
          AddPressureX =  -coef(-1) * headOld(0,jy,jz)
        ELSE
          AddPressureX = 0.0d0
        END IF
        IF (activecellPressure(jx,0,jz) == 0) THEN
          AddPressureY =  -coef(-2) * headOld(jx,0,jz)
        ELSE
          AddPressureY = 0.0d0
        END IF
        IF (activecellPressure(jx,jy,0) == 0) THEN
          AddPressureZ =  -coef(-3) * headOld(jx,jy,0)
        ELSE
          AddPressureZ = 0.0d0
        END IF
        BvecCrunchP(j) = (Ch(jx,jy,jz)+Ss*wc(jx,jy,jz)/wcs(jx,jy,jz))*headOld(jx,jy,jz) - (dt/dzz(jx,jy,jz))*(Kfacz(jx,jy,jz)-Kfacz(jx,jy,jz-1)) + AddPressureX + AddPressureY + AddPressureZ
        XvecCrunchP(j) = head(jx,jy,jz)
        !!            CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j-nx*ny,coef(-3),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx*ny,coef(3),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)

    END IF

!!! ip-jm-km corner
    jx = nx
    jy = 1
    jz = 1
    j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1
    IF (activecellPressure(jx,jy,jz) == 0) THEN
        CALL MatSetValues(amatP,1,j,1,j-1,zero,INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j+1,zero,INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j-nx,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx,zero,INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j-nx*ny,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx*ny,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j,big,INSERT_VALUES,ierr)
        BvecCrunchP(j) = head(jx,jy,jz)*big
        XvecCrunchP(j) = head(jx,jy,jz)
    ELSE
        COORDINATE = 'X'
        coef(-1) = -2.0d0*dt*Kfacx(jx-1,jy,jz) / (dxx(jx)*(dxx(jx)+dxx(jx-1)))
        IF (activecellPressure(nx+1,jy,jz) == 0) THEN
            coef(1) = -2.0d0*dt*Kfacx(jx,jy,jz) / (dxx(jx)*(dxx(jx+1)+dxx(jx)))
        ELSE
            coef(1) = 0.0d0
        END IF

        COORDINATE = 'Y'
        coef(2)  = -2.0d0*dt*Kfacy(jx,jy,jz) / (dyy(jy)*(dyy(jy+1)+dyy(jy) ))
        IF (activecellPressure(jx,0,jz) == 0) THEN
            coef(-2) = -2.0d0*dt*Kfacy(jx,jy-1,jz) / (dyy(jy)*(dyy(jy)+dyy(jy-1)))
        ELSE
            coef(-2) = 0.0d0
        END IF

        COORDINATE = 'Z'
        coef(3)  = -2.0d0*dt*Kfacz(jx,jy,jz) / (dzz(jx,jy,jz)*(dzz(jx,jy,jz+1)+dzz(jx,jy,jz) ))
        IF (activecellPressure(jx,jy,0) == 0) THEN
            coef(-3) = -2.0d0*dt*Kfacz(jx,jy,jz-1) / (dzz(jx,jy,jz)*(dzz(jx,jy,jz)+dzz(jx,jy,jz-1)))
        ELSE
            coef(-3) = 0.0d0
        END IF

        coef(0) = (Ch(jx,jy,jz)+Ss*wc(jx,jy,jz)/wcs(jx,jy,jz)) - coef(-3) - coef(-2) - coef(-1) - coef(1) - coef(2) - coef(3)

        IF (activecellPressure(nx+1,jy,jz) == 0) THEN
          AddPressureX =  -coef(1) * headOld(nx+1,jy,jz)
        ELSE
          AddPressureX = 0.0d0
        END IF
        IF (activecellPressure(jx,0,jz) == 0) THEN
          AddPressureY =  -coef(-2) * headOld(jx,0,jz)
        ELSE
          AddPressureY = 0.0d0
        END IF
        IF (activecellPressure(jx,jy,0) == 0) THEN
          AddPressureZ =  -coef(-3) * headOld(jx,jy,0)
        ELSE
          AddPressureZ = 0.0d0
        END IF
        BvecCrunchP(j) = (Ch(jx,jy,jz)+Ss*wc(jx,jy,jz)/wcs(jx,jy,jz))*headOld(jx,jy,jz) - (dt/dzz(jx,jy,jz))*(Kfacz(jx,jy,jz)-Kfacz(jx,jy,jz-1)) + AddPressureX + AddPressureY + AddPressureZ
        XvecCrunchP(j) = head(jx,jy,jz)

        CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j-nx*ny,coef(-3),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx*ny,coef(3),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)
    END IF

!! im-jp-km corner
    jx = 1
    jy = ny
    jz = 1
    j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1
    IF (activecellPressure(jx,jy,jz) == 0) THEN
        !!            CALL MatSetValues(amatP,1,j,1,j-1,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+1,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx,zero,INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j+nx,zero,INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j-nx*ny,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx*ny,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j,big,INSERT_VALUES,ierr)
        BvecCrunchP(j) = head(jx,jy,jz)*big
        XvecCrunchP(j) = head(jx,jy,jz)
    ELSE
        COORDINATE = 'X'
        coef(1) = -2.0d0*dt*Kfacx(jx,jy,jz) / (dxx(jx)*(dxx(jx+1)+dxx(jx)))
        IF (activecellPressure(0,jy,jz) == 0) THEN
            coef(-1) = -2.0d0*dt*Kfacx(jx-1,jy,jz) / (dxx(jx)*(dxx(jx)+dxx(jx-1)))
        ELSE
            coef(-1) = 0.0d0
        END IF

        COORDINATE = 'Y'
        coef(-2) = -2.0d0*dt*Kfacy(jx,jy-1,jz) / (dyy(jy)*(dyy(jy)+dyy(jy-1)))
        IF (activecellPressure(jx,ny+1,jz) == 0) THEN
            coef(2)  = -2.0d0*dt*Kfacy(jx,jy,jz) / (dyy(jy)*(dyy(jy+1)+dyy(jy) ))
        ELSE
            coef(2) = 0.0d0
        END IF

        COORDINATE = 'Z'
        coef(3)  = -2.0d0*dt*Kfacz(jx,jy,jz) / (dzz(jx,jy,jz)*(dzz(jx,jy,jz+1)+dzz(jx,jy,jz) ))
        IF (activecellPressure(jx,jy,0) == 0) THEN
            coef(-3) = -2.0d0*dt*Kfacz(jx,jy,jz-1) / (dzz(jx,jy,jz)*(dzz(jx,jy,jz)+dzz(jx,jy,jz-1)))
        ELSE
            coef(-3) = 0.0d0
        END IF

        coef(0) = (Ch(jx,jy,jz)+Ss*wc(jx,jy,jz)/wcs(jx,jy,jz)) - coef(-3) - coef(-2) - coef(-1) - coef(1) - coef(2) - coef(3)

        IF (activecellPressure(0,jy,jz) == 0) THEN
          AddPressureX =  -coef(-1) * headOld(0,jy,jz)
        ELSE
          AddPressureX = 0.0d0
        END IF
        IF (activecellPressure(jx,ny+1,jz) == 0) THEN
          AddPressureY =  -coef(2) * headOld(jx,ny+1,jz)
        ELSE
          AddPressureY = 0.0d0
        END IF
        IF (activecellPressure(jx,jy,0) == 0) THEN
          AddPressureZ =  -coef(-3) * headOld(jx,jy,0)
        ELSE
          AddPressureZ = 0.0d0
        END IF
        BvecCrunchP(j) = (Ch(jx,jy,jz)+Ss*wc(jx,jy,jz)/wcs(jx,jy,jz))*headOld(jx,jy,jz) - (dt/dzz(jx,jy,jz))*(Kfacz(jx,jy,jz)-Kfacz(jx,jy,jz-1)) + AddPressureX + AddPressureY + AddPressureZ
        XvecCrunchP(j) = head(jx,jy,jz)
        !!            CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j-nx*ny,coef(-3),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx*ny,coef(3),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)

    END IF

!! ip-jp-km corner
    jx = nx
    jy = ny
    jz = 1
    j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1
    IF (activecellPressure(jx,jy,jz) == 0) THEN
        CALL MatSetValues(amatP,1,j,1,j-1,zero,INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j+1,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx,zero,INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j+nx,zero,INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j-nx*ny,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx*ny,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j,big,INSERT_VALUES,ierr)
        BvecCrunchP(j) = head(jx,jy,jz)*big
        XvecCrunchP(j) = head(jx,jy,jz)
    ELSE
        COORDINATE = 'X'
        coef(-1) = -2.0d0*dt*Kfacx(jx-1,jy,jz) / (dxx(jx)*(dxx(jx)+dxx(jx-1)))
        IF (activecellPressure(nx+1,jy,jz) == 0) THEN
            coef(1) = -2.0d0*dt*Kfacx(jx,jy,jz) / (dxx(jx)*(dxx(jx+1)+dxx(jx)))
        ELSE
            coef(1) = 0.0d0
        END IF

        COORDINATE = 'Y'
        coef(-2) = -2.0d0*dt*Kfacy(jx,jy-1,jz) / (dyy(jy)*(dyy(jy)+dyy(jy-1)))
        IF (activecellPressure(jx,ny+1,jz) == 0) THEN
            coef(2)  = -2.0d0*dt*Kfacy(jx,jy,jz) / (dyy(jy)*(dyy(jy+1)+dyy(jy) ))
        ELSE
            coef(2) = 0.0d0
        END IF

        COORDINATE = 'Z'
        coef(3)  = -2.0d0*dt*Kfacz(jx,jy,jz) / (dzz(jx,jy,jz)*(dzz(jx,jy,jz+1)+dzz(jx,jy,jz) ))
        IF (activecellPressure(jx,jy,0) == 0) THEN
            coef(-3) = -2.0d0*dt*Kfacz(jx,jy,jz-1) / (dzz(jx,jy,jz)*(dzz(jx,jy,jz)+dzz(jx,jy,jz-1)))
        ELSE
            coef(-3) = 0.0d0
        END IF

        coef(0) = (Ch(jx,jy,jz)+Ss*wc(jx,jy,jz)/wcs(jx,jy,jz)) - coef(-3) - coef(-2) - coef(-1) - coef(1) - coef(2) - coef(3)

        IF (activecellPressure(nx+1,jy,jz) == 0) THEN
          AddPressureX =  -coef(1) * headOld(nx+1,jy,jz)
        ELSE
          AddPressureX = 0.0d0
        END IF
        IF (activecellPressure(jx,ny+1,jz) == 0) THEN
          AddPressureY =  -coef(2) * headOld(jx,ny+1,jz)
        ELSE
          AddPressureY = 0.0d0
        END IF
        IF (activecellPressure(jx,jy,0) == 0) THEN
          AddPressureZ =  -coef(-3) * headOld(jx,jy,0)
        ELSE
          AddPressureZ = 0.0d0
        END IF
        BvecCrunchP(j) = (Ch(jx,jy,jz)+Ss*wc(jx,jy,jz)/wcs(jx,jy,jz))*headOld(jx,jy,jz) - (dt/dzz(jx,jy,jz))*(Kfacz(jx,jy,jz)-Kfacz(jx,jy,jz-1)) + AddPressureX + AddPressureY + AddPressureZ
        XvecCrunchP(j) = head(jx,jy,jz)

        CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j-nx*ny,coef(-3),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx*ny,coef(3),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)

    END IF

!! im-jm-kp corner
    jx = 1
    jy = 1
    jz = nz
    j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1
    IF (activecellPressure(jx,jy,jz) == 0) THEN
        !!            CALL MatSetValues(amatP,1,j,1,j-1,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+1,zero,INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j-nx,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx*ny,zero,INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j+nx*ny,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j,big,INSERT_VALUES,ierr)
        BvecCrunchP(j) = head(jx,jy,jz)*big
        XvecCrunchP(j) = head(jx,jy,jz)
    ELSE
        COORDINATE = 'X'
        coef(1) = -2.0d0*dt*Kfacx(jx,jy,jz) / (dxx(jx)*(dxx(jx+1)+dxx(jx)))
        IF (activecellPressure(0,jy,jz) == 0) THEN
            coef(-1) = -2.0d0*dt*Kfacx(jx-1,jy,jz) / (dxx(jx)*(dxx(jx)+dxx(jx-1)))
        ELSE
            coef(-1) = 0.0d0
        END IF

        COORDINATE = 'Y'
        coef(2)  = -2.0d0*dt*Kfacy(jx,jy,jz) / (dyy(jy)*(dyy(jy+1)+dyy(jy) ))
        IF (activecellPressure(jx,0,jz) == 0) THEN
            coef(-2) = -2.0d0*dt*Kfacy(jx,jy-1,jz) / (dyy(jy)*(dyy(jy)+dyy(jy-1)))
        ELSE
            coef(-2) = 0.0d0
        END IF

        COORDINATE = 'Z'
        coef(-3) = -2.0d0*dt*Kfacz(jx,jy,jz-1) / (dzz(jx,jy,jz)*(dzz(jx,jy,jz)+dzz(jx,jy,jz-1)))
        IF (activecellPressure(jx,jy,nz+1) == 0) THEN
            coef(3)  = -2.0d0*dt*Kfacz(jx,jy,jz) / (dzz(jx,jy,jz)*(dzz(jx,jy,jz+1)+dzz(jx,jy,jz) ))
        ELSE
            coef(3) = 0.0d0
        END IF

        coef(0) = (Ch(jx,jy,jz)+Ss*wc(jx,jy,jz)/wcs(jx,jy,jz)) - coef(-3) - coef(-2) - coef(-1) - coef(1) - coef(2) - coef(3)

        IF (activecellPressure(0,jy,jz) == 0) THEN
          AddPressureX =  -coef(-1) * headOld(0,jy,jz)
        ELSE
          AddPressureX = 0.0d0
        END IF
        IF (activecellPressure(jx,0,jz) == 0) THEN
          AddPressureY =  -coef(-2) * headOld(jx,0,jz)
        ELSE
          AddPressureY = 0.0d0
        END IF
        IF (activecellPressure(jx,jy,nz+1) == 0) THEN
          AddPressureZ =  -coef(3) * headOld(jx,jy,nz+1)
        ELSE
          AddPressureZ = 0.0d0
        END IF
        BvecCrunchP(j) = (Ch(jx,jy,jz)+Ss*wc(jx,jy,jz)/wcs(jx,jy,jz))*headOld(jx,jy,jz) - (dt/dzz(jx,jy,jz))*(Kfacz(jx,jy,jz)-Kfacz(jx,jy,jz-1)) + AddPressureX + AddPressureY + AddPressureZ
        XvecCrunchP(j) = head(jx,jy,jz)

        !!            CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx*ny,coef(-3),INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j+nx*ny,coef(3),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)

    END IF

!! ip-jm-kp corner
    jx = nx
    jy = 1
    jz = nz
    j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1
    IF (activecellPressure(jx,jy,jz) == 0) THEN
        CALL MatSetValues(amatP,1,j,1,j-1,zero,INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j+1,zero,INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j-nx,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx*ny,zero,INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j+nx*ny,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j,big,INSERT_VALUES,ierr)
        BvecCrunchP(j) = head(jx,jy,jz)*big
        XvecCrunchP(j) = head(jx,jy,jz)
    ELSE
        COORDINATE = 'X'
        coef(-1) = -2.0d0*dt*Kfacx(jx-1,jy,jz) / (dxx(jx)*(dxx(jx)+dxx(jx-1)))
        IF (activecellPressure(nx+1,jy,jz) == 0) THEN
            coef(1) = -2.0d0*dt*Kfacx(jx,jy,jz) / (dxx(jx)*(dxx(jx+1)+dxx(jx)))
        ELSE
            coef(1) = 0.0d0
        END IF

        COORDINATE = 'Y'
        coef(2)  = -2.0d0*dt*Kfacy(jx,jy,jz) / (dyy(jy)*(dyy(jy+1)+dyy(jy) ))
        IF (activecellPressure(jx,0,jz) == 0) THEN
            coef(-2) = -2.0d0*dt*Kfacy(jx,jy-1,jz) / (dyy(jy)*(dyy(jy)+dyy(jy-1)))
        ELSE
            coef(-2) = 0.0d0
        END IF

        COORDINATE = 'Z'
        coef(-3) = -2.0d0*dt*Kfacz(jx,jy,jz-1) / (dzz(jx,jy,jz)*(dzz(jx,jy,jz)+dzz(jx,jy,jz-1)))
        IF (activecellPressure(jx,jy,nz+1) == 0) THEN
            coef(3)  = -2.0d0*dt*Kfacz(jx,jy,jz) / (dzz(jx,jy,jz)*(dzz(jx,jy,jz+1)+dzz(jx,jy,jz) ))
        ELSE
            coef(3) = 0.0d0
        END IF

        coef(0) = (Ch(jx,jy,jz)+Ss*wc(jx,jy,jz)/wcs(jx,jy,jz)) - coef(-3) - coef(-2) - coef(-1) - coef(1) - coef(2) - coef(3)

        IF (activecellPressure(nx+1,jy,jz) == 0) THEN
          AddPressureX =  -coef(1) * headOld(nx+1,jy,jz)
        ELSE
          AddPressureX = 0.0d0
        END IF
        IF (activecellPressure(jx,0,jz) == 0) THEN
          AddPressureY =  -coef(-2) * headOld(jx,0,jz)
        ELSE
          AddPressureY = 0.0d0
        END IF
        IF (activecellPressure(jx,jy,nz+1) == 0) THEN
          AddPressureZ =  -coef(3) * headOld(jx,jy,nz+1)
        ELSE
          AddPressureZ = 0.0d0
        END IF
        BvecCrunchP(j) = (Ch(jx,jy,jz)+Ss*wc(jx,jy,jz)/wcs(jx,jy,jz))*headOld(jx,jy,jz) - (dt/dzz(jx,jy,jz))*(Kfacz(jx,jy,jz)-Kfacz(jx,jy,jz-1)) + AddPressureX + AddPressureY + AddPressureZ
        XvecCrunchP(j) = head(jx,jy,jz)

        CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx*ny,coef(-3),INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j+nx*ny,coef(3),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)

    END IF

!! im-jp-kp corner
    jx = 1
    jy = ny
    jz = nz
    j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1
    IF (activecellPressure(jx,jy,jz) == 0) THEN
        !!            CALL MatSetValues(amatP,1,j,1,j-1,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+1,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx,zero,INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j+nx,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx*ny,zero,INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j+nx*ny,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j,big,INSERT_VALUES,ierr)
        BvecCrunchP(j) = head(jx,jy,jz)*big
        XvecCrunchP(j) = head(jx,jy,jz)
    ELSE
        COORDINATE = 'X'
        coef(1) = -2.0d0*dt*Kfacx(jx,jy,jz) / (dxx(jx)*(dxx(jx+1)+dxx(jx)))
        IF (activecellPressure(0,jy,jz) == 0) THEN
            coef(-1) = -2.0d0*dt*Kfacx(jx-1,jy,jz) / (dxx(jx)*(dxx(jx)+dxx(jx-1)))
        ELSE
            coef(-1) = 0.0d0
        END IF

        COORDINATE = 'Y'
        coef(-2) = -2.0d0*dt*Kfacy(jx,jy-1,jz) / (dyy(jy)*(dyy(jy)+dyy(jy-1)))
        IF (activecellPressure(jx,ny+1,jz) == 0) THEN
            coef(2)  = -2.0d0*dt*Kfacy(jx,jy,jz) / (dyy(jy)*(dyy(jy+1)+dyy(jy) ))
        ELSE
            coef(2) = 0.0d0
        END IF

        COORDINATE = 'Z'
        coef(-3) = -2.0d0*dt*Kfacz(jx,jy,jz-1) / (dzz(jx,jy,jz)*(dzz(jx,jy,jz)+dzz(jx,jy,jz-1)))
        IF (activecellPressure(jx,jy,nz+1) == 0) THEN
            coef(3)  = -2.0d0*dt*Kfacz(jx,jy,jz) / (dzz(jx,jy,jz)*(dzz(jx,jy,jz+1)+dzz(jx,jy,jz) ))
        ELSE
            coef(3) = 0.0d0
        END IF

        coef(0) = (Ch(jx,jy,jz)+Ss*wc(jx,jy,jz)/wcs(jx,jy,jz)) - coef(-3) - coef(-2) - coef(-1) - coef(1) - coef(2) - coef(3)

        IF (activecellPressure(0,jy,jz) == 0) THEN
          AddPressureX =  -coef(-1) * headOld(0,jy,jz)
        ELSE
          AddPressureX = 0.0d0
        END IF
        IF (activecellPressure(jx,ny+1,jz) == 0) THEN
          AddPressureY =  -coef(2) * headOld(jx,ny+1,jz)
        ELSE
          AddPressureY = 0.0d0
        END IF
        IF (activecellPressure(jx,jy,nz+1) == 0) THEN
          AddPressureZ =  -coef(3) * headOld(jx,jy,nz+1)
        ELSE
          AddPressureZ = 0.0d0
        END IF
        BvecCrunchP(j) = (Ch(jx,jy,jz)+Ss*wc(jx,jy,jz)/wcs(jx,jy,jz))*headOld(jx,jy,jz) - (dt/dzz(jx,jy,jz))*(Kfacz(jx,jy,jz)-Kfacz(jx,jy,jz-1)) + AddPressureX + AddPressureY + AddPressureZ
        XvecCrunchP(j) = head(jx,jy,jz)
        !!            CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx*ny,coef(-3),INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j+nx*ny,coef(3),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)
    END IF

!! ip-jp-kp corner
    jx = nx
    jy = ny
    jz = nz
    j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1
    IF (activecellPressure(jx,jy,jz) == 0) THEN
        CALL MatSetValues(amatP,1,j,1,j-1,zero,INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j+1,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx,zero,INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j+nx,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx*ny,zero,INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j+nx*ny,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j,big,INSERT_VALUES,ierr)
        BvecCrunchP(j) = head(jx,jy,jz)*big
        XvecCrunchP(j) = head(jx,jy,jz)
    ELSE
        COORDINATE = 'X'
        coef(-1) = -2.0d0*dt*Kfacx(jx-1,jy,jz) / (dxx(jx)*(dxx(jx)+dxx(jx-1)))
        IF (activecellPressure(nx+1,jy,jz) == 0) THEN
            coef(1) = -2.0d0*dt*Kfacx(jx,jy,jz) / (dxx(jx)*(dxx(jx+1)+dxx(jx)))
        ELSE
            coef(1) = 0.0d0
        END IF

        COORDINATE = 'Y'
        coef(-2) = -2.0d0*dt*Kfacy(jx,jy-1,jz) / (dyy(jy)*(dyy(jy)+dyy(jy-1)))
        IF (activecellPressure(jx,ny+1,jz) == 0) THEN
            coef(2)  = -2.0d0*dt*Kfacy(jx,jy,jz) / (dyy(jy)*(dyy(jy+1)+dyy(jy) ))
        ELSE
            coef(2) = 0.0d0
        END IF

        COORDINATE = 'Z'
        coef(-3) = -2.0d0*dt*Kfacz(jx,jy,jz-1) / (dzz(jx,jy,jz)*(dzz(jx,jy,jz)+dzz(jx,jy,jz-1)))
        IF (activecellPressure(jx,jy,nz+1) == 0) THEN
            coef(3)  = -2.0d0*dt*Kfacz(jx,jy,jz) / (dzz(jx,jy,jz)*(dzz(jx,jy,jz+1)+dzz(jx,jy,jz) ))
        ELSE
            coef(3) = 0.0d0
        END IF

        coef(0) = (Ch(jx,jy,jz)+Ss*wc(jx,jy,jz)/wcs(jx,jy,jz)) - coef(-3) - coef(-2) - coef(-1) - coef(1) - coef(2) - coef(3)

        IF (activecellPressure(nx+1,jy,jz) == 0) THEN
          AddPressureX =  -coef(1) * headOld(nx+1,jy,jz)
        ELSE
          AddPressureX = 0.0d0
        END IF
        IF (activecellPressure(jx,ny+1,jz) == 0) THEN
          AddPressureY =  -coef(2) * headOld(jx,ny+1,jz)
        ELSE
          AddPressureY = 0.0d0
        END IF
        IF (activecellPressure(jx,jy,nz+1) == 0) THEN
          AddPressureZ =  -coef(3) * headOld(jx,jy,nz+1)
        ELSE
          AddPressureZ = 0.0d0
        END IF
        BvecCrunchP(j) = (Ch(jx,jy,jz)+Ss*wc(jx,jy,jz)/wcs(jx,jy,jz))*headOld(jx,jy,jz) - (dt/dzz(jx,jy,jz))*(Kfacz(jx,jy,jz)-Kfacz(jx,jy,jz-1)) + AddPressureX + AddPressureY + AddPressureZ
        XvecCrunchP(j) = head(jx,jy,jz)

        CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx*ny,coef(-3),INSERT_VALUES,ierr)
        !!            CALL MatSetValues(amatP,1,j,1,j+nx*ny,coef(3),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)
    END IF
DO j = 0,nx*ny-1
    IF (DABS(BvecCrunchP(j)) < 1e-40) THEN
        BvecCrunchP(j) = 0.0d0
    END IF
END DO

! DO jz = 1,nz
!     DO jy = 2,ny-1
!         DO jx = 2,nx-1
!             IF (jz == nz) THEN
!                 j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1
!                 WRITE(*,*) '>>>>>2 jz, B = ',jz,BvecCrunchP(j)
!             END IF
!         END DO
!     END DO
! END DO

!*********************end Matrix Right Hand Side ******************************

CALL MatAssemblyBegin(amatP,MAT_FINAL_ASSEMBLY,ierr)
CALL MatAssemblyEnd(amatP,MAT_FINAL_ASSEMBLY,ierr)

!----------------------------------------------------------------------------------
RETURN
END SUBROUTINE pressureRich
