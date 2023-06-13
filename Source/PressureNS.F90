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

SUBROUTINE pressureNS (nx,ny,nz,dtyr,amatP,SteadyFlow)
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
REAL(DP)                                                              :: permchar
REAL(DP)                                                              :: pumpterm
REAL(DP)                                                              :: tdepend
REAL(DP)                                                              :: ct
REAL(DP)                                                              :: RoAveLeft
REAL(DP)                                                              :: RoAveRight
REAL(DP)                                                              :: Check
REAL(DP)                                                              :: AddPressureX
REAL(DP)                                                              :: AddPressureY
REAL(DP)                                                              :: uxp
REAL(DP)                                                              :: uxm
REAL(DP)                                                              :: uyp
REAL(DP)                                                              :: uym
REAL(DP)                                                              :: vxp
REAL(DP)                                                              :: vxm
REAL(DP)                                                              :: vyp
REAL(DP)                                                              :: vym
REAL(DP)                                                              :: uppx
REAL(DP)                                                              :: ummx
REAL(DP)                                                              :: uppy
REAL(DP)                                                              :: ummy
REAL(DP)                                                              :: vpp
REAL(DP)                                                              :: vmm
REAL(DP)                                                              :: vppx
REAL(DP)                                                              :: vmmx
REAL(DP)                                                              :: vppy
REAL(DP)                                                              :: vmmy
REAL(DP)                                                              :: upp
REAL(DP)                                                              :: umm


!  ****** PARAMETERS  ****************************

REAL(DP), PARAMETER                                                   :: visc=0.000001d0
!REAL(DP), PARAMETER                                                   :: visc=31.54d0
REAL(DP), PARAMETER                                                   :: ctTransient=4.8D-07
REAL(DP), PARAMETER                                                   :: ctSteady=0.00
REAL(DP), PARAMETER                                                   :: big=100000.0D0
REAL(DP), PARAMETER                                                   :: zero=0.0d0
REAL(DP), PARAMETER                                                   :: grav=9.806d0
REAL(DP), PARAMETER                                                   :: alphaUnsat=0.001d0
REAL(DP), PARAMETER                                                   :: alphaSat=0.000d0

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

!! Only work for 2D XY domain
jz = 1

!! *********************** Calculate predictor velocties ***************************

DO jy = 0,ny-1
    DO jx = 0,nx
        qx(jx,jy,jz) = qx(jx,jy,jz) / (365 * 86400)
    END DO
END DO
DO jy = 0,ny
    DO jx = 0,nx-1
        qy(jx,jy,jz) = qy(jx,jy,jz) / (365 * 86400)
    END DO
END DO

! -------------------- Calculate us
DO jx = 0,nx
    DO jy = 0,ny-1
        IF (jx == 0) THEN
            uxp = qx(jx+1,jy,jz)
            uxm = qx(jx,jy,jz)
            vpp = 0.0d0
            vmm = 0.0d0
            IF (jy == 0) THEN
                uyp = qx(jx,jy+1,jz)
                uym = -qx(jx,jy,jz)
            ELSE IF (jy == ny-1) THEN
                uyp = -qx(jx,jy,jz)
                uym = qx(jx,jy-1,jz)
            ELSE
                uyp = qx(jx,jy+1,jz)
                uym = qx(jx,jy-1,jz)
            END IF
        ELSE IF (jx == nx) THEN
            uxp = qx(jx,jy,jz)
            uxm = qx(jx-1,jy,jz)
            vpp = 0.0d0
            vmm = 0.0d0
            IF (jy == 0) THEN
                uyp = qx(jx,jy+1,jz)
                uym = -qx(jx,jy,jz)
            ELSE IF (jy == ny-1) THEN
                uyp = -qx(jx,jy,jz)
                uym = qx(jx,jy-1,jz)
            ELSE
                uyp = qx(jx,jy+1,jz)
                uym = qx(jx,jy-1,jz)
            END IF
        ELSE
            uxp = qx(jx+1,jy,jz)
            uxm = qx(jx-1,jy,jz)
            vpp = ArithmeticMean(qy(jx,jy+1,jz),qy(jx+1,jy+1,jz))
            vmm = ArithmeticMean(qy(jx,jy,jz),qy(jx+1,jy,jz))
            IF (jy == 0) THEN
                uyp = qx(jx,jy+1,jz)
                uym = -qx(jx,jy,jz)
            ELSE IF (jy == ny-1) THEN
                uyp = -qx(jx,jy,jz)
                uym = qx(jx,jy-1,jz)
            ELSE
                uyp = qx(jx,jy+1,jz)
                uym = qx(jx,jy-1,jz)
            END IF
        END IF
        uppx = ArithmeticMean(qx(jx,jy,jz),uxp)
        ummx = ArithmeticMean(qx(jx,jy,jz),uxm)
        uppy = ArithmeticMean(qx(jx,jy,jz),uyp)
        ummy = ArithmeticMean(qx(jx,jy,jz),uym)
        ! check if face has zero permeability
!        IF (permx(jx,jy,jz) == 0.0d0) THEN
        us(jx,jy,jz) = qx(jx,jy,jz) - (dt/dxx(jx))*(uppx*uppx-ummx*ummx) - (dt/dyy(jy))*(uppy*vpp-ummy*vmm) + &
            (visc*dt/(dxx(jx)*dxx(jx)))*(uxp - 2.0d0*qx(jx,jy,jz) + uxm) + &
            (visc*dt/(dyy(jy)*dyy(jy)))*(uyp - 2.0d0*qx(jx,jy,jz) + uym)
!        us(jx,jy,jz) = qx(jx,jy,jz) + &
!            (visc*dt/(dxx(jx)*dxx(jx)))*(uxp - 2.0d0*qx(jx,jy,jz) + uxm) + &
!            (visc*dt/(dyy(jy)*dyy(jy)))*(uyp - 2.0d0*qx(jx,jy,jz) + uym)
!        us(jx,jy,jz) = qx(jx,jy,jz)
!        WRITE(*,*) ' >>> US = ',jx,jy,us(jx,jy,jz)
    END DO
END DO

! -------------------- Calculate vs
DO jx = 0,nx-1
    DO jy = 0,ny
        IF (jy == 0) THEN
            vyp = qy(jx,jy+1,jz)
            vym = qy(jx,jy,jz)
            upp = 0.0d0
            umm = 0.0d0
            IF (jx == 0) THEN
                vxp = qy(jx+1,jy,jz)
                vxm = -qy(jx,jy,jz)
            ELSE IF (jx == nx-1) THEN
                vxp = -qy(jx,jy,jz)
                vxm = qy(jx-1,jy,jz)
            ELSE
                vxp = qy(jx+1,jy,jz)
                vxm = qy(jx-1,jy,jz)
            END IF
        ELSE IF (jy == ny) THEN
            vyp = qy(jx,jy,jz)
            vym = qy(jx,jy-1,jz)
            upp = 0.0d0
            umm = 0.0d0
            IF (jx == 0) THEN
                vxp = qy(jx+1,jy,jz)
                vxm = -qy(jx,jy,jz)
            ELSE IF (jx == nx-1) THEN
                vxp = -qy(jx,jy,jz)
                vxm = qy(jx-1,jy,jz)
            ELSE
                vxp = qy(jx+1,jy,jz)
                vxm = qy(jx-1,jy,jz)
            END IF
        ELSE
            vyp = qy(jx,jy+1,jz)
            vym = qy(jx,jy-1,jz)
            upp = ArithmeticMean(qx(jx+1,jy,jz),qx(jx+1,jy+1,jz))
            umm = ArithmeticMean(qx(jx,jy,jz),qx(jx,jy+1,jz))
            IF (jx == 0) THEN
                vxp = qy(jx+1,jy,jz)
                vxm = -qy(jx,jy,jz)
            ELSE IF (jx == nx-1) THEN
                vxp = -qy(jx,jy,jz)
                vxm = qy(jx-1,jy,jz)
            ELSE
                vxp = qy(jx+1,jy,jz)
                vxm = qy(jx-1,jy,jz)
            END IF
        END IF
        vppx = ArithmeticMean(qy(jx,jy,jz),vxp)
        vmmx = ArithmeticMean(qy(jx,jy,jz),vxm)
        vppy = ArithmeticMean(qy(jx,jy,jz),vyp)
        vmmy = ArithmeticMean(qy(jx,jy,jz),vym)
        vs(jx,jy,jz) = qy(jx,jy,jz) - (dt/dxx(jx))*(upp*vppx-umm*vmmx) - (dt/dyy(jy))*(vppy*vppy-vmmy*vmmy) + &
            (visc*dt/(dxx(jx)*dxx(jx)))*(vxp - 2.0d0*qy(jx,jy,jz) + vxm) + &
            (visc*dt/(dyy(jy)*dyy(jy)))*(vyp - 2.0d0*qy(jx,jy,jz) + vym)
!        vs(jx,jy,jz) = qy(jx,jy,jz) + &
!            (visc*dt/(dxx(jx)*dxx(jx)))*(vxp - 2.0d0*qy(jx,jy,jz) + vxm) + &
!            (visc*dt/(dyy(jy)*dyy(jy)))*(vyp - 2.0d0*qy(jx,jy,jz) + vym)
!        vs(jx,jy,jz) = qy(jx,jy,jz)
!        WRITE(*,*) ' >>> VS = ',jx,jy,vs(jx,jy,jz)
    END DO
END DO
                


!*********************begin Matrix Coefficients ******************************
!! --- Matrix coefficients for the interior of domain
DO jy = 2,ny-1
  DO jx = 2,nx-1
    j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1
    IF (activecellPressure(jx,jy,jz) == 0) THEN
      CALL MatSetValues(amatP,1,j,1,j-1,zero,INSERT_VALUES,ierr)
      CALL MatSetValues(amatP,1,j,1,j+1,zero,INSERT_VALUES,ierr)
      CALL MatSetValues(amatP,1,j,1,j-nx,zero,INSERT_VALUES,ierr)
      CALL MatSetValues(amatP,1,j,1,j+nx,zero,INSERT_VALUES,ierr)
      CALL MatSetValues(amatP,1,j,1,j,big,INSERT_VALUES,ierr)
    ELSE
      COORDINATE = 'X'
      coef(-1) = 2.0d0 / (dxx(jx)*(dxx(jx)+dxx(jx-1)))
      coef(1) = 2.0d0 / (dxx(jx)*(dxx(jx+1)+dxx(jx)))

      COORDINATE = 'Y'
      coef(-2) = 2.0d0 / (dyy(jy)*(dyy(jy)+dyy(jy-1)))
      coef(2)  = 2.0d0 / (dyy(jy)*(dyy(jy+1)+dyy(jy) ))

      coef(0) = - ( coef(-2) + coef(-1) + coef(1) + coef(2) )

      CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
      CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
      CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
      CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
      CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)
    END IF
       continue
  END DO
END DO
!! --- Matrix coefficients for jy=1 boundary
jy = 1
DO jx = 2,nx-1
    j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1
    IF (activecellPressure(jx,jy,jz) == 0) THEN
        CALL MatSetValues(amatP,1,j,1,j-1,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+1,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j,big,INSERT_VALUES,ierr)
    ELSE
        COORDINATE = 'X'
        coef(-1) = 2.0d0/(dxx(jx)*(dxx(jx)+dxx(jx-1)))
        coef(1) = 2.0d0/(dxx(jx)*(dxx(jx+1)+dxx(jx)))
        IF (activecellPressure(jx,0,jz) == 0) THEN   !! Ghost cell with fixed pressure
            COORDINATE = 'Y'
            coef(-2) = 1.0d0/(dyy(jy)*dyy(jy))
            coef(2)  = 2.0d0/(dyy(jy) *(dyy(jy+1)+dyy(jy) ))
        ELSE
            COORDINATE = 'Y'
            coef(2)  = 2.0d0/(dyy(jy) *(dyy(jy+1)+dyy(jy) ))
            coef(-2) = 0.0d0
        END IF

        coef(0) = - ( coef(-2) + coef(-1) + coef(1) + coef(2) )

        CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)
    END IF
END DO

!! --- Matrix coefficients for jy=ny boundary
jy = ny
DO jx = 2,nx-1
    j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1
    IF (activecellPressure(jx,jy,jz) == 0) THEN
        CALL MatSetValues(amatP,1,j,1,j-1,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+1,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j,big,INSERT_VALUES,ierr)
    ELSE
        COORDINATE = 'X'
        coef(-1) = 2.0d0/(dxx(jx)*(dxx(jx)+dxx(jx-1)))
        coef(1) = 2.0d0/(dxx(jx)*(dxx(jx+1)+dxx(jx)))
        IF (activecellPressure(jx,ny+1,jz) == 0) THEN
            COORDINATE = 'Y'
            coef(-2) = 2.0d0/(dyy(jy)*(dyy(jy)+dyy(jy-1)))
            coef(2)  = 1.0d0/(dyy(jy)*dyy(jy))
        ELSE
            COORDINATE = 'Y'
            coef(-2) = 2.0d0/(dyy(jy)*(dyy(jy)+dyy(jy-1)))
            coef(2) = 0.0d0
        END IF
        coef(0) = - ( coef(-2) + coef(-1) + coef(1) + coef(2) )

        CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)
    END IF
END DO


!! --- Matrix coefficients for jx=1 boundary
jx = 1
DO jy = 2,ny-1
    j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1
    IF (activecellPressure(jx,jy,jz) == 0) THEN
        CALL MatSetValues(amatP,1,j,1,j+1,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx,zero,INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j,big,INSERT_VALUES,ierr)
    ELSE
        COORDINATE = 'Y'
        coef(-2) = 2.0d0/(dyy(jy)*(dyy(jy)+dyy(jy-1)))
        coef(2)  = 2.0d0/(dyy(jy) *(dyy(jy+1)+dyy(jy) ))
        IF (activecellPressure(0,jy,jz) == 0) THEN
            COORDINATE = 'X'
            coef(-1) = 1.0d0/(dxx(jx)*dxx(jx))
            coef(1) = 2.0d0/(dxx(jx)*(dxx(jx+1)+dxx(jx)))
        ELSE
            COORDINATE = 'X'
            coef(1) = 2.0d0/(dxx(jx)*(dxx(jx+1)+dxx(jx)))
            coef(-1) = 0.0d0
        END IF
        coef(0) = - ( coef(-2) + coef(-1) + coef(1) + coef(2) )

        CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
        CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)
    END IF
END DO


!! --- Matrix coefficients for jx=nx boundary
jx = nx
DO jy = 2,ny-1
    j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1
    IF (activecellPressure(jx,jy,jz) == 0) THEN
      CALL MatSetValues(amatP,1,j,1,j-1,zero,INSERT_VALUES,ierr)
      CALL MatSetValues(amatP,1,j,1,j-nx,zero,INSERT_VALUES,ierr)
      CALL MatSetValues(amatP,1,j,1,j+nx,zero,INSERT_VALUES,ierr)
      CALL MatSetValues(amatP,1,j,1,j,big,INSERT_VALUES,ierr)
    ELSE
      COORDINATE = 'Y'
      coef(-2) = 2.0d0/(dyy(jy)*(dyy(jy)+dyy(jy-1)))
      coef(2)  = 2.0d0/(dyy(jy) *(dyy(jy+1)+dyy(jy) ))

      IF (activecellPressure(nx+1,jy,jz) == 0) THEN
        COORDINATE = 'X'
        coef(-1) = 2.0d0/(dxx(jx)*(dxx(jx)+dxx(jx-1)))
        coef(1) = 1.0d0/(dxx(jx)*dxx(jx))
      ELSE
        COORDINATE = 'X'
        coef(-1) = 2.0d0/(dxx(jx)*(dxx(jx)+dxx(jx-1)))
        coef(1) = 0.0d0
      END IF

      coef(0) = - ( coef(-2) + coef(-1) + coef(1) + coef(2) )

      CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
      CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
      CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
      CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)
    END IF
END DO

!! --- Matrix coefficient for jx=1, jy=1
jx = 1
jy = 1
j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1
IF (activecellPressure(jx,jy,jz) == 0) THEN
  CALL MatSetValues(amatP,1,j,1,j+1,zero,INSERT_VALUES,ierr)
  CALL MatSetValues(amatP,1,j,1,j+nx,zero,INSERT_VALUES,ierr)
  CALL MatSetValues(amatP,1,j,1,j,big,INSERT_VALUES,ierr)
ELSE
  IF (activecellPressure(0,jy,jz) == 0) THEN
    COORDINATE = 'X'
    coef(-1) = 1.0d0/(dxx(jx)*dxx(jx))
    coef(1) = 2.0d0/(dxx(jx)*(dxx(jx+1)+dxx(jx)))
  ELSE
    COORDINATE = 'X'
    coef(1) = 2.0d0/(dxx(jx)*(dxx(jx+1)+dxx(jx)))
    coef(-1) = 0.0d0
  END IF
  IF (activecellPressure(jx,0,jz) == 0) THEN
    COORDINATE = 'Y'
    coef(2)  = 2.0d0/(dyy(jy) *(dyy(jy+1)+dyy(jy) ))
    coef(-2) = 1.0d0/(dyy(jy)*dyy(jy))
  ELSE
    COORDINATE = 'Y'
    coef(2)  = 2.0d0/(dyy(jy) *(dyy(jy+1)+dyy(jy) ))
    coef(-2) = 0.0d0
  END IF
  coef(0) = - ( coef(-2) + coef(-1) + coef(1) + coef(2) )
  CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
  CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
  CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)
END IF

!! --- Matrix coefficient for jx=nx, jy=1
jx = nx
jy = 1
j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1
IF (activecellPressure(jx,jy,jz) == 0) THEN
    CALL MatSetValues(amatP,1,j,1,j-1,zero,INSERT_VALUES,ierr)
    CALL MatSetValues(amatP,1,j,1,j+nx,zero,INSERT_VALUES,ierr)
    CALL MatSetValues(amatP,1,j,1,j,big,INSERT_VALUES,ierr)
ELSE
    IF (activecellPressure(nx+1,jy,jz) == 0) THEN
        COORDINATE = 'X'
        coef(-1) = 2.0d0/(dxx(jx)*(dxx(jx)+dxx(jx-1)))
        coef(1) = 1.0d0/(dxx(jx)*dxx(jx))
    ELSE
        COORDINATE = 'X'
        coef(-1) = 2.0d0/(dxx(jx)*(dxx(jx)+dxx(jx-1)))
        coef(1) = 0.0d0
    END IF
    IF (activecellPressure(jx,0,jz) == 0) THEN
        COORDINATE = 'Y'
        coef(2)  = 2.0d0/(dyy(jy) *(dyy(jy+1)+dyy(jy) ))
        coef(-2) = 1.0d0/(dyy(jy)*dyy(jy))
    ELSE
        COORDINATE = 'Y'
        coef(2)  = 2.0d0/(dyy(jy) *(dyy(jy+1)+dyy(jy) ))
        coef(-2) = 0.0d0
    END IF
    coef(0) = - ( coef(-2) + coef(-1) + coef(1) + coef(2) )
    CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
    CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
    CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)
END IF

!! --- Matrix coefficient for jx=1, jy=ny
jx = 1
jy = ny
j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1
IF (activecellPressure(jx,jy,jz) == 0) THEN
    CALL MatSetValues(amatP,1,j,1,j+1,zero,INSERT_VALUES,ierr)
    CALL MatSetValues(amatP,1,j,1,j-nx,zero,INSERT_VALUES,ierr)
    CALL MatSetValues(amatP,1,j,1,j,big,INSERT_VALUES,ierr)
ELSE
    IF (activecellPressure(0,jy,jz) == 0) THEN
        COORDINATE = 'X'
        coef(-1) = 1.0d0/(dxx(jx)*dxx(jx))
        coef(1) = 2.0d0/(dxx(jx)*(dxx(jx+1)+dxx(jx)))
    ELSE
        COORDINATE = 'X'
        coef(1) = 2.0d0/(dxx(jx)*(dxx(jx+1)+dxx(jx)))
        coef(-1) = 0.0d0
    END IF
    IF (activecellPressure(jx,ny+1,jz) == 0) THEN
        COORDINATE = 'Y'
        coef(-2) = 2.0d0/(dyy(jy)*(dyy(jy)+dyy(jy-1)))
        coef(2)  = 1.0d0/(dyy(jy)*dyy(jy))
    ELSE
        COORDINATE = 'Y'
        coef(-2) = 2.0d0/(dyy(jy)*(dyy(jy)+dyy(jy-1)))
        coef(2) = 0.0d0
    END IF
    coef(0) = - ( coef(-2) + coef(-1) + coef(1) + coef(2) )
    CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
    CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
    CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)
END IF

!! --- Matrix coefficient for jx=nx, jy=ny
jx = nx
jy = ny
j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1
IF (activecellPressure(jx,jy,jz) == 0) THEN
    CALL MatSetValues(amatP,1,j,1,j-1,zero,INSERT_VALUES,ierr)
    CALL MatSetValues(amatP,1,j,1,j-nx,zero,INSERT_VALUES,ierr)
    CALL MatSetValues(amatP,1,j,1,j,big,INSERT_VALUES,ierr)
ELSE
    IF (activecellPressure(nx+1,jy,jz) == 0) THEN
        COORDINATE = 'X'
        coef(-1) = 2.0d0/(dxx(jx)*(dxx(jx)+dxx(jx-1)))
        coef(1) = 1.0d0/(dxx(jx)*dxx(jx))
    ELSE
        COORDINATE = 'X'
        coef(-1) = 2.0d0/(dxx(jx)*(dxx(jx)+dxx(jx-1)))
        coef(1) = 0.0d0
    END IF
    IF (activecellPressure(jx,ny+1,jz) == 0) THEN
        COORDINATE = 'Y'
        coef(-2) = 2.0d0/(dyy(jy)*(dyy(jy)+dyy(jy-1)))
        coef(2)  = 1.0d0/(dyy(jy)*dyy(jy))
    ELSE
        COORDINATE = 'Y'
        coef(-2) = 2.0d0/(dyy(jy)*(dyy(jy)+dyy(jy-1)))
        coef(2) = 0.0d0
    END IF
    coef(0) = - ( coef(-2) + coef(-1) + coef(1) + coef(2) )
    CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
    CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
    CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)
END IF
!*********************end Matrix Coefficients ******************************

!*********************begin Matrix Right Hand Side ******************************
!! Differentiate us and vs to get right-hand-side
DO jy = 2,ny-1
    DO jx = 2,nx-1
        j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1
        IF (activecellPressure(jx,jy,jz) == 0) THEN
         BvecCrunchP(j) = pres(jx,jy,jz)*big
         XvecCrunchP(j) = pres(jx,jy,jz)
        ELSE
         BvecCrunchP(j) = (ro(jx,jy,jz)/dt) * ((us(jx,jy-1,jz)-us(jx-1,jy-1,jz))/dxx(jx) + (vs(jx-1,jy,jz)-vs(jx-1,jy-1,jz))/dyy(jy))
         XvecCrunchP(j) = pres(jx,jy,jz)
!            WRITE(*,*) 'jx,jy,us,vs = ',jx,jy,us(jx,jy-1,jz)-us(jx-1,jy-1,jz),vs(jx-1,jy,jz)-vs(jx-1,jy-1,jz)
        END IF
          continue
    END DO
END DO

jy = 1
DO jx = 2,nx-1
    j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1
    IF (activecellPressure(jx,jy,jz) == 0) THEN
        BvecCrunchP(j) = pres(jx,jy,jz)*big
        XvecCrunchP(j) = pres(jx,jy,jz)
    ELSE
        IF (activecellPressure(jx,0,jz) == 0) THEN
            AddPressureY =  -pres(jx,0,jz)/(dyy(jy)*dyy(jy))
        ELSE
            AddPressureY = 0.0d0
        END IF
        BvecCrunchP(j) = (ro(jx,jy,jz)/dt) * ((us(jx,jy-1,jz)-us(jx-1,jy-1,jz))/dxx(jx) + (vs(jx-1,jy,jz)-vs(jx-1,jy-1,jz))/dyy(jy)) + AddPressureY
        XvecCrunchP(j) = pres(jx,jy,jz)
    END IF
END DO

jy = ny
DO jx = 2,nx-1
    j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1
    IF (activecellPressure(jx,jy,jz) == 0) THEN
        BvecCrunchP(j) = pres(jx,jy,jz)*big
        XvecCrunchP(j) = pres(jx,jy,jz)
    ELSE
        IF (activecellPressure(jx,ny+1,jz) == 0) THEN
          AddPressureY =  -pres(jx,ny+1,jz)/(dyy(jy)*dyy(jy))
        ELSE
          AddPressureY = 0.0d0
        END IF
        BvecCrunchP(j) = (ro(jx,jy,jz)/dt) * ((us(jx,jy-1,jz)-us(jx-1,jy-1,jz))/dxx(jx) + (vs(jx-1,jy,jz)-vs(jx-1,jy-1,jz))/dyy(jy)) + AddPressureY
        XvecCrunchP(j) = pres(jx,jy,jz)
    END IF
END DO

jx = 1
DO jy = 2,ny-1
    j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1
    IF (activecellPressure(jx,jy,jz) == 0) THEN
        BvecCrunchP(j) = pres(jx,jy,jz)*big
        XvecCrunchP(j) = pres(jx,jy,jz)
    ELSE
        IF (activecellPressure(0,jy,jz) == 0) THEN
            AddPressureX = -pres(0,jy,jz)/(dxx(jx)*dxx(jx))
        ELSE
            AddPressureX = 0.0d0
        END IF
        BvecCrunchP(j) = (ro(jx,jy,jz)/dt) * ((us(jx,jy-1,jz)-us(jx-1,jy-1,jz))/dxx(jx) + (vs(jx-1,jy,jz)-vs(jx-1,jy-1,jz))/dyy(jy)) + AddPressureX
        XvecCrunchP(j) = pres(jx,jy,jz)
    END IF
END DO

jx = nx
DO jy = 2,ny-1
    j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1
    IF (activecellPressure(jx,jy,jz) == 0) THEN
        BvecCrunchP(j) = pres(jx,jy,jz)*big
        XvecCrunchP(j) = pres(jx,jy,jz)
    ELSE
        IF (activecellPressure(nx+1,jy,jz) == 0) THEN
          AddPressureX =  -pres(nx+1,jy,jz)/(dxx(jx)*dxx(jx))
        ELSE
          AddPressureX = 0.0d0
        END IF
        BvecCrunchP(j) = (ro(jx,jy,jz)/dt) * ((us(jx,jy-1,jz)-us(jx-1,jy-1,jz))/dxx(jx) + (vs(jx-1,jy,jz)-vs(jx-1,jy-1,jz))/dyy(jy)) + AddPressureX
        XvecCrunchP(j) = pres(jx,jy,jz)
    END IF
END DO

jx = 1
jy = 1
j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1
IF (activecellPressure(jx,jy,jz) == 0) THEN
  BvecCrunchP(j) = pres(jx,jy,jz)*big
  XvecCrunchP(j) = pres(jx,jy,jz)
ELSE
  IF (activecellPressure(0,jy,jz) == 0) THEN
    AddPressureX =  -pres(0,jy,jz)/(dxx(jx)*dxx(jx))
  ELSE
    AddPressureX = 0.0d0
  END IF
  IF (activecellPressure(jx,0,jz) == 0) THEN
    AddPressureY =  -pres(jx,0,jz)/(dyy(jy)*dyy(jy))
  ELSE
    AddPressureY = 0.0d0
  END IF
  BvecCrunchP(j) = (ro(jx,jy,jz)/dt) * ((us(jx,jy-1,jz)-us(jx-1,jy-1,jz))/dxx(jx) + (vs(jx-1,jy,jz)-vs(jx-1,jy-1,jz))/dyy(jy)) + AddPressureX + AddPressureY
  XvecCrunchP(j) = pres(jx,jy,jz)
END IF

jx = nx
jy = 1
j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1
IF (activecellPressure(jx,jy,jz) == 0) THEN
    BvecCrunchP(j) = pres(jx,jy,jz)*big
    XvecCrunchP(j) = pres(jx,jy,jz)
ELSE
    IF (activecellPressure(nx+1,jy,jz) == 0) THEN
      AddPressureX =  -pres(nx+1,jy,jz)/(dxx(jx)*dxx(jx))
    ELSE
      AddPressureX = 0.0d0
    END IF
    IF (activecellPressure(jx,0,jz) == 0) THEN
      AddPressureY =  -pres(jx,0,jz)/(dyy(jy)*dyy(jy))
    ELSE
      AddPressureY = 0.0d0
    END IF
    BvecCrunchP(j) = (ro(jx,jy,jz)/dt) * ((us(jx,jy-1,jz)-us(jx-1,jy-1,jz))/dxx(jx) + (vs(jx-1,jy,jz)-vs(jx-1,jy-1,jz))/dyy(jy)) + AddPressureX + AddPressureY
    XvecCrunchP(j) = pres(jx,jy,jz)
END IF

jx = 1
jy = ny
j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1
IF (activecellPressure(jx,jy,jz) == 0) THEN
    BvecCrunchP(j) = pres(jx,jy,jz)*big
    XvecCrunchP(j) = pres(jx,jy,jz)
ELSE
    IF (activecellPressure(0,jy,jz) == 0) THEN
      AddPressureX =  -pres(0,jy,jz)/(dxx(jx)*dxx(jx))
    ELSE
      AddPressureX = 0.0d0
    END IF
    IF (activecellPressure(jx,ny+1,jz) == 0) THEN
      AddPressureY =  -pres(jx,ny+1,jz)/(dyy(jy)*dyy(jy))
    ELSE
      AddPressureY = 0.0d0
    END IF
    BvecCrunchP(j) = (ro(jx,jy,jz)/dt) * ((us(jx,jy-1,jz)-us(jx-1,jy-1,jz))/dxx(jx) + (vs(jx-1,jy,jz)-vs(jx-1,jy-1,jz))/dyy(jy)) + AddPressureX + AddPressureY
    XvecCrunchP(j) = pres(jx,jy,jz)
END IF

jx = nx
jy = ny
j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1
IF (activecellPressure(jx,jy,jz) == 0) THEN
    BvecCrunchP(j) = pres(jx,jy,jz)*big
    XvecCrunchP(j) = pres(jx,jy,jz)
ELSE
    IF (activecellPressure(nx+1,jy,jz) == 0) THEN
      AddPressureX =  -pres(nx+1,jy,jz)/(dxx(jx)*dxx(jx))
    ELSE
      AddPressureX = 0.0d0
    END IF
    IF (activecellPressure(jx,ny+1,jz) == 0) THEN
      AddPressureY =  -pres(jx,ny+1,jz)/(dyy(jy)*dyy(jy))
    ELSE
      AddPressureY = 0.0d0
    END IF
    BvecCrunchP(j) = (ro(jx,jy,jz)/dt) * ((us(jx,jy-1,jz)-us(jx-1,jy-1,jz))/dxx(jx) + (vs(jx-1,jy,jz)-vs(jx-1,jy-1,jz))/dyy(jy)) + AddPressureX + AddPressureY
    XvecCrunchP(j) = pres(jx,jy,jz)
END IF

DO j = 0,nx*ny-1
    IF (DABS(BvecCrunchP(j)) < 1e-40) THEN
        BvecCrunchP(j) = 0.0d0
    END IF
END DO

!*********************end Matrix Right Hand Side ******************************


CALL MatAssemblyBegin(amatP,MAT_FINAL_ASSEMBLY,ierr)
CALL MatAssemblyEnd(amatP,MAT_FINAL_ASSEMBLY,ierr)

!----------------------------------------------------------------------------------
RETURN
END SUBROUTINE pressureNS



