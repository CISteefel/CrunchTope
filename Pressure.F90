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

SUBROUTINE pressure (nx,ny,nz,dt,amatP,SteadyFlow)
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

REAL(DP), INTENT(IN)                                           :: dt

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
REAL(DP)                                                              :: ScaleFactor
REAL(DP)                                                              :: pumpterm
REAL(DP)                                                              :: tdepend
REAL(DP)                                                              :: ct
REAL(DP)                                                              :: BodyForceX
REAL(DP)                                                              :: BodyForceY
REAL(DP)                                                              :: BodyForceZ
REAL(DP)                                                              :: RoAveLeft
REAL(DP)                                                              :: RoAveRight
REAL(DP)                                                              :: AddPressureX
REAL(DP)                                                              :: AddPressureY
REAL(DP)                                                              :: AddPressureZ
REAL(DP)                                                              :: Check
REAL(DP)                                                              :: printcorrectionX
REAL(DP)                                                              :: printcorrectionY

!  ****** PARAMETERS  ****************************

REAL(DP), PARAMETER                                                   :: visc=0.001d0
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

! --- formulate equation matrix for each node (including inactive cells)

  printcorrectionX = COSD(x_angle)
  printcorrectionY = COSD(y_angle)

!!  write(*,*) ' X, Y: ',printcorrectionX,printcorrectionY
!!  read(*,*)

ScaleFactor = 1.0

IF (nx > 1 .AND. ny ==1 .AND. nz == 1) THEN           ! 1D problem assuming jx is coordinate

  jz = 1
  jy = 1

  Coordinate = 'X' 
  DO jx = 2,nx-1
    j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1 
    IF (activecellPressure(jx,jy,jz) == 0) THEN
      CALL MatSetValues(amatP,1,j,1,j-1,zero,INSERT_VALUES,ierr)
      CALL MatSetValues(amatP,1,j,1,j+1,zero,INSERT_VALUES,ierr)
      CALL MatSetValues(amatP,1,j,1,j,big,INSERT_VALUES,ierr)
      BvecCrunchP(j) = pres(jx,jy,jz)*big 
      XvecCrunchP(j) = pres(jx,jy,jz)
    ELSE
      CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
      coef(-1) = -2.0d0*RoAveLeft*harx(jx-1,jy,jz)/(dxx(jx)*(dxx(jx)+dxx(jx-1))) 
      coef(1) = -2.0d0*RoAveRight*harx(jx,jy,jz)/(dxx(jx)*(dxx(jx+1)+dxx(jx))) 
      tdepend = visc*ct*por(jx,jy,jz)/ (dt*secyr)
      coef(0) = tdepend - ( coef(-1) + coef(1) )
      CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
      CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
      CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)
      
      tdepend = visc*ct*por(jx,jy,jz)*pres(jx,jy,jz)/ (dt*secyr)
!!!      tdepend = (visc*ro(jx,jy,jz)*alphaBear*(1.0d0-por(jx,jy,jz)) + visc*ct*por(jx,jy,jz) ) *pres(jx,jy,jz)/ (dt*secyr)
      pumpterm = visc*ro(jx,jy,jz)*qg(jx,jy,jz)/(secyr*dxx(jx)*dyy(jy)*dzz(jx,jy,jz))
      BodyForceX = -COSD(x_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harx(jx,jy,jz) - RoAveLeft*RoAveLeft*harx(jx-1,jy,jz))/dxx(jx)
      BvecCrunchP(j) = (pumpterm + tdepend + BodyForceX) 
      XvecCrunchP(j) = pres(jx,jy,jz)
    END IF

    
  END DO

  jx = 1
  j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1 
  IF (activecellPressure(jx,jy,jz) == 0) THEN
    CALL MatSetValues(amatP,1,j,1,j+1,zero,INSERT_VALUES,ierr)
    CALL MatSetValues(amatP,1,j,1,j,big,INSERT_VALUES,ierr)
    BvecCrunchP(j) = pres(jx,jy,jz)*big 
    XvecCrunchP(j) = pres(jx,jy,jz)

  ELSE
      
    IF (activecellPressure(0,jy,jz) == 0) THEN
      COORDINATE = 'X'
      CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
      coef(-1) = -2.0d0*ro(jx,jy,jz)*harx(0,jy,jz)/(dxx(jx)*dxx(jx))
      coef(1) = -2.0d0*RoAveRight*harx(jx,jy,jz)/(dxx(jx)*(dxx(jx+1)+dxx(jx))) 
      BodyForceX = -COSD(x_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harx(jx,jy,jz) - ro(jx,jy,jz)*ro(jx,jy,jz)*harx(jx,jy,jz))/dxx(jx)
    ELSE
      COORDINATE = 'X'
      CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
      coef(1) = -2.0d0*RoAveRight*harx(jx,jy,jz)/(dxx(jx)*(dxx(jx+1)+dxx(jx))) 
      BodyForceX = -COSD(x_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harx(jx,jy,jz) )/dxx(jx)
      coef(-1) = 0.0d0
    END IF

    tdepend = ScaleFactor*visc*ct*por(jx,jy,jz)/ (dt*secyr)
!!!      tdepend = (visc*ro(jx,jy,jz)*alphaBear*(1.0d0-por(jx,jy,jz)) + visc*ct*por(jx,jy,jz) ) / (dt*secyr)
    coef(0) = tdepend - ( coef(-1) + coef(1) )

    CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
    CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)
    
    tdepend = visc*ct*por(jx,jy,jz)*pres(jx,jy,jz)/ (dt*secyr)
    pumpterm = visc*ro(jx,jy,jz)*qg(jx,jy,jz)/(secyr*dxx(jx)*dyy(jy)*dzz(jx,jy,jz)) 
    IF (activecellPressure(0,jy,jz) == 0) THEN  
      AddPressureX =  2.0d0*ro(jx,jy,jz)*harx(0,jy,jz)*pres(0,jy,jz)/(dxx(jx)*dxx(jx))   
    ELSE
      AddPressureX = 0.0d0
    END IF 
    BvecCrunchP(j) = (pumpterm + tdepend + BodyForceX + AddPressureX) 
    XvecCrunchP(j) = pres(jx,jy,jz)
    
  END IF

  jx = nx
  j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1 
  IF (activecellPressure(jx,jy,jz) == 0) THEN
    CALL MatSetValues(amatP,1,j,1,j-1,zero,INSERT_VALUES,ierr)
    CALL MatSetValues(amatP,1,j,1,j,big,INSERT_VALUES,ierr)
    BvecCrunchP(j) = pres(jx,jy,jz)*big 
    XvecCrunchP(j) = pres(jx,jy,jz)
  ELSE
      
    IF (activecellPressure(nx+1,jy,jz) == 0) THEN
      COORDINATE = 'X'
      CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
      coef(-1) = -2.0d0*RoAveLeft*harx(jx-1,jy,jz)/(dxx(jx)*(dxx(jx)+dxx(jx-1))) 
      coef(1) = -2.0d0*ro(jx,jy,jz)*harx(jx,jy,jz)/(dxx(jx)*dxx(jx)) 
      BodyForceX = -COSD(x_angle)*grav*SignGravity*(ro(jx,jy,jz)*ro(jx,jy,jz)*harx(jx,jy,jz) - RoAveLeft*RoAveLeft*harx(jx-1,jy,jz))/dxx(jx)
    ELSE
      COORDINATE = 'X'
      CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
      coef(-1) = -2.0d0*RoAveLeft*harx(jx-1,jy,jz)/(dxx(jx)*(dxx(jx)+dxx(jx-1))) 
      BodyForceX = -COSD(x_angle)*grav*SignGravity*(-RoAveLeft*RoAveLeft*harx(jx-1,jy,jz))/dxx(jx)
      coef(1) = 0.0d0
    END IF

    tdepend = visc*ct*por(jx,jy,jz)/ (dt*secyr)
!!!      tdepend = (visc*ro(jx,jy,jz)*alphaBear*(1.0d0-por(jx,jy,jz)) + visc*ct*por(jx,jy,jz) ) / (dt*secyr)
    coef(0) = tdepend - ( coef(-1) + coef(1) )

    CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
    CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)
    
    tdepend = visc*ct*por(jx,jy,jz)*pres(jx,jy,jz)/ (dt*secyr)
    pumpterm = visc*ro(jx,jy,jz)*qg(jx,jy,jz)/(secyr*dxx(jx)*dyy(jy)*dzz(jx,jy,jz))
    IF (activecellPressure(nx+1,jy,jz) == 0) THEN  
      AddPressureX =  2.0d0*ro(jx,jy,jz)*harx(jx,jy,jz)*pres(nx+1,jy,jz)/(dxx(jx)*dxx(jx))   
    ELSE
      AddPressureX = 0.0d0
    END IF 
    BvecCrunchP(j) = (pumpterm + tdepend + BodyForceX + AddPressureX) 
    XvecCrunchP(j) = pres(jx,jy,jz)
    
  END IF

ELSE                !!  2D or 3D problem                                 
!!  ***************************************************************************************************
!!                        ************************************************

  IF (nz == 1) THEN            !!  2D problem in X and Y

    jz = 1
    
    DO jy = 2,ny-1
      DO jx = 2,nx-1
        j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1 
        IF (activecellPressure(jx,jy,jz) == 0) THEN
          CALL MatSetValues(amatP,1,j,1,j-1,zero,INSERT_VALUES,ierr)
          CALL MatSetValues(amatP,1,j,1,j+1,zero,INSERT_VALUES,ierr)
          CALL MatSetValues(amatP,1,j,1,j-nx,zero,INSERT_VALUES,ierr)
          CALL MatSetValues(amatP,1,j,1,j+nx,zero,INSERT_VALUES,ierr)
          CALL MatSetValues(amatP,1,j,1,j,big,INSERT_VALUES,ierr)
          BvecCrunchP(j) = pres(jx,jy,jz)*big 
          XvecCrunchP(j) = pres(jx,jy,jz)
        ELSE

          COORDINATE = 'X'
          CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
          coef(-1) = -2.0d0*RoAveLeft*harx(jx-1,jy,jz)/(dxx(jx)*(dxx(jx)+dxx(jx-1))) 
          coef(1) = -2.0d0*RoAveRight*harx(jx,jy,jz)/(dxx(jx)*(dxx(jx+1)+dxx(jx))) 
          BodyForceX = -COSD(x_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harx(jx,jy,jz) - RoAveLeft*RoAveLeft*harx(jx-1,jy,jz))/dxx(jx)

          COORDINATE = 'Y'
          CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
          coef(-2) = -2.0d0*RoAveLeft*hary(jx,jy-1,jz)/(dyy(jy)*(dyy(jy)+dyy(jy-1))) 
          coef(2)  = -2.0d0*RoAveRight*hary(jx,jy,jz)/(dyy(jy) *(dyy(jy+1)+dyy(jy) )) 
          BodyForceY = -COSD(y_angle)*grav*SignGravity*(RoAveRight*RoAveRight*hary(jx,jy,jz) - RoAveLeft*RoAveLeft*hary(jx,jy-1,jz))/dyy(jy) 
 
 !!!         tdepend = ScaleFactor*visc*ct*por(jx,jy,jz)/(dt*secyr)
          tdepend = (visc*ro(jx,jy,jz)*alphaBear*(1.0d0-por(jx,jy,jz)) + visc*ct*por(jx,jy,jz) ) / (dt*secyr)
          coef(0) = tdepend - ( coef(-2) + coef(-1) + coef(1) + coef(2) )

          CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
          CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
          CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
          CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
          CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)

!!!          tdepend = visc*ct*por(jx,jy,jz)*pres(jx,jy,jz)/ (dt*secyr)
          tdepend = (visc*ro(jx,jy,jz)*alphaBear*(1.0d0-por(jx,jy,jz)) + visc*ct*por(jx,jy,jz) ) *pres(jx,jy,jz)/ (dt*secyr)
          pumpterm = visc*ro(jx,jy,jz)*qg(jx,jy,jz)/(secyr*dxx(jx)*dyy(jy)*dzz(jx,jy,jz))

            !Zhi Li for debugging
            BodyForceX = 0.0d0
            BodyForceY = 0.0d0
          BvecCrunchP(j) = (pumpterm + tdepend + BodyForceX + BodyForceY) 
          XvecCrunchP(j) = pres(jx,jy,jz)

        END IF

      END DO
    END DO

      jy = 1
      DO jx = 2,nx-1
        j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1 
        IF (activecellPressure(jx,jy,jz) == 0) THEN
          CALL MatSetValues(amatP,1,j,1,j-1,zero,INSERT_VALUES,ierr)
          CALL MatSetValues(amatP,1,j,1,j+1,zero,INSERT_VALUES,ierr)
!!          CALL MatSetValues(amatP,1,j,1,j-nx,zero,INSERT_VALUES,ierr)
          CALL MatSetValues(amatP,1,j,1,j+nx,zero,INSERT_VALUES,ierr)
          CALL MatSetValues(amatP,1,j,1,j,big,INSERT_VALUES,ierr)
          BvecCrunchP(j) = pres(jx,jy,jz)*big 
          XvecCrunchP(j) = pres(jx,jy,jz)
        ELSE

          COORDINATE = 'X'
          CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
          coef(-1) = -2.0d0*RoAveLeft*harx(jx-1,jy,jz)/(dxx(jx)*(dxx(jx)+dxx(jx-1))) 
          coef(1) = -2.0d0*RoAveRight*harx(jx,jy,jz)/(dxx(jx)*(dxx(jx+1)+dxx(jx))) 
          BodyForceX = -COSD(x_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harx(jx,jy,jz) - RoAveLeft*RoAveLeft*harx(jx-1,jy,jz))/dxx(jx)

          IF (activecellPressure(jx,0,jz) == 0) THEN   !! Ghost cell with fixed pressure
            COORDINATE = 'Y'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-2) = -2.0d0*ro(jx,jy,jz)*hary(jx,0,jz)/(dyy(jy)*dyy(jy))
            coef(2)  = -2.0d0*RoAveRight*hary(jx,jy,jz)/(dyy(jy) *(dyy(jy+1)+dyy(jy) )) 
            BodyForceY = -COSD(y_angle)*grav*SignGravity*(RoAveRight*RoAveRight*hary(jx,jy,jz) - ro(jx,jy,jz)*ro(jx,jy,jz)*hary(jx,jy,jz))/dyy(jy) 
 
          ELSE
            COORDINATE = 'Y'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(2)  = -2.0d0*RoAveRight*hary(jx,jy,jz)/(dyy(jy) *(dyy(jy+1)+dyy(jy) )) 
            BodyForceY = -COSD(y_angle)*grav*SignGravity*(RoAveRight*RoAveRight*hary(jx,jy,jz))/dyy(jy) 
            coef(-2) = 0.0d0
          END IF

 !!!         tdepend = ScaleFactor*visc*ct*por(jx,jy,jz)/(dt*secyr)
          tdepend = (visc*ro(jx,jy,jz)*alphaBear*(1.0d0-por(jx,jy,jz)) + visc*ct*por(jx,jy,jz) ) / (dt*secyr)
          coef(0) = tdepend - ( coef(-2) + coef(-1) + coef(1) + coef(2) )

          CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
          CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
          CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
          CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)

!!!          tdepend = visc*ct*por(jx,jy,jz)*pres(jx,jy,jz)/ (dt*secyr)
          tdepend = (visc*ro(jx,jy,jz)*alphaBear*(1.0d0-por(jx,jy,jz)) + visc*ct*por(jx,jy,jz) ) *pres(jx,jy,jz)/ (dt*secyr)
          pumpterm = visc*ro(jx,jy,jz)*qg(jx,jy,jz)/(secyr*dxx(jx)*dyy(jy)*dzz(jx,jy,jz))
          IF (activecellPressure(jx,0,jz) == 0) THEN  
            AddPressureY =  2.0d0*ro(jx,jy,jz)*hary(jx,0,jz)*pres(jx,0,jz)/(dyy(jy)*dyy(jy))   
          ELSE
            AddPressureY = 0.0d0
          END IF
            !Zhi Li for debugging
            BodyForceX = 0.0d0
            BodyForceY = 0.0d0
          BvecCrunchP(j) = pumpterm + tdepend + BodyForceX + BodyForceY + AddPressureY
          XvecCrunchP(j) = pres(jx,jy,jz)
        END IF
        
      END DO

      jy = ny
      DO jx = 2,nx-1
        j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1 
        IF (activecellPressure(jx,jy,jz) == 0) THEN
          CALL MatSetValues(amatP,1,j,1,j-1,zero,INSERT_VALUES,ierr)
          CALL MatSetValues(amatP,1,j,1,j+1,zero,INSERT_VALUES,ierr)
          CALL MatSetValues(amatP,1,j,1,j-nx,zero,INSERT_VALUES,ierr)
          CALL MatSetValues(amatP,1,j,1,j,big,INSERT_VALUES,ierr)
          BvecCrunchP(j) = pres(jx,jy,jz)*big 
          XvecCrunchP(j) = pres(jx,jy,jz)
        ELSE

          COORDINATE = 'X'
          CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
          coef(-1) = -2.0d0*RoAveLeft*harx(jx-1,jy,jz)/(dxx(jx)*(dxx(jx)+dxx(jx-1))) 
          coef(1) = -2.0d0*RoAveRight*harx(jx,jy,jz)/(dxx(jx)*(dxx(jx+1)+dxx(jx))) 
          BodyForceX = -COSD(x_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harx(jx,jy,jz) - RoAveLeft*RoAveLeft*harx(jx-1,jy,jz))/dxx(jx)

          IF (activecellPressure(jx,ny+1,jz) == 0) THEN
            COORDINATE = 'Y'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-2) = -2.0d0*RoAveLeft*hary(jx,jy-1,jz)/(dyy(jy)*(dyy(jy)+dyy(jy-1))) 
            coef(2)  = -2.0d0*ro(jx,jy,jz)*hary(jx,jy,jz)/(dyy(jy)*dyy(jy)) 
            BodyForceY = -COSD(y_angle)*grav*SignGravity*(ro(jx,jy,jz)*ro(jx,jy,jz)*hary(jx,jy,jz) - RoAveLeft*RoAveLeft*hary(jx,jy-1,jz))/dyy(jy) 
          ELSE
            COORDINATE = 'Y'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-2) = -2.0d0*RoAveLeft*hary(jx,jy-1,jz)/(dyy(jy)*(dyy(jy)+dyy(jy-1))) 
            BodyForceY = -COSD(y_angle)*grav*SignGravity*(-RoAveLeft*RoAveLeft*hary(jx,jy-1,jz))/dyy(jy)    
            coef(2) = 0.0d0     
          END IF
          
 !!!         tdepend = ScaleFactor*visc*ct*por(jx,jy,jz)/(dt*secyr)
          tdepend = (visc*ro(jx,jy,jz)*alphaBear*(1.0d0-por(jx,jy,jz)) + visc*ct*por(jx,jy,jz) ) / (dt*secyr)
          coef(0) = tdepend - ( coef(-2) + coef(-1) + coef(1) + coef(2) )

          CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
          CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
          CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
          CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)

!!!          tdepend = visc*ct*por(jx,jy,jz)*pres(jx,jy,jz)/ (dt*secyr)
          tdepend = (visc*ro(jx,jy,jz)*alphaBear*(1.0d0-por(jx,jy,jz)) + visc*ct*por(jx,jy,jz) ) *pres(jx,jy,jz)/ (dt*secyr)
          pumpterm = visc*ro(jx,jy,jz)*qg(jx,jy,jz)/(secyr*dxx(jx)*dyy(jy)*dzz(jx,jy,jz))

          IF (activecellPressure(jx,ny+1,jz) == 0) THEN  
            AddPressureY =  2.0d0*ro(jx,jy,jz)*hary(jx,jy,jz)*pres(jx,ny+1,jz)/(dyy(jy)*dyy(jy))   
          ELSE
            AddPressureY = 0.0d0
          END IF
            !Zhi Li for debugging
            BodyForceX = 0.0d0
            BodyForceY = 0.0d0
          BvecCrunchP(j) = (pumpterm + tdepend + BodyForceX + BodyForceY + AddPressureY) 
          XvecCrunchP(j) = pres(jx,jy,jz)
        END IF
      END DO

      jx = 1
      DO jy = 2,ny-1
        j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1 
        IF (activecellPressure(jx,jy,jz) == 0) THEN
!!          CALL MatSetValues(amatP,1,j,1,j-1,zero,INSERT_VALUES,ierr)
          CALL MatSetValues(amatP,1,j,1,j+1,zero,INSERT_VALUES,ierr)
          CALL MatSetValues(amatP,1,j,1,j-nx,zero,INSERT_VALUES,ierr)
          CALL MatSetValues(amatP,1,j,1,j+nx,zero,INSERT_VALUES,ierr)
          CALL MatSetValues(amatP,1,j,1,j,big,INSERT_VALUES,ierr)
		  
          BvecCrunchP(j) = pres(jx,jy,jz)*big 
          XvecCrunchP(j) = pres(jx,jy,jz)
        ELSE

          COORDINATE = 'Y'
          CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
          coef(-2) = -2.0d0*RoAveLeft*hary(jx,jy-1,jz)/(dyy(jy)*(dyy(jy)+dyy(jy-1))) 
          coef(2)  = -2.0d0*RoAveRight*hary(jx,jy,jz)/(dyy(jy) *(dyy(jy+1)+dyy(jy) )) 
          BodyForceY = -COSD(y_angle)*grav*SignGravity*(RoAveRight*RoAveRight*hary(jx,jy,jz) - RoAveLeft*RoAveLeft*hary(jx,jy-1,jz))/dyy(jy) 

          IF (activecellPressure(0,jy,jz) == 0) THEN
            COORDINATE = 'X'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-1) = -2.0d0*ro(jx,jy,jz)*harx(0,jy,jz)/(dxx(jx)*dxx(jx))
            coef(1) = -2.0d0*RoAveRight*harx(jx,jy,jz)/(dxx(jx)*(dxx(jx+1)+dxx(jx))) 
            BodyForceX = -COSD(x_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harx(jx,jy,jz) - ro(jx,jy,jz)*ro(jx,jy,jz)*harx(jx,jy,jz))/dxx(jx)
          ELSE
            COORDINATE = 'X'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(1) = -2.0d0*RoAveRight*harx(jx,jy,jz)/(dxx(jx)*(dxx(jx+1)+dxx(jx))) 
            BodyForceX = -COSD(x_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harx(jx,jy,jz) )/dxx(jx)
            coef(-1) = 0.0d0
          END IF

 !!!         tdepend = ScaleFactor*visc*ct*por(jx,jy,jz)/(dt*secyr)
          tdepend = (visc*ro(jx,jy,jz)*alphaBear*(1.0d0-por(jx,jy,jz)) + visc*ct*por(jx,jy,jz) ) / (dt*secyr)
          coef(0) = tdepend - ( coef(-2) + coef(-1) + coef(1) + coef(2) )

          CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
          CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
          CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
          CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)

!!!          tdepend = visc*ct*por(jx,jy,jz)*pres(jx,jy,jz)/ (dt*secyr)
          tdepend = (visc*ro(jx,jy,jz)*alphaBear*(1.0d0-por(jx,jy,jz)) + visc*ct*por(jx,jy,jz) ) *pres(jx,jy,jz)/ (dt*secyr)
          pumpterm = visc*ro(jx,jy,jz)*qg(jx,jy,jz)/(secyr*dxx(jx)*dyy(jy)*dzz(jx,jy,jz))
          IF (activecellPressure(0,jy,jz) == 0) THEN  
            AddPressureX =  2.0d0*ro(jx,jy,jz)*harx(0,jy,jz)*pres(0,jy,jz)/(dxx(jx)*dxx(jx))  
          ELSE
            AddPressureX = 0.0d0
          END IF
            !Zhi Li for debugging
            BodyForceX = 0.0d0
            BodyForceY = 0.0d0
          BvecCrunchP(j) = (pumpterm + tdepend + BodyForceX + BodyForceY + AddPressureX) 
          XvecCrunchP(j) = pres(jx,jy,jz)
        END IF
      END DO

      jx = nx
      DO jy = 2,ny-1
        j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1 
        IF (activecellPressure(jx,jy,jz) == 0) THEN
          CALL MatSetValues(amatP,1,j,1,j-1,zero,INSERT_VALUES,ierr)
          CALL MatSetValues(amatP,1,j,1,j-nx,zero,INSERT_VALUES,ierr)
          CALL MatSetValues(amatP,1,j,1,j+nx,zero,INSERT_VALUES,ierr)
          CALL MatSetValues(amatP,1,j,1,j,big,INSERT_VALUES,ierr)
          BvecCrunchP(j) = pres(jx,jy,jz)*big 
          XvecCrunchP(j) = pres(jx,jy,jz)
        ELSE

          COORDINATE = 'Y'
          CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
          coef(-2) = -2.0d0*RoAveLeft*hary(jx,jy-1,jz)/(dyy(jy)*(dyy(jy)+dyy(jy-1))) 
          coef(2)  = -2.0d0*RoAveRight*hary(jx,jy,jz)/(dyy(jy) *(dyy(jy+1)+dyy(jy) )) 
          BodyForceY = -COSD(y_angle)*grav*SignGravity*(RoAveRight*RoAveRight*hary(jx,jy,jz) - RoAveLeft*RoAveLeft*hary(jx,jy-1,jz))/dyy(jy) 
 
          IF (activecellPressure(nx+1,jy,jz) == 0) THEN
            COORDINATE = 'X'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-1) = -2.0d0*RoAveLeft*harx(jx-1,jy,jz)/(dxx(jx)*(dxx(jx)+dxx(jx-1))) 
            coef(1) = -2.0d0*ro(jx,jy,jz)*harx(jx,jy,jz)/(dxx(jx)*dxx(jx)) 
            BodyForceX = -COSD(x_angle)*grav*SignGravity*(ro(jx,jy,jz)*ro(jx,jy,jz)*harx(jx,jy,jz) - RoAveLeft*RoAveLeft*harx(jx-1,jy,jz))/dxx(jx)
          ELSE
            COORDINATE = 'X'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-1) = -2.0d0*RoAveLeft*harx(jx-1,jy,jz)/(dxx(jx)*(dxx(jx)+dxx(jx-1))) 
            BodyForceX = -COSD(x_angle)*grav*SignGravity*(-RoAveLeft*RoAveLeft*harx(jx-1,jy,jz))/dxx(jx)
            coef(1) = 0.0d0
          END IF

 !!!         tdepend = ScaleFactor*visc*ct*por(jx,jy,jz)/(dt*secyr)
          tdepend = (visc*ro(jx,jy,jz)*alphaBear*(1.0d0-por(jx,jy,jz)) + visc*ct*por(jx,jy,jz) ) / (dt*secyr)
          coef(0) = tdepend - ( coef(-2) + coef(-1) + coef(1) + coef(2) )

          CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
          CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
          CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
          CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)

!!!          tdepend = visc*ct*por(jx,jy,jz)*pres(jx,jy,jz)/ (dt*secyr)
          tdepend = (visc*ro(jx,jy,jz)*alphaBear*(1.0d0-por(jx,jy,jz)) + visc*ct*por(jx,jy,jz) ) *pres(jx,jy,jz)/ (dt*secyr)
          pumpterm = visc*ro(jx,jy,jz)*qg(jx,jy,jz)/(secyr*dxx(jx)*dyy(jy)*dzz(jx,jy,jz))
          IF (activecellPressure(nx+1,jy,jz) == 0) THEN  
            AddPressureX =  2.0d0*ro(jx,jy,jz)*harx(jx,jy,jz)*pres(nx+1,jy,jz)/(dxx(jx)*dxx(jx))   
          ELSE
            AddPressureX = 0.0d0
          END IF
            !Zhi Li for debugging
            BodyForceX = 0.0d0
            BodyForceY = 0.0d0
          BvecCrunchP(j) = (pumpterm + tdepend + BodyForceX + BodyForceY + AddPressureX) 
          XvecCrunchP(j) = pres(jx,jy,jz)
        END IF
      END DO

      jx = 1
      jy = 1
        j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1 
        IF (activecellPressure(jx,jy,jz) == 0) THEN
!!          CALL MatSetValues(amatP,1,j,1,j-1,zero,INSERT_VALUES,ierr)
          CALL MatSetValues(amatP,1,j,1,j+1,zero,INSERT_VALUES,ierr)
!!          CALL MatSetValues(amatP,1,j,1,j-nx,zero,INSERT_VALUES,ierr)
          CALL MatSetValues(amatP,1,j,1,j+nx,zero,INSERT_VALUES,ierr)
          CALL MatSetValues(amatP,1,j,1,j,big,INSERT_VALUES,ierr)
          BvecCrunchP(j) = pres(jx,jy,jz)*big 
          XvecCrunchP(j) = pres(jx,jy,jz)
        ELSE

          IF (activecellPressure(0,jy,jz) == 0) THEN
            COORDINATE = 'X'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-1) = -2.0d0*ro(jx,jy,jz)*harx(0,jy,jz)/(dxx(jx)*dxx(jx))
            coef(1) = -2.0d0*RoAveRight*harx(jx,jy,jz)/(dxx(jx)*(dxx(jx+1)+dxx(jx))) 
            BodyForceX = -COSD(x_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harx(jx,jy,jz) - ro(jx,jy,jz)*ro(jx,jy,jz)*harx(jx,jy,jz))/dxx(jx)
          ELSE
            COORDINATE = 'X'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(1) = -2.0d0*RoAveRight*harx(jx,jy,jz)/(dxx(jx)*(dxx(jx+1)+dxx(jx))) 
            BodyForceX = -COSD(x_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harx(jx,jy,jz) )/dxx(jx)
            coef(-1) = 0.0d0
          END IF
           
          IF (activecellPressure(jx,0,jz) == 0) THEN
            COORDINATE = 'Y'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(2)  = -2.0d0*RoAveRight*hary(jx,jy,jz)/(dyy(jy) *(dyy(jy+1)+dyy(jy) )) 
            coef(-2) = -2.0d0*ro(jx,jy,jz)*hary(jx,0,jz)/(dyy(jy)*dyy(jy))
            BodyForceY = -COSD(y_angle)*grav*SignGravity*(RoAveRight*RoAveRight*hary(jx,jy,jz) - ro(jx,jy,jz)*ro(jx,jy,jz)*hary(jx,jy,jz))/dyy(jy) 
          ELSE
            COORDINATE = 'Y'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(2)  = -2.0d0*RoAveRight*hary(jx,jy,jz)/(dyy(jy) *(dyy(jy+1)+dyy(jy) )) 
            BodyForceY = -COSD(y_angle)*grav*SignGravity*(RoAveRight*RoAveRight*hary(jx,jy,jz))/dyy(jy) 
            coef(-2) = 0.0d0
          END IF

 !!!         tdepend = ScaleFactor*visc*ct*por(jx,jy,jz)/(dt*secyr)
          tdepend = (visc*ro(jx,jy,jz)*alphaBear*(1.0d0-por(jx,jy,jz)) + visc*ct*por(jx,jy,jz) ) / (dt*secyr)
          coef(0) = tdepend - ( coef(-2) + coef(-1) + coef(1) + coef(2) )

          CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
          CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
          CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)

!!!          tdepend = visc*ct*por(jx,jy,jz)*pres(jx,jy,jz)/ (dt*secyr)
          tdepend = (visc*ro(jx,jy,jz)*alphaBear*(1.0d0-por(jx,jy,jz)) + visc*ct*por(jx,jy,jz) ) *pres(jx,jy,jz)/ (dt*secyr)
          pumpterm = visc*ro(jx,jy,jz)*qg(jx,jy,jz)/(secyr*dxx(jx)*dyy(jy)*dzz(jx,jy,jz))
          IF (activecellPressure(0,jy,jz) == 0) THEN  
            AddPressureX =  2.0d0*ro(jx,jy,jz)*harx(0,jy,jz)*pres(0,jy,jz)/(dxx(jx)*dxx(jx))   
          ELSE
            AddPressureX = 0.0d0
          END IF 
          IF (activecellPressure(jx,0,jz) == 0) THEN  
            AddPressureY =  2.0d0*ro(jx,jy,jz)*hary(jx,0,jz)*pres(jx,0,jz)/(dyy(jy)*dyy(jy))   
          ELSE
            AddPressureY = 0.0d0
          END IF
            !Zhi Li for debugging
            BodyForceX = 0.0d0
            BodyForceY = 0.0d0
          BvecCrunchP(j) = (pumpterm + tdepend + BodyForceX + BodyForceY + AddPressureX + AddPressureY) 
          XvecCrunchP(j) = pres(jx,jy,jz)
        END IF

      jx = nx
      jy = 1
        j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1 
        IF (activecellPressure(jx,jy,jz) == 0) THEN
          CALL MatSetValues(amatP,1,j,1,j-1,zero,INSERT_VALUES,ierr)
!!          CALL MatSetValues(amatP,1,j,1,j+1,zero,INSERT_VALUES,ierr)
!!          CALL MatSetValues(amatP,1,j,1,j-nx,zero,INSERT_VALUES,ierr)
          CALL MatSetValues(amatP,1,j,1,j+nx,zero,INSERT_VALUES,ierr)
          CALL MatSetValues(amatP,1,j,1,j,big,INSERT_VALUES,ierr)
          BvecCrunchP(j) = pres(jx,jy,jz)*big 
          XvecCrunchP(j) = pres(jx,jy,jz)
        ELSE

          IF (activecellPressure(nx+1,jy,jz) == 0) THEN
            COORDINATE = 'X'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-1) = -2.0d0*RoAveLeft*harx(jx-1,jy,jz)/(dxx(jx)*(dxx(jx)+dxx(jx-1))) 
            coef(1) = -2.0d0*ro(jx,jy,jz)*harx(jx,jy,jz)/(dxx(jx)*dxx(jx)) 
            BodyForceX = -COSD(x_angle)*grav*SignGravity*(ro(jx,jy,jz)*ro(jx,jy,jz)*harx(jx,jy,jz) - RoAveLeft*RoAveLeft*harx(jx-1,jy,jz))/dxx(jx)
          ELSE
            COORDINATE = 'X'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-1) = -2.0d0*RoAveLeft*harx(jx-1,jy,jz)/(dxx(jx)*(dxx(jx)+dxx(jx-1))) 
            BodyForceX = -COSD(x_angle)*grav*SignGravity*(-RoAveLeft*RoAveLeft*harx(jx-1,jy,jz))/dxx(jx)
            coef(1) = 0.0d0
          END IF

          IF (activecellPressure(jx,0,jz) == 0) THEN
            COORDINATE = 'Y'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(2)  = -2.0d0*RoAveRight*hary(jx,jy,jz)/(dyy(jy) *(dyy(jy+1)+dyy(jy) )) 
            coef(-2) = -2.0d0*ro(jx,jy,jz)*hary(jx,0,jz)/(dyy(jy)*dyy(jy))
            BodyForceY = -COSD(y_angle)*grav*SignGravity*(RoAveRight*RoAveRight*hary(jx,jy,jz) - ro(jx,jy,jz)*ro(jx,jy,jz)*hary(jx,jy,jz))/dyy(jy) 
          ELSE
            COORDINATE = 'Y'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(2)  = -2.0d0*RoAveRight*hary(jx,jy,jz)/(dyy(jy) *(dyy(jy+1)+dyy(jy) )) 
            BodyForceY = -COSD(y_angle)*grav*SignGravity*(RoAveRight*RoAveRight*hary(jx,jy,jz))/dyy(jy) 
            coef(-2) = 0.0d0
          END IF
 
 !!!         tdepend = ScaleFactor*visc*ct*por(jx,jy,jz)/(dt*secyr)
          tdepend = (visc*ro(jx,jy,jz)*alphaBear*(1.0d0-por(jx,jy,jz)) + visc*ct*por(jx,jy,jz) ) / (dt*secyr)
          coef(0) = tdepend - ( coef(-2) + coef(-1) + coef(1) + coef(2) )

          CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
          CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
          CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)

!!!          tdepend = visc*ct*por(jx,jy,jz)*pres(jx,jy,jz)/ (dt*secyr)
          tdepend = (visc*ro(jx,jy,jz)*alphaBear*(1.0d0-por(jx,jy,jz)) + visc*ct*por(jx,jy,jz) ) *pres(jx,jy,jz)/ (dt*secyr)
          pumpterm = visc*ro(jx,jy,jz)*qg(jx,jy,jz)/(secyr*dxx(jx)*dyy(jy)*dzz(jx,jy,jz))
          IF (activecellPressure(nx+1,jy,jz) == 0) THEN  
            AddPressureX =  2.0d0*ro(jx,jy,jz)*harx(jx,jy,jz)*pres(nx+1,jy,jz)/(dxx(jx)*dxx(jx))   
          ELSE
            AddPressureX = 0.0d0
          END IF 
          IF (activecellPressure(jx,0,jz) == 0) THEN  
            AddPressureY =  2.0d0*ro(jx,jy,jz)*hary(jx,0,jz)*pres(jx,0,jz)/(dyy(jy)*dyy(jy))   
          ELSE
            AddPressureY = 0.0d0
          END IF
            !Zhi Li for debugging
            BodyForceX = 0.0d0
            BodyForceY = 0.0d0
          BvecCrunchP(j) = (pumpterm + tdepend + BodyForceX + BodyForceY + AddPressureX + AddPressureY) 
          XvecCrunchP(j) = pres(jx,jy,jz)
        END IF

      jx = 1
      jy = ny
        j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1 
        IF (activecellPressure(jx,jy,jz) == 0) THEN
!!          CALL MatSetValues(amatP,1,j,1,j-1,zero,INSERT_VALUES,ierr)
          CALL MatSetValues(amatP,1,j,1,j+1,zero,INSERT_VALUES,ierr)
          CALL MatSetValues(amatP,1,j,1,j-nx,zero,INSERT_VALUES,ierr)
!!          CALL MatSetValues(amatP,1,j,1,j+nx,zero,INSERT_VALUES,ierr)
          CALL MatSetValues(amatP,1,j,1,j,big,INSERT_VALUES,ierr)
          BvecCrunchP(j) = pres(jx,jy,jz)*big 
          XvecCrunchP(j) = pres(jx,jy,jz)
        ELSE

          IF (activecellPressure(0,jy,jz) == 0) THEN
            COORDINATE = 'X'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-1) = -2.0d0*ro(jx,jy,jz)*harx(0,jy,jz)/(dxx(jx)*dxx(jx))
            coef(1) = -2.0d0*RoAveRight*harx(jx,jy,jz)/(dxx(jx)*(dxx(jx+1)+dxx(jx))) 
            BodyForceX = -COSD(x_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harx(jx,jy,jz) - ro(jx,jy,jz)*ro(jx,jy,jz)*harx(jx,jy,jz))/dxx(jx)
          ELSE
            COORDINATE = 'X'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(1) = -2.0d0*RoAveRight*harx(jx,jy,jz)/(dxx(jx)*(dxx(jx+1)+dxx(jx))) 
            BodyForceX = -COSD(x_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harx(jx,jy,jz) )/dxx(jx)
            coef(-1) = 0.0d0
          END IF

          IF (activecellPressure(jx,ny+1,jz) == 0) THEN
            COORDINATE = 'Y'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-2) = -2.0d0*RoAveLeft*hary(jx,jy-1,jz)/(dyy(jy)*(dyy(jy)+dyy(jy-1))) 
            coef(2)  = -2.0d0*ro(jx,jy,jz)*hary(jx,jy,jz)/(dyy(jy)*dyy(jy)) 
            BodyForceY = -COSD(y_angle)*grav*SignGravity*(ro(jx,jy,jz)*ro(jx,jy,jz)*hary(jx,jy,jz) - RoAveLeft*RoAveLeft*hary(jx,jy-1,jz))/dyy(jy) 
          ELSE
            COORDINATE = 'Y'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-2) = -2.0d0*RoAveLeft*hary(jx,jy-1,jz)/(dyy(jy)*(dyy(jy)+dyy(jy-1))) 
            BodyForceY = -COSD(y_angle)*grav*SignGravity*(-RoAveLeft*RoAveLeft*hary(jx,jy-1,jz))/dyy(jy)    
            coef(2) = 0.0d0       
          END IF

 !!!         tdepend = ScaleFactor*visc*ct*por(jx,jy,jz)/(dt*secyr)
          tdepend = (visc*ro(jx,jy,jz)*alphaBear*(1.0d0-por(jx,jy,jz)) + visc*ct*por(jx,jy,jz) ) / (dt*secyr)
          coef(0) = tdepend - ( coef(-2) + coef(-1) + coef(1) + coef(2) )

          CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
          CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
          CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)

!!!          tdepend = visc*ct*por(jx,jy,jz)*pres(jx,jy,jz)/ (dt*secyr)
          tdepend = (visc*ro(jx,jy,jz)*alphaBear*(1.0d0-por(jx,jy,jz)) + visc*ct*por(jx,jy,jz) ) *pres(jx,jy,jz)/ (dt*secyr)
          pumpterm = visc*ro(jx,jy,jz)*qg(jx,jy,jz)/(secyr*dxx(jx)*dyy(jy)*dzz(jx,jy,jz))
          IF (activecellPressure(0,jy,jz) == 0) THEN  
            AddPressureX =  2.0d0*ro(jx,jy,jz)*harx(0,jy,jz)*pres(0,jy,jz)/(dxx(jx)*dxx(jx))   
          ELSE
            AddPressureX = 0.0d0
          END IF 
          IF (activecellPressure(jx,ny+1,jz) == 0) THEN  
            AddPressureY =  2.0d0*ro(jx,jy,jz)*hary(jx,jy,jz)*pres(jx,ny+1,jz)/(dyy(jy)*dyy(jy))   
          ELSE
            AddPressureY = 0.0d0
          END IF
            !Zhi Li for debugging
            BodyForceX = 0.0d0
            BodyForceY = 0.0d0
          BvecCrunchP(j) = (pumpterm + tdepend + BodyForceX + BodyForceY + AddPressureX + AddPressureY) 
          XvecCrunchP(j) = pres(jx,jy,jz)
        END IF

      jx = nx
      jy = ny
        j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1 
        IF (activecellPressure(jx,jy,jz) == 0) THEN
          CALL MatSetValues(amatP,1,j,1,j-1,zero,INSERT_VALUES,ierr)
!!          CALL MatSetValues(amatP,1,j,1,j+1,zero,INSERT_VALUES,ierr)
          CALL MatSetValues(amatP,1,j,1,j-nx,zero,INSERT_VALUES,ierr)
!!          CALL MatSetValues(amatP,1,j,1,j+nx,zero,INSERT_VALUES,ierr)
          CALL MatSetValues(amatP,1,j,1,j,big,INSERT_VALUES,ierr)
          BvecCrunchP(j) = pres(jx,jy,jz)*big 
          XvecCrunchP(j) = pres(jx,jy,jz)
        ELSE

          IF (activecellPressure(nx+1,jy,jz) == 0) THEN
            COORDINATE = 'X'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-1) = -2.0d0*RoAveLeft*harx(jx-1,jy,jz)/(dxx(jx)*(dxx(jx)+dxx(jx-1))) 
            coef(1) = -2.0d0*ro(jx,jy,jz)*harx(jx,jy,jz)/(dxx(jx)*dxx(jx)) 
            BodyForceX = -COSD(x_angle)*grav*SignGravity*(ro(jx,jy,jz)*ro(jx,jy,jz)*harx(jx,jy,jz) - RoAveLeft*RoAveLeft*harx(jx-1,jy,jz))/dxx(jx)
          ELSE
            COORDINATE = 'X'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-1) = -2.0d0*RoAveLeft*harx(jx-1,jy,jz)/(dxx(jx)*(dxx(jx)+dxx(jx-1))) 
            BodyForceX = -COSD(x_angle)*grav*SignGravity*(-RoAveLeft*RoAveLeft*harx(jx-1,jy,jz))/dxx(jx)
            coef(1) = 0.0d0
          END IF

          IF (activecellPressure(jx,ny+1,jz) == 0) THEN
            COORDINATE = 'Y'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-2) = -2.0d0*RoAveLeft*hary(jx,jy-1,jz)/(dyy(jy)*(dyy(jy)+dyy(jy-1))) 
            coef(2)  = -2.0d0*ro(jx,jy,jz)*hary(jx,jy,jz)/(dyy(jy)*dyy(jy)) 
            BodyForceY = -COSD(y_angle)*grav*SignGravity*(ro(jx,jy,jz)*ro(jx,jy,jz)*hary(jx,jy,jz) - RoAveLeft*RoAveLeft*hary(jx,jy-1,jz))/dyy(jy) 
          ELSE
            COORDINATE = 'Y'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-2) = -2.0d0*RoAveLeft*hary(jx,jy-1,jz)/(dyy(jy)*(dyy(jy)+dyy(jy-1))) 
            BodyForceY = -COSD(y_angle)*grav*SignGravity*(-RoAveLeft*RoAveLeft*hary(jx,jy-1,jz))/dyy(jy)    
            coef(2) = 0.0d0       
          END IF
 
 !!!         tdepend = ScaleFactor*visc*ct*por(jx,jy,jz)/(dt*secyr)
          tdepend = (visc*ro(jx,jy,jz)*alphaBear*(1.0d0-por(jx,jy,jz)) + visc*ct*por(jx,jy,jz) ) / (dt*secyr)
          coef(0) = tdepend - ( coef(-2) + coef(-1) + coef(1) + coef(2) )

          CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
          CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
          CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)

!!!          tdepend = visc*ct*por(jx,jy,jz)*pres(jx,jy,jz)/ (dt*secyr)
          tdepend = (visc*ro(jx,jy,jz)*alphaBear*(1.0d0-por(jx,jy,jz)) + visc*ct*por(jx,jy,jz) ) *pres(jx,jy,jz)/ (dt*secyr)
          pumpterm = visc*ro(jx,jy,jz)*qg(jx,jy,jz)/(secyr*dxx(jx)*dyy(jy)*dzz(jx,jy,jz))
          IF (activecellPressure(nx+1,jy,jz) == 0) THEN  
            AddPressureX =  2.0d0*ro(jx,jy,jz)*harx(jx,jy,jz)*pres(nx+1,jy,jz)/(dxx(jx)*dxx(jx))  
          ELSE
            AddPressureX = 0.0d0
          END IF 
          IF (activecellPressure(jx,ny+1,jz) == 0) THEN  
            AddPressureY =  2.0d0*ro(jx,jy,jz)*hary(jx,jy,jz)*pres(jx,ny+1,jz)/(dyy(jy)*dyy(jy))   
          ELSE
            AddPressureY = 0.0d0
          END IF
            !Zhi Li for debugging
            BodyForceX = 0.0d0
            BodyForceY = 0.0d0
          BvecCrunchP(j) = (pumpterm + tdepend + BodyForceX + BodyForceY + AddPressureX + AddPressureY) 
          XvecCrunchP(j) = pres(jx,jy,jz)

        END IF

!!********************************************************************************************************************************************
  ELSE           !  3D problem
!!                  ****************************************************************************************************

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
            BvecCrunchP(j) = pres(jx,jy,jz)*big 
            XvecCrunchP(j) = pres(jx,jy,jz)
          ELSE

            COORDINATE = 'X'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-1) = -2.0d0*RoAveLeft*harx(jx-1,jy,jz)/(dxx(jx)*(dxx(jx)+dxx(jx-1))) 
            coef(1) = -2.0d0*RoAveRight*harx(jx,jy,jz)/(dxx(jx)*(dxx(jx+1)+dxx(jx))) 
            BodyForceX = -COSD(x_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harx(jx,jy,jz) - RoAveLeft*RoAveLeft*harx(jx-1,jy,jz))/dxx(jx)

            COORDINATE = 'Y'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-2) = -2.0d0*RoAveLeft*hary(jx,jy-1,jz)/(dyy(jy)*(dyy(jy)+dyy(jy-1))) 
            coef(2)  = -2.0d0*RoAveRight*hary(jx,jy,jz)/(dyy(jy) *(dyy(jy+1)+dyy(jy) )) 
            BodyForceY = -COSD(y_angle)*grav*SignGravity*(RoAveRight*RoAveRight*hary(jx,jy,jz) - RoAveLeft*RoAveLeft*hary(jx,jy-1,jz))/dyy(jy) 

            COORDINATE = 'Z'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-3) = -2.0d0*RoAveLeft*harz(jx,jy,jz-1)/(dzz(jx,jy,jz)*(dzz(jx,jy,jz)+dzz(jx,jy,jz-1))) 
            coef(3)  = -2.0d0*RoAveRight*harz(jx,jy,jz)/(dzz(jx,jy,jz) *(dzz(jx,jy,jz+1)+dzz(jx,jy,jz) )) 
            BodyForceZ = -COSD(z_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harz(jx,jy,jz) - RoAveLeft*RoAveLeft*harz(jx,jy,jz-1))/dzz(jx,jy,jz) 
 
            tdepend = ScaleFactor*visc*ct*por(jx,jy,jz)/(dt*secyr)
            coef(0) = tdepend -( coef(-3) + coef(-2) + coef(-1) + coef(1) + coef(2) + coef(3) )

            CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j-nx*ny,coef(-3),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j+nx*ny,coef(3),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)

            tdepend = visc*ct*por(jx,jy,jz)*pres(jx,jy,jz)/ (dt*secyr)
            pumpterm = visc*ro(jx,jy,jz)*qg(jx,jy,jz)/(secyr*dxx(jx)*dyy(jy)*dzz(jx,jy,jz))
            BvecCrunchP(j) = (pumpterm + tdepend + BodyForceX + BodyForceY + BodyForceZ) 
            XvecCrunchP(j) = pres(jx,jy,jz)

          END IF
        END DO
      END DO
    END DO

!!   ******************************************************************************************************************************************
!!    Now the faces

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
            BvecCrunchP(j) = pres(jx,jy,jz)*big 
            XvecCrunchP(j) = pres(jx,jy,jz)
          ELSE

          IF (activecellPressure(0,jy,jz) == 0) THEN
            COORDINATE = 'X'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-1) = -2.0d0*ro(jx,jy,jz)*harx(jx,jy,jz)/(dxx(jx)*dxx(jx))
            coef(1) = -2.0d0*RoAveRight*harx(jx,jy,jz)/(dxx(jx)*(dxx(jx+1)+dxx(jx))) 
            BodyForceX = -COSD(x_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harx(jx,jy,jz) - ro(jx,jy,jz)*ro(jx,jy,jz)*harx(jx,jy,jz))/dxx(jx)
          ELSE
            COORDINATE = 'X'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(1) = -2.0d0*RoAveRight*harx(jx,jy,jz)/(dxx(jx)*(dxx(jx+1)+dxx(jx))) 
            BodyForceX = -COSD(x_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harx(jx,jy,jz) )/dxx(jx)
            coef(-1) = 0.0d0
          END IF

            COORDINATE = 'Y'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-2) = -2.0d0*RoAveLeft*hary(jx,jy-1,jz)/(dyy(jy)*(dyy(jy)+dyy(jy-1))) 
            coef(2)  = -2.0d0*RoAveRight*hary(jx,jy,jz)/(dyy(jy) *(dyy(jy+1)+dyy(jy) )) 
            BodyForceY = -COSD(y_angle)*grav*SignGravity*(RoAveRight*RoAveRight*hary(jx,jy,jz) - RoAveLeft*RoAveLeft*hary(jx,jy-1,jz))/dyy(jy) 

            COORDINATE = 'Z'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-3) = -2.0d0*RoAveLeft*harz(jx,jy,jz-1)/(dzz(jx,jy,jz)*(dzz(jx,jy,jz)+dzz(jx,jy,jz-1))) 
            coef(3)  = -2.0d0*RoAveRight*harz(jx,jy,jz)/(dzz(jx,jy,jz) *(dzz(jx,jy,jz+1)+dzz(jx,jy,jz) )) 
            BodyForceZ = -COSD(z_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harz(jx,jy,jz) - RoAveLeft*RoAveLeft*harz(jx,jy,jz-1))/dzz(jx,jy,jz) 
 
            tdepend = ScaleFactor*visc*ct*por(jx,jy,jz)/(dt*secyr)
            coef(0) = tdepend -( coef(-3) + coef(-2) + coef(-1) + coef(1) + coef(2) + coef(3) )

!!            CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j-nx*ny,coef(-3),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j+nx*ny,coef(3),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)

            tdepend = visc*ct*por(jx,jy,jz)*pres(jx,jy,jz)/ (dt*secyr)
            pumpterm = visc*ro(jx,jy,jz)*qg(jx,jy,jz)/(secyr*dxx(jx)*dyy(jy)*dzz(jx,jy,jz))

            IF (activecellPressure(0,jy,jz) == 0) THEN  
              AddPressureX =  2.0d0*ro(jx,jy,jz)*harx(jx-1,jy,jz)*pres(0,jy,jz)/(dxx(jx)*dxx(jx)) 
!! Check This              AddPressureX =  ro(jx,jy,jz)*harx(jx,jy,jz)*pres(0,jy,jz)/(dxx(jx)*dxx(jx))    
            ELSE
              AddPressureX = 0.0d0
            END IF 

            BvecCrunchP(j) = (pumpterm + tdepend + BodyForceX + BodyForceY + BodyForceZ + AddPressureX) 
            XvecCrunchP(j) = pres(jx,jy,jz)

          END IF
        END DO
      END DO

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
            BvecCrunchP(j) = pres(jx,jy,jz)*big 
            XvecCrunchP(j) = pres(jx,jy,jz)
          ELSE

          IF (activecellPressure(nx+1,jy,jz) == 0) THEN
            COORDINATE = 'X'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-1) = -2.0d0*RoAveLeft*harx(jx-1,jy,jz)/(dxx(jx)*(dxx(jx)+dxx(jx-1))) 
            coef(1) = -2.0d0*ro(jx,jy,jz)*harx(jx,jy,jz)/(dxx(jx)*dxx(jx)) 
            BodyForceX = -COSD(x_angle)*grav*SignGravity*(ro(jx,jy,jz)*ro(jx,jy,jz)*harx(jx,jy,jz) - RoAveLeft*RoAveLeft*harx(jx-1,jy,jz))/dxx(jx)
          ELSE
            COORDINATE = 'X'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-1) = -2.0d0*RoAveLeft*harx(jx-1,jy,jz)/(dxx(jx)*(dxx(jx)+dxx(jx-1))) 
            BodyForceX = -COSD(x_angle)*grav*SignGravity*(-RoAveLeft*RoAveLeft*harx(jx-1,jy,jz))/dxx(jx)
            coef(1) = 0.0d0
          END IF

            COORDINATE = 'Y'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-2) = -2.0d0*RoAveLeft*hary(jx,jy-1,jz)/(dyy(jy)*(dyy(jy)+dyy(jy-1))) 
            coef(2)  = -2.0d0*RoAveRight*hary(jx,jy,jz)/(dyy(jy) *(dyy(jy+1)+dyy(jy) )) 
            BodyForceY = -COSD(y_angle)*grav*SignGravity*(RoAveRight*RoAveRight*hary(jx,jy,jz) - RoAveLeft*RoAveLeft*hary(jx,jy-1,jz))/dyy(jy) 

            COORDINATE = 'Z'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-3) = -2.0d0*RoAveLeft*harz(jx,jy,jz-1)/(dzz(jx,jy,jz)*(dzz(jx,jy,jz)+dzz(jx,jy,jz-1))) 
            coef(3)  = -2.0d0*RoAveRight*harz(jx,jy,jz)/(dzz(jx,jy,jz) *(dzz(jx,jy,jz+1)+dzz(jx,jy,jz) )) 
            BodyForceZ = -COSD(z_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harz(jx,jy,jz) - RoAveLeft*RoAveLeft*harz(jx,jy,jz-1))/dzz(jx,jy,jz) 
 
            tdepend = ScaleFactor*visc*ct*por(jx,jy,jz)/(dt*secyr)
            coef(0) = tdepend -( coef(-3) + coef(-2) + coef(-1) + coef(1) + coef(2) + coef(3) )

            CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
!!            CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j-nx*ny,coef(-3),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j+nx*ny,coef(3),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)

            tdepend = visc*ct*por(jx,jy,jz)*pres(jx,jy,jz)/ (dt*secyr)
            pumpterm = visc*ro(jx,jy,jz)*qg(jx,jy,jz)/(secyr*dxx(jx)*dyy(jy)*dzz(jx,jy,jz))

            IF (activecellPressure(nx+1,jy,jz) == 0) THEN  
              AddPressureX =  2.0d0*ro(jx,jy,jz)*harx(jx,jy,jz)*pres(nx+1,jy,jz)/(dxx(jx)*dxx(jx))   
            ELSE
              AddPressureX = 0.0d0
            END IF 

            BvecCrunchP(j) = (pumpterm + tdepend + BodyForceX + BodyForceY + BodyForceZ + AddPressureX) 
            XvecCrunchP(j) = pres(jx,jy,jz)

          END IF
        END DO
      END DO

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
            BvecCrunchP(j) = pres(jx,jy,jz)*big 
            XvecCrunchP(j) = pres(jx,jy,jz)
          ELSE

            COORDINATE = 'X'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-1) = -2.0d0*RoAveLeft*harx(jx-1,jy,jz)/(dxx(jx)*(dxx(jx)+dxx(jx-1))) 
            coef(1) = -2.0d0*RoAveRight*harx(jx,jy,jz)/(dxx(jx)*(dxx(jx+1)+dxx(jx))) 
            BodyForceX = -COSD(x_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harx(jx,jy,jz) - RoAveLeft*RoAveLeft*harx(jx-1,jy,jz))/dxx(jx)

          IF (activecellPressure(jx,0,jz) == 0) THEN
            COORDINATE = 'Y'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(2)  = -2.0d0*RoAveRight*hary(jx,jy,jz)/(dyy(jy) *(dyy(jy+1)+dyy(jy) )) 
            coef(-2) = -2.0d0*ro(jx,jy,jz)*hary(jx,jy,jz)/(dyy(jy)*dyy(jy))
            BodyForceY = -COSD(y_angle)*grav*SignGravity*(RoAveRight*RoAveRight*hary(jx,jy,jz) - ro(jx,jy,jz)*ro(jx,jy,jz)*hary(jx,jy,jz))/dyy(jy) 
          ELSE
            COORDINATE = 'Y'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(2)  = -2.0d0*RoAveRight*hary(jx,jy,jz)/(dyy(jy) *(dyy(jy+1)+dyy(jy) )) 
            BodyForceY = -COSD(y_angle)*grav*SignGravity*(RoAveRight*RoAveRight*hary(jx,jy,jz))/dyy(jy) 
            coef(-2) = 0.0d0
          END IF

            COORDINATE = 'Z'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-3) = -2.0d0*RoAveLeft*harz(jx,jy,jz-1)/(dzz(jx,jy,jz)*(dzz(jx,jy,jz)+dzz(jx,jy,jz-1))) 
            coef(3)  = -2.0d0*RoAveRight*harz(jx,jy,jz)/(dzz(jx,jy,jz) *(dzz(jx,jy,jz+1)+dzz(jx,jy,jz) )) 
            BodyForceZ = -COSD(z_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harz(jx,jy,jz) - RoAveLeft*RoAveLeft*harz(jx,jy,jz-1))/dzz(jx,jy,jz) 
 
            tdepend = ScaleFactor*visc*ct*por(jx,jy,jz)/(dt*secyr)
            coef(0) = tdepend -( coef(-3) + coef(-2) + coef(-1) + coef(1) + coef(2) + coef(3) )

            CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
!!            CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j-nx*ny,coef(-3),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j+nx*ny,coef(3),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)

            tdepend = visc*ct*por(jx,jy,jz)*pres(jx,jy,jz)/ (dt*secyr)
            pumpterm = visc*ro(jx,jy,jz)*qg(jx,jy,jz)/(secyr*dxx(jx)*dyy(jy)*dzz(jx,jy,jz))

            IF (activecellPressure(jx,0,jz) == 0) THEN  
              AddPressureY =  2.0d0*ro(jx,jy,jz)*hary(jx,jy-1,jz)*pres(jx,0,jz)/(dyy(jy)*dyy(jy))   
            ELSE
              AddPressureY = 0.0d0
            END IF 
  
            BvecCrunchP(j) = (pumpterm + tdepend + BodyForceX + BodyForceY + BodyForceZ + AddPressureY) 
            XvecCrunchP(j) = pres(jx,jy,jz)

          END IF
        END DO
      END DO

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
            BvecCrunchP(j) = pres(jx,jy,jz)*big 
            XvecCrunchP(j) = pres(jx,jy,jz)
          ELSE

            COORDINATE = 'X'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-1) = -2.0d0*RoAveLeft*harx(jx-1,jy,jz)/(dxx(jx)*(dxx(jx)+dxx(jx-1))) 
            coef(1) = -2.0d0*RoAveRight*harx(jx,jy,jz)/(dxx(jx)*(dxx(jx+1)+dxx(jx))) 
            BodyForceX = -COSD(x_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harx(jx,jy,jz) - RoAveLeft*RoAveLeft*harx(jx-1,jy,jz))/dxx(jx)

          IF (activecellPressure(jx,ny+1,jz) == 0) THEN
            COORDINATE = 'Y'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-2) = -2.0d0*RoAveLeft*hary(jx,jy-1,jz)/(dyy(jy)*(dyy(jy)+dyy(jy-1))) 
            coef(2)  = -2.0d0*ro(jx,jy,jz)*hary(jx,jy,jz)/(dyy(jy)*dyy(jy)) 
            BodyForceY = -COSD(y_angle)*grav*SignGravity*(ro(jx,jy,jz)*ro(jx,jy,jz)*hary(jx,jy,jz) - RoAveLeft*RoAveLeft*hary(jx,jy-1,jz))/dyy(jy) 
          ELSE
            COORDINATE = 'Y'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-2) = -2.0d0*RoAveLeft*hary(jx,jy-1,jz)/(dyy(jy)*(dyy(jy)+dyy(jy-1))) 
            BodyForceY = -COSD(y_angle)*grav*SignGravity*(-RoAveLeft*RoAveLeft*hary(jx,jy-1,jz))/dyy(jy)    
            coef(2) = 0.0d0     
          END IF

            COORDINATE = 'Z'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-3) = -2.0d0*RoAveLeft*harz(jx,jy,jz-1)/(dzz(jx,jy,jz)*(dzz(jx,jy,jz)+dzz(jx,jy,jz-1))) 
            coef(3)  = -2.0d0*RoAveRight*harz(jx,jy,jz)/(dzz(jx,jy,jz) *(dzz(jx,jy,jz+1)+dzz(jx,jy,jz) )) 
            BodyForceZ = -COSD(z_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harz(jx,jy,jz) - RoAveLeft*RoAveLeft*harz(jx,jy,jz-1))/dzz(jx,jy,jz) 
 
            tdepend = ScaleFactor*visc*ct*por(jx,jy,jz)/(dt*secyr)
            coef(0) = tdepend -( coef(-3) + coef(-2) + coef(-1) + coef(1) + coef(2) + coef(3) )

            CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
!!            CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j-nx*ny,coef(-3),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j+nx*ny,coef(3),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)

            tdepend = visc*ct*por(jx,jy,jz)*pres(jx,jy,jz)/ (dt*secyr)
            pumpterm = visc*ro(jx,jy,jz)*qg(jx,jy,jz)/(secyr*dxx(jx)*dyy(jy)*dzz(jx,jy,jz))

            IF (activecellPressure(jx,ny+1,jz) == 0) THEN  
              AddPressureY =  2.0d0*ro(jx,jy,jz)*hary(jx,jy,jz)*pres(jx,ny+1,jz)/(dyy(jy)*dyy(jy))   
            ELSE
              AddPressureY = 0.0d0
            END IF 

            BvecCrunchP(j) = (pumpterm + tdepend + BodyForceX + BodyForceY + BodyForceZ + AddPressureY) 
            XvecCrunchP(j) = pres(jx,jy,jz)

          END IF
        END DO
      END DO

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
            BvecCrunchP(j) = pres(jx,jy,jz)*big 
            XvecCrunchP(j) = pres(jx,jy,jz)
          ELSE

            COORDINATE = 'X'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-1) = -2.0d0*RoAveLeft*harx(jx-1,jy,jz)/(dxx(jx)*(dxx(jx)+dxx(jx-1))) 
            coef(1) = -2.0d0*RoAveRight*harx(jx,jy,jz)/(dxx(jx)*(dxx(jx+1)+dxx(jx))) 
            BodyForceX = -COSD(x_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harx(jx,jy,jz) - RoAveLeft*RoAveLeft*harx(jx-1,jy,jz))/dxx(jx)

            COORDINATE = 'Y'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-2) = -2.0d0*RoAveLeft*hary(jx,jy-1,jz)/(dyy(jy)*(dyy(jy)+dyy(jy-1))) 
            coef(2)  = -2.0d0*RoAveRight*hary(jx,jy,jz)/(dyy(jy) *(dyy(jy+1)+dyy(jy) )) 
            BodyForceY = -COSD(y_angle)*grav*SignGravity*(RoAveRight*RoAveRight*hary(jx,jy,jz) - RoAveLeft*RoAveLeft*hary(jx,jy-1,jz))/dyy(jy) 

          IF (activecellPressure(jx,jy,0) == 0) THEN
            COORDINATE = 'Z'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-3) = -2.0d0*ro(jx,jy,jz)*harz(jx,jy,jz)/(dzz(jx,jy,jz)*dzz(jx,jy,jz))
            coef(3)  = -2.0d0*RoAveRight*harz(jx,jy,jz)/(dzz(jx,jy,jz) *(dzz(jx,jy,jz+1)+dzz(jx,jy,jz) )) 
            BodyForceZ = -COSD(z_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harz(jx,jy,jz) - ro(jx,jy,jz)*ro(jx,jy,jz)*harz(jx,jy,jz))/dzz(jx,jy,jz) 
          ELSE
            COORDINATE = 'Z'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(3)  = -2.0d0*RoAveRight*harz(jx,jy,jz)/(dzz(jx,jy,jz) *(dzz(jx,jy,jz+1)+dzz(jx,jy,jz) )) 
            BodyForceZ = -COSD(z_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harz(jx,jy,jz))/dzz(jx,jy,jz) 
            coef(-3) = 0.0d0
          END IF 
 
            tdepend = ScaleFactor*visc*ct*por(jx,jy,jz)/(dt*secyr)
            coef(0) = tdepend -( coef(-3) + coef(-2) + coef(-1) + coef(1) + coef(2) + coef(3) )

            CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
!!            CALL MatSetValues(amatP,1,j,1,j-nx*ny,coef(-3),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j+nx*ny,coef(3),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)

            tdepend = visc*ct*por(jx,jy,jz)*pres(jx,jy,jz)/ (dt*secyr)
            pumpterm = visc*ro(jx,jy,jz)*qg(jx,jy,jz)/(secyr*dxx(jx)*dyy(jy)*dzz(jx,jy,jz))

            IF (activecellPressure(jx,jy,0) == 0) THEN  
              AddPressureZ =  2.0d0*ro(jx,jy,jz)*harz(jx,jy,jz-1)*pres(jx,jy,0)/(dzz(jx,jy,jz)*dzz(jx,jy,jz))   
            ELSE
              AddPressureZ = 0.0d0
            END IF 

            BvecCrunchP(j) = (pumpterm + tdepend + BodyForceX + BodyForceY + BodyForceZ + AddPressureZ) 
            XvecCrunchP(j) = pres(jx,jy,jz)

          END IF
        END DO
      END DO

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
            BvecCrunchP(j) = pres(jx,jy,jz)*big 
            XvecCrunchP(j) = pres(jx,jy,jz)
          ELSE

            COORDINATE = 'X'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-1) = -2.0d0*RoAveLeft*harx(jx-1,jy,jz)/(dxx(jx)*(dxx(jx)+dxx(jx-1))) 
            coef(1) = -2.0d0*RoAveRight*harx(jx,jy,jz)/(dxx(jx)*(dxx(jx+1)+dxx(jx))) 
            BodyForceX = -COSD(x_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harx(jx,jy,jz) - RoAveLeft*RoAveLeft*harx(jx-1,jy,jz))/dxx(jx)

            COORDINATE = 'Y'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-2) = -2.0d0*RoAveLeft*hary(jx,jy-1,jz)/(dyy(jy)*(dyy(jy)+dyy(jy-1))) 
            coef(2)  = -2.0d0*RoAveRight*hary(jx,jy,jz)/(dyy(jy) *(dyy(jy+1)+dyy(jy) )) 
            BodyForceY = -COSD(y_angle)*grav*SignGravity*(RoAveRight*RoAveRight*hary(jx,jy,jz) - RoAveLeft*RoAveLeft*hary(jx,jy-1,jz))/dyy(jy) 

          IF (activecellPressure(jx,jy,nz+1) == 0) THEN
            COORDINATE = 'Z'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-3) = -2.0d0*RoAveLeft*harz(jx,jy,jz-1)/(dzz(jx,jy,jz)*(dzz(jx,jy,jz)+dzz(jx,jy,jz-1))) 
            coef(3) = -2.0d0*ro(jx,jy,jz)*harz(jx,jy,jz)/(dzz(jx,jy,jz)*dzz(jx,jy,jz)) 
            BodyForceZ = -COSD(z_angle)*grav*SignGravity*(ro(jx,jy,jz)*ro(jx,jy,jz)*harz(jx,jy,jz) - RoAveLeft*RoAveLeft*harz(jx,jy,jz-1))/dzz(jx,jy,jz)
          ELSE
            COORDINATE = 'Z'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-3) = -2.0d0*RoAveLeft*harz(jx,jy,jz-1)/(dzz(jx,jy,jz)*(dzz(jx,jy,jz)+dzz(jx,jy,jz-1))) 
            BodyForceZ = -COSD(z_angle)*grav*SignGravity*(-RoAveLeft*RoAveLeft*harz(jx,jy,jz-1))/dzz(jx,jy,jz)
            coef(3) = 0.0d0
          END IF
 
            tdepend = ScaleFactor*visc*ct*por(jx,jy,jz)/(dt*secyr)
            coef(0) = tdepend -( coef(-3) + coef(-2) + coef(-1) + coef(1) + coef(2) + coef(3) )

            CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j-nx*ny,coef(-3),INSERT_VALUES,ierr)
!!            CALL MatSetValues(amatP,1,j,1,j+nx*ny,coef(3),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)

            tdepend = visc*ct*por(jx,jy,jz)*pres(jx,jy,jz)/ (dt*secyr)
            pumpterm = visc*ro(jx,jy,jz)*qg(jx,jy,jz)/(secyr*dxx(jx)*dyy(jy)*dzz(jx,jy,jz))

            IF (activecellPressure(jx,jy,nz+1) == 0) THEN  
              AddPressureZ =  2.0d0*ro(jx,jy,jz)*harz(jx,jy,jz)*pres(jx,jy,nz+1)/(dzz(jx,jy,jz)*dzz(jx,jy,jz))   
            ELSE
              AddPressureZ = 0.0d0
            END IF 

            BvecCrunchP(j) = (pumpterm + tdepend + BodyForceX + BodyForceY + BodyForceZ + AddPressureZ) 
            XvecCrunchP(j) = pres(jx,jy,jz)

          END IF
        END DO
      END DO

!!  ************************************************************************************************************************************************
!!  Now the edges

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
            BvecCrunchP(j) = pres(jx,jy,jz)*big 
            XvecCrunchP(j) = pres(jx,jy,jz)
          ELSE

          IF (activecellPressure(0,jy,jz) == 0) THEN
            COORDINATE = 'X'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-1) = -2.0d0*ro(jx,jy,jz)*harx(jx,jy,jz)/(dxx(jx)*dxx(jx))
            coef(1) = -2.0d0*RoAveRight*harx(jx,jy,jz)/(dxx(jx)*(dxx(jx+1)+dxx(jx))) 
            BodyForceX = -COSD(x_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harx(jx,jy,jz) - ro(jx,jy,jz)*ro(jx,jy,jz)*harx(jx,jy,jz))/dxx(jx)
          ELSE
            COORDINATE = 'X'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(1) = -2.0d0*RoAveRight*harx(jx,jy,jz)/(dxx(jx)*(dxx(jx+1)+dxx(jx))) 
            BodyForceX = -COSD(x_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harx(jx,jy,jz) )/dxx(jx)
            coef(-1) = 0.0d0
          END IF

          IF (activecellPressure(jx,0,jz) == 0) THEN
            COORDINATE = 'Y'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(2)  = -2.0d0*RoAveRight*hary(jx,jy,jz)/(dyy(jy) *(dyy(jy+1)+dyy(jy) )) 
            coef(-2) = -2.0d0*ro(jx,jy,jz)*hary(jx,jy,jz)/(dyy(jy)*dyy(jy))
            BodyForceY = -COSD(y_angle)*grav*SignGravity*(RoAveRight*RoAveRight*hary(jx,jy,jz) - ro(jx,jy,jz)*ro(jx,jy,jz)*hary(jx,jy,jz))/dyy(jy) 
          ELSE
            COORDINATE = 'Y'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(2)  = -2.0d0*RoAveRight*hary(jx,jy,jz)/(dyy(jy) *(dyy(jy+1)+dyy(jy) )) 
            BodyForceY = -COSD(y_angle)*grav*SignGravity*(RoAveRight*RoAveRight*hary(jx,jy,jz))/dyy(jy) 
            coef(-2) = 0.0d0
          END IF
            COORDINATE = 'Z'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-3) = -2.0d0*RoAveLeft*harz(jx,jy,jz-1)/(dzz(jx,jy,jz)*(dzz(jx,jy,jz)+dzz(jx,jy,jz-1))) 
            coef(3)  = -2.0d0*RoAveRight*harz(jx,jy,jz)/(dzz(jx,jy,jz) *(dzz(jx,jy,jz+1)+dzz(jx,jy,jz) )) 
            BodyForceZ = -COSD(z_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harz(jx,jy,jz) - RoAveLeft*RoAveLeft*harz(jx,jy,jz-1))/dzz(jx,jy,jz) 
 
            tdepend = ScaleFactor*visc*ct*por(jx,jy,jz)/(dt*secyr)
            coef(0) = tdepend -( coef(-3) + coef(-2) + coef(-1) + coef(1) + coef(2) + coef(3) )

!!            CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
!!            CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j-nx*ny,coef(-3),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j+nx*ny,coef(3),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)

            tdepend = visc*ct*por(jx,jy,jz)*pres(jx,jy,jz)/ (dt*secyr)
            pumpterm = visc*ro(jx,jy,jz)*qg(jx,jy,jz)/(secyr*dxx(jx)*dyy(jy)*dzz(jx,jy,jz))

          IF (activecellPressure(0,jy,jz) == 0) THEN  
            AddPressureX =  2.0d0*ro(jx,jy,jz)*harx(jx-1,jy,jz)*pres(0,jy,jz)/(dxx(jx)*dxx(jx))   
          ELSE
            AddPressureX = 0.0d0
          END IF 
          IF (activecellPressure(jx,0,jz) == 0) THEN  
            AddPressureY =  2.0d0*ro(jx,jy,jz)*hary(jx,jy-1,jz)*pres(jx,0,jz)/(dyy(jy)*dyy(jy))   
          ELSE
            AddPressureY = 0.0d0
          END IF   

            BvecCrunchP(j) = (pumpterm + tdepend + BodyForceX + BodyForceY + BodyForceZ + AddPressureX + AddPressureY) 
            XvecCrunchP(j) = pres(jx,jy,jz)

          END IF
        END DO

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
            BvecCrunchP(j) = pres(jx,jy,jz)*big 
            XvecCrunchP(j) = pres(jx,jy,jz)
          ELSE

          IF (activecellPressure(nx+1,jy,jz) == 0) THEN
            COORDINATE = 'X'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-1) = -2.0d0*RoAveLeft*harx(jx-1,jy,jz)/(dxx(jx)*(dxx(jx)+dxx(jx-1))) 
            coef(1) = -2.0d0*ro(jx,jy,jz)*harx(jx,jy,jz)/(dxx(jx)*dxx(jx)) 
            BodyForceX = -COSD(x_angle)*grav*SignGravity*(ro(jx,jy,jz)*ro(jx,jy,jz)*harx(jx,jy,jz) - RoAveLeft*RoAveLeft*harx(jx-1,jy,jz))/dxx(jx)
          ELSE
            COORDINATE = 'X'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-1) = -2.0d0*RoAveLeft*harx(jx-1,jy,jz)/(dxx(jx)*(dxx(jx)+dxx(jx-1))) 
            BodyForceX = -COSD(x_angle)*grav*SignGravity*(-RoAveLeft*RoAveLeft*harx(jx-1,jy,jz))/dxx(jx)
            coef(1) = 0.0d0
          END IF

          IF (activecellPressure(jx,0,jz) == 0) THEN
            COORDINATE = 'Y'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(2)  = -2.0d0*RoAveRight*hary(jx,jy,jz)/(dyy(jy) *(dyy(jy+1)+dyy(jy) )) 
            coef(-2) = -2.0d0*ro(jx,jy,jz)*hary(jx,jy,jz)/(dyy(jy)*dyy(jy))
            BodyForceY = -COSD(y_angle)*grav*SignGravity*(RoAveRight*RoAveRight*hary(jx,jy,jz) - ro(jx,jy,jz)*ro(jx,jy,jz)*hary(jx,jy,jz))/dyy(jy) 
          ELSE
            COORDINATE = 'Y'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(2)  = -2.0d0*RoAveRight*hary(jx,jy,jz)/(dyy(jy) *(dyy(jy+1)+dyy(jy) )) 
            BodyForceY = -COSD(y_angle)*grav*SignGravity*(RoAveRight*RoAveRight*hary(jx,jy,jz))/dyy(jy) 
            coef(-2) = 0.0d0
          END IF

            COORDINATE = 'Z'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-3) = -2.0d0*RoAveLeft*harz(jx,jy,jz-1)/(dzz(jx,jy,jz)*(dzz(jx,jy,jz)+dzz(jx,jy,jz-1))) 
            coef(3)  = -2.0d0*RoAveRight*harz(jx,jy,jz)/(dzz(jx,jy,jz) *(dzz(jx,jy,jz+1)+dzz(jx,jy,jz) )) 
            BodyForceZ = -COSD(z_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harz(jx,jy,jz) - RoAveLeft*RoAveLeft*harz(jx,jy,jz-1))/dzz(jx,jy,jz) 
 
            tdepend = ScaleFactor*visc*ct*por(jx,jy,jz)/(dt*secyr)
            coef(0) = tdepend -( coef(-3) + coef(-2) + coef(-1) + coef(1) + coef(2) + coef(3) )

            CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
!!            CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
!!            CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j-nx*ny,coef(-3),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j+nx*ny,coef(3),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)

            tdepend = visc*ct*por(jx,jy,jz)*pres(jx,jy,jz)/ (dt*secyr)
            pumpterm = visc*ro(jx,jy,jz)*qg(jx,jy,jz)/(secyr*dxx(jx)*dyy(jy)*dzz(jx,jy,jz))

          IF (activecellPressure(nx+1,jy,jz) == 0) THEN  
            AddPressureX =  2.0d0*ro(jx,jy,jz)*harx(jx,jy,jz)*pres(nx+1,jy,jz)/(dxx(jx)*dxx(jx))   
          ELSE
            AddPressureX = 0.0d0
          END IF 
          IF (activecellPressure(jx,0,jz) == 0) THEN  
            AddPressureY =  2.0d0*ro(jx,jy,jz)*hary(jx,jy-1,jz)*pres(jx,0,jz)/(dyy(jy)*dyy(jy))   
          ELSE
            AddPressureY = 0.0d0
          END IF 

            BvecCrunchP(j) = (pumpterm + tdepend + BodyForceX + BodyForceY + BodyForceZ + AddPressureX + AddPressureY) 
            XvecCrunchP(j) = pres(jx,jy,jz)

          END IF
        END DO

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
            BvecCrunchP(j) = pres(jx,jy,jz)*big 
            XvecCrunchP(j) = pres(jx,jy,jz)
          ELSE

         IF (activecellPressure(nx+1,jy,jz) == 0) THEN
            COORDINATE = 'X'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-1) = -2.0d0*RoAveLeft*harx(jx-1,jy,jz)/(dxx(jx)*(dxx(jx)+dxx(jx-1))) 
            coef(1) = -2.0d0*ro(jx,jy,jz)*harx(jx,jy,jz)/(dxx(jx)*dxx(jx)) 
            BodyForceX = -COSD(x_angle)*grav*SignGravity*(ro(jx,jy,jz)*ro(jx,jy,jz)*harx(jx,jy,jz) - RoAveLeft*RoAveLeft*harx(jx-1,jy,jz))/dxx(jx)
          ELSE
            COORDINATE = 'X'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-1) = -2.0d0*RoAveLeft*harx(jx-1,jy,jz)/(dxx(jx)*(dxx(jx)+dxx(jx-1))) 
            BodyForceX = -COSD(x_angle)*grav*SignGravity*(-RoAveLeft*RoAveLeft*harx(jx-1,jy,jz))/dxx(jx)
            coef(1) = 0.0d0
          END IF

          IF (activecellPressure(jx,ny+1,jz) == 0) THEN
            COORDINATE = 'Y'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-2) = -2.0d0*RoAveLeft*hary(jx,jy-1,jz)/(dyy(jy)*(dyy(jy)+dyy(jy-1))) 
            coef(2)  = -2.0d0*ro(jx,jy,jz)*hary(jx,jy,jz)/(dyy(jy)*dyy(jy)) 
            BodyForceY = -COSD(y_angle)*grav*SignGravity*(ro(jx,jy,jz)*ro(jx,jy,jz)*hary(jx,jy,jz) - RoAveLeft*RoAveLeft*hary(jx,jy-1,jz))/dyy(jy) 
          ELSE
            COORDINATE = 'Y'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-2) = -2.0d0*RoAveLeft*hary(jx,jy-1,jz)/(dyy(jy)*(dyy(jy)+dyy(jy-1))) 
            BodyForceY = -COSD(y_angle)*grav*SignGravity*(-RoAveLeft*RoAveLeft*hary(jx,jy-1,jz))/dyy(jy)    
            coef(2) = 0.0d0     
          END IF

            COORDINATE = 'Z'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-3) = -2.0d0*RoAveLeft*harz(jx,jy,jz-1)/(dzz(jx,jy,jz)*(dzz(jx,jy,jz)+dzz(jx,jy,jz-1))) 
            coef(3)  = -2.0d0*RoAveRight*harz(jx,jy,jz)/(dzz(jx,jy,jz) *(dzz(jx,jy,jz+1)+dzz(jx,jy,jz) )) 
            BodyForceZ = -COSD(z_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harz(jx,jy,jz) - RoAveLeft*RoAveLeft*harz(jx,jy,jz-1))/dzz(jx,jy,jz) 
 
            tdepend = ScaleFactor*visc*ct*por(jx,jy,jz)/(dt*secyr)
            coef(0) = tdepend -( coef(-3) + coef(-2) + coef(-1) + coef(1) + coef(2) + coef(3) )

            CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
!!            CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
!!            CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j-nx*ny,coef(-3),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j+nx*ny,coef(3),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)

            tdepend = visc*ct*por(jx,jy,jz)*pres(jx,jy,jz)/ (dt*secyr)
            pumpterm = visc*ro(jx,jy,jz)*qg(jx,jy,jz)/(secyr*dxx(jx)*dyy(jy)*dzz(jx,jy,jz))

          IF (activecellPressure(nx+1,jy,jz) == 0) THEN  
            AddPressureX =  2.0d0*ro(jx,jy,jz)*harx(jx,jy,jz)*pres(nx+1,jy,jz)/(dxx(jx)*dxx(jx))  
          ELSE
            AddPressureX = 0.0d0
          END IF 
          IF (activecellPressure(jx,ny+1,jz) == 0) THEN  
            AddPressureY =  2.0d0*ro(jx,jy,jz)*hary(jx,jy,jz)*pres(jx,ny+1,jz)/(dyy(jy)*dyy(jy))   
          ELSE
            AddPressureY = 0.0d0
          END IF 

            BvecCrunchP(j) = (pumpterm + tdepend + BodyForceX + BodyForceY + BodyForceZ + AddPressureX + AddPressureY) 
            XvecCrunchP(j) = pres(jx,jy,jz)

          END IF
        END DO

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
            BvecCrunchP(j) = pres(jx,jy,jz)*big 
            XvecCrunchP(j) = pres(jx,jy,jz)
          ELSE

          IF (activecellPressure(0,jy,jz) == 0) THEN
            COORDINATE = 'X'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-1) = -2.0d0*ro(jx,jy,jz)*harx(jx,jy,jz)/(dxx(jx)*dxx(jx))
            coef(1) = -2.0d0*RoAveRight*harx(jx,jy,jz)/(dxx(jx)*(dxx(jx+1)+dxx(jx))) 
            BodyForceX = -COSD(x_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harx(jx,jy,jz) - ro(jx,jy,jz)*ro(jx,jy,jz)*harx(jx,jy,jz))/dxx(jx)
          ELSE
            COORDINATE = 'X'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(1) = -2.0d0*RoAveRight*harx(jx,jy,jz)/(dxx(jx)*(dxx(jx+1)+dxx(jx))) 
            BodyForceX = -COSD(x_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harx(jx,jy,jz) )/dxx(jx)
            coef(-1) = 0.0d0
          END IF

          IF (activecellPressure(jx,ny+1,jz) == 0) THEN
            COORDINATE = 'Y'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-2) = -2.0d0*RoAveLeft*hary(jx,jy-1,jz)/(dyy(jy)*(dyy(jy)+dyy(jy-1))) 
            coef(2)  = -2.0d0*ro(jx,jy,jz)*hary(jx,jy,jz)/(dyy(jy)*dyy(jy)) 
            BodyForceY = -COSD(y_angle)*grav*SignGravity*(ro(jx,jy,jz)*ro(jx,jy,jz)*hary(jx,jy,jz) - RoAveLeft*RoAveLeft*hary(jx,jy-1,jz))/dyy(jy) 
          ELSE
            COORDINATE = 'Y'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-2) = -2.0d0*RoAveLeft*hary(jx,jy-1,jz)/(dyy(jy)*(dyy(jy)+dyy(jy-1))) 
            BodyForceY = -COSD(y_angle)*grav*SignGravity*(-RoAveLeft*RoAveLeft*hary(jx,jy-1,jz))/dyy(jy)    
            coef(2) = 0.0d0     
          END IF

            COORDINATE = 'Z'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-3) = -2.0d0*RoAveLeft*harz(jx,jy,jz-1)/(dzz(jx,jy,jz)*(dzz(jx,jy,jz)+dzz(jx,jy,jz-1))) 
            coef(3)  = -2.0d0*RoAveRight*harz(jx,jy,jz)/(dzz(jx,jy,jz) *(dzz(jx,jy,jz+1)+dzz(jx,jy,jz) )) 
            BodyForceZ = -COSD(z_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harz(jx,jy,jz) - RoAveLeft*RoAveLeft*harz(jx,jy,jz-1))/dzz(jx,jy,jz) 
 
            tdepend = ScaleFactor*visc*ct*por(jx,jy,jz)/(dt*secyr)
            coef(0) = tdepend -( coef(-3) + coef(-2) + coef(-1) + coef(1) + coef(2) + coef(3) )

!!            CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
!!            CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j-nx*ny,coef(-3),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j+nx*ny,coef(3),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)

            tdepend = visc*ct*por(jx,jy,jz)*pres(jx,jy,jz)/ (dt*secyr)
            pumpterm = visc*ro(jx,jy,jz)*qg(jx,jy,jz)/(secyr*dxx(jx)*dyy(jy)*dzz(jx,jy,jz))

          IF (activecellPressure(0,jy,jz) == 0) THEN  
            AddPressureX =  2.0d0*ro(jx,jy,jz)*harx(jx-1,jy,jz)*pres(0,jy,jz)/(dxx(jx)*dxx(jx))   
          ELSE
            AddPressureX = 0.0d0
          END IF 
          IF (activecellPressure(jx,ny+1,jz) == 0) THEN  
            AddPressureY =  2.0d0*ro(jx,jy,jz)*hary(jx,jy,jz)*pres(jx,ny+1,jz)/(dyy(jy)*dyy(jy))   
          ELSE
            AddPressureY = 0.0d0
          END IF 

            BvecCrunchP(j) = (pumpterm + tdepend + BodyForceX + BodyForceY + BodyForceZ + AddPressureX + AddPressureY) 
            XvecCrunchP(j) = pres(jx,jy,jz)

          END IF
        END DO


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
            BvecCrunchP(j) = pres(jx,jy,jz)*big 
            XvecCrunchP(j) = pres(jx,jy,jz)
          ELSE

          IF (activecellPressure(0,jy,jz) == 0) THEN
            COORDINATE = 'X'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-1) = -2.0d0*ro(jx,jy,jz)*harx(jx,jy,jz)/(dxx(jx)*dxx(jx))
            coef(1) = -2.0d0*RoAveRight*harx(jx,jy,jz)/(dxx(jx)*(dxx(jx+1)+dxx(jx))) 
            BodyForceX = -COSD(x_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harx(jx,jy,jz) - ro(jx,jy,jz)*ro(jx,jy,jz)*harx(jx,jy,jz))/dxx(jx)
          ELSE
            COORDINATE = 'X'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(1) = -2.0d0*RoAveRight*harx(jx,jy,jz)/(dxx(jx)*(dxx(jx+1)+dxx(jx))) 
            BodyForceX = -COSD(x_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harx(jx,jy,jz) )/dxx(jx)
            coef(-1) = 0.0d0
          END IF

            COORDINATE = 'Y'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-2) = -2.0d0*RoAveLeft*hary(jx,jy-1,jz)/(dyy(jy)*(dyy(jy)+dyy(jy-1))) 
            coef(2)  = -2.0d0*RoAveRight*hary(jx,jy,jz)/(dyy(jy) *(dyy(jy+1)+dyy(jy) )) 
            BodyForceY = -COSD(y_angle)*grav*SignGravity*(RoAveRight*RoAveRight*hary(jx,jy,jz) - RoAveLeft*RoAveLeft*hary(jx,jy-1,jz))/dyy(jy) 

          IF (activecellPressure(jx,jy,0) == 0) THEN
            COORDINATE = 'Z'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-3) = -2.0d0*ro(jx,jy,jz)*harz(jx,jy,jz)/(dzz(jx,jy,jz)*dzz(jx,jy,jz))
            coef(3)  = -2.0d0*RoAveRight*harz(jx,jy,jz)/(dzz(jx,jy,jz) *(dzz(jx,jy,jz+1)+dzz(jx,jy,jz) )) 
            BodyForceZ = -COSD(z_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harz(jx,jy,jz) - ro(jx,jy,jz)*ro(jx,jy,jz)*harz(jx,jy,jz))/dzz(jx,jy,jz) 
          ELSE
            COORDINATE = 'Z'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(3)  = -2.0d0*RoAveRight*harz(jx,jy,jz)/(dzz(jx,jy,jz) *(dzz(jx,jy,jz+1)+dzz(jx,jy,jz) )) 
            BodyForceZ = -COSD(z_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harz(jx,jy,jz))/dzz(jx,jy,jz) 
            coef(-3) = 0.0d0
          END IF
 
            tdepend = ScaleFactor*visc*ct*por(jx,jy,jz)/(dt*secyr)
            coef(0) = tdepend -( coef(-3) + coef(-2) + coef(-1) + coef(1) + coef(2) + coef(3) )

!!            CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
!!            CALL MatSetValues(amatP,1,j,1,j-nx*ny,coef(-3),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j+nx*ny,coef(3),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)

            tdepend = visc*ct*por(jx,jy,jz)*pres(jx,jy,jz)/ (dt*secyr)
            pumpterm = visc*ro(jx,jy,jz)*qg(jx,jy,jz)/(secyr*dxx(jx)*dyy(jy)*dzz(jx,jy,jz))

          IF (activecellPressure(0,jy,jz) == 0) THEN  
            AddPressureX =  2.0d0*ro(jx,jy,jz)*harx(jx-1,jy,jz)*pres(0,jy,jz)/(dxx(jx)*dxx(jx))   
          ELSE
            AddPressureX = 0.0d0
          END IF 
          IF (activecellPressure(jx,jy,0) == 0) THEN  
            AddPressureZ =  2.0d0*ro(jx,jy,jz)*harz(jx,jy,jz-1)*pres(jx,jy,0)/(dzz(jx,jy,jz)*dzz(jx,jy,jz))   
          ELSE
            AddPressureZ = 0.0d0
          END IF 

            BvecCrunchP(j) = (pumpterm + tdepend + BodyForceX + BodyForceY + BodyForceZ + AddPressureX + AddPressureZ) 
            XvecCrunchP(j) = pres(jx,jy,jz)

          END IF
        END DO

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
            BvecCrunchP(j) = pres(jx,jy,jz)*big 
            XvecCrunchP(j) = pres(jx,jy,jz)
          ELSE

          IF (activecellPressure(nx+1,jy,jz) == 0) THEN
            COORDINATE = 'X'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-1) = -2.0d0*RoAveLeft*harx(jx-1,jy,jz)/(dxx(jx)*(dxx(jx)+dxx(jx-1))) 
            coef(1) = -2.0d0*ro(jx,jy,jz)*harx(jx,jy,jz)/(dxx(jx)*dxx(jx)) 
            BodyForceX = -COSD(x_angle)*grav*SignGravity*(ro(jx,jy,jz)*ro(jx,jy,jz)*harx(jx,jy,jz) - RoAveLeft*RoAveLeft*harx(jx-1,jy,jz))/dxx(jx)
          ELSE
            COORDINATE = 'X'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-1) = -2.0d0*RoAveLeft*harx(jx-1,jy,jz)/(dxx(jx)*(dxx(jx)+dxx(jx-1))) 
            BodyForceX = -COSD(x_angle)*grav*SignGravity*(-RoAveLeft*RoAveLeft*harx(jx-1,jy,jz))/dxx(jx)
            coef(1) = 0.0d0
          END IF

            COORDINATE = 'Y'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-2) = -2.0d0*RoAveLeft*hary(jx,jy-1,jz)/(dyy(jy)*(dyy(jy)+dyy(jy-1))) 
            coef(2)  = -2.0d0*RoAveRight*hary(jx,jy,jz)/(dyy(jy) *(dyy(jy+1)+dyy(jy) )) 
            BodyForceY = -COSD(y_angle)*grav*SignGravity*(RoAveRight*RoAveRight*hary(jx,jy,jz) - RoAveLeft*RoAveLeft*hary(jx,jy-1,jz))/dyy(jy) 

          IF (activecellPressure(jx,jy,0) == 0) THEN
            COORDINATE = 'Z'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-3) = -2.0d0*ro(jx,jy,jz)*harz(jx,jy,jz)/(dzz(jx,jy,jz)*dzz(jx,jy,jz))
            coef(3)  = -2.0d0*RoAveRight*harz(jx,jy,jz)/(dzz(jx,jy,jz) *(dzz(jx,jy,jz+1)+dzz(jx,jy,jz) )) 
            BodyForceZ = -COSD(z_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harz(jx,jy,jz) - ro(jx,jy,jz)*ro(jx,jy,jz)*harz(jx,jy,jz))/dzz(jx,jy,jz) 
          ELSE
            COORDINATE = 'Z'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(3)  = -2.0d0*RoAveRight*harz(jx,jy,jz)/(dzz(jx,jy,jz) *(dzz(jx,jy,jz+1)+dzz(jx,jy,jz) )) 
            BodyForceZ = -COSD(z_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harz(jx,jy,jz))/dzz(jx,jy,jz) 
            coef(-3) = 0.0d0
          END IF
 
            tdepend = ScaleFactor*visc*ct*por(jx,jy,jz)/(dt*secyr)
            coef(0) = tdepend -( coef(-3) + coef(-2) + coef(-1) + coef(1) + coef(2) + coef(3) )

            CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
!!            CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
!!            CALL MatSetValues(amatP,1,j,1,j-nx*ny,coef(-3),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j+nx*ny,coef(3),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)

            tdepend = visc*ct*por(jx,jy,jz)*pres(jx,jy,jz)/ (dt*secyr)
            pumpterm = visc*ro(jx,jy,jz)*qg(jx,jy,jz)/(secyr*dxx(jx)*dyy(jy)*dzz(jx,jy,jz))

          IF (activecellPressure(nx+1,jy,jz) == 0) THEN  
            AddPressureX =  2.0d0*ro(jx,jy,jz)*harx(jx,jy,jz)*pres(nx+1,jy,jz)/(dxx(jx)*dxx(jx))   
          ELSE
            AddPressureX = 0.0d0
          END IF 
          IF (activecellPressure(jx,jy,0) == 0) THEN  
            AddPressureZ =  2.0d0*ro(jx,jy,jz)*harz(jx,jy,jz-1)*pres(jx,jy,0)/(dzz(jx,jy,jz)*dzz(jx,jy,jz))   
          ELSE
            AddPressureZ = 0.0d0
          END IF 

            BvecCrunchP(j) = (pumpterm + tdepend + BodyForceX + BodyForceY + BodyForceZ + AddPressureX + AddPressureZ) 
            XvecCrunchP(j) = pres(jx,jy,jz)

          END IF
        END DO

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
            BvecCrunchP(j) = pres(jx,jy,jz)*big 
            XvecCrunchP(j) = pres(jx,jy,jz)
          ELSE

          IF (activecellPressure(0,jy,jz) == 0) THEN
            COORDINATE = 'X'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-1) = -2.0d0*ro(jx,jy,jz)*harx(jx,jy,jz)/(dxx(jx)*dxx(jx))
            coef(1) = -2.0d0*RoAveRight*harx(jx,jy,jz)/(dxx(jx)*(dxx(jx+1)+dxx(jx))) 
            BodyForceX = -COSD(x_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harx(jx,jy,jz) - ro(jx,jy,jz)*ro(jx,jy,jz)*harx(jx,jy,jz))/dxx(jx)
          ELSE
            COORDINATE = 'X'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(1) = -2.0d0*RoAveRight*harx(jx,jy,jz)/(dxx(jx)*(dxx(jx+1)+dxx(jx))) 
            BodyForceX = -COSD(x_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harx(jx,jy,jz) )/dxx(jx)
            coef(-1) = 0.0d0
          END IF

            COORDINATE = 'Y'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-2) = -2.0d0*RoAveLeft*hary(jx,jy-1,jz)/(dyy(jy)*(dyy(jy)+dyy(jy-1))) 
            coef(2)  = -2.0d0*RoAveRight*hary(jx,jy,jz)/(dyy(jy) *(dyy(jy+1)+dyy(jy) )) 
            BodyForceY = -COSD(y_angle)*grav*SignGravity*(RoAveRight*RoAveRight*hary(jx,jy,jz) - RoAveLeft*RoAveLeft*hary(jx,jy-1,jz))/dyy(jy) 

          IF (activecellPressure(jx,jy,nz+1) == 0) THEN
            COORDINATE = 'Z'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-3) = -2.0d0*RoAveLeft*harz(jx,jy,jz-1)/(dzz(jx,jy,jz)*(dzz(jx,jy,jz)+dzz(jx,jy,jz-1))) 
            coef(3) = -2.0d0*ro(jx,jy,jz)*harz(jx,jy,jz)/(dzz(jx,jy,jz)*dzz(jx,jy,jz)) 
            BodyForceZ = -COSD(z_angle)*grav*SignGravity*(ro(jx,jy,jz)*ro(jx,jy,jz)*harz(jx,jy,jz) - RoAveLeft*RoAveLeft*harz(jx,jy,jz-1))/dzz(jx,jy,jz)
          ELSE
            COORDINATE = 'Z'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-3) = -2.0d0*RoAveLeft*harz(jx,jy,jz-1)/(dzz(jx,jy,jz)*(dzz(jx,jy,jz)+dzz(jx,jy,jz-1))) 
            BodyForceZ = -COSD(z_angle)*grav*SignGravity*(-RoAveLeft*RoAveLeft*harz(jx,jy,jz-1))/dzz(jx,jy,jz)
            coef(3) = 0.0d0
          END IF
 
            tdepend = ScaleFactor*visc*ct*por(jx,jy,jz)/(dt*secyr)
            coef(0) = tdepend -( coef(-3) + coef(-2) + coef(-1) + coef(1) + coef(2) + coef(3) )

!!            CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j-nx*ny,coef(-3),INSERT_VALUES,ierr)
!!            CALL MatSetValues(amatP,1,j,1,j+nx*ny,coef(3),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)

            tdepend = visc*ct*por(jx,jy,jz)*pres(jx,jy,jz)/ (dt*secyr)
            pumpterm = visc*ro(jx,jy,jz)*qg(jx,jy,jz)/(secyr*dxx(jx)*dyy(jy)*dzz(jx,jy,jz))

          IF (activecellPressure(0,jy,jz) == 0) THEN  
            AddPressureX =  2.0d0*ro(jx,jy,jz)*harx(jx-1,jy,jz)*pres(0,jy,jz)/(dxx(jx)*dxx(jx))   
          ELSE
            AddPressureX = 0.0d0
          END IF 
          IF (activecellPressure(jx,jy,nz+1) == 0) THEN  
            AddPressureZ =  2.0d0*ro(jx,jy,jz)*harz(jx,jy,jz)*pres(jx,jy,nz+1)/(dzz(jx,jy,jz)*dzz(jx,jy,jz))   
          ELSE
            AddPressureZ = 0.0d0
          END IF 

            BvecCrunchP(j) = (pumpterm + tdepend + BodyForceX + BodyForceY + BodyForceZ + AddPressureX + AddPressureZ) 
            XvecCrunchP(j) = pres(jx,jy,jz)

          END IF
        END DO

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
            BvecCrunchP(j) = pres(jx,jy,jz)*big 
            XvecCrunchP(j) = pres(jx,jy,jz)
          ELSE

          IF (activecellPressure(nx+1,jy,jz) == 0) THEN
            COORDINATE = 'X'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-1) = -2.0d0*RoAveLeft*harx(jx-1,jy,jz)/(dxx(jx)*(dxx(jx)+dxx(jx-1))) 
            coef(1) = -2.0d0*ro(jx,jy,jz)*harx(jx,jy,jz)/(dxx(jx)*dxx(jx)) 
            BodyForceX = -COSD(x_angle)*grav*SignGravity*(ro(jx,jy,jz)*ro(jx,jy,jz)*harx(jx,jy,jz) - RoAveLeft*RoAveLeft*harx(jx-1,jy,jz))/dxx(jx)
          ELSE
            COORDINATE = 'X'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-1) = -2.0d0*RoAveLeft*harx(jx-1,jy,jz)/(dxx(jx)*(dxx(jx)+dxx(jx-1))) 
            BodyForceX = -COSD(x_angle)*grav*SignGravity*(-RoAveLeft*RoAveLeft*harx(jx-1,jy,jz))/dxx(jx)
            coef(1) = 0.0d0
          END IF

            COORDINATE = 'Y'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-2) = -2.0d0*RoAveLeft*hary(jx,jy-1,jz)/(dyy(jy)*(dyy(jy)+dyy(jy-1))) 
            coef(2)  = -2.0d0*RoAveRight*hary(jx,jy,jz)/(dyy(jy) *(dyy(jy+1)+dyy(jy) )) 
            BodyForceY = -COSD(y_angle)*grav*SignGravity*(RoAveRight*RoAveRight*hary(jx,jy,jz) - RoAveLeft*RoAveLeft*hary(jx,jy-1,jz))/dyy(jy) 

          IF (activecellPressure(jx,jy,nz+1) == 0) THEN
            COORDINATE = 'Z'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-3) = -2.0d0*RoAveLeft*harz(jx,jy,jz-1)/(dzz(jx,jy,jz)*(dzz(jx,jy,jz)+dzz(jx,jy,jz-1))) 
            coef(3) = -2.0d0*ro(jx,jy,jz)*harz(jx,jy,jz)/(dzz(jx,jy,jz)*dzz(jx,jy,jz)) 
            BodyForceZ = -COSD(z_angle)*grav*SignGravity*(ro(jx,jy,jz)*ro(jx,jy,jz)*harz(jx,jy,jz) - RoAveLeft*RoAveLeft*harz(jx,jy,jz-1))/dzz(jx,jy,jz)
          ELSE
            COORDINATE = 'Z'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-3) = -2.0d0*RoAveLeft*harz(jx,jy,jz-1)/(dzz(jx,jy,jz)*(dzz(jx,jy,jz)+dzz(jx,jy,jz-1))) 
            BodyForceZ = -COSD(z_angle)*grav*SignGravity*(-RoAveLeft*RoAveLeft*harz(jx,jy,jz-1))/dzz(jx,jy,jz)
            coef(3) = 0.0d0
          END IF
 
            tdepend = ScaleFactor*visc*ct*por(jx,jy,jz)/(dt*secyr)
            coef(0) = tdepend -( coef(-3) + coef(-2) + coef(-1) + coef(1) + coef(2) + coef(3) )

            CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
!!            CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j-nx*ny,coef(-3),INSERT_VALUES,ierr)
!!            CALL MatSetValues(amatP,1,j,1,j+nx*ny,coef(3),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)

            tdepend = visc*ct*por(jx,jy,jz)*pres(jx,jy,jz)/ (dt*secyr)
            pumpterm = visc*ro(jx,jy,jz)*qg(jx,jy,jz)/(secyr*dxx(jx)*dyy(jy)*dzz(jx,jy,jz))

          IF (activecellPressure(nx+1,jy,jz) == 0) THEN  
            AddPressureX =  2.0d0*ro(jx,jy,jz)*harx(jx,jy,jz)*pres(nx+1,jy,jz)/(dxx(jx)*dxx(jx))   
          ELSE
            AddPressureX = 0.0d0
          END IF 
          IF (activecellPressure(jx,jy,nz+1) == 0) THEN  
            AddPressureZ =  2.0d0*ro(jx,jy,jz)*harz(jx,jy,jz)*pres(jx,jy,nz+1)/(dzz(jx,jy,jz)*dzz(jx,jy,jz))   
          ELSE
            AddPressureZ = 0.0d0
          END IF 

            BvecCrunchP(j) = (pumpterm + tdepend + BodyForceX + BodyForceY + BodyForceZ + AddPressureX + AddPressureZ) 
            XvecCrunchP(j) = pres(jx,jy,jz)

          END IF
        END DO

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
            BvecCrunchP(j) = pres(jx,jy,jz)*big 
            XvecCrunchP(j) = pres(jx,jy,jz)
          ELSE

            COORDINATE = 'X'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-1) = -2.0d0*RoAveLeft*harx(jx-1,jy,jz)/(dxx(jx)*(dxx(jx)+dxx(jx-1))) 
            coef(1) = -2.0d0*RoAveRight*harx(jx,jy,jz)/(dxx(jx)*(dxx(jx+1)+dxx(jx))) 
            BodyForceX = -COSD(x_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harx(jx,jy,jz) - RoAveLeft*RoAveLeft*harx(jx-1,jy,jz))/dxx(jx)

          IF (activecellPressure(jx,0,jz) == 0) THEN
            COORDINATE = 'Y'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(2)  = -2.0d0*RoAveRight*hary(jx,jy,jz)/(dyy(jy) *(dyy(jy+1)+dyy(jy) )) 
            coef(-2) = -2.0d0*ro(jx,jy,jz)*hary(jx,jy,jz)/(dyy(jy)*dyy(jy))
            BodyForceY = -COSD(y_angle)*grav*SignGravity*(RoAveRight*RoAveRight*hary(jx,jy,jz) - ro(jx,jy,jz)*ro(jx,jy,jz)*hary(jx,jy,jz))/dyy(jy) 
          ELSE
            COORDINATE = 'Y'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(2)  = -2.0d0*RoAveRight*hary(jx,jy,jz)/(dyy(jy) *(dyy(jy+1)+dyy(jy) )) 
            BodyForceY = -COSD(y_angle)*grav*SignGravity*(RoAveRight*RoAveRight*hary(jx,jy,jz))/dyy(jy) 
            coef(-2) = 0.0d0
          END IF

          IF (activecellPressure(jx,jy,0) == 0) THEN
            COORDINATE = 'Z'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-3) = -2.0d0*ro(jx,jy,jz)*harz(jx,jy,jz)/(dzz(jx,jy,jz)*dzz(jx,jy,jz))
            coef(3)  = -2.0d0*RoAveRight*harz(jx,jy,jz)/(dzz(jx,jy,jz) *(dzz(jx,jy,jz+1)+dzz(jx,jy,jz) )) 
            BodyForceZ = -COSD(z_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harz(jx,jy,jz) - ro(jx,jy,jz)*ro(jx,jy,jz)*harz(jx,jy,jz))/dzz(jx,jy,jz) 
          ELSE
            COORDINATE = 'Z'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(3)  = -2.0d0*RoAveRight*harz(jx,jy,jz)/(dzz(jx,jy,jz) *(dzz(jx,jy,jz+1)+dzz(jx,jy,jz) )) 
            BodyForceZ = -COSD(z_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harz(jx,jy,jz))/dzz(jx,jy,jz) 
            coef(-3) = 0.0d0
          END IF
 
            tdepend = ScaleFactor*visc*ct*por(jx,jy,jz)/(dt*secyr)
            coef(0) = tdepend -( coef(-3) + coef(-2) + coef(-1) + coef(1) + coef(2) + coef(3) )

            CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
!!            CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
!!            CALL MatSetValues(amatP,1,j,1,j-nx*ny,coef(-3),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j+nx*ny,coef(3),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)

            tdepend = visc*ct*por(jx,jy,jz)*pres(jx,jy,jz)/ (dt*secyr)
            pumpterm = visc*ro(jx,jy,jz)*qg(jx,jy,jz)/(secyr*dxx(jx)*dyy(jy)*dzz(jx,jy,jz))

          IF (activecellPressure(jx,0,jz) == 0) THEN  
            AddPressureY =  2.0d0*ro(jx,jy,jz)*hary(jx,jy-1,jz)*pres(jx,0,jz)/(dyy(jy)*dyy(jy))   
          ELSE
            AddPressureY = 0.0d0
          END IF 
          IF (activecellPressure(jx,jy,0) == 0) THEN  
            AddPressureZ =  2.0d0*ro(jx,jy,jz)*harz(jx,jy,jz-1)*pres(jx,jy,0)/(dzz(jx,jy,jz)*dzz(jx,jy,jz))   
          ELSE
            AddPressureZ = 0.0d0
          END IF 

            BvecCrunchP(j) = (pumpterm + tdepend + BodyForceX + BodyForceY + BodyForceZ + AddPressureY + AddPressureZ) 
            XvecCrunchP(j) = pres(jx,jy,jz)

          END IF
        END DO

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
            BvecCrunchP(j) = pres(jx,jy,jz)*big 
            XvecCrunchP(j) = pres(jx,jy,jz)
          ELSE

            COORDINATE = 'X'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-1) = -2.0d0*RoAveLeft*harx(jx-1,jy,jz)/(dxx(jx)*(dxx(jx)+dxx(jx-1))) 
            coef(1) = -2.0d0*RoAveRight*harx(jx,jy,jz)/(dxx(jx)*(dxx(jx+1)+dxx(jx))) 
            BodyForceX = -COSD(x_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harx(jx,jy,jz) - RoAveLeft*RoAveLeft*harx(jx-1,jy,jz))/dxx(jx)

          IF (activecellPressure(jx,ny+1,jz) == 0) THEN
            COORDINATE = 'Y'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-2) = -2.0d0*RoAveLeft*hary(jx,jy-1,jz)/(dyy(jy)*(dyy(jy)+dyy(jy-1))) 
            coef(2)  = -2.0d0*ro(jx,jy,jz)*hary(jx,jy,jz)/(dyy(jy)*dyy(jy)) 
            BodyForceY = -COSD(y_angle)*grav*SignGravity*(ro(jx,jy,jz)*ro(jx,jy,jz)*hary(jx,jy,jz) - RoAveLeft*RoAveLeft*hary(jx,jy-1,jz))/dyy(jy) 
          ELSE
            COORDINATE = 'Y'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-2) = -2.0d0*RoAveLeft*hary(jx,jy-1,jz)/(dyy(jy)*(dyy(jy)+dyy(jy-1))) 
            BodyForceY = -COSD(y_angle)*grav*SignGravity*(-RoAveLeft*RoAveLeft*hary(jx,jy-1,jz))/dyy(jy)    
            coef(2) = 0.0d0     
          END IF

          IF (activecellPressure(jx,jy,0) == 0) THEN
            COORDINATE = 'Z'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-3) = -2.0d0*ro(jx,jy,jz)*harz(jx,jy,jz)/(dzz(jx,jy,jz)*dzz(jx,jy,jz))
            coef(3)  = -2.0d0*RoAveRight*harz(jx,jy,jz)/(dzz(jx,jy,jz) *(dzz(jx,jy,jz+1)+dzz(jx,jy,jz) )) 
            BodyForceZ = -COSD(z_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harz(jx,jy,jz) - ro(jx,jy,jz)*ro(jx,jy,jz)*harz(jx,jy,jz))/dzz(jx,jy,jz) 
          ELSE
            COORDINATE = 'Z'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(3)  = -2.0d0*RoAveRight*harz(jx,jy,jz)/(dzz(jx,jy,jz) *(dzz(jx,jy,jz+1)+dzz(jx,jy,jz) )) 
            BodyForceZ = -COSD(z_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harz(jx,jy,jz))/dzz(jx,jy,jz) 
            coef(-3) = 0.0d0
          END IF
 
            tdepend = ScaleFactor*visc*ct*por(jx,jy,jz)/(dt*secyr)
            coef(0) = tdepend -( coef(-3) + coef(-2) + coef(-1) + coef(1) + coef(2) + coef(3) )

            CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
!!            CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
!!            CALL MatSetValues(amatP,1,j,1,j-nx*ny,coef(-3),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j+nx*ny,coef(3),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)

            tdepend = visc*ct*por(jx,jy,jz)*pres(jx,jy,jz)/ (dt*secyr)
            pumpterm = visc*ro(jx,jy,jz)*qg(jx,jy,jz)/(secyr*dxx(jx)*dyy(jy)*dzz(jx,jy,jz))

          IF (activecellPressure(jx,ny+1,jz) == 0) THEN  
            AddPressureY =  2.0d0*ro(jx,jy,jz)*hary(jx,jy,jz)*pres(jx,ny+1,jz)/(dyy(jy)*dyy(jy))   
          ELSE
            AddPressureY = 0.0d0
          END IF 
          IF (activecellPressure(jx,jy,0) == 0) THEN  
            AddPressureZ =  2.0d0*ro(jx,jy,jz)*harz(jx,jy,jz-1)*pres(jx,jy,0)/(dzz(jx,jy,jz)*dzz(jx,jy,jz))   
          ELSE
            AddPressureZ = 0.0d0
          END IF 

            BvecCrunchP(j) = (pumpterm + tdepend + BodyForceX + BodyForceY + BodyForceZ + AddPressureY + AddPressureZ) 
            XvecCrunchP(j) = pres(jx,jy,jz)

          END IF
        END DO

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
            BvecCrunchP(j) = pres(jx,jy,jz)*big 
            XvecCrunchP(j) = pres(jx,jy,jz)
          ELSE

            COORDINATE = 'X'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-1) = -2.0d0*RoAveLeft*harx(jx-1,jy,jz)/(dxx(jx)*(dxx(jx)+dxx(jx-1))) 
            coef(1) = -2.0d0*RoAveRight*harx(jx,jy,jz)/(dxx(jx)*(dxx(jx+1)+dxx(jx))) 
            BodyForceX = -COSD(x_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harx(jx,jy,jz) - RoAveLeft*RoAveLeft*harx(jx-1,jy,jz))/dxx(jx)

          IF (activecellPressure(jx,0,jz) == 0) THEN
            COORDINATE = 'Y'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(2)  = -2.0d0*RoAveRight*hary(jx,jy,jz)/(dyy(jy) *(dyy(jy+1)+dyy(jy) )) 
            coef(-2) = -2.0d0*ro(jx,jy,jz)*hary(jx,jy,jz)/(dyy(jy)*dyy(jy))
            BodyForceY = -COSD(y_angle)*grav*SignGravity*(RoAveRight*RoAveRight*hary(jx,jy,jz) - ro(jx,jy,jz)*ro(jx,jy,jz)*hary(jx,jy,jz))/dyy(jy) 
          ELSE
            COORDINATE = 'Y'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(2)  = -2.0d0*RoAveRight*hary(jx,jy,jz)/(dyy(jy) *(dyy(jy+1)+dyy(jy) )) 
            BodyForceY = -COSD(y_angle)*grav*SignGravity*(RoAveRight*RoAveRight*hary(jx,jy,jz))/dyy(jy) 
            coef(-2) = 0.0d0
          END IF

          IF (activecellPressure(jx,jy,nz+1) == 0) THEN
            COORDINATE = 'Z'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-3) = -2.0d0*RoAveLeft*harz(jx,jy,jz-1)/(dzz(jx,jy,jz)*(dzz(jx,jy,jz)+dzz(jx,jy,jz-1))) 
            coef(3) = -2.0d0*ro(jx,jy,jz)*harz(jx,jy,jz)/(dzz(jx,jy,jz)*dzz(jx,jy,jz)) 
            BodyForceZ = -COSD(z_angle)*grav*SignGravity*(ro(jx,jy,jz)*ro(jx,jy,jz)*harz(jx,jy,jz) - RoAveLeft*RoAveLeft*harz(jx,jy,jz-1))/dzz(jx,jy,jz)
          ELSE
            COORDINATE = 'Z'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-3) = -2.0d0*RoAveLeft*harz(jx,jy,jz-1)/(dzz(jx,jy,jz)*(dzz(jx,jy,jz)+dzz(jx,jy,jz-1))) 
            BodyForceZ = -COSD(z_angle)*grav*SignGravity*(-RoAveLeft*RoAveLeft*harz(jx,jy,jz-1))/dzz(jx,jy,jz)
            coef(3) = 0.0d0
          END IF
 
            tdepend = ScaleFactor*visc*ct*por(jx,jy,jz)/(dt*secyr)
            coef(0) = tdepend -( coef(-3) + coef(-2) + coef(-1) + coef(1) + coef(2) + coef(3) )

            CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
!!            CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j-nx*ny,coef(-3),INSERT_VALUES,ierr)
!!            CALL MatSetValues(amatP,1,j,1,j+nx*ny,coef(3),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)

            tdepend = visc*ct*por(jx,jy,jz)*pres(jx,jy,jz)/ (dt*secyr)
            pumpterm = visc*ro(jx,jy,jz)*qg(jx,jy,jz)/(secyr*dxx(jx)*dyy(jy)*dzz(jx,jy,jz))

          IF (activecellPressure(jx,0,jz) == 0) THEN  
            AddPressureY =  2.0d0*ro(jx,jy,jz)*hary(jx,jy-1,jz)*pres(jx,0,jz)/(dyy(jy)*dyy(jy))   
          ELSE
            AddPressureY = 0.0d0
          END IF 
          IF (activecellPressure(jx,jy,nz+1) == 0) THEN  
            AddPressureZ =  2.0d0*ro(jx,jy,jz)*harz(jx,jy,jz)*pres(jx,jy,nz+1)/(dzz(jx,jy,jz)*dzz(jx,jy,jz))   
          ELSE
            AddPressureZ = 0.0d0
          END IF 

            BvecCrunchP(j) = (pumpterm + tdepend + BodyForceX + BodyForceY + BodyForceZ + AddPressureY + AddPressureZ) 
            XvecCrunchP(j) = pres(jx,jy,jz)

          END IF
        END DO

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
            BvecCrunchP(j) = pres(jx,jy,jz)*big 
            XvecCrunchP(j) = pres(jx,jy,jz)
          ELSE

            COORDINATE = 'X'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-1) = -2.0d0*RoAveLeft*harx(jx-1,jy,jz)/(dxx(jx)*(dxx(jx)+dxx(jx-1))) 
            coef(1) = -2.0d0*RoAveRight*harx(jx,jy,jz)/(dxx(jx)*(dxx(jx+1)+dxx(jx))) 
            BodyForceX = -COSD(x_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harx(jx,jy,jz) - RoAveLeft*RoAveLeft*harx(jx-1,jy,jz))/dxx(jx)

          IF (activecellPressure(jx,ny+1,jz) == 0) THEN
            COORDINATE = 'Y'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-2) = -2.0d0*RoAveLeft*hary(jx,jy-1,jz)/(dyy(jy)*(dyy(jy)+dyy(jy-1))) 
            coef(2)  = -2.0d0*ro(jx,jy,jz)*hary(jx,jy,jz)/(dyy(jy)*dyy(jy)) 
            BodyForceY = -COSD(y_angle)*grav*SignGravity*(ro(jx,jy,jz)*ro(jx,jy,jz)*hary(jx,jy,jz) - RoAveLeft*RoAveLeft*hary(jx,jy-1,jz))/dyy(jy) 
          ELSE
            COORDINATE = 'Y'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-2) = -2.0d0*RoAveLeft*hary(jx,jy-1,jz)/(dyy(jy)*(dyy(jy)+dyy(jy-1))) 
            BodyForceY = -COSD(y_angle)*grav*SignGravity*(-RoAveLeft*RoAveLeft*hary(jx,jy-1,jz))/dyy(jy)    
            coef(2) = 0.0d0     
          END IF

          IF (activecellPressure(jx,jy,nz+1) == 0) THEN
            COORDINATE = 'Z'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-3) = -2.0d0*RoAveLeft*harz(jx,jy,jz-1)/(dzz(jx,jy,jz)*(dzz(jx,jy,jz)+dzz(jx,jy,jz-1))) 
            coef(3) = -2.0d0*ro(jx,jy,jz)*harz(jx,jy,jz)/(dzz(jx,jy,jz)*dzz(jx,jy,jz)) 
            BodyForceZ = -COSD(z_angle)*grav*SignGravity*(ro(jx,jy,jz)*ro(jx,jy,jz)*harz(jx,jy,jz) - RoAveLeft*RoAveLeft*harz(jx,jy,jz-1))/dzz(jx,jy,jz)
          ELSE
            COORDINATE = 'Z'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-3) = -2.0d0*RoAveLeft*harz(jx,jy,jz-1)/(dzz(jx,jy,jz)*(dzz(jx,jy,jz)+dzz(jx,jy,jz-1))) 
            BodyForceZ = -COSD(z_angle)*grav*SignGravity*(-RoAveLeft*RoAveLeft*harz(jx,jy,jz-1))/dzz(jx,jy,jz)
            coef(3) = 0.0d0
          END IF
 
            tdepend = ScaleFactor*visc*ct*por(jx,jy,jz)/(dt*secyr)
            coef(0) = tdepend -( coef(-3) + coef(-2) + coef(-1) + coef(1) + coef(2) + coef(3) )

            CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
!!            CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j-nx*ny,coef(-3),INSERT_VALUES,ierr)
!!            CALL MatSetValues(amatP,1,j,1,j+nx*ny,coef(3),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)

            tdepend = visc*ct*por(jx,jy,jz)*pres(jx,jy,jz)/ (dt*secyr)
            pumpterm = visc*ro(jx,jy,jz)*qg(jx,jy,jz)/(secyr*dxx(jx)*dyy(jy)*dzz(jx,jy,jz))

          IF (activecellPressure(jx,ny+1,jz) == 0) THEN  
            AddPressureY =  2.0d0*ro(jx,jy,jz)*hary(jx,jy,jz)*pres(jx,ny+1,jz)/(dyy(jy)*dyy(jy))   
          ELSE
            AddPressureY = 0.0d0
          END IF 
          IF (activecellPressure(jx,jy,nz+1) == 0) THEN  
            AddPressureZ =  2.0d0*ro(jx,jy,jz)*harz(jx,jy,jz)*pres(jx,jy,nz+1)/(dzz(jx,jy,jz)*dzz(jx,jy,jz))   
          ELSE
            AddPressureZ = 0.0d0
          END IF 

            BvecCrunchP(j) = (pumpterm + tdepend + BodyForceX + BodyForceY + BodyForceZ + AddPressureY + AddPressureZ) 
            XvecCrunchP(j) = pres(jx,jy,jz)

          END IF
        END DO



!!  **************************************************************************************************************************************************
!!  Now the corners

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
            BvecCrunchP(j) = pres(jx,jy,jz)*big 
            XvecCrunchP(j) = pres(jx,jy,jz)
          ELSE

          IF (activecellPressure(0,jy,jz) == 0) THEN
            COORDINATE = 'X'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-1) = -2.0d0*ro(jx,jy,jz)*harx(jx,jy,jz)/(dxx(jx)*dxx(jx))
            coef(1) = -2.0d0*RoAveRight*harx(jx,jy,jz)/(dxx(jx)*(dxx(jx+1)+dxx(jx))) 
            BodyForceX = -COSD(x_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harx(jx,jy,jz) - ro(jx,jy,jz)*ro(jx,jy,jz)*harx(jx,jy,jz))/dxx(jx)
          ELSE
            COORDINATE = 'X'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(1) = -2.0d0*RoAveRight*harx(jx,jy,jz)/(dxx(jx)*(dxx(jx+1)+dxx(jx))) 
            BodyForceX = -COSD(x_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harx(jx,jy,jz) )/dxx(jx)
            coef(-1) = 0.0d0
          END IF

          IF (activecellPressure(jx,0,jz) == 0) THEN
            COORDINATE = 'Y'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(2)  = -2.0d0*RoAveRight*hary(jx,jy,jz)/(dyy(jy) *(dyy(jy+1)+dyy(jy) )) 
            coef(-2) = -2.0d0*ro(jx,jy,jz)*hary(jx,jy,jz)/(dyy(jy)*dyy(jy))
            BodyForceY = -COSD(y_angle)*grav*SignGravity*(RoAveRight*RoAveRight*hary(jx,jy,jz) - ro(jx,jy,jz)*ro(jx,jy,jz)*hary(jx,jy,jz))/dyy(jy) 
          ELSE
            COORDINATE = 'Y'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(2)  = -2.0d0*RoAveRight*hary(jx,jy,jz)/(dyy(jy) *(dyy(jy+1)+dyy(jy) )) 
            BodyForceY = -COSD(y_angle)*grav*SignGravity*(RoAveRight*RoAveRight*hary(jx,jy,jz))/dyy(jy) 
            coef(-2) = 0.0d0
          END IF

          IF (activecellPressure(jx,jy,0) == 0) THEN
            COORDINATE = 'Z'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-3) = -2.0d0*ro(jx,jy,jz)*harz(jx,jy,jz)/(dzz(jx,jy,jz)*dzz(jx,jy,jz))
            coef(3)  = -2.0d0*RoAveRight*harz(jx,jy,jz)/(dzz(jx,jy,jz) *(dzz(jx,jy,jz+1)+dzz(jx,jy,jz) )) 
            BodyForceZ = -COSD(z_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harz(jx,jy,jz) - ro(jx,jy,jz)*ro(jx,jy,jz)*harz(jx,jy,jz))/dzz(jx,jy,jz) 
          ELSE
            COORDINATE = 'Z'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(3)  = -2.0d0*RoAveRight*harz(jx,jy,jz)/(dzz(jx,jy,jz) *(dzz(jx,jy,jz+1)+dzz(jx,jy,jz) )) 
            BodyForceZ = -COSD(z_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harz(jx,jy,jz))/dzz(jx,jy,jz) 
            coef(-3) = 0.0d0
          END IF
 
            tdepend = ScaleFactor*visc*ct*por(jx,jy,jz)/(dt*secyr)
            coef(0) = tdepend -( coef(-3) + coef(-2) + coef(-1) + coef(1) + coef(2) + coef(3) )

!!            CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
!!            CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
!!            CALL MatSetValues(amatP,1,j,1,j-nx*ny,coef(-3),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j+nx*ny,coef(3),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)

            tdepend = visc*ct*por(jx,jy,jz)*pres(jx,jy,jz)/ (dt*secyr)
            pumpterm = visc*ro(jx,jy,jz)*qg(jx,jy,jz)/(secyr*dxx(jx)*dyy(jy)*dzz(jx,jy,jz))

          IF (activecellPressure(0,jy,jz) == 0) THEN  
            AddPressureX =  2.0d0*ro(jx,jy,jz)*harx(jx-1,jy,jz)*pres(0,jy,jz)/(dxx(jx)*dxx(jx))   
          ELSE
            AddPressureX = 0.0d0
          END IF 
          IF (activecellPressure(jx,0,jz) == 0) THEN  
            AddPressureY =  2.0d0*ro(jx,jy,jz)*hary(jx,jy-1,jz)*pres(jx,0,jz)/(dyy(jy)*dyy(jy))   
          ELSE
            AddPressureY = 0.0d0
          END IF 
          IF (activecellPressure(jx,jy,0) == 0) THEN  
            AddPressureZ =  2.0d0*ro(jx,jy,jz)*harz(jx,jy,jz-1)*pres(jx,jy,0)/(dzz(jx,jy,jz)*dzz(jx,jy,jz))   
          ELSE
            AddPressureZ = 0.0d0
          END IF 

            BvecCrunchP(j) = (pumpterm + tdepend + BodyForceX + BodyForceY + BodyForceZ + AddPressureX + AddPressureY + AddPressureZ) 
            XvecCrunchP(j) = pres(jx,jy,jz)

          END IF

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
            BvecCrunchP(j) = pres(jx,jy,jz)*big 
            XvecCrunchP(j) = pres(jx,jy,jz)
          ELSE

          IF (activecellPressure(nx+1,jy,jz) == 0) THEN
            COORDINATE = 'X'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-1) = -2.0d0*RoAveLeft*harx(jx-1,jy,jz)/(dxx(jx)*(dxx(jx)+dxx(jx-1))) 
            coef(1) = -2.0d0*ro(jx,jy,jz)*harx(jx,jy,jz)/(dxx(jx)*dxx(jx)) 
            BodyForceX = -COSD(x_angle)*grav*SignGravity*(ro(jx,jy,jz)*ro(jx,jy,jz)*harx(jx,jy,jz) - RoAveLeft*RoAveLeft*harx(jx-1,jy,jz))/dxx(jx)
          ELSE
            COORDINATE = 'X'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-1) = -2.0d0*RoAveLeft*harx(jx-1,jy,jz)/(dxx(jx)*(dxx(jx)+dxx(jx-1))) 
            BodyForceX = -COSD(x_angle)*grav*SignGravity*(-RoAveLeft*RoAveLeft*harx(jx-1,jy,jz))/dxx(jx)
            coef(1) = 0.0d0
          END IF

          IF (activecellPressure(jx,0,jz) == 0) THEN
            COORDINATE = 'Y'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(2)  = -2.0d0*RoAveRight*hary(jx,jy,jz)/(dyy(jy) *(dyy(jy+1)+dyy(jy) )) 
            coef(-2) = -2.0d0*ro(jx,jy,jz)*hary(jx,jy,jz)/(dyy(jy)*dyy(jy))
            BodyForceY = -COSD(y_angle)*grav*SignGravity*(RoAveRight*RoAveRight*hary(jx,jy,jz) - ro(jx,jy,jz)*ro(jx,jy,jz)*hary(jx,jy,jz))/dyy(jy) 
          ELSE
            COORDINATE = 'Y'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(2)  = -2.0d0*RoAveRight*hary(jx,jy,jz)/(dyy(jy) *(dyy(jy+1)+dyy(jy) )) 
            BodyForceY = -COSD(y_angle)*grav*SignGravity*(RoAveRight*RoAveRight*hary(jx,jy,jz))/dyy(jy) 
            coef(-2) = 0.0d0
          END IF 

          IF (activecellPressure(jx,jy,0) == 0) THEN
            COORDINATE = 'Z'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-3) = -2.0d0*ro(jx,jy,jz)*harz(jx,jy,jz)/(dzz(jx,jy,jz)*dzz(jx,jy,jz))
            coef(3)  = -2.0d0*RoAveRight*harz(jx,jy,jz)/(dzz(jx,jy,jz) *(dzz(jx,jy,jz+1)+dzz(jx,jy,jz) )) 
            BodyForceZ = -COSD(z_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harz(jx,jy,jz) - ro(jx,jy,jz)*ro(jx,jy,jz)*harz(jx,jy,jz))/dzz(jx,jy,jz) 
          ELSE
            COORDINATE = 'Z'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(3)  = -2.0d0*RoAveRight*harz(jx,jy,jz)/(dzz(jx,jy,jz) *(dzz(jx,jy,jz+1)+dzz(jx,jy,jz) )) 
            BodyForceZ = -COSD(z_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harz(jx,jy,jz))/dzz(jx,jy,jz) 
            coef(-3) = 0.0d0
          END IF
 
            tdepend = ScaleFactor*visc*ct*por(jx,jy,jz)/(dt*secyr)
            coef(0) = tdepend -( coef(-3) + coef(-2) + coef(-1) + coef(1) + coef(2) + coef(3) )

            CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
!!            CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
!!            CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
!!            CALL MatSetValues(amatP,1,j,1,j-nx*ny,coef(-3),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j+nx*ny,coef(3),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)

            tdepend = visc*ct*por(jx,jy,jz)*pres(jx,jy,jz)/ (dt*secyr)
            pumpterm = visc*ro(jx,jy,jz)*qg(jx,jy,jz)/(secyr*dxx(jx)*dyy(jy)*dzz(jx,jy,jz))

          IF (activecellPressure(nx+1,jy,jz) == 0) THEN  
            AddPressureX =  2.0d0*ro(jx,jy,jz)*harx(jx,jy,jz)*pres(nx+1,jy,jz)/(dxx(jx)*dxx(jx))   
          ELSE
            AddPressureX = 0.0d0
          END IF 
          IF (activecellPressure(jx,0,jz) == 0) THEN  
            AddPressureY =  2.0d0*ro(jx,jy,jz)*hary(jx,jy-1,jz)*pres(jx,0,jz)/(dyy(jy)*dyy(jy))   
          ELSE
            AddPressureY = 0.0d0
          END IF 
          IF (activecellPressure(jx,jy,0) == 0) THEN  
            AddPressureZ =  2.0d0*ro(jx,jy,jz)*harz(jx,jy,jz-1)*pres(jx,jy,0)/(dzz(jx,jy,jz)*dzz(jx,jy,jz))   
          ELSE
            AddPressureZ = 0.0d0
          END IF 

            BvecCrunchP(j) = (pumpterm + tdepend + BodyForceX + BodyForceY + BodyForceZ + AddPressureX + AddPressureY + AddPressureZ) 
            XvecCrunchP(j) = pres(jx,jy,jz)

          END IF

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
            BvecCrunchP(j) = pres(jx,jy,jz)*big 
            XvecCrunchP(j) = pres(jx,jy,jz)
          ELSE

          IF (activecellPressure(0,jy,jz) == 0) THEN
            COORDINATE = 'X'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-1) = -2.0d0*ro(jx,jy,jz)*harx(jx,jy,jz)/(dxx(jx)*dxx(jx))
            coef(1) = -2.0d0*RoAveRight*harx(jx,jy,jz)/(dxx(jx)*(dxx(jx+1)+dxx(jx))) 
            BodyForceX = -COSD(x_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harx(jx,jy,jz) - ro(jx,jy,jz)*ro(jx,jy,jz)*harx(jx,jy,jz))/dxx(jx)
          ELSE
            COORDINATE = 'X'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(1) = -2.0d0*RoAveRight*harx(jx,jy,jz)/(dxx(jx)*(dxx(jx+1)+dxx(jx))) 
            BodyForceX = -COSD(x_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harx(jx,jy,jz) )/dxx(jx)
            coef(-1) = 0.0d0
          END IF

         IF (activecellPressure(jx,ny+1,jz) == 0) THEN
            COORDINATE = 'Y'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-2) = -2.0d0*RoAveLeft*hary(jx,jy-1,jz)/(dyy(jy)*(dyy(jy)+dyy(jy-1))) 
            coef(2)  = -2.0d0*ro(jx,jy,jz)*hary(jx,jy,jz)/(dyy(jy)*dyy(jy)) 
            BodyForceY = -COSD(y_angle)*grav*SignGravity*(ro(jx,jy,jz)*ro(jx,jy,jz)*hary(jx,jy,jz) - RoAveLeft*RoAveLeft*hary(jx,jy-1,jz))/dyy(jy) 
          ELSE
            COORDINATE = 'Y'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-2) = -2.0d0*RoAveLeft*hary(jx,jy-1,jz)/(dyy(jy)*(dyy(jy)+dyy(jy-1))) 
            BodyForceY = -COSD(y_angle)*grav*SignGravity*(-RoAveLeft*RoAveLeft*hary(jx,jy-1,jz))/dyy(jy)    
            coef(2) = 0.0d0     
          END IF

          IF (activecellPressure(jx,jy,0) == 0) THEN
            COORDINATE = 'Z'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-3) = -2.0d0*ro(jx,jy,jz)*harz(jx,jy,jz)/(dzz(jx,jy,jz)*dzz(jx,jy,jz))
            coef(3)  = -2.0d0*RoAveRight*harz(jx,jy,jz)/(dzz(jx,jy,jz) *(dzz(jx,jy,jz+1)+dzz(jx,jy,jz) )) 
            BodyForceZ = -COSD(z_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harz(jx,jy,jz) - ro(jx,jy,jz)*ro(jx,jy,jz)*harz(jx,jy,jz))/dzz(jx,jy,jz) 
          ELSE
            COORDINATE = 'Z'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(3)  = -2.0d0*RoAveRight*harz(jx,jy,jz)/(dzz(jx,jy,jz) *(dzz(jx,jy,jz+1)+dzz(jx,jy,jz) )) 
            BodyForceZ = -COSD(z_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harz(jx,jy,jz))/dzz(jx,jy,jz) 
            coef(-3) = 0.0d0
          END IF
 
            tdepend = ScaleFactor*visc*ct*por(jx,jy,jz)/(dt*secyr)
            coef(0) = tdepend -( coef(-3) + coef(-2) + coef(-1) + coef(1) + coef(2) + coef(3) )

!!            CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
!!            CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
!!            CALL MatSetValues(amatP,1,j,1,j-nx*ny,coef(-3),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j+nx*ny,coef(3),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)

            tdepend = visc*ct*por(jx,jy,jz)*pres(jx,jy,jz)/ (dt*secyr)
            pumpterm = visc*ro(jx,jy,jz)*qg(jx,jy,jz)/(secyr*dxx(jx)*dyy(jy)*dzz(jx,jy,jz))

          IF (activecellPressure(0,jy,jz) == 0) THEN  
            AddPressureX =  2.0d0*ro(jx,jy,jz)*harx(jx-1,jy,jz)*pres(0,jy,jz)/(dxx(jx)*dxx(jx))   
          ELSE
            AddPressureX = 0.0d0
          END IF 
          IF (activecellPressure(jx,ny+1,jz) == 0) THEN  
            AddPressureY =  2.0d0*ro(jx,jy,jz)*hary(jx,jy-1,jz)*pres(jx,ny+1,jz)/(dyy(jy)*dyy(jy))   
          ELSE
            AddPressureY = 0.0d0
          END IF 
          IF (activecellPressure(jx,jy,0) == 0) THEN  
            AddPressureZ =  2.0d0*ro(jx,jy,jz)*harz(jx,jy,jz-1)*pres(jx,jy,0)/(dzz(jx,jy,jz)*dzz(jx,jy,jz))   
          ELSE
            AddPressureZ = 0.0d0
          END IF 

            BvecCrunchP(j) = (pumpterm + tdepend + BodyForceX + BodyForceY + BodyForceZ + AddPressureX + AddPressureY + AddPressureZ) 
            XvecCrunchP(j) = pres(jx,jy,jz)

          END IF

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
            BvecCrunchP(j) = pres(jx,jy,jz)*big 
            XvecCrunchP(j) = pres(jx,jy,jz)
          ELSE

          IF (activecellPressure(nx+1,jy,jz) == 0) THEN
            COORDINATE = 'X'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-1) = -2.0d0*RoAveLeft*harx(jx-1,jy,jz)/(dxx(jx)*(dxx(jx)+dxx(jx-1))) 
            coef(1) = -2.0d0*ro(jx,jy,jz)*harx(jx,jy,jz)/(dxx(jx)*dxx(jx)) 
            BodyForceX = -COSD(x_angle)*grav*SignGravity*(ro(jx,jy,jz)*ro(jx,jy,jz)*harx(jx,jy,jz) - RoAveLeft*RoAveLeft*harx(jx-1,jy,jz))/dxx(jx)
          ELSE
            COORDINATE = 'X'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-1) = -2.0d0*RoAveLeft*harx(jx-1,jy,jz)/(dxx(jx)*(dxx(jx)+dxx(jx-1))) 
            BodyForceX = -COSD(x_angle)*grav*SignGravity*(-RoAveLeft*RoAveLeft*harx(jx-1,jy,jz))/dxx(jx)
            coef(1) = 0.0d0
          END IF

          IF (activecellPressure(jx,ny+1,jz) == 0) THEN
            COORDINATE = 'Y'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-2) = -2.0d0*RoAveLeft*hary(jx,jy-1,jz)/(dyy(jy)*(dyy(jy)+dyy(jy-1))) 
            coef(2)  = -2.0d0*ro(jx,jy,jz)*hary(jx,jy,jz)/(dyy(jy)*dyy(jy)) 
            BodyForceY = -COSD(y_angle)*grav*SignGravity*(ro(jx,jy,jz)*ro(jx,jy,jz)*hary(jx,jy,jz) - RoAveLeft*RoAveLeft*hary(jx,jy-1,jz))/dyy(jy) 
          ELSE
            COORDINATE = 'Y'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-2) = -2.0d0*RoAveLeft*hary(jx,jy-1,jz)/(dyy(jy)*(dyy(jy)+dyy(jy-1))) 
            BodyForceY = -COSD(y_angle)*grav*SignGravity*(-RoAveLeft*RoAveLeft*hary(jx,jy-1,jz))/dyy(jy)    
            coef(2) = 0.0d0     
          END IF

          IF (activecellPressure(jx,jy,0) == 0) THEN
            COORDINATE = 'Z'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-3) = -2.0d0*ro(jx,jy,jz)*harz(jx,jy,jz)/(dzz(jx,jy,jz)*dzz(jx,jy,jz))
            coef(3)  = -2.0d0*RoAveRight*harz(jx,jy,jz)/(dzz(jx,jy,jz) *(dzz(jx,jy,jz+1)+dzz(jx,jy,jz) )) 
            BodyForceZ = -COSD(z_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harz(jx,jy,jz) - ro(jx,jy,jz)*ro(jx,jy,jz)*harz(jx,jy,jz))/dzz(jx,jy,jz) 
          ELSE
            COORDINATE = 'Z'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(3)  = -2.0d0*RoAveRight*harz(jx,jy,jz)/(dzz(jx,jy,jz) *(dzz(jx,jy,jz+1)+dzz(jx,jy,jz) )) 
            BodyForceZ = -COSD(z_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harz(jx,jy,jz))/dzz(jx,jy,jz) 
            coef(-3) = 0.0d0
          END IF
 
            tdepend = ScaleFactor*visc*ct*por(jx,jy,jz)/(dt*secyr)
            coef(0) = tdepend -( coef(-3) + coef(-2) + coef(-1) + coef(1) + coef(2) + coef(3) )

            CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
!!            CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
!!            CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
!!            CALL MatSetValues(amatP,1,j,1,j-nx*ny,coef(-3),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j+nx*ny,coef(3),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)

            tdepend = visc*ct*por(jx,jy,jz)*pres(jx,jy,jz)/ (dt*secyr)
            pumpterm = visc*ro(jx,jy,jz)*qg(jx,jy,jz)/(secyr*dxx(jx)*dyy(jy)*dzz(jx,jy,jz))

          IF (activecellPressure(nx+1,jy,jz) == 0) THEN  
            AddPressureX =  2.0d0*ro(jx,jy,jz)*harx(jx,jy,jz)*pres(nx+1,jy,jz)/(dxx(jx)*dxx(jx))   
          ELSE
            AddPressureX = 0.0d0
          END IF 
          IF (activecellPressure(jx,ny+1,jz) == 0) THEN  
            AddPressureY =  2.0d0*ro(jx,jy,jz)*hary(jx,jy,jz)*pres(jx,ny+1,jz)/(dyy(jy)*dyy(jy))   
          ELSE
            AddPressureY = 0.0d0
          END IF 
          IF (activecellPressure(jx,jy,0) == 0) THEN  
            AddPressureZ =  2.0d0*ro(jx,jy,jz)*harz(jx,jy,jz-1)*pres(jx,jy,0)/(dzz(jx,jy,jz)*dzz(jx,jy,jz))   
          ELSE
            AddPressureZ = 0.0d0
          END IF 

            BvecCrunchP(j) = (pumpterm + tdepend + BodyForceX + BodyForceY + BodyForceZ + AddPressureX + AddPressureY + AddPressureZ) 
            XvecCrunchP(j) = pres(jx,jy,jz)

          END IF

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
            BvecCrunchP(j) = pres(jx,jy,jz)*big 
            XvecCrunchP(j) = pres(jx,jy,jz)
          ELSE

          IF (activecellPressure(0,jy,jz) == 0) THEN
            COORDINATE = 'X'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-1) = -2.0d0*ro(jx,jy,jz)*harx(jx,jy,jz)/(dxx(jx)*dxx(jx))
            coef(1) = -2.0d0*RoAveRight*harx(jx,jy,jz)/(dxx(jx)*(dxx(jx+1)+dxx(jx))) 
            BodyForceX = -COSD(x_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harx(jx,jy,jz) - ro(jx,jy,jz)*ro(jx,jy,jz)*harx(jx,jy,jz))/dxx(jx)
          ELSE
            COORDINATE = 'X'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(1) = -2.0d0*RoAveRight*harx(jx,jy,jz)/(dxx(jx)*(dxx(jx+1)+dxx(jx))) 
            BodyForceX = -COSD(x_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harx(jx,jy,jz) )/dxx(jx)
            coef(-1) = 0.0d0
          END IF

          IF (activecellPressure(jx,0,jz) == 0) THEN
            COORDINATE = 'Y'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(2)  = -2.0d0*RoAveRight*hary(jx,jy,jz)/(dyy(jy) *(dyy(jy+1)+dyy(jy) )) 
            coef(-2) = -2.0d0*ro(jx,jy,jz)*hary(jx,jy,jz)/(dyy(jy)*dyy(jy))
            BodyForceY = -COSD(y_angle)*grav*SignGravity*(RoAveRight*RoAveRight*hary(jx,jy,jz) - ro(jx,jy,jz)*ro(jx,jy,jz)*hary(jx,jy,jz))/dyy(jy) 
          ELSE
            COORDINATE = 'Y'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(2)  = -2.0d0*RoAveRight*hary(jx,jy,jz)/(dyy(jy) *(dyy(jy+1)+dyy(jy) )) 
            BodyForceY = -COSD(y_angle)*grav*SignGravity*(RoAveRight*RoAveRight*hary(jx,jy,jz))/dyy(jy) 
            coef(-2) = 0.0d0
          END IF

          IF (activecellPressure(jx,jy,nz+1) == 0) THEN
            COORDINATE = 'Z'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-3) = -2.0d0*RoAveLeft*harz(jx,jy,jz-1)/(dzz(jx,jy,jz)*(dzz(jx,jy,jz)+dzz(jx,jy,jz-1))) 
            coef(3) = -2.0d0*ro(jx,jy,jz)*harz(jx,jy,jz)/(dzz(jx,jy,jz)*dzz(jx,jy,jz)) 
            BodyForceZ = -COSD(z_angle)*grav*SignGravity*(ro(jx,jy,jz)*ro(jx,jy,jz)*harz(jx,jy,jz) - RoAveLeft*RoAveLeft*harz(jx,jy,jz-1))/dzz(jx,jy,jz)
          ELSE
            COORDINATE = 'Z'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-3) = -2.0d0*RoAveLeft*harz(jx,jy,jz-1)/(dzz(jx,jy,jz)*(dzz(jx,jy,jz)+dzz(jx,jy,jz-1))) 
            BodyForceZ = -COSD(z_angle)*grav*SignGravity*(-RoAveLeft*RoAveLeft*harz(jx,jy,jz-1))/dzz(jx,jy,jz)
            coef(3) = 0.0d0
          END IF
 
            tdepend = ScaleFactor*visc*ct*por(jx,jy,jz)/(dt*secyr)
            coef(0) = tdepend -( coef(-3) + coef(-2) + coef(-1) + coef(1) + coef(2) + coef(3) )

!!            CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
!!            CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j-nx*ny,coef(-3),INSERT_VALUES,ierr)
!!            CALL MatSetValues(amatP,1,j,1,j+nx*ny,coef(3),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)

            tdepend = visc*ct*por(jx,jy,jz)*pres(jx,jy,jz)/ (dt*secyr)
            pumpterm = visc*ro(jx,jy,jz)*qg(jx,jy,jz)/(secyr*dxx(jx)*dyy(jy)*dzz(jx,jy,jz))

          IF (activecellPressure(0,jy,jz) == 0) THEN  
            AddPressureX =  2.0d0*ro(jx,jy,jz)*harx(jx-1,jy,jz)*pres(0,jy,jz)/(dxx(jx)*dxx(jx))   
          ELSE
            AddPressureX = 0.0d0
          END IF 
          IF (activecellPressure(jx,0,jz) == 0) THEN  
            AddPressureY =  2.0d0*ro(jx,jy,jz)*hary(jx,jy-1,jz)*pres(jx,0,jz)/(dyy(jy)*dyy(jy))   
          ELSE
            AddPressureY = 0.0d0
          END IF 
          IF (activecellPressure(jx,jy,nz+1) == 0) THEN  
            AddPressureZ =  2.0d0*ro(jx,jy,jz)*harz(jx,jy,jz)*pres(jx,jy,nz+1)/(dzz(jx,jy,jz)*dzz(jx,jy,jz))   
          ELSE
            AddPressureZ = 0.0d0
          END IF 

            BvecCrunchP(j) = (pumpterm + tdepend + BodyForceX + BodyForceY + BodyForceZ + AddPressureX + AddPressureY + AddPressureZ) 
            XvecCrunchP(j) = pres(jx,jy,jz)

          END IF

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
            BvecCrunchP(j) = pres(jx,jy,jz)*big 
            XvecCrunchP(j) = pres(jx,jy,jz)
          ELSE

          IF (activecellPressure(nx+1,jy,jz) == 0) THEN
            COORDINATE = 'X'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-1) = -2.0d0*RoAveLeft*harx(jx-1,jy,jz)/(dxx(jx)*(dxx(jx)+dxx(jx-1))) 
            coef(1) = -2.0d0*ro(jx,jy,jz)*harx(jx,jy,jz)/(dxx(jx)*dxx(jx)) 
            BodyForceX = -COSD(x_angle)*grav*SignGravity*(ro(jx,jy,jz)*ro(jx,jy,jz)*harx(jx,jy,jz) - RoAveLeft*RoAveLeft*harx(jx-1,jy,jz))/dxx(jx)
          ELSE
            COORDINATE = 'X'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-1) = -2.0d0*RoAveLeft*harx(jx-1,jy,jz)/(dxx(jx)*(dxx(jx)+dxx(jx-1))) 
            BodyForceX = -COSD(x_angle)*grav*SignGravity*(-RoAveLeft*RoAveLeft*harx(jx-1,jy,jz))/dxx(jx)
            coef(1) = 0.0d0
          END IF

          IF (activecellPressure(jx,0,jz) == 0) THEN
            COORDINATE = 'Y'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(2)  = -2.0d0*RoAveRight*hary(jx,jy,jz)/(dyy(jy) *(dyy(jy+1)+dyy(jy) )) 
            coef(-2) = -2.0d0*ro(jx,jy,jz)*hary(jx,jy,jz)/(dyy(jy)*dyy(jy))
            BodyForceY = -COSD(y_angle)*grav*SignGravity*(RoAveRight*RoAveRight*hary(jx,jy,jz) - ro(jx,jy,jz)*ro(jx,jy,jz)*hary(jx,jy,jz))/dyy(jy) 
          ELSE
            COORDINATE = 'Y'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(2)  = -2.0d0*RoAveRight*hary(jx,jy,jz)/(dyy(jy) *(dyy(jy+1)+dyy(jy) )) 
            BodyForceY = -COSD(y_angle)*grav*SignGravity*(RoAveRight*RoAveRight*hary(jx,jy,jz))/dyy(jy) 
            coef(-2) = 0.0d0
          END IF

          IF (activecellPressure(jx,jy,nz+1) == 0) THEN
            COORDINATE = 'Z'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-3) = -2.0d0*RoAveLeft*harz(jx,jy,jz-1)/(dzz(jx,jy,jz)*(dzz(jx,jy,jz)+dzz(jx,jy,jz-1))) 
            coef(3) = -2.0d0*ro(jx,jy,jz)*harz(jx,jy,jz)/(dzz(jx,jy,jz)*dzz(jx,jy,jz)) 
            BodyForceZ = -COSD(z_angle)*grav*SignGravity*(ro(jx,jy,jz)*ro(jx,jy,jz)*harz(jx,jy,jz) - RoAveLeft*RoAveLeft*harz(jx,jy,jz-1))/dzz(jx,jy,jz)
          ELSE
            COORDINATE = 'Z'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-3) = -2.0d0*RoAveLeft*harz(jx,jy,jz-1)/(dzz(jx,jy,jz)*(dzz(jx,jy,jz)+dzz(jx,jy,jz-1))) 
            BodyForceZ = -COSD(z_angle)*grav*SignGravity*(-RoAveLeft*RoAveLeft*harz(jx,jy,jz-1))/dzz(jx,jy,jz)
            coef(3) = 0.0d0
          END IF
 
            tdepend = ScaleFactor*visc*ct*por(jx,jy,jz)/(dt*secyr)
            coef(0) = tdepend -( coef(-3) + coef(-2) + coef(-1) + coef(1) + coef(2) + coef(3) )

            CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
!!            CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
!!            CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j-nx*ny,coef(-3),INSERT_VALUES,ierr)
!!            CALL MatSetValues(amatP,1,j,1,j+nx*ny,coef(3),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)

            tdepend = visc*ct*por(jx,jy,jz)*pres(jx,jy,jz)/ (dt*secyr)
            pumpterm = visc*ro(jx,jy,jz)*qg(jx,jy,jz)/(secyr*dxx(jx)*dyy(jy)*dzz(jx,jy,jz))

          IF (activecellPressure(nx+1,jy,jz) == 0) THEN  
            AddPressureX =  2.0d0*ro(jx,jy,jz)*harx(jx,jy,jz)*pres(nx+1,jy,jz)/(dxx(jx)*dxx(jx))   
          ELSE
            AddPressureX = 0.0d0
          END IF 
          IF (activecellPressure(jx,0,jz) == 0) THEN  
            AddPressureY =  2.0d0*ro(jx,jy,jz)*hary(jx,jy-1,jz)*pres(jx,0,jz)/(dyy(jy)*dyy(jy))   
          ELSE
            AddPressureY = 0.0d0
          END IF 
          IF (activecellPressure(jx,jy,nz+1) == 0) THEN  
            AddPressureZ =  2.0d0*ro(jx,jy,jz)*harz(jx,jy,jz)*pres(jx,jy,nz+1)/(dzz(jx,jy,jz)*dzz(jx,jy,jz))   
          ELSE
            AddPressureZ = 0.0d0
          END IF 

            BvecCrunchP(j) = (pumpterm + tdepend + BodyForceX + BodyForceY + BodyForceZ + AddPressureX + AddPressureY + AddPressureZ) 
            XvecCrunchP(j) = pres(jx,jy,jz)

          END IF

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
            BvecCrunchP(j) = pres(jx,jy,jz)*big 
            XvecCrunchP(j) = pres(jx,jy,jz)
          ELSE

          IF (activecellPressure(0,jy,jz) == 0) THEN
            COORDINATE = 'X'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-1) = -2.0d0*ro(jx,jy,jz)*harx(jx,jy,jz)/(dxx(jx)*dxx(jx))
            coef(1) = -2.0d0*RoAveRight*harx(jx,jy,jz)/(dxx(jx)*(dxx(jx+1)+dxx(jx))) 
            BodyForceX = -COSD(x_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harx(jx,jy,jz) - ro(jx,jy,jz)*ro(jx,jy,jz)*harx(jx,jy,jz))/dxx(jx)
          ELSE
            COORDINATE = 'X'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(1) = -2.0d0*RoAveRight*harx(jx,jy,jz)/(dxx(jx)*(dxx(jx+1)+dxx(jx))) 
            BodyForceX = -COSD(x_angle)*grav*SignGravity*(RoAveRight*RoAveRight*harx(jx,jy,jz) )/dxx(jx)
            coef(-1) = 0.0d0
          END IF

          IF (activecellPressure(jx,ny+1,jz) == 0) THEN
            COORDINATE = 'Y'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-2) = -2.0d0*RoAveLeft*hary(jx,jy-1,jz)/(dyy(jy)*(dyy(jy)+dyy(jy-1))) 
            coef(2)  = -2.0d0*ro(jx,jy,jz)*hary(jx,jy,jz)/(dyy(jy)*dyy(jy)) 
            BodyForceY = -COSD(y_angle)*grav*SignGravity*(ro(jx,jy,jz)*ro(jx,jy,jz)*hary(jx,jy,jz) - RoAveLeft*RoAveLeft*hary(jx,jy-1,jz))/dyy(jy) 
          ELSE
            COORDINATE = 'Y'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-2) = -2.0d0*RoAveLeft*hary(jx,jy-1,jz)/(dyy(jy)*(dyy(jy)+dyy(jy-1))) 
            BodyForceY = -COSD(y_angle)*grav*SignGravity*(-RoAveLeft*RoAveLeft*hary(jx,jy-1,jz))/dyy(jy)    
            coef(2) = 0.0d0     
          END IF

          IF (activecellPressure(jx,jy,nz+1) == 0) THEN
            COORDINATE = 'Z'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-3) = -2.0d0*RoAveLeft*harz(jx,jy,jz-1)/(dzz(jx,jy,jz)*(dzz(jx,jy,jz)+dzz(jx,jy,jz-1))) 
            coef(3) = -2.0d0*ro(jx,jy,jz)*harz(jx,jy,jz)/(dzz(jx,jy,jz)*dzz(jx,jy,jz)) 
            BodyForceZ = -COSD(z_angle)*grav*SignGravity*(ro(jx,jy,jz)*ro(jx,jy,jz)*harz(jx,jy,jz) - RoAveLeft*RoAveLeft*harz(jx,jy,jz-1))/dzz(jx,jy,jz)
          ELSE
            COORDINATE = 'Z'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-3) = -2.0d0*RoAveLeft*harz(jx,jy,jz-1)/(dzz(jx,jy,jz)*(dzz(jx,jy,jz)+dzz(jx,jy,jz-1))) 
            BodyForceZ = -COSD(z_angle)*grav*SignGravity*(-RoAveLeft*RoAveLeft*harz(jx,jy,jz-1))/dzz(jx,jy,jz)
            coef(3) = 0.0d0
          END IF
 
            tdepend = ScaleFactor*visc*ct*por(jx,jy,jz)/(dt*secyr)
            coef(0) = tdepend -( coef(-3) + coef(-2) + coef(-1) + coef(1) + coef(2) + coef(3) )

!!            CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
!!            CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j-nx*ny,coef(-3),INSERT_VALUES,ierr)
!!            CALL MatSetValues(amatP,1,j,1,j+nx*ny,coef(3),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)

            tdepend = visc*ct*por(jx,jy,jz)*pres(jx,jy,jz)/ (dt*secyr)
            pumpterm = visc*ro(jx,jy,jz)*qg(jx,jy,jz)/(secyr*dxx(jx)*dyy(jy)*dzz(jx,jy,jz))

          IF (activecellPressure(0,jy,jz) == 0) THEN  
            AddPressureX =  2.0d0*ro(jx,jy,jz)*harx(jx-1,jy,jz)*pres(0,jy,jz)/(dxx(jx)*dxx(jx))   
          ELSE
            AddPressureX = 0.0d0
          END IF 
          IF (activecellPressure(jx,ny+1,jz) == 0) THEN  
            AddPressureY =  2.0d0*ro(jx,jy,jz)*hary(jx,jy,jz)*pres(jx,ny+1,jz)/(dyy(jy)*dyy(jy))   
          ELSE
            AddPressureY = 0.0d0
          END IF 
          IF (activecellPressure(jx,jy,nz+1) == 0) THEN  
            AddPressureZ =  2.0d0*ro(jx,jy,jz)*harz(jx,jy,jz)*pres(jx,jy,nz+1)/(dzz(jx,jy,jz)*dzz(jx,jy,jz))   
          ELSE
            AddPressureZ = 0.0d0
          END IF 

            BvecCrunchP(j) = (pumpterm + tdepend + BodyForceX + BodyForceY + BodyForceZ + AddPressureX + AddPressureY + AddPressureZ) 
            XvecCrunchP(j) = pres(jx,jy,jz)

          END IF

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
            BvecCrunchP(j) = pres(jx,jy,jz)*big 
            XvecCrunchP(j) = pres(jx,jy,jz)
          ELSE

          IF (activecellPressure(nx+1,jy,jz) == 0) THEN
            COORDINATE = 'X'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-1) = -2.0d0*RoAveLeft*harx(jx-1,jy,jz)/(dxx(jx)*(dxx(jx)+dxx(jx-1))) 
            coef(1) = -2.0d0*ro(jx,jy,jz)*harx(jx,jy,jz)/(dxx(jx)*dxx(jx)) 
            BodyForceX = -COSD(x_angle)*grav*SignGravity*(ro(jx,jy,jz)*ro(jx,jy,jz)*harx(jx,jy,jz) - RoAveLeft*RoAveLeft*harx(jx-1,jy,jz))/dxx(jx)
          ELSE
            COORDINATE = 'X'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-1) = -2.0d0*RoAveLeft*harx(jx-1,jy,jz)/(dxx(jx)*(dxx(jx)+dxx(jx-1))) 
            BodyForceX = -COSD(x_angle)*grav*SignGravity*(-RoAveLeft*RoAveLeft*harx(jx-1,jy,jz))/dxx(jx)
            coef(1) = 0.0d0
          END IF

          IF (activecellPressure(jx,ny+1,jz) == 0) THEN
            COORDINATE = 'Y'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-2) = -2.0d0*RoAveLeft*hary(jx,jy-1,jz)/(dyy(jy)*(dyy(jy)+dyy(jy-1))) 
            coef(2)  = -2.0d0*ro(jx,jy,jz)*hary(jx,jy,jz)/(dyy(jy)*dyy(jy)) 
            BodyForceY = -COSD(y_angle)*grav*SignGravity*(ro(jx,jy,jz)*ro(jx,jy,jz)*hary(jx,jy,jz) - RoAveLeft*RoAveLeft*hary(jx,jy-1,jz))/dyy(jy) 
          ELSE
            COORDINATE = 'Y'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-2) = -2.0d0*RoAveLeft*hary(jx,jy-1,jz)/(dyy(jy)*(dyy(jy)+dyy(jy-1))) 
            BodyForceY = -COSD(y_angle)*grav*SignGravity*(-RoAveLeft*RoAveLeft*hary(jx,jy-1,jz))/dyy(jy)    
            coef(2) = 0.0d0     
          END IF

          IF (activecellPressure(jx,jy,nz+1) == 0) THEN
            COORDINATE = 'Z'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-3) = -2.0d0*RoAveLeft*harz(jx,jy,jz-1)/(dzz(jx,jy,jz)*(dzz(jx,jy,jz)+dzz(jx,jy,jz-1))) 
            coef(3) = -2.0d0*ro(jx,jy,jz)*harz(jx,jy,jz)/(dzz(jx,jy,jz)*dzz(jx,jy,jz)) 
            BodyForceZ = -COSD(z_angle)*grav*SignGravity*(ro(jx,jy,jz)*ro(jx,jy,jz)*harz(jx,jy,jz) - RoAveLeft*RoAveLeft*harz(jx,jy,jz-1))/dzz(jx,jy,jz)
          ELSE
            COORDINATE = 'Z'
            CALL AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            coef(-3) = -2.0d0*RoAveLeft*harz(jx,jy,jz-1)/(dzz(jx,jy,jz)*(dzz(jx,jy,jz)+dzz(jx,jy,jz-1))) 
            BodyForceZ = -COSD(z_angle)*grav*SignGravity*(-RoAveLeft*RoAveLeft*harz(jx,jy,jz-1))/dzz(jx,jy,jz)
            coef(3) = 0.0d0
          END IF
 
            tdepend = ScaleFactor*visc*ct*por(jx,jy,jz)/(dt*secyr)
            coef(0) = tdepend -( coef(-3) + coef(-2) + coef(-1) + coef(1) + coef(2) + coef(3) )

            CALL MatSetValues(amatP,1,j,1,j-1,coef(-1),INSERT_VALUES,ierr)
!!            CALL MatSetValues(amatP,1,j,1,j+1,coef(1),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j-nx,coef(-2),INSERT_VALUES,ierr)
!!            CALL MatSetValues(amatP,1,j,1,j+nx,coef(2),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j-nx*ny,coef(-3),INSERT_VALUES,ierr)
!!            CALL MatSetValues(amatP,1,j,1,j+nx*ny,coef(3),INSERT_VALUES,ierr)
            CALL MatSetValues(amatP,1,j,1,j,coef(0),INSERT_VALUES,ierr)

            tdepend = visc*ct*por(jx,jy,jz)*pres(jx,jy,jz)/ (dt*secyr)
            pumpterm = visc*ro(jx,jy,jz)*qg(jx,jy,jz)/(secyr*dxx(jx)*dyy(jy)*dzz(jx,jy,jz))

          IF (activecellPressure(nx+1,jy,jz) == 0) THEN  
            AddPressureX =  2.0d0*ro(jx,jy,jz)*harx(jx,jy,jz)*pres(nx+1,jy,jz)/(dxx(jx)*dxx(jx))   
          ELSE
            AddPressureX = 0.0d0
          END IF 
          IF (activecellPressure(jx,ny+1,jz) == 0) THEN  
            AddPressureY =  2.0d0*ro(jx,jy,jz)*hary(jx,jy,jz)*pres(jx,ny+1,jz)/(dyy(jy)*dyy(jy))   
          ELSE
            AddPressureY = 0.0d0
          END IF 
          IF (activecellPressure(jx,jy,nz+1) == 0) THEN  
            AddPressureZ =  2.0d0*ro(jx,jy,jz)*harz(jx,jy,jz)*pres(jx,jy,nz+1)/(dzz(jx,jy,jz)*dzz(jx,jy,jz))   
          ELSE
            AddPressureZ = 0.0d0
          END IF 

            BvecCrunchP(j) = (pumpterm + tdepend + BodyForceX + BodyForceY + BodyForceZ + AddPressureX + AddPressureY + AddPressureZ) 
            XvecCrunchP(j) = pres(jx,jy,jz)

          END IF


  END IF

END IF

CALL MatAssemblyBegin(amatP,MAT_FINAL_ASSEMBLY,ierr)
CALL MatAssemblyEnd(amatP,MAT_FINAL_ASSEMBLY,ierr)

!----------------------------------------------------------------------------------
RETURN
END SUBROUTINE pressure


