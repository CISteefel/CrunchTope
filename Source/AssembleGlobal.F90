
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
    
SUBROUTINE AssembleGlobal(nx,ny,nz,ncomp,nspec,nkin,nrct,ngas,ikin,  &
    nexchange,nexch_sec,nsurf,nsurf_sec,npot,ndecay,nn,delt,time,&
    user,amatpetsc,nBoundaryConditionZone)
USE crunchtype
USE params
USE runtime
USE concentration
USE mineral
USE solver
USE medium
USE transport
USE temperature
USE ReadFlow
USE flow

#include "petsc/finclude/petscmat.h"
USE petscmat

#include "petsc/finclude/petsc.h"

!**************************** End PETSc include statements **************
IMPLICIT NONE

INTEGER(I4B), INTENT(IN)                      :: nx
INTEGER(I4B), INTENT(IN)                      :: ny
INTEGER(I4B), INTENT(IN)                      :: nz
INTEGER(I4B), INTENT(IN)                      :: ncomp
INTEGER(I4B), INTENT(IN)                      :: nspec
INTEGER(I4B), INTENT(IN)                      :: nkin
INTEGER(I4B), INTENT(IN)                      :: nrct
INTEGER(I4B), INTENT(IN)                      :: ngas
INTEGER(I4B), INTENT(IN)                      :: ikin
INTEGER(I4B), INTENT(IN)                      :: nexchange
INTEGER(I4B), INTENT(IN)                      :: nexch_sec
INTEGER(I4B), INTENT(IN)                      :: nsurf
INTEGER(I4B), INTENT(IN)                      :: nsurf_sec
INTEGER(I4B), INTENT(IN)                      :: npot
INTEGER(I4B), INTENT(IN)                      :: ndecay
INTEGER(I4B), INTENT(IN)                      :: nn

REAL(DP), INTENT(IN)                          :: delt
REAL(DP), INTENT(IN)                          :: time

!********************** Add in PETSC declarations for f90 variables *****
INTEGER(I4B)                                                     :: ipetsc
INTEGER(I4B)                                                     :: ierr
INTEGER(I4B)                                                     :: irow
INTEGER(I4B)                                                     :: jcol
!************************* End PETSc declarations ****************************

!  Internal arrays

!  Internal real variables

REAL(DP)                                      :: r
REAL(DP)                                      :: wtfactor
REAL(DP)                                      :: satgas
REAL(DP)                                      :: df

REAL(DP)                                      :: sumchgbd
REAL(DP)                                      :: fgradsum
REAL(DP)                                      :: fjwtchg
REAL(DP)                                      :: term1
REAL(DP)                                      :: term2
REAL(DP)                                      :: term3
REAL(DP)                                      :: term3a
REAL(DP)                                      :: term3b
REAL(DP)                                      :: term3c
REAL(DP)                                      :: term3d
REAL(DP)                                      :: source
REAL(DP)                                      :: fgradsume
REAL(DP)                                      :: fgradsumw
REAL(DP)                                      :: fgradsumn
REAL(DP)                                      :: fgradsums
REAL(DP)                                      :: sumrct
REAL(DP)                                      :: sumkin
REAL(DP)                                      :: rxnmin
REAL(DP)                                      :: rxnaq
REAL(DP)                                      :: aq_accum
REAL(DP)                                      :: source_jac
REAL(DP)                                      :: ex_accum
REAL(DP)                                      :: gas_accum
REAL(DP)                                      :: ex_transport
REAL(DP)                                      :: gas_transport
REAL(DP)                                      :: surf_accum
REAL(DP)                                      :: surf_transport
REAL(DP)                                      :: pot_accum
REAL(DP)                                      :: pot_transport
REAL(DP)                                      :: correct
REAL(DP)                                      :: AqueousToBulk
REAL(DP)                                      :: pi

!  Internal integers

INTEGER(I4B)                                  :: neqn
INTEGER(I4B)                                  :: nxyz
INTEGER(I4B)                                  :: ntotal
INTEGER(I4B)                                  :: jx
INTEGER(I4B)                                  :: jy
INTEGER(I4B)                                  :: jz
INTEGER(I4B)                                  :: j
INTEGER(I4B)                                  :: i
INTEGER(I4B)                                  :: jdum
INTEGER(I4B)                                  :: i2
INTEGER(I4B)                                  :: ix
INTEGER(I4B)                                  :: ix2
INTEGER(I4B)                                  :: is
INTEGER(I4B)                                  :: is2
INTEGER(I4B)                                  :: k
INTEGER(I4B)                                  :: ir
INTEGER(I4B)                                  :: np
INTEGER(I4B)                                  :: ind2
INTEGER(I4B)                                  :: ind
INTEGER(I4B)                                  :: npt
INTEGER(I4B)                                  :: npt2
INTEGER(I4B)                                  :: ik
INTEGER(I4B)                                  :: nrow
INTEGER(I4B)                                  :: ncol
INTEGER(I4B)                                  :: ns
INTEGER(I4B)                                  :: nex
INTEGER(I4B)                                  :: jpoint
INTEGER(I4B)                                  :: round

INTEGER(I4B)                                  :: nbnd
INTEGER(I4B)                                  :: nBoundaryConditionZone

INTEGER(I4B)                                  :: npz

REAL(DP)                                      :: sqrt_sion
REAL(DP)                                      :: sum
REAL(DP)                                      :: rotemp
REAL(DP)                                      :: portemp
REAL(DP)                                      :: satl
REAL(DP)                                      :: faraday
REAL(DP)                                      :: gramsperL
REAL(DP)                                      :: delta_z
REAL(DP)                                      :: CellVolume
REAL(DP)                                      :: MultiplyCell
REAL(DP)                                      :: Retardation
REAL(DP)                                      :: HyperbolicSine
REAL(DP)                                      :: volMinimum

REAL(DP)        :: term2_deriv_bulk
REAL(DP)        :: term3_deriv_bulk
REAL(DP)        :: term4_deriv_bulk
REAL(DP)        :: fanalyt_w
REAL(DP)        :: term1_w
REAL(DP)        :: check1
REAL(DP)        :: check2
REAL(DP)        :: check3
REAL(DP)        :: check4
REAL(DP)        :: qgdum
REAL(DP)                                                  :: A_transpi
REAL(DP)                                                  :: coeff_immo

!! Time normalized used if time series only defined for 1 representative year:
REAL(DP)        :: time_norm
REAL(DP), DIMENSION(:), ALLOCATABLE                   :: temp_dum
!********************* PETSc declarations ********************************
PetscFortranAddr                                                    user(6)
Mat                                                                 amatpetsc
!******************end PETSc declarations *********************************

coeff_immo = 1.0

jz = 1
if (delt < 1.0e-20) then
 !! write(*,*) ' Delt = ', delt
  !!!read(*,*)
end if
r = 1.0/delt      
wtfactor = r
pi = DACOS(-1.0d0)

faraday = 96485.0

!***************insert PETSc index initialization ******************

  ipetsc = 0

!*******************end PETSc index initialization ***************


neqn = ncomp + nsurf + nexchange + npot
nxyz = nx*ny*nz
ntotal = 0
surf_accum = 0.0d0
ex_accum = 0.0d0

!!!xn = 0.0d0

IF (nxyz == nx .AND. ihindmarsh == 1 .AND. nxyz /= 1) THEN
  indexx = 0
  bbh = 0.0
  cch = 0.0
  aah = 0.0
  yh  = 0.0
END IF

fxx = 0.0

IF (ierode /= 1) THEN
  IF (nsurf > 0) THEN
    CALL SurfaceComplex(ncomp,nsurf,nsurf_sec,nx,ny,nz)
  END IF
END IF

IF (npot > 0) THEN
  CALL jacpotential(ncomp,nsurf,nsurf_sec,npot,nx,ny,nz)
END IF

TempFlux = 0.0d0

! IF (RunTempts) THEN
!   IF (ALLOCATED(temp_dum)) THEN
!     DEALLOCATE(temp_dum)
!   END IF
!   ALLOCATE(temp_dum(nb_temp_ts))

!   IF (TS_1year) THEN
!     time_norm=time-floor(time)
!     DO i=1,nb_temp_ts
!   CALL  interp3(time_norm,delt,t_temp_ts,temp_ts(i,:),temp_dum(i),size(temp_ts(i,:)))
!     END DO
!     END IF

!       DO jz = 1,nz
!       DO jy = 1,ny
!       DO jx = 1,nx
!       DO i = 1,nb_temp_ts
!           IF (temp_region(jx,jy,jz) == reg_temp_ts(i)) THEN
!             t(jx,jy,jz) = temp_dum(i)
!           ENDIF
!       END DO
!       END DO
!       END DO
!       END DO
! ENDIF

!!!  Do the boundaries first

!!!  ***** Nernst-Planck ******
IF (species_diffusion) THEN

  IF (nBoundaryConditionZone > 0) THEN

!!!  Sweep the boundaries, BC by grid cells

    jx = 0
    jz = 1
    DO jy = 1,ny
      CALL bd_diffuse_by_grid(ncomp,nspec,jx,jy,jz,sumchgbd)
      DO i = 1,ncomp
        s_dsp(i,jx,jy,jz) = sdsp(i)
        s_chg(i,jx,jy,jz) = schg(i)
      END DO
      sumwtchg(jx,jy,jz) = sumchgbd
    END DO

    jx = nx+1
    jz = 1
    DO jy = 1,ny
      CALL bd_diffuse_by_grid(ncomp,nspec,jx,jy,jz,sumchgbd)
      DO i = 1,ncomp
        s_dsp(i,jx,jy,jz) = sdsp(i)
        s_chg(i,jx,jy,jz) = schg(i)
      END DO
      sumwtchg(jx,jy,jz) = sumchgbd
    END DO

    jy = 0
    jz = 1
    DO jx = 1,nx
      CALL bd_diffuse_by_grid(ncomp,nspec,jx,jy,jz,sumchgbd)
      DO i = 1,ncomp
        s_dsp(i,jx,jy,jz) = sdsp(i)
        s_chg(i,jx,jy,jz) = schg(i)
      END DO
      sumwtchg(jx,jy,jz) = sumchgbd
    END DO

    jy = ny+1
    jz = 1
    DO jx = 1,nx
      CALL bd_diffuse_by_grid(ncomp,nspec,jx,jy,jz,sumchgbd)
      DO i = 1,ncomp
        s_dsp(i,jx,jy,jz) = sdsp(i)
        s_chg(i,jx,jy,jz) = schg(i)
      END DO
      sumwtchg(jx,jy,jz) = sumchgbd
    END DO
    
  ELSE    !! Conventional BC by face

    nbnd = 1
    CALL bd_diffuse(ncomp,nspec,nbnd,sumchgbd)
    DO jy = 1,ny
       DO i = 1,ncomp
          s_dsp(i,0,jy,jz) = sdsp(i)
          s_chg(i,0,jy,jz) = schg(i)
       END DO
       sumwtchg(0,jy,jz) = sumchgbd
    END DO
 
    nbnd = 2
    CALL bd_diffuse(ncomp,nspec,nbnd,sumchgbd)
    DO jy = 1,ny
       DO i = 1,ncomp
          s_dsp(i,nx+1,jy,jz) = sdsp(i)
          s_chg(i,nx+1,jy,jz) = schg(i)
       END DO
       sumwtchg(nx+1,jy,jz) = sumchgbd
    END DO
 
    nbnd = 3
    CALL bd_diffuse(ncomp,nspec,nbnd,sumchgbd)
    DO jx = 1,nx
       DO i = 1,ncomp
          s_dsp(i,jx,0,jz) = sdsp(i)
          s_chg(i,jx,0,jz) = schg(i)
       END DO
       sumwtchg(jx,0,jz) = sumchgbd
    END DO

    nbnd = 4
    CALL bd_diffuse(ncomp,nspec,nbnd,sumchgbd)
    DO jx = 1,nx
       DO i = 1,ncomp
          s_dsp(i,jx,ny+1,jz) = sdsp(i)
          s_chg(i,jx,ny+1,jz) = schg(i)
       END DO
       sumwtchg(jx,ny+1,jz) = sumchgbd
    END DO
    
  END IF

END IF
!!!  ***** Nernst-Planck ******

jz = 1
DO jy = 1,ny
  DO jx = 1,nx
    j = (jz-1)*nx*ny + (jy-1)*nx + jx
   
    jpoint = jinit(jx,jy,jz)
    portemp = por(jx,jy,jz)
    rotemp = ro(jx,jy,jz)

    Retardation = 0.001d0*SolidDensity(jinit(jx,jy,jz))*(1.0-por(jx,jy,jz))/por(jx,jy,jz)

    CALL AqueousToBulkConvert(jx,jy,jz,AqueousToBulk)

    satl = 0.5*(satliq(jx,jy,jz) + satliqold(jx,jy,jz) )
    satgas = 1.0d0 - satl
    IF (cylindrical) THEN
      CellVolume = dyy(jy)*pi*( (x(jx)+dxx(jx)/2.0d0 )**2.0d0 - ( x(jx)-dxx(jx)/2.0d0 )**2.0d0  )
      df = 1.0d0
      MultiplyCell = CellVolume
      IF (CylindricalDivideVolume) THEN
        df = 1.0/CellVolume
        MultiplyCell = 1.0
      END IF
    ELSE IF (spherical) THEN
      CellVolume = (4.0d0/3.0d0)*pi*( (x(jx) + 0.5d0*dxx(jx))**3.0d0 - (x(jx) - 0.5d0*dxx(jx))**3.0d0  )
       df = 1.0d0
       MultiplyCell = CellVolume
       IF (CylindricalDivideVolume) THEN
        df = 1.0/CellVolume
        MultiplyCell = 1.0
       END IF
!!      df = 1.0/CellVolume
!!      MultiplyCell = 1.0
    ELSE
      CellVolume = dxx(jx)*dyy(jy)*dzz(jx,jy,jz)
       df = 1.0d0
       MultiplyCell = CellVolume
!!      df = 1.0/CellVolume
!!      MultiplyCell = 1.0
    END IF
    
    alf = 0.0d0
    
    IF (ierode /= 1) THEN
      IF (nexchange > 0) THEN
        CALL exchange(ncomp,nexchange,nexch_sec,jx,jy,jz)
      END IF
    END IF
    
    CALL FxTopeGlobal(nx,ny,nz,ncomp,nexchange,nexch_sec,nsurf,nsurf_sec,nrct,nspec,  &
        ngas,neqn,delt,jx,jy,jz,nBoundaryConditionZone)
    
    IF (ierode == 1) THEN
      CALL ex_activity(ncomp,nexchange,nexch_sec,jx,jy,jz)
      CALL jacexchange(ncomp,nexchange,nexch_sec,neqn,jx,jy,jz)
    END IF
    
!!    IF (nexchange > 0) THEN
!!      CALL tot_ex(ncomp,nexchange,nexch_sec,jx,jy,jz)
!!    END IF
    
    IF (ierode /= 1) THEN
      CALL jac_exchange_local(ncomp,nexchange,nexch_sec,nsurf,nsurf_sec,neqn,jx,jy,jz)
      IF (nsurf > 0) THEN
        CALL jacsurf_local(ncomp,nsurf,nsurf_sec,jx,jy,jz)
      END IF
    END IF
        
    IF (nrct > 0) THEN
      CALL reaction(ncomp,nkin,nrct,nspec,nexchange,nsurf,ndecay,jx,jy,jz,delt,time)
      round = 1
      CALL jacmin(ncomp,nspec,nexchange,nsurf,nkin,nrct,jx,jy,jz,time,round)
    END IF
    
    IF (ikin > 0) THEN
      CALL reactkin(ncomp,nspec,nrct,ikin,jx,jy,jz,AqueousToBulk,time)
      CALL jacrkin(ncomp,nspec,nrct,ikin,jx,jy,jz,AqueousToBulk)
    END IF

    IF (activecell(jx,jy,jz) == 0) GO TO 800
    
    IF (nx == 1) GO TO 100
    
!    IF (species_diffusion) THEN
!      dgradw1 = 0.0
!      dgradw2 = 0.0
!      dgrade1 = 0.0
!      dgrade2 = 0.0
!      DO i = 1,ncomp
!        dgradw1 = dgradw1 + chg(i)*a_d(jx,jy,jz)*(s_dsp(i,jx-1,jy,jz))
!        dgradw2 = dgradw2 - chg(i)*a_d(jx,jy,jz)*(s_dsp(i,jx,jy,jz))
!        dgrade1 = dgrade1 + chg(i)*c_d(jx,jy,jz)*(s_dsp(i,jx+1,jy,jz))
!        dgrade2 = dgrade2 - chg(i)*c_d(jx,jy,jz)*(s_dsp(i,jx,jy,jz))
!      END DO
!      dgradw = dgradw1 + dgradw2
!      dgrade = dgrade1 + dgrade2
!      sumsigma_w = 0.5*(dstar(jx,jy,jz)*sumwtchg(jx,jy,jz) + dstar(jx-1,jy,jz)*sumwtchg(jx-1,jy,jz))
!      sumsigma_e = 0.5*(dstar(jx,jy,jz)*sumwtchg(jx,jy,jz) + dstar(jx+1,jy,jz)*sumwtchg(jx+1,jy,jz))
      
!      DO i = 1,ncomp
!        sigma_w(i) = 0.5*(dstar(jx,jy,jz)*s_chg(i,jx,jy,jz) + dstar(jx-1,jy,jz)*s_chg(i,jx-1,jy,jz))
!        sigma_e(i) = 0.5*(dstar(jx,jy,jz)*s_chg(i,jx,jy,jz) + dstar(jx+1,jy,jz)*s_chg(i,jx+1,jy,jz))
!      END DO
!    END IF
    
    IF (nxyz == nx .AND. ihindmarsh == 1 .AND. nxyz /= 1) THEN  ! Use Hindmarsh routine
  
      IF (jx /= 1) THEN
        jdum=jx-1
        
        IF (species_diffusion) THEN
          DO i2 = 1,ncomp
            fgradsum = 0.0
            fjwtchg = 0.0
            DO i = 1,ncomp
              fgradsum = fgradsum + a_d(jx,jy,jz)*chg(i)*fjac_d(i2,i,jdum,jy,jz)
              fjwtchg = fjwtchg + 0.5*dstar(jdum,jy,jz)*chg(i)*fjac_chg(i2,i,jdum,jy,jz)
            END DO
            DO i = 1,ncomp
              term1 = a_d(jx,jy,jz)*fjac_d(i2,i,jdum,jy,jz)
              term2 = 0.5*dstar(jdum,jy,jz)*fjac_chg(i2,i,jdum,jy,jz)
              term3 = -1.0/(sumsigma_w*sumsigma_w) * fjwtchg
              
              term2_deriv_bulk = term2*dgradw/sumsigma_w
              term3_deriv_bulk = term3*sigma_w(i)*dgradw
              term4_deriv_bulk = fgradsum*sigma_w(i)/sumsigma_w
              
              fanalyt(i,i2) = df*(term1     &
               - ( term2*dgradw/sumsigma_w + term3*sigma_w(i)*dgradw + fgradsum*sigma_w(i)/sumsigma_w ))
!              fanalyt(i,i2) = df*term1  
             continue
            END DO
          END DO
        END IF

        
        DO i = 1,ncomp
          IF (immobile_species(i) == 1) THEN
          coeff_immo = 0.0
          ELSE
          coeff_immo = 1.0
          ENDIF

          DO i2 = 1,ncomp
            cch(i,i2,jx) = xgram(jdum,jy,jz)*df*a(jx,jy,jz)*fjac(i2,i,jdum,jy,jz)!*coeff_immo
          END DO

          IF (ierode == 1) THEN
            DO i2 = 1,ncomp
              cch(i,i2,jx) = cch(i,i2,jx) + df*abu(jx,jy,jz)*fch(i,i2,jdum,jy,jz)
            END DO
          END IF

          IF (isaturate == 1) THEN
            DO i2 = 1,ncomp
              cch(i,i2,jx) = cch(i,i2,jx) + df*ag(jx,jy,jz)*fgas(i2,i,jdum,jy,jz)
            END DO           
          END IF

          IF (species_diffusion) THEN
            DO i2 = 1,ncomp
              cch(i,i2,jx) = cch(i,i2,jx) + fanalyt(i,i2)
            END DO
          END IF
          
          IF (ierode == 1) THEN
            DO ix2 = 1,nexchange
              cch(i,ix2+ncomp,jx) = df*abu(jx,jy,jz)*fch(i,ix2+ncomp,jdum,jy,jz)
            END DO
            DO is2 = 1,nsurf
              cch(i,is2+ncomp+nexchange,jx) = df*abu(jx,jy,jz)*  &
                  fch(i,is2+ncomp+nexchange,jdum,jy,jz)
            END DO
            DO npt2 = 1,npot
              cch(i,npt2+ncomp+nexchange+nsurf,jx) = df*abu(jx,jy,jz)*  &
                  fjpotncomp(npt2,i,jdum,jy,jz)
            END DO
          END IF
          
        END DO
        
        IF (ierode == 1) THEN
          DO is = 1,nsurf
            DO i2 = 1,ncomp
              cch(is+ncomp+nexchange,i2,jx) = df*abu(jx,jy,jz)*fsurf(is,i2,jdum,jy,jz)
            END DO
            DO is2 = 1,nsurf
              cch(is+ncomp+nexchange,is2+ncomp+nexchange,jx) = df*abu(jx,jy,jz)*  &
                  fsurf(is,is2+ncomp,jdum,jy,jz)
            END DO
            DO npt2 = 1,npot
              cch(is+ncomp+nexchange,npt2+ncomp+nexchange+nsurf,jx) =    &
                 df*abu(jx,jy,jz)*fjpotnsurf(npt2,is,jdum,jy,jz)
            END DO
          END DO
        END IF

      END IF
      
      IF (jx /= nx) THEN
        jdum=jx+1
        
        IF (species_diffusion) THEN
          DO i2 = 1,ncomp
            fgradsum = 0.0
            fjwtchg = 0.0
            DO i = 1,ncomp
              fgradsum = fgradsum + c_d(jx,jy,jz)*chg(i)*fjac_d(i2,i,jdum,jy,jz)
              fjwtchg = fjwtchg + 0.5*dstar(jdum,jy,jz)*chg(i)*fjac_chg(i2,i,jdum,jy,jz)
            END DO
            DO i = 1,ncomp
              term1 = c_d(jx,jy,jz)*fjac_d(i2,i,jdum,jy,jz)
              term2 = 0.5*dstar(jdum,jy,jz)*fjac_chg(i2,i,jdum,jy,jz)
              term3 = -1.0/(sumsigma_e*sumsigma_e) * fjwtchg
              fanalyt(i,i2) = df*(term1  &
                - (term2*dgrade/sumsigma_e + term3*sigma_e(i)*dgrade + fgradsum*sigma_e(i)/sumsigma_e ))
!              fanalyt(i,i2) = df*term1 
            END DO
          END DO
        END IF
        
        DO i = 1,ncomp
        IF (immobile_species(i) == 1) THEN
        coeff_immo = 0.0
        ELSE
        coeff_immo = 1.0
        ENDIF

          
          DO i2 = 1,ncomp
            bbh(i,i2,jx) = xgram(jdum,jy,jz)*df*c(jx,jy,jz)*fjac(i2,i,jdum,jy,jz)!*coeff_immo
          END DO

          IF (ierode == 1) THEN
            DO i2 = 1,ncomp
              bbh(i,i2,jx) = bbh(i,i2,jx) + df*cbu(jx,jy,jz)*fch(i,i2,jdum,jy,jz)
            END DO
          END IF

          IF (isaturate == 1) THEN
            DO i2 = 1,ncomp
              bbh(i,i2,jx) = bbh(i,i2,jx) + df*cg(jx,jy,jz)*fgas(i2,i,jdum,jy,jz)
            END DO
          END IF

          IF (species_diffusion) THEN
            DO i2 = 1,ncomp
              bbh(i,i2,jx) = bbh(i,i2,jx) + fanalyt(i,i2)
            END DO
          END IF
          
          IF (ierode == 1) THEN
            DO ix2 = 1,nexchange
              bbh(i,ix2+ncomp,jx) = df*cbu(jx,jy,jz)*fch(i,ix2+ncomp,jdum,jy,jz)
            END DO
            DO is2 = 1,nsurf
              bbh(i,is2+ncomp+nexchange,jx) = df*cbu(jx,jy,jz)*  &
                  fch(i,is2+ncomp+nexchange,jdum,jy,jz)
            END DO
            DO npt2 = 1,npot
              bbh(i,npt2+ncomp+nexchange+nsurf,jx) = df*cbu(jx,jy,jz)*  &
                  fjpotncomp(npt2,i,jdum,jy,jz)
            END DO
          END IF
          
        END DO
        
        IF (ierode == 1) THEN
          DO is = 1,nsurf
            DO i2 = 1,ncomp
              bbh(is+ncomp+nexchange,i2,jx) = df*cbu(jx,jy,jz)*fsurf(is,i2,jdum,jy,jz)
            END DO
            DO is2 = 1,nsurf
              bbh(is+ncomp+nexchange,is2+ncomp+nexchange,jx) = df*cbu(jx,jy,jz)*  &
                  fsurf(is,is2+ncomp,jdum,jy,jz)
            END DO
            DO npt2 = 1,npot
              bbh(is+ncomp+nexchange,npt2+ncomp+nexchange+nsurf,jx) = df*cbu(jx,jy,jz)*  &
                  fjpotnsurf(npt2,is,jdum,jy,jz)
            END DO
          END DO
        END IF
        
      END IF
      
    ELSE     !  Case where NY .ne. 1 (2D problem)
      
      IF (jx /= 1) THEN
        jdum=jx-1

        IF (species_diffusion) THEN
          DO i2 = 1,ncomp
            fgradsum = 0.0
            fjwtchg = 0.0
            DO i = 1,ncomp
              fgradsum = fgradsum + a_d(jx,jy,jz)*chg(i)*fjac_d(i2,i,jdum,jy,jz)
              fjwtchg = fjwtchg + 0.5*dstar(jdum,jy,jz)*chg(i)*fjac_chg(i2,i,jdum,jy,jz)
            END DO
            DO i = 1,ncomp
              term1 = a_d(jx,jy,jz)*fjac_d(i2,i,jdum,jy,jz)
              term2 = 0.5*dstar(jdum,jy,jz)*fjac_chg(i2,i,jdum,jy,jz)
              term3 = -1.0/(sumsigma_w*sumsigma_w) * fjwtchg
              fanalyt(i,i2) = df*(term1      &
                - (term2*dgradw/sumsigma_w + term3*sigma_w(i)*dgradw + fgradsum*sigma_w(i)/sumsigma_w ))
!              fanalyt(i,i2) = df*term1
            END DO
          END DO
        END IF
        
        DO i = 1,ncomp
          ind = (j-1)*(neqn) + i
          
          DO i2 = 1,ncomp
            ind2 = i2
            alf(ind2,i,1) = xgram(jdum,jy,jz)*df*a(jx,jy,jz)*fjac(i2,i,jdum,jy,jz)!*coeff_immo
          END DO

          IF (ierode == 1) THEN
            DO i2 = 1,ncomp
              ind2 = i2
              alf(ind2,i,1) = alf(ind2,i,1) + df*abu(jx,jy,jz)*fch(i,i2,jdum,jy,jz)
            END DO
          END IF

          IF (isaturate == 1) THEN
            DO i2 = 1,ncomp
              ind2 = i2
              alf(ind2,i,1) = alf(ind2,i,1) + df*ag(jx,jy,jz)*fgas(i2,i,jdum,jy,jz)
            END DO
          END IF

          IF (species_diffusion) THEN
            DO i2 = 1,ncomp
              ind2 = i2
              alf(ind2,i,1) = alf(ind2,i,1) + fanalyt(i,i2)
            END DO
          END IF
          
          IF (ierode == 1) THEN
            DO ix2 = 1,nexchange
              ind2 = ix2+ncomp
              alf(ind2,i,1) = df*abu(jx,jy,jz)*fch(i,ix2+ncomp,jdum,jy,jz)
            END DO
            DO is2 = 1,nsurf
              ind2 = is2+ncomp+nexchange
              alf(ind2,i,1) = df*abu(jx,jy,jz)*         &
                 fch(i,is2+ncomp+nexchange,jdum,jy,jz)
            END DO
            DO npt2 = 1,npot
              ind2 = npt2+ncomp+nexchange+nsurf
              alf(ind2,i,1) = df*abu(jx,jy,jz)*         &
                 fjpotncomp(npt2,i,jdum,jy,jz)
            END DO
          END IF

        END DO
        
        IF (ierode == 1) THEN
          DO is = 1,nsurf
            ind = (j-1)*(neqn) + is+ncomp+nexchange
            DO i2 = 1,ncomp
              ind2 = i2
              alf(ind2,is+ncomp+nexchange,1) = df*abu(jx,jy,jz)*   &
                 fsurf(is,i2,jdum,jy,jz)
            END DO
            DO is2 = 1,nsurf
              ind2 = is2+ncomp+nexchange
              alf(ind2,is+ncomp+nexchange,1) = df*abu(jx,jy,jz)*   &
                 fsurf(is,is2+ncomp,jdum,jy,jz)
            END DO
            DO npt2 = 1,npot
              ind2 = npt2+ncomp+nexchange+nsurf
              alf(ind2,is+ncomp+nexchange,1) = df*abu(jx,jy,jz)*   &
                 fjpotnsurf(npt2,is,jdum,jy,jz)
            END DO
          END DO
        END IF
        
      END IF
      
      IF (jx /= nx) THEN
        jdum=jx+1

        IF (species_diffusion) THEN
          DO i2 = 1,ncomp
            fgradsum = 0.0
            fjwtchg = 0.0
            DO i = 1,ncomp
              fgradsum = fgradsum + c_d(jx,jy,jz)*chg(i)*fjac_d(i2,i,jdum,jy,jz)
              fjwtchg = fjwtchg + 0.5*dstar(jdum,jy,jz)*chg(i)*fjac_chg(i2,i,jdum,jy,jz)
            END DO
            DO i = 1,ncomp
              term1 = c_d(jx,jy,jz)*fjac_d(i2,i,jdum,jy,jz)
              term2 = 0.5*dstar(jdum,jy,jz)*fjac_chg(i2,i,jdum,jy,jz)
              term3 = -1.0/(sumsigma_e*sumsigma_e) * fjwtchg
              fanalyt(i,i2) = df*(term1        &
                - (term2*dgrade/sumsigma_e + term3*sigma_e(i)*dgrade + fgradsum*sigma_e(i)/sumsigma_e ))
!              fanalyt(i,i2) = df*term1
            END DO
          END DO
        END IF

        DO i = 1,ncomp
        IF (immobile_species(i) == 1) THEN
        coeff_immo = 0.0
        ELSE
        coeff_immo = 1.0
        ENDIF
          ind = (j-1)*(neqn) + i

          DO i2 = 1,ncomp
            ind2 = i2
            alf(ind2,i,3) = xgram(jdum,jy,jz)*df*c(jx,jy,jz)*fjac(i2,i,jdum,jy,jz)!*coeff_immo
          END DO

          IF (ierode == 1) THEN
            DO i2 = 1,ncomp
              ind2 = i2
              alf(ind2,i,3) = alf(ind2,i,3) + df*cbu(jx,jy,jz)*      &
                 fch(i,i2,jdum,jy,jz)
            END DO
          END IF

          IF (isaturate == 1) THEN
            DO i2 = 1,ncomp
              ind2 = i2
              alf(ind2,i,3) = alf(ind2,i,3) + df*cg(jx,jy,jz)*fgas(i2,i,jdum,jy,jz)
            END DO
          END IF

          IF (species_diffusion) THEN
            DO i2 = 1,ncomp
              ind2 = i2
              alf(ind2,i,3) = alf(ind2,i,3) + fanalyt(i,i2)
            END DO
          END IF
          
          IF (ierode == 1) THEN
            DO ix2 = 1,nexchange
              ind2 =ix2+ncomp
              alf(ind2,i,3) = df*cbu(jx,jy,jz)*fch(i,ix2+ncomp,jdum,jy,jz)
            END DO
            DO is2 = 1,nsurf
              ind2 = is2+ncomp+nexchange
              alf(ind2,i,3) = df*cbu(jx,jy,jz)*              &
                 fch(i,is2+ncomp+nexchange,jdum,jy,jz)
            END DO
            DO npt2 = 1,npot
              ind2 = npt2+ncomp+nexchange+nsurf
              alf(ind2,i,3) = df*cbu(jx,jy,jz)*              &
                 fjpotncomp(npt2,i,jdum,jy,jz)
            END DO
          END IF
          
        END DO
        
        IF (ierode == 1) THEN
          DO is = 1,nsurf
            ind = (j-1)*(neqn) + is+ncomp+nexchange
            DO i2 = 1,ncomp
              ind2 = i2
              alf(ind2,is+ncomp+nexchange,3) = df*cbu(jx,jy,jz)*     &
                 fsurf(is,i2,jdum,jy,jz)
            END DO
            DO is2 = 1,nsurf
              ind2 = is2+ncomp+nexchange
              alf(ind2,is+ncomp+nexchange,3) = df*cbu(jx,jy,jz)*     & 
                 fsurf(is,is2+ncomp,jdum,jy,jz)
            END DO
            DO npt2 = 1,npot
              ind2 = npt2+ncomp+nexchange+nsurf
              alf(ind2,is+ncomp+nexchange,3) = df*cbu(jx,jy,jz)*     & 
                 fjpotnsurf(npt2,is,jdum,jy,jz)
            END DO
          END DO
        END IF
        
      END IF
    END IF
    
    100     CONTINUE
    IF (ny == 1) GO TO 200
    
    IF (jy /= ny) THEN
      jdum=jy+1

      IF (species_diffusion) THEN
        DO i2 = 1,ncomp
          fgradsum = 0.0
          fjwtchg = 0.0
          DO i = 1,ncomp
            fgradsum = fgradsum + d_d(jx,jy,jz)*chg(i)*fjac_d(i2,i,jx,jdum,jz)
            fjwtchg = fjwtchg + 0.5*dstar(jx,jdum,jz)*chg(i)*fjac_chg(i2,i,jx,jdum,jz)
          END DO
          DO i = 1,ncomp
            term1 = d_d(jx,jy,jz)*fjac_d(i2,i,jx,jdum,jz)
            term2 = 0.5*dstar(jx,jdum,jz)*fjac_chg(i2,i,jx,jdum,jz)
            term3 = -1.0/(sumsigma_n*sumsigma_n) * fjwtchg
            fanalyt(i,i2) = df*(term1      &
              - (term2*dgradn/sumsigma_n + term3*sigma_n(i)*dgradn + fgradsum*sigma_n(i)/sumsigma_n ))
!            fanalyt(i,i2) = df*term1
          END DO
        END DO
      END IF

      DO i = 1,ncomp
      IF (immobile_species(i) == 1) THEN
        coeff_immo = 0.0
        ELSE
        coeff_immo = 1.0
        ENDIF
        ind = (j-1)*(neqn) + i

        DO i2 = 1,ncomp
          ind2 = i2
          alf(ind2,i,5) = xgram(jx,jdum,jz)*df*d(jx,jy,jz)*fjac(i2,i,jx,jdum,jz)!*coeff_immo
        END DO
        
        IF (ierode == 1) THEN
          DO i2 = 1,ncomp
            ind2 = i2
            alf(ind2,i,5) = alf(ind2,i,5) + df*dbu(jx,jy,jz)*fch(i,i2,jx,jdum,jz)
          END DO
        END IF

        IF (isaturate == 1) THEN
          DO i2 = 1,ncomp
            ind2 = i2
            alf(ind2,i,5) = alf(ind2,i,5) + df*dg(jx,jy,jz)*fgas(i2,i,jx,jdum,jz)
          END DO
        END IF

        IF (species_diffusion) THEN
          DO i2 = 1,ncomp
            ind2 = i2
            alf(ind2,i,5) = alf(ind2,i,5) + fanalyt(i,i2)
          END DO
        END IF
        
        IF (ierode == 1) THEN
          DO ix2 = 1,nexchange
            ind2 = ix2+ncomp
            alf(ind2,i,5) = df*dbu(jx,jy,jz)*fch(i,ix2+ncomp,jx,jdum,jz)
          END DO
          DO is2 = 1,nsurf
            ind2 = is2+ncomp+nexchange
            alf(ind2,i,5) = df*dbu(jx,jy,jz)*fch(i,is2+ncomp+nexchange,jx,jdum,jz)
          END DO
          DO npt2 = 1,npot
            ind2 = npt2+ncomp+nexchange+nsurf
            alf(ind2,i,5) = df*dbu(jx,jy,jz)*   &
               fjpotncomp(npt2,i,jx,jdum,jz)
          END DO
        END IF
        
      END DO
      
      IF (ierode == 1) THEN
        DO is = 1,nsurf
          ind = (j-1)*neqn + is+ncomp+nexchange
          DO i2 = 1,ncomp
            ind2 = i2
            alf(ind2,is+ncomp+nexchange,5) = df*dbu(jx,jy,jz)*fsurf(is,i2,jx,jdum,jz)
          END DO
          DO is2 = 1,nsurf
            ind2 = is2+ncomp+nexchange
            alf(ind2,is+ncomp+nexchange,5) = df*dbu(jx,jy,jz)*fsurf(is,is2+ncomp,jx,jdum,jz)
          END DO
          DO npt2 = 1,npot
            ind2 = npt2+ncomp+nexchange+nsurf
            alf(ind2,is+ncomp+nexchange,5) = df*dbu(jx,jy,jz)*   &
                 fjpotnsurf(npt2,is,jx,jdum,jz)
          END DO
        END DO
      END IF
      
    END IF

    IF (jy /= 1) THEN
      jdum=jy-1

      IF (species_diffusion) THEN
        DO i2 = 1,ncomp
          fgradsum = 0.0
          fjwtchg = 0.0
          DO i = 1,ncomp
            fgradsum = fgradsum + f_d(jx,jy,jz)*chg(i)*fjac_d(i2,i,jx,jdum,jz)
            fjwtchg = fjwtchg + 0.5*dstar(jx,jdum,jz)*chg(i)*fjac_chg(i2,i,jx,jdum,jz)
          END DO
          DO i = 1,ncomp
            term1 = f_d(jx,jy,jz)*fjac_d(i2,i,jx,jdum,jz)
            term2 = 0.5*dstar(jx,jdum,jz)*fjac_chg(i2,i,jx,jdum,jz)
            term3 = -1.0/(sumsigma_s*sumsigma_s) * fjwtchg
            fanalyt(i,i2) = df*(term1      &
              - (term2*dgrads/sumsigma_s + term3*sigma_s(i)*dgrads + fgradsum*sigma_s(i)/sumsigma_s ))
!            fanalyt(i,i2) = df*term1
          END DO
        END DO
      END IF

      DO i = 1,ncomp
      IF (immobile_species(i) == 1) THEN
        coeff_immo = 0.0
        ELSE
        coeff_immo = 1.0
        ENDIF
        ind = (j-1)*neqn + i

        DO i2 = 1,ncomp
          ind2 = i2
          alf(ind2,i,4) = xgram(jx,jdum,jz)*df*f(jx,jy,jz)*fjac(i2,i,jx,jdum,jz)!*coeff_immo
        END DO

        IF (ierode == 1) THEN
          DO i2 = 1,ncomp
            ind2 = i2
            alf(ind2,i,4) = alf(ind2,i,4) + df*fbu(jx,jy,jz)*fch(i,i2,jx,jdum,jz)
          END DO
        END IF

        IF (isaturate == 1) THEN
          DO i2 = 1,ncomp
            ind2 = i2
            alf(ind2,i,4) = alf(ind2,i,4) + df*fg(jx,jy,jz)*fgas(i2,i,jx,jdum,jz)
          END DO
        END IF

        IF (species_diffusion) THEN
          DO i2 = 1,ncomp
            ind2 = i2
            alf(ind2,i,4) = alf(ind2,i,4) + fanalyt(i,i2)
          END DO
        END IF
        
        IF (ierode == 1) THEN

          DO ix2 = 1,nexchange
            ind2= ix2+ncomp
            alf(ind2,i,4) = df*fbu(jx,jy,jz)*fch(i,ix2+ncomp,jx,jdum,jz)
          END DO
          DO is2 = 1,nsurf
            ind2 =is2+ncomp+nexchange
            alf(ind2,i,4) = df*fbu(jx,jy,jz)*fch(i,is2+ncomp+nexchange,jx,jdum,jz)
          END DO
          DO npt2 = 1,npot
            ind2 =npt2+ncomp+nexchange+nsurf
            alf(ind2,i,4) = df*fbu(jx,jy,jz)*    &
                fjpotncomp(npt2,i,jx,jdum,jz)
          END DO
        END IF
      END DO
      
      IF (ierode == 1) THEN
        DO is = 1,nsurf
          ind = (j-1)*neqn + is+ncomp+nexchange
          DO i2 = 1,ncomp
            ind2 = i2
            alf(ind2,is+ncomp+nexchange,4) = df*fbu(jx,jy,jz)*fsurf(is,i2,jx,jdum,jz)
          END DO
          DO is2 = 1,nsurf
            ind2 = is2+ncomp+nexchange
            alf(ind2,is+ncomp+nexchange,4) = df*fbu(jx,jy,jz)*fsurf(is,is2+ncomp,jx,jdum,jz)
          END DO
          DO npt2 = 1,npot
            ind2 = npt2+ncomp+nexchange+nsurf
            alf(ind2,is+ncomp+nexchange,4) = df*fbu(jx,jy,jz)*   &
                fjpotnsurf(npt2,is,jx,jdum,jz)
          END DO
        END DO
      END IF
      
    END IF
    
200 CONTINUE
    
    source = 0.0d0

    IF ((transpifix .OR. transpitimeseries) .AND. Richards) THEN
      if (ny == 1 .AND. nz == 1) THEN
      A_transpi = dyy(jy) * dzz(jx,jy,jz)
      source = source - xgram(jx,jy,jz)*transpirate_cell(jx)*A_transpi*rotemp/CellVolume
    ENDIF
    ENDIF

    IF (wells) THEN
   
      DO npz = 1,npump(jx,jy,jz)
        IF (qg(npz,jx,jy,jz) > 0.0) THEN       !  Injection well
          CONTINUE                ! Source term on R.H.S.
        ELSE IF (qg(npz,jx,jy,jz) < 0.0) THEN  ! Pumping well, S(i,j) unknown
          source = source + xgram(jx,jy,jz)*qg(npz,jx,jy,jz)*rotemp/CellVolume   !!  GIMRT source term in m^3/year
        ELSE
          CONTINUE
        END IF
      END DO

    ELSE IF (pumptimeseries) THEN
      IF (npump(jx,jy,jz)>0) THEN
        CALL interp3(time,delt,tpump,qgt(:),qg(1,jx,jy,jz),size(qgt(:)))
        ELSE
        qg(1,jx,jy,jz)=0
        END IF
      IF (qg(1,jx,jy,jz) > 0.0) THEN       !  Injection well
        CONTINUE                ! Source term on R.H.S.
      ELSE IF (qg(1,jx,jy,jz) < 0.0) THEN  ! Pumping well, S(i,j) unknown
        source = source + xgram(jx,jy,jz)*qg(1,jx,jy,jz)*rotemp/CellVolume   !!  GIMRT source term in m^3/year
      ELSE
        CONTINUE
      END IF
    END IF
    
    IF (species_diffusion) THEN

      IF (nx > 1 .AND. ny == 1) THEN

        DO i2 = 1,ncomp
          fgradsume = 0.0
          fgradsumw = 0.0
          fjwtchg = 0.0
          DO i = 1,ncomp
            fgradsume = fgradsume - c_d(jx,jy,jz)*chg(i)*fjac_d(i2,i,jx,jy,jz)
            fgradsumw = fgradsumw - a_d(jx,jy,jz)*chg(i)*fjac_d(i2,i,jx,jy,jz)
            fjwtchg = fjwtchg + 0.5*dstar(jx,jy,jz)*chg(i)*fjac_chg(i2,i,jx,jy,jz)
          END DO
          DO i = 1,ncomp
            term1 = -fjac_d(i2,i,jx,jy,jz)*( a_d(jx,jy,jz)+c_d(jx,jy,jz) )
            term2 = 0.5*dstar(jx,jy,jz)*fjac_chg(i2,i,jx,jy,jz)
            term3a = -1.0/(sumsigma_w*sumsigma_w) * fjwtchg
            term3b = -1.0/(sumsigma_e*sumsigma_e) * fjwtchg
            term3 = term3a + term3b
            fgradsum = fgradsumw + fgradsume
            term1_w = -fjac_d(i2,i,jx,jy,jz)*( a_d(jx,jy,jz) )
            
              check1 = term1_w
              check2 = term2*dgradw/sumsigma_w
              check3 = term3a*sigma_w(i)*dgradw
              check4 = fgradsumw*sigma_w(i)/sumsigma_w
              
            fanalyt(i,i2) = df*(  term1  &
              - (term2*dgradw/sumsigma_w + term3a*sigma_w(i)*dgradw + fgradsumw*sigma_w(i)/sumsigma_w)  &
              - (term2*dgrade/sumsigma_e + term3b*sigma_e(i)*dgrade + fgradsume*sigma_e(i)/sumsigma_e)  )

            fanalyt_w = df*(  check1  &
              - (term2*dgradw/sumsigma_w + term3a*sigma_w(i)*dgradw + fgradsumw*sigma_w(i)/sumsigma_w)  )
              
!            fanalyt(i,i2) = df* term1  

          END DO
        END DO

      ELSE IF (nx > 1 .AND. ny > 1) THEN             !  2d PROBLEM
        DO i2 = 1,ncomp
          fgradsume = 0.0
          fgradsumw = 0.0
          fgradsumn = 0.0
          fgradsums = 0.0
          fjwtchg = 0.0
          DO i = 1,ncomp
            fgradsume = fgradsume - c_d(jx,jy,jz)*chg(i)*fjac_d(i2,i,jx,jy,jz)
            fgradsumw = fgradsumw - a_d(jx,jy,jz)*chg(i)*fjac_d(i2,i,jx,jy,jz)
            fgradsumn = fgradsumn - d_d(jx,jy,jz)*chg(i)*fjac_d(i2,i,jx,jy,jz)
            fgradsums = fgradsums - f_d(jx,jy,jz)*chg(i)*fjac_d(i2,i,jx,jy,jz)
            fjwtchg = fjwtchg + 0.5*dstar(jx,jy,jz)*chg(i)*fjac_chg(i2,i,jx,jy,jz)
          END DO
          DO i = 1,ncomp
            term1 = -fjac_d(i2,i,jx,jy,jz)*( a_d(jx,jy,jz)+c_d(jx,jy,jz)+d_d(jx,jy,jz)+f_d(jx,jy,jz) )
            term2 = 0.5*dstar(jx,jy,jz)*fjac_chg(i2,i,jx,jy,jz)
            term3a = -1.0/(sumsigma_w*sumsigma_w) * fjwtchg
            term3b = -1.0/(sumsigma_e*sumsigma_e) * fjwtchg
            term3c = -1.0/(sumsigma_s*sumsigma_s) * fjwtchg
            term3d = -1.0/(sumsigma_n*sumsigma_n) * fjwtchg
            fanalyt(i,i2) = df*(  term1  &
              - (term2*dgradw/sumsigma_w + term3a*sigma_w(i)*dgradw + fgradsumw*sigma_w(i)/sumsigma_w)  &
              - (term2*dgrade/sumsigma_e + term3b*sigma_e(i)*dgrade + fgradsume*sigma_e(i)/sumsigma_e)  &
              - (term2*dgrads/sumsigma_s + term3c*sigma_s(i)*dgrads + fgradsums*sigma_s(i)/sumsigma_s)  &
              - (term2*dgradn/sumsigma_n + term3d*sigma_n(i)*dgradn + fgradsumn*sigma_n(i)/sumsigma_n)  )

            
!            fanalyt(i,i2) = df* term1  
            continue
          END DO
        END DO
      END IF

    END IF

!  Surface charge calculation

    CALL SurfaceCharge(ncomp,nspec,nsurf,nsurf_sec,npot,jx,jy,jz,time)

    sqrt_sion = SQRT(sion(jx,jy,jz))

    DO npt = 1,npot
      ind = (j-1)*(neqn) + npt+ncomp+nexchange+nsurf
      k = kpot(npt)
      HyperbolicSine = SINH(LogPotential(npt,jx,jy,jz) )
!!      fxx(ind) = 0.1174*sqrt_sion*HyperbolicSine - surfcharge(ksurf(ispot(npt)))
      fxx(ind) = 0.1174*sqrt_sion*HyperbolicSine - surfcharge( k )
    END DO
    
    DO i = 1,ncomp
      ind = (j-1)*(neqn) + i
      
      IF (nradmax > 0) THEN
        sumrct = 0.0
        DO k = 1,nkin
          DO np = 1,nreactmin(k)
            IF (rmin(np,k) >= 0.0) then
              sumrct = sumrct + decay_correct(i,k)*mumin(np,k,i)*rmin(np,k)
            ELSE
              sumrct = sumrct + decay_correct(i,k)*mumin_decay(np,k,i,jx,jy,jz)*rmin(np,k)
            END IF
          END DO
        END DO
      ELSE
        sumrct = 0.0d0
        DO k = 1,nkin
          DO np = 1,nreactmin(k)
            sumrct = sumrct + decay_correct(i,k)*mumin(np,k,i)*rmin(np,k)
          END DO
        END DO
      ENDIF
      
      sumkin = 0.0
      DO ir = 1,ikin
        sumkin = sumkin - mukin(ir,i)*raq_tot(ir,jx,jy,jz)
      END DO
      
!  Update the residual, adding reaction terms and exchange terms
      
      fxx(ind) = fxx(ind) + MultiplyCell*(sumrct + satl*xgram(jx,jy,jz)*portemp*rotemp*sumkin)
 
      sumrd = 0.0d0
      sumjackin = 0.0d0
      
      IF (nradmax > 0) THEN
        DO k = 1,nkin
          DO np = 1,nreactmin(k)
            IF (mumin(np,k,i) /= 0.0) THEN
              IF (rmin(np,k) >= 0.0) THEN
                DO i2 = 1,ncomp+nexchange+nsurf
                  sumrd(i2) = sumrd(i2) + decay_correct(i,k)*mumin(np,k,i)*jac_rmin(i2,np,k)
                END DO
              ELSE
                DO i2 = 1,ncomp+nexchange+nsurf
                  sumrd(i2) = sumrd(i2) + decay_correct(i,k)*mumin_decay(np,k,i,jx,jy,jz)*jac_rmin(i2,np,k)
                END DO
              END IF 
            END IF
          END DO
        END DO
      ELSE
        DO k = 1,nkin
          DO np = 1,nreactmin(k)
            IF (mumin(np,k,i) /= 0.0) THEN
              DO i2 = 1,ncomp+nexchange+nsurf
                sumrd(i2) = sumrd(i2) + decay_correct(i,k)*mumin(np,k,i)*jac_rmin(i2,np,k)         
              END DO
            END IF
          END DO
        END DO
      END IF

      DO ir = 1,ikin
        IF (mukin(ir,i) /= 0.0) THEN
          DO i2 = 1,ncomp
            sumjackin(i2) = sumjackin(i2) - mukin(ir,i)*rdkin(ir,i2)
          END DO
        END IF
      END DO


      IF (ierode == 1) THEN 
        DO i2 = 1,ncomp        
          ind2 = i2                
          rxnmin = sumrd(i2)
          rxnaq = satl*xgram(jx,jy,jz)*portemp*rotemp*sumjackin(i2)
          IF (i /= ikh2o) THEN
            aq_accum = H2Oreacted(jx,jy,jz)*satl*xgram(jx,jy,jz)*r*portemp*rotemp*fjac(i2,i,jx,jy,jz)  &
               *(1.0 + Retardation*distrib(i) )
          ELSE
            aq_accum = satl*xgram(jx,jy,jz)*r*portemp*rotemp*fjac(i2,i,jx,jy,jz)  &
               *(1.0 + Retardation*distrib(i) )
          END IF
          source_jac = source*fjac(i2,i,jx,jy,jz)  
          ex_accum = r*fch(i,i2,jx,jy,jz)
          ex_transport = df*(ebu(jx,jy,jz)+bbu(jx,jy,jz))*fch(i,i2,jx,jy,jz)
          alf(ind2,i,2) = MultiplyCell*(rxnmin + rxnaq + aq_accum + ex_accum - source_jac)         &
               + xgram(jx,jy,jz)*df*(e(jx,jy,jz)+b(jx,jy,jz))*fjac(i2,i,jx,jy,jz) + ex_transport  

        END DO   ! end of I2 loop

        DO ix2 = 1,nexchange
          ind2 = ix2+ncomp
          rxnmin = sumrd(ix2+ncomp)
          ex_accum =r*fch(i,ix2+ncomp,jx,jy,jz)
          ex_transport = df*( ebu(jx,jy,jz) + bbu(jx,jy,jz) )*fch(i,ix2+ncomp,jx,jy,jz)
          alf(ind2,i,2) = MultiplyCell*(ex_accum + rxnmin) + ex_transport
        END DO

        DO is2 = 1,nsurf
          ind2 = is2+ncomp+nexchange
          rxnmin = sumrd(is2+ncomp+nexchange)
          surf_accum = r*fch(i,is2+ncomp+nexchange,jx,jy,jz)
          surf_transport = df*      &
             ( ebu(jx,jy,jz) + bbu(jx,jy,jz) )*fch(i,is2+ncomp+nexchange,jx,jy,jz)
          alf(ind2,i,2) = MultiplyCell*(surf_accum + rxnmin) + surf_transport 
        END DO

!!  Dependence of the total aqueous concentration on the potential
        DO npt2 = 1,npot
          ind2 = npt2+ncomp+nexchange+nsurf
!          rxnmin = sumrd(is2+ncomp+nexchange)
          pot_accum = r*fjpotncomp(npt2,i,jx,jy,jz)
          pot_transport = df*      &
             ( ebu(jx,jy,jz) + bbu(jx,jy,jz) )*fjpotncomp(npt2,i,jx,jy,jz)
          alf(ind2,i,2) = MultiplyCell*pot_accum + pot_transport !  + rxnmin

!  NOTE:  Need to add dependence of reaction rate on potentials (if surface complex is in 
!         reaction rate

        END DO
!!  ***************************************************************
      
      ELSE    !! Following is for no burial/erosion

        DO i2 = 1,ncomp        
        IF (immobile_species(i2) == 1) THEN
        coeff_immo = 0.0
        ENDIF
          ind2 = i2                
          rxnmin = sumrd(i2)

          rxnaq = satl*xgram(jx,jy,jz)*portemp*rotemp*sumjackin(i2)

          IF (i /= ikh2o) THEN
            aq_accum = H2Oreacted(jx,jy,jz)*satl*xgram(jx,jy,jz)*r*portemp*rotemp*fjac(i2,i,jx,jy,jz)  &
                *(1.0 + Retardation*distrib(i) ) 

          ELSE
            aq_accum = satl*xgram(jx,jy,jz)*r*portemp*rotemp*fjac(i2,i,jx,jy,jz)  &
                *(1.0 + Retardation*distrib(i) ) 
          END IF

          source_jac = source*fjac(i2,i,jx,jy,jz) 
          ex_accum = r*fch_local(i,i2) 
          alf(ind2,i,2) = MultiplyCell*(rxnmin + rxnaq + aq_accum + ex_accum - source_jac)   &
               + xgram(jx,jy,jz)*df*(e(jx,jy,jz)+b(jx,jy,jz))*fjac(i2,i,jx,jy,jz)!*coeff_immo 

        END DO   ! end of I2 loop

        DO ix2 = 1,nexchange
          ind2 = ix2+ncomp
          rxnmin = sumrd(ix2+ncomp)
          ex_accum = r*fch_local(i,ix2+ncomp)
          alf(ind2,i,2) = MultiplyCell*(ex_accum + rxnmin)
        END DO

        DO is2 = 1,nsurf
          ind2 = is2+ncomp+nexchange
          rxnmin = sumrd(is2+ncomp+nexchange)
          surf_accum = r*fch_local(i,is2+ncomp+nexchange)
          alf(ind2,i,2) = MultiplyCell*(surf_accum + rxnmin)
        END DO

!!  Dependence of the total aqueous concentration on the potential
        DO npt2 = 1,npot
          ind2 = npt2+ncomp+nexchange+nsurf
!          rxnmin = sumrd(is2+ncomp+nexchange)
          pot_accum = fjpotncomp(npt2,i,jx,jy,jz)/delt
          alf(ind2,i,2) = MultiplyCell*pot_accum  !  + rxnmin
!!!          alf(ind2,i,2) = 0.0d0
!  NOTE:  Need to add dependence of reaction rate on potentials (if surface complex is in 
!         reaction rate
        END DO 
!!  ***************************************************************

      END IF

      IF (isaturate == 1) THEN
        DO i2 = 1,ncomp     
          ind2 = i2      
          alf(ind2,i,2) = alf(ind2,i,2) + MultiplyCell*satgas*portemp*r*fgas(i2,i,jx,jy,jz)    &
             + df*(bg(jx,jy,jz)+eg(jx,jy,jz))*fgas(i2,i,jx,jy,jz)
        END DO   ! end of I2 loop
      END IF

      IF (species_diffusion) THEN
        DO i2 = 1,ncomp
          ind2 = i2
          alf(ind2,i,2) = alf(ind2,i,2) + fanalyt(i,i2)
        END DO
      END IF

      IF (ihindmarsh == 1 .AND. nxyz == nx .AND. nx /= 1) THEN
        DO i2 = 1,neqn
          aah(i,i2,jx) = alf(i2,i,2)
        END DO

      END IF
      
    END DO     ! end of I loop

    do ix = 1,nexchange
      sum = 0.0
      do nex = 1,nexch_sec
        sum = sum + muexc(nex,ix+ncomp)*spex10(nex+nexchange,jx,jy,jz)
      end do
      totex(ix) = sum
    end do
    
    DO ix = 1,nexchange

      ind = (j-1)*(neqn) + ix+ncomp

! Update the residual

      IF (iexc == 1 .OR. iexc == 3) THEN
!!        fxx(ind) = (totex(ix) - exchangesites(ix,jx,jy,jz))
        fxx(ind) = (sumactivity(ix) - 1.0)
      ELSE
        fxx(ind) = sumactivity(ix) - 1.0
      END IF

      DO i2 = 1,ncomp+nexchange
        ind2 = i2
        alf(ind2,ix+ncomp,2) = fexch(ix+ncomp,i2)
!!        alf(ind2,ix+ncomp,2) = fexch(ix+ncomp,i2)
      END DO

      IF (nxyz == nx .AND. ihindmarsh == 1 .AND. nxyz /= 1) THEN
        DO i2 = 1,ncomp+nexchange
          ind2 = i2
          aah(ix+ncomp,i2,jx) = alf(ind2,ix+ncomp,2)
        END DO
      END IF

    END DO
    
    DO is = 1,nsurf

      ind = (j-1)*(neqn) + is+ncomp+nexchange
     
      IF (ierode == 1) THEN
      
        DO i2 = 1,ncomp
          ind2 = i2
          surf_accum = r*fsurf(is,i2,jx,jy,jz)
          surf_transport = df*(bbu(jx,jy,jz)+ebu(jx,jy,jz))*fsurf(is,i2,jx,jy,jz)
          alf(ind2,is+ncomp+nexchange,2) = MultiplyCell*(surf_accum) + surf_transport
        END DO

        DO is2 = 1,nsurf
          ind2 = is2+ncomp+nexchange
          surf_accum = r* fsurf(is,is2+ncomp,jx,jy,jz)
          surf_transport = df*(bbu(jx,jy,jz)+ebu(jx,jy,jz))*fsurf(is,is2+ncomp,jx,jy,jz)
          alf(ind2,is+ncomp+nexchange,2) = MultiplyCell*(surf_accum) + surf_transport
        END DO

!!  Dependence of the total surface complex concentration on the potential
        DO npt2 = 1,npot
          ind2 = npt2+ncomp+nexchange+nsurf
          pot_accum = r*fjpotnsurf(npt2,is,jx,jy,jz)
          pot_transport = df*(bbu(jx,jy,jz)+ebu(jx,jy,jz))*     &
               fjpotnsurf(npt2,is,jx,jy,jz)
          alf(ind2,is+ncomp+nexchange,2) = MultiplyCell*(pot_accum) + pot_transport
        END DO
!!  ***********************************************************************
      
      ELSE

        DO i2 = 1,ncomp
          ind2 = i2
          surf_accum = r*fsurf_local(is,i2)
          alf(ind2,is+ncomp+nexchange,2) = MultiplyCell*(surf_accum)
        END DO

        DO is2 = 1,nsurf
          ind2 = is2+ncomp+nexchange
          surf_accum = r*fsurf_local(is,is2+ncomp)
          alf(ind2,is+ncomp+nexchange,2) = MultiplyCell*(surf_accum)
        END DO

!!  Dependence of the total surface complex concentration on the potential
        DO npt2 = 1,npot
          ind2 = npt2+ncomp+nexchange+nsurf
          pot_accum = fjpotnsurf(npt2,is,jx,jy,jz)/delt
          alf(ind2,is+ncomp+nexchange,2) = MultiplyCell*(pot_accum) 
        END DO
!!  ***********************************************************************

      END IF
      
      IF (nxyz == nx .AND. ihindmarsh == 1 .AND. nxyz /= 1) THEN
        DO i2 = 1,neqn
          aah(is+ncomp+nexchange,i2,jx) = alf(i2,is+ncomp+nexchange,2)
        END DO
      END IF
      
    END DO           !  End of IS loop

!  Jacobian entries:  Dependence of residual of surface potential equation on 
!     other unknowns (ncomp,nsurf,npot)
    

    DO npt = 1,npot

      ncol = npt + ncomp + nexchange + nsurf
      k = kpot(npt)       !!  One to one correspondence between potential and mineral surface (k)

      IF (volinByGrid(k,jx,jy,jz) == 0.0d0 .AND. volfx(k,jx,jy,jz) < voltemp(k,jinit(jx,jy,jz)) ) THEN
        correct = wtmin(k)*specificByGrid(k,jx,jy,jz)*voltemp(k,jinit(jx,jy,jz))/volmol(k)   !!  m^2 mineral/m^3 BV
      ELSE
        volMinimum = volfx(k,jx,jy,jz)
        if (volMinimum < 1.0D-15) then
          volMinimum = 1.0D-15
        end if
        correct = wtmin(k)*specificByGrid(k,jx,jy,jz)*volMinimum/volmol(k)   !!  m^2 mineral/m^3 BV
      END IF
 
!!  Dependence of the potential on primary species concentrations
      DO i2 = 1,ncomp
        ind2 = i2
        sum = 0.0
        DO ns = 1,nsurf_sec
          IF (ksurf(islink(ns)) == kpot(npt)) THEN
            sum = sum - musurf(ns,i2)*zsurf(ns+nsurf)*spsurf10(ns+nsurf,jx,jy,jz)*faraday/correct
          END IF
        END DO
        alf(ind2,ncol,2) = sum
      END DO

!!  Dependence of the potential on surface complex concentrations
      DO is2 = 1,nsurf
        ind2 = is2+ncomp+nexchange
        sum = 0.0d0
        DO ns = 1,nsurf_sec
!!        Add on secondary surface complex if it is associated with the npt potential
          IF (ksurf(islink(ns)) == kpot(npt)) THEN
            sum = sum - musurf(ns,is2+ncomp)*zsurf(ns+nsurf)*spsurf10(ns+nsurf,jx,jy,jz)*faraday/correct
          END IF
        END DO
        alf(ind2,ncol,2) = sum
!!      Add on primary surface complex if it is associated with the npt potential
        IF (ksurf(is2) == kpot(npt)) THEN
          alf(ind2,ncol,2) =  alf(ind2,ncol,2) - zsurf(is2)*spsurf10(is2,jx,jy,jz)*faraday/correct
        END IF

      END DO        !!  End of is2 loop

!!    Dependence of the potential equation on the potential (through the surface charge)
      DO npt2 = 1,npot
        sum = 0.0
        nrow = npt2 + ncomp + nexchange + nsurf
          DO ns = 1,nsurf_sec
            delta_z = zsurf(ns+nsurf) - zsurf(islink(ns))
            IF (ksurf(islink(ns)) == kpot(npt) .AND. kpot(npt2) == kpot(npt)) THEN
              sum = sum - zsurf(ns+nsurf)*spsurf10(ns+nsurf,jx,jy,jz)*delta_z*2.0
            END IF
          END DO

        IF (nrow == ncol) THEN
           alf(nrow,ncol,2) = 0.1174d0*sqrt_sion*COSH(LogPotential(npt,jx,jy,jz)) - sum*faraday/correct
        ELSE
!!!           alf(nrow,ncol,2) =  -sum*faraday/correct
          alf(nrow,ncol,2) = 0.0d0
        END IF

    END DO          !!  End of npt (potential) loop


!!   *********************************************************************************

      IF (nxyz == nx .AND. ihindmarsh == 1 .AND. nxyz /= 1) THEN
        DO i2 = 1,neqn
          aah(ncol,i2,jx) = alf(i2,ncol,2)
        END DO
      END IF

    END DO
    
    GO TO 850
    
!  Case where nxyz = 1
    
    850     IF (nxyz == 1) THEN
      DO i = 1,neqn
        DO i2 = 1,neqn
          ind2 = i2
          aaa(i,i2) = alf(ind2,i,2)
        END DO
      END DO
      GO TO 900
    END IF
    
!  For fixed concentration boundary conditions
    

    800     IF (activecell(jx,jy,jz) == 0) THEN
      DO ix = 1,nexchange
        ind = (j-1)*(neqn) + ix+ncomp
        fxx(ind) = spex10(ix,jx,jy,jz) - spcondex10(ix,jinit(jx,jy,jz))*   &
               xgram(jx,jy,jz)*por(jx,jy,jz)*ro(jx,jy,jz)*satliq(jx,jy,jz)
      END DO
      DO i = 1,ncomp
        ind2 = i
        alf(ind2,i,2) = sp10(i,jx,jy,jz)
        IF (nxyz == nx .AND. ihindmarsh == 1 .AND. nxyz /= 1) THEN
          aah(i,i,jx) = alf(ind2,i,2)
        END IF
      END DO
      DO ix = 1,nexchange
        ind2 =  ix+ncomp
        alf(ind2,ix+ncomp,2) = spex10(ix,jx,jy,jz)
        IF (nxyz == nx .AND. ihindmarsh == 1 .AND. nxyz /= 1) THEN
          aah(ix+ncomp,ix+ncomp,jx) = alf(ind2,ix+ncomp,2)
        END IF
      END DO
      DO is = 1,nsurf
        ind2 = is+ncomp+nexchange
        alf(ind2,is+ncomp+nexchange,2) = spsurf10(is,jx,jy,jz)
        IF (nxyz == nx .AND. ihindmarsh == 1 .AND. nxyz /= 1) THEN
          aah(is+ncomp+nexchange,is+ncomp+nexchange,jx) = alf(ind2,is+ncomp+nexchange,2)
        END IF
      END DO
    END IF

    IF (ikh2o /= 0) THEN

      alf(:,ikh2o,:) = 0.0d0

      ind = (j-1)*(neqn) + ikh2o
      fxx(ind) = sp10(ikh2o,jx,jy,jz) - 55.40d0

      ind2 = ikh2o
      alf(ind2,ikh2o,2) = sp10(ikh2o,jx,jy,jz)
      IF (nxyz == nx .AND. ihindmarsh == 1 .AND. nxyz /= 1) THEN

        aah(ikh2o,:,jx) = 0.0d0
        bbh(ikh2o,:,jx) = 0.0d0
        cch(ikh2o,:,jx) = 0.0d0
        aah(:,ikh2o,jx) = 0.0d0
        bbh(:,ikh2o,jx) = 0.0d0
        cch(:,ikh2o,jx) = 0.0d0
        aah(ikh2o,ikh2o,jx) = alf(ikh2o,ikh2o,2)

      END IF

    END IF

    IF (ihindmarsh == 0 .OR. nxyz /= nx) THEN

     IF ( petscon) THEN

! *********************** insert filling matrix for PETSc case ***************************

       if(ny.eq.1) then  !   1-D case for PETSc

         if(jx.eq.1) then
             do irow =1,neqn
             do jcol = 1,neqn
               blockm(irow,jcol)=alf(jcol,irow,2)
               blockr(irow,jcol)=alf(jcol,irow,3)
             end do
             end do
             call MatSetValuesBlocked(amatpetsc,1,ipetsc,1,ipetsc,blockm,INSERT_VALUES,ierr)
             call MatSetValuesBlocked(amatpetsc,1,ipetsc,1,ipetsc+1,blockr,INSERT_VALUES,ierr)
         elseif(jx.eq.nx) then
             do irow =1,neqn
             do jcol = 1,neqn
              blockl(irow,jcol)=alf(jcol,irow,1)
              blockm(irow,jcol)=alf(jcol,irow,2)    
             end do
             end do
             call MatSetValuesBlocked(amatpetsc,1,ipetsc,1,ipetsc-1,blockl,INSERT_VALUES,ierr)
             call MatSetValuesBlocked(amatpetsc,1,ipetsc,1,ipetsc,blockm,INSERT_VALUES,ierr)
         else
             do irow =1,neqn
             do jcol = 1,neqn
              blockl(irow,jcol)=alf(jcol,irow,1)
              blockm(irow,jcol)=alf(jcol,irow,2)
              blockr(irow,jcol)=alf(jcol,irow,3)
             end do
             end do
             call MatSetValuesBlocked(amatpetsc,1,ipetsc,1,ipetsc-1,blockl,INSERT_VALUES,ierr)
             call MatSetValuesBlocked(amatpetsc,1,ipetsc,1,ipetsc,blockm,INSERT_VALUES,ierr)
             call MatSetValuesBlocked(amatpetsc,1,ipetsc,1,ipetsc+1,blockr,INSERT_VALUES,ierr)
          endif

          IF (Switcheroo) THEN
            do ind  = 1,neqn
              do ind2 = 1,neqn
                cch(ind,ind2,jx) = alf(ind2,ind,1)
                aah(ind,ind2,jx) = alf(ind2,ind,2)
                bbh(ind,ind2,jx) = alf(ind2,ind,3)
              end do
            end do
          END IF

       else   !  2-D case for PETSc

        if (jx == 1) then

         if (jy ==1) then
            do irow =1,neqn
            do jcol = 1,neqn
             blockm(irow,jcol)=alf(jcol,irow,2)
             blockr(irow,jcol)=alf(jcol,irow,3)
             blockfr(irow,jcol)=alf(jcol,irow,5)
            end do
            end do
            call MatSetValuesBlocked(amatpetsc,1,ipetsc,1,ipetsc,blockm,INSERT_VALUES,ierr)
            call MatSetValuesBlocked(amatpetsc,1,ipetsc,1,ipetsc+1,blockr,INSERT_VALUES,ierr)
            call MatSetValuesBlocked(amatpetsc,1,ipetsc,1,ipetsc+nx,blockfr,INSERT_VALUES,ierr)
         elseif (jy == ny) then
            do irow =1,neqn
            do jcol = 1,neqn
             blockfl(irow,jcol)=alf(jcol,irow,4)
             blockm(irow,jcol)=alf(jcol,irow,2)
             blockr(irow,jcol)=alf(jcol,irow,3)
            end do
            end do
            call MatSetValuesBlocked(amatpetsc,1,ipetsc,1,ipetsc-nx,blockfl,INSERT_VALUES,ierr)
            call MatSetValuesBlocked(amatpetsc,1,ipetsc,1,ipetsc,blockm,INSERT_VALUES,ierr)
            call MatSetValuesBlocked(amatpetsc,1,ipetsc,1,ipetsc+1,blockr,INSERT_VALUES,ierr)
         else
            do irow =1,neqn
            do jcol = 1,neqn
             blockfl(irow,jcol)=alf(jcol,irow,4)
             blockm(irow,jcol)=alf(jcol,irow,2)
             blockr(irow,jcol)=alf(jcol,irow,3)
             blockfr(irow,jcol)=alf(jcol,irow,5)
            end do
            end do
            call MatSetValuesBlocked(amatpetsc,1,ipetsc,1,ipetsc-nx,blockfl,INSERT_VALUES,ierr)
            call MatSetValuesBlocked(amatpetsc,1,ipetsc,1,ipetsc,blockm,INSERT_VALUES,ierr)
            call MatSetValuesBlocked(amatpetsc,1,ipetsc,1,ipetsc+1,blockr,INSERT_VALUES,ierr)
            call MatSetValuesBlocked(amatpetsc,1,ipetsc,1,ipetsc+nx,blockfr,INSERT_VALUES,ierr)
         endif

        elseif (jx == nx ) then                             

         if (jy ==1) then
            do irow =1,neqn
            do jcol = 1,neqn
             blockl(irow,jcol)=alf(jcol,irow,1)
             blockm(irow,jcol)=alf(jcol,irow,2)
             blockfr(irow,jcol)=alf(jcol,irow,5)
            end do
            end do
            call MatSetValuesBlocked(amatpetsc,1,ipetsc,1,ipetsc-1,blockl,INSERT_VALUES,ierr)
            call MatSetValuesBlocked(amatpetsc,1,ipetsc,1,ipetsc,blockm,INSERT_VALUES,ierr)
            call MatSetValuesBlocked(amatpetsc,1,ipetsc,1,ipetsc+nx,blockfr,INSERT_VALUES,ierr)
         elseif (jy == ny) then
            do irow =1,neqn
            do jcol = 1,neqn
             blockfl(irow,jcol)=alf(jcol,irow,4)
             blockl(irow,jcol)=alf(jcol,irow,1)
             blockm(irow,jcol)=alf(jcol,irow,2)
            end do
            end do
            call MatSetValuesBlocked(amatpetsc,1,ipetsc,1,ipetsc-nx,blockfl,INSERT_VALUES,ierr)
            call MatSetValuesBlocked(amatpetsc,1,ipetsc,1,ipetsc-1,blockl,INSERT_VALUES,ierr)
            call MatSetValuesBlocked(amatpetsc,1,ipetsc,1,ipetsc,blockm,INSERT_VALUES,ierr)
         else
            do irow =1,neqn
            do jcol = 1,neqn
             blockfl(irow,jcol)=alf(jcol,irow,4)
             blockl(irow,jcol)=alf(jcol,irow,1)
             blockm(irow,jcol)=alf(jcol,irow,2)
             blockfr(irow,jcol)=alf(jcol,irow,5)
            end do
            end do
            call MatSetValuesBlocked(amatpetsc,1,ipetsc,1,ipetsc-nx,blockfl,INSERT_VALUES,ierr)
            call MatSetValuesBlocked(amatpetsc,1,ipetsc,1,ipetsc-1,blockl,INSERT_VALUES,ierr)
            call MatSetValuesBlocked(amatpetsc,1,ipetsc,1,ipetsc,blockm,INSERT_VALUES,ierr)
            call MatSetValuesBlocked(amatpetsc,1,ipetsc,1,ipetsc+nx,blockfr,INSERT_VALUES,ierr)
         endif

        else
         
         if (jy ==1) then
            do irow =1,neqn
            do jcol = 1,neqn
             blockl(irow,jcol)=alf(jcol,irow,1)
             blockm(irow,jcol)=alf(jcol,irow,2)
             blockr(irow,jcol)=alf(jcol,irow,3)
             blockfr(irow,jcol)=alf(jcol,irow,5)
            end do
            end do
            call MatSetValuesBlocked(amatpetsc,1,ipetsc,1,ipetsc-1,blockl,INSERT_VALUES,ierr)
            call MatSetValuesBlocked(amatpetsc,1,ipetsc,1,ipetsc,blockm,INSERT_VALUES,ierr)
            call MatSetValuesBlocked(amatpetsc,1,ipetsc,1,ipetsc+1,blockr,INSERT_VALUES,ierr)
            call MatSetValuesBlocked(amatpetsc,1,ipetsc,1,ipetsc+nx,blockfr,INSERT_VALUES,ierr)
         elseif (jy == ny) then
            do irow =1,neqn
            do jcol = 1,neqn
             blockfl(irow,jcol)=alf(jcol,irow,4)
             blockl(irow,jcol)=alf(jcol,irow,1)
             blockm(irow,jcol)=alf(jcol,irow,2)
             blockr(irow,jcol)=alf(jcol,irow,3)
            end do
            end do
            call MatSetValuesBlocked(amatpetsc,1,ipetsc,1,ipetsc-nx,blockfl,INSERT_VALUES,ierr)
            call MatSetValuesBlocked(amatpetsc,1,ipetsc,1,ipetsc-1,blockl,INSERT_VALUES,ierr)
            call MatSetValuesBlocked(amatpetsc,1,ipetsc,1,ipetsc,blockm,INSERT_VALUES,ierr)
            call MatSetValuesBlocked(amatpetsc,1,ipetsc,1,ipetsc+1,blockr,INSERT_VALUES,ierr)
         else
            do irow =1,neqn
            do jcol = 1,neqn
             blockfl(irow,jcol)=alf(jcol,irow,4)
             blockl(irow,jcol)=alf(jcol,irow,1)
             blockm(irow,jcol)=alf(jcol,irow,2)
             blockr(irow,jcol)=alf(jcol,irow,3)
             blockfr(irow,jcol)=alf(jcol,irow,5)
            end do
            end do
            call MatSetValuesBlocked(amatpetsc,1,ipetsc,1,ipetsc-nx,blockfl,INSERT_VALUES,ierr)
            call MatSetValuesBlocked(amatpetsc,1,ipetsc,1,ipetsc-1,blockl,INSERT_VALUES,ierr)
            call MatSetValuesBlocked(amatpetsc,1,ipetsc,1,ipetsc,blockm,INSERT_VALUES,ierr)
            call MatSetValuesBlocked(amatpetsc,1,ipetsc,1,ipetsc+1,blockr,INSERT_VALUES,ierr)
            call MatSetValuesBlocked(amatpetsc,1,ipetsc,1,ipetsc+nx,blockfr,INSERT_VALUES,ierr)
         endif

        endif ! test on jx in 2-D fill

      endif ! 1-D or 2-D cases for PETSc

! ********************* end insert of PETSc matrix fill **********************8

     endif ! Watsolv or PETSc dichotomy/toggle   

    END IF

!***************************insert update of petsc index **********************    
    ipetsc = ipetsc + 1
!***************************end update of petsc index ************************
  END DO  ! End of JX loop
END DO    ! End of JY loop

900 CONTINUE

1077 FORMAT(12(1X,1PE10.3))

! ***************************insert construction for PETSc matrix ***********
 if (petscon) then
          call MatAssemblyBegin(amatpetsc,MAT_FINAL_ASSEMBLY,ierr)
          call MatAssemblyEnd(amatpetsc,MAT_FINAL_ASSEMBLY,ierr)
 endif
! ******************end construction for PETSc matrix ********************




RETURN
END SUBROUTINE AssembleGlobal
