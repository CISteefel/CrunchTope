
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
        amatpetsc,nBoundaryConditionZone,igamma)
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
USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: INT64, REAL64

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
INTEGER(I4B), INTENT(IN)                      :: nBoundaryConditionZone
INTEGER(I4B), INTENT(IN)                      :: igamma


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

INTEGER(I4B)                                  :: npz

REAL(DP)                                      :: sqrt_sion
REAL(DP)                                      :: sum
REAL(DP)                                      :: rotemp
REAL(DP)                                      :: portemp
REAL(DP)                                      :: satl
REAL(DP)                                      :: satlold
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

INTEGER(I4B)                                   :: pos_IonS
INTEGER(I4B)                                   :: pos_GammaWater
INTEGER(I4B)                                   :: pos_der
INTEGER(I4B)                                   :: nk
INTEGER(I4B)                                   :: icomp
INTEGER(I4B)                                   :: ksp
REAL(DP)                                       :: IonicStrength
REAL(DP)                                       :: ChargeSum
REAL(DP)                                       :: TotalMoles
REAL(DP)                                       :: sumder
REAL(DP)                                       :: gammaH2O
REAL(DP)                                       :: CheckH2O
INTEGER(I4B)                                   :: icomp2

REAL(DP), DIMENSION(ncomp + nsurf + nexchange + npot + 1 + 1,ncomp + nsurf + nexchange + npot + 1 + 1)  :: der_residuals

CHARACTER (LEN=3)                              :: ulabprint

!! Time normalized used if time series only defined for 1 representative year:
REAL(DP)                                       :: time_norm
REAL(DP), DIMENSION(:), ALLOCATABLE            :: temp_dum

!********************* PETSc declarations ********************************
!!! PetscFortranAddr                                                    user(6)
Mat                                                                 amatpetsc
!******************end PETSc declarations *********************************

INTEGER(INT64) :: count_start, count_end, count_rate, phase_start, phase_end, call_start, call_end
REAL(REAL64) :: dt_local, dt_insert
INTEGER(INT64) :: fx_start, fx_end, chemistry_start, chemistry_end
INTEGER(INT64) :: reaction_start, reaction_end, neighbor_start, neighbor_end, primary_start, primary_end

REAL(REAL64), SAVE :: assemble_wall_total = 0.0_REAL64
REAL(REAL64)       :: elapsed_time

CALL SYSTEM_CLOCK(COUNT_RATE=count_rate)
CALL SYSTEM_CLOCK(COUNT=count_start)
call_start = count_start

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

!!! Adding two additional unknowns, one for ionic strength, one for water activity
neqn = ncomp + nsurf + nexchange + npot + 1 + 1
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
    CALL SurfaceComplex(ncomp,nsurf,nsurf_sec,nx,ny,nz,igamma)
  END IF
END IF

IF (npot > 0) THEN
  CALL jacpotential(ncomp,nsurf,nsurf_sec,npot,nx,ny,nz)
END IF

TempFlux = 0.0d0

!!!  Do the boundaries first

!!!  ***** Nernst-Planck ******

jz = 1
DO jy = 1,ny
  DO jx = 1,nx
    CALL SYSTEM_CLOCK(COUNT=phase_start)
    j = (jz-1)*nx*ny + (jy-1)*nx + jx
   
    jpoint = jinit(jx,jy,jz)
    portemp = por(jx,jy,jz)
    rotemp = ro(jx,jy,jz)

    Retardation = 0.001d0*SolidDensity(jinit(jx,jy,jz))*(1.0-por(jx,jy,jz))/por(jx,jy,jz)

    CALL AqueousToBulkConvert(jx,jy,jz,AqueousToBulk)

    satl = satliq(jx,jy,jz)
    satgas = 1.0d0 - satl
    satlold = satliqold(jx,jy,jz)
    
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
    IF (ihindmarsh == 1 .AND. nxyz == nx .AND. nxyz /= 1) THEN
      aah(:,:,jx) = 0.0d0
      bbh(:,:,jx) = 0.0d0
      cch(:,:,jx) = 0.0d0
    END IF
    
    IF (ierode /= 1) THEN
      IF (nexchange > 0) THEN
        CALL exchange(ncomp,nexchange,nexch_sec,jx,jy,jz)
      END IF
    END IF
    
    CALL SYSTEM_CLOCK(COUNT=fx_start)
    CALL FxTopeGlobal(nx,ny,nz,ncomp,nexchange,nexch_sec,nsurf,nsurf_sec,nrct,nspec,  &
               ngas,neqn,delt,jx,jy,jz,nBoundaryConditionZone)
    CALL SYSTEM_CLOCK(COUNT=fx_end)
    fx_wall_total = fx_wall_total + REAL(fx_end-fx_start,REAL64)/REAL(count_rate,REAL64)
      
        
    ! Time reaction and kinetic routines together (not FxTopeGlobal).
    CALL SYSTEM_CLOCK(COUNT=reaction_start)
    IF (nrct > 0) THEN
      CALL reaction(ncomp,nkin,nrct,nspec,nexchange,nsurf,ndecay,jx,jy,jz,delt,time)
      round = 1
      CALL jacmin(ncomp,nspec,nexchange,nsurf,nkin,nrct,jx,jy,jz,time,round)
    END IF
    
    IF (ikin > 0) THEN
      CALL reactkin(ncomp,nspec,nrct,ikin,jx,jy,jz,AqueousToBulk,time)
      CALL jacrkin(ncomp,nspec,nrct,ikin,jx,jy,jz,AqueousToBulk)
    END IF

    CALL SYSTEM_CLOCK(COUNT=reaction_end)
    reaction_calls_wall_total = reaction_calls_wall_total + &
      REAL(reaction_end-reaction_start,REAL64)/REAL(count_rate,REAL64)

    IF (activecell(jx,jy,jz) == 0) GO TO 800
    
    ! Start before both dimension-dependent branches; stop after label 200.
    CALL SYSTEM_CLOCK(COUNT=neighbor_start)
    IF (nx == 1) GO TO 100
    
    
    IF (nxyz == nx .AND. ihindmarsh == 1 .AND. nxyz /= 1) THEN  ! Use Hindmarsh routine
  
      
    ELSE     !  Case where NY .ne. 1 (2D problem)
      
      IF (jx /= 1) THEN
        jdum=jx-1

        
        DO i = 1,ncomp
          ind = (j-1)*(neqn) + i
          
          DO i2 = 1,ncomp
            ind2 = i2
            alf(ind2,i,1) = xgram(jdum,jy,jz)*df*a(jx,jy,jz)*fjac(i2,i,jdum,jy,jz)
          END DO




          

        END DO
        

        
      END IF
      
      IF (jx /= nx) THEN
        jdum=jx+1



        DO i = 1,ncomp
          ind = (j-1)*(neqn) + i

          DO i2 = 1,ncomp
            ind2 = i2
            alf(ind2,i,3) = xgram(jdum,jy,jz)*df*c(jx,jy,jz)*fjac(i2,i,jdum,jy,jz)
          END DO




          
          
        END DO
        
        
      END IF
    END IF
    
    100     CONTINUE
    IF (ny == 1) GO TO 200
    
    IF (jy /= ny) THEN
      jdum=jy+1


      DO i = 1,ncomp
        ind = (j-1)*(neqn) + i

        DO i2 = 1,ncomp
          ind2 = i2
          alf(ind2,i,5) = xgram(jx,jdum,jz)*df*d(jx,jy,jz)*fjac(i2,i,jx,jdum,jz)
        END DO
        


        
        
      END DO
      
      
    END IF

    IF (jy /= 1) THEN
      jdum=jy-1


      DO i = 1,ncomp
        ind = (j-1)*neqn + i

        DO i2 = 1,ncomp
          ind2 = i2
          alf(ind2,i,4) = xgram(jx,jdum,jz)*df*f(jx,jy,jz)*fjac(i2,i,jx,jdum,jz)
        END DO



        
      END DO
      
      
    END IF
    
 200 CONTINUE
    CALL SYSTEM_CLOCK(COUNT=neighbor_end)
    neighbor_alf_wall_total = neighbor_alf_wall_total + &
      REAL(neighbor_end-neighbor_start,REAL64)/REAL(count_rate,REAL64)
    
    source = 0.0d0

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
      
    END IF
    

!  Surface charge calculation

    !!! CALL SurfaceCharge(ncomp,nspec,nsurf,nsurf_sec,npot,jx,jy,jz,time)

    sqrt_sion = SQRT(sion(jx,jy,jz))
    CALL SYSTEM_CLOCK(COUNT=primary_start)

    
    DO i = 1,ncomp
      ind = (j-1)*(neqn) + i
      
      IF (nradmax > 0) THEN

      ELSE
        sumrct = 0.0d0
        DO k = 1,nkin
          DO np = 1,nreactmin(k)
            sumrct = sumrct + mumin(np,k,i)*rmin(np,k)
!!!            sumrct = sumrct + decay_correct(i,k)*mumin(np,k,i)*rmin(np,k)
          END DO
        END DO
      ENDIF
      
      sumkin = 0.0
      DO ir = 1,ikin
        sumkin = sumkin - mukin(ir,i)*raq_tot(ir,jx,jy,jz)
      END DO
      
!  Update the residual, adding reaction terms and exchange terms
      
      fxx(ind) = fxx(ind) + MultiplyCell*(sumrct + 0.5*(satl+satlold)*xgram(jx,jy,jz)*portemp*rotemp*sumkin)
      
      sumrd = 0.0d0
      sumjackin = 0.0d0
      
      IF (nradmax > 0) THEN

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



      IF (ierode == 1) THEN 
        

!!  ***************************************************************
      
      ELSE    !! Following is for no burial/erosion

        DO i2 = 1,ncomp        
          ind2 = i2                
          rxnmin = sumrd(i2)

          rxnaq = satl*xgram(jx,jy,jz)*portemp*rotemp*sumjackin(i2)

          IF (i /= ikh2o) THEN
            aq_accum = H2Oreacted(jx,jy,jz)*satl*xgram(jx,jy,jz) * portemp * rotemp * fjac(i2,i,jx,jy,jz) *  &
                (1.0 + Retardation*distrib(i) )  /delt
          ELSE
            aq_accum = satl * xgram(jx,jy,jz) * portemp * rotemp * fjac(i2,i,jx,jy,jz)  /delt
          END IF

          source_jac = source*fjac(i2,i,jx,jy,jz) 
          ex_accum = r*fch_local(i,i2) 
          alf(ind2,i,2) = MultiplyCell*(rxnmin + rxnaq + aq_accum + ex_accum - source_jac)   &
               + xgram(jx,jy,jz)*df*(e(jx,jy,jz)+b(jx,jy,jz))*fjac(i2,i,jx,jy,jz) 
          continue

        END DO   ! end of I2 loop


!!  ***************************************************************

      END IF

      IF (ihindmarsh == 1 .AND. nxyz == nx .AND. nx /= 1) THEN
        DO i2 = 1,neqn
          aah(i,i2,jx) = alf(i2,i,2)
        END DO
      END IF
      
    END DO     ! end of I primary species loop
    CALL SYSTEM_CLOCK(COUNT=primary_end)
    primary_alf_wall_total = primary_alf_wall_total + &
      REAL(primary_end-primary_start,REAL64)/REAL(count_rate,REAL64)

!!! Start of IX exchange loop
    

    

    
!  Jacobian entries:  Dependence of residual of surface potential equation on 
!     other unknowns (ncomp,nsurf,npot)
    
!!! Start of NPT electrical potential loop
    

    
!!! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
    CALL SYSTEM_CLOCK(COUNT=chemistry_start)
!!! --> Start of new functionality section (Ionic Strength and Activity of Water)
    
      pos_IonS       = ncomp + nexchange + nsurf + npot + 1
      pos_gammawater = ncomp + nexchange + nsurf + npot + 1 + 1
      
!!! --> Add in residual for Ionic Strength
       
      der_residuals(pos_IonS,:) = 0D0
      
      ChargeSum = 0.0d0
      DO ik = 1,ncomp+nspec
        ChargeSum  = ChargeSum + sp10(ik,jx,jy,jz)*chg(ik)*chg(ik)
      END DO
      
!!!      sion(jx,jy,jz) = 0.5d0 * ChargeSum
      IonicStrength = sion(jx,jy,jz)
      
!!! Residual for ionic strength
      ind = (j-1)*(neqn) + pos_IonS
      fxx(ind) =  LOG( 0.5D0 * ChargeSum / IonicStrength)
 
!!!  NOTE: First subscript points to residual, second to variable differentiated with respect to...
      
      !!! Derivative of ionic strength with respect to "icomp" primary species
      DO icomp = 1, ncomp
        sumder = 0D0
        ulabPrint = ulab(icomp)
        IF (ulabPrint(1:3) /= 'H2O' .and. ulabPrint(1:3) /= 'HHO') THEN
          sumder = sp10(icomp,jx,jy,jz)*chg(icomp)*chg(icomp)
          DO ksp = 1, nspec
            nk = ksp + ncomp
            sumder = sumder + chg(nk)*chg(nk) * deriv_conc(nk,icomp,jx,jy,jz)
          END DO
        END IF

        der_residuals(pos_IonS,icomp) = sumder / ChargeSum
        continue
        
      END DO  
      
!!! Derivative of ionic strength with respect to ionic strength (diagonal)
      pos_der = ncomp + nexchange + nsurf + npot + 1     !! Diagonal element 
      sumder = 0D0
      DO ksp = 1, nspec
        nk = ksp + ncomp
        sumder = sumder + chg(nk)*chg(nk) * deriv_conc(nk,pos_der,jx,jy,jz)
      END DO
      der_residuals(pos_IonS,pos_der) = sumder / ChargeSum - 1.0D0


!!! Derivative of ionic strength with respect to gamma water
      pos_der = ncomp + nexchange + nsurf + npot + 1 + 1
      sumder = 0D0
      DO ksp = 1, nspec
        nk = ksp + ncomp
        sumder = sumder + chg(nk)*chg(nk) * deriv_conc(nk,pos_der,jx,jy,jz)
      END DO
      der_residuals(pos_IonS,pos_der) = sumder / ChargeSum
      
!!!             +++++++++++++++++++++++++++++++++
         
!!! --> Add in residual for activity of water
      
      der_residuals(pos_gammawater,:) = 0D0
      
      gammawater(jx,jy,jz) = EXP( lngammawater(jx,jy,jz) )
      gammaH2O = gammawater(jx,jy,jz)
      
      TotalMoles = 0.0d0
      DO ik = 1,ncomp+nspec
        ulabPrint = ulab(ik)
        IF (ulabPrint(1:3) /= 'H2O' .and. ulabPrint(1:3) /= 'HHO') THEN
          TotalMoles = TotalMoles + sp10(ik,jx,jy,jz)
        END IF
      END DO
      
      ind = (j-1)*(neqn) + pos_gammawater
      
      CheckH2O = (1.0d0 - 0.017d0*TotalMoles)
      
      if (CheckH2O <= 0.d0) then
        write(*,*)
        write(*,*) ' Blew Sky High'
        write(*,*) ' No convergence:  !! FAIL !! '
        write(*,*) 'CheckH2O', CheckH2O
        STOP
      end IF
      
      fxx(ind) = LOG( CheckH2O/ gammaH2O )
        
!!!   Derivative of activity of water (gammawater) with respect to primary species 

      sumder = 0D0
      DO icomp = 1, ncomp
            sumder = 0D0

        ulabPrint = ulab(icomp)
        IF (ulabPrint(1:3) /= 'H2O' .and. ulabPrint(1:3) /= 'HHO') THEN
          
          sumder = sp10(icomp,jx,jy,jz)
          DO ksp = 1, nspec
            nk = ncomp + ksp
            sumder = sumder + deriv_conc(nk,icomp,jx,jy,jz)
          END DO

        END IF

        der_residuals(pos_gammawater,icomp) = - 0.017D0 * sumder / (1.0D0 - 0.017D0 * TotalMoles) 

      END DO  

!!! Derivative of activity of water (gammawater) with respect to ionic strength
      
      pos_der = ncomp+nexchange+nsurf+npot + 1
      sumder = 0D0
      DO ksp = 1, nspec
        nk = ncomp + ksp
        sumder = sumder + deriv_conc(nk,pos_der,jx,jy,jz)
      END DO
      der_residuals(pos_gammawater,pos_der) = - 0.017D0 * sumder / (1.0D0 - 0.017D0 * TotalMoles)
!!!      der_residuals(pos_gammawater,pos_der) = 0.0d0
      
!!! Derivative of activity of water (gammawater) with respect to gamma water (diagonal)
      
      pos_der = ncomp+nexchange+nsurf+npot + 1 + 1
      sumder = 0D0
      DO ksp = 1, nspec
        nk = ncomp + ksp
        sumder = sumder + deriv_conc(nk,pos_der,jx,jy,jz)
      END DO
      der_residuals(pos_gammawater,pos_der) = - 0.017D0 * sumder / (1.0D0 - 0.017D0 * TotalMoles) - 1.0d0

      
!!! -->  der_residuals(pos_IonS,pos_der)
!!! -->  der_residuals(pos_gammawater,pos_der)
    
      
!!! ---> First, residual for ionic strength
      
      do icomp2 = 1,ncomp
        alf(icomp2,pos_IonS,2)  = der_residuals(pos_IonS,icomp2)
        IF (nxyz == nx .AND. ihindmarsh == 1 .AND. nxyz /= 1) THEN
          aah(pos_IonS,icomp2,jx) = der_residuals(pos_IonS,icomp2)
        endif 
      end do
      
      alf(pos_IonS,pos_IonS,2)  = der_residuals(pos_IonS,pos_IonS)
      IF (nxyz == nx .AND. ihindmarsh == 1 .AND. nxyz /= 1) THEN
        aah(pos_IonS,pos_IonS,jx) = der_residuals(pos_IonS,pos_IonS)
      endif
      
      alf(pos_gammawater,pos_IonS,2)  = der_residuals(pos_IonS,pos_gammawater)
      IF (nxyz == nx .AND. ihindmarsh == 1 .AND. nxyz /= 1) THEN
        aah(pos_IonS,pos_gammawater,jx) = der_residuals(pos_IonS,pos_gammawater) 
      endif
      
!!! Then, residual for gammawater
      
      do icomp2 = 1,ncomp
        alf(icomp2,pos_gammawater,2)  = der_residuals(pos_gammawater,icomp2)
        IF (nxyz == nx .AND. ihindmarsh == 1 .AND. nxyz /= 1) THEN
          aah(pos_gammawater,icomp2,jx) = der_residuals(pos_gammawater,icomp2)
        endif
      end do
      
      alf(pos_gammawater,pos_gammawater,2)  = der_residuals(pos_gammawater,pos_gammawater)
      IF (nxyz == nx .AND. ihindmarsh == 1 .AND. nxyz /= 1) THEN
        aah(pos_gammawater,pos_gammawater,jx) = der_residuals(pos_gammawater,pos_gammawater)
      endif
      
      alf(pos_IonS,pos_gammawater,2)  = der_residuals(pos_gammawater,pos_IonS)
      IF (nxyz == nx .AND. ihindmarsh == 1 .AND. nxyz /= 1) THEN
        aah(pos_gammawater,pos_IonS,jx) = der_residuals(pos_gammawater,pos_IonS)
      endif
    
    CALL SYSTEM_CLOCK(COUNT=chemistry_end)
    ion_water_wall_total = ion_water_wall_total + &
      REAL(chemistry_end-chemistry_start,REAL64)/REAL(count_rate,REAL64)

!  Case where nxyz = 1
    
    IF (nxyz == 1) THEN
      DO i = 1,neqn
        DO i2 = 1,neqn
          ind2 = i2
          aaa(i,i2) = alf(ind2,i,2)
        END DO
      END DO
      GO TO 900
    END IF
    
!  For fixed concentration boundary conditions
    
800 IF (activecell(jx,jy,jz) == 0) THEN
      
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
    
    CALL SYSTEM_CLOCK(COUNT=phase_end)
    dt_local = REAL(phase_end-phase_start,REAL64)/REAL(count_rate,REAL64)
    local_wall_total = local_wall_total + dt_local
    CALL SYSTEM_CLOCK(COUNT=phase_start)

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

!!!    END IF

!***************************insert update of petsc index **********************    
    END IF  ! PETSc path versus Hindmarsh path
    ipetsc = ipetsc + 1
!***************************end update of petsc index ************************
    CALL SYSTEM_CLOCK(COUNT=phase_end)
    dt_insert = REAL(phase_end-phase_start,REAL64)/REAL(count_rate,REAL64)
    insert_wall_total = insert_wall_total + dt_insert
  END DO  ! End of JX loop
END DO    ! End of JY loop

900 CONTINUE

1077 FORMAT(12(1X,1PE10.3))

! ***************************insert construction for PETSc matrix ***********
 CALL SYSTEM_CLOCK(COUNT=phase_start)
 if (petscon) then
          call MatAssemblyBegin(amatpetsc,MAT_FINAL_ASSEMBLY,ierr)
          call MatAssemblyEnd(amatpetsc,MAT_FINAL_ASSEMBLY,ierr)
 endif
 
 CALL SYSTEM_CLOCK(COUNT=count_end)
 final_wall_total = final_wall_total + REAL(count_end-phase_start,REAL64)/REAL(count_rate,REAL64)
 assemble_wall_total = assemble_wall_total + REAL(count_end-call_start,REAL64)/REAL(count_rate,REAL64)
 assemble_wall_time = assemble_wall_time + REAL(count_end-call_start,REAL64)/REAL(count_rate,REAL64)
 timing_calls = timing_calls + 1
 
 !!! IF (timing_calls == 1 .OR. MOD(timing_calls,100) == 0) THEN
 !!!   WRITE(*,*) "AssembleGlobal cumulative seconds (calls):", timing_calls
 !!!   WRITE(*,*) "  Local compute :", local_wall_total
 !!!    WRITE(*,*) "  PETSc insertion:", insert_wall_total
 !!!    WRITE(*,*) "  Final assembly :", final_wall_total
 !!!    WRITE(*,*) "  Other overhead :", assemble_wall_total-local_wall_total-insert_wall_total-final_wall_total
 !!!  END IF
! ******************end construction for PETSc matrix ********************

RETURN
END SUBROUTINE AssembleGlobal
