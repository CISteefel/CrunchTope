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
    
SUBROUTINE SurfaceComplex(ncomp,nspec,nsurf,nsurf_sec,nx,ny,nz,igamma)
USE crunchtype
USE concentration
USE mineral
USE medium
USE transport
USE temperature
USE RunTime

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: ncomp
INTEGER(I4B), INTENT(IN)                                    :: nspec
INTEGER(I4B), INTENT(IN)                                    :: nsurf
INTEGER(I4B), INTENT(IN)                                    :: nsurf_sec
INTEGER(I4B), INTENT(IN)                                    :: nx
INTEGER(I4B), INTENT(IN)                                    :: ny
INTEGER(I4B), INTENT(IN)                                    :: nz
INTEGER(I4B), INTENT(IN)                                    :: igamma

!  Internal variables and arrays

INTEGER(I4B)                                                :: ns
INTEGER(I4B)                                                :: i
INTEGER(I4B)                                                :: is
INTEGER(I4B)                                                :: jx
INTEGER(I4B)                                                :: jy
INTEGER(I4B)                                                :: jz
INTEGER(I4B)                                                :: ik

REAL(DP)                                                    :: sum
REAL(DP)                                                    :: delta_z
REAL(DP)                                                    :: AqueousToBulk
REAL(DP)                                                    :: LogAqueousToBulk
REAL(DP)                                                    :: MeanSaltConcentration
REAL(DP)                                                    :: MassFraction
REAL(DP)                                                    :: activity
REAL(DP)                                                    :: LogTotalSites
REAL(DP)                                                    :: LogTotalEquivalents
REAL(DP)                                                    :: check

REAL(DP)                                                    :: sum_dlnI
REAL(DP)                                                    :: sion_tmp
REAL(DP)                                                    :: aa1
REAL(DP)                                                    :: ah
REAL(DP)                                                    :: sqrt_ions
REAL(DP)                                                    :: bdt
REAL(DP)                                                    :: bdotpar
REAL(DP)                                                    :: bh
REAL(DP)                                                    :: tmp1
REAL(DP)                                                    :: tmp2
REAL(DP)                                                    :: k

LOGICAL(LGT)                                              :: Davies
LOGICAL(LGT)                                              :: Wateq_Extended_DH
LOGICAL(LGT)                                              :: Helgeson
LOGICAL(LGT)                                              :: Unity

REAL(DP)                                                        :: lnActivity
CHARACTER (LEN=3)                                               :: ulabPrint

DO jz = 1,nz
  DO jy = 1,ny
    DO jx = 1,nx

!!      CALL AqueousToBulkConvert(jx,jy,jz,AqueousToBulk)
!!      LogAqueousToBulk = DLOG(AqueousToBulk)

      DO ns = 1,nsurf_sec

        IF (nptlink(ns) /= 0) THEN                            !  Electrostatic correction

          delta_z = zsurf(ns+nsurf) - zsurf(islink(ns))       
          sum = 0.0
          DO i = 1,ncomp
              
            ulabPrint = ulab(i)
            IF (ulabPrint(1:3) == 'H2O' .or. ulabPrint(1:3) == 'HHO') THEN
              lnActivity = lngamma(i,jx,jy,jz)
            ELSE
              lnActivity = sp(i,jx,jy,jz) + lngamma(i,jx,jy,jz)
            END IF
            sum = sum + musurf(ns,i)*lnActivity
            continue
          END DO
          
          LogTotalSites = LogTotalSurface(islink(ns),jx,jy,jz) 

!!  Surface complexes
          DO is = 1,nsurf
            activity = spsurf(is,jx,jy,jz)
            sum = sum + musurf(ns,is+ncomp)*Activity
          END DO

          spsurf(ns+nsurf,jx,jy,jz) = keqsurf(ns,jx,jy,jz) + sum                                &    
            - 2.0d0*musurf(ns,islink(ns)+ncomp) * delta_z * LogPotential(nptlink(ns),jx,jy,jz)  &
            - (musurf(ns,islink(ns)+ncomp)-1.0d0) * LogTotalSites                               &
            - DLOG(musurf(ns,islink(ns)+ncomp)) 
          
          spsurf10(ns+nsurf,jx,jy,jz) = DEXP( spsurf(ns+nsurf,jx,jy,jz) )
          
!!!  Surface Complexation Cheat Sheet    
!!!    kPotential(k) --> Logical to EDL potential
!!!    ksurf(is) --> pointer for primary nsurf complex to mineral (initialized in read_surface.F90)
!!!    iedl(is) --> 0 for electrostatic, 1 for -no_edl
!!!    npot --> number of potentials
!!!    kpot(npt) --> pointer to mineral upon which the potential is developed
!!!    islink(ns) --> pointer from secondary surface complex (ns) to primary surface complex (is)
!!!    ksurf(islink(ns)) --> This would point from a secondary surface complex (ns) to a primary (islink(ns)) complex to a mineral
!!!    nptlink(ns) --> pointer of surface complex (primary or secondary) to potential (npt)
          

        ELSE                                                  !  Non-electrostatic 

          sum = 0.0
          DO i = 1,ncomp
              
            ulabPrint = ulab(i)
            IF (ulabPrint(1:3) == 'H2O' .or. ulabPrint(1:3) == 'HHO') THEN
              lnActivity = lngamma(i,jx,jy,jz)
            ELSE
              lnActivity = sp(i,jx,jy,jz) + lngamma(i,jx,jy,jz)
            END IF
          
            sum = sum + musurf(ns,i)*lnActivity
            
          END DO
          
          LogTotalSites = LogTotalSurface(islink(ns),jx,jy,jz) 

          DO is = 1,nsurf
!!!            activity = spsurf(is,jx,jy,jz) - LogTotalSites
            activity = spsurf(is,jx,jy,jz)
            sum = sum + musurf(ns,is+ncomp)*activity
          END DO
          
          spsurf(ns+nsurf,jx,jy,jz) = keqsurf(ns,jx,jy,jz) + sum                  &
               - (musurf(ns,islink(ns)+ncomp)-1.0d0)*LogTotalSites                   &
               - DLOG(musurf(ns,islink(ns)+ncomp)) 
          
          spsurf10(ns+nsurf,jx,jy,jz) = DEXP( spsurf(ns+nsurf,jx,jy,jz) )

        END IF
      END DO
      
      
!!! Derivatives
      
!-----------------------------------------------------------------------
! Build dlngamma_dlnI(i,jx,jy,jz) = d ln gamma(i) / d ln(I)
! Mirrors the four branches in your activity-coefficient block.
! Water: 0.  Unity model: 0.  Neutral: clg*0.1*I.
! Davies / Helgeson / Wateq_Extended_DH: closed forms below.
!-----------------------------------------------------------------------
  dlngamma_dlnI(:,jx,jy,jz) = 0.0d0

  IF (igamma /= 0) THEN

    DO ik = 1,ncomp+nspec
      ulabPrint = ulab(ik)

      IF (ulabPrint(1:3) == 'H2O' .OR. ulabPrint(1:3) == 'HHO') THEN
        dlngamma_dlnI(ik,jx,jy,jz) = 0.0d0          ! water decoupled from I

      ELSE IF (chg(ik) == 0.0d0) THEN               ! neutral non-water
        dlngamma_dlnI(ik,jx,jy,jz) = clg * 0.1d0 * sion_tmp

      ELSE                                          ! charged species
        IF (Davies) THEN
          aa1 = - clg * ah * chg(ik) * chg(ik)
          dlngamma_dlnI(ik,jx,jy,jz) = aa1 *                              &
               ( 0.5d0 * sqrt_IonS / ((1.0d0 + sqrt_IonS)**2)             &
                 - 0.3d0 * sion_tmp )

        ELSE IF (Helgeson .OR. Wateq_Extended_DH) THEN
          IF (Helgeson) THEN
            bdotpar = bdt
          ELSE
            bdotpar = bdotparameter(ik)
          END IF
          tmp1 = 1.0d0 + bh * azero(ik) * sqrt_IonS
          tmp2 = ah * chg(ik) * chg(ik)
          dlngamma_dlnI(ik,jx,jy,jz) = clg *                              &
               ( - 0.5d0 * tmp2 * sqrt_IonS / (tmp1*tmp1)                 &
                 + bdotpar * sion_tmp )
        END IF
      END IF
    END DO

  END IF
  
!-----------------------------------------------------------------------
! Jacobian of nsurf_sec secondary surface species (non-electrostatic).
! Unknowns differentiated:
!   sp(k)         k = 1..ncomp        (log primary aqueous concentrations)
!   lngammawater                       (= ln a_H2O,  at pos_gammawater)
!   ln(I)         ionic strength       (at pos_IonS, log-space convention)
!
! Frozen-gamma w.r.t. composition: d lngamma(i)/d sp(k) = 0.
! All ionic-strength dependence of lngamma(i) carried in dlngamma_dlnI.
!-----------------------------------------------------------------------
      
 DO ns = 1,nsurf_sec

    !---------------- d / d sp(k) -----------------------------------
    DO k = 1,ncomp
      ulabPrint = ulab(k)
      IF (ulabPrint(1:3) == 'H2O' .OR. ulabPrint(1:3) == 'HHO') THEN
        dspsurf_dsp  (ns+nsurf,k,jx,jy,jz) = 0.0d0
        dspsurf10_dsp(ns+nsurf,k,jx,jy,jz) = 0.0d0
      ELSE
        dspsurf_dsp  (ns+nsurf,k,jx,jy,jz) = musurf(ns,k)
        dspsurf10_dsp(ns+nsurf,k,jx,jy,jz) =                            &
             spsurf10(ns+nsurf,jx,jy,jz) * musurf(ns,k)
      END IF
    END DO

    !---------------- d / d lngammawater  (= d / d ln a_H2O) --------
    !  Water enters the residual through lnActivity(ikh2o) = lngamma(ikh2o)
    !  and lngamma(ikh2o) is exactly the unknown lngammawater.
    dspsurf_dlnaH2O(ns+nsurf,jx,jy,jz) = musurf(ns,ikh2o)
    dspsurf10_dlnaH2O(ns+nsurf,jx,jy,jz) =                              &
         spsurf10(ns+nsurf,jx,jy,jz) * musurf(ns,ikh2o)

    !---------------- d / d ln(I) -----------------------------------
    !  Only non-water lngamma(i) carry I dependence in this code.
    sum_dlnI = 0.0d0
    DO i = 1,ncomp
      ulabPrint = ulab(i)
      IF (ulabPrint(1:3) == 'H2O' .OR. ulabPrint(1:3) == 'HHO') CYCLE
      sum_dlnI = sum_dlnI + musurf(ns,i) * dlngamma_dlnI(i,jx,jy,jz)
    END DO
    dspsurf_dlnI(ns+nsurf,jx,jy,jz) = sum_dlnI
    dspsurf10_dlnI(ns+nsurf,jx,jy,jz) =                                 &
    spsurf10(ns+nsurf,jx,jy,jz) * sum_dlnI

 END DO
 
 

    END DO
  END DO
END DO

RETURN
END SUBROUTINE SurfaceComplex
!  **************************************************************
