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
    
SUBROUTINE SurfaceComplex(ncomp,nsurf,nsurf_sec,nx,ny,nz)
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
INTEGER(I4B), INTENT(IN)                                    :: nsurf
INTEGER(I4B), INTENT(IN)                                    :: nsurf_sec
INTEGER(I4B), INTENT(IN)                                    :: nx
INTEGER(I4B), INTENT(IN)                                    :: ny
INTEGER(I4B), INTENT(IN)                                    :: nz

!  Internal variables and arrays

INTEGER(I4B)                                                :: ns
INTEGER(I4B)                                                :: i
INTEGER(I4B)                                                :: is
INTEGER(I4B)                                                :: jx
INTEGER(I4B)                                                :: jy
INTEGER(I4B)                                                :: jz

REAL(DP)                                                    :: sum
REAL(DP)                                                    :: delta_z
REAL(DP)                                                    :: AqueousToBulk
REAL(DP)                                                    :: LogAqueousToBulk
REAL(DP)                                                    :: MeanSaltConcentration
REAL(DP)                                                    :: MassFraction
REAL(DP)                                                    :: activity
REAL(DP)                                                    :: LogTotalSites
REAL(DP)                                                    :: LogTotalEquivalents

DO jz = 1,nz
  DO jy = 1,ny
    DO jx = 1,nx

!!      CALL AqueousToBulkConvert(jx,jy,jz,AqueousToBulk)
!!      LogAqueousToBulk = DLOG(AqueousToBulk)

      DO ns = 1,nsurf_sec

        IF (nptlink(ns) /= 0) THEN                            !  Electrostatic correction

          delta_z = zsurf(ns+nsurf) - zsurf(islink(ns))

          IF (ikh2o /= 0) THEN

            sum = 0.0
            DO i = 1,ncomp
              IF (ulab(i) == 'H2O') THEN
                sum = sum + musurf(ns,i)*(gam(i,jx,jy,jz))
              ELSE
                sum = sum + musurf(ns,i)*(sp(i,jx,jy,jz) + gam(i,jx,jy,jz))
              END IF
            END DO

          ELSE

            sum = 0.0
            DO i = 1,ncomp
              sum = sum + musurf(ns,i)*(sp(i,jx,jy,jz) + gam(i,jx,jy,jz))
            END DO

          END IF

          LogTotalSites = LogTotalSurface(islink(ns),jx,jy,jz) 

!!  Surface complexes
          DO is = 1,nsurf

            activity = spsurf(is,jx,jy,jz)
            sum = sum + musurf(ns,is+ncomp)*activity

          END DO

!! NOTE:  Below is the LOg concentration of sites in units of mol/kgw
          spsurf(ns+nsurf,jx,jy,jz) = keqsurf(ns,jx,jy,jz) + sum -                     &    
        + 2.0d0*musurf(ns,islink(ns)+ncomp)*delta_z*LogPotential(nptlink(ns),jx,jy,jz)   &
        - (musurf(ns,islink(ns)+ncomp)-1.0d0)*LogTotalSites                   &
        - DLOG(musurf(ns,islink(ns)+ncomp)) 
          
          spsurf10(ns+nsurf,jx,jy,jz) = DEXP( spsurf(ns+nsurf,jx,jy,jz) )

        ELSE                                                  !  Non-electrostatic 

          IF (ikh2o /= 0) THEN

            sum = 0.0
            DO i = 1,ncomp
              IF (ulab(i) == 'H2O') THEN
                sum = sum + musurf(ns,i)*(gam(i,jx,jy,jz))
              ELSE
                sum = sum + musurf(ns,i)*(sp(i,jx,jy,jz) + gam(i,jx,jy,jz))
              END IF
            END DO

          ELSE

            sum = 0.0
            DO i = 1,ncomp
              sum = sum + musurf(ns,i)*(sp(i,jx,jy,jz) + gam(i,jx,jy,jz))
            END DO

          END IF
          
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

    END DO
  END DO
END DO

RETURN
END SUBROUTINE SurfaceComplex
!  **************************************************************
