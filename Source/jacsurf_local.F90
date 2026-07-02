!!! *** Copyright Notice ***
!!! “CrunchFlow”, Copyright (c) 2016, The Regents of the University of California, through Lawrence Berkeley National Laboratory 
!!! (subject to receipt of any required approvals from the U.S. Dept. of Energy).  All rights reserved.
!!! 
!!! If you have questions about your rights to use or distribute this software, please contact 
!!! Berkeley Lab's Innovation & Partnerships Office at  IPO@lbl.gov.
!!! 
!!! NOTICE.  This Software was developed under funding from the U.S. Department of Energy and the U.S. Government 
!!! consequently retains certain rights. As such, the U.S. Government has been granted for itself and others acting 
!!! on its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the Software to reproduce, distribute copies to the public, 
!!! prepare derivative works, and perform publicly and display publicly, and to permit other to do so.
!!!
!!! *** License Agreement ***
!!! “CrunchFlow”, Copyright (c) 2016, The Regents of the University of California, through Lawrence Berkeley National Laboratory)
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
    
SUBROUTINE jacsurf_local(ncomp,nsurf,nsurf_sec,jx,jy,jz)
USE crunchtype
USE params
USE concentration
USE solver

IMPLICIT NONE
!  External variables

INTEGER(I4B), INTENT(IN)                             :: ncomp
INTEGER(I4B), INTENT(IN)                             :: nsurf
INTEGER(I4B), INTENT(IN)                             :: nsurf_sec
INTEGER(I4B), INTENT(IN)                             :: jx
INTEGER(I4B), INTENT(IN)                             :: jy
INTEGER(I4B), INTENT(IN)                             :: jz

!  Internal variables

REAL(DP)                                             :: sum
REAL(DP)                                             :: mutemp
REAL(DP)                                             :: surftemp
REAL(DP)                                             :: delta_z


INTEGER(I4B)                                         :: i2
INTEGER(I4B)                                         :: iss
INTEGER(I4B)                                         :: ns
INTEGER(I4B)                                         :: is
INTEGER(I4B)                                         :: i
INTEGER(I4B)                                         :: l

CHARACTER (LEN=3)                                    :: ulabprint

fsurf_local = 0.0

DO ns = 1,nsurf_sec
  surftemp = spsurf10(ns+nsurf,jx,jy,jz)
  DO is = 1,nsurf
    iss = ncomp + is
    IF (musurf(ns,iss) /= 0.0) THEN
      mutemp = musurf(ns,iss)
      DO i2 = 1,ncomp+nsurf
        fsurf_local(is,i2) = fsurf_local(is,i2) +    &
           mutemp*musurf(ns,i2)*surftemp
      END DO
    END IF
  END DO
END DO

DO is = 1,nsurf
  iss = ncomp + is
  fsurf_local(is,iss) = fsurf_local(is,iss) + spsurf10(is,jx,jy,jz)
END DO

DO ns = 1,nsurf_sec

    !-----------------------------------------------------------------
    ! d spsurf(ns+nsurf) / d sp(i)         for i = 1..ncomp
    !-----------------------------------------------------------------
    DO i = 1,ncomp
     
      ulabPrint = ulab(i)
      IF (ulabPrint(1:3) == 'H2O' .OR. ulabPrint(1:3) == 'HHO') THEN
        ! Water: lnActivity = lngamma(i) only -> no sp(i) dependence
        dspsurf_dsp(ns+nsurf,i,jx,jy,jz) = 0.0d0
      ELSE
        dspsurf_dsp(ns+nsurf,i,jx,jy,jz) = musurf(ns,i)
      END IF
    END DO

    !-----------------------------------------------------------------
    ! d spsurf(ns+nsurf) / d spsurf(is)    for is = 1..nsurf
    ! (couples secondary site species to the other surface species)
    !-----------------------------------------------------------------
    DO is = 1,nsurf
      dspsurf_dspsurf(ns+nsurf,is,jx,jy,jz) = musurf(ns,is+ncomp)
    END DO
    

    !-----------------------------------------------------------------
    ! Electrostatic-only contributions:
    !   d/d LogPotential   and   d/d LogTotalSites
    !-----------------------------------------------------------------
    IF (nptlink(ns) /= 0) THEN
      delta_z = zsurf(ns+nsurf) - zsurf(islink(ns))

      dspsurf_dpot(ns+nsurf,jx,jy,jz)    = &
            -2.0d0 * musurf(ns,islink(ns)+ncomp) * delta_z

      !!!dspsurf_dlogtot(ns+nsurf,jx,jy,jz) = &
      !!!     -( musurf(ns,islink(ns)+ncomp) - 1.0d0 )
    ELSE
      dspsurf_dpot(ns+nsurf,jx,jy,jz)    = 0.0d0
     !!! dspsurf_dlogtot(ns+nsurf,jx,jy,jz) = &
     !!!       -( musurf(ns,islink(ns)+ncomp) - 1.0d0 )
    END IF

    !-----------------------------------------------------------------
    ! Chain rule for the linear-space variable spsurf10 = exp(spsurf):
    !   d spsurf10 / d sp(i) = spsurf10 * d spsurf / d sp(i)
    ! Use these forms wherever the mass-balance equations are written
    ! in concentration (rather than log-concentration) form.
    !-----------------------------------------------------------------
    DO i = 1,ncomp
      dspsurf10_dsp(ns+nsurf,i,jx,jy,jz) =                         &
           spsurf10(ns+nsurf,jx,jy,jz) *                           &
           dspsurf_dsp(ns+nsurf,i,jx,jy,jz)
    END DO

    DO is = 1,nsurf
      dspsurf10_dspsurf(ns+nsurf,is,jx,jy,jz) =                    &
           spsurf10(ns+nsurf,jx,jy,jz) *                           &
           dspsurf_dspsurf(ns+nsurf,is,jx,jy,jz)
    END DO

    IF (nptlink(ns) /= 0) THEN
      dspsurf10_dpot(ns+nsurf,jx,jy,jz) =                          &
           spsurf10(ns+nsurf,jx,jy,jz) *                           &
           dspsurf_dpot(ns+nsurf,jx,jy,jz)
    ELSE
      dspsurf10_dpot(ns+nsurf,jx,jy,jz) = 0.0d0
    END IF

!!!    dspsurf10_dlogtot(ns+nsurf,jx,jy,jz) =                         &
!!!           spsurf10(ns+nsurf,jx,jy,jz) *                           &
!!!           dspsurf_dlogtot(ns+nsurf,jx,jy,jz)

      END DO

RETURN
END SUBROUTINE jacsurf_local
!***********************************************************
