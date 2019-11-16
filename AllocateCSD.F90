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
    
SUBROUTINE AllocateCSD(ncomp,nrct,nx,ny,nz)
USE crunchtype
USE params
USE mineral
USE NanoCrystal

IMPLICIT NONE

!!  External arrays and variables

INTEGER(I4B), INTENT(IN)                           :: ncomp
INTEGER(I4B), INTENT(IN)                           :: nrct
INTEGER(I4B), INTENT(IN)                           :: nx
INTEGER(I4B), INTENT(IN)                           :: ny
INTEGER(I4B), INTENT(IN)                           :: nz

!! Internal arrays and variables

REAL(DP)                                           :: dxCSD
INTEGER(I4B)                                       :: l

nCSD = 500

IF (ALLOCATED(CrystalSizeDistribution)) THEN
  DEALLOCATE(CrystalSizeDistribution)
END IF
ALLOCATE(CrystalSizeDistribution(nrct))
CrystalSizeDistribution = .FALSE.

CrystalSizeDistribution(3) = .TRUE.
CrystalSizeDistribution(4) = .TRUE.

IF (ALLOCATED(sigma)) THEN
  DEALLOCATE(sigma)
END IF
ALLOCATE(sigma(nrct))
sigma = 0.0d0

IF (ALLOCATED(radius)) THEN
  DEALLOCATE(radius)
END IF
ALLOCATE(radius(nCSD))
radius = 0.0d0

    dxCSD = 2.0E-06
    radius(1) = 0.0d0
    DO l = 2,nCSD
      radius(l) = radius(l-1) + dxCSD
    END DO

IF (ALLOCATED(xtvd)) THEN
  DEALLOCATE(xtvd)
END IF
ALLOCATE(xtvd(nCSD))
xtvd = 0.0d0

IF (ALLOCATED(aaCSD)) THEN
  DEALLOCATE(aaCSD)
END IF
ALLOCATE(aaCSD(nCSD))
aaCSD = 0.0d0

IF (ALLOCATED(bbCSD)) THEN
  DEALLOCATE(bbCSD)
END IF
ALLOCATE(bbCSD(nCSD))
bbCSD = 0.0d0

IF (ALLOCATED(ccCSD)) THEN
  DEALLOCATE(ccCSD)
END IF
ALLOCATE(ccCSD(nCSD))
ccCSD = 0.0d0

IF (ALLOCATED(uuCSD)) THEN
  DEALLOCATE(uuCSD)
END IF
ALLOCATE(uuCSD(nCSD))
uuCSD = 0.0d0

IF (ALLOCATED(rrCSD)) THEN
  DEALLOCATE(rrCSD)
END IF
ALLOCATE(rrCSD(nCSD))
rrCSD = 0.0d0

IF (ALLOCATED(AdjustKeqLog)) THEN
  DEALLOCATE(AdjustKeqLog)
END IF
ALLOCATE(AdjustKeqLog(nCSD,nrct))
AdjustKeqLog = 0.0d0

IF (ALLOCATED(rateCrystalSize)) THEN
  DEALLOCATE(rateCrystalSize)
END IF
ALLOCATE(rateCrystalSize(nCSD,nreactmax,nrct))
rateCrystalSize = 0.0d0

IF (ALLOCATED(silogCrystalSize)) THEN
  DEALLOCATE(silogCrystalSize)
END IF
ALLOCATE(silogCrystalSize(nCSD,nreactmax,nrct))
silogCrystalSize = 0.0d0

IF (ALLOCATED(siCrystalSize)) THEN
  DEALLOCATE(siCrystalSize)
END IF
ALLOCATE(siCrystalSize(nCSD,nreactmax,nrct))
siCrystalSize = 0.0d0

IF (ALLOCATED(AreaCrystalSize)) THEN
  DEALLOCATE(AreaCrystalSize)
END IF
ALLOCATE(AreaCrystalSize(nCSD,nrct))
AreaCrystalSize = 0.0d0

IF (ALLOCATED(nCrystal)) THEN
  DEALLOCATE(nCrystal)
END IF
ALLOCATE(nCrystal(nCSD,nrct,nx,ny,nz))
nCrystal = 0

IF (ALLOCATED(LinearGrowthRate)) THEN
  DEALLOCATE(LinearGrowthRate)
END IF
ALLOCATE(LinearGrowthRate(nCSD,nrct,nx,ny,nz))
LinearGrowthRate = 0.0d0

RETURN
END SUBROUTINE AllocateCSD