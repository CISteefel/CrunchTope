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
    
SUBROUTINE SurfaceCharge(ncomp,nspec,nsurf,nsurf_sec,npot,jx,jy,jz,time)
USE crunchtype
USE concentration
USE medium
USE temperature
USE mineral

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: ncomp
INTEGER(I4B), INTENT(IN)                                    :: nspec
INTEGER(I4B), INTENT(IN)                                    :: nsurf
INTEGER(I4B), INTENT(IN)                                    :: nsurf_sec
INTEGER(I4B), INTENT(IN)                                    :: npot
INTEGER(I4B), INTENT(IN)                                    :: jx
INTEGER(I4B), INTENT(IN)                                    :: jy
INTEGER(I4B), INTENT(IN)                                    :: jz
REAL(DP), INTENT(IN)                                        :: time

!  Internal variables and arrays

INTEGER(I4B)                                                :: ns
INTEGER(I4B)                                                :: is
INTEGER(I4B)                                                :: ik
INTEGER(I4B)                                                :: k
INTEGER(I4B)                                                :: npt

REAL(DP)                                                    :: sum
REAL(DP)                                                    :: faraday
REAL(DP)                                                    :: gramsperL
REAL(DP)                                                    :: portemp
REAL(DP)                                                    :: convert
REAL(DP)                                                    :: correct
REAL(DP)                                                    :: volMinimum
REAL(DP)                                                    :: term1

faraday = 96485.0d0     
surfcharge = 0.0d0

DO npt = 1,npot

  k = kpot(npt)

  IF (volinByGrid(k,jx,jy,jz) == 0.0d0 .AND. volfx(k,jx,jy,jz) < voltemp(k,jinit(jx,jy,jz)) ) THEN
    volMinimum = voltemp(k,jinit(jx,jy,jz))
  ELSE
    volMinimum = volfx(k,jx,jy,jz)
    IF (volMinimum < 1.0D-15) volMinimum = 1.0D-15
  END IF
  
  correct = wtmin(k)*specificByGrid(k,jx,jy,jz)*volMinimum/volmol(k) 

  if (correct == 0.0) then
    write(*,*) ' Divide by zero in surface charge'
    read(*,*)
  end if
  
  term1 = faraday/correct

  DO is = 1,nsurf
    IF ( ksurf(is) == kpot(npt) ) THEN
      surfcharge(k) = surfcharge(k) + zsurf(is)*spsurf10(is,jx,jy,jz)*term1
    END IF
  END DO

  DO ns = 1,nsurf_sec
    IF (ksurf(islink(ns)) == kpot(npt)) THEN
      surfcharge(k) = surfcharge(k) + zsurf(ns+nsurf)*spsurf10(ns+nsurf,jx,jy,jz)*term1   !!  Equation 2.1 in Dzombak
    END IF
  END DO

!!!  Surface Complexation Cheat Sheet    
!!!    kPotential(k) --> Logical to EDL potential
!!!    ksurf(is) --> pointer for primary nsurf complex to mineral (initialized in read_surface.F90)
!!!    iedl(is) --> 0 for electrostatic, 1 for -no_edl
!!!    npot --> number of potentials
!!!    kpot(npt) --> pointer to mineral upon which the potential is developed
!!!    islink(ns) --> pointer from secondary surface complex (ns) to primary surface complex (is)
!!!    ksurf(islink(ns)) --> This would point from a secondary surface complex (ns) to a primary (islink(ns)) complex to a mineral
!!!    nptlink(ns) --> pointer of surface complex (primary or secondary) to potential (npt)
          

END DO


RETURN
END SUBROUTINE SurfaceCharge
