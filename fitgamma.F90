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

SUBROUTINE fitgamma(nbasis,ntemp,adh,bvec,vec)
USE crunchtype
USE params

IMPLICIT NONE

interface
SUBROUTINE ludcmp90(a,indx,d,n)
USE crunchtype
REAL(DP), DIMENSION(:,:), intent(in out)                   :: a
INTEGER(I4B), DIMENSION(:), intent(out)                    :: indx
INTEGER(I4B), INTENT(IN)                                   :: n
REAL(DP), intent(out)                                      :: d
END SUBROUTINE ludcmp90
END interface

interface
SUBROUTINE lubksb90(a,indx,b,n)
USE crunchtype
REAL(DP), DIMENSION(:,:), intent(in)                       :: a
INTEGER(I4B), DIMENSION(:), intent(in)                     :: indx
REAL(DP), DIMENSION(:), intent(inout)                      :: b
INTEGER(I4B), INTENT(IN)                                   :: n
END SUBROUTINE lubksb90
END interface

!  External variables needed for dimensions

INTEGER(I4B)                                               :: nbasis
INTEGER(I4B)                                               :: ntemp

!  External variables and arrays

REAL(DP), DIMENSION(ntemp), INTENT(IN)                     :: adh
REAL(DP), DIMENSION(nbasis), INTENT(OUT)                   :: bvec
REAL(DP), DIMENSION(nbasis,ntemp), INTENT(IN)              :: vec

!  Internal variables and arrays

REAL(DP), DIMENSION(nbasis,nbasis)                         :: w
REAL(DP)                                                   :: dd

INTEGER(I4B), DIMENSION(nbasis)                            :: indx
INTEGER(I4B)                                               :: i
INTEGER(I4B)                                               :: k
INTEGER(I4B)                                               :: j
CHARACTER (LEN=1)                                          :: trans
INTEGER(I4B)                                               :: info
INTEGER(I4B), PARAMETER                                    :: ione=1

trans = 'N'

DO j = 1, nbasis
  bvec(j) = 0.d0
  DO i = 1, ntemp
    bvec(j) = bvec(j) + adh(i)*vec(j,i)
  END DO
END DO

DO j = 1, nbasis
  DO k = j, nbasis
    w(j,k) = 0.d0
    DO i = 1, ntemp
      w(j,k) = w(j,k) + vec(j,i)*vec(k,i)
    END DO
    IF (j /= k) w(k,j) = w(j,k)
  END DO
END DO

!!CALL ludcmp90(w,indx,dd,nbasis)
!!CALL lubksb90(w,indx,bvec,nbasis)


!!  call ludcmp90(fj,indx,det,neqn)
!!  call lubksb90(fj,indx,beta,neqn)

CALL dgetrf(nbasis,nbasis,w,nbasis,indx,info)
CALL dgetrs(trans,nbasis,ione,w,nbasis,indx,bvec,nbasis,info)

!!  CALL dgetrf(neqn,neqn,fj,neqn,indx,info)
!!  CALL dgetrs(trans,neqn,ione,fj,neqn,indx,beta,neqn,info)

RETURN
END SUBROUTINE fitgamma
!**************************************************************
