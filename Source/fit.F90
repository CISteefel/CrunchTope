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
    
SUBROUTINE fit(nbasis,ntemp,alogk0,bvec,vec,iflgint,int,ntt,NameTransfer)
USE crunchtype
USE params

IMPLICIT NONE

!  *****************  Beginning of interface blocks  *************************

INTERFACE
  SUBROUTINE ludcmp90(a,indx,d,n)
  USE crunchtype
  REAL(DP), DIMENSION(:,:), intent(in out)                   :: a
  INTEGER(I4B), DIMENSION(:), intent(out)                    :: indx
  INTEGER(I4B), INTENT(IN)                                   :: n
  REAL(DP), intent(out)                                      :: d
  END SUBROUTINE ludcmp90
END INTERFACE

INTERFACE
  SUBROUTINE lubksb90(a,indx,b,n)
  USE crunchtype
  IMPLICIT NONE
  REAL(DP),  DIMENSION(:,:), INTENT(IN)                          :: a
  REAL(DP),  DIMENSION(:), INTENT(IN OUT)                        :: b
  INTEGER(I4B),  DIMENSION(:),INTENT(IN)                         :: indx
  INTEGER(I4B), INTENT(IN)                                       :: n
  END SUBROUTINE lubksb90
END INTERFACE

!  *****************  End of interface blocks  *************************

!  External variables needed for dimensions

INTEGER(I4B), INTENT(IN)                                   :: nbasis
INTEGER(I4B), INTENT(IN)                                   :: ntemp

!  Arrays or variables passed in (external)

REAL(DP), DIMENSION(ntemp), INTENT(IN)                     :: alogk0
REAL(DP), DIMENSION(nbasis), INTENT(IN OUT)                :: bvec
REAL(DP), DIMENSION(nbasis,ntemp), INTENT(IN)              :: vec

INTEGER(I4B), INTENT(OUT)                                  :: iflgint
INTEGER(I4B), DIMENSION(8), INTENT(IN OUT)                 :: int
INTEGER(I4B), INTENT(OUT)                                  :: ntt
CHARACTER(30), INTENT(IN)                                  :: NameTransfer


! Internal arrays

INTEGER(I4B), PARAMETER                                    :: ione=1

REAL(DP), DIMENSION(:), ALLOCATABLE                        :: bvectmp2
REAL(DP), DIMENSION(:,:),ALLOCATABLE                       :: wscale

REAL(DP), DIMENSION(nbasis,8)                              :: vectmp
REAL(DP), DIMENSION(nbasis,nbasis)                         :: w
REAL(DP), DIMENSION(nbasis)                                :: bvectmp
REAL(DP)                                                   :: det

LOGICAL(LGT)                                               :: NonZeroLogK

INTEGER(I4B)                                               :: i
INTEGER(I4B)                                               :: i2
INTEGER(I4B)                                               :: nbcalc
INTEGER(I4B)                                               :: j
INTEGER(I4B)                                               :: k
INTEGER(I4B), DIMENSION(nbasis)                            :: indx

CHARACTER (LEN=1)                                          :: trans
INTEGER(I4B)                                               :: info

trans = 'N'

iflgint = 0
ntt = 0
DO i = 1,ntemp
  IF (alogk0(i) /= 500.0) THEN
    INT(i) = 1
    ntt = ntt + 1
  ELSE
    iflgint = 1
    INT(i) = 0
  END IF
END DO

NonZeroLogK = .FALSE.
DO i = 1,ntemp
  IF (int(i) /= 0) THEN
    NonZeroLogK = .TRUE.
    EXIT
  END IF
END DO

IF (.NOT. NonZeroLogK) THEN
  WRITE(*,*)
  WRITE(*,*) ' Row of equilibrium constants all with 500'
  WRITE(*,*) ' Species or mineral: ',NameTransfer
  WRITE(*,*)
  STOP
END IF

vectmp = 0.0
IF (ntt == 1) THEN
  vectmp(1,:) = vec(2,:)
  nbcalc = 1
ELSE IF (ntt == 2) THEN
  vectmp(1,:) = vec(2,:)
  vectmp(2,:) = vec(3,:)
  nbcalc = 2
ELSE
  vectmp(1:5,1:ntemp) = vec(1:5,1:ntemp)
  nbcalc = 5
END IF

bvec = 0.0
bvectmp = 0.0
DO j = 1, nbcalc
  DO i = 1, ntemp
    IF (INT(i) == 1) THEN
      bvectmp(j) = bvectmp(j) + alogk0(i)*vectmp(j,i)
    END IF
  END DO
END DO

w = 0.0
DO j = 1, nbcalc
  DO k = j, nbcalc
    DO i = 1, ntemp
      IF (INT(i) == 1) THEN
        w(j,k) = w(j,k) + vectmp(j,i)*vectmp(k,i)
      END IF
    END DO
    IF (j /= k) w(k,j) = w(j,k)
  END DO
END DO

ALLOCATE(wscale(nbcalc,nbcalc))
ALLOCATE(bvectmp2(nbcalc))

wscale(1:nbcalc,1:nbcalc) = w(1:nbcalc,1:nbcalc)
bvectmp2(1:nbcalc) = bvectmp(1:nbcalc)

!!CALL ludcmp90(wscale,indx,det,nbcalc)
!!CALL lubksb90(wscale,indx,bvectmp,nbcalc)

!!  call ludcmp90(fj,indx,det,neqn)
!!  call lubksb90(fj,indx,beta,neqn)

CALL dgetrf(nbcalc,nbcalc,wscale,nbcalc,indx,info)
CALL dgetrs(trans,nbcalc,ione,wscale,nbcalc,indx,bvectmp,nbcalc,info)

!!  CALL dgetrf(neqn,neqn,fj,neqn,indx,info)
!!  CALL dgetrs(trans,neqn,ione,fj,neqn,indx,beta,neqn,info)

IF (nbcalc == 1) THEN
  bvec(2) = bvectmp(1)
ELSE IF (nbcalc == 2) THEN
  bvec(2) = bvectmp(1)
  bvec(3) = bvectmp(2)
ELSE
  bvec(1:5) = bvectmp(1:5)
END IF

DEALLOCATE(wscale)
DEALLOCATE(bvectmp2)

RETURN
END SUBROUTINE fit
!**************************************************************
