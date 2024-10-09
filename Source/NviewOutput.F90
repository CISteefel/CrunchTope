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
    

SUBROUTINE NviewOutput(ncomp,nrct,nkin,nspec,nexchange,nexch_sec,nsurf,  &
    nsurf_sec,ndecay,ikin,nx,ny,nz,realtime,nn,nint,ikmast,ikph,ikO2,master,delt)
USE crunchtype
USE params
USE concentration
USE mineral
USE solver
USE medium
USE transport
USE flow
USE temperature
USE strings

IMPLICIT NONE

!  External variables and arrays

CHARACTER (LEN=2*mls)                              :: ltitle

REAL(DP), INTENT(IN)                               :: realtime
REAL(DP), INTENT(IN)                               :: delt

INTEGER(I4B), INTENT(IN)                           :: ncomp
INTEGER(I4B), INTENT(IN)                           :: nrct
INTEGER(I4B), INTENT(IN)                           :: nspec
INTEGER(I4B), INTENT(IN)                           :: ndecay
INTEGER(I4B), INTENT(IN)                           :: nsurf
INTEGER(I4B), INTENT(IN)                           :: nsurf_sec
INTEGER(I4B), INTENT(IN)                           :: ikin
INTEGER(I4B), INTENT(IN)                           :: nkin
INTEGER(I4B), INTENT(IN)                           :: nexchange
INTEGER(I4B), INTENT(IN)                           :: nexch_sec
INTEGER(I4B), INTENT(IN)                           :: nx
INTEGER(I4B), INTENT(IN)                           :: ny
INTEGER(I4B), INTENT(IN)                           :: nz
INTEGER(I4B), INTENT(IN)                           :: nn
INTEGER(I4B), INTENT(IN)                           :: nint
INTEGER(I4B), INTENT(IN)                           :: ikmast
INTEGER(I4B), INTENT(IN)                           :: ikph
INTEGER(I4B), INTENT(IN)                           :: ikO2

CHARACTER (LEN=mls), INTENT(IN)                    :: master

!  Internal variables and arrays

CHARACTER (LEN=13), DIMENSION(nrct)                :: uminprnt
CHARACTER (LEN=13), DIMENSION(ncomp+nspec)         :: ulabprnt
CHARACTER (LEN=20)                                 :: fn
CHARACTER (LEN=4)                                  :: suf
CHARACTER (LEN=4)                                  :: suf1
CHARACTER (LEN=20)                                 :: fnv
CHARACTER (LEN=1)                                  :: tab
CHARACTER (LEN=mls), DIMENSION(nsurf+nsurf_sec)    :: prtsurf
CHARACTER (LEN=mls)                                :: tempstring
CHARACTER (LEN=mls)                                :: tempstring2
 
INTEGER(I4B), DIMENSION(ncomp+nspec)               :: len_sp
INTEGER(I4B), DIMENSION(nrct)                      :: len_min
INTEGER(I4B)                                       :: j
INTEGER(I4B)                                       :: jx
INTEGER(I4B)                                       :: jy
INTEGER(I4B)                                       :: jz
INTEGER(I4B)                                       :: ilength
INTEGER(I4B)                                       :: ik
INTEGER(I4B)                                       :: k
INTEGER(I4B)                                       :: ks
INTEGER(I4B)                                       :: ns
INTEGER(I4B)                                       :: i
INTEGER(I4B)                                       :: nex
INTEGER(I4B)                                       :: ir
INTEGER(I4B)                                       :: lsjx
INTEGER(I4B)                                       :: nlen
INTEGER(I4B)                                       :: nxyz
INTEGER(I4B)                                       :: noutput

REAL(DP), DIMENSION(ncomp)                         :: totex_bas
REAL(DP), DIMENSION(nrct)                          :: dptprt
REAL(DP), DIMENSION(nrct)                          :: dsat
REAL(DP), DIMENSION(nrct)                          :: dvolpr
REAL(DP)                                           :: sum
REAL(DP)                                           :: porprt
REAL(DP)                                           :: phprt
REAL(DP)                                           :: porcalc
REAL(DP)                                           :: times

nxyz = nx*ny*nz
times = realtime*31557600.
noutput = 0
432   FORMAT(e14.7)
465   FORMAT('       ',I3)

 
! pH
IF (ikph /= 0) THEN
   noutput = noutput + 1
   WRITE(32,432) times
   WRITE(32,465) noutput
   WRITE(32,465) nxyz
   DO jx=1,nx
      DO jy=1,ny
         j = (jy-1)*nx + jx
         DO jz=1,nz
            WRITE(32,432) -(sp(ikph,jx,jy,jz)+lngamma(ikph,jx,jy,jz))/clg
         ENDDO
      ENDDO
   ENDDO
ENDIF   ! END OF pH

! Total concentration
DO i=1,ncomp
   noutput = noutput + 1
   WRITE(32,432) times
   WRITE(32,465) noutput
   WRITE(32,465) nxyz
   DO jx=1,nx
      DO jy=1,ny
         DO jz=1,nz
            WRITE(32,432) s(i,jx,jy,jz)
         ENDDO
      ENDDO
   ENDDO
ENDDO   ! END OF ncomp

! Mineral volume fraction
DO k=1,nrct
   noutput = noutput + 1
   WRITE(32,432) times
   WRITE(32,465) noutput
   WRITE(32,465) nxyz
   DO jx=1,nx
      DO jy=1,ny
         DO jz=1,nz
            WRITE(32,432) volfx(k,jx,jy,jz)*100.0
         ENDDO
      ENDDO
   ENDDO
ENDDO   ! END OF nrct

! Porosity
IF (nrct > 1) THEN
   noutput = noutput + 1
   WRITE(32,432) times
   WRITE(32,*) noutput
   WRITE(32,*) nxyz
   DO jx=1,nx
      DO jy=1,ny
         DO jz=1,nz
            sum = 0.0
            DO k = 1,nrct
              if (mintype(k) == 0) sum = sum + volfx(k,jx,jy,jz)
            END DO
            WRITE(32,432) (1.0-sum)*100.0
         ENDDO
      ENDDO
   ENDDO
ENDIF

! Adsorbates
DO i=1,nexch_sec
   noutput = noutput + 1
   WRITE(32,432) times
   WRITE(32,465) noutput
   WRITE(32,465) nxyz
   DO jx=1,nx
      DO jy=1,ny
         DO jz=1,nz
            WRITE(32,432) spex10(i+nexchange,jx,jy,jz)
         ENDDO
      ENDDO
   ENDDO
ENDDO   ! END OF nexch_sec

RETURN
END SUBROUTINE NviewOutput

