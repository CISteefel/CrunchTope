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
    
    
SUBROUTINE NviewInit(ncomp,nrct,nkin,nspec,nexchange,nexch_sec,nsurf,ltitle,  &
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
INTEGER(I4B)                                       :: noutput

REAL(DP), DIMENSION(ncomp)                         :: totex_bas
REAL(DP), DIMENSION(nrct)                          :: dptprt
REAL(DP), DIMENSION(nrct)                          :: dsat
REAL(DP), DIMENSION(nrct)                          :: dvolpr
REAL(DP)                                           :: sum
REAL(DP)                                           :: porprt
REAL(DP)                                           :: phprt
REAL(DP)                                           :: porcalc
REAL(DP)                                           :: dzzz


432   FORMAT(e14.7)
465   FORMAT('      ',I3)

!        File Header
!
OPEN(32,FILE='crunch.ext',STATUS='UNKNOWN')
WRITE(32,401) 
401   FORMAT('         1')
WRITE(32,406) ltitle
406   FORMAT(a132)
WRITE(32,402) 
402   FORMAT('         10')
WRITE(32,400)
400   FORMAT('internal')
WRITE(32,403) 
403   FORMAT('Data from CRUNCH')
WRITE(32,404) 
404   FORMAT('0')
WRITE(32,405)  
405   FORMAT('VARIABLE')

WRITE(32,407) 
407   FORMAT('$gdef')
WRITE(32,408) 
408   FORMAT('$type rect')
WRITE(32,409) nx
409   FORMAT('$nx ',I5)
WRITE(32,410) ny
410   FORMAT('$ny ',I5)
WRITE(32,411) nz
411   FORMAT('$nz ',I5)
WRITE(32,412) 
412   FORMAT('$order xyz')

WRITE(32,421) 
421   FORMAT('$dx ')
DO jx=1, nx
   WRITE(32,432) dxx(jx)
ENDDO
WRITE(32,422) 
422   FORMAT('$dy ')

DO jy=1, ny
   WRITE(32,432) dyy(jy)
ENDDO

WRITE(32,423) 
423   FORMAT('$dz ')

dzzz=1
DO jz=1, nz
   WRITE(32,432) dzzz
ENDDO

WRITE(32,424) 
424   FORMAT('$end_internal_grid')
WRITE(32,425) 
425   FORMAT('$OperatingSystem ')
WRITE(32,426) 
426   FORMAT('$C-Compiler  ')
WRITE(32,427) 
427   FORMAT('$FortranCompiler  ')
WRITE(32,428) 
428   FORMAT('$RunID  0')
WRITE(32,429) 
429   FORMAT('$RunDate ')

! number of output variables (total_conc+minerals+porosity+pH)
noutput = ncomp + nrct + nexch_sec
IF (ikph /= 0) THEN
   noutput = noutput + 1
ENDIF ! for pH
IF (nrct > 1) THEN
   noutput = noutput + 1
ENDIF ! for porosity
WRITE(32,465) noutput 

230 FORMAT(a14)
IF (ikph /= 0) THEN
   WRITE(32,230) "pH            "
ENDIF ! for pH
DO i=1,ncomp
   WRITE(32,230) ulab(i)
ENDDO
DO k=1,nrct
   WRITE(32,230) umin(k)
ENDDO
IF (nrct > 1) THEN
   WRITE(32,230) "Porosity      "
ENDIF ! for porosity
DO i=1,nexch_sec
   WRITE(32,230) nam_exchsec(i)
ENDDO

WRITE(32,430) 
430   FORMAT('        0')
WRITE(32,465) nx*ny*nz 

DO jx=1,nx
   DO jy=1,ny
      DO jz=1,nz
!         WRITE(32, 431) jx, jy, jz
         WRITE(32, *) "rock#",jx,":",jy,":",jz
      ENDDO
   ENDDO
ENDDO
!431   FORMAT(I3,":",I3,":",I3)

231 FORMAT("THIST")
232 FORMAT("EVAR ",a14)
233 FORMAT("T *")
234 FORMAT("E *")
235 FORMAT(" ")
WRITE(32,430) 
IF (ikph /= 0) THEN
   WRITE(32,231)
   WRITE(32,232) "pH            "
   WRITE(32,233)
   WRITE(32,234)
   WRITE(32,235)
ENDIF ! for pH
DO i=1,ncomp
   WRITE(32,231)
   WRITE(32,232) ulab(i)
   WRITE(32,233)
   WRITE(32,234)
   WRITE(32,235)
ENDDO
DO k=1,nrct
   WRITE(32,231)
   WRITE(32,232) umin(k)
   WRITE(32,233)
   WRITE(32,234)
   WRITE(32,235)
ENDDO
IF (nrct > 1) THEN
   WRITE(32,231)
   WRITE(32,232) "Porosity      "
   WRITE(32,233)
   WRITE(32,234)
   WRITE(32,235)
ENDIF ! for porosity
DO i=1,nexch_sec
   WRITE(32,231)
   WRITE(32,232) nam_exchsec(i)
   WRITE(32,233)
   WRITE(32,234)
   WRITE(32,235)
ENDDO


RETURN
END SUBROUTINE NviewInit

