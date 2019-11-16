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
    
subroutine MineralBreakthroughInitialize(ncomp,nkin,nrct,nx,ny,nz,minseries)

USE crunchtype
USE params
USE runtime
USE mineral, ONLY: umin
USE concentration, ONLY: jxminseries,jyminseries,jzminseries,ulab

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                     :: ncomp
INTEGER(I4B), INTENT(IN)                                     :: nkin
INTEGER(I4B), INTENT(IN)                                     :: nrct
INTEGER(I4B), INTENT(IN)                                     :: nx
INTEGER(I4B), INTENT(IN)                                     :: ny
INTEGER(I4B), INTENT(IN)                                     :: nz
INTEGER(I4B), INTENT(IN)                                     :: minseries



!! INTERNAL VARIABLES

INTEGER(I4B)                                                  :: i
INTEGER(I4B)                                                  :: k
INTEGER(I4B)                                                  :: ksp
INTEGER(I4B)                                                  :: ik
INTEGER(I4B)                                                  :: ls
INTEGER(I4B)                                                  :: ns
INTEGER(I4B)                                                  :: nex

INTEGER(I4B)                                                  :: j
INTEGER(I4B)                                                  :: nxyz
INTEGER(I4B)                                                  :: jx
INTEGER(I4B)                                                  :: jy
INTEGER(I4B)                                                  :: jz
INTEGER(I4B)                                                  :: ll
INTEGER(I4B)                                                  :: intfile

CHARACTER (LEN=mls)                                           :: filename
CHARACTER (LEN=mls)                                           :: dummy
CHARACTER (LEN=mls)                                           :: dumstring
CHARACTER (LEN=12)                                            :: writeph

LOGICAL(LGT)                                                  :: ext


CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE                :: MineralBasis

IF (minseries >= 1) THEN
  DO ll = 1,minseries
    intfile = 150+ll  
    IF (irestart == 1 .AND. AppendRestart) THEN    !  Open breakthrough file and go to end of file
      filename = TimeSeriesMineralFile(ll)
      INQUIRE(FILE=filename,EXIST=ext)
      IF (.NOT. ext) THEN
        CALL stringlen(filename,ls)
        WRITE(*,*) 
        WRITE(*,*) ' Cannot find mineral time series file: ', filename(1:ls)
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
      OPEN(UNIT=intfile,FILE=filename,STATUS='unknown',ERR=702,POSITION='append')
    ELSE
      filename = TimeSeriesMineralFile(ll)
      OPEN(UNIT=intfile,FILE=filename,STATUS='unknown',ERR=702) 
!!CIS: Moved to individual plotting options      WRITE(intfile,2283) jxseries(ll),jyseries(ll),jzseries(ll)
2283 FORMAT('# Mineral time series at grid cell: 'i3,1x,i3,1x,i3)      

      IF (OutputTimeUnits == 'years') THEN

        IF (tecplot) THEN
          WRITE(intfile,2283) jxminseries(ll),jyminseries(ll),jzminseries(ll)
!!!          WRITE(intfile,3001) (ulab(i),i=1,ncomp) 
          WRITE(intfile,3701) (umin(k),k=1,nrct) 

        ELSE
          WRITE(intfile,2283) jxminseries(ll),jyminseries(ll),jzminseries(ll)
!!!          WRITE(intfile,3001) (ulab(i),i=1,ncomp) 
          WRITE(intfile,3701) (umin(k),k=1,nrct) 
          
        ENDIF
        
      ELSE IF (OutputTimeUnits == 'days') THEN
          
!!!  Need 3002 and 3702
          
        IF (tecplot) THEN
          WRITE(intfile,2283) jxminseries(ll),jyminseries(ll),jzminseries(ll)
!!!          WRITE(intfile,3702) (ulab(i),i=1,ncomp) 
          WRITE(intfile,3702) (umin(k),k=1,nrct)   

        ELSE
          WRITE(intfile,2283) jxminseries(ll),jyminseries(ll),jzminseries(ll)
!!!          WRITE(intfile,3702) (ulab(i),i=1,ncomp) 
          WRITE(intfile,3702) (umin(k),k=1,nrct)  
          
        ENDIF

      ELSE IF (OutputTimeUnits == 'hours') THEN
          
!!!  Need 3003 and 3703
          
        IF (tecplot) THEN
          WRITE(intfile,2283) jxminseries(ll),jyminseries(ll),jzminseries(ll)
!!!          WRITE(intfile,3702) (ulab(i),i=1,ncomp) 
          WRITE(intfile,3703) (umin(k),k=1,nrct)

        ELSE
          WRITE(intfile,2283) jxminseries(ll),jyminseries(ll),jzminseries(ll)
!!!          WRITE(intfile,3702) (ulab(i),i=1,ncomp) 
          WRITE(intfile,3703) (umin(k),k=1,nrct)
          
        ENDIF


      ELSE IF (OutputTimeUnits == 'minutes') THEN
          
        IF (tecplot) THEN
          WRITE(intfile,2283) jxminseries(ll),jyminseries(ll),jzminseries(ll)
!!!          WRITE(intfile,3702) (ulab(i),i=1,ncomp) 
          WRITE(intfile,3704) (umin(k),k=1,nrct)

        ELSE
          WRITE(intfile,2283) jxminseries(ll),jyminseries(ll),jzminseries(ll)
!!!          WRITE(intfile,3702) (ulab(i),i=1,ncomp) 
          WRITE(intfile,3704) (umin(k),k=1,nrct)
          
        ENDIF
 

      ELSE IF (OutputTimeUnits == 'seconds') THEN
          
        IF (tecplot) THEN
          WRITE(intfile,2283) jxminseries(ll),jyminseries(ll),jzminseries(ll)
!!!          WRITE(intfile,3702) (ulab(i),i=1,ncomp) 
          WRITE(intfile,3705) (umin(k),k=1,nrct)

        ELSE
          WRITE(intfile,2283) jxminseries(ll),jyminseries(ll),jzminseries(ll)
!!!          WRITE(intfile,3702) (ulab(i),i=1,ncomp) 
          WRITE(intfile,3705) (umin(k),k=1,nrct)
          
        ENDIF
        
      ELSE
        WRITE(*,*) 
        WRITE(*,*) ' Output time units not recognized'
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF

    END IF
  END DO
END IF

IF (ALLOCATED(MineralBasis)) THEN
  DEALLOCATE(MineralBasis)
END IF


!! Tecplot formats
3001 FORMAT('VARIABLES = "Time (yrs) " ',                   100(', "',A19,'"'))
3002 FORMAT('VARIABLES = "Time (days)" ',                   100(', "',A19,'"'))
3003 FORMAT('VARIABLES = "Time (hrs) " ',                   100(', "',A19,'"'))
3004 FORMAT('VARIABLES = "Time (min) " ',                   100(', "',A19,'"'))
3005 FORMAT('VARIABLES = "Time (sec) " ',                   100(', "',A19,'"'))

3006 FORMAT('  Time          ',                                       100(A17) )
3007 FORMAT('  Yrs           ',                                        100(A17) )
3008 FORMAT('  Days          ',                                        100(A17) )
3009 FORMAT('  Hrs           ',                                        100(A17) )
3010 FORMAT('  Min           ',                                        100(A17) )
3011 FORMAT('  Sec           ',                                        100(A17) )


!! Kaleidagraph formats
3701 FORMAT('  Time(yrs)',17x,150(1X,A22))
3702 FORMAT('  Time(day)',17x,150(1X,A22))
3703 FORMAT('  Time(hrs)',17x,150(1X,A22))
3704 FORMAT('  Time(min)',17x,150(1X,A22))
3705 FORMAT('  Time(sec)',17x,150(1X,A22))

RETURN

702 WRITE(*,*)
IF (ALLOCATED(MineralBasis)) THEN
  DEALLOCATE(MineralBasis)
END IF

WRITE(*,*) ' Error opening mineral breakthrough file'
WRITE(*,*)
READ(*,*)
STOP


END subroutine MineralBreakthroughInitialize