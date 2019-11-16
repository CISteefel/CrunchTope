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

SUBROUTINE ModScan(nx,ny,nz,cnhIn,wellIn,riverIn,ndimdummy,mxwell,mxrivr,mxdrn)

USE crunchtype
USE medium
USE transport
USE flow
USE modflowModule

IMPLICIT NONE
 
!  *********************  INTERFACE BLOCKS  ********************

INTERFACE
  SUBROUTINE readtype1(ndim,iii,jjj,kkk,xxx,num)
  USE crunchtype
  IMPLICIT NONE
  INTEGER(I4B), INTENT(IN)                               :: ndim
  INTEGER(I4B), DIMENSION(:), INTENT(INOUT)              :: iii
  INTEGER(I4B), DIMENSION(:), INTENT(INOUT)              :: jjj
  INTEGER(I4B), DIMENSION(:), INTENT(INOUT)              :: kkk
  REAL(DP), DIMENSION(:), INTENT(OUT)                    :: xxx
  INTEGER(I4B), INTENT(OUT)                              :: num
  END SUBROUTINE readtype1
END INTERFACE

!  ****************** END INTERFACE BLOCKS  ********************

! External arrays and variables

INTEGER(I4B), INTENT(IN)                                 :: nx
INTEGER(I4B), INTENT(IN)                                 :: ny
INTEGER(I4B), INTENT(IN)                                 :: nz

LOGICAL(LGT), DIMENSION(:), INTENT(OUT)                  :: cnhIn
LOGICAL(LGT), DIMENSION(:), INTENT(OUT)                  :: riverIn
LOGICAL(LGT), DIMENSION(:), INTENT(OUT)                  :: wellIn

INTEGER(I4B), INTENT(IN)                                 :: ndimdummy
INTEGER(I4B), INTENT(IN)                                 :: mxwell
INTEGER(I4B), INTENT(IN)                                 :: mxrivr
INTEGER(I4B), INTENT(IN)                                 :: mxdrn

!  Internal arrays and variables

INTEGER(I4B)                                             :: ndum

INTEGER(I4B)                                             :: kper
INTEGER(I4B)                                             :: kstp
INTEGER(I4B)                                             :: ncol
INTEGER(I4B)                                             :: nrow
INTEGER(I4B)                                             :: nlay
INTEGER(I4B)                                             :: kperold
INTEGER(I4B)                                             :: kstpold
INTEGER(I4B)                                             :: i          

CHARACTER (LEN=16)                                       :: text

REAL(DP), DIMENSION(nx,ny,nz)                           :: xdum
REAL(DP), DIMENSION(nx,ny)                              :: xydum
INTEGER(I4B), DIMENSION(nx,ny)                           :: lxydum

INTEGER(I4B)                                             :: ncnhmax
INTEGER(I4B)                                             :: nwellsmax
INTEGER(I4B)                                             :: nriversmax
INTEGER(I4B)                                             :: ndrainsmax


WellIn  = .FALSE.
CnhIn  = .FALSE.
RiverIn = .FALSE.

! Start reading thickness, X,Y,Z flows, constant head,
!   and source/sink information

10 READ(1)
READ(1,END = 1000) kper,kstp,ncol,nrow,nlay,text
kperold=kper
kstpold=kstp

DO
  select CASE (text)

!   read thickness
    CASE ("THKSAT")
      CALL readtype3(nx,ny,nz,text,xdum(1:nx,1:ny,1:nz))
!   read X flows
    CASE ("QXX")
      CALL readtype3(nx,ny,nz,text,xdum(1:nx,1:ny,1:nz))
!   read Y flows
    CASE ("QYY")
      CALL readtype3(nx,ny,nz,text,xdum(1:nx,1:ny,1:nz))
!   read Z flows
    CASE ("QZZ")
      CALL readtype3(nx,ny,nz,text,xdum(1:nx,1:ny,1:nz))

!   read constant head locations and rates
    CASE ("CNH")
      CALL readtype1(ndimdummy,jxHeadLoc,jyHeadLoc,jzHeadLoc,qcnh,ncnh)
      IF (ncnh /= ncnhmax .AND. kper /= 1) THEN
        WRITE(*,*) 
        WRITE(*,*) ' There should be a constant number of constant head conditions'
        WRITE(*,*) ' Stress period: ', kper
        WRITE(*,*) ' Number of constant head conditions in current stress period:   ', ncnh
        WRITE(*,*) ' Number of constant head conditions in previous stress periods: ', ncnhmax
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
      ncnhmax = MAX(ncnh,ncnhmax)
      DO i = 1, ncnh
        IF (Qcnh(i) > 1.E-30) THEN
          CnhIn(i)  = .TRUE.
        ELSE
          CONTINUE
        END IF
      END DO

!   read well locations and rates
    CASE ("WEL")
      CALL readtype1(mxwell,jxWellLoc,jyWellLoc,jzWellLoc,q,nwells)
      IF (nwells /= nwellsmax .AND. kper /= 1) THEN
        WRITE(*,*) 
        WRITE(*,*) ' There should be a constant number of wells'
        WRITE(*,*) ' Stress period: ', kper
        WRITE(*,*) ' Number of wells in current stress period:   ', nwells
        WRITE(*,*) ' Number of wells in previous stress periods: ', nwellsmax
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
      nwellsmax = MAX(nwells,nwellsmax)
      DO i = 1, nwells
        IF (Q(i) > 1.E-30) THEN
          WellIn(i) = .TRUE.
        ELSE
          CONTINUE
        END IF
      END DO

!   read drain locations and rates

    CASE ("DRN")
      CALL readtype1(mxdrn,jxDrainLoc,jyDrainLoc,jzDrainLoc,qdrain,ndrains)
      IF (ndrains /= ndrainsmax .AND. kper /= 1) THEN
        WRITE(*,*) 
        WRITE(*,*) ' There should be a constant number of wells'
        WRITE(*,*) ' Stress period: ', kper
        WRITE(*,*) ' Number of drains in current stress period:   ', ndrains
        WRITE(*,*) ' Number of drains in previous stress periods: ', ndrainsmax
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
      ndrainsmax = MAX(ndrains,ndrainsmax)

!   read river locations and rates
    CASE ("RIV")
      CALL readtype1(mxrivr,jxRiverLoc,jyRiverLoc,jzRiverLoc,qriver,nrivers)
      IF (nrivers /= nriversmax .AND. kper /= 1) THEN
        WRITE(*,*) 
        WRITE(*,*) ' There should be a constant number of wells'
        WRITE(*,*) ' Stress period: ', kper
        WRITE(*,*) ' Number of rivers in current stress period:   ', nrivers
        WRITE(*,*) ' Number of rivers in previous stress periods: ', nriversmax
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
      nriversmax = MAX(nrivers,nriversmax)
      DO i=1,nrivers
        IF (Qriver(i) > 1.E-30) THEN
          RiverIn(i) = .TRUE.
        ELSE
          CONTINUE
        END IF
      END DO

!   read recharge locations and rates
    CASE ("RCH")
      CALL readtype2(nx,ny,lxydum,xydum)
!   read evapotranspiration locations and rates
    CASE ("EVT")
      CALL readtype2(nx,ny,lxydum,xydum)
!   read head dependent boundary and rates
    CASE ("GHB")
      CALL readtypescan1(ndum)
    CASE default
!     GO TO 1000
  END select
  READ(1,END = 1000) kper,kstp,ncol,nrow,nlay,text

  IF(kper /= kperold.OR.kstp /= kstpold)THEN

    nwells = nwellsmax
    nrivers = nriversmax
    ncnh = ncnhmax
    ndrains = ndrainsmax

    BACKSPACE 1
    BACKSPACE 1
    GO TO 10
  END IF

END DO

1000 CONTINUE
nwells = nwellsmax
nrivers = nriversmax
ncnh = ncnhmax
ndrains = ndrainsmax


RETURN
END SUBROUTINE ModScan





