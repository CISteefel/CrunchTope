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
    
subroutine BreakthroughInitialize(ncomp,nspec,nkin,nrct,ngas,npot, &
    nx,ny,nz,nseries,nexchange,nexch_sec,nsurf,nsurf_sec,ikpH,nplotsurface,nplotexchange)

USE crunchtype
USE params
USE runtime
USE concentration
USE mineral
USE solver
USE medium
USE transport
USE flow
USE temperature
USE io
USE strings
USE ReadFlow
USE modflowModule
USE NanoCrystal

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                     :: ncomp
INTEGER(I4B), INTENT(IN)                                     :: nspec
INTEGER(I4B), INTENT(IN)                                     :: nkin
INTEGER(I4B), INTENT(IN)                                     :: nrct
INTEGER(I4B), INTENT(IN)                                     :: ngas
INTEGER(I4B), INTENT(IN)                                     :: npot
INTEGER(I4B), INTENT(IN)                                     :: nx
INTEGER(I4B), INTENT(IN)                                     :: ny
INTEGER(I4B), INTENT(IN)                                     :: nz
INTEGER(I4B), INTENT(IN)                                     :: nseries
INTEGER(I4B), INTENT(IN)                                     :: nexchange
INTEGER(I4B), INTENT(IN)                                     :: nexch_sec
INTEGER(I4B), INTENT(IN)                                     :: nsurf
INTEGER(I4B), INTENT(IN)                                     :: nsurf_sec
INTEGER(I4B), INTENT(IN)                                     :: ikpH
INTEGER(I4B), INTENT(OUT)                                    :: nplotsurface
INTEGER(I4B), INTENT(OUT)                                    :: nplotexchange


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
INTEGER(I4B)                                                  :: kk
INTEGER(I4B)                                                  :: intfile

CHARACTER (LEN=mls)                                           :: filename
CHARACTER (LEN=mls)                                           :: dummy
CHARACTER (LEN=mls)                                           :: dumstring
CHARACTER (LEN=12)                                            :: writeph

LOGICAL(LGT)                                                  :: ext


CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE                :: ExchangeBasis
CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE                :: SurfaceBasis

writeph = 'pH           '

IF (nexchange > 0) THEN
  nplotexchange = ncomp
  ALLOCATE(ExchangeBasis(ncomp))
  DO i = 1,ncomp
    Dummy = 'EXCH-'
    dumstring = ulab(i)
    Dummy(6:mls) = dumstring(1:mls-5)
    ExchangeBasis(i) = dummy
  END DO
ELSE
  nplotexchange = 0
END IF

IF (nsurf > 0) THEN
  nplotsurface = ncomp
  ALLOCATE(SurfaceBasis(ncomp))
  DO i = 1,ncomp
    Dummy = 'SURF-'
    dumstring = ulab(i)
    Dummy(6:mls) = dumstring(1:mls-5)
    SurfaceBasis(i) = dummy
  END DO
ELSE
  nplotsurface = 0
END IF

IF (nseries >= 1) THEN
  DO ll = 1,nseries
    intfile = 100+ll  
    IF (irestart == 1 .AND. AppendRestart) THEN    !  Open breakthrough file and go to end of file
      filename = TimeSeriesFile(ll)
      INQUIRE(FILE=filename,EXIST=ext)
      IF (.NOT. ext) THEN
        CALL stringlen(filename,ls)
        WRITE(*,*) 
        WRITE(*,*) ' Cannot find time series file: ', filename(1:ls)
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
      OPEN(UNIT=intfile,FILE=filename,STATUS='unknown',ERR=702,POSITION='append')
    ELSE
      filename = TimeSeriesFile(ll)
      OPEN(UNIT=intfile,FILE=filename,STATUS='unknown',ERR=702) 
!!CIS: Moved to individual plotting options      WRITE(intfile,2283) jxseries(ll),jyseries(ll),jzseries(ll)
2283 FORMAT('# Time series at grid cell: 'i3,1x,i3,1x,i3)      

      IF (OutputTimeUnits == 'years') THEN

        IF (tecplot) THEN
          WRITE(intfile,2283) jxseries(ll),jyseries(ll),jzseries(ll)
          IF (iplotall == 1) THEN
            IF (ikph /= 0) THEN
              WRITE(intfile,3001) writeph,(ulab(iplot(i)),i=1,nplot),        &
                (SurfaceBasis(i),i=1,nplotsurface),                          &
                (ExchangeBasis(i),i=1,nplotexchange),                        &
                (ulab(ik),ik=1,ncomp+nspec),                                 &
                (namsurf(ns),ns=1,nsurf),(namsurf_sec(ns),ns=1,nsurf_sec),   &
                (nam_exchsec(nex),nex=1,nexch_sec),                          &
                (namg(kk),kk=1,ngas)
            ELSE
              WRITE(intfile,3001) (ulab(iplot(i)),i=1,nplot),                &
                (SurfaceBasis(i),i=1,nplotsurface),                          &
                (ExchangeBasis(i),i=1,nplotexchange),                        &
                (ulab(ik),ik=1,ncomp+nspec),                                 &
                (namsurf(ns),ns=1,nsurf),(namsurf_sec(ns),ns=1,nsurf_sec),   &
                (nam_exchsec(nex),nex=1,nexch_sec),                          &
                (namg(kk),kk=1,ngas)
            END IF
          ELSE
            WRITE(intfile,3001) (TimeSeriesSpecies(i),i=1,nplot) 
          END IF
        ELSE IF (originlab) THEN
          WRITE(intfile,3006) (TimeSeriesSpecies(i),i=1,nplot) 
          WRITE(intfile,3007) (TimeSeriesUnits(i),i=1,nplot)
        ELSE
          WRITE(intfile,2283) jxseries(ll),jyseries(ll),jzseries(ll)
          IF (iplotall == 1) THEN
            IF (ikph /= 0) THEN
              WRITE(intfile,3701) writeph,(ulab(iplot(i)),i=1,nplot),        &
                (SurfaceBasis(i),i=1,nplotsurface),                          &
                (ExchangeBasis(i),i=1,nplotexchange),                        &
                (ulab(ik),ik=1,ncomp+nspec),                                 &
                (namsurf(ns),ns=1,nsurf),(namsurf_sec(ns),ns=1,nsurf_sec),   &
                (nam_exchsec(nex),nex=1,nexch_sec),                          &
                (namg(kk),kk=1,ngas)
            ELSE
              WRITE(intfile,3701) (ulab(iplot(i)),i=1,nplot),                &
                (SurfaceBasis(i),i=1,nplotsurface),                          &
                (ExchangeBasis(i),i=1,nplotexchange),                        &
                (ulab(ik),ik=1,ncomp+nspec),                                 &
                (namsurf(ns),ns=1,nsurf),(namsurf_sec(ns),ns=1,nsurf_sec),   &
                (nam_exchsec(nex),nex=1,nexch_sec),                          &
                (namg(kk),kk=1,ngas)
            END IF
          ELSE
            WRITE(intfile,3701) (TimeSeriesSpecies(i),i=1,nplot) 
          END IF
        ENDIF
        
      ELSE IF (OutputTimeUnits == 'days') THEN

        IF (tecplot) THEN
          WRITE(intfile,2283) jxseries(ll),jyseries(ll),jzseries(ll)

          IF (iplotall == 1) THEN
            IF (ikph /= 0) THEN
              WRITE(intfile,3002) writeph,(ulab(iplot(i)),i=1,nplot),        &
                (SurfaceBasis(i),i=1,nplotsurface),                          &
                (ExchangeBasis(i),i=1,nplotexchange),                        &
                (ulab(ik),ik=1,ncomp+nspec),                                 &
                (namsurf(ns),ns=1,nsurf),(namsurf_sec(ns),ns=1,nsurf_sec),   &
                (nam_exchsec(nex),nex=1,nexch_sec),                          &
                (namg(kk),kk=1,ngas)
            ELSE
              WRITE(intfile,3002) (ulab(iplot(i)),i=1,nplot),                &
                (SurfaceBasis(i),i=1,nplotsurface),                          &
                (ExchangeBasis(i),i=1,nplotexchange),                        &
                (ulab(ik),ik=1,ncomp+nspec),                                 &
                (namsurf(ns),ns=1,nsurf),(namsurf_sec(ns),ns=1,nsurf_sec),   &
                (nam_exchsec(nex),nex=1,nexch_sec),                          &
                (namg(kk),kk=1,ngas)
            END IF

          ELSE
            WRITE(intfile,3002) (TimeSeriesSpecies(i),i=1,nplot) 
          END IF

        ELSE IF (originlab) THEN

          IF (JennyDruhan) THEN
            WRITE(intfile,3036) (TimeSeriesSpecies(i),i=1,nplot) 
            WRITE(intfile,3008) (TimeSeriesUnits(i),i=1,nplot) 
          ELSE
            WRITE(intfile,3006) (TimeSeriesSpecies(i),i=1,nplot) 
            WRITE(intfile,3008) (TimeSeriesUnits(i),i=1,nplot) 
          END IF
        ELSE
          WRITE(intfile,2283) jxseries(ll),jyseries(ll),jzseries(ll)

          IF (iplotall == 1) THEN
            IF (ikph /= 0) THEN
              WRITE(intfile,3702) writeph,(ulab(iplot(i)),i=1,nplot),        &
                (SurfaceBasis(i),i=1,nplotsurface),                          &
                (ExchangeBasis(i),i=1,nplotexchange),                        &
                (ulab(ik),ik=1,ncomp+nspec),                                 &
                (namsurf(ns),ns=1,nsurf),(namsurf_sec(ns),ns=1,nsurf_sec),   &
                (nam_exchsec(nex),nex=1,nexch_sec),                          &
                (namg(kk),kk=1,ngas)
            ELSE
              WRITE(intfile,3702) (ulab(iplot(i)),i=1,nplot),                &
                (SurfaceBasis(i),i=1,nplotsurface),                          &
                (ExchangeBasis(i),i=1,nplotexchange),                        &
                (ulab(ik),ik=1,ncomp+nspec),                                 &
                (namsurf(ns),ns=1,nsurf),(namsurf_sec(ns),ns=1,nsurf_sec),   &
                (nam_exchsec(nex),nex=1,nexch_sec),                          &
                (namg(kk),kk=1,ngas)
            END IF    

          ELSE

!! Kaleidagraph Jenny Druhan
            if (JennyDruhan) then
!!              WRITE(intfile,3702) (TimeSeriesSpecies(i),i=1,nplot), 'delCaAqueous' , 'delCaMineralBulk', 'delCaMineralInstant'
              WRITE(intfile,3712) (TimeSeriesSpecies(i),i=1,nplot) 
            else
              WRITE(intfile,3702) (TimeSeriesSpecies(i),i=1,nplot) 
            end if

          END IF
        ENDIF

      ELSE IF (OutputTimeUnits == 'hours') THEN

        IF (tecplot) THEN
          WRITE(intfile,2283) jxseries(ll),jyseries(ll),jzseries(ll)
          IF (iplotall == 1) THEN
            IF (ikph /= 0) THEN
              WRITE(intfile,3003) writeph,(ulab(iplot(i)),i=1,nplot),        &
                (SurfaceBasis(i),i=1,nplotsurface),                          &
                (ExchangeBasis(i),i=1,nplotexchange),                        &
                (ulab(ik),ik=1,ncomp+nspec),                                 &
                (namsurf(ns),ns=1,nsurf),(namsurf_sec(ns),ns=1,nsurf_sec),   &
                (nam_exchsec(nex),nex=1,nexch_sec),                          &
                (namg(kk),kk=1,ngas)
            ELSE
              WRITE(intfile,3003) (ulab(iplot(i)),i=1,nplot),                &
                (SurfaceBasis(i),i=1,nplotsurface),                          &
                (ExchangeBasis(i),i=1,nplotexchange),                        &
                (ulab(ik),ik=1,ncomp+nspec),                                 &
                (namsurf(ns),ns=1,nsurf),(namsurf_sec(ns),ns=1,nsurf_sec),   &
                (nam_exchsec(nex),nex=1,nexch_sec),                          &
                (namg(kk),kk=1,ngas)
            END IF
          ELSE
            WRITE(intfile,3003) (TimeSeriesSpecies(i),i=1,nplot) 
          END IF
        ELSE IF (originlab) THEN
            WRITE(intfile,3006) (TimeSeriesSpecies(i),i=1,nplot) 
            WRITE(intfile,3009) (TimeSeriesUnits(i),i=1,nplot)  
        ELSE
          WRITE(intfile,2283) jxseries(ll),jyseries(ll),jzseries(ll)
          IF (iplotall == 1) THEN
            IF (ikph /= 0) THEN
              WRITE(intfile,3703) writeph,(ulab(iplot(i)),i=1,nplot),        &
                (SurfaceBasis(i),i=1,nplotsurface),                          &
                (ExchangeBasis(i),i=1,nplotexchange),                        &
                (ulab(ik),ik=1,ncomp+nspec),                                 &
                (namsurf(ns),ns=1,nsurf),(namsurf_sec(ns),ns=1,nsurf_sec),   &
                (nam_exchsec(nex),nex=1,nexch_sec),                          &
                (namg(kk),kk=1,ngas)
            ELSE
              WRITE(intfile,3703) (ulab(iplot(i)),i=1,nplot),                &
                (SurfaceBasis(i),i=1,nplotsurface),                          &
                (ExchangeBasis(i),i=1,nplotexchange),                        &
                (ulab(ik),ik=1,ncomp+nspec),                                 &
                (namsurf(ns),ns=1,nsurf),(namsurf_sec(ns),ns=1,nsurf_sec),   &
                (nam_exchsec(nex),nex=1,nexch_sec),                          &
                (namg(kk),kk=1,ngas)
            END IF
          ELSE
            WRITE(intfile,3703) (TimeSeriesSpecies(i),i=1,nplot) 
          END IF
        ENDIF

      ELSE IF (OutputTimeUnits == 'minutes') THEN
 
        IF (tecplot) THEN
          WRITE(intfile,2283) jxseries(ll),jyseries(ll),jzseries(ll)
          IF (iplotall == 1) THEN
            IF (ikph /= 0) THEN
              WRITE(intfile,3004) writeph,(ulab(iplot(i)),i=1,nplot),        &
                (SurfaceBasis(i),i=1,nplotsurface),                          &
                (ExchangeBasis(i),i=1,nplotexchange),                        &
                (ulab(ik),ik=1,ncomp+nspec),                                 &
                (namsurf(ns),ns=1,nsurf),(namsurf_sec(ns),ns=1,nsurf_sec),   &
                (nam_exchsec(nex),nex=1,nexch_sec),                          &
                (namg(kk),kk=1,ngas)
            ELSE
              WRITE(intfile,3004) (ulab(iplot(i)),i=1,nplot),                &
               (ulab(i),i=1,nplotsurface),(ulab(i),i=1,nplotexchange),       &
                (ulab(ik),ik=1,ncomp+nspec),                                 &
                (namsurf(ns),ns=1,nsurf),(namsurf_sec(ns),ns=1,nsurf_sec),   &
                (nam_exchsec(nex),nex=1,nexch_sec),                          &
                (namg(kk),kk=1,ngas)
            END IF
          ELSE
            WRITE(intfile,3004) (TimeSeriesSpecies(i),i=1,nplot) 
          END IF
        ELSE IF (originlab) THEN
          WRITE(intfile,3006) (TimeSeriesSpecies(i),i=1,nplot) 
          WRITE(intfile,3010) (TimeSeriesUnits(i),i=1,nplot)  
        ELSE
          WRITE(intfile,2283) jxseries(ll),jyseries(ll),jzseries(ll)
          IF (iplotall == 1) THEN
            IF (ikph /= 0) THEN
              WRITE(intfile,3704) writeph,(ulab(iplot(i)),i=1,nplot),        &
                (SurfaceBasis(i),i=1,nplotsurface),                          &
                (ExchangeBasis(i),i=1,nplotexchange),                        &
                (ulab(ik),ik=1,ncomp+nspec),                                 &
                (namsurf(ns),ns=1,nsurf),(namsurf_sec(ns),ns=1,nsurf_sec),   &
                (nam_exchsec(nex),nex=1,nexch_sec),                          &
                (namg(kk),kk=1,ngas)
            ELSE
              WRITE(intfile,3704) (ulab(iplot(i)),i=1,nplot),                &
                (SurfaceBasis(i),i=1,nplotsurface),                          &
                (ExchangeBasis(i),i=1,nplotexchange),                        &
                (ulab(ik),ik=1,ncomp+nspec),                                 &
                (namsurf(ns),ns=1,nsurf),(namsurf_sec(ns),ns=1,nsurf_sec),   &
                (nam_exchsec(nex),nex=1,nexch_sec),                          &
                (namg(kk),kk=1,ngas)
            END IF
          ELSE
            WRITE(intfile,3704) (TimeSeriesSpecies(i),i=1,nplot) 
          END IF
        ENDIF

      ELSE IF (OutputTimeUnits == 'seconds') THEN
        IF (tecplot) THEN
          WRITE(intfile,2283) jxseries(ll),jyseries(ll),jzseries(ll)
          IF (iplotall == 1) THEN
            IF (ikph /= 0) THEN
              WRITE(intfile,3005) writeph,(ulab(iplot(i)),i=1,nplot),        &
                (SurfaceBasis(i),i=1,nplotsurface),                          &
                (ExchangeBasis(i),i=1,nplotexchange),                        &
                (ulab(ik),ik=1,ncomp+nspec),                                 &
                (namsurf(ns),ns=1,nsurf),(namsurf_sec(ns),ns=1,nsurf_sec),   &
                (nam_exchsec(nex),nex=1,nexch_sec),                          &
                (namg(kk),kk=1,ngas)
            ELSE
              WRITE(intfile,3005) (ulab(iplot(i)),i=1,nplot),                &
               (ulab(i),i=1,nplotsurface),(ulab(i),i=1,nplotexchange),       &
                (ulab(ik),ik=1,ncomp+nspec),                                 &
                (namsurf(ns),ns=1,nsurf),(namsurf_sec(ns),ns=1,nsurf_sec),   &
                (nam_exchsec(nex),nex=1,nexch_sec),                          &
                (namg(kk),kk=1,ngas)
            END IF
          ELSE
            WRITE(intfile,3005) (TimeSeriesSpecies(i),i=1,nplot) 
          END IF
        ELSE IF (originlab) THEN
          WRITE(intfile,3006) (TimeSeriesSpecies(i),i=1,nplot) 
          WRITE(intfile,3011) (TimeSeriesUnits(i),i=1,nplot)  
        ELSE
          WRITE(intfile,2283) jxseries(ll),jyseries(ll),jzseries(ll)
          IF (iplotall == 1) THEN
            IF (ikph /= 0) THEN
              WRITE(intfile,3705) writeph,(ulab(iplot(i)),i=1,nplot),        &
                (SurfaceBasis(i),i=1,nplotsurface),                          &
                (ExchangeBasis(i),i=1,nplotexchange),                        &
                (ulab(ik),ik=1,ncomp+nspec),                                 &
                (namsurf(ns),ns=1,nsurf),(namsurf_sec(ns),ns=1,nsurf_sec),   &
                (nam_exchsec(nex),nex=1,nexch_sec),                          &
                (namg(kk),kk=1,ngas)
            ELSE
              WRITE(intfile,3705) (ulab(iplot(i)),i=1,nplot),                &
                (SurfaceBasis(i),i=1,nplotsurface),                          &
                (ExchangeBasis(i),i=1,nplotexchange),                        &
                (ulab(ik),ik=1,ncomp+nspec),                                 &
                (namsurf(ns),ns=1,nsurf),(namsurf_sec(ns),ns=1,nsurf_sec),   &
                (nam_exchsec(nex),nex=1,nexch_sec),                          &
                (namg(kk),kk=1,ngas)
            END IF
          ELSE
            WRITE(intfile,3705) (TimeSeriesSpecies(i),i=1,nplot) 
          END IF
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

IF (ALLOCATED(SurfaceBasis)) THEN
  DEALLOCATE(SurfaceBasis)
END IF
IF (ALLOCATED(ExchangeBasis)) THEN
  DEALLOCATE(ExchangeBasis)
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

3036 FORMAT('  Time     ', 2x,11(3X,a20),3x,'delCaAqueous        ',3x,     &
                                          'delCaMineral        ',3x,     &
                                          'del34S_Sulfate      ',3x,     &
                                          'del34S_Sulfide      ',3x,     &
                                          'del34S_Mineral      '     )  
3712 FORMAT('  Time(day)',2x,11(3X,a20),3x,'delCaAqueous        ',3x,     &
                                          'delCaMineral        ',3x,     &
                                          'del34S_Sulfate      ',3x,     &
                                          'del34S_Sulfide      ',3x,     &
                                          'del34S_Mineral      '     )
!! Kaleidagraph formats
3701 FORMAT('  Time(yrs)',17x,150(1X,A22))
3702 FORMAT('  Time(day)',17x,150(1X,A22))
3703 FORMAT('  Time(hrs)',17x,150(1X,A22))
3704 FORMAT('  Time(min)',17x,150(1X,A22))
3705 FORMAT('  Time(sec)',17x,150(1X,A22))

RETURN

702 WRITE(*,*)
IF (ALLOCATED(SurfaceBasis)) THEN
  DEALLOCATE(SurfaceBasis)
END IF
IF (ALLOCATED(ExchangeBasis)) THEN
  DEALLOCATE(ExchangeBasis)
END IF
WRITE(*,*) ' Error opening breakthrough file'
WRITE(*,*)
READ(*,*)
STOP


END subroutine BreakthroughInitialize