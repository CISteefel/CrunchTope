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
    
subroutine BreakthroughWrite(ncomp,nspec,nkin,nrct,ngas,npot, &
    nx,ny,nz,nseries,nexchange,nexch_sec,nsurf,nsurf_sec,nn,ikpH,time,nplotsurface,nplotexchange )


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
INTEGER(I4B), INTENT(IN)                                     :: nn
INTEGER(I4B), INTENT(IN)                                     :: ikpH

REAL(DP), INTENT(IN)                                         :: time
INTEGER(I4B), INTENT(IN)                                    :: nplotsurface
INTEGER(I4B), INTENT(IN)                                    :: nplotexchange


!! INTERNAL VARIABLES

INTEGER(I4B)                                                  :: i
INTEGER(I4B)                                                  :: k
INTEGER(I4B)                                                  :: ksp
INTEGER(I4B)                                                  :: ik
INTEGER(I4B)                                                  :: ls
INTEGER(I4B)                                                  :: ns
INTEGER(I4B)                                                  :: is
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

REAL(DP)                                                      :: PrintTime
REAL(DP)                                                      :: AqueousToBulk
REAL(DP)                                                      :: pHwrite
REAL(DP)                                                      :: tK
REAL(DP)                                                      :: qxsum

REAL(DP), PARAMETER                                           :: eps=1.D-12


!!  Write to file at fixed intervals
  
  IF (nplot > 0) THEN

    IF (interval == 0) THEN              !  Specific output times specified 
     IF (OutputTimeCounter <= NumOutputTimes) THEN
      IF (time-eps < OutputTime(OutputTimeCounter) .AND.  time+eps > OutputTime(OutputTimeCounter)  ) THEN
        OutputTimeCounter = OutputTimeCounter + 1
        PrintTime = OutputTimeScale*time
        DO ll = 1,nseries
          intfile = 100+ll
          jx = jxseries(ll)
          jy = jyseries(ll)
          jz = jzseries(ll)
          CALL AqueousToBulkConvert(jx,jy,jz,AqueousToBulk)
          IF (AqueousToBulk == 0.0) THEN    !  Dry cell
            AqueousToBulk = 1.0               !  So, don't use it
          END IF
          IF (nsurf > 0) THEN
            CALL SurfaceConcentration(jx,jy,jz,ncomp,nsurf,nsurf_sec)
          END IF
          IF (nexchange > 0) THEN
            CALL ExchangeConcentration(jx,jy,jz,ncomp,nexchange,nexch_sec)
          END IF
          IF (ikph /= 0) THEN
            phwrite = -(sp(ikph,jx,jy,jz)+gam(ikph,jx,jy,jz) )/clg
          ELSE
            phwrite = 0.0
          END IF
          IF (tecplot) THEN
            IF (iplotall == 1) THEN
              IF (ikph /= 0) THEN            
                WRITE(intfile,705) PrintTime,phwrite,(s(iplot(i),jx,jy,jz),i=1,nplot),   &
                  (SurfaceCon(i),i=1,nplotsurface),(ExchangeCon(i),i=1,nplotexchange),   &
                  (sp(ik,jx,jy,jz)/clg,ik=1,ncomp+nspec),                                &
                  (spsurf10(is,jx,jy,jz)/AqueousToBulk,is=1,nsurf+nsurf_sec),            &
                  (spex10(nex+nexchange,jx,jy,jz)/AqueousToBulk,nex=1,nexch_sec)  
              ELSE
                WRITE(intfile,705) PrintTime,(s(iplot(i),jx,jy,jz),i=1,nplot),         &
                  (SurfaceCon(i),i=1,nplotsurface),(ExchangeCon(i),i=1,nplotexchange),   &
                  (sp(ik,jx,jy,jz)/clg,ik=1,ncomp+nspec),                                &
                  (spsurf10(is,jx,jy,jz)/AqueousToBulk,is=1,nsurf+nsurf_sec),            &
                  (spex10(nex+nexchange,jx,jy,jz)/AqueousToBulk,nex=1,nexch_sec)  
              END IF
            ELSE
              IF (iplotph /= 0) THEN
                IF (iplotph == 1 .AND. nplot == 1) THEN
                  WRITE(intfile,185) PrintTime,phwrite
                ELSE IF (iplotph == 1 .AND. nplot > 1) THEN
                  WRITE(intfile,185) PrintTime,phwrite,(s(iplot(i),jx,jy,jz),i=2,nplot)  !!& !added next line, hardwired sergi - 2011-08-25
                                           !! ,(volfx(k,jx,jy,jz),k=1,nkin)
                ELSE IF (iplotph == nplot) THEN
                  WRITE(intfile,185) PrintTime,(s(iplot(i),jx,jy,jz),i=1,nplot-1),phwrite
                ELSE
                  WRITE(intfile,185) PrintTime,(s(iplot(i),jx,jy,jz),i=1,iplotph-1),  &
                    phwrite,(s(iplot(i),jx,jy,jz),i=iplotph+1,nplot)
                END IF
              ELSE
                WRITE(intfile,185) PrintTime,(s(iplot(i),jx,jy,jz),i=1,nplot)  & !added next line, hardwired sergi - 2011-08-25
                                            ,(volfx(k,jx,jy,jz),k=1,nkin)
              END IF
            END IF
          ELSE

            IF (iplotall == 1) THEN
              IF (ikph /= 0) THEN            
                WRITE(intfile,705) PrintTime,phwrite,(s(iplot(i),jx,jy,jz),i=1,nplot),   &
                  (SurfaceCon(i),i=1,nplotsurface),(ExchangeCon(i),i=1,nplotexchange),   &
                  (sp(ik,jx,jy,jz)/clg,ik=1,ncomp+nspec),                                &
                  (spsurf10(is,jx,jy,jz)/AqueousToBulk,is=1,nsurf+nsurf_sec),            &
                  (spex10(nex+nexchange,jx,jy,jz)/AqueousToBulk,nex=1,nexch_sec)  
              ELSE
                WRITE(intfile,705) PrintTime,(s(iplot(i),jx,jy,jz),i=1,nplot),   &
                  (SurfaceCon(i),i=1,nplotsurface),(ExchangeCon(i),i=1,nplotexchange),   &
                  (sp(ik,jx,jy,jz)/clg,ik=1,ncomp+nspec),                                &
                  (spsurf10(is,jx,jy,jz)/AqueousToBulk,is=1,nsurf+nsurf_sec),            &
                  (spex10(nex+nexchange,jx,jy,jz)/AqueousToBulk,nex=1,nexch_sec)  
              END IF
            ELSE
              IF (iplotph /= 0) THEN
                IF (iplotph == 1 .AND. nplot == 1) THEN
                  WRITE(intfile,705) PrintTime,phwrite
                ELSE IF (iplotph == 1 .AND. nplot > 1) THEN
                  WRITE(intfile,705) PrintTime,phwrite,(s(iplot(i),jx,jy,jz),i=2,nplot)
                ELSE IF (iplotph == nplot) THEN
                  WRITE(intfile,705) PrintTime,(s(iplot(i),jx,jy,jz),i=1,nplot-1),phwrite
                ELSE
                  WRITE(intfile,705) PrintTime,(s(iplot(i),jx,jy,jz),i=1,iplotph-1),  &
                    phwrite,(s(iplot(i),jx,jy,jz),i=iplotph+1,nplot)
                END IF
              ELSE
                WRITE(intfile,705) PrintTime,(s(iplot(i),jx,jy,jz)/PestScale(i),i=1,nplot)
              END IF
            END IF
          END IF
        END DO
      END IF
     END IF

    ELSE      !! End of write for specific output times (as opposed to interval)

      IF (MOD(nn,interval) == 0) THEN

        PrintTime = OutputTimeScale*time
        DO ll = 1,nseries
          intfile = 100+ll
          jx = jxseries(ll)
          jy = jyseries(ll)
          jz = jzseries(ll)
          CALL AqueousToBulkConvert(jx,jy,jz,AqueousToBulk)
          IF (AqueousToBulk == 0.0) THEN    !  Dry cell
            AqueousToBulk = 1.0               !  So, don't use it
          END IF
          IF (nsurf > 0) THEN
            CALL SurfaceConcentration(jx,jy,jz,ncomp,nsurf,nsurf_sec)
          END IF
          IF (nexchange > 0) THEN
            CALL ExchangeConcentration(jx,jy,jz,ncomp,nexchange,nexch_sec)
          END IF

          IF (ikph /= 0) THEN
            phwrite = -(sp(ikph,jx,jy,jz)+gam(ikph,jx,jy,jz) )/clg
          ELSE
            phwrite = 0.0
          END IF

          IF (tecplot) THEN
            IF (iplotall == 1) THEN
              IF (ikph /= 0) THEN            
                WRITE(intfile,705) PrintTime,phwrite,(s(iplot(i),jx,jy,jz),i=1,nplot),  &
                  (SurfaceCon(i),i=1,nplotsurface),(ExchangeCon(i),i=1,nplotexchange),   &
                  (sp(ik,jx,jy,jz)/clg,ik=1,ncomp+nspec),                                &
                  (spsurf10(is,jx,jy,jz)/AqueousToBulk,is=1,nsurf+nsurf_sec),            &
                  (spex10(nex+nexchange,jx,jy,jz)/AqueousToBulk,nex=1,nexch_sec)  
              ELSE
                WRITE(intfile,705) PrintTime,(s(iplot(i),jx,jy,jz),i=1,nplot),  &
                  (SurfaceCon(i),i=1,nplotsurface),(ExchangeCon(i),i=1,nplotexchange),   &
                  (sp(ik,jx,jy,jz)/clg,ik=1,ncomp+nspec),                                &
                  (spsurf10(is,jx,jy,jz)/AqueousToBulk,is=1,nsurf+nsurf_sec),            &
                  (spex10(nex+nexchange,jx,jy,jz)/AqueousToBulk,nex=1,nexch_sec)  
              END IF
            ELSE
              IF (iplotph /= 0) THEN
                IF (iplotph == 1 .AND. nplot == 1) THEN
                  WRITE(intfile,185) PrintTime,phwrite
                ELSE IF (iplotph == 1 .AND. nplot > 1) THEN

                  WRITE(intfile,705) PrintTime,phwrite,(s(iplot(i),jx,jy,jz),i=2,nplot)
   
                ELSE IF (iplotph == nplot) THEN

                  qxsum = 0.0
                  IF (Benchmark) THEN
                    qxSum = 0.0d0
                    do jy = 1,ny
                      qxSum = qxSum + qx(nx,jy,1)*dyy(jy)
                    end do
                  END IF

                  WRITE(intfile,705) PrintTime,(s(iplot(i),jx,jy,jz),i=1,nplot-1),phwrite,qxsum/365.0

                ELSE
                  WRITE(intfile,185) PrintTime,(s(iplot(i),jx,jy,jz),i=1,iplotph-1),  &
                    phwrite,(s(iplot(i),jx,jy,jz),i=iplotph+1,nplot)
                END IF
              ELSE
                WRITE(intfile,185) PrintTime,(s(iplot(i),jx,jy,jz),i=1,nplot)
              END IF
            END IF

          ELSE    !! Kaleidagraph (and OriginLab?)

            IF (iplotall == 1) THEN
              IF (ikph /= 0) THEN            
                WRITE(intfile,705) PrintTime,phwrite,(s(iplot(i),jx,jy,jz),i=1,nplot),  &
                  (SurfaceCon(i),i=1,nplotsurface),(ExchangeCon(i),i=1,nplotexchange),   &
                  (sp(ik,jx,jy,jz)/clg,ik=1,ncomp+nspec),                                &
                  (spsurf10(is,jx,jy,jz)/AqueousToBulk,is=1,nsurf+nsurf_sec),            &
                  (spex10(nex+nexchange,jx,jy,jz)/AqueousToBulk,nex=1,nexch_sec)  
              ELSE
                WRITE(intfile,705) PrintTime,(s(iplot(i),jx,jy,jz),i=1,nplot),  &
                  (SurfaceCon(i),i=1,nplotsurface),(ExchangeCon(i),i=1,nplotexchange),   &
                  (sp(ik,jx,jy,jz)/clg,ik=1,ncomp+nspec),                                &
                  (spsurf10(is,jx,jy,jz)/AqueousToBulk,is=1,nsurf+nsurf_sec),            &
                  (spex10(nex+nexchange,jx,jy,jz)/AqueousToBulk,nex=1,nexch_sec)  
              END IF
            ELSE
              IF (iplotph /= 0) THEN
                IF (iplotph == 1 .AND. nplot == 1) THEN
                  WRITE(intfile,705) PrintTime,phwrite
                ELSE IF (iplotph == 1 .AND. nplot > 1) THEN
                  Tk = t(jx,jy,jz) + 273.15
                  WRITE(intfile,705) PrintTime,phwrite,(s(iplot(i),jx,jy,jz),i=2,nplot)
              ELSE IF (iplotph == nplot) THEN
                  Tk = 25.0 + 273.15
                  WRITE(intfile,705) PrintTime,(s(iplot(i),jx,jy,jz),i=1,nplot-1),phwrite

                ELSE
                  WRITE(intfile,705) PrintTime,(s(iplot(i),jx,jy,jz),i=1,iplotph-1),  &
                    phwrite,(s(iplot(i),jx,jy,jz),i=iplotph+1,nplot)
                END IF
              ELSE
                WRITE(intfile,705) PrintTime,(s(iplot(i),jx,jy,jz),i=1,nplot)
              END IF
            END IF
          END IF
        END DO
      END IF

    END IF
  END IF

705 FORMAT(1X,1PE12.5,13x,150(1X,1PE22.14))
185 FORMAT(1X,1PE12.5,13x,100(1X,1PE22.14))


RETURN

END subroutine BreakthroughWrite