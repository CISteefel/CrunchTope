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

SUBROUTINE tvd(nx,ny,nz,delt,icomp)

USE crunchtype
USE runtime
USE concentration
USE medium
USE transport
USE temperature
USE flow
USE modflowModule

IMPLICIT NONE

INTEGER(I4B), INTENT(IN)                                               :: nx
INTEGER(I4B), INTENT(IN)                                               :: ny
INTEGER(I4B), INTENT(IN)                                               :: nz
REAL(DP), INTENT(IN)                                                   :: delt
INTEGER(I4B), INTENT(IN)                                               :: icomp

LOGICAL(LGT)                                                           :: TvdOption

REAL(DP) :: height,cfl,phi,tmp,sdum,cnv
REAL(DP) :: dtdx,dtdy,dtdz,area,cellvolume,porsatro,PorSorp,sourceTerm,AccumulationTerm


INTEGER(I4B) :: i,j,k,jx,jy,jz

!   Well Declarations
INTEGER(I4B) :: well

!   Constant Head Declarations
INTEGER(I4B) :: cnh

!   Recharge Declarations

!   Evapotranspiration Declarations

!   River Declarations
INTEGER(I4B) :: river

!   Drain Declarations
INTEGER(I4B) :: drain

INTEGER(I4B) :: npz

TvdOption = .TRUE.

ftvd = 0.0
gtvd = 0.0
htvd = 0.0

IF (modflow) THEN
  cnv = ModFlowCnv
ELSE
  cnv = 1.0
END IF

!   x-direction fluxes

DO jz=0,nz
  DO jy=0,ny
    DO jx=0,nx

      IF(jy /= 0 .AND. jz /= 0 .AND. xflow)THEN
        IF(qx(jx,jy,jz) >= 0)THEN !!! +ve velocity case
          area = dzz(jx,jy,jz) * dyy(jy)
          porsatro = por(jx,jy,jz)*satliq(jx,jy,jz)
          IF (modflow) THEN
            cfl= cnv*(ABS(qx(jx,jy,jz))*delt/dxx(jx) )/(area*porsatro) 
          ELSE
            cfl= cnv*(ABS(qx(jx,jy,jz))*delt/dxx(jx) )/(porsatro) 
          ENDIF
          ftvd(jx,jy,jz)= ctvd(jx+1,jy,jz)- ctvd(jx,jy,jz)
          IF(ABS(ftvd(jx,jy,jz)) > 1.e-30 .AND. TvdOption)THEN
            sdum = (ctvd(jx,jy,jz)-ctvd(jx-1,jy,jz))/ftvd(jx,jy,jz)
            tmp = MIN(2.d0,2.d0*sdum,0.33333333333333333D0*  &
                (2.d0-cfl+(1.d0+cfl)*sdum))
            phi = MAX(0.0D0,tmp)
          ELSE
            phi=0
          END IF
          ftvd(jx,jy,jz)=cnv*ro(jx,jy,jz)*qx(jx,jy,jz)* (ctvd(jx,jy,jz) + 0.5*(1.0-cfl)*phi*ftvd(jx,jy,jz) )
          
! write(*,*)'f+ ',jx,jy,jz,qx(jx,jy,jz),CFL,phi,ftvd(jx,jy,jz)
          
          IF (jx == 0) THEN
            IF (modflow) THEN
              cflowin(icomp) = cflowin(icomp) + delt * ftvd(jx,jy,jz)  !  Flux already in units of mole/yr
            ELSE
              cflowin(icomp) = cflowin(icomp) + area * delt * ftvd(jx,jy,jz)
            END IF
          ELSE IF (jx == nx) THEN
            IF (modflow) THEN
              cflowout(icomp) = cflowout(icomp) + delt * ftvd(jx,jy,jz)  !  Flux already in units of mole/yr
            ELSE
              cflowout(icomp) = cflowout(icomp) + area * delt * ftvd(jx,jy,jz)
            END IF
          END IF
          
          
        ELSE IF (xflow) THEN !!! -ve velocity case
          area = dzz(jx+1,jy,jz) * dyy(jy)
          porsatro = por(jx+1,jy,jz)*satliq(jx+1,jy,jz)
          IF (modflow) THEN
            cfl= cnv*(ABS(qx(jx,jy,jz))*delt/dxx(jx+1) )/(area*porsatro) 
          ELSE
            cfl= cnv*(ABS(qx(jx,jy,jz))*delt/dxx(jx+1) )/(porsatro) 
          ENDIF
          ftvd(jx,jy,jz)=ctvd(jx,jy,jz)-ctvd(jx+1,jy,jz)
          IF(ABS(ftvd(jx,jy,jz)) > 1.e-30 .AND. TvdOption)THEN
            sdum=(ctvd(jx+1,jy,jz)-ctvd(jx+2,jy,jz))/ftvd(jx,jy,jz)
            phi=MAX(0.d0,MIN(2.d0,2.d0*sdum, 0.33333333333333333D0*  &
                (2.d0-cfl+(1.d0+cfl)*sdum)))
          ELSE
            phi=0
          END IF
          ftvd(jx,jy,jz)= cnv*ro(jx,jy,jz)*qx(jx,jy,jz)*( ctvd(jx+1,jy,jz) +  &
              0.5*(1.0-cfl)*phi*ftvd(jx,jy,jz) )
                
          IF (jx == 0) THEN
            IF (modflow) THEN
              cflowout(icomp) = cflowout(icomp) + delt * DABS(ftvd(jx,jy,jz))  !  Flux already in units of mole/yr
            ELSE
              cflowout(icomp) = cflowout(icomp) + area * delt * DABS(ftvd(jx,jy,jz))
            END IF
          ELSE IF (jx == nx) THEN
            IF (modflow) THEN
              cflowin(icomp) = cflowin(icomp) + delt * DABS(ftvd(jx,jy,jz))  !  Flux already in units of mole/yr
            ELSE
              cflowin(icomp) = cflowin(icomp) + area * delt * DABS(ftvd(jx,jy,jz))
            END IF
          END IF
          
! write(*,*)'f- ',qx(jx,jy,jz),CFL,phi,ftvd(jx,jy,jz)
        ELSE
          CONTINUE         
        END IF
      END IF
      
!   y-direction fluxes
      
      IF(jx /= 0 .AND. jz /= 0 .AND. yflow)THEN
        IF(qy(jx,jy,jz) >= 0)THEN !!! +ve velocity case
          area = dzz(jx,jy,jz) * dxx(jx)
          porsatro = por(jx,jy,jz)*satliq(jx,jy,jz)
          IF (modflow) THEN
            cfl= cnv*(ABS(qy(jx,jy,jz))*delt/dyy(jy) )/(area*porsatro) 
          ELSE
            cfl= cnv*(ABS(qy(jx,jy,jz))*delt/dyy(jy) )/(porsatro) 
          ENDIF
          gtvd(jx,jy,jz)=ctvd(jx,jy+1,jz)-ctvd(jx,jy,jz)
          IF(ABS(gtvd(jx,jy,jz)) > 1.e-30 .AND. TvdOption)THEN
            sdum=(ctvd(jx,jy,jz)-ctvd(jx,jy-1,jz))/gtvd(jx,jy,jz)
            phi=MAX(0.d0,MIN(2.d0,2.d0*sdum, 0.33333333333333333D0*  &
                (2.d0-cfl+(1.d0+cfl)*sdum)))
          ELSE
            phi=0
          END IF
          gtvd(jx,jy,jz)= cnv*ro(jx,jy,jz)*qy(jx,jy,jz)*( ctvd(jx,jy,jz) +  &
              0.5*(1.0-cfl)*phi*gtvd(jx,jy,jz) )
! write(*,*)'g+ ',qy(jx,jy,jz),CFL,phi,gtvd(jx,jy,jz)
          
        ELSE IF (yflow) THEN    !!! -ve velocity case
          area = dxx(jx) * dzz(jx,jy+1,jz)
          porsatro = por(jx,jy-1,jz)*satliq(jx,jy-1,jz)
          IF (modflow) THEN
            cfl= cnv*(ABS(qy(jx,jy,jz))*delt/dyy(jy+1) )/(area*porsatro) 
          ELSE
            cfl= cnv*(ABS(qy(jx,jy,jz))*delt/dyy(jy+1) )/(porsatro) 
          ENDIF
          gtvd(jx,jy,jz)=ctvd(jx,jy,jz)-ctvd(jx,jy+1,jz)
          IF(ABS(gtvd(jx,jy,jz)) > 1.e-30 .AND. TvdOption)THEN
            sdum= (ctvd(jx,jy+1,jz)-ctvd(jx,jy+2,jz))/gtvd(jx,jy,jz)
            phi=MAX(0.d0,MIN(2.d0,2.d0*sdum, 0.33333333333333333D0*  &
                (2.d0-cfl+(1.d0+cfl)*sdum)))
          ELSE
            phi=0
          END IF
          gtvd(jx,jy,jz)= cnv*ro(jx,jy,jz)*qy(jx,jy,jz)*( ctvd(jx,jy+1,jz) +  &
              0.5*(1.0-cfl)*phi*gtvd(jx,jy,jz) )
! write(*,*)'g- ',qy(jx,jy,jz),CFL,phi,gtvd(jx,jy,jz)

        ELSE
          CONTINUE
        END IF
      END IF
      
!   z-direction fluxes
      
      IF(jx /= 0 .AND. jy /= 0 .AND. zflow)THEN
        height = 0.50D0*(dzz(jx,jy,jz)+dzz(jx,jy,jz+1))
        IF(qz(jx,jy,jz) >= 0)THEN !!! +ve velocitz case
          area = dxx(jx) * dyy(jy)
          porsatro = por(jx,jy,jz)*satliq(jx,jy,jz)
          IF (modflow) THEN
            cfl= cnv*(ABS(qz(jx,jy,jz))*delt/dzz(jx,jy,jz) )/(area*porsatro) 
          ELSE
            cfl= cnv*(ABS(qz(jx,jy,jz))*delt/dzz(jx,jy,jz) )/(porsatro) 
          ENDIF
! write(*,*)'h+ ',jx,jy,jz,dzz(jx,jy,JZ),CFL
          htvd(jx,jy,jz)=ctvd(jx,jy,jz+1)-ctvd(jx,jy,jz)
          IF(ABS(htvd(jx,jy,jz)) > 1.e-30 .AND. TvdOption)THEN
            sdum= (ctvd(jx,jy,jz)-ctvd(jx,jy,jz-1))/htvd(jx,jy,jz)
            phi=MAX(0.d0,MIN(2.d0,2.d0*sdum, 0.33333333333333333D0*  &
                (2.d0-cfl+(1.d0+cfl)*sdum)))
          ELSE
            phi=0
          END IF
          htvd(jx,jy,jz)= cnv*ro(jx,jy,jz)*qz(jx,jy,jz)*( ctvd(jx,jy,jz) +  &
              0.5*(1.0-cfl)*phi*htvd(jx,jy,jz) )
! write(*,*)'h+ ',jx,jy,jz,qz(jx,jy,jz),CFL,phi,htvd(jx,jy,jz)
          
        ELSE IF (zflow) THEN     !!! -ve velocity case
          height = 0.50D0*(dzz(jx,jy,jz)+dzz(jx,jy,jz+1))
          area = dxx(jx) * dyy(jy)
          porsatro = por(jx,jy,jz-1)*satliq(jx,jy,jz-1)
          IF (modflow) THEN
            cfl= cnv*(ABS(qz(jx,jy,jz))*delt/dzz(jx,jy,jz+1) )/(area*porsatro) 
          ELSE
            cfl= cnv*(ABS(qz(jx,jy,jz))*delt/dzz(jx,jy,jz+1) )/(porsatro) 
          ENDIF
          htvd(jx,jy,jz)=ctvd(jx,jy,jz)-ctvd(jx,jy,jz+1)
          IF(ABS(htvd(jx,jy,jz)) > 1.e-30 .AND. TvdOption)THEN
            sdum= (ctvd(jx,jy,jz+1)-ctvd(jx,jy,jz+2))/htvd(jx,jy,jz)
            phi=MAX(0.d0,MIN(2.d0,2.d0*sdum, 0.33333333333333333D0*  &
                (2.d0-cfl+(1.d0+cfl)*sdum)))
          ELSE
            phi=0
          END IF
          htvd(jx,jy,jz)= cnv*ro(jx,jy,jz)*qz(jx,jy,jz)*( ctvd(jx,jy,jz+1) +  &
              0.5*(1.0-cfl)*phi*htvd(jx,jy,jz) )
! write(*,*)'h- ',jx,jy,jz,qz(jx,jy,jz),CFL,phi,htvd(jx,jy,jz)

        ELSE
          CONTINUE
        END IF
      END IF
    END DO
  END DO
END DO

!  Accumulation terms 

DO jz=1,nz
  DO jy=1,ny
    DO jx=1,nx
      cellvolume = dzz(jx,jy,jz)*dxx(jx)*dyy(jy)
      dtdx = delt/dxx(jx)
      dtdy = delt/dyy(jy)
      dtdz = delt/dzz(jx,jy,jz)
      PorSorp = por(jx,jy,jz)*sorp(jx,jy,jz)
      IF (activecell(jx,jy,jz) /= 0) THEN

        sourceTerm = 0.d0
        IF (wells) THEN
          
          DO npz = 1,npump(jx,jy,jz)
            IF (qg(npz,jx,jy,jz) > 0.0) THEN
              sourceTerm = sourceTerm + ro(jx,jy,jz)*qg(npz,jx,jy,jz)*scond(icomp,intbnd(npz,jx,jy,jz))*delt/cellvolume
            ELSE IF (qg(npz,jx,jy,jz) < 0.0d0) THEN
              sourceTerm = sourceTerm + ro(jx,jy,jz)*qg(npz,jx,jy,jz)*ctvd(jx,jy,jz)*delt/cellvolume
            ELSE
              CONTINUE
            END IF      
          END DO
        ELSEIF (pumptimeseries) THEN
          IF (qg(1,jx,jy,jz) > 0.0) THEN
            sourceTerm = sourceTerm + ro(jx,jy,jz)*qg(1,jx,jy,jz)*scond(icomp,intbnd(1,jx,jy,jz))*delt/cellvolume
          ELSEIF (qg(1,jx,jy,jz) < 0.0d0) THEN
            sourceTerm = sourceTerm + ro(jx,jy,jz)*qg(1,jx,jy,jz)*ctvd(jx,jy,jz)*delt/cellvolume
          ELSE
            CONTINUE
          END IF
        END IF
            
        accumulationTerm = satliq(jx,jy,jz)*ro(jx,jy,jz)*PorSorp*ctvd(jx,jy,jz)
        ctvd(jx,jy,jz)= ( sourceTerm + accumulationTerm    &
              - ( dtdx*(ftvd(jx,jy,jz) - ftvd(jx-1,jy,jz)) +                               &
                  dtdy*(gtvd(jx,jy,jz) - gtvd(jx,jy-1,jz)) +                               &
                  dtdz*(htvd(jx,jy,jz) - htvd(jx,jy,jz-1)) )   )                            &
                      /(satliq(jx,jy,jz)*ro(jx,jy,jz)*PorSorp)
          

        
      END IF
            
    END DO
  END DO
END DO

! ********* Do the following only if MODFLOW option is selected *******************************

IF (modflow) THEN

! * Handle the recharge (currently assumes all recharge is source through layer 1)

  IF (RechargeCondition /= 0) THEN
    jz=1
    DO jy=1,ny
      DO jx=1,nx
        cellvolume = dzz(jx,jy,jz)*dxx(jx)*dyy(jy)
        ctvd(jx,jy,jz)= ctvd(jx,jy,jz) + cnv*qrecharge(jx,jy)*scond(icomp,RechargeCondition)*delt/  &
          (sorp(jx,jy,jz)*por(jx,jy,jz)*cellvolume)
        cflowin(icomp) = cflowin(icomp) + cnv*qrecharge(jx,jy) * scond(icomp,RechargeCondition) *delt
      END DO
    END DO
  END IF

!   Handle the evapotranspiration (currently assumes all evaporation is sink through layer 1)

  jz=1
  DO jy=1,ny
    DO jx=1,nx
      cellvolume = dzz(jx,jy,jz)*dxx(jx)*dyy(jy)
      ctvd(jx,jy,jz) = ctvd(jx,jy,jz) + cnv*qevt(jx,jy) * sn(icomp,jx,jy,jz) *delt/  &
        (sorp(jx,jy,jz)*por(jx,jy,jz)*cellvolume)   
      cflowout(icomp) = cflowout(icomp) - cnv*qevt(jx,jy) * sn(icomp,jx,jy,jz) *delt
    END DO
  END DO

! * Handle the wells

  DO well = 1, nwells
    jx = jxWellLoc(well)
    jy = jyWellLoc(well)
    jz = jzWellLoc(well)
    cellvolume = dzz(jx,jy,jz)*dxx(jx)*dyy(jy)
    IF (wtype(well) == 'Injection') THEN
      ctvd(jx,jy,jz) = ctvd(jx,jy,jz) + cnv*q(well) * scond(icomp,WellCondition(well)) *delt/  &
        (sorp(jx,jy,jz)*por(jx,jy,jz)*cellvolume)
      cflowin(icomp) = cflowin(icomp) + cnv*q(well) * scond(icomp,WellCondition(well)) *delt
    ELSE
      ctvd(jx,jy,jz) = ctvd(jx,jy,jz) + cnv*q(well) * sn(icomp,jx,jy,jz) *delt/  &
        (sorp(jx,jy,jz)*por(jx,jy,jz)*cellvolume)  
      cflowout(icomp) = cflowout(icomp) - cnv*q(well) * sn(icomp,jx,jy,jz) *delt   
    END IF
  END DO

!   Handle the constant heads

  DO cnh = 1, ncnh
    jx = jxHeadLoc(cnh)
    jy = jyHeadLoc(cnh)
    jz = jzHeadLoc(cnh)
    cellvolume = dzz(jx,jy,jz)*dxx(jx)*dyy(jy)
    IF (htype(cnh) == 'Injection') THEN
      ctvd(jx,jy,jz) = ctvd(jx,jy,jz) + cnv*qcnh(cnh) * scond(icomp,HeadCondition(cnh)) *delt/  &
        (sorp(jx,jy,jz)*por(jx,jy,jz)*cellvolume)
      cflowin(icomp) = cflowin(icomp) + cnv*qcnh(cnh) * scond(icomp,HeadCondition(cnh)) *delt
    ELSE
      ctvd(jx,jy,jz) = ctvd(jx,jy,jz) + cnv*qcnh(cnh) * sn(icomp,jx,jy,jz) *delt/  &
        (sorp(jx,jy,jz)*por(jx,jy,jz)*cellvolume)
      cflowout(icomp) = cflowout(icomp) - cnv*qcnh(cnh) * sn(icomp,jx,jy,jz) *delt
    END IF
  END DO

! * Handle the river

  DO river = 1, nrivers
    jx = jxRiverLoc(river)
    jy = jyRiverLoc(river)
    jz = jzRiverLoc(river)
    cellvolume = dzz(jx,jy,jz)*dxx(jx)*dyy(jy)
    IF (rtype(river) == 'Injection') THEN
      ctvd(jx,jy,jz) = ctvd(jx,jy,jz) + cnv*q(river) * scond(icomp,RiverCondition(river)) *delt/  &
        (sorp(jx,jy,jz)*por(jx,jy,jz)*cellvolume)
      cflowin(icomp) = cflowin(icomp) + cnv*q(river) * scond(icomp,RiverCondition(river)) *delt
    ELSE
      ctvd(jx,jy,jz) = ctvd(jx,jy,jz) + cnv*qriver(river) * sn(icomp,jx,jy,jz) *delt/  &
        (sorp(jx,jy,jz)*por(jx,jy,jz)*cellvolume)   
      cflowout(icomp) = cflowout(icomp) - cnv*qriver(river) * sn(icomp,jx,jy,jz) *delt   
    END IF
  END DO

!   Handle the drains

  DO drain = 1, ndrains
    jx = jxDrainLoc(drain)
    jy = jyDrainLoc(drain)
    jz = jzDrainLoc(drain)
    cellvolume = dzz(jx,jy,jz)*dxx(jx)*dyy(jy)
    ctvd(jx,jy,jz) = ctvd(jx,jy,jz) + cnv*qdrain(drain) * sn(icomp,jx,jy,jz) *delt/  &
       (sorp(jx,jy,jz)*por(jx,jy,jz)*cellvolume)
    cflowout(icomp) = cflowout(icomp) - cnv*qdrain(drain) * sn(icomp,jx,jy,jz) *delt 
  END DO

END IF   

! ***************************  End of MODFLOW block **********************************

RETURN
END SUBROUTINE tvd
