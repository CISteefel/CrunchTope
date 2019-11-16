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

SUBROUTINE tvdNuft(nx,ny,nz,delt,icomp)

USE crunchtype
USE concentration
USE medium
USE transport
USE temperature
USE flow
USE runtime

IMPLICIT NONE

INTEGER(I4B), INTENT(IN)                                               :: nx
INTEGER(I4B), INTENT(IN)                                               :: ny
INTEGER(I4B), INTENT(IN)                                               :: nz
INTEGER(I4B), INTENT(IN)                                               :: icomp
REAL(DP), INTENT(IN)                                                   :: delt
LOGICAL(LGT)                                                           :: TvdOption

REAL(DP) :: height,cfl,phi,tmp,sfx
REAL(DP) :: dtdx,dtdy,dtdz,area,cellvolume,CellArea,porsatro,PorRoSorp

INTEGER(I4B) :: i,j,k,jx,jy,jz,nco

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

TvdOption = .TRUE.

ftvd = 0.0
gtvd = 0.0
htvd = 0.0

!   x-direction fluxes

!fp! emit({#exec_list# &
!fp!            call update_borders(ctvd_info,ctvd)
!fp!            call update_borders(por_info,por)
!fp!            call update_borders(dzz_info,dzz)
!fp!           #});
!
!fp! data_sched({#expr# ctvd(jx,jy,jz) #});
DO jz=0,nz
  DO jy=0,ny
    DO jx=0,nx

      IF(jy /= 0 .AND. jz /= 0 .AND. xflow)THEN
        IF(qx(jx,jy,jz) >= 0)THEN !!! +ve velocity case
          area = dzz(jx,jy,jz) * dyy(jy)
          porsatro = por(jx,jy,jz)*satliq(jx,jy,jz)*ro(jx,jy,jz)
          IF (modflow) THEN
            cfl=(ABS(qx(jx,jy,jz))*delt/dxx(jx) )/(area*porsatro) 
          ELSE
            cfl=(ABS(qx(jx,jy,jz))*delt/dxx(jx) )/(porsatro) 
          ENDIF
!fp! compare_elem("ctvdshift",{#expr# ctvd(jx+1,jy,jz)#},{#expr#(nx+1)*ny*nz#},0,{#expr#jx+jy*(nx+1)+jz*(nx+1)*(ny)#},real8,20);
          ftvd(jx,jy,jz)=ctvd(jx+1,jy,jz)-ctvd(jx,jy,jz)
          IF(ABS(ftvd(jx,jy,jz)) > 1.e-30 .AND. TvdOption)THEN
            sfx = (ctvd(jx,jy,jz)-ctvd(jx-1,jy,jz))/ftvd(jx,jy,jz)
            tmp = MIN(2.d0,2.d0*sfx,0.33333333333333333D0*  &
                (2.d0-cfl+(1.d0+cfl)*sfx))
            phi = MAX(0.0D0,tmp)
          ELSE
            phi=0
          END IF
          ftvd(jx,jy,jz)=qx(jx,jy,jz)*  &
              ( ctvd(jx,jy,jz) + 0.5*(1.0-cfl)*phi*ftvd(jx,jy,jz) )
          
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
          porsatro = por(jx+1,jy,jz)*satliq(jx+1,jy,jz)*ro(jx+1,jy,jz)
          IF (modflow) THEN
            cfl=(ABS(qx(jx,jy,jz))*delt/dxx(jx+1) )/(area*porsatro) 
          ELSE
            cfl=(ABS(qx(jx,jy,jz))*delt/dxx(jx+1) )/(porsatro) 
          ENDIF
          ftvd(jx,jy,jz)=ctvd(jx,jy,jz)-ctvd(jx+1,jy,jz)
          IF(ABS(ftvd(jx,jy,jz)) > 1.e-30 .AND. TvdOption)THEN
            sfx=(ctvd(jx+1,jy,jz)-ctvd(jx+2,jy,jz))/ftvd(jx,jy,jz)
            phi=MAX(0.d0,MIN(2.d0,2.d0*sfx, 0.33333333333333333D0*  &
                (2.d0-cfl+(1.d0+cfl)*sfx)))
          ELSE
            phi=0
          END IF
          ftvd(jx,jy,jz)=qx(jx,jy,jz)*( ctvd(jx+1,jy,jz) +  &
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
          porsatro = por(jx,jy,jz)*satliq(jx,jy,jz)*ro(jx,jy,jz)
          IF (modflow) THEN
            cfl=(ABS(qy(jx,jy,jz))*delt/dyy(jy) )/(area*porsatro) 
          ELSE
            cfl=(ABS(qy(jx,jy,jz))*delt/dyy(jy) )/(porsatro) 
          ENDIF
          gtvd(jx,jy,jz)=ctvd(jx,jy+1,jz)-ctvd(jx,jy,jz)
!fp! compare_elem("gtvd 1",{#expr# gtvd(jx(1:nx),jy(0:ny),jz(1:nz))#},0,0,0,real8,20);
          IF(ABS(gtvd(jx,jy,jz)) > 1.e-30 .AND. TvdOption)THEN
            sfx=(ctvd(jx,jy,jz)-ctvd(jx,jy-1,jz))/gtvd(jx,jy,jz)
            phi=MAX(0.d0,MIN(2.d0,2.d0*sfx, 0.33333333333333333D0*  &
                (2.d0-cfl+(1.d0+cfl)*sfx)))
          ELSE
            phi=0
          END IF
          gtvd(jx,jy,jz)=qy(jx,jy,jz)*( ctvd(jx,jy,jz) +  &
              0.5*(1.0-cfl)*phi*gtvd(jx,jy,jz) )
!fp! compare_elem("gtvd 2",{#expr# gtvd(jx(1:nx),jy(0:ny),jz(1:nz))#},0,0,0,real8,20);
! write(*,*)'g+ ',qy(jx,jy,jz),CFL,phi,gtvd(jx,jy,jz)
          
        ELSE IF (yflow) THEN    !!! -ve velocity case
          area = dxx(jx) * dzz(jx,jy+1,jz)
          porsatro = por(jx,jy-1,jz)*satliq(jx,jy-1,jz)*ro(jx,jy-1,jz)
          IF (modflow) THEN
            cfl=(ABS(qy(jx,jy,jz))*delt/dyy(jy+1) )/(area*porsatro) 
          ELSE
            cfl=(ABS(qy(jx,jy,jz))*delt/dyy(jy+1) )/(porsatro) 
          ENDIF
          gtvd(jx,jy,jz)=ctvd(jx,jy,jz)-ctvd(jx,jy+1,jz)
!fp! compare_elem("gtvd 3",{#expr# gtvd(jx(1:nx),jy(0:ny),jz(1:nz))#},0,0,0,real8,20);
          IF(ABS(gtvd(jx,jy,jz)) > 1.e-30 .AND. TvdOption)THEN
            sfx=(ctvd(jx,jy+1,jz)-ctvd(jx,jy+2,jz))/gtvd(jx,jy,jz)
            phi=MAX(0.d0,MIN(2.d0,2.d0*sfx, 0.33333333333333333D0*  &
                (2.d0-cfl+(1.d0+cfl)*sfx)))
          ELSE
            phi=0
          END IF
          gtvd(jx,jy,jz)=qy(jx,jy,jz)*( ctvd(jx,jy+1,jz) +  &
              0.5*(1.0-cfl)*phi*gtvd(jx,jy,jz) )
!fp! compare_elem("gtvd 4",{#expr# gtvd(jx(1:nx),jy(0:ny),jz(1:nz))#},0,0,0,real8,20);
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
          porsatro = por(jx,jy,jz)*satliq(jx,jy,jz)*ro(jx,jy,jz)
          IF (modflow) THEN
            cfl=(ABS(qz(jx,jy,jz))*delt/dzz(jx,jy,jz) )/(area*porsatro) 
          ELSE
            cfl=(ABS(qz(jx,jy,jz))*delt/dzz(jx,jy,jz) )/(porsatro) 
          ENDIF
! write(*,*)'h+ ',jx,jy,jz,dzz(jx,jy,JZ),CFL
          htvd(jx,jy,jz)=ctvd(jx,jy,jz+1)-ctvd(jx,jy,jz)
          IF(ABS(htvd(jx,jy,jz)) > 1.e-30 .AND. TvdOption)THEN
            sfx=(ctvd(jx,jy,jz)-ctvd(jx,jy,jz-1))/htvd(jx,jy,jz)
            phi=MAX(0.d0,MIN(2.d0,2.d0*sfx, 0.33333333333333333D0*  &
                (2.d0-cfl+(1.d0+cfl)*sfx)))
          ELSE
            phi=0
          END IF
          htvd(jx,jy,jz)=qz(jx,jy,jz)*( ctvd(jx,jy,jz) +  &
              0.5*(1.0-cfl)*phi*htvd(jx,jy,jz) )
! write(*,*)'h+ ',jx,jy,jz,qz(jx,jy,jz),CFL,phi,htvd(jx,jy,jz)
          
        ELSE IF (zflow) THEN     !!! -ve velocity case
          height = 0.50D0*(dzz(jx,jy,jz)+dzz(jx,jy,jz+1))
          area = dxx(jx) * dyy(jy)
          porsatro = por(jx,jy,jz-1)*satliq(jx,jy,jz-1)*ro(jx,jy,jz-1)
          IF (modflow) THEN
            cfl=(ABS(qz(jx,jy,jz))*delt/dzz(jx,jy,jz+1) )/(area*porsatro) 
          ELSE
            cfl=(ABS(qz(jx,jy,jz))*delt/dzz(jx,jy,jz+1) )/(porsatro) 
          ENDIF
          htvd(jx,jy,jz)=ctvd(jx,jy,jz)-ctvd(jx,jy,jz+1)
          IF(ABS(htvd(jx,jy,jz)) > 1.e-30 .AND. TvdOption)THEN
            sfx=(ctvd(jx,jy,jz+1)-ctvd(jx,jy,jz+2))/htvd(jx,jy,jz)
            phi=MAX(0.d0,MIN(2.d0,2.d0*sfx, 0.33333333333333333D0*  &
                (2.d0-cfl+(1.d0+cfl)*sfx)))
          ELSE
            phi=0
          END IF
          htvd(jx,jy,jz)=qz(jx,jy,jz)*( ctvd(jx,jy,jz+1) +  &
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

!fp! emit({#exec_list# &
!fp!            call update_borders(ftvd_info,ftvd)
!fp!            call update_borders(gtvd_info,gtvd)
!fp!            call update_borders(htvd_info,htvd)
!fp!           #});
DO jz=1,nz
  DO jy=1,ny
    DO jx=1,nx
      dtdx = delt/dxx(jx)
      dtdy = delt/dyy(jy)
      dtdz = delt/dzz(jx,jy,jz)
      PorRoSorp = por(jx,jy,jz)*sorp(jx,jy,jz)
      IF (activecell(jx,jy,jz) /= 0) THEN
        ctvd(jx,jy,jz) = ( satliqOld(jx,jy,jz)*roOld(jx,jy,jz)*PorRoSorp*ctvd(jx,jy,jz)       & 
                        - ( dtdx*(ftvd(jx,jy,jz) - ftvd(jx-1,jy,jz)) +         &
                           dtdy*(gtvd(jx,jy,jz) - gtvd(jx,jy-1,jz)) +          &
                           dtdz*(htvd(jx,jy,jz) - htvd(jx,jy,jz-1)) )    )     &
                         /( satliq(jx,jy,jz)*ro(jx,jy,jz)*PorRoSorp)
      END IF
    END DO
  END DO
END DO

! * Handle the recharge (currently assumes all recharge is source through layer 1)

!fp! compare_flush();
!fp! compare("ctvd after cubes",{#expr# ctvd(1:nx,1:ny,1:nz)#}, real8, 0, 20);

!WRITE(*,*)' before TVD recharge: cflowin(icomp)', cflowin(icomp)

IF (irecharge == 2) THEN
  !!  Assumes for now a radial symmetry in cylindrical coordinates case (no 3D)
  jy = 1
  DO jz=1,nz
    DO jx=1,nx
!fp! set_index({#ident# jy #});
      IF (rectangular) THEN
        cellvolume = dzz(jx,jy,jz)*dxx(jx)*dyy(jy)
        CellArea = dzz(jz,jy,jz)*dxx(jx)
      ELSE 
        CellVolume = 3.1416*dyy(jy)*                                     &
          ( (x(jx)+0.5*dxx(jx))**2 - (x(jx)-0.5*dxx(jx))**2  )
        CellArea = 3.1416*( (x(jx)+0.5*dxx(jx))**2 - (x(jx)-0.5*dxx(jx))**2  )*dzz(jx,jy,jz)
      END IF
      PorRoSorp = por(jx,jy,jz)*sorp(jx,jy,jz)*ro(jx,jy,jz)
      ctvd(jx,jy,jz)= ( satliq(jx,jy,jz)*PorRoSorp*ctvd(jx,jy,jz)                  & 
                    + qrecharge(jx,jy)*scond(icomp,infiltration)*CellArea*delt/cellvolume  ) &
                   /(satliq(jx,jy,jz)*PorRoSorp)
      cflowin(icomp) = cflowin(icomp) + qrecharge(jx,jy) * scond(icomp,infiltration) *delt
    END DO
  END DO
ELSE IF (irecharge == 3) THEN
  jz = 1
  DO jy=1,ny
    DO jx=1,nx
!fp! set_index({#ident# jz #});
      IF (rectangular) THEN
        cellvolume = dzz(jx,jy,jz)*dxx(jx)*dyy(jy)
        CellArea = dyy(jy)*dxx(jx)
      ELSE 
        CellVolume = 3.1416*dyy(jy)*                                     &
          ( (x(jx)+0.5*dxx(jx))**2 - (x(jx)-0.5*dxx(jx))**2  )*dzz(jx,jy,jz)
        CellArea = 3.1416*( (x(jx)+0.5*dxx(jx))**2 - (x(jx)-0.5*dxx(jx))**2  )*dyy(jy)
      END IF
      PorRoSorp = por(jx,jy,jz)*sorp(jx,jy,jz)*ro(jx,jy,jz)
      ctvd(jx,jy,jz)= ( satliq(jx,jy,jz)*PorRoSorp*ctvd(jx,jy,jz)                  & 
                    + qrecharge(jx,jy)*scond(icomp,infiltration)*CellArea*delt/cellvolume  ) &
                   /(satliq(jx,jy,jz)*PorRoSorp)
      cflowin(icomp) = cflowin(icomp) + qrecharge(jx,jy) * scond(icomp,infiltration) *delt
    END DO
  END DO
END IF

DO jz = 1,nz
  DO jy=1,ny
    DO jx=1,nx
      IF (intbnd(jx,jy,jz) > 0) THEN
        nco = intbnd(jx,jy,jz)
        IF (rectangular) THEN
          cellvolume = dzz(jx,jy,jz)*dxx(jx)*dyy(jy)
        ELSE 
          CellVolume = 3.1416*dyy(jy)*                                     &
             ( (x(jx)+0.5*dxx(jx))**2 - (x(jx)-0.5*dxx(jx))**2  )
        END IF
!!      if (qg(jx,jy,jz) > 0.0) then
!!      write(*,*) ' Cellvolume inside TVD ',cellvolume
!!      write(*,*) ' qg = ',qg(jx,jy,jz)
!!      write(*,*) ' Nodes: ',jx,jy,jz
!!      write(*,*) ' satliq(jx,jy,jz)   : ', satliq(jx,jy,jz)
!!      write(*,*) ' satliqold(jx,jy,jz): ', satliqold(jx,jy,jz)
!!      write(*,*) ' PorRoSorp: ', PorRoSorp
!!      write(*,*) ulab(icomp)(1:15),scond(icomp,nco)
!!      write(*,*)
!!      write(*,*) ' ro(jx,jy,jz)   : ', ro(jx,jy,jz)
!!      write(*,*) ' roOld(jx,jy,jz): ', roOld(jx,jy,jz)
!!      end if
!!        PorRoSorp = por(jx,jy,jz)*sorp(jx,jy,jz)*ro(jx,jy,jz)
!!        ctvd(jx,jy,jz)= ( satliq(jx,jy,jz)*PorRoSorp*ctvd(jx,jy,jz)         &
!!                          + 0.5*(qg(jx,jy,jz)+qg(jx,jy,jz))*scond(icomp,nco)*delt/cellvolume  )  &
!!                        /(satliq(jx,jy,jz)*PorRoSorp)
!!        cflowin(icomp) = cflowin(icomp) + qg(jx,jy,jz) * scond(icomp,nco) *delt
       END IF
    END DO
  END DO
END DO

RETURN
END SUBROUTINE tvdNuft

