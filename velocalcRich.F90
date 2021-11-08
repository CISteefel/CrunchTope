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

SUBROUTINE velocalcRich(nx,ny,nz,dtyr)
USE crunchtype
USE params
USE medium
USE transport
USE temperature, ONLY: ro
USE flow
USE CrunchFunctions

IMPLICIT NONE

INTEGER(I4B), INTENT(IN)                                       :: nx
INTEGER(I4B), INTENT(IN)                                       :: ny
INTEGER(I4B), INTENT(IN)                                       :: nz

!  Internal variables and arrays

INTEGER(I4B)                                                          :: jx
INTEGER(I4B)                                                          :: jy
INTEGER(I4B)                                                          :: jz
INTEGER(I4B)                                                          :: npz

REAL(DP)                                                       :: dt,pumpterm
REAL(DP), INTENT(IN)                                           :: dtyr


REAL(DP)                                                  :: qgtinterp

!  ****** PARAMETERS  ****************************


CHARACTER (LEN=1)                                                     :: Coordinate

dt = dtyr * 365 * 86400

!   calculate darcy fluxes

DO jz = 1,nz
  DO jy = 1,ny
    DO jx = 1,nx

!!  *****    Qz and Qx    **************************

      IF (jx /= nx) THEN
        Coordinate = 'X'
        qx(jx,jy,jz) = -2.0d0 * Kfacx(jx,jy,jz) * (head(jx+1,jy,jz) - head(jx,jy,jz)) / (dxx(jx)+dxx(jx+1))
      END IF

      IF (jz /= nz) THEN
        Coordinate = 'Z'
        qz(jx,jy,jz) = -2.0d0 * Kfacz(jx,jy,jz) * (head(jx,jy,jz+1) - head(jx,jy,jz)) / (dzz(jx,jy,jz)+dzz(jx,jy,jz+1)) &
                        + Kfacz(jx,jy,jz)
        IF (activecellPressure(jx,jy,jz) == 0 .AND. activecellPressure(jx,jy,jz+1) == 1) THEN
            IF (head(jx,jy,jz) == 0.0d0) THEN
                ! no infiltration on dry surface, but exfiltration is ok
                IF (qz(jx,jy,jz) > 0.0d0) THEN
                    qz(jx,jy,jz) = 0.0d0
                END IF
            END IF
        END IF
      END IF

!!  *******             **************************

!!  *******      Qy     **************************

      IF (jy /= ny) THEN
        Coordinate = 'Y'
        if (jy == ny-1 .and. jx == nx-1) then
          continue
        end if
        qy(jx,jy,jz) = -2.0d0 * Kfacy(jx,jy,jz) * (head(jx,jy+1,jz) - head(jx,jy,jz)) / (dyy(jy)+dyy(jy+1))
        IF (y_is_vertical) THEN
            qz(jx,jy,jz) = 0.0d0
            qy(jx,jy,jz) = qy(jx,jy,jz) + Kfacy(jx,jy,jz)


!! Land surface, "fixed" cell at JY, active cell at JY+1  -- velocity = 0, or set to "pump" term

            IF (activecellPressure(jx,jy,jz) == 0 .AND. activecellPressure(jx,jy+1,jz) == 1) THEN
                qy(jx,jy,jz) = 2.0 * (qy(jx,jy,jz) - Kfacy(jx,jy,jz)) + Kfacy(jx,jy,jz)

                IF (head(jx,jy,jz) == 0.0d0) THEN
                    ! no infiltration on dry surface, but exfiltration is ok
                    IF (qy(jx,jy,jz) > 0.0d0) THEN
                        qy(jx,jy,jz) = 0.0d0
                    END IF
                END IF

                ! flux boundary condition
                pumpterm = 0.0d0
                IF (wells) THEN
                  DO npz = 1,npump(jx,jy+1,jz)
                    pumpterm = pumpterm + qg(npz,jx,jy+1,jz)/(secyr*dxx(jx)*dzz(jx,jy+1,jz))
                  END DO

                  IF (qy(jx,jy,jz) >= 0.0 .AND. pumpterm > 0.0) THEN
                      IF (npump(jx,jy+1,jz) > 0) THEN
                          qy(jx,jy,jz) = pumpterm
                      END IF
                  END IF


                ELSEIF (pumptimeseries) THEN
                
                  pumpterm = pumpterm + qg(1,jx,jy+1,jz)/(secyr*dxx(jx)*dzz(jx,jy+1,jz))
                  qy(jx,jy,jz) = pumpterm
                END IF

                ! check if there is enough room available
                IF (qy(jx,jy,jz) > 0.0 .AND. room(jx,jy+1,jz) < qy(jx,jy,jz) * dt * dxx(jx)* dzz(jx,jy+1,jz)) THEN
                    qy(jx,jy,jz) = room(jx,jy+1,jz) / dt * dxx(jx)* dzz(jx,jy+1,jz)
                END IF

            END IF
        ELSE
            WRITE(*,*) ' WARNING : Richards solver only works for 2D x-y now!'
        END IF
      END IF

    END DO
  END DO
END DO

!!  ********   End of Darcy flux   ************************

DO jz = 1,nz
  DO jx = 1,nx
      ! jy = 1

!!  ***  Calculate qy(j,0,jz)  ******************

    qy(jx,0,jz) = -2.0d0 * Kfacy(jx,0,jz) * (head(jx,1,jz) - head(jx,0,jz)) / (dyy(1))
    IF (y_is_vertical) THEN
        qy(jx,0,jz) = qy(jx,0,jz) + Kfacy(jx,0,jz)
        IF (head(jx,0,jz) == 0.0d0) THEN
            ! no infiltration on dry surface, but exfiltration is ok
            IF (qy(jx,0,jz) > 0.0d0) THEN
                qy(jx,0,jz) = 0.0d0
            END IF
        END IF
        ! Switch to flux BC if pressure BC is not defined
        IF (activecellPressure(jx,0,jz) == 1) THEN
            pumpterm = 0.0d0
            IF (wells) THEN
                DO npz = 1,npump(jx,1,jz)
                    pumpterm = pumpterm + qg(npz,jx,1,jz)/(secyr*dxx(jx)*dzz(jx,0,jz))
                END DO
                IF (npump(jx,1,jz) > 0) THEN
                    qy(jx,0,jz) = pumpterm
                ELSE
                    qy(jx,0,jz) = 0.0d0
                END IF
            
            ELSEIF (pumptimeseries) THEN
              pumpterm = pumpterm + qg(1,jx,1,jz)/(secyr*dxx(jx)*dzz(jx,0,jz))
              qy(jx,0,jz) = pumpterm
            ELSE
              ! If neither pump nor pressure is specified, qy = 0
              qy(jx,0,jz) = 0.0d0
            END IF
        END IF
    ELSE
        WRITE(*,*) ' WARNING : Richards solver only works for 2D x-y now!'
    END IF
    ! check if there is enough room available
    IF (qy(jx,0,jz) > 0.0 .AND. room(jx,1,jz) < qy(jx,0,jz) * dt * dxx(jx)* dzz(jx,1,jz)) THEN
        qy(jx,0,jz) = room(jx,1,jz) / dt * dxx(jx)* dzz(jx,1,jz)
    END IF

!!  ****   End of qy(jx,0,jz)   ******************

!!  ***  Calculate qy(j,NY,jz)  ******************

    IF (head(jx,ny+1,jz)>=head(jx,ny,jz) .AND. back_flow_closed) then
      qy(jx,ny,jz)=0
      Kfacy(jx,ny,jz)=0
    END IF

    ! jy = ny
    qy(jx,ny,jz) = -2.0d0 * Kfacy(jx,ny,jz) * (head(jx,ny+1,jz) - head(jx,ny,jz)) / (dyy(ny))
    ! gravity term
    IF (y_is_vertical) THEN
      qy(jx,ny,jz) = qy(jx,ny,jz) + Kfacy(jx,ny,jz)
      if (qy(jx,ny,jz) > 0.0) THEN
        continue
      END IF

      !forbid back flow


      ! pump source term
      IF (activecellPressure(jx,ny+1,jz) == 1) THEN

        IF (wells .or. pumptimeseries) THEN

          IF (qg(1,jx,ny,jz) /= 0.0) THEN
              ! free drainage
              qy(jx,ny,jz) = Kfacy(jx,ny,jz)
              ! qy(jx,ny,jz) = qg(1,jx,jy,jz)/(secyr*dxx(jx)*dzz(jx,jy,jz))
              IF (qy(jx,ny,jz) * dt > wc(jx,ny,jz) * dyy(ny)) THEN
                  qy(jx,ny,jz) = wc(jx,ny,jz) * dyy(ny) / dt
              END IF
          END IF

        END IF

      END IF

    ELSE
      WRITE(*,*) ' WARNING : Richards solver only works for 2D x-y now!'
    END IF
    ! check if there is enough room available
    IF (qy(jx,ny,jz) < 0.0 .AND. room(jx,ny,jz) < -qy(jx,ny,jz) * dt * dxx(jx)* dzz(jx,ny,jz)) THEN
      qy(jx,ny,jz) = -room(jx,ny,jz) / dt * dxx(jx)* dzz(jx,ny,jz)
    END IF

!!  ****   End of qy(jx,NY,jz)   ******************

  END DO
END DO

!!  ***  Calculate qx(0,jy,jz) and qx(NX,jy,jz)  ******************

DO jz = 1,nz
  DO jy = 1,ny
    IF (activecellPressure(0,jy,jz) == 0) THEN
      qx(0,jy,jz) = -2.0d0 * Kfacx(0,jy,jz) * (head(1,jy,jz) - head(0,jy,jz)) / (dxx(1))
    ELSE
      qx(0,jy,jz) = 0.0d0
    END IF

    IF (activecellPressure(nx+1,jy,jz) == 0) THEN
      qx(nx,jy,jz) = -2.0d0 * Kfacx(nx,jy,jz) * (head(nx+1,jy,jz) - head(nx,jy,jz)) / (dxx(nx))
    ELSE
      qx(nx,jy,jz) = 0.0d0
    END IF
  END DO
END DO

!!  ***  Calculate qz(jx,jy,0) and qz(jx,jy,NZ)  ******************

DO jy = 1,ny
  DO jx = 1,nx
    IF (activecellPressure(jx,jy,0) == 0) THEN
      qz(jx,jy,0) = -2.0d0 * Kfacz(jx,jy,0) * (head(jx,jy,1) - head(jx,jy,0)) / (dzz(jx,jy,1)) + Kfacz(jx,jy,0)
    ELSE
      qz(jx,jy,0) = 0.0d0
    END IF
    ! IF (activecellPressure(jx,jy,nz+1) == 0) THEN
    !   qz(jx,jy,nz) = -2.0d0 * Kfacz(jx,jy,nz) * (head(jx,jy,nz+1) - head(jx,jy,nz)) / (dzz(jx,jy,nz)) + Kfacz(jx,jy,nz)
    ! ELSE
    !   qz(jx,jy,nz) = 0.0d0
    ! END IF
    qz(jx,jy,nz) = 0.0d0
  END DO
END DO

!! Convert to m/yr
qx = qx*secyr
qy = qy*secyr
qz = qz*secyr

!!if (time>0.9) then
 !! WRITE(*,*) 'I am stopping'
 !! stop
!! end if
RETURN


END SUBROUTINE velocalcRich
