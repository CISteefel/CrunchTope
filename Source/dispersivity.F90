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
    
SUBROUTINE dispersivity(nx,ny,nz)
USE crunchtype
USE transport
USE medium, ONLY: por
USE concentration, ONLY: jinit

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                           :: nx
INTEGER(I4B), INTENT(IN)                                           :: ny
INTEGER(I4B), INTENT(IN)                                           :: nz

!  Internal variables and arrays

REAL(DP)                                                           :: qbar
REAL(DP)                                                           :: vx
REAL(DP)                                                           :: vy
REAL(DP)                                                           :: vz
REAL(DP)                                                           :: MeanVtransverse
REAL(DP)                                                           :: vx_T
REAL(DP)                                                           :: qbar_T

REAL(DP)                                                           :: porAverage1
REAL(DP)                                                           :: porAverage2
REAL(DP)                                                           :: porAverage3
REAL(DP)                                                           :: porAverage4
REAL(DP)                                                           :: vx1
REAL(DP)                                                           :: vx2
REAL(DP)                                                           :: vx3
REAL(DP)                                                           :: vx4


INTEGER(I4B)                                                       :: jx
INTEGER(I4B)                                                       :: jy
INTEGER(I4B)                                                       :: jz
REAL(DP)                                                           :: satliq_x
REAL(DP)                                                           :: satliq_xx
REAL(DP)                                                           :: satliq_yy
REAL(DP)                                                           :: satliq_zz

IF (ny==1 .AND. nz==1) THEN      !!  1D case

  DO jz = 1,nz
    DO jy = 1,ny
      DO jx = 1,nx

        satliq_x = 0.5d0*(satliqold(jx,jy,jz)+satliq(jx,jy,jz))
        IF (jx == nx) THEN
        satliq_xx = satliq_x
        ELSE
        satliq_xx = 0.5d0*(satliqold(jx+1,jy,jz)+satliq(jx+1,jy,jz))
        ENDIF
        IF (jy == ny) THEN
        satliq_yy = satliq_x
        ELSE
        satliq_yy = 0.5d0*(satliqold(jx,jy+1,jz)+satliq(jx,jy+1,jz))
        ENDIF
        IF (jz == nz) THEN
        satliq_zz = satliq_x
        ELSE
        satliq_zz = 0.5d0*(satliqold(jx,jy,jz+1)+satliq(jx,jy,jz+1))
        ENDIF

        vx = qx(jx,jy,jz)/( (0.5d0*(por(jx,jy,jz)+por(jx+1,jy,jz))) * (0.5d0*(satliq_x + satliq_xx)))
        vy = qy(jx,jy,jz)/( (0.5d0*(por(jx,jy,jz)+por(jx,jy+1,jz))) * (0.5d0*(satliq_x + satliq_yy)))
        vz = qz(jx,jy,jz)/( (0.5d0*(por(jx,jy,jz)+por(jx,jy,jz+1))) * (0.5d0*(satliq_x + satliq_zz)))

        qbar = DSQRT( vx*vx + vy*vy + vz*vz )

        IF (qbar /= 0.0) THEN

          dspx(jx,jy,jz) = alft*( 0.5d0*(por(jx,jy,jz)+por(jx+1,jy,jz)) )*qbar +  &
            (alfl-alft)*( 0.5d0*(por(jx,jy,jz)+por(jx+1,jy,jz)) )*vx*vx/qbar
          
!!!        dspx(jx,jy,jz) = alft*qbar +  &
!!!            (alfl-alft)*vx*vx/qbar
        
!!!        dspx(jx,jy,jz) = alfl*vx*vx/qbar + alft*vy*vy/qbar
        
        ELSE
          dspx(jx,jy,jz) = 0.0
        END IF
        
      END DO
    END DO
  END DO

ELSE IF (ny>1 .and. nz==1) THEN             !!  2D case

  DO jz = 1,nz
    DO jy = 1,ny
      DO jx = 1,nx

        satliq_x = 0.5d0*(satliqold(jx,jy,jz)+satliq(jx,jy,jz))
        IF (jx == nx) THEN
        satliq_xx = satliq_x
        ELSE
        satliq_xx = 0.5d0*(satliqold(jx+1,jy,jz)+satliq(jx+1,jy,jz))
        ENDIF
        IF (jy == ny) THEN
        satliq_yy = satliq_x
        ELSE
        satliq_yy = 0.5d0*(satliqold(jx,jy+1,jz)+satliq(jx,jy+1,jz))
        ENDIF
        IF (jz == nz) THEN
        satliq_zz = satliq_x
        ELSE
        satliq_zz = 0.5d0*(satliqold(jx,jy,jz+1)+satliq(jx,jy,jz+1))
        ENDIF

        vx = qx(jx,jy,jz)/( (0.5d0*(por(jx,jy,jz)+por(jx+1,jy,jz))) * (0.5d0*(satliq_x + satliq_xx)))
        vy = qy(jx,jy,jz)/( (0.5d0*(por(jx,jy,jz)+por(jx,jy+1,jz))) * (0.5d0*(satliq_x + satliq_yy)))
        vz = qz(jx,jy,jz)/( (0.5d0*(por(jx,jy,jz)+por(jx,jy,jz+1))) * (0.5d0*(satliq_x + satliq_zz)))

        qbar = DSQRT( vx*vx + vy*vy + vz*vz )

        IF (qbar /= 0.0) THEN
          dspx(jx,jy,jz) = alft*qbar +  &
            (alfl-alft)*vx*vx/qbar   
        ELSE
          dspx(jx,jy,jz) = 0.0
        END IF
      END DO
    END DO
  END DO

!!!DO jz = 1,nz
  jz = 1
  DO jy = 1,ny
    DO jx = 1,nx

      IF (jy==ny) THEN

!! Average the suckah
      porAverage1 = 0.5d0*(por(jx,jy,jz) + por(jx+1,jy,jz) )
      porAverage2 = 0.5d0*(por(jx,jy+1,jz) + por(jx+1,jy+1,jz) )
      porAverage3 = 0.5d0*(por(jx-1,jy+1,jz) + por(jx,jy+1,jz) )
      porAverage4 = 0.5d0*(por(jx-1,jy,jz) + por(jx,jy,jz) )

      vx1 = qx(jx,jy,jz)/porAverage1
      vx2 = 0.0 
      vx3 = 0.0
      vx4 = qx(jx-1,jy,jz)/porAverage4
        
     ELSE

!! Average the suckah
      porAverage1 = 0.5d0*(por(jx,jy,jz) + por(jx+1,jy,jz) )
      porAverage2 = 0.5d0*(por(jx,jy+1,jz) + por(jx+1,jy+1,jz) )
      porAverage3 = 0.5d0*(por(jx-1,jy+1,jz) + por(jx,jy+1,jz) )
      porAverage4 = 0.5d0*(por(jx-1,jy,jz) + por(jx,jy,jz) )

      vx1 = qx(jx,jy,jz)/porAverage1
      vx2 = qx(jx,jy+1,jz)/porAverage2
      vx3 = qx(jx-1,jy+1,jz)/porAverage3
      vx4 = qx(jx-1,jy,jz)/porAverage4

     END IF
 
      vx_T = 0.25*( vx1 + vx2 + vx3 + vx4 )
      
      MeanVtransverse = vx_T  

      qbar_T = DSQRT( MeanVtransverse*MeanVtransverse)

      vy = 0.0d0
      vz = 0.0d0
  
      qbar = DSQRT( vx*vx + vy*vy + vz*vz )
      
      IF (qbar /= 0.0) THEN

!!        dspy(jx,jy,jz) =  alft*qbar +  &
!!            (alfl-alft)*qy(jx,jy,jz)*qy(jx,jy,jz)/qbar
          
!!            dspy(jx,jy,jz) = alft*( 0.5d0*(por(jx,jy,jz)+por(jx,jy+1,jz)) )*qbar +  &
!!            (alfl-alft)*( 0.5d0*(por(jx,jy,jz)+por(jx,jy+1,jz)) )*vy*vy/qbar
          
!!             dspy(jx,jy,jz) = alft*qbar +  &
!!            (alfl-alft)*vy*vy/qbar

             dspy(jx,jy,jz) = alfT*qbar_T 

      ELSE
        dspy(jx,jy,jz) = 0.0
      END IF

    END DO
  END DO
!!! END DO  !! JZ

ELSE

DO jz = 1,nz
  DO jy = 1,ny
    DO jx = 1,nx

      satliq_x = 0.5d0*(satliqold(jx,jy,jz)+satliq(jx,jy,jz))
      IF (jx == nx) THEN
      satliq_xx = satliq_x
      ELSE
      satliq_xx = 0.5d0*(satliqold(jx+1,jy,jz)+satliq(jx+1,jy,jz))
      ENDIF
      IF (jy == ny) THEN
      satliq_yy = satliq_x
      ELSE
      satliq_yy = 0.5d0*(satliqold(jx,jy+1,jz)+satliq(jx,jy+1,jz))
      ENDIF
      IF (jz == nz) THEN
      satliq_zz = satliq_x
      ELSE
      satliq_zz = 0.5d0*(satliqold(jx,jy,jz+1)+satliq(jx,jy,jz+1))
      ENDIF
     
      vx = qx(jx,jy,jz)/( (0.5d0*(por(jx,jy,jz)+por(jx+1,jy,jz))) * (0.5d0*(satliq_x + satliq_xx)))
      vy = qy(jx,jy,jz)/( (0.5d0*(por(jx,jy,jz)+por(jx,jy+1,jz))) * (0.5d0*(satliq_x + satliq_yy)))
      vz = qz(jx,jy,jz)/( (0.5d0*(por(jx,jy,jz)+por(jx,jy,jz+1))) * (0.5d0*(satliq_x + satliq_zz)))

      qbar = DSQRT( vx*vx + vy*vy + vz*vz )

      IF (qbar /= 0.0) THEN
        dspx(jx,jy,jz) = alft*( 0.5d0*(por(jx,jy,jz)+por(jx+1,jy,jz)) )*qbar +  &
            (alfl-alft)*( 0.5d0*(por(jx,jy,jz)+por(jx+1,jy,jz)) )*vx*vx/qbar
          
!!!        dspx(jx,jy,jz) = alft*qbar +  &
!!!            (alfl-alft)*vx*vx/qbar
        
!!!        dspx(jx,jy,jz) = alfl*vx*vx/qbar + alft*vy*vy/qbar
        
      ELSE
        dspx(jx,jy,jz) = 0.0
      END IF
    END DO
  END DO
END DO

DO jz = 1,nz
  DO jy = 1,ny
    DO jx = 1,nx
        
      vx = qx(jx,jy,jz) /( 0.5d0*(por(jx,jy,jz)+por(jx+1,jy,jz)) )

!! Average the suckah
      porAverage1 = 0.5d0*(por(jx,jy,jz) + por(jx+1,jy,jz) )
      porAverage2 = 0.5d0*(por(jx,jy+1,jz) + por(jx+1,jy+1,jz) )
      porAverage3 = 0.5d0*(por(jx-1,jy+1,jz) + por(jx,jy+1,jz) )
      porAverage4 = 0.5d0*(por(jx-1,jy,jz) + por(jx,jy,jz) )

      vx1 = qx(jx,jy,jz)
      vx2 = qx(jx,jy+1,jz)
      vx3 = qx(jx-1,jy+1,jz)
      vx4 = qx(jx-1,jy,jz)
 
      vx_T = 0.25*( vx1/porAverage1 + vx2/porAverage2 + vx3/porAverage3 + vx4/porAverage4 )
      
      MeanVtransverse = vx_T  

      qbar_T = DSQRT( MeanVtransverse*MeanVtransverse)

      vy = 0.0d0
      vz = 0.0d0
  
      qbar = DSQRT( vx*vx + vy*vy + vz*vz )
      
      IF (qbar /= 0.0) THEN

!!        dspy(jx,jy,jz) =  alft*qbar +  &
!!            (alfl-alft)*qy(jx,jy,jz)*qy(jx,jy,jz)/qbar
          
!!            dspy(jx,jy,jz) = alft*( 0.5d0*(por(jx,jy,jz)+por(jx,jy+1,jz)) )*qbar +  &
!!            (alfl-alft)*( 0.5d0*(por(jx,jy,jz)+por(jx,jy+1,jz)) )*vy*vy/qbar
          
!!             dspy(jx,jy,jz) = alft*qbar +  &
!!            (alfl-alft)*vy*vy/qbar

             dspy(jx,jy,jz) = alfT*qbar_T 

      ELSE
        dspy(jx,jy,jz) = 0.0
      END IF

    END DO
  END DO
END DO

DO jz = 1,nz
  DO jy = 1,ny
    DO jx = 1,nx
      qbar = SQRT( qx(jx,jy,jz)*qx(jx,jy,jz) + qy(jx,jy,jz)*qy(jx,jy,jz)  &
          + qz(jx,jy,jz)*qz(jx,jy,jz) )
      IF (qbar /= 0.0) THEN
        dspz(jx,jy,jz) = alft*qbar +  &
            (alfl-alft)*qz(jx,jy,jz)*qz(jx,jy,jz)/qbar
      ELSE
        dspz(jx,jy,jz) = 0.0
      END IF
    END DO
  END DO
END DO

END IF

RETURN
END SUBROUTINE dispersivity
