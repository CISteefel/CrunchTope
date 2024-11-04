SUBROUTINE RichardsArrayAllocation(nx,ny,nz)
USE crunchtype
USE flow
USE medium

IMPLICIT NONE

INTEGER(I4B), INTENT(IN)   :: nx, ny, nz
INTEGER(I4B)               :: jx, jy, jz

! get the discretization for the Ricahrds for one-dimension
IF (ALLOCATED(x_2)) THEN
  DEALLOCATE(x_2)
  ALLOCATE(x_2(0:nx+1))
ELSE
  ALLOCATE(x_2(0:nx+1))
END IF
  
IF (ALLOCATED(dxx_2)) THEN
  DEALLOCATE(dxx_2)
  ALLOCATE(dxx_2(0:nx+1))
ELSE
  ALLOCATE(dxx_2(0:nx+1))
END IF
  
DO jx = 1, nx
  x_2(jx) = x(jx)
  dxx_2(jx) = dxx(jx)
END DO
    
! fill the ghost cells
dxx_2(0) = dxx_2(1)
dxx_2(nx+1) = dxx_2(nx)
  
x_2(0) = -0.5d0*dxx_2(0)
x_2(nx+1) = x_2(nx)+0.5d0*(dxx_2(nx) + dxx_2(nx+1))

IF (nx > 1 .AND. ny > 1) THEN ! two-dimensional problem
  ! discretization in y direction  
  IF (ALLOCATED(y_2)) THEN
    DEALLOCATE(y_2)
    ALLOCATE(y_2(0:ny+1))
  ELSE
    ALLOCATE(y_2(0:ny+1))
  END IF
  
  IF (ALLOCATED(dyy_2)) THEN
    DEALLOCATE(dyy_2)
    ALLOCATE(dyy_2(0:ny+1))
  ELSE
    ALLOCATE(dyy_2(0:ny+1))
  END IF
  
  DO jy = 1, ny
    y_2(jy) = y(jy)
    dyy_2(jy) = dyy(jy)
  END DO
    
  ! fill the ghost cells
  dyy_2(0) = dyy_2(1)
  dyy_2(nx+1) = dyy_2(ny)
  
  y_2(0) = -0.5d0*dyy_2(0)
  y_2(nx+1) = y_2(ny)+0.5d0*(dyy_2(ny) + dyy_2(ny+1))

END IF


IF (nx > 1 .AND. ny > 1 .AND. nz > 1) THEN
  WRITE(*,*)
  WRITE(*,*) ' Currently, three-dimensional Richards solver is supported.'
  WRITE(*,*)
  READ(*,*)
  STOP
END IF
  
    
IF (ny == 1 .AND. nz == 1) THEN ! one-dimensional problem
  IF (ALLOCATED(mu_water)) THEN
    DEALLOCATE(mu_water)
    ALLOCATE(mu_water(0:nx+1,ny,nz))
  ELSE
    ALLOCATE(mu_water(0:nx+1,ny,nz))
  END IF

  IF (ALLOCATED(rho_water_2)) THEN
    DEALLOCATE(rho_water_2)
    ALLOCATE(rho_water_2(0:nx+1,ny,nz))
  ELSE
    ALLOCATE(rho_water_2(0:nx+1,ny,nz))
  END IF
    
  ! allocate state variable for the Richards equation
  ! water potential psi
  IF (ALLOCATED(psi)) THEN
    DEALLOCATE(psi)
    ALLOCATE(psi(0:nx+1, ny, nz))
  ELSE
    ALLOCATE(psi(0:nx+1, ny, nz))
  END IF
     
  ! volumetric water content theta
  IF (ALLOCATED(theta)) THEN
    DEALLOCATE(theta)
    ALLOCATE(theta(0:nx+1, ny, nz))
  ELSE
    ALLOCATE(theta(0:nx+1, ny, nz))
  END IF
    
  IF (ALLOCATED(theta_prev)) THEN
    DEALLOCATE(theta_prev)
    ALLOCATE(theta_prev(0:nx+1, ny, nz))
  ELSE
    ALLOCATE(theta_prev(0:nx+1, ny, nz))
  END IF

  ! derivative of volumetric water content theta
  IF (ALLOCATED(dtheta)) THEN
    DEALLOCATE(dtheta)
    ALLOCATE(dtheta(0:nx+1, ny, nz))
  ELSE
    ALLOCATE(dtheta(0:nx+1, ny, nz))
  END IF

  ! head values
  IF (ALLOCATED(head)) THEN
    DEALLOCATE(head)
    ALLOCATE(head(0:nx+1, ny, nz))
  ELSE
    ALLOCATE(head(0:nx+1, ny, nz))
  END IF

  ! relative permeability
  IF (ALLOCATED(kr)) THEN
    DEALLOCATE(kr)
    ALLOCATE(kr(0:nx+1, ny, nz))
  ELSE
    ALLOCATE(kr(0:nx+1, ny, nz))
  END IF

  ! derivative of relative permeability
  IF (ALLOCATED(dkr)) THEN
    DEALLOCATE(dkr)
    ALLOCATE(dkr(0:nx+1, ny, nz))
  ELSE
    ALLOCATE(dkr(0:nx+1, ny, nz))
  END IF
     
  ! relative permeability at cell faces
  IF (ALLOCATED(kr_faces)) THEN
    DEALLOCATE(kr_faces)
    ALLOCATE(kr_faces(0:nx, ny, nz))
  ELSE
    ALLOCATE(kr_faces(0:nx, ny, nz))
  END IF
     
  ! physical constant used in the Richards equation
  IF (ALLOCATED(xi_2)) THEN
    DEALLOCATE(xi_2)
    ALLOCATE(xi_2(0:nx+1, ny, nz))
  ELSE
    ALLOCATE(xi_2(0:nx+1, ny, nz))
  END IF
    
  ! physical constant used in the Richards equation at cell faces
  IF (ALLOCATED(xi_2_faces)) THEN
    DEALLOCATE(xi_2_faces)
    ALLOCATE(xi_2_faces(0:nx, ny, nz))
  ELSE
    ALLOCATE(xi_2_faces(0:nx, ny, nz))
  END IF
    
  ! allocate van Genuchten parameters
  IF (ALLOCATED(theta_r)) THEN
    DEALLOCATE(theta_r)
    ALLOCATE(theta_r(0:nx+1, ny, nz))
  ELSE
    ALLOCATE(theta_r(0:nx+1, ny, nz))
  END IF
            
  IF (ALLOCATED(theta_s)) THEN
      DEALLOCATE(theta_s)
      ALLOCATE(theta_s(0:nx+1, ny, nz))
  ELSE
    ALLOCATE(theta_s(0:nx+1, ny, nz))
  END IF
    
  IF (ALLOCATED(VG_alpha)) THEN
    DEALLOCATE(VG_alpha)
    ALLOCATE(VG_alpha(0:nx+1, ny, nz))
  ELSE
    ALLOCATE(VG_alpha(0:nx+1, ny, nz))
  END IF
    
  IF (ALLOCATED(VG_n)) THEN
    DEALLOCATE(VG_n)
    ALLOCATE(VG_n(0:nx+1, ny, nz))
  ELSE
    ALLOCATE(VG_n(0:nx+1, ny, nz))
  END IF
  
  ! allocate permeability at faces
  IF (ALLOCATED(K_faces)) THEN
    DEALLOCATE(K_faces)
    ALLOCATE(K_faces(0:nx))
  ELSE
    ALLOCATE(K_faces(0:nx))
  END IF
  
  
ELSE IF (nx > 1 .AND. ny > 1 .AND. nz == 1) THEN ! two-dimensional problem
  IF (ALLOCATED(mu_water)) THEN
    DEALLOCATE(mu_water)
    ALLOCATE(mu_water(0:nx+1,0:ny+1,nz))
  ELSE
    ALLOCATE(mu_water(0:nx+1,0:ny+1,nz))
  END IF

  IF (ALLOCATED(rho_water_2)) THEN
    DEALLOCATE(rho_water_2)
    ALLOCATE(rho_water_2(0:nx+1,0:ny+1,nz))
  ELSE
    ALLOCATE(rho_water_2(0:nx+1,0:ny+1,nz))
  END IF
    
  ! allocate state variable for the Richards equation
  ! water potential psi
  IF (ALLOCATED(psi)) THEN
    DEALLOCATE(psi)
    ALLOCATE(psi(0:nx+1, 0:ny+1, nz))
  ELSE
    ALLOCATE(psi(0:nx+1, 0:ny+1, nz))
  END IF
     
  ! volumetric water content theta
  IF (ALLOCATED(theta)) THEN
    DEALLOCATE(theta)
    ALLOCATE(theta(0:nx+1, 0:ny+1, nz))
  ELSE
    ALLOCATE(theta(0:nx+1, 0:ny+1, nz))
  END IF
    
  IF (ALLOCATED(theta_prev)) THEN
    DEALLOCATE(theta_prev)
    ALLOCATE(theta_prev(0:nx+1, 0:ny+1, nz))
  ELSE
    ALLOCATE(theta_prev(0:nx+1, 0:ny+1, nz))
  END IF

  ! derivative of volumetric water content theta
  IF (ALLOCATED(dtheta)) THEN
    DEALLOCATE(dtheta)
    ALLOCATE(dtheta(0:nx+1, 0:ny+1, nz))
  ELSE
    ALLOCATE(dtheta(0:nx+1, 0:ny+1, nz))
  END IF

  ! head values
  IF (ALLOCATED(head)) THEN
    DEALLOCATE(head)
    ALLOCATE(head(0:nx+1, 0:ny+1, nz))
  ELSE
    ALLOCATE(head(0:nx+1, 0:ny+1, nz))
  END IF

  ! relative permeability
  IF (ALLOCATED(kr)) THEN
    DEALLOCATE(kr)
    ALLOCATE(kr(0:nx+1, 0:ny+1, nz))
  ELSE
    ALLOCATE(kr(0:nx+1, 0:ny+1, nz))
  END IF

  ! derivative of relative permeability
  IF (ALLOCATED(dkr)) THEN
    DEALLOCATE(dkr)
    ALLOCATE(dkr(0:nx+1, 0:ny+1, nz))
  ELSE
    ALLOCATE(dkr(0:nx+1, 0:ny+1, nz))
  END IF
     
  ! relative permeability at cell faces
  IF (ALLOCATED(kr_faces)) THEN
    DEALLOCATE(kr_faces)
    ALLOCATE(kr_faces(0:nx, 0:ny+1, nz))
  ELSE
    ALLOCATE(kr_faces(0:nx, 0:ny+1, nz))
  END IF
     
  ! physical constant used in the Richards equation
  IF (ALLOCATED(xi_2)) THEN
    DEALLOCATE(xi_2)
    ALLOCATE(xi_2(0:nx+1, 0:ny+1, nz))
  ELSE
    ALLOCATE(xi_2(0:nx+1, 0:ny+1, nz))
  END IF
    
  ! physical constant used in the Richards equation at cell faces
  IF (ALLOCATED(xi_2_faces)) THEN
    DEALLOCATE(xi_2_faces)
    ALLOCATE(xi_2_faces(0:nx, 0:ny+1, nz))
  ELSE
    ALLOCATE(xi_2_faces(0:nx, 0:ny+1, nz))
  END IF
    
  ! allocate van Genuchten parameters
  IF (ALLOCATED(theta_r)) THEN
    DEALLOCATE(theta_r)
    ALLOCATE(theta_r(0:nx+1, 0:ny+1, nz))
  ELSE
    ALLOCATE(theta_r(0:nx+1, 0:ny+1, nz))
  END IF
            
  IF (ALLOCATED(theta_s)) THEN
      DEALLOCATE(theta_s)
      ALLOCATE(theta_s(0:nx+1, 0:ny+1, nz))
  ELSE
    ALLOCATE(theta_s(0:nx+1, 0:ny+1, nz))
  END IF
    
  IF (ALLOCATED(VG_alpha)) THEN
    DEALLOCATE(VG_alpha)
    ALLOCATE(VG_alpha(0:nx+1, 0:ny+1, nz))
  ELSE
    ALLOCATE(VG_alpha(0:nx+1, 0:ny+1, nz))
  END IF
    
  IF (ALLOCATED(VG_n)) THEN
    DEALLOCATE(VG_n)
    ALLOCATE(VG_n(0:nx+1, 0:ny+1, nz))
  ELSE
    ALLOCATE(VG_n(0:nx+1, 0:ny+1, nz))
  END IF
  
  ! allocate permeability at faces
  IF (ALLOCATED(K_faces_x)) THEN
    DEALLOCATE(K_faces_x)
    ALLOCATE(K_faces_x(0:nx, 0:ny+1, nz))
  ELSE
    ALLOCATE(K_faces_x(0:nx, 0:ny+1, nz))
  END IF
ELSE IF (nx > 1 .AND. ny > 1 .AND. nz > 1) THEN
  WRITE(*,*)
  WRITE(*,*) ' Currently, three-dimensional Richards solver is supported.'
  WRITE(*,*)
  READ(*,*)
  STOP
END IF
  
END SUBROUTINE RichardsArrayAllocation