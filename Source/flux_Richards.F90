SUBROUTINE flux_Richards(nx, ny, nz)
! This subroutine computes water flux based on the Buckingham Darcy's law
  USE crunchtype
  USE io
  USE params
  USE runtime
  USE concentration
  USE medium
  USE flow
  USE transport

  IMPLICIT NONE

  INTEGER(I4B), INTENT(IN)                                   :: nx
  INTEGER(I4B), INTENT(IN)                                   :: ny
  INTEGER(I4B), INTENT(IN)                                   :: nz

  INTEGER(I4B)                                               :: i
  INTEGER(I4B)                                               :: jx
  INTEGER(I4B)                                               :: jy
  INTEGER(I4B)                                               :: jz
  INTEGER(I4B)                                               :: n_xfaces
  INTEGER(I4B)                                               :: n_yfaces

  REAL(DP)                                                   :: psi_grad ! gradient of water potential
  REAL(DP)                                                   :: q_diff ! diffusion flow
  REAL(DP)                                                   :: q_grav ! gravitational flow

  REAL(DP)                                                   :: S_1
  REAL(DP)                                                   :: S_2
  REAL(DP)                                                   :: dx_1
  REAL(DP)                                                   :: dx_2
  REAL(DP)                                                   :: gravity_effect

  INTEGER(I4B)                                               :: jx_1
  INTEGER(I4B)                                               :: jx_2
  INTEGER(I4B)                                               :: jy_1
  INTEGER(I4B)                                               :: jy_2

  INTEGER(I4B), DIMENSION(2)                   :: coordinate_cell_1 ! coordinate of cell 1 for 2D problem
  INTEGER(I4B), DIMENSION(2)                   :: coordinate_cell_2 ! coordinate of cell 2 for 2D problem

  xi_2 = rho_water_2*g/mu_water

!*********************************************************************
  IF (nx > 1 .AND. ny == 1 .AND. nz == 1) THEN ! one-dimensional problem
    ! apply van Genuchten model to all grid cells
    jy = 1
    jz = 1

    DO jx = 0, nx + 1
      CALL vanGenuchten_model(psi(jx, jy, jz), theta_r(jx, jy, jz), theta_s(jx, jy, jz), VG_alpha(jx, jy, jz), VG_n(jx, jy, jz), &
                              theta(jx, jy, jz), kr(jx, jy, jz), dtheta(jx, jy, jz), dkr(jx, jy, jz))
      head(jx, jy, jz) = psi(jx, jy, jz) + SignGravity*COSD(x_angle)*x_2(jx)
    END DO

    ! compute xi and kr at faces
    DO jx = 0, nx
      xi_2_faces(jx) = (xi_2(jx, jy, jz)*dxx_2(jx + 1) + xi_2(jx + 1, jy, jz)*dxx_2(jx))/(dxx_2(jx + 1) + dxx_2(jx))
      kr_faces(jx) = (kr(jx, jy, jz)*dxx_2(jx + 1) + kr(jx + 1, jy, jz)*dxx_2(jx))/(dxx_2(jx + 1) + dxx_2(jx))
    END DO

    ! flux at faces (only gravitational flow is upwinded)
    DO jx = 0, nx
      psi_grad = (psi(jx + 1, jy, jz) - psi(jx, jy, jz))/(x_2(jx + 1) - x_2(jx))
      q_diff = -xi_2_faces(jx)*K_faces(jx)*kr_faces(jx)*psi_grad
      q_grav = -SignGravity*COSD(x_angle)*K_faces(jx)* &
               MERGE(kr(jx + 1, jy, jz)*xi_2(jx + 1, jy, jz), kr(jx, jy, jz)*xi_2(jx, jy, jz), &
                     head(jx + 1, jy, jz) - head(jx, jy, jz) >= 0.0D0)
      qx(jx, jy, jz) = q_diff + q_grav
    END DO

!*********************************************************************
  ELSE IF (nx > 1 .AND. ny > 1 .AND. nz == 1) THEN ! two-dimensional problem

    ! apply van Genuchten model to all grid cells
    jz = 1
    DO jy = 0, ny + 1
      DO jx = 0, nx + 1
        CALL vanGenuchten_model(psi(jx, jy, jz), theta_r(jx, jy, jz), theta_s(jx, jy, jz), VG_alpha(jx, jy, jz), VG_n(jx, jy, jz), &
                                theta(jx, jy, jz), kr(jx, jy, jz), dtheta(jx, jy, jz), dkr(jx, jy, jz))
        head(jx, jy, jz) = psi(jx, jy, jz) + SignGravity*COSD(x_angle)*x_2(jx) + SignGravity*COSD(y_angle)*y_2(jy)
      END DO
    END DO

    IF (domain_shape_flow == 'regular') THEN
      n_xfaces = nx*(ny + 1) ! number of faces
      n_yfaces = (nx + 1)*ny

    ELSE
      WRITE (*, *)
      WRITE (*, *) ' Currently, only regular spatial domain is supported.'
      WRITE (*, *)
      READ (*, *)
      STOP

    END IF

    !compute xi and kr at faces by looping over faces
    DO i = 1, n_xfaces + n_yfaces
      coordinate_cell_1 = cell_to_coordinate(face_to_cell(i, 1), :)
      coordinate_cell_2 = cell_to_coordinate(face_to_cell(i, 2), :)

      jx_1 = coordinate_cell_1(1)
      jx_2 = coordinate_cell_2(1)
      jy_1 = coordinate_cell_1(2)
      jy_2 = coordinate_cell_2(2)

      IF (i < n_xfaces + 1) THEN
        ! faces paralell to x-axis
        dx_1 = dyy_2(jy_1)
        dx_2 = dyy_2(jy_2)
      ELSE
        ! faces paralell to y-axis
        dx_1 = dxx_2(jx_1)
        dx_2 = dxx_2(jx_2)
      END IF

      ! distance weighted harmonic mean
      S_1 = xi_2(jx_1, jy_1, 1)
      S_2 = xi_2(jx_2, jy_2, 1)

      xi_2_faces(i) = (S_1*dx_2 + S_2*dx_1)/(dx_1 + dx_2)

      S_1 = kr(jx_1, jy_1, 1)
      S_2 = kr(jx_2, jy_2, 1)

      kr_faces(i) = (S_1*dx_2 + S_2*dx_1)/(dx_1 + dx_2)

      ! compute fluxes
      IF (i < n_xfaces + 1) THEN
        ! faces paralell to x-axis
        dx_1 = y_2(jy_2) - y_2(jy_1)
        gravity_effect = SignGravity*COSD(y_angle)
      ELSE
        ! faces paralell to y-axis
        dx_1 = x_2(jx_2) - x_2(jx_1)
        gravity_effect = SignGravity*COSD(x_angle)
      END IF

      psi_grad = (psi(jx_2, jy_2, jz) - psi(jx_1, jy_1, jz))/dx_1
      q_diff = -xi_2_faces(i)*K_faces(i)*kr_faces(i)*psi_grad
      q_grav = -gravity_effect*K_faces(i)* &
               MERGE(kr(jx_2, jy_2, jz)*xi_2(jx_2, jy_2, jz), kr(jx_1, jy_1, jz)*xi_2(jx_1, jy_1, jz), &
                     head(jx_2, jy_2, jz) - head(jx_1, jy_1, jz) >= 0.0D0)
      q_Richards(i) = q_diff + q_grav

    END DO

!*********************************************************************
  ELSE IF (nx > 1 .AND. ny > 1 .AND. nz > 1) THEN
    WRITE (*, *)
    WRITE (*, *) ' Currently, three-dimensional Richards solver is supported.'
    WRITE (*, *)
    READ (*, *)
    STOP
  END IF

END SUBROUTINE flux_Richards
