subroutine flux_Richards(nx, ny, nz)
! This subroutine computes water flux based on the Buckingham Darcy's law
  use crunchtype
  use io
  use params
  use runtime
  use concentration
  use medium
  use flow
  use transport

  implicit none

  integer(I4B), intent(IN)                                   :: nx
  integer(I4B), intent(IN)                                   :: ny
  integer(I4B), intent(IN)                                   :: nz

  integer(I4B)                                               :: i
  integer(I4B)                                               :: jx
  integer(I4B)                                               :: jy
  integer(I4B)                                               :: jz
  integer(I4B)                                               :: n_xfaces
  integer(I4B)                                               :: n_yfaces

  real(DP)                                                   :: psi_grad ! gradient of water potential
  real(DP)                                                   :: q_diff ! diffusion flow
  real(DP)                                                   :: q_grav ! gravitational flow

  real(DP)                                                   :: S_1
  real(DP)                                                   :: S_2
  real(DP)                                                   :: dx_1
  real(DP)                                                   :: dx_2
  real(DP)                                                   :: gravity_effect

  integer(I4B)                                               :: jx_1
  integer(I4B)                                               :: jx_2
  integer(I4B)                                               :: jy_1
  integer(I4B)                                               :: jy_2

  integer(I4B), dimension(2)                   :: coordinate_cell_1 ! coordinate of cell 1 for 2D problem
  integer(I4B), dimension(2)                   :: coordinate_cell_2 ! coordinate of cell 2 for 2D problem

  xi_2 = rho_water_2*g/mu_water

!*********************************************************************
  if (nx > 1 .and. ny == 1 .and. nz == 1) then ! one-dimensional problem
    ! apply van Genuchten model to all grid cells
    jy = 1
    jz = 1

    do jx = 0, nx + 1
      call vanGenuchten_model(psi(jx, jy, jz), theta_r(jx, jy, jz), theta_s(jx, jy, jz), VG_alpha(jx, jy, jz), VG_n(jx, jy, jz), &
                              theta(jx, jy, jz), kr(jx, jy, jz), dtheta(jx, jy, jz), dkr(jx, jy, jz))
      head(jx, jy, jz) = psi(jx, jy, jz) + SignGravity*COSD(x_angle)*x_2(jx)
    end do

    ! compute xi and kr at faces
    do jx = 0, nx
      xi_2_faces(jx) = (xi_2(jx, jy, jz)*dxx_2(jx + 1) + xi_2(jx + 1, jy, jz)*dxx_2(jx))/(dxx_2(jx + 1) + dxx_2(jx))
      kr_faces(jx) = (kr(jx, jy, jz)*dxx_2(jx + 1) + kr(jx + 1, jy, jz)*dxx_2(jx))/(dxx_2(jx + 1) + dxx_2(jx))
    end do

    ! flux at faces (only gravitational flow is upwinded)
    do jx = 0, nx
      psi_grad = (psi(jx + 1, jy, jz) - psi(jx, jy, jz))/(x_2(jx + 1) - x_2(jx))
      q_diff = -xi_2_faces(jx)*K_faces(jx)*kr_faces(jx)*psi_grad
      q_grav = -SignGravity*COSD(x_angle)*K_faces(jx)* &
               merge(kr(jx + 1, jy, jz)*xi_2(jx + 1, jy, jz), kr(jx, jy, jz)*xi_2(jx, jy, jz), &
                     head(jx + 1, jy, jz) - head(jx, jy, jz) >= 0.0d0)
      qx(jx, jy, jz) = q_diff + q_grav
    end do

!*********************************************************************
  else if (nx > 1 .and. ny > 1 .and. nz == 1) then ! two-dimensional problem

    ! apply van Genuchten model to all grid cells
    jz = 1
    do jy = 0, ny + 1
      do jx = 0, nx + 1
        call vanGenuchten_model(psi(jx, jy, jz), theta_r(jx, jy, jz), theta_s(jx, jy, jz), VG_alpha(jx, jy, jz), VG_n(jx, jy, jz), &
                                theta(jx, jy, jz), kr(jx, jy, jz), dtheta(jx, jy, jz), dkr(jx, jy, jz))
        head(jx, jy, jz) = psi(jx, jy, jz) + SignGravity*COSD(x_angle)*x_2(jx) + SignGravity*COSD(y_angle)*y_2(jy)
      end do
    end do

    if (domain_shape_flow == 'regular') then
      n_xfaces = nx*(ny + 1) ! number of faces
      n_yfaces = (nx + 1)*ny

    else
      write (*, *)
      write (*, *) ' Currently, only regular spatial domain is supported.'
      write (*, *)
      read (*, *)
      stop

    end if

    !compute xi and kr at faces by looping over faces
    do i = 1, n_xfaces + n_yfaces
      coordinate_cell_1 = cell_to_coordinate(face_to_cell(i, 1), :)
      coordinate_cell_2 = cell_to_coordinate(face_to_cell(i, 2), :)

      jx_1 = coordinate_cell_1(1)
      jx_2 = coordinate_cell_2(1)
      jy_1 = coordinate_cell_1(2)
      jy_2 = coordinate_cell_2(2)

      if (i < n_xfaces + 1) then
        ! faces paralell to x-axis
        dx_1 = dyy_2(jy_1)
        dx_2 = dyy_2(jy_2)
      else
        ! faces paralell to y-axis
        dx_1 = dxx_2(jx_1)
        dx_2 = dxx_2(jx_2)
      end if

      ! distance weighted harmonic mean
      S_1 = xi_2(jx_1, jy_1, 1)
      S_2 = xi_2(jx_2, jy_2, 1)

      xi_2_faces(i) = (S_1*dx_2 + S_2*dx_1)/(dx_1 + dx_2)

      S_1 = kr(jx_1, jy_1, 1)
      S_2 = kr(jx_2, jy_2, 1)

      kr_faces(i) = (S_1*dx_2 + S_2*dx_1)/(dx_1 + dx_2)

      ! compute fluxes
      if (i < n_xfaces + 1) then
        ! faces paralell to x-axis
        dx_1 = y_2(jy_2) - y_2(jy_1)
        gravity_effect = SignGravity*COSD(y_angle)
      else
        ! faces paralell to y-axis
        dx_1 = x_2(jx_2) - x_2(jx_1)
        gravity_effect = SignGravity*COSD(x_angle)
      end if

      psi_grad = (psi(jx_2, jy_2, jz) - psi(jx_1, jy_1, jz))/dx_1
      q_diff = -xi_2_faces(i)*K_faces(i)*kr_faces(i)*psi_grad
      q_grav = -gravity_effect*K_faces(i)* &
               merge(kr(jx_2, jy_2, jz)*xi_2(jx_2, jy_2, jz), kr(jx_1, jy_1, jz)*xi_2(jx_1, jy_1, jz), &
                     head(jx_2, jy_2, jz) - head(jx_1, jy_1, jz) >= 0.0d0)
      q_Richards(i) = q_diff + q_grav

    end do

!*********************************************************************
  else if (nx > 1 .and. ny > 1 .and. nz > 1) then
    write (*, *)
    write (*, *) ' Currently, three-dimensional Richards solver is supported.'
    write (*, *)
    read (*, *)
    stop
  end if

end subroutine flux_Richards
