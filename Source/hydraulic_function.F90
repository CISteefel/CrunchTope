MODULE hydraulic_function_module

USE crunchtype
  
IMPLICIT NONE
  
PRIVATE


TYPE :: VGMparameters ! store parameters for van Genuchten Mualem model
  REAL(DP), ALLOCATABLE :: theta_r(:,:,:) ! residual volumetric water content [L3 L-3]
  REAL(DP), ALLOCATABLE :: theta_s(:,:,:) ! saturated volumetric water content [L3 L-3]
  REAL(DP), ALLOCATABLE :: alpha(:,:,:) ! van Genuchten alpha parameter [L-1]
  REAL(DP), ALLOCATABLE :: n(:,:,:) ! van Genuchten n parameter [-]
END TYPE VGMparameters

TYPE(VGMparameters), PUBLIC :: VGM_parameters

!SUBROUTINE vanGenuchtenMualem(psi, theta_r, theta_s, alpha, n, &
!                 theta, kr, dtheta, dkr)
!USE crunchtype
!IMPLICIT NONE
!
!REAL(DP), INTENT(IN) :: psi, theta_r, theta_s, alpha, n
!REAL(DP), INTENT(OUT) :: theta, kr, dtheta, dkr
!REAL(DP) :: m, l, S_e, S_e_epsilon, kr_epsilon, slope
!REAL(DP) :: term_1, term_2, term_3, term_4, term_5, term_6, outer, inner
!REAL(DP) :: dkr_dS_e, dS_e_dtheta
!  
!IF (psi >= 0.0) THEN
!  theta = theta_s
!  dtheta = 0.0d0
!  dkr = 0.0d0
!  kr = 1.0d0
!ELSE
!  m = 1.0d0 - 1.0d0/n
!  l = 0.5d0
!  
!  S_e_epsilon = 0.999d0 ! interpolation point in terms of saturation
!  
!  ! evaluate theta
!  term_1 = ABS(alpha*psi)
!  term_2 = (term_1)**n
!  S_e = (1.0d0/(1.0d0 + term_2))**(m)
!  theta =  S_e * (theta_s - theta_r) + theta_r
!  
!  ! evaluate dtheta/dpsi
!  outer = m*(1.0d0 + term_2)**(-m - 1.0d0)
!  inner = n*alpha*(term_1)**(n-1.0d0)
!  dtheta = (theta_s - theta_r)*outer*inner
!  
!  ! evaluate kr (= S_e**l*(1-(1-S_e**(1/m))**m)**2) and dkr/dS_e
!  IF (S_e < S_e_epsilon) THEN
!    term_3 = 1.0d0/(1.0d0 + term_2) ! S_e**(1/m)
!    term_4 = (1.0d0 - term_3)**m
!    kr = S_e**l*(1.0d0 - term_4)**2
!    kr = MAX(kr, 1.0d-33)
!    
!    term_5 = term_4 * ((l + 2.0d0)*term_3 - l) - l * term_3 + l 
!    dkr_dS_e = S_e**(l - 1.0d0)*((term_4 - 1.0d0) * term_5)/(term_3 - 1.0d0)
!  ELSE
!    kr_epsilon = S_e_epsilon**l*(1.0d0-(1.0d0-S_e_epsilon**(1.0d0/m))**m)**2
!    slope = (1.0d0 - kr_epsilon)/(1.0d0 - S_e_epsilon)
!    kr = kr_epsilon + slope*(S_e - S_e_epsilon)
!    dkr_dS_e = slope
!    
!  END IF
!  
!  ! evaluate dkr/dpsi
!  ! dS_e/dtheta
!  dS_e_dtheta = 1.0d0/(theta_s - theta_r)
!  
!  dkr = dkr_dS_e * dS_e_dtheta * dtheta
!
!END IF
!
!END SUBROUTINE vanGenuchten_model


END MODULE hydraulic_function_module