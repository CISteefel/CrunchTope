SUBROUTINE modified_vanGenuchten_model(psi, theta_r, theta_s, alpha, n, psi_s, &
                 theta, kr, dtheta, dkr)
USE crunchtype
IMPLICIT NONE

REAL(DP), INTENT(IN) :: psi, theta_r, theta_s, alpha, n, psi_s
REAL(DP), INTENT(OUT) :: theta, kr, dtheta, dkr
REAL(DP) :: m, l, theta_m, denominator, S_e, S_e_star, F, F_1, ratio
REAL(DP) :: term_1, term_2, term_3, dkr_dS_e, dS_e_dtheta, outer, inner
  
IF (psi >= psi_s) THEN
  theta = theta_s
  dtheta = 0.0d0
  dkr = 0.0d0
  kr = 1.0d0
ELSE
  m = 1.0d0 - 1.0d0/n
  l = 0.5d0
  
  ! constants for modified VG model
  theta_m = theta_r + (theta_s - theta_r)*(1 + (ABS(alpha*psi_s))**n)**(m)
  ratio = (theta_s - theta_r)/(theta_m - theta_r)
  F_1 = (1 - ratio**(1/m))**m
  
  ! evaluate theta
  denominator = (1 + (ABS(alpha*psi))**n)**(m)
  theta = theta_r + (theta_m - theta_r)/denominator
  
  ! evaluate dtheta/dpsi
  outer = m*(1 + (ABS(alpha * psi))**n)**(-m - 1)
  inner = n*alpha*(ABS(alpha * psi))**(n-1)
  dtheta = (theta_m - theta_r)*outer*inner
  
  ! evaluate kr
  S_e = (theta - theta_r)/(theta_s - theta_r)
  S_e_star = ratio*S_e
  F = (1-S_e_star**(1/m))**m
  kr = S_e**l*((1 - F)/(1 - F_1))**2
  kr = MAX(kr, 1.0d-33)
  
  ! evaluate dkr/dpsi
  ! dkr/dS_e
  term_1 = l*S_e**(l-1)*((1 - F)/(1 - F_1))**2
  term_2 = 2 * S_e**l * (1 - F)/(1 - F_1)**2
  term_3 = (1-(S_e_star)**(1/m))**(m-1) * (S_e_star)**(1/m-1) * ratio
  dkr_dS_e = term_1 + term_2*term_3
  
  ! dS_e/dtheta
  dS_e_dtheta = 1/(theta_s - theta_r)
  
  dkr = dkr_dS_e * dS_e_dtheta * dtheta

END IF

END SUBROUTINE modified_vanGenuchten_model