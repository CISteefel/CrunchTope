SUBROUTINE vanGenuchten_model_kr(psi, theta_r, theta_s, alpha, n, psi_s, &
                 kr)
USE crunchtype
IMPLICIT NONE

REAL(DP), INTENT(IN) :: psi, theta_r, theta_s, alpha, n, psi_s
REAL(DP), INTENT(OUT) :: kr
REAL(DP) :: m, l, theta, theta_m, denominator, S_e, S_e_star, F, F_1, ratio
  
IF (psi >= psi_s) THEN
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
  
  ! evaluate kr
  S_e = (theta - theta_r)/(theta_s - theta_r)
  S_e_star = ratio*S_e
  F = (1-S_e_star**(1/m))**m
  kr = S_e**l*((1 - F)/(1 - F_1))**2
  kr = MAX(kr, 1.0d-33)


END IF

END SUBROUTINE vanGenuchten_model_kr