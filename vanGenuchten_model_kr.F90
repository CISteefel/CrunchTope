SUBROUTINE vanGenuchten_model_kr(psi, theta_r, theta_s, alpha, n, &
                 kr)
USE crunchtype
IMPLICIT NONE

REAL(DP), INTENT(IN) :: psi, theta_r, theta_s, alpha, n
REAL(DP), INTENT(OUT) :: kr
REAL(DP) :: m, l, S_e, S_e_epsilon, kr_epsilon, slope
REAL(DP) :: term_1, term_2, term_3, term_4, term_5
  
IF (psi >= 0.0) THEN
  kr = 1.0d0
ELSE
  m = 1.0d0 - 1.0d0/n
  l = 0.5d0
  
  S_e_epsilon = 0.99d0 ! interpolation point in terms of saturation
  
  ! evaluate theta
  term_1 = ABS(alpha*psi)
  term_2 = (term_1)**n
  S_e = (1.0d0/(1.0d0 + term_2))**(m)
  
  ! evaluate kr (= S_e**l*(1-(1-S_e**(1/m))**m)**2) and dkr/dS_e
  IF (S_e < S_e_epsilon) THEN
    term_3 = 1.0d0/(1.0d0 + term_2) ! S_e**(1/m)
    term_4 = (1.0d0 - term_3)**m
    kr = S_e**l*(1.0d0 - term_4)**2
    kr = MAX(kr, 1.0d-33)
    
  ELSE
    kr_epsilon = S_e_epsilon**l*(1.0d0-(1.0d0-S_e_epsilon**(1.0d0/m))**m)**2
    slope = (1.0d0 - kr_epsilon)/(1.0d0 - S_e_epsilon)
    kr = kr_epsilon + slope*(S_e - S_e_epsilon)
    
  END IF
  
END IF

END SUBROUTINE vanGenuchten_model_kr