MODULE Richards_module
  
  USE crunchtype
  
  IMPLICIT NONE
  
  PRIVATE
  
  
  
  TYPE :: Richards_Options
    LOGICAL(LGT) :: is_steady ! True when solving the steady-state Richards equation to get the initial condition
    LOGICAL(LGT) :: is_print ! True if you want print statements from the Richards solver
    LOGICAL(LGT) :: vg_is_n ! True if the input to vg_n is the n parameter in the van Genuchten model, otherwise, the input value is interpreted as the m parameter
    LOGICAL(LGT) :: psi_is_head ! True if the primary variable psi in the Richards equation is pressure head [L] or not. If false, the input values for the initial and boundary conditions, and vg_alpha are interpreted as in terms of pressure [Pa].  
    LOGICAL(LGT) :: theta_s_is_porosity ! True if the input to theta_s is the same as the porosity
    LOGICAL(LGT) :: theta_r_is_S_r ! True if the input to theta_r is the residual saturation
  END TYPE Richards_Options
  
  PUBLIC Richards_Options!, &
  
  !CONTAINS
  
! ************************************************************************** !
  
  !SUBROUTINE
  !
  !END SUBROUTINE
  
END MODULE Richards_module