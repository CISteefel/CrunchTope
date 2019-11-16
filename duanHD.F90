!-----------------------------------------------------------------------------------
!
! Lambda coefficients - Pitzer interaction parameters for CO2 with ions (Na)
!
!-----------------------------------------------------------------------------------
subroutine calc_lambda2(pg,tk,lambda)

  use crunchtype

  implicit none
  
! arguments
  real(dp), intent(in)    :: tk       ! temperature in Kelvin
  real(dp), intent(in)    :: pg       ! gas pressure in bar  
  real(dp), intent(out)   :: lambda   ! coefficients

! internal
  integer(i4b)            :: i
  real(dp)                :: tr_inv
  real(dp)                :: R630_inv
  real(dp), dimension(11) :: L

! Duan and Sun 2003 Chem. Geol. Table 2  
  real(dp), dimension(11), parameter :: &
   
  c = (/ -0.411370585d0, &  ! c(1)
          6.07632013d-4, &  ! c(2)
         97.5347708d0,   &  ! c(3)
          0.0d0,         &  ! c(4)
          0.0d0,         &  ! c(5)
          0.0d0,         &  ! c(6)
          0.0d0,         &  ! c(7)
         -0.0237622469d0,&  ! c(8)
          0.0170656236d0,&  ! c(9)
          0.0d0,         &  ! c(10)
          1.41335834d-5  &  ! c(11)
      /)
  
  tr_inv   = 1.0d0 / tk
  r630_inv = 1.0d0 / (630.d0 - tk)
  
  L(:)  = 0.0d0
    
  L(1)  = c(1)
  L(2)  = c(2) * tk
  L(3)  = c(3) * tr_inv
  L(4)  = c(4) * tk * tk
  L(5)  = c(5) * r630_inv
  L(6)  = c(6) * pg
  L(7)  = c(7) * pg * dlog(tk)
  L(8)  = c(8) * pg * tr_inv
  L(9)  = c(9) * pg * r630_inv
  L(10) = c(10)* pg * pg * r630_inv * r630_inv
  L(11) = c(11)* tk * dlog(pg)
  
  lambda = 0.0d0
  do i=1,11
   lambda = lambda + L(i)
  end do

  return

end subroutine calc_lambda2

!-----------------------------------------------------------------------------------
!
! Mu coefficients - standard chemical potential
!                   standard Gibbs energy change for the reaction
!
!-----------------------------------------------------------------------------------
subroutine calc_mu2(pg,tk,mu)

  use crunchtype

  implicit none
  
! arguments  
  real(dp), intent(in)    :: tk ! temperature in Kelvin
  real(dp), intent(in)    :: pg ! gas pressure in bar  
  real(dp), intent(out)   :: mu

! internal
  integer(i4b)            :: i
  real(dp)                :: tk_inv
  real(dp)                :: R630_inv 
  real(dp), dimension(11) :: L 

! Duan and Sun 2003 Chem. Geol. Table 2
  real(dp), dimension(11), parameter :: & 
  
  c = (/ 28.9447706d0,   &  ! c(1)
         -0.0354581768d0,&  ! c(2)
         -4770.67077d0,  &  ! c(3)
          1.02782768d-5, &  ! c(4)
         33.8126098d0,   &  ! c(5)
          9.04037140d-3, &  ! c(6)
         -1.14934031d-3, &  ! c(7)
         -0.307405726d0, &  ! c(8)
         -0.0907301486d0,&  ! c(9)
          9.32713393d-4, &  ! c(10)
          0.0d0          &  ! c(11)
      /)
  
  tk_inv   = 1.0d0 / tk
  r630_inv = 1.0d0 / (630.d0 - tk)
  
  L(:)  = 0.0d0
    
  L(1)  = c(1)
  L(2)  = c(2) * tk
  L(3)  = c(3) * tk_inv
  L(4)  = c(4) * tk * tk
  L(5)  = c(5) * r630_inv
  L(6)  = c(6) * pg
  L(7)  = c(7) * pg * dlog(tk)
  L(8)  = c(8) * pg * tk_inv
  L(9)  = c(9) * pg * r630_inv
  L(10) = c(10)* pg * pg * r630_inv * r630_inv
  L(11) = c(11)* tk * dlog(pg)
  
  mu = 0.0d0
  do i=1,11
   mu = mu + L(i)
  end do

  return

end subroutine calc_mu2
!-----------------------------------------------------------------------------------
!
! Xi coefficients - Pitzer interaction parameters for CO2 with Na and Cl
!
!-----------------------------------------------------------------------------------
subroutine calc_xi2(pg,tk,xi)

  use crunchtype

  implicit none

! arguments  
  real(dp), intent(in)    :: tk ! temperature in Kelvin
  real(dp), intent(in)    :: pg ! gas pressure in bar  
  real(dp), intent(out)   :: xi

! internal  
  integer(i4b) :: i
  real(dp) :: tk_inv
  real(dp) :: R630_inv
  real(dp), dimension(11) :: L 

! Duan and Sun 2003 Chem. Geol. Table 2
  real(dp), dimension(11), parameter :: & 
  
  c = (/  3.36389723d-4, &  ! c(1)
         -1.98298980d-5, &  ! c(2)
          0.0d0,         &  ! c(3)
          0.0d0,         &  ! c(4)
          0.0d0,         &  ! c(5)
          0.0d0,         &  ! c(6)
          0.0d0,         &  ! c(7)
          2.12220830d-3, &  ! c(8)
         -5.24873303d-3, &  ! c(9)
          0.0d0,         &  ! c(10)
          0.0d0          &  ! c(11)
      /)
  
  tk_inv   = 1.0d0 / tk
  r630_inv = 1.0d0 / (630.d0 - tk)
  
  L(:)  = 0.0d0
    
  L(1)  = c(1)
  L(2)  = c(2) * tk
  L(3)  = c(3) * tk_inv
  L(4)  = c(4) * tk * tk
  L(5)  = c(5) * r630_inv
  L(6)  = c(6) * pg
  L(7)  = c(7) * pg * dlog(tk)
  L(8)  = c(8) * pg * tk_inv
  L(9)  = c(9) * pg * r630_inv
  L(10) = c(10)* pg * pg * r630_inv * r630_inv
  L(11) = c(11)* tk * dlog(pg)
  
  xi = 0.0d0
  do i=1,11
   xi = xi + L(i)
  end do
  
  return
  
end subroutine calc_xi2
!-----------------------------------------------------------------------------------
!
! CO2 Fugacity - Duan et al 2006 Marine Chemistry
!
!-----------------------------------------------------------------------------------
subroutine fugacity_co22(pg,tk,ln_fco2,fco2)

  use crunchtype
  use params

  implicit none

! arguments
  real(dp), intent(in)    :: pg
  real(dp), intent(in)    :: tk
  real(dp), intent(out)   :: ln_fco2
  real(dp), intent(out)   :: fco2

! internal

  real(dp) :: p1
  real(dp), dimension(15)  :: a

  integer(i4b) :: Class

  real(dp) :: term1,   term2,   &
              term3,   term4,   &
              term5,   term6

IF (tk < 305) THEN
  P1 = 10 * ( -863.6 / tk +4.705)
    IF (pg < P1) THEN
      Class = 1;
    ELSE IF (pg < 1000) THEN
      Class = 2;
    ELSE
      Class = 3;
    END IF

ELSE IF (tk < 405) THEN
  P1 = 75 + (tk - 305) * 1.25;
    IF (pg < P1) THEN
      Class = 1;
    ELSE IF (pg < 1000) THEN
      IF (tk < 340) THEN
        Class = 2;
      ELSE 
        Class = 4;
      END IF
    ELSE
      IF (tk < 340) THEN
        Class = 3;
      ELSE 
        Class = 5;
      END IF
    END IF

ELSE
  P1 = 200
    IF (pg < P1) THEN
      Class = 1;
    ELSE IF (pg < 1000) THEN
      IF (tk < 435) THEN
        Class = 4;
      ELSE
        Class = 6;
      END IF
    ELSE
      IF (tk < 435) THEN
        Class = 5;
      ELSE
        Class = 6;
      END IF
    END IF
END IF
  

SELECT CASE (Class)
  CASE (1)
	
a = &

                                        (/ 1.0d0,  &  ! a(1)  
                                           4.7586835d-3,  &  ! a(2)
                                          -3.3569963d-6,  &  ! a(3)
                                           0.0d0,         &  ! a(4)
                                          -1.3179396d0,   &  ! a(5)
                                          -3.8389101d-6,  &  ! a(6)
                                           0.0d0,         &  ! a(7)
                                           2.2815104d-3,  &  ! a(8) 
                                           0.0d0,  &  ! a(9) 
                                           0.0d0,  &  ! a(10)
                                           0.0d0,  &  ! a(11)
                                           0.0d0,  &  ! a(12)
                                           0.0d0,  &  ! a(13)
                                           0.0d0,  &  ! a(14)
                                           0.0d0   /) ! a(15)
  CASE (2)
	
a = &

                                        (/ -7.1734882d-1,  &  ! a(1)  
                                           1.5985379d-4,  &  ! a(2)
                                          -4.9286471d-7,  &  ! a(3)
                                           0.0d0,         &  ! a(4)
                                           0.0d0,   &  ! a(5)
                                          -2.7855285d-7,  &  ! a(6)
                                           1.1877015d-9,         &  ! a(7)
                                           0.0d0,  &  ! a(8) 
                                           0.0d0,  &  ! a(9) 
                                           0.0d0,  &  ! a(10)
                                           0.0d0,  &  ! a(11)
                                           -9.6539512d1,  &  ! a(12)
                                           4.4774938d-1,  &  ! a(13)
                                           1.0181078d2,  &  ! a(14)
                                           5.3783879d-6   /) ! a(15)
  CASE (3)

 a = &

                                        (/ -6.5129019d-2,  &  ! a(1)  
                                          -2.1429977d-4,  &  ! a(2)
                                          -1.1444930d-6,  &  ! a(3)
                                           0.0d0,         &  ! a(4)
                                           0.0d0,   &  ! a(5)
                                          -1.1558081d-7,  &  ! a(6)
                                           1.1952370d-9,         &  ! a(7)
                                           0.0d0,  &  ! a(8) 
                                           0.0d0,  &  ! a(9) 
                                           0.0d0,  &  ! a(10)
                                           0.0d0,  &  ! a(11)
                                           -2.2134306d2,  &  ! a(12)
                                           0.0d0,  &  ! a(13)
                                           7.1820393d1,  &  ! a(14)
                                           6.6089246d-6   /) ! a(15)
  CASE (4)
	
 a = &

                                        (/ 5.0383896d0,  &  ! a(1)  
                                          -4.4257744d-3,  &  ! a(2)
                                           0.0d0,  &  ! a(3)
                                           1.9572733d0,         &  ! a(4)
                                           0.0d0,   &  ! a(5)
                                           2.4223436d-6,  &  ! a(6)
                                           0.0d0,         &  ! a(7)
                                          -9.3796135d-4,  &  ! a(8) 
                                          -1.5026030d0,  &  ! a(9) 
                                           3.0272240d-3,  &  ! a(10)
                                          -3.1377342d1,  &  ! a(11)
                                          -1.2847063d1,  &  ! a(12)
                                           0.0d0,  &  ! a(13)
                                           0.0d0,  &  ! a(14)
                                          -1.5056648d-5   /) ! a(15)
  CASE (5)
	
a = &

                                        (/ -1.6063152d1,  &  ! a(1)  
                                          -2.7057990d-3,  &  ! a(2)
                                           0.0d0,  &  ! a(3)
                                           1.4119239d-1,         &  ! a(4)
                                           0.0d0,   &  ! a(5)
                                           8.1132965d-7,  &  ! a(6)
                                           0.0d0,         &  ! a(7)
                                          -1.1453082d-4,  &  ! a(8) 
                                           2.3895671d0,  &  ! a(9) 
                                           5.0527457d-4,  &  ! a(10)
                                          -1.7763460d1,  &  ! a(11)
                                           9.8592232d2,  &  ! a(12)
                                           0.0d0,  &  ! a(13)
                                           0.0d0,  &  ! a(14)
                                          -5.4965256d-7   /) ! a(15)
  CASE (6)
	
a = &

                                        (/ -1.5693490d-1,  &  ! a(1)  
                                           4.4621407d-4,  &  ! a(2)
                                          -9.1080591d-7,  &  ! a(3)
                                           0.0d0,         &  ! a(4)
                                           0.0d0,   &  ! a(5)
                                           1.0647399d-7,  &  ! a(6)
                                           2.4273357d-10,         &  ! a(7)
                                           0.0d0,  &  ! a(8) 
                                           3.5874255d-1,  &  ! a(9) 
                                           6.3319710d-5,  &  ! a(10)
                                          -2.4989661d2,  &  ! a(11)
                                           0.0d0,  &  ! a(12)
                                           0.0d0,  &  ! a(13)
                                           8.8876800d2,  &  ! a(14)
                                          -6.6348003d-7   /) ! a(15)


END SELECT

  
  term1 = ( a(2) + a(3) * tk + a(4) / tk + a(5) /(tk - 150) ) * pg
  term2 = ( a(6) + a(7) * tk + a(8) / tk ) * pg **2.0d0
  term3 = ( a(9) + a(10) * tk + a(11) / tk ) * dlog(pg)
  term4 = ( a(12) + a(13) * tk ) / pg
  term5 = a(14) / tk
  term6 = a(15) * tk ** 2.0d0
  
  fco2 = a(1)      &
       + term1     &
       + term2     &
       + term3     &
       + term4     &
       + term5     &
       + term6  
  ln_fco2 = LOG(fco2)


  write(*,*)
  write(*,*)'P(bar)= ',pg
  write(*,*)'T(K)  = ',tk
  write(*,*)'ln(fugacity coeff CO2)= ',ln_fco2
  write(*,*)'fugacity coeff CO2= ',fco2
  write(*,*)
  write(*,*)'-------------------------------------------------------'

  
  return

end subroutine fugacity_co22

!-----------------------------------------------------------------------------------
!
! Lambda coefficients - Pitzer interaction parameters for CO2 with ions (Na)
!
!-----------------------------------------------------------------------------------
subroutine calc_ph2o2(tk,ph2o)

  use crunchtype

  implicit none
  
! arguments
  real(dp), intent(in)    :: tk       ! temperature in Kelvin
  !!real(dp), intent(in)    :: pg       ! gas pressure in bar  
  real(dp), intent(out)   :: ph2o     ! coefficients

! internal
  integer(i4b)            :: i
  real(dp)                :: tau

  real(dp) :: term2,   term3,   &
              term4,   term5,   &
              term6,   term0

! parameters: critical pressure and temperature for H2O
  real(dp), parameter     :: pc = 220.85d0, &  !bar
                             tc = 647.29d0     !Kelvin

! Duan and Sun 2003 Chem. Geol. Table B1  
  real(dp), dimension(5), parameter :: &
   
  c = (/ -38.640844d0,    &  ! c(1)
           5.8948420d0,   &  ! c(2)
          59.876516d0,    &  ! c(3)
          26.654627d0,    &  ! c(4)
          10.637097d0     &  ! c(5)
      /)
  
  tau    = (tk - tc) / tc

  term3  = tau
  term4  = term3 * tau
  term5  = term4 * tau
  term6  = term5 * tau

  term0  = pc * tk / tc    
  term2  = c(1) * (-tau)**1.9d0
  term3  = c(2) * term3
  term4  = c(3) * term4
  term5  = c(4) * term5
  term6  = c(5) * term6

  ph2o = term0 * (1.0d0 + &
                  term2 + &
                  term3 + &
                  term4 + &
                  term5 + &
                  term6 )  

  return

end subroutine calc_ph2o2

