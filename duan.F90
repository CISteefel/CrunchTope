!!! *** Copyright Notice ***
!!! “CrunchFlow”, Copyright (c) 2016, The Regents of the University of California, through Lawrence Berkeley National Laboratory 
!!! (subject to receipt of any required approvals from the U.S. Dept. of Energy).  All rights reserved.
!!! 
!!! If you have questions about your rights to use or distribute this software, please contact 
!!! Berkeley Lab's Innovation & Partnerships Office at  IPO@lbl.gov.
!!! 
!!! NOTICE.  This Software was developed under funding from the U.S. Department of Energy and the U.S. Government 
!!! consequently retains certain rights. As such, the U.S. Government has been granted for itself and others acting 
!!! on its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the Software to reproduce, distribute copies to the public, 
!!! prepare derivative works, and perform publicly and display publicly, and to permit other to do so.
!!!
!!! *** License Agreement ***
!!! “CrunchFlow”, Copyright (c) 2016, The Regents of the University of California, through Lawrence Berkeley National Laboratory)
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
!-----------------------------------------------------------------------------------
!
! Lambda coefficients - Pitzer interaction parameters for CO2 with ions (Na)
! Written by Sergi Molins
!
!-----------------------------------------------------------------------------------

!
! Compare to Duan's model online: http://models.kl-edi.ac.cn/models/h2o_co2/index.htm
!
!-----------------------------------------------------------------------------------
!
! Lambda coefficients - Pitzer interaction parameters for CO2 with ions (Na)
!
!-----------------------------------------------------------------------------------
subroutine calc_lambda(pg,tk,lambda)

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

end subroutine calc_lambda

!-----------------------------------------------------------------------------------
!
! Mu coefficients - standard chemical potential
!                   standard Gibbs energy change for the reaction
!
!-----------------------------------------------------------------------------------
subroutine calc_mu(pg,tk,mu)

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

end subroutine calc_mu
!-----------------------------------------------------------------------------------
!
! Xi coefficients - Pitzer interaction parameters for CO2 with Na and Cl
!
!-----------------------------------------------------------------------------------
subroutine calc_xi(pg,tk,xi)

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
  
end subroutine calc_xi
!-----------------------------------------------------------------------------------
!
! CO2 Fugacity - Duan and Sun 2003 Chem. Geol. - from Duan et al 1992 GCA 56
!
!-----------------------------------------------------------------------------------
subroutine fugacity_co2(pg,tk,ln_fco2,vrInOut)

  use crunchtype
  use params

  implicit none

! arguments
  real(dp), intent(in)    :: pg
  real(dp), intent(in)    :: tk
  real(dp), intent(out)   :: ln_fco2
  real(dp), intent(inout)  :: vrInOut

! internal
  logical(lgt) :: RunVerbose
  integer(i4b) :: icvg
  integer(i4b) :: nr

  real(dp) :: pr, tr, vr
  real(dp) :: vc

  real(dp) :: pvt
  real(dp) :: ZZ

  real(dp) :: tr_inv,   &
              tr2_inv,  &
              tr3_inv,  &
              vr_inv,   &
              vr2_inv,  &
              vr4_inv,  &
              vr5_inv

  real(dp) :: term4,   term5,   &
              term6,   term7,   &
              term8,   term8_1, &
              term8_2, term8_3

! solve non-linear problem
  real(dp) :: jac,    &
              fxx,    &
              dxx,    &
              fxxinc, &
              vrinc

  real(dp) :: rgas_barL_Kmol
  
! parameters
  integer(i4b), parameter :: max=100
  real(dp), parameter     :: inc=1.0d-10
  
! parameters: critical pressure and temperature for CO2 (from Duan et al 1992b GCA)
  real(dp), parameter     :: pc = 73.825d0, &  !bar
                             tc = 304.20d0     !Kelvin  31.05d0 !!

! parameters: coefficients for equation A1 Duan and Sun (2005) Chem. Geol.  
  real(dp), dimension(15), parameter :: a = &

                                        (/ 8.99288497d-2,  &  ! a(1)  
                                          -4.94783127d-1,  &  ! a(2)
                                           4.77922245d-2,  &  ! a(3)
                                           1.03808883d-2,  &  ! a(4)
                                          -2.82516861d-2,  &  ! a(5)
                                           9.49887563d-2,  &  ! a(6)
                                           5.20600880d-4,  &  ! a(7)
                                          -2.93540971d-4,  &  ! a(8) 
                                          -1.77265112d-3,  &  ! a(9) 
                                          -2.51101973d-5,  &  ! a(10)
                                           8.93353441d-5,  &  ! a(11)
                                           7.88998563d-5,  &  ! a(12)
                                          -1.66727022d-2,  &  ! a(13)
                                           1.39800000d+0,  &  ! a(14)
                                           2.96000000d-2   /) ! a(15)

RunVerbose= .FALSE.

IF (RunVerbose) THEN

  write(*,*)'-------------------------------------------------------'
  write(*,*)
  write(*,*)'solving for Vr'
  write(*,*)
  write(*,*)'-------------------------------------------------------'
END IF

  ! solve non-linear problem for vr, and calculate ZZ = pr * vr / tr
  ! requires tk and pg, also: tc and pc

  rgas_barL_Kmol = rgas * 10.d0

  pr = pg / pc
  tr = tk / tc
  vc = rgas_barL_Kmol * tc / pc
  vr = vrInOut
  !!vr = 1.0d0 ! initial approx.

  icvg = 1
  do_NR: do nr=1,max
    
    fxx = ZZ(pr,tr,vr,a)
      
!    vol = vr * vc 
!    pvt = 
    
    fxx = fxx - pr * vr / tr
    
    IF (RunVerbose) THEN
      write(*,*)'NR iter ',nr,' reduced volume: ',vr, ' abs(residual): ',dabs(fxx)
    END IF

    if (dabs(fxx) < 1.0d-12) then
      icvg = 0
      vrInOut = vr
      exit do_NR
    end if
    
    vrinc = vr + inc
    
    fxxinc = ZZ(pr,tr,vrinc,a) - pr * vrinc / tr
    
    jac = (fxxinc - fxx) / inc
    
    if (dabs(jac) < 1.d-12 .or. jac /= jac) then
    
      write(*,*)'derivative problems in N-R for Vr ',jac
      stop
    
    end if
    
    dxx = - fxx / jac
        
    vr = vr + dxx

  end do do_NR 
  
  if (icvg == 1 .or. vr < 0.0d0) then
    write(*,*)'could not converge to find Vr, or Vr negative'
    write(*,*)'Vr= ',vr
    stop
  end if
  
  pvt = pr * vr / tr
  
  ! calculate fugacity
  ! requires z = pr vr /tr, vr and tr
  
  tr_inv = 1.0d0 / tr
  tr2_inv = tr_inv * tr_inv
  tr3_inv = tr2_inv * tr_inv
  
  vr_inv = 1.0d0 / vr
  vr2_inv = vr_inv * vr_inv
  vr4_inv = vr2_inv * vr2_inv
  vr5_inv = vr4_inv * vr_inv
  
  term4 = ( a(1) + a(2) * tr2_inv + a(3) * tr3_inv ) * vr_inv
  term5 = ( a(4) + a(5) * tr2_inv + a(6) * tr3_inv ) * vr2_inv / 2.0d0
  term6 = ( a(7) + a(8) * tr2_inv + a(9) * tr3_inv ) * vr4_inv / 4.0d0
  term7 = ( a(10) + a(11) * tr2_inv + a(12) * tr3_inv ) * vr5_inv / 5.0d0
  
  term8_1 = a(13) * tr3_inv / (2.0d0 * a(15))
  term8_2 = a(14) + 1.0d0 - ( a(14) + 1.0d0 + a(15) * vr2_inv ) * exp(-a(15) * vr2_inv) 
  term8_3 = exp(-a(15) * vr2_inv) ! only for checking NaN
  
  
  if (term8_1 /= term8_1 .or. &
      term8_2 /= term8_2 .or. &
      term8_3 /= term8_3) then
    write(*,*)'NaN or infinity in term8 of EOS'
    stop
  else   
    term8 = term8_1 * term8_2
  end if
  
  ln_fco2 = pvt       &
          - 1.0d0     &
          - dlog(pvt) &
          + term4     &
          + term5     &
          + term6     &
          + term7     &
          + term8

IF (RunVerbose) THEN
  write(*,*)
  write(*,*)
  write(*,*)'Vr    = ',vr
  write(*,*)'P(bar)= ',pg
  write(*,*)'T(K)  = ',tk
  write(*,*)'ln(fugacity coeff CO2)= ',ln_fco2
  write(*,*)'fugacity coeff CO2= ',exp(ln_fco2)
  write(*,*)
  write(*,*)'-------------------------------------------------------'
END IF
  
  return

end subroutine fugacity_co2

!-----------------------------------------------------------------------------------
!
! CO2 Fugacity - Duan and Sun 2006 Chem. Geol. 
!
!-----------------------------------------------------------------------------------
subroutine fugacity_co24(pg,tk,ln_fco2,vrInOut)

  use crunchtype
  use params

  implicit none

! arguments
  real(dp), intent(in)    :: pg
  real(dp), intent(in)    :: tk
  real(dp), intent(out)   :: ln_fco2
  real(dp), intent(inout)  :: vrInOut

! internal
  logical(lgt) :: RunVerbose

  real(dp)     :: fco2
  real(dp)     :: p1
  real(dp)     :: calc_pv
  real(dp)     :: phi
 
! parameters: critical pressure and temperature for CO2 (from Duan et al 1992b GCA)
  real(dp), parameter     :: pc = 73.825d0, &  !bar
                             tc = 304.20d0     !Kelvin  31.05d0 !!

! parameters: coefficients for equation 2 Duan et al (2006) Marine Chem.
  real(dp), dimension(15), parameter :: c1 = &

                                        (/ 1.0d0        ,  &  ! c(1)  
                                           4.75868350d-3,  &  ! c(2)
                                          -3.35699630d-6,  &  ! c(3)
                                           0.0d0        ,  &  ! c(4)
                                          -1.31793960d+0   ,  &  ! c(5)
                                          -3.83891010d-6,  &  ! c(6)
                                           0.0d0        ,  &  ! c(7)
                                           2.28151040d-3,  &  ! c(8) 
                                           0.0d0        ,  &  ! c(9) 
                                           0.0d0        ,  &  ! c(10)
                                           0.0d0        ,  &  ! c(11)
                                           0.0d0        ,  &  ! c(12)
                                           0.0d0        ,  &  ! c(13)
                                           0.0d0        ,  &  ! c(14)
                                           0.0d0           /) ! c(15)

  real(dp), dimension(15), parameter :: c2 = &

                                        (/-7.17348820d-1,  &  ! c(1)  
                                           1.59853790d-4,  &  ! c(2)
                                          -4.92864710d-7,  &  ! c(3)
                                           0.0d0        ,  &  ! c(4)
                                           0.0d0        ,  &  ! c(5)
                                          -2.78552850d-7,  &  ! c(6)
                                           1.18770150d-9,  &  ! c(7)
                                           0.0d0        ,  &  ! c(8) 
                                           0.0d0        ,  &  ! c(9) 
                                           0.0d0        ,  &  ! c(10)
                                           0.0d0        ,  &  ! c(11)
                                         -96.53951200d+0,  &  ! c(12)
                                           4.47749380d-1,  &  ! c(13)
                                         101.81078d+0   ,  &  ! c(14)
                                           5.37838790d-6   /) ! c(15)

  real(dp), dimension(15), parameter :: c3 = &

                                        (/-6.51290190d-2,  &  ! c(1)  
                                          -2.14299770d-4,  &  ! c(2)
                                          -1.14449300d-6,  &  ! c(3)
                                           0.0d0        ,  &  ! c(4)
                                           0.0d0        ,  &  ! c(5)
                                          -1.15580810d-7,  &  ! c(6)
                                           1.19523700d-9,  &  ! c(7)
                                           0.0d0        ,  &  ! c(8) 
                                           0.0d0        ,  &  ! c(9) 
                                           0.0d0        ,  &  ! c(10)
                                           0.0d0        ,  &  ! c(11)
                                        -221.34306d0    ,  &  ! c(12)
                                           0.0d0        ,  &  ! c(13)
                                          71.820393d0   ,  &  ! c(14)
                                           6.60892460d-6   /) ! c(15)

  real(dp), dimension(15), parameter :: c4 = &

                                        (/ 5.0383896d0  ,  &  ! c(1)  
                                          -4.42577440d-3,  &  ! c(2)
                                           0.0d0        ,  &  ! c(3)
                                           1.9572733d0  ,  &  ! c(4)
                                           0.0d0        ,  &  ! c(5)
                                           2.42234360d-6,  &  ! c(6)
                                           0.0d0        ,  &  ! c(7)
                                          -9.37961350d-4,  &  ! c(8) 
                                          -1.50260300d+0,  &  ! c(9) 
                                           3.02722400d-3,  &  ! c(10)
                                         -31.377342d0   ,  &  ! c(11)
                                         -12.847063d0   ,  &  ! c(12)
                                           0.0d0        ,  &  ! c(13)
                                           0.0d0        ,  &  ! c(14)
                                          -1.50566480d-5   /) ! c(15)

  real(dp), dimension(15), parameter :: c5 = &

                                        (/-1.60631520d+1,  &  ! c(1)  
                                          -2.70579900d-3,  &  ! c(2)
                                           0.0d0        ,  &  ! c(3)
                                           1.41192390d-1,  &  ! c(4)
                                           0.0d0        ,  &  ! c(5)
                                           8.11329650d-7,  &  ! c(6)
                                           0.0d0        ,  &  ! c(7)
                                          -1.14530820d-4,  &  ! c(8) 
                                           2.38956710d+0,  &  ! c(9) 
                                           5.05274570d-4,  &  ! c(10)
                                         -17.76346000d+0,  &  ! c(11)
                                         985.92232d0    ,  &  ! c(12)
                                           0.0d0        ,  &  ! c(13)
                                           0.0d0        ,  &  ! c(14)
                                          -5.49652560d-7   /) ! c(15)

  real(dp), dimension(15), parameter :: c6 = &

                                        (/-1.56934900d-1,  &  ! c(1)  
                                           4.46214070d-4,  &  ! c(2)
                                          -9.10805910d-7,  &  ! c(3)
                                           0.0d0        ,  &  ! c(4)
                                           0.0d0        ,  &  ! c(5)
                                           1.06473990d-7,  &  ! c(6)
                                           2.4273357d-10,  &  ! c(7)
                                           0.0d0        ,  &  ! c(8) 
                                           3.58742550d-1,  &  ! c(9) 
                                           6.33197100d-5,  &  ! c(10)
                                        -249.89661d0    ,  &  ! c(11)
                                           0.0d0        ,  &  ! c(12)
                                           0.0d0        ,  &  ! c(13)
                                         888.768d0      ,  &  ! c(14)
                                          -6.63480030d-7   /) ! c(15)
  
  RunVerbose= .FALSE.

  ! calculate p1
  if (tk < 305.d0) then
    
    p1 = calc_pv(tk) ! saturation pressure

  else if (tk < 405.d0) then

    p1 = 75.d0 + (tk -305.d0) * 1.25d0

  else 

    p1 = 200.d0 ! bar
  
  end if
  

  if (pg < p1 .and. tk > 273.d0 .and. tk < 573.d0 ) then
    fco2 = phi(pg,tk,c1)

! pg > p1 
  else if (pg < 1000.d0 .and. tk > 273.d0 .and. tk < 340.d0 ) then

    fco2 = phi(pg,tk,c2)

  else if (pg < 1000.d0 .and. tk > 340.d0 .and. tk < 435.d0 ) then

    fco2 = phi(pg,tk,c4)

! pg > 1000 bar
  else if (pg > 1000.d0 .and. tk > 273.d0 .and. tk < 340.d0 ) then

    fco2 = phi(pg,tk,c3)

  else if (pg > 1000.d0 .and. tk > 340.d0 .and. tk < 435.d0 ) then

    fco2 = phi(pg,tk,c5)
  
  else if (pg > 1000.d0 .and. tk < 573.d0 ) then

    fco2 = phi(pg,tk,c6)

  else 

    write(*,*)'fugacity model error: pg or tk exceeds range: ', pg, tk
    stop

  end if  

   ln_fco2 = log(fco2)

end subroutine fugacity_co24


!-----------------------------------------------------------------------------------
!
! saturated pv(CO2)  -option 1  http://en.wikipedia.org/wiki/Carbon_dioxide_data & CHERIC
!                    -option 2  http://ddbonline.ddbst.com/AntoineCalculation/AntoineCalculationCGI.exe
!
!-----------------------------------------------------------------------------------
function calc_pv(tt) 

  use crunchtype
  
  implicit none
  
  real(dp) :: calc_pv
  real(dp) :: tt ! kelvin
  real(dp) :: ln_pv

  real(dp) :: term1,   term2,   &
              term3,   term4,   &
              term5,   term6,   &
              term7
  
  real(dp), parameter :: mmhg_to_bar = 0.00133322368
    
  !term1 = log( 760.0d0 / 101.325d0 ) 
  !term2 = 24.03761d0 * log(tt)
  !term3 = - 7062.404d0 / tt
  !term4 = 166.3861d0
  !term5 = 3.368548d-5 * tt * tt
  
  !ln_pv = term1 + &
  !        term2 + &
  !        term3 + &
  !        term4 + &
  !        term5
  
  !calc_pv = mmhg_to_bar * exp(ln_pv) 

  calc_pv = mmhg_to_bar * 10**( 7.5322 - 835.06 / ( 268.223 + tt - 273.15) ) 


  return
       
end function calc_pv

!-----------------------------------------------------------------------------------
!
! EOS for fugacity - Duan and Sun 2006 Marine Chemistry
!
!-----------------------------------------------------------------------------------
function phi(pp,tt,cc) 

  use crunchtype
  
  implicit none
  
  real(dp) :: phi
  real(dp) :: tt, pp ! kelvin, bar
  real(dp), dimension(15) :: cc

  real(dp) :: term1,   term2,   &
              term3,   term4,   &
              term5,   term6,   &
              term7
  
  term1 = cc(1)
  term2 = ( cc(2) + cc(3) * tt + cc(4) / tt + cc(5) / (tt - 150.0d0) ) * pp
  term3 = ( cc(6) + cc(7) * tt + cc(8) / tt ) * pp * pp 
  term4 = ( cc(9) + cc(10) * tt + cc(11) / tt ) * log(pp) 
  term5 = ( cc(12) + cc(13) * tt ) / pp
  term6 = cc(14) / tt 
  term7 = cc(15) * tt * tt
  
  phi = term1 + &
        term2 + &
        term3 + &
        term4 + &
        term5 + &
        term6 + &
        term7

  return
       
end function phi

!-----------------------------------------------------------------------------------
!
! EOS for supercritical CO2 - Duan and Sun 2003 Chem. Geol. - from Duan et al 1992 GCA
!
!-----------------------------------------------------------------------------------
function ZZ(pr,tr,vr,a) 

  use crunchtype
  
  implicit none
  
  real(dp) :: ZZ
  real(dp) :: pr, tr, vr
  real(dp), dimension(15) :: a
  
  real(dp) :: tr_inv,   &
              tr2_inv,  &
              tr3_inv,  &
              vr_inv,   &
              vr2_inv,  &
              vr4_inv,  &
              vr5_inv

  real(dp) :: term2,   term3,   &
              term4,   term5,   &
              term6,   term6_1, &
              term6_2, term6_3
  
  ! calculate Z = pr vr /tr
  
  tr_inv = 1.0d0 / tr
  tr2_inv = tr_inv * tr_inv
  tr3_inv = tr2_inv * tr_inv
  
  vr_inv = 1.0d0 / vr
  vr2_inv = vr_inv * vr_inv
  vr4_inv = vr2_inv * vr2_inv
  vr5_inv = vr4_inv * vr_inv
  
  term2 = ( a(1) + a(2) * tr2_inv + a(3) * tr3_inv ) * vr_inv
  term3 = ( a(4) + a(5) * tr2_inv + a(6) * tr3_inv ) * vr2_inv
  term4 = ( a(7) + a(8) * tr2_inv + a(9) * tr3_inv ) * vr4_inv
  term5 = ( a(10) + a(11) * tr2_inv + a(12) * tr3_inv ) * vr5_inv
  
  term6_1 = a(13) * tr3_inv * vr2_inv
  term6_2 = a(14) + a(15) * vr2_inv
  term6_3 = exp(-a(15) * vr2_inv)
  
  if (term6_1 /= term6_1 .or. &
      term6_2 /= term6_2 .or. &
      term6_3 /= term6_3) then
    write(*,*)'NaN or infinity in term6 of EOS'
    stop
  else   
    term6 = term6_1 * term6_2 * term6_3
  end if
  
  ZZ = 1.0d0 + &
       term2 + &
       term3 + &
       term4 + &
       term5 + &
       term6

   return
       
end function ZZ

!-----------------------------------------------------------------------------------
!
! Lambda coefficients - Pitzer interaction parameters for CO2 with ions (Na)
!
!-----------------------------------------------------------------------------------
subroutine calc_ph2o(tk,ph2o)

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

end subroutine calc_ph2o

