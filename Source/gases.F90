!!! *** Copyright Notice ***
!!! �CrunchFlow�, Copyright (c) 2016, The Regents of the University of California, through Lawrence Berkeley National Laboratory 
!!! (subject to receipt of any required approvals from the U.S. Dept. of Energy).� All rights reserved.
!!!�
!!! If you have questions about your rights to use or distribute this software, please contact 
!!! Berkeley Lab's Innovation & Partnerships Office at��IPO@lbl.gov.
!!!�
!!! NOTICE.� This Software was developed under funding from the U.S. Department of Energy and the U.S. Government 
!!! consequently retains certain rights. As such, the U.S. Government has been granted for itself and others acting 
!!! on its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the Software to reproduce, distribute copies to the public, 
!!! prepare derivative works, and perform publicly and display publicly, and to permit other to do so.
!!!
!!! *** License Agreement ***
!!! �CrunchFlow�, Copyright (c) 2016, The Regents of the University of California, through Lawrence Berkeley National Laboratory)
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

SUBROUTINE gases(ncomp,ngas,jx,jy,jz)
USE crunchtype
USE params
USE concentration
USE temperature
USE runtime, ONLY: Duan,Duan2006

IMPLICIT NONE

!  External variables

INTEGER(I4B), INTENT(IN)                                   :: ncomp
INTEGER(I4B), INTENT(IN)                                   :: ngas
INTEGER(I4B), INTENT(IN)                                   :: jx 
INTEGER(I4B), INTENT(IN)                                   :: jy
INTEGER(I4B), INTENT(IN)                                   :: jz
!  Internal variables

REAL(DP)                                                   :: tempk
REAL(DP)                                                   :: denmol
REAL(DP)                                                   :: sum
REAL(DP)                                                   :: pg
REAL(DP)                                                   :: ln_fco2
REAL(DP)                                                   :: vrInOut

INTEGER(I4B)                                               :: i
INTEGER(I4B)                                               :: kk

  !!! Added July 17 by Carl (hopefully not stomped on)
  
ln_fco2 = 0.0d0

IF (Duan .OR. Duan2006) THEN
  pg = GasPressureTotal(jx,jy,jz)
END IF

tempk = t(jx,jy,jz) + 273.15

!!denmol = LOG(1.e05/(8.314*tempk))   ! P/RT = n/V, with pressure converted from bars to Pascals
denmol = DLOG( (1.0E05) /(8.314d0*tempk) )   ! P/RT = n/V, with pressure converted from bars to Pascals
!!!denmol = LOG( (1.0E05) /(8.314d0*283.15 ) )  ! P/RT = n/V, with pressure converted from bars to Pascals

!!  NOTE:  The "denmol" should convert to mol/m*3 (n/V)

DO kk = 1,ngas

    sum = 0.0
    DO i = 1,ncomp
      sum = sum + mugas(kk,i)*(sp(i,jx,jy,jz) + gam(i,jx,jy,jz))
    END DO

  ln_fco2 = 0.0d0  ! fugacity coefficient for CO2(g)
  IF (Duan) THEN
    ln_fco2 = 0.0d0  ! fugacity coefficient for CO2(g)
    if (namg(kk) == 'CO2(g)') then
      vrInOut = vrSave(jx,jy,jz)
      call fugacity_co2(pg,tempk,ln_fco2,vrInOut)
      vrSave(jx,jy,jz) = vrInOut
    end if
  ELSE IF (Duan2006) THEN
    ln_fco2 = 0.0d0  ! fugacity coefficient for CO2(g)
    if (namg(kk) == 'CO2(g)') then
      vrInOut = vrSave(jx,jy,jz)
      call fugacity_co24(pg,tempk,ln_fco2,vrInOut)
      vrSave(jx,jy,jz) = vrInOut
    end if
  END IF

!! Basically, first two terms on RHS give you the mole fraction, then multipled by n/V to gives mol/m^3

  spgas(kk,jx,jy,jz) = keqgas(kk,jx,jy,jz) + sum + denmol - ln_fco2
  spgas10(kk,jx,jy,jz) = DEXP(spgas(kk,jx,jy,jz))  ! mol/m**3

END DO

RETURN
END SUBROUTINE gases
!  **************************************************************
