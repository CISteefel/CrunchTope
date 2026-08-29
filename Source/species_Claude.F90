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


SUBROUTINE species(ncomp,nspec,nsurf,nexchange,npot,nx,ny,nz)
USE crunchtype
USE params
USE concentration
USE runtime

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: ncomp
INTEGER(I4B), INTENT(IN)                                    :: nspec
INTEGER(I4B), INTENT(IN)                                    :: nsurf
INTEGER(I4B), INTENT(IN)                                    :: nexchange
INTEGER(I4B), INTENT(IN)                                    :: npot
INTEGER(I4B), INTENT(IN)                                    :: nx
INTEGER(I4B), INTENT(IN)                                    :: ny
INTEGER(I4B), INTENT(IN)                                    :: nz

!  Internal variables and arrays

INTEGER(I4B)                                                :: ksp
INTEGER(I4B)                                                :: icomp
INTEGER(I4B)                                                :: nk
INTEGER(I4B)                                                :: jx
INTEGER(I4B)                                                :: jy
INTEGER(I4B)                                                :: jz


INTEGER(I4B)                                                :: ik
INTEGER(I4B)                                                :: neqn

CHARACTER (LEN=3)                                           :: ulabPrint

INTEGER(I4B)                                                :: kk
REAL(DP)                                                    :: sum
REAL(DP)                                                    :: sumderiv_gamma
INTEGER(I4B)                                                :: pos_der

neqn = ncomp + nexchange + nsurf +npot + 1 + 1
posLNval = neqn + 1
posval   = neqn + 2
posgammawater = ncomp + nexchange + nsurf +npot + 1 + 1
posIonS = ncomp + nexchange + nsurf + npot + 1 
      
DO jz = 1,nz
  DO jy = 1,ny
    DO jx = 1,nx
      
      deriv_conc(:,:,jx,jy,jz) = 0.d0

      ! Water (ik = 1)
      deriv_conc(1,1,jx,jy,jz) = 1D0

      ! Primary species other than H20
      DO ik = 2,ncomp
      ! derivatives
        deriv_conc(ik,ik,jx,jy,jz) = sp10(ik,jx,jy,jz)
      END DO

      ! Secondary species
      !!!DO ik = ncomp+1, ncomp+nspec
      DO ksp = 1,nspec
        ik = ksp + ncomp
        ! Contribution of water to equilibrium, icomp = 1
  
        !!!   f_ConcAqueous(ik,posLNval) = f_ConcAqueous(ik,posLNval) + GloPar%AqSpecies(ik,1)%stoichio * ln_primary(scalar%posGammaWater) 
        !!! First, to ln concentration of secondary species (contribution of water)
        sp(ik,jx,jy,jz) = sp(ik,jx,jy,jz) + muaq(ksp,1) * lnGammaWater(jx,jy,jz)  
  
        !!! Now derivative with respect to GammaWater 
        !!! f_ConcAqueous(ik,scalar%posGammaWater) = f_ConcAqueous(ik,scalar%posGammaWater) + GloPar%AqSpecies(ik,1)%stoichio 
        deriv_conc(ik, posGammaWater,jx,jy,jz) = deriv_conc(ik, posGammaWater,jx,jy,jz) + muaq(ksp,1)
  
        !!! Ionic strength
        !!! sion(jx,jy,jz)   !!! defined as primary variable ncomp+nexchange+nsurf+npot+1
        !!! Definition of ionic strength
        !!! SumCharge = 0.d0
        !!! DO ik = 1,ncomp+nspec
        !!!   SumCharge = SumCharge + chg(ik)*chg(ik)*sp10(ik,jx,jy,jz)
        !!! END DO
        !!! sion(jx,jy,jz) = 0.5*SumCharge
  
        !!! Now derivative with respect to ionic strength (posIonS)
        !!! deriv_conc(ik,posIonS,jx,jy,jz) = 0.d0
  
        ! Contribution of other primary species
        DO icomp = 2,ncomp
  
         !!!   f_ConcAqueous(ik,posLNval) = f_ConcAqueous(ik,posLNval) &
         !!!                                + GloPar%AqSpecies(ik,icomp)%stoichio * ( f_ConcAqueous(icomp,posLNval) + GammaAqueous(icomp,posLNval) )
	        sp(ik,jx,jy,jz) = sp(ik,jx,jy,jz) + muaq(ksp,icomp) * ( sp(icomp,jx,jy,jz) + lngamma(icomp,jx,jy,jz) )
	
           !!! f_ConcAqueous(ik,1:scalar%neqn) = f_ConcAqueous(ik,1:scalar%neqn)                                       &
           !!!                                         + GloPar%AqSpecies(ik,icomp)%stoichio                                   &
           !!!                                         * ( f_ConcAqueous(icomp,1:scalar%neqn)/f_ConcAqueous(icomp,posval)  & 
           !!! 										   + GammaAqueous(icomp,1:scalar%neqn) / GammaAqueous(icomp,posval) )
           deriv_conc(ik,1:neqn,jx,jy,jz) = deriv_conc(ik,1:neqn,jx,jy,jz) + muaq(ksp,icomp)               &
	                                              * (deriv_conc(icomp,1:neqn,jx,jy,jz)/sp10(icomp,jx,jy,jz)     &
										          + deriv_gamma(icomp,1:neqn,jx,jy,jz) / gamma(icomp,jx,jy,jz) )
										                                                                                                                                       
        END DO

        !!! f_ConcAqueous(ik,posLNval) = f_ConcAqueous(ik,posLNval) + EquilConstAqueous(ik,posLNval) - GammaAqueous(ik,posLNval) 
        sp(ik,jx,jy,jz) = sp(ik,jx,jy,jz) + keqaq(ksp,jx,jy,jz) - lngamma(ik,jx,jy,jz)
                                    
        !!! f_ConcAqueous(ik,1:scalar%neqn) = f_ConcAqueous(ik,1:scalar%neqn) - GammaAqueous(ik,1:scalar%neqn)/GammaAqueous(ik,posval) 
        deriv_conc(ik,1:neqn,jx,jy,jz) = deriv_conc(ik,1:neqn,jx,jy,jz) - deriv_gamma(ik,1:neqn,jx,jy,jz)/gamma(ik,jx,jy,jz)
  
        !!! f_ConcAqueous(ik,posval) = EXP(f_ConcAqueous(ik,posLNval))
        sp10(ik,jx,jy,jz) = EXP( sp(ik,jx,jy,jz) )
  
        !!! f_ConcAqueous(ik,1:scalar%neqn) = f_ConcAqueous(ik,1:scalar%neqn) * f_ConcAqueous(ik,posval)
        deriv_conc(ik,1:neqn,jx,jy,jz) = deriv_conc(ik,1:neqn,jx,jy,jz) * sp10(ik,jx,jy,jz)
                                      
      END DO
      
    !!! End of loop of nx, ny, and nz
    END DO 
  END DO
END DO

RETURN
END SUBROUTINE species
!  **************************************************************
