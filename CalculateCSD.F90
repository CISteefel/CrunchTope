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
    
subroutine CalculateCSD(jx,jy,jz,nrct,ncomp,delt)
USE crunchtype
USE params
USE temperature
USE mineral
USE crunch_interface
USE NanoCrystal

IMPLICIT NONE

!  External arrays and variables

INTEGER(I4B), INTENT(IN)                      :: jx
INTEGER(I4B), INTENT(IN)                      :: jy
INTEGER(I4B), INTENT(IN)                      :: jz
INTEGER(I4B), INTENT(IN)                      :: nrct
INTEGER(I4B), INTENT(IN)                      :: ncomp
REAL(DP), INTENT(IN)                          :: delt               

! Internal arrays and variables

REAL(DP)                                      :: dxCSD
REAL(DP)                                      :: rinv
REAL(DP)                                      :: fe
REAL(DP)                                      :: fw
REAL(DP)                                      :: ae
REAL(DP)                                      :: aw
REAL(DP)                                      :: apx
REAL(DP)                                      :: pi
REAL(DP)                                      :: Tk
REAL(DP)                                      :: VolumeSingleMolecule
REAL(DP)                                      :: gCritical
REAL(DP)                                      :: qNucleation
REAL(DP)                                      :: qNucleationPrefactor
REAL(DP)                                      :: CriticalRadius
REAL(DP)                                      :: CriticalRadiusCheck
REAL(DP)                                      :: RatioVolumeBoltzmann
REAL(DP)                                      :: dt
REAL(DP),PARAMETER                                      :: AvogadroNumber=6.02214E23
REAL(DP),PARAMETER                                      :: BoltzmannConstant=1.36E-23

INTEGER(I4B)                                  :: l
INTEGER(I4B)                                  :: k
INTEGER(I4B)                                  :: nn


!! Crystal Size Distribution calculation
    
    dt = delt/1.0
    rinv = 1.0/dt
    Tk = t(jx,jy,jz) + 273.15

!!    LinearGrowthRate = 0.0d0

    pi = DACOS(-1.0d0)
    qNucleationPrefactor = 1.0E05
    sigma = 0.0d0
    sigma(3) = 200.0/1000.0
    sigma(4) = 90.0/1000.0

!!  Or calculate the volume of a single molecule from the molar volume divided by Avogradro's number
!! (then multiply by number of units in a cluster)
    
    dxCSD = 2.0E-06
    radius(1) = 4.0E-06
    DO l = 2,nCSD
      radius(l) = radius(l-1) + dxCSD
    END DO

    DO k=1,nrct
 
      IF (CrystalSizeDistribution(k)) THEN

!!       nCrystal(20,k,jx,jy,jz) = 100

!!        LinearGrowthRate = 1.D-8

        VolumeSingleMolecule = volmol(k)/AvogadroNumber
        CALL satcalc(ncomp,nrct,jx,jy,jz)
        RatioVolumeBoltzmann = VolumeSingleMolecule/BoltzmannConstant
        gCritical = (pi*sigma(k)**3 * RatioVolumeBoltzmann*RatioVolumeBoltzmann )/          &
                ( 3.0d0* (Tk*silog(1,k)*clg)**2 )
        qNucleation = qNucleationPrefactor*DEXP (-gCritical/(BoltzmannConstant*Tk) )
        NucleationScaleFactor = qNucleationPrefactor*100000.0
        NucleationScaleFactor = 1.0d0
        qNucleation = qNucleation/NucleationScaleFactor

        DO l = 1,nCSD
          xtvd(l) = DFLOAT(nCrystal(l,k,jx,jy,jz))
        END DO

        fw = LinearGrowthRate(1,k,jx,jy,jz)
        fe = 0.5d0*(LinearGrowthRate(2,k,jx,jy,jz) + LinearGrowthRate(1,k,jx,jy,jz))
        aw = DMAX1(fw,0.0D0)
        ae = DMAX1(-fe,0.0D0)
        apx = DMAX1(-fw,0.0D0) + DMAX1(fe,0.0D0)
        aaCSD(1) = -aw*dt
        bbCSD(1) = apx*dt + dxCSD
        ccCSD(1) = -ae*dt
        uuCSD(1) = xtvd(1)
        rrCSD(1) = xtvd(1)*dxCSD + dt*qNucleation*dxCSD

        fw = 0.5d0*(LinearGrowthRate(nCSD-1,k,jx,jy,jz)+LinearGrowthRate(nCSD,k,jx,jy,jz))
        fe = LinearGrowthRate(nCSD,k,jx,jy,jz)
        aw = DMAX1(fw,0.0D0)
        ae = DMAX1(-fe,0.0D0)
        apx = DMAX1(-fw,0.0D0) + DMAX1(fe,0.0D0)
        aaCSD(nCSD) = -aw*dt
        bbCSD(nCSD) = apx*dt + dxCSD
        ccCSD(nCSD) = -ae*dt
        uuCSD(nCSD) = xtvd(nCSD)
        rrCSD(nCSD) = xtvd(nCSD)*dxCSD

        DO l = 2,nCSD-1

!!        Define fluxes at boundaries
          fw = 0.5d0*(LinearGrowthRate(l-1,k,jx,jy,jz)+LinearGrowthRate(l,k,jx,jy,jz))
          fe = 0.5d0*(LinearGrowthRate(l+1,k,jx,jy,jz)+LinearGrowthRate(l,k,jx,jy,jz))
          aw = DMAX1(fw,0.0D0)
          ae = DMAX1(-fe,0.0D0)
          apx = DMAX1(-fw,0.0D0) + DMAX1(fe,0.0D0)

          aaCSD(l) = -aw*dt
          bbCSD(l) = apx*dt + dxCSD
          ccCSD(l) = -ae*dt
          uuCSD(l) = xtvd(l)
          rrCSD(l) = xtvd(l)*dxCSD

        END DO
        
        CALL tridag_ser(aaCSD,bbCSD,ccCSD,rrCSD,uuCSD)      
        
        DO l = 1,nCSD
          nCrystal(l,k,jx,jy,jz) = NINT(uuCSD(l))
        END DO

      END IF
          
    END DO   ! End of loop for CSD calculation


RETURN
END SUBROUTINE CalculateCSD