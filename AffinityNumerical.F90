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
    
SUBROUTINE AffinityNumerical(ncomp,nrct,jx,jy,jz,np,k,sppTMP,AffinityTerm,time)
USE crunchtype
USE params
USE runtime, ONLY: JennyDruhan,UseBulkMineral,LagSurface
USE concentration, ONLY: gam,sn,sppTMP10,ulab,ikh2o
USE mineral
USE medium
USE temperature
USE isotope

IMPLICIT NONE

INTEGER(I4B), INTENT(IN)                                        :: ncomp
INTEGER(I4B), INTENT(IN)                                        :: nrct
INTEGER(I4B), INTENT(IN)                                        :: jx
INTEGER(I4B), INTENT(IN)                                        :: jy
INTEGER(I4B), INTENT(IN)                                        :: jz
INTEGER(I4B), INTENT(IN)                                        :: np
INTEGER(I4B), INTENT(IN)                                        :: k

REAL(DP),DIMENSION(:),INTENT(IN)                                :: sppTMP
REAL(DP), INTENT(OUT)                                           :: AffinityTerm
REAL(DP), INTENT(IN)                                            :: time

!  Internal variables and arrays

REAL(DP)                                                        :: sumiap
REAL(DP)                                                        :: silogTMP
REAL(DP)                                                        :: silnTMP
REAL(DP)                                                        :: snormTMP
REAL(DP)                                                        :: siTMP
REAL(DP)                                                        :: power
REAL(DP)                                                        :: term1
REAL(DP)                                                        :: sign
!!!REAL(DP)                                                        :: MoleFraction40
!!!REAL(DP)                                                        :: MoleFraction44
!!!REAL(DP)                                                        :: MoleFraction32
!!!REAL(DP)                                                        :: MoleFraction34
!!!REAL(DP)                                                        :: MoleFraction44Mineral
!!!REAL(DP)                                                        :: MoleFraction40Mineral


INTEGER(I4B)                                                    :: i
INTEGER(I4B)                                                    :: kk
INTEGER(I4B)                                                    :: npp

REAL(DP)                                                        :: MoleFractionMineral
REAL(DP)                                                        :: Denominator

INTEGER(I4B)                                                    :: kIsotopologue
INTEGER(I4B)                                                    :: Isotopologue
INTEGER(I4B)                                                    :: kMineralRare
INTEGER(I4B)                                                    :: kMineralCommon
INTEGER(I4B)                                                    :: iPrimaryRare
INTEGER(I4B)                                                    :: iPrimaryCommon

!! DummyComment

IF (kcrossaff(np,k) /= 0) THEN
  kk = kcrossaff(np,k)
  npp = 1
ELSE
  kk = k
  npp = np
END IF


  sumiap = 0.0D0
  DO i = 1,ncomp
    sumiap = sumiap + mumin(npp,kk,i)*(sppTMP(i)+gam(i,jx,jy,jz))
  END DO


IF (nIsotopeMineral > 0) THEN

  IF (IsotopeMineralRare(k)) THEN

    kIsotopologue = kPointerIsotope(k)

    kMineralRare = kIsotopeRare(kIsotopologue)
    KMineralCommon = kIsotopeCommon(kIsotopologue)
    isotopologue = PointerToPrimaryIsotope(kIsotopologue)
    iPrimaryCommon = isotopeCommon(Isotopologue)

    IF (isotopeBackReactionOption(kIsotopologue) == 'none' .OR. UseAqueousMoleFraction(kIsotopologue)) THEN
      isotopologue = PointerToPrimaryIsotope(kIsotopologue)
      iPrimaryRare = isotopeRare(Isotopologue)
      iPrimaryCommon = isotopeCommon(Isotopologue)
      denominator = sppTMP10(iPrimaryRare) + sppTMP10(iPrimaryCommon)
      MoleFractionMineral = (sppTMP10(iPrimaryRare)/denominator)
    ELSE
      MoleFractionMineral = MoleFractionMineralRare(kPointerIsotope(k))
    END IF
    sumiap = sumiap - (mumin(1,kMineralCommon,iPrimaryCommon))*DLOG(MoleFractionMineral)

  ELSE IF (IsotopeMineralCommon(k)) THEN

    kIsotopologue = kPointerIsotope(k)

    kMineralRare = kIsotopeRare(kIsotopologue)
    KMineralCommon = kIsotopeCommon(kIsotopologue)
    isotopologue = PointerToPrimaryIsotope(kIsotopologue)
    iPrimaryCommon = isotopeCommon(Isotopologue)

    IF (isotopeBackReactionOption(kIsotopologue) == 'none' .OR. UseAqueousMoleFraction(kIsotopologue)) THEN
      isotopologue = PointerToPrimaryIsotope(kIsotopologue)
      iPrimaryRare = isotopeRare(Isotopologue)
      iPrimaryCommon = isotopeCommon(Isotopologue)
      denominator = sppTMP10(iPrimaryRare) + sppTMP10(iPrimaryCommon)
      MoleFractionMineral = (sppTMP10(iPrimaryCommon)/denominator)
    ELSE
      MoleFractionMineral = MoleFractionMineralCommon(kPointerIsotope(k))
    END IF
    sumiap = sumiap - (mumin(1,kMineralCommon,iPrimaryCommon))*DLOG(MoleFractionMineral)

  ELSE
    CONTINUE
  END IF

END IF   !! Block where nIsotopeMineral > 0

silogTMP = (sumiap - keqmin(npp,kk,jx,jy,jz))/clg
siTMP = 10d0**(silogTMP)
silnTMP = silogTMP*clg

IF (siTMP > 1.0D0) THEN
  sign = 1.0D0
ELSE
  sign = -1.0D0
END IF
  
!!CIS 03/19/07  IF (AffinityDepend2(np,k) /= 1.0d0) THEN
!!CIS 03/19/07    snorm = siTMP**AffinityDepend2(np,k)
!!CIS 03/19/07  ELSE    
!!CIS 03/19/07    snorm = siTMP
!!CIS 03/19/07  END IF     

!!CIS ******* 03/19/07  ******************  
IF (AffinityDepend2(npp,kk) == 1.0D0 .and. AffinityDepend3(npp,kk) == 1.0d0) THEN
  snormTMP = siTMP
ELSE IF (AffinityDepend2(npp,kk) /= 1.0d0 .AND. AffinityDepend3(npp,kk) == 1.0d0) THEN
  snormTMP = siTMP**AffinityDepend2(npp,kk)
ELSE
  snormTMP = DEXP( -AffinityDepend2(npp,kk)*(DABS(silnTMP))**AffinityDepend3(npp,kk) )
END IF
!!CIS ******* 03/19/07  ****************** 

IF (AffinityDepend1(npp,kk) == 1.0D0) THEN
  term1 = sign*DABS(snormTMP - 1.0D0)
ELSE
  term1 = sign*DABS(snormTMP - 1.0D0)**(AffinityDepend1(npp,kk))
END IF

IF (imintype(npp,kk) == 4) THEN                               !! Precipitation only
  AffinityTerm = MAX(0.0d0,term1)
ELSE IF (imintype(npp,kk) == 5) THEN                          !! Dissolution only
  AffinityTerm = MIN(0.0d0,term1)
ELSE
  AffinityTerm = term1
END IF
   
RETURN
END SUBROUTINE AffinityNumerical
!********************************************************
