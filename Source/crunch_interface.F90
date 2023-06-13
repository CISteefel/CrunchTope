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
    
MODULE crunch_interface

INTERFACE
  SUBROUTINE CourantStepAlt(nx,ny,nz,dtmaxcour)
  USE crunchtype
  IMPLICIT NONE    
  INTEGER(I4B), INTENT(IN)                                 :: nx
  INTEGER(I4B), INTENT(IN)                                 :: ny
  INTEGER(I4B), INTENT(IN)                                 :: nz
  REAL(DP), INTENT(OUT)                                    :: dtmaxcour
  END SUBROUTINE CourantStepAlt
END INTERFACE

INTERFACE
  SUBROUTINE ludcmp90(a,indx,d,n)
  USE crunchtype
  REAL(DP), DIMENSION(:,:), intent(in out)                   :: a
  INTEGER(I4B), DIMENSION(:), intent(out)                    :: indx
  INTEGER(I4B), INTENT(IN)                                   :: n
  REAL(DP), intent(out)                                      :: d
  END SUBROUTINE ludcmp90
END INTERFACE

INTERFACE
  SUBROUTINE lubksb90(a,indx,b,n)
  USE crunchtype
  IMPLICIT NONE
  REAL(DP),  DIMENSION(:,:), INTENT(IN)                          :: a
  REAL(DP),  DIMENSION(:), INTENT(IN OUT)                        :: b
  INTEGER(I4B),  DIMENSION(:),INTENT(IN)                         :: indx
  INTEGER(I4B), INTENT(IN)                                       :: n
  END SUBROUTINE lubksb90
END INTERFACE

INTERFACE
  SUBROUTINE tridag_ser(a,b,c,r,u)
  USE crunchtype
  IMPLICIT NONE
  REAL(DP), DIMENSION(:), INTENT(IN) :: a,b,c,r
  REAL(DP), DIMENSION(:), INTENT(OUT) :: u
  END SUBROUTINE tridag_ser
END INTERFACE

INTERFACE
  FUNCTION locate(xx,x)
  USE crunchtype
  IMPLICIT NONE
  REAL(DP), DIMENSION(:), INTENT(IN) :: xx
  REAL(DP), INTENT(IN) :: x
  INTEGER(I4B) :: locate
  END FUNCTION locate
END INTERFACE

!!INTERFACE
!!SUBROUTINE assemble(nx,ny,nz,ncomp,nspec,nkin,nrct,ngas,ikin,  &
!!    nexchange,nexch_sec,nsurf,nsurf_sec,npot,ndecay,nn,delt,   &
!!    user,amatpetsc,blockfl,blockl,blockm,blockr,blockfr)
!!  USE crunchtype
!!  IMPLICIT NONE
!!  INTEGER(I4B), INTENT(IN)                      :: nx
!!  INTEGER(I4B), INTENT(IN)                      :: ny
!!  INTEGER(I4B), INTENT(IN)                      :: nz
!!  INTEGER(I4B), INTENT(IN)                      :: ncomp
!!  INTEGER(I4B), INTENT(IN)                      :: nspec
!!  INTEGER(I4B), INTENT(IN)                      :: nkin
!!  INTEGER(I4B), INTENT(IN)                      :: nrct
!!  INTEGER(I4B), INTENT(IN)                      :: ngas
!!  INTEGER(I4B), INTENT(IN)                      :: ikin
!!  INTEGER(I4B), INTENT(IN)                      :: nexchange
!!  INTEGER(I4B), INTENT(IN)                      :: nexch_sec
!!  INTEGER(I4B), INTENT(IN)                      :: nsurf
!!  INTEGER(I4B), INTENT(IN)                      :: nsurf_sec
!!  INTEGER(I4B), INTENT(IN)                      :: npot
!!  INTEGER(I4B), INTENT(IN)                      :: ndecay
!!  INTEGER(I4B), INTENT(IN)                      :: nn
!!  REAL(DP), INTENT(IN)                          :: delt

!!!********************** Add in PETSC declarations for f90 variables *****
!!REAL(DP), DIMENSION(ncomp+nexchange+nsurf+npot,ncomp+nexchange+nsurf+npot) :: blockfl
!!REAL(DP), DIMENSION(ncomp+nexchange+nsurf+npot,ncomp+nexchange+nsurf+npot) :: blockl
!!REAL(DP), DIMENSION(ncomp+nexchange+nsurf+npot,ncomp+nexchange+nsurf+npot) :: blockm
!!REAL(DP), DIMENSION(ncomp+nexchange+nsurf+npot,ncomp+nexchange+nsurf+npot) :: blockr
!!REAL(DP), DIMENSION(ncomp+nexchange+nsurf+npot,ncomp+nexchange+nsurf+npot) :: blockfr
!!!************************* End PETSc declarations ****************************
!!  END SUBROUTINE assemble
!!END INTERFACE

INTERFACE
  SUBROUTINE readModFlowParameters(nout,nchem,nparams,jxTemp,jyTemp,jzTemp,  &
              conditionNum,modflowstring,lenstring)
  USE crunchtype
  USE strings
  IMPLICIT NONE
  INTEGER(I4B), INTENT(IN)                                    :: nout
  INTEGER(I4B), INTENT(IN)                                    :: nchem
  INTEGER(I4B), INTENT(OUT)                                   :: nparams
  INTEGER(I4B), DIMENSION(:),INTENT(OUT)                      :: jxTemp
  INTEGER(I4B), DIMENSION(:),INTENT(OUT)                      :: jyTemp
  INTEGER(I4B), DIMENSION(:),INTENT(OUT)                      :: jzTemp
  INTEGER(I4B), DIMENSION(:),INTENT(OUT)                      :: conditionNum
  CHARACTER (LEN=mls), INTENT(IN)                             :: modflowstring
  INTEGER(I4B), INTENT(IN)                                    :: lenstring
  END SUBROUTINE readModFlowParameters
END INTERFACE

INTERFACE
  FUNCTION imaxloc(arr)
  USE crunchtype
  REAL(DP), DIMENSION(:), INTENT(IN) :: arr
  INTEGER(I4B) :: imaxloc
  END FUNCTION imaxloc
END INTERFACE

INTERFACE
  SUBROUTINE SurfaceConcentration(jx,jy,jz,ncomp,nsurf,nsurf_sec)
  USE crunchtype
  IMPLICIT NONE
  INTEGER(I4B), INTENT(IN)                                    :: jx
  INTEGER(I4B), INTENT(IN)                                    :: jy
  INTEGER(I4B), INTENT(IN)                                    :: jz
  INTEGER(I4B), INTENT(IN)                                    :: ncomp
  INTEGER(I4B), INTENT(IN)                                    :: nsurf
  INTEGER(I4B), INTENT(IN)                                    :: nsurf_sec
  END SUBROUTINE SurfaceConcentration
END INTERFACE

INTERFACE
  SUBROUTINE ExchangeConcentration(jx,jy,jz,ncomp,nexchange,nexch_sec)
  USE crunchtype
  IMPLICIT NONE
  INTEGER(I4B), INTENT(IN)                                    :: jx
  INTEGER(I4B), INTENT(IN)                                    :: jy
  INTEGER(I4B), INTENT(IN)                                    :: jz
  INTEGER(I4B), INTENT(IN)                                    :: ncomp
  INTEGER(I4B), INTENT(IN)                                    :: nexchange
  INTEGER(I4B), INTENT(IN)                                    :: nexch_sec
  END SUBROUTINE ExchangeConcentration
END INTERFACE

END MODULE crunch_interface
