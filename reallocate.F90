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

    
SUBROUTINE REALLOCATE(ncomp,nspec,nrct,nkin,ngas,nsurf,nexchange,ikin,nexch_sec,nsurf_sec)
USE crunchtype
USE params
USE runtime
USE concentration
USE mineral
USE solver
USE medium
USE transport
USE flow
USE temperature
USE io
USE strings
USE ReadFlow
USE modflowModule

IMPLICIT NONE

!! EXTERNAL VARIABLES AND ARRAYS

INTEGER(I4B), INTENT(IN)                             :: ncomp
INTEGER(I4B), INTENT(IN)                             :: nspec
INTEGER(I4B), INTENT(IN)                             :: nrct
INTEGER(I4B), INTENT(IN)                             :: nkin
INTEGER(I4B), INTENT(IN)                             :: ngas
INTEGER(I4B), INTENT(IN)                             :: nsurf
INTEGER(I4B), INTENT(IN)                             :: nexchange
INTEGER(I4B), INTENT(IN)                             :: ikin
INTEGER(I4B), INTENT(IN)                             :: nexch_sec
INTEGER(I4B), INTENT(IN)                             :: nsurf_sec
!!INTEGER(I4B), INTENT(IN)                             :: nchem
!!INTEGER(I4B), INTENT(IN)                             :: ninhibitaqmax
!!INTEGER(I4B), INTENT(IN)                             :: nmonodaqmax
!!INTEGER(I4B), INTENT(IN)                             :: nreactmax
!!INTEGER(I4B), INTENT(IN)                             :: ntot
!!INTEGER(I4B), INTENT(IN)                             :: nreactkinmax
!!INTEGER(I4B), INTENT(IN)                             :: ninhibitmax

!! Internal variables and arrays

REAL(DP), DIMENSION(:), ALLOCATABLE                           :: work1
REAL(DP), DIMENSION(:,:), ALLOCATABLE                         :: work2
REAL(DP), DIMENSION(:,:,:), ALLOCATABLE                       :: work3
REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE                     :: work4

INTEGER(I4B), DIMENSION(:), ALLOCATABLE                       :: workint1
INTEGER(I4B), DIMENSION(:,:), ALLOCATABLE                     :: workint2
INTEGER(I4B), DIMENSION(:,:,:), ALLOCATABLE                   :: workint3

CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE                :: workchar1
CHARACTER (LEN=mls), DIMENSION(:,:), ALLOCATABLE              :: workchar2
CHARACTER (LEN=mls), DIMENSION(:,:,:), ALLOCATABLE            :: workchar3

LOGICAL(LGT), DIMENSION(:), ALLOCATABLE                       :: worklogical1
LOGICAL(LGT), DIMENSION(:,:), ALLOCATABLE                     :: worklogical2
LOGICAL(LGT), DIMENSION(:,:,:), ALLOCATABLE                   :: worklogical3

INTEGER(I4B)                                                  :: ndim1
INTEGER(I4B)                                                  :: ndim2
INTEGER(I4B)                                                  :: ndim3
INTEGER(I4B)                                                  :: ndim4
INTEGER(I4B)                                                  :: i,j,k,nbig,ntot,ncount,np


!  **************  ReALLOCATE 1D arrays  *****************************

!  **********  Logical arrays  *****************

i = size(LocalEquilibrium,1)
ALLOCATE(worklogical1(i))
worklogical1 = LocalEquilibrium
DEALLOCATE(LocalEquilibrium)
ALLOCATE(LocalEquilibrium(nrct))
IF (nrct /= 0) LocalEquilibrium(1:nrct) = worklogical1(1:nrct)
DEALLOCATE(worklogical1)

!  **********  Integer arrays  *****************

i = size(iedl,1)
ALLOCATE(workint1(i))
workint1 = iedl
DEALLOCATE(iedl)
ALLOCATE(iedl(nsurf))
IF (nsurf /= 0) iedl(1:nsurf) = workint1(1:nsurf)
DEALLOCATE(workint1)

i = size(ivolume,1)
ALLOCATE(workint1(i))
workint1 = ivolume
DEALLOCATE(ivolume)
ALLOCATE(ivolume(nrct))
IF (nrct /= 0) ivolume(1:nrct) = workint1(1:nrct)
DEALLOCATE(workint1)

i = size(iplot,1)
ALLOCATE(workint1(i))
workint1 = iplot
DEALLOCATE(iplot)
ALLOCATE(iplot(ncomp+nspec))
IF(ncomp+nspec /= 0) iplot(1:ncomp+nspec) = workint1(1:ncomp+nspec)
DEALLOCATE(workint1)

i = size(nreactmin,1)
ALLOCATE(workint1(i))
workint1 = nreactmin
DEALLOCATE(nreactmin)
ALLOCATE(nreactmin(nrct))
IF(nrct /= 0)nreactmin(1:nrct) = workint1(1:nrct)
DEALLOCATE(workint1)

i = size(lenmin,1)
ALLOCATE(workint1(i))
workint1 = lenmin
DEALLOCATE(lenmin)
ALLOCATE(lenmin(nrct))
IF(nrct /= 0)lenmin(1:nrct) = workint1(1:nrct)
DEALLOCATE(workint1)

i = size(ksurf,1)
ALLOCATE(workint1(i))
workint1 = ksurf
DEALLOCATE(ksurf)
ALLOCATE(ksurf(nsurf))
IF(nsurf /= 0)ksurf(1:nsurf) = workint1(1:nsurf)
DEALLOCATE(workint1)

i = size(kexch,1)
ALLOCATE(workint1(i))
workint1 = kexch
DEALLOCATE(kexch)
ALLOCATE(kexch(nexchange))
IF(nexchange /= 0)kexch(1:nexchange) = workint1(1:nexchange)
DEALLOCATE(workint1)

i = size(nreactkin,1)
ALLOCATE(workint1(i))
workint1 = nreactkin
DEALLOCATE(nreactkin)
ALLOCATE(nreactkin(ikin))
IF(ikin /= 0)nreactkin(1:ikin) = workint1(1:ikin)
DEALLOCATE(workint1)

i = size(ixlink,1)
ALLOCATE(workint1(i))
workint1 = ixlink
DEALLOCATE(ixlink)
ALLOCATE(ixlink(nexch_sec))
IF(nexch_sec /= 0)ixlink(1:nexch_sec) = workint1(1:nexch_sec)
DEALLOCATE(workint1)

i = size(nclink,1)
ALLOCATE(workint1(i))
workint1 = nclink
DEALLOCATE(nclink)
ALLOCATE(nclink(nexch_sec))
IF(nexch_sec /= 0)nclink(1:nexch_sec) = workint1(1:nexch_sec)
DEALLOCATE(workint1)

i = size(iexchange,1)
ALLOCATE(workint1(i))
workint1 = iexchange
DEALLOCATE(iexchange)
ALLOCATE(iexchange(nexchange))
IF(nexchange /= 0)iexchange(1:nexchange) = workint1(1:nexchange)
DEALLOCATE(workint1)

i = size(iaqtype,1)
ALLOCATE(workint1(i))
workint1 = iaqtype
DEALLOCATE(iaqtype)
ALLOCATE(iaqtype(ikin))
IF(ikin /= 0)iaqtype(1:ikin) = workint1(1:ikin)
DEALLOCATE(workint1)

i = size(nmonodaq,1)
ALLOCATE(workint1(i))
workint1 = nmonodaq
DEALLOCATE(nmonodaq)
ALLOCATE(nmonodaq(ikin))
IF(ikin /= 0)nmonodaq(1:ikin) = workint1(1:ikin)
DEALLOCATE(workint1)

! biomass
ndim1 = nreactmax
ndim2 = nrct
i = size(chi_min,1)
j = size(chi_min,2)
ALLOCATE(workint2(i,j))
workint2 = chi_min
DEALLOCATE(chi_min)
ALLOCATE(chi_min(ndim1,ndim2))
IF(ndim1 /= 0 .AND. ndim2 /= 0) chi_min(1:ndim1,1:ndim2) = workint2(1:ndim1,1:ndim2)
DEALLOCATE(workint2)

ndim1 = nreactmax
ndim2 = nrct
i = size(direction_min,1)
j = size(direction_min,2)
ALLOCATE(workint2(i,j))
workint2 = direction_min
DEALLOCATE(direction_min)
ALLOCATE(direction_min(ndim1,ndim2))
IF(ndim1 /= 0 .AND. ndim2 /= 0) direction_min(1:ndim1,1:ndim2) = workint2(1:ndim1,1:ndim2)
DEALLOCATE(workint2)

ndim1 = ikin
!!ndim2 = nrct
i = size(chi_kin,1)
!!j = size(chi_kin,2)
ALLOCATE(workint1(i))
workint1 = chi_kin
DEALLOCATE(chi_kin)
ALLOCATE(chi_kin(ndim1))
IF(ndim1 /= 0) chi_kin(1:ndim1) = workint1(1:ndim1)
DEALLOCATE(workint1)

ndim1 = ikin
!!ndim2 = nrct
i = size(direction_kin,1)
!!j = size(direction_kin,2)
ALLOCATE(workint1(i))
workint1 = direction_kin
DEALLOCATE(direction_kin)
ALLOCATE(direction_kin(ndim1))
IF(ndim1 /= 0) direction_kin(1:ndim1) = workint1(1:ndim1)
DEALLOCATE(workint1)
! biomass end

i = size(ninhibitaq,1)
ALLOCATE(workint1(i))
workint1 = ninhibitaq
DEALLOCATE(ninhibitaq)
ALLOCATE(ninhibitaq(ikin))
IF(ikin /= 0)ninhibitaq(1:ikin) = workint1(1:ikin)
DEALLOCATE(workint1)


!  *********  Real arrays  ********************

i = size(volmol,1)
ALLOCATE(work1(i))
work1 = volmol
DEALLOCATE(volmol)
ALLOCATE(volmol(nrct))
IF(nrct /= 0)volmol(1:nrct) = work1(1:nrct)
DEALLOCATE(work1)

i = size(wtmin,1)
ALLOCATE(work1(i))
work1 = wtmin
DEALLOCATE(wtmin)
ALLOCATE(wtmin(nrct))
IF(nrct /= 0)wtmin(1:nrct) = work1(1:nrct)
DEALLOCATE(work1)

i = size(wtcomp,1)
ALLOCATE(work1(i))
work1 = wtcomp
DEALLOCATE(wtcomp)
ALLOCATE(wtcomp(ncomp))
wtcomp(1:ncomp) = work1(1:ncomp)
DEALLOCATE(work1)

i = size(alnk,1)
ALLOCATE(work1(i))
work1 = alnk
DEALLOCATE(alnk)
ALLOCATE(alnk(nreactmax*nrct))
IF(nreactmax*nrct /= 0)alnk(1:nreactmax*nrct) = work1(1:nreactmax*nrct)
DEALLOCATE(work1)

i = size(acmp,1)
ALLOCATE(work1(i))
work1 = acmp
DEALLOCATE(acmp)
ALLOCATE(acmp(ncomp+nspec))
IF(ncomp+nspec /= 0)acmp(1:ncomp+nspec) = work1(1:ncomp+nspec)
DEALLOCATE(work1)

i = size(chg,1)
ALLOCATE(work1(i))
work1 = chg
DEALLOCATE(chg)
ALLOCATE(chg(ncomp+nspec))
IF(ncomp+nspec /= 0)chg(1:ncomp+nspec) = work1(1:ncomp+nspec)
DEALLOCATE(work1)

i = size(bdotParameter,1)
ALLOCATE(work1(i))
work1 = bdotParameter
DEALLOCATE(bdotParameter)
ALLOCATE(bdotParameter(ncomp+nspec))
IF(ncomp+nspec /= 0) bdotParameter(1:ncomp+nspec) = work1(1:ncomp+nspec)
DEALLOCATE(work1)

i = size(wtaq,1)
ALLOCATE(work1(i))
work1 = wtaq
DEALLOCATE(wtaq)
ALLOCATE(wtaq(ncomp+nspec))
IF(ncomp+nspec /= 0)wtaq(1:ncomp+nspec) = work1(1:ncomp+nspec)
DEALLOCATE(work1)

nbig = ncomp+nspec+ngas+nrct+nsurf_sec

i = size(alogkp,1)
ALLOCATE(work1(i))
work1 = alogkp
DEALLOCATE(alogkp)
ALLOCATE(alogkp(nbig))
IF(nbig /= 0)alogkp(1:nbig) = work1(1:nbig)
DEALLOCATE(work1)

i = size(satkin,1)
ALLOCATE(work1(i))
work1 = satkin
DEALLOCATE(satkin)
ALLOCATE(satkin(ikin))
IF(ikin /= 0)satkin(1:ikin) = work1(1:ikin)
DEALLOCATE(work1)

i = size(keqexc,1)
ALLOCATE(work1(i))
work1 = keqexc
DEALLOCATE(keqexc)
ALLOCATE(keqexc(nexch_sec))
IF(nexch_sec /= 0)keqexc(1:nexch_sec) = work1(1:nexch_sec)
DEALLOCATE(work1)

i = size(bfit,1)
ALLOCATE(work1(i))
work1 = bfit
DEALLOCATE(bfit)
ALLOCATE(bfit(nexch_sec))
IF(nexch_sec /= 0)bfit(1:nexch_sec) = work1(1:nexch_sec)
DEALLOCATE(work1)

i = size(keqkin,1)
ALLOCATE(work1(i))
work1 = keqkin
DEALLOCATE(keqkin)
ALLOCATE(keqkin(ikin))
IF(ikin /= 0)keqkin(1:ikin) = work1(1:ikin)
DEALLOCATE(work1)

i = size(TotChargeSave,1)
ALLOCATE(work1(i))
work1 = TotChargeSave
DEALLOCATE(TotChargeSave)
ALLOCATE(TotChargeSave(nchem))
IF(nchem /= 0) TotChargeSave(1:nchem) = work1(1:nchem)
DEALLOCATE(work1)

i = size(SolidDensity,1)
ALLOCATE(work1(i))
work1 = SolidDensity
DEALLOCATE(SolidDensity)
ALLOCATE(SolidDensity(nchem))
IF(nchem /= 0) SolidDensity(1:nchem) = work1(1:nchem)
DEALLOCATE(work1)

!  ***********  Allocate 1D character arrays  ****************

i = size(condlabel,1)
ALLOCATE(workchar1(i))
workchar1 = condlabel
DEALLOCATE(condlabel)
ALLOCATE(condlabel(nchem))
IF(nchem /= 0) condlabel(1:nchem) = workchar1(1:nchem)
DEALLOCATE(workchar1)

i = size(namg,1)
ALLOCATE(workchar1(i))
workchar1 = namg
DEALLOCATE(namg)
ALLOCATE(namg(ngas))
IF(ngas /= 0)namg(1:ngas) = workchar1(1:ngas)
DEALLOCATE(workchar1)

i = size(umin,1)
ALLOCATE(workchar1(i))
workchar1 = umin
DEALLOCATE(umin)
ALLOCATE(umin(nrct))
IF(nrct /= 0)umin(1:nrct) = workchar1(1:nrct)
DEALLOCATE(workchar1)

i = size(ulab,1)
ALLOCATE(workchar1(i))
workchar1 = ulab
DEALLOCATE(ulab)
ALLOCATE(ulab(ncomp+nspec))
IF(ncomp+nspec /= 0)ulab(1:ncomp+nspec) = workchar1(1:ncomp+nspec)
DEALLOCATE(workchar1)

i = size(namkin,1)
ALLOCATE(workchar1(i))
workchar1 = namkin
DEALLOCATE(namkin)
ALLOCATE(namkin(ikin))
IF(ikin /= 0)namkin(1:ikin) = workchar1(1:ikin)
DEALLOCATE(workchar1)

i = size(namexc,1)
ALLOCATE(workchar1(i))
workchar1 = namexc
DEALLOCATE(namexc)
ALLOCATE(namexc(nexchange))
IF(nexchange /= 0)namexc(1:nexchange) = workchar1(1:nexchange)
DEALLOCATE(workchar1)

i = size(nam_exchsec,1)
ALLOCATE(workchar1(i))
workchar1 = nam_exchsec
DEALLOCATE(nam_exchsec)
ALLOCATE(nam_exchsec(nexch_sec))
IF(nexch_sec /= 0)nam_exchsec(1:nexch_sec) = workchar1(1:nexch_sec)
DEALLOCATE(workchar1)

i = size(namsurf,1)
ALLOCATE(workchar1(i))
workchar1 = namsurf
DEALLOCATE(namsurf)
ALLOCATE(namsurf(nsurf))
IF(nsurf /= 0)namsurf(1:nsurf) = workchar1(1:nsurf)
DEALLOCATE(workchar1)

i = size(namsurf_sec,1)
ALLOCATE(workchar1(i))
workchar1 = namsurf_sec
DEALLOCATE(namsurf_sec)
ALLOCATE(namsurf_sec(nsurf_sec))
IF(nsurf_sec /= 0)namsurf_sec(1:nsurf_sec) = workchar1(1:nsurf_sec)
DEALLOCATE(workchar1)


!  ****************  End of 1D reallocation  ********************************

!  **************  ReALLOCATE 2D arrays  *****************************

!  *********  Integer arrays  **************

ntot = ncomp+nspec

ndim1 = ninhibitaqmax
ndim2 = ikin
i = size(itot_inhibitaq,1)
j = size(itot_inhibitaq,2)
ALLOCATE(workint2(i,j))
workint2 = itot_inhibitaq
DEALLOCATE(itot_inhibitaq)
ALLOCATE(itot_inhibitaq(ndim1,ndim2))
IF(ndim1 /= 0 .AND. ndim2 /= 0)  &
    itot_inhibitaq(1:ndim1,1:ndim2) = workint2(1:ndim1,1:ndim2)
DEALLOCATE(workint2)

ndim1 = nmonodaqmax
ndim2 = ikin
i = size(itot_monodaq,1)
j = size(itot_monodaq,2)
ALLOCATE(workint2(i,j))
workint2 = itot_monodaq
DEALLOCATE(itot_monodaq)
ALLOCATE(itot_monodaq(ndim1,ndim2))
IF(ndim1 /= 0 .AND. ndim2 /= 0)itot_monodaq(1:ndim1,1:ndim2) = workint2(1:ndim1,1:ndim2)
DEALLOCATE(workint2)

ndim1 = ninhibitaqmax
ndim2 = ikin
i = size(inhibitaq,1)
j = size(inhibitaq,2)
ALLOCATE(workint2(i,j))
workint2 = inhibitaq
DEALLOCATE(inhibitaq)
ALLOCATE(inhibitaq(ndim1,ndim2))
IF(ndim1 /= 0 .AND. ndim2 /= 0)inhibitaq(1:ndim1,1:ndim2) = workint2(1:ndim1,1:ndim2)
DEALLOCATE(workint2)

ndim1 = nmonodaqmax
ndim2 = ikin
i = size(imonodaq,1)
j = size(imonodaq,2)
ALLOCATE(workint2(i,j))
workint2 = imonodaq
DEALLOCATE(imonodaq)
ALLOCATE(imonodaq(ndim1,ndim2))
IF(ndim1 /= 0 .AND. ndim2 /= 0)imonodaq(1:ndim1,1:ndim2) = workint2(1:ndim1,1:ndim2)
DEALLOCATE(workint2)

ndim1 = nreactmax
ndim2 = nrct
i = size(imintype,1)
j = size(imintype,2)
ALLOCATE(workint2(i,j))
workint2 = imintype
DEALLOCATE(imintype)
ALLOCATE(imintype(ndim1,ndim2))
IF(ndim1 /= 0 .AND. ndim2 /= 0)imintype(1:ndim1,1:ndim2) = workint2(1:ndim1,1:ndim2)
DEALLOCATE(workint2)

ndim1 = nreactmax
ndim2 = nrct
i = size(ndepend,1)
j = size(ndepend,2)
ALLOCATE(workint2(i,j))
workint2 = ndepend
DEALLOCATE(ndepend)
ALLOCATE(ndepend(ndim1,ndim2))
IF(ndim1 /= 0 .AND. ndim2 /= 0)ndepend(1:ndim1,1:ndim2) = workint2(1:ndim1,1:ndim2)
DEALLOCATE(workint2)

!!  ********** Hyperbolic Inhibition ******************
ndim1 = nreactmax
ndim2 = nrct
i = size(HyperbolicInhibition,1)
j = size(HyperbolicInhibition,2)
ALLOCATE(worklogical2(i,j))
worklogical2 = HyperbolicInhibition
DEALLOCATE(HyperbolicInhibition)
ALLOCATE(HyperbolicInhibition(ndim1,ndim2))
IF(ndim1 /= 0 .AND. ndim2 /= 0) HyperbolicInhibition(1:ndim1,1:ndim2) = worklogical2(1:ndim1,1:ndim2)
DEALLOCATE(worklogical2)

ndim1 = nreactmax
ndim2 = nrct
i = size(HyperbolicInhibitionPointer,1)
j = size(HyperbolicInhibitionPointer,2)
ALLOCATE(workint2(i,j))
workint2 = HyperbolicInhibitionPointer
DEALLOCATE(HyperbolicInhibitionPointer)
ALLOCATE(HyperbolicInhibitionPointer(ndim1,ndim2))
IF(ndim1 /= 0 .AND. ndim2 /= 0) HyperbolicInhibitionPointer(1:ndim1,1:ndim2) = workint2(1:ndim1,1:ndim2)
DEALLOCATE(workint2)

ndim1 = nreactmax
ndim2 = nrct
i = size(Kformation,1)
j = size(Kformation,2)
ALLOCATE(work2(i,j))
work2 = Kformation
DEALLOCATE(Kformation)
ALLOCATE(Kformation(ndim1,ndim2))
IF(ndim1 /= 0 .AND. ndim2 /= 0) Kformation(1:ndim1,1:ndim2) = work2(1:ndim1,1:ndim2)
DEALLOCATE(work2)

ndim1 = nreactmax
ndim2 = nrct
i = size(HyperbolicInhibitionDepend,1)
j = size(HyperbolicInhibitionDepend,2)
ALLOCATE(work2(i,j))
work2 = HyperbolicInhibitionDepend
DEALLOCATE(HyperbolicInhibitionDepend)
ALLOCATE(HyperbolicInhibitionDepend(ndim1,ndim2))
IF(ndim1 /= 0 .AND. ndim2 /= 0) HyperbolicInhibitionDepend(1:ndim1,1:ndim2) = work2(1:ndim1,1:ndim2)
DEALLOCATE(work2)

ndim1 = nreactmax
ndim2 = nrct
i = size(HyperbolicInhibitionName,1)
j = size(HyperbolicInhibitionName,2)
ALLOCATE(workchar2(i,j))
workchar2 = HyperbolicInhibitionName
DEALLOCATE(HyperbolicInhibitionName)
ALLOCATE(HyperbolicInhibitionName(ndim1,ndim2))
IF(ndim1 /= 0 .AND. ndim2 /= 0) HyperbolicInhibitionName(1:ndim1,1:ndim2) = workchar2(1:ndim1,1:ndim2)
DEALLOCATE(workchar2)
!!  ********** Hyperbolic Inhibition ******************

ndim1 = nreactmax
ndim2 = nrct
i = size(nmonod,1)
j = size(nmonod,2)
ALLOCATE(workint2(i,j))
workint2 = nmonod
DEALLOCATE(nmonod)
ALLOCATE(nmonod(ndim1,ndim2))
IF(ndim1 /= 0 .AND. ndim2 /= 0)nmonod(1:ndim1,1:ndim2) = workint2(1:ndim1,1:ndim2)
DEALLOCATE(workint2)

ndim1 = nreactmax
ndim2 = nrct
i = size(ninhibit,1)
j = size(ninhibit,2)
ALLOCATE(workint2(i,j))
workint2 = ninhibit
DEALLOCATE(ninhibit)
ALLOCATE(ninhibit(ndim1,ndim2))
IF(ndim1 /= 0 .AND. ndim2 /= 0)ninhibit(1:ndim1,1:ndim2) = workint2(1:ndim1,1:ndim2)
DEALLOCATE(workint2)

ndim1 = nreactmax
ndim2 = nrct
i = size(kcrossaff,1)
j = size(kcrossaff,2)
ALLOCATE(workint2(i,j))
workint2 = kcrossaff
DEALLOCATE(kcrossaff)
ALLOCATE(kcrossaff(ndim1,ndim2))
IF(ndim1 /= 0 .AND. ndim2 /= 0)kcrossaff(1:ndim1,1:ndim2) = workint2(1:ndim1,1:ndim2)
DEALLOCATE(workint2)

ndim1 = nrct
ndim2 = nchem
i = size(iarea,1)
j = size(iarea,2)
ALLOCATE(workint2(i,j))
workint2 = iarea
DEALLOCATE(iarea)
ALLOCATE(iarea(ndim1,ndim2))
IF(ndim1 /= 0 .AND. ndim2 /= 0) iarea(1:ndim1,1:ndim2) = workint2(1:ndim1,1:ndim2)
DEALLOCATE(workint2)


!  *********  Real arrays  **************

ndim1 = nrct
ndim2 = nchem
i = size(volin,1)
j = size(volin,2)
ALLOCATE(work2(i,j))
work2 = volin
DEALLOCATE(volin)
ALLOCATE(volin(ndim1,ndim2))
IF(ndim1 /= 0 .AND. ndim2 /= 0)volin(1:ndim1,1:ndim2) = work2(1:ndim1,1:ndim2)
DEALLOCATE(work2)

ndim1 = nrct
ndim2 = nchem
i = size(MineralMoles,1)
j = size(MineralMoles,2)
ALLOCATE(work2(i,j))
work2 = MineralMoles
DEALLOCATE(MineralMoles)
ALLOCATE(MineralMoles(ndim1,ndim2))
IF(ndim1 /= 0 .AND. ndim2 /= 0)MineralMoles(1:ndim1,1:ndim2) = work2(1:ndim1,1:ndim2)
DEALLOCATE(work2)

ndim1 = nrct
ndim2 = nchem
i = size(areain,1)
j = size(areain,2)
ALLOCATE(work2(i,j))
work2 = areain
DEALLOCATE(areain)
ALLOCATE(areain(ndim1,ndim2))
IF(ndim1 /= 0 .AND. ndim2 /= 0)areain(1:ndim1,1:ndim2) = work2(1:ndim1,1:ndim2)
DEALLOCATE(work2)

ndim1 = nsurf
ndim2 = nchem
i = size(site_density,1)
j = size(site_density,2)
ALLOCATE(work2(i,j))
work2 = site_density
DEALLOCATE(site_density)
ALLOCATE(site_density(ndim1,ndim2))
IF(ndim1 /= 0 .AND. ndim2 /= 0) site_density(1:ndim1,1:ndim2) = work2(1:ndim1,1:ndim2)
DEALLOCATE(work2)

ndim1 = nrct
ndim2 = nchem
i = size(specific,1)
j = size(specific,2)
ALLOCATE(work2(i,j))
work2 = specific
DEALLOCATE(specific)
ALLOCATE(specific(ndim1,ndim2))
IF(ndim1 /= 0 .AND. ndim2 /= 0) specific(1:ndim1,1:ndim2) = work2(1:ndim1,1:ndim2)
DEALLOCATE(work2)

ndim1 = nrct
ndim2 = nchem
i = size(voltemp,1)
j = size(voltemp,2)
ALLOCATE(work2(i,j))
work2 = voltemp
DEALLOCATE(voltemp)
ALLOCATE(voltemp(ndim1,ndim2))
IF(ndim1 /= 0 .AND. ndim2 /= 0) voltemp(1:ndim1,1:ndim2) = work2(1:ndim1,1:ndim2)
DEALLOCATE(work2)

ndim1 = ninhibitaqmax
ndim2 = ikin
i = size(rinhibitaq,1)
j = size(rinhibitaq,2)
ALLOCATE(work2(i,j))
work2 = rinhibitaq
DEALLOCATE(rinhibitaq)
ALLOCATE(rinhibitaq(ndim1,ndim2))
IF(ndim1 /= 0 .AND. ndim2 /= 0)rinhibitaq(1:ndim1,1:ndim2) = work2(1:ndim1,1:ndim2)
DEALLOCATE(work2)

ndim1 = ikin
ndim2 = ncomp
i = size(mukin,1)
j = size(mukin,2)
ALLOCATE(work2(i,j))
work2 = mukin
DEALLOCATE(mukin)
ALLOCATE(mukin(ikin,ndim2))
IF(ndim1 /= 0 .AND. ndim2 /= 0)mukin(1:ndim1,1:ndim2) = work2(1:ndim1,1:ndim2)
DEALLOCATE(work2)

ndim1 = nmonodaqmax
ndim2 = ikin
i = size(halfsataq,1)
j = size(halfsataq,2)
ALLOCATE(work2(i,j))
work2 = halfsataq
DEALLOCATE(halfsataq)
ALLOCATE(halfsataq(ndim1,ndim2))
IF(ndim1 /= 0 .AND. ndim2 /= 0) halfsataq(1:ndim1,1:ndim2) = work2(1:ndim1,1:ndim2)
DEALLOCATE(work2)

ndim1 = ntot
ndim2 = nchem
i = size(spcond,1)
j = size(spcond,2)
ALLOCATE(work2(i,j))
work2 = spcond
DEALLOCATE(spcond)
ALLOCATE(spcond(ndim1,ndim2))
IF(ndim1 /= 0 .AND. ndim2 /= 0)spcond(1:ndim1,1:ndim2) = work2(1:ndim1,1:ndim2)
DEALLOCATE(work2)

ndim1 = ncomp
ndim2 = nchem
i = size(scond,1)
j = size(scond,2)
ALLOCATE(work2(i,j))
work2 = scond
DEALLOCATE(scond)
ALLOCATE(scond(ndim1,ndim2))
IF(ndim1 /= 0 .AND. ndim2 /= 0)scond(1:ndim1,1:ndim2) = work2(1:ndim1,1:ndim2)
DEALLOCATE(work2)

ndim1 = ntot
ndim2 = nchem
i = size(spcond10,1)
j = size(spcond10,2)
ALLOCATE(work2(i,j))
work2 = spcond10
DEALLOCATE(spcond10)
ALLOCATE(spcond10(ndim1,ndim2))
IF(ndim1 /= 0 .AND. ndim2 /= 0)spcond10(1:ndim1,1:ndim2) = work2(1:ndim1,1:ndim2)
DEALLOCATE(work2)

ndim1 = ngas
ndim2 = nchem
i = size(spcondgas,1)
j = size(spcondgas,2)
ALLOCATE(work2(i,j))
work2 = spcondgas
DEALLOCATE(spcondgas)
ALLOCATE(spcondgas(ndim1,ndim2))
IF(ndim1 /= 0 .AND. ndim2 /= 0)spcondgas(1:ndim1,1:ndim2) = work2(1:ndim1,1:ndim2)
DEALLOCATE(work2)

ndim1 = ngas
ndim2 = nchem
i = size(spcondgas10,1)
j = size(spcondgas10,2)
ALLOCATE(work2(i,j))
work2 = spcondgas10
DEALLOCATE(spcondgas10)
ALLOCATE(spcondgas10(ndim1,ndim2))
IF(ndim1 /= 0 .AND. ndim2 /= 0)spcondgas10(1:ndim1,1:ndim2) = work2(1:ndim1,1:ndim2)
DEALLOCATE(work2)

ndim1 = nexchange+nexch_sec
ndim2 = nchem
i = size(spcondex,1)
j = size(spcondex,2)
ALLOCATE(work2(i,j))
work2 = spcondex
DEALLOCATE(spcondex)
ALLOCATE(spcondex(ndim1,ndim2))
IF(ndim1 /= 0 .AND. ndim2 /= 0)spcondex(1:ndim1,1:ndim2) = work2(1:ndim1,1:ndim2)
DEALLOCATE(work2)

ndim1 = nexchange+nexch_sec
ndim2 = nchem
i = size(spcondex10,1)
j = size(spcondex10,2)
ALLOCATE(work2(i,j))
work2 = spcondex10
DEALLOCATE(spcondex10)
ALLOCATE(spcondex10(ndim1,ndim2))
IF(ndim1 /= 0 .AND. ndim2 /= 0)spcondex10(1:ndim1,1:ndim2) = work2(1:ndim1,1:ndim2)
DEALLOCATE(work2)

ndim1 = nsurf+nsurf_sec
ndim2 = nchem
i = size(spcondsurf,1)
j = size(spcondsurf,2)
ALLOCATE(work2(i,j))
work2 = spcondsurf
DEALLOCATE(spcondsurf)
ALLOCATE(spcondsurf(ndim1,ndim2))
IF(ndim1 /= 0 .AND. ndim2 /= 0)spcondsurf(1:ndim1,1:ndim2) = work2(1:ndim1,1:ndim2)
DEALLOCATE(work2)

ndim1 = nsurf
ndim2 = nchem
i = size(LogPotentialInit,1)
j = size(LogPotentialInit,2)
ALLOCATE(work2(i,j))
work2 = LogPotentialInit
DEALLOCATE(LogPotentialInit)
ALLOCATE(LogPotentialInit(ndim1,ndim2))
IF(ndim1 /= 0 .AND. ndim2 /= 0)LogPotentialInit(1:ndim1,1:ndim2) = work2(1:ndim1,1:ndim2)
DEALLOCATE(work2)


ndim1 = nsurf+nsurf_sec
ndim2 = nchem
i = size(spcondsurf10,1)
j = size(spcondsurf10,2)
ALLOCATE(work2(i,j))
work2 = spcondsurf10
DEALLOCATE(spcondsurf10)
ALLOCATE(spcondsurf10(ndim1,ndim2))
IF(ndim1 /= 0 .AND. ndim2 /= 0)spcondsurf10(1:ndim1,1:ndim2) = work2(1:ndim1,1:ndim2)
DEALLOCATE(work2)

ndim1 = nspec
ndim2 = ncomp
i = size(muaq,1)
j = size(muaq,2)
ALLOCATE(work2(i,j))
work2 = muaq
DEALLOCATE(muaq)
ALLOCATE(muaq(ndim1,ndim2))
IF(ndim1 /= 0 .AND. ndim2 /= 0)muaq(1:ndim1,1:ndim2) = work2(1:ndim1,1:ndim2)
DEALLOCATE(work2)

ndim1 = ngas
ndim2 = ncomp
i = size(mugas,1)
j = size(mugas,2)
ALLOCATE(work2(i,j))
work2 = mugas
DEALLOCATE(mugas)
ALLOCATE(mugas(ndim1,ndim2))
IF(ndim1 /= 0 .AND. ndim2 /= 0)mugas(1:ndim1,1:ndim2) = work2(1:ndim1,1:ndim2)
DEALLOCATE(work2)

ndim1 = nsurf_sec
ndim2 = ncomp+nsurf
i = size(musurf,1)
j = size(musurf,2)
ALLOCATE(work2(i,j))
work2 = musurf
DEALLOCATE(musurf)
ALLOCATE(musurf(ndim1,ndim2))
IF(ndim1 /= 0 .AND. ndim2 /= 0)musurf(1:ndim1,1:ndim2) = work2(1:ndim1,1:ndim2)
DEALLOCATE(work2)

ndim1 = nexch_sec
ndim2 = ncomp+nexchange
i = size(muexc,1)
j = size(muexc,2)
ALLOCATE(work2(i,j))
work2 = muexc
DEALLOCATE(muexc)
ALLOCATE(muexc(ndim1,ndim2))
IF(ndim1 /= 0 .AND. ndim2 /= 0)muexc(1:ndim1,1:ndim2) = work2(1:ndim1,1:ndim2)
DEALLOCATE(work2)

ndim1 = ikin
ndim2 = ncomp
i = size(mukin,1)
j = size(mukin,2)
ALLOCATE(work2(i,j))
work2 = mukin
DEALLOCATE(mukin)
ALLOCATE(mukin(ndim1,ndim2))
IF(ndim1 /= 0 .AND. ndim2 /= 0)mukin(1:ndim1,1:ndim2) = work2(1:ndim1,1:ndim2)
DEALLOCATE(work2)

ndim1 = nrct
ndim2 = 6
i = size(volb,1)
j = size(volb,2)
ALLOCATE(work2(i,j))
work2 = volb
DEALLOCATE(volb)
ALLOCATE(volb(ndim1,6))
IF(ndim1 /= 0 .AND. ndim2 /= 0)volb(1:ndim1,1:6) = work2(1:ndim1,1:6)
DEALLOCATE(work2)

nbig = ncomp+nspec+ngas+nrct+nsurf_sec
!!  Since the equilibrium constant coefficients cycle over parallel reactions for minerals, count the total number needed
ncount = 0
DO k = 1,nrct
  do np = 1,nreactmin(k)
    ncount = ncount + 1
  END DO
END DO

!!ndim1 = nbig
ndim1 = ncomp+nspec+ngas+nsurf_sec+ncount
ndim2 = 5
i = size(as1,1)
j = size(as1,2)
ALLOCATE(work2(i,j))
work2 = as1
DEALLOCATE(as1)
ALLOCATE(as1(ndim1,5))
IF(ndim1 /= 0 .AND. ndim2 /= 0)as1(1:ndim1,1:5) = work2(1:ndim1,1:5)
DEALLOCATE(work2)

!!ndim1 = nbig
ndim1 = ncomp+nspec+ngas+nsurf_sec+ncount
ndim2 = 5
i = size(as2,1)
j = size(as2,2)
ALLOCATE(work2(i,j))
work2 = as2
DEALLOCATE(as2)
ALLOCATE(as2(ndim1,5))
IF(ndim1 /= 0 .AND. ndim2 /= 0)as2(1:ndim1,1:5) = work2(1:ndim1,1:5)
DEALLOCATE(work2)

ndim1 = nreactmax
ndim2 = nrct
i = size(ea,1)
j = size(ea,2)
ALLOCATE(work2(i,j))
work2 = ea
DEALLOCATE(ea)
ALLOCATE(ea(ndim1,ndim2))
IF(ndim1 /= 0 .AND. ndim2 /= 0)ea(1:ndim1,1:ndim2) = work2(1:ndim1,1:ndim2)
DEALLOCATE(work2)

ndim1 = nreactmax
ndim2 = nrct
i = size(sat1,1)
j = size(sat1,2)
ALLOCATE(work2(i,j))
work2 = sat1
DEALLOCATE(sat1)
ALLOCATE(sat1(ndim1,ndim2))
IF(ndim1 /= 0 .AND. ndim2 /= 0)sat1(1:ndim1,1:ndim2) = work2(1:ndim1,1:ndim2)
DEALLOCATE(work2)

ndim1 = nreactmax
ndim2 = nrct
i = size(sat2,1)
j = size(sat2,2)
ALLOCATE(work2(i,j))
work2 = sat2
DEALLOCATE(sat2)
ALLOCATE(sat2(ndim1,ndim2))
IF(ndim1 /= 0 .AND. ndim2 /= 0)sat2(1:ndim1,1:ndim2) = work2(1:ndim1,1:ndim2)
DEALLOCATE(work2)

ndim1 = nreactmax
ndim2 = nrct
i = size(rate0,1)
j = size(rate0,2)
ALLOCATE(work2(i,j))
work2 = rate0
DEALLOCATE(rate0)
ALLOCATE(rate0(ndim1,ndim2))
IF(ndim1 /= 0 .AND. ndim2 /= 0)rate0(1:ndim1,1:ndim2) = work2(1:ndim1,1:ndim2)
DEALLOCATE(work2)

ndim1 = nreactmax
ndim2 = nrct
i = size(ssa,1)
j = size(ssa,2)
ALLOCATE(work2(i,j))
work2 = ssa
DEALLOCATE(ssa)
ALLOCATE(ssa(ndim1,ndim2))
IF(ndim1 /= 0 .AND. ndim2 /= 0)ssa(1:ndim1,1:ndim2) = work2(1:ndim1,1:ndim2)
DEALLOCATE(work2)

ndim1 = nreactmax
ndim2 = nrct
i = size(AffinityDepend1,1)
j = size(AffinityDepend1,2)
ALLOCATE(work2(i,j))
work2 = AffinityDepend1
DEALLOCATE(AffinityDepend1)
ALLOCATE(AffinityDepend1(ndim1,ndim2))
IF(ndim1 /= 0 .AND. ndim2 /= 0)AffinityDepend1(1:ndim1,1:ndim2) = work2(1:ndim1,1:ndim2)
DEALLOCATE(work2)

ndim1 = nreactmax
ndim2 = nrct
i = size(AffinityDepend2,1)
j = size(AffinityDepend2,2)
ALLOCATE(work2(i,j))
work2 = AffinityDepend2
DEALLOCATE(AffinityDepend2)
ALLOCATE(AffinityDepend2(ndim1,ndim2))
IF(ndim1 /= 0 .AND. ndim2 /= 0)AffinityDepend2(1:ndim1,1:ndim2) = work2(1:ndim1,1:ndim2)
DEALLOCATE(work2)

ndim1 = nreactmax
ndim2 = nrct
i = size(AffinityDepend3,1)
j = size(AffinityDepend3,2)
ALLOCATE(work2(i,j))
work2 = AffinityDepend3
DEALLOCATE(AffinityDepend3)
ALLOCATE(AffinityDepend3(ndim1,ndim2))
IF(ndim1 /= 0 .AND. ndim2 /= 0)AffinityDepend3(1:ndim1,1:ndim2) = work2(1:ndim1,1:ndim2)
DEALLOCATE(work2)

! biomass
ndim1 = nreactmax
ndim2 = nrct
i = size(BQ_min,1)
j = size(BQ_min,2)
ALLOCATE(work2(i,j))
work2 = BQ_min
DEALLOCATE(BQ_min)
ALLOCATE(BQ_min(ndim1,ndim2))
IF(ndim1 /= 0 .AND. ndim2 /= 0) BQ_min(1:ndim1,1:ndim2) = work2(1:ndim1,1:ndim2)
DEALLOCATE(work2)

ndim1 = ikin
!!ndim2 = nrct
i = size(BQ_kin,1)
!!j = size(BQ_kin,2)
ALLOCATE(work1(i))
work1 = BQ_kin
DEALLOCATE(BQ_kin)
ALLOCATE(BQ_kin(ndim1))
IF(ndim1 /= 0) BQ_kin(1:ndim1) = work1(1:ndim1)
DEALLOCATE(work1)
! biomass end

ndim1 = nreactkinmax
ndim2 = ikin
i = size(ratek,1)
j = size(ratek,2)
ALLOCATE(work2(i,j))
work2 = ratek
DEALLOCATE(ratek)
ALLOCATE(ratek(ndim1,ndim2))
IF(ndim1 /= 0 .AND. ndim2 /= 0)ratek(1:ndim1,1:ndim2) = work2(1:ndim1,1:ndim2)
DEALLOCATE(work2)

!  *********  Character arrays  **************

ndim1 = nreactmax
ndim2 = nrct
i = size(rlabel,1)
j = size(rlabel,2)
ALLOCATE(workchar2(i,j))
workchar2 = rlabel
DEALLOCATE(rlabel)
ALLOCATE(rlabel(ndim1,ndim2))
IF(ndim1 /= 0.AND.ndim2 /= 0)rlabel(1:ndim1,1:ndim2) = workchar2(1:ndim1,1:ndim2)
DEALLOCATE(workchar2)

ndim1 = nreactmax
ndim2 = nrct
i = size(crossaff,1)
j = size(crossaff,2)
ALLOCATE(workchar2(i,j))
workchar2 = crossaff
DEALLOCATE(crossaff)
ALLOCATE(crossaff(ndim1,ndim2))
IF(ndim1 /= 0.AND.ndim2 /= 0)crossaff(1:ndim1,1:ndim2) = workchar2(1:ndim1,1:ndim2)
DEALLOCATE(workchar2)



!  ************  End of 2D array reallocation  ***********************

!  **************  ReALLOCATE 3D arrays  *****************************

!  *******  Integer arrays  ****************

ntot = ncomp + nspec

ndim1 = ncomp
ndim2 = nreactkinmax
ndim3 = ikin
i = size(itot,1)
j = size(itot,2)
k = size(itot,3)
ALLOCATE(workint3(i,j,k))
workint3 = itot
DEALLOCATE(itot)
ALLOCATE(itot(ndim1,ndim2,ndim3))
IF(ndim1 /= 0.AND.ndim2 /= 0.AND.ndim3 /= 0)  &
    itot(1:ndim1,1:ndim2,1:ndim3) = workint3(1:ndim1,1:ndim2,1:ndim3)
DEALLOCATE(workint3)

ndim1 = ntot
ndim2 = nreactmax
ndim3 = nrct
i = size(idepend,1)
j = size(idepend,2)
k = size(idepend,3)
ALLOCATE(workint3(i,j,k))
workint3 = idepend
DEALLOCATE(idepend)
ALLOCATE(idepend(ndim1,ndim2,ndim3))
IF(ndim1 /= 0.AND.ndim2 /= 0.AND.ndim3 /= 0)  &
    idepend(1:ndim1,1:ndim2,1:ndim3) = workint3(1:ndim1,1:ndim2,1:ndim3)
DEALLOCATE(workint3)

ndim1 = nmonodmax
ndim2 = nreactmax
ndim3 = nrct
i = size(imonod,1)
j = size(imonod,2)
k = size(imonod,3)
ALLOCATE(workint3(i,j,k))
workint3 = imonod
DEALLOCATE(imonod)
ALLOCATE(imonod(ndim1,ndim2,ndim3))
IF(ndim1 /= 0.AND.ndim2 /= 0.AND.ndim3 /= 0)  &
    imonod(1:ndim1,1:ndim2,1:ndim3) = workint3(1:ndim1,1:ndim2,1:ndim3)
DEALLOCATE(workint3)

ndim1 = nmonodmax
ndim2 = nreactmax
ndim3 = nrct
i = size(kmonod,1)
j = size(kmonod,2)
k = size(kmonod,3)
ALLOCATE(workint3(i,j,k))
workint3 = kmonod
DEALLOCATE(kmonod)
ALLOCATE(kmonod(ndim1,ndim2,ndim3))
IF(ndim1 /= 0.AND.ndim2 /= 0.AND.ndim3 /= 0)  &
    kmonod(1:ndim1,1:ndim2,1:ndim3) = workint3(1:ndim1,1:ndim2,1:ndim3)
DEALLOCATE(workint3)

ndim1 = ncomp
ndim2 = nreactmax
ndim3 = nrct
i = size(itot_min,1)
j = size(itot_min,2)
k = size(itot_min,3)
ALLOCATE(workint3(i,j,k))
workint3 = itot_min
DEALLOCATE(itot_min)
ALLOCATE(itot_min(ndim1,ndim2,ndim3))
IF(ndim1 /= 0.AND.ndim2 /= 0.AND.ndim3 /= 0)  &
    itot_min(1:ndim1,1:ndim2,1:ndim3) = workint3(1:ndim1,1:ndim2,1:ndim3)
DEALLOCATE(workint3)

ndim1 = ninhibitmax
ndim2 = nreactmax
ndim3 = nrct
i = size(itot_inhibit,1)
j = size(itot_inhibit,2)
k = size(itot_inhibit,3)
ALLOCATE(workint3(i,j,k))
workint3 = itot_inhibit
DEALLOCATE(itot_inhibit)
ALLOCATE(itot_inhibit(ndim1,ndim2,ndim3))
IF(ndim1 /= 0.AND.ndim2 /= 0.AND.ndim3 /= 0)  &
    itot_inhibit(1:ndim1,1:ndim2,1:ndim3) = workint3(1:ndim1,1:ndim2,1:ndim3)
DEALLOCATE(workint3)

ndim1 = ninhibitmax
ndim2 = nreactmax
ndim3 = nrct
i = size(inhibit,1)
j = size(inhibit,2)
k = size(inhibit,3)
ALLOCATE(workint3(i,j,k))
workint3 = inhibit
DEALLOCATE(inhibit)
ALLOCATE(inhibit(ndim1,ndim2,ndim3))
IF(ndim1 /= 0.AND.ndim2 /= 0.AND.ndim3 /= 0)  &
    inhibit(1:ndim1,1:ndim2,1:ndim3) = workint3(1:ndim1,1:ndim2,1:ndim3)
DEALLOCATE(workint3)


!  *******  Real arrays  ****************

ndim1 = nreactmax
ndim2 = nrct
ndim3 = ncomp
i = size(mumin,1)
j = size(mumin,2)
k = size(mumin,3)
ALLOCATE(work3(i,j,k))
work3 = mumin
DEALLOCATE(mumin)
ALLOCATE(mumin(ndim1,ndim2,ndim3))
IF(ndim1 /= 0.AND.ndim2 /= 0.AND.ndim3 /= 0)  &
    mumin(1:ndim1,1:ndim2,1:ndim3) = work3(1:ndim1,1:ndim2,1:ndim3)
DEALLOCATE(work3)

ndim1 = ntot
ndim2 = nreactmax
ndim3 = nrct
i = size(depend,1)
j = size(depend,2)
k = size(depend,3)
ALLOCATE(work3(i,j,k))
work3 = depend
DEALLOCATE(depend)
ALLOCATE(depend(ndim1,ndim2,ndim3))
IF(ndim1 /= 0.AND.ndim2 /= 0.AND.ndim3 /= 0)  &
    depend(1:ndim1,1:ndim2,1:ndim3) = work3(1:ndim1,1:ndim2,1:ndim3)
DEALLOCATE(work3)

ndim1 = nmonodmax
ndim2 = nreactmax
ndim3 = nrct
i = size(halfsat,1)
j = size(halfsat,2)
k = size(halfsat,3)
ALLOCATE(work3(i,j,k))
work3 = halfsat
DEALLOCATE(halfsat)
ALLOCATE(halfsat(ndim1,ndim2,ndim3))
IF(ndim1 /= 0.AND.ndim2 /= 0.AND.ndim3 /= 0)  &
    halfsat(1:ndim1,1:ndim2,1:ndim3) = work3(1:ndim1,1:ndim2,1:ndim3)
DEALLOCATE(work3)

ndim1 = ninhibitmax
ndim2 = nreactmax
ndim3 = nrct
i = size(rinhibit,1)
j = size(rinhibit,2)
k = size(rinhibit,3)
ALLOCATE(work3(i,j,k))
work3 = rinhibit
DEALLOCATE(rinhibit)
ALLOCATE(rinhibit(ndim1,ndim2,ndim3))
IF(ndim1 /= 0.AND.ndim2 /= 0.AND.ndim3 /= 0)  &
    rinhibit(1:ndim1,1:ndim2,1:ndim3) = work3(1:ndim1,1:ndim2,1:ndim3)
DEALLOCATE(work3)

ndim1 = ncomp
ndim2 = nreactkinmax
ndim3 = ikin
i = size(dependk,1)
j = size(dependk,2)
k = size(dependk,3)
ALLOCATE(work3(i,j,k))
work3 = dependk
DEALLOCATE(dependk)
ALLOCATE(dependk(ndim1,ndim2,ndim3))
IF(ndim1 /= 0.AND.ndim2 /= 0.AND.ndim3 /= 0)  &
    dependk(1:ndim1,1:ndim2,1:ndim3) = work3(1:ndim1,1:ndim2,1:ndim3)
DEALLOCATE(work3)

!  ************  End of 3D array reallocation  ***********************


END SUBROUTINE reallocate