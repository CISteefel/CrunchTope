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
    
SUBROUTINE FirstAllocation()

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
USE isotope

IMPLICIT NONE

!! External variables and arrays

ALLOCATE(keqmin_tmp(mreact,nm))
ALLOCATE(keqaq_tmp(mcmplx))
ALLOCATE(keqgas_tmp(ng))
ALLOCATE(keqsurf_tmp(msurf_sec))
ALLOCATE(sptmp(nc+mcmplx))
ALLOCATE(sptmp10(nc+mcmplx))
ALLOCATE(gamtmp(nc+mcmplx))
ALLOCATE(stmp(nc))
ALLOCATE(nbasin(nc))
ALLOCATE(nbkin(mcmplx))
ALLOCATE(dxxt(mzone))
ALLOCATE(dyyt(mzone))
ALLOCATE(dzzt(mzone))
ALLOCATE(nvx(mzone))
ALLOCATE(nvy(mzone))
ALLOCATE(nvz(mzone))
ALLOCATE(sexch(nc+mexch))
ALLOCATE(spextmp(mexch_sec+mexch))
ALLOCATE(spextmp10(mexch_sec+mexch))
ALLOCATE(totextmp(mexch))
ALLOCATE(spgastmp(ng))
ALLOCATE(spgastmp10(ng))
ALLOCATE(sgastmp(nc))
ALLOCATE(spsurftmp(msurf+msurf_sec))
ALLOCATE(spsurftmp10(msurf+msurf_sec))
ALLOCATE(ssurftmp(msurf))
ALLOCATE(dpsi(meqn,meqn))
ALLOCATE(namdep_nyf(5,10,mrct))


! **********  End of temporary arrays *************

!  **********  1D permanent arrays  ****************

IF (ALLOCATED(eqgas)) THEN
  DEALLOCATE(eqgas)
  ALLOCATE(eqgas(ng))
ELSE
  ALLOCATE(eqgas(ng))
END IF
IF (ALLOCATED(eqhom)) THEN
  DEALLOCATE(eqhom)
  ALLOCATE(eqhom(mcmplx))
ELSE
  ALLOCATE(eqhom(mcmplx))
END IF
IF (ALLOCATED(eqsurf)) THEN
  DEALLOCATE(eqsurf)
  ALLOCATE(eqsurf(msurf_sec))
ELSE
  ALLOCATE(eqsurf(msurf_sec))
END IF
IF (ALLOCATED(volmol)) THEN
  DEALLOCATE(volmol)
  ALLOCATE(volmol(nm))
ELSE
  ALLOCATE(volmol(nm))
END IF
IF (ALLOCATED(wtmin)) THEN
  DEALLOCATE(wtmin)
  ALLOCATE(wtmin(nm))
ELSE
  ALLOCATE(wtmin(nm))
END IF
IF (ALLOCATED(wtcomp)) THEN
  DEALLOCATE(wtcomp)
  ALLOCATE(wtcomp(nc))
ELSE
  ALLOCATE(wtcomp(nc))
END IF
IF (ALLOCATED(ivolume)) THEN
  DEALLOCATE(ivolume)
  ALLOCATE(ivolume(mrct))
ELSE
  ALLOCATE(ivolume(mrct))
END IF
IF (ALLOCATED(alnk)) THEN
  DEALLOCATE(alnk)
  ALLOCATE(alnk(mreact*nm))
ELSE
  ALLOCATE(alnk(mreact*nm))
END IF
IF (ALLOCATED(lenmin)) THEN
  DEALLOCATE(lenmin)
  ALLOCATE(lenmin(mrct))
ELSE
  ALLOCATE(lenmin(mrct))
END IF
IF (ALLOCATED(ulab)) THEN
  DEALLOCATE(ulab)
  ALLOCATE(ulab(nc+mcmplx))
ELSE
  ALLOCATE(ulab(nc+mcmplx))
END IF
IF (ALLOCATED(umin)) THEN
  DEALLOCATE(umin)
  ALLOCATE(umin(nm))
ELSE
  ALLOCATE(umin(nm))
END IF
IF (ALLOCATED(namg)) THEN
  DEALLOCATE(namg)
  ALLOCATE(namg(ng))
ELSE
  ALLOCATE(namg(ng))
END IF
IF (ALLOCATED(namkin)) THEN
  DEALLOCATE(namkin)
  ALLOCATE(namkin(maqkin))
ELSE
  ALLOCATE(namkin(maqkin))
END IF
IF (ALLOCATED(namexc)) THEN
  DEALLOCATE(namexc)
  ALLOCATE(namexc(mexch))
ELSE
  ALLOCATE(namexc(mexch))
END IF
IF (ALLOCATED(nam_exchsec)) THEN
  DEALLOCATE(nam_exchsec)
  ALLOCATE(nam_exchsec(mexch_sec))
ELSE
  ALLOCATE(nam_exchsec(mexch_sec))
END IF
IF (ALLOCATED(namsurf)) THEN
  DEALLOCATE(namsurf)
  ALLOCATE(namsurf(msurf))
ELSE
  ALLOCATE(namsurf(msurf))
END IF
IF (ALLOCATED(namsurf_sec)) THEN
  DEALLOCATE(namsurf_sec)
  ALLOCATE(namsurf_sec(msurf_sec))
ELSE
  ALLOCATE(namsurf_sec(msurf_sec))
END IF
IF (ALLOCATED(acmp)) THEN
  DEALLOCATE(acmp)
  ALLOCATE(acmp(mcmplx+nc))
ELSE
  ALLOCATE(acmp(mcmplx+nc))
END IF
IF (ALLOCATED(chg)) THEN
  DEALLOCATE(chg)
  ALLOCATE(chg(mcmplx+nc))
ELSE
  ALLOCATE(chg(mcmplx+nc))
END IF

IF (ALLOCATED(bdotParameter)) THEN
  DEALLOCATE(bdotParameter)
  ALLOCATE(bdotParameter(mcmplx+nc))
ELSE
  ALLOCATE(bdotParameter(mcmplx+nc))
END IF
bdotParameter = 0.0d0

IF (ALLOCATED(adh)) THEN
  DEALLOCATE(adh)
  ALLOCATE(adh(8))
ELSE
  ALLOCATE(adh(8))
END IF
IF (ALLOCATED(bdh)) THEN
  DEALLOCATE(bdh)
  ALLOCATE(bdh(8))
ELSE
  ALLOCATE(bdh(8))
END IF
IF (ALLOCATED(bdot)) THEN
  DEALLOCATE(bdot)
  ALLOCATE(bdot(8))
ELSE
  ALLOCATE(bdot(8))
END IF
IF (ALLOCATED(alogkp)) THEN
  DEALLOCATE(alogkp)
  ALLOCATE(alogkp(ndim))
ELSE
  ALLOCATE(alogkp(ndim))
END IF
IF (ALLOCATED(wtaq)) THEN
  DEALLOCATE(wtaq)
  ALLOCATE(wtaq(nc+mcmplx))
ELSE
  ALLOCATE(wtaq(nc+mcmplx))
END IF
IF (ALLOCATED(satkin)) THEN
  DEALLOCATE(satkin)
  ALLOCATE(satkin(maqkin))
ELSE
  ALLOCATE(satkin(maqkin))
END IF
IF (ALLOCATED(ksurf)) THEN
  DEALLOCATE(ksurf)
  ALLOCATE(ksurf(msurf))
ELSE
  ALLOCATE(ksurf(msurf))
END IF
IF (ALLOCATED(kexch)) THEN
  DEALLOCATE(kexch)
  ALLOCATE(kexch(mexch))
ELSE
  ALLOCATE(kexch(mexch))
END IF
IF (ALLOCATED(nreactkin)) THEN
  DEALLOCATE(nreactkin)
  ALLOCATE(nreactkin(maqkin))
ELSE
  ALLOCATE(nreactkin(maqkin))
END IF
IF (ALLOCATED(iplot)) THEN
  DEALLOCATE(iplot)
  ALLOCATE(iplot(mcomp+mspec))
ELSE
  ALLOCATE(iplot(mcomp+mspec))
END IF

IF (ALLOCATED(jxseries)) THEN
  DEALLOCATE(jxseries)
  ALLOCATE(jxseries(500))
ELSE
  ALLOCATE(jxseries(500))
END IF
IF (ALLOCATED(jyseries)) THEN
  DEALLOCATE(jyseries)
  ALLOCATE(jyseries(500))
ELSE
  ALLOCATE(jyseries(500))
END IF
IF (ALLOCATED(jzseries)) THEN
  DEALLOCATE(jzseries)
  ALLOCATE(jzseries(500))
ELSE
  ALLOCATE(jzseries(500))
END IF

!!!  Add jxminseries arrays
IF (ALLOCATED(jxminseries)) THEN
  DEALLOCATE(jxminseries)
  ALLOCATE(jxminseries(500))
ELSE
  ALLOCATE(jxminseries(500))
END IF
IF (ALLOCATED(jyminseries)) THEN
  DEALLOCATE(jyminseries)
  ALLOCATE(jyminseries(500))
ELSE
  ALLOCATE(jyminseries(500))
END IF
IF (ALLOCATED(jzminseries)) THEN
  DEALLOCATE(jzminseries)
  ALLOCATE(jzminseries(500))
ELSE
  ALLOCATE(jzminseries(500))
END IF


IF (ALLOCATED(jxisotopeseries)) THEN
  DEALLOCATE(jxisotopeseries)
  ALLOCATE(jxisotopeseries(500))
ELSE
  ALLOCATE(jxisotopeseries(500))
END IF
IF (ALLOCATED(jyisotopeseries)) THEN
  DEALLOCATE(jyisotopeseries)
  ALLOCATE(jyisotopeseries(500))
ELSE
  ALLOCATE(jyisotopeseries(500))
END IF
IF (ALLOCATED(jzisotopeseries)) THEN
  DEALLOCATE(jzisotopeseries)
  ALLOCATE(jzisotopeseries(500))
ELSE
  ALLOCATE(jzisotopeseries(500))
END IF

!!!IF (ALLOCATED(jxAqueousFluxSeries)) THEN
!!!  DEALLOCATE(jxAqueousFluxSeries)
!!!  ALLOCATE(jxAqueousFluxSeries(500))
!!!ELSE
!!!  ALLOCATE(jxAqueousFluxSeries(500))
!!!END IF
!!!IF (ALLOCATED(jyAqueousFluxSeries)) THEN
!!!  DEALLOCATE(jyAqueousFluxSeries)
!!!  ALLOCATE(jyAqueousFluxSeries(500))
!!!ELSE
!!!  ALLOCATE(jyAqueousFluxSeries(500))
!!!END IF
!!!IF (ALLOCATED(jzAqueousFluxSeries)) THEN
!!!  DEALLOCATE(jzAqueousFluxSeries)
!!!  ALLOCATE(jzAqueousFluxSeries(500))
!!!ELSE
!!!  ALLOCATE(jzAqueousFluxSeries(500))
!!!END IF

IF (ALLOCATED(iexchange)) THEN
  DEALLOCATE(iexchange)
  ALLOCATE(iexchange(mexch))
ELSE
  ALLOCATE(iexchange(mexch))
END IF
IF (ALLOCATED(nreactmin)) THEN
  DEALLOCATE(nreactmin)
  ALLOCATE(nreactmin(mrct))
ELSE
  ALLOCATE(nreactmin(mrct))
END IF
IF (ALLOCATED(keqexc)) THEN
  DEALLOCATE(keqexc)
  ALLOCATE(keqexc(mexch_sec))
ELSE
  ALLOCATE(keqexc(mexch_sec))
END IF
IF (ALLOCATED(bfit)) THEN
  DEALLOCATE(bfit)
  ALLOCATE(bfit(mexch_sec))
ELSE
  ALLOCATE(bfit(mexch_sec))
END IF
IF (ALLOCATED(ixlink)) THEN
  DEALLOCATE(ixlink)
  ALLOCATE(ixlink(mexch_sec))
ELSE
  ALLOCATE(ixlink(mexch_sec))
END IF
IF (ALLOCATED(nclink)) THEN
  DEALLOCATE(nclink)
  ALLOCATE(nclink(mexch_sec))
ELSE
  ALLOCATE(nclink(mexch_sec))
END IF
IF (ALLOCATED(iedl)) THEN
  DEALLOCATE(iedl)
  ALLOCATE(iedl(msurf))
ELSE
  ALLOCATE(iedl(msurf))
END IF
IF (ALLOCATED(LocalEquilibrium)) THEN
  DEALLOCATE(LocalEquilibrium)
  ALLOCATE(LocalEquilibrium(mrct))
ELSE
  ALLOCATE(LocalEquilibrium(mrct))
END IF

!  **********  2D permanent arrays  ****************

IF (ALLOCATED(volb)) THEN
  DEALLOCATE(volb)
  ALLOCATE(volb(mrct,6))
ELSE
  ALLOCATE(volb(mrct,6))
END IF
IF (ALLOCATED(as1)) THEN
  DEALLOCATE(as1)
  ALLOCATE(as1(ndim,5))
ELSE
  ALLOCATE(as1(ndim,5))
END IF
IF (ALLOCATED(as2)) THEN
  DEALLOCATE(as2)
  ALLOCATE(as2(ndim,5))
ELSE
  ALLOCATE(as2(ndim,5))
END IF
IF (ALLOCATED(sat1)) THEN
  DEALLOCATE(sat1)
  ALLOCATE(sat1(mreact,mrct))
ELSE
  ALLOCATE(sat1(mreact,mrct))
END IF
IF (ALLOCATED(sat2)) THEN
  DEALLOCATE(sat2)
  ALLOCATE(sat2(mreact,mrct))
ELSE
  ALLOCATE(sat2(mreact,mrct))
END IF
IF (ALLOCATED(ea)) THEN
  DEALLOCATE(ea)
  ALLOCATE(ea(mreact,mrct))
ELSE
  ALLOCATE(ea(mreact,mrct))
END IF
IF (ALLOCATED(rate0)) THEN
  DEALLOCATE(rate0)
  ALLOCATE(rate0(mreact,mrct))
ELSE
  ALLOCATE(rate0(mreact,mrct))
END IF
IF (ALLOCATED(ssa)) THEN
  DEALLOCATE(ssa)
  ALLOCATE(ssa(mreact,mrct))
ELSE
  ALLOCATE(ssa(mreact,mrct))
END IF
IF (ALLOCATED(imintype)) THEN
  DEALLOCATE(imintype)
  ALLOCATE(imintype(mreact,mrct))
ELSE
  ALLOCATE(imintype(mreact,mrct))
END IF
IF (ALLOCATED(ndepend)) THEN
  DEALLOCATE(ndepend)
  ALLOCATE(ndepend(mreact,mrct))
ELSE
  ALLOCATE(ndepend(mreact,mrct))
END IF

!!  ********** Hyperbolic Inhibition ******************
IF (ALLOCATED(HyperbolicInhibition)) THEN
  DEALLOCATE(HyperbolicInhibition)
  ALLOCATE(HyperbolicInhibition(mreact,mrct))
ELSE
  ALLOCATE(HyperbolicInhibition(mreact,mrct))
END IF
HyperbolicInhibition = .FALSE.
IF (ALLOCATED(HyperbolicInhibitionPointer)) THEN
  DEALLOCATE(HyperbolicInhibitionPointer)
  ALLOCATE(HyperbolicInhibitionPointer(mreact,mrct))
ELSE
  ALLOCATE(HyperbolicInhibitionPointer(mreact,mrct))
END IF
IF (ALLOCATED(Kformation)) THEN
  DEALLOCATE(Kformation)
  ALLOCATE(Kformation(mreact,mrct))
ELSE
  ALLOCATE(Kformation(mreact,mrct))
END IF
IF (ALLOCATED(HyperbolicInhibitionDepend)) THEN
  DEALLOCATE(HyperbolicInhibitionDepend)
  ALLOCATE(HyperbolicInhibitionDepend(mreact,mrct))
ELSE
  ALLOCATE(HyperbolicInhibitionDepend(mreact,mrct))
END IF
IF (ALLOCATED(HyperbolicInhibitionName)) THEN
  DEALLOCATE(HyperbolicInhibitionName)
  ALLOCATE(HyperbolicInhibitionName(mreact,mrct))
ELSE
  ALLOCATE(HyperbolicInhibitionName(mreact,mrct))
END IF
!!  ********** Hyperbolic Inhibition ******************

IF (ALLOCATED(nmonod)) THEN
  DEALLOCATE(nmonod)
  ALLOCATE(nmonod(mreact,mrct))
ELSE
  ALLOCATE(nmonod(mreact,mrct))
END IF
IF (ALLOCATED(ninhibit)) THEN
  DEALLOCATE(ninhibit)
  ALLOCATE(ninhibit(mreact,mrct))
ELSE
  ALLOCATE(ninhibit(mreact,mrct))
END IF
IF (ALLOCATED(kcrossaff)) THEN
  DEALLOCATE(kcrossaff)
  ALLOCATE(kcrossaff(mreact,mrct))
ELSE
  ALLOCATE(kcrossaff(mreact,mrct))
END IF
IF (ALLOCATED(rlabel)) THEN
  DEALLOCATE(rlabel)
  ALLOCATE(rlabel(mreact,mrct))
ELSE
  ALLOCATE(rlabel(mreact,mrct))
END IF
IF (ALLOCATED(crossaff)) THEN
  DEALLOCATE(crossaff)
  ALLOCATE(crossaff(mreact,mrct))
ELSE
  ALLOCATE(crossaff(mreact,mrct))
END IF
crossaff = 'none'
IF (ALLOCATED(muaq)) THEN
  DEALLOCATE(muaq)
  ALLOCATE(muaq(mcmplx,nc))
ELSE
  ALLOCATE(muaq(mcmplx,nc))
END IF
IF (ALLOCATED(mugas)) THEN
  DEALLOCATE(mugas)
  ALLOCATE(mugas(ng,nc))
ELSE
  ALLOCATE(mugas(ng,nc))
END IF
IF (ALLOCATED(musurf)) THEN
  DEALLOCATE(musurf)
  ALLOCATE(musurf(msurf_sec,nc+msurf))
ELSE
  ALLOCATE(musurf(msurf_sec,nc+msurf))
END IF
IF (ALLOCATED(muexc)) THEN
  DEALLOCATE(muexc)
  ALLOCATE(muexc(mexch_sec,nc+msurf))
ELSE
  ALLOCATE(muexc(mexch_sec,nc+msurf))
END IF
IF (ALLOCATED(AffinityDepend1)) THEN
  DEALLOCATE(AffinityDepend1)
  ALLOCATE(AffinityDepend1(mreact,mrct))
ELSE
  ALLOCATE(AffinityDepend1(mreact,mrct))
END IF
IF (ALLOCATED(AffinityDepend2)) THEN
  DEALLOCATE(AffinityDepend2)
  ALLOCATE(AffinityDepend2(mreact,mrct))
ELSE
  ALLOCATE(AffinityDepend2(mreact,mrct))
END IF
IF (ALLOCATED(AffinityDepend3)) THEN
  DEALLOCATE(AffinityDepend3)
  ALLOCATE(AffinityDepend3(mreact,mrct))
ELSE
  ALLOCATE(AffinityDepend3(mreact,mrct))
END IF

! biomass
IF (ALLOCATED(chi_min)) THEN
  DEALLOCATE(chi_min)
  ALLOCATE(chi_min(mreact,mrct))
ELSE
  ALLOCATE(chi_min(mreact,mrct))
END IF
IF (ALLOCATED(direction_min)) THEN
  DEALLOCATE(direction_min)
  ALLOCATE(direction_min(mreact,mrct))
ELSE
  ALLOCATE(direction_min(mreact,mrct))
END IF
IF (ALLOCATED(BQ_min)) THEN
  DEALLOCATE(BQ_min)
  ALLOCATE(BQ_min(mreact,mrct))
ELSE
  ALLOCATE(BQ_min(mreact,mrct))
END IF

! Aqueous
IF (ALLOCATED(chi_kin)) THEN
  DEALLOCATE(chi_kin)
  ALLOCATE(chi_kin(mrct))
ELSE
  ALLOCATE(chi_kin(mrct))
END IF
IF (ALLOCATED(direction_kin)) THEN
  DEALLOCATE(direction_kin)
  ALLOCATE(direction_kin(mrct))
ELSE
  ALLOCATE(direction_kin(mrct))
END IF
IF (ALLOCATED(BQ_kin)) THEN
  DEALLOCATE(BQ_kin)
  ALLOCATE(BQ_kin(mrct))
ELSE
  ALLOCATE(BQ_kin(mrct))
END IF

! biomass end


!  **********  3D permanent arrays  ****************

IF (ALLOCATED(idepend)) THEN
  DEALLOCATE(idepend)
  ALLOCATE(idepend(mtot+msurf+msurf_sec+mexch_sec,mreact,mrct))
ELSE
  ALLOCATE(idepend(mtot+msurf+msurf_sec+mexch_sec,mreact,mrct))
END IF
IF (ALLOCATED(depend)) THEN
  DEALLOCATE(depend)
  ALLOCATE(depend(mtot+msurf+msurf_sec+mexch_sec,mreact,mrct))
ELSE
  ALLOCATE(depend(mtot+msurf+msurf_sec+mexch_sec,mreact,mrct))
END IF
IF (ALLOCATED(imonod)) THEN
  DEALLOCATE(imonod)
  ALLOCATE(imonod(mtot,mreact,mrct))
ELSE
  ALLOCATE(imonod(mtot,mreact,mrct))
END IF
IF (ALLOCATED(kmonod)) THEN
  DEALLOCATE(kmonod)
  ALLOCATE(kmonod(mtot,mreact,mrct))
ELSE
  ALLOCATE(kmonod(mtot,mreact,mrct))
END IF
IF (ALLOCATED(itot_min)) THEN
  DEALLOCATE(itot_min)
  ALLOCATE(itot_min(mcomp,mreact,mrct))
ELSE
  ALLOCATE(itot_min(mcomp,mreact,mrct))
END IF
IF (ALLOCATED(itot_monod)) THEN
  DEALLOCATE(itot_monod)
  ALLOCATE(itot_monod(mcomp,mreact,mrct))
ELSE
  ALLOCATE(itot_monod(mcomp,mreact,mrct))
END IF
IF (ALLOCATED(itot_inhibit)) THEN
  DEALLOCATE(itot_inhibit)
  ALLOCATE(itot_inhibit(mcomp,mreact,mrct))
ELSE
  ALLOCATE(itot_inhibit(mcomp,mreact,mrct))
END IF
IF (ALLOCATED(inhibit)) THEN
  DEALLOCATE(inhibit)
  ALLOCATE(inhibit(mtot,mreact,mrct))
ELSE
  ALLOCATE(inhibit(mtot,mreact,mrct))
END IF
IF (ALLOCATED(halfsat)) THEN
  DEALLOCATE(halfsat)
  ALLOCATE(halfsat(mtot,mreact,mrct))
ELSE
  ALLOCATE(halfsat(mtot,mreact,mrct))
END IF
IF (ALLOCATED(rinhibit)) THEN
  DEALLOCATE(rinhibit)
  ALLOCATE(rinhibit(mtot,mreact,mrct))
ELSE
  ALLOCATE(rinhibit(mtot,mreact,mrct))
END IF
IF (ALLOCATED(mumin)) THEN
  DEALLOCATE(mumin)
  ALLOCATE(mumin(mreact,nm,nc))
ELSE
  ALLOCATE(mumin(mreact,nm,nc))
END IF

IF (ALLOCATED(itot)) THEN
  DEALLOCATE(itot)
  ALLOCATE(itot(mcomp,mreact,maqkin))
ELSE
  ALLOCATE(itot(mcomp,mreact,maqkin))
END IF
IF (ALLOCATED(ratek)) THEN
  DEALLOCATE(ratek)
  ALLOCATE(ratek(mreact,maqkin))
ELSE
  ALLOCATE(ratek(mreact,maqkin))
END IF
IF (ALLOCATED(dependk)) THEN
  DEALLOCATE(dependk)
  ALLOCATE(dependk(mcomp,mreact,maqkin))
ELSE
  ALLOCATE(dependk(mcomp,mreact,maqkin))
END IF

END SUBROUTINE FirstAllocation