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
    

SUBROUTINE speciation(ncomp,nrct,nkin,nspec,ngas,nexchange,nexch_sec,nsurf,nsurf_sec,npot,  &
    ndecay,ikin,nx,ny,nz,realtime,nn,nint,ikmast,ikph,delt)
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
USE strings

IMPLICIT NONE

!  *********************  INTERFACE BLOCKS  *****************************
INTERFACE
  SUBROUTINE GasPartialPressure(ncomp,ngas,gastmp10,jx,jy,jz)
    USE crunchtype
    INTEGER(I4B), INTENT(IN)                                   :: ncomp
    INTEGER(I4B), INTENT(IN)                                   :: ngas
    REAL(DP), DIMENSION(:)                                     :: gastmp10
    INTEGER(I4B), INTENT(IN)                                   :: jx
    INTEGER(I4B), INTENT(IN)                                   :: jy
    INTEGER(I4B), INTENT(IN)                                   :: jz
  END SUBROUTINE GasPartialPressure
END INTERFACE
!  **********************************************************************

!  External variables and arrays

REAL(DP), INTENT(IN)                               :: realtime
REAL(DP), INTENT(IN)                               :: delt

INTEGER(I4B), INTENT(IN)                           :: ncomp
INTEGER(I4B), INTENT(IN)                           :: nrct
INTEGER(I4B), INTENT(IN)                           :: nspec
INTEGER(I4B), INTENT(IN)                           :: ndecay
INTEGER(I4B), INTENT(IN)                           :: nsurf
INTEGER(I4B), INTENT(IN)                           :: nsurf_sec
INTEGER(I4B), INTENT(IN)                           :: npot
INTEGER(I4B), INTENT(IN)                           :: ikin
INTEGER(I4B), INTENT(IN)                           :: nkin
INTEGER(I4B), INTENT(IN)                           :: ngas
INTEGER(I4B), INTENT(IN)                           :: nexchange
INTEGER(I4B), INTENT(IN)                           :: nexch_sec
INTEGER(I4B), INTENT(IN)                           :: nx
INTEGER(I4B), INTENT(IN)                           :: ny
INTEGER(I4B), INTENT(IN)                           :: nz
INTEGER(I4B), INTENT(IN)                           :: nn
INTEGER(I4B), INTENT(IN)                           :: nint
INTEGER(I4B), INTENT(IN)                           :: ikmast
INTEGER(I4B), INTENT(IN)                           :: ikph

!  Internal variables and arrays

CHARACTER (LEN=13), DIMENSION(nrct)                :: uminprnt
CHARACTER (LEN=13), DIMENSION(ncomp+nspec)         :: ulabprnt
CHARACTER (LEN=mls)                                 :: fn
CHARACTER (LEN=mls)                                  :: suf
CHARACTER (LEN=mls)                                  :: suf1
CHARACTER (LEN=mls)                                 :: fnv
CHARACTER (LEN=1)                                  :: tab
CHARACTER (LEN=mls), DIMENSION(nsurf+nsurf_sec)    :: prtsurf
 
INTEGER(I4B), DIMENSION(nrct)                      :: len_min
INTEGER(I4B)                                       :: j
INTEGER(I4B)                                       :: jx
INTEGER(I4B)                                       :: jy
INTEGER(I4B)                                       :: jz
INTEGER(I4B)                                       :: ilength
INTEGER(I4B)                                       :: ik
INTEGER(I4B)                                       :: k
INTEGER(I4B)                                       :: ks
INTEGER(I4B)                                       :: ns
INTEGER(I4B)                                       :: i
INTEGER(I4B)                                       :: nex
INTEGER(I4B)                                       :: ir
INTEGER(I4B)                                       :: lsjx
INTEGER(I4B)                                       :: nlen

REAL(DP), DIMENSION(ncomp)                         :: totex_bas
REAL(DP), DIMENSION(nrct)                          :: dptprt
REAL(DP), DIMENSION(nrct)                          :: dsat
REAL(DP), DIMENSION(nrct)                          :: dvolpr
REAL(DP)                                           :: sum
REAL(DP)                                           :: porprt
REAL(DP)                                           :: phprt
REAL(DP)                                           :: porcalc
REAL(DP)                                           :: SurfacePerKgWat
REAL(DP)                                           :: EquivalentsPerKgWat
REAL(DP)                                           :: AqueousToBulk

REAL(DP)                                                   :: sumiap
REAL(DP)                                                   :: pHprint
REAL(DP)                                                   :: peprint
REAL(DP)                                                   :: Ehprint
REAL(DP)                                                   :: spprint
REAL(DP)                                                   :: totcharge
REAL(DP)                                                   :: silnTMP
REAL(DP)                                                   :: siprnt
REAL(DP)                                                   :: actprint
REAL(DP)                                                   :: actprint10
REAL(DP)                                                   :: spbase
REAL(DP)                                                   :: rone
REAL(DP)                                                   :: PrintTime
REAL(DP)                                                   :: alk
REAL(DP)                                                   :: tflux_top
REAL(DP)                                                   :: tflux_bot
REAL(DP)                                                   :: top_norm
REAL(DP)                                                   :: bot_norm
REAL(DP)                                                   :: aflux_net
REAL(DP)                                                   :: ad_net_bot
REAL(DP)                                                   :: SolidSolutionRatioTemp

CHARACTER (LEN=mls)                                        :: namtemp

INTEGER(I4B)                                               :: ix
INTEGER(I4B)                                               :: is
INTEGER(I4B)                                               :: npt

REAL(DP)                                                   :: MeanSaltConcentration
REAL(DP)                                                   :: MassFraction

REAL(DP), DIMENSION(ngas)                                  :: gastmp10
REAL(DP)                                                   :: denmol
REAL(DP)                                                   :: tk
REAL(DP)                                                   :: wtt
REAL(DP)                                                   :: pg

REAL(DP)                                                   :: lnActivity
CHARACTER (LEN=3)                                          :: ulabPrint

CHARACTER (LEN=17)                                         :: EDLoption

PrintTime = realtime*OutputTimeScale
rone = 1.0d0

suf='.out'
suf1 ='.out'
tab = CHAR(9)

DO k = 1,nrct
  uminprnt(k) = umin(k)
END DO
DO ik = 1,ncomp+nspec
  ulabprnt(ik) = ulab(ik)
END DO
DO ks = 1,nsurf
  prtsurf(ks) = namsurf(ks)
END DO
DO ns = 1,nsurf_sec
  prtsurf(ns+nsurf) = namsurf_sec(ns)
END DO


!  Write out speciation at each grid point

fn='speciation'
CALL stringlen(fn,ilength)
CALL newfile(fn,suf1,fnv,nint,ilength)
OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
WRITE(8,2283) PrintTime
WRITE(8,*)
DO jz = 1,nz
  DO jy = 1,ny
    DO jx = 1,nx

      tk = 273.15d0 + t(jx,jy,jz)
      denmol = 1.e05/(8.314*tk)   ! P/RT = n/V, with pressure converted from bars to Pascals

      IF (DensityModule /= 'temperature') THEN
!       Calculate the correction for the mass fraction of water:  kg_solution/kg_water
        MeanSaltConcentration = 0.001*( wtaq(MeanSalt(1))*s(MeanSalt(1),jx,jy,jz) +   &
            wtaq(MeanSalt(2))*s(MeanSalt(2),jx,jy,jz) ) 
        MassFraction = 1.0/(1.0 + MeanSaltConcentration)
      ELSE
        MassFraction = 1.0
      END IF
      
      AqueousToBulk = MassFraction*ro(jx,jy,jz)*satliq(jx,jy,jz)*por(jx,jy,jz)
      
      WRITE(8,*) '**************************************************************'
      WRITE(8,FMT=102) jx,jy,jz
      WRITE(8,*)

555 FORMAT(2X, 'Temperature (C)    = ',f10.3)
558 FORMAT(2X, 'Porosity           = ',f10.3)
556 FORMAT(2X, 'Liquid Saturation  = ',f10.3)
557 FORMAT(2X, 'Liquid Density     = ',f10.3)
563 FORMAT(2X, 'Activity of water  = ',f10.3)
564 FORMAT(2X, 'Pressure (bars)    = ',f10.3)
411 FORMAT(2X, 'Ionic Strength     = ',f10.3)
5022 FORMAT(2X,'Solution pH        = ',f10.3)
5023 FORMAT(2X,'Solution pe        = ',f10.3)
5024 FORMAT(2X,'Solution Eh        = ',f10.3)
205  FORMAT(2X,'Total Charge       = ',1pe10.3)

      WRITE(8,555) t(jx,jy,jz)
      WRITE(8,564) PressureCond(jinit(jx,jy,jz))
      WRITE(8,558) por(jx,jy,jz)
      WRITE(8,556) satliq(jx,jy,jz)
      WRITE(8,557) ro(jx,jy,jz)
      WRITE(8,563) gammawater(jx,jy,jz)
      WRITE(8,411) sion(jx,jy,jz)
      IF (ikph /= 0) THEN
        pHprint = -(sp(ikph,jx,jy,jz)+lngamma(ikph,jx,jy,jz))/clg
        WRITE(8,5022) pHprint
      END IF
      IF (O2Found) THEN     !  Calculate pe based on O2-H2O couple
        IF (npointH2gas == 0) THEN
          WRITE(8,5025)
          5025 FORMAT(2x,'Solution pe and Eh require H2(g) be included')
        ELSE
          peprint = (-0.5*keqgas(npointH2gas,jx,jy,jz) +             &
            0.25*(sp(npointO2aq,jx,jy,jz)+lngamma(npointO2aq,jx,jy,jz)) -     &  
            pHprint*clg + 0.25*keqgas(npointO2gas,jx,jy,jz) )/clg
          WRITE(8,5023) peprint
          Ehprint = peprint*clg*rgas*(t(jx,jy,jz)+273.15)/23.06
          WRITE(8,5024) Ehprint   
        END IF 
      END IF
      totcharge = 0.0
      DO i = 1,ncomp
        totcharge = totcharge + chg(i)*s(i,jx,jy,jz)
      END DO  
      WRITE(8,205) totcharge

      WRITE(8,*)
      WRITE(8,*) 'Total Aqueous Concentrations of Primary Species '
      WRITE(8,*) '--------------------------------------- '
      WRITE(8,*)
      WRITE(8,*) '  Species              Molality       Type '
      namtemp = 'Aqueous'
      DO i = 1,ncomp
        WRITE(8,511) ulab(i),s(i,jx,jy,jz),namtemp
      END DO

      namtemp = 'Exchange'
      IF (nexchange > 0) THEN
        SolidSolutionRatioTemp = 1000.d0*SolidDensity(jinit(jx,jy,jz))*(1.0d0-por(jx,jy,jz))
        WRITE(8,*)
        WRITE(8,*) 'Exchangers           Equiv/kgw        Equiv/g solid    Equiv/m^3 bulk'
        DO ix = 1,nexchange
          EquivalentsPerKgWat = exchangesites(ix,jx,jy,jz)/AqueousToBulk     !!! Conversion ExchangeSites/m_pm^3 to ExchangeSites/kgw
          WRITE(8,515) namexc(ix),EquivalentsPerKgWat,exchangesites(ix,jx,jy,jz)/SolidSolutionRatioTemp,exchangesites(ix,jx,jy,jz)
        END DO
      END IF

!!CIS  The following call does not work with erosion (need to call "totsurf.F90")
      IF (ierode /= 1) THEN
        call SurfaceComplex(ncomp,nsurf,nsurf_sec,nx,ny,nz)
        call totsurf_local(ncomp,nsurf,nsurf_sec,jx,jy,jz)
      END IF

      namtemp = 'Surface Complex'
      IF (nsurf > 0) THEN
        
        SolidSolutionRatioTemp = 1000.0*SolidDensity(jinit(jx,jy,jz))*(1.0d0-por(jx,jy,jz))
        
        WRITE(8,*)
        WRITE(8,*) 'Surface complex   Sites/kgw      Mol/g solid    Mol/m^3 bulk  Option'
        DO is = 1,nsurf
          IF (iedl(is) == 0) THEN
            EDLoption = 'Electrostatic'
          ELSE
            EDLoption = 'Non-Electrostatic'
          END IF
          IF (gimrt .and. ierode == 1) THEN
            SurfacePerKgWat = ssurf(is,jx,jy,jz)/AqueousToBulk         !!! Conversion SurfaceSites/m_pm^3 to SurfaceSites/kgw
            WRITE(8,515) namsurf(is),SurfacePerKgWat,ssurf(is,jy,jy,jz)/   & 
                  SolidSolutionRatioTemp,ssurf(is,jx,jy,jz),EDLoption
          ELSE
            SurfacePerKgWat = ssurf_local(is)/AqueousToBulk
            WRITE(8,515) namsurf(is),SurfacePerKgWat,ssurf_local(is)/      &
                    SolidSolutionRatioTemp,ssurf_local(is),EDLoption
          END IF
        END DO
      
      END IF
    
515   FORMAT(a14,4X,1PE11.4,4X,1PE11.4,4x,1PE11.4,4x,a17)
      
      namtemp = 'Electrostatic Potentials'
      IF (npot > 0) THEN
      
        WRITE(8,*)
        WRITE(8,*) 'Electrostatic Potentials'
        DO npt = 1,npot
          WRITE(8,518) npt,umin(kpot(npt)),LogPotential(npt,jx,jy,jz)
          WRITE(8,*) '    Primary surface complexes '
          DO is = 1,nsurf
            IF (nptPrimary(is) == npt) THEN
              WRITE(8,517) namsurf(is), zsurf(is)
            END IF
          END DO
          WRITE(8,*) '    Secondary surface complexes '
          DO ns = 1,nsurf_sec
            IF (nptlink(ns) == npt) THEN
              WRITE(8,517) namsurf_sec(ns), zsurf(nsurf+ns)
            END IF
          END DO
        END DO  
      
      END IF
    
!!!  Surface Complexation Cheat Sheet    
!!!    kPotential(k) --> Logical to EDL potential
!!!    ksurf(is) --> pointer for primary nsurf complex to mineral (initialized in read_surface.F90)
!!!    iedl(is) --> 0 for electrostatic, 1 for -no_edl
!!!    npot --> number of potentials
!!!    kpot(npt) --> pointer to mineral upon which the potential is developed
!!!    islink(ns) --> pointer from secondary surface complex (ns) to primary surface complex (is)
!!!    ksurf(islink(ns)) --> This would point from a secondary surface complex (ns) to a primary (islink(ns)) complex to a mineral
!!!    nptlink(ns) --> pointer of surface complex (primary or secondary) to potential (npt)    

518     FORMAT( '-->LogPotential', i2, ' on ', a11,2x,1PE11.4 )
517     FORMAT(7x, a13, 1x, f6.1)

      IF (nsurf > 0) THEN
        namtemp = 'Surface'
        WRITE(8,*)
        WRITE(8,*) 'Total Concentrations on Surface Hydroxyl Sites '
        WRITE(8,*) '---------------------------------------------- '
        SolidSolutionRatioTemp = 1000.0d0*SolidDensity(jinit(jx,jy,jz))*(1.0d0-por(jx,jy,jz))
        WRITE(8,*)
!!        WRITE(8,*) ' Surface complexes    Sites/kgw     Sites/BulkVolume(m^3)'
        WRITE(8,*) 'Surface complex   Sites/kgw      Sites/g solid  Sites/BulkVolume(m^3)'
        
        totex_bas = 0.0
        DO i = 1,ncomp  
          DO ns = 1,nsurf_sec
            totex_bas(i) = totex_bas(i) + musurf(ns,i)*spsurf10(ns+nsurf,jx,jy,jz)
          END DO

          IF (totex_bas(i) /= 0.0) THEN
            WRITE(8,515) ulab(i),totex_bas(i)/AqueousToBulk,totex_bas(i)/SolidSolutionRatioTemp,totex_bas(i)
          END IF
        END DO
      END IF 

      CALL exchange(ncomp,nexchange,nexch_sec,jx,jy,jz)
      IF (nexchange > 0) THEN
        namtemp = 'Exchange'
        WRITE(8,*)
        WRITE(8,*) 'Total Concentrations in Exchange Sites '
        WRITE(8,*) '-------------------------------------- '
        SolidSolutionRatioTemp = 1000.0d0*SolidDensity(jinit(jx,jy,jz))*(1.0d0-por(jx,jy,jz))
        WRITE(8,*)
        WRITE(8,*) 'Exchangers           Equiv/kgw        Equiv/g solid    Equiv/m^3 bulk'
        totex_bas = 0.0
        DO i = 1,ncomp  
          DO nex = 1,nexch_sec
            totex_bas(i) = totex_bas(i) + muexc(nex,i)*spex10(nex+nexchange,jx,jy,jz)
          END DO
          IF (totex_bas(i) /= 0.0) THEN
            WRITE(8,515) ulab(i),totex_bas(i)/AqueousToBulk,totex_bas(i)/SolidSolutionRatioTemp,totex_bas(i)
          END IF
        END DO
      END IF 

      WRITE(8,*)
      WRITE(8,*) 'Concentrations of Individual Species, Exchangers, and Surface Complexes '
      WRITE(8,*) '----------------------------------------------------------------------- '
      WRITE(8,*)
      
      208 FORMAT('                 ','          Log','       Log',1x,               '                                Activity   ', '     ')
      206 FORMAT('  Species        ','     Molality','  Activity',1x,               '     Molality ','   Conc (m^3) ', '  Coefficient', '     Type')

      WRITE(8,208)
      WRITE(8,206) 
      namtemp = 'Exchange'
      DO nex = 1,nexch_sec
        ix = nex + nexchange
        spprint = DLOG10(spex10(ix,jx,jy,jz)/AqueousToBulk)
        actprint10 = aexch(nex) 
        actprint = DLOG10(actprint10)
        WRITE(8,212) nam_exchsec(nex),spprint,actprint,spex10(ix,jx,jy,jz)/AqueousToBulk,spex10(ix,jx,jy,jz),actprint10,namtemp
      END DO

      namtemp = 'Surface'
      DO is = 1,nsurf
        spprint = DLOG10(spsurf10(is,jx,jy,jz)/AqueousToBulk)

        IF (LogTotalSurface(is,jx,jy,jz) /= 0.0d0) THEN
          actprint = spsurf(is,jx,jy,jz) - LogTotalSurface(is,jx,jy,jz)
          actprint10 = DEXP(actprint)
        ELSE
          actprint = 0.0d0
          actprint10 = 0.0d0 
        END IF

        WRITE(8,212) namsurf(is),spprint,actprint,spsurf10(is,jx,jy,jz)/AqueousToBulk,spsurf10(is,jx,jy,jz),actprint10,namtemp
      END DO

      DO ns = 1,nsurf_sec
        is = ns + nsurf
        spprint = DLOG10(spsurf10(is,jx,jy,jz))

        IF (LogTotalSurface(islink(ns),jx,jy,jz) /= 0.0d0) THEN
          actprint = spsurf(is,jx,jy,jz) - LogTotalSurface(islink(ns),jx,jy,jz)
          actprint10 = DEXP(actprint)
        ELSE
          actprint = 0.0d0
          actprint10 = 0.0d0 
        END IF

        WRITE(8,212) namsurf_sec(ns),spprint,actprint,spsurf10(is,jx,jy,jz)/AqueousToBulk,spsurf10(is,jx,jy,jz),actprint10,namtemp
      END DO
      
      207 FORMAT('                 ','            Log',1X,'         Log',  '              ', '              ', '               Activity')
      204 FORMAT(' Species         ','       Molality',1X,'    Activity', '         Molality ', '       Conc(m^3) ', '     Coefficient','    Type')
    
      WRITE(8,*)
      WRITE(8,207)
      WRITE(8,204) 
      namtemp = 'Aqueous'
      DO ik = 1,ncomp+nspec
        spbase = DLOG10(sp10(ik,jx,jy,jz))

          actprint = (sp(ik,jx,jy,jz)+lngamma(ik,jx,jy,jz))/clg

        actprint10 = 10**(actprint)
        WRITE(8,202) ulab(ik),spbase,actprint,sp10(ik,jx,jy,jz),actprint10,EXP(lngamma(ik,jx,jy,jz)),namtemp
      END DO  

    WRITE(8,*)
    WRITE(8,*) '*** Partial pressure gas (bars)   moles gas/m^3  g gas/m^3'
    WRITE(8,*)

    CALL GasPartialPressure(ncomp,ngas,gastmp10,jx,jy,jz)
    CALL gases(ncomp,ngas,jx,jy,jz)
    pg = PressureCond(jinit(jx,jy,jz))

    
    DO i = 1,ngas
      
      IF (namg(i) == 'H2(g)') THEN
        wtt = 2.016*spgas10(i,jx,jy,jz)
      ELSE IF (namg(i) == 'O2(g)') THEN
        wtt = 31.999*spgas10(i,jx,jy,jz)
      ELSE
        wtt = 0.0d0
      END IF
        
      WRITE(8,513)  namg(i),gastmp10(i),spgas10(i,jx,jy,jz), wtt
      
    END DO

      WRITE(8,*)
      WRITE(8,*) ' ***** Saturation state of minerals (log[Q/K] *****'
      WRITE(8,*)
      DO k = 1,nrct
        sumiap = 0.0D0
        DO i = 1,ncomp
          
          ulabPrint = ulab(i)
          IF (ulabPrint(1:3) == 'H2O' .or. ulabPrint(1:3) == 'HHO') THEN
            lnActivity = lngamma(i,jx,jy,jz)
          ELSE
            lnActivity = sp(i,jx,jy,jz) + lngamma(i,jx,jy,jz)
          END IF

          sumiap = sumiap + mumin(1,k,i)* lnActivity

        END DO
        silnTMP = sumiap - keqmin(1,k,jx,jy,jz)
        siprnt = silnTMP/clg
        WRITE(8,509) umin(k),siprnt
      END DO
      WRITE(8,*)
      WRITE(8,*)
      WRITE(8,*) ' ***** Volume of minerals added or subtracted *****'
      WRITE(8,*)
      DO k = 1,nrct
        WRITE(8,516) umin(k),volSave(k,jx,jy,jz)
      END DO
      WRITE(8,*)

    END DO
  END DO
END DO

CLOSE(UNIT=8,STATUS='keep')

502 FORMAT('temperature    ' ,f8.2)
503 FORMAT(a20,4X,1PE12.4)
504 FORMAT('END')
182 FORMAT(80(1X,1PE12.4))
183 FORMAT(1PE12.4,2X,1PE12.4)
184 FORMAT(100(1X,1PE14.6))

!2283 FORMAT('# Time (yrs) ',2X,1PE12.4)
2283 FORMAT('# Time      ',2X,1PE12.4)
2284 FORMAT('    Distance ',a18)
2282 FORMAT('    Distance ','        pH')
2281 FORMAT('    Distance ',4X,a18)
2285 FORMAT('    Distance    ',100(1X,a14))
2286 FORMAT('    Distance    ',100(1X,a14))
514 FORMAT(1X,a10,1X,1PE11.4)
513 FORMAT(1X,a18,1X,1PE12.4,2X,1PE12.4,2X,1PE12.4)
!!!515 FORMAT(1X,a18,1X,1PE12.4,5X,1PE12.4,5X,1PE12.4)

600 FORMAT(2X,f10.2,2X,a15)
201 FORMAT(2X,a18,2X,f8.2)
202 FORMAT(2X,a18,2X,f10.4,3X,f10.4,2X,1PE15.6,2X,1PE15.6,2X,1PE15.6,2x,a8)
211 FORMAT(2X,a18,2X,f8.3,3X,f8.3,2X,1PE12.3,2X,1PE12.3,2X,'            ',2x,a8)
203 FORMAT(2X,a18)

509 FORMAT(2X,a18,2X,f12.4)
510 FORMAT(2X,'GEOCHEMICAL CONDITION NUMBER',i3)

512 FORMAT(' Basis species    ','     Molality  ', '  Constraint type')
511 FORMAT(1X,a18,1X,1PE12.4,5X,a13,1X,a18)
516 FORMAT(1X,a18,1x,1PE16.9)

542  FORMAT(5X,'Primary species ',a18,1X,'at grid pt. ',i5)
544  FORMAT(5X,'Check geochem. condition # ',i2)
543  FORMAT(5X,'Constraint mineral ',2X,a18)
643  FORMAT(5X,'Constraint gas',2X,a18)
102 FORMAT(' ------> GRID LOCATION:   ',I3,':',I3,':',I3)
2287 FORMAT('#   Distance ',' alkalinity (eq/kg)')
2290 FORMAT('#  Component',5X,12(a7,7X))
2291 FORMAT('#           ',5X,6('mol/m2/yr',5X))
2292 FORMAT(a15, 9(1X,1PE13.6), 2(1X,1I2), 2(1X,1PE13.6))
2293 FORMAT('#           ',5X,3('mol/m2/yr',5X))
2294 FORMAT(a15, 9(1X,1PE13.6))
2296 FORMAT('# Net flow at top: ',1X,1PE13.6)
2297 FORMAT('# Net flow at top:    ',1X,1PE13.6)
2298 FORMAT('# Net flow at bottom: ',1X,1PE13.6)
807 FORMAT('                 ','          Log',1X,'       Log',1x,  &
'              ',1x,'          ', '       Activity')
804 FORMAT(' SPECIES         ','    Equiv/kgw',1X,'  Activity',1x,  &
'    Equiv/kgw ',1x,'    Activity ', '  Coefficient','    Type')




212 FORMAT(2X,a18,2X,f8.3,2X,f8.3,2X,1PE12.3,2X,1PE12.3,2X,1PE12.3,2X,a8)
!!!  WRITE(8,212) namsurf(is),spprint,actprint,spsurf10(is,jx,jy,jz)/AqueousToBulk,spsurf10(is,jx,jy,jz),actprint10,namtemp


RETURN
END SUBROUTINE speciation
!  *******************************************************
