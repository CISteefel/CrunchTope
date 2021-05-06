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
    
SUBROUTINE database(ncomp,ncmplx,mnrl,nrct,ngas,nsurf,nsurf_sec,ntemp,  &
  iprint,icomplete,jtemp,data1,namc,namcx,namrl,  &
  vbar,wtminfull,wtbas,coef,temp,z,bdotPrimary,zx,bdotSecondary,a0,ax0)
USE crunchtype
USE params
USE concentration
USE mineral
USE temperature, ONLY: RunIsothermal, TPointer, tinit
USE io

IMPLICIT NONE

!!EXTERNAL dgetrf
!!EXTERNAL dgetrs

!  ***************  Beginning of interface blocks  *******************************
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
  SUBROUTINE fit(nbasis,ntemp,alogk0,bvec,vec,iflgint,int,ntt,NameTransfer)
  USE crunchtype
  USE params
  IMPLICIT NONE
  INTEGER(I4B), INTENT(IN)                                   :: nbasis
  INTEGER(I4B), INTENT(IN)                                   :: ntemp
  REAL(DP), DIMENSION(ntemp), INTENT(IN)                     :: alogk0
  REAL(DP), DIMENSION(nbasis), INTENT(IN OUT)                :: bvec
  REAL(DP), DIMENSION(nbasis,ntemp), INTENT(IN)              :: vec
  INTEGER(I4B), INTENT(OUT)                                  :: iflgint
  INTEGER(I4B), DIMENSION(8), INTENT(IN OUT)                 :: int
  INTEGER(I4B), INTENT(OUT)                                  :: ntt
  CHARACTER(30), INTENT(IN)                                  :: NameTransfer
  END SUBROUTINE fit
END INTERFACE
!  ***************  End of interface blocks  *******************************
!  External arrays and variables

INTEGER(I4B), INTENT(IN)                            :: ncomp
INTEGER(I4B), INTENT(IN OUT)                        :: ncmplx
INTEGER(I4B), INTENT(IN OUT)                        :: mnrl
INTEGER(I4B), INTENT(IN)                            :: nrct
INTEGER(I4B), INTENT(IN OUT)                        :: ngas
INTEGER(I4B), INTENT(IN)                            :: nsurf
INTEGER(I4B), INTENT(IN OUT)                        :: nsurf_sec
INTEGER(I4B), INTENT(OUT)                           :: ntemp
INTEGER(I4B), INTENT(IN)                            :: iprint
INTEGER(I4B), INTENT(IN)                            :: icomplete
INTEGER(I4B), INTENT(IN)                            :: jtemp

CHARACTER (LEN=mls), INTENT(IN)                     :: data1

REAL(DP), INTENT(IN)                                :: temp

CHARACTER (LEN=mls), DIMENSION(:), INTENT(IN)       :: namc
CHARACTER (LEN=mls), DIMENSION(:), INTENT(IN OUT)   :: namcx
CHARACTER (LEN=mls), DIMENSION(:), INTENT(IN OUT)   :: namrl
REAL(DP), DIMENSION(:), INTENT(OUT)                 :: vbar
REAL(DP), DIMENSION(:), INTENT(OUT)                 :: wtminfull
REAL(DP), DIMENSION(:), INTENT(OUT)                 :: wtbas
REAL(DP), DIMENSION(:,:), INTENT(OUT)               :: coef
REAL(DP), DIMENSION(:), INTENT(OUT)                 :: z
REAL(DP), DIMENSION(:), INTENT(OUT)                 :: bdotPrimary
REAL(DP), DIMENSION(:), INTENT(OUT)                 :: zx
REAL(DP), DIMENSION(:), INTENT(OUT)                 :: bdotSecondary
REAL(DP), DIMENSION(:), INTENT(OUT)                 :: a0
REAL(DP), DIMENSION(:), INTENT(OUT)                 :: ax0

!  Internal arrays and variables

INTEGER(I4B), DIMENSION(8)                           :: int

REAL(DP), DIMENSION(8)                               :: temptmp(8)

INTEGER(I4B), PARAMETER                              :: ione=1


CHARACTER (LEN=mls)                                  :: dummy
CHARACTER (LEN=mls)                                  :: namdummy
CHARACTER (LEN=mls)                                  :: namdum
CHARACTER (LEN=mls)                                  :: namtemp
CHARACTER (LEN=mls)                                  :: namsec_tmp
CHARACTER (LEN=mls)                                  :: dumstring
CHARACTER (LEN=mls)                                  :: NAME
CHARACTER (LEN=mls)                                  :: namsave
CHARACTER (LEN=30)                                   :: NameTransfer
CHARACTER (LEN=mls),DIMENSION(:), ALLOCATABLE        :: namedummy

LOGICAL(LGT)                                         :: ok
LOGICAL(LGT)                                         :: speciesfound

!  Allocatable arrays to be deallocated upon exiting from database.f90

CHARACTER (LEN=mls), DIMENSION(:), allocatable       :: min_found
CHARACTER (LEN=mls), DIMENSION(:), allocatable       :: gas_found
CHARACTER (LEN=mls), DIMENSION(:), allocatable       :: aq_found
CHARACTER (LEN=mls), DIMENSION(:), allocatable       :: namprnt
CHARACTER (LEN=mls), DIMENSION(:), allocatable       :: namsec
CHARACTER (LEN=mls), DIMENSION(:), allocatable       :: nam

REAL(DP), DIMENSION(:,:), allocatable                :: smnrl
REAL(DP), DIMENSION(:,:), allocatable                :: stogas
REAL(DP), DIMENSION(:,:), allocatable                :: shom
REAL(DP), DIMENSION(:,:), allocatable                :: xsurf
REAL(DP), DIMENSION(:), allocatable                  :: alogk0
REAL(DP), DIMENSION(:), allocatable                  :: tempc
REAL(DP), DIMENSION(:,:), allocatable                :: vec
REAL(DP), DIMENSION(:,:), allocatable                :: vecgam
REAL(DP), DIMENSION(:), allocatable                  :: bvec
REAL(DP), DIMENSION(:), allocatable                  :: sto
REAL(DP), DIMENSION(:,:), allocatable                :: dummyreal
REAL(DP), DIMENSION(:), allocatable                  :: alogkeh
REAL(DP), DIMENSION(:), allocatable                  :: coefeh
REAL(DP), DIMENSION(:), allocatable                  :: alogk
REAL(DP), DIMENSION(:,:), allocatable                :: a
REAL(DP), DIMENSION(:,:), allocatable                :: b
REAL(DP), DIMENSION(:,:), allocatable                :: smat
REAL(DP), DIMENSION(:,:), allocatable                :: ainv
REAL(DP), DIMENSION(:), allocatable                  :: y
REAL(DP), DIMENSION(:,:), allocatable                :: w
REAL(DP), DIMENSION(:,:), allocatable                :: w_orig
REAL(DP), DIMENSION(:), allocatable                  :: work
REAL(DP), DIMENSION(:,:), allocatable                :: wgam
REAL(DP), DIMENSION(:), allocatable                  :: bb
REAL(DP), DIMENSION(:,:), allocatable                :: coef0

REAL(DP), DIMENSION(:), allocatable                  :: wtcmplx

INTEGER(I4B), DIMENSION(:), allocatable              :: iwork
INTEGER(I4B), DIMENSION(:), allocatable              :: iload
INTEGER(I4B), DIMENSION(:), allocatable              :: indx
INTEGER(I4B), DIMENSION(:), allocatable              :: indxpri
INTEGER(I4B), DIMENSION(:), allocatable              :: indxsec
INTEGER(I4B), DIMENSION(:), allocatable              :: indxmin
INTEGER(I4B), DIMENSION(:), allocatable              :: indxgas
INTEGER(I4B), DIMENSION(:), allocatable              :: ndxkin

INTEGER(I4B)                                         :: lwork
INTEGER(I4B)                                         :: n0
INTEGER(I4B)                                         :: isweep
INTEGER(I4B)                                         :: ls
INTEGER(I4B)                                         :: nsec
INTEGER(I4B)                                         :: nbasis
INTEGER(I4B)                                         :: iflgint
INTEGER(I4B)                                         :: iflgck
INTEGER(I4B)                                         :: j
INTEGER(I4B)                                         :: ncmplx_new
INTEGER(I4B)                                         :: ngas_new
INTEGER(I4B)                                         :: mnrl_new
INTEGER(I4B)                                         :: nreac
INTEGER(I4B)                                         :: nreacmin
INTEGER(I4B)                                         :: nreacgas
INTEGER(I4B)                                         :: nreacaq
INTEGER(I4B)                                         :: i
INTEGER(I4B)                                         :: l
INTEGER(I4B)                                         :: k
INTEGER(I4B)                                         :: nprimary
INTEGER(I4B)                                         :: jj
INTEGER(I4B)                                         :: ii
INTEGER(I4B)                                         :: imiss
INTEGER(I4B)                                         :: n
INTEGER(I4B)                                         :: iflg
INTEGER(I4B)                                         :: nct
INTEGER(I4B)                                         :: isecondary
INTEGER(I4B)                                         :: nsecond
INTEGER(I4B)                                         :: ngs
INTEGER(I4B)                                         :: kk
INTEGER(I4B)                                         :: ntt
INTEGER(I4B)                                         :: nlength
INTEGER(I4B)                                         :: nncnt
INTEGER(I4B)                                         :: iabove
INTEGER(I4B)                                         :: ibelow
INTEGER(I4B)                                         :: nmiss
INTEGER(I4B)                                         :: ix
INTEGER(I4B)                                         :: ixx
INTEGER(I4B)                                         :: nmat0
INTEGER(I4B)                                         :: m
INTEGER(I4B)                                         :: mm
INTEGER(I4B)                                         :: nmat
INTEGER(I4B)                                         :: lengthmin
INTEGER(I4B)                                         :: lenspec
INTEGER(I4B)                                         :: isweep_surf
INTEGER(I4B)                                         :: icantfind
INTEGER(I4B)                                         :: ifind
INTEGER(I4B)                                         :: is
INTEGER(I4B)                                         :: ns
INTEGER(I4B)                                         :: nreacsurf
INTEGER(I4B)                                         :: nn0
INTEGER(I4B)                                         :: ngss
INTEGER(I4B)                                         :: nss
INTEGER(I4B)                                         :: ik
INTEGER(I4B)                                         :: msub
INTEGER(I4B)                                         :: np
INTEGER(I4B)                                         :: iflgstop
INTEGER(I4B)                                         :: icheck
INTEGER(I4B)                                         :: info

REAL(DP)                                             :: flogk
REAL(DP)                                             :: vbargas
REAL(DP)                                             :: sum
REAL(DP)                                             :: dd
REAL(DP)                                             :: AA0
REAL(DP)                                             :: zz
REAL(DP)                                             :: wtt
REAL(DP)                                             :: vbar0
REAL(DP)                                             :: alogsum
REAL(DP)                                             :: ferr
REAL(DP)                                             :: berr
REAL(DP)                                             :: charge

REAL(DP), PARAMETER                                  :: eps=1.E-12
REAL(DP), PARAMETER                                  :: tiny=1.E-12
REAL(DP), PARAMETER                                  :: small=1.E-05
REAL(DP), PARAMETER                                  :: tk=273.15

CHARACTER (LEN=1)                                    :: trans

CHARACTER (LEN=mls)                                  :: parchar
CHARACTER (LEN=mls)                                  :: parfind
INTEGER(I4B)                                         :: lchar
CHARACTER (LEN=mls)                                  :: section
REAL(DP)                                             :: bdotDummy               

trans = 'N'

!     subroutine database orders an independent set of chemical
!     reactions in terms of an arbitrary set of components.

n0 = ndim

ALLOCATE(min_found(nm))
ALLOCATE(gas_found(ng))
ALLOCATE(aq_found(mcmplx+nc))
ALLOCATE(namprnt(nc)) 
ALLOCATE(namsec(ndim+1))
ALLOCATE(nam(50))
ALLOCATE(smnrl(nc,nm))
ALLOCATE(stogas(nc,ng))
ALLOCATE(shom(nc,mcmplx))
ALLOCATE(xsurf(nc,msurf_sec))
ALLOCATE(alogk0(ntmp))
ALLOCATE(tempc(ntmp))
ALLOCATE(bvec(5))
ALLOCATE(sto(50))
ALLOCATE(dummyreal(ndim,5))
ALLOCATE(alogkeh(ntmp))
ALLOCATE(coefeh(5))
ALLOCATE(alogk(ndim))
ALLOCATE(a(ndim,ndim))
ALLOCATE(b(ndim,nc))
ALLOCATE(smat(ndim,ndim))
!!ALLOCATE(ainv(ndim,ndim))
!!ALLOCATE(y(ndim))
ALLOCATE(vec(5,8))
ALLOCATE(vecgam(5,8))

ALLOCATE(wgam(5,5))
ALLOCATE(bb(ndim))
ALLOCATE(coef0(ndim,5))
ALLOCATE(wtcmplx(mcmplx))
ALLOCATE(iload(mcmplx))
ALLOCATE(indxpri(nc))
ALLOCATE(indxsec(mcmplx))
ALLOCATE(indxmin(nm))
ALLOCATE(indxgas(ng))
ALLOCATE(ndxkin(nm))


min_found = ' '
gas_found = ' '
aq_found = ' '
namprnt = ' '
namsec = ' '
nam = ' '

smnrl = 0.0
stogas = 0.0
shom = 0.0
xsurf = 0.0
alogk0 = 0.0
tempc = 0.0
vec = 0.0
vecgam = 0.0
bvec = 0.0
sto = 0.0
dummyreal = 0.0
alogkeh = 0.0
coefeh = 0.0
alogk = 0.0
a = 0.0
b = 0.0
smat = 0.0
!!ainv = 0.0
!!y = 0.0
wgam = 0.0
bb = 0.0
coef0 = 0.0
wtbas = 0.0
wtcmplx = 0.0

iload = 0
indxpri = 0
indxsec = 0
indxmin = 0
indxgas = 0
ndxkin = 0
icantfind = 0

iunit5 = 8
isweep = 0

!     open appropriate database file

WRITE(iunit2,*)
IF (data1 /= ' ') THEN
  CALL stringlen(data1,ls)
!!  WRITE(*,*) ' Using database file: ', data1(1:ls)
  WRITE(iunit2,*) ' Using database file: ', data1(1:ls)
  
  OPEN(UNIT=iunit5,FILE=data1,STATUS='old',ERR=334)
  
  
!        open(unit=iunit5,file=data1,status='old',err=334,access='sequential',recl=1024)
  
!      else if (temp .eq. 25.d0 .and. jtemp.eq.0) then
! write(*,*) '        --> using database: master25'
! write(iunit2,*)'        --> using database: master25'
!        open(iunit5, file=
!     .  '/users/steefel/threedata/master25.data',
!     .  status='old',err=334)
ELSE
  WRITE(*,*)
  WRITE(*,*) ' No default database:  must be specified in input file'
  WRITE(*,*)
  READ(*,*)
  STOP
END IF
WRITE(iunit2,*)

!!CIS 8/7/07 nsec = ncmplx + mnrl + ngas 
nsec = ncmplx + mnrl + ngas 
nbasis = 5
iflgint = 0
iflgck = 0

DO j = 1, ncomp
  indxpri(j)=0
END DO
DO j = 1, ncmplx
  indxsec(j)=0
END DO
DO j = 1, mnrl
  indxmin(j)=0
END DO
DO j = 1, ngas
  indxgas(j)=0
END DO

ncmplx_new = ncmplx
ngas_new = ngas
mnrl_new = mnrl
nreac = 0
nreacsurf = 0
nreacmin = 0
nreacgas = 0
nreacaq = 0


!     initialize matrices for log K fit (coef) and
!     for stoichiometric coefficients (a)

coef = 0.d0
a = 0.d0

!-----read data base

section = 'ReadBdot'
IncludeBdot = .FALSE.
parchar = 'includebdot'
parfind = ' '
CALL readSingleString(iunit5,parchar,parfind,section)
IF (parfind == 'includebdot') THEN  
  IncludeBdot = .TRUE.    
ELSE
  IncludeBdot = .FALSE.
  REWIND(iunit5)
END IF

READ(iunit5,*,ERR = 6001) NAME,ntemp,(tempc(l),l=1,ntemp)
IF (ntemp > ntmp) THEN
  WRITE(*,*) 'too many temperature points in database!'
  READ(*,*)
  STOP
END IF

IF (ALLOCATED(DatabaseTemperature)) THEN
  DEALLOCATE(DatabaseTemperature)
END IF
ALLOCATE(DatabaseTemperature(ntemp))

DatabaseTemperature(1:ntemp) = tempc(1:ntemp)

!!WRITE(*,*) ' Successfully read temperatures'

READ(iunit5,*,ERR=7001) NAME,(adh(i),i=1,ntemp)
READ(iunit5,*,ERR=7001) NAME,(bdh(i),i=1,ntemp)
READ(iunit5,*,ERR=7001) NAME,(bdot(i),i=1,ntemp)

!!WRITE(*,*) ' Successfully read Debye-Huckel parameters'

!  If there is more than a single temperature, fit the log K's to the
!  polynomial which follows

IF (ntemp > 1) THEN

  DO i = 1, ntemp
    vec(1,i) = DLOG(tempc(i) + tk)
    vec(2,i) = 1.d0
    vec(3,i) = tempc(i) + tk
    vec(4,i) = 1.d0/(tempc(i) + tk)
    vec(5,i) = 1.d0/((tempc(i) + tk)*(tempc(i) + tk))
  END DO
  
!  ALLOCATE(w(nbasis,nbasis))
!  ALLOCATE(indx(nbasis))

!  DO j = 1, nbasis
!    DO k = j, nbasis
!      w(j,k) = 0.d0
!      DO i = 1, ntemp
!        w(j,k) = w(j,k) + vec(j,i)*vec(k,i)
!      END DO
!      IF (j /= k) w(k,j) = w(j,k)
!    END DO
!  END DO
  
!  CALL ludcmp90(w,indx,dd,nbasis)
!  CALL ludcmp(w,nbasis,n0,indx,dd)
 
  
!  ***********  Calculate coefficients for Debye-Huckel parameters  *****************
  
  DO i = 1, ntemp
    vecgam(1,i) = 1.0
    vecgam(2,i) = tempc(i)
    vecgam(3,i) = (tempc(i))**2
    vecgam(4,i) = (tempc(i))**3
    vecgam(5,i) = (tempc(i))**4
  END DO
  
  !DO j = 1, nbasis
  !  DO k = j, nbasis
  !    wgam(j,k) = 0.d0
  !    DO i = 1, ntemp
  !      wgam(j,k) = wgam(j,k) + vecgam(j,i)*vecgam(k,i)
  !    END DO
  !    IF (j /= k) wgam(k,j) = wgam(j,k)
  !  END DO
  !END DO
  
  !CALL ludcmp(wgam,nbasis,n0,indx,dd)
  
  CALL fitgamma(nbasis,ntemp,adh,bvec,vecgam)
  DO i = 1, nbasis
    adhcoeff(i) = bvec(i)
  END DO
  
  CALL fitgamma(nbasis,ntemp,bdh,bvec,vecgam)
  DO i = 1, nbasis
    bdhcoeff(i) = bvec(i)
  END DO
  
  CALL fitgamma(nbasis,ntemp,bdot,bvec,vecgam)
  DO i = 1, nbasis
    bdtcoeff(i) = bvec(i)
  END DO

  
!  ********   End of Debye-Huckel coefficient calculation  ************
  
END IF

!*********************************************************
!  Start reading database here

!-----read block #1   Reading primary species in EQ3 database

!-----read and store aqueous species properties
5001 nprimary = 0
500 CONTINUE

IF (IncludeBdot) THEN
  READ(iunit5,*,ERR=6002) NAME,aa0,zz,wtt,bdotDummy
ELSE
  READ(iunit5,*,ERR=6002) NAME,aa0,zz,wtt
END IF

IF (NAME == 'End of primary') GO TO 501   !  End of primary species block
DO jj = 1, ncomp   ! Search through list of primary species in input file
  j = jj
  IF (NAME == namc(jj)) GO TO 530
END DO
DO ii = 1, ncmplx   ! Search through list of secondary species in input file
  i = ii
  IF (NAME == namcx(i)) GO TO 532
END DO
GO TO 500
530 CONTINUE
nprimary = nprimary + 1
aq_found(nprimary) = NAME
z(j) = zz
a0(j) = aa0
wtbas(j) = wtt
IF (IncludeBdot) THEN
  bdotPrimary(j) = bdotDummy
END IF
GO TO 500
532 CONTINUE
nprimary = nprimary + 1
aq_found(nprimary) = NAME
zx(i) = zz
ax0(i) = aa0
wtcmplx(i) = wtt
IF (IncludeBdot) THEN
  bdotSecondary(i) = bdotDummy
END IF
GO TO 500

501 CONTINUE   !  End of primary species block encountered, so get out

!*************************************************************************
!-----read block #2   --Now read through the EQ3 "secondary" species list

!-----read and store aqueous reactions


30 CONTINUE

imiss = 0
IF (IncludeBdot) THEN
  READ(iunit5,*,ERR=6003) nam(1),n,(sto(i+1),nam(i+1), i = 1, n),  &
    (alogk0(l),l=1,ntemp),aa0,zz,wtt,bdotDummy
ELSE
  READ(iunit5,*,ERR=6003) nam(1),n,(sto(i+1),nam(i+1), i = 1, n),  &
    (alogk0(l),l=1,ntemp),aa0,zz,wtt
END IF

IF (nam(1) == 'End of secondary' ) GO TO 29

sto(1) = -1    ! Secondary species in EQ3 list has stoichiometric coefficient of -1

!     ----decide if reaction occurs in chosen system

iflg = 0
nct = 0
isecondary = 0
nsecond = 0
DO  k = n+1,1,-1

  speciesfound = .FALSE.
  
  IF (nam(k) == 'h2o' .OR. nam(k) == 'H2O') THEN
    nct = nct + 1
    speciesfound = .TRUE.
  END IF

  IF (speciesfound) then
    CYCLE
  END IF

  DO j = 1, ncomp
    IF (nam(k) == namc(j)) THEN
      nct = nct + 1
      indxpri(j) = 1
      IF (k == 1) THEN  ! Secondary species in EQ3, primary species in input file
        iflg = 1
        jj = j
      END IF
      speciesfound = .TRUE.
      EXIT
    END IF
  END DO

  IF (speciesfound) then
    CYCLE
  END IF
  
  DO i = 1, ncmplx
    IF (nam(k) == namcx(i)) THEN
      nsecond = nsecond + 1
      namsec_tmp = nam(k)
      isecondary = 1
      nct = nct + 1
      indxsec(i) = 1
      IF (k == 1) THEN  ! Secondary species in EQ3,secondary species in input file
        iflg = 2
        ii = i
      END IF
      speciesfound = .TRUE.
      EXIT
    END IF
  END DO

  IF (speciesfound) then
    CYCLE
  END IF
  
  DO ngs = 1, ngas
    IF (nam(k) == namg(ngs)) THEN
      nct = nct + 1
      indxgas(ngs) = 1
      speciesfound = .TRUE.
      EXIT
    END IF
  END DO

  IF (speciesfound) then
    CYCLE
  END IF
  
  imiss = imiss + 1
  namsave = nam(k) 

END DO

!  Check to see if only one species in the reaction is missing.  If so, add
!  this species to the secondary species list.
!  Do this only if "icomplete" = 1, otherwise rely on the
!  list in the input file

IF (icomplete == 1) THEN
  IF (isweep == 3) THEN
    IF (imiss == 0) THEN
      CONTINUE
    ELSE
      GO TO 30
    END IF
  ELSE
    
    IF (imiss == 1 .AND. namsave /= 'O2(g)' .AND. namsave /= 'o2(g)') THEN
      
!            if (ntemp.gt.1) then
!              do k = 1,ntemp
!                if (alogk0(k).eq.500.0) then
!                  goto 30
!                endif
!              end do
!            endif
      
      ncmplx_new = ncmplx_new + 1
      IF (ncmplx_new > mcmplx) THEN
        WRITE(*,*)
        WRITE(*,*) ' MCMPLX not dimensioned large enough'
        WRITE(*,*) ' In params.inc'
        WRITE(*,*) ' Number of secondary species = ',ncmplx_new
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
!!      WRITE(*,*)
!!      WRITE(*,*) ' Adding a secondary species from the EQ3 database'
!!      WRITE(*,*) nam(1)
!!      WRITE(*,*) ' Number of secondary species = ',ncmplx_new
      namcx(ncmplx_new) = namsave
!!      WRITE(*,*)
      GO TO 30
    ELSE IF (imiss >= 2) THEN
      GO TO 30
    ELSE
      GO TO 30
    END IF
  END IF
ELSE
  IF (imiss == 0) THEN
    CONTINUE
  ELSE
    GO TO 30
  END IF
END IF

29  IF (icomplete == 1 .AND. isweep < 3) GO TO 41

IF (nam(1) == 'End of secondary') GO TO 39


IF (nsecond == 1) THEN
  
!  Only one secondary species in reaction--check to see if it is already loaded
  
  DO kk = 1,ncmplx
    IF (namcx(kk) == namsec_tmp) THEN
      IF (iload(kk) == 1) THEN    ! Species already loaded--skip reaction
!!        WRITE(*,*) ' Secondary species already loaded'
!!        WRITE(*,*) ' Reaction involves the species ',namsec_tmp
!!        WRITE(*,*) ' Skipping reaction'
!!        WRITE(*,*)
        GO TO 30
      ELSE
!!        WRITE(*,*) ' Secondary species found for first time'
!!        WRITE(*,*) ' Reaction involves the species ',namsec_tmp
!!        WRITE(*,*) ' Loading reaction'
!!        WRITE(*,*)
        iload(kk) = 1
        GO TO 3333
      END IF
    END IF
  END DO
END IF

3333   IF (isecondary == 1) THEN
  nreac = nreac + 1
ELSE
!!  WRITE(*,*)
!!  WRITE(*,*) ' Reaction made up only of primary species'
!!  WRITE(*,*) ' Reaction involves the species ',nam(1)
!!  WRITE(*,*)
  z(jj) = zz
  a0(jj) = aa0
  wtbas(jj) = wtt
  IF (IncludeBdot) THEN
    bdotPrimary(jj) = bdotDummy
  END IF
  ALLOCATE(namedummy(nreacaq))
  DO i = 1,nreacaq
    namedummy(i) = aq_found(nprimary+i)
  END DO
  nprimary = nprimary + 1
  aq_found(nprimary) = nam(1)
  DO i = 1,nreacaq
    aq_found(nprimary+i) = namedummy(i) 
  END DO  
  DEALLOCATE(namedummy)
  GO TO 30
END IF

namsec(nreac) = nam(1)
nreacaq = nreacaq + 1
aq_found(nreacaq+nprimary) = nam(1)

IF (nct /= n+1) THEN
  WRITE(iunit2,*) 'species missing in reaction: ',nreac,namsec(nreac)
  WRITE(*,*) 'species missing in reaction',nreac,namsec(nreac)
  READ(*,*)
  STOP
END IF

IF (ntemp > 1) THEN
  
!  Check to see if there are any logK = 500 in advance of call to subroutine fit
  
!        do i = 1,ntemp
!          if (ntemp.gt.1 .and. alogk0(i).eq.500.0) then
!            iflgint = 1
!            write(iunit2,*)
!            write(*,*)
!     write(iunit2,*) ' LogK = 500 encountered for species '
!     &    ,   namsec(nreac)
!     write(*,*) ' LogK = 500 encountered for species '
!     &    ,   namsec(nreac)
!            write(*,*) ' Non-isothermal problem--logKs missing'
!            write(iunit2,*) ' Non-isothermal problem--logKs missing'
!            write(*,*) ' Aborting run'
!            write(iunit2,*) ' Aborting run'
!            write(iunit2,*)
!            write(*,*)
!            stop
!          endif
!        end do
  
  IF (RunIsothermal) THEN
      
!!  Find the right logK value in the database for the specified temperature
    TPointer = 0
    DO i = 1,ntemp
        IF (tempc(i) == tinit) THEN
            TPointer = i
            alogk(nreac) = alogk0(TPointer)
        END IF
    END DO
    IF (Tpointer == 0) THEN
       WRITE(*,*)
       WRITE(*,*) ' Single temperature value for isothermal run not found among temperature points in database'
       WRITE(*,*)
       READ(*,*)
       STOP
    END IF
    
  ELSE

    NameTransfer = namsec(nreac)
    CALL fit(nbasis,ntemp,alogk0,bvec,vec,iflgint,INT,ntt,NameTransfer)
  
    IF (iflgint == 1) THEN
      dumstring = namsec(nreac)
      CALL stringlen(dumstring,nlength)
!    WRITE(iunit2,*)
!    WRITE(iunit2,*) '            WARNING  '
!    WRITE(iunit2,*) ' Some logKs missing for: ',dumstring(1:nlength)
!    WRITE(iunit2,*)
!    WRITE(*,*)
!    WRITE(*,*) '              WARNING:  '
!    WRITE(*,*) ' Some logKs missing for: ',dumstring(1:nlength)
      nncnt = 0
      nncnt = 0
      DO i = 1,ntemp
        IF (INT(i) == 1) THEN
          nncnt = nncnt + 1
          temptmp(nncnt) = tempc(i)
        END IF
      END DO
      IF (ntt == 1) THEN
        IF (temp /= temptmp(1)) THEN
          WRITE(*,*) ' Only one logK in database at T(C): ',temptmp(1)
          WRITE(*,*) ' Temperature of condition:          ',temp
          WRITE(*,*) ' Species missing log Ks: ', dumstring(1:nlength)
          WRITE(*,*)
          STOP
        END IF
      ELSE
        iabove = 0
        ibelow = 0
        DO i = 1,ntt
!  See if TEMP falls above at least one of the temperatures in the database
          IF (temp >= temptmp(i)) THEN
            ibelow = 1
          END IF
          IF (temp <= temptmp(i)) THEN
            iabove = 1
          END IF
        END DO
        IF (iabove == 0) THEN
          namdummy = namsec(nreac)
          CALL stringlen(namdummy,ls)
          WRITE(*,*)
          WRITE(*,*) ' Temperature for condition outside temperature range for: ',namdummy(1:ls)
          WRITE(*,*) ' Temperature of condition:  ',temp
          WRITE(*,*)
          WRITE(*,*) ' Thermodynamic data for following temperatures available '
          DO i = 1,ntt
            WRITE(*,909) temptmp(i)
          END DO
          WRITE(*,*)
          READ(*,*)
          STOP
        END IF
      909  FORMAT(1X,f10.3)
        IF (ibelow == 0) THEN
          namdummy = namsec(nreac)
          CALL stringlen(namdummy,ls)
          WRITE(*,*)
          WRITE(*,*) ' Temperature for condition outside temperature range for: ',namdummy(1:ls)
          WRITE(*,*) ' Temperature of condition:  ',temp
          WRITE(*,*)
          WRITE(*,*) ' Thermodynamic data for following temperatures available '
          DO i = 1,ntt
            WRITE(*,909) temptmp(i)
          END DO
          WRITE(*,*)
          READ(*,*)
          STOP
        END IF
      END IF
     
    END IF
  
    DO j = 1, nbasis
      coef0(nreac,j) = bvec(j)
    END DO
  
    alogk(nreac) = flogk(bvec,temp)
    
  END IF    !! End of isothermal/nonisothermal block for case where ntemp > 1
  
ELSE IF (ntemp == 1) THEN
    
  alogk(nreac) = alogk0(1)
  
END IF
  

IF (iflg == 1) THEN   ! First species in reaction designated as primary species in
!                               input file
  z(jj) = zz
  a0(jj) = aa0
  wtbas(jj) = wtt
  IF (IncludeBdot) THEN
    bdotPrimary(jj) = bdotDummy
  END IF
  
ELSE IF (iflg == 2) THEN  ! First species in reaction designated as a
!                                     secondary species in input file
  zx(ii) = zz
  ax0(ii) = aa0
  wtcmplx(ii) = wtt
  IF (IncludeBdot) THEN
    bdotSecondary(ii) = bdotDummy
  END IF
END IF

!     ----rearrange reactions - check for primary species

DO  i = 1, n+1        ! list of species in reaction
  speciesfound = .FALSE.
  DO j = 1, ncomp    ! list of primary species in input file
    IF (nam(i) == namc(j)) THEN
      b(nreac,j) = sto(i)    ! stoichiometric coefficient for primary species
      speciesfound = .TRUE.
      EXIT
    END IF
  END DO
  
  IF (speciesfound) then
    CYCLE
  END IF
  
!     ----check for secondary species
  
  DO k = 1, ncmplx   ! list of secondary species in input file
    IF (nam(i) == namcx(k)) THEN
      a(nreac,k) = -sto(i)    ! stoichiometric coefficient for secondary species
      speciesfound = .TRUE.
      EXIT
    END IF
  END DO

  IF (speciesfound) then
    CYCLE
  END IF
  
  DO ngs = 1, ngas   ! list of gases in input file
    IF (nam(i) == namg(ngs)) THEN
!       a(nreac,ncmplx+mnrl+ngs) = -sto(i)  ! stoichiometric coefficient for gas
      a(nreac,ncmplx+ngs) = -sto(i)  ! stoichiometric coefficient for gas
      speciesfound = .TRUE.
      EXIT
    END IF
  END DO

  IF (speciesfound) then
    CYCLE
  END IF

END DO

GO TO 30

39 CONTINUE

IF (nreac /= ncmplx) iflgck = 1
IF (iflgck == 1) THEN
  nmiss = ncmplx - nreac
  WRITE(*,*) ' Finished read of secondary species block'
  WRITE(*,*)
  WRITE(*,*) ' Number of reactions =            ',nreac
  WRITE(*,*) ' Number of secondary species =    ',ncmplx
  WRITE(*,*)
  WRITE(*,*)
  WRITE(*,*) ' Could not find one or more secondary species reactions'
  WRITE(*,*) ' Number of reactions missing =    ',nmiss
  WRITE(*,*)
  DO k = 1,ncmplx
    DO kk = 1,nreacaq+nprimary
      IF (aq_found(kk) == namcx(k) ) THEN
        GO TO 5113
      END IF
    END DO
    namtemp = namcx(k)
    CALL stringlen(namtemp,ls)
    WRITE(*,*) ' Secondary species in input file not found in database: ',namtemp(1:ls)
    5113     CONTINUE
  END DO
  WRITE(*,*)
  DO k = 1,ncomp
    DO kk = 1,nreacaq+nprimary
      IF (aq_found(kk) == namc(k) ) THEN
        GO TO 5114
      END IF
    END DO
    namtemp = namc(k)
    CALL stringlen(namtemp,ls)
    WRITE(*,*) ' Primary species in input file not found in database: ',namtemp(1:ls)
    5114     CONTINUE
  END DO
  WRITE(*,*)
  READ(*,*)
  STOP
END IF

!**********************************************************************
!-----read block #4
!-----read and store gas reactions

41 CONTINUE

READ(iunit5,*,ERR=6005) nam(1),vbargas,n,(sto(i+1),nam(i+1),i=1,n),  &
    (alogk0(l),l=1,ntemp),wtt

IF (nam(1) == 'End of gases') GO TO 52

nct = 0

sto(1) = -1.

!     ----decide if reaction occurs in chosen system

imiss = 0
DO  i = n+1,1,-1
  speciesfound = .FALSE.
  
!       check for water
  IF (nam(i) == 'h2o' .OR. nam(i) == 'H2O') THEN
    nct = nct + 1
    speciesfound = .TRUE.
  END IF

  IF (speciesfound) then
    CYCLE
  END IF
  
  DO  ngs = 1, ngas
    IF (nam(i) == namg(ngs)) THEN
      nct = nct + 1
      indxgas(ngs) = 1
!              write(*,*) ' Gas = ',namg(ngs)
      speciesfound = .TRUE.
      EXIT
    END IF
  END DO

  IF (speciesfound) then
    CYCLE
  END IF
  
  DO  j = 1, ncomp
    IF (nam(i) == namc(j)) THEN
      nct = nct + 1
!              write(*,*) ' Primary species = ',namc(j)
      indxpri(j) = 1
      speciesfound = .TRUE.
      EXIT
    END IF
  END DO

  IF (speciesfound) then
    CYCLE
  END IF
  
  DO  ix = 1, ncmplx
    IF (nam(i) == namcx(ix)) THEN
      ixx = ix
      nct = nct + 1
!              write(*,*) ' Secondary species = ',namcx(ix)
      indxsec(ix) = 1
      speciesfound = .TRUE.
      EXIT
    END IF
  END DO

  IF (speciesfound) then
    CYCLE
  END IF
  
  403 CONTINUE
  
  imiss = imiss + 1
  namsave = nam(i)
  
END DO 

IF (icomplete == 1) THEN
  IF (isweep == 3) THEN
    IF (imiss == 0) THEN
      CONTINUE
    ELSE
      GO TO 41
    END IF
  ELSE
    
    IF (imiss == 1) THEN
      
!            if (ntemp.gt.1) then
!              do k = 1,ntemp
!                if (alogk0(k).eq.500.0) then
!                  goto 41
!                endif
!              end do
!            endif
      
      ngas_new = ngas_new + 1
      IF (ngas_new > ng) THEN
        WRITE(*,*)
        WRITE(*,*) ' NG not dimensioned large enough'
        WRITE(*,*) ' In params.inc'
        WRITE(*,*) ' Number of gases = ',ngas_new
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
      WRITE(*,*)
      WRITE(*,*) ' Adding a gas from the EQ3 database'
      WRITE(*,*) nam(1)
      WRITE(*,*) ' Number of gases = ',ngas_new
      namg(ngas_new) = namsave
      WRITE(*,*)
      GO TO 41
    ELSE IF (imiss >= 2) THEN
      GO TO 41
    ELSE
      GO TO 41
    END IF
  END IF
ELSE
  IF (imiss == 0) THEN
    CONTINUE
  ELSE
    GO TO 41
  END IF
END IF

52  IF (icomplete == 1 .AND. isweep < 3) THEN
  IF (isweep == 0) THEN
    WRITE(*,*)
    WRITE(*,*) ' Finished first sweep of EQ3 database'
    WRITE(*,*) ' Now carrying out a second sweep to add additional species'
    WRITE(*,*)
    ncmplx = ncmplx_new
    ngas = ngas_new
    nsec = ncmplx + mnrl + ngas
    isweep = 1
    nreac = 0
    REWIND(iunit5)
    READ(iunit5,*) dummy
    GO TO 5001
  ELSE IF (isweep == 1) THEN
    WRITE(*,*)
    WRITE(*,*) ' Finished second sweep of EQ3 database'
    WRITE(*,*) ' Now repeating sweep to add additional species'
    WRITE(*,*)
    ncmplx = ncmplx_new
    ngas = ngas_new
    nsec = ncmplx + mnrl + ngas
    isweep = 2
    nreac = 0
    REWIND(iunit5)
    READ(iunit5,*) dummy
    GO TO 5001
  ELSE
    WRITE(*,*)
    WRITE(*,*) ' Finished second sweep of EQ3 database'
    WRITE(*,*) ' Now repeating sweep to add additional species'
    WRITE(*,*)
    ncmplx = ncmplx_new
    ngas = ngas_new
    nsec = ncmplx + mnrl + ngas
    isweep = 3
    nreac = 0
    REWIND(iunit5)
    READ(iunit5,*) dummy
    GO TO 5001
  END IF
END IF

IF (nam(1) == 'End of gases') GO TO 42

IF (nct /= n+1) GO TO 41

nreac = nreac + 1
nreacgas = nreacgas + 1
gas_found(nreacgas) = nam(1)

namsec(nreac) = nam(1)

IF (ntemp > 1) THEN
    
  IF (RunIsothermal) THEN
      
!!  Find the right logK value in the database for the specified temperature
    TPointer = 0
    DO i = 1,ntemp
        IF (tempc(i) == tinit) THEN
            TPointer = i
            alogk(nreac) = alogk0(TPointer)
        END IF
    END DO
    IF (Tpointer == 0) THEN
       WRITE(*,*)
       WRITE(*,*) ' Single temperature value for isothermal run not found among temperature points in database'
       WRITE(*,*)
       READ(*,*)
       STOP
    END IF
    
  ELSE    
    
    NameTransfer = namsec(nreac)
    CALL fit(nbasis,ntemp,alogk0,bvec,vec,iflgint,INT,ntt,NameTransfer)

    IF (iflgint == 1) THEN
!    WRITE(iunit2,*)
!    WRITE(*,*)
!    WRITE(iunit2,*) 'WARNING! LogK = 500 encountered for species '  &
!        ,namsec(nreac)
!    WRITE(*,*) 'WARNING! LogK = 500 encountered for species ' ,namsec(nreac)
    END IF
  
    DO  j = 1, nbasis
      coef0(nreac,j) = bvec(j)
    END DO
  
    alogk(nreac) = flogk(bvec,temp)
    
  END IF
  
ELSE IF (ntemp == 1) THEN
    
  alogk(nreac) = alogk0(1)
  
END IF

DO  i = 1, n+1

  speciesfound = .FALSE.

!     ----rearrange reactions - check for primary species
  DO  j = 1, ncomp
    IF (nam(i) == namc(j)) THEN
      b(nreac,j) = sto(i)
      speciesfound = .TRUE.
      EXIT
    END IF
  END DO
  
  IF (speciesfound) then
    CYCLE
  END IF
!     ----check for secondary species
  
  DO  k = 1, ncmplx
    IF (nam(i) == namcx(k)) THEN
      a(nreac,k) = -sto(i)
      speciesfound = .TRUE.
      EXIT
    END IF
  END DO

  IF (speciesfound) then
    CYCLE
  END IF
  
  DO  ngs = 1, ngas
    IF (nam(i) == namg(ngs)) THEN
!       a(nreac,ncmplx+mnrl+ngs) = -sto(i)
      a(nreac,ncmplx+ngs) = -sto(i)
      speciesfound = .TRUE.
      EXIT
    END IF
  END DO

  IF (speciesfound) then
    CYCLE
  END IF
  
END DO 

GO TO 41

42 CONTINUE

IF (nreac /= ncmplx+ngas) iflgck = 1
IF (iflgck == 1) THEN
  WRITE(*,*)
  WRITE(*,*) ' Mismatch between gases specified and reactions found'
  WRITE(*,*) ' Number of gas reactions =         ',nreacgas
  WRITE(*,*) ' Number of gases =                 ',ngas
  WRITE(*,*)
  DO k = 1,ngas
    DO kk = 1,nreacgas
      IF (gas_found(kk) == namg(k)) THEN
        GO TO 5112
      END IF
    END DO
    namtemp = namg(k)
    call STRINGLEN(namtemp,ls)
    WRITE(*,*)
    WRITE(*,*) ' Gas reaction not found in thermo database: ',namtemp(1:ls)
    WRITE(*,*)
    5112     CONTINUE
  END DO
  WRITE(*,*)
  READ(*,*)
  STOP
END IF

isweep = 0
!****************************************************************
!                       ***MINERALS****

!-----read block #3   ! Read minerals block in the EQ3 database
!-----read and store mineral reactions

nmat0 = nreac

IF (icomplete == 1) THEN
  
  31   CONTINUE
  
  READ(iunit5,*,ERR=6004) nam(1),vbar0,n,(sto(i+1),nam(i+1),i=1,n),  &
      (alogk0(l),l=1,ntemp),wtt
  
  sto(1) = -1.
  nct = 0
  
  IF (nam(1) == 'End of minerals') GO TO 32
  
!     ----decide if reaction occurs in chosen system
  
  DO  i = n+1,1,-1

    speciesfound = .FALSE.

!       check for water
    IF (nam(i) == 'h2o' .OR. nam(i) == 'H2O') THEN
      nct = nct + 1
     speciesfound = .TRUE.
    END IF

    IF (speciesfound) THEN
      CYCLE
    END IF
    
    DO m = 1, mnrl
      IF (nam(i) == namrl(m)) THEN
        nct = nct + 1
        indxmin(m) = 1
        mm = m
        speciesfound = .TRUE.
        EXIT
      END IF
    END DO

    IF (speciesfound) THEN
      CYCLE
    END IF
    
    DO j = 1, ncomp
      IF (nam(i) == namc(j)) THEN
        nct = nct + 1
        indxpri(j) = 1
        speciesfound = .TRUE.
        EXIT
      END IF
    END DO

    IF (speciesfound) THEN
      CYCLE
    END IF
    
    DO ix = 1, ncmplx
      IF (nam(i) == namcx(ix)) THEN
        nct = nct + 1
        indxsec(ix) = 1
        speciesfound = .TRUE.
        EXIT
      END IF
    END DO

    IF (speciesfound) THEN
      CYCLE
    END IF
    
    DO ngs = 1, ngas
      IF (nam(i) == namg(ngs)) THEN
        nct = nct + 1
        indxgas(ngs) = 1
        speciesfound = .TRUE.
        EXIT
      END IF
    END DO

    IF (speciesfound) THEN
      CYCLE
    END IF
    
    IF (i == 1) THEN
      
      mnrl_new = mnrl_new + 1
      IF (mnrl_new > nm) THEN
        WRITE(*,*) ' NM not dimensioned large enough in params.inc'
        WRITE(*,*) ' Number of minerals = ',mnrl_new
        READ(*,*)
        STOP
      END IF
      WRITE(*,*)
      WRITE(*,*) ' Adding a mineral from the EQ3 database'
      WRITE(*,*) nam(1)
      WRITE(*,*) ' Number of minerals = ',mnrl_new
      namrl(mnrl_new) = nam(1)
      WRITE(*,*)
    END IF
    
    503     CONTINUE
    
!     ----mineral not found in first pass - skip reaction for now
    
    GO TO 31
    

  END DO
  
!  Mineral found, add to reaction list
  
  nreac = nreac + 1
  nreacmin = nreacmin + 1
  min_found(nreacmin) = nam(1)
  IF (iprint == 2) WRITE(iunit2,*) nreac,nam(1)
  
  nmat = nmat0 + mm
  
  namsec(nmat) = nam(1)
  
  IF (nct /= n+1) THEN
    WRITE(iunit2,*) 'species missing in reaction: ',nmat, namsec(nmat)
    WRITE(*,*) 'species missing in reaction',nmat, namsec(nmat)
    READ(*,*)
    STOP
  END IF
  
  IF (icomplete /= 1) THEN
    if (mintype(mm) == 0) then
      vbar(mm) = vbar0*1.d-6  ! Convert to m**3/mole (Changed by Steefel)
    else if (mintype(mm) == 1) then
      vbar(mm) = 1.0d0 ! mole biomass/mole biomass  ! obsolete: vbar0 cells/mole
    end if
  ELSE
    vbar(mm) = vbar0*1.d-6  ! Convert to m**3/mole (Changed by Steefel)
  END IF
  wtminfull(mm) = wtt
  
!  Check to see if there are any logK = 500 in advance of call to subroutine fit
  
!        do i = 1,ntemp
!          if (ntemp.gt.1 .and. alogk0(i).eq.500.0) then
!            iflgint = 1
!            write(iunit2,*)
!            write(*,*)
!     write(iunit2,*) ' LogK = 500 encountered for species '
!     &    ,   namsec(nreac)
!     write(*,*) ' LogK = 500 encountered for species '
!     &    ,   namsec(nreac)
!            write(*,*) ' Non-isothermal problem--logKs missing'
!            write(iunit2,*) ' Non-isothermal problem--logKs missing'
!            write(*,*) ' Aborting run'
!            write(iunit2,*) ' Aborting run'
!            write(iunit2,*)
!            write(*,*)
!            stop
!          endif
!        end do
  
  IF (ntemp > 1) THEN
      
    IF (RunIsothermal) THEN
      
!!  Find the right logK value in the database for the specified temperature
      TPointer = 0
      DO i = 1,ntemp
        IF (tempc(i) == tinit) THEN
            TPointer = i
            alogk(nreac) = alogk0(TPointer)
        END IF
      END DO
      IF (Tpointer == 0) THEN
        WRITE(*,*)
        WRITE(*,*) ' Single temperature value for isothermal run not found among temperature points in database'
        WRITE(*,*)
        READ(*,*)
        STOP
     END IF
    
  ELSE 

    NameTransfer = namsec(nmat)
    CALL fit(nbasis,ntemp,alogk0,bvec,vec,iflgint,INT,ntt,NameTransfer)
    
!    IF (iflgint == 1) THEN
!      WRITE(iunit2,*)
!      WRITE(*,*)
!      WRITE(iunit2,*) 'WARNING! LogK = 500 encountered for species '  &
!          ,namsec(nmat)
!      WRITE(*,*) 'WARNING! LogK = 500 encountered for species ' ,namsec(nmat)
!    END IF
    
    DO  j = 1, nbasis
      coef0(nmat,j) = bvec(j)
    END DO
    
    alogk(nmat) = flogk(bvec,temp)
    
  END IF
  
  ELSE IF (ntemp == 1) THEN
    alogk(nmat) = alogk0(1)
  END IF
  
  a(nmat,nmat) = 1.
  
!     ----rearrange reactions - check for primary species
DO  i = 1, n
    
    speciesfound = .FALSE.

    DO  j = 1, ncomp
      IF (nam(i+1) == namc(j)) THEN
        b(nmat,j) = sto(i+1)
        speciesfound = .TRUE.
        EXIT
      END IF
    END DO

    IF (speciesfound) THEN
      CYCLE
    END IF
    
!     ----check for secondary species
    
    DO  k = 1, ncmplx
      IF (nam(i+1) == namcx(k)) THEN
        a(nmat,k) = -sto(i+1)
        speciesfound = .TRUE.
        EXIT
      END IF
    END DO

    IF (speciesfound) THEN
      CYCLE
    END IF
    
    DO  ngs = 1, ngas
      IF (nam(i+1) == namg(ngs)) THEN
!       a(nmat,ncmplx+mnrl+ngs) = -sto(i+1)
        a(nmat,ncmplx+ngs) = -sto(i+1)
        speciesfound = .TRUE.
        EXIT
      END IF
    END DO

    IF (speciesfound) THEN
      CYCLE
    END IF

  END DO 
  
  GO TO 31
  
  32 CONTINUE
  
  IF (mnrl_new /= mnrl) THEN
    WRITE(*,*)
    WRITE(*,*) ' Minerals added from the EQ3 database'
    WRITE(*,*)
    REWIND iunit5
    READ(iunit5,*) dummy
    5005   READ(iunit5,*) nam(1)
    IF (nam(1) == 'End of gases') THEN
      GO TO 5006
    END IF
    GO TO 5005
    5006   icheck = mnrl_new - mnrl
    WRITE(*,*)
    WRITE(*,*) ' Number of new minerals added = ',icheck
    WRITE(*,*)
    WRITE(*,*) ' *************************************'
    WRITE(*,*) ' End of first sweep through minerals'
    WRITE(*,*)
    WRITE(*,*) ' nreac = ',nreac
    WRITE(*,*) ' ncmplx = ',ncmplx
    WRITE(*,*) ' mnrl = ',mnrl
    WRITE(*,*) ' ngas = ',ngas
    WRITE(*,*) ' *************************************'
    WRITE(*,*)
    mnrl = mnrl_new
    nsec = ncmplx + ngas + mnrl
    nreac = ncmplx + ngas
    nmat0 = nreac
    isweep = 1
    GO TO 31
  END IF
  
  WRITE(*,*)
  WRITE(*,*) ' *************************************'
  IF (icomplete == 1) THEN
    WRITE(*,*) ' End of second sweep through minerals'
  ELSE
    WRITE(*,*) ' Finished reading database'
  END IF
  WRITE(*,*)
  WRITE(*,*) ' Number of reactions = ',nreac
  WRITE(*,*) ' Number of complexes = ',ncmplx
  WRITE(*,*) ' Number of minerals  = ',mnrl
  WRITE(*,*) ' Number of gases     = ',ngas
  
  IF (nreac /= ncmplx+ngas+mnrl) iflgck = 1
  IF (iflgck == 1) THEN
    WRITE(*,*) ' Could not find a mineral reaction'
    DO k = 1,mnrl
      DO kk = 1,nreacmin
        IF (min_found(kk) == namrl(k)) THEN
          GO TO 5111
        END IF
      END DO
      WRITE(*,*)
      WRITE(*,*) ' Mineral reaction not found in thermo database'
      WRITE(*,*) namrl(k)
      WRITE(*,*)
      5111     CONTINUE
    END DO
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
  
ELSE    ! Case where icomplete .ne. 1 (NO DATABASE SWEEP)
  
  DO k = 1,mnrl
    REWIND iunit5
    CALL mineralfind(iunit5)
    namtemp = namrl(k)
    CALL stringlen(namtemp,lengthmin)
    1111     READ(iunit5,*) namdummy                               !! Loop back to here to read another database entry
    IF (namdummy == 'End of minerals') THEN
      WRITE(*,*)
      WRITE(*,*) ' Mineral not found in database'
      WRITE(*,*) ' Mineral = ',namtemp(1:lengthmin)
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
    IF (namdummy == namtemp) THEN   ! Mineral found
      
!  Check to see that the mineral can be described with the
!    available primary/secondary species, or gases
      
      BACKSPACE iunit5
      
      READ(iunit5,*,ERR=6004) nam(1),vbar0,n,(sto(i+1),nam(i+1),i=1,n),  &
          (alogk0(l),l=1,ntemp),wtt
      
      DO  i = n+1,1,-1
        
        speciesfound = .FALSE.
        namdum = nam(i)
        CALL stringlen(namdum,lenspec)
        
!       check for water
        
        IF (nam(i) == 'h2o' .OR. nam(i) == 'H2O') THEN
          nct = nct + 1
          speciesfound = .TRUE.
        END IF

        IF (speciesfound) THEN
          CYCLE
        END IF
        
        DO m = 1, mnrl
          IF (nam(i) == namrl(m)) THEN
            nct = nct + 1
            indxmin(m) = 1
            mm = m
            speciesfound = .TRUE.
            EXIT
          END IF
        END DO

        IF (speciesfound) THEN
          CYCLE
        END IF
        
        DO j = 1, ncomp
          IF (nam(i) == namc(j)) THEN
            nct = nct + 1
            indxpri(j) = 1
            speciesfound = .TRUE.
            EXIT
          END IF
        END DO

        IF (speciesfound) THEN
          CYCLE
        END IF
        
        DO ix = 1, ncmplx
          IF (nam(i) == namcx(ix)) THEN
            nct = nct + 1
            indxsec(ix) = 1
            speciesfound = .TRUE.
            EXIT
          END IF
        END DO

        IF (speciesfound) THEN
          CYCLE
        END IF
        
        DO ngs = 1, ngas
          IF (nam(i) == namg(ngs)) THEN
            nct = nct + 1
            indxgas(ngs) = 1
            speciesfound = .TRUE.
            EXIT
          END IF
        END DO

        IF (speciesfound) THEN
          CYCLE
        END IF
        
!     ----species or gas in mineral reaction not found--ABORT
        
        WRITE(*,*)
        WRITE(*,*) ' Species not found in mineral reaction'
        WRITE(*,*) ' Species = ',namdum(1:lenspec)
        WRITE(*,*) ' In mineral: ',namtemp(1:lengthmin)
        WRITE(*,*)
        READ(*,*)
        STOP

      END DO
    ELSE                                            !! Keep reading database if mineral cannot be loaded
      GO TO 1111
    END IF
    
!  Mineral found, add to reaction list
    
    nreac = nreac + 1
    nreacmin = nreacmin + 1
    min_found(nreacmin) = nam(1)
        
!          nmat = nmat0 + mm
    nmat = nmat0 + k
    namsec(nmat) = nam(1)
    if (mintype(nmat-nmat0) == 0) then
      vbar(nmat-nmat0) = vbar0*1.d-6  ! Convert to m**3/mole
    else if (mintype(nmat-nmat0) == 1) then
      vbar(nmat-nmat0) = 1.0d0 ! mole biomass/mole biomass  ! obsolete: vbar0 cells/mole
    end if
    wtminfull(nmat-nmat0) = wtt
    
    
    IF (ntemp > 1) THEN
        
  IF (RunIsothermal) THEN
      
!!  Find the right logK value in the database for the specified temperature
    TPointer = 0
    DO i = 1,ntemp
        IF (tempc(i) == tinit) THEN
            TPointer = i
            alogk(nreac) = alogk0(TPointer)
        END IF
    END DO
    IF (Tpointer == 0) THEN
       WRITE(*,*)
       WRITE(*,*) ' Single temperature value for isothermal run not found among temperature points in database'
       WRITE(*,*)
       READ(*,*)
       STOP
    END IF
    
  ELSE 
      
      NameTransfer = namsec(nmat)
      CALL fit(nbasis,ntemp,alogk0,bvec,vec,iflgint,INT,ntt,NameTransfer)
      
      DO j = 1, nbasis
        coef0(nmat,j) = bvec(j)
      END DO
      
      alogk(nmat) = flogk(bvec,temp)
      
    END IF   !! End of RunIsothermal block
  
    ELSE IF (ntemp == 1) THEN
      alogk(nmat) = alogk0(1)
    END IF
    
    a(nmat,nmat) = 1.0d0
    
!     ----rearrange reactions - check for primary species
    
    DO  i = 1, n
      
      speciesfound = .FALSE.

      DO j = 1, ncomp
        IF (nam(i+1) == namc(j)) THEN
          b(nmat,j) = sto(i+1)
          speciesfound = .TRUE.
          EXIT
        END IF
      END DO

      IF (speciesfound) THEN
        CYCLE
      END IF
      
!     ----check for secondary species
      
      DO kk = 1, ncmplx
        IF (nam(i+1) == namcx(kk)) THEN
          a(nmat,kk) = -sto(i+1)
          speciesfound = .TRUE.
          EXIT
        END IF
      END DO

      IF (speciesfound) THEN
        CYCLE
      END IF
      
      DO ngs = 1, ngas
        IF (nam(i+1) == namg(ngs)) THEN
          a(nmat,ncmplx+ngs) = -sto(i+1)
          speciesfound = .TRUE.
          EXIT
        END IF
      END DO

      IF (speciesfound) THEN
        CYCLE
      END IF
      
    END DO
    
    
  END DO   !  End of mineral loop
  
!!  WRITE(*,*) ' END OF MINERAL READ'
  WRITE(*,*)
  WRITE(*,*) ' Number of reactions:                   ',nreac
  WRITE(*,*) ' Number of complexes:                   ',ncmplx
  WRITE(*,*) ' Number of minerals:                    ',mnrl
  WRITE(*,*) ' Number of gases:                       ',ngas

END IF

!  Check to see that all the aqueous species have been found

IF (icomplete /= 1) THEN
  icantfind = 0
  DO k = 1,ncmplx
    DO kk = 1,nreacaq+nprimary
      IF (aq_found(kk) == namcx(k) ) THEN
        GO TO 6113
      END IF
    END DO
    namtemp = namcx(k)
    CALL stringlen(namtemp,ls)
    WRITE(*,*) ' Secondary species in input file not found in database: ',namtemp(1:ls)
    icantfind = 1
    6113     CONTINUE
  END DO
  DO k = 1,ncomp
    DO kk = 1,nreacaq+nprimary
      IF (aq_found(kk) == namc(k) ) THEN
        GO TO 6114
      END IF
    END DO
    namtemp = namc(k)
    CALL stringlen(namtemp,ls)
    WRITE(*,*) ' Primary species in input file not found in database: ',namtemp(1:ls)
    icantfind = 1
    6114     CONTINUE
  END DO
ELSE
  DO k = 1,ncomp
    DO kk = 1,nreacaq+nprimary
      IF (aq_found(kk) == namc(k) ) THEN
        GO TO 7114
      END IF
    END DO
    namtemp = namc(k)
    CALL stringlen(namtemp,ls)
    WRITE(*,*) ' Primary species in input file not found in database: ',namtemp(1:ls)
    icantfind = 1
    7114     CONTINUE
  END DO
END IF

IF (icantfind == 1) THEN
  WRITE(*,*)
  WRITE(*,*) ' Aborting run'
  WRITE(*,*)
  READ(*,*)
  STOP
END IF

!****************************************************************
!                       ***SURFACE COMPLEXATION****

IF (nsurf == 0) GO TO 5555

isweep_surf = 0

65 REWIND (iunit5)

dumstring = 'Begin surface complexation'
CALL find_string(iunit5,dumstring,ifind)
IF (ifind == 0) THEN
  WRITE(*,*)
  WRITE(*,*) ' Beginning of "Surface Complexation" section not found'
  WRITE(*,*)
  READ(*,*)
  STOP
END IF

61 CONTINUE

READ(iunit5,'(a)') dumstring

IF (dumstring == 'End of surface complexation') GO TO 62

BACKSPACE (iunit5)

READ(iunit5,*,ERR=6006) nam(1),n,(sto(i+1),nam(i+1),i=1,n),  &
    (alogk0(l),l=1,ntemp)

sto(1) = -1.

!     ----decide if reaction occurs in chosen system

imiss = 0
DO  i = n+1,1,-1
  
  speciesfound = .FALSE.

!       check for water
  IF (nam(i) == 'h2o' .OR. nam(i) == 'H2O') THEN
    speciesfound = .TRUE.
  END IF

  IF (speciesfound) THEN
    CYCLE
  ENDIF
  
  DO ngs = 1, ngas
    IF (nam(i) == namg(ngs)) THEN
      speciesfound = .TRUE.
      EXIT
    END IF
  END DO

  IF (speciesfound) THEN
    CYCLE
  ENDIF
  
  DO j = 1, ncomp
    IF (nam(i) == namc(j)) THEN
      speciesfound = .TRUE.
      EXIT
    END IF
  END DO

  IF (speciesfound) THEN
    CYCLE
  ENDIF
  
  DO ix = 1, ncmplx
    IF (nam(i) == namcx(ix)) THEN
      ixx = ix
      speciesfound = .TRUE.
      EXIT
    END IF
  END DO

  IF (speciesfound) THEN
    CYCLE
  ENDIF
  
  DO is = 1, nsurf
    IF (nam(i) == namsurf(is)) THEN
      speciesfound = .TRUE.
      EXIT
    END IF
  END DO

  IF (speciesfound) THEN
    CYCLE
  ENDIF
  
  DO ns = 1,nsurf_sec
    IF (nam(i) == namsurf_sec(ns)) THEN
      speciesfound = .TRUE.
      EXIT
    END IF
  END DO

  IF (speciesfound) THEN
    CYCLE
  ENDIF
  
  imiss = imiss + 1
  namsave = nam(i)
  
END DO


IF (imiss == 1) THEN
  
!        if (ntemp.gt.1) then
!          do k = 1,ntemp
!            if (alogk0(k).eq.500.0) then
!              goto 61
!            endif
!          end do
!        endif
  
  nsurf_sec = nsurf_sec + 1
  IF (nsurf_sec > msurf_sec) THEN
    WRITE(*,*)
    WRITE(*,*) ' MSURF_SEC not dimensioned large enough in "params.inc"'
    WRITE(*,*) ' Number of secondary surface complexes = ',nsurf_sec
    WRITE(*,*) ' MSURF_SEC dimensioned at: ',msurf_sec
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
  namsurf_sec(nsurf_sec) = namsave
ELSE
  GO TO 61
END IF

nreac = nreac + 1
nreacsurf = nreacsurf + 1

namsec(nreac) = nam(1)

IF (ntemp > 1) THEN
    
      IF (RunIsothermal) THEN
      
!!  Find the right logK value in the database for the specified temperature
    TPointer = 0
    DO i = 1,ntemp
        IF (tempc(i) == tinit) THEN
            TPointer = i
            alogk(nreac) = alogk0(TPointer)
        END IF
    END DO
    IF (Tpointer == 0) THEN
       WRITE(*,*)
       WRITE(*,*) ' Single temperature value for isothermal run not found among temperature points in database'
       WRITE(*,*)
       READ(*,*)
       STOP
    END IF
    
  ELSE 
  NameTransfer = namsec(nreac)
  CALL fit(nbasis,ntemp,alogk0,bvec,vec,iflgint,INT,ntt,NameTransfer)
  
  DO j = 1, nbasis
    coef0(nreac,j) = bvec(j)
  END DO
  
  alogk(nreac) = flogk(bvec,temp)
  
  END IF   !! End of RunIsothermal block
  
ELSE IF (ntemp == 1) THEN
  alogk(nreac) = alogk0(1)
END IF

DO  i = 1, n+1
  
  speciesfound = .FALSE.

!     ----rearrange reactions - check for primary species
  DO j = 1, ncomp
    IF (nam(i) == namc(j)) THEN
      b(nreac,j) = sto(i)
      speciesfound = .TRUE.
      EXIT
    END IF
  END DO

  IF (speciesfound) THEN
    CYCLE
  END IF
  
  DO is = 1,nsurf
    IF (nam(i) == namsurf(is)) THEN
      b(nreac,is+ncomp) = sto(i)
      speciesfound = .TRUE.
      EXIT
    END IF
  END DO

  IF (speciesfound) THEN
    CYCLE
  END IF
  
!     ----check for secondary species
  
  DO k = 1, ncmplx
    IF (nam(i) == namcx(k)) THEN
      a(nreac,k) = -sto(i)
      speciesfound = .TRUE.
      EXIT
    END IF
  END DO

  IF (speciesfound) THEN
    CYCLE
  END IF
  
  DO ngs = 1, ngas
    IF (nam(i) == namg(ngs)) THEN
      a(nreac,ncmplx+ngs) = -sto(i)
      speciesfound = .TRUE.
      EXIT
    END IF
  END DO

  IF (speciesfound) THEN
    CYCLE
  END IF
  
  DO ns = 1, nsurf_sec
    IF (nam(i) == namsurf_sec(ns)) THEN
      a(nreac,ncmplx+ngas+mnrl+ns) = -sto(i)
      speciesfound = .TRUE.
      EXIT
    END IF
  END DO

  IF (speciesfound) THEN
    CYCLE
  END IF
  
END DO

GO TO 61

62 CONTINUE

nsec = ncmplx + ngas + mnrl + nsurf_sec

IF (isweep_surf == 0) THEN
WRITE(*,*)
WRITE(*,*) ' END OF FIRST SURFACE COMPLEX READ'
WRITE(*,*)
WRITE(*,*) ' Number of surface complexes:           ',nsurf_sec
isweep_surf = 1
GO TO 65
ELSE
WRITE(*,*)
WRITE(*,*) ' END OF SECOND SURFACE COMPLEX READ'
WRITE(*,*)
WRITE(*,*) ' Number of surface complexes:           ',nsurf_sec
END IF

5555 CONTINUE

!!  **********  END OF SURFACE COMPLEXATION READ  ********************************



100 FORMAT(' ',a20,2(1PE12.4))


IF (nreac /= nsec) THEN
  WRITE(iunit2,101) nreac,nsec
  101 FORMAT(' nreac =',i4,' nsec = ',i4)
  nn0 = nsec
  IF (nreac > nsec) nn0 = nreac
  DO  i = 1, nn0
    WRITE(iunit2,1114) i,namsec(i)
  END DO
  WRITE(*,*) ' nreac = ',nreac
  WRITE(*,*) ' nsec =', nsec
  GO TO 555
END IF

1114 FORMAT(' ',i3,2X,a20)

IF (iflgck == 1) GO TO 555

!     ----compute inverse matrix to [a]

lwork = 3*nsec
ALLOCATE(w(nsec,nsec))
ALLOCATE(w_orig(nsec,nsec))
ALLOCATE(indx(nsec))
ALLOCATE(ainv(nsec,nsec))
ALLOCATE(y(nsec))
ALLOCATE(work(lwork))
ALLOCATE(iwork(lwork))

ainv = 0.0
y = 0.0

DO  k = 1, nsec
  DO  l = 1, nsec
    w(l,k) = a(l,k)
  END DO
END DO

w_orig = w

!       if (iprint .eq. 2) then
DO  k = 1, nsec
!          write(iunit2,9999) namsec(k),(w(k,l),l=1,nsec)
  9999     FORMAT(' ',a8,(' ',25F6.2))
END DO
DO  k = 1, nsec
!          write(iunit2,9999) namsec(k),(b(k,j),j=1,ncomp)
END DO
! endif

IF (nsec > 0) THEN
  CALL dgetrf(nsec,nsec,w,nsec,indx,info)
!!  CALL ludcmp90(w,indx,dd,nsec)
!!  CALL ludcmp(w,nsec,n0,indx,dd)

  DO  i = 1, nsec
    DO  k = 1, nsec
      y(k) = 0.
      IF (k == i) y(k) = 1.
      bb(k) = y(k)
    END DO
  
!!  CALL lubksb(w,nsec,n0,indx,y)
!!  CALL lubksb90(w,indx,y,nsec)

    CALL dgetrs(trans,nsec,ione,w,nsec,indx,y,nsec,info)

    CALL dgerfs(trans,nsec,ione,w_orig,nsec,w,nsec,indx,bb,nsec,y,nsec,ferr,berr,work,iwork,info)
    CALL dgerfs(trans,nsec,ione,w_orig,nsec,w,nsec,indx,bb,nsec,y,nsec,ferr,berr,work,iwork,info) 
    CALL dgerfs(trans,nsec,ione,w_orig,nsec,w,nsec,indx,bb,nsec,y,nsec,ferr,berr,work,iwork,info)

    DO  l = 1, nsec
      ainv(l,i) = y(l)
    END DO
  END DO

!     write(iunit2,106)
!     do 210 i = 1, nsec
!       write(iunit2,105) i,(ainv(i,l), l = 1, nsec)
!210  continue

  105  FORMAT(' ',i3,10(1PE13.6))
  106  FORMAT(/,' transformation matrix')

  DEALLOCATE(work)

END IF

!     ----check inverse
!!iflg = 0
DO  i = 1, nsec
  DO  j = 1, nsec
    sum = 0.
    DO  k = 1, nsec
      sum = sum + ainv(i,k)*a(k,j)
    END DO
    IF (i == j) THEN
      IF (sum < 1.d0-eps .OR. sum > 1.d0+eps) THEN
        WRITE(*,*) 'identity check failed',i,j,sum
        iflg = 1
        stop
      END IF
    ELSE IF (i /= j) THEN
      IF (ABS(sum) > eps) THEN
        WRITE(*,*) 'identity check failed',i,j,sum
        iflg = 1
        stop
      END IF
    END IF
    smat(j,i) = sum
  END DO
END DO

IF (iflg == 1) THEN
  WRITE(iunit2,107)
  DO  i = 1, nsec
    WRITE(iunit2,105) i,(smat(j,i), j = 1, nsec)
  END DO
  107    FORMAT(/,' identity check')
  
! stop
END IF

!     ----compute new reaction coefficients and logK's
DO  i = 1, nsec
  DO  j = 1, ncomp+nsurf
    sum = 0.
    DO  k = 1, nsec
      sum = sum + ainv(i,k)*b(k,j)
    END DO
    smat(j,i) = sum
  END DO
  alogsum = 0.
  DO  ii = 1, nsec
    alogsum = alogsum + ainv(i,ii)*alogk(ii)
    IF (ntemp > 1) THEN
      DO  j = 1, nbasis
        coef(i,j) = coef(i,j) + ainv(i,ii)*coef0(ii,j)
      END DO
    END IF
  END DO
  alogkp(i) = alogsum
END DO

!************************************************************

!     ----construct submatrices shom, smnrl, sol.

DO  j = 1, ncomp
  DO  i = 1, ncmplx
    shom(j,i) = smat(j,i)
  END DO
  DO  m = 1, mnrl
    mm = ncmplx + m + ngas
    smnrl(j,m) = smat(j,mm)
  END DO
  DO  ngs = 1, ngas
    ngss = ncmplx + ngs
    stogas(j,ngs) = smat(j,ngss)
  END DO
END DO

DO i = 1, ncomp+nsurf
  DO ns = 1,nsurf_sec
    nss = ns + ncmplx + ngas + mnrl
    xsurf(i,ns) = smat(i,nss)
  END DO
END DO

DO ik = 1,ncmplx
  DO i = 1,ncomp
    muaq(ik,i) = shom(i,ik)
    IF (ABS(muaq(ik,i)) < 1.E-35) THEN
      muaq(ik,i) = 0.0
    END IF
  END DO
END DO

IF (icomplete == 1) THEN
  DO k = 1,mnrl
    DO i = 1,ncomp
      mumin(1,k,i) = smnrl(i,k)
      IF (ABS(mumin(1,k,i)) < 1.E-35) THEN
        mumin(1,k,i) = 0.0
      END IF
    END DO
    nreactmin(k) = 1
  END DO
ELSE
  msub = 0
  DO k = 1,nrct
    DO np = 1,nreactmin(k)
      msub = msub + 1
      DO i = 1,ncomp
        mumin(np,k,i) = smnrl(i,msub)
        IF (ABS(mumin(np,k,i)) < 1.E-35) THEN
          mumin(np,k,i) = 0.0
        END IF
      END DO
!             write(*,4444) namrl(msub),(mumin(np,k,i),i=1,ncomp)
    END DO
  END DO
  
END IF

DO ngs = 1,ngas
  DO i = 1,ncomp
    mugas(ngs,i) = stogas(i,ngs)
    IF (ABS(mugas(ngs,i)) < 1.E-35) THEN
      mugas(ngs,i) = 0.0
    END IF
  END DO
END DO

DO ns = 1,nsurf_sec
  DO i = 1,ncomp+nsurf
    musurf(ns,i) = xsurf(i,ns)
    IF (ABS(musurf(ns,i)) < 1.E-35) THEN
      musurf(ns,i) = 0.0
    END IF
  END DO
END DO

DO i = 1,ncomp
  wtaq(i) = wtbas(i)
END DO

!     ----store logK's
DO  i = 1, ncmplx
  ik = ncomp + i
  eqhom(i) = alogkp(i)
  wtaq(ik) = wtcmplx(i)
END DO
DO  m = 1, mnrl
! mm = ncmplx + m
  mm = ncmplx + m + ngas
  alnk(m) = alogkp(mm)
END DO
DO  ngs = 1, ngas
! ngss = ncmplx + mnrl + ngs
  ngss = ncmplx + ngs
  eqgas(ngs) = alogkp(ngss)
END DO
DO ns = 1, nsurf_sec
  ngss = ncmplx + ngas + mnrl + ns
  eqsurf(ns) = alogkp(ngss)
END DO

!  Fix up order of equilibrium constants so that it is back to
!  ncmplx, then mnrl, then ngas

DO m = 1,mnrl
  mm = m + ncmplx
  alogkp(mm) = alnk(m)
END DO
DO ngs = 1,ngas
  ngss = ncmplx + mnrl + ngs
  alogkp(ngss) = eqgas(ngs)
END DO

!      do ik = 1,ncmplx+ngas+mnrl
!        write(*,*) alogkp(ik)
!      end do
!      pause

!     write out temperature interpolation coefficients

IF (ntemp > 1) THEN
  WRITE(iunit2,*) ' temperature interpolation coefficients'
  WRITE(iunit2,7776)
  7776 FORMAT(' ',15X,'  ln(T+Tk)        1         T+Tk      (T+Tk)**-1   (t+tk)**-2')
!      if (iprint.eq.2) then
!      do 777 i = 1, nsec
!   write(iunit2,7777) namsec(i),(coef(i,j),j=1,nbasis)
!  777 continue
!      end if
  7777 FORMAT(' ',a20,5(1PE12.4))
END IF

!  599 format(' ',27x,<ncomp>a7)
!  600 format(' ',a12,1pe12.4,<ncomp>(0pf7.2))
!  601 format(' ',a12,1pe12.4,<ncomp+nsurf>(0pf7.2))
!  602 format(' ',27x,<ncomp+nsurf>a7)
599 FORMAT(' ',27X,100A7)
600 FORMAT(' ',a12,1PE12.4,100(0PF7.2))
601 FORMAT(' ',a12,1PE12.4,100(0PF7.2))
602 FORMAT(' ',27X,100A7)

IF (ncmplx > 0) THEN
  IF (iprint == 1) THEN
    WRITE(iunit2,*)
    DO i = 1,ncomp
      namtemp = namc(i)
      namdummy = namtemp(1:6)
      namdummy(7:7) = ' '
      namprnt(i) = namdummy(1:7)
    END DO
    DO is = 1,nsurf
      namtemp = namsurf(is)
      namdummy = namtemp(1:6)
      namdummy(7:7) = ' '
      namprnt(is+ncomp) = namdummy(1:7)
    END DO
    WRITE(iunit2,*) '   HOMOGENEOUS REACTIONS'
    WRITE(iunit2,*)
    WRITE(iunit2,*) 'Species         log K         Stoichiometric Coefficients'
    WRITE(iunit2,599) (namprnt(j),j=1,ncomp)
    DO  i = 1, ncmplx
      WRITE(iunit2,600) namcx(i),eqhom(i),(shom(j,i),j=1,ncomp)
!!  Check electrical balance on reactions
      charge = 0.0
      DO j = 1,ncomp
        charge = charge + shom(j,i)*z(j)
      END DO
      charge = charge - zx(i)
      IF (ABS(charge) > tiny) THEN
        WRITE(*,*)  
        WRITE(*,*) ' Homogeneous reaction not charge balanced'
        WRITE(*,603) namcx(i), charge
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
 603 FORMAT(1x,a15,1x,f10.2)
    END DO
  ELSE
    WRITE(iunit2,*)
    WRITE(iunit2,*) '   Homogeneous Reactions'
    WRITE(iunit2,*) 'Species         log K '
    DO i = 1, ncmplx
      WRITE(iunit2,600) namcx(i),eqhom(i)
    END DO
  END IF
END IF

IF (mnrl > 0) THEN
  IF (iprint == 1) THEN
    WRITE(iunit2,*)
    WRITE(iunit2,*) '     MINERAL REACTIONS'
    WRITE(iunit2,*)
    WRITE(iunit2,*) 'Mineral         log K         Stoichiometric Coefficients'
    WRITE(iunit2,599) (namprnt(j),j=1,ncomp)
    iflgstop = 0
    DO  m = 1, mnrl
      WRITE(iunit2,600) namrl(m),alnk(m),(smnrl(j,m),j=1,ncomp)

!!  Check electrical balance on reactions
      charge = 0.0
      DO j = 1,ncomp
        charge = charge + smnrl(j,m)*z(j)
      END DO
      IF (ABS(charge) > 1.E-10 .AND. namrl(m) /= 'H(c)') THEN
        WRITE(*,*)  
        WRITE(*,*) ' Mineral reaction not charge balanced'
        WRITE(*,604) namrl(m), charge
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
 604 FORMAT(1x,a15,1x,1pe12.2)

      IF (vbar(m) == 0.) THEN
        iflgstop = 1
        dumstring = namrl(m)
        CALL stringlen(dumstring,ls)
        WRITE(*,*) '     WARNING  '
        WRITE(*,*) ' Zero molar volume for mineral ',dumstring(1:ls)
        WRITE(*,*) '   ABORTING RUN'
        WRITE(*,*)
        WRITE(iunit2,*) '     WARNING  '
        WRITE(iunit2,*) ' Zero molar volume for mineral ',dumstring(1:ls)
        WRITE(iunit2,*) '   ABORTING RUN'
        WRITE(iunit2,*)
!        STOP
      END IF
    END DO
  ELSE
    WRITE(iunit2,*)
    WRITE(iunit2,*) '     MINERAL REACTIONS'
    WRITE(iunit2,*)
    WRITE(iunit2,*) 'Mineral        log K '
    iflgstop = 0
    DO m = 1, mnrl
      WRITE(iunit2,600) namrl(m),alnk(m)
      IF (vbar(m) == 0.) THEN
        iflgstop = 1
        dumstring = namrl(m)
        CALL stringlen(dumstring,ls)
        WRITE(*,*) '     WARNING  '
        WRITE(*,*) ' Zero molar volume for mineral ',dumstring(1:ls)
        WRITE(*,*) '   ABORTING RUN'
        WRITE(*,*)
        WRITE(iunit2,*) '     WARNING  '
        WRITE(iunit2,*) ' Zero molar volume for mineral ',dumstring(1:ls)
        WRITE(iunit2,*) '   ABORTING RUN'
        WRITE(iunit2,*)
        READ(*,*)
        STOP
      END IF
    END DO
  END IF
END IF

IF (ngas > 0) THEN
  IF (iprint == 1) THEN
    WRITE(iunit2,*)
    WRITE(iunit2,*) '    GAS REACTIONS'
    WRITE(iunit2,*)
    WRITE(iunit2,*) 'Gases           log K         Stoichiometric Coefficients'
    WRITE(iunit2,599) (namprnt(j),j=1,ncomp)
    DO  i = 1, ngas
      WRITE(iunit2,600) namg(i),eqgas(i),(stogas(j,i),j=1,ncomp)
    END DO
  ELSE
    WRITE(iunit2,*)
    WRITE(iunit2,*) 'Gas Reactions'
    WRITE(iunit2,*) 'Gas            log K'
    DO i = 1, ngas
      WRITE(iunit2,600) namg(i),eqgas(i)
    END DO
  END IF
END IF


IF (nsurf_sec > 0) THEN
  WRITE(iunit2,*)
  WRITE(iunit2,*) '    SURFACE COMPLEXATION REACTIONS'
  WRITE(iunit2,*)
  WRITE(iunit2,*) 'Surface complexes log K       Stoichiometric Coefficients'
  WRITE(iunit2,602) (namprnt(j),j=1,ncomp+nsurf)
  DO is = 1, nsurf_sec
    WRITE(iunit2,601) namsurf_sec(is),eqsurf(is),(xsurf(j,is),j=1,ncomp+nsurf)
  END DO
END IF

!     close data file
CLOSE(iunit5,STATUS='keep')

DEALLOCATE(min_found)
DEALLOCATE(gas_found)
DEALLOCATE(aq_found)
DEALLOCATE(namprnt) 
DEALLOCATE(namsec)
DEALLOCATE(nam)
DEALLOCATE(smnrl)
DEALLOCATE(stogas)
DEALLOCATE(shom)
DEALLOCATE(xsurf)
DEALLOCATE(alogk0)
DEALLOCATE(tempc)
DEALLOCATE(vec)
DEALLOCATE(vecgam)
DEALLOCATE(bvec)
DEALLOCATE(sto)
DEALLOCATE(dummyreal)
DEALLOCATE(alogkeh)
DEALLOCATE(coefeh)
DEALLOCATE(alogk)
DEALLOCATE(a)
DEALLOCATE(b)
DEALLOCATE(smat)
DEALLOCATE(ainv)
DEALLOCATE(y)
DEALLOCATE(w)
DEALLOCATE(wgam)
DEALLOCATE(bb)
DEALLOCATE(coef0)
DEALLOCATE(wtcmplx)
DEALLOCATE(iload)
DEALLOCATE(indx)
DEALLOCATE(indxpri)
DEALLOCATE(indxsec)
DEALLOCATE(indxmin)
DEALLOCATE(indxgas)
DEALLOCATE(ndxkin)

RETURN

6001 WRITE(*,*) 'Error reading temperature data: stop'
READ(*,*)
STOP
7001 WRITE(*,*) 'Error reading Debye-Huckel parameters: stop'
READ(*,*)
STOP
6002 WRITE(*,*) 'Error reading primary species: stop'
WRITE(*,*) NAME
READ(*,*)
STOP
6003 WRITE(*,*) 'Error reading secondary species: stop'
WRITE(*,*) nam(1)
READ(*,*)
STOP
6004 WRITE(*,*) 'Error reading minerals: stop'
WRITE(*,*) nam(1)
READ(*,*)
STOP
6005 WRITE(*,*) 'Error reading gases: stop'
WRITE(*,*) nam(1)
READ(*,*)
STOP
6006 WRITE(*,*) 'Error reading surface complexes: stop'
WRITE(*,*) nam(1)
READ(*,*)
STOP
334 WRITE(*,*) 'Error in opening database file: stop'
READ(*,*)
STOP

555 CONTINUE

WRITE(iunit2,*) 'STOP -- number of reactions not equal to ',  &
    'number of secondary species!'
WRITE(*,*) 'STOP -- number of reactions not equal to ',  &
    'number of secondary species!'
DO  j = 1, ncomp
  IF (indxpri(j) /= 1) THEN
    WRITE(iunit2,*) 'species not found: ',namc(j)
    WRITE(*,*) 'species not found: ',namc(j)
  END IF
END DO
DO  j = 1, ncmplx
  IF (indxsec(j) /= 1) THEN
    WRITE(iunit2,*) 'species not found: ',namcx(j)
    WRITE(*,*) 'species not found: ',namcx(j)
  END IF
END DO
DO  j = 1, mnrl
  IF (indxmin(j) /= 1) THEN
    WRITE(iunit2,*) 'species not found: ',namrl(j)
    WRITE(*,*) 'species not found: ',namrl(j)
  END IF
END DO
DO  j = 1, ngas
  IF (indxgas(j) /= 1) THEN
    WRITE(iunit2,*) 'species not found: ',namg(j)
    WRITE(*,*) 'species not found: ',namg(j)
  END IF
END DO

4444 FORMAT(a15,15(1X,f6.2))
READ(*,*)
STOP
END SUBROUTINE database
!******************************************************************
FUNCTION flogk(b,t)
USE crunchtype
REAL(DP), INTENT(IN)                          :: t
REAL(DP)                                      :: temp
REAL(DP)                                      :: flogk
REAL(DP), DIMENSION(5)                        :: b

temp = t + 273.15
flogk = b(1)*LOG(temp) + b(2)  &
    + b(3)*temp + b(4)/temp  &
    + b(5)/(temp*temp)

RETURN
END FUNCTION flogk
!**************************************************************
