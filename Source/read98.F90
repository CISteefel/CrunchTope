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

SUBROUTINE read98(ncomp,ncmplx,nkin,nrct,ngas,nsurf,nsurf_sec,data1,  &
    icomplete,ispeciate,igenericrates,GMSsecondary,GMSmineral,RateGeneric)
USE crunchtype
USE params
USE concentration
USE mineral
USE temperature
USE temperature
USE io
USE strings

IMPLICIT NONE

INTERFACE
  SUBROUTINE read_minkin(nout,ncomp,nspec,nkin,ngas,data1,namrl,  &
    ispeciate,igenericrates,nammin_tmp)
  USE crunchtype
  USE CrunchFunctions
  USE params
  USE concentration
  USE mineral
  USE strings
  IMPLICIT NONE
  INTEGER(I4B), INTENT(IN)                                    :: nout
  INTEGER(I4B), INTENT(IN)                                    :: ncomp
  INTEGER(I4B), INTENT(IN)                                    :: nspec
  INTEGER(I4B), INTENT(OUT)                                   :: nkin
  INTEGER(I4B), INTENT(IN)                                    :: ngas
  INTEGER(I4B), INTENT(IN)                                    :: ispeciate
  INTEGER(I4B), INTENT(IN)                                    :: igenericrates
  CHARACTER (LEN=mls), DIMENSION(:), INTENT(IN OUT)           :: nammin_tmp
  CHARACTER (LEN=mls), INTENT(IN)                             :: data1
  CHARACTER (LEN=mls), DIMENSION(:), INTENT(IN OUT)           :: namrl
  END SUBROUTINE read_minkin
END INTERFACE
INTERFACE
  SUBROUTINE database(ncomp,ncmplx,mnrl,nrct,ngas,nsurf,nsurf_sec,ntemp,  &
    iprint,icomplete,jtemp,data1,namc,namcx,namrl,  &
    vbar,wtminfull,wtbas,coef,temp,z,bdotPrimary,zx,bdotSecondary,a0,ax0)
  USE crunchtype
  USE params
  USE concentration
  USE mineral
  USE io
  IMPLICIT NONE
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
  END SUBROUTINE database
END INTERFACE

!  External variables and arrays

INTEGER(I4B), INTENT(IN OUT)                                           :: ncomp
INTEGER(I4B), INTENT(IN OUT)                                           :: ncmplx
INTEGER(I4B), INTENT(IN OUT)                                           :: nkin
INTEGER(I4B), INTENT(IN OUT)                                           :: nrct
INTEGER(I4B), INTENT(IN OUT)                                           :: ngas
INTEGER(I4B), INTENT(IN OUT)                                           :: nsurf
INTEGER(I4B), INTENT(IN OUT)                                           :: nsurf_sec
INTEGER(I4B), INTENT(IN)                                               :: icomplete
INTEGER(I4B), INTENT(IN)                                               :: ispeciate
INTEGER(I4B), INTENT(IN)                                               :: igenericrates

CHARACTER (LEN=mls), INTENT(IN)                                        :: data1

LOGICAL(LGT), INTENT(IN)                                               :: GMSsecondary
LOGICAL(LGT), INTENT(IN)                                               :: GMSmineral

REAL(DP), INTENT(IN)                                                   :: RateGeneric

!  Internal variables and arrays

CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE                     :: namc
CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE                     :: namcx
CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE                     :: namrl
CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE                     :: nammin_tmp
CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE                     :: dumlabel

REAL(DP), DIMENSION(:), ALLOCATABLE                                :: vbar
REAL(DP), DIMENSION(:), ALLOCATABLE                                :: wtminfull
REAL(DP), DIMENSION(:), ALLOCATABLE                                :: wtbas
REAL(DP), DIMENSION(:), ALLOCATABLE                                :: a0
REAL(DP), DIMENSION(:), ALLOCATABLE                                :: z
REAL(DP), DIMENSION(:), ALLOCATABLE                                :: bdotPrimary

REAL(DP), DIMENSION(:), ALLOCATABLE                                :: ax0
REAL(DP), DIMENSION(:), ALLOCATABLE                                :: zx
REAL(DP), DIMENSION(:), ALLOCATABLE                                :: bdotSecondary
REAL(DP), DIMENSION(:), ALLOCATABLE                                :: vbarkin
!!!REAL(DP), DIMENSION(:), ALLOCATABLE                                :: azero
REAL(DP), DIMENSION(:,:), ALLOCATABLE                              :: coef

CHARACTER (LEN=mls)                                                :: parchar
CHARACTER (LEN=mls)                                                :: parfind
REAL(DP)                                                           :: realjunk
INTEGER(I4B)                                                       :: lchar

CHARACTER (LEN=50)                                                 :: label
CHARACTER (LEN=mls)                                                :: string1
CHARACTER (LEN=mls)                                                :: namtemp
CHARACTER (LEN=mls)                                                :: namtemp2
CHARACTER (LEN=mls)                                                :: section

LOGICAL(LGT)                                                       :: found
LOGICAL(LGT)                                                       :: O2gasFound
LOGICAL(LGT)                                                       :: H2gasFound

INTEGER(I4B)                                                       :: laffinity
INTEGER(I4B)                                                       :: iprint
INTEGER(I4B)                                                       :: nin
INTEGER(I4B)                                                       :: nout
INTEGER(I4B)                                                       :: i
INTEGER(I4B)                                                       :: ksp
INTEGER(I4B)                                                       :: ncmplx_old
INTEGER(I4B)                                                       :: ik
INTEGER(I4B)                                                       :: njunk
INTEGER(I4B)                                                       :: k
INTEGER(I4B)                                                       :: ls
INTEGER(I4B)                                                       :: np
INTEGER(I4B)                                                       :: kk
INTEGER(I4B)                                                       :: lmin
INTEGER(I4B)                                                       :: mnrl
INTEGER(I4B)                                                       :: mm
INTEGER(I4B)                                                       :: ltrim
INTEGER(I4B)                                                       :: jh2o
INTEGER(I4B)                                                       :: jph
INTEGER(I4B)                                                       :: joh
INTEGER(I4B)                                                       :: jo2
INTEGER(I4B)                                                       :: jco2
INTEGER(I4B)                                                       :: j
INTEGER(I4B)                                                       :: loop
INTEGER(I4B)                                                       :: ncount
INTEGER(I4B)                                                       :: ks
INTEGER(I4B)                                                       :: lsurf
INTEGER(I4B)                                                       :: lssec
INTEGER(I4B)                                                       :: lsmin
INTEGER(I4B)                                                       :: m
INTEGER(I4B)                                                       :: msub
INTEGER(I4B)                                                       :: ll
INTEGER(I4B)                                                       :: id
INTEGER(I4B)                                                       :: iff
INTEGER(I4B)                                                       :: ids
INTEGER(I4B)                                                       :: lzs
INTEGER(I4B)                                                       :: ifind
INTEGER(I4B)                                                       :: is
INTEGER(I4B)                                                       :: ns

REAL(DP)                                                           :: yrsec
REAL(DP)                                                           :: temp

LOGICAL(LGT)                                                       :: NameListFormat
LOGICAL(LGT)                                                       :: IonsOK

ALLOCATE(namc(nc))
ALLOCATE(namcx(mcmplx))
ALLOCATE(namrl(nm*mreact))
ALLOCATE(nammin_tmp(nm))
ALLOCATE(dumlabel(50))
ALLOCATE(vbar(nm))
ALLOCATE(wtminfull(nm))
ALLOCATE(wtbas(nc))
ALLOCATE(a0(nc))
ALLOCATE(z(nc))
ALLOCATE(bdotPrimary(nc))
ALLOCATE(ax0(mcmplx))
ALLOCATE(zx(mcmplx))
ALLOCATE(bdotSecondary(mcmplx))
ALLOCATE(vbarkin(nm))
IF (ALLOCATED(azero)) THEN
  DEALLOCATE(azero)
END IF
ALLOCATE(azero(nc+mcmplx))
ALLOCATE(coef(ndim,5))

!********************************
!  INPUT DATA
!********************************

ncmplx_old = 0
npointO2aq = 0
npointO2gas = 0
npointH2gas = 0
O2Found = .false.
O2gasFound = .false.
H2gasFound = .false.
yrsec = 1.0/secyr
temp = tinit
iprint = 1

nin = iunit1
nout = 4

!******************************************
!  Primary (component) species
!******************************************
section = 'primary_species'
CALL readblock(nin,nout,section,found,ncomp)

IF (found) THEN
!!  WRITE(*,*)
!!  WRITE(*,*) ' Primary species block found'
!!  WRITE(*,*)
  
ELSE
  
  WRITE(*,*)
  WRITE(*,*) ' No primary species given'
  WRITE(*,*) ' This is going to be a pretty boring simulation'
  WRITE(*,*) ' Lets bag it!'
  READ(*,*)
  STOP
  
END IF

REWIND nout

ikh2o = 0
DO i = 1,ncomp
  IF (ncomp > nc) THEN
    WRITE(*,*)
    WRITE(*,*) ' Primary species dimension too small! ',ncomp
    WRITE(*,*) ' Redimension NC in params.inc'
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
  IF (icomplete /= 1 .AND. ncomp > mcomp) THEN
    WRITE(*,*)
    WRITE(*,*) ' Primary species dimension too small! ',ncomp
    WRITE(*,*) ' Redimension MCOMP in params.inc'
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
  READ(nout,'(a)') namc(i)
  IF (namc(i) == 'h2o' .OR. namc(i) == 'H2O') ikh2o = i
  IF (namc(i) == 'end') THEN
    GO TO 33
  END IF
END DO

33  CONTINUE

WRITE(iunit2,*)
WRITE(*,*)
WRITE(*,*)      ' Number of components:               ',ncomp
WRITE(iunit2,*) 'Number of components =          ',ncomp
WRITE(iunit2,*)

DO i = 1,ncomp
!!  WRITE(*,1025) i,namc(i)
  WRITE(iunit2,1025) i,namc(i)
END DO

1025 FORMAT(i3,1X,a18)

DO i = 1,ncomp
  ulab(i) = namc(i)
END DO

!*********************************
!  Read reversible Secondary Species
!*********************************

! If a database sweep is specified, don't read the
! secondary species or the gases

IF (icomplete /= 1) THEN
  section = 'secondary_species'
  CALL readblock(nin,nout,section,found,ncmplx)
  
  IF (found) THEN
!!    WRITE(*,*)
!!    WRITE(*,*) ' Secondary species block found'
!!    WRITE(*,*)
  ELSE
    WRITE(*,*) ' No secondary species given'
    WRITE(*,*)
    ncmplx = 0
  END IF
  
  REWIND nout
  
  DO ksp = 1,ncmplx
    IF (ncmplx > mcmplx) THEN
      WRITE(*,*)
      WRITE(*,*) ' Second. species dimension too small! ',ncmplx
      WRITE(*,*) ' Redimension MCMPLX in params.inc'
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
    IF (ncmplx > mspec .AND. icomplete /= 1) THEN
      WRITE(*,*)
      WRITE(*,*) ' Second. species dimension too small! ',ncmplx
      WRITE(*,*) ' Redimension MSPEC in params.inc'
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
    READ(nout,'(a)') namcx(ksp)
    IF (namcx(ksp) == 'end') THEN
      GO TO 43
    END IF
  END DO
  
  43    CONTINUE
  
  ncmplx_old = ncmplx
  
  WRITE(iunit2,*)
  WRITE(*,*)
  WRITE(*,*)      ' Number of secondary species:        ', ncmplx
  WRITE(iunit2,*) 'Number of secondary species in input file = ', ncmplx
  WRITE(iunit2,*)
  
  DO ksp = 1,ncmplx
!!    WRITE(*,1025) ksp,namcx(ksp)
    WRITE(iunit2,1025) ksp,namcx(ksp)
  END DO
  
  DO ksp = 1,ncmplx
    ik = ksp + ncomp
    ulab(ik) = namcx(ksp)
  END DO
  
END IF

!  End of secondary species read

!********************************
!  Read gases
!********************************

!  If database sweep is specified, don't read the gases

IF (icomplete /= 1) THEN
  section = 'gases'
  CALL readblock(nin,nout,section,found,ngas)
  
  IF (found) THEN
!!    WRITE(*,*) ' Gases block found'
!!    WRITE(*,*)
    
    REWIND nout
    
    DO i = 1, ngas
      IF (i > ng) THEN
        WRITE(*,*)
        WRITE(*,*) '  NG dimensioned too small',ngas
        WRITE(*,*) '  Redimension NG in params.inc'
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
      IF (i > mgas .AND. icomplete /= 1) THEN
        WRITE(*,*)
        WRITE(*,*) '  MGAS dimensioned too small',ngas
        WRITE(*,*) '  Redimension MGAS in params.inc'
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
      READ(nout,'(a)',ERR=333) namg(i)
    END DO
    WRITE(iunit2,*)
    WRITE(*,*)
    WRITE(iunit2,*) 'Number of gases in input file =     ',ngas
    WRITE(*,*)      ' Number of gases:                    ',ngas
    WRITE(iunit2,*)
    DO i = 1,ngas
      WRITE(iunit2,1040) i,namg(i)
!!      WRITE(*,1040) i,namg(i)
    END DO
    WRITE(iunit2,*)
!!    WRITE(*,*)
    DO ik = 1,ncomp+ncmplx
      IF (ulab(ik) == 'O2(aq)') THEN
        O2Found = .true.
        npointO2aq = ik
      END IF
    END DO
    IF (O2Found) THEN
!     Check to see if O2(g) has been included--if not, add it
      DO i = 1,ngas
        IF (namg(i) == 'O2(g)') THEN
          O2gasFound = .true.
          npointO2gas = i
        END IF
        IF (namg(i) == 'H2(g)') THEN
          H2gasFound = .true.
          npointH2gas = i
        END IF
      END DO

      IF (.not. O2gasfound) THEN
        ngas = ngas + 1
        namg(ngas) = 'O2(g)'
        npointO2gas = ngas
      END IF
!!      IF (.not. H2gasfound) THEN
!!        ngas = ngas + 1
!!       namg(ngas) = 'H2(g)'
!!       npointH2gas = ngas
!!      END IF

    END IF

  ELSE
    ngas = 0
    DO ik = 1,ncomp+ncmplx
      IF (ulab(ik) == 'O2(aq)') THEN
        O2Found = .true.
        npointO2aq = ik
      END IF
    END DO
    IF (O2Found) THEN
      namg(1) = 'O2(g)'
!!      namg(2) = 'H2(g)'
!!      ngas = 2
      ngas = 1
      npointO2gas = 1
!!      npointH2gas = 2
    END IF
  END IF
  
END IF

!  *****End of gas read****

!********************************
!  Read kinetically reacting minerals
!********************************

!  Don't read minerals if database sweep is specified

IF (icomplete /= 1) THEN
  section = 'minerals'
  CALL readblock(nin,nout,section,found,njunk)
  
  IF (found) THEN
!!    WRITE(*,*) ' Minerals block found'
!!    WRITE(*,*)
    
    CALL read_minkin(nout,ncomp,ncmplx,nkin,ngas,data1,namrl,  &
        ispeciate,igenericrates,nammin_tmp)
    
    nreactmax = 0
    nmonodmax = 0
    ninhibitmax = 0
    DO k = 1,nkin
      umin(k) = nammin_tmp(k)
      nreactmax = MAX0(nreactmin(k),nreactmax)
       DO np = 1,nreactmin(k)
         nmonodmax = MAX0(nmonod(np,k),nmonodmax)
         ninhibitmax = MAX0(ninhibit(np,k),ninhibitmax)
       END DO
      
    END DO
    
    DO k = 1,nkin
      CALL stringlen(umin(k),ls)
      lenmin(k) = ls
    END DO
    
! Check for cross affinity terms--make sure minerals listed are present 
    
    IF (ispeciate /= 1 .AND. igenericrates /= 1) THEN
      DO k = 1,nkin
        DO np = 1,nreactmin(k)
          IF (crossaff(np,k) /= 'none') THEN
            DO kk = 1,nkin
              IF (crossaff(np,k) == umin(kk)) THEN
                kcrossaff(np,k) = kk
                GO TO 555
              END IF
            END DO
            namtemp = umin(k)
            CALL stringlen(namtemp,lmin)
            namtemp2 = crossaff(np,k)
            CALL stringlen(namtemp2,laffinity)
            
            WRITE(*,*)
            WRITE(*,*) ' Mineral for cross affinity not found in mineral list'
            WRITE(*,*) ' Phase: ', namtemp(1:lmin)
            WRITE(*,*) ' Cross-affinity phase not found: ',namtemp2(1:laffinity)
            WRITE(*,*)
            READ(*,*)
            STOP
            555          CONTINUE
          END IF
        END DO
      END DO
      
    END IF     !  Dont do this IF block if only speciating
    
    
    nrct = nkin
    mnrl = 0
    DO kk = 1,nkin
      DO np = 1,nreactmin(kk)
        mnrl = mnrl + 1
      END DO
    END DO
    
    DO mm = 1,mnrl
      CALL stringlen(namrl(mm),lmin)
      ltrim = lmin - 7    ! trim "_default" off string
      string1 = namrl(mm)
      IF (ltrim >= 1 .AND. lmin <= mls) THEN
        IF (string1(ltrim:lmin) == '_default') THEN
          namrl(mm) = string1(1:lmin-8)
        END IF
      END IF
    END DO
    
    IF (igenericrates == 1) THEN
      rate0 = RateGeneric
      ea = 0.0
      imintype = 1
      ndepend = 0
      WRITE(*,*)
      WRITE(*,*) ' Using generic rates--kinetic database skipped'
      WRITE(*,*)
    END IF
    
!!!  Now, check for nucleation (or not)
    
  ELSE
    WRITE(*,*) ' No kinetic reactions involving minerals specified'
    WRITE(*,*)
    nkin = 0
    nrct = 0
    mnrl = 0
  END IF
  
ELSE
  
  nrct = 0
  nkin = 0
  mnrl = 0
  DO k = 1,mrct
    nreactmin(k) = 1
  END DO
  
END IF

WRITE(iunit2,*)
WRITE(*,*)
WRITE(*,*)      ' Number of kinetic minerals:         ',nkin
WRITE(iunit2,*) 'Number of kinetic minerals = ',nkin
WRITE(iunit2,*)

!      do k = 1,nkin
!        write(*,1025) k,namk(k)
!        write(iunit2,1025) k,namk(k)
!      end do
!      if (nkin.eq.0) then
!        write(*,*)
!        write(*,*) ' No kinetically reacting minerals in input file'
!        write(*,*)
!      endif

1040 FORMAT(i3,1X,a18)


!******************************************
jh2o = 0
jph = 0
joh = 0
jo2 = 0
jco2 = 0
DO j=1,ncomp
  IF(namc(j) == 'h2o') jh2o=j
  IF(namc(j) == 'H2O') jh2o=j
  IF(namc(j) == 'H2o') jh2o=j
END DO

loop = 1
IF (jh2o /= 0 ) THEN
  loop = 2
  WRITE(iunit2,*) 'water is present: loop = ',loop
END IF

!***********************************************************************
!         *********** SURFACE COMPLEXATION SECTION **************

section='surface_complexation'
CALL readblock(nin,nout,section,found,ncount)

IF (found) THEN
  
!!  CALL FormatForNameList(nout,NameListFormat)

  CALL read_surface(nout,ncomp,nkin,nsurf)
  
  WRITE(*,*)
  WRITE(*,*)      ' Number of surface complexes:        ',nsurf
  DO ks = 1,nsurf
    namtemp = namsurf(ks)
    string1 = umin(ksurf(ks))
    CALL stringlen(namtemp,lsurf)
    CALL stringlen(string1,lmin)
!!    WRITE(*,*) namtemp(1:lsurf), ' on: ',string1(1:lmin)
  END DO
ELSE
  WRITE(*,*)
  WRITE(*,*) ' No surface complexation block found'
END IF


!     ****

WRITE(*,*) 
WRITE(*,*) '    --> Calling database subroutine'
WRITE(*,*)

IF (icomplete == 1) THEN
  CALL database(ncomp,ncmplx,mnrl,nrct,ngas,nsurf,nsurf_sec,ntemp,  &
      iprint,icomplete,jtemp,data1,namc,namcx,namrl,  &
      vbar,wtminfull,wtbas,coef,temp,z,bdotPrimary,zx,bdotSecondary,a0,ax0)
  nrct = mnrl
ELSE
  CALL database(ncomp,ncmplx,mnrl,nrct,ngas,nsurf,nsurf_sec,ntemp,  &
      iprint,icomplete,jtemp,data1,namc,namcx,namrl,  &
      vbar,wtminfull,wtbas,coef,temp,z,bdotPrimary,zx,bdotSecondary,a0,ax0)
END IF

WRITE(*,*) 
WRITE(*,*) '    --> Finished reading database'

!  Update secondary species names added during sweep of EQ3 database

DO ksp = 1,ncmplx
  ulab(ncomp+ksp) = namcx(ksp)
END DO
IF (icomplete == 1) THEN
  DO k = 1,mnrl
    umin(k) = namrl(k)
  END DO
END IF

IF (GMSsecondary) THEN
  OPEN(UNIT=51,FILE='GMSsecondary.out',STATUS='unknown')
  DO ksp = 1,ncmplx
    namtemp = ulab(ksp+ncomp)
    CALL stringlen(namtemp,lssec)
    WRITE(51,*) namtemp(1:lssec)
  END DO
  CLOSE(UNIT=51,STATUS='KEEP')
  STOP
END IF

IF (GMSmineral) THEN
  OPEN(UNIT=51,FILE='GMSmineral.out',STATUS='unknown')
  DO k = 1,mnrl
    namtemp = umin(k)
    CALL stringlen(namtemp,lsmin)
    WRITE(51,*) namtemp(1:lsmin)
  END DO
  CLOSE(UNIT=51,STATUS='KEEP')
  STOP
END IF

IF (nsurf > 0) THEN
  IF (ALLOCATED(zsurf)) THEN
    DEALLOCATE(zsurf)
  END IF
  ALLOCATE(zsurf(nsurf+nsurf_sec))
  !!!CLOSE(nout,STATUS='delete')
  OPEN(UNIT=iunit5,FILE=data1,STATUS='old',ERR=334)
  !!!OPEN(UNIT=nout,FILE='CrunchJunk2.out',STATUS='unknown') 
  REWIND nout
  string1 = 'Begin surface complexation parameters'
  CALL find_string(iunit5,string1,ifind)
  IF (ifind == 0) THEN
    WRITE(*,*)
    WRITE(*,*) '   Beginning of surface complexation parameters section not found in database'
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
  602 READ(IUNIT5,'(a)') string1

  IF (string1 == 'End surface complexation parameters') THEN
    CONTINUE
  ELSE
    WRITE(nout,*) string1
    GO TO 602
  END IF

  DO is = 1,nsurf
    namtemp = namsurf(is)
    CALL stringlen(namtemp,ls)
    REWIND nout
    parchar = namtemp
    parfind = ' '
    realjunk = 0.0
    CALL readCaseSensitivePar(nout,lchar,parchar,parfind,realjunk,section)
    IF (parfind == ' ') THEN  ! Parameter timestep_max not found
      WRITE(*,*)
      WRITE(*,*) '    Surface complex parameters not found in database '
      WRITE(*,*) '    Looking in database file'
      WRITE(*,*) '  Searching for: ', namtemp(1:ls)
      WRITE(*,*)
      READ(*,*)
      STOP
    ELSE
      zsurf(is) = realjunk
    END IF
  END DO

  DO ns = 1,nsurf_sec
    namtemp = namsurf_sec(ns)
    CALL stringlen(namtemp,ls)
    REWIND nout
    parchar = namtemp
    parfind = ' '
    realjunk = 0.0
    CALL readCaseSensitivePar(nout,lchar,parchar,parfind,realjunk,section)
    IF (parfind == ' ') THEN  ! Parameter timestep_max not found
      WRITE(*,*)
      WRITE(*,*) '    Surface complex parameters not found in database '
      WRITE(*,*) '    Looking in database file'
      WRITE(*,*) '  Searching for: ', namtemp(1:ls)
      WRITE(*,*)
      READ(*,*)
      STOP
    ELSE
      zsurf(ns+nsurf) = realjunk
    END IF

  END DO

  CLOSE(iunit5,STATUS='keep')

END IF

!  If ncmplx > ncmplx_old (additional secondary species added from
!  the sweep of the EQ3 database), then write out the new species list

IF (ncmplx_old < ncmplx) THEN
  WRITE(iunit2,*)
  WRITE(iunit2,*) '   NEW SECONDARY SPECIES LIST '
  WRITE(iunit2,*)
  DO ksp = 1,ncmplx
    namtemp = ulab(ksp+ncomp)
    CALL stringlen(namtemp,lssec)
    WRITE(iunit2,1050) namtemp
  END DO
  WRITE(iunit2,*)
1050 FORMAT(a72)
  
  WRITE(iunit2,*)
  WRITE(iunit2,*) '   NEW MINERAL LIST '
  WRITE(iunit2,*)
  DO k = 1,mnrl
    namtemp = umin(k)
    CALL stringlen(namtemp,lsmin)
    WRITE(iunit2,*) namtemp(1:lsmin)
  END DO
  WRITE(iunit2,*)
END IF

!  951 format(a<lssec>)
!  952 format(a<lsmin>)
951 FORMAT(a35)
952 FORMAT(a35)


!  Check to see if mnrl > nkin and add the names of the minerals to the
!  mineral list

!      if (mnrl.gt.nkin) then
!        write(*,*)
!        write(*,*) ' Adding additional minerals from EQ3 database'
!        write(*,*) ' Minerals will not react, however'
!        write(*,*)
!        write(*,*) '   NEW MINERAL LIST '
!        write(*,*)
!        do k = 1,mnrl
!          write(*,*) umin(k)
!        end do
!        write(*,*)
!        write(iunit2,*)
!        write(iunit2,*) ' Adding additional minerals from EQ3 database'
!        write(iunit2,*) ' Minerals will not react, however'
!        write(iunit2,*)
!        write(iunit2,*) '   NEW MINERAL LIST '
!        write(iunit2,*)
!        do k = 1,mnrl
!          write(iunit2,*) umin(k)
!        end do
!        write(iunit2,*)
!      endif


!************************
!  Make Lichtnerian variables compatible with those of Steefel
!************************

DO ik = 1,ncomp
  nbasin(ik) = ik
!         ulab(ik) = namc(ik)
  chg(ik) = z(ik)
  azero(ik) = a0(ik)
  bdotParameter(ik) = bdotPrimary(ik)
END DO

DO ik = 1,ncmplx
  nbkin(ik) = ik+ncomp
!         ulab(ik+ncomp) = namcx(ik)
  chg(ik+ncomp) = zx(ik)
  azero(ik+ncomp) = ax0(ik)
  bdotParameter(ik+ncomp) = bdotSecondary(ik)
END DO

IonsOK = .FALSE.
!!! Check to see there are ions present, otherwise it blow up in the activity of water calculation
DO ik = 1,ncomp+ncmplx
  IF (chg(ik) /= 0.0d0) THEN
    IonsOK = .TRUE.
  END IF
END DO

IF (.NOT. IonsOK) THEN
  write(*,*)
  write(*,*) ' You need ions in solution for ionic strength master variable'
  write(*,*) ' You can add to CONDITIONS (for example) '
  write(*,*)
  write(*,*) ' Na+    0.001 '
  write(*,*) ' Cl-    0.001 '
  write(*,*)
  stop
END IF

WRITE(iunit2,*)
DO m = 1, mnrl
  vbarkin(m)=vbar(m)
END DO

msub = 0
DO k = 1,nkin
  DO ll = 1,nreactmin(k)
    msub = msub + 1
    IF (ll == 1) THEN
      volmol(k) = vbarkin(msub)
      wtmin(k) = wtminfull(msub)
    END IF
    rate0(ll,k) = 10**(rate0(ll,k))/yrsec
  END DO
END DO

DO i = 1,ncomp
  wtcomp(i) = wtbas(i)
END DO

DO ik = 1,ndim
  as1(ik,1) = coef(ik,1)
  as1(ik,2) = coef(ik,2)
  as1(ik,3) = coef(ik,3)
  as1(ik,4) = coef(ik,4)
  as1(ik,5) = coef(ik,5)
END DO

DO ik = 1,ncomp+ncmplx
  acmp(ik) = azero(ik)
END DO

DEALLOCATE(namc)
DEALLOCATE(namcx)
DEALLOCATE(namrl)
DEALLOCATE(nammin_tmp)
DEALLOCATE(dumlabel)
DEALLOCATE(vbar)
DEALLOCATE(wtminfull)
DEALLOCATE(wtbas)
DEALLOCATE(a0)
DEALLOCATE(z)
DEALLOCATE(ax0)
DEALLOCATE(zx)
DEALLOCATE(bdotPrimary)
DEALLOCATE(bdotSecondary)
DEALLOCATE(vbarkin)
!!!DEALLOCATE(azero)
DEALLOCATE(coef)

RETURN

1111 FORMAT('  # ','  Mineral    ', '      idep', ' sat1 ','sat2  ',  &
    '  rk_25  ','     Ea  ','    thresh ')

333 WRITE(iunit2,*) 'error from readat: ',label
WRITE(*,*) 'error from readat: ',label

1110 FORMAT(i3,1X,a18,2X,i2,2X,f4.1,2X,f4.1,2X, 1PE10.2,2X,1PE10.2,2X,f4.1)
2260 FORMAT(' Revised logK for mineral: ',a20,' logK = ',1PE12.4)
1  FORMAT(1(/))
2  FORMAT(2(/))
3  FORMAT(3(/))
4  FORMAT(4(/))
5  FORMAT(5(/))
6  FORMAT(6(/))
7  FORMAT(7(/))
8  FORMAT(8(/))
9  FORMAT(9(/))
10  FORMAT(10(/))
11  FORMAT(11(/))
12  FORMAT(12(/))
13  FORMAT(13(/))

READ(*,*)
STOP
334 WRITE(*,*) 'Error in opening database file to read surface complex parameters: stop'
READ(*,*)
STOP

END SUBROUTINE read98

!*********************************************************************


