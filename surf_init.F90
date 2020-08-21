SUBROUTINE surf_init(ncomp,nspec,nsurf,nsurf_sec,nco)
USE crunchtype
USE concentration
USE mineral

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: ncomp
INTEGER(I4B), INTENT(IN)                                    :: nspec
INTEGER(I4B), INTENT(IN)                                    :: nsurf
INTEGER(I4B), INTENT(IN)                                    :: nsurf_sec
INTEGER(I4B), INTENT(IN)                                    :: nco


!  Internal variables and arrays

INTEGER(I4B)                                                :: ns
INTEGER(I4B)                                                :: i
INTEGER(I4B)                                                :: is

REAL(DP)                                                    :: sum
REAL(DP)                                                    :: delta_z
REAL(DP)                                                    :: activity
REAL(DP)                                                    :: LogTotalSites
REAL(DP)                                                    :: LogTotalEquivalents
REAL(DP)                                                    :: LogDependence

DO ns = 1,nsurf_sec

  IF (nptlink(ns) /= 0) THEN
      
!!!  Primary surface complex charge minus secondary surface complex charge 

    delta_z = zsurf(islink(ns)) - zsurf(ns+nsurf)

!!  Aqueous species
    sum = 0.0d0
    DO i = 1,ncomp
      IF (ulab(i) == 'H2O') THEN
      sum = sum + musurf(ns,i)*(gamtmp(i))
      ELSE
      sum = sum + musurf(ns,i)*(sptmp(i) + gamtmp(i))
      END IF
    END DO

    LogTotalSites = DLOG(c_surf(islink(ns),nco))

!!  Surface complexes
    DO is = 1,nsurf
!!!      activity = spsurftmp(is) - LogTotalSites
      activity = spsurftmp(is)
      sum = sum + musurf(ns,is+ncomp)*activity
    END DO

    spsurftmp(ns+nsurf) = keqsurf_tmp(ns) + sum                               &
        + 2.0d0*musurf(ns,islink(ns)+ncomp)*delta_z*LogPotential_tmp(nptlink(ns))   &
        - (musurf(ns,islink(ns)+ncomp)-1.0d0)*LogTotalSites                   &
        - DLOG(musurf(ns,islink(ns)+ncomp)) 

    spsurftmp10(ns+nsurf) = DEXP(spsurftmp(ns+nsurf))
    continue

  ELSE                                  !! Non-electrostatic option

!!  Aqueous species
    sum = 0.0
    DO i = 1,ncomp
      IF (ulab(i) == 'H2O') THEN
        sum = sum + musurf(ns,i)*(gamtmp(i))
      ELSE
        sum = sum + musurf(ns,i)*(sptmp(i) + gamtmp(i))
      END IF
    END DO

    LogTotalSites = DLOG(c_surf(islink(ns),nco))

!!  Surface complexes
    DO is = 1,nsurf
      activity = spsurftmp(is)
      sum = sum + musurf(ns,is+ncomp)*activity
    END DO
    
    spsurftmp(ns+nsurf) = keqsurf_tmp(ns) + sum                               &
        - (musurf(ns,islink(ns)+ncomp)-1.0d0)*LogTotalSites                   &
        - DLOG(musurf(ns,islink(ns)+ncomp)) 
    
    spsurftmp10(ns+nsurf) = DEXP(spsurftmp(ns+nsurf))
    
  END IF

END DO

RETURN
END SUBROUTINE surf_init
!  **************************************************************
