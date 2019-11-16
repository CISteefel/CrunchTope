SUBROUTINE AqueousToBulkConvert(jx,jy,jz,AqueousToBulk)
USE crunchtype
USE concentration
USE RunTime
USE temperature
USE transport
USE medium

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: jx
INTEGER(I4B), INTENT(IN)                                    :: jy
INTEGER(I4B), INTENT(IN)                                    :: jz
REAL(DP), INTENT(OUT)                                       :: AqueousToBulk

!  Internal variables and arrays

REAL(DP)                                                    :: MeanSaltConcentration
REAL(DP)                                                    :: MassFraction

IF (DensityModule /= 'temperature') THEN     !  Want kgw/kg_solution
  MeanSaltConcentration = 0.001*(wtaq(MeanSalt(1))*sn(MeanSalt(1),jx,jy,jz) +   &
     wtaq(MeanSalt(2))*sn(MeanSalt(2),jx,jy,jz)) 
  MassFraction = 1.0/(1.0 + MeanSaltConcentration)
ELSE
  MassFraction = 1.0
END IF
AqueousToBulk = MassFraction*por(jx,jy,jz)*satliq(jx,jy,jz)*ro(jx,jy,jz)

RETURN
END SUBROUTINE AqueousToBulkConvert
