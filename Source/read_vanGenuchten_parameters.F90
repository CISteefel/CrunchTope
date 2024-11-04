SUBROUTINE read_vanGenuchten_parameters(nout,parchar,parfind, nx,ny,nz,nzones_VG_params,VG_error)
USE crunchtype
USE CrunchFunctions
USE params
USE flow
USE strings

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nout
INTEGER(I4B), INTENT(IN)                                    :: nx
INTEGER(I4B), INTENT(IN)                                    :: ny
INTEGER(I4B), INTENT(IN)                                    :: nz
INTEGER(I4B), INTENT(OUT)                                   :: nzones_VG_params
INTEGER(I4B), INTENT(INOUT)                                 :: VG_error ! error flag for reading van Genuchten parameters
CHARACTER (LEN=mls), INTENT(IN)                             :: parchar ! character for each van Genuchten parameter in the input file
CHARACTER (LEN=mls), INTENT(IN OUT)                         :: parfind

!  Internal variables and arrays

INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs
INTEGER(I4B)                                                :: nxyz
INTEGER(I4B)                                                :: nlen1
INTEGER(I4B)                                                :: ls_a
INTEGER(I4B)                                                :: ls_b
INTEGER(I4B)                                                :: l

REAL(DP)                                                    :: VG_params_tmp

nxyz = nx*ny*nz

VG_params_zone(0) = 0.0
REWIND nout

nzones_VG_params = 0
10 READ(nout,'(a)',END=500) zone
nlen1 = LEN(zone)
CALL majuscules(zone,nlen1)
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  CALL convan(ssch,lzs,res)
  IF (ssch == parchar) THEN
    parfind = parchar
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF(ls /= 0) THEN
      lzs=ls
      CALL convan(ssch,lzs,res)
      IF (res == 'n') THEN
        VG_params_tmp = DNUM(ssch)
      ELSE                !  An ascii string--so bag it.
        WRITE(*,*)
        WRITE(*,*) ' Cant interpret string following a VG parameter', parchar
        WRITE(*,*) ' Looking for numerical value'
        WRITE(*,*)
        VG_error = 1
        RETURN
      END IF
      
! Now look for ASCII string indicating location of a VG parameter
      
      id = ids + ls
      CALL sschaine(zone,id,iff,ssch,ids,ls)
      IF(ls /= 0) THEN
        lzs=ls
        CALL convan(ssch,lzs,res)
        IF (res == 'a') THEN
          IF (ssch == 'default' .OR. ssch == 'all') THEN
            VG_params_zone(0) = VG_params_tmp
          ELSE IF (ssch == 'zone') THEN
            
!  "Zone" specified, so look for locations
            
            nzones_VG_params = nzones_VG_params + 1
            IF (nzones_VG_params > mperm) THEN
              WRITE(*,*)
              WRITE(*,*)  ' Number of VG parameters zones dimensioned too small'
              WRITE(*,*)  ' Number of VG parameters zones = ',nzones_VG_params
              WRITE(*,*)  ' Dimension of VG parameters zones = ',mperm
              WRITE(*,*)  ' Contact the code developer at CISteefel@lbl.gov '
              WRITE(*,*)
              VG_error = 1
              RETURN
            END IF
            
            VG_params_zone(nzones_VG_params) = VG_params_tmp
            
            id = ids + ls
            CALL sschaine_hyph(zone,id,iff,ssch_a,ssch_b,ids,ls_a,ls_b,ls)
            IF(ls /= 0) THEN
              lzs=ls_a
              CALL convan(ssch_a,lzs,res)
              IF (res == 'n') THEN
                jxx_VG_params_lo(nzones_VG_params) = JNUM(ssch_a)
              ELSE                !  An ascii string--so bag it.
                WRITE(*,*)
                WRITE(*,*) ' A grid location should follow zone specification'
                WRITE(*,*) ' Dont know what to do with this string', parchar
                WRITE(*,*)
                VG_error = 1
                RETURN
              END IF
              IF (ls_b /= 0) THEN
                lzs=ls_b
                CALL convan(ssch_b,lzs,res)
                IF (res == 'n') THEN
                  jxx_VG_params_hi(nzones_VG_params) = JNUM(ssch_b)

                ELSE                !  An ascii string--so bag it.
                  WRITE(*,*)
                  WRITE(*,*) ' A grid location should follow zone specification'
                  WRITE(*,*) ' Dont know what to do with this string', parchar
                  WRITE(*,*)
                  VG_error = 1
                  RETURN
                END IF
              ELSE
                jxx_VG_params_hi(nzones_VG_params) = jxx_VG_params_lo(nzones_VG_params)   !  Assume jxx_VG_params_hi=jxx_VG_params_lo
              END IF
            ELSE                  ! Zero length trailing string
              WRITE(*,*)
              WRITE(*,*) ' No X or Y grid location given for a VG parameter', parchar
              WRITE(*,*) ' VG parameter zone ',nzones_VG_params
              WRITE(*,*)
              VG_error = 1
              RETURN
            END IF
            
            WRITE(*,*)
            WRITE(*,*) ' VG parameter zone number ',nzones_VG_params
            WRITE(*,*) ' jxx_VG_params_lo = ', jxx_VG_params_lo(nzones_VG_params)
            WRITE(*,*) ' jxx_VG_params_hi = ',jxx_VG_params_hi(nzones_VG_params)
            WRITE(*,*)
            
            id = ids + ls
            CALL sschaine_hyph(zone,id,iff,ssch_a,ssch_b,ids,ls_a,ls_b,ls)
            IF(ls /= 0) THEN
              lzs=ls_a
              CALL convan(ssch_a,lzs,res)
              IF (res == 'n') THEN
                jyy_VG_params_lo(nzones_VG_params) = JNUM(ssch_a)
              ELSE                !  An ascii string--so bag it.
                WRITE(*,*)
                WRITE(*,*) ' No Y location for a VG parameter', parchar
                WRITE(*,*)
                VG_error = 1
                RETURN
              END IF
              IF (ls_b /= 0) THEN
                lzs=ls_b
                CALL convan(ssch_b,lzs,res)
                IF (res == 'n') THEN
                  jyy_VG_params_hi(nzones_VG_params) = JNUM(ssch_b)
                ELSE                !  An ascii string--so bag it.
                  WRITE(*,*)
                  WRITE(*,*) ' A grid location should follow zone specification'
                  WRITE(*,*) ' Dont know what to do with this string after a VG parameter', parchar
                  WRITE(*,*)
                  VG_error = 1
                  RETURN
                END IF
              ELSE
                jyy_VG_params_hi(nzones_VG_params) = jyy_VG_params_lo(nzones_VG_params)   !  Assume jyy_VG_params_hi=jyy_VG_params_lo
              END IF
            ELSE                  ! Zero length trailing string
              WRITE(*,*)
              WRITE(*,*) ' No Y location given for a VG parameter zone'
              WRITE(*,*) ' VG parameter zone number ',nzones_VG_params
              WRITE(*,*)
              VG_error = 1
              RETURN
            END IF
              
            WRITE(*,*)
            WRITE(*,*) ' jyy_VG_params_lo = ',jyy_VG_params_lo(nzones_VG_params)
            WRITE(*,*) ' jyy_VG_params_hi = ',jyy_VG_params_hi(nzones_VG_params)
            WRITE(*,*)    
            
            id = ids + ls
            CALL sschaine_hyph(zone,id,iff,ssch_a,ssch_b,ids,ls_a,ls_b,ls)
            IF(ls /= 0) THEN
              lzs=ls_a
              CALL convan(ssch_a,lzs,res)
              IF (res == 'n') THEN
                jzz_VG_params_lo(nzones_VG_params) = JNUM(ssch_a)
              ELSE                !  An ascii string--so bag it.
                WRITE(*,*)
                WRITE(*,*) ' No Z location for a VG parameter ', parchar
                WRITE(*,*)
                VG_error = 1
                RETURN
              END IF
              IF (ls_b /= 0) THEN
                lzs=ls_b
                CALL convan(ssch_b,lzs,res)
                IF (res == 'n') THEN
                  jzz_VG_params_hi(nzones_VG_params) = JNUM(ssch_b)
                ELSE                !  An ascii string--so bag it.
                  WRITE(*,*)
                  WRITE(*,*) ' A grid location should follow zone specification'
                  WRITE(*,*) ' Dont know what to do with this string after a VG parameter', parchar
                  WRITE(*,*)
                  VG_error = 1
                  RETURN  
                END IF
              ELSE
                jzz_VG_params_hi(nzones_VG_params) = jzz_VG_params_lo(nzones_VG_params)   !  Assume jzz_VG_params_hi=jzz_VG_params_lo
              END IF
            ELSE                  ! Zero length trailing string
              WRITE(*,*)
              WRITE(*,*) ' No Z location given for a VG parameter zone'
              WRITE(*,*) ' VG parameter zone number ',nzones_VG_params
              WRITE(*,*)
              VG_error = 1
              RETURN
            END IF
              
            WRITE(*,*)
            WRITE(*,*) ' jzz_VG_params_lo = ',jzz_VG_params_lo(nzones_VG_params)
            WRITE(*,*) ' jzz_VG_params_hi = ',jzz_VG_params_hi(nzones_VG_params)
            WRITE(*,*)    

          ELSE
            WRITE(*,*)
            WRITE(*,*) ' Dont understand string following a VG parameter value'
            WRITE(*,*) ssch(1:ls)
            WRITE(*,*)
            VG_error = 1
            RETURN
          END IF
          
        ELSE                !  A number--so bag it.
          WRITE(*,*)
          WRITE(*,*) ' Cant interpret string following a VG parameter value'
          WRITE(*,*) ' Looking for ASCII string'
          WRITE(*,*)
          VG_error = 1
          RETURN
        END IF
      ELSE   ! Assume this is default if nothing else given
        VG_params_zone(0) = VG_params_tmp
      END IF
    ELSE
      WRITE(*,*)
      WRITE(*,*) ' No value given for a VG parameter', parchar
      VG_error = 1
      RETURN
    END IF
  ELSE
    GO TO 10
  END IF
  
END IF

GO TO 10

500 DO l = 1,nzones_VG_params
  IF (jxx_VG_params_hi(l) > nx) THEN
    WRITE(*,*)
    WRITE(*,*) 'You have specified a VG parameter at JX > NX'
    WRITE(*,*)
    VG_error = 1
    RETURN
  END IF
  IF (jyy_VG_params_hi(l) > ny) THEN
    WRITE(*,*)
    WRITE(*,*) 'You have specified a VG parameter at JY > NY'
    WRITE(*,*)
    VG_error = 1
    RETURN
  END IF
  IF (jzz_VG_params_hi(l) > nz) THEN
    WRITE(*,*)
    WRITE(*,*) 'You have specified a VG parameter at JZ > NZ'
    WRITE(*,*)
    VG_error = 1
    RETURN
  END IF
  IF (jxx_VG_params_lo(l) < 0) THEN
    WRITE(*,*)
    WRITE(*,*) 'You have specified a VG parameter at JX < 0'
    WRITE(*,*)
    VG_error = 1
    RETURN
  END IF
  IF (jyy_VG_params_lo(l) < 1) THEN
    WRITE(*,*)
    WRITE(*,*) 'You have specified a VG parameter at JY < 1'
    WRITE(*,*)
    VG_error = 1
    RETURN
  END IF
  IF (jzz_VG_params_lo(l) < 1) THEN
    WRITE(*,*)
    WRITE(*,*) 'You have specified a VG parameter at JZ < 1'
    VG_error = 1
    RETURN
  END IF
END DO

RETURN
END SUBROUTINE read_vanGenuchten_parameters