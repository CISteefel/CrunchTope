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


SUBROUTINE read_SolidDensity(nout,nchem)
USE crunchtype
USE CrunchFunctions
USE params
USE strings
USE concentration

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nout
INTEGER(I4B), INTENT(IN)                                    :: nchem


!  Internal variables and arrays

INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs
INTEGER(I4B)                                                :: nlen1
INTEGER(I4B)                                                :: lsave

!  This routine reads the solid density, or instructions to calculate this from mineral volume fractions

SolidDensityFrom(nchem) = 3               !! Density from mineral volume fractions

REWIND nout

100 READ(nout,'(a)',END=300) zone
nlen1 = LEN(zone)
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  CALL convan(ssch,lzs,res)
  lsave = lzs
  IF (ssch == 'soliddensity' .OR. ssch == 'solid_density') THEN            !! SolidDensity keyword found
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF (ls /= 0) THEN
      lzs=ls
      CALL convan(ssch,lzs,res)
      lsave = lzs
      IF (res == 'n') THEN
        SolidDensity(nchem) = DNUM(ssch)
        SolidDensityFrom(nchem) = 1
      ELSE
        IF (ssch == '-ss_ratio' .OR. ssch == 'solidsolutionratio' .OR. ssch == 'solidsolution_ratio' .OR. ssch == 'solid_solution__ratio') THEN
          id = ids + ls
          CALL sschaine(zone,id,iff,ssch,ids,ls)
          IF (ls /= 0) THEN
            lzs=ls
            CALL convan(ssch,lzs,res)
            lsave = lzs
            IF (res == 'n') THEN
              SolidSolutionRatio(nchem) = DNUM(ssch)
              SolidDensityFrom(nchem) = 2
            ELSE
              WRITE(*,*)
              WRITE(*,*) ' A numerical value for the solid:solution ratio'
              WRITE(*,*) '   should follow the parameter flag "SolidSolutionRatio" '
              WRITE(*,*) ssch(1:lsave)
              WRITE(*,5050) nchem
              WRITE(*,*)
              READ(*,*) 
              STOP
            END IF
          ELSE
            WRITE(*,*)
            WRITE(*,*) ' No values given for solid solution ratio'
            WRITE(*,*) ssch(1:lsave)
            WRITE(*,5050) nchem
            WRITE(*,*)
            READ(*,*) 
            STOP
          END IF
        ELSE IF (ssch == 'calculate' .OR. ssch == 'calculatefromminerals' .OR. ssch == 'calculatefrommineral') THEN     !!  Solid density (and solid:solution ratio) to be calculated from mineral volume fractions
          SolidDensityFrom(nchem) = 3
        ELSE
          SolidDensityFrom(nchem) = 3                   !! Default
        END IF
      END IF
    ELSE                                                !! Nothing following the keyword "SolidDensity"
      WRITE(*,*)
      WRITE(*,*) ' Nothing following the keyword "SolidDensity"'
      WRITE(*,5050) nchem
      WRITE(*,*)
      READ(*,*) 
      STOP
    END IF
    GO TO 300                                           !! Exit if the "SolidDensity" keyword has been found
  ELSE
    GO TO 100
  END IF
ELSE
  GO TO 100
END IF

5050 FORMAT(1X,'Condition number ',i2,' in input file')

300  RETURN
END SUBROUTINE read_SolidDensity
