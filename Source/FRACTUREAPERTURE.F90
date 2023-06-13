! --------------------------------------------------
!!!MODULE FRACTUREAPERTURE         ! global declarations
! --------------------------------------------------
!!!INCLUDE 'C2F.FD' 
!!!INCLUDE 'C2F.FI' 

!!!END MODULE
 
! --------------------------------------------------
SUBROUTINE rmesh51(nx,ny)
! --------------------------------------------------
!!!USE FRACTUREAPERTURE
USE CrunchType
USE concentration, ONLY: jinit
USE medium, ONLY: dxx,dyy,x,y
USE flow, ONLY: permx,permy
USE medium, ONLY: por
USE mineral, ONLY: area


IMPLICIT NONE

INTEGER(I4B), INTENT(IN)                                     :: nx
INTEGER(I4B), INTENT(IN)                                     :: ny
  
REAL(DP), DIMENSION(:,:), ALLOCATABLE                         :: fracseg
REAL(DP), DIMENSION(:,:), ALLOCATABLE                         :: cc
REAL(DP), DIMENSION(:,:), ALLOCATABLE                         :: mt
\

INTEGER(I4B)                                                  :: jx
INTEGER(I4B)                                                  :: jy
INTEGER(I4B)                                                  :: idummy

INTEGER(I4B)                                                  :: nfe
INTEGER(I4B)                                                  :: nfxe
INTEGER(I4B)                                                  :: nfxw
INTEGER(I4B)                                                  :: nfye
INTEGER(I4B)                                                  :: nfyw
INTEGER(I4B)                                                  :: nfs
INTEGER(I4B)                                                  :: ne1
INTEGER(I4B)                                                  :: ne2

INTEGER(I4B)                                                  :: i
INTEGER(I4B)                                                  :: j
INTEGER(I4B)                                                  :: jinitPrint

INTEGER(I4B)                                                  :: i00
INTEGER(I4B)                                                  :: i01
INTEGER(I4B)                                                  :: i02
INTEGER(I4B)                                                  :: i03
INTEGER(I4B)                                                  :: j00
INTEGER(I4B)                                                  :: j01
INTEGER(I4B)                                                  :: j02
INTEGER(I4B)                                                  :: j03

REAL(DP)                                                      :: t1

!!!ALLOCATE(fracseg(8999,40))    
ALLOCATE(cc(nx*ny,10))
cc = 0.0

!!!  So 2nd index = 8 is fracture segment number
!!!  Index = 9 permeability
!!!  

! - - - begin - - -

!!!               /*------------------------------------------------*/
!!!               for (i = 1; i <= 8999; i++)
!!!               {
!!!                   for (j = 0; j <= 10; j++)
!!!                   {
!!!                       fracseg[i][j] = 0;
!!!                   }
!!!               }//initialize the cc array
!!!              /*------------------------------------------------*/

  !!!fracseg = 0
  nfs = 0
  
  !!!WRITE (fl1,901) "number of fracture segments: ",nfs
  !!!               for (i = 1; i <= nc; i++)
  !!!             {
  !!!DO i = 1,nx*ny
  
  do jy = 1,ny
    do jx = 1,nx
	  i = (jy-1) * nx + jx
	  cc(i,1) = x(jx)
  end do
  end do 
  
  do jy = 1,ny
    do jx = 1,nx
	  i = (jy-1) * nx + jx
	  cc(i,2) = y(jy)
	end do
  end do 
  
   do jy = 1,ny
    do jx = 1,nx
	  i = (jy-1) * nx + jx
	  cc(i,4) = dxx(jx)
	end do
  end do 
  
  do jy = 1,ny
    do jx = 1,nx
	  i = (jy-1) * nx + jx
	  cc(i,3) = dyy(jy)
	end do
  end do 

  do jy = 1,ny
    do jx = 1,nx
	  i = (jy-1) * nx + jx
    jinitPrint = jinit(jx,jy,1)
    if (jinitPrint == 3) then
      jinitPrint = 1    !! Call it rockmatrix
    end if
	  cc(i,5) = DFLOAT(jinitPrint)
	end do
  end do
  
  !!!cc(:,8) = 0
  
  do jy = 1,ny
    do jx = 1,nx
	    i = (jy-1) * nx + jx
	    cc(i,7) = dxx(jx)*dyy(jy)
	  end do
  end do
  
  do jy = 1,ny
    do jx = 1,nx
	    i = (jy-1) * nx + jx
	    cc(i,10) = por(jx,jy,1)
	  end do
  end do
   
  ne1 = nx
  ne2 = ny
  
  DO jy = 1,ny
    DO jx = 1,nx
	
	  i = (jy-1) * nx + jx
    i00 = i
  
    nfe = 0
    nfxe = 0
    nfye = 0
	
	!!  NOTE: The are the C indices, so Fortran would increment by 1
	!!  cc(#,0)
	!!  cc(#,1) x coordinate (cell center)
	!!  cc(#,2) y coordinate (cell center)
	!!  cc(#,3) Delta y
	!!  cc(#,4) Delta x
	!!  cc(#,5) material type
	!!  cc(#,6)
	!!  cc(#,7) volume of cell, dx*dy
	!!  cc(#,8) fracture segment number
	!!  cc(#,9) permeability
	
	!!! NOTE: Indices switched in transforming from C to Fortran
	!!! NOTE: Adding 1 to adjust array from C (starting at 0) to Fortran (starting at 1)
	

    if (dabs(cc(i,5) - 1) < 0.0001)                    GO TO 10  !!! Material type
      if (dabs(cc(i,8)) > 0.0001)                      GO TO 10  !! This element is already identified to one fracture segment  
        IF (MOD(i,ne1) == 0 .AND. i > (ne2 - 1) * ne1) GO TO 10  !!! top right corner, so exit
	  
          nfs = nfs + 1
          nfxe = nfxe + 1
          nfye = nfye + 1
          i00 = i
	
	      !!!  ne1 is nx, ne2 is ny
          
        IF (MOD(i00,ne1) == 0)                                GO TO 20    !!! right column
	
          IF (dABS(cc(i00,5) - cc(i00+1,5) ) > 0.0001)        GO TO 20     !!! Case where MT are different at i and i+1
      
    DO j = 0,ne1
      IF ((i00 + j) / ne1 /= (i00 + j + 1) / ne1)             GO TO 20    !!! Different rows
        IF (dABS(cc(i00 + j,5) - cc(i00 + j+1,5)) > 0.0001)      GO TO 20    !!! Exit when different MT
      nfxe = nfxe + 1
    END DO
20  CONTINUE
    
  
    IF (i00 > (ne2 - 1) * ne1)                                GO TO 30    !!!  Top row, so don't look up
      IF (dABS(cc(i00,5) - cc(i00+ne1,5)) > 0.0001)            GO TO 30    !!!  Different MT
    DO j = 0,ne2
      j00 = i00 + j * ne1
      j01 = i00 + (j + 1) * ne1
      IF (j00 > (ne2 - 1) * ne1)                              GO TO 30    !!!  Top row, so exit
        IF (dABS(cc(j00,5) - cc(j01,5)) > 0.0001)              GO TO 30    !!!  If different MT, then exit
      nfye = nfye + 1
    END DO
30  CONTINUE
    
	if (i == 514) then
    continue
  end if
  
    IF (nfxe > 1 .and. cc(i00,4) * nfxe <= cc(i00,3) * nfye) THEN            !!  Choose shortest aperture
      nfe = nfxe                                              !!! Number of elements in fracture segment
      t1 = 0.0d0
      do idummy = 0,nfxe-1
        t1 = t1 + cc(i00+idummy,10)*cc(i00+idummy,4)
        t1 = t1 + cc(i00+idummy,4)
      end do
!!!      t1 = cc(i00,4) * nfxe                                   !! dy*number of elements --> Width in X direction
      i02 = 1
      j02 = 0
	  
	  !!!  Nfxe is number of cells in X direction
	  
      GO TO 40
    END IF
	
    IF (nfye > 1 .and. cc(i00,4) * nfxe > cc(i00,3) * nfye) THEN            !!  Choose shortest aperture
      nfe = nfye                                             !!! Number of elements in fracture segment
      t1 = 0.0d0
      do idummy = 0,nfye-1
        t1 = t1 + cc(i00+idummy,10)*cc(i00+idummy,3)
        t1 = t1 + cc(i00+idummy,3)
      end do
!!!      t1 = cc(i00,3) * nfye                                  !! dy*number of elements --> Width in Y direction
      i02 = ne1
      j02 = 1

    END IF
    40 CONTINUE
    
    j01 = nfe - 1
	
	  DO i01 = 0,j01
    
      !!! i00: Starting point
      !!! j01: End point
    
      i03 = i00 + i01 * i02
	                             
      cc(i03,8) = nfs;             !!! Fracture segment number
      cc(i03,9) = (0.01*t1 * 0.01*t1 / 12.0)   !!! Permeability		
!!!      cc(i03,9) = 1.0E-13

	    !!!apert = apert + sqrt(poro[i03]*width of the cell i03)/12.0
	  
	  END DO
	  
!!!      fracseg(nfs,1) = i00;
!!!      fracseg(nfs,2) = i00 + (nfe - 1) * i02;
!!!      fracseg(nfs,3) = j02                      !!(0 means in x direction, 1 means in y direction)
!!!      fracseg(nfs,4) = nfe
!!!      fracseg(nfs,5) = t1
  
     10 CONTINUE
	
    END DO
  end do
  
  do jy = 1,ny
    do jx = 1,nx
      i = (jy-1) * nx + jx
      permx(jx,jy,1) = cc(i,9)*1.0E-08
      permy(jx,jy,1) = cc(i,9)*1.0E-08
      permx(jx,jy,1) = cc(i,9)*0.01
      permy(jx,jy,1) = cc(i,9)*0.01
      permx(jx,jy,1) = cc(i,9)*0.01
      permy(jx,jy,1) = cc(i,9)*0.01
!!!      area(2,jx,jy,1) = 0.01*2.0/cc(i,9)
      if (permx(jx,jy,1) < 1.0E-18) then
        permx(jx,jy,1) = 1.0E-18
      end if
      if (permy(jx,jy,1) < 1.0E-18) then
        permy(jx,jy,1) = 1.0E-18
      end if
    end do
  end do
  
  !!! Add soil layer
  do jy = 1,ny
    do jx = 1,nx
      if (jy==1 .or. jy==2) then
        permx(jx,jy,1) = 4.15E-13
        permy(jx,jy,1) = 4.15E-13
      end if
    end do
  end do
  
  continue
  
!!!  WRITE (fl1,901) "number of fracture segments: ",nfs

!!!901 FORMAT (A,I0)
!!!902 FORMAT (A,I0,A,I0)
!!!903 FORMAT (A,I0,A,I0,A,I0)
  
!!!DEALLOCATE(fracseg)    
DEALLOCATE(cc)


RETURN
END SUBROUTINE

