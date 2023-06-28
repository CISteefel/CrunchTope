SUBROUTINE interp3(tstart, tint, tseries, ydat, yout,info)
 

    
  USE crunchtype
  USE flow


      implicit none
    
      !  External variables and arrays
      real(DP),intent(in) :: tstart
      real(DP),intent(in) :: tint
      INTEGER(I4B),intent(in) :: info
      real(DP),intent(in) :: ydat(info)
      real(DP),intent(in) :: tseries(info)
      real(DP),intent(out) :: yout

      !  Internal variables and arrays
      INTEGER(I4B) :: index
      INTEGER(I4B) :: countval
      INTEGER(I4B) :: tpos_low
      INTEGER(I4B) :: tpos_high
      REAL(DP)  :: slope
      REAL(DP)  :: slope1
      REAL(DP)  :: slope2
      REAL(DP)  :: ymin
      REAL(DP)  :: ymax
      real(DP),dimension(:), ALLOCATABLE     :: tseries_select
      real(DP),dimension(:), ALLOCATABLE     :: tseries_within
      real(DP),dimension(:), ALLOCATABLE     :: ydat_within
      real(DP),dimension(:), ALLOCATABLE     :: weight

      real(DP),dimension(:), ALLOCATABLE     :: ydatdum
      real(DP),dimension(:), ALLOCATABLE     :: ydatdum2

      real(DP),dimension(:), ALLOCATABLE     :: tseries2
      real(DP),dimension(:), ALLOCATABLE     :: ydat2

!! Glue an additional time series to the actual time series to ensure continuity

IF (ALLOCATED(tseries2)) THEN
    DEALLOCATE(tseries2)
    ALLOCATE(tseries2(1:size(tseries)*2-1))
  ELSE
    ALLOCATE(tseries2(1:size(tseries)*2-1))
  END IF

IF (ALLOCATED(ydat2)) THEN
    DEALLOCATE(ydat2)
    ALLOCATE(ydat2(1:size(tseries)*2-1))
  ELSE
    ALLOCATE(ydat2(1:size(tseries)*2-1))
  END IF

tseries2(1:size(tseries)) = tseries
tseries2(size(tseries)+1 : size(tseries)*2-1) = tseries(2:size(tseries))+tseries(size(tseries))

ydat2(1:size(ydat)) = ydat
ydat2(size(ydat)+1 : size(ydat)*2-1) = ydat(2:size(ydat))

!!Select the bounds and part of the time series corresponding to the interval

      tpos_low = minloc(abs(tseries2 - tstart),1)
      if (((tseries2(tpos_low) - tstart)>0) .AND. (tpos_low>1)) then
        tpos_low = tpos_low-1
      end if
    
      tpos_high = minloc(abs(tseries2 - (tstart + tint)),1)
      if (((tseries2(tpos_high) - (tstart + tint))<0) .AND. (tpos_high<size(tseries2))) then
        tpos_high = tpos_high+1
      end if
    
      if (tpos_low == tpos_high) then
        tpos_high = tpos_low+1
      end if

tseries_select = tseries2(tpos_low:tpos_high)


!! Check if time series values are within the interval

countval=count(tseries_select>tstart .AND. tseries_select<(tstart+tint))


IF (countval<1) THEN
!! 1) if no values of the time series are within the desired time interval: interpolate by using the nearest values of the the time series

slope=((ydat2(tpos_high)-ydat2(tpos_low))/(tseries2(tpos_high)-tseries2(tpos_low)))  

yout=((tstart+(tint)/2)-tseries2(tpos_low))*slope + ydat2(tpos_low)



ELSE
  !! 2) Weighted average if values of the time series are within time interval

  tseries_within=tseries2((tpos_low+1):(tpos_high-1))
  ydat_within=ydat2((tpos_low+1):(tpos_high-1))

  slope1=((ydat2(tpos_low+1)-ydat2(tpos_low))/(tseries2(tpos_low+1)-tseries2(tpos_low)))  
  ymin=(tstart-tseries2(tpos_low))*slope1 + ydat2(tpos_low)

  slope2=((ydat2(tpos_high)-ydat2(tpos_high-1))/(tseries2(tpos_high)-tseries2(tpos_high-1))) 
  ymax=(tstart+tint-tseries2(tpos_high-1))*slope2 + ydat2(tpos_high-1)

  IF (ALLOCATED(weight)) THEN
    DEALLOCATE(weight)
    ALLOCATE(weight(1:(size(tseries_within,dim=1)+1)))
  ELSE
    ALLOCATE(weight(1:(size(tseries_within,dim=1)+1)))
  END IF

  IF (ALLOCATED(ydatdum)) THEN
    DEALLOCATE(ydatdum)
    ALLOCATE(ydatdum(1:(size(tseries_within,dim=1)+1)))
  ELSE
    ALLOCATE(ydatdum(1:(size(tseries_within,dim=1)+1)))
  END IF

  IF (ALLOCATED(ydatdum2)) THEN
    DEALLOCATE(ydatdum2)
    ALLOCATE(ydatdum2(1:(size(tseries_within,dim=1)+1)))
  ELSE
    ALLOCATE(ydatdum2(1:(size(tseries_within,dim=1)+1)))
  END IF

DO index = 1,(size(tseries_within,dim=1)+1)

  if (index==1) THEN

  
  weight(index)=(tseries_within(index)-tstart)/(tint)
  
  ydatdum(index)=(ydat_within(index)+ymin)/2 

  elseif (index==((size(tseries_within,dim=1)+1))) THEN
  
  
  weight(index)=((tstart+tint)-tseries_within(index-1))/(tint)  
  ydatdum(index)=(ymax+ydat_within(index-1))/2 
  else
  
  weight(index)=(tseries_within(index)-tseries_within(index-1))/(tint)  
  ydatdum(index)=(ydat_within(index)+ydat_within(index-1))/2 
  
  END IF

  ydatdum2(index)=weight(index)*ydatdum(index)  
  
END DO
yout=sum(ydatdum2)


    END IF      

!!    if (tstart>2) then
!!      WRITE(*,*) 'I am stopping'
!!      stop
!!     end if
    end subroutine interp3