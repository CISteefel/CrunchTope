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

!!Select the pump series corresponding to the interval

      tpos_low=minloc(abs(tseries-tstart),1)
      if (((tseries(tpos_low)-tstart)>0) .AND. (tpos_low>1)) then
        tpos_low=tpos_low-1
      end if
    
      tpos_high=minloc(abs(tseries-(tstart+tint)),1)
      if (((tseries(tpos_high)-(tstart+tint))<0) .AND. (tpos_high<size(tseries))) then
        tpos_high=tpos_high+1
      end if
    
      if (tpos_low==tpos_high) then
        tpos_high=tpos_low+1
      end if

tseries_select=tseries(tpos_low:tpos_high)


!! Check if pump values are within the interval

countval=count(tseries_select>tstart .AND. tseries_select<(tstart+tint))


if (countval<1) THEN
!! 1) interpolate if they are outside

slope=((ydat(tpos_high)-ydat(tpos_low))/(tseries(tpos_high)-tseries(tpos_low)))  

yout=((tstart+(tint)/2)-tseries(tpos_low))*slope + ydat(tpos_low)



ELSE
  !! 2) Weighted average if pump time series within time interval

  tseries_within=tseries((tpos_low+1):(tpos_high-1))
  ydat_within=ydat((tpos_low+1):(tpos_high-1))

  slope1=((ydat(tpos_low+1)-ydat(tpos_low))/(tseries(tpos_low+1)-tseries(tpos_low)))  
  ymin=(tstart-tseries(tpos_low))*slope1 + ydat(tpos_low)

  slope2=((ydat(tpos_high)-ydat(tpos_high-1))/(tseries(tpos_high)-tseries(tpos_high-1))) 
  ymax=(tstart+tint-tseries(tpos_high-1))*slope2 + ydat(tpos_high-1)

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