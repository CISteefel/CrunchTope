SUBROUTINE interp2(tstart, tint, ydat, yout,info)
    ! Inputs: xData = a vector of the x-values of the data to be interpolated
    !         yData = a vector of the y-values of the data to be interpolated
    !         xVal  = a vector of the x-values where interpolation should be performed
    ! Output: yVal  = a vector of the resulting interpolated values
    
  USE crunchtype
  USE flow


      implicit none
    
      !  External variables and arrays
      real(DP),intent(in) :: tstart
      real(DP),intent(in) :: tint
      INTEGER(I4B),intent(in) :: info
      real(DP),intent(in) :: ydat(info)
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
      real(DP),dimension(:), ALLOCATABLE     :: tpump_select
      real(DP),dimension(:), ALLOCATABLE     :: tpump_within
      real(DP),dimension(:), ALLOCATABLE     :: ydat_within
      real(DP),dimension(:), ALLOCATABLE     :: weight
      real(DP),dimension(:), ALLOCATABLE     :: ydatdum
      real(DP),dimension(:), ALLOCATABLE     :: ydatdum2
      real(DP),dimension(:,:,:,:), ALLOCATABLE     :: checkdummy
      real(DP),dimension(:), ALLOCATABLE     :: checkdummy2
      REAL(DP)  :: dummy1
      REAL(DP)  :: dummy2
      REAL(DP)  :: dummy3
      INTEGER(I4B) :: dummycheck
      ! Possible checks on inputs could go here
      ! Things you may want to check:
      !   monotonically increasing xData
      !   size(xData) == size(yData)
      !   size(xVal) == size(yVal)

!!Select the pump series corresponding to the interval

      tpos_low=minloc(abs(tpump-tstart),1)
      if (((tpump(tpos_low)-tstart)>0) .AND. (tpos_low>1)) then
        tpos_low=tpos_low-1
      end if
    
      tpos_high=minloc(abs(tpump-(tstart+tint)),1)
      if (((tpump(tpos_high)-(tstart+tint))<0) .AND. (tpos_high<size(tpump))) then
        tpos_high=tpos_high+1
      end if
    
      if (tpos_low==tpos_high) then
        tpos_high=tpos_low+1
      end if

      checkdummy=qgt
      checkdummy2=tpump
tpump_select=tpump(tpos_low:tpos_high)


!! Check if pump values are within the interval

countval=count(tpump_select>tstart .AND. tpump_select<(tstart+tint))


if (countval<1) THEN
!! 1) interpolate if they are outside

slope=((ydat(tpos_high)-ydat(tpos_low))/(tpump(tpos_high)-tpump(tpos_low)))  

yout=((tstart+(tint)/2)-tpump(tpos_low))*slope + ydat(tpos_low)



ELSE
  !! 2) Weighted average if pump time series within time interval

  tpump_within=tpump((tpos_low+1):(tpos_high-1))
  ydat_within=ydat((tpos_low+1):(tpos_high-1))

  slope1=((ydat(tpos_low+1)-ydat(tpos_low))/(tpump(tpos_low+1)-tpump(tpos_low)))  
  ymin=(tstart-tpump(tpos_low))*slope1 + ydat(tpos_low)

  slope2=((ydat(tpos_high)-ydat(tpos_high-1))/(tpump(tpos_high)-tpump(tpos_high-1))) 
  ymax=(tstart+tint-tpump(tpos_high-1))*slope2 + ydat(tpos_high-1)

  IF (ALLOCATED(weight)) THEN
    DEALLOCATE(weight)
    ALLOCATE(weight(1:(size(tpump_within,dim=1)+1)))
  ELSE
    ALLOCATE(weight(1:(size(tpump_within,dim=1)+1)))
  END IF

  IF (ALLOCATED(ydatdum)) THEN
    DEALLOCATE(ydatdum)
    ALLOCATE(ydatdum(1:(size(tpump_within,dim=1)+1)))
  ELSE
    ALLOCATE(ydatdum(1:(size(tpump_within,dim=1)+1)))
  END IF

  IF (ALLOCATED(ydatdum2)) THEN
    DEALLOCATE(ydatdum2)
    ALLOCATE(ydatdum2(1:(size(tpump_within,dim=1)+1)))
  ELSE
    ALLOCATE(ydatdum2(1:(size(tpump_within,dim=1)+1)))
  END IF

  dummycheck=(size(tpump_within,dim=1)+1)

DO index = 1,(size(tpump_within,dim=1)+1)

  if (index==1) THEN

  
  weight(index)=(tpump_within(index)-tstart)/(tint)
  
  ydatdum(index)=(ydat_within(index)+ymin)/2 
  dummy1=index
  elseif (index==((size(tpump_within,dim=1)+1))) THEN
  
  
  weight(index)=((tstart+tint)-tpump_within(index-1))/(tint)  
  ydatdum(index)=(ymax+ydat_within(index-1))/2 
  dummy2=index
  else
  
  weight(index)=(tpump_within(index)-tpump_within(index-1))/(tint)  
  ydatdum(index)=(ydat_within(index)+ydat_within(index-1))/2 
  
  END IF

  ydatdum2(index)=weight(index)*ydatdum(index)  
  
END DO
dummy3=index   
yout=sum(ydatdum2)


    END IF      

    if (tstart>2) then
      WRITE(*,*) 'I am stopping'
      stop
      end if
    end subroutine interp2