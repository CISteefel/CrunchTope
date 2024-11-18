      subroutine ldbar3(ipc,ipcmax,nacdpr,narxmx,narxt,ndbmax,
     $ ntprmx,ntprt,xdbval,zdbval)
c
c     This subroutine loads data on the temperature grid from the
c     one-dimensional holding array (xdbval) into the ipc-th part
c     of the three-dimensional proper data array (represented here by
c     the dummy name zdbval).
c
c     This subroutine is called by:
c
c       EQPT/rdpar.f
c       EQPT/pcraq.f
c       EQPT/pcrsg.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       ipc    = the fixed dimension of the proper array
c       narxt  = array of numbers of coefficients in the temperature
c                  ranges
c       ntprt  = the number of temperature ranges on the standard
c                  temperature grid
c       xdbval = holding array (1-dimensional)
c
c     Principal output:
c
c       nacdpr = array containging the number of actual data points
c                  by range on the "log K" temperature grid; excludes
c       zdbval = dummy name for the proper data array (3-dimensional)
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer ipcmax,narxmx,ndbmax,ntprmx
c
      integer nacdpr(ntprmx),narxt(ntprmx)
c
      integer ipc,ntprt
c
      real*8 xdbval(ndbmax),zdbval(narxmx,ntprmx,ipcmax)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,k,n,ntpr
c
c-----------------------------------------------------------------------
c
c     Load the data from the xdbval array into the zdbval array.
c     Note that a point which comprises the boundary between two
c     ranges is a member of both ranges.
c
      ntpr = 1
      k = 0
      do n = 1,narxt(1)
        i = n
        zdbval(n,ntpr,ipc) = xdbval(i)
        if (xdbval(i) .lt. 9999999.) k = k + 1
      enddo
      nacdpr(ntpr) = k
c
      do ntpr = 2,ntprt
        k = 0
        zdbval(1,ntpr,ipc) = xdbval(i)
        if (xdbval(i) .lt. 9999999.) k = k + 1
        do n = 2,narxt(ntpr)
          i = i + 1
          zdbval(n,ntpr,ipc) = xdbval(i)
          if (xdbval(i) .lt. 9999999.) k = k + 1
        enddo
        nacdpr(ntpr) = k
      enddo
c
      end
